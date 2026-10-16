"""Prototype-level and framework-level removal arms for the T180 ablation.

Level 2 ("proto"):  remove train/val structures whose prototype label equals the label of any target of the family.
    label = anonymous formula | Pearson symbol | space group | Wyckoff letters per species (AFLOW-style, spglib symprec 0.1;
    species ordered by (reduced count, Wyckoff string)).
Level 3 ("frame"):  proto removal plus every train/val structure whose species-erased structure matches (StructureMatcher, ltol .2,
    stol .3, angle 5, primitive, scaled) the species-erased structure of any target of the family, across stoichiometries.
Original MP-20 membership kept (27,136 train / 9,047 val); removal from train and val, test never filtered.
"""
import csv, hashlib, json, math, sys, warnings
from collections import Counter, defaultdict
from functools import reduce
from multiprocessing import Pool
from pathlib import Path
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
HERE = Path(__file__).resolve().parent; MP = Path("/Users/jason/Downloads/project/2026_CSP_diagnosis_benchmark/mp_20")
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")
OUT = HERE / "data_prototype"; ARMS_FAMILIES = ["AB2", "ABC2"]  # splits are written for these; counts for all six
_M = None
def matcher():
    global _M
    if _M is None: _M = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False)
    return _M
def describe(arg):
    key, cif = arg
    try:
        s = Structure.from_str(cif, fmt="cif"); sga = SpacegroupAnalyzer(s, symprec=0.1); d = sga.get_symmetry_dataset()
        wy = defaultdict(list); seen = set()
        for i, site in enumerate(s):
            orb = d.equivalent_atoms[i]
            if orb in seen: continue
            seen.add(orb); wy[site.specie.symbol].append(d.wyckoffs[i])
        cnt = Counter(site.specie.symbol for site in s); g = reduce(math.gcd, cnt.values())
        parts = sorted((cnt[e] // g, "".join(sorted(wy[e]))) for e in cnt)
        conv = sga.get_conventional_standard_structure(); cs = sga.get_crystal_system()
        letter = {"triclinic": "a", "monoclinic": "m", "orthorhombic": "o", "tetragonal": "t", "trigonal": "h", "hexagonal": "h", "cubic": "c"}[cs]
        cent = sga.get_space_group_symbol()[0]; pearson = f"{letter}{cent}{len(conv)}"
        label = f"{''.join(f'{n}' for n,_ in parts)}|{pearson}|{d.number}|{'_'.join(w for _, w in parts)}"
        er = Structure(s.lattice, ["Fe"] * len(s), s.frac_coords); prim = er.get_primitive_structure(tolerance=0.1)
        fsg = SpacegroupAnalyzer(prim, symprec=0.1).get_space_group_number()
        return key, dict(label=label, red="".join(map(str, sorted(n for n, _ in parts))), fkey=f"{len(prim)}|{fsg}", natoms=len(s))
    except Exception as e:
        return key, dict(label=None, error=str(e)[:80])
def frame_fit(arg):
    tkey, tcif, cands = arg
    t = Structure.from_str(tcif, fmt="cif"); te = Structure(t.lattice, ["Fe"] * len(t), t.frac_coords).get_primitive_structure(tolerance=0.1)
    out = []
    for cid, ccif in cands:
        try:
            c = Structure.from_str(ccif, fmt="cif"); ce = Structure(c.lattice, ["Fe"] * len(c), c.frac_coords).get_primitive_structure(tolerance=0.1)
            if matcher().fit(te, ce): out.append(cid)
        except Exception: pass
    return tkey, out
def main():
    from pymatgen.core import Composition
    SIG = {"AB2": (1, 2), "ABC2": (1, 1, 2), "AB": (1, 1), "ABC2D6": (1, 1, 2, 6), "ABC2D4": (1, 1, 2, 4), "AB3": (1, 3)}
    sig = lambda f: tuple(sorted(int(v) // reduce(math.gcd, [int(x) for x in Composition(f).get_el_amt_dict().values()]) for v in Composition(f).get_el_amt_dict().values()))
    qs = json.loads((CAMP / "queries.json").read_text()); fam = {k: [q["material_id"] for q in qs if sig(q["primitive_cell_formula"]) == s] for k, s in SIG.items()}  # v3.1 targets
    rows = {}; fields = None
    for sp in ("train", "val", "test"):
        with open(MP / f"{sp}.csv", newline="") as fh:
            r = csv.DictReader(fh); fields = r.fieldnames; rows[sp] = list(r)
    refs = {x["material_id"]: x["cif"] for x in json.loads((CAMP / "references.json").read_text())}
    items = [(f"{sp}:{r['material_id']}", r["cif"]) for sp in ("train", "val") for r in rows[sp]] + [(f"target:{m}", c) for m, c in refs.items()]
    with Pool(8) as p: info = dict(p.map(describe, items, chunksize=64))
    bad = [k for k, v in info.items() if v["label"] is None]; print("labelled", len(info) - len(bad), "failed", len(bad), flush=True)
    lab = {k: v["label"] for k, v in info.items() if v["label"]}
    json.dump({k: v for k, v in info.items()}, open(OUT / "labels.json", "w")) if OUT.exists() else (OUT.mkdir(), json.dump(info, open(OUT / "labels.json", "w")))
    trainval = {k: v for k, v in info.items() if k.startswith(("train:", "val:")) and v["label"]}
    byproto = defaultdict(set); byframe = defaultdict(list)
    for k, v in trainval.items(): byproto[v["label"]].add(k); byframe[v["fkey"]].append(k)
    cifs = {f"{sp}:{r['material_id']}": r["cif"] for sp in ("train", "val") for r in rows[sp]}
    tasks = []
    for m, c in refs.items():
        fk = info[f"target:{m}"].get("fkey")
        if fk: tasks.append((m, c, [(k, cifs[k]) for k in byframe.get(fk, [])]))
    print("framework candidate pairs", sum(len(t[2]) for t in tasks), flush=True)
    with Pool(8) as p: fits = dict(p.map(frame_fit, tasks, chunksize=2))
    report = {}
    for f, targets in fam.items():
        proto = set().union(*[byproto.get(info[f"target:{m}"]["label"], set()) for m in targets if info.get(f"target:{m}", {}).get("label")]) if targets else set()
        frame = set(proto).union(*[set(fits.get(m, [])) for m in targets])
        cross = {k for k in frame if info[k]["red"] != info[f"target:{targets[0]}"]["red"] and all(info[k]["red"] != info[f"target:{m}"]["red"] for m in targets)}
        report[f] = dict(targets=len(targets), distinct_target_prototypes=len({info[f"target:{m}"]["label"] for m in targets}),
                         proto=dict(train=sum(k.startswith("train:") for k in proto), val=sum(k.startswith("val:") for k in proto)),
                         frame=dict(train=sum(k.startswith("train:") for k in frame), val=sum(k.startswith("val:") for k in frame)),
                         frame_cross_stoichiometry=len(cross),
                         targets_with_zero_proto=sum(1 for m in targets if not byproto.get(info[f"target:{m}"]["label"])),
                         targets_with_zero_frame=sum(1 for m in targets if not fits.get(m) and not byproto.get(info[f"target:{m}"]["label"])))
        if f in ARMS_FAMILIES:
            for level, rem in (("proto", proto), ("frame", frame)):
                arm = f"{f}_{level}"
                for sp in ("train", "val", "test"):
                    keep = [r for r in rows[sp] if sp == "test" or f"{sp}:{r['material_id']}" not in rem]
                    d = OUT / "splits" / arm; d.mkdir(parents=True, exist_ok=True)
                    with open(d / f"{sp}.csv", "w", newline="") as fh:
                        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(keep)
                    report[f].setdefault("arms", {}).setdefault(arm, {})[sp] = dict(total=len(rows[sp]), remaining=len(keep), removed=len(rows[sp]) - len(keep), sha256=hashlib.sha256((d / f"{sp}.csv").read_bytes()).hexdigest())
    json.dump(fits, open(OUT / "framework_fits.json", "w"))
    # nested arms: whole family removed AND every same-framework entry of the family's targets removed (family + prototype remains)
    for f in ("AB",):
        kept = {sp: {r["material_id"] for r in csv.DictReader(open(HERE / "data_selected/splits" / f / f"{sp}.csv", newline=""))} for sp in ("train", "val")}
        frame_f = set().union(*[set(fits.get(m, [])) for m in fam[f]]) | set().union(*[byproto.get(info[f"target:{m}"]["label"], set()) for m in fam[f] if info[f"target:{m}"].get("label")])
        arm = f"{f}_famframe"; report[f].setdefault("arms", {})[arm] = {}
        for sp in ("train", "val", "test"):
            keep = [r for r in rows[sp] if sp == "test" or (r["material_id"] in kept[sp] and f"{sp}:{r['material_id']}" not in frame_f)]
            d = OUT / "splits" / arm; d.mkdir(parents=True, exist_ok=True)
            with open(d / f"{sp}.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(keep)
            report[f]["arms"][arm][sp] = dict(total=len(rows[sp]), remaining=len(keep), removed=len(rows[sp]) - len(keep), sha256=hashlib.sha256((d / f"{sp}.csv").read_bytes()).hexdigest())
    (OUT / "removal_report.json").write_text(json.dumps(report, indent=1)); print(json.dumps({f: {k: v for k, v in r.items() if k != "arms"} for f, r in report.items()}, indent=1))
    # label-collision check: sample same-label train pairs, fraction that fit_anonymous
    import random; random.seed(0); m = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False)
    groups = [v for v in byproto.values() if len(v) >= 2]; ok = n = 0
    for _ in range(300):
        a, b = random.sample(sorted(random.choice(groups)), 2); n += 1
        try: ok += bool(m.fit_anonymous(Structure.from_str(cifs[a], fmt="cif"), Structure.from_str(cifs[b], fmt="cif")))
        except Exception: pass
    print(f"same-label pairs that fit_anonymous: {ok}/{n}")
if __name__ == "__main__": main()
