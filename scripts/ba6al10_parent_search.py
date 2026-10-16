"""Do retained MP-20 structures with more species than Ba6Al10 (mp-261, P6_3/mmc) share its framework (species ignored, supercells <=4x)?"""
import json, sys, warnings
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ablation_degeneracy_test_relaxed as D
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
def info(s):
    d = SpacegroupAnalyzer(s, symprec=.01).get_symmetry_dataset(); g = lambda k: d[k] if isinstance(d, dict) else getattr(d, k)
    return int(g("number")), "".join(sorted(g("wyckoffs")))
def work(args):
    tp, cands = args; tsg, tw = info(tp); out = []
    for c in cands:
        try:
            s = Structure.from_str(c["cif"], fmt="cif")
            if D.M.fit(tp, s):
                sg, w = info(s); out.append(dict(parent_formula=s.composition.reduced_formula, nsp=len(s.composition.elements), parent_sg=sg, parent_wyckoff=w, parent_nsites=len(s), same_sg=sg == tsg, exact=(sg == tsg and w == tw), target_sg=tsg, target_wyckoff=tw, target_nsites=len(tp)))
        except Exception: pass
    return out
if __name__ == "__main__":
    refs = {r["material_id"]: r["cif"] for r in json.loads((D.CAMP / "references.json").read_text())}
    tp = Structure.from_str(refs["mp-261"], fmt="cif").get_primitive_structure(tolerance=.01); nsp = len(tp.composition.elements)
    df = pd.concat([pd.read_csv(D.ROOT / "mp_20" / f) for f in ("train.csv", "val.csv")])
    with Pool(8) as p: prepped = [x for x in p.map(D.prep, df.cif.tolist(), chunksize=200) if x]
    cands = [c for c in prepped if c["nsp"] != nsp or True]
    cands = [c for c in cands if (c["nsites"] % len(tp) == 0 or len(tp) % c["nsites"] == 0) and max(c["nsites"], len(tp)) / min(c["nsites"], len(tp)) <= 4]
    chunks = [cands[i::8] for i in range(8)]
    with Pool(8) as p: res = [r for x in p.map(work, [(tp, ch) for ch in chunks]) for r in x]
    d = pd.DataFrame(res); d.to_csv(D.ROOT / "results/new_t180/ba6al10_framework_matches.csv", index=False)
    print("target", tp.composition.reduced_formula, len(tp), "species", nsp, "| candidates", len(cands), "| framework matches", len(d))
    if len(d):
        print(d.groupby("nsp").size().to_dict()); print("same SG", int(d.same_sg.sum()), "exact (SG+Wyckoff)", int(d.exact.sum()))
        print(d[d.nsp > nsp].drop_duplicates("parent_formula").head(15)[["parent_formula", "nsp", "parent_sg", "parent_nsites"]].to_string(index=False))
        print(d[d.nsp <= nsp].drop_duplicates("parent_formula").head(10)[["parent_formula", "nsp", "parent_sg", "parent_nsites"]].to_string(index=False))
