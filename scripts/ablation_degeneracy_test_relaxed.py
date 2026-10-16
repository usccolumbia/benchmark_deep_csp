"""Test the degeneracy explanation of the prototype ablation.
For every T180 target of an ablated family, find structures left in the MP20 train+val data (everything except the target's own reduced-stoichiometry family)
whose framework (sites only, species ignored; StructureMatcher with FrameworkComparator) matches the target's primitive cell, and record how many distinct
species those structures have relative to the target. 'Parents' (more species than the target) can be reduced to the target by merging species; 'children'
(fewer species) would have to be split."""
import json, sys, warnings
from collections import Counter
from math import gcd
from functools import reduce
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
ROOT = Path(__file__).resolve().parent.parent
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")
FAM = json.loads((ROOT / "results/new_t180/ablation_family_targets.json").read_text())
SIG = {"AB2": (1, 2), "ABC2": (1, 1, 2), "AB": (1, 1), "ABC2D6": (1, 1, 2, 6)}
def sig(comp):
    v = [int(round(x)) for x in Composition(comp).get_el_amt_dict().values()]; g = reduce(gcd, v); return tuple(sorted(x // g for x in v))
M = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=True, comparator=FrameworkComparator())
def prep(cif):
    try:
        s = Structure.from_str(cif, fmt="cif"); p = s.get_primitive_structure(tolerance=.01)
        return dict(sig=sig(s.composition), nsp=len(s.composition.elements), nsites=len(p), sg=SpacegroupAnalyzer(p, symprec=.01).get_space_group_number(), cif=p.to(fmt="cif"))
    except Exception: return None
def work(args):
    mid, ref_cif, cands = args
    t = Structure.from_str(ref_cif, fmt="cif"); tp = t.get_primitive_structure(tolerance=.01); nsp = len(tp.composition.elements)
    hits = Counter()
    for c in cands:
        try:
            if M.fit(tp, Structure.from_str(c["cif"], fmt="cif")): hits["more" if c["nsp"] > nsp else "equal" if c["nsp"] == nsp else "fewer"] += 1
        except Exception: pass
    return mid, dict(n_cands=len(cands), **{k: hits.get(k, 0) for k in ("more", "equal", "fewer")}, nsp=nsp)
if __name__ == "__main__":
    refs = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
    df = pd.concat([pd.read_csv(ROOT / "mp_20" / f) for f in ("train.csv", "val.csv")])
    with Pool(8) as p: prepped = [x for x in p.map(prep, df.cif.tolist(), chunksize=200) if x]
    print("prepped", len(prepped), flush=True)
    jobs = []; fam_of = {}
    for arm, mids in FAM.items():
        for mid in mids:
            fam_of[mid] = arm; tp = Structure.from_str(refs[mid], fmt="cif").get_primitive_structure(tolerance=.01)
            sg = SpacegroupAnalyzer(tp, symprec=.01).get_space_group_number()
            cands = [c for c in prepped if c["sig"] != SIG[arm] and (c["nsites"] % len(tp) == 0 or len(tp) % c["nsites"] == 0) and max(c["nsites"], len(tp)) / min(c["nsites"], len(tp)) <= 4]
            jobs.append((mid, refs[mid], cands))
    with Pool(8) as p: res = dict(p.map(work, jobs, chunksize=1))
    rows = [dict(material_id=m, arm=fam_of[m], **r) for m, r in res.items()]
    pd.DataFrame(rows).to_csv(ROOT / "results/new_t180/ablation_degeneracy_relaxed.csv", index=False); print("done", len(rows))
