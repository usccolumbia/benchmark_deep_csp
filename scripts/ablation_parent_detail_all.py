"""For every AB2 target: framework-matching parents (more species, same procedure as ablation_degeneracy_test_relaxed.py) with space group, Wyckoff letters and cell size,
to see whether parents share the target's space group and Wyckoff positions."""
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
def wy(s):
    d = SpacegroupAnalyzer(s, symprec=.01).get_symmetry_dataset()
    g = lambda k: d[k] if isinstance(d, dict) else getattr(d, k)
    return int(g("number")), "".join(sorted(g("wyckoffs")))
ARM = {}
def work(args):
    mid, arm, ref_cif, cands = args
    tp = Structure.from_str(ref_cif, fmt="cif").get_primitive_structure(tolerance=.01); nsp = len(tp.composition.elements); tsg, tw = wy(tp); out = []
    for c in cands:
        if c["nsp"] <= nsp: continue
        try:
            s = Structure.from_str(c["cif"], fmt="cif")
            if D.M.fit(tp, s):
                sg, w = wy(s); out.append(dict(material_id=mid, arm=arm, target_formula=tp.composition.reduced_formula, target_sg=tsg, target_wyckoff=tw, target_nsites=len(tp), parent_formula=s.composition.reduced_formula, parent_sg=sg, parent_wyckoff=w, parent_nsites=len(s), same_sg=sg == tsg, same_sg_and_wyckoff=(sg == tsg and w == tw), same_nsites=len(s) == len(tp)))
        except Exception: pass
    return out
if __name__ == "__main__":
    refs = {r["material_id"]: r["cif"] for r in json.loads((D.CAMP / "references.json").read_text())}
    df = pd.concat([pd.read_csv(D.ROOT / "mp_20" / f) for f in ("train.csv", "val.csv")])
    with Pool(8) as p: prepped = [x for x in p.map(D.prep, df.cif.tolist(), chunksize=200) if x]
    jobs = []
    for arm in D.FAM:
      for mid in D.FAM[arm]:
        tp = Structure.from_str(refs[mid], fmt="cif").get_primitive_structure(tolerance=.01)
        cands = [c for c in prepped if c["sig"] != D.SIG[arm] and (c["nsites"] % len(tp) == 0 or len(tp) % c["nsites"] == 0) and max(c["nsites"], len(tp)) / min(c["nsites"], len(tp)) <= 4]
        jobs.append((mid, arm, refs[mid], cands))
    with Pool(8) as p: res = p.map(work, jobs, chunksize=1)
    pd.DataFrame([r for x in res for r in x]).to_csv(D.ROOT / "results/new_t180/ablation_parent_detail_all.csv", index=False); print("done")
