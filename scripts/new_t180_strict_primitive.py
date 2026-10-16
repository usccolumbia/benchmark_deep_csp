"""Strict-primitive sensitivity on T180: reference and candidates are reduced to their primitive cells at 0.01 A and the
matcher's own primitive reduction is disabled (ltol .2, stol .3, angle 5, scale, no supercell). Writes results/new_t180/t180_strict_primitive_per_case.csv."""
import csv, json, sys, warnings
from multiprocessing import Pool
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import new_t180_topk_sg_all as T
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
M = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=False, scale=True, attempt_supercell=False)
REFS = {r["material_id"]: r["cif"] for r in json.loads((T.CAMP / "references.json").read_text())}
def work(t):
    a, mid = t
    gt = Structure.from_str(REFS[mid], fmt="cif").get_primitive_structure(tolerance=.01); rank = ""; n = 0
    for root in T.ALGOS[a]:
        d = root / a / mid
        if d.exists() and list(d.glob("*.cif")):
            files = sorted(d.glob("*.cif"), key=T.nat)[:20]; n = len(files)
            for i, p in enumerate(files, 1):
                try:
                    if M.fit(gt, Structure.from_file(p).get_primitive_structure(tolerance=.01)): rank = i; break
                except Exception: pass
            break
    return dict(algorithm=a, material_id=mid, n_cifs=n, strict_rank=rank)
if __name__ == "__main__":
    algos = [a for a in T.ALGOS if a != "cspml"]; mids = list(REFS)
    with Pool(8) as p: rows = p.map(work, [(a, m) for a in algos for m in mids], chunksize=8)
    out = Path(__file__).resolve().parent.parent / "results/new_t180/t180_strict_primitive_per_case.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("done", len(rows))
