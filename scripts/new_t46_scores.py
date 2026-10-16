"""T46 (v3.1) match rates for the paper's 16 methods: the 33 targets shared with T180 come from the T180 campaign scores,
the 13 T46-only targets from the new runs (t46x13_*), scored with the same matcher and space-group rule as new_t180_topk_sg_all.py."""
import json, re, warnings
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
ROOT = Path(__file__).resolve().parent.parent
N = Path("/private/tmp/claude-501/t46")
P, X = N / "t46x13_primary_20261009/outputs", N / "t46x13_extra4_20261009/outputs"
OUTS = {a: P for a in ["tcspv3","crystalflow","crystalite","omatg","uni3dar","uflowcsp","equicsp","mattergen","diffcsp","gemsdiff","crystalgrw","crystalformer"]}
OUTS.update({a: X for a in ["tgdmat","symmcd","cdvae","cond-cdvae"]})
_C = {}
nat = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]
def sg(s, sp):
    try: return SpacegroupAnalyzer(s, symprec=sp).get_space_group_number()
    except Exception: return None
def init():
    _C["refs"] = {r["material_id"]: r["cif"] for r in json.loads((N / "t46x13_primary_20261009/references.json").read_text())}
    _C["m"] = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False)
def work(t):
    a, mid = t; gt = Structure.from_str(_C["refs"][mid], fmt="cif"); d = OUTS[a] / a / mid
    files = sorted(d.glob("*.cif"), key=nat)[:20] if d.exists() else []
    row = dict(algorithm=a, material_id=mid, n_cifs=len(files), sm_rank="", top1_sm="", top1_sg01="", top1_sg_pred01="", top1_sg1="", top1_sg_pred1="", gt_sg01=sg(gt, .01), gt_sg1=sg(gt, .1))
    for i, p in enumerate(files, 1):
        try:
            s = Structure.from_file(p)
            if i == 1:
                row["top1_sg_pred01"], row["top1_sg_pred1"] = sg(s, .01), sg(s, .1)
                row["top1_sg01"] = row["top1_sg_pred01"] == row["gt_sg01"]; row["top1_sg1"] = row["top1_sg_pred1"] == row["gt_sg1"]
            if row["sm_rank"] == "" and _C["m"].fit(gt, s): row["sm_rank"] = i
        except Exception: pass
    if files: row["top1_sm"] = row["sm_rank"] == 1
    return row
if __name__ == "__main__":
    ids_new = [q["material_id"] for q in json.loads((N / "t46x13_primary_20261009/queries.json").read_text())]
    ids46 = list(pd.read_csv(ROOT / "testset_4benchmark_versions/v3_1_20261008_t180_replacement/T46/metadata.csv").material_id)
    shared = [m for m in ids46 if m not in ids_new]; print(len(ids46), len(ids_new), len(shared), flush=True)
    with Pool(8, initializer=init) as pool: rows = pool.map(work, [(a, m) for a in OUTS for m in ids_new], chunksize=2)
    new = pd.DataFrame(rows)
    old = pd.read_csv(ROOT / "algorithms_benchmark/T180/current_v3_1/metrics/t180_topk_sg_per_case_all.csv"); old = old[old.material_id.isin(shared) & old.algorithm.isin(OUTS)]
    d = pd.concat([old, new[old.columns]]); d.to_csv(ROOT / "results/new_t180/t46_topk_sg_per_case.csv", index=False); print(len(d), d.algorithm.nunique(), d.material_id.nunique())
