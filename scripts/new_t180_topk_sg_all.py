"""Recompute StructureMatcher ranks and top-1 space-group agreement for the new T180 campaign.

Reads exported CIFs (finished_predictions_*/<algo>/<mid>/*.cif), ranks them by natural
file-name order exactly like the campaign evaluator, and records the first matching
rank (<=20) under the primary matcher, plus top-1 space-group agreement at symprec 0.1
and 0.01 for the top-1 candidate. Used for the manuscript Results (Sec. 3).
"""
import csv, json, re, sys, warnings
from multiprocessing import Pool
from pathlib import Path
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

E = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_family_ablation_hyperion_20261006/results_from_clusters/enigma_T180")
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")
MAIN = E / "finished_predictions_20261007T042743Z"
GRW = E / "finished_predictions_20261007T143216Z_crystalgrw"
X = sorted((E / "extra4").glob("finished_predictions_*extra4*"))
ALGOS = {a: [MAIN] for a in ["tcspv3","crystalflow","diffcsp","equicsp","gemsdiff","mattergen","omatg","uni3dar","crystalite","uflowcsp","cspml","crystalformer"]}
ALGOS["crystalgrw"] = [GRW]
ALGOS["uflowcsp"] = [E / "finished_predictions_20261008_uflow_reorb"]  # re-relaxed with the common ORB-v3 settings (500 steps, fmax 0.01)
for a in ["cdvae","cond-cdvae","symmcd","tgdmat"]: ALGOS[a] = X
REPL = E / "finished_predictions_20261008_replacement"  # mp-3576 (replaces mp-1183076)
for a in ALGOS: ALGOS[a] = list(ALGOS[a]) + [REPL]
OUT = Path(__file__).resolve().parent.parent / "results/new_t180/t180_topk_sg_per_case_all.csv"
_C = {}

def nat(p):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]

def sg(s, sp):
    try:
        return SpacegroupAnalyzer(s, symprec=sp).get_space_group_number()
    except Exception:
        return None

def init():
    _C["refs"] = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
    _C["m"] = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False)

def work(t):
    a, mid = t
    gt = Structure.from_str(_C["refs"][mid], fmt="cif")
    files = []
    for root in ALGOS[a]:
        d = root / a / mid
        if d.exists() and list(d.glob("*.cif")):
            files = sorted(d.glob("*.cif"), key=nat)[:20]; break
    row = dict(algorithm=a, material_id=mid, n_cifs=len(files), sm_rank="", top1_sm="", top1_sg01="", top1_sg_pred01="",
               top1_sg1="", top1_sg_pred1="", gt_sg01=sg(gt, .01), gt_sg1=sg(gt, .1))
    for i, p in enumerate(files, 1):
        try:
            s = Structure.from_file(p)
            if i == 1:
                row["top1_sg_pred01"], row["top1_sg_pred1"] = sg(s, .01), sg(s, .1)
                row["top1_sg01"] = row["top1_sg_pred01"] == row["gt_sg01"]
                row["top1_sg1"] = row["top1_sg_pred1"] == row["gt_sg1"]
            if row["sm_rank"] == "" and _C["m"].fit(gt, s):
                row["sm_rank"] = i
        except Exception:
            pass
    if files:
        row["top1_sm"] = row["sm_rank"] == 1
    return row

if __name__ == "__main__":
    ids = [q["material_id"] for q in json.loads((CAMP / "queries.json").read_text())]
    tasks = [(a, m) for a in ALGOS for m in ids]
    fields = ["algorithm", "material_id", "n_cifs", "sm_rank", "top1_sm", "top1_sg01", "top1_sg_pred01", "top1_sg1", "top1_sg_pred1", "gt_sg01", "gt_sg1"]
    with Pool(8, initializer=init) as pool, open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for i, r in enumerate(pool.imap_unordered(work, tasks, chunksize=2), 1):
            w.writerow({k: r[k] for k in fields}); fh.flush()
            if i % 100 == 0: print(i, len(tasks), flush=True)
    print("done", flush=True)
