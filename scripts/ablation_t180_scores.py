"""Score the T180 family-ablation arms (CrystalFlow, EquiCSP) against the reused baselines.

For each (method, arm, target in the removed family) the 20 exported candidates of the ablated model and of the baseline
(new T180 campaign) are compared with the primitive reference with the campaign matcher (ltol .2, stol .3, angle 5,
primitive_cell, scale, no supercell).  Fit = some candidate matches; consensus = some candidate matches AND has the
reference space group (spglib via pymatgen, symprec 0.1); degraded = baseline fit and ablated no fit.
Usage: python scripts/ablation_t180_scores.py <ablation_outputs_root>   (root contains <algo>/<arm>/<mid>/*.cif)"""
import json, re, sys, warnings
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import new_t180_topk_sg_all as T  # baseline export roots (T.ALGOS)

OUT = Path(sys.argv[1]); CAMP = T.CAMP
refs = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
fam = json.loads((ROOT / "results/new_t180/ablation_family_targets.json").read_text())
m = StructureMatcher(ltol=.2, stol=.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False)
nat = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]
def sg(s, sp=0.1):
    try: return SpacegroupAnalyzer(s, symprec=sp).get_space_group_number()
    except Exception: return None
def score(files, gt, gsg):
    fit = cons = False; rank = None
    for i, p in enumerate(files[:20], 1):
        try: s = Structure.from_file(p)
        except Exception: continue
        if m.fit(gt, s):
            if not fit: fit = True; rank = i
            if sg(s) == gsg: cons = True
    return fit, cons, rank
def base_files(algo, mid):
    for root in T.ALGOS[algo]:
        d = root / algo / mid
        if d.exists() and list(d.glob("*.cif")): return sorted(d.glob("*.cif"), key=nat)
    return []
rows = []
for algo in ("crystalflow", "equicsp"):
    for arm_dir in sorted((OUT / algo).glob("*")):
        arm = arm_dir.name
        for mid in fam.get(arm, []):
            gt = Structure.from_str(refs[mid], fmt="cif").get_primitive_structure(tolerance=0.01); gsg = sg(gt)
            af = sorted((arm_dir / mid).glob("*.cif"), key=nat); bf = base_files(algo, mid)
            a_fit, a_cons, a_rank = score(af, gt, gsg); b_fit, b_cons, b_rank = score(bf, gt, gsg)
            rows.append(dict(algo=algo, arm=arm, material_id=mid, n_abl=len(af), n_base=len(bf), base_fit=b_fit, base_cons=b_cons, base_rank=b_rank,
                             abl_fit=a_fit, abl_cons=a_cons, abl_rank=a_rank, degraded=(b_fit and not a_fit), degraded_cons=(b_cons and not a_cons)))
            print(rows[-1], flush=True)
d = pd.DataFrame(rows); d.to_csv(ROOT / "results/new_t180/ablation_per_target.csv", index=False)
g = d.groupby(["algo", "arm"]).agg(n=("material_id", "size"), base_fit=("base_fit", "sum"), base_cons=("base_cons", "sum"), abl_fit=("abl_fit", "sum"),
                                   abl_cons=("abl_cons", "sum"), degraded=("degraded", "sum"), degraded_cons=("degraded_cons", "sum")).reset_index()
g.to_csv(ROOT / "results/new_t180/ablation_summary.csv", index=False); print(g.to_string(index=False))
