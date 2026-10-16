"""Distance metrics (CSPBenchMetrics definitions) between the reference and each method's top-1
candidate for the two T180 case-study targets. Energy distance uses MatterSim (frozen benchmark MLIP)
on relaxed structures; Chamfer, ICP ('superpose') RMSD, CrystalNN-fingerprint and OFM distances are
taken verbatim from github.com/usccolumbia/CSPBenchMetrics (code/distance_single.py)."""
import json, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
STAGE = sys.argv[1] if len(sys.argv) > 1 else 'geometry'
if STAGE == 'geometry':
    import open3d as o3d, point_cloud_utils as pcu  # open3d first: the reverse order segfaults
from pymatgen.core import Structure
if STAGE == 'geometry':
    from matminer.featurizers.site import CrystalNNFingerprint
    from matminer.featurizers.structure import SiteStatsFingerprint, OrbitalFieldMatrix
ROOT = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(ROOT))
if STAGE == 'energy':
    from csp_bench.mlip import MatterSimMLIP
E = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_family_ablation_hyperion_20261006/results_from_clusters/enigma_T180")
MAIN = E / "finished_predictions_20261007T042743Z"; GRW = E / "finished_predictions_20261007T143216Z_crystalgrw"
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")
if STAGE == 'energy':
    sys.path.insert(0, str(ROOT / "scripts"))
    import new_t180_topk_sg_all as T  # same export roots as the top-k results
ALGOS = ["tcspv3","crystalflow","crystalite","omatg","uni3dar","uflowcsp","equicsp","mattergen","diffcsp","cond-cdvae","gemsdiff","tgdmat","symmcd","crystalgrw","crystalformer","cdvae"]
TARGETS = ["mp-1105290", "mp-261"]
if STAGE == 'geometry':
  ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('ops', distance_cutoffs=None, x_diff_weight=0), stats=('mean','std_dev','minimum','maximum'))
  ofm = OrbitalFieldMatrix()
nat = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]

def finger(s1, s2):
    try: return float(np.linalg.norm(np.array(ssf.featurize(s1)) - np.array(ssf.featurize(s2))))
    except Exception: return 9999.0

def calc_rmsd(source, target, matching):
    matching = np.array(matching)
    pa = [source.points[i] for i in matching[:, 0]]; pb = [target.points[i] for i in matching[:, 1]]
    if matching.shape[0] < len(source.points):
        for i in set(range(len(source.points))) - set(matching[:, 0]): pa.append(source.points[i])
        for i in set(range(len(target.points))) - set(matching[:, 1]): pb.append(target.points[i])
    source.points = o3d.utility.Vector3dVector(pa); target.points = o3d.utility.Vector3dVector(pb)
    return float(np.sqrt(np.mean(np.sum(np.square(np.asarray(source.points) - np.asarray(target.points)), axis=1))))

def icp_rmsd(s1, s2):
    src, tgt = o3d.geometry.PointCloud(), o3d.geometry.PointCloud()
    src.points = o3d.utility.Vector3dVector(np.array([x.coords for x in s1.sites])); tgt.points = o3d.utility.Vector3dVector(np.array([y.coords for y in s2.sites]))
    reg = o3d.pipelines.registration.registration_icp(src, tgt, 5, np.identity(4), o3d.pipelines.registration.TransformationEstimationPointToPoint(), o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=1000))
    src.transform(reg.transformation); rmsd = reg.inlier_rmse
    if reg.fitness < 1:
        try: rmsd = calc_rmsd(src, tgt, np.asarray(reg.correspondence_set))
        except Exception: rmsd = 9999.0
    return float(rmsd)

def chamfer(s1, s2):
    return float(pcu.chamfer_distance(np.array([x.coords for x in s1.sites]), np.array([y.coords for y in s2.sites])))

def ofm_dist(s1, s2):
    return float(np.linalg.norm(ofm.featurize(s1) - ofm.featurize(s2)))

if __name__ == "__main__":
    OUTE = ROOT / "results/new_t180/case_study_energy.csv"
    if STAGE == "energy":   # run in the env that has MatterSim; records the chosen top-1 file for each pair
        refs = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
        mlip = MatterSimMLIP(model=str(ROOT / "mattersim-v1.0.0-5M.pth"), device="cpu"); rows = []
        for mid in TARGETS:
            gt = Structure.from_str(refs[mid], fmt="cif").get_primitive_structure(tolerance=0.01)
            egt = mlip.relax(gt, fmax=0.02, steps=500, cell_relax=True, constrain_symmetry=False).energy_per_atom
            for a in ALGOS:
                fs = []
                for root in T.ALGOS[a]:
                    d = root / a / mid
                    if d.exists() and list(d.glob("*.cif")): fs = sorted(d.glob("*.cif"), key=nat); break
                if not fs: rows.append(dict(material_id=mid, algorithm=a, error="no candidate")); continue
                st = Structure.from_file(fs[0]); r = dict(material_id=mid, algorithm=a, path=str(fs[0]), file=fs[0].name, n_sites=len(st))
                try: r["energy_distance"] = abs(mlip.relax(st, fmax=0.02, steps=500, cell_relax=True, constrain_symmetry=False).energy_per_atom - egt)
                except Exception as e: r["energy_distance"] = None; r["error"] = f"energy: {e}"
                rows.append(r); print(r, flush=True)
        pd.DataFrame(rows).to_csv(OUTE, index=False); print("energy done")
    else:                   # geometry distances for the recorded candidates (no MatterSim needed)
        refs = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
        e = pd.read_csv(OUTE); rows = []
        for _, x in e.iterrows():
            r = x.to_dict()
            if isinstance(r.get("path"), str):
                gt = Structure.from_str(refs[r["material_id"]], fmt="cif").get_primitive_structure(tolerance=0.01); s_ = Structure.from_file(r["path"])
                r["chamfer"] = chamfer(gt, s_); r["superpose_rmsd"] = icp_rmsd(gt, s_); r["fingerprint"] = finger(gt, s_)
                try: r["ofm"] = ofm_dist(gt, s_)
                except Exception: r["ofm"] = None
            rows.append(r); print(r, flush=True)
        pd.DataFrame(rows).drop(columns=["path"], errors="ignore").to_csv(ROOT / "results/new_t180/case_study_distances.csv", index=False); print("geometry done")
