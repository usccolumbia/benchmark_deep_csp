"""Collect the CIFs behind the case-study structure figures into P623_Deep_CSP_benchmark/figures/case_study/ (reference + each method's rank-1 candidate)."""
import json, re, shutil, hashlib, warnings
from pathlib import Path
import pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
E = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_family_ablation_hyperion_20261006/results_from_clusters/enigma_T180"); MAIN = E / "finished_predictions_20261007T042743Z"
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31"); OUT = Path("P623_Deep_CSP_benchmark/figures/case_study"); OUT.mkdir(exist_ok=True)
PLAN = {"mp-1105290": ("Co3Sb4O6F6", ["tcspv3", "crystalflow", "omatg", "equicsp", "mattergen"]), "mp-261": ("Ba6Al10", ["tcspv3", "crystalflow", "omatg", "uni3dar", "equicsp"])}
nat = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]
refs = {r["material_id"]: r["cif"] for r in json.load(open(CAMP / "references.json"))}
res = pd.read_csv("results/new_t180/t180_topk_sg_per_case.csv"); rows = []
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
for mid, (stem, algos) in PLAN.items():
    d = OUT / stem; d.mkdir(exist_ok=True)
    ref = d / f"reference_{stem}_{mid}.cif"; ref.write_text(refs[mid])
    rows.append({"target": mid, "panel": "(a)", "role": "reference", "method": "Materials Project", "candidate_id": mid, "rank": "", "file": str(ref.relative_to(OUT)),
                 "source": "T180 v3.1 references.json", "rank_basis": "", "top1_match": "", "pred_sg_symprec0.1": "", "sha256": sha(ref)})
    for letter, a in zip("bcdef", algos):
        src = sorted((MAIN / a / mid).glob("*.cif"), key=nat)[0]; cid = src.stem
        formula = Structure.from_file(src).composition.formula.replace(" ", "")
        dst = d / f"rank01_{formula}_{a}_{cid}.cif"; shutil.copy2(src, dst)
        r = res[(res.algorithm == a) & (res.material_id == mid)].iloc[0]
        basis = "TCSPv3 native chemical score (lower is better)" if a == "tcspv3" else "ORB-v3 relaxation, energy per atom ascending"
        rows.append({"target": mid, "panel": f"({letter})", "role": "top-1 candidate", "method": a, "candidate_id": cid, "rank": 1, "file": str(dst.relative_to(OUT)),
                     "source": str(src), "rank_basis": basis, "top1_match": bool(r.top1_sm), "pred_sg_symprec0.1": int(r.top1_sg_pred1), "sha256": sha(dst)})
pd.DataFrame(rows).to_csv(OUT / "manifest.csv", index=False)
(OUT / "README.md").write_text("CIFs behind the case-study structure figures (case_Co3Sb4O6F6, case_Ba6Al10). Panel (a) is the Materials Project reference; panels (b)-(f) are each method's rank-1 exported candidate from the T180 v3.1 campaign. File names carry rank, formula, method and candidate ID. Use manifest.csv for paths: it records the ranking basis, source path, match status and SHA256. The figures were drawn by scripts/new_t180_case_study_figures.py, which reduces the reference to its 0.01 A primitive cell.\n")
print(pd.DataFrame(rows)[["target", "panel", "method", "file", "top1_match"]].to_string(index=False))
