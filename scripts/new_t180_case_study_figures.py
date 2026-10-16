"""Structure figures for the two T180 case-study targets (reference + top-1 candidates)."""
import json, re, warnings
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from ase.visualize.plot import plot_atoms
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
E = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_family_ablation_hyperion_20261006/results_from_clusters/enigma_T180")
MAIN = E / "finished_predictions_20261007T042743Z"
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")
FIG = ROOT / "P623_Deep_CSP_benchmark/figures"
NAME = {"tcspv3": "TCSPv3", "crystalflow": "CrystalFlow", "omatg": "OMatG", "uni3dar": "Uni-3DAR", "equicsp": "EquiCSP", "mattergen": "MatterGen", "diffcsp": "DiffCSP"}
PLAN = {"mp-1105290": ("Co3Sb4O6F6", ["tcspv3", "crystalflow", "omatg", "equicsp", "mattergen"], "-80x,25y,0z"),
        "mp-261": ("Ba6Al10", ["tcspv3", "crystalflow", "omatg", "uni3dar", "equicsp"], "-80x,20y,0z")}
nat = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]
refs = {r["material_id"]: r["cif"] for r in json.loads((CAMP / "references.json").read_text())}
res = pd.read_csv(ROOT / "results/new_t180/t180_topk_sg_per_case.csv")
for mid, (stem, algos, rot) in PLAN.items():
    panels = [("Reference", Structure.from_str(refs[mid], fmt="cif").get_primitive_structure(tolerance=0.01), None)]
    for a in algos:
        s = Structure.from_file(sorted((MAIN / a / mid).glob("*.cif"), key=nat)[0])
        r = res[(res.algorithm == a) & (res.material_id == mid)].iloc[0]
        panels.append((NAME[a], s, (bool(r.top1_sm), int(r.top1_sg_pred1))))
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 6.6))
    for ax, (lab, s, info), letter in zip(axes.ravel(), panels, "abcdef"):
        at = AseAtomsAdaptor.get_atoms(s)
        plot_atoms(at, ax, rotation=rot, radii=0.45, show_unit_cell=2)
        ax.set_axis_off()
        sub = "" if info is None else f"\nspace group {info[1]}, " + ("match" if info[0] else "no match")
        ax.set_title(f"({letter}) {lab}{sub}", fontsize=10)
    elems = sorted({e.symbol for e in panels[0][1].composition.elements}, key=lambda x: atomic_numbers[x])
    handles = [plt.Line2D([], [], marker="o", ls="", ms=9, color=jmol_colors[atomic_numbers[e]], markeredgecolor="k", label=e) for e in elems]
    fig.legend(handles=handles, loc="lower center", ncol=len(elems), frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIG / f"case_{stem}.png", dpi=200); fig.savefig(FIG / f"case_{stem}.pdf"); plt.close(fig)
    print("saved", stem)
