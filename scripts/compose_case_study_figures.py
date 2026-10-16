"""Compose the case-study figures (Figs. 7-8) from the VESTA screenshots in figures/case_study/primitive/<target>/<method>.png.
Missing TCSPv3 screenshot for Ba6Al10 is replaced by a text panel until primitive/Ba6Al10/tcspv3.png exists."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
FIG = Path("P623_Deep_CSP_benchmark/figures"); CS = FIG / "case_study"
man = pd.read_csv(CS / "manifest.csv")
NAME = {"tcspv3": "TCSPv3", "crystalflow": "CrystalFlow", "omatg": "OMatG", "uni3dar": "Uni-3DAR", "equicsp": "EquiCSP", "mattergen": "MatterGen"}
PLAN = {"mp-1105290": ("Co3Sb4O6F6", ["tcspv3", "crystalflow", "omatg", "equicsp", "mattergen"]), "mp-261": ("Ba6Al10", ["tcspv3", "crystalflow", "omatg", "uni3dar", "equicsp"])}
def load(p):
    im = Image.open(p).convert("RGBA"); bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im); rgb = np.array(bg.convert("RGB"))
    mask = (rgb < 245).any(axis=2); ys, xs = np.where(mask); pad = 12
    return rgb[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
def pad_to(a, H, W):
    out = np.full((H, W, 3), 255, dtype=np.uint8); y = (H - a.shape[0]) // 2; x = (W - a.shape[1]) // 2; out[y:y + a.shape[0], x:x + a.shape[1]] = a; return out
RULES = {  # element -> (name, rgb predicate on the screenshot pixels), used only to sample legend colours
    "Co": lambda r, g, b: (b > 110) & (r < 60) & (g < 60), "Sb": lambda r, g, b: (r > 170) & (g > 110) & (g < 160) & (b < 110), "O": lambda r, g, b: (r > 200) & (g < 40) & (b < 40),
    "F": lambda r, g, b: (r > 150) & (r < 200) & (g > 160) & (g < 205) & (b > 200), "Ba": lambda r, g, b: (g > 170) & (r < 70) & (b < 70), "Al": lambda r, g, b: (r > 100) & (r < 150) & (g > 135) & (g < 175) & (b > 165) & (b < 200)}
def colour(a, el):
    r, g, b = [a[..., i].astype(int) for i in range(3)]; m = RULES[el](r, g, b)
    return tuple(np.median(a[m], axis=0) / 255) if m.sum() > 50 else (0.5, 0.5, 0.5)
def find(d, name):
    for p in d.iterdir():
        if p.suffix.lower() == ".png" and p.stem.lower() == name: return p
for mid, (stem, algos) in PLAN.items():
    d = CS / "primitive" / stem; panels = [("Reference", find(d, "reference"), None)]
    for a in algos:
        r = man[(man.target == mid) & (man.method == a)].iloc[0]
        info = f"\nspace group {int(r['pred_sg_symprec0.1'])}, " + ("match" if bool(r.top1_match) else "no match")
        panels.append((NAME[a], find(d, a), info, a, r))
    ims = {lab: load(path) for lab, path, *_ in panels if path is not None}
    H = max(a.shape[0] for a in ims.values()); W = max(a.shape[1] for a in ims.values()); ref_arr = ims["Reference"]
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 7.6)); axes = axes.ravel(); shown, miss = [], []
    for k, (lab, path, info, *rest) in enumerate(panels):
        ax = axes[k]; ax.set_axis_off(); ax.set_title(f"({'abcdef'[k]}) {lab}{info or ''}", fontsize=10)
        if path is not None: ax.imshow(pad_to(ims[lab], H, W)); shown.append(lab)
        else: ax.imshow(np.full((H, W, 3), 255, dtype=np.uint8)); ax.text(0.5, 0.5, "structure view pending", ha="center", va="center", fontsize=9, color="0.4", transform=ax.transAxes); miss.append(lab)
    els = [e for e in ("Co", "Sb", "O", "F") if stem.startswith("Co")] + [e for e in ("Ba", "Al") if stem.startswith("Ba")]
    handles = [plt.Line2D([], [], marker="o", ls="", ms=9, color=colour(ref_arr, e), markeredgecolor="k", label=e) for e in els]
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.legend(handles=handles, loc="lower center", ncol=len(els), frameon=False, fontsize=11); fig.savefig(FIG / f"case_{stem}.png", dpi=200); fig.savefig(FIG / f"case_{stem}.pdf"); plt.close(fig); print("saved", stem, shown, "missing:", miss)
