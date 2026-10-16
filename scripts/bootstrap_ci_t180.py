"""Percentile bootstrap (over the 180 targets) 95% intervals for the exported top-1/top-5/top-20 match rates.

Targets without a validated candidate count as failures (denominator 180 for every algorithm). The same resampled
target indices are used for all algorithms, 10,000 resamples, seed 2026."""
import numpy as np, pandas as pd
from pathlib import Path
R = Path(__file__).resolve().parent.parent / "results/new_t180"
KEY = {"TCSPv3": "tcspv3", "CrystalFlow": "crystalflow", "Crystalite": "crystalite", "OMatG": "omatg", "Uni-3DAR": "uni3dar", "uFlowCSP": "uflowcsp",
       "EquiCSP": "equicsp", "MatterGen": "mattergen", "DiffCSP": "diffcsp", "cond-CDVAE": "cond-cdvae", "GemsDiff": "gemsdiff", "TGDMat": "tgdmat",
       "SymmCD": "symmcd", "CrystalGRW": "crystalgrw", "CrystalFormer": "crystalformer", "CDVAE": "cdvae"}
def main(B=10000):
    t = pd.read_csv(R / "t180_topk_sg_per_case_all.csv"); ids = sorted(t.material_id.unique()); n = len(ids)
    idx = np.random.default_rng(2026).integers(0, n, size=(B, n)); rows = []
    for name, a in KEY.items():
        g = t[t.algorithm == a].set_index("material_id").reindex(ids); r = g.sm_rank.to_numpy(dtype=float)
        for k in (1, 5, 20):
            hit = (r <= k).astype(float); p = hit.mean() * 100; bs = hit[idx].mean(axis=1) * 100
            lo, hi = np.percentile(bs, [2.5, 97.5]); rows.append(dict(algo=name, k=k, hits=int(hit.sum()), rate=p, lo=lo, hi=hi))
    d = pd.DataFrame(rows); d.to_csv(R / "bootstrap_ci_topk.csv", index=False); return d
if __name__ == "__main__":
    print(main().pivot(index="algo", columns="k", values=["rate", "lo", "hi"]).round(1).to_string())
