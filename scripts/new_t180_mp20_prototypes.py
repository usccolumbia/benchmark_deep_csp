"""For each T180 target, count MP20 train+val structures with the same sorted
reduced stoichiometric coefficients AND same space group (spglib symprec 0.01)."""
import json, re, sys, warnings
import pandas as pd
from multiprocessing import Pool
from pathlib import Path
warnings.filterwarnings("ignore")
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
ROOT = Path(__file__).resolve().parent.parent
CAMP = Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31")

def sig(comp): return tuple(sorted(float(v) for v in Composition(comp).reduced_composition.values()))

def sgn(cif):
    try: return SpacegroupAnalyzer(Structure.from_str(cif, fmt="cif"), symprec=.01).get_space_group_number()
    except Exception: return None

if __name__ == "__main__":
    qs = json.loads((CAMP / "queries.json").read_text())
    tsg = json.loads((ROOT / "results/tcspv3_prototype_excluded_t180/target_sg.json").read_text())
    tsig = {q["material_id"]: sig(q["primitive_cell_formula"]) for q in qs}
    need = set(tsig.values())
    df = pd.concat([pd.read_csv(ROOT / "mp_20" / f) for f in ("train.csv", "val.csv")])
    df["sig"] = df.cif.map(lambda c: sig(Structure.from_str(c, fmt="cif").composition.reduced_formula))
    sub = df[df.sig.isin(need)].copy()
    print(len(df), len(sub), flush=True)
    with Pool(8) as p: sub["sg"] = p.map(sgn, sub.cif.tolist(), chunksize=20)
    out = {}
    for mid, s in tsig.items():
        m = sub[(sub.sig == s) & (sub.sg == tsg.get(mid)) & (sub.material_id != mid)]
        out[mid] = dict(n_sig=int((sub.sig == s).sum()), n_same_sg=int(len(m)))
    (ROOT / "results/new_t180/mp20_prototype_counts.json").write_text(json.dumps(out, indent=1))
    print(sum(1 for v in out.values() if v["n_same_sg"] > 0), "of", len(out), "targets have an MP20 train/val prototype")
