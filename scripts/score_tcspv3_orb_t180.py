import json,warnings,sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np,pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
OUT=Path("/private/tmp/claude-501/tcspv3_orb/out");REF=Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31/references.json")
refs={r["material_id"]:r["cif"] for r in json.loads(REF.read_text())}
M=StructureMatcher(ltol=.2,stol=.3,angle_tol=5,primitive_cell=True,scale=True,attempt_supercell=False)
def sg(s,sp=.1):
    try:return SpacegroupAnalyzer(s,symprec=sp).get_space_group_number()
    except Exception:return None
def first_ranks(gt,gsg,files):
    sm=sgr=cons=None
    for i,p in enumerate(files,1):
        try:
            s=Structure.from_file(p);fit=M.fit(gt,s);ok=sg(s)==gsg
        except Exception:continue
        if fit and sm is None:sm=i
        if ok and sgr is None:sgr=i
        if fit and ok and cons is None:cons=i
    return sm,sgr,cons
def work(mid):
    gt=Structure.from_str(refs[mid],fmt="cif");gsg=sg(gt)
    rep=json.loads((OUT/mid/"orb_report.json").read_text());c=rep["candidates"]
    by_orb=[OUT/mid/"relaxed"/x["relaxed_file"] for x in sorted(c,key=lambda x:x["orb_v3_rank"])]
    by_native=[OUT/mid/"relaxed"/x["relaxed_file"] for x in sorted(c,key=lambda x:x["source_file"])]
    o=first_ranks(gt,gsg,by_orb);n=first_ranks(gt,gsg,by_native)
    return dict(material_id=mid,n_cifs=len(c),orb_sm=o[0],orb_sg=o[1],orb_cons=o[2],nat_sm=n[0],nat_sg=n[1],nat_cons=n[2])
if __name__=="__main__":
    ids=sorted(refs)
    with Pool(8) as p:rows=p.map(work,ids,chunksize=3)
    d=pd.DataFrame(rows);d.to_csv("/private/tmp/claude-501/tcspv3_orb/per_case.csv",index=False)
    rng=np.random.default_rng(2026);idx=rng.integers(0,len(d),size=(10000,len(d)))
    for tag,lab in (("orb","TCSPv3 + ORB-v3 relax + ORB rank"),("nat","TCSPv3 + ORB-v3 relax, native order")):
        print(lab)
        for m,name in (("sm","structure"),("sg","space group"),("cons","consensus")):
            r=pd.to_numeric(d[f"{tag}_{m}"],errors="coerce")
            for k in (1,5,20):
                h=(r<=k).fillna(False).to_numpy().astype(float);lo,hi=np.percentile(h[idx].mean(1)*100,[2.5,97.5])
                print(f"  {name:12s} top-{k:<2d} {int(h.sum()):3d}/180 {h.mean()*100:5.1f}% [{lo:4.1f},{hi:4.1f}]")
    print("targets with no candidates:",int((d.n_cifs==0).sum()),"| median n:",d.n_cifs.median())
