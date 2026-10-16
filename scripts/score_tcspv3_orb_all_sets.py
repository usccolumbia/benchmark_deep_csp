import json,re,warnings,sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np,pandas as pd
warnings.filterwarnings("ignore")
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
S=Path("/private/tmp/claude-501/tcspv3_orb_stage");S180=Path("/private/tmp/claude-501/tcspv3_orb")
REFS={"T180":Path("/Users/jason/Downloads/project/2026-OOD-ML-LTC/T180_campaign_20261008_v31/references.json"),"T46x13":S/"refs/T46x13_references.json","T240x60":S/"refs/T240x60_references.json","T800":S/"refs/T800_references.json"}
M=StructureMatcher(ltol=.2,stol=.3,angle_tol=5,primitive_cell=True,scale=True,attempt_supercell=False)
nat=lambda p:[int(x) if x.isdigit() else x for x in re.split(r"(\d+)",p.name)]
def sg(s):
    try:return SpacegroupAnalyzer(s,symprec=.1).get_space_group_number()
    except Exception:return None
def ranks(gt,gsg,files):
    sm=sgr=cons=None
    for i,p in enumerate(files,1):
        try:
            s=Structure.from_file(p);fit=M.fit(gt,s);ok=sg(s)==gsg
        except Exception:continue
        if fit and sm is None:sm=i
        if ok and sgr is None:sgr=i
        if fit and ok and cons is None:cons=i
    return sm,sgr,cons
def work(t):
    st,mid,cif=t
    inp=(S180/"inputs"/mid) if st=="T180" else (S/"inputs"/st/mid);out=(S180/"out"/mid) if st=="T180" else (S/"out"/st/mid)
    gt=Structure.from_str(cif,fmt="cif");g=sg(gt)
    if not inp.exists():return dict(set=st,material_id=mid,n=0)
    rep=json.loads((out/"orb_report.json").read_text())["candidates"]
    nativ=sorted(inp.glob("*.cif"),key=nat)
    rel=[out/"relaxed"/c["relaxed_file"] for c in sorted(rep,key=lambda c:c["source_file"])]
    orb=[out/"relaxed"/c["relaxed_file"] for c in sorted(rep,key=lambda c:c["orb_v3_rank"])]
    r=dict(set=st,material_id=mid,n=len(nativ),n_conv=sum(1 for c in rep if c.get("converged")))
    for tag,fs in (("nat",nativ),("rel",rel),("orb",orb)):
        a=ranks(gt,g,fs[:20]);r[tag+"_sm"],r[tag+"_sg"],r[tag+"_cons"]=a
    return r
if __name__=="__main__":
    tasks=[]
    for st in REFS:
        refs=json.loads(REFS[st].read_text());refs={r["material_id"]:r["cif"] for r in refs} if isinstance(refs,list) else refs
        for mid,cif in refs.items():tasks.append((st,mid,cif))
    print(len(tasks),"targets",flush=True)
    with Pool(8) as p:rows=p.map(work,tasks,chunksize=4)
    pd.DataFrame(rows).to_csv(S/"per_target_all.csv",index=False);print("done",flush=True)
