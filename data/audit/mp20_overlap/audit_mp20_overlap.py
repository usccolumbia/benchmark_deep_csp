from pathlib import Path
import csv,json,hashlib,sys,time,pickle,warnings,collections,concurrent.futures,importlib.metadata,shutil
from pymatgen.core import Structure, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher, SpeciesComparator
BASE=Path('/Users/jason/Downloads/project/2026_CSP_diagnosis_benchmark')
ROOT=BASE/'testset_4benchmark'
OUT=BASE/'results/mp20_overlap_audit_20260929'
SETS=['T46','T180','T240','T800','T2400']
PARAMS=dict(ltol=.2,stol=.3,angle_tol=5,primitive_cell=True,scale=True,attempt_supercell=False,allow_subset=False)
STRICT=dict(ltol=.01,stol=.02,angle_tol=.5,primitive_cell=True,scale=False,attempt_supercell=False,allow_subset=False)
def sha(b):return hashlib.sha256(b).hexdigest()
def key(comp):return tuple(sorted((e,round(float(v),10)) for e,v in comp.element_composition.fractional_composition.get_el_amt_dict().items()))
def parse(job):
 k,cif=job
 try:
  with warnings.catch_warnings(record=True) as ws:
   warnings.simplefilter('always');s=Structure.from_str(cif,fmt='cif')
  return k,s,[str(w.message) for w in ws],None
 except Exception as e:return k,None,[],repr(e)
def writecsv(name,rows,cols=None):
 if cols is None:cols=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/name).open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=cols);w.writeheader();w.writerows(rows)
def match(job):
 h,s,refs=job;sm=StructureMatcher(**PARAMS,comparator=SpeciesComparator());strict=StructureMatcher(**STRICT,comparator=SpeciesComparator());hits=[];errors=[]
 for ref in refs:
  try:
   rms=sm.get_rms_dist(s,ref['structure']) if sm.fit(s,ref['structure']) else None
   if rms is not None:
    st=strict.get_rms_dist(s,ref['structure']) if strict.fit(s,ref['structure']) else None
    hits.append(dict(mp20_material_id=ref['material_id'],mp20_split=ref['split'],mp20_row=ref['row'],mp20_formula=ref['structure'].composition.reduced_formula,mp20_nsites=len(ref['structure']),rms_normalized=rms[0],max_displacement_normalized=rms[1],strict_no_scaling_match=st is not None,strict_rms_normalized=st[0] if st else '',volume_per_atom_ratio=(s.volume/len(s))/(ref['structure'].volume/len(ref['structure']))))
  except Exception as e:errors.append(dict(test_cif_sha256=h,mp20_material_id=ref['material_id'],error=repr(e)))
 return h,hits,errors,len(refs)
def main():
 t=time.time();OUT.mkdir(parents=True,exist_ok=True);shutil.copy2(__file__,OUT/'audit_mp20_overlap.py');csv.field_size_limit(10000000)
 source=[];tests=[];mp=[];issues=list(csv.DictReader((BASE/'testset_issues.csv').open()));text_by_hash={};integrity=[]
 for name in SETS:
  folder=ROOT/name;records=json.loads((folder/'records.json').read_text());meta=list(csv.DictReader((folder/'metadata.csv').open()));meta_by_id={r['material_id']:r for r in meta};actual=set(p.stem for p in (folder/'cifs').glob('*.cif'))
  for fname in ['records.json','metadata.csv','manifest.json']:
   p=folder/fname;source.append(dict(path=str(p),sha256=sha(p.read_bytes()),bytes=p.stat().st_size))
  integrity.append(dict(test_set=name,records=len(records),metadata_rows=len(meta),cif_files=len(actual),unique_ids=len(set(r['material_id'] for r in records)),missing_cif_ids=sorted(set(r['material_id'] for r in records)-actual),extra_cif_ids=sorted(actual-set(r['material_id'] for r in records))))
  for i,r in enumerate(records):
   p=folder/'cifs'/(r['material_id']+'.cif');b=p.read_bytes();h=sha(b);text_by_hash[h]=b.decode();md=meta_by_id.get(r['material_id'],{});diff=[k for k in md if str(r.get(k,''))!=md[k]]
   tests.append(dict(test_set=name,record_row=i,material_id=r['material_id'],cif_path=str(p),cif_sha256=h,embedded_cif_byte_match=b==r['cif'].encode(),metadata_json_differing_fields=';'.join(diff),metadata_formula=r.get('primitive_formula') or r.get('formula') or r.get('pretty_formula'),metadata_spacegroup=r.get('spacegroup'),record=r))
 for split in ['train','val']:
  p=BASE/'mp_20'/f'{split}.csv';source.append(dict(path=str(p),sha256=sha(p.read_bytes()),bytes=p.stat().st_size))
  for i,r in enumerate(csv.DictReader(p.open())):
   h=sha(r['cif'].encode());text_by_hash[h]=r['cif'];mp.append(dict(split=split,row=i,material_id=r['material_id'],cif_sha256=h,metadata_formula=r['pretty_formula']))
 source.append(dict(path=str(BASE/'testset_issues.csv'),sha256=sha((BASE/'testset_issues.csv').read_bytes())))
 (OUT/'source_manifest.json').write_text(json.dumps(source,indent=2))
 print('Inputs:',len(tests),'test rows;',len(mp),'MP20 rows;',len(text_by_hash),'unique CIF texts',flush=True)
 parsed={};parse_errors=[];parser_warnings=[]
 cache=OUT/'parsed_structures.pkl';fingerprint=sha(json.dumps(source,sort_keys=True).encode())
 if cache.exists():
  data=pickle.loads(cache.read_bytes())
  if data['fingerprint']==fingerprint:parsed=data['parsed'];parse_errors=data['errors'];parser_warnings=data['warnings'];print('Resumed parsed structures',len(parsed),flush=True)
 if not parsed:
  with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
   for i,(h,s,ws,err) in enumerate(pool.map(parse,text_by_hash.items(),chunksize=32),1):
    if err:parse_errors.append(dict(cif_sha256=h,error=err))
    else:parsed[h]=s
    if ws:parser_warnings.append(dict(cif_sha256=h,warnings=' | '.join(ws)))
    if i%3000==0:print('Parsed',i,'/',len(text_by_hash),'elapsed',round(time.time()-t),flush=True)
  cache.write_bytes(pickle.dumps(dict(fingerprint=fingerprint,parsed=parsed,errors=parse_errors,warnings=parser_warnings)))
 writecsv('parse_errors.csv',parse_errors,['cif_sha256','error']);writecsv('parser_warnings.csv',parser_warnings,['cif_sha256','warnings'])
 by_comp=collections.defaultdict(list);by_id=collections.defaultdict(list)
 for r in mp:
  by_id[r['material_id']].append(r)
  if r['cif_sha256'] in parsed:
   r['structure']=parsed[r['cif_sha256']];by_comp[key(r['structure'].composition)].append(r)
 metadata_issues=[];inventory=[];id_hits=[];jobs={}
 for r in tests:
  h=r['cif_sha256'];s=parsed.get(h);row={k:v for k,v in r.items() if k!='record'}
  if s is not None:
   row.update(actual_formula=s.composition.reduced_formula,actual_nsites=len(s),same_composition_mp20_rows=len(by_comp[key(s.composition)]))
   try:row['metadata_formula_matches_cif']=key(Composition(r['metadata_formula']))==key(s.composition)
   except Exception as e:row['metadata_formula_matches_cif']=False;row['metadata_formula_error']=repr(e)
   if not row['metadata_formula_matches_cif']:
    with warnings.catch_warnings():
     warnings.simplefilter('ignore'); sg=s.get_space_group_info(symprec=.01)
    metadata_issues.append(dict(test_set=r['test_set'],material_id=r['material_id'],metadata_formula=r['metadata_formula'],actual_formula=s.composition.reduced_formula,metadata_spacegroup=r['metadata_spacegroup'],actual_spacegroup_symprec_0p01=sg[1],cif_path=r['cif_path'],cif_sha256=h,embedded_cif_byte_match=r['embedded_cif_byte_match']))
   if h not in jobs:jobs[h]=(h,s,by_comp[key(s.composition)])
  inventory.append(row)
  for other in by_id.get(r['material_id'],[]):id_hits.append(dict(test_set=r['test_set'],test_material_id=r['material_id'],mp20_split=other['split'],mp20_material_id=other['material_id']))
 writecsv('test_inventory.csv',inventory);writecsv('metadata_mismatches.csv',metadata_issues);writecsv('id_overlaps.csv',id_hits,['test_set','test_material_id','mp20_split','mp20_material_id'])
 print('ID overlaps',len(id_hits),'metadata formula mismatches',len(metadata_issues),'structure candidate pairs',sum(len(j[2]) for j in jobs.values()),flush=True)
 all_hits={};matching_errors=[];tested=0
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for i,(h,hits,errs,n) in enumerate(pool.map(match,jobs.values(),chunksize=4),1):
   all_hits[h]=hits;matching_errors.extend(errs);tested+=n
   if i%400==0:print('Matched test CIFs',i,'/',len(jobs),'candidate pairs',tested,'elapsed',round(time.time()-t),flush=True)
 overlap=[]
 for r in tests:
  for hit in all_hits.get(r['cif_sha256'],[]):
   overlap.append(dict(test_set=r['test_set'],test_material_id=r['material_id'],test_formula=parsed[r['cif_sha256']].composition.reduced_formula,test_nsites=len(parsed[r['cif_sha256']]),same_material_id=r['material_id']==hit['mp20_material_id'],test_cif_sha256=r['cif_sha256'],**hit))
 writecsv('structure_overlaps.csv',overlap);writecsv('matching_errors.csv',matching_errors,['test_cif_sha256','mp20_material_id','error'])
 reconciled=[]
 for row in issues:
  for ts in row['test_set'].split(';'):
   if row['issue']=='duplicate_of_mp20_structure':
    found=[x for x in overlap if x['test_set']==ts and x['test_material_id']==row['test_material_id'] and x['mp20_material_id']==row['matched_mp20_material_id'] and x['mp20_split']==row['mp20_split']]
   else:found=[x for x in metadata_issues if x['test_set']==ts and x['material_id']==row['test_material_id']]
   reconciled.append(dict(**{**row,'test_set':ts},reproduced=bool(found),rms_normalized=found[0].get('rms_normalized','') if found else '',strict_no_scaling_match=found[0].get('strict_no_scaling_match','') if found else ''))
 writecsv('supplied_report_reconciliation.csv',reconciled)
 summary=[]
 for ts in SETS:
  rows=[r for r in tests if r['test_set']==ts];ov=[r for r in overlap if r['test_set']==ts]
  summary.append(dict(test_set=ts,test_records=len(rows),id_overlap_records=len({r['test_material_id'] for r in id_hits if r['test_set']==ts}),structure_overlap_records=len({r['test_material_id'] for r in ov}),train_overlap_records=len({r['test_material_id'] for r in ov if r['mp20_split']=='train'}),val_overlap_records=len({r['test_material_id'] for r in ov if r['mp20_split']=='val'}),strict_structure_overlap_records=len({r['test_material_id'] for r in ov if r['strict_no_scaling_match']}),matched_pairs=len(ov),metadata_formula_mismatches=len([r for r in metadata_issues if r['test_set']==ts])))
 writecsv('summary.csv',summary)
 result=dict(status='complete' if not parse_errors and not matching_errors else 'complete_with_errors',testsets=summary,mp20_rows=dict(collections.Counter(r['split'] for r in mp)),mp20_unique_material_ids=len({r['material_id'] for r in mp}),mp20_train_val_shared_ids=len({r['material_id'] for r in mp if r['split']=='train'}&{r['material_id'] for r in mp if r['split']=='val'}),test_records=len(tests),unique_test_cif_texts=len(jobs),same_composition_pairs_tested=tested,parse_errors=len(parse_errors),matching_errors=len(matching_errors),supplied_report_rows=len(issues),supplied_report_expanded_rows=len(reconciled),supplied_report_reproduced_rows=sum(r['reproduced'] for r in reconciled),structure_overlap_unique_test_ids=len({r['test_material_id'] for r in overlap}),input_integrity=integrity,all_cifs_equal_json=all(r['embedded_cif_byte_match'] for r in tests),all_metadata_equal_json=all(not r['metadata_json_differing_fields'] for r in tests),matcher=PARAMS,strict_matcher=STRICT,comparator='SpeciesComparator',decision='StructureMatcher.fit(test, mp20), symmetric=False; get_rms_dist used only for diagnostics after fit passes',composition_preselection='Actual CIF elemental fractional composition, rounded to 10 decimals; no metadata formula prefilter',mp20_cif_column='cif',pymatgen_version=importlib.metadata.version('pymatgen'),python=sys.version,elapsed_seconds=round(time.time()-t,2))
 (OUT/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
