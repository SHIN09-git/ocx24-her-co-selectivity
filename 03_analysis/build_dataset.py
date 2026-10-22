"""Build electrode-level paired records; never print held-out outcomes."""
from pathlib import Path
import pandas as pd, numpy as np, json, re, hashlib, datetime, ast
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'01_sources/ocx/data/experimental_data'
OUT=ROOT/'02_data'; OUT.mkdir(exist_ok=True)
SEED=20261012
def dump(obj,path):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x: int(x) if isinstance(x,np.integer) else str(x)),encoding='utf-8')
def formula(s):
    if pd.isna(s): return None
    p=str(s).split('-')
    if len(p)%2: raise ValueError(f'Invalid composition: {s}')
    out={}
    for e,x in zip(p[::2],p[1::2]):
        assert re.fullmatch('[A-Z][a-z]?',e),s
        out[e]=out.get(e,0)+float(x)
    assert min(out.values())>=0 and abs(sum(out.values())-1)<.01,s
    return {e:x/sum(out.values()) for e,x in out.items() if x>0}
def prep(s): return re.sub(r'_rep\d+$','',s)
raw=pd.read_csv(SRC/'ExpDataDump_241113.csv')
raw['raw_csv_line']=np.arange(len(raw))+2
numeric=raw[pd.to_numeric(raw.rep,errors='coerce').notna()].copy()
numeric['rep']=numeric.rep.astype(int)
assert not numeric.duplicated(['sample id','reaction','current density']).any()
h=numeric[(numeric.reaction=='HER') & (numeric['current density']==50)].copy()
c=numeric[(numeric.reaction=='CO2R') & (numeric['current density']==200)].copy()
assert h['sample id'].is_unique and c['sample id'].is_unique
assert set(h['sample id'])==set(c['sample id'])
hc=c.merge(h,on='sample id',suffixes=('_co2','_her'),validate='one_to_one')
for col in ['source','batch number','batch date','composition','post processing id','rep']:
    assert (hc[col+'_co2']==hc[col+'_her']).all(),col
cohc=hc.copy()
excluded=cohc[cohc.voltage_her.isna() | cohc.fe_co_co2.isna()].copy()
excluded[['sample id','raw_csv_line_her','raw_csv_line_co2']].assign(reason='missing_HER_voltage_or_CO_FE').to_csv(OUT/'exclusions.csv',index=False,lineterminator='\r\n')
hc=hc[hc.voltage_her.notna() & hc.fe_co_co2.notna()].copy()
f=pd.DataFrame({'electrode_id':hc['sample id'],'prep_group':hc['sample id'].map(prep),'replicate':hc.rep_co2,'source':hc.source_co2,'batch_number':hc['batch number_co2'].astype(int),'batch_date_raw':hc['batch date_co2'].astype(int),'nominal_composition':hc.composition_co2,'her_voltage_raw_V':hc.voltage_her,'her_voltage_V':hc.voltage_her.abs(),'her_csv_line':hc.raw_csv_line_her,'co2_csv_line':hc.raw_csv_line_co2})
f['batch_id']=f.source+'_'+f.batch_number.astype(str)
f['nominal_parsed']=f.nominal_composition.map(formula)
f['element_group']=f.nominal_parsed.map(lambda d:'-'.join(sorted(d)))
els=sorted({e for d in f.nominal_parsed for e in d})
for e in els: f['nominal_'+e]=f.nominal_parsed.map(lambda d:d.get(e,0.))
xrf=pd.read_csv(SRC/'XRFDataDump-241113.csv'); xrf['xrf_csv_line']=np.arange(len(xrf))+2
xrf['mapping_key']=[sid if s=='vsp' else prep(sid) for sid,s in zip(xrf['sample id'],xrf.source)]
xrf['chosen_formula']=np.where(xrf.source=='vsp',xrf['xrf-gde_formula'],xrf['xrf-wafer_formula'])
xrf_mapping={}; ambiguities=[]
for key,g in xrf.groupby('mapping_key'):
    vals=g.chosen_formula.dropna().unique()
    if len(vals)==1:
        r=g[g.chosen_formula==vals[0]].iloc[0]
        xrf_mapping[key]={'formula':vals[0],'line':r.xrf_csv_line,'record_id':r['sample id'],'source':r.source}
    elif len(vals)>1: ambiguities.append({'mapping_key':key,'formulas':vals.tolist(),'type':'conflicting_XRF'})
f['xrf_mapping_key']=np.where(f.source=='vsp',f.electrode_id,f.prep_group)
f['xrf_composition']=f.xrf_mapping_key.map(lambda k:xrf_mapping.get(k,{}).get('formula'))
f['xrf_csv_line']=f.xrf_mapping_key.map(lambda k:xrf_mapping.get(k,{}).get('line'))
f['xrf_record_id']=f.xrf_mapping_key.map(lambda k:xrf_mapping.get(k,{}).get('record_id'))
f['xrf_substrate']=np.where(f.source=='vsp','GDE','wafer_proxy')
f['xrf_parsed']=f.xrf_composition.map(formula)
for e in els: f['xrf_'+e]=f.xrf_parsed.map(lambda d: np.nan if d is None else d.get(e,0.))
f['xrf_available']=f.xrf_composition.notna()
clean=pd.read_csv(SRC/'ExpDataDump_241113_clean.csv')
cl=clean[(clean.reaction=='HER')&(clean['current density']==50)].set_index('sample id')
f['in_official_clean']=f.electrode_id.isin(cl.index)
disagree=[]
for _,r in f[f.in_official_clean].iterrows():
    a=formula(cl.loc[r.electrode_id,'xrf comp']); b=r.xrf_parsed
    if a!=b: disagree.append({'electrode_id':r.electrode_id,'official':a,'reconstructed':b})
# Preserve ambiguous mappings instead of reproducing an implicit last-row overwrite.
dump(disagree,OUT/'xrf_official_disagreements.json')
# XRD uses same wafer characterization across uoft repeat electrodes.
xrd=pd.read_csv(SRC/'XRDDataDump-241113.csv'); xrd['xrd_csv_line']=np.arange(len(xrd))+2
xrd['mapping_key']=[sid if s=='vsp' else prep(sid) for sid,s in zip(xrd['sample id'],xrd.source)]
xrd_map={}
for key,g in xrd.groupby('mapping_key'):
    if len(g)==1: xrd_map[key]=g.iloc[0]
    else: ambiguities.append({'mapping_key':key,'count':len(g),'type':'multiple_XRD'})
def getxrd(k,col):
    return xrd_map[k][col] if k in xrd_map else np.nan
f['xrd_csv_line']=f.xrf_mapping_key.map(lambda k:getxrd(k,'xrd_csv_line'))
f['xrd_available']=f.xrd_csv_line.notna()
f['xrd_n_solutions']=f.xrf_mapping_key.map(lambda k:getxrd(k,'total no of solutions'))
def minlist(v):
    if pd.isna(v): return np.nan
    a=ast.literal_eval(v)
    return min(a) if len(a) else np.nan
f['xrd_min_rwp']=f.xrf_mapping_key.map(lambda k:minlist(getxrd(k,'allsolutions_rwp')))
f['xrd_target_qscore']=f.xrf_mapping_key.map(lambda k:getxrd(k,'solutions_target_qscore'))
batches=pd.read_csv(SRC/'supporting_data/batches_info.csv').dropna(subset=['source','batch number'])
batches['batch_id']=batches.source+'_'+batches['batch number'].astype(int).astype(str)
assert batches.batch_id.is_unique
bmap=batches.set_index('batch_id')['batch date'].astype(int).to_dict()
f['batch_date_metadata']=f.batch_id.map(bmap)
f['batch_date_conflict']=f.batch_date_metadata.notna() & (f.batch_date_metadata!=f.batch_date_raw)
f['batch_metadata_available']=f.batch_date_metadata.notna()
# Existing small DFT summaries, strict exact XRF/source match only; no outcomes used.
processed=pd.read_csv(ROOT/'01_sources/ocx/data/processed_data/CO2R_40_70_all.csv')
energy_cols=[e+'_mean_energy' for e in ['C','CO','CHO','COCOH','H','OH']]
def canon(d): return None if d is None else tuple((e,round(d[e],8)) for e in sorted(d))
processed['xrf_key']=processed['xrf comp'].map(formula).map(canon)
emap={}
for key,g in processed.groupby(['source','xrf_key']):
    if len(g)==1: emap[key]=g.iloc[0]
for col in energy_cols:
    f['dft_'+col]=[emap.get((r.source,canon(r.xrf_parsed)),{}).get(col,np.nan) for _,r in f.iterrows()]
f['dft_summary_available']=f[['dft_'+x for x in energy_cols]].notna().all(axis=1)
f['dft_summary_row']=[emap.get((r.source,canon(r.xrf_parsed)),{}).get('sample id',None) for _,r in f.iterrows()]
f['extended_available']=f.xrf_available & f.xrd_min_rwp.notna() & f.dft_summary_available
# No outcome-dependent group allocation; selection uses only the set of nominal elements.
families=sorted(f.element_group.unique())
rng=np.random.default_rng(SEED)
test_families=sorted(rng.choice(families,size=int(np.ceil(.2*len(families))),replace=False).tolist())
f['partition']=np.where(f.element_group.isin(test_families),'holdout','development')
assert f.groupby('prep_group').element_group.nunique().max()==1
assert f.groupby('prep_group').partition.nunique().max()==1
f=f.drop(columns=['nominal_parsed','xrf_parsed']).sort_values('electrode_id').reset_index(drop=True)
f.to_csv(OUT/'features.csv',index=False,lineterminator='\r\n')
targets=hc.set_index('sample id')[['fe_co_co2','fe_h2_co2','voltage_co2','fe_gas_total_co2','fe_liquid_co2']].rename(columns={'fe_co_co2':'co_fe_pct','fe_h2_co2':'h2_fe_pct','voltage_co2':'co2_voltage_V','fe_gas_total_co2':'gas_total_pct','fe_liquid_co2':'liquid_residual_pct'})
targets['co2_voltage_V']=targets.co2_voltage_V.abs()
targets.index.name='electrode_id'
targets=targets.reset_index().merge(f[['electrode_id','partition']],on='electrode_id',validate='one_to_one')
for part in ['development','holdout']:
    targets[targets.partition==part].drop(columns='partition').to_csv(OUT/(part+'_targets.csv'),index=False,lineterminator='\r\n')
# Secondary outcomes kept in separate files, same family partition at every current.
long=numeric[numeric.reaction=='CO2R'].merge(f[['electrode_id','partition']],left_on='sample id',right_on='electrode_id',validate='many_to_one')
for part in ['development','holdout']:
    long[long.partition==part][['electrode_id','current density','fe_co','fe_h2','voltage','raw_csv_line']].to_csv(OUT/(part+'_currents.csv'),index=False,lineterminator='\r\n')
counts=lambda d: {'electrodes':len(d),'prep_groups':d.prep_group.nunique(),'element_groups':d.element_group.nunique(),'replicated_groups':int((d.groupby('prep_group').size()>1).sum())}
audit={'raw_rows':len(raw),'numeric_rows':len(numeric),'raw_pairs':len(cohc),'excluded':len(excluded),'usable':counts(f),'official_clean':counts(f[f.in_official_clean]),'xrf':counts(f[f.xrf_available]),'xrd':counts(f[f.xrd_available]),'extended':counts(f[f.extended_available]),'development':counts(f[f.partition=='development']),'holdout':counts(f[f.partition=='holdout']),'sources':f.source.value_counts().to_dict(),'batches':f.batch_id.nunique(),'batch_date_conflict_electrodes':int(f.batch_date_conflict.sum()),'batch_metadata_missing_electrodes':int((~f.batch_metadata_available).sum()),'xrf_official_disagreements':len(disagree),'ambiguous_mapping_count':len(ambiguities),'elements':els,'holdout_element_groups':test_families,'seed':SEED,'quality_scope':'global data-quality counts only, no held-out outcome summaries inspected','gas_FE_gt100_rows_all_currents':int((numeric[numeric.reaction=='CO2R'].fe_gas_total>100).sum()),'liquid_residual_lt0_rows_all_currents':int((numeric[numeric.reaction=='CO2R'].fe_liquid<0).sum())}
dump(audit,ROOT/'04_results/input_audit.json')
dump(ambiguities,OUT/'mapping_ambiguities.json')
f[f.batch_date_conflict|~f.batch_metadata_available][['electrode_id','batch_id','batch_date_raw','batch_date_metadata']].to_csv(OUT/'batch_discrepancies.csv',index=False,lineterminator='\r\n')
f.groupby(['partition','element_group','source']).agg(electrodes=('electrode_id','size'),prep_groups=('prep_group','nunique')).reset_index().to_csv(OUT/'group_dictionary.csv',index=False,lineterminator='\r\n')
split={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seed':SEED,'rule':'10 randomly selected sorted nominal element sets; no target-dependent adjustment','holdout_element_groups':test_families,'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'features.csv',OUT/'development_targets.csv',OUT/'holdout_targets.csv']}}
dump(split,ROOT/'00_admin/split_manifest.json')
# A reproducible 20-record source check packet restricted to development data.
packet=f[f.partition=='development'].sample(n=20,random_state=SEED).merge(targets.drop(columns='partition'),on='electrode_id',validate='one_to_one')
packet['automated_identity_check']='pass_exact_id_source_rep_composition_current'
packet['human_review_status']='pending_not_claimed'
packet.to_csv(OUT/'key_records_review_packet.csv',index=False,lineterminator='\r\n')
print(json.dumps(audit,ensure_ascii=False,indent=2,default=int))
