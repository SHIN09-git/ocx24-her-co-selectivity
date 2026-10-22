"""Independent raw-record and leakage/metric checks, no holdout targets opened."""
from analysis_core import *
import csv
checks=[]
def check(name,condition):
    if not condition: raise AssertionError(name)
    checks.append({'check':name,'passed':True})
verify_lock(); check('frozen protocol and split hashes',True)
d=load_development()
f=pd.read_csv(ROOT/'02_data/features.csv')
check('one row per electrode',f.electrode_id.is_unique)
check('all preparation groups stay in one partition',f.groupby('prep_group').partition.nunique().max()==1)
check('all element groups stay in one partition',f.groupby('element_group').partition.nunique().max()==1)
check('HER voltage magnitude positive',(f.her_voltage_V>0).all())
for prefix in ['nominal','xrf']:
    cols=[c for c in f if re.fullmatch(prefix+r'_[A-Z][a-z]?',c)]
    z=f if prefix=='nominal' else f[f.xrf_available]
    check(prefix+' composition sums to unity',np.allclose(z[cols].sum(axis=1),1))
for grouping in ['element_group','prep_group','batch_id','source']:
    for ti,vi in cv_splits(d,grouping):
        a=d.iloc[ti]; b=d.iloc[vi]
        check(grouping+' outer fold prep isolation',not (set(a.prep_group)&set(b.prep_group)))
        check(grouping+' outer fold group isolation',not (set(a[grouping])&set(b[grouping])))
with open(ROOT/'01_sources/ocx/data/experimental_data/ExpDataDump_241113.csv',encoding='utf-8',newline='') as stream:
    rawrows=list(csv.DictReader(stream))
packet=pd.read_csv(ROOT/'02_data/key_records_review_packet.csv')
for _,r in packet.iterrows():
    hr=rawrows[int(r.her_csv_line)-2]; cr=rawrows[int(r.co2_csv_line)-2]
    check('raw CSV independent check '+r.electrode_id,
          hr['sample id']==cr['sample id']==r.electrode_id and hr['reaction']=='HER' and cr['reaction']=='CO2R'
          and float(hr['current density'])==50 and float(cr['current density'])==200
          and int(hr['rep'])==int(cr['rep'])==r.replicate
          and hr['source']==cr['source']==r.source
          and hr['composition']==cr['composition']==r.nominal_composition
          and np.isclose(abs(float(hr['voltage'])),r.her_voltage_V)
          and np.isclose(float(cr['fe_co']),r.co_fe_pct))
toy=pd.DataFrame({'prep_group':['a','a','a','b']})
check('group-equal metric independent of number of repeats',np.isclose(group_mae(toy,[0,0,0,0],[1,1,1,9]),5))
ti,vi=cv_splits(d,'element_group')[0]; tr=d.iloc[ti]; te=d.iloc[vi]
for var in ['M0','M1','source_mean','batch_mean','M0B','M1B']:
    a,b,num,cat=prepare_inputs(tr,te,'nominal',var)
    altered=te.copy(); altered['co_fe_pct']=1e9; altered['h2_fe_pct']=-1e9; altered['co2_voltage_V']=1e9
    aa,bb,_,_=prepare_inputs(tr,altered,'nominal',var)
    check(var+' features independent of evaluation outcomes',a.equals(aa) and b.equals(bb))
    check(var+' input column whitelist',set(a)==set(num+cat) and not any(c in a for c in ['co_fe_pct','h2_fe_pct','co2_voltage_V']))
    if var in ['source_mean','batch_mean']:
        altered.her_voltage_V=1e9
        _,cc,_,_=prepare_inputs(tr,altered,'nominal',var)
        check(var+' fallback uses training HER only',b.equals(cc))
# Replication of metrics from saved predictions with independent Python loops.
p=pd.read_csv(RESULTS/'nominal_random_forest_element_group_co_fe_pct/predictions.csv')
s=json.loads((RESULTS/'nominal_random_forest_element_group_co_fe_pct/summary.json').read_text())
for col in ['M0','M1']:
    errors={}
    for r in p.to_dict('records'): errors.setdefault(r['prep_group'],[]).append(abs(r['y']-r[col]))
    manual=sum(sum(v)/len(v) for v in errors.values())/len(errors)
    check(col+' reported MAE independently recomputed',np.isclose(manual,s[col+'_MAE']))
check('retained test has 10 families and 14 prep groups',f[f.partition=='holdout'].element_group.nunique()==10 and f[f.partition=='holdout'].prep_group.nunique()==14)
jsave({'checks_passed':len(checks),'checks':checks,'human_review_completed':False,'sealed_outcomes_read':False},RESULTS/'validation_checks.json')
print(json.dumps({'checks_passed':len(checks),'human_review_completed':False,'sealed_outcomes_read':False}))
