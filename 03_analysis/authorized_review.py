"""Automated replay of source and saved-prediction checks.

Does not import the modeling code, retrain, or modify frozen files. The checks
support an AI-assisted audit; they are not an independent human or laboratory review.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv, json, math, statistics, re, hashlib, datetime

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'01_sources/ocx/data/experimental_data'
OUT=ROOT/'04_results/authorized_review'
OUT.mkdir(exist_ok=True)

def readcsv(path):
    with path.open(encoding='utf-8-sig',newline='') as fh:
        reader=csv.DictReader(fh)
        rows=[]
        for record in reader:
            record['_csv_line']=reader.line_num
            rows.append(record)
        return rows

def save(obj,name):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')

def close(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-10)

def prep(s): return s.rsplit('_rep',1)[0]
def key(r): return r['sample id'] if r['source']=='vsp' else prep(r['sample id'])
def missing(s): return s.strip().lower() in ('','nan')
def boolean(s): return s=='True'

raw=readcsv(SRC/'ExpDataDump_241113.csv')
numeric=[r for r in raw if r['rep'].isdigit()]
her=[r for r in numeric if r['reaction']=='HER' and float(r['current density'])==50]
co=[r for r in numeric if r['reaction']=='CO2R' and float(r['current density'])==200]
assert len({r['sample id'] for r in her})==len(her)
assert len({r['sample id'] for r in co})==len(co)
H={r['sample id']:r for r in her}; C={r['sample id']:r for r in co}
assert H.keys()==C.keys()
features=readcsv(ROOT/'02_data/features.csv'); F={r['electrode_id']:r for r in features}
packet=readcsv(ROOT/'02_data/key_records_review_packet.csv')
targets=[]
for partition in ('development','holdout'):
    for r in readcsv(ROOT/f'02_data/{partition}_targets.csv'):
        r['_partition']=partition; targets.append(r)
T={r['electrode_id']:r for r in targets}
assert len(features)==len(F)==len(targets)==len(T)==441
assert len(packet)==len({r['electrode_id'] for r in packet})==20
checks=[]
identity=['sample id','source','batch number','batch date','composition','post processing id','rep']
for eid,f in F.items():
    h,c,t=H[eid],C[eid],T[eid]
    assert all(h[col]==c[col] for col in identity),eid
    assert c['sample id']==f['electrode_id']
    assert f['prep_group']==prep(eid) and f['source']==c['source']
    assert int(f['replicate'])==int(c['rep']) and f['nominal_composition']==c['composition']
    assert int(f['batch_number'])==int(c['batch number']) and int(f['batch_date_raw'])==int(c['batch date'])
    assert f['batch_id']==c['source']+'_'+str(int(c['batch number']))
    assert eid==f"{c['source']}{int(c['batch number'])}_{int(c['batch date'])}_{c['composition']}_{c['post processing id']}_rep{int(c['rep'])}"
    assert int(f['her_csv_line'])==h['_csv_line'] and int(f['co2_csv_line'])==c['_csv_line']
    assert close(f['her_voltage_raw_V'],h['voltage']) and close(f['her_voltage_V'],abs(float(h['voltage'])))
    for dst,src in [('co_fe_pct','fe_co'),('h2_fe_pct','fe_h2'),('gas_total_pct','fe_gas_total'),('liquid_residual_pct','fe_liquid')]:
        assert close(t[dst],c[src]),(eid,dst)
    assert close(t['co2_voltage_V'],abs(float(c['voltage'])))
    assert f['partition']==t['_partition']
    tokens=c['composition'].split('-'); composition={tokens[i]:float(tokens[i+1]) for i in range(0,len(tokens),2)}
    total=sum(composition.values())
    for col in f:
        if re.fullmatch(r'nominal_[A-Z][a-z]?',col):
            assert close(f[col],composition.get(col[8:],0)/total),(eid,col)
    assert f['element_group']=='-'.join(sorted(e for e,v in composition.items() if v>0))
    checks.append({'electrode_id':eid,'source':c['source'],'replicate':int(c['rep']),
                   'her_csv_line':h['_csv_line'],'co2_csv_line':c['_csv_line'],
                   'HER_current_mA_cm2':float(h['current density']),'CO2_current_mA_cm2':float(c['current density']),
                   'HER_raw_V':h['voltage'],'HER_absolute_V':f['her_voltage_V'],'CO_FE_pct':c['fe_co'],
                   'identity_values_composition_and_lines':'pass'})
review_by_id={r['electrode_id']:r for r in checks}
key_checks=[]
for i,p in enumerate(packet,1):
    eid=p['electrode_id']; f,t=F[eid],T[eid]
    assert f['partition']=='development'
    for col in ('source','prep_group','nominal_composition','replicate','her_csv_line','co2_csv_line'):
        assert p[col]==f[col],(eid,col)
    assert close(p['her_voltage_V'],f['her_voltage_V']) and close(p['co_fe_pct'],t['co_fe_pct'])
    key_checks.append({'number':i,**review_by_id[eid]})
save(checks,'all_441_pair_checks.json'); save(key_checks,'key_20_record_checks.json')

# Reconstruct the ambiguity inventory from both source tables, not the existing
# ambiguity report, and check every feature row's chosen or missing mapping.
xrf=readcsv(SRC/'XRFDataDump-241113.csv'); xrd=readcsv(SRC/'XRDDataDump-241113.csv')
XR=defaultdict(list); XD=defaultdict(list)
for r in xrf: XR[key(r)].append(r)
for r in xrd: XD[key(r)].append(r)
detected=[]
for k,rs in sorted(XR.items()):
    vals=sorted({r['xrf-gde_formula'] if r['source']=='vsp' else r['xrf-wafer_formula'] for r in rs}-{''})
    if len(vals)>1:
        affected=[f['electrode_id'] for f in features if f['xrf_mapping_key']==k]
        detected.append({'mapping_key':k,'type':'conflicting_XRF','formulas':vals,'affected_electrodes':affected,
                         'source_records':[{'line':r['_csv_line'],'id':r['sample id'],'source':r['source'],'wafer':r['xrf-wafer_formula'],'GDE':r['xrf-gde_formula']} for r in rs],
                         'decision':'exclude ambiguous XRF, preserve nominal analysis'})
for k,rs in sorted(XD.items()):
    if len(rs)>1:
        detected.append({'mapping_key':k,'type':'multiple_XRD','count':len(rs),
                         'affected_electrodes':[f['electrode_id'] for f in features if f['xrf_mapping_key']==k],
                         'source_records':[{'line':r['_csv_line'],'id':r['sample id'],'source':r['source'],'target':r['target comp'],'XRF':r['xrf comp'],'solutions':r['total no of solutions']} for r in rs],
                         'decision':'exclude ambiguous XRD; no XRD model was fitted'})
expected=json.loads((ROOT/'02_data/mapping_ambiguities.json').read_text(encoding='utf-8'))
assert {(r['mapping_key'],r['type']) for r in detected}=={(r['mapping_key'],r['type']) for r in expected}
for a in expected:
    d=next(r for r in detected if (r['mapping_key'],r['type'])==(a['mapping_key'],a['type']))
    if 'formulas' in a: assert sorted(a['formulas'])==d['formulas']
    if 'count' in a: assert a['count']==d['count']
for f in features:
    k=f['electrode_id'] if f['source']=='vsp' else f['prep_group']
    assert k==f['xrf_mapping_key']
    rs=XR.get(k,[]); field='xrf-gde_formula' if f['source']=='vsp' else 'xrf-wafer_formula'
    vals={r[field] for r in rs}-{''}
    assert boolean(f['xrf_available'])==(len(vals)==1)
    if len(vals)==1:
        assert f['xrf_composition']==next(iter(vals))
        chosen=next(r for r in rs if r['_csv_line']==int(float(f['xrf_csv_line'])))
        assert chosen['sample id']==f['xrf_record_id'] and chosen['source']==f['source'] and chosen[field]==f['xrf_composition']
    else: assert missing(f['xrf_composition']) and missing(f['xrf_csv_line'])
    rows=XD.get(k,[])
    assert boolean(f['xrd_available'])==(len(rows)==1)
    if len(rows)==1: assert int(float(f['xrd_csv_line']))==rows[0]['_csv_line'] and rows[0]['source']==f['source']
    else: assert missing(f['xrd_csv_line'])
save(detected,'mapping_checks.json')

batches=[r for r in readcsv(SRC/'supporting_data/batches_info.csv') if r['source'] and r['batch number']]
B={r['source']+'_'+str(int(float(r['batch number']))):r for r in batches}
assert len(B)==len(batches)
conflicts=[]
for f in features:
    b=B[f['batch_id']]
    assert int(float(b['batch date']))==int(float(f['batch_date_metadata']))
    conflict=int(float(b['batch date']))!=int(f['batch_date_raw'])
    assert boolean(f['batch_date_conflict'])==conflict
    if conflict: conflicts.append((f['electrode_id'],f['batch_id'],int(f['batch_date_raw']),int(float(b['batch date']))))
registered=readcsv(ROOT/'02_data/batch_discrepancies.csv')
assert set(conflicts)=={(r['electrode_id'],r['batch_id'],int(r['batch_date_raw']),int(float(r['batch_date_metadata']))) for r in registered}
counts=Counter((b,d,m) for eid,b,d,m in conflicts)
batch_checks=[{'batch_id':b,'raw_date':d,'metadata_date':m,'electrodes':n,'decision':'retain both dates; model batch identity only'} for (b,d,m),n in sorted(counts.items())]
save(batch_checks,'batch_checks.json')
excluded=readcsv(ROOT/'02_data/exclusions.csv')
assert {r['sample id'] for r in excluded}==H.keys()-F.keys()
for e in excluded:
    eid=e['sample id']; h,c=H[eid],C[eid]
    assert missing(h['voltage']) or missing(c['fe_co'])
    assert h['_csv_line']==int(e['raw_csv_line_her']) and c['_csv_line']==int(e['raw_csv_line_co2'])
save([{'id':e['sample id'],'HER_voltage':H[e['sample id']]['voltage'],'CO_FE_pct':C[e['sample id']]['fe_co'],'decision':'retain original exclusion'} for e in excluded],'exclusion_checks.json')

# Independently calculate group-equal errors and screening from saved rows.
# No model fitting, new endpoint selection, or new holdout evaluation is done.
def score(rows):
    groups=defaultdict(list)
    for r in rows: groups[r['prep_group']].append(r)
    means={p:{k:statistics.fmean(float(r[k]) for r in rs) for k in ('y','M0','M1')} for p,rs in groups.items()}
    result={v+'_MAE':statistics.fmean(statistics.fmean(abs(float(r['y'])-float(r[v])) for r in rs) for rs in groups.values()) for v in ('M0','M1')}
    result['gain_pp_or_V']=result['M0_MAE']-result['M1_MAE']; result['prep_groups']=len(groups)
    k=math.ceil(.2*len(groups)); true=set(sorted(means,key=lambda p:(-means[p]['y'],p))[:k])
    result['screening']={v:{'top_k':k,'hits':len(true&set(sorted(means,key=lambda p:(-means[p][v],p))[:k]))} for v in ('M0','M1')}
    return result

metric_checks=[]
for path in sorted((ROOT/'04_results').glob('*/predictions.csv')):
    summary=path.parent/'summary.json'
    if not summary.exists(): continue
    rows=readcsv(path); target=json.loads((path.parent/'config.json').read_text())['target']
    for r in rows: assert close(r['y'],T[r['electrode_id']][target])
    s=score(rows); expected=json.loads(summary.read_text())
    for k in ('M0_MAE','M1_MAE','gain_pp_or_V','prep_groups'): assert close(s[k],expected[k]),(str(path),k)
    for v in ('M0','M1'):
        for k in ('top_k','hits'): assert s['screening'][v][k]==expected['screening'][v][k]
    metric_checks.append({'task':path.parent.name,**s,'status':'pass'})
ho=readcsv(ROOT/'04_results/holdout_evaluated_predictions.csv')
groups=defaultdict(list)
for r in ho:
    assert close(r['y'],T[r['electrode_id']][r['target']])
    groups[(r['budget'],r['algorithm'],r['target'])].append(r)
expected=json.loads((ROOT/'04_results/holdout_results.json').read_text())
for k,rows in groups.items():
    s=score(rows); ref=next(r for r in expected if (r['budget'],r['algorithm'],r['target'])==k)
    for field in ('M0_MAE','M1_MAE','gain_pp_or_V','prep_groups'): assert close(s[field],ref[field]),(k,field)
    if 'screening' in ref:
        for v in ('M0','M1'):
            for field in ('top_k','hits'): assert s['screening'][v][field]==ref['screening'][v][field]
    else:
        assert k[2]!='co_fe_pct'
        s.pop('screening')
    metric_checks.append({'task':'holdout_'+'_'.join(k),**s,'status':'pass'})
save(metric_checks,'metric_checks.json')

protocol=json.loads((ROOT/'00_admin/protocol_lock.json').read_text())
for p,digest in protocol['sha256'].items():
    assert hashlib.sha256((ROOT/Path(p.replace('\\','/'))).read_bytes()).hexdigest()==digest
frozen=json.loads((ROOT/'00_admin/predictions_frozen_before_unseal.json').read_text())
assert hashlib.sha256((ROOT/'04_results/holdout_predictions_before_unseal.csv').read_bytes()).hexdigest()==frozen['predictions_sha256']
summary={'review_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Automated replay of published checks; no new human authorization or review is asserted',
         'reviewer':'Automated Python checks originally prepared with Codex','independent_human_review':False,'raw_pair_records_checked':len(checks),
         'designated_key_records_checked':len(key_checks),'mapping_anomalies_checked':len(detected),
         'batch_date_conflict_electrodes_checked':len(conflicts),'batch_conflict_groups':len(batch_checks),
         'exclusions_checked':len(excluded),'saved_metric_tasks_checked':len(metric_checks),
         'frozen_inputs_and_predictions_unchanged':True,'discrepancies_in_checks':0,
         'models_refitted':0,'status':'source_and_calculation_checks_passed',
         'review_disclosure':'See AUTHORSHIP.md in the release; this script generates computational checks only'}
save(summary,'summary.json')
print(json.dumps({'summary':summary,'key_records':key_checks,'mapping_records':detected,'batch_records':batch_checks},ensure_ascii=False,indent=2))
