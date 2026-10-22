"""Leakage-controlled small-sample prediction and grouped statistics."""
from pathlib import Path
import os
for _k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:
    os.environ.setdefault(_k,'1')
import json, hashlib, itertools, re, time
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
import joblib

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'04_results'
SEED=20261012
def jsave(obj,path):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)),encoding='utf-8')
def verify_lock():
    lock=json.loads((ROOT/'00_admin/protocol_lock.json').read_text())
    for p,h in lock['sha256'].items():
        assert hashlib.sha256((ROOT/Path(p.replace('\\','/'))).read_bytes()).hexdigest()==h, f'Frozen file changed: {p}'
def load_development():
    verify_lock()
    f=pd.read_csv(ROOT/'02_data/features.csv')
    y=pd.read_csv(ROOT/'02_data/development_targets.csv')
    d=f[f.partition=='development'].merge(y,on='electrode_id',validate='one_to_one')
    assert len(d)==404 and d.element_group.nunique()==38
    return d.sort_values('electrode_id').reset_index(drop=True)
def subset(d,budget):
    if budget=='xrf': d=d[d.xrf_available]
    elif budget=='extended': d=d[d.extended_available]
    return d.reset_index(drop=True)
def base_features(d,budget,variant):
    prefix='nominal' if budget=='nominal' else 'xrf'
    numeric=[c for c in d.columns if re.fullmatch(prefix+r'_[A-Z][a-z]?',c)]
    if budget=='extended':
        numeric += ['xrd_n_solutions','xrd_min_rwp','xrd_target_qscore']+[c for c in d if c.startswith('dft_') and c.endswith('_mean_energy')]
    cat=['source']+(['batch_id'] if variant in ['M0B','M1B'] else [])
    if variant not in ['M0','M0B']: numeric+=['her_input']
    assert not any(x in numeric+cat for x in ['co_fe_pct','h2_fe_pct','co2_voltage_V','prep_group','electrode_id','gas_total_pct'])
    return numeric,cat
def prepare_inputs(tr,te,budget,variant):
    tr=tr.copy(); te=te.copy()
    if variant in ['M1','M1B']:
        tr['her_input']=tr.her_voltage_V; te['her_input']=te.her_voltage_V
    elif variant in ['source_mean','batch_mean']:
        # Equal preparation weights, using training HER only, never target outcomes.
        a=tr.groupby(['prep_group','source','batch_id'],as_index=False).her_voltage_V.mean()
        overall=a.her_voltage_V.mean()
        source=a.groupby('source').her_voltage_V.mean()
        batch=a.groupby('batch_id').her_voltage_V.mean()
        for z in [tr,te]:
            v=z.source.map(source).fillna(overall)
            if variant=='batch_mean': v=z.batch_id.map(batch).fillna(v)
            z['her_input']=v
    num,cat=base_features(tr,budget,variant)
    return tr[num+cat],te[num+cat],num,cat
def weights(d):
    w=1/d.groupby('prep_group').prep_group.transform('size').to_numpy(dtype=float)
    return w/w.mean()
def make_model(algorithm,params,num,cat):
    prep=ColumnTransformer([
        ('numeric',Pipeline([('imputer',SimpleImputer(strategy='median',keep_empty_features=True)),('scale',StandardScaler())]),num),
        ('categorical',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cat)],remainder='drop')
    if algorithm=='ridge': estimator=Ridge(**params)
    elif algorithm=='random_forest':
        estimator=RandomForestRegressor(n_estimators=128,max_features=1.0,random_state=SEED,n_jobs=1,**params)
    else: raise ValueError(algorithm)
    return Pipeline([('prepare',prep),('model',estimator)])
def clip_pred(p,target): return np.clip(p,0,100) if target.endswith('_fe_pct') else p
def fit_predict(tr,te,budget,variant,algorithm,params,target):
    a,b,num,cat=prepare_inputs(tr,te,budget,variant)
    m=make_model(algorithm,params,num,cat)
    m.fit(a,tr[target].to_numpy(),model__sample_weight=weights(tr))
    return clip_pred(m.predict(b),target),m
def group_mae(d,y,p):
    return pd.DataFrame({'group':d.prep_group.to_numpy(),'e':np.abs(np.asarray(y)-np.asarray(p))}).groupby('group').e.mean().mean()
def cv_splits(d,group,n=5,seed=SEED):
    if group=='source': return list(LeaveOneGroupOut().split(d,groups=d[group]))
    return list(GroupKFold(n_splits=min(n,d[group].nunique()),shuffle=True,random_state=seed).split(d,groups=d[group]))
def grid(algorithm):
    if algorithm=='ridge': return [{'alpha':x} for x in [.1,1.,10.,100.]]
    return [{'max_depth':a,'min_samples_leaf':b} for a,b in itertools.product([3,6],[3,8])]
def tune(tr,budget,variant,algorithm,target,inner_group):
    spl=cv_splits(tr,inner_group,n=4,seed=SEED+1)
    scores=[]
    for params in grid(algorithm):
        preds=np.full(len(tr),np.nan)
        for ti,vi in spl:
            a=tr.iloc[ti]; b=tr.iloc[vi]
            assert not set(a[inner_group])&set(b[inner_group])
            assert not set(a.prep_group)&set(b.prep_group)
            preds[vi],_=fit_predict(a,b,budget,variant,algorithm,params,target)
        assert np.isfinite(preds).all()
        scores.append({'params':params,'inner_group_mae':float(group_mae(tr,tr[target],preds))})
    return min(scores,key=lambda x:x['inner_group_mae'])['params'],scores
def paired_stats(predictions,n_boot=10000,cluster='element_group'):
    d=predictions.copy()
    d['e0']=abs(d.y-d.M0); d['e1']=abs(d.y-d.M1)
    g=d.groupby(['prep_group',cluster],as_index=False)[['e0','e1']].mean()
    cl=g.groupby(cluster).agg(e0=('e0','sum'),e1=('e1','sum'),n=('prep_group','size'))
    rng=np.random.default_rng(SEED+90)
    idx=rng.integers(0,len(cl),size=(n_boot,len(cl)))
    den=cl.n.to_numpy()[idx].sum(axis=1)
    delta=(cl.e0.to_numpy()[idx].sum(axis=1)-cl.e1.to_numpy()[idx].sum(axis=1))/den
    m0=g.e0.mean(); m1=g.e1.mean()
    ci=np.quantile(delta,[.025,.975])
    return {'electrodes':len(d),'prep_groups':len(g),'element_groups':len(cl),'M0_MAE':float(m0),'M1_MAE':float(m1),'gain_pp_or_V':float(m0-m1),'relative_gain_pct':float(100*(m0-m1)/m0),'gain_CI95_low':float(ci[0]),'gain_CI95_high':float(ci[1]),'bootstrap_positive_fraction_descriptive':float((delta>0).mean())}
def screening(d,col):
    g=d.groupby('prep_group',as_index=False).agg(y=('y','mean'),p=(col,'mean'))
    k=int(np.ceil(.2*len(g)))
    truth=set(g.sort_values(['y','prep_group'],ascending=[False,True]).head(k).prep_group)
    pred=set(g.sort_values(['p','prep_group'],ascending=[False,True]).head(k).prep_group)
    hits=len(truth&pred)
    return {'candidate_groups':len(g),'top_k':k,'hits':hits,'recall':hits/k,'enrichment':(hits/k)/(k/len(g))}
def run_cv(d,budget='nominal',algorithm='random_forest',group='element_group',target='co_fe_pct',variants=('M0','M1'),tag=None):
    d=subset(d,budget).dropna(subset=[target]).reset_index(drop=True)
    tag=tag or f'{budget}_{algorithm}_{group}_{target}'
    dest=RESULTS/tag; dest.mkdir(parents=True,exist_ok=True)
    config={'budget':budget,'algorithm':algorithm,'group':group,'target':target,'variants':list(variants),'ids':d.electrode_id.tolist()}
    signature=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    if (dest/'complete.json').exists():
        old=json.loads((dest/'complete.json').read_text())
        assert old['signature']==signature,'Task config changed'
        print('CACHED '+tag,flush=True)
        return pd.read_csv(dest/'predictions.csv')
    outer=cv_splits(d,group)
    out=d[['electrode_id','prep_group','element_group','batch_id','source']].copy(); out['y']=d[target]; out['fold']=-1
    tuning=[]; memberships=[]
    for fold,(ti,vi) in enumerate(outer):
        tr=d.iloc[ti]; te=d.iloc[vi]
        assert not set(tr.prep_group)&set(te.prep_group)
        assert not set(tr[group])&set(te[group])
        out.loc[vi,'fold']=fold
        for eid in te.electrode_id: memberships.append({'electrode_id':eid,'fold':fold})
        inner='element_group' if group in ['element_group','source'] else group
        for var in variants:
            best,scores=tune(tr,budget,var,algorithm,target,inner)
            pred,model=fit_predict(tr,te,budget,var,algorithm,best,target)
            out.loc[vi,var]=pred
            tuning.append({'fold':fold,'variant':var,'best':best,'scores':scores,'train_groups':tr[group].nunique(),'validation_groups':te[group].nunique()})
            joblib.dump({'model':model,'train_ids':tr.electrode_id.tolist(),'validation_ids':te.electrode_id.tolist(),'params':best,'budget':budget,'variant':var,'target':target},dest/f'fold{fold}_{var}.joblib')
        print(json.dumps({'task':tag,'fold':fold+1,'folds':len(outer),'status':'completed'},ensure_ascii=False),flush=True)
    assert out[list(variants)].notna().all().all()
    out.to_csv(dest/'predictions.csv',index=False)
    pd.DataFrame(memberships).to_csv(dest/'fold_assignments.csv',index=False)
    jsave(tuning,dest/'tuning.json'); jsave(config,dest/'config.json')
    if 'M0' in variants and 'M1' in variants:
        summary=paired_stats(out)
        summary['screening']={v:screening(out,v) for v in ['M0','M1']}
        jsave(summary,dest/'summary.json')
        print(json.dumps({'task':tag,'summary':summary},ensure_ascii=False),flush=True)
    jsave({'signature':signature,'finished_unix':time.time()},dest/'complete.json')
    return out
