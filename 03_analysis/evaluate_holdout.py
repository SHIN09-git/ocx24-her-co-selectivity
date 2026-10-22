"""Replay the original prediction-before-label ordering; this is not a new blind test."""
from analysis_core import *
import datetime

if __name__=='__main__':
    verify_lock()
    if (ROOT/'00_admin/holdout_unsealed.json').exists():
        raise SystemExit('Already unsealed. Preserve original predictions; do not retune on these outcomes.')
    d=load_development()
    f=pd.read_csv(ROOT/'02_data/features.csv')
    test=f[f.partition=='holdout'].reset_index(drop=True)
    specs=[]
    for budget in ['nominal','xrf']:
        for alg in ['random_forest','ridge']:
            specs.append((budget,alg,'co_fe_pct'))
    specs.extend([('nominal','random_forest','h2_fe_pct'),('nominal','random_forest','co2_voltage_V')])
    rows=[]; tuning=[]
    for budget,alg,target in specs:
        tr=subset(d,budget); te=subset(test,budget)
        result=te[['electrode_id','prep_group','element_group','source','batch_id']].copy()
        result['budget']=budget; result['algorithm']=alg; result['target']=target
        for var in ['M0','M1']:
            best,scores=tune(tr,budget,var,alg,target,'element_group')
            pred,model=fit_predict(tr,te,budget,var,alg,best,target)
            result[var]=pred
            tuning.append({'budget':budget,'algorithm':alg,'target':target,'variant':var,'best':best,'scores':scores})
            joblib.dump(model,RESULTS/f'final_{budget}_{alg}_{target}_{var}.joblib')
        rows.append(result)
        print(json.dumps({'holdout_input_only_predictions_completed':[budget,alg,target]},ensure_ascii=False),flush=True)
    pending=pd.concat(rows,ignore_index=True)
    pred_path=RESULTS/'holdout_predictions_before_unseal.csv'; pending.to_csv(pred_path,index=False)
    jsave(tuning,RESULTS/'holdout_model_selection_development_only.json')
    stamp={'scope':'Public reproduction replay, not a new sealed evaluation','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'predictions_sha256':hashlib.sha256(pred_path.read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((ROOT/'00_admin/preregistered_analysis_plan.json').read_bytes()).hexdigest(),'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'03_analysis').glob('*.py')},'decision_before_unseal':'Close positive-screening route: development main comparisons do not meet 3 pp and 10% gain. Holdout is a final prespecified evaluation, not a chance to retune.'}
    jsave(stamp,ROOT/'00_admin/predictions_frozen_before_unseal.json')
    # Replay the historical ordering; the public outcomes are already known.
    y=pd.read_csv(ROOT/'02_data/holdout_targets.csv')
    records=[]; evaluated=[]
    for (budget,alg,target),a in pending.groupby(['budget','algorithm','target'],sort=False):
        z=a.merge(y[['electrode_id',target]],on='electrode_id',validate='one_to_one').rename(columns={target:'y'})
        stats=paired_stats(z); stats.update(budget=budget,algorithm=alg,target=target)
        if target=='co_fe_pct': stats['screening']={v:screening(z,v) for v in ['M0','M1']}
        records.append(stats); evaluated.append(z)
    pd.concat(evaluated,ignore_index=True).to_csv(RESULTS/'holdout_evaluated_predictions.csv',index=False)
    jsave(records,RESULTS/'holdout_results.json')
    coverage={'electrodes':len(test),'prep_groups':test.prep_group.nunique(),'element_groups':test.element_group.nunique(),'electrode_sources':test.source.value_counts().to_dict(),'prep_sources':test.groupby('source').prep_group.nunique().to_dict(),'note':'Input-selected families are uniform by element-set identity, not balanced by electrode count; one held-out source may dominate. No redraw.'}
    jsave(coverage,RESULTS/'holdout_coverage.json')
    jsave({'scope':'Public reproduction replay, not a new sealed evaluation','unsealed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'prediction_sha256':stamp['predictions_sha256'],'protocol_unchanged':True,'no_further_model_selection_allowed':True},ROOT/'00_admin/holdout_unsealed.json')
    print(json.dumps({'holdout_results':records,'coverage':coverage},ensure_ascii=False,indent=2),flush=True)
