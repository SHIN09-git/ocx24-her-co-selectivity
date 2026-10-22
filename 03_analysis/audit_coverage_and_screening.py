"""Complete prespecified reporting from stored predictions; fit no new models."""
from analysis_core import *

verify_lock()
F=pd.read_csv(ROOT/'02_data/features.csv').set_index('electrode_id')
num=[c for c in F if re.fullmatch(r'nominal_[A-Z][a-z]?',c)]
def coverage(tr,te,task,fold):
    a={c[8:] for c in num if (tr[c]>0).any()}; b={c[8:] for c in num if (te[c]>0).any()}
    unseen=sorted(b-a)
    affected=(te[['nominal_'+e for e in unseen]]>0).any(axis=1) if unseen else pd.Series(False,index=te.index)
    hmin,hmax=tr.her_voltage_V.min(),tr.her_voltage_V.max()
    return {'task':task,'fold':fold,'training_electrodes':len(tr),'validation_electrodes':len(te),
            'training_prep_groups':tr.prep_group.nunique(),'validation_prep_groups':te.prep_group.nunique(),
            'training_element_sets':tr.element_group.nunique(),'validation_element_sets':te.element_group.nunique(),
            'shared_element_sets':len(set(tr.element_group)&set(te.element_group)),
            'training_elements':' '.join(sorted(a)),'validation_elements':' '.join(sorted(b)),
            'unseen_nominal_elements':' '.join(unseen),'unseen_element_count':len(unseen),
            'affected_electrodes':int(affected.sum()),'affected_prep_groups':te.loc[affected,'prep_group'].nunique(),
            'unseen_sources':' '.join(sorted(set(te.source)-set(tr.source))),
            'unseen_batch_count':len(set(te.batch_id)-set(tr.batch_id)),
            'training_HER_min_V':hmin,'training_HER_max_V':hmax,
            'validation_HER_outside_training_range':int(((te.her_voltage_V<hmin)|(te.her_voltage_V>hmax)).sum())}

cov=[]; scores=[]; cluster_stats=[]
for p in sorted(RESULTS.glob('*/predictions.csv')):
    cfg=json.loads((p.parent/'config.json').read_text())
    z=pd.read_csv(p)
    first=cfg['variants'][0]
    for fold,t in z.groupby('fold'):
        rec=joblib.load(p.parent/f'fold{fold}_{first}.joblib')
        tr=F.loc[rec['train_ids']]; te=F.loc[rec['validation_ids']]
        assert set(te.index)==set(t.electrode_id)
        assert not set(tr.prep_group)&set(te.prep_group)
        cov.append(coverage(tr,te,p.parent.name,int(fold)))
        if cfg['target']=='co_fe_pct':
            for v in cfg['variants']:
                scores.append({'task':p.parent.name,'fold':int(fold),'model':v,**screening(t,v)})
    if cfg['group']=='batch_id' and set(['M0','M1'])<=set(z):
        # Original element-cluster intervals are retained; this is the requested
        # batch-level stress-test uncertainty, with explicit cluster labels.
        s=paired_stats(z,cluster='batch_id'); s['cluster_count']=s.pop('element_groups')
        cluster_stats.append({'task':p.parent.name,'cluster':'batch_id',**s})
pd.DataFrame(cov).to_csv(RESULTS/'validation_input_coverage.csv',index=False)
pd.DataFrame(scores).to_csv(RESULTS/'fold_screening_metrics.csv',index=False)
pd.DataFrame(cluster_stats).to_csv(RESULTS/'batch_clustered_intervals.csv',index=False)
tr=F[F.partition=='development']; te=F[F.partition=='holdout']
jsave(coverage(tr,te,'final_holdout','holdout'),RESULTS/'holdout_element_coverage.json')

# Per-family and per-source errors from exactly the existing predictions.
ho=pd.read_csv(RESULTS/'holdout_evaluated_predictions.csv')
families=[]
for (budget,alg,target),q in ho.groupby(['budget','algorithm','target']):
    q=q.copy(); q['e0']=abs(q.y-q.M0); q['e1']=abs(q.y-q.M1)
    g=q.groupby(['prep_group','element_group'],as_index=False)[['e0','e1']].mean()
    t=g.groupby('element_group').agg(prep_groups=('prep_group','size'),M0_MAE=('e0','mean'),M1_MAE=('e1','mean')).reset_index()
    t['gain']=t.M0_MAE-t.M1_MAE; t['budget']=budget; t['algorithm']=alg; t['target']=target; families.append(t)
pd.concat(families,ignore_index=True).to_csv(RESULTS/'holdout_family_errors.csv',index=False)
base=pd.read_csv(RESULTS/'nominal_random_forest_element_group_co_fe_pct/predictions.csv')
base['e0']=abs(base.y-base.M0); base['e1']=abs(base.y-base.M1)
g=base.groupby(['prep_group','element_group','source'],as_index=False)[['e0','e1']].mean()
s=g.groupby('source').agg(prep_groups=('prep_group','size'),M0_MAE=('e0','mean'),M1_MAE=('e1','mean')).reset_index()
s['gain']=s.M0_MAE-s.M1_MAE; s.to_csv(RESULTS/'primary_crossfit_source_errors.csv',index=False)
influence=[]
for fam in sorted(g.element_group.unique()):
    remain=g[g.element_group!=fam]
    influence.append({'omitted_evaluation_family':fam,'retained_prep_groups':len(remain),'M0_MAE':remain.e0.mean(),'M1_MAE':remain.e1.mean(),'gain':remain.e0.mean()-remain.e1.mean(),'refit':False})
pd.DataFrame(influence).to_csv(RESULTS/'leave_one_family_evaluation_influence.csv',index=False)

cv=pd.DataFrame(cov)
focus=cv[cv.task.isin(['nominal_random_forest_element_group_co_fe_pct','nominal_random_forest_batch_id_co_fe_pct','nominal_random_forest_source_co_fe_pct'])]
main_scores=pd.DataFrame(scores).query("task == 'nominal_random_forest_element_group_co_fe_pct'")
summary={'models_refit':0,'fold_coverage_rows':len(cov),'fold_screening_rows':len(scores),
         'main_and_stress_coverage':focus[['task','fold','unseen_nominal_elements','affected_electrodes','validation_electrodes']].to_dict('records'),
         'primary_fold_screening':main_scores.to_dict('records'),
         'batch_clustered_intervals':cluster_stats,'scope':'Coverage uses only stored inputs; all metrics use stored predictions. Supplemental reporting after unsealing; no new selection or independent validation.'}
jsave(summary,RESULTS/'reporting_completeness_audit.json')
print(json.dumps(summary,ensure_ascii=False,indent=2))
