from analysis_core import *

def permutation_diagnostics(d):
    tag='nominal_random_forest_element_group_co_fe_pct'
    base=pd.read_csv(RESULTS/tag/'predictions.csv')
    allrows=[]; overall=[]
    for kind in ['prep_group','batch_id']:
        errors={}; reference={}; moved={}
        for fold in sorted(base.fold.unique()):
            obj=joblib.load(RESULTS/tag/f'fold{fold}_M1.joblib')
            tr=d.set_index('electrode_id').loc[obj['train_ids']].reset_index()
            te=d.set_index('electrode_id').loc[obj['validation_ids']].reset_index()
            groups=list(te.groupby(kind,sort=True).indices.values())
            mask=te.groupby(kind).electrode_id.transform('size').to_numpy()>1
            _,x,_,_=prepare_inputs(tr,te,'nominal','M1')
            p=clip_pred(obj['model'].predict(x),'co_fe_pct')
            rng=np.random.default_rng(SEED+1000+int(fold)+(0 if kind=='prep_group' else 100))
            batches=[]
            for i in range(200):
                xx=x.copy()
                hv=xx.her_input.to_numpy().copy()
                for ids in groups: hv[ids]=rng.permutation(hv[ids])
                xx['her_input']=hv; batches.append(xx)
            # Predict in one batch to avoid repeated forest overhead.
            pp=clip_pred(obj['model'].predict(pd.concat(batches,ignore_index=True)),'co_fe_pct').reshape(200,len(te))
            errs=abs(pp-te.co_fe_pct.to_numpy()[None,:])
            for k,eid in enumerate(te.electrode_id):
                if mask[k]:
                    errors[eid]=errs[:,k]; reference[eid]=abs(p[k]-te.co_fe_pct.iloc[k]); moved[eid]=True
        ids=sorted(errors)
        er=np.array([errors[e] for e in ids]); ref=np.array([reference[e] for e in ids])
        delta=er.mean(axis=1)-ref
        rows=d.set_index('electrode_id').loc[ids][['prep_group','element_group','source']].reset_index()
        rows['correct_pair_error']=ref; rows['mean_shuffled_error']=er.mean(axis=1); rows['shuffled_minus_correct']=delta; rows['permutation_group']=kind
        allrows.append(rows)
        # Electrode-equal loss, resample whole element families.
        cl=rows.groupby('element_group').agg(total=('shuffled_minus_correct','sum'),n=('electrode_id','size'))
        rng=np.random.default_rng(SEED+212)
        idx=rng.integers(0,len(cl),size=(10000,len(cl)))
        vals=cl.total.to_numpy()[idx].sum(1)/cl.n.to_numpy()[idx].sum(1)
        overall.append({'permutation_group':kind,'electrodes_with_permutable_group':len(ids),'prep_groups':rows.prep_group.nunique(),'element_groups':len(cl),'correct_pair_electrode_MAE':float(ref.mean()),'mean_shuffled_electrode_MAE':float(er.mean()),'shuffled_minus_correct_MAE':float(delta.mean()),'CI95_low':float(np.quantile(vals,.025)),'CI95_high':float(np.quantile(vals,.975)),'permutation_iterations':200,'interpretation':'Validation-only reliance diagnostic; positive means correct pairing helps; not a causal or exact-null test'})
    pd.concat(allrows).to_csv(RESULTS/'permutation_electrode_errors.csv',index=False)
    jsave(overall,RESULTS/'permutation_summary.json')
    print(json.dumps({'permutations':overall},ensure_ascii=False),flush=True)

def noise_and_units(d):
    rows=[]
    for g,z in d.groupby('prep_group'):
        if len(z)>1:
            co=z.co_fe_pct.to_numpy(); hv=z.her_voltage_V.to_numpy()
            diffs=[abs(co[i]-co[j]) for i in range(len(co)) for j in range(i)]
            loo=[abs(co[i]-np.delete(co,i).mean()) for i in range(len(co))]
            rows.append({'prep_group':g,'element_group':z.element_group.iloc[0],'source':z.source.iloc[0],'n':len(z),'co_FE_sd_pp':float(np.std(co,ddof=1)),'her_voltage_sd_V':float(np.std(hv,ddof=1)),'co_FE_pairwise_abs_difference_pp':float(np.mean(diffs)),'within_group_LOO_reference_MAE_pp':float(np.mean(loo))})
    r=pd.DataFrame(rows); r.to_csv(RESULTS/'replicate_noise.csv',index=False)
    jsave({'development_replicated_groups':len(r),'median_CO_SD_pp':r.co_FE_sd_pp.median(),'mean_CO_SD_pp':r.co_FE_sd_pp.mean(),'median_pairwise_CO_difference_pp':r.co_FE_pairwise_abs_difference_pp.median(),'mean_within_group_LOO_reference_MAE_pp':r.within_group_LOO_reference_MAE_pp.mean(),'median_HER_SD_V':r.her_voltage_sd_V.median(),'interpretation':'Repeat dispersion at fixed preparation ID; not a formal irreducible-error lower bound and not all independent synthesis variability'},RESULTS/'noise_summary.json')
    src=ROOT/'01_sources/ocx/core/data_handling.py'
    source=src.read_text()
    assert 'get_production_rate(row.fe_c2h4, 2,' in source
    raw=pd.read_csv(ROOT/'01_sources/ocx/data/experimental_data/ExpDataDump_241113.csv')
    q=raw[(raw.reaction=='CO2R') & pd.to_numeric(raw.rep,errors='coerce').notna()].copy()
    q['c2h4_rate_correct_mol_s_cm2']=abs(q['current density'])/1000*(q.fe_c2h4/100)/(12*96485)
    q['c2h4_rate_if_n2_mol_s_cm2']=abs(q['current density'])/1000*(q.fe_c2h4/100)/(2*96485)
    nonzero=q.fe_c2h4>0
    ratio=q.loc[nonzero,'c2h4_rate_if_n2_mol_s_cm2']/q.loc[nonzero,'c2h4_rate_correct_mol_s_cm2']
    assert np.allclose(ratio,6)
    # Does not read or use sealed CO FE. Unit correction is ancillary, not a model input.
    q[['sample id','current density','fe_c2h4','c2h4_rate_correct_mol_s_cm2','c2h4_rate_if_n2_mol_s_cm2']].to_csv(RESULTS/'ethylene_native_point_unit_audit.csv',index=False)
    jsave({'electrons_CO2_to_C2H4':12,'upstream_code_electrons':2,'nonzero_native_rows':int(nonzero.sum()),'ratio_min':float(ratio.min()),'ratio_max':float(ratio.max()),'scope':'Native-point rate comparison; FE unchanged, primary CO FE modeling unaffected; interpolation ratios from processed table are a different calculation'},RESULTS/'unit_audit.json')

if __name__=='__main__':
    d=load_development()
    noise_and_units(d)
    # All comparisons below were specified before the first training run.
    p=run_cv(d,variants=('M0B','M1B','source_mean','batch_mean'),tag='nominal_RF_batch_and_mean_controls')
    base=pd.read_csv(RESULTS/'nominal_random_forest_element_group_co_fe_pct/predictions.csv')
    combined=base.merge(p.drop(columns=['prep_group','element_group','batch_id','source','y','fold']),on='electrode_id',validate='one_to_one')
    combined.to_csv(RESULTS/'control_predictions.csv',index=False)
    cs={}
    for col in ['source_mean','batch_mean']:
        z=combined.drop(columns='M1').rename(columns={col:'M1'})
        cs['M0_vs_'+col]=paired_stats(z)
    z=combined.drop(columns=['M0','M1']).rename(columns={'M0B':'M0','M1B':'M1'})
    cs['M0B_vs_M1B']=paired_stats(z)
    jsave(cs,RESULTS/'batch_control_summary.json')
    permutation_diagnostics(d)
    for group in ['batch_id','source']:
        run_cv(d,group=group)
    for target in ['h2_fe_pct','co2_voltage_V']:
        run_cv(d,target=target)
    quality=(d.co_fe_pct.between(0,100)&d.h2_fe_pct.between(0,100)&d.gas_total_pct.le(100)&d.liquid_residual_pct.ge(0))
    run_cv(d[quality],tag='nominal_RF_quality_filtered')
    run_cv(d[d.in_official_clean],tag='nominal_RF_official_clean')
    # Largest family is chosen by preparation counts, independent of model outcomes.
    largest=d.groupby('element_group').prep_group.nunique().sort_values(ascending=False,kind='stable').index[0]
    run_cv(d[d.element_group!=largest],tag='nominal_RF_without_largest_family')
    jsave({'largest_family_by_preparation_count':largest},RESULTS/'family_removal.json')
    for src in sorted(d.source.unique()):
        run_cv(d[d.source==src],tag='nominal_RF_within_source_'+src)
    rows=[]
    for f in RESULTS.glob('*/summary.json'):
        a=json.loads(f.read_text()); a.pop('screening',None); rows.append({'task':f.parent.name,**a})
    pd.DataFrame(rows).to_csv(RESULTS/'all_development_comparisons.csv',index=False)
