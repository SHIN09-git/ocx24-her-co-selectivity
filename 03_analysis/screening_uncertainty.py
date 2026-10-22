from analysis_core import *
def metrics(y,p,ids):
    n=len(y); k=int(np.ceil(.2*n))
    a=np.lexsort((ids,-y))[:k]; b=np.lexsort((ids,-p))[:k]
    hit=len(set(a)&set(b)); return hit/k,(hit/k)/(k/n)
def compute(z,label):
    g=z.groupby(['prep_group','element_group'],as_index=False)[['y','M0','M1']].mean()
    rng=np.random.default_rng(SEED+311)
    clusters=[x.index.to_numpy() for _,x in g.groupby('element_group')]
    samples={c:[] for c in ['M0','M1']}
    for _ in range(2000):
        take=np.concatenate([clusters[i] for i in rng.integers(0,len(clusters),len(clusters))])
        b=g.loc[take].reset_index(drop=True)
        for model in samples: samples[model].append(metrics(b.y.to_numpy(),b[model].to_numpy(),np.arange(len(b))))
    out=[]
    for model in samples:
        v=np.array(samples[model]); actual=screening(z,model)
        out.append({'partition':label,'model':model,**actual,'recall_CI95_low':float(np.quantile(v[:,0],.025)),'recall_CI95_high':float(np.quantile(v[:,0],.975)),'enrichment_CI95_low':float(np.quantile(v[:,1],.025)),'enrichment_CI95_high':float(np.quantile(v[:,1],.975))})
    return out
dev=pd.read_csv(RESULTS/'nominal_random_forest_element_group_co_fe_pct/predictions.csv')
ho=pd.read_csv(RESULTS/'holdout_evaluated_predictions.csv')
ho=ho[(ho.budget=='nominal')&(ho.algorithm=='random_forest')&(ho.target=='co_fe_pct')]
rows=compute(dev,'development_crossfit')+compute(ho,'holdout')
pd.DataFrame(rows).to_csv(RESULTS/'screening_uncertainty.csv',index=False)
# Recorded-repeat bootstrap only. Singletons retain their observed value.
rng=np.random.default_rng(SEED+412); groups=list(dev.groupby('prep_group',sort=True))
ids=np.array([g for g,_ in groups]); k=int(np.ceil(.2*len(ids)))
baseline={c:np.array([z[c].mean() for _,z in groups]) for c in ['y','M0','M1']}
baseline_sets={c:set(np.lexsort((ids,-v))[:k]) for c,v in baseline.items()}
stability={c:[] for c in baseline}
for _ in range(2000):
    vals={c:[] for c in baseline}
    for _,z in groups:
        selected=rng.integers(0,len(z),len(z))
        for c in vals: vals[c].append(z[c].to_numpy()[selected].mean())
    for c,v in vals.items():
        chosen=set(np.lexsort((ids,-np.array(v)))[:k]); old=baseline_sets[c]
        stability[c].append(len(chosen&old)/len(chosen|old))
jsave({'bootstrap_iterations':2000,'candidate_groups':len(ids),'groups_with_repeats':int(sum(len(z)>1 for _,z in groups)),'top_k':k,'jaccard_mean':{c:float(np.mean(v)) for c,v in stability.items()},'interpretation':'Conditional stability using observed repeat electrodes only; singleton uncertainty and model-training variation are not estimated'},RESULTS/'ranking_stability.json')
print(json.dumps({'screening':rows,'ranking_stability':{c:float(np.mean(v)) for c,v in stability.items()}},ensure_ascii=False))
