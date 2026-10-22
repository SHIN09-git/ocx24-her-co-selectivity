from analysis_core import *
if __name__=='__main__':
    d=load_development()
    for group in ['element_group','prep_group']:
        for budget in ['nominal','xrf']:
            for alg in ['random_forest','ridge']:
                run_cv(d,budget=budget,algorithm=alg,group=group)
    rows=[]
    for f in RESULTS.glob('*/summary.json'):
        a=json.loads(f.read_text()); a.pop('screening',None); rows.append({'task':f.parent.name,**a})
    pd.DataFrame(rows).to_csv(RESULTS/'primary_comparisons.csv',index=False)
