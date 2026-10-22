"""Export the numerical table behind the archived family-error figure."""
from analysis_core import *
p=pd.read_csv(RESULTS/'nominal_random_forest_element_group_co_fe_pct/predictions.csv')
p['e0']=abs(p.y-p.M0); p['e1']=abs(p.y-p.M1)
pg=p.groupby(['prep_group','element_group'],as_index=False)[['e0','e1']].mean()
fam=pg.groupby('element_group').agg(e0=('e0','mean'),e1=('e1','mean'),n=('prep_group','size'))
fam['gain']=fam.e0-fam.e1
fam.sort_values('gain').to_csv(RESULTS/'primary_family_errors.csv')
