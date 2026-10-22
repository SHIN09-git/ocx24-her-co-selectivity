from analysis_core import *
raw=pd.read_csv(ROOT/'01_sources/ocx/data/experimental_data/ExpDataDump_241113.csv')
raw['raw_csv_line']=np.arange(len(raw))+2
n=raw[pd.to_numeric(raw.rep,errors='coerce').notna()].copy()
c=n[n.reaction=='CO2R'].copy()
for col in ['fe_co','fe_h2','fe_ch4','fe_c2h4']:
    c[col+'_outside_0_100']=~c[col].between(0,100)
c['gas_total_above_100']=c.fe_gas_total>100; c['liquid_residual_below_0']=c.fe_liquid<0
flagcols=[c for c in c if c.endswith('_outside_0_100')]+['gas_total_above_100','liquid_residual_below_0']
c[['sample id','current density','raw_csv_line']+flagcols].to_csv(ROOT/'02_data/quality_flags_all_currents.csv',index=False)
coverage=c.groupby('sample id').voltage.agg(['min','max','count'])
# Voltage signs are negative in some sources; compute on absolute magnitude explicitly.
c['absolute_voltage']=c.voltage.abs(); coverage=c.groupby('sample id').absolute_voltage.agg(['min','max','count'])
valid=coverage[coverage['count']==5]
covered=(valid['min']<=3.3)&(valid['max']>=3.3)
p=pd.read_csv(ROOT/'01_sources/ocx/data/processed_data/CO2R_40_70_all.csv')
co_theory=p['current density'].abs()/1000*p.fe_co/100/(2*96485)
eth_theory=p['current density'].abs()/1000*p.fe_c2h4/100/(12*96485)
mask=(p.fe_co>1e-6)&(p.fe_c2h4>1e-6)&(p.CO_pr>0)&(p.C2H4_pr>0)
rat=(p.C2H4_pr[mask]/eth_theory[mask])/(p.CO_pr[mask]/co_theory[mask])
audit={'complete_voltage_sequences':len(valid),'3_3V_inside_native_voltage_range':int(covered.sum()),'processed_rows':len(p),'processed_loco_groups':int(p.loco_cv.nunique()),'quality_counts':{x:int(c[x].sum()) for x in flagcols},'ethylene_processed_normalized_ratio':{'eligible_rows':int(mask.sum()),'min':float(rat.min()),'median':float(rat.median()),'max':float(rat.max()),'explanation':'(listed ethylene rate / 12-electron rate at displayed FE,j) divided by analogous CO ratio; common interpolation factor cancels'}}
jsave(audit,RESULTS/'supplemental_data_audit.json')
print(json.dumps(audit,ensure_ascii=False,indent=2))
