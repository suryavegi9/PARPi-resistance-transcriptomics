"""Publication assembly only. Reads frozen CSVs; no expression analysis or enrichment."""
from pathlib import Path
import pandas as pd
import numpy as np
import json, hashlib
import argparse
parser=argparse.ArgumentParser(description="Assemble tables from frozen outputs only; no DE or GSEA.")
parser.add_argument('--output',required=True)
args=parser.parse_args()
P=Path(__file__).resolve().parents[2]
O=Path(args.output).resolve()
if O.exists() or O==P or P in O.parents: parser.error('Use a new directory outside the repository.')
if not O.parent.is_dir(): parser.error('Output parent must exist.')
O.mkdir()
sources=[]
def read(rel):
 p=P/rel
 sources.append({'source':rel,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 d=pd.read_csv(p)
 if 'pathway' not in d and 'ID' in d: d=d.rename(columns={'ID':'pathway'})
 if 'pathway' in d:
  assert len(d)==50 and d.pathway.nunique()==50,rel
  d=d.set_index('pathway').sort_index()
 return d
def unique(name):
 hits=list((P/'results').rglob(name)); assert len(hits)==1,(name,hits)
 return read(str(hits[0].relative_to(P)))
old=read('results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv')
a=unique('D7A_M468_Hallmark_GSEA.csv'); b=unique('D7B_SUM1315_Hallmark_GSEA.csv'); d8=unique('D8_primary_Hallmark_GSEA.csv'); d5=unique('D5_H3_P3_Hallmark_master_comparison.csv')
f8='results/FINAL_D1_D8_integration/Figure8_FINAL/'
acc=read(f8+'Figure8_accession_source_table.csv'); mod=read(f8+'Figure8_model_level_source_table.csv'); tr=read(f8+'Figure8_treatment_state_source_table.csv')
for x in [a,b,d8,d5,acc,mod,tr]: assert old.index.equals(x.index)
same=np.sign(a.NES)==np.sign(b.NES)
d7nes=((a.NES+b.NES)/2).where(same)
d7dir=np.sign(d7nes).fillna(0).astype(int)
d7sig=same & (a['p.adjust']<.05) & (b['p.adjust']<.05)
def classify(up,dn,su,sd,dt,st,specific):
 return np.select([(su>=1)&(sd>=1),((up>=dt)&(su>=st)&(sd==0))|((dn>=dt)&(sd>=st)&(su==0)),((up>=dt)&(sd==0))|((dn>=dt)&(su==0)),su+sd==1],['Context-dependent','Strong recurrent','Directional',specific],default='Inconsistent')
S1=pd.DataFrame(index=old.index)
S1['classification']=acc.final_class
for c in ['n_UP','n_DOWN','n_sig_UP','n_sig_DOWN']:
 # Stored source uses final_ prefixes; do not infer summary values from labels.
 S1[c]=acc['final_'+c]
dirs=pd.DataFrame({f'D{i}':old[f'D{i}_direction'] for i in range(1,7)})
sigs=pd.DataFrame({f'D{i}':old[f'D{i}_sig_direction'] for i in range(1,7)})
dirs['D7']=d7dir; dirs['D8']=np.sign(d8.NES).astype(int)
sigs['D7']=d7dir*d7sig; sigs['D8']=dirs.D8*(d8['p.adjust']<.05)
for key,calc in [('n_UP',(dirs>0).sum(axis=1)),('n_DOWN',(dirs<0).sum(axis=1)),('n_sig_UP',(sigs>0).sum(axis=1)),('n_sig_DOWN',(sigs<0).sum(axis=1))]: assert (S1[key]==calc).all(),key
assert (classify(S1.n_UP,S1.n_DOWN,S1.n_sig_UP,S1.n_sig_DOWN,6,4,'Dataset-specific')==S1.classification).all()
S1['n_no_formal_direction']=(dirs==0).sum(axis=1)
for i in range(1,9):
 name=f'D{i}'
 S1[name+'_NES']=old[name+'_NES'] if i<=6 else d7nes if i==7 else d8.NES
 S1[name+'_direction']=dirs[name]
 S1[name+'_sig_direction']=sigs[name].astype(int)
S1['D7_models_same_direction']=same
for name,x in [('D7A_M468',a),('D7B_SUM1315',b),('D8_primary',d8)]:
 S1[name+'_NES']=x.NES; S1[name+'_padj']=x['p.adjust']
S1['D1_D6_historical_classification']=acc.old_class
S2=pd.DataFrame(index=old.index)
S2['classification']=mod.model12_class
for c in ['n_UP','n_DOWN','n_sig_UP','n_sig_DOWN']: S2[c]=mod['model12_'+c]
units={'D1_def':'D1_def','D1_prof':'D1_prof','D2':'D2','D3':'D3','D4':'D4','HCC_H1':'H1','HCC_H2':'H2','PSN_P1':'P1','D6':'D6'}
NES=pd.DataFrame(index=old.index); SIG=pd.DataFrame(index=old.index)
for unit,prefix in units.items():
 NES[unit]=old[prefix+'_NES']; SIG[unit]=old[prefix+'_padj']<.05
 if unit=='PSN_P1': SIG[unit]&=(old.P1_sens_padj<.05)&(np.sign(old.P1_NES)==np.sign(old.P1_sens_NES))
 S2[unit+'_NES']=NES[unit]; S2[unit+'_padj']=old[prefix+'_padj']; S2[unit+'_significant_vote']=SIG[unit]
for unit,x in [('D7_M468',a),('D7_SUM1315',b),('D8_PH039',d8)]:
 NES[unit]=x.NES; SIG[unit]=x['p.adjust']<.05
 S2[unit+'_NES']=NES[unit]; S2[unit+'_padj']=x['p.adjust']; S2[unit+'_significant_vote']=SIG[unit]
for key,calc in [('n_UP',(NES>0).sum(axis=1)),('n_DOWN',(NES<0).sum(axis=1)),('n_sig_UP',((NES>0)&SIG).sum(axis=1)),('n_sig_DOWN',((NES<0)&SIG).sum(axis=1))]: assert (S2[key]==calc).all(),key
assert (classify(S2.n_UP,S2.n_DOWN,S2.n_sig_UP,S2.n_sig_DOWN,8,6,'Unit-specific')==S2.classification).all()
S2['PSN_P1_sensitivity_NES']=old.P1_sens_NES; S2['PSN_P1_sensitivity_padj']=old.P1_sens_padj
S2['D6_sensitivity_NES']=old.D6_sens_NES; S2['D6_sensitivity_padj']=old.D6_sens_padj
S2['D6_robust_significant_context_only']=old.D6_robust_significant
S2['D1_D6_historical_model9_classification']=mod.model9_class
off=old[['D1_NES','D4_NES','D6_NES']]; on=old[['D2_NES','D3_NES']].copy(); on['D7_NES']=d7nes
offmean=off.mean(axis=1); onmean=on.mean(axis=1); delta=onmean-offmean; direct=d5[['H3_NES','P3_NES']].mean(axis=1)
agree=np.sign(delta)==np.sign(direct)
clean=same & (((delta>0)&(on.min(axis=1)>off.max(axis=1)))|((delta<0)&(on.max(axis=1)<off.min(axis=1))))
robust=d5.robust_shared_exposure_candidate
classes=np.select([robust&agree&clean,robust&agree&~clean,robust&~agree,~robust&clean],['Strong convergent','Convergent','D5-specific exposure','Observational-only'],default='Weak/mixed')
S3=pd.DataFrame({'classification':classes,'OFF_mean_NES':offmean,'ON_mean_NES':onmean,'observational_delta':delta,'D5_exposure_mean':direct,'robust_D5':robust,'same_exposure_direction':agree,'complete_ON_OFF_separation':clean,'ON_accessions_contributing':on.notna().sum(axis=1),'D7_models_same_direction':same},index=old.index)
for c in ['OFF_mean_NES','ON_mean_NES','observational_delta','D5_exposure_mean']: assert np.allclose(S3[c],tr[c],rtol=1e-13,atol=1e-13),c
for c in ['robust_D5','same_exposure_direction']: assert (S3[c]==tr[c]).all(),c
for c in off: S3[c]=off[c]
for c in on: S3[c]=on[c]
for c in ['H3_NES','H3_padj','P3_NES','P3_padj','P3_sensitivity_NES','P3_sensitivity_padj']: S3[c]=d5[c]
S3['OFF_min_NES']=off.min(axis=1); S3['OFF_max_NES']=off.max(axis=1)
S3['ON_min_NES']=on.min(axis=1); S3['ON_max_NES']=on.max(axis=1)
levels=['Strong recurrent','Directional','Context-dependent','Dataset-specific','Inconsistent']
counts1=S1.classification.value_counts().reindex(levels,fill_value=0).tolist()
counts2=S2.classification.value_counts().reindex([x.replace('Dataset-specific','Unit-specific') for x in levels],fill_value=0).tolist()
counts3=S3.classification.value_counts().reindex(['Strong convergent','Convergent','D5-specific exposure','Observational-only','Weak/mixed'],fill_value=0).tolist()
assert counts1==[0,4,24,12,10]; assert counts2==[0,5,35,5,5]; assert counts3==[0,9,9,2,30]
rho=delta.rank(method='average').corr(direct.rank(method='average')); assert abs(rho-(-.2557022809123649))<1e-14; assert agree.sum()==24
payload={}
for name,df in [('Table_S1',S1),('Table_S2',S2),('Table_S3',S3)]:
 df=df.reset_index(); df.to_csv(O/(name+'.csv'),index=False,na_rep='NA')
 payload[name]={'columns':list(df.columns),'rows':json.loads(df.to_json(orient='values',double_precision=15))}
 # JSON double_precision caps decimal digits; use standard json for full float precision instead.
 payload[name]['rows']=[[('NA' if pd.isna(v) else v.item() if isinstance(v,np.generic) else v) for v in row] for row in df.itertuples(index=False,name=None)]
pd.DataFrame(sources).to_csv(O/'Source_manifest.csv',index=False)
check={'accession_counts':counts1,'model_counts':counts2,'treatment_counts':counts3,'treatment_rho':rho,'treatment_direction_agreement':int(agree.sum()),'D7_concordant_pathways':int(same.sum()),'D7_significant_pathways':int(d7sig.sum()),'D8_significant_pathways':int((d8['p.adjust']<.05).sum()),'rows_each':50,'pathway_checks':'All 50 rows matched to saved Figure 8 count/classification/treatment fields; all per-unit values copied from locked sources.','DE_GSEA_executed':False}
(O/'Validation_report.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
(O/'tables.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False),encoding='utf-8')
print(json.dumps(check,indent=2)); print('Treatment classes:'); print(S3.groupby('classification').apply(lambda x: ', '.join(x.index),include_groups=False).to_string())
