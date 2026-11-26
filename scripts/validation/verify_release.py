"""Read-only integrity and locked-checkpoint verification; no biological analysis."""
from pathlib import Path
import argparse,csv,hashlib,json,math,re
import pandas as pd
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-release-manifest',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]; checks=[]
    def check(name,condition,detail):
        checks.append(dict(check=name,passed=bool(condition),detail=detail))
    def read(rel):return pd.read_csv(root/rel)
    fig='results/FINAL_D1_D8_integration/Figure8_FINAL/'
    a=read(fig+'Figure8_accession_source_table.csv').set_index('pathway').sort_index()
    m=read(fig+'Figure8_model_level_source_table.csv').set_index('pathway').sort_index()
    t=read(fig+'Figure8_treatment_state_source_table.csv').set_index('pathway').sort_index()
    levels=['Strong recurrent','Directional','Context-dependent','Dataset-specific','Inconsistent']
    check('eight accession classes',a.final_class.value_counts().reindex(levels,fill_value=0).tolist()==[0,4,24,12,10],'0/4/24/12/10')
    check('twelve unit classes',m.model12_class.value_counts().reindex([x.replace('Dataset','Unit') for x in levels],fill_value=0).tolist()==[0,5,35,5,5],'0/5/35/5/5')
    d7a=read('results/GSE165914/D7_analysis/D7A_M468_Hallmark_GSEA.csv').set_index('pathway').sort_index()
    d7b=read('results/GSE165914/D7_analysis/D7B_SUM1315_Hallmark_GSEA.csv').set_index('pathway').sort_index()
    same=np.sign(d7a.NES)==np.sign(d7b.NES)
    check('D7 direction and significance',same.sum()==33 and (same&(d7a['p.adjust']<.05)&(d7b['p.adjust']<.05)).sum()==6,'33/50; 6/50')
    d8=read('results/GSE163854/D8_analysis/D8_primary_Hallmark_GSEA.csv')
    robust=read('results/GSE163854/D8_analysis/D8_Hallmark_limma_vs_edgeR.csv')
    check('D8 primary and sensitivity',(d8['p.adjust']<.05).sum()==15 and robust.robust_significant.sum()==14,'15 significant; 14 robust')
    check('classification changes',(a.old_class!=a.final_class).sum()==18,'18/50')
    rho=t.observational_delta.rank().corr(t.D5_exposure_mean.rank())
    check('treatment checkpoints',abs(rho+.255702280912365)<1e-12 and t.same_exposure_direction.sum()==24 and t.robust_D5.sum()==18 and (t.same_exposure_direction&t.robust_D5).sum()==9,'rho=-0.255702280912365; 24/50; 18; 9/18')
    s=[read('supplementary/Table_S'+str(i)+'.csv').set_index('pathway').sort_index() for i in (1,2,3)]
    check('supplement dimensions and identifiers',all(len(x)==50 and x.index.is_unique and x.index.equals(a.index) for x in s),'50 unique matched Hallmarks each')
    check('S1/S2 final classification',s[0].classification.equals(a.final_class) and s[1].classification.equals(m.model12_class),'Exact final classes')
    classes=['Strong convergent','Convergent','D5-specific exposure','Observational-only','Weak/mixed']
    check('treatment classes',s[2].classification.value_counts().reindex(classes,fill_value=0).tolist()==[0,9,9,2,30],'0/9/9/2/30')
    hypothesis=read(fig+'Figure8_frozen_hypothesis_heatmap_source.csv').pathway.unique()
    check('frozen hypotheses',len(hypothesis)==10 and (a.loc[hypothesis,'final_class']=='Context-dependent').sum()==9,'9 of 10 context-dependent')
    check('glycolysis',a.loc['HALLMARK_GLYCOLYSIS','final_class']=='Inconsistent' and m.loc['HALLMARK_GLYCOLYSIS','model12_class']=='Context-dependent','Inconsistent accession; context-dependent model')
    d4=read('provenance/D4_reconstruction/20261125_132322/reconstructed_DE.csv')
    filt=read('provenance/D4_reconstruction/20261125_132322/reconstructed_filter.csv')
    check('D4 locked counts',len(filt)==15623 and len(d4)==14277 and (d4['adj.P.Val']<.05).sum()==1984 and ((d4['adj.P.Val']<.05)&(d4.logFC>=1)).sum()==588 and ((d4['adj.P.Val']<.05)&(d4.logFC<=-1)).sum()==453,'15623/14277/1984/588/453')
    d4report=(root/'provenance/D4_reconstruction/20261125_132322/D4_reconstruction_report.txt').read_text()
    check('D4 correlation evidence','0.39477637734058' in d4report,'Archived report; no DE rerun')
    for directory,file,n in [('D1_D6_reconstruction/20261125_133657','column_comparison.csv',49),('final_synthesis/20261125_165053','final_comparison.csv',46)]:
        x=read('provenance/'+directory+'/'+file)
        check(file,len(x)==n and x['pass'].all(),str(n)+' passing fields')
    for rel in ['provenance/SOURCE_TRANSFER_MANIFEST.csv','provenance/publication_tables/Source_manifest.csv']:
        rows=list(csv.DictReader((root/rel).open(encoding='utf-8')))
        if rel.endswith('publication_tables/Source_manifest.csv'):
            rows=[r for r in rows if r['source'].replace('\\','/').startswith('results/')]
        ok=all(hashlib.sha256((root/r.get('path',r.get('source','')).replace('\\','/')).read_bytes()).hexdigest()==r['sha256'] for r in rows)
        check(rel,ok,'Transferred files / frozen result sources; original historical script hashes are preserved separately')
    if not args.skip_release_manifest:
        rows=list(csv.DictReader((root/'RELEASE_MANIFEST.csv').open(encoding='utf-8')))
        check('release manifest',all((root/r['path']).stat().st_size==int(r['bytes']) and hashlib.sha256((root/r['path']).read_bytes()).hexdigest()==r['sha256'] for r in rows),str(len(rows))+' files')
        actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '.git' not in p.parts and '__pycache__' not in p.parts and p.name!='RELEASE_MANIFEST.csv'}
        check('manifest coverage',actual=={r['path'] for r in rows},'All release files except manifest itself')
    missing=[]
    for p in root.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-z]+:',target) or target.startswith('#'):continue
            target=target.split('#')[0]
            if target and not (p.parent/target).exists():missing.append(p.relative_to(root).as_posix()+': '+target)
    check('relative Markdown links',not missing,missing)
    print(json.dumps({'checks':checks,'passed':all(x['passed'] for x in checks),'no_DE_GSEA':True},indent=2))
    if not all(x['passed'] for x in checks):raise SystemExit(1)

if __name__=='__main__':main()
