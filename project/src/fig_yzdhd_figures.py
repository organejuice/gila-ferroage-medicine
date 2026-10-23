"""Render frozen YZDHD tables; no fit, new significance test or upstream rerun.

Style-only inheritance from the nature-figure Python defaults. Existing project
table loader/provenance patterns are retained; the panels use their actual units.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
CFG=json.loads((ROOT/'configs/fig_yzdhd_figures.json').read_text())
OUT=ROOT/CFG['output_dir']; OUT.mkdir(parents=True,exist_ok=True)
DATA=OUT/'source_data'; DATA.mkdir(exist_ok=True)
font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,
 'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.65,'savefig.facecolor':'white'})
SOURCES={}; OUTPUTS=[]; CAPTIONS={}
BLUE='#326B8B'; RED='#AA5C50'; GRAY='#B7BDC1'
width_mm=174

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(rel):
 p=ROOT/rel; SOURCES[rel]=sha(p)
 return pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t')
def save_data(d,name):d.to_csv(DATA/name,index=False,sep='\t')
def panel(ax,s):ax.text(-.14,1.06,s,transform=ax.transAxes,fontsize=11,fontweight='bold',va='top')
def save(fig,n):
 fig.canvas.draw()
 for ext in ['.pdf','.svg','.png','.tif']:
  p=OUT/f'Fig{n}{ext}'
  if ext=='.tif':fig.savefig(p,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
  else:fig.savefig(p,dpi=300)
  if ext=='.tif':
   with Image.open(p) as raw: rgb=raw.convert('RGB')
   rgb.save(p,compression='tiff_lzw',dpi=(600,600))
  OUTPUTS.append({'file':p.name,'sha256':sha(p)})
 plt.close(fig)

G1='results/archive/2026-09-20-disease-side-v3.1-g2-stop-historical/meta/random_effects_meta_summary.tsv'
ST='results/2026-10-13-manuscript/supplementary_tables/'
meta=read(G1); sig=read(ST+'Table_S1_G1_244_signature.tsv.gz')
assert len(meta)==19738 and len(sig)==244
save_data(meta,'Fig1a_all_meta_genes.tsv');save_data(sig,'Fig1b_244_signature.tsv')
fig,axes=plt.subplots(1,2,figsize=(width_mm/25.4,80/25.4),gridspec_kw={'width_ratios':[2.3,1]})
fig.subplots_adjust(left=.10,right=.98,bottom=.23,top=.89,wspace=.46)
ax=axes[0]; colors=np.where(meta.gene.isin(sig.gene),np.where(meta.meta_logFC>0,RED,BLUE),GRAY)
if not (meta.meta_FDR>0).all():raise ValueError('Non-positive FDR cannot be logarithmically displayed')
ax.scatter(meta.meta_logFC,-np.log10(meta.meta_FDR),s=3,c=colors,alpha=.65,rasterized=False)
ax.axhline(-np.log10(.05),color='#444444',lw=.7,ls='--');ax.axvline(0,color='#777777',lw=.5)
ax.set(xlabel='Meta-analysis log2 fold change (AD − control)',ylabel='−log10(FDR)');panel(ax,'a')
ax=axes[1];counts=[int((sig.meta_logFC>0).sum()),int((sig.meta_logFC<0).sum())]
assert counts==[128,116]
ax.bar(['Positive','Negative'],counts,color=[RED,BLUE],width=.65)
for i,c in enumerate(counts):ax.text(i,c+3,str(c),ha='center',fontsize=9)
ax.set(ylim=(0,150),ylabel='Genes in complete-region subset');panel(ax,'b')
save(fig,1)
CAPTIONS['Fig1']='Human microglial AD association and the direction-preserved subset. a All 19,738 genes from the random-effects meta-analysis across five region units in three cohorts. Dashed line indicates BH FDR 0.05; colored points mark the 244 genes estimable in every region unit with preserved direction after each cohort omission. Gray points include other significant genes. The region covariance was not explicitly estimated. b The 244-gene subset contains 128 positive and 116 negative AD-minus-control associations. Direction preservation does not establish preservation of statistical significance or therapeutic direction. No genes are omitted.'

sen=read(ST+'Table_S2_iMGL_senescence.tsv');iron=read(ST+'Table_S3_iMGL_iron_LPO.tsv');astro=read(ST+'Table_S4_astrocyte_modules.tsv')
NAMES={'M3_SASP_Senescence':'SASP','inflammation_control':'Inflammation control','senescence_core':'Senescence core','M5_Iron_Uptake':'Iron uptake','M5_Iron_Storage':'Iron storage','M5_Iron_Export':'Iron export','M5_Ferritinophagy':'Ferritinophagy','M6_LPO_ACSL4':'LPO/ACSL4','M7_FerroAging':'Ferro-aging reference','M7_FerroAging_Author':'Ferro-aging reference'}
required=CFG['fig2c_missing_policy']['required']
valid_c=np.isfinite(astro[required].to_numpy(float)).all(axis=1)
excluded_c=astro.loc[~valid_c,'module'].tolist()
assert len(astro)==7 and valid_c.sum()==6 and excluded_c==['M5_Iron_Storage']
(DATA/'Fig2c_exclusion_log.json').write_text(json.dumps({'before':7,'displayed':6,'excluded_modules':excluded_c,'rule':'All plotted coefficients, confidence limits and recorded BH FDR must be finite.','reason':'Iron storage had only one of its two module genes represented; no module effect could be estimated.','FDR_family_preserved':7,'imputation':False},indent=2),encoding='utf-8')
astro=astro.loc[valid_c].reset_index(drop=True)
save_data(sen,'Fig2a_iMGL_senescence.tsv');save_data(iron,'Fig2b_iMGL_iron_LPO.tsv');save_data(astro[['dataset','analysis','module','n_present','n_total','n_case','n_control']+required],'Fig2c_astrocyte.tsv')
fig=plt.figure(figsize=(width_mm/25.4,160/25.4));gs=fig.add_gridspec(2,2,height_ratios=[1.2,1],left=.28,right=.96,bottom=.09,top=.91,wspace=.52,hspace=.56)
def heat(ax,d,effect):
 order=list(d.module.drop_duplicates()); z=d.pivot(index='module',columns='time_hours',values=effect).loc[order]
 q=d.pivot(index='module',columns='time_hours',values='FDR').loc[order]
 im=ax.imshow(z,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto')
 ax.set_yticks(range(len(z)),[NAMES[x] for x in z.index]);ax.set_xticks(range(3),['4','24','96']);ax.set_xlabel('Time after Aβ (h)')
 for i in range(len(z)):
  for j in range(3):
   if q.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=10,color='white' if abs(z.iloc[i,j])>1.5 else 'black')
 return im
ax=fig.add_subplot(gs[0,0]);heat(ax,sen,'logFC_like_mean_difference');panel(ax,'a')
ax=fig.add_subplot(gs[0,1]);im=heat(ax,iron,'effect_treated_minus_control');panel(ax,'b')
cb=fig.colorbar(im,ax=[fig.axes[0],ax],location='top',fraction=.07,pad=.12);cb.set_label('Treated − baseline module mean (expression-score units)',fontsize=8)
ax=fig.add_subplot(gs[1,:]); y=np.arange(len(astro));
for i,r in astro.iterrows():
 if np.isfinite(r.logFC):ax.errorbar(r.logFC,i,xerr=[[r.logFC-r.CI_low],[r.CI_high-r.logFC]],fmt='o',color=BLUE if r.FDR_BH_7<.05 else '#777777',markersize=4,capsize=2)
ax.axvline(0,color='#777777',ls='--',lw=.7);ax.set_yticks(y,[NAMES[x] for x in astro.module]);ax.invert_yaxis();ax.set_xlabel('Adjusted disease coefficient (APP/PS1-21 − wild type)');panel(ax,'c');save(fig,2)
CAPTIONS['Fig2']='Model-specific glial transcriptional programs. a,b All nine senescence/control and eighteen iron/LPO/reference module-time comparisons in the iMGL time-course analysis. Color is the recorded treated-minus-baseline expression-score mean difference; asterisks denote the source BH FDR below 0.05 within the respective analysis families. Each comparison has three treated and three baseline culture replicates from one donor, not three biological donors. These RNA scores do not directly measure senescence, iron flux, lipid peroxidation or lactylation. c Six estimable modules in the primary mouse astrocyte model, adjusted for circadian terms, batch and sex (24 APP/PS1-21 and 24 wild-type animal units). Points are disease coefficients; horizontal bars are the source 95% confidence intervals. Blue indicates the original seven-test BH FDR below 0.05, gray an estimable nonsignificant result. Iron storage was omitted from this panel because only one of its two genes was present, preventing estimation of the module effect. The prespecified seven-module correction family was unchanged. Units and models differ between a,b and c and their effect magnitudes should not be pooled.'

pairs=read(ST+'Table_S5_E1_compound_target_pairs.csv');network=read('results/2026-10-03-stage-a-mechanism-030/network_proximity_experiments400.tsv');network=network[network.interface=='G1_robust244'].copy()
cnt=pairs.groupby('靶点').size().sort_values();assert len(cnt)==20 and cnt.sum()==80 and len(network)==8
cnt.rename_axis('target').reset_index(name='unique_compound_target_pairs').to_csv(DATA/'Fig3a_target_pair_counts.tsv',sep='\t',index=False);save_data(network,'Fig3b_G1_primary_network.tsv')
fig,axes=plt.subplots(1,2,figsize=(width_mm/25.4,130/25.4),gridspec_kw={'width_ratios':[1,1.1]});fig.subplots_adjust(left=.12,right=.98,bottom=.16,top=.91,wspace=.65)
ax=axes[0];ax.barh(cnt.index,cnt.values,color=BLUE);ax.set_xlabel('Unique accepted compound–target pairs');ax.set_xticks([0,5,10,15,20,25]);panel(ax,'a')
ax=axes[1];lab=[f'{r.tier} / {r.overlap_mode}' for r in network.itertuples()];ax.scatter(network.FDR_network,range(8),color=BLUE,s=25);ax.axvline(.05,color=RED,ls='--',lw=.8);ax.set_yticks(range(8),lab);ax.invert_yaxis();ax.set(xlim=(0,1.08),xlabel='Primary-network BH FDR');panel(ax,'b');save(fig,3)
CAPTIONS['Fig3']='Experimental-record connections do not provide positive global human network proximity. a All 80 deduplicated accepted experimental-record compound–target pairs across twenty program-connected targets. Counts represent source-reviewed records, not administered composition, target engagement, equal assay quality or therapeutic efficacy. b All eight human-signature comparisons in the primary experimental-channel network at confidence 400: evidence tiers E1, E2, H and P, with shared drug/disease members included or removed. Each point is the recorded network BH FDR from 10,000 degree-matched randomizations; dashed line marks 0.05. No comparison is significant. Evidence-tier definitions and assay restrictions are provided in the methods. Incomplete target/network coverage limits inference from this negative proximity result.'

review=read('results/2026-10-11-virtual-ko-scientific-review-035/target_verification_conclusions.tsv')
ref=read('results/2026-10-11-virtual-ko-matching-correction-036/functional_response_comparison.tsv')
storage=ref[(ref.variant=='corrected_full_CP10K') & ref.function.str.contains('Storage')].copy()
assert len(review)==29 and len(storage)==58
fields=CFG['fig4_missing_policy']['panel_a_required']
valid_a=np.isfinite(review[fields].to_numpy(float)).all(axis=1)
valid_b=np.isfinite(storage[CFG['fig4_missing_policy']['panel_b_required']].to_numpy(float)).all(axis=1)
assert (review.loc[~valid_a,'original_shared_weighted_terms']==0).all()
assert (storage.loc[~valid_b,'n_members']==1).all()
keys=['analysis','context','target']
excluded_a=review.loc[~valid_a,keys].to_dict('records')
excluded_b=storage.loc[~valid_b,keys+['seed']].to_dict('records')
exclusions={'panel_a':{'before':29,'displayed':int(valid_a.sum()),'excluded':excluded_a,'reason':'No original shared weighted pathway terms; reproduction-rate denominator is zero.'},'panel_b':{'before':58,'displayed':int(valid_b.sum()),'excluded':excluded_b,'reason':'Only one iron-storage module gene represented; module reference comparison not estimable.'},'policy':'Exclude nonfinite plotted statistics only; retain finite zero and nonsignificant values. Original source tables are unchanged.'}
(DATA/'Fig4_exclusion_log.json').write_text(json.dumps(exclusions,indent=2),encoding='utf-8')
review=review.loc[valid_a].reset_index(drop=True)
storage=storage.loc[valid_b].copy()
assert len(review)==20 and len(storage)==42
save_data(review[keys+['original_shared_weighted_terms','lodo_planned_pairs','lodo_estimable_pairs','bootstrap_planned_pairs','bootstrap_estimable_pairs']+fields],'Fig4a_estimable_reproduction_rates.tsv')
save_data(storage[keys+['seed','function','n_members','genes','n_selected_controls','n_completed_controls','reference_tail_fraction']],'Fig4b_estimable_storage_comparisons.tsv')
# These obsolete files were generated figure exports, never raw inputs.
for name in ['Fig4a_all_29_KO_reviews.tsv','Fig4b_all_58_corrected_storage_comparisons.tsv']:
 (DATA/name).unlink(missing_ok=True)
fig,axes=plt.subplots(1,2,figsize=(width_mm/25.4,180/25.4),gridspec_kw={'width_ratios':[.8,1.15]});fig.subplots_adjust(left=.24,right=.98,bottom=.10,top=.92,wspace=.92)
def context_label(r):
 return f'{r.target} / '+('A' if r.context.startswith('Astro') else 'M')+' / '+('AD' if r.context.endswith('_AD') else 'Ctrl')+('†' if r.analysis=='target_retention' else '')
labels=[context_label(r) for r in review.itertuples()]
ax=axes[0]; z=review[fields].to_numpy(float);im=ax.imshow(z,vmin=0,vmax=1,cmap='Blues',aspect='auto');ax.set_yticks(range(len(review)),labels);ax.set_xticks([0,1],['LODO','Bootstrap']);
for i in range(len(review)):
 for j in range(2):ax.text(j,i,f'{z[i,j]:.2f}',ha='center',va='center',color='white' if z[i,j]>.65 else 'black',fontsize=8)
panel(ax,'a');ax.set_xlabel('Median unconditional reproduction rate')
ax=axes[1]
contexts=storage[keys].drop_duplicates().reset_index(drop=True)
for i,r in enumerate(contexts.itertuples()):
 sub=storage[(storage.target==r.target)&(storage.context==r.context)&(storage.analysis==r.analysis)]
 for _,s in sub.iterrows():
  if pd.notna(s.reference_tail_fraction):ax.scatter(s.reference_tail_fraction,i+(-.13 if s.seed==42 else .13),marker='o' if s.seed==42 else 's',s=20,facecolors=BLUE if s.seed==42 else 'none',edgecolors=BLUE)
ax.set_yticks(range(len(contexts)),[context_label(r) for r in contexts.itertuples()]);ax.set(ylim=(len(contexts)-.3,-.7),xlim=(-.03,1.03),xlabel='Corrected reference tail fraction');ax.scatter([],[],s=20,color=BLUE,label='Seed 42');ax.scatter([],[],s=20,facecolors='none',edgecolors=BLUE,marker='s',label='Seed 43');ax.legend(loc='upper right',bbox_to_anchor=(1,1.08),ncol=2,columnspacing=.8,handletextpad=.3);panel(ax,'b');save(fig,4)
CAPTIONS['Fig4']='Virtual knockout stability and limitations of reference specificity. a Twenty target/background combinations with estimable median unconditional reproduction proportions of original shared weighted pathway terms under leave-one-donor-out (LODO) and donor bootstrap resampling. Nine of the twenty-nine assessed combinations had no original shared weighted terms and were excluded because the reproduction-rate denominator was zero. A, astrocyte; M, microglia; Ctrl, control; dagger, target-retention sensitivity analysis. These summaries do not measure replication in independent donors. b Forty-two estimable corrected full-CP10K reference comparisons across twenty-one combinations for the historical two-gene iron-storage module (two seeds per combination). Sixteen of the fifty-eight assessed records represented only one module gene and were excluded. The panels therefore have separate row labels. Finite zero and nonsignificant values are retained. Reference tail fractions are descriptive comparisons against twenty matched controls, not formal specificity P values or therapeutic effects. This historical module was subsequently retired from candidate selection; other modules and the complete assessment records remain in the source tables.'

audit={'configuration':CFG,'sources_sha256':SOURCES,'outputs':OUTPUTS,'captions':CAPTIONS,
 'upstream_scientific_tables_unchanged':all(sha(ROOT/p)==v for p,v in SOURCES.items()),
 'runtime':{'python':__import__('sys').version,'matplotlib':matplotlib.__version__,'pandas':pd.__version__},
 'image_integrity':'Plots from recorded tables only; no artificial observations, fit or new P value. Fig2c and Fig4 exclude nonfinite plotted statistics only, with counts and reasons in their source_data exclusion logs. Source tables and multiple-testing families unchanged.',
 'visual_review':'PENDING'}
(OUT/'figure_provenance.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'figure_legends.md').write_text('\n\n'.join(f'**{k}.** {v}' for k,v in CAPTIONS.items())+'\n',encoding='utf-8')
(OUT/'figure_captions.json').write_text(json.dumps([{'caption':f'Fig. {n} '+CAPTIONS[f'Fig{n}'],'pdf':f'Fig{n}.pdf','tiff':f'Fig{n}.tif'} for n in range(1,5)],ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'figures':4,'sources':len(SOURCES),'source_preservation':audit['upstream_scientific_tables_unchanged']}))

