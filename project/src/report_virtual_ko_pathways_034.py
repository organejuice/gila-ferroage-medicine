"""Audit KO-specific backgrounds, replay ORA independently, and report pathway contributors."""
from pathlib import Path
import hashlib, importlib.util, json, sys
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, beta
import yaml

ROOT=Path(__file__).resolve().parents[1]
CFG=yaml.safe_load((ROOT/'configs/virtual_ko_pathways_034.yaml').read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir'];ALPHA=CFG['fdr']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(values):
    v=np.asarray(values,float);out=np.full(len(v),np.nan);ix=np.flatnonzero(np.isfinite(v));order=ix[np.argsort(v[ix],kind='stable')]
    out[order]=np.minimum(1,np.minimum.accumulate((v[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1]);return out
def genes(x):return set() if pd.isna(x) or x=='' else set(str(x).split(';'))
def main():
    frozen=json.loads((OUT/'frozen_input_hashes.json').read_text(encoding='utf-8'))
    for p,h in frozen.items():assert sha(ROOT/p)==h,p
    jobs=pd.read_csv(OUT/'job_registry.tsv',sep='\t');assert len(jobs)==84
    registry=pd.read_csv(OUT/'assessment_registry.tsv',sep='\t');assert len(registry)==198
    spec=importlib.util.spec_from_file_location('ko_original_layers',ROOT/'src/133_summarize_virtual_ko.py')
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    memberships=[]
    for row in registry.to_dict('records'):
        for layer,targets in original.LAYERS.items():
            if row['target'] in targets and ((row['cell_type']=='Astrocyte')==layer.startswith('G5')):
                memberships.append(dict(layer=layer,research_category=original.PRIORITY[layer][row['target']],**row))
    pd.DataFrame(memberships).to_csv(OUT/'stratified_assessment_registry.tsv',sep='\t',index=False)
    assert registry.loc[registry.cell_type=='Microglia','target'].nunique()==20 and registry.loc[registry.cell_type=='Astrocyte','target'].nunique()==13
    libraries={}
    for lib in CFG['libraries']:
        rows=[s.split('\t') for s in (ROOT/CFG['library_dir']/f'{lib}.gmt').read_text(encoding='utf-8').splitlines()]
        libraries[lib]={s[0]:(s[1],set(s[2:])) for s in rows}
    all_o=[];all_g=[];all_q=[];effects=[];audit=[]
    for job in jobs.to_dict('records'):
        d=OUT/'jobs'/job['job_id'];a=json.loads((d/'completion.json').read_text(encoding='utf-8'))
        assert a['status']=='PASS' and a['source_sha256']==job['sha256'] and a['job_id']==job['job_id']
        dr=pd.read_csv(ROOT/job['path'],sep='\t');dr=dr[dr.gene!=job['target']].copy()
        assert dr.gene.is_unique and np.isfinite(dr.distance).all() and (dr.distance>=0).all()
        bg=set(dr.gene);sig=set(dr.loc[dr['p.adj']<ALPHA,'gene'])
        for n in ['job_id','analysis','context','cell_type','state','target','seed']:dr[n]=job[n]
        dr['network_significant']=dr['p.adj']<ALPHA;dr['expression_direction']='NOT_ESTIMATED'
        effects.append(dr)
        o=pd.read_csv(d/'ORA_all.tsv.gz',sep='\t');g=pd.read_csv(d/'GSEA_all.tsv.gz',sep='\t');q=pd.read_csv(d/'QC.tsv',sep='\t')
        assert len(o)==a['ora_tests'] and len(g)==a['gsea_tests']
        assert not o.duplicated(['library','term']).any() and not g.duplicated(['library','mode','pathway']).any()
        for df in [o,g,q]:
            for n in ['job_id','analysis','context','target','seed']:assert set(df[n])=={job[n]},(job,n)
        assert set(zip(q.library,q['mode']))=={(lib,mode) for lib in CFG['libraries'] for mode in CFG['gsea']['modes']}
        assert q.missing_P.sum()==0,(job,q.to_dict('records'))
        max_error=0
        mapping=pd.read_csv(d/'gene_mapping.tsv.gz',sep='\t')
        for lib,terms in libraries.items():
            annotated=set.union(*(v[1] for v in terms.values()));universe=bg&annotated;query=sig&universe
            expected_o={t for t,(_,s) in terms.items() if CFG['ora']['min_size']<=len(s&universe)<=CFG['ora']['max_size']}
            x=o[o.library==lib].copy();assert set(x.term)==expected_o
            m=mapping[mapping.library==lib]
            assert set(m.gene)==bg and len(m)==len(bg)
            assert set(m.loc[m.annotated,'gene'])==universe and set(m.loc[m.significant,'gene'])==sig
            expected_p=[]
            for row in x.itertuples():
                term_genes=terms[row.term][1]&universe;hits=term_genes&query
                assert row.background_size==len(universe) and row.query_size==len(query) and row.term_size==len(term_genes) and row.overlap==len(hits)
                assert genes(row.genes)==hits and job['target'] not in hits
                assert row.description==terms[row.term][0]
                expected_p.append(hypergeom.sf(len(hits)-1,len(universe),len(term_genes),len(query)) if query else 1.)
            err=float(np.max(np.abs(x.P.to_numpy()-np.array(expected_p))))
            max_error=max(max_error,err);assert err<1e-10
            assert np.allclose(bh(x.P),x.FDR_library,atol=1e-10)
            expected_g={t for t,(_,s) in terms.items() if CFG['gsea']['min_size']<=len(s&bg)<=CFG['gsea']['max_size']}
            for mode in CFG['gsea']['modes']:
                y=g[(g.library==lib)&(g['mode']==mode)]
                assert set(y.pathway)==expected_g and np.isfinite(y[['pval','NES','ES']]).all().all()
                assert ((y.pval>=0)&(y.pval<=1)).all() and ((y.ES>=-1e-12)&(y.ES<=1+1e-12)).all()
                assert np.allclose(bh(y.pval),y.FDR_library,atol=1e-10)
                for row in y.itertuples():
                    assert row.description==terms[row.pathway][0] and row.size==len(terms[row.pathway][1]&bg)
                    assert genes(row.leadingEdge)<=terms[row.pathway][1]&bg
        assert np.allclose(bh(o.P),o.FDR_job_all_libraries,atol=1e-10)
        for mode in CFG['gsea']['modes']:
            x=g[g['mode']==mode];assert np.allclose(bh(x.pval),x.FDR_job_all_libraries,atol=1e-10)
        all_o.append(o);all_g.append(g);all_q.append(q)
        audit.append(dict(job_id=job['job_id'],status='PASS',background_genes=len(bg),significant_genes=len(sig),ora_tests=len(o),gsea_tests=len(g),max_hypergeom_error=max_error))
        print('AUDITED',len(audit),'/84',job['job_id'],flush=True)
    o=pd.concat(all_o,ignore_index=True);g=pd.concat(all_g,ignore_index=True);qc=pd.concat(all_q,ignore_index=True);effect=pd.concat(effects,ignore_index=True)
    o['FDR_run_all_tests']=bh(o.P)
    g['FDR_run_all_tests']=g.groupby('mode',sort=False).pval.transform(lambda x:bh(x))
    n=CFG['gsea']['nperm'];k=g.nMoreExtreme.to_numpy()
    assert np.allclose(g.pval,(k+1)/(n+1),atol=1e-12)
    g['P_MC_lower95']=np.where(k==0,0,beta.ppf(.025,k,n-k+1))
    g['P_MC_upper95']=np.where(k==n,1,beta.ppf(.975,k+1,n-k))
    g['FDR_job_MC_upper95']=g.groupby(['job_id','mode'],sort=False).P_MC_upper95.transform(lambda x:bh(x))
    g['significant_MC_upper95']=g.FDR_job_MC_upper95<ALPHA
    o['interpretation_eligible']=(o.FDR_job_all_libraries<ALPHA)&(o.overlap>=CFG['ora']['min_overlap'])
    g['interpretation_eligible']=g.FDR_job_all_libraries<ALPHA
    for df,name in [(o,'ORA_all_tests.tsv.gz'),(g,'GSEA_all_tests.tsv.gz'),(effect,'gene_perturbation_all.tsv.gz')]:df.to_csv(OUT/name,sep='\t',index=False,compression='gzip')
    o[o.interpretation_eligible].to_csv(OUT/'ORA_significant.tsv',sep='\t',index=False)
    g[g.interpretation_eligible].to_csv(OUT/'GSEA_significant.tsv',sep='\t',index=False)
    effect[effect.network_significant].to_csv(OUT/'gene_perturbation_significant.tsv',sep='\t',index=False)
    qc.to_csv(OUT/'QC_all_jobs.tsv',sep='\t',index=False)
    pd.DataFrame(audit).to_csv(OUT/'audit_jobs.tsv',sep='\t',index=False)
    # Two seeds are stability checks, not independent biological replicates.
    stable=[];summary=[]
    for keys,j in jobs.groupby(['analysis','context','target'],sort=True):
        analysis,context,target=keys;assert set(j.seed)=={42,43}
        eo=o[(o.analysis==analysis)&(o.context==context)&(o.target==target)]
        eg=g[(g.analysis==analysis)&(g.context==context)&(g.target==target)]
        ds=effect[(effect.analysis==analysis)&(effect.context==context)&(effect.target==target)]
        sig42=set(ds.loc[(ds.seed==42)&ds.network_significant,'gene']);sig43=set(ds.loc[(ds.seed==43)&ds.network_significant,'gene'])
        stats=dict(analysis=analysis,context=context,target=target,n_sig_42=len(sig42),n_sig_43=len(sig43),n_shared_genes=len(sig42&sig43),shared_genes=';'.join(sorted(sig42&sig43)))
        layers=[layer for layer,targets in original.LAYERS.items() if target in targets and (context.startswith('Astrocyte')==layer.startswith('G5'))]
        stats['layers']=';'.join(layers);stats['research_categories']=';'.join(f'{layer}:{original.PRIORITY[layer][target]}' for layer in layers)
        for method,mode,table,term_col,contrib in [('ORA','significant_gene_set',eo,'term','genes')]+[('GSEA',mode,eg[eg['mode']==mode],'pathway','leadingEdge') for mode in CFG['gsea']['modes']]:
            a=table[(table.seed==42)&table.interpretation_eligible];b=table[(table.seed==43)&table.interpretation_eligible]
            sa=set(zip(a.library,a[term_col]));sb=set(zip(b.library,b[term_col]));stats[f'{method}_{mode}_terms42']=len(sa);stats[f'{method}_{mode}_terms43']=len(sb);stats[f'{method}_{mode}_shared']=len(sa&sb)
            merged=a.merge(b,on=['library',term_col],suffixes=('_42','_43'))
            for r in merged.to_dict('records'):
                ga=genes(r[contrib+'_42']);gb=genes(r[contrib+'_43'])
                stable.append(dict(analysis=analysis,context=context,target=target,method=method,mode=mode,library=r['library'],term=r[term_col],description=r['description_42'],FDR42=r['FDR_job_all_libraries_42'],FDR43=r['FDR_job_all_libraries_43'],NES42=r.get('NES_42'),NES43=r.get('NES_43'),MC_upper_FDR42=r.get('FDR_job_MC_upper95_42'),MC_upper_FDR43=r.get('FDR_job_MC_upper95_43'),contributors42=';'.join(sorted(ga)),contributors43=';'.join(sorted(gb)),shared_contributors=';'.join(sorted(ga&gb)),contributor_jaccard=len(ga&gb)/len(ga|gb) if ga|gb else np.nan))
        summary.append(stats)
    ss=pd.DataFrame(summary);st=pd.DataFrame(stable,columns=['analysis','context','target','method','mode','library','term','description','FDR42','FDR43','NES42','NES43','MC_upper_FDR42','MC_upper_FDR43','contributors42','contributors43','shared_contributors','contributor_jaccard'])
    st.to_csv(OUT/'two_seed_pathways.tsv',sep='\t',index=False);ss.to_csv(OUT/'target_pathway_summary.tsv',sep='\t',index=False)
    wg=st[(st.method=='GSEA')&(st['mode']=='weighted')];ug=st[(st.method=='GSEA')&(st['mode']=='unweighted_sensitivity')]
    robust=wg.merge(ug,on=['analysis','context','target','library','term'],suffixes=('_weighted','_unweighted'))
    robust['all_four_MC_upper_FDR_pass']=np.logical_and.reduce([robust[n]<ALPHA for n in ['MC_upper_FDR42_weighted','MC_upper_FDR43_weighted','MC_upper_FDR42_unweighted','MC_upper_FDR43_unweighted']])
    robust.to_csv(OUT/'two_seed_weighting_robust_GSEA.tsv',sep='\t',index=False)
    reports=OUT/'target_reports';reports.mkdir(exist_ok=True)
    for row in ss.to_dict('records'):
        analysis,context,target=[row[k] for k in ['analysis','context','target']]
        ds=effect[(effect.analysis==analysis)&(effect.context==context)&(effect.target==target)]
        ts=st[(st.analysis==analysis)&(st.context==context)&(st.target==target)]
        lines=[f'# {target}：{analysis} / {context}','','结果只表示网络调控扰动，不能判断表达上/下调、通路激活/抑制或治疗获益。两个种子不是独立队列。', '',f"原研究定位：{row['research_categories']}。",'',f"显著受扰动基因：seed42={row['n_sig_42']}，seed43={row['n_sig_43']}；共同={row['n_shared_genes']}。",'','共同显著基因：'+row['shared_genes'],'','## 两种子均显著的通路','', '|方法|模式|库|通路|FDR42|FDR43|共同贡献基因|','|---|---|---|---|---:|---:|---|']
        for r in ts.sort_values(['method','mode','FDR42']).head(45).to_dict('records'):
            lines.append(f"|{r['method']}|{r['mode']}|{r['library']}|{r['description']}|{r['FDR42']:.3g}|{r['FDR43']:.3g}|{r['shared_contributors']}|")
        if not len(ts):lines.append('|—|—|—|未检出两种子共同显著通路|—|—|—|')
        lines+=['','展示最多45行，完整共同通路在two_seed_pathways.tsv；显著性主标准为每个KO/模式跨三个库合并BH<0.05，ORA还要求至少2个重叠基因。','', '## 每个种子的主要受扰动基因','', '|种子|基因|distance|网络FC|网络FDR|','|---|---|---:|---:|---:|']
        for seed in [42,43]:
            x=ds[(ds.seed==seed)&ds.network_significant].sort_values('distance',ascending=False).head(15)
            for r in x.to_dict('records'):lines.append(f"|{seed}|{r['gene']}|{r['distance']:.4g}|{r['FC']:.4g}|{r['p.adj']:.3g}|")
        lines+=['','FC为相对背景的平方扰动距离比，不是表达倍数。所有显著基因见gene_perturbation_significant.tsv，全部受检基因见gene_perturbation_all.tsv.gz。','', '## 单种子富集详情','']
        for method,mode,table,termcol,contrib in [('ORA','significant_gene_set',o,'term','genes')]+[('GSEA',mode,g[g['mode']==mode],'pathway','leadingEdge') for mode in CFG['gsea']['modes']]:
            x=table[(table.analysis==analysis)&(table.context==context)&(table.target==target)&table.interpretation_eligible]
            lines+=[f'### {method} / {mode}','','|种子|库|通路|FDR|贡献基因|','|---|---|---|---:|---|']
            for seed in [42,43]:
                for lib in CFG['libraries']:
                    subset=x[(x.seed==seed)&(x.library==lib)].sort_values('FDR_job_all_libraries').head(5)
                    for r in subset.to_dict('records'):lines.append(f"|{seed}|{lib}|{r['description']}|{r['FDR_job_all_libraries']:.3g}|{r[contrib]}|")
            lines+=['','每种子每库最多展示5条；完整结果包括所有未显著/零重叠检验，见ORA_all_tests.tsv.gz、GSEA_all_tests.tsv.gz。']
        (reports/f'{analysis}__{context}__{target}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    lines=['# 虚拟敲除后GO、KEGG、Reactome与GSEA报告','',f'已完成84个有效KO输出的三个库ORA与两种GSEA分析（加权主分析/不加权敏感性）；全部198项原评估保留，低表达不可估计不被解释为阴性。未重算网络或KO。','',
        '## 方法与解释','', 'GO仅复用既有GO_BP（生物过程），不含GO细胞组分或分子功能。GO/KEGG/Reactome来自阶段A冻结GMT快照，不联网更新。ORA用每个KO实际受检且该库有注释的基因作背景，显著网络扰动基因为前景，排除KO自身；大小5–500，跨库BH<0.05且至少2重叠基因。GSEA保留全部受检基因，按非负扰动距离降序，scoreType=pos，大小15–500；gseaParam1主分析、0敏感性，fgseaSimple10万次固定基因集置换，seed34034+job序号。四个本地单核R工作进程写入互不重叠的任务目录。GSEA中正NES表示该集合集中在高扰动端，不是通路激活。','',
        '主FDR是每个KO/模式内跨三个库BH。另附库内和全运行FDR；主FDR不控制挑选整个运行任意结果的错误率，跨靶点比较可参考FDR_run_all_tests。跨种子和加权敏感性检查不是独立生物学验证。相同通路高度重叠、共同贡献基因和少数极端扰动会造成重复富集，不能据通路数量排序靶点。AD与对照分别分析，不将显著/不显著差异当成正式疾病交互检验。','',
        '## 逐靶点通路结果','', '|分析|细胞背景|靶点|显著基因42/43|共同基因|ORA共同通路|加权GSEA共同通路|不加权GSEA共同通路|','|---|---|---|---:|---:|---:|---:|---:|']
    for r in ss.to_dict('records'):
        link=f"target_reports/{r['analysis']}__{r['context']}__{r['target']}.md"
        lines.append(f"|{r['analysis']}|{r['context']}|[{r['target']}]({link})|{r['n_sig_42']}/{r['n_sig_43']}|{r['n_shared_genes']}|{r['ORA_significant_gene_set_shared']}|{r['GSEA_weighted_shared']}|{r['GSEA_unweighted_sensitivity_shared']}|")
    lines+=['','表中的通路数为带库标识的注释条目，未去除GO/Reactome嵌套冗余，不能当独立机制数量。详细报告列出受扰动基因、强度、FDR、通路及贡献基因；全结果在配套TSV/GZ中。','',
        '## 审计与边界','',f'84/84输出逐个核验输入哈希、靶点排除、受检背景、库注释、所有通路覆盖、两种GSEA模式、leading edge与前景贡献基因、BH校正。ORA全检验P值与独立scipy超几何尾概率最大误差{max(a["max_hypergeom_error"] for a in audit):.3g}。',
        f'ORA检验共{len(o):,}项；GSEA检验共{len(g):,}项；GSEA中零次超越的记录共{int(qc.p_at_resolution.sum()):,}，最小P值为1/100001，不能排序比此分辨率更小的尾概率。P_MC_lower95/upper95是逐项二项分布置换误差区间，FDR_job_MC_upper95为上限敏感性（不是同时置信区间或生物学效应区间）。警告保留在QC_all_jobs.tsv；自适应初试只作性能记录，不纳入结果。','',
        '细胞标签、固定供体/脑区抽样、网络随机性和G5人类补充背景的限制继承033报告。维持原候选类别及对接限制。既有数据库、R环境和项目方法复用，无新增依赖；本机计算，Git/DVC本地缓存保存正式资产。未配置远端，当前为可追溯探索性验证，不称严格跨机器或真实生物学复现。未删除原始数据或正式KO结果；没有新建一次性转换文件。','',
        '方法依据：[scTenifoldKnk论文](https://pubmed.ncbi.nlm.nih.gov/35510185/) DOI10.1016/j.patter.2022.100434；[GSEA原始论文](https://pubmed.ncbi.nlm.nih.gov/16199517/) DOI10.1073/pnas.0506580102；[fgsea方法预印本](https://www.biorxiv.org/content/10.1101/060012v3) DOI10.1101/060012。参数选择及无符号排名适配详见本轮Decision Log。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    audit_result=dict(status='PASS',audited_jobs=84,expected_jobs=84,assessment_records=198,ORA_tests=len(o),GSEA_tests=len(g),input_hashes_checked=len(frozen),max_hypergeom_error=max(a['max_hypergeom_error'] for a in audit),no_expression_direction_claim=True)
    (OUT/'acceptance_audit.json').write_text(json.dumps(audit_result,indent=2)+'\n',encoding='utf-8')
    outputs=[p for p in OUT.iterdir() if p.is_file() and p.name!='result_hashes.json']+list(reports.glob('*.md'))
    (OUT/'result_hashes.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():sha(p) for p in outputs},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(audit_result),flush=True)
if __name__=='__main__':main()
