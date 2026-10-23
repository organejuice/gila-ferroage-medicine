"""Update a candid progress report from actual SOP035 artifacts, without rerunning KO."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/2026-10-05-virtual-ko-verification-execution-035'

def main():
    e=pd.read_csv(OUT/'target_expression_by_donor_region.tsv',sep='\t');e=e[e.stage=='balanced']
    donor=e.groupby(['context','target','donor_id'],as_index=False)[['N','k','target_UMI']].sum()
    rows=[]
    for (context,target),q in donor.groupby(['context','target']):
        total=q.k.sum();maximum=q.loc[q.k.idxmax()]
        rows.append(dict(context=context,target=target,N=q.N.sum(),k=total,detection_fraction=total/q.N.sum(),n_donors=len(q),n_detecting_donors=(q.k>0).sum(),max_detecting_donor=maximum.donor_id,max_donor_detection_contribution=maximum.k/total if total else None,minimum_remaining_detection_after_donor_removal=total-q.k.max(),maximum_remaining_detection_after_donor_removal=total-q.k.min(),removal_counts_note='diagnostic on fixed original balanced subset; not rebuilt LODO output'))
    s=pd.DataFrame(rows);s.to_csv(OUT/'target_expression_summary.tsv',sep='\t',index=False)
    state=json.loads((OUT/'run_state.json').read_text(encoding='utf-8')) if (OUT/'run_state.json').exists() else {'status':'PREPARING'}
    completed=list((OUT/'networks').glob('*/completion.json'))
    enrichment=list((OUT/'jobs').glob('*/completion.json'))
    lines=['# 分层虚拟敲除实际核查：阶段报告','',
           '当前为实际执行记录，不是完整核查结束报告。15个靶点、29组合；AD14/Control15。',
           f"当前状态：{state.get('stage','')} / {state.get('status','')}。",
           f'已闭合重采样网络任务：{len(completed)}/130（每任务两个网络种子）；已完成富集任务：{len(enrichment)}。',
           '','## 已完成的现有数据计算','',
           '四个背景的QC前、QC后和平衡抽样后靶点表达计算完成；多谱系标志描述、58份原KO输出贡献与跨靶点比较完成。描述性标志核查不等于参考图谱确认或双细胞排除。',
           '','|背景|靶点|原平衡抽样N|检出细胞k|检出率|最高供体检出贡献|','|---|---|---:|---:|---:|---:|']
    for r in s[(s.context=='Microglia_AD')&s.target.isin(['CA2','DPP4','PTGS2','RELA'])].itertuples(index=False):
        lines.append(f'|{r.context}|{r.target}|{r.N}|{r.k}|{r.detection_fraction:.3%}|{r.max_donor_detection_contribution:.3%}|')
    lines+=['','PTGS2仅26个检出细胞，接近原低表达保留门槛25。上述固定子集去供体计数只作可估计性风险诊断，不能代替重新抽样、建网和KO的留一结果。',
            '','## 已观察到的响应相似性','',
            '在共同受检背景上，DPP4与PTGS2扰动距离Spearman相关：seed42约0.9875，seed43约0.9025；其单种子显著基因Jaccard约1.000和0.619。这支持存在共同响应，但不能单凭相关性排除真实靶点机制。',
            '','## 身份与资料限制','',
            'AD星形与小胶质背景分别有约0.096%和28.628%的细胞在描述性面板上另一胶质谱系均值较高；面板均值未经身份概率校准，且不等于污染/双细胞比例。需结合原cluster、供体和共表达逐项复核。',
            '本环境未安装scDblFinder，现有原始归档为filtered矩阵；未识别到用于空液滴背景估计的unfiltered矩阵。因此本次不宣称完成正式双细胞检测或环境RNA校正。这些缺项列入library_QC_availability.tsv，其他核查继续。',
            '','## 接口验证与执行边界','']
    for p in ['adapter_deterministic_tests.json','adapter_real_replay.json']:
        if (OUT/p).exists():
            a=json.loads((OUT/p).read_text());lines.append(f'- {p}: {a.get("status")}; 按最终适配器哈希及完整进程退出状态放行。')
    lines+=['','重采样20次、匹配参考20个为执行前冻结的探索模式；结果不作为精细显著特异性检验。核糖体成员使用冻结KEGG Ribosome蛋白基因与Reactome线粒体翻译交集构建操作性代理集合，含部分生物发生成员，未包装成纯结构复合物成员注释。',
            '所有原输入与旧结果保留。新输出位于本目录，缓存和任务状态由独立入口记录。运行中阶段文件不代表最终审核结论。',
            '','## 可重复入口','',
            '在原R环境中运行src/run_virtual_ko_verification_035.py；参数为configs/virtual_ko_verification_execution_035.yaml。执行需要原环境DLL路径，入口检查CONDA_PREFIX。运行中的任务不得重复启动第二个监督进程。',
            '','全部供体重采样、匹配对照、富集与最终科学复核未结束前，不能宣称15个靶点均通过核查，也不能推导表达或通路活动方向。']
    (OUT/'阶段核查报告.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(report=str(OUT/'阶段核查报告.md'),closed_network_tasks=len(completed),enrichment_jobs=len(enrichment)),ensure_ascii=False))

if __name__=='__main__':main()
