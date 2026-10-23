"""Build the reviewable delivery archive after document and reference audits."""
from pathlib import Path
import hashlib,json,re,zipfile,platform,importlib.metadata
from pypdf import PdfReader
import build_fig_yzdhd_submission as build
ROOT=build.ROOT; PROC=build.PROC; OUT=build.OUT; QA=build.QA
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    admin=OUT/'administrative'; auditdir=OUT/'audit'; auditdir.mkdir(exist_ok=True)
    source=(build.RAW/'manuscript_en.md').read_text(encoding='utf-8')
    manuscript=(OUT/'manuscript/main_manuscript.md').read_text(encoding='utf-8')
    wordcount=len(manuscript.split('## Abstract\n',1)[1].split('Keywords:',1)[0].split())
    (admin/'Chinese_abstract.md').write_text('# '+build.TITLE+'\n\n## 中文摘要\n\n'+source.split('## 中文摘要\n',1)[1].strip()+'\n',encoding='utf-8')
    (admin/'submission_checklist_CN.md').write_text('''# 投稿核对清单

本文为本项目益智地黄汤与阿尔茨海默病胶质细胞研究。ADNI文件仅用于作者及贡献元数据；此前ADNI正文投稿包不适用于本次投稿。

## 已完成

- Research文章；单盲稿保留作者信息。
- 233词摘要、6关键词、66篇作者—年份引用及字母序参考文献；所有文献均引用。
- 2张可编辑正文表；4张主图嵌入Word并提供独立174mm、RGB、600dpi TIFF及嵌入字体PDF。
- 标题页含9名作者、单位、通讯邮箱、贡献、致谢、资助与利益冲突；用户补充的Shi Jian与Yu Zhang邮箱已录入。
- 用户确认无资助、无利益冲突；未挪用ADNI资助或致谢。
- 原8张补充表字节一致保留，另提供标准CSV；Online Resource 1–3有元数据及正文引用、描述。
- 研究结果、关闭范围和C0证据上限保持；未重跑统计或创造实验结论。
- 文献独立终审、原稿保存审计与图像导出校验；逐页视觉检查记录。

## 作者实际提交前确认

1. 全体作者同意最终正文、顺序及本项目贡献，确认原创性与未同时向其他期刊提交。作者贡献与利益冲突仍须在系统内录入。
2. 确认本机构对公开去标识组学二次分析的伦理判定要求；本文仅说明公共数据及未新增实验，不虚构委员会或豁免号。
3. 原公开研究的数据使用条款、引用、代码和第三方软件许可符合此次共享用途。整套源分子数据及跨机器重建尚未独立恢复，不能据代码ZIP声称全流程严格复现。
4. 可选ORCID、推荐审稿人与出版选项仅据真实信息填写；上传后检查系统合成PDF。
5. 本文具有转录组与网络功能分析，但药理候选仍需实验验证。按研究文章提交的范围与新发现门槛由编辑判断。

本包完成材料制作，未向期刊投稿系统提交。
''',encoding='utf-8')
    for name in ['original_reference_audit.json','original_reference_audit.md','final_reference_audit.json','final_reference_audit.md','package_audit.json','source_manifest.json','author_metadata.json','editorial_insertions.json','title_revision_audit.json']:
        p=PROC/name
        assert p.exists(),p
        (auditdir/name).write_bytes(p.read_bytes())
    documents=[]
    for stem in ['main_manuscript','title_page','cover_letter','ESM_1']:
        folder=QA/('render_'+stem); pdf=folder/(stem+'.pdf'); count=len(PdfReader(pdf).pages)
        pages=sorted(folder.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[-1]))
        assert len(pages)==count
        documents.append({'name':stem,'pages':count,'pdf_sha256':sha(pdf),'page_image_hashes':[{'page':p.name,'sha256':sha(p)} for p in pages]})
    visual={'status':'Passed','method':'All final page images visually reviewed by the model; final PDFs and page-image hashes retained. Figures reviewed separately and source/export checks recorded.','documents':documents}
    (auditdir/'visual_QA.json').write_text(json.dumps(visual,ensure_ascii=False,indent=2),encoding='utf-8')
    environment={'bundle':'26.1013.11914','document_python':'C:/Users/administer/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe','document_python_version':platform.python_version(),'figure_python':'D:/conda/envs/adni_ai/python.exe','office':'D:/Libreoffice/program/soffice.exe','office_reason':'User explicitly specified native local office; bundle contains no Windows LibreOffice','document_packages':{n:importlib.metadata.version(n) for n in ['python-docx','pypdf','Pillow']},'no_dependencies_installed':True,'no_new_conda_environment':True,'scope':'Editorial continuation of existing project'}
    (auditdir/'environment.lock.json').write_text(json.dumps(environment,ensure_ascii=False,indent=2),encoding='utf-8')
    (auditdir/'submission_parameters.json').write_bytes((ROOT/'configs/fig_yzdhd_submission_20261014.json').read_bytes())
    (auditdir/'decision_log.md').write_bytes((ROOT/'docs/decisions/2026-10-14-FIG-manuscript-source-correction.md').read_bytes())
    (auditdir/'title_decision_log.md').write_bytes((ROOT/'docs/decisions/2026-10-14-FIG-YZDHD-title-revision.md').read_bytes())
    (OUT/'README_CN.md').write_text('# '+build.TITLE+'\n\n'+'''# 本项目论文投稿包

目标期刊：Functional & Integrative Genomics。正文为益智地黄汤相关候选成分与阿尔茨海默病胶质细胞转录组、模型乳酸化及虚拟敲除研究，基于2026-10-13本项目论文包调整。ADNI论文只提供作者和贡献信息。此前ADNI正文转换包已被本包取代，不应用于本次投稿。

上传正文manuscript/main_manuscript.docx、标题页和投稿信；figures/Fig1–4的PDF/TIFF为出版原图。图源数据、SVG和预览保留用于编辑复核，无需全部重复上传。supplementary/ESM_1.pdf为补充说明，ESM_2.zip包含标准CSV和原始八表，ESM_3.zip为代码与配置。administrative/提供作者录入表、中文摘要和投稿核对清单。audit/保留文献、来源、参数、环境与视觉审计。

正文66篇参考文献均引用，摘要233词、6关键词；2张可编辑表格，4张主图。无资助、无利益冲突已按作者最新确认写入。只改稿件与展示，没有修改原始数据、重跑模型或夸大C0药理证据。公众数据来源与人类/模型推断边界见正文。

本机复用了已有项目源码、冻结结果和实验技术库的引用治理条目，使用已有Python与绘图环境，没有新增安装、环境或远程计算。WikiSkill沿用本项目架构，仅记录纠错和执行轨迹；未演化或晋级新技能，baseline/best/test不适用，接受和拒绝迭代均0。支持可追溯文稿构建与架构复现；不宣称完成全研究跨机器生物学严格复现。原研究代码、环境、配置及科学结果保持。

复建顺序：src/prepare_fig_yzdhd.py、既有figure builder、src/build_fig_yzdhd_submission.py、src/render_fig_yzdhd_submission.py、src/fig_independent_final_reference_audit.py、src/audit_fig_yzdhd_submission.py；最终归档为src/finalize_fig_yzdhd_submission.py。须在完整本地项目中执行，使用记录环境及冻结输入。代码归档没有分发分子原始计数矩阵，不能视为自包含复现包。

文件归档可校验SHA-256，DVC管理输入与成品。临时转换预览在审核完成后清理；来源快照、正式脚本、表格、决策记录、锁文件和最终结果保留。作者仍需确认本机构伦理判定、全部作者批准、数据使用条款和系统合成PDF。本任务没有实际向期刊提交。
''',encoding='utf-8')
    for name in ['language_style_audit.json','language_review_baseline.md','ancillary_language_edits.json']:
        (auditdir/name).write_bytes((PROC/name).read_bytes())
    for name in ['submission_checklist_CN.md']:
        p=admin/name; p.write_text(p.read_text(encoding='utf-8').replace('233词',str(wordcount)+'词'),encoding='utf-8')
    p=OUT/'README_CN.md'; p.write_text(p.read_text(encoding='utf-8').replace('摘要233词','摘要'+str(wordcount)+'词'),encoding='utf-8')
    (auditdir/'language_decision_log.md').write_bytes((ROOT/'docs/decisions/2026-10-14-FIG-YZDHD-language-review.md').read_bytes())
    (auditdir/'language_style_review_CN.md').write_bytes((ROOT/'docs/reports/2026-10-14-FIG-YZDHD-language-review.md').read_bytes())
    for name in ['ingredient_nomenclature_sources_20261015.json']:
        (auditdir/name).write_bytes((PROC/name).read_bytes())
    (auditdir/'ingredient_nomenclature_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-ingredient-nomenclature.md').read_bytes())
    (auditdir/'table_format_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-three-line-tables.md').read_bytes())
    (auditdir/'biological_conclusion_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-biological-conclusion.md').read_bytes())
    (auditdir/'conclusion_template_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-conclusion-template.md').read_bytes())
    (auditdir/'discussion_template_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-discussion-template.md').read_bytes())
    (auditdir/'reference_style_profile_20261015.json').write_bytes((PROC/'reference_style_profile_20261015.json').read_bytes())
    (auditdir/'reference_style_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-reference-style-edit.md').read_bytes())
    (auditdir/'Fig4_missing_values_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-Fig4-estimable-only.md').read_bytes())
    (auditdir/'Fig2c_missing_values_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-Fig2c-estimable-only.md').read_bytes())
    (auditdir/'editorial_technical_check_decision.md').write_bytes((ROOT/'docs/decisions/2026-10-15-FIG-YZDHD-editorial-technical-check.md').read_bytes())
    (admin/'response_to_technical_check.md').write_text('''# Response to the technical check

1. Clinical trial registration

This study is a secondary analysis of public molecular datasets and published experimental records and is not a clinical trial. We have added the exact statement “Clinical trial number: not applicable.” under Clinical trial registration in Statements and Declarations.

2. Registration date and NLRP3

Trial registry details and a registration date are not applicable. NLRP3 denotes a gene investigated in this study, rather than a clinical trial registration number. The manuscript now explicitly states “Clinical trial number: not applicable.”

3. Table S2 caption

We have added the following caption to the Supplementary table captions section of the main manuscript: “Table S2. Senescence-associated and inflammatory transcriptional module responses in amyloid-β-exposed induced pluripotent stem cell-derived microglia.” Captions for all eight supplementary tables are now listed in the main manuscript, with a matching caption index in Online Resource 2. The table data remain unchanged.
''',encoding='utf-8')
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.suffix!='.png')
    manifest=[{'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]
    (PROC/'delivery_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    target=ROOT/'results/Functional_Integrative_Genomics_YZDHD_complete_submission.zip'
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files: z.write(p,'FIG_YZDHD_submission/'+p.relative_to(OUT).as_posix())
        z.writestr('FIG_YZDHD_submission/delivery_manifest.json',json.dumps(manifest,ensure_ascii=False,indent=2))
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        assert all(hashlib.sha256(z.read('FIG_YZDHD_submission/'+p['path'])).hexdigest()==p['sha256'] for p in manifest)
    report={'status':'Passed','archive':str(target),'bytes':target.stat().st_size,'sha256':sha(target),'file_count':len(manifest),'all_members_match_sha256':True,'source':'Correct YZDHD manuscript; ADNI author metadata only'}
    (PROC/'delivery_archive_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__': main()
