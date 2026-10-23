"""Freeze the correct manuscript and author-only metadata for journal conversion."""
from pathlib import Path
import hashlib, json, shutil, csv
from docx import Document

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/raw/FIG_YZDHD_20261014'
OUT=ROOT/'results/2026-10-14-FIG-YZDHD'
PROC=ROOT/'data/processed/FIG_YZDHD_20261014'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    for p in [RAW,OUT,PROC,ROOT/'workspaces/2026-10-14-FIG-YZDHD/raw/traces']:
        p.mkdir(parents=True,exist_ok=True)
    base=ROOT/'results/2026-10-13-manuscript'
    sources=[base/p for p in ['manuscript_en.md','claim_evidence.tsv','content_audit.json','provenance.json','editorial_audit_zh.md']]
    sources+=sorted((base/'supplementary_tables').glob('*'))
    manifest=[]
    for p in sources:
        dest=RAW/p.relative_to(base)
        dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists(): shutil.copy2(p,dest)
        assert sha(p)==sha(dest),str(p)
        manifest.append({'source':str(p),'frozen_copy':str(dest.relative_to(ROOT)),'sha256':sha(p)})
    author_source=Path('E:/我的文档/论文'+chr(36)*4+'/生信论文/ADNI数据库/2026-09-17_alzheimers_research_and_therapy_submission_ready/manuscript/main_manuscript - 1.docx')
    paragraphs=[p.text for p in Document(author_source).paragraphs]
    contributions=paragraphs[paragraphs.index("Authors' contributions")+1]
    title='Glial lactylation, senescence and iron homeostasis in Alzheimer’s disease: an integrative transcriptomic and network pharmacology analysis of Yizhi Dihuang Decoction'
    # Reuse verified author metadata, independently of the deleted ADNI package.
    authors=json.loads((PROC/'author_metadata.json').read_text(encoding='utf-8'))['authors']
    assert authors[1]['email']=='author2@example.invalid'
    assert authors[7]['email']=='author8@example.invalid'
    metadata={'journal':'Functional & Integrative Genomics','article_type':'Research','title':title,'authors':authors,'source_author_paragraphs':paragraphs[1:5],'source_contribution_text':contributions,'contributions':'Hengjun Liu: Conceptualization, Methodology, Software, Formal analysis, Data curation, Visualization and Writing – original draft. Shi Jian, TingTing Zhang, Shangyi Cai, Yunze Shi, Zhengtong Xu, Yude Wang and Yu Zhang: Validation, Investigation and Writing – review and editing. Rui Qie: Conceptualization, Supervision, Project administration and Writing – review and editing.','funding':'The authors received no funding for this work.','competing_interests':'The authors declare no competing interests.','funding_conflicts_source':'User confirmation on 2026-10-14: 无资助、无利益冲突','author_source':str(author_source),'author_source_sha256':sha(author_source),'author_source_boundary':'Only author names, affiliations, correspondence and contribution roles are reused. No ADNI scientific content, cohort data, acknowledgements, grants or ethics declarations are transferred.'}
    (PROC/'author_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    (PROC/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    settings={'source_manuscript':'data/raw/FIG_YZDHD_20261014/manuscript_en.md','output_dir':str(OUT.relative_to(ROOT)),'article_type':'Research','references_minimum':61,'source_boundary':'Project manuscript for scientific content; designated ADNI document for author metadata only','biological_reanalysis':False,'new_environment':False,'existing_results_unchanged':True,'font':'Times New Roman','body_points':11,'spacing':1.25,'page_size':'A4 inherited manuscript design','author_year_citations':True,'keywords':6,'figures_separate':True,'submission_to_portal':False}
    settings['title']=title
    (ROOT/'configs/fig_yzdhd_submission_20261014.json').write_text(json.dumps(settings,ensure_ascii=False,indent=2),encoding='utf-8')
    trace={'run_id':'2026-10-14-FIG-YZDHD','backend':'Current Codex model-led editorial work, existing project WikiSkill architecture reused','correction':'Prior package wrongly used ADNI body; superseded for this user task. Correct source is this project YZDHD manuscript.','baseline_validation':None,'best_validation':None,'test':None,'accepted_iterations':0,'rejected_iterations':0,'active_skills_added':[],'cross_project_promotion':False,'biological_reanalysis':False,'deterministic_checks':'Source hashes and author extraction passed; final package audit pending'}
    tracepath=ROOT/'workspaces/2026-10-14-FIG-YZDHD/raw/traces/source_scope_correction.json'
    if not tracepath.exists(): tracepath.write_text(json.dumps(trace,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'source_files_frozen':len(manifest),'author_count':len(authors),'correct_project':True},ensure_ascii=False))
if __name__=='__main__': main()
