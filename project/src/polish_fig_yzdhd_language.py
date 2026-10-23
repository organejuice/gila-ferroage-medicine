"""Reviewed editorial substitutions; all scientific inputs remain unchanged."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]
PROC=ROOT/'data/processed/FIG_YZDHD_20261014'
CONFIG=ROOT/'configs/fig_yzdhd_language_edits_20261014.json'

def apply(text):
    baseline=PROC/'language_review_baseline.md'
    if not baseline.exists(): baseline.write_text(text,encoding='utf-8')
    assert baseline.read_text(encoding='utf-8')==text,'Unreviewed upstream manuscript change'
    edits=json.loads(CONFIG.read_text(encoding='utf-8'))
    body,bib=text.split('## References\n',1)
    records=[]
    for e in edits:
        count=body.count(e['before'])
        assert count==e['occurrences'],(e['before'],count,e['occurrences'])
        body=body.replace(e['before'],e['after'])
        records.append(e)
    result=body+'## References\n'+bib
    citations=lambda s: re.findall(r'\((?:[A-Za-zÀ-ÿ’\' .-]+\d{4}[a-z]?[; ]*)+\)',s)
    numeric=lambda s: re.findall(r'(?<![A-Za-z])[-−]?\d[\d,.]*(?:%|\^[-−]?\d+)?',re.sub(r'G[1-5]|D[378]|M3|GSE\d+|Online Resource [1-3]','',s))
    report={'status':'Passed','scope':'Reviewed paper-style substitutions; no new data, methods, numerical results or reference entries',
            'edit_rules':len(records),'edits':records,'baseline_sha256':hashlib.sha256(text.encode()).hexdigest(),
            'polished_sha256':hashlib.sha256(result.encode()).hexdigest(),
            'bibliography_identical':bib==result.split('## References\n',1)[1],
            'citation_sequence_identical':citations(text)==citations(result),
            'numerical_token_sequence_identical_excluding_internal_codes_accessions_resource_labels':numeric(text)==numeric(result)}
    assert all(report[k] for k in ['bibliography_identical','citation_sequence_identical','numerical_token_sequence_identical_excluding_internal_codes_accessions_resource_labels']),{k:v for k,v in report.items() if k!='edits'}
    (PROC/'language_style_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return result

def audited_baseline(polished):
    """Validate the complete approved transform before auditing inherited source prose."""
    baseline=(PROC/'language_review_baseline.md').read_text(encoding='utf-8')
    assert apply(baseline)==polished,'Unexpected unreviewed manuscript edits'
    return baseline
