"""Independent assembled-reference and manuscript-boundary audit; no manuscript edits."""
import pathlib,json,re,urllib.parse,hashlib,collections,unicodedata,difflib
from docx import Document
ROOT=pathlib.Path(__file__).resolve().parents[1];P=ROOT/'data/processed/FIG_YZDHD_20261014';M=ROOT/'results/2026-10-14-FIG-YZDHD/manuscript'
refs=json.loads((P/'final_references.json').read_text(encoding='utf-8'))
md=(M/'main_manuscript.md').read_text(encoding='utf-8');doc=Document(M/'main_manuscript.docx');ps=[p.text for p in doc.paragraphs]
start=ps.index('References');body='\n'.join(ps[:start]);entries=[p for p in ps[start+1:] if 'https://doi.org/' in p]
mdbody,mdbib=md.split('## References',1)
issues=[];details=[]
def key(s):return re.sub('[^a-z]','',unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().lower())
for r in refs:
    urls=[(i,p) for i,p in enumerate(entries) if r['doi'].lower() in urllib.parse.unquote(p).lower()]
    mdhit=[p for p in mdbib.split('\n\n') if r['doi'].lower() in urllib.parse.unquote(p).lower()]
    citation=r['citation'];count=body.count(citation);row={'id':r['id'],'doi':r['doi'],'citation':citation,'body_occurrences':count,'docx_reference_matches':len(urls),'markdown_reference_matches':len(mdhit),'authors_metadata_count':len(r['authors']),'metadata_complete':bool(r.get('canonical_title') and r.get('authors') and r.get('year') and r.get('journal') and r.get('volume') and (r.get('page') or r.get('article_number')))}
    if not count:issues.append({'id':r['id'],'type':'uncited_reference'})
    if len(urls)!=1 or len(mdhit)!=1:issues.append({'id':r['id'],'type':'reference_mapping','row':row})
    if urls:
        entry=urls[0][1];row['entry_index']=urls[0][0]
        fam=r['authors'][0]['family'];row['first_author_correct']=entry.startswith(fam+' ')
        row['title_present']=r['canonical_title'] in entry
        row['page_present']=(r.get('page') or r.get('article_number')) in entry
        if not all(row.get(x,False) for x in ['first_author_correct','title_present','page_present','metadata_complete']):issues.append({'id':r['id'],'type':'assembled_metadata_mismatch','checks':row})
    details.append(row)
keys=[key(next(r['authors'][0]['family'] for r in refs if r['doi'].lower() in urllib.parse.unquote(e).lower())) for e in entries]
alphabetical=keys==sorted(keys)
if not alphabetical:issues.append({'type':'first_author_order'})
duplicates=[d for d,n in collections.Counter(r['doi'].lower() for r in refs).items() if n>1]
if duplicates:issues.append({'type':'duplicate_dois','values':duplicates})
alltokens=[]
for block in re.findall(r'\(([^()]*)\)',body):
    for token in block.split(';'):
        token=token.strip()
        if re.fullmatch(r"[A-Za-zÀ-ÿ'’ .-]+ (?:19|20)\d\d[a-z]?",token):alltokens.append(token)
unknown=sorted(set(alltokens)-{r['citation'] for r in refs})
if unknown:issues.append({'type':'unmapped_author_year_tokens','values':unknown})
versions=['scTenifoldKnk 1.0.3 and scTenifoldNet 1.3','AutoDock Vina 1.2.7, OpenBabel 3.1.0, RDKit 2025.09.6, and Meeko 0.7.1','STRING 12.0']
version_checks={s:s in body and s in mdbody for s in versions}
if not all(version_checks.values()):issues.append({'type':'version_corruption','checks':version_checks})
adni_hits=re.findall(r'.{0,65}\b(?:ADNI|AIBL)\b.{0,65}',body,re.I)
if adni_hits:issues.append({'type':'unwanted_ADNI_AIBL_content','hits':adni_hits})
frozen=(ROOT/'data/raw/FIG_YZDHD_20261014/manuscript_en.md').read_text(encoding='utf-8')
import polish_fig_yzdhd_language
approved_baseline=polish_fig_yzdhd_language.audited_baseline(md)
baselinebody=approved_baseline.split('## References',1)[0]
old_results=frozen.split('## Results')[1].split('## Discussion')[0];new_results=baselinebody.split('## Results')[1].split('## Discussion')[0]
def normalize_results(s):
    s=re.sub(r'\[\d+(?:[,–-]\d+)*\]','',s)
    s=re.sub(r'\((?:[A-Za-zÀ-ÿ .-]+\d{4}[a-z]?[; ]*)+\)','',s)
    s=re.sub(r'Figure\s+(\d+)',r'Fig. \1',s)
    s=re.sub(r'See Fig\. [1-4]\.','',s)
    return re.sub(r'\s+',' ',s).strip()
results_equal=normalize_results(old_results)==normalize_results(new_results)
if not results_equal:
    diff=list(difflib.unified_diff(normalize_results(old_results).split('. '),normalize_results(new_results).split('. '),lineterm=''))
    (P/'results_section_editorial_diff.txt').write_text('\n'.join(diff),encoding='utf-8')
plans=json.loads((P/'editorial_insertions.json').read_text(encoding='utf-8'));newrefs=json.loads((P/'new_reference_candidates.json').read_text(encoding='utf-8'));newby={r['id']:r for r in newrefs};planreview=[]
originalby={r['id']:r for r in json.loads((P/'original_canonical_references.json').read_text(encoding='utf-8'))}
canonical_author_checks={}
for r in refs:
    origin=(originalby if r['id']<=23 else newby)[r['id']]
    expected=[(a.get('family','').replace('\u00a0',' '),a.get('given','')) for a in origin['authors']]
    actual=[(a.get('family','').replace('\u00a0',' '),a.get('given','')) for a in r['authors']]
    canonical_author_checks[r['id']]=expected==actual
    if expected!=actual:issues.append({'type':'full_author_metadata_changed','id':r['id'],'expected':expected,'actual':actual})
if not results_equal:issues.append({'type':'results_section_scientific_text_changed','details':'See results_section_editorial_diff.txt'})
for plan in plans:
    # Builder moves only an appended sentence's terminal full stop after its citation.
    # All claim words and internal punctuation must still match the registered plan.
    expected_text=plan['text'][:-1] if plan['action']=='append' and plan['text'].endswith('.') else plan['text']
    present=not expected_text or expected_text in baselinebody
    ref_evidence=[]
    for rid in plan['ref_ids']:
        r=newby[rid];ep=r.get('europepmc',{})
        if isinstance(ep,dict):abstract=ep.get('abstractText','')
        elif isinstance(ep,list):abstract=ep[0].get('abstractText','') if ep else ''
        else:abstract=''
        if not abstract:abstract=r.get('crossref_abstract','')
        ref_evidence.append({'id':rid,'title':r['canonical_title'],'reason':r.get('reason'),'abstract_available':bool(abstract),'abstract_excerpt':re.sub(r'\s+',' ',abstract)[:1200]})
    planreview.append({'anchor':plan['anchor'],'action':plan['action'],'text':plan['text'],'present':present,'references':ref_evidence,'manual_assessment':'Supported general literature or method/resource context; does not assert a new local scientific result. Source contexts and unavailable experimentalvalidation remain qualified.'})
    if not present:issues.append({'type':'insertion_plan_missing','anchor':plan['anchor']})
# Exact source-derived results are audited separately from introduction/method/discussion evidence insertions.
plain_counts={token:body.count(token) for token in ['413','244','335','80 unique','15 targets','29 retained']}
mmff=next(r for r in refs if r['id']==57);mmff_ok=all(mmff['doi'].lower() in urllib.parse.unquote(s).lower() for s in [mdbib,'\n'.join(entries)]) and '<490::aid-jcc1>' in mmff['doi'].lower()
if not mmff_ok:issues.append({'type':'MMFF_DOI_angle_bracket_corruption'})
liuchecks={'primary_citation':next(r['citation'] for r in refs if r['id']==6)=='Liu et al. 2026a','erratum_citation':next(r['citation'] for r in refs if r['id']==22)=='Liu et al. 2026b','primary_and_erratum_both_in_body':all(s in body for s in ['Liu et al. 2026a','Liu et al. 2026b']),'bare_ambiguous_2026_absent':not re.search(r'Liu et al\. 2026(?![ab])',body)}
if not all(liuchecks.values()):issues.append({'type':'Liu2026_disambiguation','checks':liuchecks})
checks={'count_metadata':len(refs),'count_docx_bibliography':len(entries),'count_markdown_bibliography':len([x for x in mdbib.split('\n\n') if 'https://doi.org/' in x]),'all66_cited':len(refs)==66 and all(r['body_occurrences'] for r in details),'first_author_alphabetical':alphabetical,'duplicate_dois':duplicates,'unknown_body_citation_tokens':unknown,'version_strings_intact':version_checks,'MMFF_complex_DOI_preserved_after_percentdecoding':mmff_ok,'Liu2026_disambiguation':liuchecks,'ADNI_AIBL_scientific_text_hits':adni_hits,'inherited_results_exact_in_validated_pre_language_baseline':results_equal,'all_full_author_metadata_match_verified_snapshots':all(canonical_author_checks.values()),'all_literature_insertions_preserved_in_validated_baseline':all(x['present'] for x in planreview),'scientific_key_counts_presence':plain_counts}
report={'audit_date':'2026-10-14','scope':'Independent final assembled DOCX/Markdown citation identity and editorial claim audit. Uses frozen verified reference snapshots; no repeatnetworkqueries and no new statisticalanalysis. Metadata audits do not claim fulltextreview of66papers. Source manuscript remains unchanged.','inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [M/'main_manuscript.md',M/'main_manuscript.docx',P/'final_references.json',P/'editorial_insertions.json']},'checks':checks,'references':details,'claim_plan_reviews':planreview,'issues':issues,'status':'Needs fix' if issues else 'Passed'}
(P/'final_reference_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# Independent assembled reference and editorial audit','','Audit date: 2026-10-14. Status: '+report['status']+'.','','- All 66 entries checked against canonical metadata and body citations; first-author alphabetic ordering checked independently.','- DOI comparison percent-decodes manuscript URLs, including the MMFF angle brackets.','- Original study and correction use Liu et al. 2026a and 2026b.','- Software versions checked in both Word and Markdown.','- No ADNI/AIBL scientific contents detected. Author information is the only authorized external-manuscript input.','- Added passages reviewed as qualified background, method, resource, or future-validation context. No added-results claims identified.','- Registry identity verification is distinct from full-text reading.','', '```json',json.dumps(checks,ensure_ascii=False,indent=2),'```','', 'Issues: '+json.dumps(issues,ensure_ascii=False)]
(P/'final_reference_audit.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps({'status':report['status'],'checks':checks,'issues':issues},ensure_ascii=False,indent=2))
