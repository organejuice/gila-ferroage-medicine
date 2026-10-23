"""Finalize field decisions from the original-reference metadata snapshot."""
import json,pathlib,re,urllib.request,urllib.parse,hashlib
ROOT=pathlib.Path(__file__).resolve().parents[1];out=ROOT/'data/processed/FIG_YZDHD_20261014'
rows=json.loads((out/'original_reference_metadata_snapshot.json').read_text(encoding='utf-8'))
def norm(s):return re.sub(r'[^a-z0-9]','',s.lower())
def initials(given):return ''.join(x[0] for x in re.findall(r'[A-Za-zÀ-ÿ]+',given))
final=[]
for row in rows:
    r={k:row[k] for k in ['id','original_entry','doi','canonical_title','authors','year','journal','volume','issue','page','article_number','dates']}
    ep=row['europepmc'][0] if row['europepmc'] else None
    r['sources']=[{'source':'Crossref','url':row['crossref_url'],'status':'resolved'}]
    if ep:r['sources'].append({'source':'PubMed via EuropePMC','url':'https://pubmed.ncbi.nlm.nih.gov/'+ep['id']+'/','retrieval_url':row['europepmc_url'],'pmid':ep['id'],'title':ep['title'],'author_string':ep.get('authorString'),'status':'doi matched'})
    r['issues']=[]
    original_title=row['original_entry'].split('. ',1)[1].split('. '+row['journal'])[0]
    if row['id']==17:original_title=original_title.split('. Proceedings')[0]
    if row['id']==22:original_title=original_title.removeprefix('Correction: ')
    if norm(original_title)!=norm(row['canonical_title']):r['issues'].append({'severity':'check','field':'title','message':'Title parser differs; compare manually.','original_parsed':original_title})
    if row['issue'] and not re.search(re.escape(str(row['volume']))+r'\('+re.escape(row['issue'])+r'\)',row['original_entry']):r['issues'].append({'severity':'info','field':'issue','message':'Add canonical issue for completeness; identity unaffected.'})
    if ep:
        ji=ep.get('journalInfo',{})
        r['registry_comparison']={'doi_match':ep.get('doi','').lower()==row['doi'].lower(),'title_match':norm(ep['title'])==norm(row['canonical_title']),'year_match':str(ep.get('pubYear'))==str(row['year']),'volume_match':ji.get('volume')==row['volume'],'issue_match':ji.get('issue')==row['issue'],'page_match':ep.get('pageInfo') in [row['page'],row['article_number'],row['page'].split('-')[0]]}
        r['comment_corrections']=ep.get('commentCorrectionList',{}).get('commentCorrection',[])
    if row['id']==7:r['page']='e47';r['article_number']='e47';r['issues'].append({'severity':'info','field':'page','message':'Crossref e47-e47 normalized to article number e47, concordant with PubMed.'})
    if row['id']==11:
        r['sources'].append({'source':'Elsevier ScienceDirect original1997 publisher record (web indexed preview)','url':'https://www.sciencedirect.com/science/article/pii/S0169409X96004231','status':'title authors year volume issue pages DOI concordant; direct fetch403'})
        r['version_note']='Retain original1997 DOI. PubMed-indexed2001 reprint DOI10.1016/S0169-409X(00)00129-0 has volume46/pages3-26; do not mix metadata.'
    if row['id']==6:
        r['publication_type']='Journal Article';r['version_note']='Real2026 publication, before audit2026-10-14. Linked erratum22 concerns mismatched FigureS2P source data, not publishedfigure or analyses/conclusions.'
    if row['id']==18:
        r['issues'].append({'severity':'context','field':'database_version','message':'Identity is valid2025 database paper, but frozen STRING12.0 is explicitly covered by2023 paper10.1093/nar/gkac1000. Prefer replace18 with2023 source; retain version12.0 methods.'})
    if row['id']==22:
        r['publication_type']='Published Erratum';r['display_title']='Correction: '+r['canonical_title']
        r['sources'].append({'source':'Live PubMed publisher-corrected issue record','url':'https://pubmed.ncbi.nlm.nih.gov/42617592/','status':'2026Oct6;38(10):2146; online2026Aug19'})
        r['sources'].append({'source':'Elsevier ScienceDirect correction notice (webindexed)','url':'https://www.sciencedirect.com/science/article/abs/pii/S1550413126003360','status':'source data FigureS2P corrected; publishedfigure analyses conclusions unchanged','corrected_source_data':'https://doi.org/10.5281/zenodo.21806776'})
        r['version_note']='Canonical registries omit prefix Correction:. Prefix or bracketed [Correction] is a transparent publication-type label, not title identity error. Cached EuropePMC lacks finalissue/pages; currentCrossref and livePubMed concordant.'
        r['issues'].append({'severity':'info','field':'registry_update','message':'Use final38(10):2146 from currentCrossref+livePubMed, not cachedEuropePMC inpress PII.'})
    r['formatted_author_year_entry']=', '.join(a['family'].replace('\u00a0',' ')+' '+initials(a.get('given','')) for a in r['authors'])+' ('+str(r['year'])+') '+r.get('display_title',r['canonical_title'])+'. '+r['journal']+' '+r['volume']+('('+r['issue']+')' if r['issue'] else '')+':'+(r['page'] or r['article_number'])+'. https://doi.org/'+urllib.parse.quote(r['doi'],safe='/():')
    r['status']='Verified; context correction recommended' if row['id']==18 else 'Verified'
    final.append(r)
# DOI-specific replacement; never substitute newestversion for frozenanalysis.
doi='10.1093/nar/gkac1000';url='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
with urllib.request.urlopen(url,timeout=50) as res:m=json.load(res)['message']
replacement={'replaces_original_id':18,'doi':doi,'canonical_title':m['title'][0],'authors':m['author'],'year':2023,'journal':m['container-title'][0],'volume':m['volume'],'issue':m['issue'],'page':m['page'],'article_number':m.get('article-number',''),'sources':[url,'https://pubmed.ncbi.nlm.nih.gov/36370105/','https://version-12.string-db.org/cgi/access'],'rationale':'2023 primarypaper explicitly describes STRING12.0 used in frozenanalysis; no networkversionchange.'}
(out/'string12_reference_replacement.json').write_text(json.dumps(replacement,ensure_ascii=False,indent=2),encoding='utf-8')
report={'audit_date':'2026-10-14','scope':'23 originalYZDHD references; ADNI scientificcontents excluded. Registry metadata audit; not fulltextreviewofall23. Correction22 notice and original1997publisherpreview reviewed.','input_sha256':hashlib.sha256((ROOT/'docs/manuscript/manuscript_en.md').read_bytes()).hexdigest(),'reuse':'Local referenceidentity/retraction/versiongovernance; existing Crossref+EuropePMC workflow. No environmentinstallation or analysisrerun.','summary':{'original_entries':23,'crossref_resolved':23,'doi_matched_pubmed':22,'publisher_fallback':1,'identity_critical_errors':0,'context_correction_recommended':[18],'known_erratum_pair':[6,22]},'references':final,'recommended_replacement':replacement}
(out/'original_reference_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['# Original YZDHD reference verification','', 'Auditdate2026-10-14. All23 DOIidentities verified: Crossref23/23, PubMed/EuropePMC22/23, original1997Lipinski publisherrecord1/23. No DOI points to anotherpaper. Truncated etal bibliographies expanded with fullcanonical authors.','', 'Metadataaudit does not claim fulltextreading ofall23 papers. Frozen analysis versions and results unchanged.','', '## Decisions','', '- Ref18: replace2025general STRINGpaper with2023 STRING12.0paper DOI10.1093/nar/gkac1000; methods retain12.0.', '- Refs6/22: valid2026article+published erratum. Finalerratum38(10):2146; cachedEuropePMC remains inpress. CorrectionconcernsFigureS2Psource-data compilation; publisherstates figure/analyses/conclusions unchanged.', '- Ref11: retain1997originalDOI and23(1–3):3–25. Do not combine with2001reprint metadata.', '- Completeissue fields throughout; normalize limmae47-e47 toarticle e47. Correctionprefix is descriptivepublicationtypelabel.','', '## Entries','', '|OriginalID|Status|Canonicalfirstauthor/year|Volume(issue):pages/article|DOI|','|---|---|---|---|---|']
for r in final:lines.append(f"|{r['id']}|{r['status']}|{r['authors'][0]['family']} {r['year']}|{r['volume']}({r['issue']}):{r['page'] or r['article_number']}|[{r['doi']}](https://doi.org/{urllib.parse.quote(r['doi'],safe='/():')})|")
(out/'original_reference_audit.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(report['summary'],indent=2));print('Unexpectedtitleissues',[(r['id'],r['issues']) for r in final if any(x['field']=='title' for x in r['issues'])])
