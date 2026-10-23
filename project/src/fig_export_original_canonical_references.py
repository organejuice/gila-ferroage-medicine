"""Export corrected author-year bibliography metadata, keeping original audit separate."""
import json,pathlib
root=pathlib.Path(__file__).resolve().parents[1];out=root/'data/processed/FIG_YZDHD_20261014'
audit=json.loads((out/'original_reference_audit.json').read_text(encoding='utf-8'))
refs=[]
for item in audit['references']:
    x=audit['recommended_replacement'] if item['id']==18 else item
    row={k:x.get(k,'') for k in ['canonical_title','year','journal','volume','issue','page','article_number','doi','sources']}
    row['id']=item['id'];row['authors']=[{'family':a.get('family','').replace('\u00a0',' '),'given':a.get('given',''),'collective':a.get('name','')} for a in x['authors']];row['status']='Verified'
    if item['id']==22:row['publication_type']='Published Erratum';row['display_title']='Correction: '+row['canonical_title']
    if item['id']==18:row['change_reason']=audit['recommended_replacement']['rationale']
    refs.append(row)
(out/'original_canonical_references.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8')
assert len(refs)==23 and len({r['doi'].lower() for r in refs})==23
print('23 corrected canonical references exported; STRING18 now2023 version12.0 source.')
