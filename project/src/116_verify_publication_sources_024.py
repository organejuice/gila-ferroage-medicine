"""Legal public supplement attempts and exact historical activity rechecks."""
from pathlib import Path
import importlib.util
import json

s=importlib.util.spec_from_file_location('s024',Path(__file__).with_name('113_supplement_formula_targets_024.py'))
m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
cfg=json.loads(m.CONFIG.read_text(encoding='utf-8')); out=m.ROOT/'results'/cfg['run_id']; raw=m.ROOT/'data/raw'/cfg['run_id']
records=[]; client=m.Client('publication',cfg,raw)
urls=['https://journals.sagepub.com/doi/10.1177/15491684261460820']+[f'https://journals.sagepub.com/doi/suppl/10.1177/15491684261460820/suppl_file/sj-docx-{i}-rej-10.1177_15491684261460820.docx' for i in range(1,4)]
for url in urls:
    try:
        body,meta=client.get(url)
        records.append({**meta,'review_status':'DOWNLOADED_NOT_YET_FORMULA_VERIFIED'})
    except m.SourceError as exc:
        records.append({'url':url,'retrieved_at':m.now(),'status':'SOURCE_BLOCKED','detail':str(exc),'review_status':'NOT_OBTAINED'})
    m.writetable(out/'publication_supplement_access.tsv',records)
genes={r['gene'] for r in m.readtable(m.ROOT/cfg['inputs']['baseline_023'])}
hist=[r for r in m.readtable(m.ROOT/cfg['inputs']['historical_assays']) if r['Target_Symbol'] in genes and r['Evidence_Source']=='ChEMBL']
client=m.Client('historical_ChEMBL',cfg,raw); result=[]
for aid in sorted({r['Activity_ID'] for r in hist if r['Activity_ID']}):
    try:
        data,meta=client.json(cfg['sources']['ChEMBL']+'activity/'+aid+'.json')
        molecule,mm=client.json(cfg['sources']['ChEMBL']+'molecule/'+data['molecule_chembl_id']+'.json')
        target,tm=client.json(cfg['sources']['ChEMBL']+'target/'+data['target_chembl_id']+'.json')
        actual=m.canonical(molecule.get('molecule_structures',{}).get('canonical_smiles',''))
        queries=m.readtable(out/'pubchem_structure_resolution.tsv')
        expected={r['PubChem_CID']:r.get('canonical_smiles','') for r in queries}
        candidate_ids=sorted({r['PubChem_CID'] for r in hist if r['Activity_ID']==aid})
        result.append({'activity_id':aid,'actual_chembl_molecule':data.get('molecule_chembl_id',''),'actual_target':data.get('target_chembl_id',''),
                       'actual_canonical_smiles':actual,'candidate_CIDs':';'.join(candidate_ids),
                       'exact_candidate_CIDs':';'.join(cid for cid in candidate_ids if expected.get(cid)==actual),
                       'molecule_response_sha256':mm['sha256'],'target_response_sha256':tm['sha256'],
                       'target_accessions':';'.join(x.get('accession','') for x in target.get('target_components',[])),
                       'actual_assay':data.get('assay_chembl_id',''),'source_status':'VERIFIED_PUBLIC_ACTIVITY','response_sha256':meta['sha256'],'source_url':meta['url']})
    except m.SourceError as exc: result.append({'activity_id':aid,'source_status':'SOURCE_BLOCKED','detail':str(exc)})
    m.writetable(out/'historical_ChEMBL_original_activity_verification.tsv',result)
print(json.dumps({'publication':len(records),'historical_activity_records':len(result)},ensure_ascii=False))
