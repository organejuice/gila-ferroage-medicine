"""Audit immutable snapshots, input provenance, offline replay and deliverables."""
from pathlib import Path
import importlib.util
import json
import subprocess
import sys
import yaml

spec=importlib.util.spec_from_file_location('s024',Path(__file__).with_name('113_supplement_formula_targets_024.py'))
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
cfg=json.loads(m.CONFIG.read_text(encoding='utf-8')); OUT=m.ROOT/'results'/cfg['run_id']; RAW=m.ROOT/'data/raw'/cfg['run_id']

def main():
    rawrows=[]
    for p in sorted(RAW.rglob('*.json')):
        data=json.loads(p.read_text(encoding='utf-8'))
        if 'body_path' not in data: continue
        body=m.ROOT/data['body_path']
        actual=m.filehash(body) if body.exists() else ''
        rawrows.append({'metadata_path':str(p.relative_to(m.ROOT)),'source':data['source'],'url':data['url'],
          'retrieved_at':data['retrieved_at'],'HTTP_status':data['status'],'body_path':data['body_path'],
          'declared_sha256':data['sha256'],'actual_sha256':actual,'match':actual==data['sha256'],'bytes':data['bytes']})
    m.writetable(OUT/'raw_source_verification.tsv',rawrows)
    extras=['configs/g2_v3.2_execution.yaml','configs/g3_v3.2_execution.yaml','configs/g4_v3.2_execution.yaml',
            'src/32_run_g4_v3_2.py','results/archive/2026-09-20-pre-v3-freeze/results/tables/D5_ferroaging_reference.tsv']
    manifest=OUT/'supplementary_disease_input_sha256.tsv'
    if not manifest.exists(): m.writetable(manifest,[{'path':p,'sha256':m.filehash(m.ROOT/p),'bytes':(m.ROOT/p).stat().st_size} for p in extras])
    extras_match=all(m.filehash(m.ROOT/r['path'])==r['sha256'] for r in m.readtable(manifest))
    result=subprocess.run([sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_formula_target_supplement_024.py','-v'],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8')
    (OUT/'deterministic_tests.txt').write_text(result.stdout+result.stderr,encoding='utf-8')
    outputs=['evidence_layered.tsv','G1_244_layered_connections.tsv','exploratory_disease_connections.tsv','G1_layer_summary.tsv','historical_8_and_PTGS2_disposition.tsv']
    hashes=[]
    for _ in range(2):
        subprocess.run([sys.executable,'-X','utf8','src/114_finalize_formula_targets_024.py'],cwd=m.ROOT,check=True,stdout=subprocess.DEVNULL)
        hashes.append({p:m.filehash(OUT/p) for p in outputs})
    replay=[{'path':p,'first_sha256':hashes[0][p],'second_sha256':hashes[1][p],'match':hashes[0][p]==hashes[1][p]} for p in outputs]
    m.writetable(OUT/'offline_replay_verification.tsv',replay)
    state=json.loads((OUT/'machine_audit.json').read_text(encoding='utf-8'))
    checks=[
      {'criterion':'all_query_rows_have_status','status':'PASS' if state['query_rows']==state['query_rows_with_explicit_status'] else 'FAIL'},
      {'criterion':'all_sources_completed','status':'PARTIAL_SOURCE_BLOCKED' if state['status']!='COMPLETE' else 'PASS'},
      {'criterion':'frozen_biological_inputs_match','status':'PASS' if state['all_biological_inputs_hash_unchanged'] and extras_match else 'FAIL'},
      {'criterion':'raw_snapshot_hashes_match','status':'PASS' if all(r['match'] for r in rawrows) else 'FAIL'},
      {'criterion':'deterministic_rule_tests','status':'PASS' if result.returncode==0 else 'FAIL'},
      {'criterion':'offline_target_and_intersection_replay','status':'PASS' if all(r['match'] for r in replay) else 'FAIL'},
      {'criterion':'new_E1_G1_original_assay_reviews','status':'PASS' if state['E1_G1_assays_review_pending']==0 else 'FAIL'},
      {'criterion':'historical_8_and_PTGS2_disposition','status':'PASS' if state['historical_disposition_rows']>8 else 'FAIL'},
      {'criterion':'C0_and_017_unchanged','status':'PASS' if state['C0'] and not state['reopen_017'] else 'FAIL'},
      {'criterion':'paper_supplement_and_formula_equivalence','status':'NOT_OBTAINED_SOURCE_BLOCKED'},
      {'criterion':'DVC_remote_restore','status':'NOT_CONFIGURED'}]
    m.writetable(OUT/'acceptance_checklist.tsv',checks)
    tracked=['src/113_supplement_formula_targets_024.py','src/114_finalize_formula_targets_024.py','src/115_fetch_formula_prior_sources_024.py','src/116_verify_publication_sources_024.py','src/117_audit_formula_supplement_024.py',
             'tests/test_formula_target_supplement_024.py','configs/fixed_formula_target_supplement_024.json','configs/fixed_formula_validated_rules_024.json','configs/fixed_formula_activity_comment_rules_024.json','configs/TCM_select_024_pip_freeze.txt',
             'configs/wikiskill_fixed_formula_supplement_024.json','docs/protocol/2026-09-30-fixed-formula-target-supplement-SOP-024.md',
             'docs/protocol/2026-09-30-fixed-formula-target-supplement-amendment-024.md','docs/decisions/2026-09-30-fixed-formula-target-supplement-024.md',
             str((OUT/'manual_assay_adjudication.tsv').relative_to(m.ROOT))]
    m.writetable(OUT/'final_implementation_sha256.tsv',[{'path':p,'sha256':m.filehash(m.ROOT/p),'bytes':(m.ROOT/p).stat().st_size} for p in tracked])
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=m.ROOT,capture_output=True,text=True).stdout.strip()
    m.writejson(OUT/'execution_provenance.json',{'Git_HEAD':head,'code_changes_uncommitted':True,'reproducibility':'local offline replay verified; remote restoration unavailable; not full strict reproducibility',
          'environment':'existing TCM_select; dependencies not installed','resource':'local Windows CPU; public data HTTP retrieval only',
          'rule_tests_passed':result.returncode==0,'snapshot_records':len(rawrows),'offline_replay_equal':all(r['match'] for r in replay),
          'temporary_workspace_test_cache_removed':not (m.ROOT/'.cache/024-restart-test').exists(),
          'outside_workspace_temporary_test_cleanup_unverified':True})
    print(json.dumps({'raw_snapshots':len(rawrows),'tests':result.returncode,'replay_equal':all(r['match'] for r in replay),'acceptance':checks},ensure_ascii=False))

if __name__=='__main__': main()
