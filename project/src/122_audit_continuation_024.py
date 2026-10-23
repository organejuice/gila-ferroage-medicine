"""Audit continuation snapshots and offline replay without touching baseline."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

spec=importlib.util.spec_from_file_location('c024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m


def main():
    if (c.RAW/'Publication_user_supplied').exists():
        subprocess.run([sys.executable,'-X','utf8','src/124_ingest_publication_supplements_024.py'],cwd=m.ROOT,check=True)
    snapshots=[]
    for path in sorted(c.RAW.rglob('*.json')):
        meta=json.loads(path.read_text(encoding='utf-8'))
        if 'body_path' not in meta: continue
        body=m.ROOT/meta['body_path']
        actual=m.filehash(body) if body.exists() else ''
        snapshots.append(dict(metadata=str(path.relative_to(m.ROOT)),source=meta['source'],url=meta['url'],declared_sha256=meta['sha256'],actual_sha256=actual,match=actual==meta['sha256'],bytes=meta['bytes']))
    m.writetable(c.OUT/'raw_source_verification.tsv',snapshots)
    tests=subprocess.run([sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_formula*024.py','-v'],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8')
    (c.OUT/'deterministic_tests.txt').write_text(tests.stdout+tests.stderr,encoding='utf-8')
    outputs=['evidence_layered.tsv','G1_244_layered_connections.tsv','exploratory_disease_connections.tsv','G1_G5_gate_layer_summary.tsv','historical_8_and_PTGS2_disposition.tsv']
    outputs += [name for name in ['publication_supplement_table1.tsv','publication_supplement_table2.tsv','publication_gene_disease_comparison.tsv','publication_G1_G5_summary.tsv'] if (c.OUT/name).exists()]
    hashes=[]
    for iteration in range(2):
        proc=subprocess.run([sys.executable,'-X','utf8','src/121_finalize_continuation_024.py'],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8',check=True)
        hashes.append({name:m.filehash(c.OUT/name) for name in outputs})
    replay=[dict(path=name,first_sha256=hashes[0][name],second_sha256=hashes[1][name],match=hashes[0][name]==hashes[1][name]) for name in outputs]
    m.writetable(c.OUT/'offline_replay_verification.tsv',replay)
    state=json.loads((c.OUT/'machine_audit.json').read_text(encoding='utf-8'))
    baseline=[]
    for row in m.readtable(c.OLD/'frozen_input_verification.tsv'):
        baseline.append(dict(path=row['path'],expected=row['frozen_sha256'],actual=m.filehash(m.ROOT/row['path'])))
    checks={'raw_snapshots_match':all(r['match'] for r in snapshots),'tests_pass':tests.returncode==0,
        'offline_replay_match':all(r['match'] for r in replay),'frozen_inputs_match':all(r['actual']==r['expected'] for r in baseline),
        'C0_preserved':state['C0'] and not state['reopen_017'],'G1_original_review_pending':state['E1_G1_assays_review_pending'],
        'all_sources_completed':state['status']=='COMPLETE','stage_status':state['status'],
        'DVC_remote_restore':'NOT_CONFIGURED','environment_rebuilt':False}
    m.writejson(c.OUT/'acceptance_audit.json',checks)
    files=['src/118_continue_formula_supplement_024.py','src/119_query_herb_candidates_024.py','src/120_resolve_legacy_structures_024.py','src/121_finalize_continuation_024.py','src/122_audit_continuation_024.py','src/123_supervise_continuation_024.py','tests/test_formula_continuation_024.py','configs/fixed_formula_supplement_024_continuation.json','docs/decisions/2026-10-03-supplement-024-continuation.md']
    files.append('src/124_ingest_publication_supplements_024.py')
    m.writetable(c.OUT/'implementation_sha256.tsv',[dict(path=name,sha256=m.filehash(m.ROOT/name)) for name in files])
    print(json.dumps(checks,ensure_ascii=False),flush=True)
    if not all(checks[k] for k in ['raw_snapshots_match','tests_pass','offline_replay_match','frozen_inputs_match','C0_preserved']): raise SystemExit(1)


if __name__=='__main__': main()
