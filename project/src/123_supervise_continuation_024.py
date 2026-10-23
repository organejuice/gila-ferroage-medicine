"""Keep local 024 jobs observable and finalize/audit their outputs on completion."""
import importlib.util
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

spec=importlib.util.spec_from_file_location('c024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m
PID_FILE=c.OUT/'worker_processes.json'


def running(ids):
    if not ids: return []
    cmd="Get-CimInstance Win32_Process -Filter \"Name = 'python.exe'\" | Where-Object { $_.ProcessId -in @("+','.join(str(i) for i in ids)+") -and $_.CommandLine -match 'src[\\\\/](118_continue_formula_supplement_024|119_query_herb_candidates_024|120_resolve_legacy_structures_024)' } | Select-Object -ExpandProperty ProcessId"
    result=subprocess.run(['powershell','-NoProfile','-Command',cmd],capture_output=True,text=True,check=True)
    return [int(line.strip()) for line in result.stdout.splitlines() if line.strip().isdigit()]


def ledger_progress(name,expected):
    path=c.OUT/name
    rows=m.readtable(path) if path.exists() else []
    return dict(expected=expected,recorded=len(rows),success=sum(r['status'].startswith('SUCCESS') for r in rows),statuses=dict(Counter(r['status'] for r in rows)))


def snapshot(active,phase='QUERYING'):
    pp=c.OUT/'publication_supplement_status.json'
    paper=json.loads(pp.read_text(encoding='utf-8'))['status'] if pp.exists() else 'SOURCE_BLOCKED'
    state=dict(run_id=c.RUN,phase=phase,updated_at=m.now(),active_worker_PIDs=active,
        BindingDB=ledger_progress('ledger_BindingDB.tsv',2918),HERB=ledger_progress('HERB_CID_query_ledger.tsv',2276),
        legacy_identity=ledger_progress('legacy_structure_identity_followup.tsv',642),paper_supplements=paper,C0=True,reopen_017=False)
    m.writejson(c.OUT/'progress.json',state)
    lines=['# 024续跑实时进度','',f'状态：{phase}；更新时间（UTC）：{state["updated_at"]}。此文件由本机流程自动更新。','',
        '|任务|已成功查询/核验|预期总数|','|---|---:|---:|']
    for key in ['BindingDB','HERB','legacy_identity']:
        value=state[key]; lines.append(f'|{key}|{value["success"]}|{value["expected"]}|')
    lines+=['','论文三个补充文件已由用户提供并解析，处方等价性与候选成分—靶点对仍未建立。' if paper=='USER_SUPPLIED_PARSED' else '论文3个补表原链接仍受阻，Figshare DOI检索未命中。失败、未查与身份待定不记成阴性。','',
        '查询结束后自动执行122审计：合并主G1/探索G2–G5连接、校验原始快照、规则测试及两次离线重放。最终文件为report.md、acceptance_audit.json和G1_G5_gate_layer_summary.tsv。',
        '运行中report.md可能是较早快照，以本实时进度为准；只有phase=FINISHED时才进入最终验收。药物仍C0，017关闭范围保持。']
    (c.OUT/'progress.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return state


def main():
    jobs=json.loads(PID_FILE.read_text(encoding='utf-8'))
    if isinstance(jobs,dict): jobs=[jobs]
    ids=[int(r['ProcessId']) for r in jobs]
    while True:
        active=running(ids)
        try: snapshot(active)
        except (ValueError,KeyError,json.JSONDecodeError):
            time.sleep(2); continue
        if not active: break
        time.sleep(30)
    snapshot([],phase='AUDITING')
    result=subprocess.run([sys.executable,'-X','utf8','src/122_audit_continuation_024.py'],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8')
    (c.OUT/'final_audit_execution.txt').write_text(result.stdout+result.stderr,encoding='utf-8')
    if result.returncode:
        snapshot([],phase='AUDIT_FAILED'); raise SystemExit(result.returncode)
    # Save completed new source snapshots and large evidence inputs locally in DVC.
    dvc=Path('C:/Users/administer/.local/bin/dvc.exe')
    names=[str(c.RAW.relative_to(m.ROOT))]+[str((c.OUT/name).relative_to(m.ROOT)) for name in ['evidence_BindingDB.json','evidence_HERB.json','evidence_ChEMBL.json','evidence_PubChem.json']]
    dvcresult=subprocess.run([str(dvc),'add',*names],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8')
    (c.OUT/'dvc_snapshot_execution.txt').write_text(dvcresult.stdout+dvcresult.stderr,encoding='utf-8')
    dvcreplay=subprocess.run([str(dvc),'repro','--single-item','fixed_formula_supplement_024_continuation_replay'],cwd=m.ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace') if dvcresult.returncode==0 else None
    if dvcreplay:
        (c.OUT/'dvc_replay_execution.txt').write_text(dvcreplay.stdout+dvcreplay.stderr,encoding='utf-8')
    cleaned=[]
    for folder,patterns in [(m.ROOT/'src/__pycache__',[f'{prefix}_*.pyc' for prefix in ['118','119','120','121','122','123','124']]),(m.ROOT/'tests/__pycache__',['test_formula_continuation_024*.pyc'])]:
        if not folder.exists(): continue
        for pattern in patterns:
            for path in folder.glob(pattern):
                if path.resolve().parent!=folder.resolve(): raise ValueError('cache cleanup escaped workspace')
                path.unlink(); cleaned.append(str(path.relative_to(m.ROOT)))
    m.writejson(c.OUT/'delivery_status.json',dict(final_audit_returncode=result.returncode,dvc_add_returncode=dvcresult.returncode,dvc_replay_returncode=dvcreplay.returncode if dvcreplay else None,DVC_remote='NOT_CONFIGURED',Git_mixed_changes_preserved=True,completed_at=m.now(),temporary_bytecode_removed=cleaned))
    snapshot([],phase='FINISHED' if dvcresult.returncode==0 and dvcreplay and dvcreplay.returncode==0 else 'DVC_SNAPSHOT_FAILED')


if __name__=='__main__': main()
