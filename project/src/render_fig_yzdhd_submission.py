"""Canonical document rendering with the user-specified native Windows office."""
from pathlib import Path
import os,sys,runpy,subprocess,shutil,json
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/2026-10-14-FIG-YZDHD'
QA=ROOT/'tmp/FIG_YZDHD_20261014'
RENDERER=Path('C:/Users/administer/.codex/plugins/cache/openai-primary-runtime/documents/26.1013.11914/skills/documents/render_docx.py')
os.environ['PATH']='D:/Libreoffice/program;'+os.environ['PATH']
native_run=subprocess.run
def tolerant_run(*args,**kwargs):
    if kwargs.get('text') or kwargs.get('universal_newlines'): kwargs.setdefault('errors','replace')
    return native_run(*args,**kwargs)
subprocess.run=tolerant_run
files=[OUT/'manuscript/main_manuscript.docx',OUT/'administrative/title_page.docx',OUT/'administrative/cover_letter.docx',QA/'ESM_1.docx']
selected=sys.argv[1:]
if selected: files=[f for f in files if f.stem in selected]
counts={}
for f in files:
    destination=QA/('render_'+f.stem)
    sys.argv=['render_docx.py',str(f),'--output_dir',str(destination),'--emit_pdf','--dpi','120']
    runpy.run_path(str(RENDERER),run_name='__main__')
    n=len(PdfReader(destination/(f.stem+'.pdf')).pages); counts[f.stem]=n
    for p in destination.glob('page-*.png'):
        if int(p.stem.split('-')[1])>n: p.unlink()
    if f.stem=='ESM_1': shutil.copy2(destination/'ESM_1.pdf',OUT/'supplementary/ESM_1.pdf')
print(json.dumps(counts))
