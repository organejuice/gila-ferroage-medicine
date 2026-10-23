"""Offline reproducible entry for the final AI evidence review of SOP035."""
from pathlib import Path
import os, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]

def main():
    prefix=Path('D:/conda/envs/zhang2026-nicotine-copd-r')
    env=os.environ.copy();env['CONDA_PREFIX']=str(prefix)
    env['PATH']=os.pathsep.join(str(prefix/p) for p in ['Library/bin','lib/R/bin','lib/R/bin/x64'])+os.pathsep+env.get('PATH','')
    subprocess.run([str(prefix/'lib/R/bin/Rscript.exe'),'src/review_virtual_ko_identity_035.R'],cwd=ROOT,env=env,check=True)
    subprocess.run([sys.executable,'-X','utf8','src/scientific_review_virtual_ko_035.py'],cwd=ROOT,check=True)
    subprocess.run([sys.executable,'-X','utf8','src/finalize_scientific_review_035.py'],cwd=ROOT,check=True)
    subprocess.run([sys.executable,'-X','utf8','src/record_scientific_review_provenance_035.py'],cwd=ROOT,check=True)

if __name__=='__main__':main()
