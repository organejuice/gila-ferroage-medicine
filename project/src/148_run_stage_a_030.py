"""Offline stage-A pipeline; resources must already be frozen by 144."""
import subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    subprocess.run([sys.executable,'-X','utf8','src/144_prepare_stage_a_030.py','--inputs-only'],cwd=ROOT,check=True)
    subprocess.run(['D:/R/R-4.3.3/bin/Rscript.exe','src/145_enrich_stage_a_030.R'],cwd=ROOT,check=True)
    # Force network recomputation only when requested; default uses stored deterministic null results.
    if '--recompute-network' in sys.argv:
        run=json.loads((ROOT/'configs/stage_a_030.json').read_text())['run_id'];out=ROOT/'results'/run
        for p in out.glob('network_proximity_*.tsv'):p.unlink()
    subprocess.run([sys.executable,'-X','utf8','src/146_network_stage_a_030.py'],cwd=ROOT,check=True)
    subprocess.run([sys.executable,'-X','utf8','src/147_finalize_stage_a_030.py'],cwd=ROOT,check=True)
if __name__=='__main__':main()
