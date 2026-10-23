"""Offline local replay; --all also regenerates workbook and expression intermediates."""
import gzip,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';INT=ROOT/'data/interim'/RUN
R='D:/R/R-4.3.3/bin/Rscript.exe';PY='D:/conda/envs/TCM_select/python.exe'
def call(args):subprocess.run(args,cwd=ROOT,check=True)
def main():
    if '--all' in sys.argv:
        work=INT/'workbooks';work.mkdir(parents=True,exist_ok=True)
        old=ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources';raw=ROOT/'data/raw'/RUN
        for source,name in [(raw/'GSE228927_FPKM_anno.xls.gz','GSE228927_FPKM_anno.xls'),(old/'GSE275834_BLEODMY_vs_BLEO.xls.gz','GSE275834_BLEODMY_vs_BLEO.xls'),(raw/'GSE275834_BLEO_vs_CTRL.deseq2.annotation_HMC3_.xls.gz','GSE275834_BLEO_vs_CTRL.deseq2.annotation_HMC3_.xls')]:
            (work/name).write_bytes(gzip.decompress(source.read_bytes()))
        call([R,'src/126_export_workbooks_026.R']);call([PY,'src/127_prepare_expression_026.py']);call([PY,'src/129_audit_TQ_reference_026.py']);call([PY,'src/130_prepare_receiver_extension_026.py']);call([PY,'src/132_prepare_protein_sensitivity_026.py'])
    call([PY,'src/129_audit_TQ_reference_026.py'])
    for args in [[],['receiver_dataset_manifest.json','receiver'],['protein_dataset_manifest.json','protein'],['sample_sensitivity_manifest.json','sample_sensitivity']]:call([R,'src/128_fit_mechanism_effects_026.R']+args)
    call([PY,'src/131_review_samples_sources_026.py']);call([PY,'src/133_report_audit_mechanism_026.py']);call([PY,'src/134_finalize_provenance_026.py']);call([PY,'-m','unittest','discover','-s','tests','-p','test_mechanism_026.py'])
if __name__=='__main__':main()
