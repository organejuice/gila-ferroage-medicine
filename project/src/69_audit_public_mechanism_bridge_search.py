"""Compile discovery metadata and audit downloaded public inputs; no effect fitting."""
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = '2026-10-02-public-mechanism-bridge-search-025'
RAW = ROOT / 'workspaces' / RUN / 'raw' / 'sources'
OUT = ROOT / 'results' / RUN
spec = importlib.util.spec_from_file_location('existing_geo_summary', ROOT / 'src/42_summarize_g5_ncbi_candidates.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Curated decisions are discovery judgements, not reproduced biological results.
CURATED = [
('GSE325931','A','10.1002/2211-5463.70182','human iPSC tri-culture; microglial pseudobulk','iron/RSL3 and culture-context perturbation','CPM_ONLY_DOWNLOADED','17 GEO libraries vs 16 CPM columns: missing a31 co-culture rep4; tri vehicle 6 and iron/RSL3 4 are present; batch/well independence requires audit; no astrocyte matrix or raw cell counts located'),
('GSE290755','A','10.1038/s41467-026-76131-0','human iPSC microglia','ARID5A CRISPRi; iron/inflammation/ferroptosis','RNA_RPKM_DOWNLOADED','9 GEO libraries: RNA control/KD 3 each; RPF control 2/KD 1; downloaded RNA RPKM has KD 3/control 2, not 3/3; file named mrna_reps_condition1and2 actually contains translation efficiency statistics; phenotype co-culture belongs to publisher source data'),
('GSE228927','A','','mouse primary astrocytes','resting/LPS/hemin microglial conditioned medium','OFFICIAL_FILES_LISTED','3 groups x 3 samples; FPKM; stroke context; residual LPS/hemin effects require paper controls'),
('GSE175485','A','','mouse astrocytes','activated microglial vs endothelial conditioned medium','OFFICIAL_FILES_LISTED','processed XLSX files; sample-level design needs further audit'),
('GSE315100','A','10.64898/2026.01.14.699608','human iPSC microglia','telomere/senescence induction; DLK1 secretory route','COUNTS_DOWNLOADED','28 GEO samples; matrix only bal_001-016: Ctrl/BIBR1532/Ctrl+etoposide/BIBR1532+etoposide 4 each; same WTC11 male line; full SOFT mapping saved; DLK1 recipient phenotype not included in this bulk matrix; preprint'),
('GSE275834','A','','human HMC3 immortalized microglia','BLEO senescence and BLEO+DMY intervention','OFFICIAL_FILES_LISTED','CTRL/BLEO/BLEODMY 3 each; processed DE tables and raw reads; DMY is an external reference, no confirmed current-formula identity'),
('GSE345471','A','','mouse microglia','aged vehicle vs aged neutral-pH CatB inhibitor vs young','OFFICIAL_FILES_LISTED','3 groups x 3 samples; GSM confirms microglia; processed FPKM; aging not AD'),
('GSE275256','B','','human HMC3 immortalized microglia','chemical senescence dose series','ESUMMARY_VERIFIED_SERIES_FETCH_FAILED','vehicle/SLO half/SLO full 3 each; induction benchmark only'),
('GSE240609','B','10.1016/j.celrep.2024.114982','human isogenic iPSC neuron-microglia co-culture','APOE3/Christchurch x neuronal WT/PSEN1','OFFICIAL_FILES_LISTED','4 deposited libraries across 4 conditions; cannot claim replicated RNA differential effects from this GEO alone'),
('GSE286095','B','','aged mouse sorted microglia','Atg7 conditional genotype perturbation','OFFICIAL_FILES_LISTED','gene count matrix; aging not automatically the amyloid intervention cohort'),
('GSE286094','B','','mouse single-cell','ATG7/ferroptosis/amyloid route','OFFICIAL_FILES_LISTED','scRNA-seq; per-sample genotype/amyloid and biological unit audit pending'),
('GSE333352','B','10.1038/s43587-026-01154-7','mouse GeoMx spatial profiling','spatial senescent/DAM reference','OFFICIAL_FILES_LISTED','143 GEO entries are spatial segments, not 143 mice; senolytic-treated transcriptome contrast not verified; no human AD-control claim'),
('GSE315101','B','10.64898/2026.01.14.699608','mouse snRNA-seq','DLK1-related recipient/cellular effects','OFFICIAL_FILES_LISTED','specific deposited perturbation arms require GSM audit; preprint'),
('GSE345470','B','','mouse microglial cell models','CatB/lysosome/senescence route','ESUMMARY_VERIFIED_SERIES_FETCH_FAILED','MG6/MG7 control/palmitate each 3; not inhibitor rescue transcriptome'),
('GSE341515','C','','mouse whole midbrain/colon','p53 genotype x PFF synucleinopathy','OFFICIAL_FILES_LISTED','24 libraries across organs; D+Q experiment in study does not establish D+Q RNA deposition; not AD or purified microglia'),
('GSE277435','C','','human MS Visium','spatial senescent-like glia reference','OFFICIAL_FILES_LISTED','12 slices are not necessarily 12 donors; multiple sclerosis rather than AD'),
('GSE280341','C','','organoid scRNA-seq','inflammation-induced senescent-like glia','ESUMMARY_VERIFIED_SERIES_FETCH_FAILED','donor and induction arms pending audit; MS rather than AD'),
('GSE235076','B','','mouse neurons','microglial/astrocytic apoE particles recipient response','OFFICIAL_FILES_LISTED','cross-cell secreted-material perturbation; no histone lactylation direct measurement'),
('GSE246888','B','','mouse neurons','apoE particles recipient response','OFFICIAL_FILES_LISTED','same article family as GSE235076; avoid independent-study double counting'),
('GSE18740','B','','mouse BV2','luteolin x LPS factorial intervention','OFFICIAL_FILES_LISTED','4 groups x 3 biological replicates; actual compound perturbation; current recipe chemical identity/exposure still unverified'),
('GSE161340','C','','mouse hippocampal sc/snRNA-seq','age/senescence reference','OFFICIAL_FILES_LISTED','GEO design is 4m vs 24m age comparison; do not relabel as senolytic-treated transcriptome'),
('GSE229553','C','','mouse WT/P301S microglia','tauopathy senescence reference','OFFICIAL_FILES_LISTED','WT/P301S observed state, not deposited senolytic contrast'),
('GSE188765','EXISTING','','mouse 5XFAD microglia/brain','H4K12la chromatin reference','OFFICIAL_FILES_LISTED','existing G2 reused; no new effect reproduced'),
('GSE227221','EXISTING','10.1038/s41593-022-01221-3','human iPSC microglia','iron/RSL3 time course','LOCAL_EXISTING_MATRIX','existing D2/G3/G4; broader tri-culture Synapse deposits have distinct access routes'),
('GSE296797','EXISTING','','mouse EAE microglia/macrophages','H4K12la/H3K18la and glycolysis','ESUMMARY_VERIFIED','existing alternative model; not independent AD validation'),
('GSE272754','EXISTING','','mouse BV2 OGD','A485/OGD transcriptome and H3K9la paper','ESUMMARY_VERIFIED','RNA-seq is not site-specific chromatin measurement; H3K9la is not H4K12la'),
('GSE283233','EXISTING','','mouse BV2 ischemia model','SMEK1/H3K9la','ESUMMARY_VERIFIED','RNA intervention/model reference; H3K9la is not H4K12la'),
('E-MTAB-15005','A','10.1038/s41467-025-66012-3','five human brain cell lines','BrdU induced senescence; microglia and astrocyte comparison','EBI_API_AND_MATRIX_VERIFIED','repository design BrdU vs DMSO; TFAP4 knockdown RNA deposition not established; independent cell cultures not communication'),
('PXD079947','B','10.1155/omcl/9948225','human HMC3','iron/RSL3 ferroptosis proteomics','PUBLIC_PROTEOMEXCHANGE_RECORD','proteomics public Sep13 2026; broader transcriptomic/lipidomic data stated available upon request, not automatically public'),
('FIGSHARE_30252505','B','10.1038/s42003-026-09772-y','5XFAD Srr conditional model','competing lactylation direction','PUBLISHER_DECLARED_FIGURE_SOURCE_DATA','Srr deletion associated with increased H3K18la and improved male phenotype; Lyz2-Cre specificity limits; no raw omics access confirmed'),
('YZDHD_2026','C','10.1177/15491684261460820','MCI-like mouse model','exact-name decoction publication','ABSTRACT_AND_SUPPLEMENT_LIST_VERIFIED','publisher restricted full text; current 16-herb identity not verified; GEO exact-name query returned zero; no matched wholeformula omics located'),
('TRAP1_VSMC','C','10.1093/eurheartj/ehae379','vascular smooth muscle','H4K12la to SASP intervention reference','REQUEST_ONLY_OR_SUPPLEMENT','paper reuses AD chromatin for SASP promoters; own direct validation in VSMC; cannot substitute for microglial same-system proof'),
]

def series_fields(path):
    values = {}
    if path.exists():
        for line in path.read_text(encoding='utf-8-sig').splitlines():
            if line.startswith('!Series_') and ' = ' in line:
                key, value = line.split(' = ', 1)
                values.setdefault(key.replace('!Series_', ''), []).append(value)
    return values

def main():
    records = {}
    for path in sorted(RAW.glob('geo_esummary_*.json')):
        result = json.loads(path.read_text(encoding='utf-8-sig'))['result']
        for uid in result.get('uids', []):
            rec = result[str(uid)]
            records[rec['accession']] = rec
    candidates = []
    for acc, tier, doi, context, bridge, status, limitation in CURATED:
        rec = records.get(acc, {})
        sf = series_fields(RAW / (acc + '.txt'))
        samples = rec.get('samples', [])
        row = dict(accession=acc,priority=tier,doi=doi,
                   title=' | '.join(sf.get('title', [])) or rec.get('title', '') or acc,
                   context=context,bridge=bridge,access_status=status,
                   n_deposited_entries=len(sf.get('sample_id', [])) or len(samples) or '',
                   public_status=' | '.join(sf.get('status', [])),
                   design=' | '.join(sf.get('overall_design', [])),
                   pubmed=' | '.join(sf.get('pubmed_id', [])) or module.as_text(rec.get('pubmedids')),
                   files=' | '.join(sf.get('supplementary_file', [])),
                   sample_titles=module.as_text(samples),limitations=limitation,
                   url=('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc if acc.startswith('GSE') else 'https://doi.org/'+doi))
        candidates.append(row)
    with (OUT/'candidate_registry.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(candidates[0]),delimiter='\t');w.writeheader();w.writerows(candidates)
    # Matrix existence and readable header are access checks, not a reanalysis gate.
    matrices=[]
    for path in sorted(RAW.glob('*.gz')):
        if path.name.endswith('family.soft.gz'):
            continue
        with gzip.open(path,'rb') as f:
            signature=f.read(8)
        if signature.startswith(bytes.fromhex('d0cf11e0a1b11ae1')):
            matrices.append(dict(file=path.name,bytes=path.stat().st_size,type='binary XLS compound document; downloaded, workbook contents not parsed'))
            continue
        with gzip.open(path,'rt',encoding='utf-8-sig') as f:
            header=f.readline().strip(); first=f.readline().strip()
        matrices.append(dict(file=path.name,bytes=path.stat().st_size,header=header[:4000],first_row_columns=len(first.split('\t')),header_columns=len(header.split('\t'))))
    for path in sorted(RAW.glob('*TPM_values.csv')):
        with path.open(encoding='utf-8-sig') as f: header=f.readline().strip()
        matrices.append(dict(file=path.name,bytes=path.stat().st_size,header=header[:4000]))
    sample_rows=[]
    for path in sorted(RAW.glob('*family.soft.gz')):
        current=None
        with gzip.open(path,'rt',encoding='utf-8-sig') as f:
            for line in f:
                line=line.rstrip('\n')
                if line.startswith('^SAMPLE = '):
                    current={'series':path.name.split('_')[0],'sample':line.split(' = ',1)[1]}
                    sample_rows.append(current)
                elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                    k,v=line.split(' = ',1);k=k.replace('!Sample_','')
                    if k in ['title','source_name_ch1','characteristics_ch1','description']:
                        current[k]=current.get(k,'')+' | '+v
    with (OUT/'sample_mapping.tsv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['series','sample','title','source_name_ch1','characteristics_ch1','description'],delimiter='\t');w.writeheader();w.writerows(sample_rows)
    wanted=json.loads((OUT/'geo_unique_ids.json').read_text(encoding='utf-8-sig'))
    audit=dict(run=RUN,scope='discovery_and_access_only',unique_geo_search_ids=len(wanted),official_summary_records=len(records),
               curated_rows=len(candidates),priorities={t:sum(x['priority']==t for x in candidates) for t in ['A','B','C','EXISTING']},
               no_duplicate_curated_accessions=len(candidates)==len({x['accession'] for x in candidates}),
               required_decision_fields_complete=all(x['bridge'] and x['limitations'] and x['access_status'] for x in candidates),
               downloaded_matrix_checks=matrices,biological_effects_recomputed=False,new_wet_lab=False,
               unresolved_input_discrepancies=['GSE325931: missing co-culture a31/rep4; desired tri-culture vehicle 6 vs iron/RSL3 4 mapped and present, independent-batch design still requires audit','GSE290755: 6 RNA GEO libraries vs 5 RPKM columns; control replicate 3 missing in downloaded RPKM','GSE315100: 28 GEO samples vs 16 count columns; analyzed matrix scope is bal_001-016 only'],
               limitation='Broad GEO queries capped at 100 each; metadata screen is not exhaustive full-text or sample-level review; uncertain donors/arms retained as unresolved.')
    (OUT/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    manifest=[]
    for path in sorted(RAW.rglob('*')):
        if path.is_file():manifest.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False,indent=2))

if __name__ == '__main__': main()
