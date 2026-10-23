"""Candidate file/sample audit from official cached snapshots only."""
import gzip,importlib.util,json,re
from collections import Counter
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-disease-interface-perturbation-audit-027';RAW=ROOT/'data/raw'/RUN;OUT=ROOT/'results'/RUN
spec=importlib.util.spec_from_file_location('p127',ROOT/'src/127_prepare_expression_026.py');p127=importlib.util.module_from_spec(spec);spec.loader.exec_module(p127)
def char(r,key):
    for x in r.get('characteristics_ch1','').split(' | '):
        k,sep,v=x.partition(':')
        if sep and k.strip().lower()==key:return v.strip()
    return 'UNKNOWN'
def write(d,name):d.to_csv(OUT/name,sep='\t',index=False)
def main():
    probes=json.loads((RAW/'supplementary_HEAD_snapshot.json').read_text())
    for item in probes:
        listing=pd.read_csv(RAW/(item['accession']+'_filelist.txt'),sep='\t');match=listing[listing.Name==item['url'].rsplit('/',1)[-1]]
        declared=int(match.Size.iloc[0]) if len(match)==1 else None;length=int(item['bytes']) if str(item.get('bytes','')).isdigit() else None
        item['filelist_declared_bytes']=declared;item['size_check']='MATCH' if declared is not None and length==declared else 'MISMATCH' if declared is not None and length is not None else 'UNVERIFIED_LENGTH_OR_NO_FILELIST_ENTRY'
    write(pd.DataFrame(probes),'candidate_file_access.tsv')
    samples=[];groups=[];decisions=[];files=[]
    for acc in ['GSE286094','GSE315101','GSE333352']:
        records=p127.soft(RAW/(acc+'_family.soft.gz'));listing=pd.read_csv(RAW/(acc+'_filelist.txt'),sep='\t');listing['accession']=acc;files.append(listing)
        rows=[]
        for r in records:
            if acc=='GSE286094':cohort='ATG7_5xFAD_factorial';group=char(r,'genotype');unit=r['title'];unit_status='DECLARED_BIO_REP_ANIMAL_ID_NOT_EXPLICIT'
            elif acc=='GSE315101':
                group=char(r,'genotype');cohort='human_oligodendrocyte_sDLK1' if 'oligodendrocyte' in char(r,'tissue') else 'AAV_DLK1_8mo' if 'LG180' in r['title'] else 'aged_cGAS_control_only' if 'LG107' in r['title'] else 'Terc_8mo';unit=r['title'];unit_status='NAMED_LIBRARY_ANIMAL_OR_DONOR_ID_NOT_EXPLICIT'
            else:cohort='aged_GeoMx_28mo';group=char(r,'label');unit=char(r,'mouse id');unit_status='EXPLICIT_MOUSE_ID_AOIS_NESTED'
            urls=re.findall(r'https?://[^ |]+|ftp://[^ |]+',' | '.join(v for k,v in r.items() if k.startswith('supplementary_file')))
            names=[u.rsplit('/',1)[-1] for u in urls]
            required=['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'] if acc!='GSE333352' else ['.dcc.gz']
            complete=all(any(n.endswith(s) for n in names) for s in required)
            listed=all(n in set(listing.Name) for n in names)
            relevant=[x for x in probes if x['accession']==acc and x['url'].rsplit('/',1)[-1] in names]
            row=dict(accession=acc,GSM=r['GSM'],title=r['title'],cohort=cohort,group=group,organism=r.get('organism_ch1','UNKNOWN'),sex=char(r,'sex'),age=char(r,'age'),region=char(r,'tissue region') if acc=='GSE333352' else char(r,'tissue'),independent_unit=unit,unit_status=unit_status,ROI=r['title'].split(',')[0] if acc=='GSE333352' else 'NOT_APPLICABLE',supplementary_files=';'.join(names),expected_file_roles_declared=complete,filelist_match=listed,file_urls_accessible=len(relevant)==len(names) and all(x['accessible'] for x in relevant),matrix_dimensions_verified=False,barcodes_cell_annotation_verified=False,numeric_spatial_coordinates_available=False)
            rows.append(row);samples.append(row)
        frame=pd.DataFrame(rows)
        for (cohort,group),d in frame.groupby(['cohort','group']):groups.append(dict(accession=acc,cohort=cohort,group=group,libraries_or_AOIs=len(d),distinct_unit_labels=d.independent_unit.nunique(),independence_verified=acc=='GSE333352',sex_counts=json.dumps(Counter(d.sex),sort_keys=True),age_counts=json.dumps(Counter(d.age),sort_keys=True)))
        access=[x for x in probes if x['accession']==acc]
        archive=listing[listing['#Archive/File']=='Archive']
        decisions.append(dict(accession=acc,samples_or_AOIs=len(frame),explicit_animals=frame.independent_unit.nunique() if acc=='GSE333352' else 'NOT_CONFIRMED',cohorts=frame.cohort.nunique(),declared_file_roles_complete=bool(frame.expected_file_roles_declared.all()),filelist_all_matched=bool(frame.filelist_match.all()),file_urls_checked=len(access),file_urls_accessible=sum(x['accessible'] for x in access),archive_bytes=int(archive.Size.iloc[0]),matrix_structure_verified=False,cell_type_annotation_supplied='NO' if acc!='GSE333352' else 'SEGMENT_LABELS_ONLY_NOT_SINGLE_CELL',numeric_spatial_coordinates_verified=False,decision='PREPARE_NEW_MODEL_SCOPED_ANALYSIS_WITH_ID_CONFIRMATION' if acc=='GSE286094' else 'SPLIT_COHORTS_AND_CONFIRM_UNITS_BEFORE_ANALYSIS' if acc=='GSE315101' else 'AGED_WITHIN_ANIMAL_SPATIAL_DESCRIPTION_ONLY',reason='Four genotypes declared biological libraries 3/3/2/3; animal ID/pooling/sex absent; confirm before pseudobulk' if acc=='GSE286094' else 'AAV 4Ctrl/4DLK1, human oligodendrocyte 3/3, Terc 5/5, old cGAS 3control-only; never merge cohorts; series says male but sample annotations include females' if acc=='GSE315101' else '143 AOIs from 5 mice, all aged; no young or AD/drug control; 5 animals are the independent units; segmented GeoMx, no verified coordinates or full M/A chain',next_required='Validate matrix dimensions, barcode namespace, gene assembly; confirm animal/donor pooling; annotate/QC without automatic thresholds; freeze new sample-level model',human_AD_therapeutic_eligible=False,reopen_G6=False))
    write(pd.DataFrame(samples),'candidate_sample_registry.tsv');write(pd.DataFrame(groups),'candidate_group_independence_audit.tsv');write(pd.DataFrame(decisions),'candidate_feasibility_decisions.tsv');write(pd.concat(files,ignore_index=True),'candidate_archive_filelist.tsv')
    # The study-level statement is retained separately from per-sample unknown fields.
    series=[]
    for acc in ['GSE286094','GSE315101','GSE333352']:
        text=gzip.open(RAW/(acc+'_family.soft.gz'),'rt',encoding='utf8').read()
        for line in text.splitlines():
            if line.startswith(('!Series_title =','!Series_overall_design =','!Series_last_update_date =')):
                key,value=line.split(' = ',1);series.append(dict(accession=acc,field=key,value=value))
    write(pd.DataFrame(series),'candidate_series_scope.tsv')
    print('candidate samples/AOIs',len(samples),'file URLs',len(probes),'decisions',len(decisions),flush=True)
if __name__=='__main__':main()
