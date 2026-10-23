"""Archive user-supplied supplements and compare gene-only priors offline.

No paper gene is attributed to a candidate compound. Original rows are retained.
"""
import importlib.util
import json
import shutil
import sys
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

spec=importlib.util.spec_from_file_location('c024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
DOI='10.1177/15491684261460820'
PUB=c.RAW/'Publication_user_supplied'


def text(node):
    return ''.join(t.text or '' for t in node.findall('.//w:t',NS))


def parse_docx(path):
    with ZipFile(path) as z:
        bad=z.testzip()
        if bad: raise ValueError('DOCX CRC failure: '+bad)
        root=ET.fromstring(z.read('word/document.xml'))
        if root.findall('.//w:del',NS) or root.findall('.//w:ins',NS):
            raise ValueError('tracked content requires explicit review')
        body=root.find('w:body',NS)
        paragraphs=[text(n) for n in body.findall('w:p',NS) if text(n).strip()]
        tables=[]
        for table in body.findall('w:tbl',NS):
            rows=[[text(cell) for cell in row.findall('w:tc',NS)] for row in table.findall('w:tr',NS)]
            if not rows or len(set(rows[0]))!=len(rows[0]): raise ValueError('invalid headers')
            if any(len(row)!=len(rows[0]) for row in rows): raise ValueError('table width mismatch')
            tables.append([dict(zip(rows[0],row)) for row in rows[1:]])
        media=[{'member':n,'bytes':len(z.read(n))} for n in z.namelist() if n.startswith('word/media/') and not n.endswith('/')]
    return dict(paragraphs=paragraphs,tables=tables,media=media)


def archive():
    PUB.mkdir(parents=True,exist_ok=True)
    for index in (1,2,3):
        name=f'sj-docx-{index}-rej-10.1177_15491684261460820.docx'
        src=Path('F:/BaiduNetdiskDownload')/name; dest=PUB/name
        mp=PUB/(name+'.json')
        if dest.exists() and mp.exists():
            old=json.loads(mp.read_text(encoding='utf-8'))
            if m.filehash(dest)!=old['sha256']: raise ValueError('archived file changed: '+name)
            continue
        if not dest.exists(): shutil.copyfile(src,dest)
        if m.filehash(src)!=m.filehash(dest): raise ValueError('archived file differs: '+name)
        parsed=parse_docx(dest)
        meta=dict(source='Publication_user_supplied',url=f'https://journals.sagepub.com/doi/suppl/{DOI}/suppl_file/{name}',
            body_path=str(dest.relative_to(m.ROOT)),sha256=m.filehash(dest),bytes=dest.stat().st_size,
            retrieved_at=m.now(),provenance='user supplied local file; not a successful publisher download',
            supplied_path=str(src),zip_crc_verified=True,document_content=parsed)
        if not mp.exists(): m.writejson(mp,meta)


def main():
    archive()
    parsed=[parse_docx(PUB/f'sj-docx-{i}-rej-10.1177_15491684261460820.docx') for i in (1,2,3)]
    if [len(p['tables']) for p in parsed]!=[0,1,1]: raise ValueError('unexpected supplement contents')
    compounds=parsed[1]['tables'][0]; targets=parsed[2]['tables'][0]
    if len(compounds)!=154 or len(targets)!=37: raise ValueError('unexpected row counts')
    if list(compounds[0])!=['Herb','Compound','Formula','MW','#Rotatable bonds','#H-bond acceptors','#H-bond donors','GI absorption','BBB permeant']:
        raise ValueError('compound schema changed')
    if list(targets[0])!=['Gene','Score','Average shortest path length','Closeness centrality','Degree','Betweenness centrality']:
        raise ValueError('target schema changed')
    m.writetable(c.OUT/'publication_supplement_table1.tsv',[dict(source_row=i+2,**r) for i,r in enumerate(compounds)])
    m.writetable(c.OUT/'publication_supplement_table2.tsv',[dict(source_row=i+2,**r) for i,r in enumerate(targets)])
    fs=importlib.util.spec_from_file_location('f024',Path(__file__).with_name('114_finalize_formula_targets_024.py'))
    f=importlib.util.module_from_spec(fs); fs.loader.exec_module(f)
    sets,contexts=f.disease_sets()
    approved={r['symbol'] for r in m.readtable(f.RAW/'HGNC/hgnc_complete_set.tsv') if r['status']=='Approved'}
    genes={r['Gene'].upper() for r in targets}
    if not genes<=approved: raise ValueError('unapproved paper gene symbols: '+str(genes-approved))
    comparisons=[]
    for name,disease in sets.items():
        for gene in sorted(genes):
            comparisons.append(dict(disease_set=name,gene=gene,in_disease_set=gene in disease,tier='H',
                evidence_scope='paper PPI topology gene-only prior; no compound-target pair or source assay',
                disease_context=contexts[name],document=DOI))
    m.writetable(c.OUT/'publication_gene_disease_comparison.tsv',comparisons)
    gates=[]
    for gate in ['G1','G2','G3','G4','G5']:
        matched={r['gene'] for r in comparisons if r['disease_set'].startswith(gate+'_') and r['disease_set']!='G1_significant413_exploratory' and r['in_disease_set']}
        gates.append(dict(gate=gate,genes=len(matched),gene_symbols=';'.join(sorted(matched)),scope='gene-only H comparison'))
    m.writetable(c.OUT/'publication_G1_G5_summary.tsv',gates)
    herbs=Counter(r['Herb'] for r in compounds)
    mapping={
        'Acorus tatarinowii Schott':'石菖蒲','Alpinia oxyphylla Miq.':'益智仁','Carthamus tinctorius L.':'红花',
        'Chinemys reevesii (Gray)':'龟甲','Cinnamomum cassia (L.) J.Presl':'肉桂','Cistanche deserticola Y.C.Ma':'肉苁蓉',
        'Cornus officinalis Siebold et Zucc.':'山茱萸','Dendrobium nobile Lindl.':'石斛','Morinda officinalis F.C.How':'巴戟天',
        'Ophiopogon japonicus (Thunb.) Ker Gawl.':'麦冬','Polygala tenuifolia':'远志','Polygala tenuifolia Willd.':'远志',
        'Prunus persica (L.) Batsch':'桃仁','Rehmannia glutinosa (Gaertn.) DC.':'地黄（炮制未说明）',
        'Salvia miltiorrhiza Bunge':'丹参','Schisandra chinensis (Turcz.) Baill.':'五味子'}
    if set(herbs)!=set(mapping): raise ValueError('herb mapping requires review')
    m.writetable(c.OUT/'publication_herb_comparison.tsv',[dict(paper_herb=h,paper_compound_rows=n,Chinese_correspondence=mapping[h],
        interpretation='botanical correspondence only; preparation/dose equivalence unverified') for h,n in sorted(herbs.items())])
    manifest=[]
    for mp in sorted(PUB.glob('*.docx.json')):
        meta=json.loads(mp.read_text(encoding='utf-8'))
        manifest.append({k:meta[k] for k in ['source','url','body_path','sha256','bytes','provenance']})
    m.writetable(c.OUT/'publication_supplement_manifest.tsv',manifest)
    state=dict(status='USER_SUPPLIED_PARSED',files=3,compound_rows=154,distinct_compound_names=len({r['Compound'] for r in compounds}),
        target_rows=37,approved_genes=len(genes),paper_herb_labels=len(herbs),normalized_botanical_taxa=15,
        formula_equivalence='NOT_ESTABLISHED',paper_extra_herb='益智仁',project_herbs_not_listed=['龙骨','水蛭'],
        processed_rehmannia_and_doses_verified=False,compound_identifiers_present=False,compound_target_pairs_present=False,
        E1_promoted=0,candidate_set_changed=False,C0=True,reopen_017=False,G1_G5=gates)
    m.writejson(c.OUT/'publication_supplement_status.json',state)
    lines=['# 同名益智地黄汤论文补充文件核查','',
        'DOI：'+DOI+'。三个文件由用户提供，ZIP完整性和SHA-256已核验；出版社直链仍受阻，两种来源状态分别保留。','',
        '- 文件1：补充图S1，Figure 6A对应的完整Western blot图像；仅凭图像不建立候选化合物直接结合证据。',
        '- 文件2：补充表1，154条药味—成分记录，148个不同成分名称，包含分子式、分子量、GI absorption和BBB permeant。无CID、SMILES或逐成分—靶点对；High/Yes不能替代本项目ADMET筛选或脑暴露测量。',
        '- 文件3：补充表2，37个PPI核心靶点及拓扑指标。mTOR规范为HGNC批准符号MTOR，原始写法保留。拓扑分数不是活性测量或直接结合证据。','',
        '药味有16种原始标签；两种Polygala tenuifolia写法合并后为15种植物/动物分类名。包含益智仁，未列出本项目龙骨、水蛭；地黄的炮制状态及剂量无法由补表确认。补表成分清单也不能反向证明未列药味未用于正文实验。当前不能宣称同处方复现。','',
        '|集合|与37个核心基因交集|基因|','|---|---:|---|']
    lines.extend(f'|{r["gate"]}|{r["genes"]}|{r["gene_symbols"] or "无"}|' for r in gates)
    lines+=['','上述是独立gene-only H比较，没有接入候选成分—靶点连接表，也不与E1/E2/P计数相加。G1为冻结244基因；其余集合保持原模型/同源映射限制。药物C0、017关闭范围和冻结候选不变。','',
        '复现：TCM_select环境执行src/124_ingest_publication_supplements_024.py；已归档的原始DOCX与独立SHA清单位于本次raw/Publication_user_supplied。后续离线重放不依赖网盘源文件。']
    (c.OUT/'publication_supplement_review.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(state,ensure_ascii=False),flush=True)


def audit():
    outputs=['publication_supplement_table1.tsv','publication_supplement_table2.tsv','publication_gene_disease_comparison.tsv','publication_G1_G5_summary.tsv','publication_supplement_status.json']
    hashes=[]
    for iteration in range(2):
        main()
        hashes.append({name:m.filehash(c.OUT/name) for name in outputs})
    meta=[json.loads(p.read_text(encoding='utf-8')) for p in PUB.glob('*.docx.json')]
    checks=dict(files_verified=len(meta),raw_DOCX_sha256_match=all(m.filehash(m.ROOT/r['body_path'])==r['sha256'] for r in meta),
        paper_offline_replay_match=hashes[0]==hashes[1],paper_output_sha256=hashes[1],scope='publication ingest only; final all-source audit pending')
    m.writejson(c.OUT/'publication_ingest_audit.json',checks)
    if len(meta)!=3 or not checks['raw_DOCX_sha256_match'] or not checks['paper_offline_replay_match']: raise ValueError('publication ingest audit failed')


if __name__=='__main__':
    audit() if '--audit' in sys.argv else main()
