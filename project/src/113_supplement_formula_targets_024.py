"""Frozen public-target supplement 024. Run prepare, fetch, then analyze.

Raw responses are immutable successful request snapshots. A failed request is
never a negative assay. No credentials/proxy overrides, name-derived gene IDs,
similarity-to-exact conversion, or cross-source evidence inflation are used.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import re
import sys
import time
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import urlencode

import requests
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'configs/fixed_formula_target_supplement_024.json'
VALIDATED_RULES = json.loads((ROOT/'configs/fixed_formula_validated_rules_024.json').read_text(encoding='utf-8'))
COMMENT_RULES = json.loads((ROOT/'configs/fixed_formula_activity_comment_rules_024.json').read_text(encoding='utf-8'))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def filehash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def readtable(path):
    path = Path(path)
    op = gzip.open if path.suffix == '.gz' else open
    with op(path, 'rt', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t' if '.tsv' in path.name else ','))


def writejson(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def writetable(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or sorted({k for r in rows for k in r})
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def canonical(smiles):
    m = Chem.MolFromSmiles(str(smiles or ''))
    return Chem.MolToSmiles(m, isomericSmiles=True) if m else ''


def inchikey(smiles):
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m else ''


def number(value):
    try:
        n = float(str(value).strip())
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def chembl_outcome(comment):
    label=str(comment or '').strip().lower()
    if label in COMMENT_RULES['inactive_exact_labels']: return 'Inactive'
    if label in COMMENT_RULES['inconclusive_exact_labels']: return 'Inconclusive'
    if label in COMMENT_RULES['unspecified_exact_labels'] or any(p in label for p in COMMENT_RULES['no_measured_curve_phrases']): return 'Unspecified'
    if label=='active': return 'Active'
    return 'quantitative_record'


def classify(r):
    if r.get('similarity_only'):
        return 'S', 'actual assay compound not the exact candidate structure'
    if r.get('association_only'):
        return 'H', 'association is neither verified experimental modulation nor traceable prediction'
    if r.get('predicted'):
        return ('P', 'traceable database prediction') if r.get('identity_exact') and r.get('prediction_method') else ('H', 'prediction identity/method unresolved')
    reasons = []
    if r.get('candidate_identity_verified') is False: reasons.append('catalogue CID/structure not confirmed')
    if not r.get('identity_exact'): reasons.append('chemical identity unresolved/mismatched')
    if not r.get('target_stable_id') or not r.get('target_symbol'): reasons.append('stable target ID not resolved')
    if str(r.get('species')) != '9606': reasons.append('human target not established')
    if 'unknown origin' in (r.get('assay_description') or '').lower(): reasons.append('original assay species/origin ambiguous')
    if r.get('single_protein') is not True: reasons.append('single protein not established')
    if r.get('assay_mode') != 'binding': reasons.append('direct binding/biochemical assay not established')
    if r.get('source') == 'ChEMBL' and str(r.get('confidence')) != '9': reasons.append('assay-target confidence not 9')
    if str(r.get('activity_outcome', '')).lower() in ('inactive', 'inconclusive', 'unspecified'): reasons.append('not a positive activity result')
    if number(r.get('value')) is None or number(r.get('value')) <= 0: reasons.append('valid positive quantitative value missing')
    if r.get('units') not in VALIDATED_RULES['molar_units']: reasons.append('interpretable molar unit missing')
    if r.get('activity_type') not in VALIDATED_RULES['concentration_endpoints']: reasons.append('quantitative concentration endpoint not established')
    if r.get('relation') not in ('=', '<', '<=', '>', '>='): reasons.append('measurement relation unresolved')
    if not r.get('assay_id') or not r.get('document'): reasons.append('original assay/document missing')
    if r.get('validity_comment'): reasons.append('source validity warning')
    if r.get('target_type') in ('CELL-LINE', 'ORGANISM', 'TISSUE'): reasons.append('phenotypic target')
    return ('E2', '; '.join(reasons)) if reasons else ('E1', 'quantitative exact-identity human single-protein binding record')


def evidence_key(r):
    # Source-independent key is used only when an original measurement identity
    # exists; same DOI + value alone cannot establish identical assay conditions.
    original = r.get('original_measurement_id', '')
    if original:
        return '|'.join([r.get('smiles', ''), r.get('target_hgnc_id') or r.get('target_stable_id', ''), original])
    # Without a proven common measurement ID preserve distinct annotations and
    # outcomes; AID/SID alone can describe several target chains/accessions.
    return '|'.join([r.get('source', ''), r.get('activity_id', ''), r.get('assay_id', ''), r.get('compound_key', ''), r.get('target_stable_id', ''), str(r.get('value', '')), r.get('units', ''), r.get('relation', ''),r.get('target_accession',''),r.get('activity_type',''),r.get('activity_outcome','')])


def split_batches(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


class SourceError(RuntimeError):
    pass


class Client:
    def __init__(self, source, cfg, raw):
        self.source, self.cfg, self.raw = source, cfg, Path(raw) / source
        self.raw.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.last = 0.0
        self.log = []

    def get(self, url, params=None, allow404=False):
        full = url + ('?' + urlencode(params) if params else '')
        key = digest(full.encode())
        payload = self.raw / (key + '.body')
        meta = self.raw / (key + '.json')
        if meta.exists():
            old = json.loads(meta.read_text(encoding='utf-8'))
            if old['status'] in (200, 404) and payload.exists() and filehash(payload) == old['sha256']:
                if old['status'] == 404 and not allow404: raise SourceError('cached HTTP404: ' + full)
                return payload.read_bytes(), old
        failures = []
        delays = [0] + self.cfg['http']['backoff_seconds']
        for attempt, delay in enumerate(delays):
            time.sleep(delay)
            time.sleep(max(0, self.cfg['http']['interval_seconds'] - (time.monotonic() - self.last)))
            self.last = time.monotonic()
            try:
                r = self.session.get(url, params=params, timeout=self.cfg['http']['timeout_seconds'])
                retry_after=r.headers.get('Retry-After','')
                if r.status_code in (408,429,500,502,503,504) and retry_after:
                    try: pause=float(retry_after)
                    except ValueError:
                        try: pause=(parsedate_to_datetime(retry_after)-datetime.now(timezone.utc)).total_seconds()
                        except (ValueError,TypeError): pause=0
                    time.sleep(max(0,pause))
                if r.status_code == 200 or (r.status_code == 404 and allow404):
                    data = r.content
                    payload.write_bytes(data)
                    record = {'source': self.source, 'url': full, 'retrieved_at': now(), 'status': r.status_code, 'sha256': digest(data), 'bytes': len(data), 'body_path': str(payload.relative_to(ROOT)), 'attempts': attempt + 1}
                    writejson(meta, record)
                    self.log.append(record)
                    return data, record
                if r.status_code not in (408, 429, 500, 502, 503, 504):
                    raise SourceError(f'HTTP{r.status_code}: {full}')
                failures.append(f'HTTP{r.status_code}')
            except requests.RequestException as exc:
                failures.append(type(exc).__name__ + ': ' + str(exc)[:260])
        failure = {'source': self.source, 'url': full, 'retrieved_at': now(), 'status': 'RETRY_EXHAUSTED', 'failures': failures, 'attempts': len(delays)}
        writejson(self.raw / (key + '.failure.json'), failure)
        raise SourceError('; '.join(failures))

    def json(self, url, params=None, allow404=False):
        data, meta = self.get(url, params, allow404)
        if meta['status'] == 404: return None, meta
        try:
            return json.loads(data), meta
        except ValueError as exc:
            raise SourceError('JSON parse failure: ' + meta['url']) from exc

    def pages(self, url, field, params=None):
        rows, hashes = [], []
        visited = set()
        while url:
            if url in visited: raise SourceError('pagination cycle')
            visited.add(url)
            data, meta = self.json(url, params)
            if not isinstance(data, dict) or field not in data: raise SourceError('missing paged field ' + field)
            rows.extend(data[field]); hashes.append(meta['sha256'])
            pg = data.get('page_meta', {})
            nxt = pg.get('next')
            if nxt and nxt.startswith('/'): nxt = 'https://www.ebi.ac.uk' + nxt
            url, params = nxt, None
            if not url and pg.get('total_count') is not None and len(rows) != int(pg['total_count']):
                raise SourceError('incomplete pagination')
        return rows, ';'.join(hashes)


def prepare(cfg, out):
    manifest = out / 'source_sha256.tsv'
    inputs = {k: ROOT / v for k, v in cfg['inputs'].items()}
    inputs.update({'configuration': CONFIG, 'code': Path(__file__)})
    hashes = [{'role': k, 'path': str(v.relative_to(ROOT)), 'sha256': filehash(v), 'bytes': v.stat().st_size} for k, v in inputs.items()]
    if manifest.exists():
        old = {r['role']: r['sha256'] for r in readtable(manifest)}
        if old != {r['role']: r['sha256'] for r in hashes}: raise RuntimeError('Frozen inputs/code changed; use new run or explicit amendment')
    else:
        writetable(manifest, hashes)
    candidates = readtable(inputs['candidates'])
    src = readtable(inputs['source_records'])
    models = {r['Canonical_Isomeric_SMILES']: r['Step7_Status'] for r in readtable(inputs['admet'])}
    maps = defaultdict(lambda: {'cids': set(), 'herbs': set(), 'names': set(), 'candidate_ids': set()})
    for r in src:
        s = r['Structure_Canonical_Isomeric_SMILES']
        if not s: continue
        maps[s]['cids'].add(r['PubChem_CID']) if r['PubChem_CID'] else None
        maps[s]['herbs'].add(r['Herb']); maps[s]['names'].add(r['Ingredient_Name_Raw'])
    rows = []
    for r in readtable(inputs['structures']):
        s = r['Canonical_Isomeric_SMILES']; m = maps[s]
        rows.append({'compound_key': 'STR:' + digest(s.encode())[:24], 'smiles': s, 'inchikey': inchikey(s), 'cids': ';'.join(sorted(m['cids'])), 'herbs': ';'.join(sorted(m['herbs'])), 'names': ';'.join(sorted(m['names'])), 'phase': '1' if models.get(s) == 'pass_predictive_heuristic' else '2', 'identity_status': 'local_structure_pending_source_crosscheck'})
    rows.sort(key=lambda r: (r['phase'], r['compound_key']))
    cids = {r['PubChem_CID'] for r in candidates if r['PubChem_CID']}
    checks = {'candidate_cids': len(cids), 'structures': len(rows), 'priority_structures': sum(r['phase'] == '1' for r in rows)}
    for k,v in checks.items():
        if v != cfg['expected'][k]: raise RuntimeError(f'input count {k}: {v}')
    writetable(out / 'compound_queries.tsv', rows)
    writetable(out / 'candidate_identity_manifest.tsv', candidates)
    writejson(out / 'environment.json', {'python': sys.version, 'rdkit': rdBase.rdkitVersion, 'requests': requests.__version__, 'platform': platform.platform(), 'executable': sys.executable, 'dependencies_installed': False})
    print(json.dumps({'prepare': checks}, ensure_ascii=False), flush=True)


def load_hgnc(path):
    rows = readtable(path)
    by_entrez, by_uniprot = {}, {}
    for r in rows:
        if r.get('status') != 'Approved': continue
        rec = {'symbol': r['symbol'], 'hgnc_id': r['hgnc_id'], 'entrez': r.get('entrez_id', '')}
        if r.get('entrez_id'): by_entrez[r['entrez_id']] = rec
        for a in r.get('uniprot_ids', '').split('|'):
            if a: by_uniprot[a] = rec
    return by_entrez, by_uniprot


def base_evidence(q, source):
    return {'compound_key': q['compound_key'], 'smiles': q['smiles'], 'cids': q.get('confirmed_cids') or q['cids'], 'catalogue_cids':q['cids'], 'herbs': q['herbs'], 'phase': q['phase'], 'source': source, 'candidate_identity_verified':bool(q.get('confirmed_cids')), 'action_direction': 'undetermined', 'manual_source_review': 'pending', 'formula_presence_measured': False}


def ledger(q, source, status, detail='', response_hash='', pages_complete=False):
    return {'compound_key': q['compound_key'], 'phase': q['phase'], 'source': source, 'query_key': q['inchikey'] if source == 'ChEMBL' else q['smiles'] if source == 'BindingDB' else q['cids'], 'status': status, 'detail': detail, 'response_sha256': response_hash, 'pages_complete': pages_complete, 'updated_at': now()}


def checkpoint(out, source, rows, evidence):
    writetable(out / f'ledger_{source}.tsv', rows, ['compound_key','phase','source','query_key','status','detail','response_sha256','pages_complete','updated_at'])
    writejson(out / f'evidence_{source}.json', evidence)


def chembl_fetch(cfg, out, raw, queries, uniprot):
    source = 'ChEMBL'; client = Client(source, cfg, raw)
    done, evidence = [], []
    root = cfg['sources'][source]
    try:
        status, meta = client.json(root + 'status.json')
        writejson(out / 'ChEMBL_version.json', status)
        for phase in ('1','2','3'):
            items = [q for q in queries if q['phase'] == phase]
            for batch in split_batches(items, cfg['batch_sizes']['chembl_molecules']):
                keys = ','.join(q['inchikey'] for q in batch if q['inchikey'])
                molecules, mh = client.pages(root + 'molecule.json', 'molecules', {'molecule_structures__standard_inchi_key__in': keys, 'limit': 1000})
                mol_by_key = defaultdict(list)
                for m in molecules: mol_by_key[m.get('molecule_structures', {}).get('standard_inchi_key', '')].append(m)
                q_by_mol = {}; mol_smiles = {}
                for q in batch:
                    for m in mol_by_key[q['inchikey']]:
                        mol_byid = m['molecule_chembl_id']; q_by_mol[mol_byid] = q
                        mol_smiles[mol_byid] = canonical(m.get('molecule_structures', {}).get('canonical_smiles', ''))
                activities, ah = [], ''
                for ids in split_batches(list(q_by_mol), cfg['batch_sizes']['chembl_activities']):
                    part, ph = client.pages(root + 'activity.json', 'activities', {'molecule_chembl_id__in': ','.join(ids), 'target_organism': 'Homo sapiens', 'limit': 1000})
                    activities.extend(part); ah += ';' + ph
                assays, targets, documents = {}, {}, {}
                for field, endpoint, collection, idkey in [('assay_chembl_id','assay','assays','assay_chembl_id'),('target_chembl_id','target','targets','target_chembl_id'),('document_chembl_id','document','documents','document_chembl_id')]:
                    ids = sorted({a.get(field, '') for a in activities if a.get(field)})
                    dest = {'assay': assays, 'target': targets, 'document': documents}[endpoint]
                    for chunk in split_batches(ids, cfg['batch_sizes']['chembl_metadata']):
                        data, ph = client.pages(root + endpoint + '.json', collection, {idkey + '__in': ','.join(chunk), 'limit': 1000})
                        dest.update({r[idkey]: r for r in data}); ah += ';' + ph
                count_by_key = Counter()
                for a in activities:
                    q = q_by_mol[a['molecule_chembl_id']]; t = targets.get(a.get('target_chembl_id'), {}); ass = assays.get(a.get('assay_chembl_id'), {})
                    accs = [c.get('accession', '') for c in t.get('target_components', [])]
                    mapped = [uniprot[x] for x in accs if x in uniprot]
                    symbol = mapped[0]['symbol'] if len(mapped) == 1 and len(accs) == 1 else ''
                    doc = documents.get(a.get('document_chembl_id'), {})
                    r = base_evidence(q,source)
                    r.update({'identity_exact': mol_smiles[a['molecule_chembl_id']] == q['smiles'], 'actual_compound_id': a['molecule_chembl_id'], 'actual_smiles': mol_smiles[a['molecule_chembl_id']], 'target_symbol': symbol, 'target_stable_id': 'UniProt:' + accs[0] if len(accs)==1 else '', 'species': '9606' if t.get('organism') == 'Homo sapiens' else '', 'target_name': t.get('pref_name',''), 'single_protein': t.get('target_type') == 'SINGLE PROTEIN', 'target_type': t.get('target_type',''), 'assay_mode': 'binding' if ass.get('assay_type') == 'B' else 'functional_or_phenotypic', 'confidence': ass.get('confidence_score',''), 'assay_id': a.get('assay_chembl_id',''), 'activity_id': str(a.get('activity_id','')), 'activity_type': a.get('standard_type') or '', 'value': a.get('standard_value') or '', 'units': a.get('standard_units') or '', 'relation': a.get('standard_relation') or '', 'document': doc.get('doi') or (('PMID:'+str(doc['pubmed_id'])) if doc.get('pubmed_id') else ''), 'document_id': a.get('document_chembl_id') or '', 'assay_description': ass.get('description',''), 'validity_comment': a.get('data_validity_comment') or '', 'potential_duplicate': a.get('potential_duplicate'), 'activity_outcome': 'quantitative_record', 'response_sha256': mh+ah})
                    r['source_activity_comment']=a.get('activity_comment') or ''
                    r['activity_outcome']=chembl_outcome(r['source_activity_comment'])
                    r['tier'], r['tier_reason'] = classify(r); evidence.append(r); count_by_key[q['compound_key']] += 1
                for q in batch:
                    identity_mols = mol_by_key[q['inchikey']]
                    done.append(ledger(q,source,'SUCCESS_WITH_RECORDS' if count_by_key[q['compound_key']] else 'SUCCESS_NO_RECORDS', f'exact_key_molecules={len(identity_mols)}; human_activity_rows={count_by_key[q["compound_key"]]}', mh+ah,True))
                checkpoint(out,source,done,evidence)
                print(f'{source} {len(done)}/{len(queries)} evidence={len(evidence)}',flush=True)
    except (SourceError, ValueError, KeyError) as e:
        covered={r['compound_key'] for r in done}
        done.extend(ledger(q,source,'SOURCE_UNAVAILABLE_NOT_QUERIED',str(e)[:350]) for q in queries if q['compound_key'] not in covered)
        checkpoint(out,source,done,evidence)
        print(source+' blocked: '+str(e)[:200],flush=True)


def bindingdb_fetch(cfg,out,raw,queries,uniprot):
    source='BindingDB'; client=Client(source,cfg,raw); done=[]; evidence=[]; failures=0
    previous_path=out/f'ledger_{source}.tsv'
    old={r['compound_key']:r for r in readtable(previous_path)} if previous_path.exists() else {}
    old_evidence=json.loads((out/f'evidence_{source}.json').read_text(encoding='utf-8')) if (out/f'evidence_{source}.json').exists() else []
    old_by_key=defaultdict(list)
    for r in old_evidence: old_by_key[r['compound_key']].append(r)
    for q in queries:
        if old.get(q['compound_key'],{}).get('status','').startswith('SUCCESS'):
            done.append(old[q['compound_key']]); evidence.extend(old_by_key[q['compound_key']]); continue
        if failures>=cfg['http']['circuit_breaker_consecutive_failures']:
            done.append(ledger(q,source,'SOURCE_UNAVAILABLE_NOT_QUERIED','source circuit open after exhausted requests')); continue
        try:
            data,meta=client.json(cfg['sources'][source],{'smiles':q['smiles'],'cutoff':'1.0','response':'application/json'})
            block=data.get('getLindsByUniprotResponse',{}) if isinstance(data,dict) else {}
            if not block: raise SourceError('unrecognised BindingDB response schema')
            records=block.get('bdb.affinities') or []
            if isinstance(records,dict): records=[records]
            for a in records:
                r=base_evidence(q,source); actual=canonical(a.get('bdb.smiles',''))
                accession=a.get('bdb.uniprot',''); mapped=uniprot.get(accession,{})
                val=str(a.get('bdb.affinity','')).strip(); match=re.fullmatch(r'\s*([<>]=?|=)?\s*([\d.eE+-]+)\s*',val)
                r.update({'identity_exact': actual==q['smiles'], 'actual_smiles':actual,'actual_compound_id':str(a.get('bdb.monomerid','')), 'target_symbol':mapped.get('symbol',''), 'target_stable_id':('UniProt:'+accession) if accession else '', 'target_name':a.get('bdb.target',''), 'species':'9606' if a.get('bdb.species') in ('Human','Homo sapiens') else a.get('bdb.species',''), 'single_protein':bool(accession), 'assay_mode':'binding', 'value':match.group(2) if match else val,'relation':(match.group(1) or '=') if match else '', 'units':a.get('bdb.units',''), 'activity_type':a.get('bdb.affinity_type',''), 'assay_id':a.get('bdb.reactant_set_id',''), 'activity_id':a.get('bdb.reactant_set_id',''), 'document':a.get('bdb.doi',''), 'activity_outcome':'quantitative_record','response_sha256':meta['sha256']})
                r['tier'],r['tier_reason']=classify(r); evidence.append(r)
            done.append(ledger(q,source,'SUCCESS_WITH_RECORDS' if records else 'SUCCESS_NO_RECORDS','API returns target names; missing stable IDs/units/documents remain E2',meta['sha256'],True)); failures=0
        except (SourceError,ValueError,KeyError) as e:
            failures+=1
            status='PARSE_FAILED' if 'parse failure' in str(e) or 'response schema' in str(e) else 'QUERY_FAILED'
            done.append(ledger(q,source,status,str(e)[:350]))
        if len(done)%10==0:
            checkpoint(out,source,done,evidence); print(f'{source} {len(done)}/{len(queries)} evidence={len(evidence)}',flush=True)
    checkpoint(out,source,done,evidence)


def pubchem_rows(data):
    if data is None: return []
    t=data.get('Table',{})
    columns=t.get('Columns',{}).get('Column')
    if columns is None or 'Row' not in t: raise SourceError('unrecognised PubChem assay summary schema')
    return [dict(zip(columns,r.get('Cell',[]))) for r in t['Row']]


def pubchem_fetch(cfg,out,raw,queries,entrez):
    source='PubChem'; client=Client(source,cfg,raw); done=[]; evidence=[]
    roots=cfg['sources'][source]; q_by_cid=defaultdict(list)
    for q in queries:
        for cid in q.get('confirmed_cids','').split(';'):
            if cid: q_by_cid[cid].append(q)
    cid_order=sorted(q_by_cid,key=lambda cid:(min(q['phase'] for q in q_by_cid[cid]),int(cid)))
    completed_cids=set(); statuses={}
    def get_batch(cids):
        try:
            data,meta=client.json(roots+'compound/cid/'+','.join(cids)+'/assaysummary/JSON',allow404=True)
        except SourceError:
            if len(cids)>1:
                mid=len(cids)//2
                return get_batch(cids[:mid])+get_batch(cids[mid:])
            raise
        if data is None and len(cids)>1:
            merged=[]
            for cid in cids: merged.extend(get_batch([cid]))
            return merged
        rows=pubchem_rows(data)
        for cid in cids: statuses[cid]={'status':'SUCCESS_WITH_RECORDS' if any(str(r.get('CID'))==cid for r in rows) else 'SUCCESS_NO_RECORDS','sha':meta['sha256']}; completed_cids.add(cid)
        return [(r,meta['sha256']) for r in rows]
    try:
        for batch in split_batches(cid_order,cfg['batch_sizes']['pubchem_assays']):
            for a,hs in get_batch(batch):
                cid=str(a.get('CID','')); gene=str(a.get('Target GeneID','')); mapping=entrez.get(gene,{})
                if not gene and not a.get('Target Accession'): continue
                # Summaries retain all outcomes; original descriptions are reviewed
                # separately for eligible candidate intersections.
                aid=str(a.get('AID','')); desc={}
                for q in q_by_cid.get(cid,[]):
                    r=base_evidence(q,source)
                    r.update({'identity_exact': len({x['smiles'] for x in q_by_cid[cid]})==1,'actual_compound_id':'CID:'+cid,'target_symbol':mapping.get('symbol',''),'target_stable_id':'NCBI_Gene:'+gene if gene else '', 'species':'9606' if mapping else '', 'single_protein':False, 'target_accession':a.get('Target Accession',''), 'assay_mode':'unverified_assay_context','assay_id':'AID:'+aid,'activity_id':aid+':SID:'+str(a.get('SID',''))+':panel:'+str(a.get('Panel Member ID','')),'activity_outcome':a.get('Activity Outcome',''),'value':a.get('Activity Value [uM]',''),'units':'uM' if a.get('Activity Value [uM]') else '', 'relation':'=' if a.get('Activity Value [uM]') else '', 'activity_type':a.get('Activity Name',''), 'document':'PMID:'+str(a['PubMed ID']) if a.get('PubMed ID') else '', 'assay_description':a.get('Assay Name',''), 'assay_type':a.get('Assay Type',''),'description_retrieved':bool(desc),'response_sha256':hs,'RNAi':a.get('RNAi','')})
                    r['tier'],r['tier_reason']=classify(r); evidence.append(r)
            done=[]
            for q in queries:
                cids=[x for x in q.get('confirmed_cids','').split(';') if x]
                if cids and all(c in completed_cids for c in cids):
                    done.append(ledger(q,source,'SUCCESS_WITH_RECORDS' if any(statuses[c]['status']=='SUCCESS_WITH_RECORDS' for c in cids) else 'SUCCESS_NO_RECORDS','CID assay summary; E1 requires target-specific original assay adjudication',';'.join(statuses[c]['sha'] for c in cids),True))
            checkpoint(out,source,done,evidence); print(f'{source} CIDs={len(completed_cids)}/{len(cid_order)} structures={len(done)}/{len(queries)} evidence={len(evidence)}',flush=True)
    except (SourceError,ValueError,KeyError) as e:
        print(source+' blocked: '+str(e)[:200],flush=True)
    covered={r['compound_key'] for r in done}
    done.extend(ledger(q,source,'IDENTITY_HOLD' if not q.get('confirmed_cids') else 'SOURCE_UNAVAILABLE_NOT_QUERIED','no source-confirmed CID/structure' if not q.get('confirmed_cids') else 'query/description/pagination incomplete; consult raw failure log') for q in queries if q['compound_key'] not in covered)
    checkpoint(out,source,done,evidence)


def resolve(cfg,out,raw):
    client=Client('PubChem_structure',cfg,raw)
    candidates=readtable(ROOT/cfg['inputs']['candidates'])
    cid_candidates={r['PubChem_CID']:r for r in candidates if r['PubChem_CID']}
    existing=readtable(out/'compound_queries.tsv')
    by_smiles={r['smiles']:r for r in existing}
    identity=[]; resolved_cids=set(); properties={}
    def batch_get(cids):
        data,meta=client.json(cfg['sources']['PubChem']+'compound/cid/'+','.join(cids)+'/property/InChIKey,SMILES/JSON',allow404=True)
        if data is None and len(cids)>1:
            for cid in cids: batch_get([cid])
            return
        found={str(r['CID']):r for r in (data or {}).get('PropertyTable',{}).get('Properties',[])}
        for cid in cids:
            resolved_cids.add(cid)
            p=found.get(cid,{})
            s=canonical(p.get('SMILES') or p.get('IsomericSMILES') or '')
            properties[cid]={'PubChem_CID':cid,'source_smiles':p.get('SMILES') or p.get('IsomericSMILES') or '', 'canonical_smiles':s,'inchikey':p.get('InChIKey',''),'status':'SUCCESS_WITH_STRUCTURE' if s else 'CID_NOT_FOUND_OR_NO_PARSEABLE_STRUCTURE','response_sha256':meta['sha256']}
    error=''
    try:
        for batch in split_batches(sorted(cid_candidates,key=int),cfg['batch_sizes']['pubchem_properties']):
            batch_get(batch)
            writetable(out/'pubchem_structure_resolution.tsv',list(properties.values()))
            print(f'PubChem structure CID={len(resolved_cids)}/{len(cid_candidates)}',flush=True)
    except (SourceError,ValueError,KeyError) as e:
        error=str(e)[:350]
    for cid in cid_candidates:
        if cid not in resolved_cids:
            properties[cid]={'PubChem_CID':cid,'status':'SOURCE_UNAVAILABLE_NOT_QUERIED','detail':error}
    for cid,p in properties.items():
        s=p.get('canonical_smiles','')
        if not s: continue
        candidate=cid_candidates[cid]
        if s not in by_smiles:
            q={'compound_key':'STR:'+digest(s.encode())[:24],'smiles':s,'inchikey':p['inchikey'],'cids':cid,'herbs':candidate['Herbs'],'names':candidate['Candidate_Name'],'phase':'3','identity_status':'public_CID_structure_resolved'}
            by_smiles[s]=q
        else:
            q=by_smiles[s]
            q['cids']=';'.join(sorted(set(q['cids'].split(';'))|{cid} - {''}))
            q['identity_status']='public_CID_confirmed_structure'
    all_queries=sorted(by_smiles.values(),key=lambda r:(r['phase'],r['compound_key']))
    writetable(out/'compound_queries_resolved.tsv',all_queries)
    writetable(out/'pubchem_structure_resolution.tsv',list(properties.values()))
    for q in all_queries:
        verified=[c for c in q['cids'].split(';') if properties.get(c,{}).get('canonical_smiles')==q['smiles']]
        contradicted=[c for c in q['cids'].split(';') if properties.get(c,{}).get('canonical_smiles') and properties[c]['canonical_smiles']!=q['smiles']]
        identity.append({**q,'confirmed_cids':';'.join(verified),'conflicting_cids':';'.join(contradicted)})
    writetable(out/'compound_identity_audit.tsv',identity)
    print(json.dumps({'structure_resolution':Counter(r['status'] for r in properties.values()),'query_structures_after_resolution':len(all_queries)},ensure_ascii=False),flush=True)


def fetch(cfg,out,raw,only='all'):
    client=Client('HGNC',cfg,raw)
    body,meta=client.get(cfg['sources']['HGNC'])
    hgnc=raw/'HGNC'/'hgnc_complete_set.tsv'
    if not hgnc.exists(): hgnc.write_bytes(body)
    entrez,uniprot=load_hgnc(hgnc)
    queries=readtable(out/'compound_identity_audit.tsv') if (out/'compound_identity_audit.tsv').exists() else readtable(out/'compound_queries.tsv')
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks=[]
        if only in ('all','ChEMBL'): tasks.append(pool.submit(chembl_fetch,cfg,out,raw,queries,uniprot))
        if only in ('all','BindingDB'): tasks.append(pool.submit(bindingdb_fetch,cfg,out,raw,queries,uniprot))
        if only in ('all','PubChem'): tasks.append(pool.submit(pubchem_fetch,cfg,out,raw,queries,entrez))
        for task in tasks: task.result()


def main():
    p=argparse.ArgumentParser(); p.add_argument('stage',choices=['prepare','resolve','fetch']); p.add_argument('--source',default='all',choices=['all','ChEMBL','BindingDB','PubChem']); a=p.parse_args()
    cfg=json.loads(CONFIG.read_text(encoding='utf-8')); out=ROOT/'results'/cfg['run_id']; raw=ROOT/'data/raw'/cfg['run_id']
    out.mkdir(parents=True,exist_ok=True)
    history=out/'code_revision_log.json'
    revisions=json.loads(history.read_text(encoding='utf-8')) if history.exists() else []
    revisions.append({'stage':a.stage,'time':now(),'code_sha256':filehash(Path(__file__)),'configuration_sha256':filehash(CONFIG)})
    writejson(history,revisions)
    if a.stage=='prepare': prepare(cfg,out)
    elif a.stage=='resolve': resolve(cfg,out,raw)
    elif a.stage=='fetch': fetch(cfg,out,raw,a.source)


if __name__=='__main__': main()
