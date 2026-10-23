"""Cache official PDB metadata and coordinates for the expanded docking run."""
import concurrent.futures
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/2026-10-04-expanded-docking-20'
RAW.mkdir(parents=True, exist_ok=True)
IDS = ['6ESM','4JLG','3ANR','1X70','3HS4','4NRE','6NCF','8GHB','4C1M','7ALV','5EQI','8VT4','4B2D','6M61','5IKR','8X7L','1NFI']

def fetch(pdb):
    try:
        def get(url, name):
            p = RAW / name
            if not p.exists():
                r = requests.get(url, timeout=45)
                r.raise_for_status()
                p.write_bytes(r.content)
            return p
        p = get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}', f'{pdb}_entry.json')
        obj = json.loads(p.read_text(encoding='utf-8'))
        get(f'https://files.rcsb.org/download/{pdb}.pdb', f'{pdb}.pdb')
        entities = []
        for eid in obj['rcsb_entry_container_identifiers']['polymer_entity_ids']:
            p = get(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}', f'{pdb}_entity{eid}.json')
            e = json.loads(p.read_text(encoding='utf-8'))
            entities.append({'description':e['rcsb_polymer_entity']['pdbx_description'],
                'chains':e['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'],
                'uniprot':e['rcsb_polymer_entity_container_identifiers'].get('uniprot_ids',[]),
                'mutation':e['entity_poly'].get('pdbx_mutation',''),
                'taxa':[o.get('ncbi_taxonomy_id') for o in e.get('rcsb_entity_source_organism',[])]})
        atoms = (RAW/f'{pdb}.pdb').read_text().splitlines()
        hets = sorted({(l[17:20].strip(),l[21],l[22:26].strip()) for l in atoms if l.startswith('HETATM') and l[17:20]!='HOH'})
        return {'pdb':pdb,'title':obj['struct']['title'],'resolution':obj['rcsb_entry_info'].get('resolution_combined'),
                'entities':entities,'heterogens':hets,'doi':obj.get('rcsb_primary_citation',{}).get('pdbx_database_id_DOI'),
                'publication_doi':obj.get('rcsb_primary_citation',{}).get('pdbx_database_id_DOI'),
                'citation':obj.get('rcsb_primary_citation',{})}
    except Exception as exc:
        return {'pdb':pdb,'error':str(exc)}

if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for r in pool.map(fetch, IDS):
            dest=ROOT/'data/interim/2026-10-04-expanded-docking-20'
            dest.mkdir(parents=True,exist_ok=True)
            (dest/f"{r['pdb']}_inventory.json").write_text(json.dumps(r,indent=2),encoding='utf-8')
            print(json.dumps({k:v for k,v in r.items() if k!='citation'},ensure_ascii=False),flush=True)
