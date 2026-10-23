"""Snapshot original manuscript reference identity using two public metadata registries."""
import concurrent.futures,json,pathlib,re,urllib.request,urllib.parse,time,hashlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'data/processed/FIG_YZDHD_20261014';OUT.mkdir(parents=True,exist_ok=True)
text=(ROOT/'docs/manuscript/manuscript_en.md').read_text(encoding='utf-8')
bib=text.split('## References')[1].split('## 中文摘要')[0]
entries=re.findall(r'^([0-9]+)\. (.+)$',bib,re.M)
def get(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'ResearchReferenceAudit/1.0'}),timeout=50) as r:return json.load(r)
        except Exception as exc:
            if attempt==2:return {'error':str(exc)}
            time.sleep(1+attempt)
def work(pair):
    num,entry=pair;doi=re.search(r'\[doi:(.+?)\]',entry).group(1)
    crurl='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
    epurl='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query='+urllib.parse.quote('DOI:"'+doi+'"')+'&format=json&resultType=core'
    cr=get(crurl);ep=get(epurl);m=cr.get('message',{})
    matches=ep.get('resultList',{}).get('result',[])
    dates={k:m[k] for k in ['published','published-print','published-online','issued'] if k in m}
    year=(m.get('published-print') or m.get('published') or m.get('issued') or {}).get('date-parts',[[None]])[0][0]
    return {'id':int(num),'original_entry':entry,'doi':doi,'canonical_title':m.get('title',[''])[0],'authors':m.get('author',[]),'year':year,'journal':m.get('container-title',[''])[0],'volume':m.get('volume',''),'issue':m.get('issue',''),'page':m.get('page',''),'article_number':m.get('article-number',''),'dates':dates,'crossref_url':crurl,'europepmc_url':epurl,'crossref':m,'europepmc':matches,'errors':{'crossref':cr.get('error'),'europepmc':ep.get('error')}}
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    rows=sorted(pool.map(work,entries),key=lambda x:x['id'])
(OUT/'original_reference_metadata_snapshot.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
for x in rows:
    print(json.dumps({k:x[k] for k in ['id','doi','canonical_title','year','journal','volume','issue','page','article_number','errors']},ensure_ascii=False))
    print('CRauthors',[(a.get('family'),a.get('given')) for a in x['authors']])
    print('EP',[(a.get('id'),a.get('title'),a.get('authorString'),a.get('pubYear'),a.get('journalInfo'),a.get('pageInfo'),a.get('commentCorrectionList')) for a in x['europepmc']])
