"""Capture publicly linked BATMAN bulk TTIs and HERB access evidence."""
import importlib.util
import json
from pathlib import Path

spec=importlib.util.spec_from_file_location('targets024',Path(__file__).with_name('113_supplement_formula_targets_024.py'))
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
cfg=json.loads(m.CONFIG.read_text(encoding='utf-8'))
out=m.ROOT/'results'/cfg['run_id']; raw=m.ROOT/'data/raw'/cfg['run_id']
records=[]
urls={
 'BATMAN':[
  'http://bionet.ncpsb.org.cn/batman-tcm/',
  'http://bionet.ncpsb.org.cn/batman-tcm/js/713.c42ca4e9.js',
  'http://batman2.cloudna.cn/downloadApiFile/data/browser/known_browse_by_ingredients.txt.gz',
  'http://batman2.cloudna.cn/downloadApiFile/data/browser/predicted_browse_by_ingredients.txt.gz'],
 'HERB':['http://herb.ac.cn/v2','http://herb.ac.cn/static/umi.js?version=762']}
for source,links in urls.items():
 c=m.Client(source,cfg,raw)
 for url in links:
  try:
   body,meta=c.get(url)
   records.append({**meta,'status':'SUCCESS_CAPTURED','meaning':'publicly linked BATMAN bulk file' if url.endswith('.gz') else 'public website/access provenance only; page success is not target query completion'})
   print(source,len(body),url,flush=True)
  except m.SourceError as e:
   records.append({'source':source,'url':url,'status':'SOURCE_UNAVAILABLE','detail':str(e)[:350]})
m.writetable(out/'prior_source_access.tsv',records)
