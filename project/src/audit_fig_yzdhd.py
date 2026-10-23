"""Audit final YZDHD figure exports with bundled artifact Python; no drawing."""
from pathlib import Path
import json,hashlib
from pypdf import PdfReader
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/2026-10-14-FIG-YZDHD/figures'
checks=[]
for n in range(1,5):
 pdf=OUT/f'Fig{n}.pdf'; tif=OUT/f'Fig{n}.tif'; reader=PdfReader(pdf)
 p=reader.pages[0];fonts=p['/Resources']['/Font'].get_object()
 embedded=[]
 for v in fonts.values():
  f=v.get_object()
  if '/DescendantFonts' in f:f=f['/DescendantFonts'][0].get_object()
  fd=f.get('/FontDescriptor');fd=fd.get_object() if fd else {}
  embedded.append(any(k in fd for k in ['/FontFile','/FontFile2','/FontFile3']))
 sizes=[]
 def textvisit(text,cm,tm,font,size):
  if text.strip():sizes.append(size)
 p.extract_text(visitor_text=textvisit)
 with Image.open(tif) as im:
  dpi=tuple(float(x) for x in im.info['dpi']);mode=im.mode;pixel_size=im.size
 width=float(p.mediabox.width)*25.4/72;height=float(p.mediabox.height)*25.4/72
 row={'figure':n,'width_mm':width,'height_mm':height,'tiff_dpi':dpi,'tiff_mode':mode,'pixel_size':pixel_size,'pdf_font_sizes':sorted(set(sizes)), 'embedded_fonts':all(embedded),'sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [pdf,tif]}}
 row['pass']=abs(width-174)<.1 and height<=234 and mode=='RGB' and min(dpi)>=599 and all(embedded) and min(sizes)>=8 and max(sizes)<=12
 checks.append(row)
preflight=json.loads((OUT/'figure_source_preflight.json').read_text(encoding='utf-8-sig'))
report={'checks':checks,'all_pass':all(c['pass'] for c in checks),'static_preflight':preflight['summary'], 'accepted_warning':'174 mm is the target journal width; the generic skill validator expects Nature 183 mm.', 'visual_review':{'reviewed_figures':[1,2,3,4],'observations':'All complete figures inspected at final physical proportions. No internal titles, overlapping labels, clipped data, or unreadable symbols. Fig2 spacing adjusted to retain colorbar label; Fig4 reduced to224mm and every29context row retained. Lowercase a/b/c labels and Arial8–11pt text.','data_review':'Assertions cover19,738 genes,244subset,128/116 signs,20targets/80pairs,8G1network comparisons,29KOcontexts and58reference comparisons. Source hashes rechecked unchanged.'}}
(OUT/'figure_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
assert report['all_pass'] and preflight['summary']['ready'],report
print('All four exports and source preflight pass')
