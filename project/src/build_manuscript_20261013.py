"""Build and audit the manuscript from frozen local results; never rerun science.

Run with the bundled Python recorded in configs/manuscript_20261013.json.
Only this manuscript source, its output directory and its own trace are written.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import html
import json
import math
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import reportlab
from pypdf import PdfReader
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = ROOT / "configs/manuscript_20261013.json"
CFG = json.loads(CFG_PATH.read_text(encoding="utf-8"))
OUT = ROOT / CFG["output_dir"]
SOURCE = ROOT / CFG["source"]
RUN = CFG["run_id"]
G1 = "results/archive/2026-09-20-disease-side-v3.1-g2-stop-historical"
G5 = "results/2026-09-22-disease-side-v3.4-g5"
A30 = "results/2026-10-03-stage-a-mechanism-030"
R36 = "results/2026-10-11-virtual-ko-matching-correction-036"
R35 = "results/2026-10-11-virtual-ko-scientific-review-035"
LINK = re.compile(r"\[([^\]]+)\]\((https?://\S+?)\)(?=[.,;\s]|$)")
CITE = re.compile(r"\[(\d+(?:[–-]\d+)?(?:,\s*\d+(?:[–-]\d+)?)*)\]")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read_rows(relative: str) -> list[dict]:
    p = ROOT / relative
    op = gzip.open if p.suffix == ".gz" else open
    with op(p, "rt", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="," if p.suffix == ".csv" else "\t"))


def load(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8-sig"))


def expand(numbers: str) -> list[int]:
    result = []
    for token in numbers.split(","):
        ns = re.split("[–-]", token.strip())
        result.extend(range(int(ns[0]), int(ns[-1]) + 1))
    return result


def order_references(text: str) -> tuple[str, dict]:
    body, end = text.split("## References\n", 1)
    refs, chinese = end.split("## 中文摘要\n", 1)
    entries = dict((int(n), paragraph) for n, paragraph in re.findall(
        r"(?m)^(\d+)\. (.*?)(?=\n\n\d+\. |\Z)", refs, re.S))
    order = []
    for match in CITE.finditer(body):
        for n in expand(match[1]):
            if n not in order:
                order.append(n)
    if set(order) != set(entries):
        raise ValueError(f"Cited/reference mismatch: {set(order) ^ set(entries)}")
    mapping = {old: i + 1 for i, old in enumerate(order)}
    body = CITE.sub(lambda m: "[" + ",".join(str(mapping[n]) for n in expand(m[1])) + "]", body)
    reference_text = "\n\n".join(f"{mapping[n]}. {entries[n].strip()}" for n in order)
    return body + "## References\n\n" + reference_text + "\n\n## 中文摘要\n" + chinese, mapping


def rich(text: str, pdf: bool = False) -> str:
    escaped = html.escape(text)
    escaped = LINK.sub(lambda m: f'<a href="{m[2]}">{m[1]}</a>', escaped)
    tag = "super" if pdf else "sup"
    escaped = re.sub(r"10\^(-\d+)", lambda m: f"10<{tag}>{m[1]}</{tag}>", escaped)
    if pdf:
        escaped = re.sub(r"[\u4e00-\u9fff]+", lambda m: f'<font name="ChineseLocal">{m[0]}</font>', escaped)
    return escaped


def blocks(text: str):
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [x.strip() for x in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", x) for x in row):
                    rows.append(row)
                i += 1
            yield "table", rows
            continue
        if line.startswith("#"):
            n = len(line) - len(line.lstrip("#"))
            yield "h" + str(n), line[n:].strip()
            i += 1
            continue
        if line.startswith("- "):
            yield "bullet", line[2:]
            i += 1
            continue
        group = [line]
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "|", "- ")):
            group.append(lines[i].strip())
            i += 1
        yield "p", " ".join(group)


def export_html(text: str):
    parts = []
    for kind, value in blocks(text):
        if kind == "table":
            parts.append('<div class="table-wrap"><table>')
            for i, row in enumerate(value):
                tag = "th" if i == 0 else "td"
                parts.append("<tr>" + "".join(f"<{tag}>{rich(c)}</{tag}>" for c in row) + "</tr>")
            parts.append("</table></div>")
        elif kind == "bullet":
            parts.append(f'<p class="point">• {rich(value)}</p>')
        elif kind.startswith("h"):
            parts.append(f"<{kind}>{rich(value)}</{kind}>")
        else:
            parts.append(f"<p>{rich(value)}</p>")
    css = """body{max-width:880px;margin:50px auto;padding:0 28px;color:#171717;
    font:18px/1.6 Georgia,'Times New Roman','Microsoft YaHei',serif}
    h1{font-size:32px;line-height:1.25;margin-bottom:32px}h2{font-size:25px;margin-top:36px}
    h3{font-size:20px;margin-top:24px}p{margin:0 0 16px}a{color:#244b73}
    table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.4;margin:10px 0 18px}
    th,td{border-bottom:1px solid #d5d5d5;text-align:left;padding:9px 8px;vertical-align:top}
    th{background:#f1f3f5}.table-wrap{overflow-x:auto}.point{padding-left:18px}
    @media print{body{font-size:11pt;max-width:none;margin:0}h2,h3{break-after:avoid}tr{break-inside:avoid}}
    """
    title = next(v for k, v in blocks(text) if k == "h1")
    (OUT / "manuscript_en.html").write_text(
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
        'content="width=device-width,initial-scale=1"><title>' + html.escape(title) +
        "</title><style>" + css + "</style></head><body><main>" + "\n".join(parts) +
        "</main></body></html>", encoding="utf-8")


def export_pdf(text: str):
    for name, fn in [("TimesLocal", "times.ttf"), ("TimesBoldLocal", "timesbd.ttf")]:
        pdfmetrics.registerFont(TTFont(name, "C:/Windows/Fonts/" + fn))
    pdfmetrics.registerFont(TTFont("ChineseLocal", "C:/Windows/Fonts/msyh.ttc", subfontIndex=0))
    styles = {
        "h1": ParagraphStyle("Title", fontName="TimesBoldLocal", fontSize=19, leading=23,
                             spaceAfter=20, keepWithNext=True),
        "h2": ParagraphStyle("Section", fontName="TimesBoldLocal", fontSize=14, leading=18,
                             spaceBefore=15, spaceAfter=8, keepWithNext=True),
        "h3": ParagraphStyle("Subsection", fontName="TimesBoldLocal", fontSize=12, leading=15,
                             spaceBefore=10, spaceAfter=6, keepWithNext=True),
        "p": ParagraphStyle("Body", fontName="TimesLocal", fontSize=11, leading=14,
                            spaceAfter=8, alignment=TA_JUSTIFY, splitLongWords=False),
        "table": ParagraphStyle("Cell", fontName="TimesLocal", fontSize=9, leading=11,
                                alignment=TA_LEFT),
        "ref": ParagraphStyle("Reference", fontName="TimesLocal", fontSize=10, leading=13,
                              spaceAfter=7, alignment=TA_LEFT, splitLongWords=False),
        "cn": ParagraphStyle("Chinese", fontName="ChineseLocal", fontSize=10.5, leading=16,
                             spaceAfter=8, wordWrap="CJK"),
    }
    flow = []
    in_refs = False
    in_cn = False
    for kind, value in blocks(text):
        if kind == "h2":
            if value == "中文摘要":
                flow.append(PageBreak())
            in_refs = value == "References"
            in_cn = value == "中文摘要"
        if kind == "table":
            widths = [86, 103, 52, 52, 89, 105] if len(value[0]) == 6 else [143, 63, 69, 122, 90]
            cells = [[Paragraph(rich(c, True), styles["table"]) for c in row] for row in value]
            table = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
            table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef1f4")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#66717b")),
                ("LINEBELOW", (0, 1), (-1, -1), 0.25, colors.HexColor("#c8c8c8")),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]))
            flow.extend([table, Spacer(1, 8)])
            continue
        if in_cn:
            style = styles["cn"]
        elif kind.startswith("h"):
            style = styles[kind]
        else:
            style = styles["ref"] if in_refs else styles["p"]
        p = Paragraph(rich(value, True), style, bulletText="•" if kind == "bullet" else None)
        flow.append(KeepTogether([p]) if in_refs and kind == "p" else p)

    def page(canvas, doc):
        canvas.saveState()
        canvas.setFont("TimesLocal", 9)
        canvas.setFillColor(colors.HexColor("#555555"))
        canvas.drawRightString(A4[0] - 54, 30, str(doc.page))
        canvas.restoreState()

    title = next(v for k, v in blocks(text) if k == "h1")
    doc = SimpleDocTemplate(str(OUT / "manuscript_en.pdf"), pagesize=A4,
        leftMargin=54, rightMargin=54, topMargin=48, bottomMargin=48,
        title=title, author="", subject="Human glial associations and YZDHD candidate evaluation")
    doc.build(flow, onFirstPage=page, onLaterPages=page)


def audit(text: str) -> dict:
    tests = {}
    def check(name: str, observed, expected):
        tests[name] = {"observed": observed, "expected": expected, "pass": observed == expected}

    rows = read_rows(f"{G1}/meta/random_effects_meta_summary.tsv")
    hits = [r for r in rows if float(r["meta_FDR"]) < 0.05]
    full = [r for r in hits if float(r["k"]) == 5]
    check("G1_meta_genes", len(rows), 19738)
    check("G1_FDR_hits", len(hits), 413)
    check("G1_complete_region_hits", len(full), 244)
    check("G1_positive", sum(float(r["meta_logFC"]) > 0 for r in full), 128)
    check("G1_negative", sum(float(r["meta_logFC"]) < 0 for r in full), 116)
    genes = {r["gene"] for r in full}
    loco = [r for r in read_rows(f"{G1}/meta/leave_one_cohort_out.tsv") if r["gene"] in genes]
    check("G1_LOCO_rows", len(loco), 244 * 3)
    check("G1_LOCO_direction_preserved", all(r["direction_stable"].lower() == "true" for r in loco), True)
    g5 = read_rows(f"{G5}/tables/GSE261698_primary_module_effects.tsv")
    check("G5_modules", len(g5), 7)
    check("G5_significant", sum(r["status"] != "non_estimable" and float(r["FDR_BH_7"]) < .05 for r in g5), 5)
    check("G5_non_estimable", sum(r["status"] == "non_estimable" for r in g5), 1)
    for row in g5:
        if row["status"] != "non_estimable":
            check("G5_effect_print_" + row["module"],
                  f'{float(row["logFC"]):.4f}'.replace("-", "−") in text, True)
    chem = read_rows("results/2026-09-30-fixed-formula-admet-recovery-021/formula_admet_predictions.tsv")
    check("ADMET_input", len(chem), 841)
    check("ADMET_pass", sum(r["Step7_Status"] == "pass_predictive_heuristic" for r in chem), 335)
    pairs = read_rows("results/2026-10-04-g1-g5-docking-reuse/E1成分靶点去重清单.csv")
    check("E1_unique_pair_rows", len(pairs), 80)
    a = load(f"{A30}/result_summary.json")
    for key, expected in [("GSEA_main_tests", 6039), ("GSEA_main_global_FDR_significant", 33),
                          ("GSEA_significant_all_three", 13), ("network_primary_significant", 43)]:
        check(key, a[key], expected)
    final = load(f"{R36}/final_status.json")
    for key, expected in [("combinations_CP10K_corrected", 29), ("reference_lists_changed", 2),
                          ("corrected_reference_KOs", 1116), ("reused_reference_KOs", 1112),
                          ("expanded_KO_attempts", 102), ("source_assets_unchanged", 1202),
                          ("formal_specificity_P", False), ("drug_evidence_ceiling", "C0")]:
        check(key, final[key], expected)
    before_refs, refs_cn = text.split("## References", 1)
    reference_section = refs_cn.split("## 中文摘要", 1)[0]
    check("reference_count", len(re.findall(r"(?m)^\d+\. ", reference_section)), 23)
    check("cited_reference_count", len({n for m in CITE.finditer(before_refs) for n in expand(m[1])}), 23)
    check("placeholder_free", not re.search(r"TODO|TBD|INSERT|Evidence needed|AUTHOR_INPUT_NEEDED", text), True)
    # Compare prose only; titles in references intentionally reproduce bibliographic titles.
    ref_text = (ROOT / "tmp/manuscript/reference.txt").read_text(encoding="utf-8")
    words = lambda s: re.findall(r"[a-z]+", s.lower())
    n = 12
    original = words(ref_text)
    reference_ngrams = {tuple(original[i:i+n]) for i in range(len(original)-n+1)}
    manuscript_words = words(before_refs)
    shared = sorted({" ".join(manuscript_words[i:i+n]) for i in range(len(manuscript_words)-n+1)
                     if tuple(manuscript_words[i:i+n]) in reference_ngrams})
    check("reference_prose_shared_12grams", len(shared), 0)
    result = {"status": "PASS" if all(t["pass"] for t in tests.values()) else "FAIL",
              "scope": "saved-result numeric assertions and manuscript integrity; no biological rerun",
              "tests": tests, "shared_12grams": shared,
              "English_words_before_references": len(re.findall(r"\b[A-Za-z]+\b", before_refs))}
    return result


def provenance(sources: list[str], reference_map: dict):
    git = "C:/Users/administer/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/git/cmd/git.exe"
    parent = subprocess.check_output([git, "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    all_sources = sources + [CFG["source"], "docs/manuscript/editorial_audit_zh.md", "docs/manuscript/claim_evidence.tsv",
        "configs/manuscript_20261013.json", "src/build_manuscript_20261013.py",
        "docs/decisions/2026-10-13-manuscript-framing.md", CFG["reference_pdf"]]
    manifest = []
    for relative in dict.fromkeys(all_sources):
        p = Path(relative) if Path(relative).is_absolute() else ROOT / relative
        manifest.append({"path": str(p), "bytes": p.stat().st_size, "sha256": sha(p)})
    payload = {"run_id": RUN, "time_utc": datetime.now(timezone.utc).isoformat(),
        "Git_parent": parent, "Git_new_commit": False,
        "Git_worktree_preexisting_changes": True, "source_files": manifest,
        "reference_number_map_at_build": reference_map,
        "ReportLab": reportlab.Version, "science_rerun": False,
        "compute": "local document construction; web only for reference verification",
        "strict_cross_machine_reproduction": False}
    (OUT / "provenance.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Extraction is a temporary derivative of the user-provided style reference.
    tmp = ROOT / "tmp/manuscript"
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "reference.txt").write_text("\n".join(p.extract_text() for p in PdfReader(CFG["reference_pdf"]).pages), encoding="utf-8")
    source_text = re.sub(r"10−(\d+)", r"10^-\1", SOURCE.read_text(encoding="utf-8"))
    text, mapping = order_references(source_text)
    SOURCE.write_text(text, encoding="utf-8")
    result = audit(text)
    (OUT / "content_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    if result["status"] != "PASS":
        failed = {k:v for k,v in result["tests"].items() if not v["pass"]}
        raise RuntimeError(f"Content audit failed: {failed}")
    shutil.copy2(SOURCE, OUT / SOURCE.name)
    export_html(text)
    export_pdf(text)
    supplement_sources = {
        "Table_S1_G1_244_signature.tsv.gz": "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/G1_microglia_signed_sender_244.tsv.gz",
        "Table_S2_iMGL_senescence.tsv": "results/2026-09-20-disease-side-v3.2-g3-model-scoped/tables/D2_model_temporal_evidence.tsv",
        "Table_S3_iMGL_iron_LPO.tsv": "results/2026-09-20-disease-side-v3.2-g4/tables/D2_G4_RNA_sample_level_module_effects.tsv",
        "Table_S4_astrocyte_modules.tsv": f"{G5}/tables/GSE261698_primary_module_effects.tsv",
        "Table_S5_E1_compound_target_pairs.csv": "results/2026-10-04-g1-g5-docking-reuse/E1成分靶点去重清单.csv",
        "Table_S6_candidate_structural_scope.csv": "results/2026-10-04-expanded-docking-20/20靶点分析汇总.csv",
        "Table_S7_virtual_KO_final_review.tsv": f"{R35}/target_verification_conclusions.tsv",
        "Table_S8_corrected_reference_comparison.tsv": f"{R36}/functional_response_comparison.tsv",
    }
    sdir = OUT / "supplementary_tables"
    sdir.mkdir(exist_ok=True)
    for dest, original in supplement_sources.items():
        shutil.copy2(ROOT / original, sdir / dest)
    sources = list(supplement_sources.values()) + [
        "README.md", "params.yaml", "environment.yml", "requirements.txt",
        "docs/protocol/disease_side_latest.md", f"{G1}/execution_report.md",
        f"{G1}/meta/random_effects_meta_summary.tsv", f"{G1}/meta/leave_one_cohort_out.tsv",
        "results/2026-09-20-disease-side-v3.2-g2/execution_report.md",
        "results/2026-09-20-disease-side-v3.2-g2/g2_decision.json",
        "results/2026-09-30-fixed-formula-drug-prior-020/report.md",
        "results/2026-09-30-fixed-formula-admet-recovery-021/formula_admet_predictions.tsv",
        "results/2026-10-03-fixed-formula-target-supplement-024-continuation/report.md",
        f"{A30}/report.md", f"{A30}/result_summary.json", f"{A30}/GSEA_sensitivity_comparison.tsv",
        "results/2026-10-04-stage-b-tcmsp-docking-032/README.md",
        "results/2026-10-04-expanded-docking-20/20靶点对接分析报告.md",
        f"{R35}/最终科学复核报告.md", f"{R36}/最小修正报告.md", f"{R36}/final_status.json",
        f"{R36}/completion_review.json", "docs/reports/2026-10-13-matching-correction-036-completion-review.md",
        "src/25_run_g1_v3_1.R", "src/30_run_g3_v3_2.py", "src/32_run_g4_v3_2.py",
        "src/47_run_g5_v3_4.R", "src/108_process_fixed_formula_prior_019.py",
        "configs/g5_v3.4_execution.yaml", "configs/virtual_ko_matching_correction_036.yaml",
        "F:/TCM select/config/pipeline_policy.json", "F:/TCM select/code/6_ADMET_filter.py"]
    provenance(sources, mapping)
    table_manifest = {name: {"source": p, "sha256": sha(ROOT / p)} for name,p in supplement_sources.items()}
    (sdir / "table_provenance.json").write_text(json.dumps(table_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (sdir / "README.md").write_text(
        "# Supplementary tables\n\nThese are byte-identical copies of frozen project outputs. "
        "Program membership, experimental-record evidence, computational response and biological validation "
        "must remain distinct. Table S7 records the 035 review; Table S8 supplements it with 036 matching "
        "correction and does not imply rerunning donor resampling. See table_provenance.json for original "
        "paths and SHA-256. Full upstream outputs remain in the project and are not duplicated here.\n",
        encoding="utf-8")
    trace_dir = ROOT / "workspaces" / RUN / "raw/traces"
    trace_dir.mkdir(parents=True, exist_ok=True)
    trace = {"run_id": RUN, "user_request": "Write completed project paper following bbaf623 language style; English with Chinese abstract",
        "backend": "existing project modules; current Codex writing session",
        "deterministic_content_audit": result["status"], "skill_evolution": False,
        "baseline_validation": None, "best_validation": None, "test": None,
        "accepted_iterations": 0, "rejected_iterations": 0, "active_skills_added": [],
        "cross_project_promotion": False, "biological_reanalysis": False}
    trace_path = trace_dir / "manuscript_build.json"
    if not trace_path.exists():
        trace_path.write_text(json.dumps(trace, ensure_ascii=False, indent=2), encoding="utf-8")
    audit_path = trace_dir.parent.parent / "audit.json"
    audit_path.write_text(json.dumps(trace, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"content_audit": result["status"], "checks": len(result["tests"]),
        "English_words": result["English_words_before_references"],
        "PDF_pages": len(PdfReader(OUT / "manuscript_en.pdf").pages),
        "supplementary_tables": len(supplement_sources)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
