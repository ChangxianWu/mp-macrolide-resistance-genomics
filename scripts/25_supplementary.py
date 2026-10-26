# -*- coding: utf-8 -*-
"""25_supplementary.py — Supplementary material generation (Supplementary Table S1)
Generates: MY9_Supplementary_TableS1_genome_metadata.docx in the submission package's attachments and checklist folder
      (full metadata + typing results generated dynamically from the current typing matrix, closing the loop with manuscript citations)
Dependency: the docx module is copied from the CV2 project (consistent with the manuscript generator)
"""
import csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
PACK = os.path.join(BASE, "投稿材料包")
ATT = os.path.join(PACK, "02_附件与清单")
os.makedirs(ATT, exist_ok=True)
# Supplementary tables use python-docx (consistent with the supplement convention of other projects; the JS docx is used only by the main-manuscript generator)
import docx as pydocx
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

mat = list(csv.DictReader(open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t"))
mat.sort(key=lambda r: (r["country"], r["year"], r["ena_acc"]))

doc = pydocx.Document()
# Landscape + narrow margins: an 11-column table compressed onto portrait pages produced 30 pages of garbled output (observed in ScholarOne proofs); landscape keeps every row unbroken
sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = Inches(11), Inches(8.5)
for attr in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
    setattr(sec, attr, Inches(0.5))
doc.add_heading("Supplementary Table S1. Genome panel: metadata and typing results", level=2)
doc.add_paragraph(
    "One row per included genome (n = %d). Country is the leading field of the geographic "
    "locality as submitted; \"not applicable\" entries are retained but excluded from "
    "geographic stratification. Resistance calls come from dual-anchor 23S rRNA typing at "
    "positions 2063, 2064 and 2617 (wild-type bases A, A and C). P1 type uses canonical "
    "k-mer counts against type-1 (M129) and type-2 (FH) references; sequence types come from "
    "the public eight-locus MLST scheme (nearest-neighbour fallback marked)."
    % len(mat)).style.font.size = Pt(9)

# Metadata enrichment: join from the download manifest + inclusion pool (pass through submitter text only, no fabrication; rules 11/33)
dl_map = {}
with open(os.path.join(BASE, "MY9_下载清单.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        dl_map[r["ena_acc"]] = r["assembly_accession"].split(".")[0]
pool = {}
with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        pool[r["assembly_accession"].split(".")[0]] = r

def meta_of(acc):
    p = pool.get(dl_map.get(acc, ""), {})
    host = (p.get("host") or "").strip()
    src = (p.get("isolation_source") or "").strip()
    hostsrc = host + ("; " + src if src and src != host else "")
    return (p.get("assembly_accession") or "").strip(), (p.get("year") or "").strip(), hostsrc

cols = ["ena_acc", "assembly_accession", "country", "year", "host_or_source",
        "site2063", "site2064", "site2617", "resist_23S", "P1_type", "MLST_ST"]
hdr = ["ENA accession", "Assembly", "Country", "Year", "Host / source",
       "2063", "2064", "2617", "23S genotype", "P1 type", "MLST ST"]
tbl = doc.add_table(rows=1, cols=len(hdr))
tbl.style = "Table Grid"
tbl.autofit = False
# Fixed column widths (9.8 in total, against 10 in usable landscape width): prevents automatic layout from squeezing accessions/years into wrapped lines
COLW = [1.1, 1.35, 1.15, 0.45, 2.35, 0.5, 0.5, 0.5, 0.75, 0.55, 0.6]
for i, w in enumerate(COLW):
    tbl.columns[i].width = Inches(w)
for i, h in enumerate(hdr):
    c = tbl.rows[0].cells[i]
    c.width = Inches(COLW[i])
    c.text = h
    for run in c.paragraphs[0].runs:
        run.bold = True
        run.font.size = Pt(8)
# Repeat the header across pages + keep rows from splitting across pages
trPr = tbl.rows[0]._tr.get_or_add_trPr()
trPr.append(OxmlElement("w:tblHeader"))
for r in mat:
    asm, yr, hostsrc = meta_of(r["ena_acc"])
    row = tbl.add_row()
    trPr = row._tr.get_or_add_trPr()
    trPr.append(OxmlElement("w:cantSplit"))
    for i, k in enumerate(cols):
        v = {"assembly_accession": asm or r.get("assembly_accession", ""),
             "year": yr or r.get("year", ""),
             "host_or_source": hostsrc or r.get("host_or_source", "")}.get(k, r.get(k, ""))
        row.cells[i].width = Inches(COLW[i])
        row.cells[i].text = str(v)
        for p in row.cells[i].paragraphs:
            for run in p.runs:
                run.font.size = Pt(7.5)

out = os.path.join(ATT, "MY9_Supplementary_TableS1_genome_metadata.docx")
doc.save(out)
print("written", out, len(mat), "rows")
