# -*- coding: utf-8 -*-
"""25_supplementary.py — 补充材料生成（补充表 S1）
生成: 投稿材料包/02_附件与清单/MY9_Supplementary_TableS1_genome_metadata.docx
      （按当前分型矩阵动态生成的完整元数据 + 分型结果, 与稿件引用闭环）
依赖: docx 模块从 CV2 项目复制（与稿件生成器一致）
"""
import csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
PACK = os.path.join(BASE, "投稿材料包")
ATT = os.path.join(PACK, "02_附件与清单")
os.makedirs(ATT, exist_ok=True)
# 补充表用 python-docx（补充件与其他课题同惯例; JS docx 仅用于主稿生成器）
import docx as pydocx
from docx.shared import Pt, Inches
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

mat = list(csv.DictReader(open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t"))
mat.sort(key=lambda r: (r["country"], r["year"], r["ena_acc"]))

doc = pydocx.Document()
# 横向 + 窄边距：11 列表格在纵向页被压缩成 30 页乱码（ScholarOne 校样实测），横向后单行不折行
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

# 元数据富集: 从下载清单 + 入组池 join (仅透传提交方原文, 不编造; 规则 11/33)
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
# 固定列宽（合计 9.8 in，横向可用宽度 10 in）：避免自动布局把登录号/年份挤成断行
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
# 跨页重复表头 + 行不跨页断开
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
