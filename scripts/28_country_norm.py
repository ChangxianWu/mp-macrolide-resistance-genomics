# -*- coding: utf-8 -*-
"""28_country_norm.py — 国家/地区名规范化（投稿口径）
规则（跟随 Lancet/PMCID 惯例）:
  - Taiwan → "Taiwan, China"; China → "Mainland China"（正文/图表/表格统一; NCBI 原文保留在 S1 补充表 country 列）
  - USA → "United States"
  - "not applicable" → "Not specified"
实施点: results/MY9_耐药频率_按国家.tsv 与 manuscript/tables/table1_panel.csv 的 Region 列重出,
        其余表以 table1 为锚; 稿件 config 中正文提及处一并替换。
输入: results/MY9_分型矩阵.tsv（原文 country 列不动）
"""
import csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
OUT = os.path.join(BASE, "manuscript", "tables")

NORM = {"Taiwan": "Taiwan, China", "China": "Mainland China", "USA": "United States",
        "not applicable": "Not specified"}
def norm(c):
    return NORM.get(c, c)

def disp_st(k):
    """ST 展示命名规范化（评审 P2; 与 MY9_24 同口径）: uncertain(near X)->X-related, novel(near X)->X-novel。"""
    if k.startswith("uncertain(near "):
        return k.replace("uncertain(near ", "").rstrip(")") + "-related"
    if k.startswith("novel(near "):
        return k.replace("novel(near ", "").rstrip(")") + "-novel"
    return k

mat = list(csv.DictReader(open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t"))

# ---- 重出 table1_panel.csv（规范化 Region 名） ----
from collections import defaultdict, Counter
import numpy as np
by_c = defaultdict(lambda: {"n": 0, "res": 0, "wt": 0, "years": [], "st": Counter()})
for r in mat:
    c = norm(r["country"].split(":")[0].strip() or "Not specified")
    d = by_c[c]
    d["n"] += 1
    if r["resist_23S"] == "resistant": d["res"] += 1
    elif r["resist_23S"] == "WT": d["wt"] += 1
    if r["year"].isdigit(): d["years"].append(int(r["year"]))
    d["st"][r["MLST_ST"] or "novel/uncertain"] += 1
rows = []
for c in sorted(by_c, key=lambda c: -by_c[c]["n"]):
    d = by_c[c]
    if d["n"] < 5 and c != "Not specified": continue
    if c == "Not specified" and d["n"] < 5: continue
    called = d["res"] + d["wt"]
    stx = ", ".join(f"ST{disp_st(k)} ({v})" for k, v in d["st"].most_common(3))
    rows.append([c, d["n"], d["res"], d["wt"],
                 f"{d['res']/called*100:.1f}" if called else "NA",
                 f"{np.median(d['years']):.0f}" if d["years"] else "NA", stx])
tot_res = sum(1 for r in mat if r["resist_23S"] == "resistant")
tot_wt = sum(1 for r in mat if r["resist_23S"] == "WT")
_allyears = [int(r["year"]) for r in mat if r["year"].isdigit()]
_allst = Counter(r["MLST_ST"] or "novel/uncertain" for r in mat if r["resist_23S"] in ("resistant", "WT"))
_allstx = ", ".join(f"ST{k} ({v})" for k, v in _allst.most_common(3))
rows.append(["All included genomes", len(mat), tot_res, tot_wt,
             f"{tot_res/(tot_res+tot_wt)*100:.1f}",
             f"{np.median(_allyears):.0f}" if _allyears else "NA", _allstx])
p = os.path.join(OUT, "table1_panel.csv")
with open(p, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Region", "N_genomes", "N_resistant", "N_wildtype", "Pct_resistant", "Median_year", "Top_MLST"])
    w.writerows(rows)
print("table1 rewritten with normalized names")

# ---- 重出耐药频率_按国家.tsv（规范化, 供内部一致） ----
by_c2 = defaultdict(lambda: [0, 0])
for r in mat:
    c = norm(r["country"].split(":")[0].strip() or "Not specified")
    if r["resist_23S"] == "resistant": by_c2[c][0] += 1
    elif r["resist_23S"] == "WT": by_c2[c][1] += 1
p2 = os.path.join(RES, "MY9_耐药频率_按国家.tsv")
with open(p2, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["country", "resist", "wt", "freq_resist"])
    for c in sorted(by_c2, key=lambda c: -(by_c2[c][0] + by_c2[c][1])):
        a, b = by_c2[c]
        w.writerow([c, a, b, f"{a/(a+b)*100:.3f}" if a + b else "NA"])
print("耐药频率_按国家.tsv rewritten")

# ---- 稿件 config 措辞: "in China / in Taiwan" → 规范名 ----
p3 = os.path.join(BASE, "manuscript", "manuscript_config.js")
s = open(p3, encoding="utf-8").read()
n = 0
for o, w_ in [('in Taiwan)', 'in Taiwan, China)'),
              ('in Taiwan ', 'in Taiwan, China '),
              ('China, Taiwan and South Korea', 'China, Taiwan (China), and South Korea'),
              ('across China, Taiwan and South Korea', 'across China, Taiwan (China), and South Korea'),
              ('spans China, Taiwan and South Korea', 'spans China, Taiwan (China), and South Korea'),
              ('("Taiwan")', '("Taiwan, China")')]:
    if o in s:
        s = s.replace(o, w_); n += 1
open(p3, "w", encoding="utf-8").write(s)
print(f"config wording updated: {n} sites")
print("NOTE: S1 保留 NCBI 原文 country（审计友好）; 规范化只作用于汇总层")
