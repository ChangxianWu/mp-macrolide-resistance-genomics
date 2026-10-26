# -*- coding: utf-8 -*-
"""28_country_norm.py — Country/region name normalization (submission convention)
Rules (following Lancet/PMCID conventions):
  - Taiwan -> "Taiwan, China"; China -> "Mainland China" (unified across main text/figures/tables; raw NCBI text is retained in the country column of Supplementary Table S1)
  - USA -> "United States"
  - "not applicable" -> "Not specified"
Implementation points: rewrite the Region column of the per-country resistance-frequency table under results/ and manuscript/tables/table1_panel.csv;
        other tables are anchored to table1; mentions in the main text of the manuscript config are replaced as well.
Input: the integrated per-genome typing matrix under results/ (the raw country column is left untouched)
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
    """Normalize ST display names (review P2; same convention as MY9_24): uncertain(near X) -> X-related, novel(near X) -> X-novel."""
    if k.startswith("uncertain(near "):
        return k.replace("uncertain(near ", "").rstrip(")") + "-related"
    if k.startswith("novel(near "):
        return k.replace("novel(near ", "").rstrip(")") + "-novel"
    return k

mat = list(csv.DictReader(open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t"))

# ---- Rewrite table1_panel.csv (normalized Region names) ----
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

# ---- Rewrite the per-country frequency table (normalized, for internal consistency) ----
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

# ---- Manuscript config wording: "in China / in Taiwan" -> canonical names ----
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
