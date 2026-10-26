# -*- coding: utf-8 -*-
"""06_summary_tables.py — resistance mutation frequency × time × geography summary tables
Input: the per-genome 23S typing matrix under results/ (output of step 04)
Output: the per-year resistance-frequency table under results/
        the per-country resistance-frequency table under results/
        the W2 summary Markdown under results/
Site naming: literature numbering A2063G/A2064G/C2617 (slice position = literature position + 64, see the calibration report).
Note: the WT base at 2063 and 2064 is A in both cases; the site column of the typing TSV is already split per position, so the summary is computed per column.
"""
import csv, os, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = os.path.join(BASE, "results", "MY9_23S分型.tsv")
RESDIR = os.path.join(BASE, "results")

rows = []
with open(IN, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)

def bucket(year):
    y = int(year)
    return f"{(y // 5) * 5}s"

def aggregate(rows):
    n = len(rows)
    called = [r for r in rows if r["site2063"] != "-"]
    c = collections.Counter()
    c["n_total"] = n
    c["n_called"] = len(called)
    c["n_not_called"] = n - len(called)
    c["n_2063_mut"] = sum(1 for r in called if r["site2063"] != "A")
    c["n_2064_mut"] = sum(1 for r in called if r["site2064"] != "A")
    c["n_2617_mut"] = sum(1 for r in called if r["site2617"] != "C")
    c["n_2063G"] = sum(1 for r in called if r["site2063"] == "G")
    c["n_2063T"] = sum(1 for r in called if r["site2063"] == "T")
    c["n_2063C"] = sum(1 for r in called if r["site2063"] == "C")
    c["n_2064G"] = sum(1 for r in called if r["site2064"] == "G")
    c["n_2064C"] = sum(1 for r in called if r["site2064"] == "C")
    c["n_2617A"] = sum(1 for r in called if r["site2617"] == "A")
    c["n_2617G"] = sum(1 for r in called if r["site2617"] == "G")
    c["n_resist"] = sum(1 for r in called if r["site2063"] != "A" or r["site2064"] != "A" or r["site2617"] != "C")
    c["freq_resist"] = c["n_resist"] / c["n_called"] if c["n_called"] else 0.0
    c["n_2063only"] = sum(1 for r in called
                          if r["site2063"] != "A" and r["site2064"] == "A" and r["site2617"] == "C")
    c["n_2063_2064"] = sum(1 for r in called
                           if r["site2063"] != "A" and r["site2064"] != "A")
    c["n_2063_2617"] = sum(1 for r in called
                           if r["site2063"] != "A" and r["site2617"] != "C")
    c["n_2064only"] = sum(1 for r in called
                          if r["site2063"] == "A" and r["site2064"] != "A" and r["site2617"] == "C")
    c["n_2617only"] = sum(1 for r in called
                          if r["site2063"] == "A" and r["site2064"] == "A" and r["site2617"] != "C")
    return c

def write_table(path, groups, rows_groups, title):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group", "n_total", "n_called", "n_not_called", "n_resist", "freq_resist",
                    "n_2063G", "n_2063T", "n_2063C", "n_2064G", "n_2064C", "n_2617A", "n_2617G",
                    "n_2063only", "n_2063_2064", "n_2063_2617", "n_2064only", "n_2617only"])
        for g in groups:
            c = aggregate(rows_groups[g])
            w.writerow([g, c["n_total"], c["n_called"], c["n_not_called"], c["n_resist"],
                        f"{c['freq_resist']:.3f}",
                        c["n_2063G"], c["n_2063T"], c["n_2063C"], c["n_2064G"], c["n_2064C"], c["n_2617A"], c["n_2617G"],
                        c["n_2063only"], c["n_2063_2064"], c["n_2063_2617"], c["n_2064only"], c["n_2617only"]])
    print("written:", path)

# --- by year ---
by_year = collections.defaultdict(list)
for r in rows:
    y = r["year"].strip()
    if y.isdigit() and 1950 <= int(y) <= 2026:
        by_year[y].append(r)
years = sorted(by_year)
write_table(os.path.join(RESDIR, "MY9_耐药频率_按年份.tsv"), years, by_year, "by year")

# --- by country ---
by_ctry = collections.defaultdict(list)
for r in rows:
    cc = r["country"].split(":")[0].strip()
    if cc and cc != "not applicable":
        by_ctry[cc].append(r)
ctries = sorted(by_ctry, key=lambda k: -len(by_ctry[k]))
write_table(os.path.join(RESDIR, "MY9_耐药频率_按国家.tsv"), ctries, by_ctry, "by country")

# --- summary markdown ---
tot = aggregate(rows)
top = [(y, aggregate(by_year[y])) for y in years if len(by_year[y]) >= 10]
tct = [(k, aggregate(by_ctry[k])) for k in ctries if len(by_ctry[k]) >= 10]
md = [f"# MY9 W2 摘要：23S 耐药突变分型（{tot['n_total']} 株全分型）", "",
      f"- 分型成功 {tot['n_called']}/{tot['n_total']}（anchor_not_found={tot['n_not_called']}，锚点 v2 ±2nt 间隙设计）",
      f"- 耐药突变株（任一 site 非 WT）: **{tot['n_resist']}**（freq={tot['freq_resist']:.1%}）",
      f"- 位点分布: A2063G={tot['n_2063G']}, A2063T={tot['n_2063T']}, A2063C={tot['n_2063C']}, "
      f"A2064G={tot['n_2064G']}, C2617A={tot['n_2617A']}, C2617G={tot['n_2617G']}",
      f"- 组合: 仅2063={tot['n_2063only']}, 2063+2064={tot['n_2063_2064']}, 2063+2617={tot['n_2063_2617']}, "
      f"仅2064={tot['n_2064only']}, 仅2617={tot['n_2617only']}", "",
      "## 年份趋势（n≥10）", "", "| 年份 | n | 耐药率 | A2063G |", "|---|---|---|---|"]
for y, c in top:
    md.append(f"| {y} | {c['n_called']} | {c['freq_resist']:.1%} | {c['n_2063G']} |")
md += ["", "## 国家分布（n≥10）", "", "| 国家/地区 | n | 耐药率 | A2063G |", "|---|---|---|---|"]
for k, c in tct:
    md.append(f"| {k} | {c['n_called']} | {c['freq_resist']:.1%} | {c['n_2063G']} |")
md += ["", "> 生成: MY9_06_summary_tables.py；输入 results/MY9_23S分型.tsv；位点命名见 results/MY9_位点校准报告.md。"]
with open(os.path.join(RESDIR, "MY9_W2_摘要.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(md) + "\n")

print(f"total: called {tot['n_called']}/{tot['n_total']}, resist {tot['n_resist']} ({tot['freq_resist']:.1%})")
print(f"sites: 2063G={tot['n_2063G']} 2063T={tot['n_2063T']} 2063C={tot['n_2063C']} "
      f"2064G={tot['n_2064G']} 2617A={tot['n_2617A']} 2617G={tot['n_2617G']}")
print("written: results/MY9_W2_摘要.md")
