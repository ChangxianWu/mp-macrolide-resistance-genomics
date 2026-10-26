# -*- coding: utf-8 -*-
"""34_fig5_temporal_series.py — Figure 5: Japanese time series + year-by-year composition of the 25-SNP group

Rationale: part 6 of the review asked for a visualization of "the deepest temporal series" and "the effect of sampling expansion on ancestral-state reconstruction".
      MY9_33 already quantified both: the 352-genome group remains a single connected component at within-group 10/15/20 SNP thresholds (no sub-cluster of >=5 genomes splits off);
      after pruning the 114 Japanese backfill genomes, the minimum number of changes drops from 34 to 24 (whole tree) and from 6 to 3 (ST3 subtree).
Design: (A) Japanese assemblies as stacked bars by collection year (wild type / resistance-associated), annotating the 1976-1994 all-wild-type window
      and the first resistance-associated genotypes in 2000; (B) the 25-SNP group as stacked bars by year (regional composition),
      annotating the three sampling waves and the absence of any split on within-group re-clustering.
Single source: the integrated per-genome typing matrix under results/ (genotype/year/country) + the tight-clone
      membership table under results/coresnp/ (group members).
Output: manuscript/figures/MY9_Fig5_temporal_series.png / .pdf
"""
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
FIGDIR = BASE / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

C_NORM = {"Taiwan": "Taiwan, China", "China": "Mainland China",
          "USA": "United States", "not applicable": "Not specified"}
C_WT, C_RES = "#95a5a6", "#c0392b"


def cnorm(c):
    return C_NORM.get(c.split(":")[0].strip(), c.split(":")[0].strip())


mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(RES / "MY9_分型矩阵.tsv", encoding="utf-8-sig"), delimiter="\t")}
mem = {r["ena_acc"] for r in csv.DictReader(
    open(RES / "coresnp" / "MY9_紧克隆成员.tsv", encoding="utf-8-sig"), delimiter="\t")}
yr = lambda r: int(r["year"]) if r["year"].isdigit() else None

jp = [r for r in mat.values() if cnorm(r["country"]) == "Japan"]
jp_wt = Counter(yr(r) for r in jp if r["resist_23S"] == "WT")
jp_res = Counter(yr(r) for r in jp if r["resist_23S"] == "resistant")
pre2000_wt = sum(v for y, v in jp_wt.items() if y < 2000)
first_res = min(y for y, v in jp_res.items() if v)
clone = [r for r in mat.values() if r["ena_acc"] in mem]
c_years = [yr(r) for r in clone if yr(r) is not None]
waves = {"2010-2018": sum(1 for y in c_years if y <= 2018),
         "2019-2020": sum(1 for y in c_years if 2019 <= y <= 2020),
         "2023-2025": sum(1 for y in c_years if y >= 2023)}
assert sum(waves.values()) == len(clone) == 352, "wave composition mismatch"

fig = plt.figure(figsize=(11.5, 5.1))
gs = fig.add_gridspec(1, 2, hspace=0.32, wspace=0.22,
                      left=0.055, right=0.985, top=0.82, bottom=0.27)

# ---- (A) Japanese time series ----
ax = fig.add_subplot(gs[0, 0])
years = list(range(1976, 2025))
wt_v = [jp_wt.get(y, 0) for y in years]
res_v = [jp_res.get(y, 0) for y in years]
ax.bar(years, wt_v, color=C_WT, alpha=0.95, width=0.8, label="Wild type at assayed loci")
ax.bar(years, res_v, bottom=wt_v, color=C_RES, alpha=0.95, width=0.8,
       label="Resistance-associated genotype")
ax.axvspan(1975.5, 1994.5, color="#1e8449", alpha=0.07)
ax.text(1985, 24.6, f"1976-1994: all wild type\nat assayed loci (n={pre2000_wt})",
        ha="center", fontsize=9, color="#1e8449")
ax.annotate(f"first resistance-associated\ngenotypes ({first_res}: C2617G + A2063C)",
            xy=(2000, 2.1), xytext=(1999.5, 15.5), fontsize=8.5, color="#c0392b",
            ha="center", arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.9))
for y, tot in zip(years, [a + b for a, b in zip(wt_v, res_v)]):
    if tot and y >= 2008:
        ax.text(y, tot + 0.5, str(tot), ha="center", va="bottom", fontsize=8, color="#333")
ax.set_ylabel("genomes", fontsize=10.5)
ax.set_xlim(1974.5, 2025.5)
ax.set_ylim(0, 33)
ax.set_xticks(list(range(1976, 2025, 4)))
ax.tick_params(labelsize=9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(f"(A) Japanese assemblies by collection year (n={len(jp)})",
             fontsize=11.5, loc="left")
ax.legend(fontsize=9, loc="upper left", frameon=False)

# ---- (B) 25-SNP group by year and region ----
ax = fig.add_subplot(gs[0, 1])
cyears = sorted(set(c_years))
regions = Counter(cnorm(r["country"]) for r in clone)
top_regions = [k for k, _ in regions.most_common(4)]
REGION_COL = {"Mainland China": "#117a65", "Taiwan, China": "#2874a6",
              "Japan": "#e67e22", "South Korea": "#8e44ad"}
per_year = {y: Counter() for y in cyears}
for r in clone:
    y = yr(r)
    k = cnorm(r["country"])
    per_year[y]["Other regions" if k not in top_regions else k] += 1
labels = top_regions + ["Other regions"]
cols = [REGION_COL.get(k, "#b3b6b7") for k in labels]
bottom = [0] * len(cyears)
for lab, col in zip(labels, cols):
    vals = [per_year[y][lab] for y in cyears]
    ax.bar(cyears, vals, bottom=bottom, color=col, alpha=0.92, width=0.8, label=lab)
    bottom = [a + b for a, b in zip(bottom, vals)]
for y, tot in zip(cyears, bottom):
    ax.text(y, tot + 4, str(tot), ha="center", va="bottom", fontsize=8.5, color="#333")
ax.set_ylabel("genomes", fontsize=10.5)
ax.set_xticks(cyears)
ax.set_xticklabels(cyears, fontsize=9, rotation=45)
ax.set_ylim(0, 250)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(f"(B) 25-SNP single-linkage group by collection year (n={len(clone)})",
             fontsize=11.5, loc="left")
ax.legend(fontsize=9, loc="upper left", frameon=False)
takeaway = (f"Sampled in three periods ({waves['2010-2018']} in 2010-2018, {waves['2019-2020']} in 2019-2020, "
            f"{waves['2023-2025']} in 2023-2025); one connected component when the within-group "
            f"single-linkage is repeated at 10, 15 and 20 SNPs.")
import textwrap
takeaway = (f"Sampled in three periods ({waves['2010-2018']} in 2010-2018, {waves['2019-2020']} in 2019-2020, "
            f"{waves['2023-2025']} in 2023-2025); one connected component when the within-group "
            f"single-linkage is repeated at 10, 15 and 20 SNPs.")
ax.text(0, -0.30, "\n".join(textwrap.wrap(takeaway, 62)), transform=ax.transAxes,
        fontsize=9, color="#222", va="top", ha="left")

fig.suptitle("Temporal structure: the Japanese series and the 25-SNP single-linkage group",
             fontsize=13.5, y=0.94)

out_png = FIGDIR / "MY9_Fig5_temporal_series.png"
fig.savefig(out_png, dpi=300)
fig.savefig(FIGDIR / "MY9_Fig5_temporal_series.pdf")
plt.close(fig)
print(f"Fig5 -> {out_png}")
print(f"  Japan {len(jp)} (WT {sum(jp_wt.values())} / resistant {sum(jp_res.values())}); "
      f"pre-2000 WT {pre2000_wt}; first resistant {first_res}; "
      f"clone years {min(c_years)}-{max(c_years)}; waves {waves}")
