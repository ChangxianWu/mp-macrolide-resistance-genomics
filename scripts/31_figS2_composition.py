# -*- coding: utf-8 -*-
"""31_figS2_composition.py — Supplementary Figure S2: epidemic-clone composition

Rationale: when Figure 3 was redesigned, the clone-composition dashboard was demoted from the main
      figure to a supplementary figure — the main figure became phylogeny + change events;
      composition is descriptive statistics and therefore fits better in the supplement, and the
      original four small bar charts (one of which, P1 type, had only a single bar) had too low an information density.
Improvements: the year/country/MLST panels are all annotated with counts and percentages; P1 type and the 2063 genotype are drawn as two
      100% stacked bars (a single full-width segment is itself the visual conclusion that the composition is perfectly homogeneous), and the takeaway sentence below draws its values from Table 4;
      titles moved to panel level (to avoid overlapping text in the previous version).
Single source: the tight-clone membership table under results/coresnp/ (members) + the integrated per-genome
      typing matrix under results/ (strain-level attributes) +
      manuscript/tables/table4_clone.csv (size and resistance proportion).
Output: manuscript/figures/MY9_FigS2_clone_composition.png / .pdf
"""
import csv
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
CS = RES / "coresnp"
FIGDIR = BASE / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

C_NORM = {"Taiwan": "Taiwan, China", "China": "Mainland China",
          "USA": "United States", "not applicable": "Not specified"}
C_CLONE = "#117a65"


def cnorm(c):
    return C_NORM.get(c.split(":")[0].strip(), c.split(":")[0].strip())


def disp_st(k):
    """Normalize the displayed ST naming (review P2; same convention as MY9_24)."""
    if k.startswith("uncertain(near "):
        return k.replace("uncertain(near ", "").rstrip(")") + "-related"
    if k.startswith("novel(near "):
        return k.replace("novel(near ", "").rstrip(")") + "-novel"
    return k


mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(RES / "MY9_分型矩阵.tsv", encoding="utf-8-sig"), delimiter="\t")}
mem = list(csv.DictReader(open(CS / "MY9_紧克隆成员.tsv", encoding="utf-8-sig"), delimiter="\t"))
clone = [r["ena_acc"] for r in mem]
t4 = {r["SNP_threshold"]: r for r in csv.DictReader(
    open(BASE / "manuscript" / "tables" / "table4_clone.csv", encoding="utf-8-sig"))}
row25 = t4["25"]

years = Counter(int(mat[a]["year"]) for a in clone if mat.get(a, {}).get("year", "").isdigit())
ctry = Counter(cnorm(mat[a].get("country", "")) for a in clone)
st = Counter((mat[a].get("MLST_ST") or "uncertain") for a in clone)
p1 = Counter((mat[a].get("P1_type") or "uncertain") for a in clone)
geno = Counter(mat[a].get("site2063", "?") for a in clone)

fig = plt.figure(figsize=(11.5, 7.0))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.30,
                      left=0.075, right=0.975, top=0.86, bottom=0.075)

# (A) Collection year
ax = fig.add_subplot(gs[0, 0])
ks = sorted(years)
ax.bar([str(k) for k in ks], [years[k] for k in ks], color=C_CLONE, alpha=0.9)
for i, k in enumerate(ks):
    ax.text(i, years[k] + max(years.values()) * 0.02, str(years[k]),
            ha="center", va="bottom", fontsize=8.5, color="#333")
ax.set_ylabel("genomes", fontsize=10.5)
ax.set_title(f"(A) Collection year (median {row25['Median_year']})", fontsize=11.5, loc="left")
ax.tick_params(labelsize=9.5, rotation=45)
ax.set_ylim(0, max(years.values()) * 1.14)
ax.spines[["top", "right"]].set_visible(False)

# (B) Country/region
ax = fig.add_subplot(gs[0, 1])
top = ctry.most_common(6)
ys = list(range(len(top)))[::-1]
ax.barh(ys, [v for _, v in top], color="#2874a6", alpha=0.9)
ax.set_yticks(ys)
ax.set_yticklabels([k for k, _ in top], fontsize=8.5)
for y, (k, v) in zip(ys, top):
    ax.text(v + len(clone) * 0.01, y, f"{v} ({v / len(clone) * 100:.0f}%)",
            va="center", fontsize=8.5, color="#333")
ax.set_xlabel("genomes", fontsize=10.5)
ax.set_xlim(0, max(v for _, v in top) * 1.22)
ax.set_title("(B) Country or region", fontsize=11.5, loc="left")
ax.tick_params(labelsize=9.5)
ax.spines[["top", "right"]].set_visible(False)

# (C) MLST sequence type
ax = fig.add_subplot(gs[1, 0])
tops = st.most_common(4)
xs = range(len(tops))
ax.bar(list(xs), [v for _, v in tops], color="#8e44ad", alpha=0.9)
ax.set_xticks(list(xs))
ax.set_xticklabels([f"ST{disp_st(k)}" for k, _ in tops], fontsize=9.5)
for i, (k, v) in enumerate(tops):
    ax.text(i, v + len(clone) * 0.01, f"{v} ({v / len(clone) * 100:.0f}%)",
            ha="center", va="bottom", fontsize=8.5, color="#333")
ax.set_ylabel("genomes", fontsize=10.5)
ax.set_ylim(0, max(v for _, v in tops) * 1.16)
ax.set_title("(C) MLST sequence type", fontsize=11.5, loc="left")
ax.tick_params(labelsize=9.5)
ax.spines[["top", "right"]].set_visible(False)

# (D) P1 type + 2063 genotype: two 100% stacked bars (compositional homogeneity is itself the conclusion; the segment loop is generic,
#     so any future group with more than one segment gets automatic colors/external labels)
ax = fig.add_subplot(gs[1, 1])
groups = [
    ("P1 type", [(k, v) for k, v in p1.most_common()]),
    ("23S rRNA genotype (position 2063)", [
        ("A2063G" if k == "G" else ("wild type" if k == "A" else "A2063" + str(k)), v)
        for k, v in geno.most_common()]),
]
BASE_COL, ALT = {"P1 type": C_CLONE, "23S rRNA genotype (position 2063)": "#c0392b"}, \
    ["#117a65", "#7fb3a5", "#c0392b", "#e67e22"]
for yi, (gname, segs) in enumerate(groups):
    left = 0
    for si, (lab, v) in enumerate(segs):
        col = BASE_COL[gname] if len(segs) == 1 else ALT[si % len(ALT)]
        ax.barh(yi, v, left=left, color=col, alpha=0.9, height=0.52,
                edgecolor="white", linewidth=1.2)
        txt = f"{lab}  {v} ({v / len(clone) * 100:.0f}%)"
        if v >= len(clone) * 0.45:
            ax.text(left + v / 2, yi, txt, ha="center", va="center",
                    fontsize=10.5, color="white", fontweight="bold")
        else:
            ax.text(left + v / 2, yi + 0.36, txt, ha="center", va="bottom",
                    fontsize=9, color="#333")
        left += v
ax.set_yticks([0, 1])
ax.set_yticklabels([g for g, _ in groups], fontsize=10.5)
ax.set_xlim(0, len(clone) * 1.01)
ax.set_ylim(-1.75, 1.45)
ax.set_xlabel("genomes", fontsize=10.5)
ax.tick_params(labelsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("(D) Lineage markers within the single-linkage group", fontsize=11.5, loc="left")
takeaway = (f"All {row25['Resistant_called']} group members carry a resistance-associated "
            f"genotype ({row25['Pct_resistant']}% of called; Table 4); median within-group "
            f"pairwise distance {row25['Median_within_SNP']} SNP (Table 4).")
ax.text(0, -0.78, "\n".join(textwrap.wrap(takeaway, 58)),
        fontsize=9.5, color="#222", va="top", ha="left", linespacing=1.45)

fig.suptitle(f"Composition of the 25-SNP single-linkage group ({len(clone)} assemblies, "
             f"{row25['Pct_of_panel']}% of the panel)", fontsize=13.5, y=0.955)

out_png = FIGDIR / "MY9_FigS2_clone_composition.png"
fig.savefig(out_png, dpi=300)
fig.savefig(FIGDIR / "MY9_FigS2_clone_composition.pdf")
plt.close(fig)
print(f"FigS2 -> {out_png}")
print(f"  clone {len(clone)} | years {ks[0]}-{ks[-1]} | regions {len(ctry)} | "
      f"ST {dict(st.most_common(3))} | P1 {dict(p1.most_common(2))} | 2063 {dict(geno.most_common(3))}")
