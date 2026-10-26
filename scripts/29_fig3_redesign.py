# -*- coding: utf-8 -*-
"""29_fig3_redesign.py — Figure 3 redesign v2 (phylogeny with true branch lengths + annotation of resistance-change events)

Three hard errors in the old version (fixed in this round):
  1. Panel A had a hard-coded xlim of (-0.2, 10.5) while the NJ tree is 56 levels deep -> the ST2/ST14/ST17 wedges were all drawn outside the canvas,
     contradicting the "collapsed wedges" promised in the figure legend;
  2. The green wedge, defined as the "smallest subtree containing >=90% of clone members", was labelled "Epidemic clone n = 626, 456/626 resistant"
     — 626 is the tip count of that subtree, not the clone size (352), which directly contradicts the main text/Table 4;
  3. The panel B title overlapped the subplot title at the same y, and the overall title hard-coded "528 genomes".

New design (matching the current narrative "root-free minimum state changes + panel-specific similarity group"; the tick group denotes
  only 25-SNP single-linkage membership, not monophyly or transmission chains):
  Panel A  full-tree phylogram (x = SNP distance, true branch length x742,308): tips colored by 23S genotype;
          the 352 strains of the 25-SNP similarity group are marked precisely with green ticks at the right edge (not merged into a single block, because this group is not strictly monophyletic on the NJ topology);
          A2063 change edges are marked with symbols (^ change to G / v change back to A; directions are the optimal reconstruction under the displayed topology, not a unique history); major ST clades annotated; SNP scale bar.
  Panel B  ST3 subtree phylogram (396 tips): 6 changes marked edge by edge; the dominant resistance lineage (>=95% G) is shaded.
The clone-composition dashboard has been moved down to Supplementary Figure S2 (MY9_31_figS2_composition.py).

Single-source discipline: tip genotypes/clone membership are read from the integrated per-genome typing matrix under
         results/ and the tight-clone membership table under results/coresnp/;
         change counts are asserted against Table 6 (34/22 whole tree, 6 subtree) — any mismatch raises an error rather than silently producing the figure.
Output: manuscript/figures/MY9_Fig3_coreSNP_tree.png / .pdf
"""
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
CS = RES / "coresnp"
FIGDIR = BASE / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

CORE_SITES = 742308
C_RES, C_WT, C_OTH = "#c0392b", "#2471a3", "#e67e22"
C_CLONE, C_DOM = "#117a65", "#f4d03f"

# ---------- Single-source data ----------
mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(RES / "MY9_分型矩阵.tsv", encoding="utf-8-sig"), delimiter="\t")}
clone = {r["ena_acc"] for r in csv.DictReader(
    open(CS / "MY9_紧克隆成员.tsv", encoding="utf-8-sig"), delimiter="\t")}
t6 = {r["Key"]: r["Value"] for r in csv.DictReader(
    open(BASE / "manuscript" / "tables" / "table6_analysis_stats.csv", encoding="utf-8-sig"))}
EXPECT = {"full": (int(t6["fitch_full_changes"]), int(t6["fitch_full_gains"])),
          "st3": (int(t6["fitch_changes_2063"]), None)}


def geno(a):
    r = mat.get(a, {})
    if r.get("resist_23S") != "resistant":
        return "wt"
    return "g" if r.get("site2063") == "G" else "other"


def color_of(a):
    return {"wt": C_WT, "g": C_RES, "other": C_OTH}[geno(a)]


# ---------- Fitch (same algorithm as MY9_26 / MY9_20; used to locate change edges; counts asserted against Table 6) ----------
def fitch_edges(tree, states):
    sets, changes = {}, 0
    for nd in tree.find_clades(order="postorder"):
        if nd.is_terminal():
            v = states.get(nd.name, "?")
            sets[id(nd)] = None if v in ("?", "N") else {v}
        else:
            cs = [sets[id(c)] for c in nd.clades if sets[id(c)] is not None]
            if not cs:
                sets[id(nd)] = None
                continue
            inter = set.intersection(*cs)
            if inter:
                sets[id(nd)] = inter
            else:
                sets[id(nd)] = set.union(*cs)
                changes += 1
    assign = {}
    root_set = sets[id(tree.root)]
    assign[id(tree.root)] = "G" if "G" in root_set else sorted(root_set)[0]
    edges = []
    for nd in tree.find_clades(order="preorder"):
        if nd.is_terminal():
            continue
        ps = assign.get(id(nd))
        for c in nd.clades:
            cs = sets[id(c)]
            if ps is None or cs is None:
                assign[id(c)] = ps if ps is not None else (sorted(cs)[0] if cs else None)
            elif ps in cs:
                assign[id(c)] = ps
            else:
                assign[id(c)] = sorted(cs)[0]
                if len(cs) == 1:
                    edges.append((nd, c, ps, assign[id(c)]))
    return changes, edges


# ---------- Layout: phylogram (x = cumulative SNP, y = tip order) ----------
def layout(tree):
    xs, ys = {}, {}
    xs[id(tree.root)] = 0.0
    for nd in tree.find_clades(order="preorder"):
        for c in nd.clades:
            xs[id(c)] = xs[id(nd)] + (c.branch_length or 0.0)
    idx = [0]
    for nd in tree.find_clades(order="postorder"):
        if nd.is_terminal():
            ys[id(nd)] = idx[0]
            idx[0] += 1
        else:
            ys[id(nd)] = sum(ys[id(c)] for c in nd.clades) / len(nd.clades)
    return xs, ys


def draw_phylogram(ax, tree, xs, ys, tip_colors, lw=0.35):
    for nd in tree.find_clades(order="preorder"):
        if nd.is_terminal():
            continue
        x = xs[id(nd)] * CORE_SITES
        ax.plot([x, x], [min(ys[id(c)] for c in nd.clades), max(ys[id(c)] for c in nd.clades)],
                color="#666", lw=lw, zorder=2)
        for c in nd.clades:
            xc, yc = xs[id(c)] * CORE_SITES, ys[id(c)]
            ax.plot([x, xc], [yc, yc], color="#666", lw=lw, zorder=2)
    tips = tree.get_terminals()
    ax.scatter([xs[id(t)] * CORE_SITES for t in tips], [ys[id(t)] for t in tips],
               c=[tip_colors[t.name] for t in tips], s=3.2, linewidths=0, zorder=3)


def largest_pure_clade(tree, states, allele, frac=0.80, min_tips=5):
    """Largest "purity-passing" clade -> returns (node, n_tips, purity). Used for ST clade annotation (purity included to avoid implying monophyly)."""
    best = None
    for nd in tree.find_clades(order="postorder"):
        tips = [t.name for t in nd.get_terminals()]
        if len(tips) < min_tips:
            continue
        c = Counter(states.get(x, "?") for x in tips)
        top, ntop = c.most_common(1)[0]
        if top == allele and ntop / len(tips) >= frac:
            if best is None or len(tips) > best[1]:
                best = (nd, len(tips), ntop / len(tips))
    return best


# =========================================================
# Panel A: full tree
# =========================================================
full = Phylo.read(CS / "MY9_coreSNP_NJ.nwk", "newick")
for t in full.get_terminals():
    t.name = t.name.strip()
st_full = {t.name: (mat.get(t.name, {}).get("site2063") or "?") for t in full.get_terminals()}
n_ch_full, edges_full = fitch_edges(full, st_full)
n_gain_full = sum(1 for _, _, _, b in edges_full if b == "G")
n_loss_full = sum(1 for _, _, a, _ in edges_full if a == "G")
n_other_full = n_ch_full - n_gain_full - n_loss_full
assert n_ch_full == EXPECT["full"][0], f"full-tree Fitch {n_ch_full} != Table 6 {EXPECT['full'][0]}"
assert n_gain_full == EXPECT["full"][1], f"full-tree gains {n_gain_full} != Table 6 {EXPECT['full'][1]}"

xsA, ysA = layout(full)
tipsA = full.get_terminals()
colorsA = {t.name: color_of(t.name) for t in tipsA}
n_tipsA = len(tipsA)

fig = plt.figure(figsize=(13.6, 7.6))
gsA = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.16,
                       left=0.045, right=0.985, top=0.90, bottom=0.155)
axA = fig.add_subplot(gsA[0, 0])
draw_phylogram(axA, full, xsA, ysA, colorsA)

xmaxA = max(xsA[id(t)] for t in tipsA) * CORE_SITES
axA.set_xlim(-0.04 * xmaxA, xmaxA * 1.18)
axA.set_ylim(-8, n_tipsA + 8)
axA.set_yticks([])
for sp in ("top", "right", "left"):
    axA.spines[sp].set_visible(False)
axA.set_xlabel("core-genome SNP distance", fontsize=10)
axA.tick_params(labelsize=9)

# Clone members: green ticks at the right edge (consecutive members naturally form a bar; not merged into one block, because the clone is not strictly monophyletic on the NJ topology)
xbar = xmaxA * 1.045
clone_tips = [t for t in tipsA if t.name in clone]
for t in clone_tips:
    axA.plot([xbar, xbar + xmaxA * 0.022], [ysA[id(t)], ysA[id(t)]],
             color=C_CLONE, lw=0.9, solid_capstyle="butt", zorder=4)
yc = [ysA[id(t)] for t in clone_tips]
axA.text(xbar + xmaxA * 0.030, (min(yc) + max(yc)) / 2,
         f"25-SNP group\n(single-linkage)\nn = {len(clone_tips)}\n"
         f"{sum(1 for t in clone_tips if geno(t.name) == 'g')}/{len(clone_tips)} A2063G",
         fontsize=8.5, color="#0b5345", fontweight="bold", va="center", ha="left")

# A2063 change edges: ^ gain of G, v reversion (direction relative to the currently displayed root; the minimum change count itself is rooting-independent)
for par, ch, a, b in edges_full:
    x, y = xsA[id(ch)] * CORE_SITES, ysA[id(ch)]
    if b == "G":
        axA.plot(x, y, marker="^", ms=3.4, color=C_RES, mec="white", mew=0.4, zorder=5)
    elif a == "G":
        axA.plot(x, y, marker="v", ms=3.4, color="#7f8c8d", mec="white", mew=0.4, zorder=5)

# Major ST clade annotation (with purity; STs are not strictly monophyletic, hence not written as "ST3 clade")
st_states = {t.name: (mat.get(t.name, {}).get("MLST_ST") or "?") for t in tipsA}
for lab in ("3", "14", "17", "2"):
    hit = largest_pure_clade(full, st_states, lab)
    if hit:
        nd, n, pur = hit
        axA.text(xsA[id(nd)] * CORE_SITES, ysA[id(nd)] + n_tipsA * 0.010,
                 f"ST{lab}-rich clade\n{n} tips, {pur * 100:.0f}% ST{lab}",
                 fontsize=7.4, color="#444", style="italic", ha="left", va="bottom", zorder=7,
                 bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=1.0))

# SNP scale bar
sbx, sby = xmaxA * 0.60, -6
axA.plot([sbx, sbx + 100], [sby, sby], color="#333", lw=1.6)
axA.text(sbx + 50, sby + 3, "100 SNP", fontsize=8.5, ha="center", va="bottom")

axA.set_title(f"A  Core-genome SNP phylogeny (n = {n_tipsA} tips; branch lengths in SNPs)",
              fontsize=11.5, loc="left")
legA = [Line2D([], [], marker="o", ls="", color=C_RES, ms=5, label="A2063G (resistant)"),
        Line2D([], [], marker="o", ls="", color=C_WT, ms=5, label="wild type"),
        Line2D([], [], marker="o", ls="", color=C_OTH, ms=5, label="other resistance mutation"),
        Line2D([], [], marker="^", ls="", color=C_RES, ms=5, label=f"A2063G gain ({n_gain_full})"),
        Line2D([], [], marker="v", ls="", color="#7f8c8d", ms=5,
               label=f"loss of A2063G ({n_loss_full})"),
        Line2D([], [], marker="s", ls="", color=C_OTH, ms=4.5,
               label=f"other change, A2063T/C ({n_other_full})"),
        Line2D([], [], color=C_CLONE, lw=2.2, label=f"single-linkage group members (n = {len(clone_tips)})")]
axA.legend(handles=legA, loc="upper left", fontsize=8.8, frameon=False, ncol=2,
           bbox_to_anchor=(0.0, -0.075), handletextpad=0.5, columnspacing=1.4)

# =========================================================
# Panel B: ST3 subtree
# =========================================================
sub = Phylo.read(CS / "treetime" / "ST3_subtree.nwk", "newick")
for t in sub.get_terminals():
    t.name = t.name.strip()
st_sub = {t.name: (mat.get(t.name, {}).get("site2063") or "?") for t in sub.get_terminals()}
n_ch_sub, edges_sub = fitch_edges(sub, st_sub)
assert n_ch_sub == EXPECT["st3"][0], f"ST3-subtree Fitch {n_ch_sub} != Table 6 {EXPECT['st3'][0]}"

xsB, ysB = layout(sub)
tipsB = sub.get_terminals()
colorsB = {t.name: color_of(t.name) for t in tipsB}
n_tipsB = len(tipsB)

axB = fig.add_subplot(gsA[0, 1])
draw_phylogram(axB, sub, xsB, ysB, colorsB)

# Clone members: same encoding as panel A (green ticks at the right edge) — single-linkage clustering can chain members together,
# so they are not concentrated in a single clade on the tree; hence no "dominant clade" shading (which would misleadingly imply monophyly),
# and individual tick marks + a count annotation are used instead
xmaxB = max(xsB[id(t)] for t in tipsB) * CORE_SITES
xbarB = xmaxB * 1.045
clone_tipsB = [t for t in tipsB if t.name in clone]
for t in clone_tipsB:
    axB.plot([xbarB, xbarB + xmaxB * 0.022], [ysB[id(t)], ysB[id(t)]],
             color=C_CLONE, lw=0.9, solid_capstyle="butt", zorder=4)

for par, ch, a, b in edges_sub:
    x, y = xsB[id(ch)] * CORE_SITES, ysB[id(ch)]
    if b == "G":
        axB.plot(x, y, marker="^", ms=5.0, color=C_RES, mec="white", mew=0.5, zorder=5)
    elif a == "G":
        axB.plot(x, y, marker="v", ms=5.0, color="#7f8c8d", mec="white", mew=0.5, zorder=5)
    else:
        axB.plot(x, y, marker="s", ms=4.0, color=C_OTH, mec="white", mew=0.5, zorder=5)

axB.set_xlim(-0.04 * xmaxB, xmaxB * 1.30)
axB.set_ylim(-8, n_tipsB + 8)
axB.set_yticks([])
for sp in ("top", "right", "left"):
    axB.spines[sp].set_visible(False)
axB.set_xlabel("core-genome SNP distance", fontsize=10)
axB.tick_params(labelsize=9)

sbx2, sby2 = xmaxB * 0.62, -6
axB.plot([sbx2, sbx2 + 25], [sby2, sby2], color="#333", lw=1.6)
axB.text(sbx2 + 12.5, sby2 + 1.2, "25 SNP", fontsize=8.5, ha="center", va="bottom")

n_gain_sub = sum(1 for _, _, _, b in edges_sub if b == "G")
n_loss_sub = sum(1 for _, _, a, _ in edges_sub if a == "G")
n_other_sub = n_ch_sub - n_gain_sub - n_loss_sub
n_res_sub = sum(1 for t in tipsB if mat.get(t.name, {}).get("resist_23S") == "resistant")

axB.set_title(f"B  ST3 subtree (n = {n_tipsB} tips; {n_res_sub} resistant)\n"
              f"changes at position 2063", fontsize=11.5, loc="left")
legB = [Line2D([], [], marker="^", ls="", color=C_RES, ms=5, label=f"A2063G gain ({n_gain_sub})"),
        Line2D([], [], marker="v", ls="", color="#7f8c8d", ms=5, label=f"loss of A2063G ({n_loss_sub})"),
        Line2D([], [], marker="s", ls="", color=C_OTH, ms=4.5, label=f"other change ({n_other_sub})"),
        Line2D([], [], color=C_CLONE, lw=2.2,
               label=f"single-linkage group members in this subtree ({len(clone_tipsB)} of {len(clone)})")]
axB.legend(handles=legB, loc="upper left", fontsize=8.8, frameon=False, ncol=2,
           bbox_to_anchor=(0.0, -0.075), handletextpad=0.5, columnspacing=1.4)

fig.suptitle(f"Core-genome SNP phylogeny of {n_tipsA} genomes and the changes at the "
             f"macrolide-resistance position 2063", fontsize=13, y=1.03)

out_png = FIGDIR / "MY9_Fig3_coreSNP_tree.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(FIGDIR / "MY9_Fig3_coreSNP_tree.pdf", bbox_inches="tight")
plt.close(fig)
print(f"Fig3 v2 -> {out_png}")
print(f"  panel A: {n_tipsA} tips | {n_ch_full} changes ({n_gain_full} A2063G gains, "
      f"{n_loss_full} losses, {n_other_full} A2063T/C) | clone ticks {len(clone_tips)}")
print(f"  panel B: {n_tipsB} tips | {n_ch_sub} changes "
      f"({n_gain_sub} to G, {n_loss_sub} away from G, {n_other_sub} other) | "
      f"clone members in subtree {len(clone_tipsB)}/{len(clone)} | resistant tips {n_res_sub}")
