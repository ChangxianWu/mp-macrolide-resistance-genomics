# -*- coding: utf-8 -*-
"""29_fig3_redesign.py — Figure 3 重设计 v2（真实枝长系统发生 + 耐药变化事件标注）

旧版三个硬错误（本轮修复）:
  ① 面板 A xlim 写死 (-0.2, 10.5) 而 NJ 树深 56 层 → ST2/ST14/ST17 楔形全部画在画布外,
     与图注承诺的"折叠楔形"不符;
  ② 绿色楔形按"最小包含 ≥90% 克隆成员的子树"标注为 "Epidemic clone n = 626, 456/626 resistant"
     —— 626 是该子树 tip 数而非克隆规模（352）, 与正文/Table 4 直接矛盾;
  ③ 面板 B 标题与子图标题同 y 重叠, 总标题写死 "528 genomes"。

新版设计（对应当前叙事"根自由最小状态变化 + 面板特异性相似组"; 刻度组只表示 25-SNP
  single-linkage 成员, 不表示单系或传播链）:
  面板 A  全树 phylogram（x = SNP 距离, 真实枝长 ×742,308）: tip 按 23S 基因型着色;
          25-SNP 相似组 352 株在右缘以绿色刻度精确标出（不并成整块, 因该组在 NJ 拓扑上非严格单系）;
          A2063 变化边以符号标出（▲ 变为 G / ▽ 变回 A; 方向为展示拓扑下的最优重建, 非唯一历史）; ST 主要支系标注; SNP 标尺。
  面板 B  ST3 子树 phylogram（396 tips）: 6 次变化逐边标出; 主导耐药支系（≥95% G）涂色。
克隆组成仪表板已下沉为 Supplementary Figure S2（MY9_31_figS2_composition.py）。

单源纪律: tip 基因型/克隆成员读 results/MY9_分型矩阵.tsv 与 results/coresnp/MY9_紧克隆成员.tsv;
         变化次数与 Table 6 断言一致（34/22 全树, 6 子树）——不一致直接报错, 不静默出图。
输出: manuscript/figures/MY9_Fig3_coreSNP_tree.png / .pdf
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

# ---------- 单源数据 ----------
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


# ---------- Fitch（与 MY9_26 / MY9_20 同算法; 用于定位变化边, 计数与 Table 6 断言） ----------
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


# ---------- 布局: phylogram（x = 累积 SNP, y = tip 序） ----------
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
    """最大"纯度达标"支系 → 返回 (node, n_tips, purity)。用于 ST 支系标注（含纯度, 避免暗示单系）。"""
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
# 面板 A: 全树
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

# 克隆成员: 右缘绿色刻度（连续者自然成条; 不并成整块, 因克隆在 NJ 拓扑上非严格单系）
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

# A2063 变化边: ▲ 获得 G, ▽ 回复（方向相对当前显示根; 最小变化次数本身与定根无关）
for par, ch, a, b in edges_full:
    x, y = xsA[id(ch)] * CORE_SITES, ysA[id(ch)]
    if b == "G":
        axA.plot(x, y, marker="^", ms=3.4, color=C_RES, mec="white", mew=0.4, zorder=5)
    elif a == "G":
        axA.plot(x, y, marker="v", ms=3.4, color="#7f8c8d", mec="white", mew=0.4, zorder=5)

# 主要 ST 支系标注（带纯度; ST 非严格单系, 故不写成"ST3 clade"）
st_states = {t.name: (mat.get(t.name, {}).get("MLST_ST") or "?") for t in tipsA}
for lab in ("3", "14", "17", "2"):
    hit = largest_pure_clade(full, st_states, lab)
    if hit:
        nd, n, pur = hit
        axA.text(xsA[id(nd)] * CORE_SITES, ysA[id(nd)] + n_tipsA * 0.010,
                 f"ST{lab}-rich clade\n{n} tips, {pur * 100:.0f}% ST{lab}",
                 fontsize=7.4, color="#444", style="italic", ha="left", va="bottom", zorder=7,
                 bbox=dict(facecolor="white", alpha=0.78, edgecolor="none", pad=1.0))

# SNP 标尺
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
# 面板 B: ST3 子树
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

# 克隆成员: 与面板 A 同一编码（右缘绿色刻度）——单连锁聚类可串联, 成员在树上不集中于单一支系,
# 故不做"主导支系"着色（会误导为单系）, 改用逐株刻度 + 计数标注
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
