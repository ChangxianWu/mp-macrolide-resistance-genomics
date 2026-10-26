# -*- coding: utf-8 -*-
"""30_figS1_threshold.py — 补充图 Figure S1: 阈值敏感性（数据单源 Table 4）
面板 A: 各阈值下最大簇规模 + 簇内耐药比例（左轴 n, 右轴 %）
面板 B: 最大簇规模随阈值的变化（簇生长曲线, 含 25→100 增量标注）
数据: manuscript/tables/table4_clone.csv（MY9_24 生成, 单源; 禁手抄数字）
输出: manuscript/figures/MY9_FigS1_threshold_sensitivity.png/pdf (300 dpi)
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T4 = os.path.join(BASE, "manuscript", "tables", "table4_clone.csv")
FIGDIR = os.path.join(BASE, "manuscript", "figures")
os.makedirs(FIGDIR, exist_ok=True)

rows = list(csv.DictReader(open(T4, encoding="utf-8-sig")))
rows.sort(key=lambda r: int(r["SNP_threshold"]))
thr = [int(r["SNP_threshold"]) for r in rows]
n = [int(r["N_cluster"]) for r in rows]
pct_res = [float(r["Pct_resistant"]) for r in rows]
pct_panel = [float(r["Pct_of_panel"]) for r in rows]
med = [int(r["Median_within_SNP"]) for r in rows]

C_N, C_P = "#117a65", "#c0392b"
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.6))

# ---- 面板 A: 最大簇规模 + 耐药比例 ----
axA.plot(thr, n, "o-", color=C_N, lw=2, ms=7, label="group size (left axis)")
for x, y in zip(thr, n):
    axA.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 8),
                 ha="center", fontsize=10, color=C_N, fontweight="bold")
axA.set_xlabel("Single-linkage SNP threshold", fontsize=11)
axA.set_ylabel("Genomes in largest cluster", fontsize=11, color=C_N)
axA.tick_params(axis="y", labelcolor=C_N)
axA.set_ylim(300, 400)
ax2 = axA.twinx()
ax2.plot(thr, pct_res, "s--", color=C_P, lw=1.8, ms=6, label="resistance-associated genotype within group (right axis)")
for x, y in zip(thr, pct_res):
    ax2.annotate(f"{y:.1f}%", (x, y), textcoords="offset points", xytext=(0, 10),
                 ha="center", fontsize=9.5, color=C_P)
ax2.set_ylabel("Resistance-associated genotype within group (%)", fontsize=11, color=C_P)
ax2.tick_params(axis="y", labelcolor=C_P)
ax2.set_ylim(90, 102)
axA.set_xticks(thr)
axA.spines[["top"]].set_visible(False)
axA.set_title("A  Group size vs genotype proportion", fontsize=12, loc="left")
h1, l1 = axA.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
axA.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=9.5, frameon=False)

# ---- 面板 B: 簇生长与阈值偏离 ----
base25 = n[0]
delta = [y - base25 for y in n]
axB.bar([str(t) for t in thr], n, color=["#117a65"] + ["#7fb3a5"] * (len(n) - 1), alpha=0.9)
for i, (t, y, d) in enumerate(zip(thr, n, delta)):
    axB.text(i, y + 1.5, f"+{d}" if d else "primary", ha="center", fontsize=10,
             color="#0b5345" if d == 0 else "#555", fontweight="bold" if d == 0 else "normal")
axB.set_xlabel("Single-linkage SNP threshold", fontsize=11)
axB.set_ylabel("Genomes in largest cluster", fontsize=11)
axB.set_ylim(300, 445)
axB.spines[["top", "right"]].set_visible(False)
axB.set_title("B  Growth of the largest cluster", fontsize=12, loc="left")

fig.suptitle("Sensitivity of the genomic-similarity group definition to the single-linkage SNP threshold (Table 4)",
             fontsize=13.5, y=1.02)
fig.tight_layout()
out = os.path.join(FIGDIR, "MY9_FigS1_threshold_sensitivity.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR, "MY9_FigS1_threshold_sensitivity.pdf"), bbox_inches="tight")
plt.close(fig)
print("FigS1 ->", out)
print("thr:", thr, "n:", n, "pct_res:", pct_res, "median_within:", med, "pct_panel:", pct_panel)
