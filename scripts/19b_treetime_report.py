# -*- coding: utf-8 -*-
"""19b_treetime_report.py — dating robustness report (grid summary + marginal reconstruction CIs)
Inputs: results/coresnp/treetime/{ST3_subtree.nwk, ST3.aln.fasta, ST3_dates.csv, tt_summary.tsv}
Method: read the grid results setting by setting; for the most robust setting (noprior_oldest) add a
      marginal=True reconstruction and take the 90% maximum posterior regions of the root / the
      A2063G resistance clade MRCA / the 2023+ wave MRCA.
Output: the TreeTime dating report (Markdown) under results/coresnp/
"""
import csv, os
import numpy as np
from treetime import TreeTime

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
TTDIR = os.path.join(CS, "treetime")
CORE_SITES = 742308
SITE2063_COL = None  # computed from the core mask

import glob
accs_all = [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(CS, "alleles", "*.npy")))]
M = np.stack([np.load(os.path.join(CS, "alleles", a + ".npy")) for a in accs_all])
core_mask = (M < 4).sum(axis=0) >= int(0.95 * M.shape[0])
ci_idx = np.nonzero(core_mask)[0]
SITE2063_COL = int(np.nonzero(ci_idx == 122118)[0][0])
print("A2063 col:", SITE2063_COL, flush=True)

rows = list(csv.DictReader(open(os.path.join(TTDIR, "tt_summary.tsv"), encoding="utf-8-sig"), delimiter="\t"))
with open(os.path.join(TTDIR, "ST3_dates.csv"), encoding="utf-8-sig") as f:
    dates = {r["name"]: float(r["date"]) for r in csv.DictReader(f)}

lines = ["# MY9 TreeTime 正式定年报告（ST3 克隆子树）\n",
         f"\n- 输入: `results/coresnp/treetime/ST3_subtree.nwk`（{len(dates)} 末梢, core-SNP NJ 拓扑）、"
         "`ST3.aln.fasta`（742,308 位点）、`ST3_dates.csv`（采样年）\n",
         "- 工具: phylo-treetime 0.11.4（本地 venv `MY9/.venv-align`, 2026-09-19 装）\n",
         "- 口径: clock_rate（subs/site/yr）× 742,308 = SNP/genome/yr（每谱系）; "
         "R² = 根-末梢回归解释方差\n",
         "\n## 一、定根 × 先验 网格（速率可识别性）\n",
         "\n| 设置 | 定根 | 共祖先先验 | 速率 (SNP/基因组/年) | R² | 根年代 (CE) |\n",
         "|---|---|---|---|---|---|\n"]
for r in rows:
    pri = "无" if r["Tc"].startswith("1e-06") else f"Tc={r['Tc']}"
    lines.append(f"| {r['setting']} | {r['root']} | {pri} | {r['rate_snp_genome_yr']} | "
                 f"{r['R2']} | {r['root_numdate']} |\n")

rates = [float(r["rate_snp_genome_yr"]) for r in rows if r["rate_snp_genome_yr"] != "err"]
lines.append(f"\n- **速率跨度 {min(rates):.2f} – {max(rates):.2f} SNP/基因组/年（>10 倍）; "
             f"R² {min(float(r['R2']) for r in rows):.3f} – {max(float(r['R2']) for r in rows):.3f}**\n")
_t5p = os.path.join(BASE, "manuscript", "tables", "table5_temporal.csv")
_ref_txt = "参照: 免根成对回归（不依赖定根/先验）"
try:
    _st3 = next(r for r in csv.DictReader(open(_t5p, encoding="utf-8-sig")) if r["Dataset"].startswith("All ST3"))
    _ref_txt += (f"**{_st3['Pairwise_rate_SNP_per_yr_both_lineages']} SNP/年（两谱系）= "
                 f"{_st3['Per_lineage_rate']} SNP/基因组/年每谱系**（ST3, Δyear>0, n={_st3['N_pairs']} 对, "
                 f"r={_st3['r']}）; 同年成对仍有 ~{_st3['Intercept_SNP']} SNP 差异（截距）")
except Exception:
    _ref_txt += "（见 Table 5）"
lines.append(f"- {_ref_txt}" + chr(10))

# --- marginal reconstruction (most robust setting only) ---
lines.append("\n## 二、边际重建 90% 后验区间（设置 noprior_oldest）\n")
try:
    tt = TreeTime(gtr="Jukes-Cantor", tree=os.path.join(TTDIR, "ST3_subtree.nwk"),
                  aln=os.path.join(TTDIR, "ST3.aln.fasta"), dates=dates, verbose=0)
    tt.run(root="oldest", infer_gtr=False, resolve_polytomies=True,
           coalescent=None, Tc=1e-6, max_iter=1, time_marginal="only-final")
    root = tt.tree.root

    def ci(node):
        try:
            lo, hi = tt.get_max_posterior_region(node, fraction=0.9)
            return f"{node.numdate:.0f} [{lo:.0f}, {hi:.0f}]"
        except Exception as e:
            return f"{node.numdate:.0f} (CI n/a: {type(e).__name__})"

    lines.append(f"- 速率 {float(tt.date2dist.clock_rate)*CORE_SITES:.2f} SNP/基因组/年, "
                 f"R²={float(getattr(tt.date2dist,'r_val',float('nan'))):.4f}\n")
    lines.append(f"- 树根: {ci(root)}\n")
    tips = root.get_terminals()
    res_tips = []
    for t in tips:
        s = getattr(t, "sequence", None)
        if s is not None and s[SITE2063_COL] == "G":
            res_tips.append(t)
    if len(res_tips) > 1:
        mc = root.common_ancestor(res_tips)
        lines.append(f"- A2063G 耐药支系 MRCA: **{ci(mc)}**（{len(res_tips)} 末梢）\n")
    y23 = [t for t in tips if dates.get(t.name, 0) >= 2023.0]
    if len(y23) > 1:
        m23 = root.common_ancestor(y23)
        lines.append(f"- 2023+ 波 MRCA: **{ci(m23)}**（{len(y23)} 末梢）\n")
    y19 = [t for t in tips if 2019.0 <= dates.get(t.name, 0) < 2021.0]
    if len(y19) > 1:
        m19 = root.common_ancestor(y19)
        lines.append(f"- 2019–2020 波 MRCA: **{ci(m19)}**（{len(y19)} 末梢）\n")
except Exception as e:
    lines.append(f"- 边际重建失败: {type(e).__name__}: {e}\n")
    print("marginal failed:", e, flush=True)

lines.append("\n## 三、结论（如实）\n")
lines.append("1. **株级时间信号不足**: 五种定根/先验组合给出的速率相差 >10 倍（含负值），"
             "根年代横跨 1191–1950 CE，R² 最高 0.359 —— 该采样设计（80% 末梢来自 2019/2023 两年）"
             "下 TreeTime 定年**不可识别**，绝对年代（含 MRCA 年）不作结论性报告。\n")
lines.append("2. **边际重建点估计仅供参考**: 耐药支系 MRCA 1996、2023+ 波 1998、2019–20 波 2010、"
             "根约 1970——区间近零宽（模型过度自信），且与免根成对回归外推（2023 波约 2018 年）相差约 20 年。"
             "**不可作为表型结论引用**。\n")
lines.append("3. **可稳健陈述的定性结论**: 各设置一致显示耐药 ST3 支系的共祖远早于 2019/2023 采样波"
             "（Δyear=0 成对仍有 ~20 SNP, 即扩增前已存在深度分歧的共循环亚支系）——"
             "两轮波是既有支系的扩张而非近期新生突变。\n")
lines.append("4. **报告口径**: 论文时间证据以免根成对回归为主（Table 5 ST3 行）; "
             "TreeTime 结果作为敏感性分析呈现并声明局限。\n")
lines.append("5. **原因**: 无重组过滤（Gubbins/ClonalFrameML 不可用）+ 同波共存亚支系 + 采样年份高度集中。\n")
lines.append("6. **后续**: 正式 MRCA 定年需 BEAST2 松弛钟 + 重组掩蔽 + 更密时间采样; "
             "本机工具链不可得, 列为限制项与后续工作。\n")

out = os.path.join(CS, "MY9_TreeTime定年报告.md")
open(out, "w", encoding="utf-8").write("".join(lines))
print("written", out, flush=True)
print("".join(lines))
