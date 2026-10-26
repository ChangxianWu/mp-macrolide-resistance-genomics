# -*- coding: utf-8 -*-
"""23_epidemic_clade.py — Epidemic tight-clone definition (single-linkage clustering + threshold sensitivity)
Input: the core-SNP pairwise distance matrix under results/coresnp/ (substitutions/sites x 742,308 = SNPs), and the integrated per-genome typing matrix under results/
Method: single-linkage (union-find) clustering at SNP threshold t; report the largest cluster's size/composition/resistance rate/median within-cluster distance;
      sensitivity across t = 25/50/75/100; recompute the root-free pairwise molecular clock for the largest cluster (with the full ST3 set as control).
Output: the epidemic tight-clone definition report (Markdown) under results/
"""
import csv, os
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
CS = os.path.join(RES, "coresnp")
CORE_SITES = 742308  # Authoritative core (>=95%) value for the v3 panel of 646 strains (measured in MY9_15 stage2 on 2026-09-22; do not revert)

mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
L = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
accs = [ln.split("\t")[0] for ln in L[1:]]
D = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in L[1:]])
D = (D + D.T) / 2 * CORE_SITES
n = len(accs)
st = np.array([mat.get(a, {}).get("MLST_ST", "") or "novel" for a in accs])
res = np.array([mat.get(a, {}).get("resist_23S", "") for a in accs])
yr = np.array([float(mat[a]["year"]) if mat.get(a, {}).get("year", "").isdigit() else np.nan
               for a in accs])
ctry = np.array([mat.get(a, {}).get("country", "").split(":")[0].strip() for a in accs])

def clusters(t):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    iu, ju = np.nonzero(np.triu(D <= t, 1))
    for a, b in zip(iu, ju):
        union(int(a), int(b))
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=len, reverse=True)

lines = ["# MY9 流行紧克隆定义（单连锁 SNP 阈值聚类）\n",
         "\n- 距离口径: core-SNP 成对距离（缺失对忽略后归一 × 742,308 = SNP）\n"]
result = {}
for t in (25, 50, 75, 100):
    cls = clusters(t)
    big = cls[0]
    idx = np.array(big)
    sub = D[np.ix_(idx, idx)]; iu = np.triu_indices(len(idx), 1)
    med = np.median(sub[iu])
    from collections import Counter
    stc = Counter(st[idx]); rc = res[idx]
    n_res = int((rc == "resistant").sum()); n_wt = int((rc == "WT").sum())
    yr_ok = yr[idx][~np.isnan(yr[idx])]
    cc = Counter(ctry[idx])
    lines.append(f"\n## 阈值 t = {t} SNP\n")
    lines.append(f"- 最大簇 n={len(idx)}（占 {len(idx)/n:.0%}）, 内部成对中位 {med:.0f} SNP\n")
    lines.append(f"- ST 构成: " + ", ".join(f"ST{k or 'novel'} {v}" for k, v in stc.most_common(5)) + "\n")
    lines.append(f"- 耐药 {n_res} / 敏感 {n_wt} / 未call {len(idx)-n_res-n_wt}"
                 f"（耐药率 {n_res/max(1,n_res+n_wt):.1%}）\n")
    lines.append(f"- 年份: 中位 {np.median(yr_ok):.0f}, 范围 {yr_ok.min():.0f}–{yr_ok.max():.0f}"
                 f"; ≥2019 占 {(yr_ok>=2019).mean():.0%}\n")
    lines.append(f"- 国家/地区 Top: " + ", ".join(f"{k} {v}" for k, v in cc.most_common(6)) + "\n")
    result[t] = idx
    print(f"t={t}: n={len(idx)} med={med:.0f} resist={n_res}/{n_res+n_wt}")

# Root-free pairwise molecular clock for the largest cluster (t=50) vs the full ST3 set
idx = result[50]
sub = D[np.ix_(idx, idx)]
y = yr[idx]
dy, dd = [], []
for i in range(len(idx)):
    for j in range(i + 1, len(idx)):
        if np.isnan(y[i]) or np.isnan(y[j]):
            continue
        dt = abs(y[i] - y[j])
        if dt > 0:
            dy.append(dt); dd.append(sub[i, j])
dy = np.array(dy); dd = np.array(dd)
slope, intercept = np.polyfit(dy, dd, 1)
rr = np.corrcoef(dy, dd)[0, 1]
i3 = np.nonzero(st == "3")[0]
sub3 = D[np.ix_(i3, i3)]; y3 = yr[i3]
dy3, dd3 = [], []
for i in range(len(i3)):
    for j in range(i + 1, len(i3)):
        if np.isnan(y3[i]) or np.isnan(y3[j]):
            continue
        dt = abs(y3[i] - y3[j])
        if dt > 0:
            dy3.append(dt); dd3.append(sub3[i, j])
dy3 = np.array(dy3); dd3 = np.array(dd3)
s3, b3 = np.polyfit(dy3, dd3, 1)
lines.append("\n## 免根成对分子钟对照\n")
lines.append(f"- 紧克隆（t=50 最大簇, n={len(idx)}）: 斜率 {slope:.2f} SNP/年（两谱系）, "
             f"r={rr:.3f}, 截距 {intercept:.1f} SNP, n={len(dy):,} 对\n")
lines.append(f"- 全 ST3（n={len(i3)}）: 斜率 {s3:.2f}, r={np.corrcoef(dy3, dd3)[0,1]:.3f}, "
             f"截距 {b3:.1f}, n={len(dy3):,} 对\n")
lines.append("\n## 解读\n")
lines.append("- ST3 序列型在 core-SNP 层面**非紧密单系**: 内部 p90=81、p99=242、最大 520 SNP,"
             "含法/德老株（1970–2011）与南京 2024 离群株; 与 ST17 中位距离仅 156 SNP（最小 118）\n")
lines.append("- 论文应以**紧克隆子群**（本表 t 阈值簇）作为流行单元陈述, 而非整段 ST3;"
             "ST3 作为 MLST 序列型仍适用于耐药率分层表, 但系统发生叙述需用紧克隆\n")

out = os.path.join(RES, "MY9_流行克隆定义.md")
open(out, "w", encoding="utf-8").write("".join(lines))
print("written", out)
print("".join(lines[:14]))

# Tight-clone (largest 25-SNP cluster) member list -> reused for the folded Fig3 presentation (algorithm unchanged)
_big = clusters(25)[0]
with open(os.path.join(os.path.dirname(out), "MY9_紧克隆成员.tsv"), "w", encoding="utf-8-sig", newline="") as _f:
    _f.write("ena_acc" + chr(9) + "MLST_ST" + chr(9) + "resist_23S" + chr(9) + "year" + chr(9) + "country" + chr(10))
    for _j in _big:
        _a = accs[_j]
        _mm = mat.get(_a, {})
        _f.write(_a + chr(9) + str(_mm.get("MLST_ST","")) + chr(9) + str(_mm.get("resist_23S","")) + chr(9) + str(_mm.get("year","")) + chr(9) + str(_mm.get("country","")) + chr(10))
print(f"clone members exported: {len(_big)}")
