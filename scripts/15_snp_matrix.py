# -*- coding: utf-8 -*-
"""15_snp_matrix.py — 本地 core-SNP stage 2：核心位点组装 + SNP 距离 + NJ + root-to-tip
输入: results/coresnp/alleles/*.npy（MY9_14 产物, 0=A 1=C 2=G 3=T 4=缺失/掩蔽）
      + results/MY9_分型矩阵.tsv（年份/ST/耐药注释） + results/MY9_mash距离.phylip（互验）
输出: results/coresnp/MY9_SNP距离.phylip
      results/coresnp/MY9_coreSNP_NJ.nwk（中点定根）
      results/coresnp/MY9_coreSNP_root_to_tip.tsv
      results/coresnp/MY9_coreSNP_摘要.md
方法:
  core = 非缺失株比例 ≥95% 的参考位点（indel±8nt/重复区/结构断点已在 stage1 掩蔽为 4）
  距离 = core 位点上等位不同的计数（成对缺失按不计入处理: 任一方为 4 的位点不比较）
  NJ/中点定根/root-to-tip 同 MY9_12; 另报 SNP 距离 vs Mash 距离 Spearman 相关（互验）。
"""
import csv, glob, math, os, time, collections

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")

accs, mats = [], []
for p in sorted(glob.glob(os.path.join(CS, "alleles", "*.npy"))):
    accs.append(os.path.basename(p)[:-4])
    mats.append(np.load(p))
M = np.stack(mats)
mats.clear()
n = len(accs)
print(f"genomes: {n}, ref sites: {M.shape[1]}", flush=True)

core_mask = (M < 4).sum(axis=0) >= int(0.95 * n)
C = M[:, core_mask]
core_len = C.shape[1]
print(f"core sites (>=95% taxa covered): {core_len}", flush=True)
# 低覆盖株剔除（组装质量差, 成对可比位点过少 → 距离不可比）: 树/回归用 subset
genome_cov = (M < 4).mean(axis=1)
keep_g = genome_cov >= 0.80
print(f"genomes with coverage>=0.80: {int(keep_g.sum())}/{n}（低覆盖剔除出树/回归）", flush=True)
# 变异位点数（速览, 双方可比才计）:
first = 0
v0 = (C != C[first]) & (C < 4)[first][None, :] & (C < 4)
print(f"variable vs first genome (comparable only): {int(v0.any(axis=0).sum())}", flush=True)

# --- 成对 SNP 距离（缺失不计入; 按可比位点数归一化为替换/位点） ---
print("pairwise SNP distances...", flush=True)
t0 = time.time()
idx_keep = np.nonzero(keep_g)[0]
nk = len(idx_keep)
Ck = C[idx_keep]
Vk = (Ck < 4)
D = np.zeros((nk, nk), dtype=np.float64)
NCOMP = np.zeros((nk, nk), dtype=np.float64)
CH_R, CH_C = 8, 120_000
for i0 in range(0, nk, CH_R):
    i1 = min(i0 + CH_R, nk)
    for j0 in range(0, core_len, CH_C):
        j1 = min(j0 + CH_C, core_len)
        a = Ck[i0:i1, j0:j1]
        v = Vk[i0:i1, j0:j1]
        b = Ck[:, j0:j1]
        bv = Vk[:, j0:j1]
        both = v[:, None, :] & bv[None, :, :]
        diff = (a[:, None, :] != b[None, :, :]) & both
        D[i0:i1, :] += diff.sum(axis=2)
        NCOMP[i0:i1, :] += both.sum(axis=2)
    print(f"  rows {i1}/{nk} ({time.time()-t0:.0f}s)", flush=True)
with np.errstate(invalid="ignore", divide="ignore"):
    D = np.where(NCOMP > 0, D / np.maximum(NCOMP, 1), 0.0)
print(f"normalization: comparable sites min pair = {NCOMP[NCOMP>0].min():.0f}", flush=True)

accs_k = [accs[i] for i in idx_keep]
with open(os.path.join(CS, "MY9_SNP距离.phylip"), "w", encoding="utf-8") as f:
    f.write(f"{nk}\n")
    for i in range(nk):
        f.write(accs_k[i] + "\t" + " ".join(f"{D[i,j]:.6f}" for j in range(nk)) + "\n")
print("phylip written", flush=True)

# --- NJ（与 MY9_12 同算法） ---
print("NJ tree building...", flush=True)
t0 = time.time()
active = list(range(nk))
node_name = {i: accs_k[i] for i in range(nk)}
D2 = {i: dict(zip(range(nk), D[i].tolist())) for i in range(nk)}
br = {}
next_id = n
while len(active) > 2:
    m = len(active)
    r = {a: sum(D2[a][b] for b in active) for a in active}
    best, bi, bj = 1e18, None, None
    for x, a in enumerate(active):
        Da = D2[a]
        for b in active[x + 1:]:
            v = Da[b] * (m - 2) - r[a] - r[b]
            if v < best:
                best, bi, bj = v, a, b
    dab = D2[bi][bj]
    lb = 0.5 * dab + (r[bi] - r[bj]) / (2 * (m - 2))
    lc = dab - lb
    nid = next_id; next_id += 1
    D2[nid] = {}
    br[nid] = [(bi, max(lb, 1e-6)), (bj, max(lc, 1e-6))]
    node_name[nid] = None
    for a in active:
        if a in (bi, bj):
            continue
        dnew = 0.5 * (D2[bi][a] + D2[bj][a] - dab)
        D2[nid][a] = dnew
        D2[a][nid] = dnew
    D2[nid][nid] = 0.0
    active = [a for a in active if a not in (bi, bj)] + [nid]
a, b = active
dab = D2[a][b]
root_id = next_id
br[root_id] = [(a, dab / 2), (b, dab / 2)]
print(f"NJ done in {time.time()-t0:.0f}s", flush=True)

def build(node):
    if node_name.get(node):
        return node_name[node]
    return "(" + ",".join(f"{build(c)}:{L:.6f}" for c, L in br[node]) + ")"

with open(os.path.join(CS, "MY9_coreSNP_NJ.nwk"), "w", encoding="utf-8") as f:
    f.write(build(root_id) + ";\n")

# --- 中点定根 + root-to-tip ---
import collections
adj = collections.defaultdict(list)
for node, children in br.items():
    for c, L in children:
        adj[node].append((c, L))
        adj[c].append((node, L))

def farthest(src):
    dist = {src: 0.0}
    stack = [src]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + L
                stack.append(v)
    w = max(dist, key=dist.get)
    return w, dist

u, _ = farthest(root_id)
v, du = farthest(u)
dv = {v: 0.0}
stack = [v]
while stack:
    x = stack.pop()
    for y, L in adj[x]:
        if y not in dv:
            dv[y] = dv[x] + L
            stack.append(y)
diam = du[v]
id_of = {accs[i]: i for i in range(n)}

mat = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
rows = []
for i, name in enumerate(accs_k):
    dmid = (du[i] + dv[i] - diam) / 2.0
    mt = mat.get(name, {})
    rows.append({"tip": name, "year": mt.get("year", ""), "ST": mt.get("MLST_ST", ""),
                 "resist": mt.get("resist_23S", ""), "P1": mt.get("P1_cluster", ""),
                 "root_dist": f"{dmid:.6f}"})
with open(os.path.join(CS, "MY9_coreSNP_root_to_tip.tsv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

def regress(pred):
    pts = [(int(r["year"]), float(r["root_dist"])) for r in rows
           if r["year"].isdigit() and 2005 <= int(r["year"]) <= 2026 and pred(r)]
    if len(pts) < 10:
        return None
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    syy = sum((y - my) ** 2 for _, y in pts)
    return len(pts), sxy / sxx, (sxy * sxy) / (sxx * syy)

# --- Mash 互验 ---
mash = {}
lines = open(os.path.join(RESDIR, "MY9_mash距离.phylip"), encoding="utf-8").read().strip().splitlines()
mnames = [ln.split("\t")[0] for ln in lines[1:]]
MD = [[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]]
mid = {nm: k for k, nm in enumerate(mnames)}
xs, ys = [], []
for i in range(nk):
    for j in range(i + 1, nk):
        if accs_k[i] in mid and accs_k[j] in mid:
            xs.append(D[i, j]); ys.append(MD[mid[accs_k[i]]][mid[accs_k[j]]])
from scipy.stats import spearmanr
rho = spearmanr(xs, ys).statistic if len(xs) > 100 else float("nan")

lines_out = ["# MY9 core-SNP 摘要（本地锚链法, stage2）", "",
             f"- 基因组 {n}; 参考位点 {M.shape[1]}; **core(≥95%) = {core_len}**; 树/回归用株 {nk}（覆盖≥0.80）",
             f"- SNP 距离(替换/位点) vs Mash 距离 Spearman ρ = {rho:.4f}（n={len(xs)} 对）", ""]
for label, pred in [("all (2005–2026)", lambda r: True),
                    ("ST3", lambda r: r["ST"] == "3"),
                    ("ST14", lambda r: r["ST"] == "14"),
                    ("div-A", lambda r: r["P1"].startswith("div-A")),
                    ("ST2 (WT clone)", lambda r: r["ST"] == "2")]:
    out = regress(pred)
    if out:
        nn, slope, r2 = out
        lines_out.append(f"- root-to-tip {label}: n={nn}, slope={slope:.3e} SNP/年, R²={r2:.3f}")
lines_out += ["", "> 缺失成对忽略; indel 邻域 ±8nt 掩蔽; 无重组过滤（Gubbins 类工具不可用, 论文限制项声明）。"]
with open(os.path.join(CS, "MY9_coreSNP_摘要.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines_out) + "\n")
print("\n".join(lines_out))
print("stage2 done -> results/coresnp/", flush=True)
