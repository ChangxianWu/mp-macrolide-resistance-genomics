# -*- coding: utf-8 -*-
"""15b_clock.py — SNP 级分子钟检验（子树内定根）+ 分层 ρ 诊断
输入: results/coresnp/MY9_SNP距离.phylip（MY9_15 产物, 替换/位点）+ MY9_mash距离.phylip + 分型矩阵
输出: results/coresnp/MY9_分子钟检验.md + 控制台

方法: 对每个目标群（ST3/ST14/div-A）取其 SNP 距离子矩阵 → NJ → 子树中点定根 →
      root-to-tip 对收集年份回归（正斜率+高 R² = 群体内分子钟成立）。
      Mash 互验按距离分层（Mash sketch 在 ~1e-4 量级有噪声底, 全域 ρ 被紧簇打乱属预期）。
"""
import csv, os, time, collections
import numpy as np
from scipy.stats import spearmanr

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")

lines = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
nk = int(lines[0])
accs = [ln.split("\t")[0] for ln in lines[1:]]
D = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]])
idx = {a: i for i, a in enumerate(accs)}
mat = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}

def nj(sub):
    """sub: 索引列表 → (br dict, root_id, node_name)"""
    m0 = len(sub)
    active = list(range(m0))
    name = {i: accs[sub[i]] for i in range(m0)}
    D2 = {i: dict(zip(range(m0), D[np.ix_(sub, sub)][i].tolist())) for i in range(m0)}
    br = {}
    nxt = m0
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
        nid = nxt; nxt += 1
        D2[nid] = {}
        br[nid] = [(bi, max(lb, 1e-9)), (bj, max(dab - lb, 1e-9))]
        name[nid] = None
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
    root = nxt
    br[root] = [(a, dab / 2), (b, dab / 2)]
    return br, root, name

def depths_from(br, src, root):
    adj = collections.defaultdict(list)
    for nd, ch in br.items():
        for c, L in ch:
            adj[nd].append((c, L))
            adj[c].append((nd, L))
    dist = {src: 0.0}
    stack = [src]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + L
                stack.append(v)
    return dist, adj

def farthest(adj, src):
    dist = {src: 0.0}
    stack = [src]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + L
                stack.append(v)
    w = max(dist, key=dist.get)
    return w, dist[w]

def clock(group_name, pred, min_n=15):
    sub = [idx[a] for a in accs if pred(mat.get(a, {}))]
    if len(sub) < min_n:
        return f"- {group_name}: n={len(sub)} < {min_n}, 跳过"
    br, root, name = nj(sub)
    # 中点定根: 找直径两端
    adj = collections.defaultdict(list)
    for nd, ch in br.items():
        for c, L in ch:
            adj[nd].append((c, L))
            adj[c].append((nd, L))
    u, _ = farthest(adj, root)
    v, diam = farthest(adj, u)
    du, _ = depths_from(br, u, root)
    dv, _ = depths_from(br, v, root)
    pts = []
    for i in range(len(sub)):
        nm = name[i]
        if nm is None:
            continue
        mt = mat.get(nm, {})
        y = mt.get("year", "")
        if y.isdigit() and 2005 <= int(y) <= 2026:
            dmid = (du[i] + dv[i] - diam) / 2.0
            pts.append((int(y), dmid))
    if len(pts) < 10:
        return f"- {group_name}: n={len(pts)} 年份点不足"
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    syy = sum((y - my) ** 2 for _, y in pts)
    slope = sxy / sxx
    r2 = (sxy * sxy) / (sxx * syy)
    # 换算 SNP/年: 乘 core 长度（距离单位=替换/位点）
    core_len = 742308
    return f"- {group_name}: n={len(pts)}, slope={slope*core_len:.1f} SNP/年, R²={r2:.3f}"

def _adj_of(br):
    adj = collections.defaultdict(list)
    for nd, ch in br.items():
        for c, L in ch:
            adj[nd].append((c, L))
            adj[c].append((nd, L))
    return adj

out = ["# MY9 SNP 级分子钟检验（子树内中点定根）", ""]
out.append(clock("ST3（东亚耐药克隆）", lambda m: m.get("MLST_ST") == "3"))
out.append(clock("ST14", lambda m: m.get("MLST_ST") == "14"))
out.append(clock("div-A（台湾分歧 P1 谱系）", lambda m: str(m.get("P1_cluster", "")).startswith("div-A")))
out.append(clock("ST2（美国野生型对照）", lambda m: m.get("MLST_ST") == "2"))

# --- 分层 ρ ---
lines = open(os.path.join(RESDIR, "MY9_mash距离.phylip"), encoding="utf-8").read().strip().splitlines()
mnames = [ln.split("\t")[0] for ln in lines[1:]]
MD = [[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]]
mid = {nm: k for k, nm in enumerate(mnames)}
strata = collections.defaultdict(list)
for i in range(nk):
    for j in range(i + 1, nk):
        if accs[i] in mid:
            s = D[i, j]
            k = ("<10" if s < 10 / 742308 else "10-100" if s < 100 / 742308 else
                 "100-500" if s < 500 / 742308 else "500-2000" if s < 2000 / 742308 else ">2000")
            strata[k].append((s, MD[mid[accs[i]]][mid[accs[j]]]))
out += ["", "## SNP 距离 vs Mash 距离 分层 Spearman ρ（Mash 在紧簇端分辨率不足属预期）", ""]
for k in ["<10", "10-100", "100-500", "500-2000", ">2000"]:
    pts = strata.get(k, [])
    if len(pts) >= 30:
        rho = spearmanr([x for x, _ in pts], [y for _, y in pts]).statistic
        out.append(f"- {k} SNP 层: n={len(pts)}, ρ={rho:.3f}")
out += ["", "> 单位换算: 距离矩阵为替换/位点; SNP 数 ≈ 距离 × 742308。"]
with open(os.path.join(CS, "MY9_分子钟检验.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))
