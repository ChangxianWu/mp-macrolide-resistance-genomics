# -*- coding: utf-8 -*-
"""17_datetree.py — 数据驱动定根的 root-to-tip 定年（TempEst/RRT 式, 纯本地）
输入: results/coresnp/MY9_SNP距离.phylip（替换/位点）+ MY9_分型矩阵.tsv
输出: results/coresnp/MY9_定年报告.md + 控制台

方法: 目标群 → NJ 子树（无根）→ 扫描每条边上的定根位置（粗扫 1/4·1/2·3/4 + 最优边细化）
      → 取 root-to-tip 回归 R² 最大的根 → 报告: 斜率(SNP/年), R², 根日期(回归零点 = 该群样本 MRCA 估计)。
      与免根成对回归互为印证; 真正的贝叶斯定年（BEAST2）可后续在服务器抽检。
"""
import csv, os, collections
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
CORE_LEN = 742308

lines = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
nk = int(lines[0])
accs = [ln.split("\t")[0] for ln in lines[1:]]
D = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]])
idx = {a: i for i, a in enumerate(accs)}
mat = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}

def nj(sub):
    m0 = len(sub)
    active = list(range(m0))
    name = {i: accs[sub[i]] for i in range(m0)}
    Dsub = D[np.ix_(sub, sub)]
    D2 = {i: dict(zip(range(m0), Dsub[i].tolist())) for i in range(m0)}
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
        br[nid] = [(bi, max(lb, 1e-12)), (bj, max(dab - lb, 1e-12))]
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
    return br, name

def adjacency(br):
    adj = collections.defaultdict(list)
    for nd, ch in br.items():
        for c, L in ch:
            adj[nd].append((c, L))
            adj[c].append((nd, L))
    return adj

def bfs_depths(adj, src, blocked=None):
    dist = {src: 0.0}
    stack = [src]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v == blocked or v in dist:
                continue
            dist[v] = dist[u] + L
            stack.append(v)
    return dist

def best_root_rtt(tips, years, br):
    """无根树上的最优定根 root-to-tip。tips: 叶节点 id 集; years: 叶→年份。
    返回 (best_R2, slope_SNP_per_year, root_date, best_edge, best_frac)"""
    adj = adjacency(br)
    leaves = [nd for nd in adj if nd in years]
    best = None
    seen_edges = set()
    for u in list(adj.keys()):
        for v, w in adj[u]:
            e = (min(u, v), max(u, v))
            if e in seen_edges or u == v:
                continue
            seen_edges.add(e)
            comp_u = bfs_depths(adj, u, blocked=v)
            comp_v = bfs_depths(adj, v, blocked=u)
            for f in (0.25, 0.5, 0.75):
                xs, ys = [], []
                for t in leaves:
                    y = years.get(t)
                    if y is None:
                        continue
                    if t in comp_u:
                        d = comp_u[t] + f * w
                    else:
                        d = comp_v[t] + (1 - f) * w
                    xs.append(y)
                    ys.append(d * CORE_LEN)
                if len(xs) < 10:
                    continue
                mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
                sxx = sum((x - mx) ** 2 for x in xs)
                if sxx == 0:
                    continue
                sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                syy = sum((y - my) ** 2 for y in ys)
                slope0 = sxy / sxx
                if slope0 <= 0:      # 只接受正时钟（负斜率=反时钟, 生物上不可能）
                    continue
                r2 = sxy * sxy / (sxx * syy) if syy > 0 else 0.0
                if best is None or r2 > best[0]:
                    a = my - slope0 * mx
                    root_date = -a / slope0 if slope0 != 0 else float("nan")
                    best = (r2, slope0, root_date, (u, v, w), f)
    # 细化最优边
    if best:
        r2b, _, _, (u, v, w), _ = best
        comp_u = bfs_depths(adj, u, blocked=v)
        comp_v = bfs_depths(adj, v, blocked=u)
        for f in np.arange(0.05, 0.96, 0.05):
            xs, ys = [], []
            for t in leaves:
                y = years.get(t)
                if y is None:
                    continue
                d = comp_u[t] + f * w if t in comp_u else comp_v[t] + (1 - f) * w
                xs.append(y); ys.append(d * CORE_LEN)
            mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
            sxx = sum((x - mx) ** 2 for x in xs)
            sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            syy = sum((y - my) ** 2 for y in ys)
            slope0 = sxy / sxx if sxx > 0 else 0.0
            if slope0 <= 0:
                continue
            r2 = sxy * sxy / (sxx * syy) if syy > 0 else 0.0
            if r2 > best[0]:
                a = my - slope0 * mx
                best = (r2, slope0, -a / slope0 if slope0 != 0 else float("nan"), (u, v, w), float(f))
    return best

def analyze(group_name, pred):
    sub = [idx[a] for a in accs if pred(mat.get(a, {}))]
    if len(sub) < 15:
        return f"- {group_name}: n={len(sub)}, 跳过"
    br, name = nj(sub)
    # 叶节点 → 年份
    years = {}
    for nd, nm in name.items():
        if nm is not None and br.get(nd):
            pass
    leaves = [nd for nd, nm in name.items() if nm is not None]
    years = {}
    for nd in leaves:
        y = mat.get(name[nd], {}).get("year", "")
        if y.isdigit() and 2005 <= int(y) <= 2026:
            years[nd] = int(y)
    if len(years) < 15:
        return f"- {group_name}: 年份点 {len(years)} 不足"
    res = best_root_rtt(leaves, years, br)
    if res is None:
        return (f"- **{group_name}**: n(tips with year)={len(years)}, "
                f"扫描全部定根位置均无正斜率解 → root-to-tip 定根未识别（不报告速率/MRCA）")
    r2, slope, rdate, edge, frac = res
    return (f"- **{group_name}**: n(tips with year)={len(years)}, 最优定根 R²={r2:.3f}, "
            f"速率={slope:.1f} SNP/年, 样本 MRCA 估计≈{rdate:.1f} 年")

out = ["# MY9 定年报告（数据驱动定根 root-to-tip, RRT/TempEst 式）", "",
       "> 单位: SNP/年（core 742,308 位点）; 根日期 = 最优定根下回归零点 ≈ 该群采样样本的 MRCA 年份。",
       "> 结果与免根成对回归、TreeTime 网格运行互为印证; BEAST2/TreeTime 正式定年可在可联网环境以已导出的比对+树+日期复跑。", ""]
out.append(analyze("ST3 全体（东亚耐药克隆）", lambda m: m.get("MLST_ST") == "3"))
out.append(analyze("ST3 近期波（2023+）", lambda m: m.get("MLST_ST") == "3" and m.get("year", "").isdigit() and int(m["year"]) >= 2022))
out.append(analyze("ST14（次级克隆）", lambda m: m.get("MLST_ST") == "14"))
out.append(analyze("type2/ST2+ST14 系", lambda m: m.get("P1_type") == "type2"))
with open(os.path.join(CS, "MY9_定年报告.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
print("\n".join(out))
