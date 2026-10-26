# -*- coding: utf-8 -*-
"""33_cluster_substructure.py — Sub-structure of the 352-member cluster + sampling sensitivity of ancestral-state reconstruction

Two analyses, both derived entirely from existing single-source data, without altering any upstream result:
  A. Cluster sub-structure: within the 352-genome 25-SNP single-linkage group, apply single-linkage again at 10/15/20 SNP thresholds,
     to answer "is the 2019 wave a resurgence of the same sub-lineage, or a different sub-lineage taking over?" (review part 6, path 2).
     Wave definition: <=2018 / 2019-2020 / 2023-2025 (consistent with the year histogram in FigS2).
  B. Sampling sensitivity: take the set difference between the v1_backup download manifest (included=1, 532 genomes in the pre-backfill panel) and the current panel
     to identify the ENA backfill batch, prune the backfilled genomes from the current 642-tip NJ tree, and recompute the minimum number of Fitch changes at site 2063 on the same topology
     (whole tree + ST3 subtree) — this turns "how adding the Japanese time series (all wild type before 2000, three resistance-associated waves after 2000)
     affects ancestral-state reconstruction" into reproducible numbers (an executable substitute for review part 6, path 1: same topology controlled, only the leaf set changes).
     Self-check: without pruning, the authoritative values from table6 must be reproduced (whole tree 34 / ST3 subtree 6); otherwise an error is raised.

Input: the core-SNP SNP distance matrix under results/coresnp/ (pairwise distances); the tight-clone membership
      table under results/coresnp/ (group members); the integrated per-genome typing matrix under results/
      (genotype/year/country/ST); the v1 backup download manifest (backfill-batch identification);
      the core-SNP NJ tree and the ST3 subtree under results/coresnp/ (topology).
Output: a Markdown report and a TSV table covering the cluster sub-structure and sampling sensitivity, written under results/coresnp/.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
CS = RES / "coresnp"

WAVES = [("<=2018", lambda y: y <= 2018), ("2019-2020", lambda y: 2019 <= y <= 2020),
         ("2023-2025", lambda y: y >= 2023)]

mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(RES / "MY9_分型矩阵.tsv", encoding="utf-8-sig"), delimiter="\t")}
members = [r["ena_acc"] for r in csv.DictReader(
    open(CS / "MY9_紧克隆成员.tsv", encoding="utf-8-sig"), delimiter="\t")]
assert len(members) == 352, f"clone members {len(members)} != 352"

# ---- Pairwise distances (phylip variant: "acc\td1 d2 ..." with a full row per genome, same parsing convention as MY9_24) ----
L = open(CS / "MY9_SNP距离.phylip", encoding="utf-8").read().strip().splitlines()
n = int(L[0].split()[0])
accs, D = [], {}
for ln in L[1:]:
    parts = ln.split("\t")
    if len(parts) < 2:
        continue
    acc = parts[0].strip()
    accs.append(acc)
    D[acc] = {j: float(v) for j, v in enumerate(parts[1].split())}
assert len(accs) == n, f"phylip parse {len(accs)} != {n}"
_idx = {a: i for i, a in enumerate(accs)}
def dist(a, b):
    return D[a][_idx[b]] if _idx[b] in D[a] else D[b][_idx[a]]

# ---- A. Within-group sub-structure ----
def single_linkage(items, thr):
    parent = {x: x for x in items}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    m = len(items)
    for i in range(m):
        for j in range(i + 1, m):
            if dist(items[i], items[j]) <= thr:
                parent[find(items[i])] = find(items[j])
    groups = defaultdict(list)
    for x in items:
        groups[find(x)].append(x)
    return sorted(groups.values(), key=len, reverse=True)

def year_of(a):
    y = mat[a]["year"]
    return int(y) if y.isdigit() else None

lines_md = ["# MY9_簇亚结构与采样敏感性（MY9_33, 2026-09-30）", "",
            "## A. 352 株 25-SNP 组内亚结构（single-linkage, 组内距离矩阵）", ""]
rows_tsv = [["threshold", "cluster_rank", "n", "n_le2018", "n_2019_2020", "n_2023_2025", "top_countries", "accs"]]
summary = {}
for thr in (10, 15, 20):
    cls = single_linkage(members, thr)
    big = [c for c in cls if len(c) >= 5]
    waves_big = []
    for rank, c in enumerate(big, 1):
        yc = Counter()
        for a in c:
            y = year_of(a)
            for wname, cond in WAVES:
                if y is not None and cond(y):
                    yc[wname] += 1
        cc = Counter(mat[a]["country"].split(":")[0].strip() for a in c)
        topc = ", ".join(f"{k} {v}" for k, v in cc.most_common(3))
        waves_big.append({w: yc.get(w, 0) for w, _ in WAVES})
        rows_tsv.append([thr, rank, len(c), yc.get("<=2018", 0), yc.get("2019-2020", 0),
                         yc.get("2023-2025", 0), topc, " ".join(sorted(c))])
    lines_md.append(f"## 阈值 {thr} SNP: ≥5 株亚簇 {len(big)} 个, 最大 {len(big[0]) if big else 0} 株, "
                    f"其余 {len(members) - sum(len(c) for c in big)} 株散布于 <5 株小簇")
    for rank, (c, wc) in enumerate(zip(big, waves_big), 1):
        lines_md.append(f"- 亚簇 {rank}: n={len(c)}, 波次构成 ≤2018={wc['<=2018']}, 2019-2020={wc['2019-2020']}, 2023-2025={wc['2023-2025']}")
    summary[thr] = {"big": big, "waves": waves_big}
    lines_md.append("")

# ---- B. Sampling-sensitivity Fitch ----
# Backfill-batch identification (single source): set difference between the v1_backup download manifest, included=1 (532 genomes in the pre-backfill panel), and the current panel
led = list(csv.DictReader(open(BASE / "MY9_下载清单_v1_backup.csv", encoding="utf-8-sig")))
pre_panel = {r["ena_acc"] for r in led if r["included"] == "1" and r["trio_ok"] == "1"}
recovered = {a for a in set(mat) - pre_panel}  # backfill batch = current panel − v1 included set
lines_md += ["## B. 祖先态重建采样敏感性（同拓扑剪枝对照）", "",
             f"- ENA 回填批（v1 面板 included=1 共 {len(pre_panel)} 株之外, 当前在面板的回填株）: {len(recovered)} 株"]

def tip_state(a):
    g = mat.get(a, {}).get("site2063", "")
    return g if g in ("A", "G", "T", "C") else None

def fitch_changes(tree_obj):
    """Fitch for site 2063 (same algorithm as MY9_26/MY9_29): post-order, bottom-up; when the intersection of the child sets is empty, take the union and add 1 change.
    State = the site2063 base; missing/uncalled leaves are unknown (None, skipped)."""
    state, changes = {}, 0
    for nd in tree_obj.find_clades(order="postorder"):
        if nd.is_terminal():
            s = tip_state(nd.name)
            state[id(nd)] = None if s is None else frozenset(s)
        else:
            cs = [state[id(c)] for c in nd.clades if state[id(c)] is not None]
            if not cs:
                state[id(nd)] = None
                continue
            inter = set(cs[0])
            for s in cs[1:]:
                inter &= s
            if inter:
                state[id(nd)] = frozenset(inter)
            else:
                state[id(nd)] = frozenset(set().union(*cs))
                changes += 1
    return changes

tree = Phylo.read(str(CS / "MY9_coreSNP_NJ.nwk"), "newick")
tree_tips = [t.name for t in tree.get_terminals()]
full_all = fitch_changes(tree)
sub0 = Phylo.read(str(CS / "treetime" / "ST3_subtree.nwk"), "newick")
sub_all = fitch_changes(sub0)
assert full_all == 34, f"full-tree Fitch self-check {full_all} != 34"
assert sub_all == 6, f"ST3-subtree Fitch self-check {sub_all} != 6"

sub = Phylo.read(str(CS / "treetime" / "ST3_subtree.nwk"), "newick")
sub_tips = [t.name for t in sub.get_terminals()]

pruned_n = len([t for t in tree_tips if t not in recovered])
sub_pruned_n = len([t for t in sub_tips if t not in recovered])
for t in [t for t in tree.get_terminals() if t.name in recovered]:
    tree.prune(t)
for t in [t for t in sub.get_terminals() if t.name in recovered]:
    sub.prune(t)
full_pre = fitch_changes(tree)
sub_pre = fitch_changes(sub)
assert len([t for t in tree.get_terminals()]) == pruned_n
bf_jp = sum(1 for a in recovered if "Japan" in mat[a]["country"])
bf_res = sum(1 for a in recovered if mat[a]["resist_23S"] == "resistant")
bf_pre_wt = sum(1 for a in recovered if mat[a]["year"].isdigit() and int(mat[a]["year"]) < 2000
                and mat[a]["resist_23S"] == "WT")
bf_pre = sum(1 for a in recovered if mat[a]["year"].isdigit() and int(mat[a]["year"]) < 2000)
assert bf_pre == bf_pre_wt, "pre-2000 backfill contains non-WT assemblies"
lines_md += [
    f"- 当前面板（646）: 全树 {full_all} 次变化, ST3 子树 {sub_all} 次（与 table6 权威值一致, 自检通过）",
    f"- 剪去回填批后（树 {pruned_n} tips, 其中 ST3 子树 {sub_pruned_n} tips）: 全树 {full_pre} 次, ST3 子树 {sub_pre} 次",
    f"- 回填批构成: {len(recovered)} 株（日本 {bf_jp} 株）; 2000 年前 {bf_pre} 株全部野生型; 2000-2024 含耐药相关 {bf_res} 株",
    f"- 结论: 加入日本时间序列使 ST3 子树内最小变化数从 {sub_pre} 升至 {sub_all}, 全树从 {full_pre} 升至 {full_all}",
]

out_md = CS / "MY9_簇亚结构与采样敏感性.md"
out_md.write_text("\n".join(lines_md) + "\n", newline="")
with open(CS / "MY9_簇亚结构与采样敏感性.tsv", "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerows(rows_tsv)
    f.write("\n")
    w.writerow(["metric", "value"])
    w.writerow(["recovered_batch_n", len(recovered)])
    w.writerow(["fitch_full_current", full_all])
    w.writerow(["fitch_st3_current", sub_all])
    w.writerow(["fitch_full_prerecovery", full_pre])
    w.writerow(["fitch_st3_prerecovery", sub_pre])
    w.writerow(["prerecovery_tree_tips", pruned_n])
    w.writerow(["recovered_batch_japan_n", bf_jp])
    w.writerow(["expansion_pre2000_wt", bf_pre_wt])
    w.writerow(["expansion_resistant", bf_res])
print(out_md)
print("\n".join(lines_md))
