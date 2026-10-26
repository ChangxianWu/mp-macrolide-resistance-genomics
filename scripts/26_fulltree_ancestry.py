# -*- coding: utf-8 -*-
"""26_fulltree_ancestry.py — 全树剪枝面板（528 tips）A2063 Fitch 祖先重建
此前只在 ST3 子树内重建, 覆盖不了 P1 type2 等耐药株 → 全基因组获得次数未评估
状态: 每株 site2063 (A/G/T/C/缺失), 直接取自 results/MY9_分型矩阵.tsv
树: results/coresnp/MY9_coreSNP_NJ.nwk (528 tips)
输出: results/coresnp/MY9_全树祖先重建.md + 台账 CSV
"""
import csv, os
from Bio import Phylo

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
CS = os.path.join(RES, "coresnp")

mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
tree = Phylo.read(os.path.join(CS, "MY9_coreSNP_NJ.nwk"), "newick")
tips = tree.get_terminals()
state = {}
for t in tips:
    r = mat.get(t.name)
    s = (r["site2063"] if r and r["site2063"] in "AGTC" else "?")
    state[t.name] = s
    t.name = t.name.strip()

# Fitch: 自底向上集合 + 变化计数; 变化边定位（显式栈, 无递归风险）
post = list(tree.find_clades(order="postorder"))
sets = {}
changes = 0
for nd in post:
    if nd.is_terminal():
        v = state.get(nd.name, "?")
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
# 前序赋值 + 记录变化边
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
            edges.append((nd, c, ps, assign[id(c)], sorted(cs)))

lines = [f"# MY9 全树（{len(tips)} tips）A2063 祖先重建（Fitch, 根无关）\n"]
from collections import Counter
tipc = Counter(state.values())
lines.append(f"- 末梢状态: " + ", ".join(f"{k}={v}" for k, v in sorted(tipc.items())) + "\n")
lines.append(f"- **最小变化次数 = {changes}**\n")
for par, ch, a, b, cs in edges:
    d = [x.name for x in ch.get_terminals()]
    dres = [x for x in d if state.get(x) == "G"]
    lines.append(f"- 变化边 {a}→{b}: 下游末梢 {len(d)}（其中 G {len(dres)}）\n")
    lines.append(f"  - 该支构成: ST " + str(Counter(mat.get(x, {}).get('MLST_ST', '?') or 'novel' for x in d).most_common(4)) +
                 "; P1 " + str(Counter(mat.get(x, {}).get('P1_type', '?') for x in d).most_common(3)) + "\n")
# 每条 gain 支下游 G 计入账目
gains = [(ch, a, b) for _, ch, a, b, _ in edges if b == "G"]
tot_gain = sum(len([x for x in ch.get_terminals() if state.get(x) == "G"]) for ch, _, _ in gains)
lines.append(f"\n- 获得边（→G）{len(gains)} 条; 其下游 G 末梢合计 = {tot_gain} / {tipc.get('G', 0)}\n")
lines.append(f"- fitch_full_changes,{changes}\n- fitch_full_gains,{len(gains)}\n")
out_md = os.path.join(CS, "MY9_全树祖先重建.md")
open(out_md, "w", encoding="utf-8").write("".join(lines))
print("".join(lines))
print("written", out_md)
