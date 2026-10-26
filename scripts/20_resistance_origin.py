# -*- coding: utf-8 -*-
"""20_resistance_origin.py — ancestral state reconstruction of resistance mutations (A2063G origin and timing)
Inputs: results/coresnp/treetime/{ST3_subtree.nwk, ST3.aln.fasta, ST3_dates.csv}
        results/coresnp/alleles/*.npy (rebuild of the core column mask, same convention as MY9_18)
Methods:
  1) **Root-independent Fitch parsimony** (primary evidence): count the minimum number of A/G changes
     on the NJ topology and locate the change edges and the downstream tips -> the "single origin vs
     multiple independent mutations" call does not depend on rooting or on the clock;
  2) TreeTime ML joint ancestral reconstruction (fixed Tc; coalescent='opt' is not used — it hangs in practice):
     gives rough dates for the gain nodes (with a caveat that the dating is not robust).
Sites (M129 genome, 0-based, double-checked via MY9_04 call_site geometry + the phenotype-discriminating column):
     2063 = 122118 | 2064 = 122119 | 2617 = 122672
Output: the resistance-origin report (Markdown) under results/coresnp/
"""
import csv, glob, os
import numpy as np
from Bio import Phylo

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
TTDIR = os.path.join(CS, "treetime")
CORE_SITES = 742308
SITES = {"A2063": 122118, "A2064": 122119, "C2617": 122672}
WT = {"A2063": "A", "A2064": "A", "C2617": "C"}
MUT = {"A2063": "G", "A2064": "G", "C2617": "G"}

# --- core column mask (exactly the same convention as MY9_18) ---
accs_all = [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(CS, "alleles", "*.npy")))]
M = np.stack([np.load(os.path.join(CS, "alleles", a + ".npy")) for a in accs_all])
n = M.shape[0]
core_mask = (M < 4).sum(axis=0) >= int(0.95 * n)
ci = np.nonzero(core_mask)[0]
print(f"core columns: {len(ci)}", flush=True)
col_of = {}
for name, g0 in SITES.items():
    w = np.nonzero(ci == g0)[0]
    col_of[name] = int(w[0]) if len(w) else -1
    print(f"{name} genome0={g0} -> alignment col {col_of[name]}", flush=True)

# --- tip states: read the column directly from the ST3 alignment (no TreeTime dependency) ---
seqs = {}
name = None
cur = []
for line in open(os.path.join(TTDIR, "ST3.aln.fasta")):
    line = line.strip()
    if line.startswith(">"):
        if name:
            seqs[name] = "".join(cur)
        name = line[1:]; cur = []
    elif line:
        cur.append(line)
if name:
    seqs[name] = "".join(cur)
print(f"ST3 alignment sequences: {len(seqs)}", flush=True)

tree = Phylo.read(os.path.join(TTDIR, "ST3_subtree.nwk"), "newick")
tips = tree.get_terminals()
tip_state = {}
for t in tips:
    t.name = t.name.strip()
    s = seqs.get(t.name)
    tip_state[t.name] = s[col_of["A2063"]] if s else "?"

lines = []
lines.append("# MY9 耐药突变祖先状态重建（A2063G 起源与年代）\n")
lines.append(f"- 输入: ST3 子树 {len(tips)} 株, core {len(ci):,} 位点\n")
lines.append("- 位点坐标双重核实: MY9_04 call_site 几何（up+24）与表型判别列"
             "（423 株 G / 95 株 A, M129 参考 GGACGGAAAGACC 基序 AA 双位点）一致\n")

for name in ("A2063", "A2064", "C2617"):
    col = col_of[name]
    if col < 0:
        lines.append(f"\n## {name}: 不在 core（覆盖率不足），跳过\n"); continue
    W, Mu = WT[name], MUT[name]
    st = {}
    for t in tips:
        s = seqs.get(t.name)
        st[t.name] = (s[col] if s else "?")
    n_wt = sum(1 for v in st.values() if v == W)
    n_mut = sum(1 for v in st.values() if v == Mu)
    n_unk = len(st) - n_wt - n_mut
    # ---- Fitch parsimony (iterative post-order + pre-order assignment) ----
    post = list(tree.find_clades(order="postorder"))
    sets = {}
    changes = 0
    for nd in post:
        if nd.is_terminal():
            v = st.get(nd.name, "?")
            sets[id(nd)] = None if v in ("?", "N") else {v}
        else:
            child_sets = [sets[id(c)] for c in nd.clades if sets[id(c)] is not None]
            if not child_sets:
                sets[id(nd)] = None; continue
            inter = set.intersection(*child_sets)
            if inter:
                sets[id(nd)] = inter
            else:
                sets[id(nd)] = set.union(*child_sets); changes += 1
    # pre-order assignment; record the change edges
    assign = {}
    root_set = sets[id(tree.root)]
    assign[id(tree.root)] = sorted(root_set)[0] if root_set else None
    edges = []
    for nd in tree.find_clades(order="preorder"):
        ps = assign.get(id(nd))
        if nd.is_terminal():
            continue
        for c in nd.clades:
            cs = sets[id(c)]
            if ps is None or cs is None:
                assign[id(c)] = ps if ps is not None else (sorted(cs)[0] if cs else None)
            elif ps in cs:
                assign[id(c)] = ps
            else:
                assign[id(c)] = sorted(cs)[0]
                if len(cs) == 1:
                    edges.append((nd, c, ps, sorted(cs)[0]))
                else:
                    changes += 1  # polymorphic subtree: requires an extra change
    lines.append(f"\n## {name}（WT={W}, 突变={Mu}）\n")
    lines.append(f"- 末梢状态: WT {n_wt} / 突变 {n_mut} / 未定或缺失 {n_unk}\n")
    lines.append(f"- **Fitch 简约法（根无关）: 最小变化次数 = {changes}**\n")
    for par, ch, a, b in edges:
        n_desc_mut = len([t for t in ch.get_terminals() if st.get(t.name) == Mu])
        lines.append(f"  - 变化边 {par.name}→{ch.name}: {a}→{b}, 下游末梢 {len(ch.get_terminals())}"
                     f"（其中 {Mu} {n_desc_mut}）\n")
    if changes == 1:
        lines.append(f"- **结论: {name}{Mu} 为单次起源**（树上最小 1 次变化）\n")
    elif changes > 1:
        lines.append(f"- {name}{Mu} 需 {changes} 次独立变化（或回复突变）\n")

# --- TreeTime ML ancestral reconstruction (rough dates; see the dating report for robustness) ---
try:
    from treetime import TreeTime
    with open(os.path.join(TTDIR, "ST3_dates.csv"), encoding="utf-8-sig") as f:
        dates = {r["name"]: float(r["date"]) for r in csv.DictReader(f)}
    tt = TreeTime(gtr="Jukes-Cantor", tree=os.path.join(TTDIR, "ST3_subtree.nwk"),
                  aln=os.path.join(TTDIR, "ST3.aln.fasta"), dates=dates, verbose=0)
    tt.run(root="oldest", infer_gtr=False, resolve_polytomies=True,
           coalescent=None, Tc=1e-6, max_iter=2)
    rate = float(tt.date2dist.clock_rate)
    lines.append(f"\n## TreeTime ML 祖先重建（参考年代, 定年稳健性有保留）\n")
    lines.append(f"- rate={rate:.3e} subs/site/yr = {rate*CORE_SITES:.2f} SNP/genome/yr, "
                 f"R²={float(getattr(tt.date2dist,'r_val',float('nan'))):.4f}, root={tt.tree.root.numdate:.0f}\n")
    for name in ("A2063",):
        col = col_of[name]; W, Mu = WT[name], MUT[name]
        stt = {}
        for nd in tt.tree.find_clades():
            s = getattr(nd, "sequence", None)
            stt[id(nd)] = (s[col] if s is not None else "N")
        gains = []
        for par in tt.tree.find_clades():
            for ch in par.clades:
                if stt[id(par)] == W and stt[id(ch)] == Mu:
                    gains.append(ch)
        tips_mut = [t for t in tt.tree.root.get_terminals() if stt[id(t)] == Mu]
        lines.append(f"- {name}: ML 增益事件 {len(gains)} 次; 突变末梢 {len(tips_mut)}\n")
        for g in gains:
            lines.append(f"  - 增益节点 {g.name}: 年代 {g.numdate:.0f}, "
                         f"下游末梢 {len(g.get_terminals())}\n")
        if len(tips_mut) > 1:
            mc = tt.tree.root.common_ancestor(tips_mut)
            lines.append(f"  - 突变株 MRCA 年代: {mc.numdate:.0f}"
                         f"（含 {len(mc.get_terminals())} 末梢）\n")
except Exception as e:
    lines.append(f"\n( TreeTime 段失败: {type(e).__name__}: {e} )\n")

out = os.path.join(CS, "MY9_耐药起源.md")
open(out, "w", encoding="utf-8").write("".join(lines))
print("written", out, flush=True)
print("".join(lines))
