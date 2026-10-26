# -*- coding: utf-8 -*-
"""24_manuscript_tables.py — 稿件表格 CSV（表格直读 CSV, 数字不落正文）
输入: results/MY9_分型矩阵.tsv、results/MY9_耐药频率_按年份.tsv、results/MY9_流行克隆定义.md、
      results/coresnp/MY9_SNP距离.phylip
输出: manuscript/tables/table1_panel.csv ... table5_temporal.csv
"""
import csv, os, re
import numpy as np
from collections import Counter, defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
CS = os.path.join(RES, "coresnp")
OUT = os.path.join(BASE, "manuscript", "tables")
os.makedirs(OUT, exist_ok=True)
CORE_SITES = 742308

mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
rows = list(mat.values())
N = len(rows)

def w(name, header, data):
    p = os.path.join(OUT, name)
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(header)
        cw.writerows(data)
    print(f"{name}: {len(data)} 行")

# ---- Table 1: panel by country/region (n>=5) ----
def disp_st(k):
    """ST 展示命名规范化（评审 P2）: uncertain(near X)->X-related, novel(near X)->X-novel。"""
    if k.startswith("uncertain(near "):
        return k.replace("uncertain(near ", "").rstrip(")") + "-related"
    if k.startswith("novel(near "):
        return k.replace("novel(near ", "").rstrip(")") + "-novel"
    return k

by_c = defaultdict(lambda: {"n": 0, "res": 0, "wt": 0, "years": [], "st": Counter()})
for r in rows:
    c = r["country"].split(":")[0].strip() or "not applicable"
    d = by_c[c]
    d["n"] += 1
    if r["resist_23S"] == "resistant": d["res"] += 1
    elif r["resist_23S"] == "WT": d["wt"] += 1
    if r["year"].isdigit(): d["years"].append(int(r["year"]))
    d["st"][r["MLST_ST"] or "novel/uncertain"] += 1
order = sorted(by_c, key=lambda c: -by_c[c]["n"])
t1 = []
for c in order:
    d = by_c[c]
    if d["n"] < 5: continue
    called = d["res"] + d["wt"]
    pct = f"{d['res']/called*100:.1f}" if called else "NA"
    stx = ", ".join(f"ST{disp_st(k)} ({v})" for k, v in d["st"].most_common(3))
    t1.append([c, d["n"], d["res"], d["wt"],
               f"{pct}", f"{np.median(d['years']):.0f}" if d["years"] else "NA", stx])
tot_res = sum(1 for r in rows if r["resist_23S"] == "resistant")
tot_wt = sum(1 for r in rows if r["resist_23S"] == "WT")
_allyears = [int(r["year"]) for r in rows if r["year"].isdigit()]
_allst = Counter(r["MLST_ST"] or "novel/uncertain" for r in rows if r["resist_23S"] in ("resistant", "WT"))
_allstx = ", ".join(f"ST{disp_st(k)} ({v})" for k, v in _allst.most_common(3))
t1.append(["All included genomes", N, tot_res, tot_wt,
           f"{tot_res/(tot_res+tot_wt)*100:.1f}",
           f"{np.median(_allyears):.0f}" if _allyears else "NA", _allstx])
w("table1_panel.csv", ["Region", "N_genomes", "N_resistant", "N_wildtype",
                       "Pct_resistant", "Median_year", "Top_MLST"], t1)

# ---- Table 2: by year (n>=5) ----
by_y = defaultdict(lambda: {"n": 0, "res": 0, "wt": 0, "g": 0, "t": 0, "c": 0,
                            "a2064g": 0, "a2064c": 0, "c2617": 0})
for r in rows:
    y = r["year"]
    if not y.isdigit(): continue
    y = int(y)
    if y < 2005: continue
    d = by_y[y]
    d["n"] += 1
    if r["resist_23S"] == "resistant": d["res"] += 1
    elif r["resist_23S"] == "WT": d["wt"] += 1
    if r["site2063"] == "G": d["g"] += 1
    elif r["site2063"] == "T": d["t"] += 1
    elif r["site2063"] == "C": d["c"] += 1
    if r["site2064"] == "G": d["a2064g"] += 1
    if r["site2064"] == "C": d["a2064c"] += 1
    if r["site2617"] == "G": d["c2617"] += 1
t2 = []
for y in sorted(by_y):
    d = by_y[y]
    if d["n"] < 5: continue
    cal = d["res"] + d["wt"]
    t2.append([y, d["n"], d["res"], d["wt"], f"{d['res']/cal*100:.1f}" if cal else "NA",
               d["g"], d["t"], d["c"], d["a2064g"], d["a2064c"], d["c2617"]])
w("table2_by_year.csv", ["Year", "N_genomes", "N_resistant", "N_wildtype", "Pct_resistant",
                         "A2063G", "A2063T", "A2063C", "A2064G", "A2064C", "C2617G"], t2)

# ---- Table 3: by MLST ST / P1 type (n>=10) ----
def disp_st(k):
    """ST 展示命名规范化（评审 P2）: uncertain(near X)->X-related, novel(near X)->X-novel。"""
    if k.startswith("uncertain(near "):
        return k.replace("uncertain(near ", "").rstrip(")") + "-related"
    if k.startswith("novel(near "):
        return k.replace("novel(near ", "").rstrip(")") + "-novel"
    return k

# 国家/地区展示名规范化（与 MY9_28 同口径; NCBI 原文只在 S1 保留）
NORM = {"Taiwan": "Taiwan, China", "China": "Mainland China", "USA": "United States",
        "not applicable": "Not specified"}
def norm_regions(s):
    out = []
    for part in s.split(", "):
        name, _, cnt = part.rpartition(" ")
        out.append((NORM.get(name, name) + " " + cnt) if name else part)
    return ", ".join(out)

def strat(keyfn):
    d = defaultdict(lambda: {"n": 0, "res": 0, "wt": 0})
    for r in rows:
        k = keyfn(r)
        if not k: continue
        e = d[k]; e["n"] += 1
        if r["resist_23S"] == "resistant": e["res"] += 1
        elif r["resist_23S"] == "WT": e["wt"] += 1
    return d
t3 = []
for label, keyfn in [("MLST ST", lambda r: ("ST" + disp_st(r["MLST_ST"])) if r["MLST_ST"] else None),
                     ("P1 type", lambda r: ("P1 " + r["P1_type"]) if r["P1_type"] else None)]:
    d = strat(keyfn)
    for k in sorted(d, key=lambda x: -d[x]["n"]):
        e = d[k]
        if e["n"] < 10: continue
        cal = e["res"] + e["wt"]
        t3.append([label, k, e["n"], e["res"], e["wt"],
                   f"{e['res']/cal*100:.1f}" if cal else "NA"])
w("table3_by_type.csv", ["Stratifier", "Group", "N_genomes", "N_resistant", "N_wildtype",
                         "Pct_resistant"], t3)

# ---- Table 4: epidemic clone thresholds（复用 MY9_23 结果, 直接解析其 md 表） ----
t4 = []
md = open(os.path.join(RES, "MY9_流行克隆定义.md"), encoding="utf-8").read()
for m in re.finditer(r"最大簇 n=(\d+)（占 (\d+)%）, 内部成对中位 (\d+) SNP", md):
    pass
for m in re.finditer(r"## 阈值 t = (\d+) SNP\n- 最大簇 n=(\d+)（占 (\d+)%）, 内部成对中位 (\d+) SNP\n"
                     r"- ST 构成: ([^\n]+)\n- 耐药 (\d+) / 敏感 (\d+) / 未call (\d+)（耐药率 ([\d.]+)%）\n"
                     r"- 年份: 中位 (\d+), 范围 (\d+)–(\d+); ≥2019 占 (\d+)%\n"
                     r"- 国家/地区 Top: ([^\n]+)\n", md):
    t4.append([m.group(1), m.group(2), m.group(3), m.group(4),
               f"{m.group(6)}/{int(m.group(6))+int(m.group(7))}", m.group(9),
               m.group(10), norm_regions(m.group(14))])
w("table4_clone.csv", ["SNP_threshold", "N_cluster", "Pct_of_panel", "Median_within_SNP",
                       "Resistant_called", "Pct_resistant", "Median_year", "Top_regions"], t4)

# ---- Table 5: temporal signal（成对回归口径 + TreeTime 网格跨度） ----
L = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
accs = [ln.split("\t")[0] for ln in L[1:]]
D = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in L[1:]])
D = (D + D.T) / 2 * CORE_SITES
st = np.array([mat.get(a, {}).get("MLST_ST", "") or "novel" for a in accs])
yr = np.array([float(mat[a]["year"]) if mat.get(a, {}).get("year", "").isdigit() else np.nan for a in accs])
def pairfit(mask):
    idx = np.nonzero(mask)[0]
    ys = yr[idx]
    dy, dd = [], []
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            if np.isnan(ys[i]) or np.isnan(ys[j]): continue
            dt = abs(ys[i] - ys[j])
            if dt > 0:
                dy.append(dt); dd.append(D[idx[i], idx[j]])
    dy = np.array(dy); dd = np.array(dd)
    s, b = np.polyfit(dy, dd, 1)
    return s, b, np.corrcoef(dy, dd)[0, 1], len(dy)
t5 = []
for lab, mask in [("All ST3 genomes", st == "3"),
                  ("Epidemic clone (<=25 SNP single-linkage)", None)]:
    if mask is None:
        # 复用 MY9_23 的簇定义
        p = list(range(len(accs)))
        def find(x):
            while p[x] != x:
                p[x] = p[p[x]]; x = p[x]
            return x
        iu, ju = np.nonzero(np.triu(D <= 25, 1))
        for a, b in zip(iu, ju):
            ra, rb = find(int(a)), find(int(b))
            if ra != rb: p[rb] = ra
        g = defaultdict(list)
        for i in range(len(accs)): g[find(i)].append(i)
        big = max(g.values(), key=len)
        m = np.zeros(len(accs), bool); m[big] = True
        mask = m
    s, b, r, n = pairfit(mask)
    t5.append([lab, f"{s:.2f}", f"{s/2:.2f}", f"{r:.3f}", f"{b:.1f}", n])
# TreeTime 网格跨度
tt = os.path.join(CS, "treetime", "tt_summary.tsv")
if os.path.exists(tt):
    vals = [float(r["rate_snp_genome_yr"]) for r in csv.DictReader(open(tt, encoding="utf-8-sig"), delimiter="\t")
            if r["rate_snp_genome_yr"] not in ("err",)]
    r2s = [float(r["R2"]) for r in csv.DictReader(open(tt, encoding="utf-8-sig"), delimiter="\t")
           if r["R2"] not in ("err",)]
    t5.append(["TreeTime ML grid (5 settings)", f"{min(vals):.2f} to {max(vals):.2f}", "not identifiable",
               f"{min(r2s):.3f} to {max(r2s):.3f}", "NA", len(vals)])
w("table5_temporal.csv", ["Dataset", "Pairwise_rate_SNP_per_yr_both_lineages",
                          "Per_lineage_rate", "r", "Intercept_SNP", "N_pairs"], t5)

# ---- Table 6 (统计量单源): 正文引用的非表内数值（规则 11/33） ----
import re as _re
summ = open(os.path.join(CS, "MY9_coreSNP_摘要.md"), encoding="utf-8").read()
rho = (_re.search(r"Spearman ρ = ([\d.]+)", summ) or [])[1]
epi = open(os.path.join(RES, "MY9_流行病学动态.md"), encoding="utf-8").read()
mfit = _re.search(r"K=([\d.]+)±([\d.]+).*?r=([\d.]+)±([\d.]+)/年, \*\*t50=([\d.]+)±([\d.]+)\*\*, R²=([\d.]+)", epi)
origin = open(os.path.join(CS, "MY9_耐药起源.md"), encoding="utf-8").read()
mfit2 = _re.search(r"最小变化次数 = (\d+)", origin)
mmut = _re.search(r"A2063: ML 增益事件 (\d+) 次; 突变末梢 (\d+)", origin)
ml_sec = origin.split("## TreeTime ML")[-1]
mgain = sorted([int(x) for x in _re.findall(r"增益节点 \S+: 年代 -?\d+, 下游末梢 (\d+)", ml_sec)], reverse=True)
rr = list(csv.DictReader(open(os.path.join(RES, "MY9_耐药频率_按年份.tsv"), encoding="utf-8-sig"), delimiter="\t"))
tot_g = sum(int(x["n_2063G"]) for x in rr if x["group"].strip().isdigit())
tot_t = sum(int(x["n_2063T"]) for x in rr if x["group"].strip().isdigit())
tot_c = sum(int(x["n_2063C"]) for x in rr if x["group"].strip().isdigit())
tot_2064 = sum(int(x["n_2064G"]) for x in rr if x["group"].strip().isdigit())
tot_2617 = sum(int(x["n_2617G"]) for x in rr if x["group"].strip().isdigit())
t6 = [
    ["core_sites", "742308", "core positions (>=95% of assemblies)"],
    ["tree_tips", str(len(accs)), "genomes in the core-SNP tree (coverage >= 0.80)"],
    ["snp_vs_mash_spearman", rho, "cross-check of the anchor-chain distances"],
    ["logistic_K", mfit.group(1), "weighted logistic asymptote"],
    ["logistic_r_per_year", mfit.group(3), "weighted logistic growth rate"],
    ["logistic_t50", mfit.group(5), "weighted logistic midpoint (year)"],
    ["logistic_R2", mfit.group(7), "weighted logistic fit"],
    ["fitch_changes_2063", mfit2.group(1), "minimum changes at position 2063 (Fitch, root-free)"],
    ["ml_gain_events", mmut.group(1), "maximum-likelihood gain events at position 2063 (ST3 subtree)"],
    ["resistant_tips_in_subtree", mmut.group(2), "resistant tips in the ST3 subtree"],
    ["main_gain_node_tips", str(mgain[0]) if mgain else "NA", "descendant tips of the largest ML gain node"],
    ["A2063G_total", str(tot_g), "genomes with A2063G"],
    ["A2063T_total", str(tot_t), "genomes with A2063T"],
    ["A2063C_total", str(tot_c), "genomes with A2063C"],
    ["A2064G_total", str(tot_2064), "genomes with A2064G"],
    ["A2064C_total", str(sum(1 for r in rows if r["site2064"] == "C")), "genomes with A2064C"],
    ["C2617G_total", str(tot_2617), "genomes with C2617G"],
    ["coord_2063", "122119", "1-based M129 genome position of A2063"],
    ["coord_2064", "122120", "1-based M129 genome position of A2064"],
    ["coord_2617", "122673", "1-based M129 genome position of C2617"],
    ["taxid", "2104", "NCBI taxonomy identifier (Mycoplasmoides pneumoniae)"],
    ["ec_equivalents", "2058/2059/2611", "E. coli numbering equivalents"],
]
# ST3 内部离散度与 ST3<->ST17 距离（直接由距离矩阵计算, 不解析中文散文）
_i3 = np.nonzero(st == "3")[0]
_sub3 = D[np.ix_(_i3, _i3)]
_iu3 = np.triu_indices(len(_i3), 1)
_v3 = _sub3[_iu3]
_i17 = np.nonzero(np.array([mat.get(a, {}).get("MLST_ST", "") == "17" for a in accs]))[0]
t6 += [
    ["st3_within_median", f"{np.median(_v3):.0f}", "median within-ST3 pairwise SNP distance"],
    ["st3_within_p90", f"{np.percentile(_v3, 90):.0f}", "p90 within-ST3 pairwise SNP distance"],
    ["st3_within_p99", f"{np.percentile(_v3, 99):.0f}", "p99 within-ST3 pairwise SNP distance"],
    ["st3_within_max", f"{_v3.max():.0f}", "maximum within-ST3 pairwise SNP distance"],
    ["st3_vs_st17_median", f"{np.median(D[np.ix_(_i3, _i17)]):.0f}", "median ST3-ST17 pairwise SNP distance"],
    ["st3_subtree_tips", str(len(list(__import__("Bio.Phylo", fromlist=["read"]).read(os.path.join(CS, "treetime", "ST3_subtree.nwk"), "newick").get_terminals()))), "tips in the ST3 subtree used for dating and ancestral reconstruction"],
    ["other_regions", str(N - sum(int(r[1]) for r in t1[:-1])), "genomes from regions with fewer than five genomes, including unspecified"],
    ["python_versions", "3.13/3.12", "analysis runtimes (main environment and TreeTime virtual environment)"],
]
# 全树 Fitch 结果（MY9_26 产物）并入 table6, 保持单源链路（规则 142）
_ft = open(os.path.join(BASE, "results", "coresnp", "MY9_全树祖先重建.md"), encoding="utf-8").read()
_mch = _re.search(r"fitch_full_changes,(\d+)", _ft)
_mga = _re.search(r"fitch_full_gains,(\d+)", _ft)
# 四类会计分解（评审 P0-1）: 变化边方向计数, 直接解析权威 md 的边清单
_cls = Counter(_re.findall(r"- 变化边 (A→G|G→A|A→T|A→C):", _ft))
# 克隆单次获得检验点估计（评审 P0-2）: MY9_32 扰动检验 md 的真实树行
_rb = open(os.path.join(CS, "MY9_树稳健性_扰动.md"), encoding="utf-8").read()
_mcg = _re.search(r"克隆成员置 G/其余未知时所需变为 G 的边数 (\d+)", _rb)
_mpt = _re.search(r"全树变更 (\d+)–(\d+)（中位 \d+）", _rb)
_mps = _re.search(r"ST3 子树内 (\d+)–(\d+)（中位 \d+）", _rb)
t6 += [
    ["fitch_full_changes", _mch.group(1) if _mch else "NA", "Fitch changes at 2063 over the whole panel"],
    ["fitch_full_gains", _mga.group(1) if _mga else "NA", "independent A2063G gain edges over the whole panel"],
    ["fitch_full_reversions", _cls.get("G→A", 0), "Fitch reversions to A at 2063 (G→A edges, displayed topology)"],
    ["fitch_full_A2063T_gains", _cls.get("A→T", 0), "Fitch A→T transitions at 2063"],
    ["fitch_full_A2063C_gains", _cls.get("A→C", 0), "Fitch A→C transitions at 2063"],
    ["fitch_clone_member_gains", _mcg.group(1) if _mcg else "NA",
     "additional gains required to explain clone members' shared A2063G (non-members unknown; 0 = single ancestral gain compatible)"],
    ["perturb_total_changes", (_mpt.group(1) + "-" + _mpt.group(2)) if _mpt else "NA",
     "range of total changes over 15 NNI-perturbed trees (MY9_32)"],
    ["perturb_st3_within", (_mps.group(1) + "-" + _mps.group(2)) if _mps else "NA",
     "range of within-ST3 changes over 15 NNI-perturbed trees (MY9_32)"],
]
# 簇亚结构与采样敏感性（MY9_33 产物并入 table6, 保持单源链路; 与全树重建 md/ST3 权威值对账）
_ss_rows = list(csv.reader(open(os.path.join(CS, "MY9_簇亚结构与采样敏感性.tsv"), encoding="utf-8-sig"), delimiter="\t"))
_ss_m, _ss_thr, _mode = {}, {}, None
for _r in _ss_rows:
    if not _r:
        continue
    if _r[0] == "metric" and len(_r) == 2:
        _mode = "m"
        continue
    if _mode is None and len(_r) == 8 and _r[0].isdigit():
        if _r[1] == "1":
            _ss_thr[_r[0]] = int(_r[2])
    elif _mode == "m" and len(_r) == 2:
        _ss_m[_r[0]] = _r[1]
assert _ss_m.get("fitch_full_current") == (_mch.group(1) if _mch else None), "MY9_33 fitch_full 与 MY9_26 对账失败"
assert _ss_m.get("fitch_st3_current") == "6", "MY9_33 fitch_st3 与权威值 6 对账失败"
assert _ss_m.get("recovered_batch_japan_n") == _ss_m.get("recovered_batch_n"), "回填批并非全部来自日本, 稿件描述需复核"
t6 += [
    ["panel_expansion_japan", _ss_m["recovered_batch_n"], "Japanese assemblies added via ENA BioSample metadata recovery (1976-2024)"],
    ["expansion_pre2000_wt", _ss_m["expansion_pre2000_wt"], "pre-2000 assemblies in the expansion batch (all wild type at the assayed loci)"],
    ["expansion_resistant", _ss_m["expansion_resistant"], "resistance-associated assemblies within the expansion batch (2000-2024)"],
    ["fitch_full_prerecovery", _ss_m["fitch_full_prerecovery"], "whole-panel Fitch minimum changes after pruning the expansion batch (same topology)"],
    ["fitch_st3_prerecovery", _ss_m["fitch_st3_prerecovery"], "within-ST3 Fitch minimum changes after pruning the expansion batch (same topology)"],
    ["prerecovery_tree_tips", _ss_m["prerecovery_tree_tips"], "tips retained after pruning the expansion batch (same topology)"],
    ["subthr_max_component", "/".join(str(_ss_thr[t]) for t in ("10", "15", "20")), "largest within-group component when single-linkage is repeated at 10/15/20 SNPs (352 = group not split)"],
]
# ST7 唯一耐药株与 2000 年前面板核查（第七轮评审 P1-2/P1-3; 由分型矩阵直接计算）
def _cc2(c): return c.split(":")[0].strip()
def _mut(r):
    if r["site2063"] != "A": return "A2063" + r["site2063"]
    if r["site2064"] != "A": return "A2064" + r["site2064"]
    if r["site2617"] != "C": return "C2617" + r["site2617"]
    return "other"
def _disp_c(c):
    c0, _, reg = c.partition(":")
    name = NORM.get(c0.strip(), c0.strip())
    return name + (" (" + reg.strip() + ")" if reg.strip() else "")
_st7r = [r for r in rows if (r["MLST_ST"] or "").strip() == "7" and r["resist_23S"] == "resistant"]
_pre = [r for r in rows if r["year"].isdigit() and int(r["year"]) < 2000]
_prer = sorted([r for r in _pre if r["resist_23S"] == "resistant"], key=lambda r: r["year"])
_parts = [f"{_disp_c(r['country'])} in {r['year']} (ST{disp_st(r['MLST_ST'])})" for r in _prer]
_pre_detail = (" and ".join(_parts) if len(_parts) == 2 else "; ".join(_parts)) or "none"
# pre-2000 耐药株邻接关系（评审 P1-1; 克隆阈值口径, D 为 SNP 单位）
_pre_t = [r["ena_acc"] for r in _prer]
_n25 = {}
_near = {}
for a in _pre_t:
    _i = accs.index(a)
    _row = D[_i].copy()
    _row[_i] = np.inf
    _n25[a] = int((_row <= 25).sum())
    _j = int(np.argmin(_row))
    _near[a] = float(_row[_j])
_row99 = D[accs.index(_pre_t[1])].copy()
_row99[accs.index(_pre_t[1])] = np.inf
_jacc = accs[int(np.argmin(_row99))]
_jr = mat[_jacc]
# 评审第十一轮: 两株纽约 1999 的独立性——入组池元数据核验 BioSample 登记号
# 池 CSV 键带前缀与版本(GCF_002128205.1), 矩阵键为 GCA_002128205, 统一取数字段对齐
_pool = {r["assembly_accession"].split("_")[1].split(".")[0]: r for r in csv.DictReader(
    open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig"))}
_bs99 = [_pool[a.split("_")[1]]["biosample"] for a in (mat[_pre_t[1]]["ena_acc"], _jacc)]
assert _bs99[0] and _bs99[1] and _bs99[0] != _bs99[1], f"NY1999 BioSample check failed: {_bs99}"
# same-city 升级: 完整 country 字符串逐字比对(USA: New York), 不再只比国家字段
_wt_same_city_year = int(
    _jr["resist_23S"] == "WT"
    and _jr["year"] == mat[_pre_t[1]]["year"]
    and _jr["country"] == mat[_pre_t[1]]["country"]
) if _pre_t else 0
t6 += [
    ["pre2000_neighbors_le25", "/".join(str(_n25[a]) for a in _pre_t),
     "neighbours within 25 SNPs for the pre-2000 resistance-associated assemblies (Denmark 1988 first, then New York 1999)"],
    ["pre2000_nearest_snp", "/".join(f"{_near[a]:.0f}" for a in _pre_t),
     "SNP distance to the nearest neighbour for the pre-2000 resistance-associated assemblies (Denmark 1988 first, then New York 1999)"],
    ["pre2000_1999_nearest_wt_same_city_year", str(_wt_same_city_year),
     "nearest neighbour of the 1999 New York assembly is a wild-type assembly sampled the same year in the same city (1 = yes)"],
    ["pre2000_1999_distinct_biosamples", "1",
     "the 1999 New York resistance-associated assembly and its nearest wild-type neighbour are distinct BioSample records"],
    ["pre2000_1999_biosample_ids", "/".join(_bs99),
     "BioSample accessions: the 1999 New York resistance-associated assembly first, then its nearest wild-type neighbour"],
]
t6 += [
    ["st7_res_mutation", _mut(_st7r[0]) if _st7r else "NA", "mutation of the single resistance-associated exact-ST7 assembly"],
    ["st7_res_country", _disp_c(_st7r[0]["country"]) if _st7r else "NA", "country of the single resistance-associated exact-ST7 assembly"],
    ["st7_res_year", _st7r[0]["year"] if _st7r else "NA", "sampling year of the single resistance-associated exact-ST7 assembly"],
    ["pre2000_assemblies", str(len(_pre)), "assemblies sampled before 2000"],
    ["pre2000_resistant", str(len(_prer)), "resistance-associated assemblies sampled before 2000"],
    ["pre2000_resistant_detail", _pre_detail, "pre-2000 resistance-associated assemblies (country, year, ST)"],
]
w("table6_analysis_stats.csv", ["Key", "Value", "Note"], t6)
print("tables ->", OUT)
