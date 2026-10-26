# -*- coding: utf-8 -*-
"""22_figures_en.py — 投稿版英文图件（Fig 1–4）
输入: results/MY9_分型矩阵.tsv（分型+元数据）
      results/MY9_P1聚类.tsv（P1 谱系）
      results/MY9_耐药频率_按年份.tsv（频率×年）
      results/coresnp/MY9_coreSNP_NJ.nwk（core-SNP NJ 树）
      results/coresnp/treetime/MY9_coreSNP_NJ.nwk（同拓扑副本）
输出: manuscript/figures/MY9_Fig1..4_*.{png,pdf}（英文标注, 300 dpi PNG）
数据逻辑与核心分析输出一致。
"""
import csv, os, collections, math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 300

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
FIGDIR = os.path.join(BASE, "manuscript", "figures")
os.makedirs(FIGDIR, exist_ok=True)
CORE_SITES = 742308

# 国家名规范化（投稿口径: 跟随 Lancet 惯例）
COUNTRY_NORM = {"Taiwan": "Taiwan, China", "China": "Mainland China",
                "USA": "United States", "not applicable": "Not specified"}
def cnorm(c):
    return COUNTRY_NORM.get(c, c)


mat = {r["ena_acc"]: r for r in csv.DictReader(
    open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
rows = list(mat.values())
n_total = len(rows)

# ---------------- Fig 1 ----------------
by_year = collections.defaultdict(lambda: [0, 0])
for r in rows:
    y = r["year"]
    if y.isdigit() and 2005 <= int(y) <= 2026 and r["resist_23S"] in ("resistant", "WT"):
        by_year[int(y)][0 if r["resist_23S"] == "resistant" else 1] += 1
years = [y for y in sorted(by_year) if sum(by_year[y]) >= 5]
freqs = [by_year[y][0] / sum(by_year[y]) * 100 for y in years]
ns = [sum(by_year[y]) for y in years]

by_ctry = collections.defaultdict(lambda: [0, 0])
for r in rows:
    c = cnorm(r["country"].split(":")[0].strip())
    if r["resist_23S"] in ("resistant", "WT") and c and c != "Not specified":
        by_ctry[c][0 if r["resist_23S"] == "resistant" else 1] += 1
top_c = [c for c in sorted(by_ctry, key=lambda c: -sum(by_ctry[c]))
         if sum(by_ctry[c]) >= 5][:12]

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axes[0]
ax.bar(range(len(years)), freqs, color="#c0392b", alpha=0.88)
for i, (y, f, n) in enumerate(zip(years, freqs, ns)):
    ax.text(i, f + 2.5, f"{f:.0f}", ha="center", fontsize=8.5, color="#333")
ax.set_xticks(range(len(years)))
ax.set_xticklabels(years, rotation=45, fontsize=9.5)
ax.set_ylim(0, 118)
ax.set_ylabel("Resistance-associated genotype (%)", fontsize=11)
ax.set_xlabel("Collection year", fontsize=11)
ax.set_title("A  Genotype frequency by year (years with n \u2265 5)", fontsize=11.5)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
fs = [by_ctry[c][0] / sum(by_ctry[c]) * 100 for c in top_c]
nsc = [sum(by_ctry[c]) for c in top_c]
ax.barh(range(len(top_c))[::-1], fs, color="#2874a6", alpha=0.88)
for i, (c, f, n) in enumerate(zip(top_c, fs, nsc)):
    ax.text(f + 2, len(top_c) - 1 - i, f"{f:.0f}% (n={n})", va="center", fontsize=9)
ax.set_yticks(range(len(top_c))[::-1])
ax.set_yticklabels(top_c, fontsize=10)
ax.set_xlim(0, 132)
ax.set_xlabel("Resistance-associated genotype (%)", fontsize=11)
ax.set_title("B  Genotype frequency by country/region (n \u2265 5, top 12)", fontsize=11.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "MY9_Fig1_resistance_frequency.png"))
fig.savefig(os.path.join(FIGDIR, "MY9_Fig1_resistance_frequency.pdf"))
plt.close(fig)
print("Fig1 written")

# ---------------- Fig 2 ----------------
cl = {}
with open(os.path.join(RESDIR, "MY9_P1聚类.tsv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cl[r["ena_acc"]] = r["P1_cluster"].split("\u2248")[0]
def p1_group(v):
    if v.startswith("type1"): return "type1 (M129-like)"
    if v.startswith("type2"): return "type2 (FH-like)"
    if v.startswith("div-single"): return "divergent (singleton)"
    if v.startswith("type-mixed"): return "mixed/recombinant"
    if v.startswith("div"): return "divergent lineage"
    return "not extracted"
by_py = collections.defaultdict(collections.Counter)
for r in rows:
    y = r["year"]
    if y.isdigit() and 2005 <= int(y) <= 2026:
        by_py[int(y)][p1_group(cl.get(r["ena_acc"], "?"))] += 1
years2 = sorted(by_py)
groups = ["type1 (M129-like)", "type2 (FH-like)", "divergent lineage",
          "mixed/recombinant", "divergent (singleton)"]
colors = ["#c0392b", "#2874a6", "#e67e22", "#8e44ad", "#95a5a6"]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
bottom = np.zeros(len(years2))
for g, col in zip(groups, colors):
    vals = np.array([by_py[y].get(g, 0) for y in years2])
    ax.bar(range(len(years2)), vals, bottom=bottom, color=col, label=g, alpha=0.92)
    bottom += vals
ax.set_xticks(range(len(years2))); ax.set_xticklabels(years2, rotation=45, fontsize=9.5)
ax.set_ylabel("Number of genomes", fontsize=11); ax.set_xlabel("Collection year", fontsize=11)
ax.set_title(f"P1 lineage composition by year (2005\u20132026, n={n_total}, strand-normalized typing)", fontsize=11.5)
ax.legend(fontsize=9.5, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "MY9_Fig2_P1_lineages.png"))
fig.savefig(os.path.join(FIGDIR, "MY9_Fig2_P1_lineages.pdf"))
plt.close(fig)
print("Fig2 written")

# ---------------- Fig 3: NJ tree ----------------
def parse_newick(path):
    txt = open(path).read().strip().rstrip(";")
    pos = 0
    def node():
        nonlocal pos
        if txt[pos] == "(":
            pos += 1
            children = [node()]
            while txt[pos] == ",":
                pos += 1
                children.append(node())
            assert txt[pos] == ")"
            pos += 1
            name, L = label()
            return {"children": children, "name": name, "L": L}
        name, L = label()
        return {"children": [], "name": name, "L": L}
    def label():
        nonlocal pos
        start = pos
        while pos < len(txt) and txt[pos] not in ",():":
            pos += 1
        name = txt[start:pos].strip()
        L = 0.0
        if pos < len(txt) and txt[pos] == ":":
            pos += 1
            s2 = pos
            while pos < len(txt) and (txt[pos].isdigit() or txt[pos] in ".eE-"):
                pos += 1
            L = float(txt[s2:pos])
        return name, L
    return node()

root = parse_newick(os.path.join(CS, "MY9_coreSNP_NJ.nwk"))
leaf_y, counter = {}, [0]
def assign_y(nd):
    if not nd["children"]:
        leaf_y[id(nd)] = counter[0]; counter[0] += 1
    else:
        ys = [assign_y(c) for c in nd["children"]]
        leaf_y[id(nd)] = sum(ys) / len(ys)
    return leaf_y[id(nd)]
assign_y(root)
def depth(nd, d=0.0):
    nd["depth"] = d
    for c in nd["children"]:
        depth(c, d + c["L"])
depth(root)
def walk_all(nd):
    yield nd
    for c in nd["children"]:
        yield from walk_all(c)
allnodes = list(walk_all(root))
maxd = max(n["depth"] for n in allnodes)

C_RES, C_WT, C_NC = "#c0392b", "#2471a3", "#b2babb"
# 支序图布局（忽略分支长）: 大枝长离散（ST3 内 ~1e-5 vs ST3↔ST2 ~1e-3）使任何单轴都压缩,
# 故按等长呈现拓扑, 分支长信息以定量注释框给出。
def height(nd):
    if not nd["children"]:
        nd["h"] = 0
        return 0
    nd["h"] = 1 + max(height(c) for c in nd["children"])
    return nd["h"]
MAXH = height(root)
def X(nd):
    return MAXH - nd["h"]

# ST3 支系定位（子树内 ST3 占≥85% 且规模≥200 的最大子树）
st3_tips_all = {a for a, m in mat.items() if m.get("MLST_ST") == "3"}
def collect(nd):
    if not nd["children"]:
        return [nd["name"]]
    out = []
    for c in nd["children"]:
        out.extend(collect(c))
    return out
best = {"size": 0, "node": None, "tips": set()}
for nd in allnodes:
    if not nd["children"]:
        continue
    tips = collect(nd)
    k = sum(1 for t in tips if t in st3_tips_all)
    if len(tips) >= 200 and k / len(tips) >= 0.85 and len(tips) > best["size"]:
        best = {"size": len(tips), "node": nd, "tips": set(tips)}

# 定量注释: within-ST3（全体/2023 波）与 ST3↔ST2 中位成对 SNP
_l = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
_accs = [ln.split("\t")[0] for ln in _l[1:]]
_D = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in _l[1:]])
_D = (_D + _D.T) / 2 * CORE_SITES
_st = np.array([mat.get(a, {}).get("MLST_ST", "") for a in _accs])
_yr = np.array([float(mat[a]["year"]) if mat.get(a, {}).get("year", "").isdigit() else np.nan
                for a in _accs])
def _med_within(mask):
    sub = _D[np.ix_(np.nonzero(mask)[0], np.nonzero(mask)[0])]
    iu = np.triu_indices(sub.shape[0], 1)
    return np.median(sub[iu])
def _med_between(m1, m2):
    return np.median(_D[np.ix_(np.nonzero(m1)[0], np.nonzero(m2)[0])].ravel())
i3, i2 = (_st == "3"), (_st == "2")
w23 = i3 & (_yr >= 2023)
ann_lines = [f"Median pairwise SNP distance",
             f"Within ST3 (all): {_med_within(i3):.0f}",
             f"Within ST3 2023 wave: {_med_within(w23):.0f}",
             f"ST3 vs ST2: {_med_between(i3, i2):.0f}"]

# 流行紧克隆: 单连锁 ≤25 SNP 最大簇（与 MY9_23 同口径）
_p = list(range(len(_accs)))
def _find(x):
    while _p[x] != x:
        _p[x] = _p[_p[x]]; x = _p[x]
    return x
_iu, _ju = np.nonzero(np.triu(_D <= 25, 1))
for a, b in zip(_iu, _ju):
    ra, rb = _find(int(a)), _find(int(b))
    if ra != rb:
        _p[rb] = ra
_groups = {}
for i in range(len(_accs)):
    _groups.setdefault(_find(i), []).append(i)
_big = max(_groups.values(), key=len)
epi_tips = {_accs[i] for i in _big}
ann_lines.append(f"25-SNP single-linkage group (n={len(epi_tips)}): "
                 f"{sum(1 for i in _big if mat[_accs[i]]['resist_23S']=='resistant')} resistance-associated")

fig, ax = plt.subplots(figsize=(7.5, 13))
def draw(nd):
    x, y = X(nd), leaf_y[id(nd)]
    if nd["children"]:
        ch = [(X(c), leaf_y[id(c)]) for c in nd["children"]]
        ylo, yhi = min(v for _, v in ch), max(v for _, v in ch)
        if yhi == ylo:
            yhi = ylo + 1e-9
        ax.plot([x, x], [ylo, yhi], color="#555", lw=0.45, zorder=1)
        for c in nd["children"]:
            ax.plot([x, X(c)], [y, leaf_y[id(c)]], color="#555", lw=0.45, zorder=1)
            draw(c)
    else:
        m = mat.get(nd["name"], {})
        col = {"resistant": C_RES, "WT": C_WT}.get(m.get("resist_23S", "-"), C_NC)
        ax.plot(x, y, "o", color=col, ms=2.2, zorder=2)
        if nd["name"] in epi_tips:
            ax.plot(x, y, "o", mfc="none", mec="#0e6655", ms=4.0, mew=0.7, zorder=3)
draw(root)

if best["node"] is not None:
    ys_b = [leaf_y[id(n)] for n in allnodes
            if not n["children"] and n["name"] in best["tips"]]
    if ys_b:
        xb = X(best["node"])
        ax.plot([xb - 0.6, xb - 0.6], [min(ys_b), max(ys_b)], color="#16a085", lw=2.4, zorder=3)
        ax.text(xb - 0.85, (min(ys_b) + max(ys_b)) / 2,
                f"ST3 clonal\ncomplex\n(n={best['size']})", color="#16a085",
                fontsize=9, ha="right", va="center")

ax.set_xlim(-MAXH * 0.06, MAXH * 1.02)
ax.set_ylim(-2, counter[0] + 14)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlabel("Cladogram (branch lengths not to scale)", fontsize=9)
handles = [Rectangle((0, 0), 1, 1, color=c) for c in (C_RES, C_WT, C_NC)]
handles.append(plt.Line2D([], [], ls="none", marker="o", mfc="none", mec="#0e6655", ms=5, mew=0.8))
labels = ["macrolide-resistant (23S)", "wild-type", "not called",
          "25-SNP single-linkage group members"]
ax.legend(handles, labels, fontsize=8, loc="upper left")
ax.text(0.985, 0.012, "\n".join(ann_lines), transform=ax.transAxes, fontsize=8,
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", fc="#fdfefe", ec="#95a5a6", lw=0.6))
n_tips = sum(1 for n in allnodes if not n["children"])
ax.set_title(f"Core-SNP neighbour-joining tree (n={n_tips}, midpoint-rooted); tips coloured by 23S genotype",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "MY9_Fig3_coreSNP_tree.png"))
fig.savefig(os.path.join(FIGDIR, "MY9_Fig3_coreSNP_tree.pdf"))
plt.close(fig)
print("Fig3 written")

# ---------------- Fig 4: temporal dynamics ----------------
from scipy.optimize import curve_fit

def logistic(t, K, r, t50):
    return K / (1.0 + np.exp(-r * (t - t50)))

year_rows = [r for r in csv.DictReader(
    open(os.path.join(RESDIR, "MY9_耐药频率_按年份.tsv"), encoding="utf-8-sig"), delimiter="\t")
    if r.get("group", "").strip()]
yrs, fr, wts = [], [], []
for r in year_rows:
    g = r["group"].strip()
    if g.isdigit() and int(r["n_called"]) >= 5 and int(g) >= 1995:
        yrs.append(int(g)); fr.append(float(r["freq_resist"])); wts.append(int(r["n_called"]))
yrs = np.array(yrs, float); fr = np.array(fr, float); wts = np.array(wts, float)
popt, pcov = curve_fit(logistic, yrs, fr, p0=[0.95, 0.4, 2015.0],
                       sigma=1.0 / np.sqrt(wts), absolute_sigma=False, maxfev=20000)
perr = np.sqrt(np.diag(pcov))
resid = fr - logistic(yrs, *popt)
r2 = 1 - (resid ** 2).sum() / ((fr - np.average(fr, weights=wts)) ** 2).sum()

lines = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
accs = [ln.split("\t")[0] for ln in lines[1:]]
Dm = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]])
Dm = (Dm + Dm.T) / 2
years_arr = np.array([float(mat[a]["year"]) if mat.get(a, {}).get("year", "").isdigit() else np.nan
                      for a in accs])
grp = np.array([mat.get(a, {}).get("MLST_ST", "") for a in accs])
sel = np.nonzero((grp == "3") & np.isfinite(years_arr))[0]
ys = years_arr[sel]
dy, dd = [], []
for i in range(len(sel)):
    for j in range(i + 1, len(sel)):
        dt = abs(ys[i] - ys[j])
        if dt > 0:
            dy.append(dt); dd.append(Dm[sel[i], sel[j]] * CORE_SITES)
dy = np.array(dy); dd = np.array(dd)
slope, intercept = np.polyfit(dy, dd, 1)
r_pear = np.corrcoef(dy, dd)[0, 1]
bins = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30]
bx, by_, _ = [], [], []
for k in range(len(bins) - 1):
    m = (dy > bins[k]) & (dy <= bins[k + 1])
    if m.sum() >= 20:
        bx.append((bins[k] + bins[k + 1]) / 2); by_.append(np.median(dd[m]))

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
xs = np.linspace(1995, 2026, 200)
ax.plot(xs, logistic(xs, *popt), color="#c0392b", lw=2,
        label=f"logistic fit: K={popt[0]:.2f}, r={popt[1]:.2f}/yr, t50={popt[2]:.1f}")
ax.scatter(yrs, fr, s=wts * 4, color="#2c3e50", zorder=3, label="observed (size \u221d n typed)")
for y, f, w in zip(yrs, fr, wts):
    ax.annotate(f"{int(w)}", (y, f), textcoords="offset points",
                xytext=(0, 4 + 2.6 * np.sqrt(w)),
                ha="center", fontsize=8.5, color="#555")
ax.set_xlabel("Collection year", fontsize=11); ax.set_ylabel("Resistance-associated genotype proportion", fontsize=11)
ax.set_ylim(-0.05, 1.08)
ax.set_title("A  Annual genotype proportion and descriptive logistic curve", fontsize=11.5)
ax.legend(fontsize=9.5, loc="lower right", markerscale=0.22, handletextpad=0.6,
          borderpad=0.5, labelspacing=0.8)
ax.grid(alpha=0.3); ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
ax.scatter(dy, dd, s=3, alpha=0.05, color="#7f8c8d")
ax.plot(bx, by_, "-o", color="#c0392b", ms=6, lw=2, label="binned median")
xs2 = np.linspace(0, dy.max(), 50)
ax.plot(xs2, slope * xs2 + intercept, "--", color="#2980b9", lw=2,
        label=f"descriptive slope: {slope:.2f} SNP/yr (both lineages, r={r_pear:.2f})")
ax.set_xlabel("|\u0394 sampling year|", fontsize=11); ax.set_ylabel("Pairwise SNP distance", fontsize=11)
ax.set_title("B  Pairwise SNP distance vs sampling-year difference (ST3)", fontsize=11.5)
ax.legend(fontsize=9.5); ax.grid(alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "MY9_Fig4_temporal_dynamics.png"))
fig.savefig(os.path.join(FIGDIR, "MY9_Fig4_temporal_dynamics.pdf"))
plt.close(fig)
print("Fig4 written")
print("figures ->", FIGDIR)
