# -*- coding: utf-8 -*-
"""21_epi_dynamics.py — temporal dynamics of resistance: logistic growth fit + root-free pairwise molecular clock + Fig4
Inputs: the per-year resistance frequency table under results/
        the integrated per-genome typing matrix under results/ (group/country/year)
        the core-SNP distance matrix (PHYLIP, substitutions/site, 528 strains) under results/coresnp/
Outputs: the epidemiological dynamics report (Markdown) under results/
         the Figure 4 panel (resistance rate + logistic fit; pairwise SNP vs Δyear) under results/figures/

Conventions:
- Logistic growth: freq(t) = K / (1 + exp(-r (t - t50))), weighted least squares (weights = number of strains called that year);
  only years with n_called>=5 are used; a sensitivity analysis excluding China/Taiwan is also reported
- Root-free pairwise molecular clock: pairwise comparisons within the same clone (ST3) with Δyear>0, SNP count = distance × 748,341;
  the slope is the combined rate of both lineages -> halved per lineage; this estimate does not depend on rooting
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
CS = os.path.join(RES, "coresnp")
FIG = os.path.join(RES, "figures")
os.makedirs(FIG, exist_ok=True)
CORE_SITES = 742308
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ---------- Inputs ----------
year_rows = []
with open(os.path.join(RES, "MY9_耐药频率_按年份.tsv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if not r.get("group", "").strip():
            continue
        year_rows.append(r)
print("年份行:", len(year_rows))

mat = {}
with open(os.path.join(RES, "MY9_分型矩阵.tsv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mat[r["ena_acc"]] = r

# ---------- A. Logistic growth ----------
def logistic(t, K, r, t50):
    return K / (1.0 + np.exp(-r * (t - t50)))

yrs, fr, wts = [], [], []
for r in year_rows:
    g = r["group"].strip()
    if not g.isdigit():
        continue
    y = int(g); n = int(r["n_called"])
    if n < 5 or y < 1995:
        continue
    yrs.append(y); fr.append(float(r["freq_resist"])); wts.append(n)
yrs = np.array(yrs, float); fr = np.array(fr, float); wts = np.array(wts, float)
print("拟合年份点:", dict(zip(yrs.astype(int).tolist(), wts.astype(int).tolist())))

p0 = [0.95, 0.4, 2015.0]
popt, pcov = curve_fit(logistic, yrs, fr, p0=p0, sigma=1.0 / np.sqrt(wts),
                       absolute_sigma=False, maxfev=20000)
perr = np.sqrt(np.diag(pcov))
resid = fr - logistic(yrs, *popt)
r2 = 1 - (resid ** 2).sum() / ((fr - np.average(fr, weights=wts)) ** 2).sum()
print(f"logistic K={popt[0]:.3f}±{perr[0]:.3f} r={popt[1]:.3f}±{perr[1]:.3f} "
      f"t50={popt[2]:.1f}±{perr[2]:.1f} R2={r2:.3f}")

# Sensitivity: exclude China/Taiwan
sub = {a: r for a, r in mat.items() if r["country"].split(":")[0].strip() not in ("China", "Taiwan")}
from collections import defaultdict
by_year = defaultdict(lambda: [0, 0])
for r in sub.values():
    y = r["year"].strip()
    if not y.isdigit():
        continue
    d = by_year[int(y)]
    d[1] += 1
    if r["resist_23S"] != "WT":
        d[0] += 1
ys2, fr2, wt2 = [], [], []
for y in sorted(by_year):
    nr, nt = by_year[y]
    if nt >= 5 and y >= 1995:
        ys2.append(y); fr2.append(nr / nt); wt2.append(nt)
print("剔除东亚后拟合点:", list(zip(ys2, wt2)))
try:
    popt2, pcov2 = curve_fit(logistic, np.array(ys2, float), np.array(fr2, float),
                             p0=p0, sigma=1.0 / np.sqrt(np.array(wt2, float)),
                             absolute_sigma=False, maxfev=20000)
    perr2 = np.sqrt(np.diag(pcov2))
except Exception as e:
    popt2 = perr2 = None
    print("东亚剔除拟合失败:", e)

# ---------- B. Root-free pairwise molecular clock ----------
lines = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
accs = [ln.split("\t")[0] for ln in lines[1:]]
Dm = np.array([[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]])
Dm = (Dm + Dm.T) / 2
years = np.array([float(mat[a]["year"]) if mat.get(a, {}).get("year", "").isdigit() else np.nan
                  for a in accs])
grp = np.array([mat.get(a, {}).get("MLST_ST", "") for a in accs])

sel = np.nonzero((grp == "3") & np.isfinite(years))[0]
ys = years[sel]
dy, dd = [], []
for i in range(len(sel)):
    for j in range(i + 1, len(sel)):
        dt = abs(ys[i] - ys[j])
        if dt > 0:
            dy.append(dt); dd.append(Dm[sel[i], sel[j]] * CORE_SITES)
dy = np.array(dy); dd = np.array(dd)
slope, intercept = np.polyfit(dy, dd, 1)
r_pear = np.corrcoef(dy, dd)[0, 1]
print(f"成对回归: n={len(dy)} 斜率={slope:.3f} SNP/年 r={r_pear:.3f} 截距={intercept:.1f}")
# Fraction of pairs with 0 SNP (close relatives from the same wave)
print("成对距离=0 比例: %.3f" % ((dd == 0).mean()))

bins = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30]
bx, by_, bn = [], [], []
for k in range(len(bins) - 1):
    m = (dy > bins[k]) & (dy <= bins[k + 1])
    if m.sum() >= 20:
        bx.append((bins[k] + bins[k + 1]) / 2); by_.append(np.median(dd[m])); bn.append(int(m.sum()))

# ---------- Report ----------
rep = ["# MY9 流行病学时间动态（逻辑增长 + 免根分子钟）\n",
       f"- 逻辑增长拟合（年份点 n≥5, 加权）: **K={popt[0]:.3f}±{perr[0]:.3f}**, "
       f"r={popt[1]:.3f}±{perr[1]:.3f}/年, **t50={popt[2]:.1f}±{perr[2]:.1f}**, R²={r2:.3f}\n"]
if popt2 is not None:
    rep.append(f"- 敏感性（剔除 China/Taiwan, 点 n≥5）: K={popt2[0]:.3f}±{perr2[0]:.3f}, "
               f"r={popt2[1]:.3f}±{perr2[1]:.3f}, t50={popt2[2]:.1f}±{perr2[2]:.1f}\n")
rep.append(f"- 免根成对分子钟（ST3, Δyear>0）: n={len(dy):,} 对, "
           f"斜率 **{slope:.2f} SNP/年（两谱系合计）→ 每谱系 {slope/2:.2f} SNP/年**; r={r_pear:.3f}; "
           f"截距 {intercept:.1f} SNP（同代基线差）\n")
rep.append(f"- 成对距离为 0 的比例: {(dd==0).mean():.3f}（同波紧克隆）\n")
rep.append("- 分箱中位（Δyear→SNP）: " +
           "; ".join(f"{b:.0f}y:{v:.0f}" for b, v in zip(bx, by_)) + "\n")
open(os.path.join(RES, "MY9_流行病学动态.md"), "w", encoding="utf-8").write("".join(rep))
print("".join(rep))

# ---------- Fig4 ----------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
xs = np.linspace(1995, 2026, 200)
ax.plot(xs, logistic(xs, *popt), color="#c0392b", lw=2,
        label=f"logistic fit: K={popt[0]:.2f}, r={popt[1]:.2f}/yr, t50={popt[2]:.1f}")
ax.scatter(yrs, fr, s=wts * 4, color="#2c3e50", zorder=3,
           label="observed (size ∝ n typed)")
for y, f, w in zip(yrs, fr, wts):
    ax.annotate(f"{int(w)}", (y, f), textcoords="offset points", xytext=(0, 7),
                ha="center", fontsize=7, color="#555")
ax.set_xlabel("Collection year"); ax.set_ylabel("Macrolide resistance frequency")
ax.set_ylim(-0.05, 1.08); ax.set_title("A  Resistance prevalence and logistic growth")
ax.legend(fontsize=8, loc="lower right", markerscale=0.22, handletextpad=0.6,
          borderpad=0.5, labelspacing=0.8); ax.grid(alpha=0.3)

ax = axes[1]
ax.scatter(dy, dd, s=3, alpha=0.05, color="#7f8c8d")
ax.plot(bx, by_, "-o", color="#c0392b", ms=6, lw=2, label="binned median")
xs2 = np.linspace(0, dy.max(), 50)
ax.plot(xs2, slope * xs2 + intercept, "--", color="#2980b9", lw=2,
        label=f"pairwise regression: {slope:.2f} SNP/yr (both lineages, r={r_pear:.2f})")
ax.set_xlabel("|Δ sampling year|"); ax.set_ylabel("Pairwise SNP distance")
ax.set_title("B  Root-free pairwise molecular clock (ST3 clone)")
ax.legend(fontsize=8, markerscale=1.0, borderpad=0.5, labelspacing=0.8); ax.grid(alpha=0.3)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"MY9_Fig4_时间动态.{ext}"), dpi=200)
print("Fig4 ->", FIG)
