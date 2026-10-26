# -*- coding: utf-8 -*-
"""19_treetime.py — TreeTime 正式定年（ST3 克隆子树）
输入: results/coresnp/treetime/{ST3_subtree.nwk, ST3.aln.fasta, ST3_dates.csv}
      results/MY9_分型矩阵.tsv（国家/年份, 用于分组 MRCA）
网格: 定根方式 × 共祖先先验（固定 Tc; coalescent='opt' 在本机实测卡死, 弃用）
输出: results/coresnp/treetime/tt_summary.tsv + tt_run.log
      results/coresnp/treetime/timetree_<setting>.nexus（时间树）
      results/coresnp/MY9_TreeTime定年报告.md（由 MY9_19b 汇总生成）

口径:
- clock_rate 单位 subs/site/year; × core 位点数 742,308 = SNP/genome/year（每谱系）
- R² = TreeTime date2dist.r_val（= Treg.explained_variance, 根-末梢回归解释方差）
- 90% 后验区间仅对关键节点调用 get_max_posterior_region（全节点扫描代价过高）
本机 venv: MY9/.venv-align (treetime 0.11.4)
"""
import csv, os, sys, time
from treetime import TreeTime
from Bio import Phylo

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
TTDIR = os.path.join(RESDIR, "coresnp", "treetime")
CORE_SITES = 742308  # v2 面板 654 株 core 权威值
LOG = open(os.path.join(TTDIR, "tt_run.log"), "w", encoding="utf-8")

def log(msg):
    print(msg, flush=True)
    LOG.write(msg + "\n"); LOG.flush()

with open(os.path.join(TTDIR, "ST3_dates.csv"), encoding="utf-8-sig") as f:
    dates = {r["name"]: float(r["date"]) for r in csv.DictReader(f)}
log(f"dated tips: {len(dates)}")

mat = {}
mpath = os.path.join(RESDIR, "MY9_分型矩阵.tsv")
if os.path.exists(mpath):
    with open(mpath, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            mat[r["ena_acc"]] = r

SETTINGS = [
    ("Tc0.03_best",   dict(root="best",   Tc=0.03)),
    ("Tc0.03_min_dev", dict(root="min_dev", Tc=0.03)),
    ("Tc0.03_oldest", dict(root="oldest", Tc=0.03)),
    ("noprior_best",  dict(root="best",   Tc=1e-6)),
    ("noprior_oldest", dict(root="oldest", Tc=1e-6)),
]

GROUPS = {
    "wave2023plus": lambda n: dates.get(n, 0) >= 2023.0,
    "wave2019_20":  lambda n: 2019.0 <= dates.get(n, 0) < 2021.0,
    "China2023":    lambda n: (mat.get(n, {}).get("country", "") == "China"
                               and dates.get(n, 0) >= 2023.0),
}

rows = []
for name, kw in SETTINGS:
    t0 = time.time()
    try:
        tt = TreeTime(gtr="Jukes-Cantor",
                      tree=os.path.join(TTDIR, "ST3_subtree.nwk"),
                      aln=os.path.join(TTDIR, "ST3.aln.fasta"),
                      dates=dates, verbose=0)
        log(f"[{name}] loaded ({time.time()-t0:.0f}s)")
        tt.run(root=kw["root"], infer_gtr=False, resolve_polytomies=True,
               coalescent=None, Tc=kw["Tc"], max_iter=2)
        rate = float(tt.date2dist.clock_rate)
        r2 = float(getattr(tt.date2dist, "r_val", float("nan")))
        root = tt.tree.root
        rec = {"setting": name, "root": kw["root"], "Tc": kw["Tc"],
               "rate_subs_site_yr": f"{rate:.4e}",
               "rate_snp_genome_yr": f"{rate * CORE_SITES:.3f}",
               "R2": f"{r2:.4f}", "root_numdate": f"{root.numdate:.0f}"}
        tips = root.get_terminals()
        for gname, fn in GROUPS.items():
            sel = [t for t in tips if fn(t.name)]
            if len(sel) < 2:
                rec[gname] = "n/a"; continue
            try:
                clade = root.common_ancestor(sel)
                ci = tt.get_max_posterior_region(clade, fraction=0.9)
                rec[gname] = f"{clade.numdate:.0f} [{ci[0]:.0f},{ci[1]:.0f}] n={len(sel)}"
            except Exception as e:
                rec[gname] = f"err:{type(e).__name__}"
        rec["secs"] = f"{time.time()-t0:.0f}"
        rows.append(rec)
        log(f"[{name}] rate={rec['rate_snp_genome_yr']} SNP/yr R2={rec['R2']} "
            f"root={rec['root_numdate']} 2023MRCA={rec.get('wave2023plus')} ({rec['secs']}s)")
        if name in ("Tc0.03_best", "noprior_oldest"):
            Phylo.write(tt.tree, os.path.join(TTDIR, f"timetree_{name}.nexus"), "nexus")
    except Exception as e:
        rows.append({"setting": name, "root": kw["root"], "Tc": kw["Tc"],
                     "rate_subs_site_yr": "err", "rate_snp_genome_yr": "err",
                     "R2": "err", "root_numdate": "err",
                     "wave2023plus": f"err:{type(e).__name__}:{e}", "secs": f"{time.time()-t0:.0f}"})
        log(f"[{name}] FAILED: {type(e).__name__}: {e}")

cols = ["setting", "root", "Tc", "rate_subs_site_yr", "rate_snp_genome_yr",
        "R2", "root_numdate", "wave2023plus", "wave2019_20", "China2023", "secs"]
with open(os.path.join(TTDIR, "tt_summary.tsv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
log("summary -> " + os.path.join(TTDIR, "tt_summary.tsv"))
LOG.close()
