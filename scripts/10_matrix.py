# -*- coding: utf-8 -*-
"""10_matrix.py — integrated typing matrix + cross tables (23S resistance × P1 type × MLST ST)
Input: the per-genome 23S resistance table, the P1 typing table, the MLST table
      and the W2 QC table under results/
Output: the integrated per-genome typing matrix under results/ (per-strain panorama)
      the P1 type x 23S resistance cross table under results/
      the MLST ST x 23S resistance cross table under results/
"""
import csv, os, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")

def read_tsv(name):
    with open(os.path.join(RESDIR, name), encoding="utf-8") as f:
        return {r["ena_acc"]: r for r in csv.DictReader(f, delimiter="\t")}

t23 = read_tsv("MY9_23S分型.tsv")
p1 = read_tsv("MY9_P1分型.tsv")
mlst = read_tsv("MY9_MLST.tsv")
qc = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_W2_QC.csv"), encoding="utf-8-sig"))}
pcl = {}
_cl_path = os.path.join(RESDIR, "MY9_P1聚类.tsv")
if os.path.exists(_cl_path):
    with open(_cl_path, encoding="utf-8-sig") as f:
        pcl = {r["ena_acc"]: r["P1_cluster"] for r in csv.DictReader(f, delimiter="\t")}

accs = sorted(set(t23) & set(p1) & set(mlst))
out = os.path.join(RESDIR, "MY9_分型矩阵.tsv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ena_acc", "country", "year", "assembly_status", "total_len", "len_flag",
                "site2063", "site2064", "site2617", "resist_23S",
                "P1_type", "P1_cluster", "MLST_ST", "clonal_complex"])
    for a in accs:
        t, p, m, q = t23[a], p1[a], mlst[a], qc.get(a, {})
        resist = "-" if t["site2063"] == "-" else (
            "resistant" if (t["site2063"] != "A" or t["site2064"] != "A" or t["site2617"] != "C") else "WT")
        w.writerow([a, t["country"], t["year"], q.get("assembly_status", ""), q.get("total_len", ""),
                    q.get("len_flag", ""), t["site2063"], t["site2064"], t["site2617"], resist,
                    p["P1_type"], pcl.get(a, ""), m["ST"], m.get("clonal_complex", "")])
print(f"matrix: {len(accs)} genomes -> {out}")

def cross_table(key_field, fname, order=None):
    """Rows = key_field values, columns = resistance status; also includes the 2063G count"""
    tab = collections.defaultdict(collections.Counter)
    for a in accs:
        t, src = t23[a], (p1 if key_field == "P1_type" else mlst)[a]
        k = src[key_field]
        resist = "-" if t["site2063"] == "-" else (
            "resistant" if (t["site2063"] != "A" or t["site2064"] != "A" or t["site2617"] != "C") else "WT")
        tab[k][resist] += 1
        if t["site2063"] == "G":
            tab[k]["2063G"] += 1
    keys = order or sorted(tab, key=lambda k: -sum(v for kk, v in tab[k].items() if kk in ("resistant", "WT")))
    with open(os.path.join(RESDIR, fname), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([key_field, "n_called", "resistant", "WT", "freq_resistant", "n_2063G"])
        for k in keys:
            c = tab[k]
            n = c["resistant"] + c["WT"]
            if n == 0:
                continue
            w.writerow([k, n, c["resistant"], c["WT"], f"{c['resistant']/n:.3f}", c["2063G"]])
    print(f"written {fname}")
    for k in keys[:8]:
        c = tab[k]
        n = c["resistant"] + c["WT"]
        if n:
            print(f"  {k}: n={n} resist={c['resistant']} ({c['resistant']/n:.1%}) 2063G={c['2063G']}")

cross_table("P1_type", "MY9_交叉_P1×23S耐药.tsv", order=["type1", "type2", "ambiguous", "not_called"])
cross_table("ST", "MY9_交叉_ST×23S耐药.tsv")
