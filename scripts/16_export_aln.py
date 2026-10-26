# -*- coding: utf-8 -*-
"""16_export_aln.py — 导出 TreeTime/BEAST 输入（variable-sites 比对 + 日期表）
输入: results/coresnp/alleles/*.npy + MY9_分型矩阵.tsv + reference（M129 序列）
输出: results/coresnp/MY9_variable_sites.aln.fasta（panel assemblies + M129 参考, 核心可变位点）
      results/coresnp/MY9_dates.csv（tip,year）
      results/coresnp/MY9_variable_sites_meta.tsv（位点清单: refpos, ref, alleles）
注: 本机 treetime 不可装（镜像无兼容发行）时, 该导出供任何可联网环境直接运行:
    treetime --tree MY9_coreSNP_NJ.nwk --aln MY9_variable_sites.aln.fasta --dates MY9_dates.csv
"""
import csv, glob, os, re
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
BIOT = "ACGT"

accs = [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(CS, "alleles", "*.npy")))]
M = np.stack([np.load(os.path.join(CS, "alleles", a + ".npy")) for a in accs])
n, L = M.shape

core_mask = (M < 4).sum(axis=0) >= int(0.95 * n)
C = M[:, core_mask]
Vk = C < 4
n_valid = Vk.sum(axis=0)
n_alleles = np.zeros(C.shape[1], dtype=np.int16)
for b in range(4):
    n_alleles += (C == b).sum(axis=0) > 0
var_mask = n_alleles >= 2
print(f"core={int(core_mask.sum())}, variable sites={int(var_mask.sum())}", flush=True)

def embl_seq(path):
    embl = open(path, encoding="utf-8").read()
    body = embl.split("SQ   ", 1)[1].split("//", 1)[0]
    return "".join(re.findall(r"[A-Za-z]+", " ".join(body.splitlines()[1:]))).upper()

ref = embl_seq(os.path.join(BASE, "reference", "U00089.2.embl.cache.txt"))
ref_codes = np.full(L, 4, dtype=np.uint8)
ci = np.nonzero(core_mask)[0]
ref_codes[ci] = np.array([BIOT.index(b) for b in (ref[p] for p in ci)], dtype=np.uint8)

CV = C[:, var_mask]
RV = ref_codes[ci][var_mask]
pos = ci[var_mask]

with open(os.path.join(CS, "MY9_variable_sites.aln.fasta"), "w") as f:
    def wseq(name, arr):
        s = "".join(BIOT[b] if b < 4 else "-" for b in arr)
        f.write(f">{name}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i+80] + "\n")
    wseq("M129_reference", RV)
    for i, a in enumerate(accs):
        wseq(a, CV[i])

mat = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}
with open(os.path.join(CS, "MY9_dates.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "date"])
    for i, a in enumerate(accs):
        y = mat.get(a, {}).get("year", "")
        if y.isdigit():
            w.writerow([a, y])

with open(os.path.join(CS, "MY9_variable_sites_meta.tsv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["refpos_1based", "ref", "n_alleles"])
    for p, b, na in zip(pos, RV, n_alleles[var_mask]):
        w.writerow([p + 1, BIOT[b], na])

print(f"exported: alignment {CV.shape[0]+1} seqs x {CV.shape[1]} variable sites -> results/coresnp/")
