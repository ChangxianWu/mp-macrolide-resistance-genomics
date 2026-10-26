# -*- coding: utf-8 -*-
"""01b_ena_manifest.py — 用 ENA Portal 的 assembly_name 构建下载清单
输入: ../ena_assembly_meta.tsv (accession,assembly_name,assembly_level)
      ../MY9_入组池_基因组元数据.csv
输出: ../MY9_下载清单.csv  (assembly_accession(带版本), ena_acc(无版本), trio_ok, included,
                            ena_fasta_url, ncbi_fna_gz_url, filename)
无需 NCBI 在线（NCBI URL 备用列，待其解除限流后可用）。
"""
import csv, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

names = {}
with open(os.path.join(BASE, "ena_assembly_meta.tsv"), encoding="utf-8") as f:
    next(f)
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            names[c[0]] = c[1]

pool = {}
with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        key = r["assembly_accession"].split(".")[0]
        pool[key] = r

out = os.path.join(BASE, "MY9_下载清单.csv")
missing = 0
n = 0
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["assembly_accession", "ena_acc", "biosample", "trio_ok", "included",
                "ena_fasta_url", "ncbi_fna_gz_url", "filename"])
    for key, r in pool.items():
        ver = r["assembly_accession"]
        name = names.get(key)
        if not name:
            # GCF(GenBank GCA 配对): RefSeq 与 GenBank 共享数字段
            if key.startswith("GCF"):
                gca_key = "GCA" + key[3:]
                name = names.get(gca_key)
                if name:
                    key = gca_key
        if not name:
            missing += 1
            continue
        digits = key.split("_")[1]
        ncbi_url = (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{ver[:3]}/"
                    f"{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{ver}_{name}/{ver}_{name}_genomic.fna.gz")
        w.writerow([ver, key, r["biosample"], r["trio_ok"], r["included"],
                    f"https://www.ebi.ac.uk/ena/browser/api/fasta/{key}", ncbi_url, key + ".fna"])
        n += 1
print(f"manifest rows={n}  missing_in_ENA={missing}  -> {out}")
