# -*- coding: utf-8 -*-
"""01b_ena_manifest.py — Build the download manifest using assembly_name from the ENA Portal
Input:  ../ena_assembly_meta.tsv (accession,assembly_name,assembly_level)
        the genome metadata CSV for the inclusion pool
Output: the download manifest CSV under the repository root
        (assembly_accession (with version), ena_acc (without version), trio_ok, included,
        ena_fasta_url, ncbi_fna_gz_url, filename)
No NCBI connectivity required (the NCBI URL is a fallback column, usable once rate limiting is lifted).
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
            # GCF (RefSeq) / GCA (GenBank) pairing: RefSeq and GenBank share the numeric part
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
