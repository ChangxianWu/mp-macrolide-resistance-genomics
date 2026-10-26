# -*- coding: utf-8 -*-
"""01_build_manifest.py — Build the FASTA download manifest
Input:  ../assembly_raw.tsv + the genome metadata CSV for the inclusion pool
Output: the download manifest CSV under the repository root
        (assembly_accession, biosample, trio_ok, url, filename)
Deps:   standard library only; NCBI esummary is queried to fill in assembly_name
        (required to construct the FTP path)
"""
import csv, json, sys, time, urllib.request, urllib.error, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def fetch_json(url, tries=5):
    for i in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MY9-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            wait = 5 * i
            print(f"[retry {i}] {type(e).__name__}: {e}; sleep {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"FAILED after {tries} tries: {url}")

def main():
    asm = {}
    with open(os.path.join(BASE, "assembly_raw.tsv"), encoding="utf-8") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 6:
                asm[c[0]] = c[1]  # uid -> accession
    trio = {}
    with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            trio[row["assembly_accession"]] = (row["included"], row["trio_ok"])

    uids = list(asm.keys())
    names = {}
    for i in range(0, len(uids), 150):
        chunk = ",".join(uids[i:i + 150])
        r = fetch_json(f"{EUTILS}/esummary.fcgi?db=assembly&retmode=json&id={chunk}")
        for uid in r["result"]["uids"]:
            names[uid] = r["result"][uid].get("assemblyname", "")
        print(f"names {min(i+150, len(uids))}/{len(uids)}", flush=True)
        time.sleep(0.5)

    out = os.path.join(BASE, "MY9_下载清单.csv")
    n = 0
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["assembly_accession", "biosample", "trio_ok", "included", "url", "filename"])
        for uid, acc in asm.items():
            if acc not in trio:
                continue
            name = names.get(uid, "")
            if not name:
                print(f"WARN no assemblyname for {acc}; skipped")
                continue
            digits = acc.split("_")[1].split(".")[0]
            prefix = acc[:3]
            url = (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/"
                   f"{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{acc}_{name}/{acc}_{name}_genomic.fna.gz")
            w.writerow([acc, "", trio[acc][1], trio[acc][0], url, acc + ".fna.gz"])
            n += 1
    print(f"manifest rows={n} -> {out}")

if __name__ == "__main__":
    main()
