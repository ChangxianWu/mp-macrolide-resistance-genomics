# -*- coding: utf-8 -*-
"""02_download_fasta.py — Batch download of FASTA files (idempotent and resumable, multithreaded,
exponential backoff)
Usage: python MY9_02_download_fasta.py [--included-only] [--workers 6]
Output: ../data/fasta/<acc>.fna.gz + ../download_log.tsv
"""
import argparse, concurrent.futures as cf, csv, gzip, os, sys, time, urllib.request, urllib.error

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(BASE, "data", "fasta")
LOG = os.path.join(BASE, "download_log.tsv")

def fetch(url, dest, tries=4):
    for i in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MY9-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) < 2 or data[0] != 0x1F or data[1] != 0x8B:
                raise RuntimeError("not a gzip file")
            with open(dest, "wb") as f:
                f.write(data)
            return len(data)
        except Exception as e:
            wait = 4 * i
            print(f"[retry {i}] {os.path.basename(dest)}: {type(e).__name__}; sleep {wait}s", flush=True)
            time.sleep(wait)
    return -1

def work(row):
    dest = os.path.join(OUTDIR, row["filename"])
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return (row["assembly_accession"], "skip", os.path.getsize(dest))
    n = fetch(row["url"], dest)
    return (row["assembly_accession"], "ok" if n > 0 else "fail", n)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--included-only", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []
    with open(os.path.join(BASE, "MY9_下载清单.csv"), encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if a.included_only and r["included"] != "1":
                continue
            rows.append(r)
    print(f"to download: {len(rows)} (workers={a.workers})", flush=True)
    results = []
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for res in ex.map(work, rows):
            results.append(res)
            if len(results) % 25 == 0:
                ok = sum(1 for _, s, _ in results if s in ("ok", "skip"))
                print(f"progress {ok}/{len(results)}", flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        for acc, st, n in results:
            f.write(f"{acc}\t{st}\t{n}\n")
    ok = sum(1 for _, s, _ in results if s in ("ok", "skip"))
    print(f"DONE ok/skip={ok} fail={len(results)-ok}")

if __name__ == "__main__":
    main()
