# -*- coding: utf-8 -*-
"""02b_download_ena.py — ENA Browser API 批量下载 FASTA（幂等可续传）
用法: python MY9_02b_download_ena.py [--included-only] [--workers 5]
输出: ../data/fasta/<ena_acc>.fna + ../download_log.tsv
"""
import argparse, concurrent.futures as cf, csv, os, time, urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(BASE, "data", "fasta")
LOG = os.path.join(BASE, "download_log.tsv")

def fetch(url, dest, tries=4):
    for i in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MY9-pipeline/1.0"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) < 1000 or not data.startswith(b">"):
                raise RuntimeError("not a fasta")
            with open(dest, "wb") as f:
                f.write(data)
            return len(data)
        except Exception as e:
            wait = 5 * i
            print(f"[retry {i}] {os.path.basename(dest)}: {type(e).__name__}; sleep {wait}s", flush=True)
            time.sleep(wait)
    return -1

def work(row):
    dest = os.path.join(OUTDIR, row["filename"])
    if os.path.exists(dest) and os.path.getsize(dest) > 100000:
        return (row["ena_acc"], "skip", os.path.getsize(dest))
    n = fetch(row["ena_fasta_url"], dest)
    return (row["ena_acc"], "ok" if n > 0 else "fail", n)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--included-only", action="store_true")
    ap.add_argument("--workers", type=int, default=5)
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
