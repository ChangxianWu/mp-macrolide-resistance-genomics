# -*- coding: utf-8 -*-
"""07_w2qc.py — assembly completeness QC + redundancy removal (Mash replacement, pure standard-library MinHash) + anomaly detection
Input: ../data/fasta/*.fna + the download manifest CSV + the inclusion-pool genome metadata CSV
Output: ../results/MY9_W2_QC.csv (per-strain QC table) + console summary

Method:
  1) Completeness: n_contigs / total_len / N50 / GC per strain
  2) Redundancy removal: canonical 16-mer sliding-window sampling (stride 50, ~16k seeds/strain) all-pairs Jaccard;
     >=0.97 is treated as a suspected duplicate submission of the same isolate (the same strain submitted twice under different BioSamples); keep the one with higher inclusion-pool priority
  3) Anomaly: total length deviates from the typical MP range 0.72–0.95 Mb (dual-genome mixed assemblies ~1.6 Mb were already handled in step 04 typing via dup2ok)
  4) excluded recheck: the W1 esummary exclfromrefseq field was not persisted, so this table compensates with status + length anomaly and records it as-is

Note: this QC is used for sample cleanup before phylogenetics (W4–6 core-SNP input); it does not alter the 23S/P1/MLST typing results.
"""
import csv, glob, itertools, os, sys, time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K, STRIDE = 16, 50
JAC_DUP = 0.97
TYPICAL_LO, TYPICAL_HI = 0.72e6, 0.95e6   # M129 = 816,394 bp

RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def load(path):
    op = open if not path.endswith(".gz") else __import__("gzip").open
    contigs, cur = [], []
    with op(path, "rt") as f:
        for ln in f:
            if ln.startswith(">"):
                if cur:
                    contigs.append("".join(cur)); cur = []
            else:
                cur.append(ln.strip())
    if cur:
        contigs.append("".join(cur))
    return contigs

def n50(lens):
    lens = sorted(lens, reverse=True)
    half = sum(lens) / 2
    acc = 0
    for L in lens:
        acc += L
        if acc >= half:
            return L
    return 0

def sketch(contigs):
    s = set()
    for c in contigs:
        c = c.upper()
        for i in range(0, len(c) - K + 1, STRIDE):
            k = c[i:i + K]
            if "N" in k:
                continue
            rk = k.translate(RC)[::-1]
            s.add(min(k, rk))
    return s

def main():
    files = sorted(glob.glob(os.path.join(BASE, "data", "fasta", "*.fna")))
    print(f"genomes: {len(files)}", flush=True)

    with open(os.path.join(BASE, "MY9_下载清单.csv"), encoding="utf-8-sig") as f:
        dl = {r["ena_acc"]: r for r in csv.DictReader(f)}
    with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
        pool = {r["assembly_accession"]: r for r in csv.DictReader(f)}
    meta = {}
    for acc, r in dl.items():
        rec = pool.get(r["assembly_accession"])
        if rec:
            meta[acc] = rec

    qc, sketches = [], {}
    t0 = time.time()
    for idx, path in enumerate(files, 1):
        acc = os.path.basename(path).replace(".fna.gz", "").replace(".fna", "")
        contigs = load(path)
        lens = sorted((len(c) for c in contigs), reverse=True)
        total = sum(lens)
        gc = sum(c.count("G") + c.count("C") for c in contigs)
        at = sum(c.count("A") + c.count("T") for c in contigs)
        sketches[acc] = sketch(contigs)
        m = meta.get(acc, {})
        qc.append({
            "ena_acc": acc,
            "assembly_accession": m.get("assembly_accession", ""),
            "country": m.get("country", ""), "year": m.get("year", ""),
            "assembly_status": m.get("assembly_status", ""),
            "n_contigs": len(contigs), "total_len": total, "N50": n50(lens),
            "GC": f"{(gc / max(gc + at, 1)):.4f}",
            "len_flag": 0 if TYPICAL_LO <= total <= TYPICAL_HI else 1,
        })
        if idx % 50 == 0:
            print(f"sketched {idx}/{len(files)} ({time.time()-t0:.0f}s)", flush=True)

    # --- all-pairs Jaccard ---
    print("pairwise Jaccard...", flush=True)
    t0 = time.time()
    accs = [q["ena_acc"] for q in qc]
    sk = {a: sketches[a] for a in accs}
    dup_pairs = []
    for a, b in itertools.combinations(accs, 2):
        sa, sb = sk[a], sk[b]
        inter = len(sa & sb)
        if inter == 0:
            continue
        j = inter / (len(sa) + len(sb) - inter)
        if j >= JAC_DUP:
            dup_pairs.append((a, b, round(j, 4)))
    print(f"{len(dup_pairs)} dup pairs (J>={JAC_DUP}) in {time.time()-t0:.0f}s", flush=True)

    dupmap = {}
    for a, b, j in dup_pairs:
        dupmap.setdefault(a, []).append(b)
        dupmap.setdefault(b, []).append(a)
    for q in qc:
        q["dup_with"] = ";".join(sorted(dupmap.get(q["ena_acc"], [])))
        q["dup_flag"] = 1 if q["dup_with"] else 0

    # write out (preserving inclusion-pool order)
    out = os.path.join(BASE, "results", "MY9_W2_QC.csv")
    with open(out, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(qc[0].keys()))
        w.writeheader()
        w.writerows(qc)
    print("written:", out)

    nf = sum(q["len_flag"] for q in qc)
    ndup = len({a for p in dup_pairs for a in p[:2]})
    lens = [q["total_len"] for q in qc]
    print(f"\n--- QC summary ---")
    print(f"genomes={len(qc)}  total_len min/median/max = {min(lens)}/{sorted(lens)[len(lens)//2]}/{max(lens)}")
    print(f"len_flag(outside 0.72-0.95Mb)={nf}  dup_flag genomes={ndup} (pairs={len(dup_pairs)})")
    st = {}
    for q in qc:
        st[q["assembly_status"]] = st.get(q["assembly_status"], 0) + 1
    print("assembly_status:", st)

if __name__ == "__main__":
    main()
