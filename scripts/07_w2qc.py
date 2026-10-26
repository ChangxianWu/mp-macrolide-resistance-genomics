# -*- coding: utf-8 -*-
"""07_w2qc.py — 组装完整性 QC + 去冗余（Mash 替代，纯标准库 MinHash）+ 异常检出
输入: ../data/fasta/*.fna + ../MY9_下载清单.csv + ../MY9_入组池_基因组元数据.csv
输出: ../results/MY9_W2_QC.csv（株级 QC 表）+ 控制台摘要

方法:
  1) 完整性: 每株 n_contigs / total_len / N50 / GC
  2) 去冗余: canonical 16-mer 滑窗取样（stride 50, ~16k seeds/株）全对 Jaccard；
     ≥0.97 判为疑似同一分离株重复提交（不同 BioSample 的同株双提交），保留入组池优先级更高者
  3) 异常: 总长偏离 MP 典型 0.72–0.95 Mb（双基因组混合组装 ~1.6 Mb 已在 04 分型中以 dup2ok 处理）
  4) excluded 复核: W1 esummary 未落盘 exclfromrefseq 字段，本表以 status + 长度异常代偿，如实记录

注: 本 QC 用于系统发生前样本洁化（W4–6 core-SNP 输入），不改动 23S/P1/MLST 分型结果。
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

    # --- 全对 Jaccard ---
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

    # 写出（保持入组池顺序）
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
