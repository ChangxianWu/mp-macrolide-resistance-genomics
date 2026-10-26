# -*- coding: utf-8 -*-
"""14_coresnp.py — 本地 core-SNP stage 1：参考锚链法逐株变异 calling（自研, 替代 snippy）

背景: snippy 不可用（WSL 启动挂死 / GitHub 封锁 / minimap2·mappy 均无 Windows 发行）——
      本机内存 + numpy/Biopython 足以支撑 0.8 Mb 基因组 × 646 assemblies 的 contig-vs-reference 比对。
方法（MUMmer 同款策略）:
  1) 参考唯一 25-mer 索引（正向 2-bit 打包, 参考内仅出现一次）
  2) 每株: 正/反向两组 25-mer 打包查索引（qry 内任一方向多拷贝者弃）→ (ref_pos, qry_pos, strand) 锚
  3) 同链锚链化（ref 同向递增且 qry 沿链方向单调）→ 链内相邻锚空隙:
       等长 → numpy 直比（SNP 直出）; 不等长且两侧 ≤3000 nt → Biopython NW（indel ±8 nt 掩蔽）;
       任一侧 >3000 nt → 未覆盖（结构差异/组装断点/重复区）
  4) 产物: 每株 816,394 位点等位向量（0=A,1=C,2=G,3=T,4=未覆盖/掩蔽）→ results/coresnp/alleles/<acc>.npy
     + results/coresnp/calling_summary.tsv（断点续跑: 已完成株跳过）
"""
import glob, os, re, time

import numpy as np
from Bio import Align

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results", "coresnp")
ALDIR = os.path.join(RESDIR, "alleles")
K = 25
GAP_ALIGN_MAX = 3000
INDEL_FLANK = 8
BIOT = "ACGT"
BIDX = {b: i for i, b in enumerate(BIOT)}

TBL = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate("ACGT"):
    TBL[ord(b)] = i
    TBL[ord(b.lower())] = i

def codes_from_embl(path):
    embl = open(path, encoding="utf-8").read()
    body = embl.split("SQ   ", 1)[1].split("//", 1)[0]
    seq = "".join(re.findall(r"[A-Za-z]+", " ".join(body.splitlines()[1:]))).upper()
    return TBL[np.frombuffer(seq.encode(), dtype=np.uint8)]

def codes_of_fasta(path):
    seq = []
    for ln in open(path):
        if not ln.startswith(">"):
            seq.append(ln.strip())
    return TBL[np.frombuffer("".join(seq).encode(), dtype=np.uint8)]

def pack_windows(codes, k=K):
    n = len(codes)
    W = n - k + 1
    if W <= 0:
        return (np.zeros(0, dtype=np.uint64),) * 2
    c = codes.astype(np.uint64)
    packs = np.zeros(W, dtype=np.uint64)
    rpacks = np.zeros(W, dtype=np.uint64)
    for j in range(k):
        packs += c[j:j + W] << np.uint64(2 * (k - 1 - j))
        rpacks += (np.uint64(3) - c[j:j + W]) << np.uint64(2 * j)
    return packs, rpacks

def build_ref_index(ref_codes):
    packs, _ = pack_windows(ref_codes)
    vals, idx, cnt = np.unique(packs, return_index=True, return_counts=True)
    keep = cnt == 1
    uvals, upos = vals[keep], idx[keep]
    order = np.argsort(uvals)
    return uvals[order], upos[order], len(packs)

aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                                open_gap_score=-5, extend_gap_score=-2)

def call_snps(ref_codes, qry_codes, uvals, upos):
    n_ref = len(ref_codes)
    allele = np.full(n_ref, 4, dtype=np.uint8)
    stats = {"snp": 0, "gap_aln": 0, "gap_uncovered": 0, "blocks": 0}

    qry_f, qry_r = pack_windows(qry_codes)
    if len(qry_f) == 0:
        return allele, stats
    # 拷贝数感知的唯一性阈值: 双基因组混合组装(len≈2×816kb)中每 k-mer 出现 ≈2 次,
    # 按"仅一次"过滤会把整株杀成零锚 —— 阈值 = 拷贝数（单拷贝基因组退化为 cnt==1, 与原始行为一致）
    copies = max(1, int(round(len(qry_codes) / 816394)))
    max_cnt = copies
    allv = np.concatenate([qry_f, qry_r])
    uv, cnt = np.unique(allv, return_counts=True)
    uniq_vals = uv[cnt <= max_cnt]
    ii = np.clip(np.searchsorted(uvals, uniq_vals), 0, len(uvals) - 1)
    hit = uvals[ii] == uniq_vals
    vals, rp0 = uniq_vals[hit], upos[ii[hit]]

    of = np.argsort(qry_f, kind="stable"); sf = qry_f[of]
    jf = np.clip(np.searchsorted(sf, vals), 0, len(sf) - 1); mf = sf[jf] == vals
    orr = np.argsort(qry_r, kind="stable"); sr = qry_r[orr]
    jr = np.clip(np.searchsorted(sr, vals), 0, len(sr) - 1); mr = sr[jr] == vals

    chains = []
    for strand in (1, 0):
        if strand == 1:
            sel = mf
            qp_all = of[jf][sel]
        else:
            sel = mr & ~mf
            qp_all = orr[jr][sel & mr] if False else orr[jr][sel]
        rp_all = rp0[sel]
        A = np.column_stack([rp_all, qp_all])
        if len(A) == 0:
            continue
        A = A[np.argsort(A[:, 0], kind="stable")]
        sign = 1 if strand == 1 else -1
        cur = [0]
        for i in range(1, len(A)):
            dr = int(A[i, 0] - A[i - 1, 0])
            dq = sign * int(A[i, 1] - A[i - 1, 1])
            if dr > 0 and dq > 0 and dr <= GAP_ALIGN_MAX * 4 and dq <= GAP_ALIGN_MAX * 4:
                cur.append(i)
            else:
                chains.append((strand, A[cur]))
                cur = [i]
        chains.append((strand, A[cur]))

    def fill_anchor(strand, rp, qp):
        # 锚匹配（正/反向）在"对齐链投影"下均为恒等: 反向锚 rc(qry 窗)==ref 窗
        # ⇒ 对齐碱基 = comp(qry[qp+K-1-j]) = ref[rp+j]。等位存对齐链投影碱基, 故正反锚都填参考碱基。
        # （旧版反向填 3-ref 是把"contig 反向"误当"生物学反转"——对反向 contig 基因组会产生 ~50% 假差异）
        allele[rp:rp + K] = ref_codes[rp:rp + K]

    for strand, blk in chains:
        stats["blocks"] += 1
        sign = 1 if strand == 1 else -1
        for rp, qp in blk:
            fill_anchor(strand, int(rp), int(qp))
        for i in range(len(blk) - 1):
            rp1, qp1 = int(blk[i, 0]), int(blk[i, 1])
            rp2, qp2 = int(blk[i + 1, 0]), int(blk[i + 1, 1])
            g_ref = rp2 - (rp1 + K)
            if strand == 1:
                g_qry = qp2 - (qp1 + K)
                qs = qp1 + K
            else:
                g_qry = qp1 - (qp2 + K)
                qs = qp2 + K
            if g_ref < 0 or g_qry < 0 or (g_ref == 0 and g_qry == 0):
                continue
            if g_ref > GAP_ALIGN_MAX or g_qry > GAP_ALIGN_MAX:
                stats["gap_uncovered"] += 1
                continue
            rs = rp1 + K
            rseg = ref_codes[rs:rs + g_ref]
            if strand == 1:
                qraw = qry_codes[qs:qs + g_qry]
            else:
                qraw = (3 - qry_codes[qs:qs + g_qry][::-1]) & 3
            if (rseg >= 4).any() or (qraw >= 4).any():
                stats["gap_uncovered"] += 1
                continue
            if g_ref == g_qry:
                diff = rseg != qraw
                allele[rs:rs + g_ref][~diff] = rseg[~diff]
                dp = np.nonzero(diff)[0]
                for off in dp:
                    allele[rs + off] = qraw[off]
                stats["snp"] += int(diff.sum())
            else:
                if g_ref == 0 or g_qry == 0:
                    # 纯插入/缺失: 掩蔽事件邻域（该区不进 core）
                    lo = max(0, rs - INDEL_FLANK)
                    hi = min(n_ref, rs + max(g_ref, 1) + INDEL_FLANK)
                    allele[lo:hi] = 4
                    continue
                a1, a2 = aligner.align("".join(BIOT[b] for b in rseg),
                                       "".join(BIOT[b] for b in qraw))[0]
                rp = rs
                for ca, cb in zip(a1, a2):
                    if ca != "-" and cb != "-":
                        b = BIDX.get(cb, 4)
                        allele[rp] = b
                        if ca != cb:
                            stats["snp"] += 1
                        rp += 1
                    else:
                        lo = max(0, rp - INDEL_FLANK)
                        hi = min(n_ref, rp + INDEL_FLANK + 1)
                        allele[lo:hi] = 4
                        if ca != "-":
                            rp += 1
                stats["gap_aln"] += 1
    return allele, stats

def main():
    os.makedirs(ALDIR, exist_ok=True)
    ref_codes = codes_from_embl(os.path.join(BASE, "reference", "U00089.2.embl.cache.txt"))
    print(f"reference: {len(ref_codes)} bp", flush=True)
    uvals, upos, W = build_ref_index(ref_codes)
    print(f"ref unique {K}-mer anchors: {len(uvals)}/{W}", flush=True)

    files = sorted(glob.glob(os.path.join(BASE, "data", "fasta", "*.fna")))
    sumpath = os.path.join(RESDIR, "calling_summary.tsv")
    done = set()
    if os.path.exists(sumpath):
        with open(sumpath, encoding="utf-8") as f:
            next(f, None)
            for ln in f:
                done.add(ln.split("\t")[0])
        print(f"[resume] {len(done)} already called", flush=True)
    newfile = not os.path.exists(sumpath)
    out = open(sumpath, "a", encoding="utf-8", newline="")
    if newfile:
        out.write("acc\tn_snp\tgaps_aligned\tgaps_uncovered\tcovered_frac\tsec\n")
    t0 = time.time()
    for i, path in enumerate(files, 1):
        acc = os.path.basename(path).replace(".fna", "")
        outp = os.path.join(ALDIR, acc + ".npy")
        if acc in done and os.path.exists(outp) and os.path.getsize(outp) >= len(ref_codes):
            continue
        tt = time.time()
        allele, st = call_snps(ref_codes, codes_of_fasta(path), uvals, upos)
        np.save(outp, allele)
        cov = float((allele < 4).mean())
        out.write(f"{acc}\t{st['snp']}\t{st['gap_aln']}\t{st['gap_uncovered']}\t{cov:.4f}\t{time.time()-tt:.1f}\n")
        out.flush()
        if i % 25 == 0 or i == len(files):
            print(f"[{i}/{len(files)}] {acc}: snp={st['snp']} cov={cov:.3f} ({time.time()-t0:.0f}s)", flush=True)
    out.close()
    print("stage1 done", flush=True)

if __name__ == "__main__":
    main()
