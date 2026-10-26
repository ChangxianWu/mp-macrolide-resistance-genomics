# -*- coding: utf-8 -*-
"""09_mlst.py — MLST 分型（pubMLST pubmlst_mpneumoniae_seqdef 方案 1，8 位点）
输入: ../data/fasta/*.fna + pubMLST API（等位基因 FASTA + profiles_csv，缓存至 reference/）
输出: ../results/MY9_MLST.tsv + 控制台摘要

方法:
  1) 每位点拉取全部等位基因序列（7± 条/位点，~470-680 nt）+ profiles_csv（ST ← 8 位点等位组合）
  2) 基因组内定位: 等位基因全长精确匹配（正/反链）; 未命中 → 20-mer 共享最高者记 nearest（标注 mismatch）
  3) ST: 8 位点组合精确查表; 组合未收录 → novel_ST(near_ST_x)
局限: pubMLST 未鉴权时仅含 ≤2024-12-31 提交的 profile（对定义无影响）; 新等位基因需 W4+ 更新。
"""
import csv, glob, json, os, sys, urllib.request, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFDIR = os.path.join(BASE, "reference")
RESDIR = os.path.join(BASE, "results")
DB = "pubmlst_mpneumoniae_seqdef"
LOCI = ["ppa", "pgm", "gyrB", "gmk", "glyA", "atpA", "arcC", "adk"]
CACHE = os.path.join(REFDIR, "mlst_cache.json")

def get(u, binary=False, tries=5):
    for i in range(1, tries + 1):
        try:
            req = urllib.request.Request(u, headers={"Accept": "text/plain"})
            d = urllib.request.urlopen(req, timeout=120).read()
            return d.decode("utf-8", "replace") if not binary else d
        except Exception as e:
            print(f"[retry {i}] {type(e).__name__} {u[-40:]}; sleep {5*i}s", flush=True)
            import time; time.sleep(5 * i)
    raise RuntimeError(f"FAILED {u}")

# --- 1) 等位基因 + profiles（缓存）---
if os.path.exists(CACHE):
    cache = json.load(open(CACHE))
    print("[cache] mlst_cache.json", flush=True)
else:
    cache = {"alleles": {}, "profiles": []}
    for locus in LOCI:
        t = get(f"https://rest.pubmlst.org/db/{DB}/loci/{locus}/alleles_fasta")
        alleles = {}
        cur = None
        for ln in t.splitlines():
            if ln.startswith(">"):
                cur = ln[1:].strip().split("_")[-1]
                alleles[cur] = []
            elif cur:
                alleles[cur].append(ln.strip())
        cache["alleles"][locus] = {k: "".join(v).upper() for k, v in alleles.items()}
        print(f"{locus}: {len(cache['alleles'][locus])} alleles", flush=True)
    p = get(f"https://rest.pubmlst.org/db/{DB}/schemes/1/profiles_csv")
    rows = [l.split("\t") for l in p.strip().splitlines() if l.strip()]
    cache["profiles"] = rows
    json.dump(cache, open(CACHE, "w"))
    print("[saved]", CACHE, flush=True)

alleles = cache["alleles"]
prof_rows = cache["profiles"]
header = prof_rows[0]
st_map = {}
for r in prof_rows[1:]:
    r = (r + [""] * len(header))[:len(header)]   # 尾列缺失补空
    combo = tuple(r[header.index(loc)] for loc in LOCI)
    st_map[combo] = (r[0], r[header.index("clonal_complex")])
print(f"profiles: {len(st_map)} STs", flush=True)

RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s):
    return s.translate(RC)[::-1]

def kmers(seq, k=20):
    return {seq[i:i+k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i+k]}

def call_locus(genome, kg, allele_seqs):
    """返回 (allele, note)。精确匹配优先，否则最近邻(20-mer)。"""
    best_exact, best_near, best_j = None, None, -1
    for num, seq in allele_seqs.items():
        if seq in genome or rc(seq) in genome:
            if best_exact is None:
                best_exact = num
        kg_a = kmers(seq)
        j = len(kg & kg_a) / max(len(kg_a), 1)
        if j > best_j:
            best_j, best_near = j, num
    if best_exact is not None:
        return best_exact, "exact"
    return f"near_{best_near}", f"nearest_j={best_j:.2f}"

# --- 2) 逐株分型 ---
def fseq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

pool = {}
with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        pool[r["assembly_accession"]] = r
meta = {}
with open(os.path.join(BASE, "MY9_下载清单.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rec = pool.get(r["assembly_accession"])
        if rec:
            meta[r["ena_acc"]] = rec

files = sorted(glob.glob(os.path.join(BASE, "data", "fasta", "*.fna")))
if not files:
    sys.exit("no fasta")
out = os.path.join(RESDIR, "MY9_MLST.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ena_acc", "country", "year"] + LOCI + ["ST", "clonal_complex", "note"])
    for k, path in enumerate(files, 1):
        acc = os.path.basename(path).replace(".fna", "")
        g = fseq(path)
        kg = kmers(g)
        calls, notes = [], []
        for locus in LOCI:
            num, note = call_locus(g, kg, alleles[locus])
            calls.append(num)
            if note != "exact":
                notes.append(f"{locus}:{note}")
        combo = tuple(c.replace("near_", "") for c in calls)
        exact_combo = all(not c.startswith("near_") for c in calls)
        if exact_combo and combo in st_map:
            st, cc = st_map[combo]
        elif exact_combo:
            near = min(st_map.keys(), key=lambda c2: sum(a != b for a, b in zip(c2, combo)))
            st, cc = f"novel(near {st_map[near][0]})", st_map[near][1]
        else:
            near = min(st_map.keys(), key=lambda c2: sum(a != b for a, b in zip(c2, combo)))
            st, cc = f"uncertain(near {st_map[near][0]})", st_map[near][1]
        m = meta.get(acc, {})
        w.writerow([acc, m.get("country", ""), m.get("year", "")] + calls + [st, cc, ";".join(notes)])
        if k % 50 == 0:
            print(f"typed {k}/{len(files)}", flush=True)

# --- 3) 摘要 ---
rows = list(csv.DictReader(open(out), delimiter="\t"))
c = collections.Counter(r["ST"] for r in rows)
print("--- ST distribution ---")
for kk, v in c.most_common(15):
    print(f"{kk}: {v}")
ccy = collections.Counter((r["country"].split(":")[0].strip(), r["ST"]) for r in rows)
print("--- top country×ST ---")
for (nat, st), v in ccy.most_common(12):
    print(f"{nat} {st}: {v}")
print(f"written {out}")
