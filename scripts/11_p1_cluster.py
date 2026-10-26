# -*- coding: utf-8 -*-
"""11_p1_cluster.py — self-clustering and naming of P1 divergent lineages
Input: ../data/fasta/*.fna + reference/U00089.2.embl.cache.txt (M129 SQ)
      + reference/CP010546.embl.cache.txt (FH SQ) + reference/p1_panel_ncbi.fasta (PHOL panel)
Output: the per-genome P1 cluster assignment table + the cluster summary markdown under ../results/

Method:
  1) Extraction: locate where 20-mers sampled from M129 P1 (every 25 nt, ~195 in total) hit each strain;
     hits >=10 → window = [min(pos)-500, max(pos)+20+500] (P1 neighbourhood, no exact boundaries needed; 20-mers are insensitive to assembly order)
  2) Clustering: all-pairs Jaccard over the window 20-mers, single linkage; threshold taken from a histogram valley (printed by the script, default 0.5)
  3) Naming: clusters containing M129-P1 = type1; containing FH-P1 = type2; containing a PHOL member = phol; the rest = div-A/B/…
  4) Cluster-level statistics: n / country / year range / ST / 23S resistance rate (joined with the integrated typing matrix under results/)
"""
import csv, glob, itertools, json, os, re, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
REFDIR = os.path.join(BASE, "reference")
K = 20
SAMPLE_STEP = 25
MIN_HITS = 10
JAC_T = 0.5

def embl_sq_seq(embl):
    body = embl.split("SQ   ", 1)[1].split("//", 1)[0]
    return "".join(re.findall(r"[A-Za-z]+", " ".join(body.splitlines()[1:]))).upper()

def fseq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def kmers(seq, k=K):
    """canonical k-mer set (strand-invariant)"""
    out = set()
    for i in range(len(seq) - k + 1):
        x = seq[i:i+k]
        if "N" not in x:
            out.add(min(x, x.translate(_RC)[::-1]))
    return out

# --- reference P1 ---
embl129 = open(os.path.join(REFDIR, "U00089.2.embl.cache.txt"), encoding="utf-8").read()
g129 = embl_sq_seq(embl129)
p1_m129 = g129[180857:185741]
samps = [p1_m129[i:i+K] for i in range(0, len(p1_m129) - K + 1, SAMPLE_STEP)]
samps = [s for s in samps if "N" not in s]
print(f"anchor samples from M129 P1: {len(samps)}", flush=True)

embl_fh = open(os.path.join(REFDIR, "CP010546.embl.cache.txt"), encoding="utf-8").read()
fh = embl_sq_seq(embl_fh)   # parse the EMBL SQ block (fseq would treat annotation text as sequence — do not use it)
p1_fh = fh[179323:184228]
assert p1_fh.startswith("ATG"), "FH P1 切片异常"

panel = {}
if os.path.exists(os.path.join(REFDIR, "p1_panel_ncbi.fasta")):
    name = None
    for ln in open(os.path.join(REFDIR, "p1_panel_ncbi.fasta")):
        if ln.startswith(">"):
            name = ln[1:].strip().split()[0]
            panel[name] = []
        elif name:
            panel[name].append(ln.strip())
    panel = {k: "".join(v).upper() for k, v in panel.items()}
print(f"PHOL panel members: {len(panel)}", flush=True)

# --- extraction ---
meta = {}
pool = {}
with open(os.path.join(BASE, "MY9_入组池_基因组元数据.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        pool[r["assembly_accession"]] = r
with open(os.path.join(BASE, "MY9_下载清单.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rec = pool.get(r["assembly_accession"])
        if rec:
            meta[r["ena_acc"]] = rec

files = sorted(glob.glob(os.path.join(BASE, "data", "fasta", "*.fna")))
windows = {}   # key -> (seq, n_anchor_hits)
fails = []
for k, path in enumerate(files, 1):
    acc = os.path.basename(path).replace(".fna", "")
    g = fseq(path)
    gr = g.translate(_RC)[::-1]
    if sum(1 for s in samps if s in gr) > sum(1 for s in samps if s in g):
        g = gr                            # reverse-strand contig genome: extract using the rc strand
    hits = sorted(p for p in (g.find(s) for s in samps) if p >= 0)
    # hit chain: the largest chain with adjacent spacing <=2000 nt (field strains may have extra RepMP copies causing distant scattered hits)
    best_chain = []
    if hits:
        chain = [hits[0]]
        for p in hits[1:]:
            if p - chain[-1] <= 2000:
                chain.append(p)
            else:
                if len(chain) > len(best_chain):
                    best_chain = chain
                chain = [p]
        if len(chain) > len(best_chain):
            best_chain = chain
    if len(best_chain) >= MIN_HITS:
        s0, e0 = best_chain[0], best_chain[-1] + K
        win = g[max(0, s0 - 100): e0 + 100]   # small pad: keep the P1 body itself; flanking conserved sequence would raise the cross-cluster Jaccard background noise
        if 2500 <= len(win) <= 8000:
            windows[acc] = (win, len(best_chain))
        else:
            fails.append((acc, f"window_len={len(win)}"))
    else:
        fails.append((acc, f"chain_hits={len(best_chain)}"))
    if k % 50 == 0:
        print(f"extracted {k}/{len(files)}", flush=True)
print(f"extracted: {len(windows)}, fails: {len(fails)}", flush=True)

# --- clustering (including reference members) ---
members = dict(windows)
members["REF_M129_P1"] = (p1_m129, len(samps))
members["REF_FH_P1"] = (p1_fh, sum(1 for s in samps if s in fh))
for nm, s in panel.items():
    members[f"REF_{nm}"] = (s, sum(1 for s in samps if s in s))
keys = list(members)
sk = {a: kmers(members[a][0]) for a in keys}
print(f"clustering {len(keys)} members...", flush=True)

# --- two-stage classification: nearest reference first (Jaccard + margin), then self-cluster the divergent strains ---
j_m129 = {a: (len(sk[a] & sk["REF_M129_P1"]) / len(sk[a] | sk["REF_M129_P1"])) for a in keys}
j_fh = {a: (len(sk[a] & sk["REF_FH_P1"]) / len(sk[a] | sk["REF_FH_P1"])) for a in keys}
J_REF_MIN = 0.35   # below this = unlike either reference (divergent lineage)
MARGIN = 0.2       # required lead ratio for the winner

assign = {}
div_pool = []
for a in keys:
    if a.startswith("REF_"):
        continue
    j1, j2 = j_m129[a], j_fh[a]
    hi, lo = max(j1, j2), min(j1, j2)
    if hi < J_REF_MIN:
        div_pool.append(a)
    elif (hi - lo) / hi >= MARGIN:
        assign[a] = "type1(M129-like)" if j1 > j2 else "type2(FH-like)"
    else:
        assign[a] = "type-mixed(recombinant?)"

# self-clustering of divergent strains (single linkage)
parent = {a: a for a in div_pool}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for a, b in itertools.combinations(div_pool, 2):
    sa, sb = sk[a], sk[b]
    inter = len(sa & sb)
    if inter and inter / (len(sa) + len(sb) - inter) >= JAC_T:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

groups = collections.defaultdict(list)
for a in div_pool:
    groups[find(a)].append(a)

# cluster-level annotation: mean Jaccard of each cluster against each reference member (evidence for subtype assignment; 2j/2f/2g etc. come from literature-labelled sequences)
ref_keys = [k for k in keys if k.startswith("REF_")]
def mean_j(members, ref):
    vals = []
    for m in members:
        sm, sr = sk[m], sk[ref]
        u = len(sm | sr)
        vals.append(len(sm & sr) / u if u else 0.0)
    return sum(vals) / len(vals)

n_div = 0
cluster_notes = {}
for root, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    n_div += 1
    lab = f"div-{chr(64 + n_div)}" if n_div <= 26 else f"div-{n_div}"
    if len(mem) < 2:
        lab = f"div-single({mem[0]})"
    scores = sorted(((mean_j(mem, rf), rf) for rf in ref_keys), reverse=True)
    if scores and scores[0][0] >= 0.15:
        lab += f"≈{scores[0][1].replace('REF_', '')}(J={scores[0][0]:.2f})"
    cluster_notes[lab] = scores[:3]
    for m in mem:
        assign[m] = lab
# annotate type1/type2 clusters as well
for kk, lab in list(assign.items()):
    pass  # type1/type2 labels are kept unchanged (the reference definition is unambiguous)

cluster_of = assign
print(f"assigned: type1/2={sum(1 for v in assign.values() if v.startswith('type1') or v.startswith('type2'))}, "
      f"mixed={sum(1 for v in assign.values() if v.startswith('type-mixed'))}, "
      f"divergent={len(div_pool)} ({len(groups)} clusters)", flush=True)

# --- output ---
mat = {}
with open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mat[r["ena_acc"]] = r

out = os.path.join(RESDIR, "MY9_P1聚类.tsv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ena_acc", "P1_cluster", "P1_type", "anchor_hits", "window_len",
                "country", "year", "ST", "resist_23S"])
    for acc in sorted(windows):
        m = mat.get(acc, {})
        w.writerow([acc, cluster_of.get(acc, "?"), m.get("P1_type", ""), windows[acc][1],
                    len(members[acc][0]), m.get("country", ""), m.get("year", ""),
                    m.get("MLST_ST", ""), m.get("resist_23S", "")])
print(f"written {out}")

# --- cluster-level summary ---
cl = collections.defaultdict(lambda: {"n": 0, "ctry": collections.Counter(), "yr": [],
                                      "st": collections.Counter(), "res": collections.Counter()})
for acc, lab in cluster_of.items():
    m = mat.get(acc, {})
    d = cl[lab]
    d["n"] += 1
    d["ctry"][m.get("country", "").split(":")[0].strip() or "NA"] += 1
    if m.get("year", "").isdigit():
        d["yr"].append(int(m["year"]))
    d["st"][m.get("MLST_ST", "?")] += 1
    d["res"][m.get("resist_23S", "-")] += 1

md = ["# MY9 P1 分歧谱系聚类摘要（W4）", "",
      f"- 提取成功 {len(windows)}/{len(files)}（anchor 采样 {len(samps)} 个, ≥{MIN_HITS} 命中）",
      f"- 单连锁聚类阈值 J≥{JAC_T}；簇数（不含参考）= {len(cl)}", "",
      "| 簇 | n | 主要国家 | 年份范围 | 主 ST | 23S 耐药率 |", "|---|---|---|---|---|---|"]
for lab, d in sorted(cl.items(), key=lambda kv: -kv[1]["n"]):
    yrs = f"{min(d['yr'])}–{max(d['yr'])}" if d["yr"] else "-"
    nres = d["res"].get("resistant", 0)
    ncalled = d["res"].get("resistant", 0) + d["res"].get("WT", 0)
    freq = f"{nres/ncalled:.0%}" if ncalled else "-"
    md.append(f"| {lab} | {d['n']} | {', '.join(f'{k}({v})' for k, v in d['ctry'].most_common(3))} | {yrs} | "
              f"{', '.join(f'{k}({v})' for k, v in d['st'].most_common(2))} | {freq} |")
md += ["", "> 提取= M129 P1 采样 20-mer 定位（±500 nt 窗）；聚类= 窗口 20-mer Jaccard 单连锁；",
       "> 参考成员: M129-P1（type1 锚）、FH-P1（type2 锚）、PHOL panel×5（reference/p1_panel_ncbi.fasta）。"]
with open(os.path.join(RESDIR, "MY9_P1聚类摘要.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(md) + "\n")
print("\n".join(md[7:]))
