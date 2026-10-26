# -*- coding: utf-8 -*-
"""12_mash_nj.py — Whole-genome distance phylogeny transition product (numpy + standard library)
Positioning: the formal core-SNP phylogeny (snippy→IQ-TREE→TreeTime) requires a server gate; this
      script first delivers preliminary temporal-signal evidence from Mash distance + NJ +
      root-to-tip (labelled in the paper's Methods as whole-genome k-mer distance phylogeny and
      cross-validated against the core-SNP tree).
Input: ../data/fasta/*.fna + the integrated typing matrix under ../results/ (year/ST/resistance annotations)
Output: the Mash distance matrix under results/ (full matrix, PHYLIP format, cached for reuse)
      results/MY9_NJ_tree.nwk (midpoint-rooted)
      results/MY9_root_to_tip.tsv (tip name, year, ST, resistance, root_dist)
      console: root-to-tip regression (all isolates / ST strata) R² and slope

Methods:
  sketch (v2): position-independent hash-threshold sampling — 2-bit packed 21-mers; after Knuth
    multiplicative hashing, those below the threshold enter the sketch
    (v1 fixed-stride sampling collapsed J under any 1 nt shift or contig-boundary difference:
     measured J=0.08 vs. the full-set truth of 0.9985, deprecated)
  distance: d = -ln(2J/(1+J))/21 (Mash formula)
  NJ: standard algorithm; midpoint rooting (root node inserted explicitly at the midpoint of the
    u-v path); Newick output
"""
import csv, glob, math, os, sys, time, collections

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
K = 21
SKETCH_FRAC = 0.05          # keep ~5% of k-mers → ~40k per isolate
GOLD = 0x9E3779B1
MASK32 = 0xFFFFFFFF

TBL = np.full(256, 255, dtype=np.uint8)   # default 255 = ambiguous/invalid → any k-mer containing it is dropped
for i, b in enumerate("ACGT"):
    TBL[ord(b)] = i
    TBL[ord(b.lower())] = i

def load_codes(path):
    seq = []
    for ln in open(path):
        if not ln.startswith(">"):
            seq.append(ln.strip())
    arr = np.frombuffer("".join(seq).encode(), dtype=np.uint8)
    return TBL[arr]

def sketch(codes):
    """Position-independent hash-threshold sketch (canonical 2-bit packed 21-mers, keep those with hash<threshold) → set(uint64)"""
    n = len(codes)
    W = n - K + 1
    if W <= 0:
        return set()
    c = codes.astype(np.uint64)
    packs = np.zeros(W, dtype=np.uint64)
    rpacks = np.zeros(W, dtype=np.uint64)
    for j in range(K):
        packs += c[j:j + W] << np.uint64(2 * (K - 1 - j))
        rpacks += (np.uint64(3) - c[j:j + W]) << np.uint64(2 * j)   # reverse-complement-equivalent packing
    canon = np.minimum(packs, rpacks)
    cs = np.concatenate(([0], np.cumsum((codes == 255).astype(np.int64))))
    badcnt = cs[K:] - cs[:W]          # count of ambiguous codes per window
    canon = canon[badcnt == 0]
    h = ((canon * np.uint64(GOLD)) & np.uint64(MASK32)).astype(np.float64) / 4294967296.0
    return set(canon[h < SKETCH_FRAC].tolist())

def fseq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()

def mash_d(ja, jb):
    inter = len(ja & jb)
    if inter == 0:
        return 0.25
    j = inter / (len(ja) + len(jb) - inter)
    return min(0.25, -math.log((2 * j) / (1 + j)) / K) if j > 0 else 0.25

files = sorted(glob.glob(os.path.join(BASE, "data", "fasta", "*.fna")))
names = [os.path.basename(p).replace(".fna", "") for p in files]
n = len(names)
print(f"genomes: {n}", flush=True)

PHY = os.path.join(RESDIR, "MY9_mash距离.phylip")
D = None
if os.path.exists(PHY):
    lines = open(PHY, encoding="utf-8").read().strip().splitlines()
    if lines and lines[0].strip() == str(n):
        D = [[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]]
        if len(D) == n:
            print(f"[cache] loaded distance matrix from {PHY}", flush=True)
        else:
            D = None
if D is None:
    sk = {}
    t0 = time.time()
    for i, p in enumerate(files, 1):
        sk[names[i - 1]] = sketch(load_codes(p))
        if i % 50 == 0:
            print(f"sketched {i}/{n} ({time.time()-t0:.0f}s)", flush=True)

    print("pairwise Mash distances...", flush=True)
    t0 = time.time()
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = mash_d(sk[names[i]], sk[names[j]])
            D[i][j] = D[j][i] = d
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{n} rows ({time.time()-t0:.0f}s)", flush=True)
    sk.clear()

    with open(PHY, "w", encoding="utf-8") as f:
        f.write(f"{n}\n")
        for i in range(n):
            f.write(names[i] + "\t" + " ".join(f"{D[i][j]:.6f}" for j in range(n)) + "\n")
    print("phylip written", flush=True)

# --- NJ ---
print("NJ tree building...", flush=True)
t0 = time.time()
active = list(range(n))
node_name = {i: names[i] for i in range(n)}
D2 = {i: {} for i in range(n)}
for i in range(n):
    for j in range(n):
        D2[i][j] = D[i][j]
br = {}  # internal node -> [(child, blen)]
next_id = n
while len(active) > 2:
    m = len(active)
    r = {a: sum(D2[a][b] for b in active) for a in active}
    best, bi, bj = 1e18, None, None
    for x, a in enumerate(active):
        Da = D2[a]
        for b in active[x + 1:]:
            v = Da[b] * (m - 2) - r[a] - r[b]
            if v < best:
                best, bi, bj = v, a, b
    dab = D2[bi][bj]
    lb = 0.5 * dab + (r[bi] - r[bj]) / (2 * (m - 2))
    lc = dab - lb
    nid = next_id; next_id += 1
    D2[nid] = {}
    br[nid] = [(bi, lb), (bj, lc)]
    node_name[nid] = None
    for a in active:
        if a in (bi, bj):
            continue
        dnew = 0.5 * (D2[bi][a] + D2[bj][a] - dab)
        D2[nid][a] = dnew
        D2[a][nid] = dnew
    D2[nid][nid] = 0.0
    active = [a for a in active if a not in (bi, bj)] + [nid]
# Final pair
a, b = active
dab = D2[a][b]
root_id = next_id; next_id += 1
br[root_id] = [(a, dab / 2), (b, dab / 2)]
print(f"NJ done in {time.time()-t0:.0f}s", flush=True)

# Recursively build Newick (with branch lengths)
def build(node):
    if node_name.get(node):
        return node_name[node]
    parts = []
    for c, L in br[node]:
        parts.append(f"{build(c)}:{L:.6f}")
    return "(" + ",".join(parts) + ")"

tree = build(root_id) + ";"
with open(os.path.join(RESDIR, "MY9_NJ_tree.nwk"), "w", encoding="utf-8") as f:
    f.write(tree + "\n")
print("newick written", flush=True)

# --- Midpoint rooting + root-to-tip ---
# Two-pass DFS to compute node depths (undirected tree)
adj = collections.defaultdict(list)
for node, children in br.items():
    for c, L in children:
        adj[node].append((c, L))
        adj[c].append((node, L))

def farthest(start):
    dist = {start: 0.0}
    stack = [start]
    while stack:
        u = stack.pop()
        for v, L in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + L
                stack.append(v)
    far = max(dist, key=dist.get)
    return far, dist[far]

u, _ = farthest(root_id)
v, diam = farthest(u)
# The point on the u-v path lying diam/2 from either end = midpoint; simplification: compute
# root-to-tip with u as the root, then subtract the offset
def depths_from(src):
    dist = {src: 0.0}
    stack = [src]
    while stack:
        x = stack.pop()
        for y, L in adj[x]:
            if y not in dist:
                dist[y] = dist[x] + L
                stack.append(y)
    return dist
du = depths_from(u)
dv = depths_from(v)
id_of = {names[i]: i for i in range(n)}
# root-to-tip (rooted at the u-v midpoint): d_mid(t) = (du[t] + dv[t] - diam)/2
mat = {}
with open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mat[r["ena_acc"]] = r

rows = []
for name in names:
    dmid = (du[id_of[name]] + dv[id_of[name]] - diam) / 2.0
    m = mat.get(name, {})
    rows.append({"tip": name, "year": m.get("year", ""), "ST": m.get("MLST_ST", ""),
                 "resist": m.get("resist_23S", ""), "root_dist": f"{dmid:.6f}"})
with open(os.path.join(RESDIR, "MY9_root_to_tip.tsv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tip", "year", "ST", "resist", "root_dist"])
    w.writeheader(); w.writerows(rows)

def regress(sel):
    pts = [(int(r["year"]), float(r["root_dist"])) for r in rows
           if r["year"].isdigit() and 2005 <= int(r["year"]) <= 2026 and sel(r)]
    if len(pts) < 10:
        return None
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    syy = sum((y - my) ** 2 for _, y in pts)
    slope = sxy / sxx
    r2 = (sxy * sxy) / (sxx * syy)
    return len(pts), slope, r2

for label, sel in [("all (2005-2026)", lambda r: True),
                   ("ST3 only", lambda r: r["ST"] == "3"),
                   ("ST14 only", lambda r: r["ST"] == "14"),
                   ("ST2 only (WT clone)", lambda r: r["ST"] == "2")]:
    out = regress(sel)
    if out:
        nn, slope, r2 = out
        print(f"root-to-tip {label}: n={nn} slope={slope:.2e} R2={r2:.3f}")
print("done -> results/MY9_NJ_tree.nwk, MY9_root_to_tip.tsv, MY9_mash距离.phylip")
