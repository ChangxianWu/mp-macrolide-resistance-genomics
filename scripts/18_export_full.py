# -*- coding: utf-8 -*-
"""18_export_full.py — TreeTime formal dating inputs: full core alignment + the ST3 subtree trio
Outputs: results/coresnp/treetime/full_core.aln.fasta (528 strains + M129, core 742,308 sites)
         results/coresnp/treetime/MY9_coreSNP_NJ.nwk (rebuilt, same algorithm/input as MY9_15)
         results/coresnp/treetime/full_dates.csv
         results/coresnp/treetime/ST3_subtree.nwk / ST3.aln.fasta / ST3_dates.csv
Idempotent: skips when full_core.aln.fasta and the NJ reconstruction already exist (a rerun only fills in the ST3 part).
All tree operations use explicit-stack iteration — the binary NJ tree has ~1055 nodes and recursion would overflow the stack.
"""
import csv, glob, os, re, time
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(BASE, "results")
CS = os.path.join(RESDIR, "coresnp")
TD = os.path.join(CS, "treetime")
os.makedirs(TD, exist_ok=True)
BIOT = "ACGT"

accs_all = [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(CS, "alleles", "*.npy")))]
M = np.stack([np.load(os.path.join(CS, "alleles", a + ".npy")) for a in accs_all])
n, L = M.shape
cov = (M < 4).mean(axis=1)
keep = cov >= 0.80
accs = [a for a, k in zip(accs_all, keep) if k]
C = M[keep]
print(f"genomes kept: {len(accs)}/{n}", flush=True)

core_mask = (M < 4).sum(axis=0) >= int(0.95 * n)
core_len = int(core_mask.sum())
print(f"core sites: {core_len}", flush=True)

def embl_seq(path):
    embl = open(path, encoding="utf-8").read()
    body = embl.split("SQ   ", 1)[1].split("//", 1)[0]
    return "".join(re.findall(r"[A-Za-z]+", " ".join(body.splitlines()[1:]))).upper()

ref = embl_seq(os.path.join(BASE, "reference", "U00089.2.embl.cache.txt"))
ci = np.nonzero(core_mask)[0]

mat = {r["ena_acc"]: r for r in csv.DictReader(open(os.path.join(RESDIR, "MY9_分型矩阵.tsv"), encoding="utf-8-sig"), delimiter="\t")}

if not os.path.exists(os.path.join(TD, "full_core.aln.fasta")):
    print("writing full core alignment (~400 MB)...", flush=True)
    with open(os.path.join(TD, "full_core.aln.fasta"), "w") as f:
        def emit(nm, row, ref_row=False):
            f.write(f">{nm}\n")
            for s0 in range(0, len(ci), 100_000):
                if ref_row:
                    s = "".join(ref[p] for p in ci[s0:s0 + 100_000])
                else:
                    s = "".join(BIOT[b] if b < 4 else "N" for b in row[ci[s0:s0 + 100_000]])
                for i in range(0, len(s), 80):
                    f.write(s[i:i + 80] + "\n")
        emit("M129_reference", None, ref_row=True)
        for i, a in enumerate(accs):
            emit(a, C[i])
    print("full aln written", flush=True)
else:
    print("full_core.aln.fasta already present -> skip", flush=True)

# ---- NJ tree reconstruction (same algorithm/input as MY9_15 -> identical topology) ----
lines = open(os.path.join(CS, "MY9_SNP距离.phylip"), encoding="utf-8").read().strip().splitlines()
nk = int(lines[0])
accs_d = [ln.split("\t")[0] for ln in lines[1:]]
Dfull = [[float(x) for x in ln.split("\t")[1].split()] for ln in lines[1:]]
assert set(accs_d) == set(accs), "tree taxa != kept genomes"
order = [accs_d.index(a) for a in accs]
D = np.array([[Dfull[i][j] for j in order] for i in order])

active = list(range(len(accs)))
node_name = {i: accs[i] for i in range(len(accs))}
D2 = {i: dict(zip(range(len(accs)), D[i].tolist())) for i in range(len(accs))}
br = {}
nxt = len(accs)
t0 = time.time()
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
    nid = nxt; nxt += 1
    D2[nid] = {}
    br[nid] = [(bi, max(lb, 1e-12)), (bj, max(dab - lb, 1e-12))]
    node_name[nid] = None
    for a in active:
        if a in (bi, bj):
            continue
        dnew = 0.5 * (D2[bi][a] + D2[bj][a] - dab)
        D2[nid][a] = dnew
        D2[a][nid] = dnew
    D2[nid][nid] = 0.0
    active = [a for a in active if a not in (bi, bj)] + [nid]
a, b = active
dab = D2[a][b]
root = nxt
br[root] = [(a, dab / 2), (b, dab / 2)]
print(f"NJ rebuilt in {time.time()-t0:.0f}s", flush=True)

def build_newick(node):
    # Iterative post-order: leaf -> name; internal -> (child1,child2)
    out = {}
    stack = [(node, False)]
    while stack:
        x, exp = stack.pop()
        if exp:
            out[x] = "(" + ",".join(f"{out[c]}:{l:.8f}" for c, l in br[x]) + ")"
        elif x in out:
            continue
        else:
            if node_name.get(x):
                out[x] = node_name[x]
            else:
                stack.append((x, True))
                for c, _ in br[x]:
                    stack.append((c, False))
    return out[node]

tree_newick = build_newick(root) + ";"
open(os.path.join(TD, "MY9_coreSNP_NJ.nwk"), "w").write(tree_newick + "\n")
print("NJ newick written", flush=True)

with open(os.path.join(TD, "full_dates.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "date"])
    for a in accs:
        y = mat.get(a, {}).get("year", "")
        if y.isdigit():
            w.writerow([a, y])

# ---- ST3 MRCA: one bottom-up pass to count ST3 leaves per node subtree, then walk top-down ----
st3_set = set(i for i, a in enumerate(accs) if mat.get(a, {}).get("MLST_ST") == "3")
total_st3 = len(st3_set)
print(f"ST3 leaves: {total_st3}", flush=True)

st3cnt = {}
stack = [(root, False)]
while stack:
    x, exp = stack.pop()
    if x in st3cnt:
        continue
    if node_name.get(x):
        st3cnt[x] = 1 if x in st3_set else 0
        continue
    if exp:
        st3cnt[x] = sum(st3cnt[c] for c, _ in br[x])
    else:
        stack.append((x, True))
        for c, _ in br[x]:
            stack.append((c, False))
assert st3cnt[root] == total_st3, "root subtree ST3 count mismatch"

mrca = root
while True:
    kids = [(c, l) for c, l in br[mrca] if st3cnt[c] == total_st3]
    if len(kids) != 1:
        break
    mrca = kids[0][0]
print("ST3 MRCA node:", mrca, flush=True)

# ---- ST3 subtree pruning (iterative, including single-branch collapse) ----
def prune(mr, seen=None):
    """Return (pruned newick fragment, effective terminal branch length, leaf count)"""
    seen = set()
    out = {}
    stack = [(mr, 0.0, False)]
    while stack:
        x, e, exp = stack.pop()
        if exp:
            if node_name.get(x):
                if x in st3_set:
                    out[(x, e)] = (node_name[x], e, 1)
                else:
                    out[(x, e)] = ("", 0.0, 0)
            else:
                surv = [(c, l) for c, l in br.get(x, []) if out.get((c, l), ("", 0.0, 0))[2] > 0]
                k = sum(out[(c, l)][2] for c, l in surv)
                if k == 0:
                    out[(x, e)] = ("", 0.0, 0)
                elif k == 1:
                    cc, ll = surv[0]
                    s, ee, _ = out[(cc, ll)]
                    out[(x, e)] = (s, e + ee, k)
                else:
                    inner = "(" + ",".join(f"{out[(c, l)][0]}:{out[(c, l)][1]:.8f}" for c, l in surv) + ")"
                    out[(x, e)] = (inner, e, k)
        elif (x, e) in out:
            continue
        else:
            stack.append((x, e, True))
            for c, l in br.get(x, []):
                stack.append((c, l, False))
    return out[(mr, 0.0)]

s_sub, e_sub, k_sub = prune(mrca)
if k_sub == 0:
    raise SystemExit("no ST3 leaves under MRCA — abort")
if k_sub == 1:
    st3_newick = f"{s_sub}:{e_sub:.8f};"
else:
    st3_newick = s_sub + ";"
names_in = re.findall(r"[)(,]([A-Za-z0-9_.]+):", st3_newick)
assert len(names_in) == k_sub, f"newick leaf count {len(names_in)} != {k_sub}"
open(os.path.join(TD, "ST3_subtree.nwk"), "w").write(st3_newick + "\n")

with open(os.path.join(TD, "ST3.aln.fasta"), "w") as f:
    for i, a in enumerate(accs):
        if i in st3_set:
            f.write(f">{a}\n")
            for s0 in range(0, len(ci), 100_000):
                s = "".join(BIOT[b] if b < 4 else "N" for b in C[i][ci[s0:s0 + 100_000]])
                for i2 in range(0, len(s), 80):
                    f.write(s[i2:i2 + 80] + "\n")

with open(os.path.join(TD, "ST3_dates.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "date"])
    for i, a in enumerate(accs):
        if i in st3_set:
            y = mat.get(a, {}).get("year", "")
            if y.isdigit():
                w.writerow([a, y])

print(f"ST3 subtree exports written: leaves={k_sub}", flush=True)
print("all treetime inputs ->", TD, flush=True)
