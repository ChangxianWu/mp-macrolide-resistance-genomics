# -*- coding: utf-8 -*-
"""08_p1_type.py — P1 adhesin genotyping (type 1 / type 2 + specific 20-mer classifier, v2)
Input: ../data/fasta/*.fna + reference/U00089.2.embl.cache.txt + reference/CP010546.embl.cache.txt
Output: the per-genome P1 typing table under ../results/ + console summary

Key v2 changes (relative to the v1 probe panel):
  1. Reference genomes are always parsed from the **EMBL SQ block** (single source for annotation + sequence). Reason: in this environment the ENA FASTA endpoint
     once returned content inconsistent with the annotation coordinate system for U00089.2 (816,466 nt vs the official 816,394 bp, P1 start not ATG),
     which had misaligned the classifier input; the current FASTA set has been spot-checked and is stable with no such issue.
  2. The classifier now counts over the full set of P1-specific 20-mers (the v1 divergent-window 30-mer probe was too strict for field strains: P1 evolves fast,
     and the most divergent region matched the two references too rarely → 380 ambiguous/145 not_called).
     Specific set = kmer(P1_ref) - kmer(other P1) - kmer(rest of the reference genome) (to exclude RepMP copy cross-matches).
  3. Validation: M129 genome = 816,394 bp; P1 slice starts with ATG and has length 4884/4905.
Method (literature basis):
  M129 (ATCC 29342) = P1 type 1 reference; FH (ATCC 15531, CP010546) = type 2 reference (Schwartz 2009 JCM).
  M129 MPN141 = 180858..185741 (/gene="P1"); FH P1 = F539_00815 = 179324..184228 (annotation product="adhesin",
  the only one of 4.5–5.3 kb). Subtypes 1a/1b/2a/2b/2c(+2d) require additional references and will be refined later; major types 1/2 are the W3 deliverable.
"""
import csv, glob, os, re, sys, collections

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFDIR = os.path.join(BASE, "reference")
RESDIR = os.path.join(BASE, "results")
K = 20
RATIO = 0.8
MIN_HITS = 30

def fseq_text(t):
    """Concatenate all FASTA records (the P1 of a multi-contig assembly may not be on the first contig — the initial v2 version took only the first,
    which produced a spurious 153 not_called; the k-mer classifier is insensitive to concatenation order, so concatenation is safe)"""
    return "".join(l.strip() for l in t.splitlines() if not l.startswith(">")).upper()

def embl_sq_seq(embl):
    """Parse the EMBL SQ block into a genome sequence (supports IUPAC ambiguity codes such as m/r/s)"""
    body = embl.split("SQ   ", 1)[1].split("//", 1)[0]
    lines = body.splitlines()
    seq = "".join(re.findall(r"[A-Za-z]+", " ".join(lines[1:])))
    return seq.upper()

def embl_cds(embl, want):
    for m in re.finditer(r"FT   CDS\s+([^\n]+)\n((?:FT   {17}\S[^\n]*\n)+)", embl):
        quals = m.group(2)
        if want(quals):
            rng = m.group(1).strip()
            mm = re.search(r"complement\((\d+)\.\.(\d+)\)", rng) or re.search(r"(\d+)\.\.(\d+)", rng)
            return int(mm.group(1)), int(mm.group(2)), "complement" in rng
    return None

RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s):
    return s.translate(RC)[::-1]

# --- 1) M129 (type 1 reference): EMBL SQ ---
embl129 = open(os.path.join(REFDIR, "U00089.2.embl.cache.txt"), encoding="utf-8").read()
g129 = embl_sq_seq(embl129)
assert len(g129) == 816394, f"M129 EMBL genome len {len(g129)} != 816394"
p1_cds = embl_cds(embl129, lambda q: '/gene="P1"' in q and "MPN_141" in q)
assert p1_cds, "M129 P1 CDS 未找到"
s0_m129, e_m129, comp = p1_cds
p1_type1 = g129[s0_m129-1:e_m129]
if comp:
    p1_type1 = rc(p1_type1)
assert p1_type1.startswith("ATG") and len(p1_type1) == 4884, \
    f"M129 P1 切片异常: start={p1_type1[:3]} len={len(p1_type1)}"
print(f"M129 genome={len(g129)}  P1(type1 ref) {s0_m129}..{e_m129} len={len(p1_type1)} start=ATG ✓", flush=True)

# --- 2) FH (type 2 reference): EMBL SQ ---
FA_FH_EMBL = os.path.join(REFDIR, "CP010546.embl.cache.txt")
assert os.path.exists(FA_FH_EMBL), "CP010546 EMBL 缓存缺失（先运行 08 v1 或手动抓取）"
embl_fh = open(FA_FH_EMBL, encoding="utf-8").read()
fh_g = embl_sq_seq(embl_fh)
fh_len_expect = int(re.search(r"^SQ   Sequence (\d+) BP", embl_fh, re.M).group(1))
assert len(fh_g) == fh_len_expect, f"FH genome len {len(fh_g)} != SQ 声明 {fh_len_expect}"
fh_p1 = None
for m in re.finditer(r"FT   CDS\s+([^\n]+)\n((?:FT   {17}\S[^\n]*\n)+)", embl_fh):
    quals = m.group(2)
    if '/product="adhesin"' in quals:
        rng = m.group(1).strip()
        mm = re.search(r"complement\((\d+)\.\.(\d+)\)", rng) or re.search(r"(\d+)\.\.(\d+)", rng)
        s0, s1 = int(mm.group(1)), int(mm.group(2))
        if 4500 <= s1 - s0 + 1 <= 5300:
            assert fh_p1 is None, "FH 有多个候选 P1 CDS"
            fh_p1 = (s0, s1, "complement" in rng)
assert fh_p1, "FH P1 CDS 未找到"
s0_fh, e_fh, comp = fh_p1
p1_type2 = fh_g[s0_fh-1:e_fh]
if comp:
    p1_type2 = rc(p1_type2)
assert p1_type2.startswith("ATG") and len(p1_type2) == 4905, \
    f"FH P1 切片异常: start={p1_type2[:3]} len={len(p1_type2)}"
print(f"FH genome={len(fh_g)}  P1(type2 ref) {s0_fh}..{e_fh} len={len(p1_type2)} start=ATG ✓", flush=True)

# --- 3) specific 20-mer sets ---
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def kmers(seq, k=K):
    """canonical k-mer set (min(fwd, rc)) — also valid for the P1 region of reverse-strand contigs"""
    out = set()
    for i in range(len(seq) - k + 1):
        x = seq[i:i+k]
        if "N" not in x:
            out.add(min(x, x.translate(_RC)[::-1]))
    return out

k11, k12 = kmers(p1_type1), kmers(p1_type2)
rest129 = kmers(g129[:s0_m129 - 1] + g129[e_m129:])
restfh = kmers(fh_g[:s0_fh - 1] + fh_g[e_fh:])
spec1 = k11 - k12 - rest129 - restfh
spec2 = k12 - k11 - rest129 - restfh
print(f"specific {K}-mers: type1={len(spec1)} type2={len(spec2)}", flush=True)

# --- 4) classification ---
def classify(genome):
    kg = kmers(genome)
    h1, h2 = len(kg & spec1), len(kg & spec2)
    tot = h1 + h2
    if tot == 0:
        return "not_called", h1, h2
    if max(h1, h2) < MIN_HITS:
        return "ambiguous", h1, h2
    r = h1 / tot
    if r >= RATIO:
        return "type1", h1, h2
    if r <= 1 - RATIO:
        return "type2", h1, h2
    return "ambiguous", h1, h2

# --- metadata join (same as 04) ---
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
    sys.exit("no fasta under data/fasta")
os.makedirs(RESDIR, exist_ok=True)
out = os.path.join(RESDIR, "MY9_P1分型.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ena_acc", "country", "year", "P1_type", "hits_type1", "hits_type2",
                "n_spec_type1", "n_spec_type2"])
    for k, path in enumerate(files, 1):
        acc = os.path.basename(path).replace(".fna", "")
        g = fseq_text(open(path).read())
        t, v1, v2 = classify(g)
        m = meta.get(acc, {})
        w.writerow([acc, m.get("country", ""), m.get("year", ""), t, v1, v2, len(spec1), len(spec2)])
        if k % 50 == 0:
            print(f"typed {k}/{len(files)}", flush=True)

# --- summary ---
rows = list(csv.DictReader(open(out), delimiter="\t"))
c = collections.Counter(r["P1_type"] for r in rows)
print("--- P1 type(total) ---")
for kk, v in c.most_common():
    print(f"{kk}: {v}")
print("--- by country (n>=5) ---")
cc = collections.Counter((r["country"].split(":")[0].strip(), r["P1_type"]) for r in rows
                         if r["country"] and r["country"] != "not applicable")
nations = collections.Counter(r["country"].split(":")[0].strip() for r in rows)
for nat, n in nations.most_common(12):
    if n < 5:
        continue
    t1 = cc.get((nat, "type1"), 0); t2 = cc.get((nat, "type2"), 0)
    print(f"{nat}: n={n} type1={t1} type2={t2} not_typed={n-t1-t2}")
print("--- by year (n>=5) ---")
yy = collections.Counter((r["year"], r["P1_type"]) for r in rows if r["year"].isdigit())
years = sorted({r["year"] for r in rows if r["year"].isdigit()}, key=int)
for y in years:
    n = sum(v for (a, _), v in yy.items() if a == y)
    if n < 5:
        continue
    print(f"{y}: n={n} type1={yy.get((y,'type1'),0)} type2={yy.get((y,'type2'),0)}")
print(f"written {out}")
