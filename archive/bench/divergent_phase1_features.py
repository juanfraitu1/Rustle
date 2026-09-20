#!/usr/bin/env python3
"""Phase 1: does an RNA-ONLY feature separate true divergent paralogs (Compara-confirmed) from prefilter
false-groupings, on the gold calibration set? Computes k-mer features from the de-novo transcript NUCLEOTIDE
sequences (denovo_transcripts.fa) for each labeled locus pair and reports AUC vs the Compara label.

Features (all RNA-only, no genome/protein at inference):
  core_recip   : contiguous-core coverage (the shipped baseline that FAILS to separate) -- from edges file.
  raw_shared   : # shared canonical 18-mers.
  jaccard      : |A∩B| / |A∪B| over canonical 18-mer SETS.
  containment  : |A∩B| / min(|A|,|B|).
  span_cov     : span of shared-kmer positions in the shorter tx / its length -- shared k-mers SPREAD across
                 the whole transcript (true paralog) vs concentrated in one domain (false grouping). THE hypothesis.
  len_ratio    : min(len)/max(len).
KILL if best AUC < 0.70.
"""
import sys, os, csv, bisect, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compara_fetch as cf
import pysam

ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FASTA = "/home/juanfra/winloci_scratch/denovo_transcripts.fa"
EDGES = "denovo_family_edges.tsv"
CORE_MAX = 0.20
K = 18

# --- gene -> (ensg, paralogue set) STRICTLY from the cache (no network: error/missing entries -> unmapped) ---
def gene_info_from_cache(symbols):
    cache = cf.load_cache()
    info = {}
    for sym in symbols:
        lk = cache.get(cf.cache_key("lookup", sym))
        ensg = lk["data"].get("id") if lk and lk.get("status") == "ok" and isinstance(lk.get("data"), dict) else None
        para = set()
        pa = cache.get(cf.cache_key("homology_paralogues", sym))
        if pa and pa.get("status") == "ok":
            data = pa.get("data") or {}
            d = data.get("data") if isinstance(data, dict) else None
            homs = d[0].get("homologies", []) if isinstance(d, list) and d else []
            para = {h.get("id") for h in homs if h.get("id")}
        info[sym] = (ensg, para)
    return info

# --- locus -> gene (max overlap) ---
genes = collections.defaultdict(list)
for r in csv.DictReader(open(ANNOT), delimiter="\t"):
    genes[r["chrom"]].append((int(r["start"]), int(r["end"]), r["gene"]))
for c in genes:
    genes[c].sort()
starts = {c: [g[0] for g in v] for c, v in genes.items()}
def map_locus(ch, s, e):
    if ch not in genes:
        return None
    arr, st = genes[ch], starts[ch]
    i = bisect.bisect_right(st, e) - 1
    best, bov = None, 0
    while i >= 0 and arr[i][1] >= s - 200_000:
        gs, ge, nm = arr[i]
        ov = min(e, ge) - max(s, gs)
        if ov > bov:
            bov, best = ov, nm
        if gs < s - 200_000:
            break
        i -= 1
    return best if bov > 0 else None

coord, nexon = {}, {}
for r in csv.DictReader(open(META), delimiter="\t"):
    coord[r["id"]] = (r["chrom"], int(r["start"]), int(r["end"]))
    nexon[r["id"]] = int(r["n_exon"])

# --- labeled LOCUS pairs (core<0.20, both-named diff-gene, Compara-labelable) ---
raw_pairs = []  # (locA, locB, core, geneA, geneB)
gene_set = set()
for r in csv.DictReader(open(EDGES), delimiter="\t"):
    cr = float(r["core_recip"])
    if cr >= CORE_MAX or r["a"] not in coord or r["b"] not in coord:
        continue
    ga, gb = map_locus(*coord[r["a"]]), map_locus(*coord[r["b"]])
    if not ga or not gb or ga.startswith("LOC") or gb.startswith("LOC") or ga == gb:
        continue
    raw_pairs.append((r["a"], r["b"], cr, ga, gb))
    gene_set.add(ga); gene_set.add(gb)

info = gene_info_from_cache(sorted(gene_set))
pairs = []  # (locA, locB, core, label)
for la, lb, cr, ga, gb in raw_pairs:
    ea, pa = info.get(ga, (None, set()))
    eb, pb = info.get(gb, (None, set()))
    if not ea or not eb:
        continue
    label = 1 if (eb in pa or ea in pb) else 0
    pairs.append((la, lb, cr, label))
npos = sum(p[3] for p in pairs); nneg = len(pairs) - npos
print("labeled locus pairs: %d  (pos=%d neg=%d)" % (len(pairs), npos, nneg))

# --- canonical 18-mers (code via 2-bit, min(fwd,rc)) per locus, with first positions ---
fa = pysam.FastaFile(FASTA)
B = {65: 0, 67: 1, 71: 2, 84: 3}  # A C G T
def low_complexity(window):
    """DUST-like: an 18-mer is low-complexity if it has < 8 distinct 3-mers (repeats/homopolymers)."""
    return len({window[i:i+3] for i in range(len(window) - 2)}) < 8
def canon(seq):
    """Return (full_kmers, masked_kmers, n): both code->first-pos dicts; masked drops low-complexity windows."""
    s = seq.upper().encode()
    n = len(s); km = {}; kmm = {}
    if n < K:
        return km, kmm, n
    fwd = rev = 0
    mask = (1 << (2 * K)) - 1
    valid = 0
    for i, ch in enumerate(s):
        b = B.get(ch, -1)
        if b < 0:
            valid = 0; fwd = rev = 0; continue
        fwd = ((fwd << 2) | b) & mask
        rev = (rev >> 2) | ((3 - b) << (2 * (K - 1)))
        valid += 1
        if valid >= K:
            code = fwd if fwd < rev else rev
            pos = i - K + 1
            if code not in km:
                km[code] = pos
            if code not in kmm and not low_complexity(s[pos:pos + K]):
                kmm[code] = pos
    return km, kmm, n

cache_km = {}
def get_km(loc):
    if loc not in cache_km:
        try:
            cache_km[loc] = canon(fa.fetch(loc))
        except Exception:
            cache_km[loc] = ({}, {}, 0)
    return cache_km[loc]

# --- protein-ORF: translate the transcript's own longest ORF (RNA-derived, self-contained) and take its
#     5-mer set. Divergent paralogs conserve PROTEIN even when synonymous DNA divergence kills DNA k-mers.
CODON = {}
_bases = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_aa):
    CODON["".join(_bases[(_i >> s) & 3] for s in (4, 2, 0))] = _a
PK = 5
def protein_kmers(seq):
    s = seq.upper()
    best = ""
    for f in range(3):
        aa = []
        for i in range(f, len(s) - 2, 3):
            aa.append(CODON.get(s[i:i+3], "X"))
        pep = "".join(aa)
        for run in pep.split("*"):          # longest stop-free peptide across the 3 forward frames
            if len(run) > len(best):
                best = run
    return {best[i:i+PK] for i in range(len(best) - PK + 1)} if len(best) >= PK else set()
cache_pk = {}
def get_pk(loc):
    if loc not in cache_pk:
        try:
            cache_pk[loc] = protein_kmers(fa.fetch(loc))
        except Exception:
            cache_pk[loc] = set()
    return cache_pk[loc]

# --- features per pair ---
feats = collections.defaultdict(lambda: ([], []))  # name -> (pos_vals, neg_vals)
def add(name, val, label):
    feats[name][0 if label else 1].append(val) if False else None  # placeholder
def cont_of(ka, kb):
    sa, sb = set(ka), set(kb)
    inter = sa & sb
    raw = len(inter)
    cont = raw / min(len(sa), len(sb)) if sa and sb else 0.0
    return raw, cont, inter
for la, lb, cr, label in pairs:
    ka, kam, na = get_km(la); kb, kbm, nb = get_km(lb)
    if not ka or not kb:
        continue
    sa, sb = set(ka), set(kb)
    inter = sa & sb
    raw = len(inter)
    uni = len(sa | sb)
    jac = raw / uni if uni else 0.0
    cont = raw / min(len(sa), len(sb)) if sa and sb else 0.0
    short_k, short_len = (ka, na) if na <= nb else (kb, nb)
    ipos = [short_k[c] for c in inter if c in short_k]
    span_cov = (max(ipos) - min(ipos)) / short_len if len(ipos) >= 2 and short_len else 0.0
    len_ratio = min(na, nb) / max(na, nb) if max(na, nb) else 0.0
    # MASKED (low-complexity removed) — attacks the repeat-driven false sharing in the negatives.
    m_raw, m_cont, _ = cont_of(kam, kbm)
    # EXON STRUCTURE — paralogs conserve architecture; domain-sharers / read-throughs do not.
    ea, eb = nexon.get(la, 0), nexon.get(lb, 0)
    exon_ratio = min(ea, eb) / max(ea, eb) if max(ea, eb) else 0.0
    exon_diff = -abs(ea - eb)  # negative so "more similar" is larger (AUC is direction-agnostic anyway)
    # PROTEIN-ORF — the lever for divergent paralogs (protein conserved under synonymous DNA divergence).
    pa, pb = get_pk(la), get_pk(lb)
    pint = len(pa & pb)
    prot_cont = pint / min(len(pa), len(pb)) if pa and pb else 0.0
    prot_jac = pint / len(pa | pb) if (pa or pb) else 0.0
    for name, val in (("core_recip", cr), ("raw_shared", raw), ("jaccard", jac),
                      ("containment", cont), ("span_cov", span_cov), ("len_ratio", len_ratio),
                      ("masked_raw", m_raw), ("masked_cont", m_cont),
                      ("exon_ratio", exon_ratio), ("exon_diff", exon_diff),
                      ("prot_cont", prot_cont), ("prot_jac", prot_jac)):
        feats[name][0 if label else 1].append(val)

def auc(pos, neg):
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return float("nan")
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks = [0.0] * len(allv); i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = r
        i = j
    sp = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 1)
    u = sp - n1 * (n1 + 1) / 2.0
    a = u / (n1 * n0)
    return max(a, 1 - a)  # direction-agnostic discriminative power

import statistics as st
print("\n=== Phase 1: feature AUC vs Compara label (direction-agnostic) ===")
print("%-13s %6s   %12s %12s" % ("feature", "AUC", "mean(pos)", "mean(neg)"))
best = ("", 0.0)
for name in ("core_recip", "raw_shared", "jaccard", "containment", "span_cov", "len_ratio",
             "masked_raw", "masked_cont", "exon_ratio", "exon_diff", "prot_cont", "prot_jac"):
    pos, neg = feats[name]
    a = auc(pos, neg)
    print("%-13s %6.3f   %12.4f %12.4f" % (name, a, st.mean(pos) if pos else 0, st.mean(neg) if neg else 0))
    if a > best[1]:
        best = (name, a)
print("\nBEST: %s AUC=%.3f  -> %s" % (best[0], best[1], "PASS (>=0.70)" if best[1] >= 0.70 else "KILL (<0.70)"))
