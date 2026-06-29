#!/usr/bin/env python3
"""Vet the unannotated loci with NO gorilla-protein hit — separate genuinely NOVEL gene families
(multi-copy paralogs + coherent ORF + spliced + complex sequence) from lncRNA (spliced, no ORF) and
NOISE (low-complexity / single-exon-genomic / repeat). Uses the already-computed transcripts, the
genome-mapping PAF (paralog count), the blastx hits, and the de-novo skeleton (exon structure)."""
import sys, json, math
from itertools import product
sys.path.insert(0, "bench")
import copy_assign as CA

WT = "/home/juanfra/winloci_scratch/refabsent"
COMP = str.maketrans("ACGTN", "TGCAN")
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): aas[i] for i, c in enumerate(product("TCAG", repeat=3))}
def orf_aa(s):
    best = 0
    for strand in (s, s.translate(COMP)[::-1]):
        for f in range(3):
            prot = "".join(CODON.get(strand[i:i+3], "X") for i in range(f, len(strand)-2, 3))
            for p in prot.split("*"):
                k = p.find("M")
                if k >= 0: best = max(best, len(p)-k)
    return best
def cplx4(s):  # distinct 4-mers / possible = linguistic complexity (1=complex, low=repetitive)
    if len(s) < 8: return 0
    k = {s[i:i+4] for i in range(len(s)-3)}
    return len(k) / min(256, len(s)-3)

# load transcripts + protein hits + paralog counts
tx = {}
hdr = None
for line in open(f"{WT}/unannot_tx.fa"):
    if line.startswith(">"): hdr = line[1:].strip()
    else: tx[hdr] = line.strip()
hit = set(l.split("\t")[0] for l in open(f"{WT}/unannot_blastx.tsv"))
from collections import defaultdict
hits = defaultdict(list)
for line in open(f"{WT}/unannot_tx.paf"):
    f = line.split("\t")
    if len(f) < 12: continue
    if int(f[9])/max(1,int(f[10])) >= 0.80 and int(f[10])/max(1,int(f[1])) >= 0.3:
        hits[f[0]].append((f[5], int(f[7])))
def n_loci(q):
    pts = sorted(hits.get(q, []))
    L = []
    for c, p in pts:
        if L and L[-1][0] == c and p - L[-1][1] < 200000: L[-1] = (c, p)
        else: L.append((c, p))
    return len(L)
# exon counts from skeleton
meta = CA.load_meta(); skel = CA.load_skel()
def n_exons(tid):
    if tid not in meta: return 1
    c, s, e, st, nr = meta[tid]
    return len(skel.get((c, s, e), [])) + 1

noprot = [h for h in tx if h not in hit]
rows = []
for h in noprot:
    tid = h.split("|")[1]; seq = tx[h]
    nl = n_loci(h); orf = orf_aa(seq); cx = cplx4(seq); ne = n_exons(tid); nr = meta.get(tid, (0,0,0,0,0))[4]
    # classify
    if cx < 0.35:                          c = "NOISE (low-complexity/repeat)"
    elif nl >= 2 and orf >= 100 and ne >= 2: c = "NOVEL gene-FAMILY candidate"
    elif nl >= 2 and orf >= 100:           c = "novel multi-copy (single-exon)"
    elif orf >= 100 and ne >= 2:           c = "novel single-copy gene candidate"
    elif ne >= 2 and orf < 80:             c = "lncRNA candidate"
    else:                                   c = "other/ambiguous"
    rows.append(dict(locus=h.split("|")[0], tid=tid, n_paralog_loci=nl, orf_aa=orf,
                     n_exons=ne, tx_len=len(seq), cplx=round(cx, 2), reads=nr, klass=c))
from collections import Counter
print(f"=== vetting {len(noprot)} no-protein-hit unannotated loci ===")
for k, v in Counter(r["klass"] for r in rows).most_common(): print(f"  {v:>3}  {k}")
print("\n=== NOVEL gene-FAMILY candidates (multi-copy + ORF>=100aa + spliced + complex) ===")
nov = sorted([r for r in rows if r["klass"].startswith("NOVEL")], key=lambda r: (-r["n_paralog_loci"], -r["orf_aa"]))
print(f"{'locus':28s} {'paralogs':>8} {'ORFaa':>5} {'exons':>5} {'txlen':>5} {'cplx':>4} {'reads':>5}")
for r in nov:
    print(f"{r['locus']:28s} {r['n_paralog_loci']:>8} {r['orf_aa']:>5} {r['n_exons']:>5} {r['tx_len']:>5} {r['cplx']:>4} {r['reads']:>5}")
json.dump(rows, open(f"{WT}/unannotated_novel_vetted.json", "w"), indent=1)
print(f"\nwrote {WT}/unannotated_novel_vetted.json")
