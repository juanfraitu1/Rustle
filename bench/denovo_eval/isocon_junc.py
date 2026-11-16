#!/usr/bin/env python3
"""PREREG r1764 (docs/PREREG_isocon_junctions_2026-11-16.md): junction precision of ours vs IsoCon on the target panel,
at a matched junction count, then slip classification of the wrong junctions.
  isocon_junc.py TARGETS.tsv REFSEQ.gtf BAM FASTA OURS.gtf ISOCON.gtf      (env as famtruth, gorilla)"""
import sys, os, re, collections, pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
tfile, ref, bam, fa, ours, iso = sys.argv[1:7]
G = F.load_ggo(F.load_cat()); idx = F.exon_index(G)
cls = {}
for l in list(open(tfile))[1:]:
    f = l.rstrip('\n').split('\t'); cls[f[2]] = f[1]
bounds = {g: ({e for _, e in G[g]['exons']}, {s for s, _ in G[g]['exons']}) for g in cls}
def exons(path):
    tx = collections.defaultdict(list); meta = {}
    for line in open(path):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
    return tx, meta
def scan(path):  # gene -> set of (donor0, acceptor0), as toolcompare.py
    tx, meta = exons(path); out = collections.defaultdict(set)
    for t, ex in tx.items():
        if len(ex) < 2: continue
        ex.sort(); g = F.map_copy(idx, meta[t][0], meta[t][1], ex)
        if g in cls: out[g] |= {(a[1], b[0]) for a, b in zip(ex, ex[1:])}
    return out
tx, meta = exons(ref); truth = collections.defaultdict(set)
for t, ex in tx.items():
    ex.sort(); truth[meta[t][0]] |= {(a[1], b[0]) for a, b in zip(ex, ex[1:])}
O, I = scan(ours), scan(iso)
blk = lambda R: sum(d in bounds[g][0] and a in bounds[g][1] for g in R for d, a in R[g]) / max(1, sum(map(len, R.values())))
print(f"step A, block truth, all modelled copies: ours {blk(O):.3f}  isocon {blk(I):.3f}")
B = pysam.AlignmentFile(bam); FA = pysam.FastaFile(fa)
RC = str.maketrans('ACGT', 'TGCA')
def support(g):  # (donor, acceptor) -> [primary, secondary] records with exactly that intron
    js = O[g] | I[g]; c = G[g]['chrom']; S = collections.defaultdict(lambda: [0, 0])
    for r in B.fetch(c, max(0, min(d for d, _ in js) - 1), max(a for _, a in js) + 1):
        if r.is_supplementary or r.is_unmapped: continue
        p = r.reference_start
        for op, n in r.cigartuples:
            if op == 3: S[(p, p + n)][r.is_secondary] += 1
            if op in (0, 2, 3, 7, 8): p += n
    return S
def motif(g, d, a):
    m = (FA.fetch(G[g]['chrom'], d, d + 2) + FA.fetch(G[g]['chrom'], a - 2, a)).upper()
    return m if G[g]['strand'] == '+' else m.translate(RC)[::-1]
def slip(g, j, J, S):  # annotated partner in the same tool's set, one end equal / both within 10 bp, <= 0.2 of its support
    c = G[g]['chrom']; n = sum(S[j]); best = None
    for k in J:
        if k == j or k not in truth[c]: continue
        dd, da = abs(k[0] - j[0]), abs(k[1] - j[1])
        if dd <= 10 and da <= 10 and (best is None or sum(S[k]) > sum(S[best])): best = k
    if best is None: return 'novel'
    return 'slip' if n <= 0.2 * sum(S[best]) else 'alt_site'
both = sorted(g for g in O if g in I)
T = collections.defaultdict(lambda: collections.Counter()); wrong = []
for g in both:
    c = G[g]['chrom']; S = support(g); k = cls[g]
    ni = len(I[g]); cut = set(sorted(O[g], key=lambda j: -sum(S[j]))[:ni])  # shortcut: cut per copy, not pooled; ties broken by sort order
    for lab, J in (('ours', O[g]), ('ours_cut', cut), ('isocon', I[g])):
        for j in J:
            ok = j in truth[c]; T[(k, lab)]['n'] += 1; T[(k, lab)]['ok'] += ok
            if not ok and lab != 'ours':
                kind = slip(g, j, J, S); T[(k, lab)][kind] += 1
                wrong.append((lab, k, G[g]['name'], c, j[0], j[1], S[j][0], S[j][1], motif(g, *j), kind))
print(f"copies modelled by both: {len(both)} (multi {sum(cls[g] == 'multi' for g in both)}, single {sum(cls[g] == 'single' for g in both)}); truth = RefSeq introns")
print("class\tset\tjunctions\tcorrect\tprecision\twrong\tslip\talt_site\tnovel")
for k in ('multi', 'single', 'all'):
    for lab in ('ours', 'ours_cut', 'isocon'):
        c = T[(k, lab)] if k != 'all' else T[('multi', lab)] + T[('single', lab)]
        print(f"{k}\t{lab}\t{c['n']}\t{c['ok']}\t{c['ok'] / max(1, c['n']):.3f}\t{c['n'] - c['ok']}\t{c['slip']}\t{c['alt_site']}\t{c['novel']}")
print("\nwrong junctions (ours at the matched cut, isocon): set class gene chrom donor acceptor primary secondary motif kind")
for w in sorted(wrong): print(*w, sep='\t')
