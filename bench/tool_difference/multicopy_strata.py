#!/usr/bin/env python3
"""Multi-copy strata of reference genes, from the reads alone (no assembly): docs/archive/2026-10/merged/PREREG_inbetween_operating_point_2026-10-11.md.

usage: multicopy_strata.py --bam B --molecules M.tsv --ref REF.gtf[.gz] --exclude CHROM --out PREFIX

Per reference gene (its exons merged), every mapped non-supplementary record whose aligned blocks (M/=/X) overlap an exon by
>= 1 bp is counted as: a PRIMARY; a GOOD secondary (AS >= 0.98 x the molecule's genome-wide best AS, the `good` pool's rule);
a TIED secondary (AS == the molecule's best AS: the molecule's primary was chosen among equally scored placements).
  S1 secondary-only   0 primaries and >= 2 good secondaries
  S2 coin-toss-starved >= 1 primary, >= 2 tied secondaries and at least as many tied secondaries as primaries
Writes PREFIX.genes.tsv (every gene with any record) and PREFIX.transcripts.tsv (multi-exon reference transcripts of S1/S2).
`--molecules` is `as_table`'s output (name, best_as, second_as, n_records, primary_chrom, primary_as)."""
import argparse, collections, gzip, re, subprocess, sys

BIN = 10_000

def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

ap = argparse.ArgumentParser()
for a in ('--bam', '--molecules', '--ref', '--exclude', '--out'):
    ap.add_argument(a, required=True)
args = ap.parse_args()

# reference genes: merged exons; multi-exon transcripts per gene
gene_ex = collections.defaultdict(list); tx_ex = collections.defaultdict(list); tx_gene = {}; gene_chrom = {}
for ln in opener(args.ref):
    f = ln.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon' or f[0] == args.exclude:
        continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    g = re.search(r'gene_id "([^"]*)"', f[8]).group(1) or t  # RefSeq gene-level features carry gene_id "": key them by their id
    iv = (int(f[3]) - 1, int(f[4]))
    gene_ex[(f[0], g)].append(iv); tx_ex[t].append(iv); tx_gene[t] = (f[0], g)
index = collections.defaultdict(list)  # (chrom, bin) -> [(start, end, gene_key)]
for gk, ivs in gene_ex.items():
    ivs.sort(); merged = []
    for s, e in ivs:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    gene_ex[gk] = merged
    for s, e in merged:
        for b in range(s // BIN, (e - 1) // BIN + 1):
            index[(gk[0], b)].append((s, e, gk))
print(f"[strata] {len(gene_ex):,} reference genes off {args.exclude}", file=sys.stderr)

best = {}
with open(args.molecules) as fh:
    for ln in fh:
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 4 and int(f[3]) >= 2:  # only molecules with a secondary record can have a good / tied secondary
            best[f[0]] = int(f[1])
print(f"[strata] {len(best):,} multi-record molecules", file=sys.stderr)

cig = re.compile(r'(\d+)([MIDNSHP=X])')

def count_chrom(chrom):
    """Per-gene [primary, good secondary, tied secondary] over the records of one chromosome (exact per-record rule)."""
    cnt = collections.defaultdict(lambda: [0, 0, 0]); n = 0
    p = subprocess.Popen(['samtools', 'view', '-F', '2052', args.bam, chrom], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for ln in p.stdout:
        n += 1
        f = ln.split('\t', 11)
        flag = int(f[1]); sec = flag & 256; tied = False
        if sec:
            b = best.get(f[0])
            i = ln.find('\tAS:i:')
            if b is None or i < 0:
                continue
            j = ln.find('\t', i + 6)
            a = int(ln[i + 6:j if j > 0 else len(ln)].strip())
            if a < 0.98 * b:
                continue
            tied = a == b
        pos = int(f[3]) - 1; hit = set()
        for l, op in cig.findall(f[5]):
            l = int(l)
            if op in 'M=X':
                s0, e0 = pos, pos + l
                for bn in range(s0 // BIN, (e0 - 1) // BIN + 1):
                    for xs, xe, gk in index.get((chrom, bn), ()):
                        if xs < e0 and xe > s0:
                            hit.add(gk)
                pos += l
            elif op in 'DN':
                pos += l
        for gk in hit:
            c = cnt[gk]
            if not sec:
                c[0] += 1
            else:
                c[1] += 1
                c[2] += tied
    p.wait()
    return chrom, n, dict(cnt)

import multiprocessing
chroms = sorted({gk[0] for gk in gene_ex})
cnt = {}; n = 0
with multiprocessing.get_context('fork').Pool(4) as pool:
    for chrom, k, c in pool.imap_unordered(count_chrom, chroms):
        n += k; cnt.update(c)
        print(f"[strata] {chrom}: {k:,} records", file=sys.stderr)
strat = {}
with open(args.out + '.genes.tsv', 'w') as o:
    o.write('chrom\tgene_id\tprimary\tgood_secondary\ttied_secondary\tstratum\n')
    for gk, (pr, gs, ts) in sorted(cnt.items()):
        s = 'S1' if pr == 0 and gs >= 2 else ('S2' if pr >= 1 and ts >= 2 and ts >= pr else '-')
        strat[gk] = s
        o.write(f"{gk[0]}\t{gk[1]}\t{pr}\t{gs}\t{ts}\t{s}\n")
k = collections.Counter(strat.values())
with open(args.out + '.transcripts.tsv', 'w') as o:
    o.write('transcript_id\tchrom\tgene_id\tstratum\tn_exons\n')
    nt = collections.Counter()
    for t, ex in tx_ex.items():
        gk = tx_gene[t]; s = strat.get(gk, '-')
        if s != '-' and len(ex) > 1:
            o.write(f"{t}\t{gk[0]}\t{gk[1]}\t{s}\t{len(ex)}\n"); nt[s] += 1
print(f"[strata] {n:,} records; genes S1 {k['S1']:,}, S2 {k['S2']:,}; multi-exon reference transcripts S1 {nt['S1']:,}, S2 {nt['S2']:,}", file=sys.stderr)
