#!/usr/bin/env python3
"""O2's certificate on the reads of secondary-only genes (docs/PREREG_s1_o2_certificate_2026-10-12.md).

usage: s1_o2.py --bam B --fasta F --molecules M.tsv --strata PREFIX --ref REF.gtf[.gz] --exclude CHROM
                --op-gtf GTF --op-tmap TMAP --out PREFIX

Every GOOD secondary record (AS >= 0.98 x the molecule's best) overlapping an S1 gene's exon (copy B) is paired with the molecule's
primary record (copy A). Columns = read positions aligned (=/X) in both whose reference bases differ (read orientation). The
certificate is O2's (`bench/asj_o2/asj_o2.py`: eps 0.003 / 3, alpha 1e-3, strict maximum, abstain otherwise)."""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'asj_o2'))
from asj_o2 import ALPHA, EPS, upper_tail  # noqa: E402

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
QUERY_OPS = {0, 1, 4, 5, 7, 8}  # M I S H = X consume the read (H counts toward the original read length)
OPCODE = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}


def aligned_columns(pos, cigar, reverse):
    """(0-based reference start, pysam cigartuples, reverse) -> {original read position: (reference position, matched)} for =/X columns"""
    total = sum(l for op, l in cigar if op in QUERY_OPS)
    cols, q, r = {}, 0, pos
    for op, l in cigar:
        if op in (7, 8):
            for i in range(l):
                o = q + i if not reverse else total - 1 - (q + i)
                cols[o] = (r + i, op == 7)
            q += l; r += l
        elif op == 0:  # M: match / mismatch unknown without the base; not used as a column
            q += l; r += l
        elif op in (1, 4, 5):
            q += l
        elif op in (2, 3):
            r += l
    return cols


def certify_pair(cols_a, cols_b, base_a, base_b):
    """-> (verdict 'A' | 'B' | None, k_A, k_B, columns, p). base_x(o) = copy x's reference base at read position o (read orientation)."""
    ka = kb = m = 0
    for o in cols_a.keys() & cols_b.keys():
        if base_a(o) != base_b(o):
            m += 1
            ka += cols_a[o][1]
            kb += cols_b[o][1]
    if ka == kb:
        return None, ka, kb, m, 1.0
    best, k = ('A', ka) if ka > kb else ('B', kb)
    p = upper_tail(k, m, EPS)
    return (best if p < ALPHA else None), ka, kb, m, p


def eichler_verdict(as_b, as_primary, best, second, margin=10):
    """Eichler's AS-margin rule for one (copy B record, primary) pair: the molecule goes to its best alignment iff no other alignment
    scores within `margin` (best - second >= margin; second < 0 = no other alignment), else it is discarded."""
    if second >= 0 and best - second < margin:
        return 'discard'
    if as_b == best:
        return 'B'
    if as_primary == best:
        return 'A'
    return 'other'


def gene_class(to_b, to_a):
    if to_b >= 2:
        return 'supported'
    if to_b == 0 and to_a >= 2:
        return 'contradicted'
    return 'undecidable'


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def parse_cigar(s):
    return [(OPCODE[op], int(n)) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', s)]


def read_primaries(lines):
    """{read: (chrom, 0-based start, reverse, cigar)} from SAM lines of primary records, the first record per read kept."""
    prim = {}
    for ln in lines:
        f = ln.split('\t', 6)
        prim.setdefault(f[0], (f[2], int(f[3]) - 1, bool(int(f[1]) & 16), parse_cigar(f[5])))
    return prim


def main():
    import pysam
    ap = argparse.ArgumentParser()
    for a in ('--bam', '--fasta', '--molecules', '--strata', '--ref', '--exclude', '--op-gtf', '--op-tmap', '--out'):
        ap.add_argument(a, required=True)
    ap.add_argument('--stratum', default='S1', choices=['S1', 'S2'])
    ap.add_argument('--names-only', action='store_true', help='write OUT.names.txt and stop (the primaries scan is split elsewhere)')
    ap.add_argument('--primaries', default=None, help='SAM lines of the molecules\' primary records, instead of scanning the BAM')
    args = ap.parse_args()
    s1 = set()
    for ln in open(args.strata + '.genes.tsv'):
        f = ln.rstrip('\n').split('\t')
        if f[-1] == args.stratum:
            s1.add((f[0], f[1]))
    ex = collections.defaultdict(list)
    for ln in opener(args.ref):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon' or f[0] == args.exclude:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        g = re.search(r'gene_id "([^"]*)"', f[8]).group(1) or t
        if (f[0], g) in s1:
            ex[(f[0], g)].append((int(f[3]) - 1, int(f[4])))
    for k, v in ex.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        ex[k] = m
    best = {}
    eich = {}  # name -> (best AS, second AS, primary AS)
    for ln in open(args.molecules):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 4 and int(f[3]) >= 2:
            best[f[0]] = int(f[1])
            eich[f[0]] = (int(f[1]), int(f[2]), int(f[5]) if len(f) > 5 and f[5].lstrip('-').isdigit() else -1)
    print(f"[s1_o2] {len(s1):,} {args.stratum} genes, {len(ex):,} with exons; {len(best):,} multi-record molecules", file=sys.stderr)
    bam = pysam.AlignmentFile(args.bam)
    secs = {}  # (gene, qname, ref, pos, cigar) -> (cigartuples, reverse, AS, tied)
    for gk, exons in ex.items():
        for s, e in exons:
            for rec in bam.fetch(gk[0], s, e):
                if not rec.is_secondary or rec.is_supplementary or rec.is_unmapped:
                    continue
                b = best.get(rec.query_name)
                if b is None or not rec.has_tag('AS'):
                    continue
                a = rec.get_tag('AS')
                if a < 0.98 * b:
                    continue
                if not any(bs < e and be > s for bs, be in rec.get_blocks()):
                    continue
                secs[(gk, rec.query_name, rec.reference_name, rec.reference_start, rec.cigarstring)] = (
                    rec.cigartuples, rec.is_reverse, a, a == b)
    names = sorted({k[1] for k in secs})
    print(f"[s1_o2] {len(secs):,} good secondary records of {len(names):,} molecules", file=sys.stderr)
    nf = args.out + '.names.txt'
    with open(nf, 'w') as o:
        o.write('\n'.join(names) + '\n')
    if args.names_only:
        print(f"[s1_o2] names only -> {nf}", file=sys.stderr)
        return
    if args.primaries:
        with open(args.primaries) as fh:
            prim = read_primaries(fh)
    else:
        p = subprocess.Popen(['samtools', 'view', '-F', '2308', '-N', nf, args.bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        prim = read_primaries(p.stdout)
        p.wait()
    fa = pysam.FastaFile(args.fasta)

    def bases(chrom, cols, reverse):
        lo = min(r for r, _ in cols.values()); hi = max(r for r, _ in cols.values()) + 1
        seq = fa.fetch(chrom, lo, hi).upper()
        if reverse:
            return lambda o: seq[cols[o][0] - lo].translate(COMP)
        return lambda o: seq[cols[o][0] - lo]

    rows = []
    for (gk, name, ref, pos, _), (cig, rev, a, tied) in secs.items():
        b_, s_, pa_ = eich[name]
        ev = eichler_verdict(a, pa_, b_, s_)
        pr = prim.get(name)
        if pr is None:
            rows.append((gk, name, tied, None, 0, 0, 0, 1.0, 'no_primary', ev)); continue
        ca = aligned_columns(pr[1], pr[3], pr[2]); cb = aligned_columns(pos, cig, rev)
        if not ca or not cb:
            rows.append((gk, name, tied, None, 0, 0, 0, 1.0, 'no_columns', ev)); continue
        v, ka, kb, m, pv = certify_pair(ca, cb, bases(pr[0], ca, pr[2]), bases(ref, cb, rev))
        rows.append((gk, name, tied, v, ka, kb, m, pv, 'tested', ev))
    with open(args.out + '.records.tsv', 'w') as o:
        o.write('chrom\tgene_id\tread\ttied\tverdict\tk_A\tk_B\tcolumns\tp\tstatus\teichler\n')
        for gk, name, tied, v, ka, kb, m, pv, st, ev in rows:
            o.write(f"{gk[0]}\t{gk[1]}\t{name}\t{int(tied)}\t{v or '-'}\t{ka}\t{kb}\t{m}\t{pv:.3g}\t{st}\t{ev}\n")
    per = collections.defaultdict(lambda: [0, 0, 0])  # gene -> [to B, to A, abstain]
    for gk, _, _, v, *_ in rows:
        per[gk][0 if v == 'B' else (1 if v == 'A' else 2)] += 1
    cls = {gk: gene_class(c[0], c[1]) for gk, c in per.items()}
    for gk in ex:
        cls.setdefault(gk, 'undecidable')
    eper = collections.defaultdict(lambda: [0, 0, 0])  # gene -> [to B, elsewhere (A or other), discard]
    for r in rows:
        eper[r[0]][0 if r[9] == 'B' else (1 if r[9] in ('A', 'other') else 2)] += 1
    ecls = {gk: gene_class(eper[gk][0], eper[gk][1]) if gk in eper else 'undecidable' for gk in cls}
    with open(args.out + '.genes.tsv', 'w') as o:
        o.write('chrom\tgene_id\tto_stratum_copy\tto_primary_copy\tabstain\tclass\teichler_to_stratum_copy\teichler_elsewhere\teichler_discard\teichler_class\n')
        for gk in sorted(cls):
            c = per.get(gk, [0, 0, 0]); e = eper.get(gk, [0, 0, 0])
            o.write(f"{gk[0]}\t{gk[1]}\t{c[0]}\t{c[1]}\t{c[2]}\t{cls[gk]}\t{e[0]}\t{e[1]}\t{e[2]}\t{ecls[gk]}\n")
    # OP transcripts touching S1 genes, by class
    eq = set()
    for ln in open(args.op_tmap):
        f = ln.split('\t')
        if len(f) > 4 and f[2] == '=':
            eq.add(f[4])
    bins = collections.defaultdict(list)
    for gk, exons in ex.items():
        for s, e in exons:
            for b in range(s // 10000, (e - 1) // 10000 + 1):
                bins[(gk[0], b)].append((s, e, gk))
    tx_genes = collections.defaultdict(set)
    for ln in open(args.op_gtf):
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
        for b in range(s // 10000, (e - 1) // 10000 + 1):
            for xs, xe, gk in bins.get((f[0], b), ()):
                if xs < e and xe > s:
                    tx_genes[t].add(gk)
    rank = {'supported': 0, 'contradicted': 1, 'undecidable': 2}
    tx_cls = {t: min((cls[g] for g in gs), key=rank.get) for t, gs in tx_genes.items()}
    # summary
    st = collections.Counter(r[8] for r in rows)
    print(f"records {len(rows):,}: tested {st['tested']:,}, no primary {st['no_primary']:,}, no columns {st['no_columns']:,}")
    for lab, sel in (('tied', True), ('strictly worse', False)):
        sub = [r for r in rows if r[2] == sel and r[8] == 'tested']
        c = collections.Counter(r[3] for r in sub)
        print(f"  {lab:15s} {len(sub):6,} records: to S1 copy {c['B']:6,}  to primary copy {c['A']:6,}  abstain {c[None]:6,}"
              f"  (no distinguishing column {sum(r[6] == 0 for r in sub):,})")
    ec = collections.Counter(r[9] for r in rows)
    print(f"  Eichler (AS margin 10): to {args.stratum} copy {ec['B']:,}  to primary copy {ec['A']:,}  to a third placement {ec['other']:,}  discard {ec['discard']:,}")
    both = [(r[3], r[9]) for r in rows if r[3] in ('A', 'B') and r[9] in ('A', 'B')]
    agree = sum(a == b for a, b in both)
    print(f"  records both rules call A or B: {len(both):,}, agree {agree:,} ({agree / max(1, len(both)):.3f})")
    xt = collections.Counter((cls[g], ecls[g]) for g in cls)
    print("  gene classes O2 (rows) x Eichler (columns): " + '; '.join(
        f"{a}: " + ', '.join(f"{b} {xt[(a, b)]}" for b in ('supported', 'contradicted', 'undecidable')) for a in ('supported', 'contradicted', 'undecidable')))
    with_tx = {g for gs in tx_genes.values() for g in gs}
    for lab, keep in ((f'all {args.stratum} genes', set(cls)), (f'{args.stratum} genes with an OP transcript', with_tx)):
        c = collections.Counter(cls[g] for g in keep)
        print(f"{lab:32s} {len(keep):5,}: supported {c['supported']:5,}  contradicted {c['contradicted']:5,}  undecidable {c['undecidable']:5,}")
    for k in ('supported', 'contradicted', 'undecidable'):
        ts = [t for t, c in tx_cls.items() if c == k]
        print(f"  OP transcripts in {k:12s} genes: {len(ts):5,}, reference-exact {sum(t in eq for t in ts):4,}")


if __name__ == '__main__':
    main()
