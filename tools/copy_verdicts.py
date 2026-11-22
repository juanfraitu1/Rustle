#!/usr/bin/env python3
"""A provenance verdict per locus from the reads that distinguish two copies (register 1904; O2 extension).

    copy_verdicts.py ASSEMBLY.gtf BAM OUT_PREFIX [--own 2] [--decisive 3] [--drop] [--threads N]

Loci = gene_id of the GTF; a record belongs to the locus it shares the most exon bases with on its strand (as
eq_classes.py), else to `outside`. Per read, its best minimap2 `ms` at each locus it reaches (primary and secondaries).
A read is DECISIVE for locus A over locus B when its best ms at A is higher than at B: sequencing errors cost the same on
both records, so the difference is a base where the two copies differ.
  SUPPORTED  >= --own reads whose best ms is at this locus alone (or that reach no other locus)
  REFUTED    not supported, and a SUPPORTED locus A has >= --decisive reads decisive for A over it and none the other way
  UNDECIDED  reads reach it and nothing separates it from its partners (for these reads the copies are identical)
Writes OUT_PREFIX.verdicts.tsv (locus, chrom, start, end, strand, verdict, own, reads, refuted_by, for, against, partners)
and, with --drop, OUT_PREFIX.gtf = the assembly without the transcripts of REFUTED loci, the others tagged copy_verdict.
--selftest runs the pure parts."""
import sys, os, re, collections, argparse, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
OUTSIDE = 'outside'

def verdicts(reads, own_min=2, dec_min=3):
    """reads: iterable of {locus: best ms}. -> {locus: (verdict, own, n_reads, refuted_by, n_for, n_against, partners)}"""
    own = collections.Counter(); n = collections.Counter(); dec = collections.Counter(); part = collections.defaultdict(set)
    for r in reads:
        ls = [l for l in r if l != OUTSIDE]
        for l in ls: n[l] += 1
        if not r: continue
        b = max(r.values()); top = [l for l in r if r[l] == b]
        if len(top) == 1 and top[0] != OUTSIDE: own[top[0]] += 1
        for a in ls:
            for c in ls:
                if a < c:
                    part[a].add(c); part[c].add(a)
                    if r[a] > r[c]: dec[(a, c)] += 1
                    elif r[c] > r[a]: dec[(c, a)] += 1
    sup = {l for l in n if own[l] >= own_min}; out = {}
    for l in n:
        if l in sup: out[l] = ('SUPPORTED', own[l], n[l], '.', 0, 0, sorted(part[l])); continue
        best = None
        for a in part[l]:
            if a in sup and dec[(a, l)] >= dec_min and dec[(l, a)] == 0 and (best is None or dec[(a, l)] > dec[(best, l)]): best = a
        if best is not None: out[l] = ('REFUTED', own[l], n[l], best, dec[(best, l)], 0, sorted(part[l]))
        else: out[l] = ('UNDECIDED', own[l], n[l], '.', sum(dec[(a, l)] for a in part[l]), sum(dec[(l, a)] for a in part[l]), sorted(part[l]))
    return out

G = {}
def scan(c):
    import pysam
    from gtf_names import overlap_genes
    import eq_classes as E
    gex = G['gex']; R = collections.defaultdict(dict)
    idx = {s: ([a for a, _, _ in gex[(c, s)]], max(b - a for a, b, _ in gex[(c, s)])) for s in '+-' if (c, s) in gex}
    for r in pysam.AlignmentFile(G['bam']).fetch(c):
        if r.is_unmapped or r.is_supplementary or not r.has_tag('ms'): continue
        s = E.tx_strand(r)
        ov = overlap_genes([(a + 1, b) for a, b in r.get_blocks()], gex[(c, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else OUTSIDE
        m = r.get_tag('ms'); d = R[r.query_name]
        if m > d.get(g, -1): d[g] = m
    return dict(R)

def selftest():
    rd = ([{'A': 100}] * 5 + [{'A': 100, 'B': 100}] * 6 + [{'A': 100, 'B': 95}] * 4      # B: tied reads only, 4 reads say A -> REFUTED
          + [{'C': 80, 'D': 80}] * 7                                                     # C, D identical -> UNDECIDED
          + [{'E': 90, 'F': 90}] * 3 + [{'E': 90, 'F': 88}] * 3 + [{'E': 85, 'F': 90}]   # F has one read of its own side -> not refuted
          + [{'H': 50, OUTSIDE: 60}] * 4 + [{'I': 70}])                                  # H: beaten by a place with no locus; I: one read
    v = verdicts(rd)
    assert v['A'][0] == 'SUPPORTED' and v['B'][:5] == ('REFUTED', 0, 10, 'A', 4), v['B']
    assert v['C'][0] == v['D'][0] == 'UNDECIDED' and v['C'][6] == ['D']
    assert v['E'][0] == 'SUPPORTED' and v['F'][0] == 'UNDECIDED' and v['F'][4:6] == (3, 1), v['F']
    assert v['H'][0] == 'UNDECIDED' and v['I'][0] == 'UNDECIDED' and OUTSIDE not in v
    assert verdicts(rd, dec_min=5)['B'][0] == 'UNDECIDED'
    print('copy_verdicts selftest ok')

def main():
    if '--selftest' in sys.argv: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('out'); ap.add_argument('--own', type=int, default=2)
    ap.add_argument('--decisive', type=int, default=3); ap.add_argument('--drop', action='store_true'); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam, eq_classes as E
    meta, gex = E.load_loci(a.gtf); G.update(gex=gex, bam=a.bam)
    refs = [c for c in pysam.AlignmentFile(a.bam).references]
    R = collections.defaultdict(dict)
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        for part in pool.imap_unordered(scan, refs, chunksize=1):
            for q, d in part.items():
                t = R[q]
                for g, m in d.items():
                    if m > t.get(g, -1): t[g] = m
    v = verdicts(R.values(), a.own, a.decisive)
    with open(a.out + '.verdicts.tsv', 'w') as w:
        w.write('locus\tchrom\tstart\tend\tstrand\tverdict\town\treads\trefuted_by\tfor_others\tfor_it\tpartners\n')
        for l in sorted(v, key=lambda l: (meta[l][0], meta[l][1], l)):
            x = v[l]; w.write('\t'.join(map(str, [l, *meta[l][:4], *x[:6], ','.join(x[6][:20]) or '.'])) + '\n')
    c = collections.Counter(x[0] for x in v.values())
    print(f"copy_verdicts: {len(v)} loci with reads: SUPPORTED {c['SUPPORTED']}, REFUTED {c['REFUTED']}, UNDECIDED {c['UNDECIDED']} ({sum(1 for x in v.values() if x[0] == 'UNDECIDED' and x[6])} with partners)")
    if a.drop:
        n = 0
        with open(a.out + '.gtf.tmp', 'w') as w:
            for l in open(a.gtf):
                f = l.rstrip('\n').split('\t')
                g = re.search(r'gene_id "([^"]+)"', f[8]) if len(f) > 8 else None
                x = v.get(g.group(1)) if g else None
                if x and x[0] == 'REFUTED': n += f[2] == 'transcript'; continue
                if x and f[2] == 'transcript': l = l.rstrip('\n').rstrip(';') + f'; copy_verdict "{x[0]}";\n'
                w.write(l)
        os.replace(a.out + '.gtf.tmp', a.out + '.gtf'); print(f'copy_verdicts: {n} transcripts of REFUTED loci dropped -> {a.out}.gtf')

if __name__ == '__main__':
    main()
