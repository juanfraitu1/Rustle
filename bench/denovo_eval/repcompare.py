#!/usr/bin/env python3
"""LOCUS REPRESENTATIVE vs StringTie / FLAIR (2026-10-30, user: "ensure that the assembly part is as similar as possible to stringtie
and flair? I feel that the chosen representatives are still different so it's time to measure and pick a new default").
A de novo locus = one gene_id group of our GTF (what `mcl_families --from-gtf` uses); the family node is ONE representative transcript
of it. Rules compared, each picking one transcript per locus:
  most-reads     reads attribute, tie longer span                         (shipped default)
  most-junctions junction count, tie reads, tie span                      (`--representative most-junctions`, opt-in)
  most-cov       StringTie-style per-base coverage (`cov` attribute), tie reads
  most-tpm       `TPM` attribute, tie reads
  longest        spliced length (sum of exons), tie reads
  reads-x-junc   reads x (junctions + 1), tie reads
Per locus the StringTie TOP transcript = the StringTie transcript overlapping the locus span on the same strand with the highest TPM
(tie cov); FLAIR has no abundance in its isoforms GTF, so FLAIR is a SET. Measures over loci with >= 2 transcripts (where rules differ):
  st_exact      the rule's pick has the same intron chain as the StringTie top (spliced loci only; both spliced)
  st_compat     the pick's introns are a subset or superset of the StringTie top's
  flair_member  the pick's intron chain exists among the FLAIR isoforms at the locus
  ref_exact     the pick's chain equals an annotated transcript's chain
  st_top_in_locus  the StringTie top's chain is one of OUR transcripts at the locus (ceiling for st_exact)
Also the same numbers over ALL loci (single-transcript loci included) and the agreement between rules (same pick).
  repcompare.py OUR.gtf STRINGTIE.gtf[.gz] FLAIR.gtf[.gz] REFSEQ.gtf[.gz] [--out OUT.tsv]
"""
import sys, re, gzip, collections, bisect, csv
a = sys.argv[1:]
OUT = a[a.index('--out') + 1] if '--out' in a else None
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or a[i - 1] != '--out')]
OURS, ST, FL, REF = pos[:4]
def op(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
def load(path, want=('reads', 'cov', 'TPM')):
    """transcript_id -> dict(chrom, strand, gene, exons[(s,e) 0-based], attrs)"""
    tx = {}
    for line in op(path):
        if line[0] == '#': continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if t not in tx:
            g = re.search(r'gene_id "([^"]+)"', f[8])
            tx[t] = {'chrom': f[0], 'strand': f[6] if f[6] in '+-' else '+', 'gene': g.group(1) if g else t, 'exons': [], 'at': {}}
        if f[2] == 'transcript':  # the abundance attributes live on the transcript line
            for k in want:
                m = re.search(r'(?:^|; )' + k + r' "([^"]+)"', f[8])
                if m: tx[t]['at'][k] = float(m.group(1))
        else: tx[t]['exons'].append((int(f[3]) - 1, int(f[4])))
    for d in tx.values(): d['exons'].sort()
    return tx
def introns(ex): return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1) if ex[i + 1][0] - ex[i][1] >= 50)
def junctions(ex): return len(introns(ex))
def span(ex): return ex[-1][1] - ex[0][0]
def slen(ex): return sum(e - s for s, e in ex)
ours = load(OURS); st = load(ST); fl = load(FL, ()); ref = load(REF, ())
print(f"ours {len(ours)} transcripts; StringTie {len(st)}; FLAIR {len(fl)}; RefSeq {len(ref)}", file=sys.stderr)
class Idx:
    def __init__(self, tx):
        by = collections.defaultdict(list)
        for t, d in tx.items(): by[(d['chrom'], d['strand'])].append((d['exons'][0][0], d['exons'][-1][1], t))
        self.v = {k: sorted(v) for k, v in by.items()}; self.s = {k: [x[0] for x in v] for k, v in self.v.items()}
        self.m = {k: max(x[1] - x[0] for x in v) for k, v in self.v.items()}
    def over(self, chrom, strand, s, e):
        k = (chrom, strand); v = self.v.get(k, []); i = bisect.bisect_left(self.s.get(k, []), s - self.m.get(k, 0)); out = []
        while i < len(v) and v[i][0] < e:
            if v[i][1] > s: out.append(v[i][2])
            i += 1
        return out
sti, fli, refi = Idx(st), Idx(fl), Idx(ref)
loci = collections.defaultdict(list)
for t, d in ours.items(): loci[d['gene']].append(t)
RULES = {
    'most-reads': lambda d: (d['at'].get('reads', 0), span(d['exons'])),
    'most-junctions': lambda d: (junctions(d['exons']), d['at'].get('reads', 0), span(d['exons'])),
    'most-cov': lambda d: (d['at'].get('cov', 0), d['at'].get('reads', 0)),
    'most-tpm': lambda d: (d['at'].get('TPM', 0), d['at'].get('reads', 0)),
    'longest': lambda d: (slen(d['exons']), d['at'].get('reads', 0)),
    'reads-x-junc': lambda d: (d['at'].get('reads', 0) * (junctions(d['exons']) + 1), d['at'].get('reads', 0)),
}
# StringTie-like rules need the locus: a chain is MAXIMAL when its intron set is not a proper subset of another transcript's at the locus;
#   maximal-reads   most reads among maximal chains (fragments never represent a locus), tie span
#   flow-reads      each transcript's reads are credited to every maximal chain that contains its introns (StringTie's flow in spirit);
#                   the maximal chain with the largest credited sum wins, tie own reads
def locus_rules(ts):
    chains = {t: set(introns(ours[t]['exons'])) for t in ts}
    maximal = [t for t in ts if not any(chains[t] < chains[u] for u in ts if u != t)]
    rd = lambda t: ours[t]['at'].get('reads', 0)
    mr = max(maximal, key=lambda t: (rd(t), span(ours[t]['exons']), t))
    credit = {t: sum(rd(u) for u in ts if chains[u] <= chains[t]) for t in maximal}
    fr = max(maximal, key=lambda t: (credit[t], rd(t), t))
    return {'maximal-reads': mr, 'flow-reads': fr}
ATTR = collections.Counter()
for d in ours.values():
    for k in ('reads', 'cov', 'TPM'):
        if k in d['at']: ATTR[k] += 1
print('attributes present:', dict(ATTR), 'of', len(ours), file=sys.stderr)
DIAG = collections.Counter()
C = collections.Counter(); picks = collections.defaultdict(dict)
rows = []
for g, ts in loci.items():
    d0 = ours[ts[0]]; chrom, strand = d0['chrom'], d0['strand']
    ls, le = min(ours[t]['exons'][0][0] for t in ts), max(ours[t]['exons'][-1][1] for t in ts)
    multi = len(ts) >= 2
    sts = sti.over(chrom, strand, ls, le)
    top = max(sts, key=lambda t: (st[t]['at'].get('TPM', 0), st[t]['at'].get('cov', 0), t)) if sts else None
    topc = introns(st[top]['exons']) if top else None
    flc = {introns(fl[t]['exons']) for t in fli.over(chrom, strand, ls, le)}
    refc = {introns(ref[t]['exons']) for t in refi.over(chrom, strand, ls, le)}
    ourc = {introns(ours[t]['exons']) for t in ts}
    for scope in (['all', 'multi'] if multi else ['all']):
        C[(scope, 'loci')] += 1
        if top is not None and topc: C[(scope, 'st_top_spliced')] += 1
        if top is not None and topc and topc in ourc: C[(scope, 'st_top_in_locus')] += 1
    lr = locus_rules(ts)
    if multi and top is not None and topc and topc in ourc:  # where the StringTie top is ours but most-reads picks another chain: why?
        mr = max(ts, key=lambda t: (RULES['most-reads'](ours[t]), t)); mrc = introns(ours[mr]['exons'])
        if mrc != topc:
            tt = next(t for t in ts if introns(ours[t]['exons']) == topc)
            DIAG['differ'] += 1
            if set(mrc) < set(topc): DIAG['pick_is_fragment_of_top'] += 1
            elif set(topc) < set(mrc): DIAG['top_is_fragment_of_pick'] += 1
            else: DIAG['incompatible'] += 1
            if ours[tt]['at'].get('reads', 0) == ours[mr]['at'].get('reads', 0): DIAG['tie_on_reads'] += 1
            if len(topc) > len(mrc): DIAG['top_has_more_junctions'] += 1
            if lr['flow-reads'] == tt: DIAG['flow_reads_picks_top'] += 1
            if lr['maximal-reads'] == tt: DIAG['maximal_reads_picks_top'] += 1
    for r in list(RULES) + ['maximal-reads', 'flow-reads']:
        p = lr[r] if r in lr else max(ts, key=lambda t: (RULES[r](ours[t]), t)); picks[g][r] = p
        pc = introns(ours[p]['exons'])
        for scope in (['all', 'multi'] if multi else ['all']):
            if top is not None and topc and pc:
                C[(scope, r, 'st_pairs')] += 1
                if pc == topc: C[(scope, r, 'st_exact')] += 1
                if set(pc) <= set(topc) or set(topc) <= set(pc): C[(scope, r, 'st_compat')] += 1
            if pc:
                C[(scope, r, 'spliced')] += 1
                if pc in flc: C[(scope, r, 'flair_member')] += 1
                if pc in refc: C[(scope, r, 'ref_exact')] += 1
    rows.append((g, chrom, strand, len(ts), top, picks[g]))
def pr(scope):
    n = C[(scope, 'loci')]
    print(f"\n## {scope} loci: {n}; StringTie top spliced at {C[(scope, 'st_top_spliced')]}; StringTie top chain is one of ours at {C[(scope, 'st_top_in_locus')]} ({C[(scope, 'st_top_in_locus')] / max(1, C[(scope, 'st_top_spliced')]):.3f} of spliced tops)")
    print("rule\tst_pairs\tst_exact\tst_compat\tspliced\tflair_member\tref_exact")
    for r in list(RULES) + ['maximal-reads', 'flow-reads']:
        sp = C[(scope, r, 'st_pairs')]; s = C[(scope, r, 'spliced')]
        print(f"{r}\t{sp}\t{C[(scope, r, 'st_exact')] / max(1, sp):.3f}\t{C[(scope, r, 'st_compat')] / max(1, sp):.3f}\t{s}\t{C[(scope, r, 'flair_member')] / max(1, s):.3f}\t{C[(scope, r, 'ref_exact')] / max(1, s):.3f}")
pr('all'); pr('multi')
print('\n## multi loci where the StringTie top chain is ours but most-reads picks another chain:', dict(DIAG))
names = list(RULES) + ['maximal-reads', 'flow-reads']; multi_g = [g for g, ts in loci.items() if len(ts) >= 2]
print("\n## same pick as most-reads on multi-transcript loci")
for r in names[1:]:
    same = sum(1 for g in multi_g if picks[g][r] == picks[g]['most-reads'])
    print(f"{r}\t{same}/{len(multi_g)}\t{same / max(1, len(multi_g)):.3f}")
if OUT:
    with open(OUT, 'w') as w:
        w.write('locus\tchrom\tstrand\tn_tx\tst_top\t' + '\t'.join(names) + '\n')
        for g, chrom, strand, n, top, pk in rows: w.write(f"{g}\t{chrom}\t{strand}\t{n}\t{top or ''}\t" + '\t'.join(pk[r] for r in names) + '\n')
