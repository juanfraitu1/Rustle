#!/usr/bin/env python3
"""Member-level recall: of the truth family members that are EXPRESSED by the reads, how many does the run find, and where
are the rest lost? Descriptive.

  member_funnel.py [--species ...] --gtf RUN.gtf --bam BAM --hits PAIRSCORE_HITS [--threads N] [--out TSV] COPIES.tsv
  member_funnel.py --selftest
Truth members = protein-coding genes in >= 1 Compara pair at simiiformes (famtruth). Expression is read-based, not
model-based: P = >= 2 molecules with a PRIMARY record overlapping the gene's exons, T = >= 2 with ANY record there.
A member is FINDABLE when it is P-expressed and >= 1 of its truth partners is too. Each findable member lands in the
first stage that fits:
  placed       it is in a hit pair of pairscore (PAIRSCORE_HITDUMP of the same table): beside a true partner
  wrong_family it maps to a copy of a family, but no true partner is there
  no_family    it has a model (famtruth.expressed_from_gtf: a transcript whose best gene it is) but no family copy
  no_model     no transcript has it as its best gene
Also counted: T-only members (seen only through secondary records). env: as pairscore.py."""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def stage(g, hits, in_family, modelled):
    return 'placed' if g in hits else 'wrong_family' if g in in_family else 'no_family' if g in modelled else 'no_model'

def _count(job):
    import pysam
    bam, chrom, items = job; B = pysam.AlignmentFile(bam); out = {}
    if chrom not in B.references: return {g: (0, 0) for g, _ in items}
    for g, ex in items:
        prim, anyr = set(), set()
        for r in B.fetch(chrom, ex[0][0], ex[-1][1]):
            if r.is_supplementary or r.is_unmapped: continue
            if any(s < y and x < e for s, e in r.get_blocks() for x, y in ex):
                anyr.add(r.query_name)
                if not r.is_secondary: prim.add(r.query_name)
        out[g] = (len(prim), len(anyr))
    return out

def main(a):
    import multiprocessing, famtruth as F
    opt = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    G = F.genes_from_argv(a, opt('--species', 'human')); idx = F.exon_index(G)
    truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
    truth = {tuple(sorted(p)) for p in truth if G[p[0]]['name'] != G[p[1]]['name']}
    partners = collections.defaultdict(set)
    for x, y in truth: partners[x].add(y); partners[y].add(x)
    bych = collections.defaultdict(list)
    for g in partners: bych[G[g]['chrom']].append((g, sorted(G[g]['exons'])))
    with multiprocessing.Pool(int(opt('--threads', 8))) as P:
        cnt = {}
        for d in P.imap_unordered(_count, [(opt('--bam'), c, v) for c, v in bych.items()]): cnt.update(d)
    Pexp = {g for g, (p, t) in cnt.items() if p >= 2}; Texp = {g for g, (p, t) in cnt.items() if t >= 2}
    findable = {g for g in Pexp if partners[g] & Pexp}
    modelled = F.expressed_from_gtf(opt('--gtf'), idx)
    fams, _, _ = F.load_catalog(a[-1], G, idx); in_family = {g for v in fams.values() for g in v}
    hits = {g for l in open(opt('--hits')) for g in l.split('\t')[:2]}
    st = {g: stage(g, hits, in_family, modelled) for g in findable}
    c = collections.Counter(st.values()); n = len(findable)
    print(f'truth members {len(partners)}; P-expressed {len(Pexp)}; T-only (secondary records only) {len(Texp - Pexp)}; findable (a truth partner is P-expressed too) {n}')
    for k in ('placed', 'wrong_family', 'no_family', 'no_model'): print(f'{k}\t{c[k]}\t{c[k] / max(1, n):.3f}')
    tonly = Texp - Pexp
    print(f'T-only members: {len(tonly)}; modelled {len(tonly & modelled)}; in a family {len(tonly & in_family)}; placed {len(tonly & hits)}')
    if opt('--out'):
        with open(opt('--out'), 'w') as w:
            w.write('gene\tname\tchrom\tprimary\tany\tstage\tpartners_expressed\n')
            for g in sorted(findable, key=lambda g: (st[g], G[g]['name'])):
                w.write(f"{g}\t{G[g]['name']}\t{G[g]['chrom']}\t{cnt[g][0]}\t{cnt[g][1]}\t{st[g]}\t{len(partners[g] & Pexp)}\n")

def selftest():
    assert stage('a', {'a'}, {'a'}, {'a'}) == 'placed' and stage('a', set(), {'a'}, {'a'}) == 'wrong_family'
    assert stage('a', set(), set(), {'a'}) == 'no_family' and stage('a', set(), set(), set()) == 'no_model'
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(sys.argv[1:])
