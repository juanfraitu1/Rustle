#!/usr/bin/env python3
"""Audit of gene-joining models (2026-11-14, descriptive, no bar): transcripts of one assembled GTF whose exons touch
>= 2 coding RefSeq genes that share no exon base (same strand), by class, shape and F1 bridge tag, and the joined genes
left without a model of their own in that GTF.
  class   annot_rt = a RefSeq coding gene itself spans two of the joined genes (annotated readthrough / fusion gene);
          truth    = a joined pair is a truth paralog pair (famtruth, primates level);
          named    = no truth pair, but two joined genes share a name stem (ZNF, PRM: relatives by name, not proven copies);
          loc      = the shared stem is LOC (uncharacterised neighbours: no evidence either way);
          other    = none of these (readthrough between unrelated neighbours or a merge artefact)
  junc2   = the transcript carries an annotated intron of two joined genes that share no exon (a splice-level join;
            the rest touch the second gene by exon overlap only, often a long end)
  tagged  = the transcript carries F1's fusion_of attribute
  own     = a joined gene also has a transcript that joins nothing
usage: fusion_audit.py ASSEMBLY.gtf REFSEQ.gtf[.gz] OUT.tsv [FAM.copies.tsv]   (famtruth env vars for the truth class)"""
import collections, csv, gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
from cat_vs_refseq_lib import Idx, merge, share  # noqa: E402

op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
A = lambda k, a: (re.search(k + r' "([^"]+)"', a) or [None, ''])[1]
stem = lambda n: re.sub(r'(L?\d+[A-Z]?|P\d*)$', '', n)  # as tandem_fused.py


def load(p, key):
    ex, meta = collections.defaultdict(list), {}
    for l in op(p):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'transcript': meta.setdefault(A(key, f[8]), [f[0], f[6], A('gene_id', f[8]), int(A('reads', f[8]) or 0), set(), 'fusion_of "' in f[8]])
        if f[2] != 'exon': continue
        k = A(key, f[8]); ex[k].append((int(f[3]) - 1, int(f[4])))
        meta.setdefault(k, [f[0], f[6], A('gene_name', f[8]) or A('gene_id', f[8]), 0, set(), False])[4].add(A('transcript_id', f[8]))
    return {k: (meta[k][0], meta[k][1], merge(v)) for k, v in ex.items()}, meta


def classify(hits, R, name, ppairs):
    """(kind, class, joined genes) of one transcript's RefSeq hits; kind as merge_audit.py."""
    pairs = [(a, b) for i, a in enumerate(hits) for b in hits[i + 1:]]
    dis = [(a, b) for a, b in pairs if not share(R[a][2], R[b][2])]
    if not dis: return ('shared-exon cluster' if pairs else 'single'), None, []
    kind = 'fused' if len(dis) == len(pairs) else 'mixed'
    if any(share(R[h][2], R[a][2]) and share(R[h][2], R[b][2]) for a, b in dis for h in hits if h not in (a, b)): cls = 'annot_rt'
    elif any(frozenset((name[a], name[b])) in ppairs for a, b in dis): cls = 'truth'
    elif any(stem(name[a]) == stem(name[b]) != 'LOC' for a, b in dis): cls = 'named'
    elif any(stem(name[a]) == stem(name[b]) for a, b in dis): cls = 'loc'
    else: cls = 'other'
    return kind, cls, sorted({g for p in dis for g in p})


def main(gtf, rs, out, copies=None):
    R, rmeta = load(rs, 'gene_id')
    coding = {g for g in R if any(t.startswith(('rna-NM_', 'rna-XM_')) for t in rmeta[g][4])}
    R = {g: R[g] for g in coding}; name = {g: rmeta[g][2].removeprefix('gene-') for g in R}; ri = Idx(R)
    tex = collections.defaultdict(list); gin = collections.defaultdict(set)
    for l in op(rs):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon' and A('gene_id', f[8]) in R: tex[(A('gene_id', f[8]), A('transcript_id', f[8]))].append((int(f[3]) - 1, int(f[4])))
    for (g, _), e in tex.items(): e.sort(); gin[g].update((a[1], b[0]) for a, b in zip(e, e[1:]))
    ppairs = set(); univ = set()
    try:
        import famtruth as F
        G = F.load_cat(); pairs, u = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['primates']))
        ppairs = {frozenset((G[a]['name'], G[b]['name'])) for a, b in pairs}; univ = {G[g]['name'] for g in u}
    except (Exception, SystemExit) as e:  # gorilla / no env: name stems only
        print(f'[fusion_audit] no paralog truth ({type(e).__name__}): paralog class = shared name stem only', file=sys.stderr)
    T, tm = load(gtf, 'transcript_id'); CL = ('annot_rt', 'truth', 'named', 'loc', 'other')
    n = collections.Counter(); rows = []; own = set(); joined = collections.defaultdict(set); loci = collections.defaultdict(set)
    for t, (c, s, ex) in T.items():
        hits = ri.hits(c, s, ex); kind, cls, genes = classify(hits, R, name, ppairs)
        if cls is None: own.update(hits); continue
        ti = {(a[1], b[0]) for a, b in zip(ex, ex[1:])}; jg = [g for g in genes if gin[g] & ti]
        j2 = any(not share(R[a][2], R[b][2]) for i, a in enumerate(jg) for b in jg[i + 1:]); tag = tm[t][5]
        far = max((b[0] - a[1] for a, b in zip(ex, ex[1:])), default=0) >= 500000
        for sh in ['mono' if len(ex) == 1 else 'spliced'] + ['junc2'] * j2:
            n[(cls, sh, 'tx')] += 1; n[(cls, sh, 'reads')] += tm[t][3]; n[(cls, sh, 'tagged')] += tag; n[(cls, sh, 'far')] += far
            loci[(cls, sh, tag)].add(tm[t][2])
        for g in genes: joined[g].add((cls, j2))
        rows.append((t, tm[t][2], c, ex[0][0], ex[-1][1], s, cls, kind, len(ex), tm[t][3], int(tag), int(j2), ','.join(name[g] for g in genes)))
    with open(out, 'w') as fh:
        print('transcript\tlocus\tchrom\tstart\tend\tstrand\tclass\tkind\texons\treads\ttagged\tjunc2\tgenes', file=fh)
        for r in sorted(rows, key=lambda r: (r[6], -r[9])): print(*r, sep='\t', file=fh)
    L = lambda f: set().union(*(v for k, v in loci.items() if f(k)))
    print(f'{gtf}: transcripts {len(T)} (F1-tagged {sum(m[5] for m in tm.values())}); coding RefSeq genes {len(R)}; gene-joining transcripts {len(rows)} in {len(L(lambda k: 1))} loci; junc2 {sum(r[11] for r in rows)} in {len(L(lambda k: k[1] == "junc2"))} loci; untagged junc2 loci {len(L(lambda k: k[1] == "junc2" and not k[2]))}')
    print('class\tshape\ttranscripts\treads\tF1_tagged\tintron>=500kb\tloci')
    for cls in CL:
        for sh in ('spliced', 'junc2', 'mono'): print(cls, sh, *(n[(cls, sh, k)] for k in ('tx', 'reads', 'tagged', 'far')), len(L(lambda k: k[:2] == (cls, sh))), sep='\t')
    ab = [g for g in joined if g not in own]
    print(f'joined genes {len(joined)}: with a model of their own {len(joined) - len(ab)}, without {len(ab)} (through a junc2 transcript {sum(any(j for _, j in joined[g]) for g in ab)}); of those in the paralog truth universe {sum(name[g] in univ for g in ab)}')
    print('  genes without their own model by joining class (a gene can be in two): ' + ', '.join(f'{c} {sum(any(x[0] == c for x in joined[g]) for g in ab)} (junc2 {sum(any(x == (c, True) for x in joined[g]) for g in ab)})' for c in CL))
    if copies:
        fams = collections.defaultdict(set)
        for r in csv.DictReader(open(copies), delimiter='\t'): fams[r['family_id']].add(r['gene_id'])
        print(f'catalog families {len(fams)}: with a member locus holding a joining transcript ' + ', '.join(f'{c} {sum(bool(v & L(lambda k: k[0] == c))  for v in fams.values())}' for c in CL) + f'; any {sum(bool(v & L(lambda k: 1)) for v in fams.values())}')


def _selfcheck():
    R = {'a': ('c', '+', [[0, 10]]), 'b': ('c', '+', [[20, 30]]), 'ab': ('c', '+', [[5, 25]]), 'a2': ('c', '+', [[40, 50]])}
    nm = {'a': 'X1', 'b': 'Y', 'ab': 'X1-Y', 'a2': 'X2'}
    assert classify(['a', 'b'], R, nm, set())[1] == 'other' and classify(['a', 'b', 'ab'], R, nm, set())[1] == 'annot_rt'
    assert classify(['a', 'a2'], R, nm, set())[1] == 'named' and classify(['a', 'a2'], R, {'a': 'LOC1', 'a2': 'LOC22'}, set())[1] == 'loc' and classify(['a', 'ab'], R, nm, set())[:2] == ('shared-exon cluster', None)
    assert classify(['a', 'b'], R, nm, {frozenset(('X1', 'Y'))})[1] == 'truth' and classify(['a'], R, nm, set())[0] == 'single'


if __name__ == '__main__':
    _selfcheck(); main(*sys.argv[1:5])
