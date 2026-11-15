#!/usr/bin/env python3
"""Truth audit of false family pairs (docs/PREREG_family_join_2026-11-13.md): do the two annotated genes share exon sequence?

  falseaudit.py [--species human|hg002] --fasta GENOME --hits N FALSE_TSV...     (FALSE_TSV = PAIRSCORE_FALSEDUMP file)
  falseaudit.py --selftest
Classes, first that applies: related (pairscore.related) / shared (>= 0.5 of the shorter gene covered, pooled identity
>= 0.85) / partial (one record >= 300 bp at >= 0.70) / none. env: as pairscore.py; needs minimap2 and samtools on PATH.
"""
import sys, os, subprocess, tempfile, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def classify(recs, len_a, len_b):
    """recs: [(qs, qe, nmatch, blocklen)] of one gene pair (query side) -> shared / partial / none."""
    if not recs:
        return 'none'
    iv = sorted((s, e) for s, e, _, _ in recs); cov = 0; end = 0
    for s, e in iv:
        cov += max(0, e - max(s, end)); end = max(end, e)
    ident = sum(r[2] for r in recs) / max(1, sum(r[3] for r in recs))
    if cov >= 0.5 * min(len_a, len_b) and ident >= 0.85:
        return 'shared'
    return 'partial' if any(bl >= 300 and nm / bl >= 0.70 for _, _, nm, bl in recs) else 'none'


def run(a):
    import famtruth as F, importlib.util
    opt = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    files = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or a[i - 1] not in ('--species', '--fasta', '--hits'))]
    G = F.genes_from_argv(a, opt('--species', 'human')); hits = [int(x) for x in opt('--hits').split(',')]
    pairs = {f: [tuple(l.split('\t')[:2]) for l in open(f)] for f in files}
    rel = {tuple(l.split('\t')[:2]) for f in files for l in open(f) if l.rstrip('\n').split('\t')[4] == '1'}   # pairscore's related flag
    genes = sorted({g for v in pairs.values() for p in v for g in p})
    with tempfile.TemporaryDirectory() as d:
        reg = os.path.join(d, 'r'); fa = os.path.join(d, 'g.fa'); L = {}
        with open(fa, 'w') as w:
            for g in genes:
                ex = sorted(G[g]['exons'])
                out = subprocess.run(['samtools', 'faidx', opt('--fasta')] + [f"{G[g]['chrom']}:{s + 1}-{e}" for s, e in ex],
                                     capture_output=True, text=True, check=True).stdout
                seq = ''.join(l for l in out.split('\n') if not l.startswith('>'))
                L[g] = len(seq); w.write(f'>{g}\n{seq}\n')
        paf = subprocess.run(['minimap2', '-x', 'asm20', '-c', '-X', '-N', '500', '-p', '0.01', '--secondary=yes', fa, fa],
                             capture_output=True, text=True, check=True).stdout
    recs = collections.defaultdict(list)
    for l in paf.split('\n'):
        f = l.split('\t')
        if len(f) > 10 and f[0] != f[5]:
            recs[(f[0], f[5])].append((int(f[2]), int(f[3]), int(f[9]), int(f[10])))
    print('file\thits\tfalse\trelated\tshared\tpartial\tnone\tprecision\tprec_wo_related_shared\tprec_wo_partial_too')
    for (f, v), h in zip(pairs.items(), hits):
        c = collections.Counter()
        for x, y in v:
            k = 'related' if ((x, y) in rel or (y, x) in rel) else None
            if k is None:
                ka, kb = classify(recs.get((x, y), []), L[x], L[y]), classify(recs.get((y, x), []), L[x], L[y])
                k = min(ka, kb, key=('shared', 'partial', 'none').index)
            c[k] += 1
            if '--list' in a and k in ('partial', 'none'): print(f'#  {k}\t{G[x]["name"]}\t{G[y]["name"]}')
        n = len(v); p = lambda fl: f'{h / max(1, h + fl):.3f}'
        print(f'{os.path.basename(f)}\t{h}\t{n}\t{c["related"]}\t{c["shared"]}\t{c["partial"]}\t{c["none"]}\t{p(n)}\t{p(c["partial"] + c["none"])}\t{p(c["none"])}')


def selftest():
    assert classify([], 1000, 1000) == 'none'
    assert classify([(0, 600, 570, 600)], 1000, 2000) == 'shared'          # 0.6 of the shorter at 0.95
    assert classify([(0, 400, 380, 400)], 1000, 2000) == 'partial'         # 0.4 of the shorter
    assert classify([(0, 600, 450, 600)], 1000, 2000) == 'partial'         # covered, identity 0.75
    assert classify([(0, 200, 190, 200)], 1000, 2000) == 'none'            # under 300 bp
    assert classify([(0, 300, 290, 300), (200, 600, 390, 400)], 1000, 2000) == 'shared'   # overlap counted once: 600
    print('selftest ok')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else run(sys.argv[1:])
