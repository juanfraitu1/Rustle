#!/usr/bin/env python3
"""Row 1668 step 0: features of single-exon models a polish change would restore, and a held-out rule search.

  mono_features.py FASTA FUNNEL.tsv BASE.gtf ALLOFF.gtf OUT.tsv      (env as pairscore.py; --selftest)
"""
import sys, os, re, bisect, collections, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
COMP = str.maketrans('ACGT', 'TGCA')

def rich(w): return bool(w) and (w.count('A') >= 0.6 * len(w) or 'AAAAAA' in w)

def end_feats(seq, s, e, strand):
    """(pas, primed) at the 3' end of [s, e) on `strand`; windows read in transcript orientation."""
    if strand == '+': up, dn = seq[max(0, e - 40):e], seq[e:e + 20]
    else: up, dn = seq[s:s + 40].translate(COMP)[::-1], seq[max(0, s - 20):s].translate(COMP)[::-1]
    return ('AATAAA' in up or 'ATTAAA' in up), rich(dn)

def read_tx(gtf):
    tx = collections.defaultdict(list); reads = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = f[8].split('transcript_id "')[1].split('"')[0]
        if f[2] == 'exon': tx[t].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
        else: reads[t] = int(re.search(r'reads "(\d+)"', f[8]).group(1))   # the attribute is on the transcript line only
    return tx, reads

def is_train(c):
    n = c[3:]
    return n == 'X' or (n.isdigit() and int(n) % 2 == 1)

GRID = list(itertools.product((2, 3, 5, 10, 20, 50), ('any', 'not_exon', 'free'), (0, 1), (0, 1), (0, 10000, 5000)))
def keep(r, rule):
    mr, sh, pas, noprime, ml = rule
    return (r['reads'] >= mr and (sh == 'any' or (sh == 'not_exon' and r['shadow'] != 'exon') or r['shadow'] == 'free')
            and (not pas or r['pas']) and (not noprime or not r['primed']) and (not ml or r['len'] <= ml))
def score(rows, rule):
    k = [r for r in rows if keep(r, rule)]
    return len({r['gene'] for r in k if r['label'] == 'member'}), len(k), len({r['gene'] for r in k if r['label'] == 'annotated'})

def main(fa, funnel, base, alloff, out):
    import pysam, famtruth as F
    G = F.genes_from_argv([], 'human'); idx = F.exon_index(G); FA = pysam.FastaFile(fa)
    nomodel = {l.split('\t')[0] for l in open(funnel) if l.split('\t')[5] == 'no_model'}
    btx, _ = read_tx(base); atx, areads = read_tx(alloff)
    bmono = {v[0] for v in btx.values() if len(v) == 1}
    ex = collections.defaultdict(list); span = collections.defaultdict(list)
    for v in btx.values():
        if len(v) > 1:
            for c, st, a, b in v: ex[c].append((a, b))
            span[(v[0][0], v[0][1])].append((min(x[2] for x in v), max(x[3] for x in v)))
    def hit(lst, s, e, mx):
        i = bisect.bisect_left(lst, (s - mx,))
        while i < len(lst) and lst[i][0] < e:
            if lst[i][1] > s: return True
            i += 1
        return False
    for d in (ex, span):
        for k in d: d[k].sort()
    mx = {k: max(b - a for a, b in v) for d in (ex, span) for k, v in d.items()}
    rows = []; seqc = {}
    for t, v in atx.items():
        if len(v) != 1 or v[0] in bmono: continue
        c, st, s, e = v[0]
        if c not in seqc: seqc.clear(); seqc[c] = FA.fetch(c).upper()
        pas, primed = end_feats(seqc[c], s, e, st)
        sh = 'exon' if hit(ex[c], s, e, mx.get(c, 0)) else 'span' if hit(span[(c, st)], s, e, mx.get((c, st), 0)) else 'free'
        g = F.map_copy(idx, c, st, [(s, e)]); lab = 'none'
        if g:
            gl = sum(y - x for x, y in G[g]['exons']); o = sum(max(0, min(y, e) - max(x, s)) for x, y in G[g]['exons'])
            lab = 'member' if g in nomodel else 'annotated' if 2 * o >= gl else 'none'
        rows.append(dict(t=t, chrom=c, reads=areads[t], len=e - s, pas=pas, primed=primed, shadow=sh, label=lab, gene=g or ''))
    with open(out, 'w') as w:
        w.write('transcript\tchrom\treads\tlen\tpas\tprimed\tshadow\tlabel\tgene\n')
        for r in rows: w.write('\t'.join(str(r[k]) for k in ('t', 'chrom', 'reads', 'len', 'pas', 'primed', 'shadow', 'label', 'gene')) + '\n')
    tr = [r for r in rows if is_train(r['chrom'])]; te = [r for r in rows if not is_train(r['chrom'])]
    for nm, part in (('train', tr), ('test', te)):
        print(nm, 'models', len(part), 'members', len({r['gene'] for r in part if r['label'] == 'member'}), 'annotated genes', len({r['gene'] for r in part if r['label'] == 'annotated'}))
        for f in ('shadow', 'pas', 'primed'):
            print(' ', f, {lab: dict(collections.Counter(r[f] for r in part if r['label'] == lab)) for lab in ('member', 'annotated', 'none')})
        for lab in ('member', 'annotated', 'none'):
            x = sorted(r['reads'] for r in part if r['label'] == lab); y = sorted(r['len'] for r in part if r['label'] == lab)
            if x: print('  reads/len quartiles', lab, [x[len(x) * q // 4] for q in (1, 2, 3)], [y[len(y) * q // 4] for q in (1, 2, 3)])
    sc = [(rule, score(tr, rule)) for rule in GRID]
    ok = [(m, -n, rule) for rule, (m, n, a) in sc if m and n <= 50 * m]
    print('train rules within 50 models per member:', len(ok))
    for rule, (m, n, a) in sorted(sc, key=lambda x: (-x[1][0] / max(1, x[1][1]) ** 0.5))[:8]: print('  train', rule, 'members', m, 'models', n, 'annotated genes', a)
    if ok:
        m, n, rule = max(ok); tm, tn, ta = score(te, rule)
        print('CHOSEN', rule, 'train members', m, 'models', -n, '| test members', tm, 'models', tn, 'annotated genes', ta, 'GO' if tm >= 15 and tn <= 50 * tm else 'STOP')
    else: print('no train rule within 50 models per member: STOP')

def selftest():
    seq = 'C' * 60 + 'AATAAA' + 'C' * 14 + 'A' * 20 + 'C' * 20          # + strand end at 80: PAS upstream, A tract downstream
    assert end_feats(seq, 0, 80, '+') == (True, True) and end_feats(seq, 0, 60, '+') == (False, False)
    rc = seq.translate(COMP)[::-1]; n = len(seq)
    assert end_feats(rc, n - 80, n, '-') == (True, True)
    r = dict(reads=5, shadow='span', pas=True, primed=False, len=2000, label='member', gene='g')
    assert keep(r, (5, 'not_exon', 1, 1, 5000)) and not keep(r, (10, 'any', 0, 0, 0)) and not keep(r, (2, 'free', 0, 0, 0))
    assert score([r, dict(r, label='none', gene='')], (2, 'any', 0, 0, 0)) == (1, 2, 0)
    assert is_train('chr1') and is_train('chrX') and not is_train('chr2') and not is_train('chrY')
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(*sys.argv[1:6])
