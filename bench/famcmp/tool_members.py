#!/usr/bin/env python3
"""Row 1670 (descriptive): which expressed truth family members does each assembler model?

  tool_members.py [--species human|gorilla] FUNNEL.tsv NAME=GTF ...      (env as pairscore.py; --selftest)
Members = the findable truth members of a member_funnel.py table (expressed by >= 2 primary molecules, a truth partner
expressed too). Per tool, a member is
  modelled  some transcript has it as its best gene (most shared exon bases on the transcript's strand; an unstranded
            transcript is tried on both strands),
  clean     one such transcript lies >= 80% inside the member's exons and covers >= 50% of them (not a blob or a fusion).
Prints counts per tool, the overlap between tools (clean), and writes FUNNEL.tools.tsv.
"""
import sys, os, re, gzip, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))

def clean(tex, gex):
    """transcript exons vs gene exons (sorted lists of [s, e))"""
    o = sum(max(0, min(y, d) - max(x, c)) for x, y in gex for c, d in tex)
    return o >= 0.8 * sum(d - c for c, d in tex) and o >= 0.5 * sum(y - x for x, y in gex)

def transcripts(path):
    tx = collections.defaultdict(list); meta = {}
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    return [(meta[t][0], meta[t][1], sorted(e)) for t, e in tx.items()]

def main(argv):
    import famtruth as F
    species = argv[argv.index('--species') + 1] if '--species' in argv else 'human'
    args = [a for i, a in enumerate(argv) if a != '--species' and (i == 0 or argv[i - 1] != '--species')]
    funnel, arms = args[0], [a.split('=', 1) for a in args[1:]]
    G = F.genes_from_argv(argv, species); idx = F.exon_index(G)
    rows = [l.rstrip('\n').split('\t') for l in open(funnel)][1:]; mem = {r[0] for r in rows}
    res = {}
    for name, gtf in arms:
        mod, cl = set(), set(); n = 0
        for c, st, ex in transcripts(gtf):
            n += 1; best = None
            for s in ((st,) if st in '+-' else ('+', '-')):
                g = F.map_copy(idx, c, s, ex)
                if g and g in mem:
                    o = sum(max(0, min(y, d) - max(x, a)) for x, y in G[g]['exons'] for a, d in ex)
                    if not best or o > best[0]: best = (o, g)
            if best:
                mod.add(best[1])
                if clean(ex, G[best[1]]['exons']): cl.add(best[1])
        res[name] = (mod, cl); print(f'{name}\ttranscripts {n}\tmembers {len(mem)}\tmodelled {len(mod)} ({len(mod) / len(mem):.3f})\tclean {len(cl)} ({len(cl) / len(mem):.3f})')
    names = [a[0] for a in arms]; pat = collections.Counter(tuple(n for n in names if g in res[n][1]) for g in mem)
    for k, v in pat.most_common(): print('clean in', '+'.join(k) or 'none', v)
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        for n in names: by[r[5]][n] += r[0] in res[n][1]
        by[r[5]]['n'] += 1
    for st, c in by.items(): print('our stage', st, 'n', c['n'], {n: c[n] for n in names})
    with open(funnel.rsplit('.tsv', 1)[0] + '.tools.tsv', 'w') as w:
        w.write('gene\tname\tprimary\tour_stage\t' + '\t'.join(f'{n}_modelled\t{n}_clean' for n in names) + '\n')
        for r in rows: w.write(f'{r[0]}\t{r[1]}\t{r[3]}\t{r[5]}\t' + '\t'.join(f'{int(r[0] in res[n][0])}\t{int(r[0] in res[n][1])}' for n in names) + '\n')

def selftest():
    g = [(0, 100), (200, 300)]
    assert clean([(0, 100)], g) and not clean([(0, 100), (200, 300), (1000, 3000)], g) and not clean([(0, 40)], g)
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(sys.argv[1:])
