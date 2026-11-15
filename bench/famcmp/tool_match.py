#!/usr/bin/env python3
"""Step 0 for the tool-matching goal (2026-11-11): how close is an assembly's multi-exon intron-chain set to a tool's?

  tool_match.py REFSEQ_GTF TOOL=GTF [TOOL=GTF ...] -- ARM=GTF [ARM=GTF ...]      (selftest: tool_match.py selftest)

Per arm x tool: shared chains, arm-only, tool-only, recall of the tool's set, share of the arm in the tool's set, Jaccard;
per arm: chains, RefSeq-correct, precision. Chain = (contig, strand, introns)."""
import sys, re, gzip, collections

def chains(lines):
    ex = collections.defaultdict(list)
    for l in lines:
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    out = set()
    for (c, s, _), e in ex.items():
        e.sort()
        if len(e) > 1: out.add((c, s, tuple((a[1], b[0]) for a, b in zip(e, e[1:]))))
    return out
def cmp(a, t): s = len(a & t); return s, len(a) - s, len(t) - s, s / max(len(t), 1), s / max(len(a), 1), s / max(len(a | t), 1)
def read(p): return chains((gzip.open if p.endswith('.gz') else open)(p, 'rt'))

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        g = lambda t, exons, s='+': [f'c\tx\texon\t{a + 1}\t{b}\t.\t{s}\t.\ttranscript_id "{t}";' for a, b in exons]
        a = chains(g('1', [(0, 10), (20, 30)]) + g('2', [(0, 10), (40, 50)]) + g('m', [(0, 99)]))
        t = chains(g('x', [(5, 10), (20, 25)]) + g('y', [(0, 10), (20, 30)], '-'))
        assert cmp(a, t) == (1, 1, 1, 0.5, 0.5, 1 / 3), cmp(a, t)
        print('ok'); sys.exit()
    i = sys.argv.index('--'); ref = read(sys.argv[1])
    tools = [(a.split('=')[0], read(a.split('=', 1)[1])) for a in sys.argv[2:i]]
    print('arm\tchains\tcorrect\tprecision\t' + '\t'.join(f'{n}:shared\t{n}:arm_only\t{n}:tool_only\t{n}:recall\t{n}:arm_in_tool\t{n}:jaccard' for n, _ in tools))
    for a in sys.argv[i + 1:]:
        n, p = a.split('=', 1); A = read(p); ok = len(A & ref)
        print(f'{n}\t{len(A)}\t{ok}\t{ok / max(len(A), 1):.3f}\t' + '\t'.join('%d\t%d\t%d\t%.3f\t%.3f\t%.3f' % cmp(A, t) for _, t in tools), flush=True)
