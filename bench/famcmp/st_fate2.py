#!/usr/bin/env python3
"""r1636b: read-level fate, inside StringTie, of chains another assembly has and StringTie lacks (and the reverse).

  st_fate2.py RUN1_JSONL CONTIG OURS_GTF STRINGTIE_GTF REFSEQ_GTF        (log from -p 1 with read_in, read_to_trf, graph_node_new, pred_*)

Ours-only chain: its reads = read_in records whose exons give exactly the chain. Each read's transfrag chain is rebuilt from
its graph nodes (an intron between two consecutive nodes that do not touch). Read outcome: equal / end_lost (the transfrag
chain is a contiguous run of the read's) / junction_moved (same count, every junction within 25 bp) / other / no_transfrag.
The chain takes its reads' most common outcome; when equal, the predictions with that chain say what happened next."""
import sys, json, collections
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read
from tool_match_drop import classify
log, contig = sys.argv[1:3]; ours, st, ref = (read(p) for p in sys.argv[3:6])
def chain_of(exons):
    e = [tuple(map(int, x.split('-'))) for x in exons.split(',') if '-' in x]
    return tuple((a[1], b[0] - 1) for a, b in zip(e, e[1:]))
nodes = {}; raw = collections.defaultdict(list); trf = {}; pred = {}; fate = collections.defaultdict(list); out = set(); bad = 0
for l in open(log):
    try: e = json.loads(l, strict=False)
    except ValueError: bad += 1; continue
    p = e['payload']; b = e.get('bundle'); s = e['step']
    if s == 'graph_node_new': nodes[(b, p['s'], p['g'], p['node'])] = (e['start'], e['end'])
    elif s == 'read_in' and p.get('action') in ('new_entry', 'collapsed_into_existing') and p.get('njunc', 0) > 0: raw[chain_of(p['exons'])].append((b, p['idx']))
    elif s == 'read_to_trf' and len(p['nodes']) > len(trf.get((b, p['read']), (0, 0, ''))[2]): trf[(b, p['read'])] = (p['s'], p['g'], p['nodes'])  # the read's longest transfrag
    elif s == 'pred_create': pred[(b, p['pid'])] = chain_of(p.get('exons', ''))
    elif s == 'pred_kill': fate[(b, p['pid'])].append('killed:' + p['reason'])
    elif s == 'pred_flag_off': fate[(b, p['pid'])].append('flag_off:' + p.get('func', '?'))
    elif s == 'pred_delete': fate[(b, p['pid'])].append('deleted:' + p.get('func', '?'))
    elif s == 'pred_output': out.add((b, p['pid']))
# a node is logged when created; a later split shortens it, so its end is the base before the next node of its graph
_g = collections.defaultdict(list)
for (b, s_, g, n), (a, e) in nodes.items(): _g[(b, s_, g)].append((a, e, n))
for k, L in _g.items():
    L.sort()
    for (a, e, n), nx in zip(L, L[1:] + [None]): nodes[k + (n,)] = (a, min(e, nx[0] - 1) if nx and nx[0] <= e else e)
def trf_chain(b, s, g, ns):
    c = [nodes.get((b, s, g, int(n))) for n in ns.split(',') if n]
    if None in c: return None
    return tuple((x[1], y[0] - 1) for x, y in zip(c, c[1:]) if y[0] > x[1] + 1)
run = lambda a, b: any(b[k:k + len(a)] == a for k in range(len(b) - len(a) + 1))
def read_outcome(ch, key):
    if key not in trf: return 'no_transfrag'
    t = trf_chain(key[0], *trf[key])
    if t is None: return 'no_transfrag'
    if t == ch: return 'equal'
    if len(t) < len(ch) and run(t, ch): return 'end_lost' if t else 'all_introns_lost'
    if len(t) == len(ch) and all(abs(x[0] - y[0]) <= 25 and abs(x[1] - y[1]) <= 25 for x, y in zip(t, ch)): return 'junction_moved'
    return 'other'
pred_by_chain = collections.defaultdict(list)
for k, c in pred.items(): pred_by_chain[c].append(k)
trf_reads = collections.Counter()
for k, (s, g, ns) in trf.items():
    t = trf_chain(k[0], s, g, ns)
    if t: trf_reads[t] += 1
idx = collections.defaultdict(lambda: (collections.defaultdict(list), collections.defaultdict(set)))
for c, s, ch in st:
    bf, bi = idx[(c, s)]; bf[(len(ch), ch[0][0])].append(ch)
    for i in ch: bi[i].add(ch)
tab = collections.defaultdict(lambda: [0, 0]); cross = collections.Counter()
for c, s, ch in ours - st:
    if c != contig: continue
    rs = raw.get(ch, [])
    if not rs: f = 'no_read_with_chain_in_log'
    else:
        f = collections.Counter(read_outcome(ch, k) for k in rs).most_common(1)[0][0]
        if f == 'equal':
            ks = pred_by_chain.get(ch, [])
            if not ks: f = 'equal:never_predicted'
            else:
                r = collections.Counter(x.split(':')[0] + ':' + x.split(':')[1] for k in ks for x in (fate.get(k) or ['no:event'])[:1])
                f = 'equal:' + r.most_common(1)[0][0]
    tab[f][0] += 1; tab[f][1] += (c, s, ch) in ref; cross[(f, classify(ch, None, *idx[(c, s)]))] += 1
n = sum(v[0] for v in tab.values()); print(f'unparsable lines {bad}; ours-only chains on {contig}: {n}\nfate\tchains\tshare\tcorrect_share\tnear\tsub\tsuper\tother_shared\tnovel_locus')
for f, (k, ok) in sorted(tab.items(), key=lambda x: -x[1][0]):
    print(f'{f}\t{k}\t{k / n:.3f}\t{ok / k:.3f}\t' + '\t'.join(str(cross[(f, x)]) for x in ('near', 'sub', 'super', 'other_shared', 'novel_locus')))
# the add side: StringTie's own exclusive chains, by how many reads' transfrags are exactly that chain
oc = {(c, s, ch) for c, s, ch in st - ours if c == contig}; t2 = collections.defaultdict(lambda: [0, 0])
for k in oc:
    r = trf_reads.get(k[2], 0); b = '0 (built by flow / extension)' if r == 0 else '1' if r == 1 else '2+'
    t2[b][0] += 1; t2[b][1] += k in ref
print(f'\nStringTie-only chains on {contig}: {len(oc)}\nreads whose transfrag is the chain\tchains\tcorrect_share')
for b, (k, ok) in sorted(t2.items()): print(f'{b}\t{k}\t{ok / k:.3f}')
# r1639 oracle (STEP0_tool_match §8): chains of the reads' transfrags with >= k reads, against StringTie's output on the contig
stc = {ch for c, s, ch in st if c == contig}; oursc = {ch for c, s, ch in ours if c == contig}
print('\noracle set\tchains\tshared\treached\tinside\tjaccard')
for name, S in [('default', oursc)] + [(f'transfrag chains >= {k} reads', {t for t, n in trf_reads.items() if n >= k}) for k in (1, 2, 3)] + [('default + transfrag >= 2', oursc | {t for t, n in trf_reads.items() if n >= 2})]:
    i = len(S & stc); print(f'{name}\t{len(S)}\t{i}\t{i / len(stc):.3f}\t{i / max(len(S), 1):.3f}\t{i / len(S | stc):.3f}')
