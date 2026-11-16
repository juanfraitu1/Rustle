"""summary.py TAG : who outputs each RefSeq chain, by read evidence, and what rules S1 / S2 recover (TAG.t1.tsv, TAG.rules.keys)."""
import sys, collections
n = sys.argv[1]; add = collections.defaultdict(set)
for l in open(n + '.rules.keys'): r, c, ch = l.rstrip('\n').split('\t'); add[(c, ch)].add(r)
H = None; T = collections.Counter(); per = collections.Counter(); rec = collections.Counter()
tools = ['st', 'flair', 'iq', 'bambu', 'bambu_ndr1', 'cl2']
for l in open(n + '.t1.tsv'):
    f = l.rstrip('\n').split('\t')
    if H is None: H = f; continue
    r = dict(zip(H, f)); o = r['ours'] == '1'; t = any(r[x] == '1' for x in tools); tier = r['tier']
    T[(tier, 'ours' if o else '-', 'tool' if t else '-')] += 1; T[('all', 'ours' if o else '-', 'tool' if t else '-')] += 1
    for x in ['ours', 'final'] + tools: per[x] += r[x] == '1'
    if not o and not t and tier in ('full', 'linked'):
        rec[(tier, 'all')] += 1; rs = add.get((r['chrom'], r['chain']), ())
        for x in rs: rec[(tier, x)] += 1
        rec[(tier, 'S1 or S2')] += bool(rs); rec[(tier, 'floor1')] += r['floor1'] == '1'
print(n, 'RefSeq chains output per arm:', ', '.join('%s %d' % (x, per[x]) for x in ['ours', 'final'] + tools))
print('tier\tours only\ttools only\tboth\tnobody')
for tier in ('full', 'linked', 'unlinked', 'nodata', 'all'): print(tier, T[(tier, 'ours', '-')], T[(tier, '-', 'tool')], T[(tier, 'ours', 'tool')], T[(tier, '-', '-')], sep='\t')
for tier in ('full', 'linked'): print('nobody, %s: %d; in the 1-read-floor run %d; S1 %d, S2 %d, either %d' % (tier, rec[(tier, 'all')], rec[(tier, 'floor1')], rec[(tier, 'S1')], rec[(tier, 'S2')], rec[(tier, 'S1 or S2')]))
