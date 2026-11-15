#!/usr/bin/env python3
"""PREREG row 1721: early-end models of nested_sub.py --models against annotated polyA sites.

  polya_ends.py MODELS_GTF SITES_BED6 FASTA [SKIP_CONTIG,...]
Supported = a site on the strand within 25 bp of the model's 3' end; background = 200 bp upstream in transcript sense.
Signature = peaked end (id field), polyA signal upstream, no A-rich run after (alt_ends.flags)."""
import sys, re, bisect, collections, pysam
from alt_ends import flags
S = collections.defaultdict(list)
for l in open(sys.argv[2]):
    f = l.split('\t'); S[f[0], f[5].strip()].append(int(f[1]))
for v in S.values(): v.sort()
def near(c, s, p):
    v = S.get((c, s), ()); i = bisect.bisect_left(v, p - 25); return int(i < len(v) and v[i] <= p + 25)
fa = pysam.FastaFile(sys.argv[3]); skip = set(sys.argv[4].split(',')) if len(sys.argv) > 4 else set(); T = collections.Counter()
for l in open(sys.argv[1]):
    f = l.split('\t')
    if len(f) < 9 or f[2] != 'transcript' or f[0] in skip: continue
    cls, grp, peaked, reads, _ = re.search(r'transcript_id "([^"]+)"', f[8]).group(1).rsplit('_', 5)[1:]
    if cls != 'NEST3': continue
    plus = f[6] == '+'; p = int(f[4]) if plus else int(f[3]) - 1; sig, arich = flags(fa.fetch, f[0], plus, p)
    k = (grp, 'signature' if peaked == '1' and sig and not arich else 'no signature')
    T[k] += 1; T[k + ('s',)] += near(f[0], f[6], p); T[k + ('b',)] += near(f[0], f[6], p - 200 if plus else p + 200)
print('group\tset\tmodels\tsupported\tbackground')
for grp, name in (('EA', 'emitted annotated'), ('NA', 'not emitted annotated'), ('EU', 'emitted unannotated'), ('NU', 'not emitted unannotated')):
    for sg in ('signature', 'no signature', 'all'):
        n, s, b = (sum(T[(grp, x) + t] for x in (('signature', 'no signature') if sg == 'all' else (sg,))) for t in ((), ('s',), ('b',)))
        print(f'{name}\t{sg}\t{n}\t{s / max(1, n):.3f}\t{b / max(1, n):.3f}')
