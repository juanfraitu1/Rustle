#!/usr/bin/env python3
"""Where the ranked list puts annotated isoforms the main assembly lacks (register 1966, step 1).

    ranked_thin.py ANNOTATION.gff RANKED.gtf MAIN.gtf
    ranked_thin.py --tags ANNOTATION.gff RANKED.gtf MAIN.gtf TAGGED.gtf ...   (the nascent_of tag among each tier's chains)

RANKED.gtf = PREFIX.ranked.gtf (attributes list_tier, list_rank). Per tier, in tenths of the tier's own order: spliced
chains, annotated among them, and the cumulative share of the annotated isoforms of that tier reached so far. Tier 1 is
judged whole. `MAIN.gtf` (PREFIX.gtf or the mode's GTF) only says how many annotated isoforms the list adds to it."""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thin_isoforms as T

def tiers(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9: continue
        if f[2] == 'transcript':
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); a = re.search(r'list_tier "(\d+)"', f[8]); b = re.search(r'list_rank "(\d+)"', f[8])
            meta[t] = (f[0], f[6], int(a.group(1)) if a else 0, int(b.group(1)) if b else 0)
        elif f[2] == 'exon': ex[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    out = collections.defaultdict(list)
    for t, e in ex.items():
        e.sort()
        if len(e) > 1 and t in meta: c, s, tier, rank = meta[t]; out[tier].append((rank, (c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))))
    return {k: [c for _, c in sorted(v)] for k, v in out.items()}

def tagged(path, attr='nascent_of'):
    """intron chains of the transcripts of a GTF that carry the attribute"""
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.split('\t')
        if len(f) < 9: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if f[2] == 'transcript' and attr + ' "' in f[8]: meta[t] = (f[0], f[6])
        elif f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4])))
    out = set()
    for t, (c, s) in meta.items():
        e = sorted(ex[t])
        if len(e) > 1: out.add((c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))))
    return out

def tags(gff, ranked, main_gtf, *gtfs):
    ann = T.chains_gff(gff); have = T.chains_gtf(main_gtf) & ann; R = tiers(ranked); tg = set()
    for g in gtfs: tg |= tagged(g)
    print(f'chains with the nascent tag in the tagged assemblies: {len(tg)}')
    for k in sorted(R):
        if k == 1: continue
        v = R[k]; new = [c for c in v if c in ann and c not in have]; un = [c for c in v if c not in ann]
        print(f'## tier {k}: annotated and new {len(new)}, tagged {sum(c in tg for c in new)} ({sum(c in tg for c in new) / max(1, len(new)):.3f}); unannotated {len(un)}, tagged {sum(c in tg for c in un)} ({sum(c in tg for c in un) / max(1, len(un)):.3f})')
        n = len(v) // 10; re_ = [c for c in v if c not in tg] + [c for c in v if c in tg]
        for lab, o in (('as shipped', v), ('tagged chains last', re_)): print(f'   first tenth, {lab}: {sum(c in ann for c in o[:n])} annotated of {n} ({sum(c in ann for c in o[:n]) / max(1, n):.3f})')

def main(gff, ranked, main_gtf):
    ann = T.chains_gff(gff); have = T.chains_gtf(main_gtf) & ann; R = tiers(ranked)
    print(f'annotated isoforms in the main assembly: {len(have)}')
    for k in sorted(R):
        v = R[k]; hit = [c in ann for c in v]; new = [c in ann and c not in have for c in v]; n = len(v); tot = sum(new)
        print(f'## tier {k}: {n} spliced chains, {sum(hit)} annotated ({sum(hit) / max(1, n):.3f}); annotated ones the main assembly lacks: {tot}')
        if k == 1 or n < 10: continue
        print('tenth of the tier\tchains\tannotated\tshare annotated\tannotated and new, cumulative\tshare of the tier\'s new annotated')
        cum = 0
        for i in range(10):
            a, b = n * i // 10, n * (i + 1) // 10; cum += sum(new[a:b]); print(f'{i + 1}\t{b - a}\t{sum(hit[a:b])}\t{sum(hit[a:b]) / max(1, b - a):.3f}\t{cum}\t{cum / max(1, tot):.3f}')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        import tempfile
        d = tempfile.mkdtemp(); p = os.path.join(d, 'r.gtf'); L = []
        for t, tier, rank, s2 in (('a', 2, 2, 21), ('b', 2, 1, 41)):
            L += [f'c\tx\ttranscript\t1\t60\t.\t+\t.\ttranscript_id "{t}"; list_tier "{tier}"; list_rank "{rank}";', f'c\tx\texon\t1\t10\t.\t+\t.\ttranscript_id "{t}";', f'c\tx\texon\t{s2}\t60\t.\t+\t.\ttranscript_id "{t}";']
        L[0] += ' nascent_of "b";'; open(p, 'w').write('\n'.join(L) + '\n'); assert tiers(p) == {2: [('c', '+', ((10, 40),)), ('c', '+', ((10, 20),))]}; assert tagged(p) == {('c', '+', ((10, 20),))}; print('ranked_thin selftest ok')
    elif sys.argv[1] == '--tags': tags(*sys.argv[2:])
    else: main(*sys.argv[1:])
