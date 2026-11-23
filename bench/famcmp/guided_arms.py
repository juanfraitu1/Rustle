#!/usr/bin/env python3
"""guided_arms.py REF GUIDE_PREFIX PREFIX [name=tool.gtf ...]  (register 1934, docs/PREREG_guided_final_oneread_2026-11-23.md)

Two candidate rules for the guided mode, scored against today's PREFIX.guided.gtf (rule G1d, register 1622) on a guide with hidden
transcripts (bench/famcmp/guide_proto.py split, GUIDE_HIDE=tx). Multi-exon intron chains throughout.
  arm A: G1d's verdicts applied to PREFIX.final.gtf (chain score filter and layers first) instead of PREFIX.gtf; floors unchanged
         (10 reads with a new intron beside guide introns, 3 reads otherwise); guide chains with a full-chain read as today.
  arm B: arm A + every chain of the floor-1 assembly (PREFIX.oneread.d/f1.gtf) that PREFIX.gtf lacks, whose introns are ALL guide
         introns, and that is neither a contiguous run nor an intron-retaining version of a guide chain.
Prints one block per arm: chains, guide found, hidden found, precision outside the guide; new and lost chains against today
(counted separately, never the net); the bars of the pre-registration. Extra name=gtf arguments are scored as reference rows."""
import sys, re, gzip, collections
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
NEW_INTRON_READS, NOVEL_READS = 10, 3
def load(path):  # (chrom, strand, chain) -> (max reads, a transcript id)
    ex = collections.defaultdict(list); rd = {}
    for l in op(path):
        f = l.split('\t')
        if len(f) < 9: continue
        if f[2] not in ('exon', 'transcript'): continue   # gene lines carry no transcript_id (IsoQuant)
        t = re.search(r'transcript_id "([^"]*)"', f[8]).group(1)
        if f[2] == 'exon': ex[(f[0], f[6], t)].append((int(f[3]), int(f[4])))
        elif f[2] == 'transcript':
            m = re.search(r'\breads "([^"]+)"', f[8]); rd[(f[0], f[6], t)] = float(m.group(1)) if m else 0.0
    out = {}
    for (c, s, t), e in ex.items():
        if len(e) < 2: continue
        e.sort(); ch = (c, s, tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)))
        out[ch] = max(out.get(ch, (0.0, t)), (rd.get((c, s, t), 0.0), t))
    return out
def within(ch, g, contiguous):
    if contiguous: return any(g[i:i + len(ch)] == ch for i in range(len(g) - len(ch) + 1))
    it = iter(g); return all(any(j == i for j in it) for i in ch)
def main():
    ref_p, gpre, pre = sys.argv[1:4]; extra = [a.split('=', 1) for a in sys.argv[4:]]
    ref = set(load(ref_p)); gui = set(load(gpre + '.guide.gtf')); hid = ref - gui
    byi = collections.defaultdict(list)
    for c, s, ch in gui:
        for i in ch: byi[(c, s, i)].append(ch)
    def cls(k):
        c, s, ch = k
        if k in gui: return 'guide'
        n = sum((c, s, i) in byi for i in ch)
        if n == len(ch):
            cand = byi[(c, s, ch[0])]
            if any(within(ch, g, True) for g in cand): return 'fragment'
            if any(within(ch, g, False) for g in cand): return 'retained'
            return 'combination'
        return 'new_intron' if n else 'outside'
    raw = load(pre + '.gtf'); fin = load(pre + '.final.gtf'); today = set(load(pre + '.guided.gtf')); f1 = load(pre + '.oneread.d/f1.gtf')
    resc = today & gui
    keep = lambda k, r: (c := cls(k)) not in ('guide', 'fragment', 'retained') and r >= (NEW_INTRON_READS if c == 'new_intron' else NOVEL_READS)
    A = {k for k, v in fin.items() if keep(k, v[0])} | resc
    pool = {k for k in f1 if k not in raw}; add = {k for k in pool if cls(k) == 'combination'} - A
    B = A | add
    prec = lambda S: len(S & hid) / max(1, len(S - gui))
    row = lambda n, S: print('%s\t%d\t%d\t%d\t%.4f' % (n, len(S), len(S & gui), len(S & hid), prec(S)))
    print('# guide chains %d, hidden chains %d; floor-1 chains absent from PREFIX.gtf %d (hidden %d)' % (len(gui), len(hid), len(pool), len(pool & hid)))
    print('arm\tchains\tguide found\thidden found\tprecision outside the guide')
    row('de novo final.gtf', set(fin)); row('guided today', today); row('arm A', A); row('arm B', B)
    for n, p in extra: row(n, set(load(p)))
    lostA, newA = today - A, A - today
    print('A vs today\tlost %d (hidden %d, share %.4f)\tnew %d (hidden %d)' % (len(lostA), len(lostA & hid), len(lostA & hid) / max(1, len(lostA)), len(newA), len(newA & hid)))
    a1 = len(lostA & hid) <= 0.02 * len(today & hid); a2 = len(lostA & hid) / max(1, len(lostA)) <= 0.5 * prec(today)
    print('bar A1 hidden lost <= 2%% of today\'s %d: %s | bar A2 hidden share of the lost chains <= half of today\'s precision %.4f: %s | arm A %s' % (len(today & hid), a1, prec(today), a2, 'PASS' if a1 and a2 else 'FAIL'))
    print('B vs A\tnew %d (hidden %d, share %.4f)\tlost 0' % (len(add), len(add & hid), len(add & hid) / max(1, len(add))))
    b1 = len(add & hid) / max(1, len(add)) >= prec(A); b2 = len(add & hid) >= 0.02 * len(A & hid)
    print('bar B1 hidden share of the added chains >= arm A\'s precision %.4f: %s | bar B2 hidden added >= 2%% of arm A\'s %d: %s | arm B %s' % (prec(A), b1, len(A & hid), b2, 'PASS' if b1 and b2 else 'FAIL'))
    cnt = collections.Counter(); h = collections.Counter()
    for k in pool: c = cls(k); cnt[c] += 1; h[c] += k in hid
    print('floor-1 pool by class (chains / hidden): ' + '; '.join('%s %d / %d' % (c, cnt[c], h[c]) for c in sorted(cnt)))
if __name__ == '__main__': main()
