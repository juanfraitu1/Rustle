#!/usr/bin/env python3
"""Difference model: why a multi-exon intron chain is in one assembler's output and not the other's (read evidence only).
usage: diff_model.py primary.sam(cols 1-6) ref.gtf ours.gtf other.gtf other_name
Read chains use the assembler's default de-duplication key (start, end, chain). Every chain exclusive to one side gets the
FIRST class that applies (order below); `ref` = share of the class whose chain a reference transcript has exactly.
Other-only classes:   ours_gate   >= 2 reads carry it exactly (ours would admit at the floor; lost in gate / polish / motif)
                      single_shared  1 read carries it, every junction also carried by another read
                      single_novel   1 read carries it, a junction only that read carries
                      corrected   no read carries it; a read chain matches it within 25 bp at every junction (junction correction)
                      stitched    no read carries it, none within 25 bp (paths built from several reads: flow / extension)
Ours-only classes:    near_stronger  a junction lies within 25 bp (donor and acceptor) of one more reads carry (StringTie -E style)
                      two_reads   exactly 2 reads carry it (FLAIR's default support is 3)
                      subset      a contiguous sub-chain of another of OUR transcripts (subset isoforms other tools drop)
                      other"""
import sys, re, collections, bisect

W = 25

def read_chains(sam):
    seen = set(); chains = collections.Counter()
    for ln in open(sam):
        f = ln.split('\t'); pos = int(f[3]) - 1; start = pos; introns = []; aligned = False
        for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
            n = int(n)
            if op in 'M=X': aligned = True; pos += n
            elif op == 'D': pos += n
            elif op == 'N':
                if aligned: introns.append((pos, pos + n))
                pos += n
        key = (start, pos, tuple(introns))
        if key in seen: continue
        seen.add(key)
        if introns: chains[tuple(introns)] += 1
    return chains

def gtf_chains(p):
    tx = collections.defaultdict(list)
    for ln in open(p):
        f = ln.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            tx[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    out = set()
    for ex in tx.values():
        ex.sort()
        if len(ex) > 1: out.add(tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))
    return out

reads = read_chains(sys.argv[1]); REF = gtf_chains(sys.argv[2]); OURS = gtf_chains(sys.argv[3]); OTH = gtf_chains(sys.argv[4]); name = sys.argv[5]
jsup = collections.Counter()
for c, n in reads.items():
    for j in c: jsup[j] += n
js = sorted(jsup); jd = [j[0] for j in js]
def stronger_near(j):
    lo = bisect.bisect_left(jd, j[0] - W); hi = bisect.bisect_right(jd, j[0] + W)
    return any(k != j and abs(k[1] - j[1]) <= W and jsup[k] > jsup[j] for k in js[lo:hi])
# read chains bucketed by (length, first donor // W) for the 25-bp chain match
buck = collections.defaultdict(list)
for c in reads: buck[(len(c), c[0][0] // W)].append(c)
def corrected(c):
    for b in (c[0][0] // W - 1, c[0][0] // W, c[0][0] // W + 1):
        for r in buck.get((len(c), b), []):
            if all(abs(a[0] - x[0]) <= W and abs(a[1] - x[1]) <= W for a, x in zip(c, r)): return True
    return False
by_first = collections.defaultdict(list)
for d in OURS:
    for i, j in enumerate(d): by_first[j].append((d, i))
def subset(c):
    return any(len(d) > len(c) and d[i:i + len(c)] == c for d, i in by_first[c[0]])
def other_cls(c):
    n = reads.get(c, 0)
    if n >= 2: return 'ours_gate'
    if n == 1: return 'single_shared' if all(jsup[j] >= 2 for j in c) else 'single_novel'
    return 'corrected' if corrected(c) else 'stitched'
def ours_cls(c):
    if any(stronger_near(j) for j in c): return 'near_stronger'
    if reads.get(c, 0) == 2: return 'two_reads'
    if subset(c): return 'subset'
    return 'other'
def table(chs, f, title, order):
    cnt = collections.Counter(); tru = collections.Counter()
    for c in chs:
        k = f(c); cnt[k] += 1; tru[k] += c in REF
    tot = sum(cnt.values())
    print(f"{title}: {tot:,} chains, reference-exact {sum(tru.values()) / max(1, tot):.3f}")
    for k in order:
        print(f"  {k:14s} {cnt[k]:6,}  share {cnt[k] / max(1, tot):.3f}  reference-exact {tru[k] / max(1, cnt[k]):.3f}")
both = OURS & OTH
print(f"ours {len(OURS):,} chains, {name} {len(OTH):,}, shared {len(both):,} (reference-exact {sum(c in REF for c in both) / max(1, len(both)):.3f})")
table(OTH - OURS, other_cls, f"{name}-only", ['ours_gate', 'single_shared', 'single_novel', 'corrected', 'stitched'])
table(OURS - OTH, ours_cls, "ours-only", ['near_stronger', 'two_reads', 'subset', 'other'])
