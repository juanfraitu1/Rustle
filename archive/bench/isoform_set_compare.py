#!/usr/bin/env python3
"""Compare two isoform sets on the SAME reads and the same family copies.

Written for `isoseq collapse` vs ours, but takes any GTF/GFF with `transcript`/`exon` rows and a
`transcript_id` attribute (StringTie, flair, isoseq collapse, ours). Reports what each set does with the
reads, not which isoforms are "right" — there is no transcript truth here.

usage: isoform_set_compare.py <A.gtf> <B.gtf> <bam> <copies.tsv> [--labels A,B]
"""
import sys, re, collections, subprocess, csv, statistics as st

a_p, b_p, bam, copies_p = sys.argv[1:5]
labels = (sys.argv[sys.argv.index('--labels') + 1].split(',') if '--labels' in sys.argv else ['A', 'B'])

cop = list(csv.DictReader(open(copies_p), delimiter='\t'))
fam = [(r['chrom'], int(r['start']), int(r['end'])) for r in cop]
def ov(a, b, c, d): return max(0, min(b, d) - max(a, c))

def load(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        m = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
        if not m: continue
        t = m.group(1)
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4])))
        elif f[2] == 'transcript': meta[t] = (f[0], int(f[3]) - 1, int(f[4]))
    out = {}
    for t, v in ex.items():
        v.sort()
        c, s, e = meta.get(t, (None, v[0][0], v[-1][1]))
        if c is None: continue
        if not any(c == fc and ov(s, e, fs, fe) > 0 for fc, fs, fe in fam): continue
        out[t] = dict(chain=tuple((x[1], y[0]) for x, y in zip(v, v[1:])), nex=len(v), c=c, s=s, e=e)
    return out

A, B = load(a_p), load(b_p)
def chains(S): return {v['chain'] for v in S.values() if v['chain']}
def juncs(S):
    j = set()
    for v in S.values(): j.update((v['c'],) + x for x in v['chain'])
    return j
ca, cb, ja, jb = chains(A), chains(B), juncs(A), juncs(B)

def introns(pos, cig):
    o = []; p = pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD': p += n
        elif op == 'N': o.append((p, p + n)); p += n
    return tuple(o)

seen = set(); tot = spl = ea = eb = 0; jr = collections.Counter()
by = collections.defaultdict(list)
for c, s, e in fam: by[c].append((s, e))
for c in by:
    iv = sorted(by[c]); m = []
    for s, e in iv:
        if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    for s, e in m:
        out = subprocess.run(f"samtools view -F 2308 {bam} {c}:{s+1}-{e}", shell=True, capture_output=True, text=True).stdout
        for line in out.split('\n'):
            f = line.split('\t')
            if len(f) < 6 or f[0] in seen: continue
            seen.add(f[0]); ch = introns(int(f[3]) - 1, f[5]); tot += 1
            if ch:
                spl += 1
                if ch in ca: ea += 1
                if ch in cb: eb += 1
                for x in ch: jr[(c,) + x] += 1
rj = {j for j, n in jr.items() if n >= 2}
def stray(S):
    return sum(1 for v in S.values() if v['e'] - v['s'] > 2 * max((ov(v['s'], v['e'], fs, fe) for fc, fs, fe in fam if fc == v['c']), default=1))
print(f"molecules in the family's copies: {tot:,} ({spl:,} spliced); read-supported junctions (>=2): {len(rj):,}\n")
print(f"{'':44s} {labels[0]:>12} {labels[1]:>12}")
rows = [
    ("transcripts overlapping a copy", len(A), len(B)),
    ("distinct spliced intron chains", len(ca), len(cb)),
    (f"  of those, shared between the two", len(ca & cb), len(ca & cb)),
    ("spliced molecules whose EXACT chain is present", ea, eb),
    ("  as a fraction of spliced molecules", round(ea / max(spl, 1), 3), round(eb / max(spl, 1), 3)),
    ("read-supported junctions recovered", len(ja & rj), len(jb & rj)),
    ("  as a fraction of read-supported junctions", round(len(ja & rj) / max(len(rj), 1), 3), round(len(jb & rj) / max(len(rj), 1), 3)),
    ("junctions asserted with <2 molecules", len(ja - rj), len(jb - rj)),
    ("transcripts >2x their best-overlapping copy", stray(A), stray(B)),
]
for name, x, y in rows: print(f"{name:44s} {str(x):>12} {str(y):>12}")
def has_copy(path):
    n = 0
    for l in open(path):
        if '\ttranscript\t' in l and ('assigned_copy' in l or 'copy_index' in l): n += 1
    return n
print(f"\n{'transcripts carrying a COPY attribute':44s} {has_copy(a_p):>12} {has_copy(b_p):>12}")
print("  (a set with 0 here cannot say which copy an isoform came from — that is the comparison's point)")
