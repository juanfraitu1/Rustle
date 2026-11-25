#!/usr/bin/env python3
"""Seed with uniquely placed primaries first, then bring in the tied reads (register 1963).

    two_step.py layers ASSEMBLY.gtf PROV.reads.tsv.gz OUT_PREFIX     arm A: OUT.unique.gtf / OUT.tied.gtf
    two_step.py unique IN.bam AS_TABLE OUT.bam                        primaries of the uniquely placed reads (step 1's input)
    two_step.py cond   IN.bam AS_TABLE STEP1.gtf OUT.bam              arm B: today's seeds minus tied alignments off step 1

A read is TIED when a rival alignment lies within 10 ms of its best (AS_TABLE column 7), else uniquely placed.
layers: a transcript is `unique` when >= 2 uniquely placed primaries carry its chain (provenance relation exact,
exact_after_replacement or unspliced, rival_within_10ms = no).
cond: seeds = primaries + secondaries at AS >= 0.98 x best AS (the `good` pool). For a tied read, among its seeds within
10 ms of its best: when at least one lies on a step-1 transcript (>= 50 exon bases), those that do not are dropped, the
primary included; a kept one is then written as the primary. When none does, the group stays whole. Only removes."""
import sys, os, re, gzip, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
BAND = 10

def layers(gtf, prov, out):
    n = collections.Counter()
    with gzip.open(prov, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); rel, tx, riv = h.index('relation'), h.index('transcript'), h.index('rival_within_10ms')
        for l in f:
            c = l.rstrip('\n').split('\t')
            if c[rel] in ('exact', 'exact_after_replacement', 'unspliced') and c[riv] == 'no': n[c[tx]] += 1
    u = open(out + '.unique.gtf', 'w'); t = open(out + '.tied.gtf', 'w'); k = collections.Counter()
    for l in open(gtf):
        m = re.search(r'transcript_id "([^"]+)"', l)
        if not m: continue
        w = u if n[m.group(1)] >= 2 else t; w.write(l)
        if '\ttranscript\t' in l: k['unique' if w is u else 'tied'] += 1
    print(f'two_step layers: {k["unique"]} unique, {k["tied"]} tied')

def band_records(bam, best):                                    # read -> [(ms, is_secondary, contig, start, blocks)] of its seeds in the band
    import pysam
    B = pysam.AlignmentFile(bam); R = collections.defaultdict(list); rivals = collections.Counter()
    for r in B:
        if r.is_unmapped or r.is_supplementary or not r.has_tag('ms'): continue
        b = best.get(r.query_name)
        if b is None or b[1] is None or r.get_tag('ms') < b[1] - BAND: continue
        rivals[r.query_name] += 1                               # any alignment in the band makes the read tied
        if r.is_secondary and b[0] is not None and r.has_tag('AS') and r.get_tag('AS') < 0.98 * b[0]: continue
        R[r.query_name].append((r.get_tag('ms'), r.is_secondary, r.reference_name, r.reference_start, r.get_blocks()))
    return R, {q for q, v in rivals.items() if v >= 2}

def unique(bam, tab, out):
    import pysam, pool_bam as PB
    best = PB.load_best(tab, bam); _, tied = band_records(bam, best); B = pysam.AlignmentFile(bam); n = 0
    with pysam.AlignmentFile(out, 'wb', template=B) as w:
        for r in B:
            if not r.is_unmapped and not r.is_secondary and r.query_name not in tied: w.write(r); n += 1
    pysam.index(out); print(f'two_step unique: {n} records of uniquely placed reads, {len(tied)} tied reads left out')

def cond(bam, tab, gtf, out):
    import pysam, pool_bam as PB
    ex = collections.defaultdict(list); BIN = 20000
    for l in open(gtf):
        f = l.split('\t', 6)
        if len(f) > 5 and f[2] == 'exon':
            s, e = int(f[3]) - 1, int(f[4])
            for k in range(s // BIN, e // BIN + 1): ex[(f[0], k)].append((s, e))
    def on_step1(c, blocks):
        n = 0; seen = set()
        for s, e in blocks:
            for k in range(s // BIN, e // BIN + 1):
                for a, b in ex.get((c, k), ()):
                    if a < e and s < b and (a, b, s) not in seen: seen.add((a, b, s)); n += min(b, e) - max(a, s)
        return n >= 50
    best = PB.load_best(tab, bam); R, tied = band_records(bam, best); drop = set(); promote = set(); n = collections.Counter()
    for q in tied:
        g = R[q]; ok = [on_step1(x[2], x[4]) for x in g]
        if not any(ok): n['tied reads, no place on step 1: group kept whole'] += 1; continue
        n['tied reads with a place on step 1'] += 1
        for x, o in zip(g, ok):
            if not o: drop.add((q, x[2], x[3], x[1])); n['alignments dropped: ' + ('secondary' if x[1] else 'PRIMARY')] += 1
        if not any(o and not x[1] for x, o in zip(g, ok)):      # the primary went: the best kept one takes its place
            x = max((x for x, o in zip(g, ok) if o), key=lambda x: x[0]); promote.add((q, x[2], x[3]))
    B = pysam.AlignmentFile(bam); lost = {q for q, c, s, sec in drop if not sec}
    with pysam.AlignmentFile(out + '.unsorted.bam', 'wb', template=B) as w:
        for r in B:
            if r.is_unmapped: continue
            q = r.query_name; b = best.get(q)
            if r.is_supplementary:
                if q not in lost: w.write(r)
                continue
            if r.is_secondary and b is not None and b[0] is not None and r.has_tag('AS') and r.get_tag('AS') < 0.98 * b[0]: n['secondaries outside the pool'] += 1; continue
            if (q, r.reference_name, r.reference_start, r.is_secondary) in drop: continue
            if r.is_secondary and (q, r.reference_name, r.reference_start) in promote:
                promote.discard((q, r.reference_name, r.reference_start)); r.is_secondary = False; r.mapping_quality = 0; n['secondaries written as the primary'] += 1
            w.write(r); n['records written'] += 1
    assert not promote, f'{len(promote)} promoted records not found'
    pysam.sort('-o', out, out + '.unsorted.bam'); os.remove(out + '.unsorted.bam'); pysam.index(out)
    print('two_step cond: ' + '; '.join(f'{k} {v}' for k, v in sorted(n.items())))

if __name__ == '__main__':
    {'layers': layers, 'unique': unique, 'cond': cond}[sys.argv[1]](*sys.argv[2:])
