#!/usr/bin/env python3
"""A BAM holding only the secondary alignments a seeding-pool rule admits (register 1942).

    pool_bam.py IN.bam AS_TABLE OUT.bam --rule as:0.98 | ms:0.98 | msgap:10 | asintron:0.98:10:30 | anchor:ms:0.9:3 | janchor:ms:0.9:3  [--threads N]

The assembler's `good` pool admits a secondary alignment when its AS is at least 0.98 x the read's best AS over the whole
BAM (the best comes from the score table). To try another rule without changing the assembler, the rule is applied here
and the result is assembled with `--seed-pool all` (every secondary that is left). Primary and supplementary records are
always kept. Rules, each against the read's genome-wide best from AS_TABLE (as_table columns: best AS = 2, best ms = 7):
  as:R      AS >= R x best AS      (as:0.98 = the `good` pool; the control that this emulation is faithful)
  ms:R      ms >= R x best ms      (the same width on the score minimap2 ranks by)
  msgap:N   ms >= best ms - N      (msgap:10 = the near-tie band of the labels)
  asintron:R:G:P  the `good` pool PLUS the secondaries that tie on ms and lose on AS only by their introns (register
            1953): kept when AS >= R x best AS, or when ms >= best ms - G and AS >= R x best AS - P x (introns of this
            alignment). AS charges about 24 per intron and ms nothing, so a spliced alignment that fits as well as an
            intronless one, base for base, is 24 x introns behind on AS and falls out of the AS pool.
  anchor:S:R:K  (register 1961) the `good` pool PLUS the secondaries at >= R x the read's best on score S (as or ms) that
            land where the sample has reads of its own: at least 50 of the alignment's bases lie where >= K primary
            alignments with mapping quality >= 1 are stacked. A place nothing is placed at without a tie-break gets none.
  janchor:S:R:K  as anchor, and every intron of the secondary must be an intron of one of those primary alignments (same
            two ends): the secondary can add support to a structure the place's own reads show, never a new one.
A secondary without the tag, or of a read missing from the table, is kept (as the assembler does). Reads without a
position (unmapped) are not copied: the assembler does not read them."""
import sys, os, re, argparse, multiprocessing
G = {}

def work(job):
    import pysam
    i, contig = job; tag = 'AS' if G['rule'][0] == 'as' else 'ms'; part = f"{G['out']}.part.{i:06d}.bam"; kept = seen = 0; kind, val = G['rule']; best = G['best']
    B = pysam.AlignmentFile(G['bam']); deep = None
    if kind in ('anchor', 'janchor'):
        tag = 'AS'; own = set() if kind == 'janchor' else None; deep = anchored_bases(B, contig, val[2], own)
    with pysam.AlignmentFile(part, 'wb', template=B) as w:
        for r in B.fetch(contig):
            if r.is_secondary:
                seen += 1; b = best.get(r.query_name)
                if kind == 'asintron':
                    if b is not None and b[0] is not None and r.has_tag('AS'):
                        a = r.get_tag('AS'); R, Gp, Pn = val
                        if a < R * b[0] and not (b[1] is not None and r.has_tag('ms') and r.get_tag('ms') >= b[1] - Gp
                                                  and a >= R * b[0] - Pn * sum(1 for o, _ in r.cigartuples if o == 3)): continue
                elif kind in ('anchor', 'janchor'):
                    if b is not None and b[0] is not None and r.has_tag('AS') and r.get_tag('AS') < 0.98 * b[0]:
                        s = val[0]; top = b[0] if s == 'as' else b[1]; t = r.get_tag('AS' if s == 'as' else 'ms') if r.has_tag('AS' if s == 'as' else 'ms') else None
                        if top is None or t is None or t < val[1] * top or sum(int(deep[y] - deep[x]) for x, y in r.get_blocks()) < 50: continue
                        if own is not None and not all(j in own for j in junctions(r)): continue
                elif b is not None and r.has_tag(tag):
                    t = r.get_tag(tag); top = b[0] if kind == 'as' else b[1]
                    if top is not None and (t < val * top if kind in ('as', 'ms') else t < top - val): continue
                kept += 1
            w.write(r)
    return part, seen, kept

def junctions(r):
    """the introns (N operations) of an alignment as (start, end)"""
    out = []; p = r.reference_start
    for o, n in r.cigartuples:
        if o == 3: out.append((p, p + n))
        if o in (0, 2, 3, 7, 8): p += n
    return out

def anchored_bases(B, contig, k, own=None):
    """prefix sums over the contig of 'at least k primary alignments with mapping quality >= 1 cover this base'"""
    import numpy as np
    d = np.zeros(B.get_reference_length(contig) + 2, dtype=np.int32)
    for r in B.fetch(contig):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 1: continue
        for x, y in r.get_blocks(): d[x] += 1; d[y] -= 1
        if own is not None: own.update(junctions(r))
    np.cumsum(d, out=d); d = (d >= k).astype(np.int32); np.cumsum(d, out=d)
    return np.concatenate(([0], d))            # deep[y] - deep[x] = anchored bases in [x, y)

def load_best(tab, bam=None):
    best = {}
    for l in open(tab):
        if l[0] == '#':
            m = re.search(r'\tbam_bytes=(\d+)', l)          # as_table's header records its BAM's size: another BAM's scores filter silently wrong
            if m and bam and int(m.group(1)) != os.path.getsize(bam): sys.exit(f'pool_bam: {tab} was made from a BAM of {m.group(1)} bytes; {bam} has {os.path.getsize(bam)}')
            continue
        f = l.rstrip('\n').split('\t')
        a = int(f[1]) if f[1].lstrip('-').isdigit() else None; m = int(f[6]) if len(f) > 6 and f[6].lstrip('-').isdigit() else None
        best[f[0]] = (a, m)
    return best

def run(bam, tab, out, rule, threads):
    import pysam
    kind, val = rule.split(':', 1); assert kind in ('as', 'ms', 'msgap', 'asintron', 'anchor', 'janchor'), rule
    if kind in ('anchor', 'janchor'): s, R, K = val.split(':'); assert s in ('as', 'ms'), rule; val = (s, float(R), int(K))
    else: val = tuple(map(float, val.split(':'))) if kind == 'asintron' else float(val)
    G.update(bam=bam, out=out, rule=(kind, val), best=load_best(tab, bam))
    if kind != 'as' and not any(m is not None for _, m in G['best'].values()):      # an AS-only table: an ms rule would keep everything and say nothing
        sys.exit(f'pool_bam: {tab} has no ms columns; --rule {rule} cannot be applied (rebuild the table with as_table)')
    refs = list(pysam.AlignmentFile(bam).references)
    with multiprocessing.get_context('fork').Pool(threads) as pool: res = pool.map(work, list(enumerate(refs)), chunksize=1)
    parts = [p for p, _, _ in res]; pysam.cat('-o', out + '.tmp', *parts)
    for p in parts: os.remove(p)
    os.replace(out + '.tmp', out); pysam.index(out)
    seen = sum(s for _, s, _ in res); kept = sum(k for _, _, k in res)
    print(f'pool_bam ({rule}): {kept} of {seen} secondary alignments kept ({kept / max(1, seen):.3f}) -> {out}')
    return seen, kept

def selftest():
    import tempfile, pysam
    d = tempfile.mkdtemp(); bam = os.path.join(d, 'x.bam'); tab = os.path.join(d, 't.tsv')
    hdr = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c', 'LN': 9000}]}
    def rd(name, pos, flag, AS, ms):
        r = pysam.AlignedSegment(); r.query_name = name; r.reference_id = 0; r.reference_start = pos; r.cigarstring = '100M'; r.flag = flag; r.mapping_quality = 60; r.query_sequence = 'A' * 100
        r.set_tag('AS', AS); r.set_tag('ms', ms); return r
    with pysam.AlignmentFile(bam, 'wb', header=hdr) as w:
        for r in (rd('r', 100, 0, 1000, 1000), rd('r', 1000, 256, 985, 1000), rd('r', 2000, 256, 970, 995), rd('r', 3000, 256, 990, 900), rd('q', 4000, 256, 10, 10)): w.write(r)
    pysam.index(bam); open(tab, 'w').write('#as_table\nr\t1000\t990\t4\tc\t1000\t1000\t1000\t1000\n')
    names = lambda p: [(r.reference_start, r.is_secondary) for r in pysam.AlignmentFile(p)]
    assert run(bam, tab, os.path.join(d, 'a.bam'), 'as:0.98', 1) == (4, 3) and [x[0] for x in names(os.path.join(d, 'a.bam'))] == [100, 1000, 3000, 4000]
    assert run(bam, tab, os.path.join(d, 'm.bam'), 'ms:0.98', 1) == (4, 3) and [x[0] for x in names(os.path.join(d, 'm.bam'))] == [100, 1000, 2000, 4000]
    assert run(bam, tab, os.path.join(d, 'g.bam'), 'msgap:10', 1) == (4, 3) and [x[0] for x in names(os.path.join(d, 'g.bam'))] == [100, 1000, 2000, 4000]
    assert run(bam, tab, os.path.join(d, 'g0.bam'), 'msgap:0', 1) == (4, 2)
    open(tab, 'w').write('#as_table\tbam=/some/other.bam\tbam_bytes=1\n')
    try: load_best(tab, bam); raise AssertionError('a table of another BAM was accepted')
    except SystemExit: pass
    open(tab, 'w').write('#as_table\nr\t1000\t990\t4\tc\t1000\t1000\t1000\t1000\n')
    def sp(name, pos, AS, ms, cigar):
        r = pysam.AlignedSegment(); r.query_name = name; r.reference_id = 0; r.reference_start = pos; r.cigarstring = cigar; r.flag = 256; r.mapping_quality = 0; r.query_sequence = 'A' * 100
        r.set_tag('AS', AS); r.set_tag('ms', ms); return r
    b2 = os.path.join(d, 'y.bam')
    with pysam.AlignmentFile(b2, 'wb', header=hdr) as w:      # a parent with 5 introns, 120 behind on AS, tied on ms: kept; the same without the ms tie, or without introns: dropped
        for r in (rd('r', 100, 0, 1000, 1000), sp('r', 1000, 880, 1000, '20M50N20M50N20M50N20M50N10M50N10M'), sp('r', 3000, 880, 950, '20M50N20M50N20M50N20M50N10M50N10M'), sp('r', 5000, 880, 1000, '100M'), sp('r', 7000, 985, 900, '100M')): w.write(r)
    pysam.index(b2)
    assert run(b2, tab, os.path.join(d, 'i.bam'), 'asintron:0.98:10:30', 1) == (4, 2) and [x[0] for x in names(os.path.join(d, 'i.bam'))] == [100, 1000, 7000]
    assert run(b2, tab, os.path.join(d, 'j.bam'), 'as:0.98', 1) == (4, 1)
    b3 = os.path.join(d, 'z.bam')           # anchor: secondaries at 0.93 of the best; one where 3 placed reads are stacked, one where 2 are, one where the 3 have mapping quality 0
    def pr(name, pos, mq):
        r = rd(name, pos, 0, 500, 500); r.mapping_quality = mq; return r
    with pysam.AlignmentFile(b3, 'wb', header=hdr) as w:
        for r in (rd('r', 100, 0, 1000, 1000), pr('a1', 1000, 60), pr('a2', 1000, 60), pr('a3', 1000, 60), rd('r', 1010, 256, 930, 930), pr('b1', 3000, 60), pr('b2', 3000, 60), rd('r', 3010, 256, 930, 930),
                  pr('c1', 5000, 0), pr('c2', 5000, 0), pr('c3', 5000, 0), rd('r', 5010, 256, 930, 880), rd('r', 7000, 256, 985, 900)): w.write(r)
    pysam.index(b3)
    assert run(b3, tab, os.path.join(d, 'k.bam'), 'anchor:as:0.9:3', 1) == (4, 2) and [x[0] for x in names(os.path.join(d, 'k.bam')) if x[1]] == [1010, 7000]
    b4 = os.path.join(d, 'w.bam')          # janchor: the place's reads show intron 1040-1100; a secondary with that intron is kept, one with another intron is not
    def spl(name, pos, flag, cigar, AS):
        r = rd(name, pos, flag, AS, AS); r.cigarstring = cigar; r.mapping_quality = 60 if flag == 0 else 0; return r
    with pysam.AlignmentFile(b4, 'wb', header=hdr) as w:
        for r in (rd('r', 100, 0, 1000, 1000), spl('a1', 1000, 0, '40M60N60M', 500), spl('a2', 1000, 0, '40M60N60M', 500), spl('a3', 1000, 0, '40M60N60M', 500), spl('r', 1000, 256, '40M60N60M', 930), spl('r', 1001, 256, '45M60N55M', 930)): w.write(r)
    pysam.index(b4)
    assert run(b4, tab, os.path.join(d, 'j1.bam'), 'janchor:as:0.9:3', 1) == (2, 1) and run(b4, tab, os.path.join(d, 'j2.bam'), 'anchor:as:0.9:3', 1) == (2, 2)
    assert run(b3, tab, os.path.join(d, 'k2.bam'), 'anchor:ms:0.9:2', 1) == (4, 3) and run(b3, tab, os.path.join(d, 'k3.bam'), 'anchor:ms:0.95:2', 1) == (4, 1)
    open(tab, 'w').write('#as_table\nr\t1000\t990\t4\tc\t1000\n')              # a table without ms columns
    try: run(bam, tab, os.path.join(d, 'n.bam'), 'msgap:10', 1); raise AssertionError('an ms rule ran on a table without ms')
    except SystemExit: pass
    assert run(bam, tab, os.path.join(d, 'n.bam'), 'as:0.98', 1) == (4, 3)
    print('pool_bam selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else:
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out'); ap.add_argument('--rule', required=True); ap.add_argument('--threads', type=int, default=8)
        a = ap.parse_args(); run(a.bam, a.tab, a.out, a.rule, a.threads)
