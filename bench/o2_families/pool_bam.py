#!/usr/bin/env python3
"""A BAM holding only the secondary alignments a seeding-pool rule admits (register 1942).

    pool_bam.py IN.bam AS_TABLE OUT.bam --rule as:0.98 | ms:0.98 | msgap:10 | asintron:0.98:10:30  [--threads N]

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
A secondary without the tag, or of a read missing from the table, is kept (as the assembler does). Reads without a
position (unmapped) are not copied: the assembler does not read them."""
import sys, os, re, argparse, multiprocessing
G = {}

def work(job):
    import pysam
    i, contig = job; tag = 'AS' if G['rule'][0] == 'as' else 'ms'; part = f"{G['out']}.part.{i:06d}.bam"; kept = seen = 0; kind, val = G['rule']; best = G['best']
    B = pysam.AlignmentFile(G['bam'])
    with pysam.AlignmentFile(part, 'wb', template=B) as w:
        for r in B.fetch(contig):
            if r.is_secondary:
                seen += 1; b = best.get(r.query_name)
                if kind == 'asintron':
                    if b is not None and b[0] is not None and r.has_tag('AS'):
                        a = r.get_tag('AS'); R, Gp, Pn = val
                        if a < R * b[0] and not (b[1] is not None and r.has_tag('ms') and r.get_tag('ms') >= b[1] - Gp
                                                  and a >= R * b[0] - Pn * sum(1 for o, _ in r.cigartuples if o == 3)): continue
                elif b is not None and r.has_tag(tag):
                    t = r.get_tag(tag); top = b[0] if kind == 'as' else b[1]
                    if top is not None and (t < val * top if kind in ('as', 'ms') else t < top - val): continue
                kept += 1
            w.write(r)
    return part, seen, kept

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
    kind, val = rule.split(':', 1); assert kind in ('as', 'ms', 'msgap', 'asintron'), rule
    val = tuple(map(float, val.split(':'))) if kind == 'asintron' else float(val)
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
