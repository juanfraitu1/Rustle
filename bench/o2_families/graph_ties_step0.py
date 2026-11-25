#!/usr/bin/env python3
"""Step 0 (register 1943, docs/STEP0_graph_ties_2026-11-25.md): are a read's tied places a property of its locus?

    graph_ties_step0.py P.gtf TX.sam BAM AS_TABLE OUT_PREFIX [--threads N] [--every 50]

P.gtf   assembly made from primary alignments only (its transcripts are the loci).
TX.sam  the spliced sequences of P.gtf's transcripts aligned to the genome with secondaries (`--eqx`): the edges of the
        family graph, one alignment per (transcript, copy place).
For each tied read (a second record within 10 `ms` of its best, as tools/tie_table.py) the read's span is carried from its
primary alignment into its transcript T and through every other alignment H of T back to the genome. Reported: how many of
the read's in-band secondary records such a projection recovers, the gap predicted from the copies' differences under the
read, and how often the projection calls a tie for a read that has none (every Nth untied read).
Writes OUT_PREFIX.summary.txt and OUT_PREFIX.records.tsv (one row per secondary record looked at)."""
import sys, os, re, zlib, bisect, collections, multiprocessing, argparse
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E

BAND = 10
G = {}  # shared with the forked workers

def load_tx(gtf):
    """-> tx {tid: (contig, strand, gene, start, end)}, index {(contig, strand): (starts, [(s, e, tid)], longest exon)}; 0-based half-open"""
    tx, ex = {}, collections.defaultdict(list)
    for l in open(gtf):
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
        if t not in tx: tx[t] = [f[0], f[6], re.search(r'gene_id "([^"]+)"', f[8]).group(1), s, e]
        else: tx[t][3] = min(tx[t][3], s); tx[t][4] = max(tx[t][4], e)
        ex[(f[0], f[6])].append((s, e, t))
    idx = {}
    for k, v in ex.items():
        v.sort(); idx[k] = ([x[0] for x in v], v, max(e - s for s, e, _ in v))
    return tx, idx

def tx_overlaps(idx, contig, strand, blocks):
    """exon bases each transcript shares with the blocks: Counter {tid: bases}"""
    c = collections.Counter(); k = idx.get((contig, strand))
    if not k: return c
    starts, v, longest = k
    for s, e in blocks:
        i = bisect.bisect_left(starts, s - longest)
        while i < len(v) and v[i][0] < e:
            o = min(e, v[i][1]) - max(s, v[i][0])
            if o > 0: c[v[i][2]] += o
            i += 1
    return c

def scan(contig):
    """records of the tied reads (in band) and of every Nth untied read: (name, secondary, contig, reverse, strand, blocks, ms, tied)"""
    import pysam
    out = []; ties, every = G['ties'], G['every']
    for r in pysam.AlignmentFile(G['bam']).fetch(contig):
        if r.is_unmapped or r.is_supplementary: continue
        q = r.query_name; best = ties.get(q)
        if best is None:
            if zlib.crc32(q.encode()) % every: continue
        elif r.is_secondary and r.get_tag('ms') < best - BAND: continue
        diff = None
        if not r.is_secondary:
            diff = set(); g = r.reference_start
            for op, n in r.cigartuples:
                if op == 8: diff.update(range(g, g + n))
                elif op == 1: diff.update((g - 1, g))
                elif op == 2: diff.update(range(g - 1, g + n + 1))
                if op in (0, 2, 3, 7, 8): g += n
        out.append((q, r.is_secondary, contig, r.is_reverse, E.tx_strand(r), r.get_blocks(), r.get_tag('ms'), best is not None, diff))
    return out

def segments(ref_start, cigar):
    """-> [(op, q0, q1, r0, r1)] over the whole query (clips included), and the query length"""
    q, r, seg = 0, ref_start, []
    for op, n in cigar:
        if op in (0, 7, 8): seg.append((op, q, q + n, r, r + n)); q += n; r += n
        elif op == 1: seg.append((1, q, q + n, r, r)); q += n
        elif op == 2: seg.append((2, q, q, r, r + n)); r += n
        elif op == 3: seg.append((3, q, q, r, r + n)); r += n
        elif op in (4, 5): q += n
    return seg, q

def query_span(seg, gs, ge):
    """query interval of the aligned bases inside the genome interval [gs, ge)"""
    qa = qb = None
    for op, q0, q1, r0, r1 in seg:
        if op in (0, 7, 8) and r0 < ge and r1 > gs:
            a = q0 + max(0, gs - r0); b = q1 - max(0, r1 - ge)
            qa = a if qa is None else min(qa, a); qb = b if qb is None else max(qb, b)
    return (qa, qb) if qa is not None and qb > qa else None

def project(seg, qa, qb):
    """the query interval through one alignment: -> genome blocks, cost of the differences inside it (see the doc)"""
    blocks, cost, covered, ev = [], 0, 0, []
    for op, q0, q1, r0, r1 in seg:
        if op in (0, 7, 8):
            a, b = max(q0, qa), min(q1, qb)
            if b <= a: continue
            covered += b - a
            if op == 8: cost += 5 * (b - a); ev += [(x, 5) for x in range(a, b)]
            s, e = r0 + (a - q0), r0 + (b - q0)
            if blocks and blocks[-1][1] == s: blocks[-1][1] = e
            else: blocks.append([s, e])
        elif op == 1:
            n = min(q1, qb) - max(q0, qa)
            if n > 0: covered += n; cost += min(6 + n, 24) + n; ev.append((max(q0, qa), min(6 + n, 24) + n))
        elif op == 2:
            if qa < q0 < qb:
                cost += min(6 + (r1 - r0), 24); ev.append((q0, min(6 + (r1 - r0), 24)))
                if blocks and blocks[-1][1] == r0: blocks[-1][1] = r1  # a deletion does not split the block
    return blocks, cost + (qb - qa - covered), ev

def q2r(seg, q):
    """genome position of a query position in one alignment (None when not aligned)"""
    for op, q0, q1, r0, r1 in seg:
        if op in (0, 7, 8) and q0 <= q < q1: return r0 + (q - q0)
    return None

def ovl(a, b):
    n = 0
    for s, e in a:
        for x, y in b:
            o = min(e, y) - max(s, x)
            if o > 0: n += o
    return n

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('gtf'); ap.add_argument('txsam'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out')
    ap.add_argument('--threads', type=int, default=16); ap.add_argument('--every', type=int, default=50)
    a = ap.parse_args()
    import pysam
    tx, idx = load_tx(a.gtf)
    G.update(bam=a.bam, every=a.every, ties=E.load_ties(a.tab, E.floor_fn(BAND), 'ms'))
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        res = pool.map(scan, list(pysam.AlignmentFile(a.bam).references), chunksize=1)
    reads = collections.defaultdict(lambda: [None, [], False])  # name -> [primary, secondaries, tied]
    for part in res:
        for q, sec, contig, rev, strand, blocks, ms, tied, diff in part:
            d = reads[q]; d[2] = tied
            if sec: d[1].append((contig, rev, strand, blocks, ms))
            else: d[0] = (contig, rev, strand, blocks, ms, diff)
    del res
    # each read's transcripts: the best by shared exon bases, and up to 5 of its locus
    choice = {}; need = set()
    for q, (p, secs, tied) in reads.items():
        if p is None: continue
        c = tx_overlaps(idx, p[0], p[2], p[3])
        if not c: continue
        best = max(c, key=lambda t: (c[t], t)); gene = tx[best][2]
        choice[q] = [best] + [t for t, _ in c.most_common(6) if t != best and tx[t][2] == gene][:4]
        need.update(choice[q])
    hits = collections.defaultdict(list)  # tid -> [(contig, reverse, segments, qlen)]
    for r in pysam.AlignmentFile(a.txsam):
        if r.is_unmapped or r.is_supplementary or r.query_name not in need: continue
        seg, qlen = segments(r.reference_start, r.cigartuples)
        hits[r.query_name].append((r.reference_name, r.is_reverse, seg, qlen, r.reference_start, r.reference_end))
    selfhit = {}
    for t, hs in hits.items():
        contig, strand, _, s, e = tx[t]
        cand = [(min(e, h[5]) - max(s, h[4]), i) for i, h in enumerate(hs) if h[0] == contig and h[1] == (strand == '-')]
        cand = [c for c in cand if c[0] > 0]
        if cand: selfhit[t] = max(cand)[1]

    def projections(t, p):
        """every other alignment of transcript t under the read's primary p: [(contig, flip, blocks, predicted gap)]"""
        if t not in selfhit: return None
        hs = hits[t]; sh = hs[selfhit[t]]
        span = query_span(sh[2], p[3][0][0], p[3][-1][1])
        if span is None: return None
        _, c0, _ = project(sh[2], *span); out = []; rd = p[5]; ql = sh[3]
        for i, h in enumerate(hs):
            if i == selfhit[t]: continue
            qa, qb = span if h[1] == sh[1] else (h[3] - span[1], h[3] - span[0])
            bl, c, ev = project(h[2], qa, qb)
            if not bl: continue
            if h[0] == p[0] and ovl(bl, p[3]) > 0.5 * sum(e - s for s, e in p[3]): continue  # the read's own place
            # read-aware gap (post hoc arm): where the copies differ and the read itself differs from its own place at that
            # base, the read is taken to carry the other copy's base, so the difference counts for the copy, not against it
            ra = c - c0
            for x, cost in ev:
                g = q2r(sh[2], x if h[1] == sh[1] else ql - 1 - x)
                if g is not None and g in rd: ra -= 2 * cost
            out.append((h[0], h[1] != sh[1], bl, c - c0, max(0, ra)))
        return out

    n = collections.Counter(); cal = collections.Counter(); rows = []
    for q, (p, secs, tied) in reads.items():
        kind = 'tied' if tied else 'control'
        if p is None: n[kind, 'no primary record'] += 1; continue
        if q not in choice: n[kind, 'primary on no transcript'] += 1; continue
        pr = projections(choice[q][0], p)
        if pr is None: n[kind, 'transcript has no alignment at its own place'] += 1; continue
        n[kind, 'reads'] += 1
        alt = None
        if any(x[3] <= BAND for x in pr): n[kind, 'reads with a predicted tie'] += 1
        if any(x[4] <= BAND for x in pr): n[kind, 'reads with a predicted tie, read-aware'] += 1
        used = set()
        for contig, rev, strand, blocks, ms in secs:
            gap = max(0, (G['ties'][q] if tied else p[4]) - ms); alen = sum(e - s for s, e in blocks); flip = rev != p[1]
            def match(pp):
                best = None
                for j, (c, f, bl, g, g2) in enumerate(pp or ()):
                    if c == contig and f == flip:
                        o = ovl(bl, blocks)
                        if o >= 0.5 * alen and (best is None or o > best[0]): best = (o, g, j, g2)
                return best
            m = match(pr); cls = lambda g: '0' if g <= 0 else ('1-10' if g <= BAND else '>10')
            if tied:
                n['tied', 'in-band secondary records'] += 1
                if m:
                    n['tied', 'recovered'] += 1; used.add(m[2])
                    if m[1] <= BAND: n['tied', 'recovered with predicted gap <= 10'] += 1
                    n['tied', f'actual {"0" if gap == 0 else "1-10"} / predicted {cls(m[1])}'] += 1
                    n['tied', f'read-aware: actual {"0" if gap == 0 else "1-10"} / predicted {cls(m[3])}'] += 1
                    if m[3] <= BAND: n['tied', 'recovered with read-aware gap <= 10'] += 1
                else:
                    if alt is None: alt = [projections(t, p) for t in choice[q][1:]]
                    if any(match(x) for x in alt): n['tied', 'recovered only by another transcript of the locus'] += 1
            elif m:
                cal[cls(gap), cls(m[1])] += 1
            rows.append((q, kind, choice[q][0], contig, blocks[0][0], blocks[-1][1], gap, m[1] if m else 'NA'))
        if tied: n['tied', 'predicted ties with no in-band record of the read'] += sum(1 for j, x in enumerate(pr) if x[3] <= BAND and j not in used)
    with open(a.out + '.records.tsv', 'w') as w:
        w.write('read\tkind\ttranscript\tcontig\tstart\tend\tactual_gap\tpredicted_gap\n'); w.writelines('\t'.join(map(str, r)) + '\n' for r in rows)
    T = lambda k: n['tied', k]; C = lambda k: n['control', k]; tied_all = sum(v for (k, _), v in n.items() if k == 'tied' and _ in ('reads', 'no primary record', 'primary on no transcript', 'transcript has no alignment at its own place'))
    fr = lambda x, y: f'{x} of {y} ({x / y:.3f})' if y else f'{x} of 0'
    L = [f'tied reads {tied_all}; control reads (every {a.every}th untied) {sum(v for (k, kk), v in n.items() if k == "control" and kk in ("reads", "no primary record", "primary on no transcript", "transcript has no alignment at its own place"))}']
    for k, v in sorted(n.items()): L.append(f'  {k[0]}\t{k[1]}\t{v}')
    ra = T('reads'); rb = T('recovered'); rd = C('reads with a predicted tie') * a.every
    L += ['', f'(a) reach: tied reads with a usable transcript {fr(ra, tied_all)}  bar >= 0.70',
          f'(b) recall: in-band secondary records recovered {fr(rb, T("in-band secondary records"))}  bar >= 0.90',
          f'    with any of up to 5 transcripts of the locus {fr(rb + T("recovered only by another transcript of the locus"), T("in-band secondary records"))}',
          f'(c) band: recovered records with predicted gap <= 10 {fr(T("recovered with predicted gap <= 10"), rb)}  bar >= 0.90',
          f'(d) false ties: control reads with a predicted tie {C("reads with a predicted tie")} of {C("reads")}, x {a.every} = {rd}; against tied reads {ra}: ratio {rd / ra if ra else float("nan"):.2f}  bar <= 0.20',
          f'    tied reads with a predicted tie {fr(T("reads with a predicted tie"), ra)}',
          f'POST HOC, read-aware gap: (c) {fr(T("recovered with read-aware gap <= 10"), rb)}; (d) control reads with a predicted tie {C("reads with a predicted tie, read-aware")}, x {a.every} = {C("reads with a predicted tie, read-aware") * a.every}, ratio {C("reads with a predicted tie, read-aware") * a.every / ra if ra else 0:.2f}',
          '', 'control reads, secondary records matched to a projection: actual gap class x predicted gap class']
    for k, v in sorted(cal.items()): L.append(f'  actual {k[0]}\tpredicted {k[1]}\t{v}')
    open(a.out + '.summary.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))

if __name__ == '__main__': main()
