#!/usr/bin/env python3
"""PREREG TANDEM CUT (register 1516, docs/PREREG_tandem_cut_2026-11-03.md): Rule T, a `--bridge-units-list` detector.
usage: tandem_cut.py PREFIX.gtf GENOME.fa OUT.units.tsv [--report OUT.report.tsv] [--threads N]
Per gene_id and strand: cut point p = an intron midpoint; LEFT / RIGHT / CROSSERS by position; MINORITY = reads(crossers) <
min(reads(LEFT), reads(RIGHT)); HOMOLOGY = LEFT vs RIGHT exon-union sequences, minimap2 -c -x asm20 best hit by matches with
identity >= 0.80 and aligned span >= 0.50 of the shorter; the admissible p with the largest min side wins, every crosser is cut
at its intron containing p, then LEFT and RIGHT recurse (depth <= 3). Output: tid chrom strand junctions label gene."""
import sys, re, os, collections, subprocess, tempfile, argparse, pysam
ap = argparse.ArgumentParser(); ap.add_argument('gtf'); ap.add_argument('genome'); ap.add_argument('out')
ap.add_argument('--report'); ap.add_argument('--threads', type=int, default=4); ap.add_argument('--max-depth', type=int, default=3)
ap.add_argument('--mode', choices=['locus', 'neighbour'], default='locus', help='locus = Rule T (r1516); neighbour = Rule T2 (r1517, docs/PREREG_tandem_cut_neighbour_2026-11-03.md)')
ap.add_argument('--debug', action='store_true', help='print every admissible minority point with its homology')
ap.add_argument('--no-homology', action='store_true', help='skip the gate (counts only; never for a registered arm)')
a = ap.parse_args()
A = lambda k, s: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]
tx = {}
for line in open(a.gtf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 9: continue
    t = A('transcript_id', f[8])
    r = tx.setdefault(t, {'g': A('gene_id', f[8]), 'chrom': f[0], 'strand': f[6], 'ex': [], 'reads': None})
    if f[2] == 'transcript':
        v = A('reads', f[8]); r['reads'] = int(v) if v is not None and v.isdigit() else None
    elif f[2] == 'exon': r['ex'].append((int(f[3]), int(f[4])))
for r in tx.values(): r['ex'].sort()
loci = collections.defaultdict(list)
for t, r in tx.items(): loci[(r['g'], r['chrom'], r['strand'])].append(t)
fa = pysam.FastaFile(a.genome)
comp = str.maketrans('ACGTacgt', 'TGCAtgca')
def union_seq(ts, chrom, strand):
    iv = sorted(e for t in ts for e in tx[t]['ex']); m = []
    for s, e in iv:
        if m and s <= m[-1][1] + 1: m[-1][1] = max(m[-1][1], e)
        else: m.append([s, e])
    seq = ''.join(fa.fetch(chrom, s - 1, e).upper() for s, e in m)
    return seq.translate(comp)[::-1] if strand == '-' else seq
hom_cache = {}
def homology(left, right, chrom, strand):
    key = (tuple(sorted(left)), tuple(sorted(right)))
    if key in hom_cache: return hom_cache[key]
    L, R = union_seq(left, chrom, strand), union_seq(right, chrom, strand)
    shorter = min(len(L), len(R)); best = (0.0, 0.0)
    if shorter >= 100:
        with tempfile.TemporaryDirectory() as d:
            q, t = os.path.join(d, 'q.fa'), os.path.join(d, 't.fa')
            open(q, 'w').write(f'>L\n{L}\n'); open(t, 'w').write(f'>R\n{R}\n')
            paf = subprocess.run(['minimap2', '-c', '-x', 'asm20', '-t', str(a.threads), t, q], capture_output=True, text=True).stdout
        bm = -1
        for line in paf.splitlines():
            f = line.split('\t'); matches, blk = int(f[9]), int(f[10]); span = int(f[3]) - int(f[2])
            if matches > bm: bm = matches; best = (matches / blk if blk else 0.0, span / shorter)
    hom_cache[key] = best; return best
def reads(ts): return sum(tx[t]['reads'] for t in ts)
cuts = collections.defaultdict(set); report = []; n_loci_tested = n_minority = n_hom = 0
def solve(ts, chrom, strand, gene, depth):
    global n_loci_tested, n_minority, n_hom
    if depth > a.max_depth or len(ts) < 2 or any(tx[t]['reads'] is None for t in ts): return
    cands = set()
    for t in ts:
        ex = tx[t]['ex']
        for (a1, a2), (b1, b2) in zip(ex, ex[1:]): cands.add((a2 + b1) // 2)
    if not cands: return
    n_loci_tested += 1; best = None
    for p in sorted(cands):
        left = [t for t in ts if tx[t]['ex'][-1][1] < p]; right = [t for t in ts if tx[t]['ex'][0][0] > p]
        cross = [t for t in ts if t not in left and t not in right]
        if not left or not right: continue
        if any(not any(x[1] < p < y[0] for x, y in zip(tx[t]['ex'], tx[t]['ex'][1:])) for t in cross): continue
        rl, rr, rc = reads(left), reads(right), reads(cross)
        if not rc < min(rl, rr): continue
        n_minority += 1
        ident, cov = (1.0, 1.0) if a.no_homology else homology(left, right, chrom, strand)
        if a.debug: print(f'[debug] {gene} p={p} left {len(left)}/{rl}r right {len(right)}/{rr}r cross {len(cross)}/{rc}r identity {ident:.3f} cov {cov:.3f}', file=sys.stderr)
        if ident < 0.80 or cov < 0.50: continue
        n_hom += 1
        score = min(rl, rr)
        if best is None or score > best[0]: best = (score, p, left, right, cross, ident, rl, rr, rc)
    if best is None: return
    score, p, left, right, cross, ident, rl, rr, rc = best
    share = rc / (rc + min(rl, rr))
    for t in cross:
        ex = tx[t]['ex']
        for x, y in zip(ex, ex[1:]):
            if x[1] < p < y[0]: cuts[t].add((x[1] + 1, y[0] - 1)); break
    report.append((gene, chrom, strand, p, len(left), rl, len(right), rr, len(cross), rc, f'{ident:.4f}', f'{share:.4f}', depth))
    solve(left, chrom, strand, gene, depth + 1); solve(right, chrom, strand, gene, depth + 1)
def neighbour_mode():
    """Rule T2: per junction J, CROSSERS = transcripts using J; UP / DOWN = same-strand non-users ending before J / starting after J
    that share an exonic base with a crosser's exons on that side (any gene_id)."""
    global n_loci_tested, n_minority, n_hom
    import bisect
    exidx = collections.defaultdict(list)   # (chrom, strand) -> sorted (start, end, tid)
    users = collections.defaultdict(list)   # (chrom, strand, s, e) -> tids
    for t, r in tx.items():
        if r['strand'] not in '+-': continue
        for x in r['ex']: exidx[(r['chrom'], r['strand'])].append((x[0], x[1], t))
        for x, y in zip(r['ex'], r['ex'][1:]): users[(r['chrom'], r['strand'], x[1] + 1, y[0] - 1)].append(t)
    for k in exidx: exidx[k].sort()
    starts = {k: [v[0] for v in lst] for k, lst in exidx.items()}
    def overlapping(chrom, strand, a, b):
        lst = exidx[(chrom, strand)]; i = bisect.bisect_right(starts[(chrom, strand)], b)
        return {tid for (s0, e0, tid) in lst[max(0, i - 20000):i] if e0 >= a}
    for (chrom, strand, s, e), cross in sorted(users.items()):
        if any(tx[t]['reads'] is None for t in cross): continue
        up = set(); down = set()
        for t in cross:
            for x in tx[t]['ex']:
                if x[1] < s: up |= overlapping(chrom, strand, x[0], x[1])
                elif x[0] > e: down |= overlapping(chrom, strand, x[0], x[1])
        cs = set(cross)
        up = [t for t in up if t not in cs and tx[t]['ex'][-1][1] < s and tx[t]['reads'] is not None]
        down = [t for t in down if t not in cs and tx[t]['ex'][0][0] > e and tx[t]['reads'] is not None]
        if not up or not down: continue
        n_loci_tested += 1
        rl, rr, rc = reads(up), reads(down), reads(cross)
        if not rc < min(rl, rr): continue
        n_minority += 1
        ident, cov = (1.0, 1.0) if a.no_homology else homology(up, down, chrom, strand)
        if a.debug: print(f'[debug] {chrom}:{s}-{e}:{strand} genes {sorted({tx[t]["g"] for t in cross})} up {len(up)}/{rl}r down {len(down)}/{rr}r cross {len(cross)}/{rc}r identity {ident:.3f} cov {cov:.3f}', file=sys.stderr)
        if ident < 0.80 or cov < 0.50: continue
        n_hom += 1
        for t in cross: cuts[t].add((s, e))
        report.append((','.join(sorted({tx[t]['g'] for t in cross})), chrom, strand, f'{s}-{e}', len(up), rl, len(down), rr, len(cross), rc, f'{ident:.4f}', f'{rc / (rc + min(rl, rr)):.4f}', 1))
if a.mode == 'neighbour': neighbour_mode()
else:
    for (g, chrom, strand), ts in sorted(loci.items()):
        if strand in '+-' and len(ts) >= 2: solve(ts, chrom, strand, g, 1)
with open(a.out, 'w') as w:
    w.write('tid\tchrom\tstrand\tjunctions\tlabel\tgene\n')
    for t in sorted(cuts):
        w.write(f"{t}\t{tx[t]['chrom']}\t{tx[t]['strand']}\t{';'.join(f'{s}-{e}' for s, e in sorted(cuts[t]))}\t{'tandem2' if a.mode == 'neighbour' else 'tandem'}\t{tx[t]['g']}\n")
if a.report:
    with open(a.report, 'w') as w:
        w.write('gene\tchrom\tstrand\tp\tn_left\treads_left\tn_right\treads_right\tn_cross\treads_cross\tidentity\tshare\tdepth\n')
        for r in report: w.write('\t'.join(map(str, r)) + '\n')
print(f'loci tested {n_loci_tested}; candidate points passing minority {n_minority}, passing homology too {n_hom}; '
      f'cuts {len(report)} on {len(cuts)} transcripts -> {a.out}')
