#!/usr/bin/env python3
"""ms ties vs O2's AS tie window on a real BAM (2026-10-30, user: "does it work for tied reads, not only AS but also ms since that is
the true tag that determines primaries from secondaries?"). minimap2 ranks hits by ms (register 1409); O2's pile window is AS
within AS_TIE (10) of the read's best AS (o3port.near). This measures, per read with >= 2 non-supplementary records, how the two
tie definitions overlap, genome-wide and restricted to reads whose records fall on family copies.
  tie_ms_vs_as.py BAM [COPIES_TSV] > report.tsv
Classes per read: ms_tie = its top two records have equal ms; as_win = >= 2 records within AS_TIE of its best AS (O2's tied set).
Reported: counts of ms_tie & as_win, ms_tie & not as_win (coin-toss reads O2 never sees), as_win & not ms_tie (reads minimap2 ranked
that O2 still treats as tied), each also with the same NM / same intron count among the tied pair; for ms ties outside the AS
window, the AS gap and intron-count difference of the two ms-tied records."""
import sys, collections, csv, bisect
import pysam
AS_TIE = 10
bam = pysam.AlignmentFile(sys.argv[1])
copies = None
if len(sys.argv) > 2:
    copies = collections.defaultdict(list)
    for r in csv.DictReader(open(sys.argv[2]), delimiter='\t'):
        copies[r['chrom']].append((int(r['start']), int(r['end']), r['family_id']))
    for c in copies: copies[c].sort()
    starts = {c: [x[0] for x in v] for c, v in copies.items()}
def fam_of(chrom, s, e):
    if copies is None: return None
    v = copies.get(chrom, [])
    i = bisect.bisect_left(starts.get(chrom, []), s - 3_000_000)
    while i < len(v) and v[i][0] < e:
        if v[i][1] > s: return v[i][2]
        i += 1
    return None
recs = collections.defaultdict(list)   # read -> [(ms, AS, NM, introns, is_secondary, family)]
n = 0
for rec in bam.fetch(until_eof=True):
    if rec.is_unmapped or rec.is_supplementary or not rec.has_tag('ms'): continue
    ni = sum(1 for op, _ in rec.cigartuples or () if op == 3)
    recs[rec.query_name].append((rec.get_tag('ms'), rec.get_tag('AS'), rec.get_tag('NM'), ni, rec.is_secondary,
                                 fam_of(rec.reference_name, rec.reference_start, rec.reference_end)))
    n += 1
    if n % 5_000_000 == 0: print(f"{n} records", file=sys.stderr)
C = collections.Counter(); gaps = collections.Counter(); dint = collections.Counter()
for name, L in recs.items():
    if len(L) < 2: continue
    scope = ['all']
    fams = {x[5] for x in L if x[5]}
    if copies is not None and len(fams) == 1 and all(x[5] for x in L): scope.append('family')  # every record on one family's copies
    bym = sorted(L, key=lambda x: (-x[0], -x[1]))
    ms_tie = bym[0][0] == bym[1][0]
    best_as = max(x[1] for x in L)
    W = [x for x in L if x[1] >= best_as - AS_TIE]
    as_win = len(W) >= 2
    byas = sorted(L, key=lambda x: (-x[1], x[2]))
    for sc in scope:
        C[(sc, 'reads>=2')] += 1
        if ms_tie:
            C[(sc, 'ms_tie')] += 1
            same_nm = bym[0][2] == bym[1][2]; same_int = bym[0][3] == bym[1][3]
            if same_nm: C[(sc, 'ms_tie.same_nm')] += 1
            if same_int: C[(sc, 'ms_tie.same_introns')] += 1
            both_in = all(x[1] >= best_as - AS_TIE for x in bym[:2])
            if both_in: C[(sc, 'ms_tie.in_as_window')] += 1
            else:
                C[(sc, 'ms_tie.outside_as_window')] += 1
                g = abs(bym[0][1] - bym[1][1]); gaps[(sc, min(g // 10 * 10, 100))] += 1
                dint[(sc, min(abs(bym[0][3] - bym[1][3]), 5))] += 1
            if bym[0][1] == bym[1][1] and same_nm: C[(sc, 'ms_tie.exact_as_nm')] += 1
        if as_win:
            C[(sc, 'as_win')] += 1
            if not ms_tie:
                C[(sc, 'as_win.not_ms_tie')] += 1
                if byas[0][1] == byas[1][1] and byas[0][2] == byas[1][2]: C[(sc, 'as_win.not_ms_tie.hopeless_by_as_nm')] += 1
                prim = [x for x in L if not x[4]]
                if prim and prim[0][0] == bym[0][0]: C[(sc, 'as_win.not_ms_tie.primary_is_max_ms')] += 1
print("scope\tclass\treads")
for (sc, k), v in sorted(C.items()): print(f"{sc}\t{k}\t{v}")
print("scope\tAS_gap_bin(ms ties outside window)\treads")
for (sc, g), v in sorted(gaps.items()): print(f"{sc}\t{g}\t{v}")
print("scope\tintron_diff(ms ties outside window)\treads")
for (sc, d), v in sorted(dint.items()): print(f"{sc}\t{d}\t{v}")
