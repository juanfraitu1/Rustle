"""Read-level node-split preprocessor, per `docs/PREREG_node_split_2026-09-21.md` (md5 `7596b35a`).

Two triggers, both operating on PRIMARY reads BEFORE assembly ever builds one locus out of a chimeric
molecule -- avoiding register 846's two documented failure modes (doubling via transitivity; coverage
inflation from re-scoring a shorter piece), since the split fragments are assembled from scratch through
the ordinary, unmodified pipeline rather than being cut out of an already-built node.

TRIGGER 1 (chimeric-bridge junction): reuses the assembler's own `is_chimeric_bridge` test (a skeleton
is a bridge if it shares a junction with two others whose spans are mutually disjoint), applied per READ
instead of per skeleton (skeletons don't exist yet at this stage) -- a read is a bridge if, among OTHER
reads sharing ANY of its own junctions, two exist whose spans are disjoint from each other. Cut point is
the one junction separating that read's "left-neighbour" group from its "right-neighbour" group.

TRIGGER 2 (read-identity turnover / coverage cliff): at each of a read's own junctions, compare the set
of OTHER reads covering the flanking exon immediately before vs immediately after. Cut where that set
turns over almost completely (>= 0.90, chosen before looking).

Writes a MODIFIED BAM: every flagged read is replaced by two sub-alignments (soft-clip the other half),
mirroring `split_mischained_reads`'s exact mechanic; everything else passes through untouched.
"""
import argparse
import collections
import re
import sys

import pysam

CIGAR_OPS = re.compile(r'(\d+)([MIDNSHP=X])')


def load_reads(bam_path, chrom):
    """[(qname, ref_start_1b, introns, blocks, record)] -- introns as [(donor,acceptor)] 1-based
    inclusive of the intron itself (matches this session's established convention); blocks as the
    exon-consuming (start,end) pairs, same coordinate system."""
    bam = pysam.AlignmentFile(bam_path, 'rb')
    out = []
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        cigar = r.cigarstring
        if not cigar:
            continue
        cur = r.reference_start + 1  # 1-based
        introns, blocks = [], []
        bs = cur
        for n, op in CIGAR_OPS.findall(cigar):
            n = int(n)
            if op in 'M=X':
                cur += n
            elif op == 'D':
                cur += n
            elif op == 'N':
                blocks.append((bs, cur - 1))
                introns.append((cur, cur + n - 1))
                cur += n
                bs = cur
        blocks.append((bs, cur - 1))
        out.append((r.query_name, r.reference_start + 1, introns, blocks, r))
    return out


def build_junction_index(reads):
    by_junc = collections.defaultdict(list)
    for idx, (qn, rs, introns, blocks, rec) in enumerate(reads):
        for j in introns:
            by_junc[j].append(idx)
    return by_junc


def read_span(reads, idx):
    return reads[idx][3][0][0], reads[idx][3][-1][1]


MAX_JUNCTION_DEGREE = 300  # a genuine two-locus bridge needs only ~2 witnesses on each side; NPIP-class
# repeat junctions run into the THOUSANDS (measured: up to 4,666 reads sharing one junction on chr16's
# NPIP/PKD1 region) and add pure cost, not signal, to the disjoint-span test below.


def trigger1_chimeric_bridge(reads, by_junc):
    """-> {read_idx: cut_junction} for reads flagged as a chimeric bridge.

    The disjoint-span test is INTERVAL-MERGE based (O(n log n) in the neighbour count), not the naive
    all-pairs O(n^2) check `is_chimeric_bridge`'s own doc-comment describes -- a read's own junctions can
    have thousands of other reads attached in a tandem-repeat region (measured: up to 4,666 on chr16's
    NPIP/PKD1 window, with reads carrying up to 44 introns), and an all-pairs check over that many
    neighbours is not just slow, it is asymptotically intractable exactly where this trigger matters
    most. Semantically identical to the pairwise check: "some two neighbour spans are disjoint" iff
    the neighbours' spans, sorted and merged, do not collapse into ONE contiguous block."""
    cuts = {}
    for idx, (qn, rs, introns, blocks, rec) in enumerate(reads):
        if len(introns) < 2:
            continue  # need >=2 introns: a cut junction plus evidence on both sides
        neighbours = []  # (junction, span_start, span_end)
        for j in introns:
            others = by_junc.get(j, ())
            if len(others) > MAX_JUNCTION_DEGREE:
                continue  # a repeat-driven junction: no genuine two-locus signal in its degree
            for oidx in others:
                if oidx != idx:
                    s, e = read_span(reads, oidx)
                    neighbours.append((j, s, e))
        if len(neighbours) < 2:
            continue
        # sort by span start; merge overlapping spans; a read is a bridge iff >1 merged group remains,
        # and the cut point is the LATEST of this read's own junctions attached to the FIRST group.
        neighbours.sort(key=lambda x: x[1])
        first_group_end = neighbours[0][2]
        cut_junction = neighbours[0][0]
        bridge_junction = None
        for j, s, e in neighbours[1:]:
            if s < first_group_end:  # still inside the first merged group
                first_group_end = max(first_group_end, e)
                if introns.index(j) > introns.index(cut_junction):
                    cut_junction = j
                continue
            # a second, disjoint group starts here
            bridge_junction = cut_junction
            break
        if bridge_junction:
            cuts[idx] = bridge_junction
    return cuts


MAX_COVER_SCAN = 200   # repeat-region circuit breaker for the backward interval scan below
MAX_EXON_SPAN = 20000   # no ordinary single exon block is this long; bounds the backward scan distance


def trigger2_turnover(reads, floor=0.90):
    """-> {read_idx: cut_junction}. Builds a coordinate-sorted interval index of all reads' blocks
    once, reused for every exon-flank lookup (O(n log n) overall, not O(n^2))."""
    all_blocks = []
    for idx, (qn, rs, introns, blocks, rec) in enumerate(reads):
        for b in blocks:
            all_blocks.append((b[0], b[1], idx))
    all_blocks.sort()
    starts = [b[0] for b in all_blocks]
    import bisect

    def covering_reads(s, e):
        # `all_blocks` is start-sorted: every block overlapping [s,e] has start <= e, so bisecting to
        # just past e and walking backward finds them all -- capped at MAX_COVER_SCAN purely as a
        # repeat-region circuit breaker (measured need: a handful of iterations in ordinary exons; NPIP-
        # class pileups can have thousands of blocks starting in a few hundred bp, same pathology as
        # trigger1's junction degree). Returns (covering_set, truncated) -- ⚠a capped scan is a BIASED
        # partial sample, not a smaller true answer: measured directly on chr16's NPIP/PKD1 window, an
        # uncapped v0 flagged 34% of local reads, almost all from truncation-induced apparent turnover
        # (whichever ~200 reads the cap happens to include differ before vs after by sampling alone, not
        # by biology) -- so a truncated comparison must ABSTAIN, matching trigger1's abstain-on-saturated-
        # junction design, not silently answer from a biased sample.
        i = bisect.bisect_right(starts, e) - 1
        out = set()
        j = i
        scanned = 0
        while j >= 0 and scanned < MAX_COVER_SCAN:
            bs, be, ridx = all_blocks[j]
            if s - bs > MAX_EXON_SPAN:
                return out, False  # exhausted genuinely, not by the cap
            if be >= s:
                out.add(ridx)
            j -= 1
            scanned += 1
        return out, (j >= 0)  # cap hit while more candidates remained -> truncated

    cuts = {}
    for idx, (qn, rs, introns, blocks, rec) in enumerate(reads):
        if len(introns) < 2:
            continue
        for k in range(len(introns) - 1):
            exon_before = blocks[k]
            exon_after = blocks[k + 1]
            before, trunc_b = covering_reads(*exon_before)
            after, trunc_a = covering_reads(*exon_after)
            before -= {idx}; after -= {idx}
            if trunc_b or trunc_a:
                continue  # abstain: a capped scan cannot tell turnover from sampling bias
            if not before or not after:
                continue
            union = before | after
            turnover = 1 - len(before & after) / len(union)
            if turnover >= floor:
                cuts[idx] = introns[k]
                break
    return cuts


def write_split_bam(in_bam, out_bam, chrom, cuts, reads):
    """cuts: {read_idx: (donor,acceptor)} -- the ONE junction each flagged read is cut at."""
    bam_in = pysam.AlignmentFile(in_bam, 'rb')
    bam_out = pysam.AlignmentFile(out_bam, 'wb', template=bam_in)
    cut_by_name = {}
    for idx, junc in cuts.items():
        cut_by_name[reads[idx][0]] = (reads[idx][4], junc)

    n_split = 0
    for r in bam_in.fetch(chrom):
        if r.query_name in cut_by_name and not r.is_secondary and not r.is_supplementary:
            _, (d, a) = cut_by_name[r.query_name]
            left, right = split_record(r, d, a)
            if left and right:
                bam_out.write(left); bam_out.write(right)
                n_split += 1
                continue
        bam_out.write(r)
    bam_in.close(); bam_out.close()
    print(f'{chrom}: {n_split} reads split', file=sys.stderr)


def split_record(r, donor, acceptor):
    """Split ONE alignment record at intron (donor,acceptor) into a left and right sub-record, each
    keeping its own genomic-consuming CIGAR and soft-clipping the rest of the query -- mirrors
    `split_mischained_reads`'s cut-and-keep-both-flanks mechanic exactly, at the SAM-record level.

    A valid CIGAR carries S/H ONLY at its two ends (never embedded), so any pre-existing clip is
    stripped from `ops` FIRST and re-attached to whichever new record owns that end -- an earlier
    version left a leading 'S' embedded in `left_ops` via the generic op branch without advancing the
    query-position counter `q`, silently undercounting `cut_q` by the clip length and writing a CIGAR
    whose query-consuming length did not match `query_sequence` (`samtools sort` refused every such
    record: "CIGAR and query sequence lengths differ")."""
    cigar = r.cigarstring
    if not cigar:
        return None, None
    ops = [(int(n), op) for n, op in CIGAR_OPS.findall(cigar)]
    orig_lead = ops[0] if ops and ops[0][1] in 'SH' else None
    orig_trail = ops[-1] if ops and ops[-1][1] in 'SH' and len(ops) > (1 if orig_lead else 0) else None
    core = ops[(1 if orig_lead else 0):(len(ops) - (1 if orig_trail else 0))]

    cur_ref = r.reference_start + 1
    q = 0
    left_ops, right_ops = [], []
    cut_q = None
    for n, op in core:
        if op in 'M=X':
            (left_ops if cut_q is None else right_ops).append((n, op))
            cur_ref += n; q += n
        elif op == 'D':
            (left_ops if cut_q is None else right_ops).append((n, op))
            cur_ref += n
        elif op == 'I':
            (left_ops if cut_q is None else right_ops).append((n, op))
            q += n
        elif op == 'N':
            if cur_ref == donor and cur_ref + n - 1 == acceptor and cut_q is None:
                cut_q = q  # this is the cut junction: stop accumulating into `left`
                cur_ref += n
                continue
            (left_ops if cut_q is None else right_ops).append((n, op))
            cur_ref += n
        else:
            (left_ops if cut_q is None else right_ops).append((n, op))
    if cut_q is None or not left_ops or not right_ops:
        return None, None

    def build(sub_ops, q_end_of_this_side, ref_start, is_left):
        a = pysam.AlignedSegment(r.header)
        a.query_name = r.query_name + ('/L' if is_left else '/R')
        a.query_sequence = r.query_sequence
        a.query_qualities = r.query_qualities
        a.flag = r.flag & ~0x900  # never secondary/supplementary
        a.reference_id = r.reference_id
        a.reference_start = ref_start - 1
        a.mapping_quality = r.mapping_quality
        # LEFT keeps the original leading clip (if any) and gains a synthetic trailing clip covering
        # the RIGHT half's query bases; RIGHT is the mirror image.
        lead = [orig_lead] if (is_left and orig_lead) else \
               ([(cut_q + (orig_lead[0] if orig_lead else 0), 'S')] if not is_left else [])
        trail = [orig_trail] if (not is_left and orig_trail) else \
                ([(len(r.query_sequence) - q_end_of_this_side, 'S')] if is_left else [])
        full = lead + sub_ops + trail
        a.cigarstring = ''.join(f'{n}{op}' for n, op in full)
        return a

    lead_len = orig_lead[0] if orig_lead else 0
    left = build(left_ops, lead_len + cut_q, r.reference_start + 1, True)
    right = build(right_ops, None, acceptor + 1, False)
    return left, right
