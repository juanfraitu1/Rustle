#!/usr/bin/env python3
"""Counting reads at a locus — the ONE correct way, and the wrong way named so it cannot be reached
by accident.

⛔⛔ **THE ERROR THIS EXISTS TO PREVENT (ledger §6cm, 2026-09-02).** Counting reads that *overlap a
locus span* instead of reads that have an *aligned block inside it*. On 2026-09-02 that inflated a
headline 3.4x and produced a whole retracted mechanism: `NPIPP1` appeared to have 1,608 reads
collapsing into a 4-read locus, when 1,151 of them — **71.6%** — merely spliced across it with no
aligned base, because `PDXDC1` (168 kb) physically CONTAINS `NPIPP1` and its transcript passes
straight through. The real figure is 457, itself an upper bound.

The project rule is old and was not enforced anywhere: **`N` in an RNA CIGAR is an intron, spliced
OUT; a read spliced OVER a locus is no evidence for it.** `samtools view -c REGION` counts the wrong
thing, silently, and reads beautifully in a script.

Use `reads_with_block_in()`. `reads_overlapping_span()` exists only so the difference is visible and
so a caller who genuinely wants span overlap has to say so.

⚠ Primary alignments only (`-F 2308`) throughout — the standing invariant before any per-read CIGAR
statistic, so one molecule is never two witnesses.
"""
import re
import subprocess

SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_blocks(pos0, cigar):
    """Reference blocks a read actually aligns to. `N` breaks the run; `M/=/X/D` extend it."""
    out, p, cur = [], pos0, None
    for ln, op in _CIG.findall(cigar):
        ln = int(ln)
        if op in "M=XD":
            cur = (cur[0] if cur else p, p + ln)
            p += ln
        elif op == "N":
            if cur:
                out.append(cur)
                cur = None
            p += ln
    if cur:
        out.append(cur)
    return out


def _fetch(bam, chrom, start, end):
    r = subprocess.run(
        [SAMTOOLS, "view", "-F", "2308", bam, f"{chrom}:{start + 1}-{end}"],
        capture_output=True, text=True,
    )
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) > 5 and f[5] != "*":
            yield f


def reads_with_block_in(bam, chrom, start, end, blocks=None):
    """⭐ THE CORRECT COUNT: primaries with >= 1 ALIGNED BASE inside [start, end).

    `blocks` optionally restricts to a locus's own exon blocks — stricter still, and what you want
    when a large gene's exons fall inside the span of a small one nested within it.
    """
    n = 0
    for f in _fetch(bam, chrom, start, end):
        bl = aligned_blocks(int(f[3]) - 1, f[5])
        if blocks is None:
            if any(bs < end and start < be for bs, be in bl):
                n += 1
        elif any(bs < e2 and s2 < be for bs, be in bl for s2, e2 in blocks):
            n += 1
    return n


def reads_overlapping_span(bam, chrom, start, end):
    """⚠ THE MISLEADING COUNT — what `samtools view -c` gives. Includes reads that splice straight
    over the locus contributing no aligned base. Returns `(n_overlapping, n_spliced_over)` so the
    caller cannot quote the first without seeing the second."""
    tot = over = 0
    for f in _fetch(bam, chrom, start, end):
        tot += 1
        bl = aligned_blocks(int(f[3]) - 1, f[5])
        if not any(bs < end and start < be for bs, be in bl):
            over += 1
    return tot, over


def spanning_genes(gff_gz, chrom, start, end, min_cover=0.60):
    """⭐ ASK THIS BEFORE BLAMING THE PIPELINE for an oversized locus (§6cm).

    Returns annotated genes covering >= `min_cover` of the node. On 2026-09-02 two nodes called
    "mis-chained giants" turned out to be `SNX29` (covers 100.0%) and `PDXDC1` (99.9%) — real genes,
    correctly assembled, with canonical junctions and majority read support because they are real
    transcripts. A rule fitted to that "blob class" would have been fitted to two real genes.
    """
    import gzip
    out = []
    with gzip.open(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "pseudogene"):
                continue
            gs, ge = int(f[3]) - 1, int(f[4])
            if gs < end and start < ge:
                cov = (min(ge, end) - max(gs, start)) / max(1, end - start)
                if cov >= min_cover:
                    m = re.search(r"Name=([^;\n]+)", f[8])
                    out.append((m.group(1) if m else "?", gs, ge, cov))
    return sorted(out, key=lambda x: -x[3])


if __name__ == "__main__":
    import sys
    bam, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    tot, over = reads_overlapping_span(bam, chrom, start, end)
    good = reads_with_block_in(bam, chrom, start, end)
    print(f"{chrom}:{start}-{end}")
    print(f"  aligned block inside (USE THIS) : {good}")
    print(f"  overlapping the span            : {tot}")
    print(f"  spliced OVER, no aligned base   : {over} = {over / max(1, tot):.1%}")
