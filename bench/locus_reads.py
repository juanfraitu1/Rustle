"""One-pass, locus-indexed read access for O1/O3 analysis.

WHY THIS EXISTS — it encodes, in code, four criteria this project has got wrong at least once each.
Every one cost a retraction or a wrong headline; the defaults here are the corrected forms.

1. ⚠ A READ BELONGS TO A LOCUS ONLY IF IT LIES MOSTLY INSIDE IT.
   `samtools view <region>` returns everything overlapping by >= 1 bp, which counts reads of an adjacent
   gene that merely clip the boundary. Measured at NPIP (`docs/o1_ledger.md` §5d): 5,544 reads by >=1 bp
   overlap against 419 inside -- only 7.6% survive -- and one locus went from 4,784 to 19. That error put
   "28/31 expressed" and "5,544 reads" into a published artifact. DEFAULT: `min_inside=0.50`.

2. ⚠ `N` IS AN INTRON, `D` IS NOT. Splitting on both silently merges exons.

3. ⚠ THE READ GATE POOLS OVER THE LOCUS, NOT THE CANDIDATE. The builder sums support across the
   connected component of the junction-incidence graph (`pool_locus_support: true`,
   `denovo_assemble.rs`), so a 2-read chain survives when its locus totals >= 3. Counting per-candidate
   under-reports; counting "reads in the interval" over-reports by ~5x. `pooled_support()` is the
   builder's rule.

4. ⚠ PRIMARY ONLY (`-F 2308`) for any per-read CIGAR statistic.

SPEED, secondarily: one indexed pass instead of a `samtools view` subprocess per locus. Measured on the
3-contig NPIP BAM, the per-locus form made 31 process spawns and timed out at 120 s; the single pass
takes ~44 s, and with pysam's own index it is faster still. The catalog run it feeds is 16-30 min, so
this is not the bottleneck -- but it is free.
"""
from __future__ import annotations
import collections
from typing import Iterable, Sequence

try:
    import pysam
except ImportError:  # pragma: no cover - environment without pysam
    pysam = None

FLAG_EXCLUDE = 0x100 | 0x200 | 0x800 | 0x4  # secondary|qcfail|supplementary|unmapped  == -F 2308


class Read:
    __slots__ = ("start", "end", "introns", "reverse", "mapq", "clip_frac", "name")

    def __init__(self, start, end, introns, reverse, mapq, clip_frac, name):
        self.start, self.end, self.introns = start, end, introns
        self.reverse, self.mapq, self.clip_frac, self.name = reverse, mapq, clip_frac, name

    @property
    def exon_blocks(self):
        """Aligned blocks with introns removed — the read's EXONIC footprint."""
        out, cur = [], self.start
        for a, b in self.introns:
            if a > cur:
                out.append((cur, a))
            cur = b
        if self.end > cur:
            out.append((cur, self.end))
        return out


def _parse(aln):
    ref = aln.reference_start
    introns, lead, trail, qlen = [], 0, 0, 0
    ops = aln.cigartuples or []
    for i, (op, n) in enumerate(ops):
        if op in (0, 7, 8):      # M = X
            ref += n; qlen += n
        elif op == 2:            # D  -- NOT an intron
            ref += n
        elif op == 3:            # N  -- intron
            introns.append((ref, ref + n)); ref += n
        elif op == 1:            # I
            qlen += n
        elif op in (4, 5):       # S H
            qlen += n
            if i == 0: lead = n
            if i == len(ops) - 1: trail = n
    return ref, introns, (lead + trail) / qlen if qlen else 0.0


class LocusReads:
    """Reads grouped by locus in ONE pass. `loci` is an iterable of (chrom, start, end), 0-based."""

    def __init__(self, bam: str, loci: Iterable[Sequence], min_inside: float = 0.50):
        if pysam is None:
            raise RuntimeError("pysam is required; install it or use the streamed fallback")
        self.loci = [(str(c), int(s), int(e)) for c, s, e in loci]
        self.min_inside = min_inside
        by_chrom = collections.defaultdict(list)
        for i, (c, s, e) in enumerate(self.loci):
            by_chrom[c].append((s, e, i))
        self.reads = collections.defaultdict(list)
        with pysam.AlignmentFile(bam, "rb") as fh:
            for c, spans in by_chrom.items():
                lo, hi = min(s for s, _, _ in spans), max(e for _, e, _ in spans)
                for aln in fh.fetch(c, lo, hi):
                    if aln.flag & FLAG_EXCLUDE:
                        continue
                    end, introns, clip = _parse(aln)
                    a = aln.reference_start
                    if end <= a:
                        continue
                    for s, e, i in spans:
                        # ⚠ the read must lie MOSTLY INSIDE — not merely overlap (criterion 1)
                        if (min(end, e) - max(a, s)) >= min_inside * (end - a):
                            self.reads[i].append(
                                Read(a, end, introns, aln.is_reverse, aln.mapping_quality, clip,
                                     aln.query_name))

    def at(self, i: int):
        return self.reads.get(i, [])

    def chains(self, i: int) -> collections.Counter:
        """Exact intron chains and their read counts (criterion 2)."""
        return collections.Counter(tuple(r.introns) for r in self.at(i) if r.introns)

    def pooled_support(self, i: int) -> int:
        """The BUILDER's gate: support summed over the junction-incidence component (criterion 3)."""
        ch = self.chains(i)
        if not ch:
            return 0
        junc = collections.defaultdict(set)
        for k in ch:
            for j in k:
                junc[j].add(k)
        best, seen = 0, set()
        for k in ch:
            if k in seen:
                continue
            comp, stack = {k}, [k]
            while stack:
                for j in stack.pop():
                    for o in junc[j]:
                        if o not in comp:
                            comp.add(o); stack.append(o)
            seen |= comp
            best = max(best, sum(ch[x] for x in comp))
        return best

    def expressed(self, i: int, min_reads: int = 3) -> bool:
        return len(self.at(i)) >= min_reads


# ---------------------------------------------------------------------------------------------------
# SELF-CHECK. Run: python3 bench/locus_reads.py
#
# Validates against numbers established INDEPENDENTLY of this file (docs/o1_ledger.md §5d and §5f), not
# against itself. A helper whose whole purpose is to encode corrected criteria has to be checkable, or it
# becomes another place for the same errors to hide.
if __name__ == "__main__":
    import json, sys, time

    BAM = "/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam"
    LOCI = "/mnt/linuxdisk/home/juanfraitu/npip/ggo_loci.json"
    try:
        loci = json.load(open(LOCI))
    except OSError:
        print(f"skip: {LOCI} not present on this machine")
        sys.exit(0)

    t0 = time.time()
    lr = LocusReads(BAM, loci)
    elapsed = time.time() - t0

    total = sum(len(lr.at(i)) for i in range(len(loci)))
    expressed = sum(1 for i in range(len(loci)) if lr.expressed(i, 3))
    # §5f: the four loci that clear the pooled read gate yet still build no node
    four = {("NC_073242.2", 29415572): 6, ("NC_073242.2", 99221059): 4,
            ("NC_073242.2", 28377934): 3, ("NC_073242.2", 31431689): 4}
    idx = {(str(c), int(s)): i for i, (c, s, _) in enumerate(loci)}

    checks = [("§5d total reads INSIDE the 31 loci", total, 419),
              ("§5d loci expressed (>=3 reads inside)", expressed, 23)]
    for (c, s), want in four.items():
        i = idx.get((c, s))
        if i is not None:
            checks.append((f"§5f pooled support at {c}:{s}", lr.pooled_support(i), want))

    bad = 0
    for name, got, want in checks:
        if got != want:
            bad += 1
        print(f"[{'ok ' if got == want else 'FAIL'}] {name}: got {got}, expected {want}")
    print(f"\none pass over {len(loci)} loci in {elapsed:.1f}s; {len(checks) - bad}/{len(checks)} checks pass")
    sys.exit(1 if bad else 0)
