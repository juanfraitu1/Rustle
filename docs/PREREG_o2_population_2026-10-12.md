# How many multi-placement reads can O2 assign, and how many does it move off the aligner's primary? 2026-10-12

**Pre-registered before the script runs.** User (2026-10-12): "how common are assignable reads to copies, do we have a lot of cases in
which it is impossible or do we have enough cases of a decent number of reads correctly assigned that StringTie and FLAIR cannot do
since they don't use secondary alignments".

**Population** (per species, benchmark BAMs `GGO_mm.bam` / `A119b.t2t.bam`, all chromosomes): molecules with another alignment within
the good window — `as_table`'s second-best AS >= 0.98 x best AS. **Sample**: 20,000 molecules per species, uniform, seed 1
(`random.Random(1).sample` over the population in table order). **Per molecule**: the primary record vs every non-supplementary
record with AS >= 0.98 x best (O2's certificate exactly as `bench/tool_difference/s1_o2.py`). Classes, first that applies:
**moved** (some rival certified over the primary: the aligner's — and so StringTie's / FLAIR's — placement is wrong and O2 can
correct it), **confirmed** (the primary certified over every rival), **identical** (no distinguishing column against any rival:
no read-level method can assign it), **undecided** (otherwise). Reported with 95% Wilson intervals and extrapolated to the
population; also split by the primary's AS rank (true tie primary = best = second; primary below the best; primary best with a close
rival), the mosaic share, and, for moved molecules, whether the winning placement is spliced where the primary is not (the
retrocopy direction). Descriptive; nothing adopted. Expected (written now): identical and confirmed dominate; moved is a small share.
