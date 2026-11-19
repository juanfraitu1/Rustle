# Goal: resolution-aware assembly of multi-copy gene families (set by the user 2026-11-19)

## The goal
For every multi-copy family, decide from the reads alone whether its copies can be separated, and assemble at that
resolution:
- one transcript per copy where the reads tell the copies apart;
- one transcript per copy set where they cannot (near-identical paralogs, tandem arrays);
- in neither case invent expression at copies the reads do not support.

**One measure.** Transcripts scored against the origin of their reads (own-haplotype truth: HG002 on CHM13, KB3781 on its
own reference; realistic simulation for the regimes real data lacks). A transcript is right when the copies it claims
equal the copies its reads fit equally
(on HG002 / KB3781: the copies of the individual's own haplotypes a read matches equally well, the truth's origin set; in a
simulation the literal source copy is a stricter truth, under which no read-based method is right on identical copies:
development ledger, step 6). Precision and recall per regime, StringTie at primaries-only and at the defaults as
the trusted baselines.

## Four questions
1. **Near ties: what can be decided, and what should be reported?** Settled: near ties are about 6% of family reads, the
   primary is right on about 92% of the decidable ones, no placement rule beats it (rows 1642-1651, 1777-1785); the
   copy-set statement is accurate (row 1796 B1: 0.82 vs 0.11). Open: a per-family call "separable / copy set / not
   measurable" read off the output, tested against origin truth.
2. **Tandem arrays: what is the right unit?** Row 1820: exact ties spread structure across near-identical units (TSPY
   11 -> 25 of 38) at a cost in false chains; some arrays get no structure (USP17L, REXO1L). Open: an array reported as
   one transcript set over N units, the reads' spread over units as the quantity, tested on testis samples; needs a truth
   not capped at 12 loci (simulation, or the individual's own assembly).
3. **When should secondary alignments seed loci and isoforms?** Settled: wider seeding bands add chains mostly at shadow
   loci (row 1767); undecidable records add almost nothing over primaries (row 1820). Open: secondaries as junction
   support only, never as seeds; the exact-tie band inside arrays only. Judged on origin truth, new and lost chains apart.
4. **What can be borrowed across paralogs?** Settled: pooling reads and lifting structure across copies add no correct
   chains (rows 1519, 1794, 1798, 1809); perfect provenance is worth x1.0-1.25 at most. Open: borrowing as evidence about
   structure, not reads: a junction seen in one copy supports the same junction column in a copy with few reads (the
   `exonmap` axis, row 1595), reported as a labelled class, never as expression.

## Frame and limits
Every result feeds the variation-graph model: a copy set is a set of paths through the family graph. Rules work on any
species with a reference genome, need no second sample at analysis time, finished and generalised by early January 2027.

## Data budget (read before any row)
- HG002 on CHM13, hash quarters of the selected reads: residues 0-24 seen (rows 1776-1791, development), 25-49 spent
  (row 1796), 50-74 spent (row 1820), **75-99 the last unseen quarter**: reserve it for one final judging row of the rules
  this goal produces, pre-registered together. Development uses 0-24 and KB3781.
- KB3781: own reference, own haplotypes; few undecidable reads (186 exact ms ties of 100,810 molecules): mechanics only.
- Arrays (question 2): HG002 is a cell line (14 reads on the TSPY array); testis samples (A119b, OR6737) have no origin
  truth. A simulation is the only judge there.

## Plan (order)
1. The measure: `bench/famwise_eval/copyset_score.py` — claims of any arm (per-copy GTF transcripts, `famwise` family
   transcripts) against read origin: precision (claim's copy set = its supporters' origin) and recall (true
   (copy set, chain) units with >= 2 supporting reads recovered). Development numbers on the seen quarter for Rustle P / D,
   the `famwise` layer, StringTie P / D.
2. Question 1: a per-family resolution call from the reads alone, checked on the seen quarter.
3. Questions 3 and 4: junction-support-only secondaries; borrowed junction columns as a labelled class.
4. Question 2: array simulation (identical and near-identical units, testis-like spread).
5. One pre-registration for the rules of 2-4, judged once on residues 75-99.

## Status (2026-11-19, end of the first PSU session)
Development ledger: workspace `goal_multicopy_dev/LEDGER.md` (steps 1-11). Judged: row 1821.
1. Near ties / the call: `resolution_call.py`. Row 1821 on the last unseen HG002 quarter: R1 fails narrowly (+1.0% over the
   uncalled famwise layer), R2 passes (precision 0.888 -> 0.909, resolved units +2.5% over the per-copy output), R3 passes.
   Grouping copies is what helps; the call adds little over the layer.
2. Arrays: simulation only (identical units: one copy set over N units is the least wrong output, the spread over units is
   not recoverable from placements; >= 0.2% divergence: primaries place every read). OPEN: no origin truth exists for real
   testis arrays (needs a testis Iso-Seq sample of an individual with assembled haplotypes).
3. Seeding: `RUSTLE_SEED_SPLICED=support` (opt-in, built): on the judging quarter it keeps about half the right new chains of
   the good pool and cuts wrong-copy new chains 4x (60 -> 15); removes silent-unit chains in arrays at 0.1-0.2%.
4. Borrowing: `borrow.py`, labelled `borrowed_structure`: structure at 5 more copy-set units, 0.17 right when read as
   expression (development only, not judged; HG002 has no unseen quarter left).
Family-member completeness (user's question): missing members are mostly silent; none lost to wrong placement; tied
seeding restores 3-7; the actionable gap is the assembler's low-support floor (StringTie finds 24-75 more expressed
members; RUSTLE_READ_ISOFORM_K=2 closes 9 of 22 on development).
