# Pre-registration: the SMN1 / SMN2 showcase — PSV-decided reads, assign-or-abstain, and an assembler that honours them (2026-10-20)

**User (2026-10-20):** "O2 ... is to be able to offer a way to find the true provenance of multi-mapping long reads where the
flag as primary or secondary is just a coin toss and meaningless ... among the [families] there must be a family in which
we can tell based on PSVs, tags ... and use that to improve the assemblies that just ignore secondary alignments." Asked
which tried families' members the PSV columns tell apart, and offered SMN1/SMN2 as the showcase: "yes lets do that".

## The family and the reads

**SMN1 / SMN2** on CHM13 (pool family `S19SMN2`: `.0` = SMN2 chr5:70,808,816-70,837,687 −; `.1` = SMN1 chr5:71,381,817-
71,409,805 +), 99.8% identical. **HG002 Kinnex** reads (`hg002proof/chm13.bam`), the one real donor whose own assembly gives
every read's copy (`hg002truth`: classes unique / unique_elsewhere). SMN2 skips exon 7 in most transcripts; SMN1 keeps it.

**Seen before, stated plainly.** On these reads: register 1344 (the corrector's picks, StringTie/FLAIR, exon-7 readout) and
register 1348 (the set rule's read-level accuracy on SMN: 2,885 tied reads with a copy in the window, 2,258 decided, 2,249
right; 627 ambiguous, on which minimap2 is right 28). This is a demonstration with pre-registered bars on a matched
donor, not a held-out test. What has not been seen: what any assembler builds from the set rule's homes.

## O2's rule for this showcase (code at this file's commit)

For every read whose minimap2 primary lies on SMN1 or SMN2:
- **tied** when a rival copy's alignment is within 30 points of its best intron-free score (`ms`; `o2window.py`,
  window 30, stitched alignments excluded);
- a tied read is **decided** when, at every column where its two alignments differ by a substitution, the read agrees
  with one copy and not the other (`pareto_set` leaves one copy); it goes to that copy;
- a tied read that is not decided is **abstained**: identical over its bases on both copies, or mixed; it is given to
  neither copy;
- a read that is not tied keeps minimap2's primary.

No thresholds beyond the window; no use of the assembly.

## Arms (one record per read, `o2loci.py arms`)

| arm | each read's primary |
|---|---|
| **mm2** | minimap2's primary (what StringTie, FLAIR and a primaries-only assembler see) |
| **O2** | decided → its copy; abstained → **left out of the family**; untied → minimap2's |
| **O2-coin** | as O2, but abstained → minimap2's primary (reported; shows what the coin toss costs) |
| **truth** | HG002's assembly copy (unique / unique_elsewhere), else minimap2's |

**Assemblers on each arm's BAM:** Rustle (`tools/rustle_pipeline.sh assemble --seed-pool primary`, the thesis's own:
it builds loci from the records it is given), StringTie 3.0.3 `-L`, FLAIR 3.0.0 `collapse --trust_ends`.

## Readout (register 1344's, unchanged): exon-7 inclusion per copy

A read or transcript on a copy **spans** when its blocks overlap exon 6 and exon 8; it is FL when a block overlaps exon 7
by >= 20 bp, else Δ7. Inclusion = FL / (FL + Δ7): reads counted one each, StringTie transcripts by `cov`, FLAIR isoforms by
mapped reads, Rustle transcripts by `reads`. Exon coordinates: `o2loci_smn/smn_spec.json`.

**The reference values** are the truth arm's reads: SMN2 0.617, SMN1 0.985 (register 1344).

## Bars

- **R1 (reads):** on each copy, the O2 arm's read inclusion is closer to the truth arm's than mm2's is.
- **R2 (loci, O2 + Rustle against the tools that ignore secondaries):** on each copy, |Rustle on O2 − truth reads| is
  smaller than each of |StringTie on mm2 − truth reads|, |FLAIR on mm2 − truth reads| and |Rustle on mm2 − truth reads|.
- **R3 (the assembler honours the homes):** on each copy, |Rustle on O2 − Rustle on truth| <= 0.03.
- **R4 (abstaining is not free of cost but is honest):** per copy, the O2 arm's spanning-read count is at least 60% of
  the truth arm's (the abstained reads are the price; report their number).

**Verdict:** PASS = R1-R4 hold; FAIL otherwise. No power bar: SMN has thousands of reads.

**A PASS would claim:** on a 99.8%-identical family in a real donor, reading each tied read's PSV columns gives most reads
their copy and leaves the undecidable ones out, and an assembler that keeps those homes recovers each copy's own exon-7
pattern where StringTie, FLAIR and a primaries-only assembly mix them. It would not claim other families or donors.

## Rules

- Code written and tested before any arm is assembled; bugs found while running go into an amendment with the rerun.
- Foreground, one assembler at a time. Human numbers only.
- The metric-auditor checks the result on a frozen copy before it is reported.

## Amendments

(none yet)
