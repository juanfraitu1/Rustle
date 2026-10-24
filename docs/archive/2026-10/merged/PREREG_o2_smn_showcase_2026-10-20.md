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

## Results (code eca1c331; `o2smn/score.log`, `bars.json`; run once): **FAIL (R2 on SMN2, R4 on SMN2)**

3,651 reads with a primary on SMN1/SMN2; 3,623 tied within 30 points; the set decides 2,861 and abstains on 762. O2 moves
8 reads against minimap2; the truth moves 616.

| exon-7 inclusion | mm2 | O2 | O2-coin | truth |
|---|---|---|---|---|
| reads SMN2 | 0.653 | 0.635 | 0.646 | **0.617** |
| reads SMN1 | 0.887 | 0.988 | 0.890 | **0.985** |
| Rustle SMN2 | 0.602 | 0.590 | 0.602 | 0.585 |
| Rustle SMN1 | 0.801 | 0.969 | 0.801 | 0.941 |
| StringTie SMN2 | 0.654 | 0.634 | 0.649 | 0.618 |
| StringTie SMN1 | 0.894 | 1.000 | 0.896 | 0.989 |
| FLAIR SMN2 | 0.763 | 0.710 | 0.763 | 0.700 |
| FLAIR SMN1 | 0.897 | 0.929 | 0.897 | 0.934 |

- **R1 holds** on both copies (O2 0.018 / 0.003 from the truth reads; minimap2 0.036 / 0.098).
- **R2 holds on SMN1** (Rustle on O2 0.016 from truth; StringTie, FLAIR and Rustle on primaries 0.091 / 0.088 / 0.184)
  and **against StringTie and FLAIR on SMN2** (0.027 vs 0.037 / 0.146), but **fails against Rustle on primaries on SMN2**
  (0.015): minimap2's SMN2 reads are mostly right, the coin toss damages SMN1.
- **R3 holds** (0.005 / 0.028): Rustle keeps the homes it is given; FLAIR does not (O2-coin = mm2 for it).
- **R4 fails on SMN2:** abstaining leaves 723 of its 1,355 spanning reads (53%); SMN1 keeps 96%.
- O2-coin (abstained reads at the primary) is minimap2's picture again: the abstained reads are where minimap2 errs.

**Reading.** The PSV columns give SMN1 its own isoforms (exon 7 kept in 98.8% of reads, as in truth, against 88.7% under
minimap2), and an assembler that keeps the homes reproduces that; StringTie and FLAIR on primaries do not. The price is
abstention: nearly half of SMN2's reads are identical to SMN1 over their span and are left out. Audit pending.

## Audit (metric-auditor, 2026-10-20): holds with caveats

- Every number reproduces. The FAIL stands.
- **The SMN1 gain is removal, not assignment:** of the 616 reads the truth moves, O2 fixes 8 and abstains on 599; the
  other 2,241 decisions confirm a primary minimap2 already had right.
- **349 of the abstained FL reads carry SMN2's c.840 T on an SMN1 primary** while matching CHM13's SMN1 at other columns:
  HG002's SMN2 differs from CHM13's (the donor-vs-reference case), and the every-column rule abstains. "Identical over
  their span" was wrong for them.
- **R2 on SMN2 is infeasible:** Rustle on the truth arm is itself 0.033 from the truth reads; Rustle credits `reads` to
  16% of the family's reads, and R3 on SMN1 turns on one read (8 vs 7). The FLAIR "O2" column re-admitted the abstained
  reads from reads.fq and is not FLAIR on the O2 arm. R4's inputs pre-dated the prereg (r1348's 627 ambiguous).
- Wording not to use: "the PSV columns give most reads their copy"; "O2 recovers SMN1's isoforms"; "FLAIR on O2".

## Tags (descriptive, `o2tags.py`; `o2smn/tags.log`, `o2smn/tags/*.log`)

Per tied read, which statistic of its best alignment on each copy picks the donor's copy:

| statistic | SMN: decides / right / wrong (2,885 tied) | on the 627 abstained |
|---|---|---|
| AS, ms | 2,863 / 2,257 / 606 | 606 / 17 / 589 |
| NM | 2,509 / 2,252 / 257 | 262 / 14 / 248 |
| de | 2,511 / 2,256 / 255 | 263 / 16 / 247 |
| **substitutions only (NM − CIGAR indel bases = MD mismatches)** | **2,604 / 2,603 / 1** | **354 / 354 / 0** |
| indel bases | 2,683 / 1,476 / 1,207 | 603 / 11 / 592 |

Substitutions alone also on THOC3 5,482 / 5; DUSP22 1,927 / 1; GTF2H2C 1,408 / 3; NOMO1 348 / 0; GOLGA8A 257 / 10
(minimap2 256 / 282). Where the donor's copies differ from CHM13 every statistic fails: TBC1D3B 0 / 75 (AS 0 / 147), PDPR
179 / 98, GTF2IRD2 157 / 37.
