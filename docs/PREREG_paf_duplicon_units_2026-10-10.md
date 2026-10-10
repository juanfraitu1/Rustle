# Our own duplicon units from the locus self-alignment (PAF units), 2026-10-10

**Pre-registered before any unit is computed.** User (2026-10-10), after "do we have our own way to find duplications?": pre-register duplicon-like units cut from the
alignments we already compute, so the expansion level of `mcl_families --emit-hierarchy` (`docs/PREREG_family_hierarchy_flag_2026-10-09.md`) no longer needs DupMasker,
which exists on disk for human CHM13 v1.0 only.

**Scope and what is not re-proposed.** Units only ANNOTATE: they replace the `--hierarchy-duplicons` track; family membership is untouched (families come from the same PAF
and the same `mcl_families` settings). Letting SD or duplicon structure decide membership failed before (register 676, 677/997, 817, 818) and is not tried. No genome-only
discovery (user 09-13): the alignments are between the loci the family definition already uses. Register 1257-1260 ("structural units") are about read-graph locus models,
not this.

## Method (fixed now; no parameter is fitted to the data)

Input: `loci.paf` (all-vs-all minimap2 of the 2,334 locus sequences, the family definition's own input) and `loci.gff3` (exons). Work in genomic coordinates (a locus
`chrom:s-e` starts at `s-1`), so loci that overlap share their sequence.

1. **Alignments.** Keep a PAF record when (a) the two loci differ and the two aligned intervals do not overlap on the genome (overlapping loci are the same DNA);
   (b) identity = matching bases / block length (columns 10/11) >= 0.90 and block length >= 1,000 bp: **the standard SD definition** (Bailey et al. 2001; the
   SDs DupMasker's duplicons were built from); (c) each aligned interval carries >= 1 exonic base of its locus: **the family definition's own `--min-exonic-bp 1`
   rule**, which also keeps intron-only repeat alignments out.
2. **Breakpoints.** Start from the ends of every kept alignment, on both sides. Project every breakpoint lying inside an alignment through that alignment's CIGAR
   (`cg:Z`) to its partner; add the projection unless a breakpoint already lies within 1,000 bp of it (the SD minimum length: two boundaries closer than an SD
   cannot both bound one). Repeat until no breakpoint is added, at most 20 rounds (the round count and whether the fixpoint was reached are reported; at the cap the
   last state is used).
3. **Atoms.** Cut every genomic stretch covered by a kept alignment at the breakpoints; atoms shorter than 1,000 bp are not units (their bases carry no unit, like
   bases outside every duplicon).
4. **Units.** Join atom a and atom b when some kept alignment maps >= 50% of a onto b and >= 50% of b onto a (reciprocal majority); a unit = a connected
   component (union-find), id `PU<n>`.
5. **Output.** BED `chrom start end unit_id`, read by the existing `--hierarchy-duplicons` loader unchanged; dominant unit per exon, ownership, signature,
   expansions, fusions as in the flag.

Code: `bench/paf_duplicons/paf_units.py` (+ `test_paf_units.py`, test-first: CIGAR projection on both strands, breakpoint snapping, the X/Y/Z chain below,
reciprocal joining, genomic self-overlap skip), `bench/paf_duplicons/score_units.py` (+ tests).
Chain the method must get right (unit test): X-Y aligned over AB, Y-Z over B only ⇒ units {A}, {B}; X and Y have signature {A,B}, Z {B} (a fragment inside the
same expansion), not one unit for everything.

## Test on human CHM13 v1.0 (truth = DupMasker, the registered final run)

Comparator: `/mnt/linuxdisk/tmp/hier_flag/final/ours.hierarchy.members.tsv` (DupMasker colors v1.0). Ours: the same command with `--hierarchy-duplicons <PAF units BED>`.
Per family (>= 2 members), the **cover** = its expansions as member sets plus its `unique` class as one more block (matched only with the other run's `unique`).

- **Family score** = maximum-weight bipartite matching of the two covers' blocks on Jaccard (members), sum of matched Jaccard / max(#blocks of each); 1 = identical.
- **Expansion sensitivity** = share of DupMasker expansions that equal one of ours exactly; **precision** = share of ours that equal one of DupMasker's
  (user reporting rule: sensitivity, precision, bipartite).
- **Baseline B1** = one block holding every member of the family (the "expansion = family" null).
- **Structured families S** = families whose DupMasker cover is not one block holding every member (B1 is wrong there by construction).

**Split.** Family id hashed: `int(sha256(id), 16) % 2 == 0` = DEV, else HELD-OUT. DEV is scored and written first; at most ONE amendment may follow, registered
before the held-out half is computed; the held-out half decides.

- **H0 (guard).** `clusters.tsv`, `loci.tsv`, `pairs.tsv` byte-identical to the DupMasker run (units only annotate).
- **P1 (primary).** Mean family score over held-out families >= **0.80**.
- **P2 (beats the null where it matters).** On held-out S, ours scores above B1 in more families than below, one-sided sign test p < 0.05 (ties dropped).
- **P3 (Soto placement replicates).** Soto families (>= 2 members) equal to exactly one expansion with PAF units >= **194** (90% of the 215 with DupMasker).
  All families; Soto is a comparator, not the objective, so P3 does not decide adoption.
- **P4 (named fusions).** PKD1P6-NPIPP1 and PDXDC2P-NPIPB14P are fusion loci `fusion_of` the NPIP family with PAF units (both needed).

**Verdict.** P1 and P2 hold on held-out ⇒ **ADOPT** PAF units as the default duplicon source of `--emit-hierarchy` when no `--hierarchy-duplicons` track is given
(Rust port, opt-in until reviewed); otherwise NOT ADOPTED, with the failing families listed. P3 and P4 are reported either way.
**Descriptive (reported, no bar):** alignments kept, rounds, atoms, units; share of members' exonic bases carrying a unit vs a DupMasker duplicon; expansions per
family under both; dev and held-out side by side; expansion sensitivity / precision; the NPIP family (MCL9) cover under both.

## Phase 2 (only after the Phase 1 verdict; descriptive, gorilla)

Gorilla has no duplicon track, so Phase 2 has no truth. Inputs (the gorilla locus set and its self-alignment, built like `loci.paf`) are fixed in an amendment
before the run. Named expectation registered now: **LRPAP1** (11 gorilla loci: 8 full copies + 3 5′ fragments, `project_lrpap1_family`) — the 8 full copies share
one maximal signature and the 3 fragments' signatures are subsets of it (one expansion, 3 shared or contained fragments). The NPIP family is reported as is.
