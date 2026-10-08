# Pre-registration — O3: what happens to the reads of copies the maternal haplotype lacks, and can a truth-free chain recover them?

**Written 2026-10-08, before any read is mapped or any number is looked at for this question.** DRAFT for user review; nothing below has been run.
Origin: advisor interest in the LRPAP1 page (UvFQf5QX): maternal and paternal haplotypes of one animal carry different copy numbers. Real-data version of
**Arm B** of `docs/archive/2026-09/PREREG_o3_reference_bias_2026-09-23.md` (its Arm A, simulation, refuted "missing copies hide as unmapped/tied reads" to 5% divergence).
Reuses, unchanged, the chain and bars of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` (Amendments 8, 9, 10). **No new thresholds.**

## 1. Questions
- **Q1 (fate).** Reads from a copy present on the paternal haplotype and absent from the maternal one: aligned to the maternal assembly, are they UNMAPPED, TIED
  (multi-mappers), or ABSORBED as untied primaries on another copy? Per locus, never pooled.
- **Q2 (O3 recovery).** Starting from the maternal alignment only (truth-free), does the net -> cluster -> flag chain recover those copies? Two clusterers,
  same inputs, same bars: IsoCon 0.3.3 and the in-house `o3_candidates` stage.

## 2. Substrate
KB3781 ("Jim", male), fibroblast Iso-Seq = the assembly's own cell line (`project_o3_matched_individual`). Reference = maternal assembly `mat`
(`GCA_028885495.2`, 225 sequences, `winloci_data/mGorGor1.mat.splice.mmi`); truth = paternal `pat` (`GCA_028885475.2`, 24 sequences). Testis (OR6737) is a
different animal: **excluded** (neither haplotype is its truth). Gorilla only; no human numbers.

## 3. Truth: MAT-ABSENT loci
A locus is MAT-ABSENT iff it exists on `pat` and has no counterpart on `mat`:
1. **Catalog loci**: `refabsent/bonly.tsv` rows with `hap = pat` (haplotype-only loci on the 9 chromosomes `_pri` took from `mat`; asm20, identity >= 0.90,
   coverage >= 0.80, overlapping no lifted copy of the family) — 35 loci / 34-family file `fams_bonly.txt`. Expressed = >= 3 reads best-placed there
   (`bonly_expressed.tsv`); LARGE = >= 20 reads. Known expressed: GWFAM175_B0 (281), GWFAM205_B0 (77), GWFAM227_B0 (13), GWFAM175_B1 (6).
2. **LRPAP1** (not in the catalog): the 11 loci of `lrpap1.copies.tsv` + `partial.copies.tsv`; mat-absent = those with no mat hit in `copies8.mat.paf` /
   `partial3.mat.paf`, read from the file and not from the artifact flags: **p12** (chr12 23.07 Mb, LOC134756368, 83 fibroblast reads) and the **chrY** copy.
3. **chrY copy = SEX CONTROL, not a copy-number locus** (male: Y is paternal only; same confound as the X-linked MAGE probes, `project_haplotype_cnv_proven`).
   Reported separately; excluded from every bar.
4. **Reverse direction** (held out; reference = `pat`, truth = `mat`): `bonly.tsv` rows with `hap = mat` (92 loci; expressed: GWFAM70_B4 24, GWFAM390_B1 20,
   GWFAM70_B1 16, GWFAM26_B2 13, GWFAM26_B3 9, ...). Primary direction is `mat` (the user's question); `pat` is the substrate held back.

SHARED-COPY CONTROL: reads whose best `pat` record falls on a pat locus that has a `mat` ortholog (lift in `copies_lift.tsv`, classes T2d/T2i).

## 4. Reads (selection bias disclosed)
- **R34**: `refabsent/scored.fa`, 32,219 net reads of the 34 families, already aligned to `mat` and `pat` with the baseline `@PG` command
  (`reads.{mat,pat}.bam`). Selection = a baseline (`_pri`) record on a copy of the family, so **every read was mapped on `_pri` before this test**.
- **R_LRP**: reads of the 11 LRPAP1 loci by the same rule (any primary/secondary record overlapping a locus body, <= 2,000 per locus, seed 1); aligned to
  `mat` and `pat` with the baseline `@PG` command (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, minimap2 2.30 vs 2.31, disclosed).
- **R_unm**: the unmapped primaries of the fibroblast BAM (959 reads, median 69 bp), aligned to `mat`/`pat`; descriptive only.
- Consequence: UNMAPPED on `mat` is observable only for reads that mapped on `_pri`; reads unmapped on `_pri` are covered by R_unm alone.

## 5. Read truth label and fate (Q1)
A read belongs to locus L iff its best `pat` record (alignment score, untied at 0.98 — the `express` rule of `refabsent_truth.py`) overlaps L.
Fate on `mat`, mutually exclusive, in this order (every constant is already registered):
1. **UNMAPPED**: no record, or primary query coverage < 0.80 (the registered hit-coverage rule; < 0.80 reported as PARTIAL, summed into this class).
2. **TIED**: primary MAPQ 0, or a secondary within 0.98 of the primary's alignment score (the registered tie rule).
3. **ABSORBED**: untied primary elsewhere; split into ON-NEAREST-PARALOG (the mat locus with the highest identity to L) and OTHER. Median `de` reported.
Per locus report: n reads, fate counts, nearest-mat-paralog identity (asm20 of L's sequence against `mat`), `de` of absorbed reads vs `de` of shared-copy reads.

**Bar (verbatim from the 09-23 registration, applied per LARGE locus, `pat`-copy reads on `mat`):**
>= 80% ABSORBED with `de` ~ divergence, <= 10% UNMAPPED, <= 10% TIED => advisor's claim **REFUTED at that locus**; 20-50% UNMAPPED+TIED => **PARTLY**;
>= 50% UNMAPPED or TIED => **SUPPORTED** at that locus.
**Prediction, written before looking:** REFUTED at every LARGE locus; GWFAM175_B0 ABSORBED at `de` ~ 0.066 (the `_pri` arm of Amendment 10 is this very
configuration: chr5 is mat-sourced); p12 ABSORBED on the chr14 ortholog of p14 (mRNA differs by 4 bp in 2,446; mat best hit NM 66 / 15.1 kb at MAPQ 60);
chrY copy ~ same fate as an autosomal absent copy of its nearest paralog. UNMAPPED > 10% is expected only where the nearest mat relative is far (see context panel).
**Context panel (no decision):** the 2026-08-14 whole-genome excision panel (`project_o3_excision_wholegenome`: ABSORBED 64.2%, ORPHANED 33.3%) re-plotted as
fate vs nearest-surviving-relative identity, so the natural cases (nearest relative >= 0.91) sit on the right of a spectrum whose left half is synthetic. Reuse of
existing tables only; if a per-family table cannot be recovered the panel is omitted, not recomputed.

## 6. Recovery chain (Q2) — reference = `mat`, truth-free until scoring
Net per family = reads with a record (primary or secondary) on any `mat` copy locus of the family (loci = `copies.mat.paf` hits at identity >= 0.90, coverage
>= 0.80; for LRPAP1 the 11 loci's `mat` hits), plus unmapped reads attributed as in A13; <= 1,000 reads per family, seed 1 (Amendments 9/10).
- **Arm I (IsoCon):** `IsoCon pipeline --nr_cores 4`, defaults, env `isocon` (as 10-01). Outputs -> flag (identity x coverage < 0.999 vs `mat`) -> link to a `mat`
  locus at whole-length d <= delta = 0.00958 (allele) -> merge (Amendment 8) -> candidates. **FLAG iff >= 2 transcripts** (floor registered in Amendment 9/10).
- **Arm H (in-house):** `o3_candidates` as shipped at main b29afa55, **opt-in, Amendment 15 fix NOT applied** (the A14 consensus defect is known and
  disclosed; this arm will probably reproduce its false flags — that is data, not a surprise). Same net, same link/merge/flag rules, same scorer.
- **Scoring (truth used only here):** a candidate RECOVERS locus L iff its best `pat` hit (identity x coverage >= 0.999) overlaps L (the Amendment 10 scorer,
  `refabsent_score.py`, with `mat` as baseline instead of `_pri`).

**Decisions (reused bars, applied per arm, mat-absent loci of Section 3):**
- **R1 (= D1)** every expressed BEYOND-delta LARGE locus is recovered by a FLAG. Smaller expressed loci reported, not judged.
- **R2 (= D2)** no expressed WITHIN-delta locus is flagged as a new copy — the delta rule's own check.
- **R3 (= D3)** families with no mat-absent expressed locus: those with a >= 2-transcript class b/c candidate <= 20%.
- **R4 (reported, not judged):** candidate-level sensitivity, precision and bipartite matching (`feedback_report_sens_prec_bipartite`); IsoCon outputs per
  recovered transcript; read-level O2 effect (reads of a recovered locus: primary on `mat` paralog vs on the candidate, `de` before/after, Amendment 10 O2 side).
- Neither arm is declared the winner; each bar is reported per arm.

**p12 is pre-registered as a DESIGNED-UNDETECTABLE case, not a failure.** Its mRNA is 4 bp (0.16%) from p14, below delta (0.96%), so the registered link rule
will classify it as an allele of the `mat` p14 locus. It is reported on its own line with three readouts: (a) do IsoCon / in-house separate p12 reads from p14
reads at the sequence level; (b) is the p12 transcript linked to p14 by delta; (c) exploratory **haplotype count**: distinct sequence types at the p12/p14
locus vs the 2 alleles one locus can carry (KB3781 has p14 on both haplotypes, so p12 + pat-p14 + mat-p14 = 3 types at one apparent locus). This rule is
**not** part of any decision here; if it looks useful it gets its own registration. Cluster-merging risk is stated up front: IsoCon merged TBC1D3 paralogs 6
edits apart in simulation (`RNA_ALLELE_ISOCON_2026-10-01.md`).

## 7. Outcomes that would change what we claim
- Q1: any LARGE locus with UNMAPPED+TIED >= 50% => the advisor's mechanism holds there; the O3 route would need an unmapped/tied-read recovery step (Arm B/C of 09-23).
- Q2: R1 fails for GWFAM175_B0 (recovered once before, 4 transcripts, id 0.9991) => the chain does not transfer from `_pri` to a haploid maternal reference.
- n is small by construction (4 expressed mat-absent catalog loci + p12): this is a demonstration per locus, not a rate. No pooled percentage will be quoted.

## 8. Not covered
Between-individual copy differences (one animal); the testis library; loci absent from both haplotypes; reads unmapped on `_pri` beyond R_unm; any threshold
tuning (none permitted); the Amendment 15 consensus fix (separate registered work).

## 9. Compute plan (CRASH RULE: foreground, serial, one heavy job at a time; `/mnt/linuxdisk` has ~25 GB free, so no new index is written)
(1) extract R_LRP + R_unm; (2) map to `mat` then `pat` (existing .mmi, ~14 GB RSS each); (3) fate tables; (4) nets; (5) IsoCon per family (resumable
`iso_batch.sh`); (6) in-house `o3_candidates` per family; (7) link/merge/score; (8) artifact. Register rows to be drafted after results, not before.

---

## Amendment 1 (2026-10-08, written before any read is mapped or any number is looked at) — paternal-reference run and side-by-side

DRAFT for user review. Requested by the user after approving the sections above: show side by side what the same copies look like where they
are PRESENT (the paternal genome) and where they are MISSING (the maternal genome), and how the recovery chain's work in the maternal run compares
with the paternal run, where the same transcripts are simply in the reference. Changes S3 item 4 and S6; everything else stands.

**A1.1 Symmetric runs.** S3-S6 are run twice with one switch (`O3_REF`): reference `mat` (primary, exactly as registered above) and reference `pat`
(truth = `mat`). Same read sets (R34, R_LRP, R_unm), same indexes, same 34 families + LRPAP1, same constants, same bars; nothing is tuned between runs.
**Order:** the `mat` run is completed and its verdicts written to `W/mat/VERDICTS.txt` before the `pat` run starts (a rule validated where it was
developed is not validated: the `pat` run is the substrate held back).
- Truth for `O3_REF=pat`: `bonly.tsv` rows with `hap = mat` (92 loci; known expressed: GWFAM70_B4 24, GWFAM390_B1 20, GWFAM70_B1 16, GWFAM26_B2 13,
  GWFAM26_B3 9, ...), plus LRPAP1 loci present on `mat` and absent from `pat` (expected none: the PAFs decide, as for `mat`).
- Supersedes S3 item 4 ("held out, not run"): the reverse direction is run once, after the `mat` run, and reported per locus.
- Not covered: the two `_pri` copies on mat-sourced chromosomes with lift < 50% that `truth_classes` never refined (GWFAM175:2, GWFAM491:1) are
  not in the truth set of either run.

**A1.2 Side-by-side readouts (descriptive; no new bar).** For every absent locus with >= 3 reads, per run:
- (a) *paired reads*: for each read, `de` and MAPQ of its primary on the reference haplotype (copy missing) and on the other haplotype (copy present);
- (b) *divergence pile*, LARGE loci only: per-bin primary-alignment coverage and mismatch counts (60 bins, positions as fractions of the interval) along
  the nearest reference paralog (where the reads land) and along the locus itself;
- (c) *the chain's view of the same copy in the two runs*: IsoCon outputs whose best hit (identity x coverage >= 0.999) overlaps the locus, counted
  against the truth haplotype in run A and against the reference in run B; how many the link step kept as new copies in each run; the flagged
  candidates recovering the locus (run A); the in-house stage candidates whose nearest locus overlaps it (run B: `d`, flagged).
(a) is partly circular on the `other` haplotype by construction (a read is labelled to a locus by its untied best placement there), so the
`other`-side fate is not reported as a result; only `de`, MAPQ and the count of reads dropped as ambiguous are.

**A1.3 Predictions, before looking.**
- **S1** For every LARGE locus the median `de` on the haplotype that has the copy is smaller than on the haplotype that lacks it; on the present
  haplotype it lies within the registered within-species range (<= delta = 0.00958).
- **S2** The IsoCon outputs of an expressed absent copy exist in both runs (they come from the reads, not from the reference): in run B they match the
  reference at >= 0.999 and none is kept as new (the flag step skips such outputs by construction — shown for consistency, not as a test); in run A they
  are kept as new and, where R1 passes, recovered.
- **S3** In the `pat` run the registered bars R1-R3 are evaluated unchanged for the mother-only copies. No outcome is predicted beyond the `mat` run's;
  it is reported even if it fails.

**A1.4 Visual contract.** The artifact shows both directions behind a toggle; per LARGE locus a left/right pair (reference where the copy is missing |
other haplotype where it is present) with the paired-read plot, the divergence pile, the fate bars and the chain's view.

**A1.5 Cost.** About doubles the compute of S6 (second IsoCon pass, second set of in-house batches, second set of contig alignments). No new index.

---

## Amendment 2 (2026-10-08, written during Task 4 after the truth loci were built and BEFORE any read was labelled or mapped to a locus) — p12 is not a mat-absent locus by the registered rule

Finding (from the PAFs the prereg names as the authority, not from the artifact flags): `truth_lift` classes p12 `T?` (lift 0.914, 112 mismatches) and lifts it to
`CM054594.2:27,426,744-27,441,619`, where `mat` carries a locus matching p12's body at identity 0.972 / coverage 1.0 (asm20 -N 50 -p 0.5). The chr12 LRPAP1 cluster
has 4 loci on `pat` (22.55, 23.07, 24.79, 30.2 Mb) and 3 on `mat` (25.87, 27.43, 34.29 Mb), but the lifts are many-to-one (c01 and c03 both lift to mat 34.29 Mb, c01 with
222 mismatches), so the registered synteny rule cannot say WHICH `pat` locus has no `mat` counterpart. Consequence for S3.2 / S5: under the registered rule the only
LRPAP1 locus absent from `mat` is the chrY copy (sex control); p12 is NOT a mat-absent locus.
Ruling (no new threshold): p12 is kept as a DESCRIPTIVE locus, kind `lrpap1_desc` = an LRPAP1 locus on a chromosome `_pri` took from the truth haplotype whose lift
class is `T?`. It gets a fate row, side-by-side panels and the p12 chain line; it carries NO bar verdict, is not counted in R1-R4's expressed set, and is labelled
in the artifact as "diverged counterpart on the reference". The S5 prediction for p12 ("absorbed on the chr14 ortholog of p14") stands as a prediction about where
its reads go; it is no longer a test of the mat-absent claim. The copy-number statement for LRPAP1 on chr12 is "4 loci on pat, 3 on mat; which one is extra is not
resolved by the registered rule".
Amendment 2, labelling note (same day, before any fate was computed): on `pat` all 83 reads whose primary record lies on p12 are TIED with p14 (score ratio within 0.98; the two
copies differ by 4 bp in 2,446), so the registered untied-placement label gives p12 zero reads. Descriptive loci (`lrpap1_desc`) are therefore labelled by the read's
PRIMARY record on the truth haplotype, ties allowed, restricted to the LRPAP1 net. No catalog locus and no bar uses this rule; the tie itself is reported as a result.
