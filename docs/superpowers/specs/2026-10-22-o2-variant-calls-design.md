# O2 from the sample's own variant calls (pileup + variant calling) — design (2026-10-22)

**Status:** approved by the user on 2026-10-22. Plan: `docs/superpowers/plans/2026-10-22-o2-variant-calls.md`.
Amended the same day while planning:
- K1 is computed on HG002 only;
- PAR families are kept out of the pooled bars.

## Why

The thesis claimed that among AS-tied multimappers a PSV or a BAM tag would always name the true copy. Measured, that fails
in one genome. An exact AS tie means the read covers no difference between the copies that favours one (r1350, r1362), and
the tags are functions of the same read-versus-reference comparison. Where minimap2 is wrong on real reads, the sample's
copies differ from the reference's: in HG002, 73% of ties are ties in the donor too, and 560 of 584 reads cannot be decided
against CHM13 (r1329).

The information that breaks those ties is in the sample's *other* reads. The reads minimap2 places confidently show which
variants this individual carries at each copy. The advisor keeps asking why we do not do a pileup and variant calling. This
design does exactly that, with an established caller, and uses the calls to place the tied reads.

**The user's constraint (2026-10-22):** there will always be one official reference, and we cannot assemble a genome per
sample. The inputs are the reference and the sample's reads. A donor's assembly is only ever the answer key.

## What the user decided (2026-10-22)

- O3 runs before O2 and feeds it: the sample's own variants, inferred from its reads against the reference, are O2's
  evidence.
- An off-the-shelf long-read RNA variant caller does the pileup and calling.
- **The deliverable is one hand-picked family shown end to end.** For each tied read it shows:
  - the candidate placements;
  - the PSVs, the BAM tags and the sample's variants, side by side;
  - our placement against minimap2's, with the truth, or a possible missing copy flagged from the variants.
- **Families:**
  - showcase GTF2IRD2 (HG002);
  - backups NBPF1 and PDPR;
  - DAZ as the honest abstention.
- Any substrate is acceptable.

## Goals

1. **A tool** that, for a named family:
   - places AS-tied reads using the sample's per-copy variant calls;
   - abstains when nothing separates the copies;
   - flags reads that no reference copy explains as a possible missing copy.
2. **A per-family case file** showing the evidence for every tied read.
3. **A pre-registered end-to-end pass** on GTF2IRD2, with HG002's real reads and truth.
4. **The same frozen tool** shown to do no harm in one genome and to hold on a held-out donor.

## Non-goals

- Assembly improvement (StringTie/FLAIR on the placed reads): a later step.
- Phasing or haplotypes as an output.
- A per-sample assembly as an input.
- Indel-based assignment. The caller calls SNVs only; indels are shown in the case file but do not decide.
- A genome-wide run, or any new default in Rustle's Rust code.
- Gorilla numbers pooled with human.
- Resolving DAZ.

## Anti-goals (a failure even if it runs)

- Using the donor assembly, or any truth, as input to the placement or the flag.
- Tuning the caller's parameters, or the rule, on the showcase families or the held-out donor.
- Swapping in a backup family after seeing GTF2IRD2 fail, without reporting the failure.
- Quoting showcase numbers as general accuracy (hand-picked families overstate: metric trap 09-18).
- Comparing O2's precision on its own placements with minimap2's on all reads (r1340). Always compare on the same reads, with
  the selective baseline included.

## Inputs

Paths are under `/mnt/linuxdisk/tmp/o3_graph/` unless absolute. Truth is the answer key only: no step but scoring reads it.

| | HG002 (panels A, C) | A119b testis (panel B) | HG00733 sim (panel E) | one genome (panel D) |
|---|---|---|---|---|
| reads on CHM13 v2.0 | `hg002proof/chm13.bam` (Kinnex FLNC candidates; `splice:hq -uf --secondary=yes -N 100 -p 0.5`) | `yag_testis/yag.bam` | `donorsim/HG00733/sim/chm13.bam` | round-23 worlds `fast/heldout23/*` |
| family copies | `o2psv_dev/pool_copies.tsv` | `r19seq/families19cal_full2_all.json` key `seq19` (as in `yag_testis/run.sh`) | as donorsim built them | as in the worlds |
| truth | `hg002proof/truth.tsv` (+ `truth_untied.tsv`) | none exists | `donorsim/HG00733/sim/truth.tsv` | the worlds' truth |
| genome | `/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa` (the alignment target) | same | same | same |

## The method

### 1. Tied reads (O2's domain, unchanged)

A read is tied when its alignments within 10 AS of its best lie on two or more copies of one family and on no locus off it
(Eichler; `o2correct.py`'s `tied()` and `collect()`). Reads with a unique best (margin ≥ 10) belong to minimap2 and are not
touched.

**One change from `tied()` (my call):** `tied()` also drops a read whose minimap2 primary is more than 9 AS below the best.
Such a read is still AS-tied, on the copies that share the best AS, so it stays in the domain. It is reported in its own row
("primary outside the tied set"). Example: in HG002, 10 of NPIPB2's 12 `no_counterpart` reads have their primary on NPIPB2.2
(AS 1734) while .5 and .7 tie at AS 1819.

### 2. Teachers and calls (the pileup and variant calling)

- **A copy's teachers:** reads whose primary alignment is on that copy and whose every alignment on another copy of the
  family is ≥ 10 AS lower. This is the complement of step 1 and the corrector's definition of confident reads. A tied read
  never teaches.
- **The caller:** the teachers' primary records go into one BAM per family, and **longcallR 2.0.0** calls SNVs on it against
  the reference over each copy's span.
- **Settings:** longcallR's defaults, with one change: `--min-mapq 0`. The teachers are already chosen by AS margin, and
  confident reads in near-identical copies carry low MAPQ.
- **Preset `hifi-masseq` on every panel (my call).** In the v2.0.0 source the two HiFi presets differ only in strand-bias
  filtering:
  - `hifi-masseq` is for "transcript strand" reads, which is what all our libraries are (oriented FLNC or simulated
    transcripts, aligned `-uf`);
  - `hifi-isoseq` assumes reads from both strands, and its strand-bias filter would silently drop real SNPs here, since
    every read of a gene sits on one strand.
- **What the HiFi presets set (v2.0.0 source, checked):**
  - min-depth 6 (the `--help` value of 10 belongs to the ONT presets);
  - phase score 11;
  - a dense window of 100 bp holding ≥ 5 candidates.
- **What the VCF holds:** every candidate site gets a record with a GT, and FILTER is one of:
  - `PASS`: phased het, or hom-alt;
  - `LowQual`: a het that could not be phased;
  - `HomRef`;
  - `Multiallelic`: 1/2;
  - `RnaEdit`;
  - `Dense`.

### 3. A copy's genotype at a column

For a tied read, a **column** is a read base aligned in every one of its tied alignments (the common window,
`common_window()`). Each tied alignment maps the column to a reference position on its copy.

- **Known:** the copy's teachers at that position meet the preset's own depth filter (6 for the HiFi presets), with depth
  counted as the caller counts it. That means the same read filters: read length, divergence, base quality and distance to
  the read end. The genotype is then the alleles of the record's GT, whatever its FILTER (my call):
  - 0/1 → {ref, alt}, unphased (`LowQual`) or phased (`PASS`);
  - 1/1 → {alt};
  - 1/2 (`Multiallelic`) → {alt1, alt2};
  - `HomRef` → {ref};
  - no record → {reference base}.

  Phasing is not needed here, and an unphased het can only make two copies overlap, never separate them falsely.
- **Unknown:** the depth is below 6, or the record is `RnaEdit` or `Dense`. A `Dense` record (≥ 5 candidates in 100 bp) is
  what a collapsed extra copy produces on its nearest reference copy, so it also feeds step 5.

### 4. Placement: assign or abstain (no tuned threshold)

- **Separating column:** a column separates copies c and d when both genotypes are known and share no allele.
- **Vote:** the read's base there, at base quality ≥ the caller's `--min-baseq` (10), votes for c when it is in c's genotype.
- **Beating:** c beats d when at least one separating column votes for c and none votes for d.
- **Assign** c when c beats every other tied copy. Otherwise the output is the set of tied copies that no copy beats (abstain).
- **Heterozygous sites need no phasing:** both of a copy's alleles count as that copy. Sites where the copies share an allele,
  including partial RNA editing, never separate.
- Every placement lists the columns that decided it.

### 5. Possible missing copy (O3: runs first and overrides)

- **Foreign allele:** at a column of a tied read where every tied copy's genotype is known, the read's base (quality ≥ 10) is
  foreign when no tied copy's genotype contains it. An unknown copy never stands in with its reference base (my call). That
  stand-in would make any ordinary SNP of a weakly expressed copy "foreign".
- **Recurrent foreign allele:** carried by ≥ 2 tied reads.
- **Linked:** two recurrent foreign alleles are linked when ≥ 2 reads carry both.
- **Flag:** a connected component of linked recurrent foreign alleles, which is the missing copy's haplotype assembled
  allele by allele. Its reads are the tied reads carrying ≥ 2 of its alleles. Changes of the A-to-I editing type (A>G on the
  copy's transcript strand) are listed but never count toward a flag.
- **A flagged read is reported as "possible missing copy" and is not placed by O2.** This follows the user's 10-18 rule that
  unexplained PSVs override.
- **The teacher side, reported only, never counted toward a bar:**
  - `Dense` records on a copy;
  - teachers carrying ≥ 2 linked alleles their own copy's genotype lacks.

  Both are signs of a collapsed copy. They are kept apart from the bars because teachers cover every untied read, including
  about 258k chrM reads that the truth labels `no_counterpart` only because the HG002-to-CHM13 chain has no chrM.
- **Known limit, stated, not solved:** if an extra copy's reads are confidently placed on a reference copy, its alleles enter
  that copy's genotype (the PRAMEF case, r1346). Step 4's no-shared-allele rule should make such copies abstain rather than
  break. Panel E measures this.

## The case file (per family)

- **`<family>.reads.tsv`:** one row per (tied read, tied copy):
  - whether it is minimap2's primary;
  - AS, NM, de, ms, substitutions (NM less indel bases) and indels (reuse `o2tags.py`);
  - reference-PSV columns for and against, and sample-call columns for and against;
  - the decision, the truth and the outcome (fix, break, both right, both wrong, abstain, flagged).
- **`<family>.columns.tsv`:** one row per (read, separating or foreign column): each copy's position, reference base, sample
  genotype and depth, plus the read's base and quality.
- **`<family>.summary.md`:**
  - the counts per arm;
  - for every fix, which evidence decided it: a sample call that differs from the reference, or one that confirms a reference
    PSV;
  - whether AS, NM and the substitution count alone would have picked right;
  - three worked reads: a fix, an abstention and a flag.

## Arms (all on the same reads)

| arm | what it is |
|---|---|
| N1 | minimap2's primary |
| N2 | minimap2's primary at MAPQ ≥ 1 (the selective baseline, r1340) |
| N3 | this rule with reference-only genotypes: same known/unknown gating, but every known genotype = {reference base}. Isolates what the variant calling adds |
| N4 | the frozen corrector (`o2correct.py donor`, v2) where its output exists |
| VC | this design |

## Panels

- **A. Showcase (HG002 real reads):** GTF2IRD2. Backups NBPF1 and PDPR are always run and reported.
- **B. DAZ (A119b testis):** no truth exists, because no testis donor has an assembly. Descriptive only. VC is expected to
  abstain on nearly all tied reads. Any placement is listed with its column: A119b's DAZ2 and DAZ3 differ at one expressed site
  (r1336).
- **C-human. Missing copy (HG002, every pool family):**
  - **Truth:** tied reads in this design's domain whose donor copy has no CHM13 counterpart (`no_counterpart` in
    `truth.tsv` ∪ `truth_untied.tsv`).
  - **What the domain holds:** the corrector's domain holds 19 of the 31 exact-tie `no_counterpart` reads plus about 90 near
    ties (mostly GTF2H2C and SPDYE1). Step 1's wider domain adds NPIPB2's 10.
  - **Excluded:** families with a chrM copy, because the chain has no chrM, so their `no_counterpart` label is an artifact.
- **C-gorilla. Missing copy (gorilla KB3781 fibroblasts on GGO):** the NPIP family, with the maternal-only unit C as truth
  (r1363–r1365; inputs and truth as in `docs/archive/2026-10/merged/PREREG_o3_npip_maternal_copy_2026-10-21.md`). Its verdict stays apart from
  human's.
- **D. No harm (one genome):** the round-23 worlds. In one genome the calls are the reference, so VC should match N3 and never
  move a right read.
- **E. Held out (HG00733, simulated from its own haplotypes, every family):** run once with the frozen tool. It is never touched
  during development.

## Pass criteria (to be pre-registered; my proposals)

**How reads are scored** (amended 2026-10-22 after the Task 5 review; memory metric trap 09-22, "a scorer that intersects
predictions with the truth universe inflates precision"):
- **Scored reads** are the tied reads with a definite truth.
  - Truth classes: unique, unique_elsewhere, elsewhere, no_counterpart.
  - **Right:** a placement is right only when the placed copy is the truth copy.
  - **Wrong:** every other placement on a scored read. That covers a truth copy outside the tied set, `elsewhere`, and
    `no_counterpart`.
- **Not scored:** `own_tie` reads (they go to K1), and poor_match, unmapped and no-truth reads (reported as counts only).
- **Same reads for every arm.** Precision (wrong placements / placements), fixes and breaks are counted over all scored
  reads, for every arm on the same reads. A break includes a read whose minimap2 primary lies outside the tied set and is
  the truth, which VC moves.
- **Reported alongside.** The decidable subset (truth copy among the tied copies) is reported too, for comparison with
  earlier registers.
- **Not evaluable.** A bar with nothing to measure is "not evaluable", and the panel's verdict is then UNDERPOWERED, not
  FAIL. That covers a zero denominator, no reads for a count bar, and no VC placements for E2.
- **Two failure modes are invisible to that score**, so calibration bars (K1, K2) cover them, as in the 10-18 prereg's T3 and
  T5:
  - **K1 (donor-tied reads):** reads whose truth copies tie in the donor too (`own_tie`), so any placement is a guess.
    VC places ≤ 2% of them. HG002 only: the simulated truths (D, E) name every read's source copy, so a guess there already
    counts as a wrong placement.
- **PAR families** (`famident/par_families.json`: X/Y alleles, not families) stay out of every pooled bar and are reported
  apart, as the user split them on 10-18.
  - **K2 (false flags):** flagged tied reads whose truth copy is present in the reference, ≤ 2% of the tied reads with a
    present copy. A flag counts as a move: it takes the read's placement away.

**Bars by panel:**
- **A (GTF2IRD2):**
  - A1: wrong placements ≤ 2% of VC's placements.
  - A2: fixes ≥ 10.
  - A3: VC fixes more reads than N3.
  - K1 and K2 apply to HG002's pool families together, since panel C runs them all anyway (HG002 is spent, so this is
    calibration, not a fresh test).
  - The backups get the same bars and are reported. If GTF2IRD2 fails, the showcase is FAIL; a backup never replaces it
    silently.
- **B:** descriptive.
- **C-human:**
  - C1: at least one flag holds ≥ 2 `no_counterpart` tied reads. Teacher-side signals do not count.
  - C2 (reported): every flag is listed with its reads' truth classes. The audit checks each counted flag's donor locus for a
    CHM13 copy, to rule out a lift artifact.
- **C-gorilla:**
  - C1g: at least one flag holds ≥ 2 tied reads from unit C.
  - C2g (reported), as C2.
- **D:**
  - wrong placements ≤ 2% of placements;
  - reads N1 places right that VC moves or flags ≤ 2% of them.
- **E:**
  - E1: wrong placements ≤ 2% of placements;
  - E2: fixes > breaks;
  - K2.

## Order of work

1. Install longcallR 2.0.0 from bioconda in its own env, `longcallr`, and record the version. On a development family,
   confirm the preset values and FILTER classes read from the source (step 2).
2. Write `bench/o3_graph/o2vc.py` with `test_o2vc.py`, test-first:
   - genotype sets, separation, beating and foreign alleles;
   - VCF parsing;
   - mapping columns across strands and through secondary records without SEQ, reusing `o2correct.py`'s machinery.
3. Develop only on spent HG002 families that are in no panel: SMN2 and MSTRG3150.
4. Write and commit `docs/archive/2026-10/merged/PREREG_o2_variant_calls_2026-10-22.md`, with the bars above and the code commit frozen.
5. Run panels A–E once each, serially, as light jobs (`tools/rlock.sh light`).
6. Audit with the `metric-auditor` agent, then add the register row (next number 1374), update memory and report.

## Practical

- **Location:** branch `o3-fast`, worktree `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset`; outputs to
  `/mnt/linuxdisk/tmp/o3_graph/o2vc/`.
- **WSL2 crash rule:** foreground, serial, one heavy job at a time. These jobs are light (per-family BAMs).
- **Before the prereg and scorer:** read the metric-trap catalogue (memory `feedback_metric_traps`).

## Calls I made (for the user to check)

1. **longcallR rather than Clair3-RNA:** a single binary on bioconda, built for long-read RNA, and needing no Docker, which this
   machine lacks.
2. **Caller settings:** defaults except `--min-mapq 0`, with preset `hifi-masseq` on every panel.
3. **Known/unknown:**
   - known when depth ≥ the preset's 6;
   - the genotype is the record's GT whatever its FILTER;
   - `RnaEdit` and `Dense` count as unknown.
4. **One separating column is enough** for one copy to beat another, provided no column contradicts it. The case file reports
   how many placements rest on a single column.
5. **The flag:**
   - built from foreign alleles found only where every tied copy is known;
   - recurrent alleles linked by ≥ 2 shared reads;
   - editing-type changes excluded;
   - teacher-side signals reported, never counted.
6. **The domain** includes AS-tied reads whose minimap2 primary lies outside the tied set.
7. **Development** happens on SMN2 and MSTRG3150 only.
8. **The pass bars** above, including the ≥ 10-fix floor for the showcase and the calibration bars K1/K2.

## Left to the builder

- Parsing the VCF.
- One longcallR run per copy span or one per family BAM.
- The case-file layout beyond the fields listed.
- Which reads become the worked examples.
- Locating the gorilla NPIP inputs from the r1363 prereg.

## Context: related work and known gaps (2026-10-22; notes, not requirements)

Recorded after the user compared this design with another agent's PSV-model suggestion.

**Prior art:**
- **Paraphase.** Chen et al., AJHG 110(2), 2023, doi:10.1016/j.ajhg.2023.01.001. It phases full-length paralog haplotypes from HiFi DNA reads and gives copy number and phased variants (SMN1/SMN2; later extended genome-wide to 316 genes). This is the DNA analogue of this design's O3 → O2 step. `bench/o3_graph/samplecopy.py` is already Paraphase-like.
- **Parascopy.** Prodanov & Bansal, Nat Commun 13, 2022, doi:10.1038/s41467-022-30930-3. It estimates paralog-specific copy number from short-read WGS jointly with each PSV's allele frequency across many samples, which amounts to PSV reliability.

**What this design already covers:**
- PSV observations taken from the read's own tied alignments, with no multiple alignment of the copies.
- A sample SNP or gene conversion at a reference PSV. These enter the copy's called genotype, and the column stops separating: `erased_ref`.
- Homopolymer indel noise. Only SNVs decide.
- Reads that no copy explains. They are flagged as a possible missing copy, not forced onto a copy.

**Known gaps:**
1. **Real indel PSVs are not used.** longcallR calls SNVs only. HG002's SMN and THOC3 copies differ from CHM13 by a 1-bp insertion. This is the first extension if the showcase falls short, and it needs its own prereg.
2. **A copy collapsed onto a reference copy is only reported, not fixed.** The signals are `Dense` records and weak hets. A Paraphase-style haplotype count per family is the natural fix.
3. **PSV reliability is learned in one pass, from confident reads only.** Tied reads never feed back. Iterative EM was tried before and did not help (r579). Rescue/EM arms were harmful on real reads (r1351–1359).
4. **Weak hets can erase a reference PSV and create placements.** On SMN2, one LowQual 0/1 call at AF 0.14 erased the counter-vote for 688 reads. Task 6 compares three genotype policies on the development families and freezes one before the prereg.

**Deliberately not used:** a facility-location solver on top. It was assessed for O2 and rejected (47–59% precise). Each read is decided on its own evidence.
