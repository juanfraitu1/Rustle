# Template-guided copy assembly: the assembler uses the family relation — design (2026-10-14)

## Intent (agreed with the user)

**The user's words (2026-10-14):** "I think I would like to ensure that the assembler is also using the family rule, we
technically have better results than stringtie and flair."

**Decisions (user, in order):**
1. Family-guided copy assembly, over running the closed loop. Register 1272 shows O2 moves few reads, and the
   missing copies already have their own primary reads.
2. Scope: **any assembled locus of the gene is a template** (family member or singleton).
3. Approach: **a Python prototype first**, ported into the Rust assembler only if it passes.
4. Design sections 1-3 approved ("yes, lets continue").

**The claim to test.** Copies that have their own reads but no full locus can be assembled from those reads, guided by
an expressed locus of the same gene. This raises Liftoff copy recall without lowering precision and without losing the
in-between operating point's lead over StringTie and FLAIR.

## What exists (verified 2026-10-14)

- **The gap.** These are read-supported Liftoff copy pairs at identity ≥ 0.95, scored against the shipped families.
  The bulk of the loss is copies with no assembled locus, not family grouping.

  | sample | pairs | shared | a copy has no assembled locus | assembled but not in a family | split |
  |---|---|---|---|---|---|
  | human A119b | 514 | 75 | 378 | 45 | 16 |
  | human testis (public, ERR13885926) | 38 | 12 | 24 | 2 | 0 |
  | gorilla KB3781 | 78 | 4 | 71 | 3 | 0 |

- **The operating point helps a little.** The in-between operating point (`--direct-transfrags-union --read-isoform-k
  2`, good pool; `/mnt/linuxdisk/tmp/asm_gap/op/human.OP_good.gtf`) fully assembles 350 of A119b's 799 read-supported
  copy loci. The shipped assembler fully assembles 318 (+43 gained, −11 lost).
- **The ceiling.** 309 extra copies still have no full locus under the operating point. Their gene's own placement is:
  - in one of our families: 61;
  - assembled but a singleton: 46;
  - not assembled: 202.

  So a template exists for at most **107** copies (35%).
- **Most missing copies are multi-exon.** In A119b, of the copy loci without a full locus, 360 are multi-exon and 121
  single-exon. 217 are partly covered by a locus and 264 have none. The single-exon read floor is not the main cause.
- **Tools and inputs.**
  - Python 3 with pysam 0.23.3, plus minimap2, samtools and gffcompare.
  - Splice indexes: human `/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi` (12.3 GB); gorilla
    `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi` (13.6 GB).
  - BAMs and genomes: as `figures/inputs.local.tsv` and `figures/samples.tsv` give them.
  - Representatives: assembled transcripts carry a `reads` attribute.
  - The guided finder (`bench/guided_pipeline.py`) already projects spliced hits (`tx_exon_blocks`, `project_exons`,
    `parse_paf`).
  - Figure 8's helpers (`figures/_liftoff.py`) give the Liftoff loci, read support and the ≥ 50% coverage rule.
- **Refuted neighbours, not re-proposed:**
  - `--tied-seed` (registers 338, 435, 980);
  - family "rescue" levers on non-identifiable copies (647, 648);
  - all secondaries seeding loci (1059);
  - SD-region nodes in semi-guided mode (1007/1009: precision .973 → .833);
  - annotated coordinates as the de novo vertex set (294).

  Here only expressed, read-assembled loci are templates, and every added structure is observed in the copy's own
  primary reads.

## 1. Method (`bench/template_copies/`)

1. **Templates.** Every assembled locus of the base GTF contributes its representative transcript: the transcript with
   the most `reads`, ties broken by the longer exon sum, then the smaller transcript id. The spliced sequence comes from
   the sample's genome.
2. **Projection.** `minimap2 -x splice -N 50` of the template sequences against the species' splice index.
   - Keep a hit when its identity is ≥ 0.80 and ≥ 50% of the template is aligned. These are the guided finder's
     registered thresholds (figure 8, G1).
   - The hit's aligned blocks are the projected exons.
3. **Candidates.** A projected copy whose projected exon union is covered < 50% by every base locus (figure 8's
   coverage rule). The template's own locus is excluded.
   - Candidates found by several templates are merged into one when their projected exon unions overlap by ≥ 50% of
     the smaller union. The highest-identity template's projection is kept and the others are listed.
4. **Evidence.** The reads whose primary alignment (`-F 2308`) has an aligned block on the candidate's projected exon
   union.
   - A candidate goes on only with **≥ 2 such reads**. That is the figure 8 / C2 read-support rule and the assembler's
     read floor.
5. **The copy's transcripts are built only from its own reads.** The projection chooses the region and the reads; the
   reads give every coordinate.
   - A read's aligned blocks (its CIGAR split at `N`) that overlap the projected exons are its exon blocks.
   - Overlapping blocks of different reads merge into one observed exon.
   - A read's `N` gap links the observed exons of the two blocks it separates.
   - Each connected component of observed exons and links is one locus.
     - If the component is a simple path, its transcript is the path.
     - Otherwise each distinct read chain is its own transcript, and chains contained in another are dropped.
   - Nothing is filled from the template, and no boundary tolerance is used.
6. **Output.**
   - `TG_*` loci in a GTF, each tagged with its template locus, identity, own-read count and observed / projected
     junctions.
   - The base GTF and the added loci merged, plus a JSON report: templates, hits, candidates by fate, and added loci.
7. **Discipline.** Reads enter through their primary alignments only. There is no O2 assignment, no EM or fractional
   reads, and no matching or facility location. Families are not changed by this stage; recomputing them on the merged
   GTF is a downstream option.

## 2. Evaluation (pre-registered before any real run)

**Samples.**
- Development: human A119b, on the shipped base and on the operating-point base.
- Held out: human testis (public, ERR13885926) and gorilla KB3781, on the shipped base.
- Tool lead: human A119b and gorilla OR6737 on `OP_primary`, the base register 1268's B1/B2 used.

**Metrics** (per sample; species never pooled):
- **Recall (reported).** Read-supported Liftoff extra copies (identity ≥ 0.95, figure 8's tables) covered ≥ 50% by a
  locus, base versus merged. Also reported: how many of the template-reachable copies were reached.
- **P1, whole set (bar).** Against the base, two quantities must not drop more than 0.2 points: the merged set's
  intron-chain precision against the sample's RefSeq (gffcompare 0.12.10, as figure 1), and the share of its loci at a
  Liftoff locus (figure 8's `located`). These are the closed loop's registered tolerances (LOOP1 b, c). The bar is weak:
  a few hundred added loci barely move whole-set rates.
- **P2, the added loci themselves (bar).** The added loci's `located` share must be at least that of the base loci in
  the same low-support regime: base loci whose representative transcript has 2-5 `reads`, where most added loci
  fall. This is a point-estimate comparison, with a 95% interval over loci reported. This bar is the informative one.
- **Tool lead (bar).** On A119b and OR6737, `OP_primary` + added loci must still pass register 1268's B1 and B2:
  intron-chain sensitivity above StringTie's and FLAIR's, and intron-chain precision above FLAIR's.

**Readings, fixed now:**

| development | held out (both) | reading |
|---|---|---|
| P1, P2 and the tool lead pass | P1 and P2 pass on both | Port into the Rust assembler as an opt-in second pass; the recall gain is the reported effect. |
| pass | fail on either | Not adopted; the held-out failures are listed. |
| any bar fails | — | Not adopted. |

## 3. Components, testing, cost

- **`bench/template_copies/`** (Python, test-first):
  - `templates.py`: representatives and spliced FASTA.
  - `candidates.py`: PAF hits to projected exon chains, the identity and aligned filters, exclusion of covered
    positions, merging of duplicates.
  - `evidence.py`: own primary reads per candidate with pysam, as blocks and junctions.
  - `build.py`: the copy's transcripts from its own reads.
  - `run_tg.py`: the CLI per sample.
  - `score_tg.py`: recall, P1, P2 and the tool lead.
- **Reuse:**
  - `figures/_liftoff.py` for loci, support, coverage and `located`;
  - the guided finder's PAF and projection helpers;
  - figure 1's gffcompare rule;
  - `bench/tool_difference/score_op.py` for B1/B2.
- **Tests:**
  - toy genome, BAM and GTF fixtures (pysam writes the BAMs), with hand-computed transcripts;
  - a minimap2 projection test on a toy FASTA;
  - the edge cases: a candidate with 1 read, an exon pair without a junction read, a candidate overlapping a base
    locus, and two templates hitting one position.
- **Cost.**
  - Projection is a heavy job, since the human index is 12.3 GB. Templates are mapped in batches, each a foreground
    call under the heavy lock that stays under 10 minutes.
  - Evidence, building and scoring are light jobs.
  - The pre-registration is committed before the first real projection.

## Out of scope

- Filling exons from the template (the gap-filling variant).
- Annotation-derived templates (that is the guided mode).
- O2 read moves (the closed loop).
- Changing O1's family rule or the shipped defaults.
- The Rust port, which follows only from reading 1.
- The 202 copies whose gene is not assembled at all; the assembler's low-expression sensitivity is a separate study.

## Results (2026-10-14)

Reading 3: not adopted, for the build tested (locus strand from the template hit, whole read blocks, no
de-duplication; a build with these fixed is untested). On A119b the stage covers 160 more read-supported Liftoff
copies (234 → 394 of 641), 119 of them only by single-exon added loci. Only 0.0813 of its 26,726 added loci (23,415
distinct exon sets) sit at a Liftoff locus, against 0.2318 for base loci with 2-5 reads (and below every base depth
group), so P1 and P2 fail. Full results and the final-review findings: `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md`;
register row 1285.
