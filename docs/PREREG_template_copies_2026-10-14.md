# Template-guided copy assembly (2026-10-14)

**Pre-registered before the first real projection.**
- Spec: `docs/superpowers/specs/2026-10-14-template-copies-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-template-copies.md`.
- Code: `bench/template_copies/` at commit `79ea8caa`.
- Every number below comes from that code, run as written in the plan's Tasks 6-9.

## Question

Can copies that have their own reads but no full locus be assembled from those reads, guided by an expressed locus of
the same gene, without lowering precision and without losing the operating point's lead over StringTie and FLAIR? This
is the assembler using the family relation (user, 2026-10-14).

## Known before running (recorded now)

- **Where read-supported Liftoff copy pairs are lost** (identity ≥ 0.95, shipped families):

  | sample | pairs | shared | a copy has no assembled locus | assembled but not in a family | split |
  |---|---|---|---|---|---|
  | human A119b | 514 | 75 | 378 | 45 | 16 |
  | human testis (public, ERR13885926) | 38 | 12 | 24 | 2 | 0 |
  | gorilla KB3781 | 78 | 4 | 71 | 3 | 0 |

- **The operating point helps a little.** Of A119b's 799 read-supported copy loci, the shipped assembler fully covers
  318 and `OP_good` 350 (+43, −11).
- **The ceiling.** 309 extra copies have no full locus under `OP_good`. Their source placement is:
  - in a family: 61;
  - an assembled singleton: 46;
  - not assembled: 202.

  So a template exists for at most 107.
- **Copy shape.** Of A119b's copy loci without a full locus, 360 are multi-exon and 121 single-exon; 217 are partly
  covered and 264 have no locus at all.
- **Calibration (planning).** One 1,000-template batch of A119b's shipped templates was projected on the human index
  (89 s, 16.3 GB; `/mnt/linuxdisk/tmp/o1_template_copies/calib/`) to size the batches. Nothing from it was scored or
  used.
- **Refuted neighbours, not re-proposed:**
  - `--tied-seed` (registers 338, 435, 980);
  - family rescue levers (647, 648);
  - all secondaries seeding loci (1059);
  - SD-region semi-guided nodes (1007/1009);
  - annotated coordinates as de novo vertices (294);
  - the closed loop is expected to be small (1272).
- **Templates per run** (inputs, from `run_tg.py templates`):

  | run | base GTF | templates | batches |
  |---|---|---|---|
  | A119b_shipped | `runs/human_A119b/human_A119b.gtf` | 94,685 | 16 |
  | A119b_OP | `asm_gap/op/human.OP_primary.gtf` | 84,734 | 15 |
  | testis | `runs/human_testis/human_testis.gtf` | 13,048 | 3 |
  | KB3781 | `runs/gorilla_KB3781/gorilla_KB3781.gtf` | 17,672 | 3 |
  | OR6737_OP | `asm_gap/op/gorilla.OP_primary.gtf` | 19,983 | 4 |

## Method (fixed)

As the spec's section 1, with the plan's rulings R1-R3:
- identity ≥ 0.80 and ≥ 50% of the template aligned (`minimap2 -x splice -c -N 50 -p 0.1`);
- a candidate is covered < 50% by every base locus, and duplicates overlapping by ≥ 50% of the smaller union merge;
- ≥ 2 own primary reads (`-F 2308`) per candidate and per emitted locus;
- transcripts come only from the own reads' blocks and junctions;
- the run fails when more than half the candidates lie on contigs absent from the BAM.

## Metrics and bars (per run; species never pooled)

- **Recall (reported).** Read-supported Liftoff extra copies (identity ≥ 0.95) covered ≥ 50% by a locus, base vs
  merged; the template-reachable ones and how many were reached.
- **P1 (bar).** Two quantities must not drop: the merged set's gffcompare intron-chain precision against the sample's
  RefSeq must be ≥ the base's − 0.2 (percentage points), and its `located` share ≥ the base's − 0.002.
- **P2 (bar).** The added loci's `located` share must be ≥ that of the base loci whose representative has 2-5 reads,
  with at least one locus added. A 95% interval over the added loci is reported.
- **Tool lead (bar; A119b_OP and OR6737_OP).** Register 1268's scorer, genome-wide minus the development chromosome
  (human chr20, gorilla NC_073244.2). The arm `OP_primary_TG` (merged) must have intron-chain sensitivity > StringTie's
  and > FLAIR's (B1), and intron-chain precision > FLAIR's (B2).

## Readings (fixed now)

| development: A119b_shipped and A119b_OP P1 + P2, tool lead on A119b_OP and OR6737_OP | held out: testis and KB3781 P1 + P2 | reading |
|---|---|---|
| all pass | all pass | Port into the Rust assembler as an opt-in second pass; the recall gain is the reported effect. |
| all pass | any fails | Not adopted; the held-out failures are listed. |
| any fails | — | Not adopted. |

## Rulings (made while planning)

- **R1.** Each emitted locus needs ≥ 2 own reads. The spec's ≥ 2 rule is applied to every component as well as to the
  candidate.
- **R2.** minimap2 also gets two flags:
  - `-c`, which gives the CIGARs; without them there are no projected exons;
  - `-p 0.1`, the guided finder's setting. minimap2's default 0.8 would drop secondary hits below 80% of the best
    score, which are exactly the diverged copies.
- **R3.** `assemble` fails when more than half the candidates lie on contigs absent from the BAM. That guards against
  an index and a BAM that name contigs differently.
- **R4.** Development's operating-point base is `OP_primary`, the base register 1268's B1/B2 judge. `OP_good` is not
  re-run (compute); its pre-design ceiling numbers are reported as context.
- **R5.** The tool lead reads register 1268's scorer's arm table, and `score_tg.lead` evaluates B1 and B2 for the TG
  arm. The scorer prints its own claims only for its four standard arms.
- **R6.** Batches hold 6,000 templates. One 1,000-template batch took 89 s and 16.3 GB on the human index, index load
  included.
- **R7.** P2 compares against base loci whose representative transcript has 2-5 `reads` (from `templates.tsv`). P2
  passes only when at least one locus was added.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
