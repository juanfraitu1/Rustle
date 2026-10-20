# Pre-registration: the second round -- loci rebuilt where O2's corrected multimappers have a say (2026-10-20)

**User (2026-10-20):** "it needs to make loci where the corrected multimappers now have a say, that is why I suggest it as
a post-assembly second round". Then, to the proposal of a held-out test on other known families: "yes please".

## The second round

This is the closed loop of `docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md`, with O2's corrector in place of the
union certificate. The certificate assigned 0 reads, so the loop never changed anything.

- **Pass 1** is the shipped `assemble` stage. Loci are seeded from primaries plus secondaries with AS >= 0.98 x the read's
  best. A read tied between two copies therefore echoes at both.
- **The home table** (`bench/loop_home.py o2 --every-read`, b15b5948) gives every read whose primary lies on a family copy one
  home:
  - O2's copy, when the corrector (`bench/o3_graph/o2correct.py donor`, version 2) decided it;
  - minimap2's primary copy, when the read is not tied within 10 AS (Eichler: minimap2's call stands);
  - no home for a tie the corrector could not break: it keeps pass 1's echo.
- **Pass 2** is the same assembly with `RUSTLE_READ_HOME_TABLE`. A homed read enters the read pool only through its records at
  home.

## Development: every number seen before this file

All of it is HG002, chr5, SMN1/SMN2, scored as exon-7 inclusion of the assembled loci against the donor truth (0.617 / 0.985),
StringTie read on `cov` (Amendment 2 of register 1344):

| arm | SMN2 | SMN1 |
|---|---|---|
| pass 2 with O2 homes | 0.585 | 0.969 |
| primaries only | 0.602 | 0.801 |
| pass 1 | 0.743 | 0.744 |
| StringTie | 0.654 | 0.894 |
| FLAIR | 0.763 | 0.897 |

SMN is excluded below. The corrector's read-level calls were seen in development for every HG002 family. What is held out
here is the **loci** of the families below: no second-round output for any of them has been looked at.

## Families (fixed by a power rule, not by results)

Every qualifying non-PAR family (`famident/qualify_nonpar.json`: every member pair at [0.90, 1.00) identity) with at least 20
of the corrector's domain reads carrying a donor truth (unique / unique_elsewhere), SMN excluded. There are 15:
AC0064533, AC1275021, DUSP22, FCGR2A, GOLGA8A, GTF2H2C, GTF2IRD2, HSPA6, NCF1, NOMO1, PDPR, SERF1B, TBC1D3B, THOC3 and
ZC3H11A.

File: `/mnt/linuxdisk/tmp/o3_graph/famident/second_round_families.json`, sha1 6a6624394d9c.

## Arms (HG002 Kinnex, `hg002proof/chm13.bam`; CHM13 v2.0)

| arm | what |
|---|---|
| **LOOP-O2** | pass 2 with the home table above |
| **P** | primaries only (`--seed-pool primary`): minimap2's one home for every read |
| **GOOD** | pass 1, the shipped assembly |
| **ST** | StringTie 3.0.3 `-L` on the primaries (`-F 2308`) over the 15 families' copy spans |
| **FL** | FLAIR 3.0.0 `bam2bed` + `collapse --trust_ends --generate_map` on the same primaries |
| ORACLE | pass 2 with the donor truth's homes (`--every-read --truth`, never O2's calls); reported, not in a bar |

- **Assembler.** The Rustle arms use `tools/rustle_pipeline.sh assemble` genome-wide with one best-AS table, binary
  `rustle_target/release`. The binary is unchanged in `src/` since it was built.
- **Corrector.** Its output is regenerated at this commit and must equal `o2psv_dev/hg002.v2.tsv` (audit of register 1344).

## Metric

**Truth shares.** Per family, the family's reads whose donor-truth class is unique / unique_elsewhere, counted by their true
copy (`truth.tsv` + `truth_untied.tsv`).

**Arm shares.** Per family, the arm's assembled transcripts on each copy (`summary.attribute`: strand, then exon overlap),
weighted:
- `reads` for the Rustle arms;
- `cov` for StringTie, read by name;
- mapped reads for FLAIR.

**Per-family score.** The total variation distance (TVD) between the arm's shares and the truth's. An arm with no transcript
on the family scores 1. Scorer: `bench/o3_graph/o2second.py score`.

## Bars

- **S1:** the pooled TVD (the mean over families, weighted by truth reads) of LOOP-O2 is below that of P, GOOD, ST and FL,
  each.
- **S2:** LOOP-O2's TVD is lower than P's in more families than it is higher, and the same against GOOD.
- **S3 (power):** at least 5 families have truth reads.

**Verdict:**
- **PASS:** S1 and S2 hold.
- **FAIL:** one of them fails.
- **UNDERPOWERED:** S3 fails.

The unweighted family mean is reported beside the pooled one (THOC3 carries the most reads).

**A PASS would claim:** on 15 near-identical families the corrector never had its loci judged on, rebuilding loci after O2
gives each copy its own reads. The loci are closer to the donor's truth than the shipped assembly, minimap2's primaries,
StringTie or FLAIR.

**It would not claim:**
- another sample: HG002 is the only one with per-read truth, and a simulation follow-up is separate;
- transcript structures, since the metric is the read share per copy;
- copies the reference lacks (O3).

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2second/run.sh`.
- A bug is recorded as an amendment before any rerun.
