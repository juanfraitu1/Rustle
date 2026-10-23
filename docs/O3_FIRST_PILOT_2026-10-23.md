EXPLORATORY (pilot, post hoc; no pre-registration, nothing claimed)

# O1 -> O3 -> O2 pilot on HG002, with Strainy and the advisor's silent-copy check (2026-10-22/23)

**Code:** `bench/o3_graph/o3first_pilot.py` (`copies`, `silent`, `strainy-inputs`, `score`). **Work dir:** `/mnt/linuxdisk/tmp/o3_graph/o3first`.
**Inputs:** HG002 Kinnex reads on CHM13 (`hg002proof/chm13.bam`); the O2-VC copy table (`o2psv_dev/pool_copies.tsv`). The table lists
every exon variant of a copy (793 of its 1,686 copies overlap), so each copy is merged into its exon union; chrM copies are dropped.
Families: the 94 with tied reads from the 2026-10-22 spike, with the same labels: positives (a) SPDYE1, RASA4B, USP17L11, POTEM,
BMS1P1; (b) = (a) + NPIPB2; 73 others; the 15 PAR families are excluded.

## 1. O3 before O2: family-level detection (every arm on the same 79 families)
- `candidate_copies` is the shipped opt-in O3 stage, run with its defaults on each family's whole read net before any tied read
  is placed: 7 batches, 4-6.5 min each, 14.4 GB peak.
- Strainy 1.2 (bioconda) ran as `-s phase -m hifi`, with the family's copies (spliced exon unions) as the collapsed
  sequences and every read with a record on them. `--max-unitig-coverage` was raised to 100000 because RNA depth exceeds the
  500 default; every other setting was left at Strainy's default, including 20x minimum coverage, 1 kb minimum length and MAPQ >= 20.
  - A copy's haplotypes = the most Strainy clusters (>= 2 reads each, unclustered excluded) that cover one of Strainy's SNP
    positions together.
  - A copy is flagged when that count is larger than its allele slots: 2, or 1 for a non-PAR X/Y copy, since HG002 is male.

| arm | (a) /5 | (b) /6 | others flagged /73 |
|---|---|---|---|
| Strainy: haplotypes > allele slots | 1 SPDYE1 | 2 (+NPIPB2) | 7: NBPF1, AC1275021, CBWD2, SMG1, GTF2H2C, ANKRD36C, DUSP22 |
| candidate_copies: flagged candidate | 1 SPDYE1 | 1 | 14 |
| candidate_copies: any candidate | 4 SPDYE1, RASA4B, POTEM, BMS1P1 | 4 | 20 |
| spike tier 1 (10-22) | 2 SPDYE1, RASA4B | 3 | 9 |
| today's flag (O2-VC, 10-22) | 1 SPDYE1 | 2 | 7 |

- **Strainy coverage.** Strainy phased >= 1 copy in 50 of the 79 families. USP17L11 and POTEM have no copy at >= 20 reads, so
  Strainy never phased them. GOLGA8A hit the 280 s limit. SPDYE1.20 and NPIPB2.16 each carry 3 concurrent haplotypes.
  RASA4B's 2 phased copies have 1 each, and BMS1P1's 3 phased copies have at most 2.
- **Labels.** 3 of Strainy's 7 "false alarms" (CBWD2, NBPF1, GTF2H2C) are families where the truth labels hundreds of reads as
  coming from donor copies with no CHM13 counterpart: 395, 243 and 212 reads. They were not among the verified positives and may
  be liftover artifacts. Every arm's count moves if they are real, so the label set must be verified before any arm is compared.
- **Candidate support.** candidate_copies' unflagged candidates in RASA4B, POTEM and BMS1P1 have 3-4 reads each, the same scale as
  its candidates in the other families.

## 2. The advisor's claim: a silent copy is found through its family's multimapping reads
- **Definitions.** A silent copy has no read with its primary record on it; an expressed copy has >= 3. A silent copy is
  reached when a read of an expressed copy of its family has a secondary record with >= 100 bp on it. Identity and coverage
  come from minimap2 -c -x asm20 against the nearest expressed copy.
- **Result.** 236 silent copies in 75 families (NUMT families excluded); 179 (76%) are reached, 166 (70%) by >= 2 reads.

| identity | silent | reached |
|---|---|---|
| >= .99 | 66 | 66 (7 of them fragments covering < 50% of the gene) |
| .97-.99 | 41 | 32 |
| .90-.97 | 104 | 64 |
| < .90, or no asm20 hit | 25 | 17 |

Expression is not the limit, divergence is. These are annotated copies; unannotated ones were not measured.

## Not done
- The read-level O2 effect: whether the donor-only copies' reads move to the candidates, and the false moves.
- Verifying the labels of CBWD2, NBPF1, GTF2H2C and WASHC1.
- A Strainy run with a lower coverage floor.
