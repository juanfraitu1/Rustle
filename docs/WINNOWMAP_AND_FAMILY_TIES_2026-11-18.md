# Winnowmap at multi-copy families, and how large the assignment problem is (2026-11-17 / 18)

One page for rows 1782, 1784, 1787, 1790 and 1791 of `NEGATIVE_RESULTS_REGISTER.md` (and the Winnowmap arm of row 1775).
Each row has its own pre-registration: `PREREG_winnowmap_hg002_2026-11-17.md`, `PREREG_winnowmap2_spliced_2026-11-17.md`,
`PREREG_winnowmap_orangutan_2026-11-18.md`, `PREREG_family_ties_2026-11-18.md`, `PREREG_tie_split_2026-11-18.md`.
Nothing here changes a default. Pages: Winnowmap https://claude.ai/artifact/K1D1JGHm9ciBJ71n6A5gnb, families
https://claude.ai/artifact/EkB8rtLd7ZhqTy77Jb2MBq (private to the user until shared).

## The two questions
1. The user, after an advisor meeting: Winnowmap was built for repeats, so why does it not help with multi-copy gene
   families?
2. The user, after the answer: the primary is right on about 90% of the decidable near ties and the rest cannot be
   decided, so is read assignment a big deal at all?

## Short answers
1. Winnowmap repairs SEEDING in satellites (weighted minimizers keep the rare k-mers of a repeat). At a gene family the
   candidate copies are already found by minimap2; what is missing is a way to RANK them. Winnowmap adds none for
   spliced reads, and when its copy-aware method is forced on it becomes more confident, not more correct.
2. On average no: wrong primaries a better choice could fix are 0.05% of family reads. For a short list of families
   yes: 14 of the 155 expressed families of HG002 have a fifth to all of their reads near-tied.

## Method common to all rows
- Truth: `bench/o2_truth/neutral_truth.py` (old test, `NT_REPAIR` unset). Reads and the reference sequence of each
  candidate locus are aligned to the individual's two haplotype assemblies. Statuses of a read: decidable (one
  candidate locus coincides with the read's best haplotype hit), tie (several do), outside (none does), no
  near-perfect hit. Strata by minimap2's records: one locus, clear lead (> 10 AS), near tie (<= 10 AS).
- Rule for an aligner arm: GAIN = at least 1 percentage point more right than minimap2's primary on the decidable near
  ties (>= 200 reads) AND at most 3 more wrong reads in the clear-lead stratum. Step 0: >= 0.95 mapped, spliced share
  within 0.05 of minimap2.
- Samples: HG002 Kinnex on CHM13 with the T2T v1.1 haplotypes, a 25% sample of the molecules on truth-family genes
  (136,074 reads; `r1776_hg002_room`); KB3781 gorilla and AG06213 Sumatran orangutan on their own assemblies.
- Flags of every alignment: `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`. Winnowmap 2.03 with
  `-W` from `meryl` k = 15, `distinct=0.9998`.

## Part 1: Winnowmap

| Row | Sample, build | Decidable near ties | minimap2 right / wrong | Winnowmap right / wrong | Clear lead, wrong: minimap2 / Winnowmap | Verdict |
|---|---|---|---|---|---|---|
| 1775 | KB3781, stock | 62 | 57 / 5 | 58 / 4 | not judged | no room (5 reads) |
| 1782 | HG002, stock | 711 | 639 / 72 | 605 / 106 | 34 / 58 | NO GAIN, worse |
| 1784 | HG002, patched | 713 | 638 / 75 | 609 / 104 | 34 / 64 | NO GAIN, worse |
| 1787 | orangutan, stock | 174 | 174 / 0 | 174 / 0 | 0 / 1 | NO GAIN, no room |
| 1787 | orangutan, patched | 174 | 174 / 0 | 173 / 1 | 0 / 4 | NO GAIN, no room |

- **Stock (row 1782).** Same primary locus as minimap2 on 98.0% of reads. On the decidable near ties it fixes 18 of
  minimap2's errors and breaks 52 right reads; 54 are wrong in both. Only 197 of 7,770 near ties get a clear lead on
  Winnowmap's own scores. My prediction "fixes fewer than 10, breaks about as many" was wrong.
- **Patched (row 1784).** The splice presets set `SVaware = false`, and the Winnowmap2 substring method only runs on
  reads of 10,000 bases or more. The patch (`r1784_winnowmap2_spliced/patch.diff`, `src/options.c`) turns it on for
  spliced reads from 1,000 bases with substrings of 500 to 4,000. Step 0 passes. It fixes 27 and breaks 56. It gives
  2,816 of the 7,770 near ties (36%) a clear lead and is wrong on 16% of the decidable ones among them; 95.5% of
  reads get MAPQ 60. More confident, not more correct.
- **Would Winnowmap's minimizers inside minimap2 help?** Not tested and advised against: Winnowmap's primary sat at a
  locus with no minimap2 record on 7 of 42,770 reads, and that locus was the origin once. The true copy is among
  minimap2's records on 71 of the 72 primary-wrong near ties. Seeding is not what fails.
- **Tandem copies (post hoc on rows 1782 / 1784; definition fixed before counting: primary and best rival on the same
  chromosome and strand within 1 Mb).** The pipeline already excludes satellites
  (`RUSTLE_EXCLUDE_BED`), so tandem genes were the remaining place for Winnowmap to help.

  | Class of near tie | Decidable | minimap2 right / wrong | stock | patched |
  |---|---|---|---|---|
  | tandem | 85 | 68 / 17 | 69 / 16 | 69 / 17 |
  | same chromosome, inverted or farther | 499 | 444 / 55 | 409 / 90 | 413 / 87 |
  | another chromosome | 127 | 127 / 0 | 127 / 0 | 127 / 0 |

  No effect at tandem copies (85 reads: weak evidence); all of Winnowmap's loss is at same-chromosome copies that are
  inverted or far apart.
- **Second ape (row 1787).** Near ties are 499 of 84,100 reads (0.6%), because the reference is a haplotype of the
  same animal. The tandem class has 90 reads and none decidable (64 "outside": not investigated, possible O3 input).
  The other apes with haplotype assemblies have the same self-reference design. HG002 on CHM13 is the only sample
  pair on hand where the reference is another individual.

## Part 2: how large the problem is (rows 1790, 1791)

HG002 on CHM13, the same 136,074 reads.

| Read group | Reads | Share of family reads |
|---|---|---|
| Near ties | 7,770 | 5.7% |
| - tie by the truth | 4,249 | 3.1% |
| - no near-perfect hit | 2,153 | 1.6% |
| - outside | 655 | 0.5% |
| - decidable | 713 | 0.5% |
| Decidable and primary wrong | 75 | 0.05% |

**Per family (row 1790).** Family = a component of the truth pairs (`famtruth`, level simiiformes); a read belongs to
the family its primary overlaps.

| Set | Families |
|---|---|
| In the annotation (1,459 genes) | 419 |
| With a read in the sample | 228 |
| With >= 50 reads (counted) | 155 |
| No near tie at all | 97 |
| >= 20% of reads near-tied | 14 |
| >= 50% of reads near-tied | 7 |

- CONCENTRATED: half of the near ties sit in 4 families (2.6% of the counted ones).
- LOCALLY LARGE: the 14 families hold 13,274 reads (about 10% of family reads) and 6,273 near ties (81% of all).
- KB3781 on its own assembly: 1.3% near ties, half in 2 families (EOLA1 / EOLA2, FAM156A), 3 families at >= 20%.

**What the ties are (row 1791).** A `tie` of the truth is two things. SEVERAL COPIES: the read fits two or more places
of one HG002 haplotype equally; no sequence difference in the read can place it. ONE PLACE: the candidate loci of the
reference all land on one best hit of the read; the truth cannot judge, and the read may still be placeable. The
registered bar (>= 80% SEVERAL COPIES) was missed: 74.9% in the 14 families, 69.4% over all. The prediction failed.

Near ties of the 14 families (6,273): 44% SEVERAL COPIES (2,785), 15% ONE PLACE (934), 34% no home (2,154: no
near-perfect hit or outside), 6% decidable (400, primary wrong on 55).

| Family | Reads | Near-tied | SEVERAL COPIES | ONE PLACE | Decidable | No home | Primary wrong |
|---|---|---|---|---|---|---|---|
| EIF3C / EIF3CL | 2,234 | 99.4% | 1,750 | 19 | 2 | 450 | 1 |
| FAM156A / B | 384 | 97.4% | 288 | 0 | 0 | 86 | 0 |
| SULT1A3 / 4 | 129 | 96.9% | 80 | 7 | 9 | 29 | 0 |
| RASA4 / RASA4B | 211 | 77.3% | 19 | 53 | 24 | 67 | 5 |
| TBC1D3 (9 genes) | 151 | 63.6% | 0 | 16 | 30 | 50 | 17 |
| SMN1 / SMN2 | 1,516 | 59.7% | 193 | 344 | 117 | 251 | 6 |
| CCZ1 / CCZ1B | 851 | 52.1% | 1 | 232 | 72 | 138 | 1 |
| RGPD (6 genes) | 121 | 46.3% | 7 | 36 | 2 | 11 | 1 |
| BOLA2 (5 genes) | 180 | 44.4% | 37 | 10 | 4 | 29 | 0 |
| CBWD (3 genes) | 951 | 25.6% | 5 | 68 | 72 | 98 | 10 |
| NPIP (23 genes) | 2,241 | 25.5% | 146 | 123 | 54 | 248 | 14 |
| EIF5A / EIF5AL1 | 2,431 | 24.2% | 0 | 0 | 11 | 578 | 0 |
| LRRC37 (4 genes) | 204 | 22.1% | 3 | 2 | 3 | 37 | 0 |
| NOMO (3 genes) | 1,670 | 21.7% | 256 | 24 | 0 | 82 | 0 |

## What follows for the project (my reading; the framing is the user's and the advisor's call)
- Choosing between copies read by read is nearly solved by minimap2's primary, and no aligner or rule tried so far
  improves it (this page and `NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md`).
- The open problem is measurability per family. EIF3C, FAM156, NOMO, SULT1A: copies identical over the reads, PSVs
  cannot help, report one count per group of indistinguishable copies (stage `eqclasses`). SMN, CCZ1, CBWD, RASA4,
  TBC1D3: mostly ONE PLACE, the reads may be placeable but this truth cannot test it.
- A third of the near ties in these families fit no candidate copy well: input for missing copies (O3).
- PSVs from one reference are only valid where the individual's copies differ at the same sites. TBC1D3 (17 of 30
  decidable near ties with a wrong primary) is the lead to follow, on 30 reads.

## Corrections made on the way
- "The true copy is among minimap2's records for all 72 primary-wrong near ties" was 71 of 72.
- Row 1790's first table depended on the run: a read on two families went to either one (set order). Now a stable
  order; verdicts and totals unchanged, RASA4 / RASA4B is 77.3% (first reported 81.9%), the 14 families hold 13,274
  reads (first 13,262).
- "59% of the near ties in the 14 families cannot be placed by any method" is replaced by 44% (row 1791).
- The scorer changed between rows 1782 and 1784 (another session; new column), so the decidable counts differ by two.

## Limits
- One library per sample, a 25% read sample for HG002; only expressed families are seen (191 of 419 have no read, 73
  more fewer than 50). Another tissue would give another list.
- The statuses carry the known defects of the old truth test (`NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md`); the repaired
  test (`NT_REPAIR=1`, row 1785) was not used here. The near-tied shares do not depend on the truth.
- SEVERAL COPIES counts a read whose PSV is cancelled by a read error. A read that fits one copy per haplotype is ONE
  PLACE.
- The tandem split is post hoc with 85 decidable reads; a full-sample HG002 run (about 340 expected) was not done.
- Not done: Winnowmap's minimizers inside minimap2 (advised against); the 64 "outside" tandem reads of the orangutan.

## Where things are (workspace `/anvme/workspace/v138da13-rustle`)
- `r1782_winnowmap_hg002/`: `wm.sbatch`, `wm_ties.py`, `repetitive_k15.txt`, `hg/reads.tsv`, `report/index.html`.
- `r1784_winnowmap2_spliced/`: patched source, `patch.diff`, `wm2.sbatch`, `tandem_split.py`, `hg/reads.tsv`.
- `r1787_winnowmap_orangutan/`: haplotypes, `orang.sbatch`, `pab/reads.tsv`.
- `r1790_family_ties/`: `family_ties.py`, `hg.families.tsv`, `hg.read_family.tsv`, `kb.families.tsv`, `make_page.py`,
  `report/index.html`.
- `r1791_tie_split/`: `tie_split.py`, `split.sbatch`, `tie_split.tsv`, `split14.txt`.
