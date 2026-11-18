# Read placement at multi-copy families on the neutral truth, and a defect in that truth (2026-11-17)

> Naming corrected 2026-11-18: "Eichler's rule" here always meant our adaptation to one reference of the Eichler lab's 10-point filter (Guitart 2024); the paper's own procedure is different and was not tested. Definition: `docs/copy_assignment_definition.md`. Rules, bars and predictions are unchanged.

One page for rows 1775-1778, 1781 and 1783 of `NEGATIVE_RESULTS_REGISTER.md`. Each row has its own pre-registration
(`PREREG_eichler_baseline`, `PREREG_hg002_room`, `PREREG_allele_graph`, `PREREG_assign_arm`, `PREREG_tie_ms_arm`,
`PREREG_own_psv`, all dated 2026-11-17). Nothing here changes a default.

## Why this was done
The advisor asked for a baseline for multi-copy gene families (the adapted Eichler filter: an alignment more than 10 AS above the rest
gets the read) and for evidence that the method works there, because the assembler alone is not the selling point. The
user's reasoning: the family gives the candidate copies, the sites where copies differ (PSVs) tell them apart.

## The truth used
`bench/o2_truth/neutral_truth.py`: a read's placements are its records merged into loci; the read and each placement's
reference sequence are aligned to the individual's two haplotype assemblies; the origin is the placement whose haplotype
locus overlaps the read's near-perfect best hit. No catalogue and no arm enters it. Strata: clear lead (best locus more than
10 AS above the next) and near tie (within 10). HG002 on CHM13 with the T2T v1.1 haplotypes, a 25% sample of the
molecules on truth-family genes: 7,770 near ties, 710 of them with one knowable origin ("decidable").

## Results on the 710 decidable near-tie reads (primary alignment: 639 right, 71 wrong)

| Row | Rule | Own calls right / wrong | Moves off the primary right / wrong | Rule, else primary: right / wrong | Verdict |
|---|---|---|---|---|---|
| 1777 | graph of the produced GTF's exons + the sample's alleles, unchecked | 117 / 65 | 7 / 48 | 608 / 102 | FAIL |
| 1778 | default assign stage (exact AS ties) | 5 / 2 | 0 / 1 | 638 / 72 | FAIL, abstains |
| 1781 | assign stage with `--tie-ms 10` | 241 / 32 | 5 / 24 | 623 / 87 | FAIL |
| 1783 | own alleles from MAPQ-60 reads, used only where the two copies' allele sets are disjoint | 229 / 16 | 5 / 12 | 632 / 78 | FAIL |

Clear-lead reads: the adapted Eichler filter read as "move to the best AS" is wrong on the reads where the best-AS locus is not the
primary (KB3781 763 of 763, HG002 479 of 488 against it; rows 1775, 1776). Guard reads: row 1777 moved 74, all wrong; rows
1778, 1781 and 1783 moved none.

What held in every row:
- The family step supplies the options: for 57 of the 71 primary-wrong reads the wrong locus and the labelled origin are
  copies of one family (row 1778).
- The default assign stage is safe: it abstains. Reference PSVs do not separate the copies for the reads it gates (17 of 24
  have two or more deciding sites and a margin of 0).
- Every rule that makes more calls makes more errors.

## Why the own-allele rule moved reads wrongly (row 1783 diagnostic)
The 12 wrong moves are three locus pairs (chr1 247.89 / 247.93 Mb 8 reads, chr7 103.98 / 103.88 Mb 3, chr9 42.1 / 80.1 Mb 1).
In all of them the read's bases are present at the labelled origin among primaries of MAPQ below 60 and absent from its
MAPQ-60 reads. MAPQ 60 selects against the reads of a copy that resemble another copy, so a teaching set filtered on it
makes two copies look more different than they are. Row 1333's proposed fix ("teach only with clearly placed reads") has
this bias.

## The defect in the truth
Found by aligning the 12 reads and the chr1 segments to both haplotypes (`r1783_own_psv/hapcheck.py`, `readhap.py`).

chr1, 8 reads. On both haplotypes the read's best hit has 0-5 mismatches and lies at the positional counterpart of the
copy the rule called, with that copy's strand and 153-kb span; at the counterpart of the labelled origin it has 8-14
mismatches. The labelled origin has one form on both haplotypes. The rule was right; the label and the primary alignment
were wrong.

Two causes in `neutral_truth.py`:
1. The origin test is a span overlap on (haplotype, chromosome): no strand, no exon overlap. A spliced placement spanning
   109 kb contains the other copy on the opposite strand.
2. A placement whose CHM13 sequence has no near-perfect haplotype hit is dropped as a candidate. Here the reference's copy
   is 2% from HG002's copy at that position, so the read's real placement was not a candidate and the other one was
   labelled the origin. The read should have been undecidable, or the placement located by position.

chr7 (3 reads): best hits with 0-2 mismatches on the maternal haplotype, at the called copy's counterpart for two reads
and the labelled origin's for one; unresolved. chr9 (1 read): best hit (1 mismatch) at a maternal locus that matches none
of the main CHM13 placements in position and strand; a candidate for a copy the reference lacks, not followed.

## What stands and what does not
- Stands: the room exists (the primary alignment is wrong for some near-tie reads); the family step supplies the options;
  the default assign stage abstains; the unchecked allele graph and the wider gate add errors on clear-lead and near-tie
  reads; MAPQ 60 is a biased teaching filter.
- Open: every count on the 710 reads. The number of labels of the chr1 kind is not measured. If only the 8 chr1 reads
  flip, row 1783's moves are 13 right / 4 wrong and the primary has 8 more errors; that is not claimed.
- Not to do before the repair: score another rule on this table, or quote the table to the advisor as final.

## The repair as planned (run the same day as row 1785: results in the next section)
1. Origin = a placement whose haplotype hit is on the same strand as the read's hit and shares exon bases with it.
2. A read is undecidable when a competitor within 10 AS of its best has no near-perfect haplotype hit.
3. Re-score every arm (primary, both adaptations of the Eichler filter, Winnowmap on KB3781, rows 1777, 1778, 1781, 1783) on the new
   table, HG002 and KB3781. The alignments exist; about half an hour of compute.
The bars stay those of the rows; the repaired table is reported next to the old one, not instead of it.

## After the repair (rows 1785, 1786; added the same day)
The repair was run under `NT_REPAIR=1` (opt-in) and every arm re-scored; the counts above this section are pre-repair.

HG002 near ties: 600 decidable (566 now `competitor unlocatable`, the 8 chr1 reads among them), primary 553 right / 47 wrong.

| Row | Rule | Own calls right / wrong | Moves right / wrong | Rule, else primary |
|---|---|---|---|---|
| 1777 | allele graph, unchecked | 108 / 44 | 5 / 28 | 536 / 64 |
| 1778 | default assign stage | 5 / 1 | 0 / 0 | 552 / 48 |
| 1781 | assign with `--tie-ms 10` | 235 / 19 | 4 / 15 | 543 / 57 |
| 1783 | own alleles, certified | 205 / 5 | 3 / 4 | 552 / 48 |

No rule gains on the primary's 47 wrong; the verdicts are unchanged. The own-allele rule's calls are now 2.4% wrong, but it
confirms the primary rather than correcting it.

The adapted Eichler filter (row 1786), wrong share of "take the best AS" against "keep the primary" by lead: HG002 6-10: 16.4 vs
2.0%; 11-20: 7.7 vs 1.0%; 21-50: 19.2 vs 0.8%; 51-200: 2.4 vs 0.07%. KB3781 11-20: 36.8 vs 0%; 21-50: 42.4 vs 0.14%.
There is no threshold at 10, moving to the best AS above 10 is not safe, and keeping the primary is. All of those errors
are reads where AS and minimap2's ranking score disagree (3% of clear-lead reads).

The filter itself was then checked against its source (2026-11-18): the Eichler lab's filter (Guitart 2024) compares primary
alignments of separate mappings to haplotype assemblies, per paralog cluster. What rows 1775 and 1786 test is that filter
adapted to one reference and to individual copies, in two forms ("keep the primary" and "move to the best AS"); neither
is the paper's procedure, which has not been run. Definition and wording: `docs/copy_assignment_definition.md`. Page
for the advisor: artifact "The Eichler lab's 10-point filter as a baseline for multi-copy gene families".

Open: `NT_REPAIR=1` is still opt-in (the user decides on the default; rows 1782 and 1784 of another session scored
without it); the chr7 and chr9 cases of row 1783 are unresolved.

## Files
Workspaces `r1775_eichler_baseline`, `r1776_hg002_room`, `r1777_allele_graph`, `r1778_assign_arm`, `r1781_tie_ms_arm`,
`r1783_own_psv` under `/anvme/workspace/v138da13-rustle`. Scripts `bench/o2_truth/{neutral_truth,allele_graph,assign_arm,own_psv}.py`.
