# PREREG — EXCISION EXTENSION: does the "missing copy" signature generalize beyond NPIP/ZNF875? (2026-09-15)

Question (advisor, RNA-only, no WGS): the original excision experiment (`docs/PREREG_excision_2026-09-05.md`,
confirmed reproduced in `docs/PREREG_o3_flag_pass_2026-09-06.md`) showed that erasing an almost-identical copy from
NPIP/ZNF875 leaves a detectable "consistent sites per kb" signature in the reads that get rejected from their
best surviving candidate. Both source families are unusual (NPIP = large, gene-conversion-heavy retrogene array;
ZNF875 = a single near-100% pair). Before trusting the signature as a general missing-copy detector, test it on
families never touched by that pre-registration: different sizes, different genes, ordinary tandem SD arrays
instead of retrogene arrays, and a broader (not just >=0.97) identity range.

Instrument: unchanged — `bench/o2_excision.py <fam_dir> <X> <out_dir>` (remove copy X from the family's
copies.tsv/copies.fa, re-run `copy_assign --families/--copies-fa/--regions`, follow X's MAPQ-60 reads), scored by
`bench/o3_flag_pass.py --dirs <dirs> --bam npip3.bam --fasta npip3_contigs.fa --gff <gff> --units catalog.units.tsv`
(the exact §6ff detector/Poisson flag rule, unchanged, no new constants). Substrate: the same 3-contig O3 dev
substrate `sweep_v13` (NC_073241.2/073242.2/073244.2), same `npip3.bam`/`npip3_contigs.fa`.

Target selection (chosen by pairwise identity + read count BEFORE any excision was run, from `sweep_v13`'s own
`copies.tsv`/`copies.fa`, restricted to `kept_full`/`kept_trimmed` members only, excluding `fam_MCL1_073242`
which is the already-used NPIP family): for each candidate family, the two `kept` copies with the highest
minimap2 asm20 (`-c --eqx -N 20 -p 0 --secondary=yes`) pairwise identity were found, and the higher-read-count
member of that pair is X, the other is the predicted nearest surviving candidate Y. Spans the family-size range
sweep_v13 actually has (3 to 17 members) and moves off zinc-finger/NPIP retrogene biology onto other gene
families where possible:

| dir name | family | size | gene(s) | X (copy_idx, n_reads) | Y (copy_idx) | identity(X,Y) |
|---|---|---|---|---|---|---|
| `znf578_x4`  | MCL4_073244  | 17 | ZNF578/ZNF808   | 4 (139 reads)  | 1  | 0.927-0.929 |
| `znf430_x2`  | MCL9_073244  |  9 | ZNF430/LOC115931567 | 2 (210 reads)  | 1  | 0.931-0.936 |
| `loc_x6`     | MCL6_073244  |  7 | LOC101129171/92/101130248/101134578/101136009 | 6 (561 reads) | 0 | 0.908-0.915 |
| `znf823_x4`  | MCL43_073244 |  5 | ZNF823          | 4 (99 reads)   | 3  | 0.969 |
| `znf91_x1`   | MCL30_073244 |  3 | ZNF91           | 1 (215 reads)  | 0  | 0.909 |

Measurements: identical to the original pre-registration — fate of X's reads (SILENT misassignment /
origin-rejected / O2 orphan / O3 orphan); consistent-site density (per kb) among origin-rejected reads vs Y,
against the un-excised control and Y's own reads; the Poisson flag rule (`alpha=0.001`, no new constant).

Predictions:
- P1 silent misassignment <= 5% of X's MAPQ-60 reads for every target (matches the original P1 bound).
- P2 the flag rule (`bench/o3_flag_pass.py`) FLAGS `missing_copy` for at least 4 of the 5 new targets — the
  identity range here (0.909-0.969) is at or below the original targets' (0.966-0.99+), so a lower hit rate than
  5/5 would not by itself contradict the original result, but 5/5 flagged is the expected outcome if the
  signature is a general property of "an almost-identical sibling exists" rather than NPIP/ZNF875-specific.
- P3 the flagged targets' consistent-site density is reported (no threshold pre-committed beyond the existing
  Poisson test) so it can be compared against the original 5 (2.9-117.6 sites/kb).
- P4 no post-hoc re-selection: all 5 targets above are run and reported, whether or not each one flags.

Fail rule: if fewer than 3 of 5 new targets flag `missing_copy`, the signature does NOT generalize past the two
original families and O3's excision-based detector is family-specific, not general — report which targets failed
and at what identity/coverage. If P1 fails for any target, report the identity/coverage at which the origin
certificate stops protecting against a missing copy, same as the original pre-registration's fail rule.
