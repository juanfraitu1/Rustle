# Report-only: O3 -> O2 on KB3781's real reference-absent copies (refabsent panel), 2026-10-24

**Status:** written before the port is run on these families. The refabsent panel (34 families with haplotype-only loci, truth
from KB3781's own haplotypes, registers 1205-1207) has been used before (O3 stage, 2026-10-01), so this is a REPORT, not a test:
no verdict, nothing tuned on it.

**Method.** `o3port.py --semi --guard` as of 11779a31 (P14's method), on `refabsent/R0.bam` (fibroblast reads on `_pri`), copies =
`refabsent/copies_genes.tsv`, read set = `refabsent/fam/<family>.fa`. Script: `bench/o3_graph/o2refabsent.py`.

**Truth.** Each read's best record on the haplotype the family's copies lift to: on a lifted copy -> that copy; on a B-only locus
-> absent; a tie between two loci -> unknown (dropped). Within-delta B-only loci (114 of 127) are near-identical to a `_pri` locus,
so their truth is close to allelic and is reported apart from beyond-delta loci.

**Reported.** (1) Per family: expressed B-only loci vs the port's certified / putative flags. (2) O2 on tied reads from reference
copies: placed, wrong, vs minimap2 MAPQ >= 1. (3) O2 on tied reads from B-only loci: placed anywhere but a certified new copy =
wrong; vs minimap2 MAPQ >= 1 (always wrong there).

## Result (run 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2refabsent`)
Truth: 16,236 reads from reference copies, 3,408 from B-only loci, 1,108 tied between loci (dropped). A first truth build mapped
haplotype contigs by their `chrN` prefix, so `chrN_*_random_utig*` overwrote `chrN` and 87 correct GWFAM64 placements scored as
wrong; fixed before the numbers below (the port's output was not rerun or changed).

**O3 flags vs expressed B-only loci** (8 families have one; reads in parentheses):

| family | expressed B-only loci | certified | putative |
|---|---|---|---|
| GWFAM175 | B0 (281), B1 (6) | 0 | 0 |
| GWFAM205 | B0 (77) | 1 | 0 |
| GWFAM227 | B0 (13) | 0 | 0 |
| GWFAM26 | B2 (13), B3 (9) | 0 | 0 |
| GWFAM382 | B0 (4) | 0 | 0 |
| GWFAM390 | B1 (20) | 1 | 0 |
| GWFAM4 | B0 (5) | 0 | 0 |
| GWFAM70 | B1 (16), B4 (24) | 0 | 0 |
| GWFAM256 | none | 0 | 1 |
| GWFAM360 | none | 1 | 0 |

- O3 flags 2 of the 8 families with an expressed haplotype-only locus (GWFAM205, GWFAM390); 2 of the 26 families without one are
  flagged (GWFAM256 putative, GWFAM360 certified). Below MIN_SUPPORT (6 reads) a locus cannot be flagged: GWFAM382, GWFAM4, GWFAM175_B1.

**O2 on tied reads:**
- From reference copies: 432 tied, **157 placed, 0 wrong**; minimap2 MAPQ >= 1 418 placed, 29 wrong (6.9%).
- From B-only loci: 18 tied, 2 placed, **2 wrong** (GWFAM175, onto copy 1: the silent free-slot case of P14, now in real data);
  minimap2 MAPQ >= 1 16 placed, 16 wrong.
- 308 tied reads have no truth (tied between haplotype loci).

## Addendum: O3 with the port's IsoCon side (the stage, `candidate_copies`, run on the 34 families; work dir `.../o2refabsent/h`)
The run above had no IsoCon side (no stage products). With it (`o3port_p4.sh` stage / concat / nets / port c cc): 10 certified
candidates + 1 putative. Each contig was aligned to KB3781's two haplotypes and to `_pri`; a contig the haplotypes hold at ~1.0
but `_pri` only far lower is a copy the reference lacks (the truth table lists only loci found from the families' copies, so it
misses some):

| contig | best haplotype identity | best `_pri` identity | reading |
|---|---|---|---|
| cand_GWFAM4_1 | 0.9978 | 0.8673 | **reference-absent copy** |
| cand_GWFAM103_1 | 0.9993 | 0.8980 | **reference-absent copy** |
| cand_GWFAM175_0 | 1.0000 | 0.9378 | **reference-absent copy** (pat chr5 40.10-40.13 Mb, in the tandem array next to B1) |
| cand_GWFAM70_0 | 0.9996 | 0.9750 | **reference-absent copy** |
| putative_GWFAM256_0 | 0.9996 | 0.9908 | allele-level gap (0.9%): undecided |
| cand_GWFAM205_0 | 0.9996 | 0.9958 | allele-level (0.4%) |
| cand_GWFAM4_0 | 0.9905 | 0.9865 | allele-level / noisy consensus |
| cand_GWFAM360_0 | 1.0000 | 0.9971 | allele-level (0.3%) |
| cand_GWFAM205_1 | 1.0000 (70% covered) | 0.9988 | not absent |
| cand_GWFAM169_0 | 0.9974 | 0.9974 | not absent (false flag) |
| cand_GWFAM103_0 | 0.7148 | 0.7148 | junk contig (false flag) |

- **4 clearly reference-absent copies found, in 4 families** (GWFAM4, 103, 175, 70), each held by a haplotype at >= 0.998 and by
  `_pri` at <= 0.975. Only GWFAM175 was in the truth table's beyond-delta expressed list; the other three were unlisted.
- 4 flags are at an allele-level gap (0.3-0.9%) and cannot be called absent or allelic from this check; 3 are false flags.
- The IsoCon side is what finds them: without it the port flagged GWFAM205 / GWFAM390 / GWFAM360 / GWFAM256 and none of the four.
