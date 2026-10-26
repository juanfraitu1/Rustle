# StringTie-gap items — dev measurements (2026-10-25, branch assembler-stringtie-gaps)

Dev only, no held-out run; nothing adopted. `polish_port.py` is a Python port of `polish_drop_set`
that reproduces the shipped polish exactly on gorilla dev (3,862 kept, same ids; drops ISM 2249 / mono 649 /
frac 446 / ret 201 as logged). Inputs: `/mnt/linuxdisk/tmp/asm_gap/chr20_none.gtf` (GGO NC_073244.2, GGO_mm.bam)
and `h20_du_nopolish.gtf` (A119b chr20, deep); gffcompare vs RefSeq. Outputs `/mnt/linuxdisk/tmp/stgaps/`.

## Item 1 — StringTie-style isoform fraction (rlink.cpp:20058-20120): local reference + ISM credit

| substrate | variant | kept | chains = | chain Sn/Pr | tx Sn/Pr |
|---|---|---|---|---|---|
| GGO dev | shipped (gene max) | 3862 | 1571 | 28.3/40.8 | 26.5/40.7 |
| | local (exon-overlap max) | 3864 | 1571 | 28.3/40.8 | 26.5/40.7 |
| | ISM fragment reads credited to container | 3771 | 1561 | 28.1/41.5 | 26.3/41.4 |
| | local + credit | 3773 | 1561 | 28.1/41.5 | 26.3/41.4 |
| | no fraction filter | 4233 | 1589 | 28.6/37.6 | 26.8/37.6 |
| A119b chr20 | shipped | 6068 | 1096 | 25.6/19.5 | 24.0/18.1 |
| | local | 6120 | 1100 | 25.7/19.4 | 24.1/18.0 |
| | credit | 5270 | 1043 | 24.3/21.7 | 22.9/19.8 |
| | local + credit | 5295 | 1050 | 24.5/21.7 | 23.0/19.8 |
| | no fraction filter | 7494 | 1153 | 26.9/16.4 | 25.3/15.4 |

Only 13/446 (GGO) and 181/1727 (A119b) fraction drops share no exonic bp with their gene's best transcript.
**Verdict: locality is inert (≤ +4 chains); credit is a Sn↔Pr trade along the same frontier as the fraction
value itself (rows 859/864/873/1083). Not worth a code change.**

## Item 2 — "exon coverage share" rule (rlink.cpp:19895)

The audit mis-described it: the rule sits in the **opposite-strand** branch of the pairwise prediction filter —
it kills a prediction overlapping an antisense prediction when none of its exons carries > ½ of the local
(both-strand) coverage. It is an antisense-artifact filter. Ceiling here: shipped output has class `x`
(antisense exonic) 23/3862 (GGO, 0.6%) and 165/6068 (A119b, 2.7%) — removing ALL of them moves tx precision
by ≤ 0.5 points on A119b, ~0.2 on GGO. **Low priority; not built.**

The precision mass is class `j` (shares ≥ 1 junction, novel chain): 1727 GGO / 3202 A119b of the shipped
transcripts.

## Item 3 — poly-A / internal priming (rlink.cpp:486 `genome_polyA_terminated`, applied at :605)

StringTie's rule applies to **unspliced** long reads only (read-sequence test: ≥ 10 A in the last ≤ 25 aligned bases at
≥ 80%, or the T mirror at the start); it never touches multi-exon 3′ ends. (The 3′-end side already exists as
`--polish-tes pas-end`, r1133/r1136, opt-in.)
- Exact port on dev A119b chr20: flags **4 / 458** single-exon transcripts — inert.
- Genomic form (≥ 60% A / A6 in 20 bp past the aligned end, or T mirror before the start; r1064/r1135 instrument):
  246 / 458 at primed share ≥ 0.8, 0 `=`; read-level primed 0.739 vs 0.066-0.081 with ends shifted ±1 kb / +5 kb.
- **Held-out PASS** (`PREREG_mono_primed_2026-10-26.md`, sha1 c87260b8; scorer `score_mono_primed.py`):

| chr | single-exon (=) | flagged / NULL | `=` lost | tx Pr BASE → DROP | chains |
|---|---|---|---|---|---|
| A119b chr11 | 1448 (1) | 838 / 42 | 0 | 19.0 → 20.5 | identical |
| A119b chr12 | 1525 (2) | 852 / 32 | 1 | 16.7 → 17.9 | identical |

- **Wired opt-in**: `copy_assign --polish-mono-primed` (driver `RUSTLE_POLISH_MONO_PRIMED=1`), last step of each
  contig's polish. Flag off: chr11/chr12 `.gtf`, `.params.tsv`, `.families.gtf` byte-identical to the frozen main
  binary (sha1 685801dd). Flag on: dropped set == the Python scorer's on both chromosomes. Unit test
  `mono_primed_read_checks_a_past_the_end_or_t_before_the_start`. Gorilla OR6737 dev: ~10 single-exon transcripts,
  inert. Default flip = user's call.

## Item 4 — coverage-drop TSS/CPAS splitting (`longtrim`): not attempted (ends only; StringTie's lead is not ends).

## Junction correction — StringTie long-read rule (ST) and sequence-checked form (SEQ)

`PREREG_junction_correct_2026-10-26.md` (written before any run); `junction_correct.py` rewrites the BAM's CIGARs, the frozen
assembler (sha1 685801dd) runs unchanged; `none` round trip byte-identical to BASE (gorilla dev). Outputs `/mnt/linuxdisk/tmp/stgaps/jc/`.
In `-L` mode StringTie marks every read "mismatched" (rlink.cpp:735), so every non-guide junction is corrected: each site
independently to the best-supported same-strand site within 25 bp (support x 0.9 > own).

| substrate | reads moved ST / SEQ | chains BASE / ST / SEQ | chain Sn/Pr BASE | ST | SEQ |
|---|---|---|---|---|---|
| GGO NC_073244.2 (dev) | 9,765 / 22 | 1571 / 1414 / 1572 | 28.3/40.8 | 25.5/39.9 | 28.3/40.8 |
| A119b chr20 (dev) | 20,397 / 163 | 1059 / 972 / 1059 | 24.7/21.0 | 22.7/20.5 | 24.7/21.0 |
| **A119b chr11** | 35,125 / 160 | 2089 / 1900 / 2090 | 21.3/21.8 | 19.4/21.0 | 21.3/21.8 |
| **A119b chr12** | 49,136 / 109 | 2033 / 1878 / 2033 | 21.7/19.1 | 20.1/18.8 | 21.7/19.1 |
| **GGO NC_073230.2** | 9,663 / 9 | 1148 / 1048 / 1148 | 24.5/35.1 | 22.4/34.3 | 24.5/35.1 |

**Verdict: both FAIL; not built.** ST fails G and N on 3/3 held-out (−7.6 to −9.0% chains, precision also down 0.3-0.8).
SEQ fails G on 3/3 (inert): 89-99% of ST's candidate moves are rejected because the read's own bases fit its original
site better (cf. r1140). Register rows to add at merge time.

## Dominance prune (PREREG_dominance_prune_2026-10-26) — human A119b dev run (not the verdict)
Outputs `/mnt/linuxdisk/tmp/stgaps/fw/hsa_dp/` (support 56,249 molecules; minimap2 map-hifi -N 100 -p 0.8 to 217,580
pooled multi-exon transcripts, 420 s / 3.2 GB). Dropped 5,668 / 12,184 pool-added candidates.

| arm | single R / P | multi R / P | multi `=` |
|---|---|---|---|
| PRIM (primary-only) | 22.1 / 19.8 | 17.7 / 14.9 | 1,645 |
| POOL (shipped) | 22.1 / 19.8 | 20.7 / 12.6 | 1,924 |
| PRUNE | 22.1 / 19.8 | 19.3 / 13.6 | 1,793 |

D1 FAIL (13.6 < 14.9), D2 FAIL (148 kept of 279 gained; needed 223), D3 PASS. Of the 131 dropped `=` transcripts, 113
have every support read scoring 1-19 AS points better on another copy's transcript, and 171 support pairs have no hit on
the candidate at all. The reads come from another copy that has the same intron chain. gffcompare `=` cannot tell a
read-observed copy from a structure the pool placed at the wrong paralog. So D2's "gain" partly rewards reference-shaped
projections. The gorilla held-out was not run: dev fails both precision and gain clauses.

## `--polish-mono-primed` default ON (user decision 2026-10-26)
`copy_assign` now applies it by default; `--no-polish-mono-primed` (driver `RUSTLE_POLISH_MONO_PRIMED=0`) restores the old
output. Checked on A119b chr12: default == the earlier `--polish-mono-primed` run and `--no-polish-mono-primed` == the
frozen main binary, byte for byte (gtf, params.tsv, families.gtf).

## Family layer FW (PREREG_family_layer_2026-10-26, sha1 718822c1) — FAIL on held-out (F1)
`bench/stgaps/fw_layer.py`; outputs `/mnt/linuxdisk/tmp/stgaps/fwl/{hsa,ggo}/`. BASE and MIN1 (`RUSTLE_MIN_CHAIN_READS=1`)
from the same binary (this branch, `--polish-mono-primed` default). FW = BASE + single-read chains at an O1 family copy
whose lifted family chain has >= 2 molecules over the family's copies.

| substrate | candidates | confirmed (admitted) | unconfirmed | unlifted | own 0 |
|---|---|---|---|---|---|
| human A119b (dev) | 22,247 | 601 | 7,846 | 8,165 | 5,634 |
| gorilla KB3781 (held-out; 326 families, 918 copies) | 2,193 | 75 | 1,044 | 769 | 305 |

| substrate | arm | single R / P | multi R / P | multi `=` |
|---|---|---|---|---|
| human (dev) | BASE | 22.1 / 19.8 | 20.7 / 12.6 | 1,924 |
| | MIN1 | 24.0 / 14.0 | 22.7 / 9.1 | 2,112 |
| | FW | 22.1 / 19.7 | 20.8 / 12.5 | 1,932 |
| gorilla (held-out) | BASE | 28.1 / 37.2 | 13.4 / 26.5 | 951 |
| | MIN1 | 29.9 / 29.2 | 14.8 / 19.3 | 1,052 |
| | FW | 28.1 / 37.2 | 13.5 / 26.3 | 956 |

`=` share among own-support-1 candidates (MIN1 tmap): human confirmed 10/601 = 1.7% vs unconfirmed 141/7,846 = 1.8%;
gorilla confirmed 10/75 = 13.3% vs unconfirmed 25/1,044 = 2.4%.
Held-out clauses: F1 FAIL (+0.1 point, 5 transcripts; bar 0.5), F2 PASS (−0.2), F3 PASS (5.5x), F4 PASS (every BASE line
kept; single-copy unchanged). Dev: F1 FAIL, F3 FAIL (no signal). **Not adopted.** Even admitting EVERY single-read chain
(MIN1) gains only 2.0 / 1.4 multi-copy recall points at a 3.5 / 7.2 point precision cost, so single-read chains are not
where the multi-copy gap is. About half of the own-support-1 candidates have a junction that does not lift to the
representative (indel at the splice site or partial copy alignment).

## Recall on observed reference transcripts (descriptive; `bench/stgaps/observed.py`)
Observed at k = >= k de-duplicated primary spliced molecules carry the reference transcript's exact intron chain on its
contig (one pass over each BAM: gorilla KB3781 12.6M spliced primaries / 324 s, human A119b 12.4M / 485 s).

| substrate | stratum | all: n, BASE R | obs>=1: n, BASE / MIN1 R | obs>=2: n, BASE / MIN1 R |
|---|---|---|---|---|
| gorilla KB3781 | single | 88,740, 28.1 | 33,756, 73.8 / 78.4 | 27,306, 91.2 / 91.1 |
| | multi | 7,093, 13.5 | 1,362 (19%), 66.5 / 72.8 | 1,025, 86.0 / 86.0 |
| human A119b | single | 166,114, 22.1 | 55,594, 66.0 / 71.6 | 42,178, 87.0 / 86.9 |
| | multi | 9,288, 20.7 | 2,521 (27%), 69.0 / 74.2 | 1,875, 87.3 / 87.2 |

Human obs>=2, other tools: StringTie single 64.7 / multi 65.7, FLAIR 77.7 / 74.9 (ours 87.0 / 87.3). Multi-copy
transcripts whose molecules are >= half MAPQ 0: human 209 (ours 92.8, StringTie 72.7, FLAIR 60.3), gorilla 45 (91.1).
**Reading.** The single- vs multi-copy recall gap is mostly a gap in what the reads show: only 19% (gorilla) / 27% (human)
of multi-copy reference transcripts have any primary read carrying their whole chain at their own locus, against 38% /
33% single-copy. Once observed twice, ours recovers multi-copy transcripts as well as single-copy ones in human (87.3
vs 87.0), and 5 points lower in gorilla (86.0 vs 91.2, 143 of 1,025 missed). Single-read chains (obs = 1) are the
remaining gap, and MIN1 shows admitting them all costs far more precision than it recovers (FW section). ⚠ "Own locus"
means the reads' PRIMARY placement: a multi-copy transcript whose reads went to a paralog counts as unobserved here.
That is the O2 question, not an assembly loss.
