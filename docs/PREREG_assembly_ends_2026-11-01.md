# PREREG ASSEMBLY ENDS (2026-11-01, register 1485 / 1486)

**Status:** written before any held-out run and before any code for E2. Code at the commit that adds this file.
**Question (user, 2026-11-01):** after the single-copy SQANTI / gffcompare panel (register 1482), emit more complete transcripts.
**Diagnosis it rests on** (chr20 panel, 370 genes, descriptive, `~/asm_chr20_diag/`): of 222 reference isoforms another tool matches
exactly and we do not, 80 are never assembled (68 with <= 1 read carrying the exact chain) and 142 are assembled then polished away;
the sub-chain (ISM) collapse alone drops 67 (median 16 coordinate-unique exact-chain reads). Keeping dropped sub-chains by read count
is not selective (5-8% of the added are reference matches at any floor 3-50, shipped precision 19.1%; registers 858, 1063, 1066).

**Substrate.** A119b on CHM13 (`A119b.t2t.bam`), `copy_assign --assemble-only --region CHR:0-LEN`, seeding exactly as the base run
(good pool, best-AS table `~/runs/human_A119b/a1.molecules.tsv`), shipped polish = `tools/rustle_pipeline.sh` POLISH (full,
fraction 0.02, mono shadow, mono quantile 0.82, ISM ratio 0.7, retained ratio 10). Truth: CAT/Liftoff (`chm13v2.0_CAT_Liftoff.slim.gff3.gz`
-> `gffread -T`, per contig), whole contig, `gffcompare -r`. Metrics: matching intron chains (count) and transcript-level precision.
Development contig chr20 (already seen); held-out contigs chr11 and chr7. Disclosure: register 1139 ran `--polish-tss rescue`
genome-wide on A119b descriptively (+1,711 transcripts / +94 chains, precision not reported), so the held-out contigs are not unseen
in aggregate for E1.

## E1 (register 1485): `--polish-tss rescue` as the assembly default
Arms: SHIP vs SHIP + `--polish-tss rescue`. Run once, nothing tuned.
- **PASS** iff on BOTH held-out contigs: matching chains >= SHIP's AND transcript precision >= SHIP's - 0.1 points.
- PASS -> propose it as the pipeline default (the user decides); FAIL -> stays opt-in.
Reported beside: chr20 numbers (dev), transcripts added, the added set's own exact-match rate.

## E2 (register 1486): rescue dropped sub-chains by end evidence
Idea: a shorter isoform y inside a longer emitted x is a real alternative-end isoform when each end where y differs from x is proven
by the run's own evidence -- a different 5' end by a read-proven TSS cluster (`--polish-tss`'s proof), a different 3' end by a
PAS-proven TES cluster (`--polish-tes`'s proof) -- and a 5'-truncation fragment when not. `--polish-tss rescue` already covers
y that are 3'-flush with x; E2 adds y whose 3' end differs (3'-short and both-short).
**Development (chr20 only):** measure, with `--polish-subchain tag`, `--polish-tss tag` and `--polish-tes tag` on the arm without the
ISM collapse, which side the dropped sub-chains miss and how often each class matches the reference with and without a proven end.
Then implement the rule (test-first) as an opt-in mode and FREEZE its exact definition in an amendment to this file **before** any
held-out run. If development shows no class whose proven-end subset beats half the shipped precision, E2 stops there (reported).
**Held-out bars (chr11 and chr7, both), arm = SHIP + E2 (and + E1 if E1 passed):**
- B1: matching chains >= SHIP + 1% of SHIP's chains;
- B2: transcript precision >= SHIP's - 0.3 points;
- B3: the transcripts E2 adds match the reference at >= half of SHIP's transcript precision.
PASS (all three on both) -> propose E2 as a default (the user decides); FAIL -> opt-in or removed.

## RESULT (2026-11-01, job 4405763, `~/asm_ends/`)
**E1 PASS on both held-out contigs** (gffcompare vs CAT, whole contig; exact transcript precision = matching transcripts / query):

| contig | arm | transcripts | matching chains | transcript precision |
|---|---|---|---|---|
| chr11 (held-out) | SHIP | 10,231 | 1,695 | 16.577 |
| chr11 (held-out) | + tss rescue | 10,318 | 1,702 | 16.505 (-0.072) |
| chr7 (held-out) | SHIP | 11,675 | 1,553 | 13.310 |
| chr7 (held-out) | + tss rescue | 11,749 | 1,562 | 13.303 (-0.007) |
| chr20 (dev) | SHIP | 5,319 | 876 | 16.507 |
| chr20 (dev) | + tss rescue | 5,366 | 882 | 16.474 |

Added sets match the reference at 8% (chr11, 7 of 87) and 12% (chr7, 9 of 74). The user had asked for the default on a pass:
`tools/rustle_pipeline.sh` now runs `--polish-tss rescue` unless `RUSTLE_POLISH_TSS` says otherwise (`off` gives the earlier output).
**E2 STOPPED at development** (chr20, the arm without the ISM collapse, `--polish-subchain/--polish-tss/--polish-tes tag`): the 3,024
sub-chains the collapse drops match the reference at 4.4% (shipped 16.5%). The class E2 targeted, 3'-short sub-chains, matches at
2.1% (97; 2.7% with a TSS proof, 2.1% without a TES proof; no 3'-short one has a PAS-proven end), so 3' end evidence cannot separate
it. The best proven-end subsets are sub-chains the end-compatible tag does not cover: + TES-proven 10.8% (622, +67 chains) and
+ TES- and TSS-proven 11.8% (68, +8). Projected on chr20 each fails a held-out bar already on development: the first B2 (precision
15.91 < 16.21), the second B1 (884 < 885 chains), the 3'-short class B1 and B3. No E2 code was written; no held-out run was made.
Disclosure: the stop condition as registered (no class above half the shipped precision) was not met by the two untagged subsets;
E2 was stopped on the projection against the held-out bars instead, which the file did not name as a stop rule.

## E3 (register 1487): read-consensus 5' ends (IsoCon / FLAIR's lesson), frozen 2026-11-01 before any code or held-out run
**Diagnosis it rests on** (chr20 panel, descriptive): on 628 multi-exon FSM transcripts with >= 2 exact-chain reads, our emitted 5'
end IS the most upstream read start (identical distributions), median 63 bp from the reference TSS, 43% within 50 bp; the median of the
exact-chain reads' starts would give 24 bp / 70%; FLAIR (one model per isoform) 21 bp, IsoCon 9 bp (best of ~4 candidates).
**Rule (new opt-in flag `--polish-tss-end off|median`, needs the `--polish-tss` evidence):** on each contig's final transcripts (after
`--polish-tss`), for every multi-exon transcript NOT rescued or split by `--polish-tss`, its own starts = the deduplicated primary
spliced records whose first intron is the transcript's first intron, whose intron chain equals the transcript's, and whose oriented
5' end lies inside the transcript's oriented first exon. With >= 2 own starts, the 5' end moves to their lower median
(sorted[(n-1)/2], oriented coordinates). Nothing else changes (chain, 3' end, `reads`, transcript set); the locus boundaries the
families stage uses are untouched (the step runs on the emitted GTF only).
**Arms:** SHIP (the pipeline default since r1485, i.e. with `--polish-tss rescue`) vs SHIP + `--polish-tss-end median`.
**Metric:** gffcompare vs CAT per contig; over query transcripts with class `=` and >= 2 exons, M = share whose 5' end lies within
50 bp of the matched reference transcript's 5' end.
**Bars (held-out chr11 and chr7, both; run once):** B1 M(E3) >= M(SHIP) + 10 points; B2 matching intron chains and matching transcripts
identical to SHIP (by construction: any difference is a bug and fails); B3 base-level sensitivity >= SHIP - 2.0 points.
PASS -> propose as the pipeline default (the user decides). Reported beside: chr20, the 3' ends (unchanged by construction),
how many transcripts moved and by how much.

### E3 RESULT (2026-11-01, job 4405816, `~/asm_ends/`, scorer `bench/assembly_ends/end_score.py`): PASS on both held-out contigs

| contig | arm | `=` multi-exon | 5' within 50 bp (M) | median 5' distance | base sensitivity | chains | matching transcripts |
|---|---|---|---|---|---|---|---|
| chr11 (held-out) | SHIP | 1,702 | 40.60% | 77 bp | 45.5 | 1,702 | 1,703 |
| chr11 (held-out) | + median | 1,702 | **64.57%** | 28 bp | 44.0 | 1,702 | 1,703 |
| chr7 (held-out) | SHIP | 1,562 | 39.31% | 79 bp | 50.6 | 1,562 | 1,563 |
| chr7 (held-out) | + median | 1,562 | **63.44%** | 30 bp | 49.3 | 1,562 | 1,563 |
| chr20 (dev) | SHIP | 882 | 46.15% | 57 bp | 53.1 | 882 | 884 |
| chr20 (dev) | + median | 882 | 67.91% | 27 bp | 51.8 | 882 | 884 |

B1 +24.0 / +24.1 points (bar +10) PASS; B2 chains and matching transcripts identical PASS; B3 base sensitivity -1.5 / -1.3 (bar -2.0)
PASS. Base precision rises (25.8 -> 26.2, 18.7 -> 19.1). Moved: 6,081 (chr11), 6,478 (chr7), 3,174 (chr20) transcripts; the rescued /
split forms (87 / 74 / 47) stay at their proven mode. Checked: same transcript set, no 3' end and no internal exon changed on any
contig; without the flag the chr20 output is byte-identical to the r1485 arm. The 3' share moves by 0.06 on chr11 only because
gffcompare picks a different representative query for a reference once starts move. Status: flag `--polish-tss-end median` (opt-in,
driver `RUSTLE_POLISH_TSS_END=median`); PASS proposes it as the pipeline default, the user decides. Tests: 3 new (`tss_end_median_*`),
copy_assign 89/89.

## E4 (register 1488): 3' UTR extension by unspliced coverage to a polyA-proven end, frozen 2026-11-01 before any code or run
**Diagnosis it rests on** (chr20 panel, descriptive): FLAIR covers 292k reference exonic bases we miss (16% of the reference), 54% of
them in annotated LAST exons; 62% of those lie past our model's 3' end inside the annotated last exon (median read depth there 16); in
the ten largest cases (3.5-8 kb short) 340 of 399 reads ending > 50 bp past our end are UNSPLICED (5'-truncated molecules that never
reach the last splice junction, so pass 1 never joins them to the spliced model).
**Rule (new opt-in flag `--polish-tes-extend off|on`; reads the `--polish-tss` evidence, plus a new stream of deduplicated primary
UNSPLICED records, strand = the read's orientation as for spliced records):** on each contig's final transcripts (after
`--polish-tss`, `--polish-tss-end`, `--polish-tes`), for every multi-exon transcript T with oriented 3' end E:
1. the unspliced records on T's strand cluster their oriented 3' ends by single linkage (gap > 21 bp starts a cluster, >= 2 reads, mode
   = most-ended position ties 3'-most), exactly `--polish-tes`'s cluster rule; a cluster is PAS-proven by `--polish-tes`'s test
   (canonical AATAAA / ATTAAA in [mode-35, mode-10], mode not internally primed);
2. candidate ends: PAS-proven cluster modes M > E (oriented) with M - E <= 20,000 bp;
3. M is admissible when every base of (E, M] is covered by >= 2 unspliced records on T's strand, and (E, M] overlaps no exon of a
   transcript of ANOTHER gene_id on the same strand (read-through guard);
4. T's 3' end moves to the most distal admissible M (none: unchanged). Nothing else changes (chain, 5' end, `reads`, the set).
**Arms:** SHIP (pipeline default, `--polish-tss rescue`) vs SHIP + `--polish-tes-extend on`. Reported beside: + `--polish-tss-end median`.
**Bars (held-out chr11 and chr7, both; run once):** B1 base-level sensitivity >= SHIP + 2.0 points; B2 base-level precision >= SHIP
- 1.0 point; B3 matching intron chains and matching transcripts identical (by construction); B4 among `=` multi-exon transcripts, the
share with the 3' end within 50 bp of the reference >= SHIP - 3 points.
PASS -> propose as default (the user decides).
