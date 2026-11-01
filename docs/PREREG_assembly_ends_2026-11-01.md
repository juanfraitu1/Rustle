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
(sorted[(n-1)/2], oriented coordinates). Nothing else changes (chain, 3' end, `reads`, transcript set). [CORRECTION 2026-11-01: the families stage reads
`PREFIX.families.gtf`, written after every polish step, so this layer DOES change the families' input ends; its effect on families is
untested.]
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

### E4 RESULT (2026-11-01, job 4405941, `~/asm_ends/`): PASS on both held-out contigs

| contig | arm | base Sn / Pr | chains | matching tx | `=` 3' within 50 bp | median 3' | 5' within 50 bp |
|---|---|---|---|---|---|---|---|
| chr11 (held-out) | SHIP | 45.5 / 25.8 | 1,702 | 1,703 | 49.00% | 62 bp | 40.60% |
| chr11 (held-out) | + extend | **47.6** / 25.9 | 1,702 | 1,703 | 48.00% | 79.5 bp | 40.60% |
| chr11 (held-out) | + extend + median 5' | 46.1 / 26.2 | 1,702 | 1,703 | 48.06% | 78 bp | 64.57% |
| chr7 (held-out) | SHIP | 50.6 / 18.7 | 1,562 | 1,563 | 43.73% | 171.5 bp | 39.31% |
| chr7 (held-out) | + extend | **53.0** / 18.7 | 1,562 | 1,563 | 42.19% | 226 bp | 39.31% |
| chr7 (held-out) | + extend + median 5' | 51.8 / 19.0 | 1,562 | 1,563 | 42.19% | 226 bp | 63.44% |
| chr20 (dev) | SHIP | 53.1 / 26.9 | 882 | 884 | 60.66% | 10 bp | 46.15% |
| chr20 (dev) | + extend | 54.4 / 26.5 | 882 | 884 | 60.66% | 9.5 bp | 46.15% |

B1 base sensitivity +2.1 / +2.4 (bar +2.0) PASS; B2 base precision +0.1 / 0.0 (bar -1.0) PASS; B3 chains and matching transcripts
identical PASS; B4 3' within 50 bp -1.0 / -1.5 points (bar -3) PASS. Moved: 2,777 (chr11), 3,273 (chr7), 1,479 (chr20) transcripts
from 300,767 / 398,711 / 120,518 unspliced records. Without the flag the chr20 output is byte-identical to the r1485 arm. The two layers
compose: with the median 5' layer base sensitivity sits between (46.1 / 51.8), precision is highest, both end shares keep their gains.
Status: opt-in `--polish-tes-extend on` (driver `RUSTLE_POLISH_TES_EXTEND=on`); PASS proposes it as default, the user decides. Tests: 5
new (`tes_extend_*`, `tss_evidence_keeps_deduplicated_oriented_unspliced_records`), copy_assign 94/94; the coverage floor was mutated to
1 to confirm a test fails.

## E5 (register 1489): single-exon models from end-clustered unspliced reads, frozen 2026-11-01 before any code or held-out run
**Diagnosis it rests on** (chr20, descriptive, `~/asm_chr20_diag/`): pass 1 builds a single-exon model as the union of overlapping
unspliced reads; where one overlaps an expressed CAT single-exon transcript it is a median 17x longer; the raw (unpolished) 5,394
single-exon models match 4 of chr20's 600 single-exon CAT transcripts exactly, the shipped 271 match 2. The Python prototype of the rule
below (W 50, 3 reads) emits 366 models, 37 exact (10%); swapped into SHIP chr20: base precision 26.9 -> 40.8, base sensitivity 53.1 ->
54.5, transcript precision 16.5 -> 16.4, matching transcripts 878 -> 888. chr20 is development; the rule is the prototype's, untuned.
**Rule (new opt-in flag `--mono-model off|cluster`; reads the `--polish-tes-extend` unspliced evidence stream, so it needs
`--polish-tss` on):** as the LAST step of each contig's final transcripts (after `--polish-tss`, `--polish-tss-end`, `--polish-tes`,
`--polish-tes-extend`, `--polish-mono-primed`):
1. every single-exon transcript of the contig is removed;
2. the contig's deduplicated primary unspliced records, per strand (the read's orientation), sorted by (start, end), are clustered
   greedily: the first unused record seeds a cluster of every unused record whose start lies within 50 bp after the seed's start and
   whose end lies within 50 bp of the seed's end;
3. a cluster of >= 3 records gives a model [lower median start, lower median end] (genomic, 1-based);
4. the model is dropped when it overlaps an exon of a multi-exon transcript of the contig on the same strand, when its 3' end is
   internally primed (`--polish-tes`'s test), or when its 3' end has no canonical PAS in [3' - 35, 3' - 10] (oriented);
5. each kept model is emitted as its own gene (`gene_id` / `transcript_id` from contig, strand and start) with `reads` = cluster size.
Multi-exon transcripts are untouched.
**Arms:** SHIP (pipeline default, `--polish-tss rescue`) vs SHIP + `--mono-model cluster` (the flag turns the unspliced stream on).
**Bars (held-out chr11 and chr7, both; run once):** B1 base precision >= SHIP + 5.0 points; B2 base sensitivity >= SHIP; B3 transcript
precision >= SHIP - 0.1 points; B4 matching intron chains identical (multi-exon models untouched). Reported beside: exact matches to
CAT single-exon transcripts, models emitted, chr20.
PASS -> propose as default (the user decides).

### E5 RESULT (2026-11-01, job 4406074, `~/asm_ends/`): PASS on both held-out contigs

| contig | arm | transcripts | base Sn / Pr | transcript Sn / Pr (exact) | chains | matching transcripts | single-exon removed / emitted |
|---|---|---|---|---|---|---|---|
| chr11 (held-out) | SHIP | 10,318 | 45.5 / 25.8 | 12.5 / 16.505 | 1,702 | 1,703 | |
| chr11 (held-out) | + mono cluster | 10,310 | 45.8 / **43.0** | 12.7 / 16.702 | 1,702 | 1,722 | 624 / 616 |
| chr7 (held-out) | SHIP | 11,749 | 50.6 / 18.7 | 14.2 / 13.303 | 1,562 | 1,563 | |
| chr7 (held-out) | + mono cluster | 11,513 | 50.7 / **35.2** | 14.3 / 13.715 | 1,562 | 1,579 | 971 / 735 |
| chr20 (dev) | SHIP | 5,366 | 53.1 / 26.9 | 16.2 / 16.474 | 882 | 884 | |
| chr20 (dev) | + mono cluster | 5,463 | 54.5 / 40.7 | 16.3 / 16.365 | 882 | 894 | 271 / 368 |

B1 base precision +17.2 / +16.5 points (bar +5.0) PASS; B2 base sensitivity +0.3 / +0.1 PASS; B3 transcript precision +0.20 / +0.41
(bar -0.1) PASS; B4 chains identical PASS. Matching transcripts +19 / +16. Without the flag the chr20 output is byte-identical to the r1485
arm. Disclosure: the first run of this job (4406056) crashed on every contig before writing a GTF (two clusters with one median start
got one transcript id; `--bridge-regroup` refused the duplicate). No score was produced or read; the id now carries start and end plus
a counter (test `mono_model_gives_clusters_with_one_start_distinct_ids`, watched failing first); rule and bars unchanged. Status: opt-in
`--mono-model cluster` (driver `RUSTLE_MONO_MODEL=cluster`); PASS proposes it as default, the user decides. Tests: 4 new, copy_assign 98/98.

### Combined E3 + E4 + E5 (descriptive, 2026-11-01, job 4406122): the three opt-in layers together vs SHIP

| contig | arm | transcripts | base Sn / Pr | transcript precision | matching transcripts | 5' within 50 bp | 3' within 50 bp |
|---|---|---|---|---|---|---|---|
| chr11 | SHIP | 10,318 | 45.5 / 25.8 | 16.505 | 1,703 | 40.60% | 49.00% |
| chr11 | all three | 10,065 | 45.8 / 43.6 | **17.119** | 1,723 | 64.57% | 48.06% |
| chr7 | SHIP | 11,749 | 50.6 / 18.7 | 13.303 | 1,563 | 39.31% | 43.73% |
| chr7 | all three | 11,271 | 51.2 / 35.2 | **14.001** | 1,578 | 63.44% | 42.19% |
| chr20 | SHIP | 5,366 | 53.1 / 26.9 | 16.474 | 884 | 46.15% | 60.66% |
| chr20 | all three | 5,400 | 54.3 / 41.3 | 16.556 | 894 | 67.91% | 60.66% |

Chains identical on every contig. The layers interact in one place: extended 3' ends overlap more unspliced clusters, so the single-exon
step keeps fewer models (chr11 10,310 -> 10,065 transcripts), which lifts transcript precision above any single layer. Base
sensitivity with all three sits between the median-5' arm (which trims) and the extension arm (which adds). On the 370-gene single-copy
panel (`docs/PREREG_sqanti_single_copy_2026-11-01.md`) the combination reads base 71.9 / 58.5, transcript 22.8 / 18.6, 5' within 50 bp
69.8%: there, inside protein-coding gene spans, the added single-exon models rarely match, so transcript precision drops 0.5 points
(19.1 -> 18.6) while whole contigs gain 0.1-0.7.

**Correction (2026-11-01, before any default change):** E3, E4 and E5 all act on the GTF that `PREFIX.families.gtf` is written from
(after every polish step, minus the bridges), so they change the families stage's input: transcript ends (E3, E4) and the single-exon
models (E5). Intron chains are identical, but locus extents and single-exon loci move; register 975 / 979 showed locus extents matter
for family graph admission. None of the three has been scored on families; that check precedes any default.
