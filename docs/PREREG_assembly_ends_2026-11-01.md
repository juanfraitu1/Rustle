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
