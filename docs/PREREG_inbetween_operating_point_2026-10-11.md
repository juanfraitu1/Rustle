# An operating point between StringTie and FLAIR that wins on multi-copy loci, 2026-10-11

**Pre-registered before any arm below is run genome-wide and before the strata are computed.** User (2026-10-11): "something
in between StringTie and FLAIR, that has the strengths of both, but ours should outperform them in secondary-only loci or in
places where O2 found the right assignment or corrected the arbitrarily selected primary when AS is equal." User choices
(same day): success = beat each tool on its weak side; multi-copy arm = the good-secondaries pool only (no O2 closed loop);
decisive samples = the benchmark samples genome-wide.

## The operating point (fixed from the difference model, register 1267; no search)

`OP` = `copy_assign --direct-transfrags-union --read-isoform-k 2` on top of the shipped assembly settings, i.e.
`RUSTLE_DIRECT_TRANSFRAGS=union RUSTLE_READ_ISOFORM_K=2 tools/rustle_pipeline.sh assemble ...`. It keeps all our chains (no
25 bp junction correction, floor 2: the classes register 1267 found right more often than the tools' additions), adds
StringTie's exclusive class (single-read chains whose every junction another read carries), and does not add FLAIR's
(`single_novel`). Two pools: `OP_primary` (`--no-seed-secondaries`) and `OP_good` (the default good pool: secondaries with
AS >= 0.98 x the molecule's genome-wide best AS). Reported only: our cached shipped default runs (`default_primary`,
`default_good`, figure run cache 2026-10-10).

## Samples and scope

Gorilla OR6737 (`GGO_mm.bam`) and human A119b (`A119b.t2t.bam`), genome-wide **minus the development chromosome**
(NC_073244.2 / chr20), against the lab's annotation-free StringTie 3.0.3 and FLAIR 3.0.1 GTFs on the same BAMs, scored by
gffcompare against `GGO.GCF_029281585.2_RefSeq.gtf.gz` / `A119b.chm13v2.0_RefSeq.gtf.gz`. Both samples were used before
(register 1264-1266; disclosed); no rule here is fitted on them.

## Multi-copy strata (from reads only; `bench/tool_difference/multicopy_strata.py`, sha1 c663dc6b)

Per reference gene, records whose aligned blocks overlap its exons: primaries; GOOD secondaries (AS >= 0.98 x the molecule's
best, from `as_table`'s `molecules.tsv`); TIED secondaries (AS == the molecule's best: its primary was one of several equally
scored placements). **S1 secondary-only**: 0 primaries and >= 2 good secondaries. **S2 coin-toss-starved**: >= 1 primary,
>= 2 tied secondaries and tied >= primaries. Evaluated: the multi-exon reference transcripts of S1 / S2 genes.

## Claims (`bench/tool_difference/score_op.py`, sha1 f65413d1; each in both species)

- **B1** `OP_primary` intron-chain sensitivity > StringTie's and > FLAIR's.
- **B2** `OP_primary` intron-chain precision > FLAIR's.
- **C1** In S1, `OP_good` matches exactly (gffcompare `=`) more reference transcripts than StringTie and than FLAIR.
- **C2** The same in S2.
- **C3** In each stratum, `OP_good`'s precision (its transcripts touching a stratum gene's exon that match a reference
  transcript exactly) >= FLAIR's genome-wide intron-chain precision (the weak-side bar of B2).

A stratum with < 20 multi-exon reference transcripts in a species is UNDERPOWERED there (no verdict for C1/C2/C3 on it).
**Verdict**: B1-B2 in both species => OP is "between StringTie and FLAIR" as defined; C1-C3 in both species => it also wins
on multi-copy loci. Each failure is named; nothing is re-tuned after a number is seen. Reported, no bar: base / intron /
transcript / locus levels, OP_primary vs StringTie precision, default arms, stratum sizes.

**Expected, written now.** B1/B2 pass in both (dev chr20: gorilla chain 31.3 / 34.1 vs StringTie 24.7 / 37.0 and FLAIR 25.1 /
23.6; human 26.7 / 17.2 vs 20.1 / 16.8 and 23.9 / 7.6). C1 should pass wherever S1 is powered, since neither tool's transcripts
rest on secondary-only coverage; C3 is the real test (do transcripts built from secondaries match the annotation?). S2 may be
small: register 875 found few AS-tied molecules in catalog families.

**Amendment 1 (2026-10-11, before any stratum or arm number).** (a) Instrument crash fix: RefSeq gene-level features carry
`gene_id ""` (gorilla and human references); both scripts now key such a gene by its `transcript_id` (`gene-LOC...`) instead of
crashing. No definition changes. (b) Execution: `--read-isoform-k` forces the buffered reader, which held 14+ GB genome-wide
with secondaries (the gorilla good-pool genome-wide run stalled; the human good-pool run was stopped by Claude Code's memory
reaper). The arms run one contig per process (`tools/rustle_pipeline.sh --contig`, documented equal to the genome-wide run on
that contig) and `--read-isoform-k` now also runs on the streaming reader (1304476b; structurally identical to the buffered
path on human chr20). Same arms, same flags.
(c) The strata script ran past the 500 s run limit on gorilla before printing any stratum count; it now counts each
chromosome in a separate worker (4 processes) with the identical per-record rule (sha1 below), before any stratum number.
multicopy_strata.py sha1 f4d7bc36798e6f89fc914b6d52fe7c3b1aa9cd78

## Result (2026-10-11): B1, B2, C1, C2, C3 PASS in both species

Scored with `bench/tool_difference/score_op.py` (sha1 45901172, the amendment-1 version); outputs
`/mnt/linuxdisk/tmp/asm_gap/op/score_{gorilla,human}.txt`. Genome-wide minus the development chromosome, sensitivity / precision:

| arm | gorilla base | intron | intron chain | transcript | locus | human base | intron | intron chain | transcript | locus |
|---|---|---|---|---|---|---|---|---|---|---|
| OP_primary | 48.6 / 63.7 | 65.7 / 78.1 | **29.3 / 30.4** | 26.3 / 30.2 | 31.7 / 70.1 | 50.6 / 18.2 | 66.7 / 56.4 | **23.6 / 15.6** | 22.0 / 14.0 | 36.1 / 23.5 |
| OP_good | 49.1 / 62.0 | 66.1 / 77.8 | 29.4 / 30.0 | 26.5 / 29.8 | 32.0 / 68.3 | 51.7 / 17.8 | 67.0 / 54.9 | 23.8 / 15.0 | 22.2 / 13.4 | 36.3 / 22.6 |
| StringTie | 50.3 / 62.6 | 64.2 / 83.1 | 24.4 / 34.4 | 21.9 / 34.1 | 30.7 / 64.1 | 59.6 / 14.4 | 63.5 / 57.8 | 18.1 / 15.0 | 16.9 / 12.7 | 33.7 / 20.1 |
| FLAIR | 48.7 / 45.2 | 62.5 / 61.3 | 25.5 / 17.3 | 23.0 / 15.7 | 28.4 / 48.4 | 66.4 / 8.8 | 67.7 / 30.1 | 22.4 / 6.2 | 21.1 / 3.6 | 31.9 / 6.3 |
| default_primary (reported) | 45.1 / 67.7 | 59.6 / 82.8 | 26.8 / 35.3 | 24.0 / 35.1 | 30.2 / 70.7 | 48.1 / 17.4 | 60.8 / 64.4 | 21.8 / 18.5 | 20.3 / 16.1 | 34.5 / 22.8 |
| default_good (reported) | 45.7 / 65.9 | 60.0 / 82.5 | 27.0 / 35.0 | 24.2 / 34.7 | 30.5 / 69.1 | 49.2 / 17.1 | 61.2 / 62.6 | 22.0 / 17.7 | 20.5 / 15.3 | 34.8 / 21.9 |

**B1 PASS** (chain sensitivity: gorilla 29.3 vs 24.4 / 25.5; human 23.6 vs 18.1 / 22.4). **B2 PASS** (chain precision: gorilla 30.4
vs FLAIR 17.3; human 15.6 vs 6.2). Reported: against StringTie's chain precision the OP is below in gorilla (30.4 vs 34.4) and above
in human (15.6 vs 15.0).

Multi-copy strata (reference multi-exon transcripts matched exactly; precision = matching / all arm transcripts touching the stratum):

| | gorilla S1 (542 tx, 305 genes) | gorilla S2 (252 tx, 117 genes) | human S1 (926 tx, 304 genes) | human S2 (592 tx, 266 genes) |
|---|---|---|---|---|
| OP_good | **34** (34 / 171 = 19.9%) | **84** (90 / 367 = 24.5%) | **36** (37 / 224 = 16.5%) | **210** (218 / 1,368 = 15.9%) |
| OP_primary | 0 | 56 (36.5%) | 0 | 100 (17.9%) |
| StringTie | 0 | 51 (39.0%) | 0 | 77 (13.3%) |
| FLAIR | 0 | 42 (36.8%) | 0 | 72 (8.0%) |
| default_good (reported) | 33 (28.7%) | 81 (32.8%) | 33 (19.3%) | 198 (18.7%) |

**C1 PASS** (S1: 34 / 36 vs 0 for both tools; neither tool emits a transcript on secondary-only genes). **C2 PASS** (S2: 84 vs 51 /
42; 210 vs 77 / 72). **C3 PASS** in every stratum (19.9 / 24.5 / 16.5 / 15.9% >= FLAIR's genome-wide 17.3 / 6.2%; the gorilla S1 margin
is 2.6 points). **Limits, stated plainly:** S1 sensitivity is low (6.3% / 3.9% of the stratum's reference transcripts); in S2 the good
pool doubles the matches but its precision there is below the primary pool's and, in gorilla, below the tools'; the gain over our own
`default_good` is small (S1 +1 / +3, S2 +3 / +12 matches). The multi-copy advantage comes from the good-secondaries pool, not from the
operating point's admission rules. Gorilla was scored by the duplicate session (the same frozen scripts), human by this one.
