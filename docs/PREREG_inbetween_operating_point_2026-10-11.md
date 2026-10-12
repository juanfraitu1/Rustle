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
