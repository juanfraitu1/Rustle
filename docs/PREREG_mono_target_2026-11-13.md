# PREREG row 1668: can a targeted rule pick the good single-exon models out of the dropped ones? (2026-11-13)

Assembly-layer step 0, offline, no code change. User request after row 1667 (all three single-exon rules off: 103 of
291 no-model members gain a model, at +236,544 single-exon models).

## Data
The single-exon models of row 1667's all-off run (`r1667_mono_attr/g_psq/a1.families.gtf`) that are absent from the
baseline (`r1666_mono_root/g_off/a1.families.gtf`). Per model, features the assembler already has:
`reads` (GTF attribute), length, `pas` (AATAAA or ATTAAA in the 40 bp before the 3' end), `primed` (20 bp past the 3'
end >= 60% A or A6), shadow class vs the baseline multi-exon transcripts (`exon` = overlaps a multi-exon exon on
either strand; `span` = inside a multi-exon span on its strand; `free`).
Labels (famtruth.map_copy, own strand): `member` = one of the 291; `annotated` = another CAT gene, model covers
>= half its exon bases; `none`.

## Rule search (frozen before running; `bench/stgaps/mono_features.py`)
Train = odd-numbered autosomes + chrX; test = even-numbered autosomes + chrY. Grid: reads >= {2, 3, 5, 10, 20, 50};
shadow in {any, not exon, free}; pas in {any, yes}; primed in {any, no}; length <= {any, 10 kb, 5 kb}. On train, take
the rule that recovers the most members with <= 50 added models per member (ties: fewer added models).
Go: that one rule, on test, recovers >= 15 members at <= 50 added models per member. Then it is registered as a
switch-guarded Rust rule with HG002 as the unseen sample (A119b is spent, gorilla has too few single-exon drops).
Otherwise the single-exon line closes: the members cannot be told from the noise with these features.
The full train/test table is reported either way.

## Result: STOP (job 4468883; the first run, 4468879, read `reads` from exon lines and saw 0 everywhere: invalid, rerun)
Train 122,326 models, 63 members; test 114,218 models, 40 members. No rule of the grid reaches 50 added models per
member on train. Best ratios on train: reads >= 10, free, PAS, not primed, <= 5 kb: 3 members in 163 models (54 per
member); reads >= 5, free, PAS, not primed: 9 members in 1,441 (160 per member).
The members do not differ from the unannotated models: reads quartiles 3 / 7 / 11 (members) vs 3 / 6 / 14 (none);
PAS 22% vs 13%; primed 53% vs 65%; free of any multi-exon model 66% vs 53%. Annotated genes that already have a model
are the well-supported ones (reads 8 / 29 / 119, 73% on a multi-exon exon). The single-exon line closes: with reads,
length, PAS, priming and shadow class, the missing members cannot be told from the noise.
