# Are the secondary-only copies expressed? O2's certificate on their reads, 2026-10-12

**Pre-registered before the script exists and before any read is tested.** User (2026-10-12), after register 1268 (our good pool
assembles transcripts on secondary-only genes S1, where StringTie and FLAIR have none): "for the loci we get from secondary-only
alignments, does SQANTI prove that they are real?" Answer given: no — SQANTI checks structure, not which copy the reads came from;
a paralog's reads can rebuild an identically structured annotated transcript on a silent copy. The test is O2's certificate.

## Population

The S1 genes of `docs/PREREG_inbetween_operating_point_2026-10-11.md` (strata `/mnt/linuxdisk/tmp/asm_gap/op/<species>.strata`,
benchmark samples gorilla OR6737 / human A119b, development chromosome excluded). For each S1 gene: every GOOD secondary record
(AS >= 0.98 x the molecule's genome-wide best AS) whose aligned blocks overlap one of the gene's exons by >= 1 bp — the records the
good pool used there. Each is paired with the same molecule's PRIMARY record (its placement elsewhere, copy A; the S1 locus is copy B).

## The certificate (O2's, unchanged: `copy_pair_significance` as ported in `bench/asj_o2/asj_o2.py`)

Columns = read positions aligned (`=` or `X`) in BOTH records whose reference bases differ between the two placements (read
orientation; bases from the FASTA): the PSVs this read spans. `k_A` / `k_B` = columns where the read matches A's / B's base (`=` in
that record). Strict maximum (`k_A != k_B`), p = P(X >= k_best), X ~ Binomial(columns, 0.003 / 3); assigned when p < 1e-3 (two
copies), else abstain. Indel columns are not used. Every (secondary record, primary record) pair is tested once.

## Per-gene class (fixed now)

- **supported**: >= 2 records certified to the S1 copy (B) — the read floor of the assembler;
- **contradicted**: 0 records certified to B and >= 2 certified to A;
- **undecidable**: otherwise (mostly: every record abstains — identical copies or too few PSVs).

## What is reported (both species; nothing is adopted or tuned)

1. Records tested; certified to B / to A / abstained; split by TIED (AS == the molecule's best) vs STRICTLY WORSE secondaries.
2. S1 genes by class, all S1 genes and those where OP_good emits a transcript.
3. OP_good's transcripts touching S1 genes, by their gene's class, and how many match a reference transcript exactly (tmap `=`).

**Expected, written now.** Strictly worse secondaries: certified to A far more often than to B (the aligner already prefers A by
>= 1 AS point). Tied secondaries: mostly abstain (identical copies). So most S1 genes will be undecidable or contradicted, and only
"supported" S1 genes may be called expressed copies; that, not register 1268's 34 / 36 matched transcripts, is the claim we would
make for the advisor.

**Instrument** (committed before running): `bench/tool_difference/s1_o2.py` sha1 6b6831e392a2349c7fb5f51b6eb9113b08484117, tests `test_s1_o2.py` (7 pass).
