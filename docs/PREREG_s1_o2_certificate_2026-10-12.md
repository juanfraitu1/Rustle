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

## Result (2026-10-12): the secondary-only copies are NOT supported by O2 — their reads belong to the other copy

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.s1o2.{txt,records.tsv,genes.tsv}`. Instrument self-check (gorilla read
m64076_221110_210557/127338300/ccs): every `=` column's reference base equals the read base on both placements (0 / 1,827 and
0 / 1,457), so the columns and bases are right.

| | gorilla | human |
|---|---|---|
| good secondary records tested (S1 genes) | 14,698 | 18,613 |
| AS == molecule best ("tied"): to S1 copy / to primary copy / abstain | 3 / 7,531 / 414 | 2 / 3,883 / 5,947 (5,859 with no distinguishing column) |
| AS < best (strictly worse): to S1 copy / to primary copy / abstain | 16 / 6,259 / 475 | 24 / 7,292 / 1,465 |
| read matches the PRIMARY copy at every distinguishing column | 9,780 (66.5%) | 7,107 (38.2%) |
| read matches the S1 copy at every distinguishing column | 13 | 16 |
| mosaic (>= 2 columns for each copy) | 1,805 (12.3%) | 3,655 (19.6%) |
| S1 genes: supported / contradicted / undecidable | 1 / 400 / 131 | 1 / 352 / 419 |
| ... with an OP_good transcript | 1 / 151 / 49 | 1 / 202 / 134 |
| OP_good transcripts there (reference-exact): supported / contradicted / undecidable | 4 (2) / 193 (38) / 45 (16) | 1 (0) / 267 (44) / 185 (3) |

(The transcript counts here use every S1 gene; register 1268's 171 / 224 used the S1 genes with a multi-exon reference transcript.)
**Reading.** In both species exactly ONE secondary-only gene has >= 2 reads certified to its own copy. Most S1 genes are
contradicted (gorilla 400 / 532) or undecidable because the copies are identical over the reads (human 419 / 772, almost all with no
distinguishing column). The transcripts our good pool builds there — including most of the reference-exact ones — are SHADOWS of the
paralog the reads came from, rebuilt on a copy with the same exon structure. Register 1268's C1 ("wins on secondary-only loci")
stands as scored but does not mean those copies are expressed. Corrections to the wording of the earlier preregistrations: (1) "tied"
(AS == the molecule's best) is not "tied with the primary" — the primary's AS is below the best in 3.8% (gorilla) / 3.9% (human) of
multi-record molecules; true ties (primary == best == second) are 0.5% / 1.6%; (2) S2 ("coin-toss-starved") used the same "tied" and
deserves the same certificate check before it is read as O2's regime. Useful side result: O2's certificate separates the shadows
(contradicted) from the undecidable — the step that should prune good-pool transcripts, as the user's 10-11 rule puts it (a primary or
secondary is not trusted as the origin until O2 says so).
