# Row 1886 (2026-11-22, user decisions): the short-exon repair is part of `assemble`; `PREFIX.ranked.gtf` is written

Two decisions of the user on 2026-11-22, after rows 1878, 1880 and 1884 (asked in the session, both answered).

## 1. Short-exon repair: default of `assemble`, folded in
User: "Default, folded in", decided against row 1878's registered verdict (TRUTH-ONLY: primates 2 of 4, bar 3 of 4; pays on 9 of
11 datasets, the two gorillas within 4 chains and paying against CAT).

`tools/rustle_pipeline.sh`, stage `assemble`:
- `RUSTLE_SHORTEXONS` unset = 1 when `assemble` is named on its own: `bench/famcmp/short_exons.py` rewrites the reads that skip an
  exon of <= 15 bp into `PREFIX.shortexons.d/reads.bam`, and the stage then runs ONCE on that BAM (the same function, called
  again with the BAM swapped). The copy is removed at the end (`RUSTLE_SHORTEXONS_KEEP=1` keeps it). Side files as before:
  `PREFIX.shortexons.units.tsv`, `.reads.tsv`, `.log`.
- `RUSTLE_SHORTEXONS=0` = the command before. A value other than 0 / 1 stops the run. A failed repair stops the stage and names
  the switch.
- **Inside `all` unset = 0.** `families` and `assign` read the original BAM, whose reads still skip the repaired exons; neither
  stage was measured on repaired models (row 1878 says so). `RUSTLE_SHORTEXONS=1` forces it there. This follows the precedent
  of the seeding pool (row 1812: `assemble` alone and `all` differ). Promoting it inside `all` needs that measurement first.
- The opt-in stage `shortexons` (a second assembly beside `PREFIX.gtf`) is unchanged; its inner run sets `RUSTLE_SHORTEXONS=0`.
- The floor-1 run of stage `oneread` follows the same rule as its caller (named alone: repaired; inside `all`: not).

**Parity** (workspace `r1885_fold_parity`, alias `r1886_shortexons_default`; jobs 4494751 heart, 4494752 WTC11 r2; patched driver
= this commit's):
| dataset | `assemble`, nothing set | `RUSTLE_SHORTEXONS=0` | wall, 40 threads against 24 |
|---|---|---|---|
| mouse heart | `x.gtf` and `x.final.gtf` byte-identical to row 1878's `se/x.shortexons.gtf` / `.final.gtf` | byte-identical to row 1878's `base` | 196 s against 74 s |
| WTC11 r2 | the same | the same | 207 s against 98 s |
The BAM copy was gone after both runs. Tests: `tests/pipeline_seed_pool.rs` (two new cases: one assembly on the rewritten BAM,
copy removed; a failed repair stops the stage, a bad value is refused), 40 of 40; `tests/pipeline_driver_dispatch.rs` 26 of 26.

**New default numbers** (`final.gtf`, chains / correct / precision; row 1878's `se` arm = this default):
fly simulated 6,938 / 6,762 / 0.9746; mouse simulated 25,826 / 25,141 / 0.9735; zebrafish simulated 13,271 / 12,147 / 0.9153;
mouse ES 41,669 / 17,093 / 0.4102; heart 21,183 / 12,950 / 0.6113; adrenal 23,046 / 14,184 / 0.6155; cortex 29,803 / 15,526 /
0.5210; OR6737 57,859 / 24,311 / 0.4202; KB3781 53,175 / 23,406 / 0.4402; chimpanzee testis 46,976 / 20,317 / 0.4325; WTC11 r2
31,345 / 20,054 / 0.6398. Not rerun under the new default: A119b, HG002 and the other seven ape samples.

**Known costs.** Wall time of the stage about 2-2.7 times on these two libraries (the rewrite is a second pass over the BAM);
1,351 s on the fly genome in row 1878 (one temporary BAM per contig, 1,870 contigs); the work directory needs room for a
second BAM. Every figure and base run from before this commit predates the default.

## 2. `PREFIX.ranked.gtf`: one ranked list, two tiers
User: "Write two tiers" (the third tier of row 1884 is not written).

- New `tools/ranked_list.py FINAL.gtf ONEREAD_RANKED.gtf OUT.gtf`: tier 1 = every transcript of `PREFIX.final.gtf`, most `reads`
  first; tier 2 = the one-read chains of `PREFIX.oneread.ranked.gtf` whose intron chain `final.gtf` lacks, best `oneread_score`
  first. Transcript lines get `list_tier`, `list_rank`, and `list_sensitive "1"` for the named cut (tier 1 + the best 10% of
  tier 2). No transcript is changed, added or removed.
- Written by the `oneread` step (so by `all`, where that step is a default): `RUSTLE_RANKED_LIST=0` = off; a failure is a
  warning and leaves the other products alone. It needs `PREFIX.final.gtf` (no file with `RUSTLE_ASSEMBLY_LAYERS=off`).
- Check: on WTC11 r2 the file's `list_sensitive` set is 36,611 chains with 21,068 annotated, the figures of row 1880's table.
- What it is for (row 1880): cut at another assembler's output size the list held at least that assembler's annotated chains
  on 65 of 67 comparisons over twenty datasets, 66 of 67 with the third tier of row 1884. Tier 2 is one read per chain: a recall tier.
- `list_rank` counts single-exon transcripts too; row 1880 counted multi-exon chains only.

## Not done
- `all` with the repair on: families / assign on repaired models.
- The public comparison page and the transfer-test page show the numbers from before this default.
- Speed of the rewrite on genomes with thousands of contigs.

## 3. Third tier of the list, opt-in (row 1891, user 2026-11-22: "Add third tier, opt-in")
- `RUSTLE_RANKED_TAIL=1` (unset = 0): the `oneread` step runs one more assembly with a one-read floor and the dropping rules off
  (isoform fraction 0, sub-chain ratio 0, retained-intron ratio 0, spliced share 0, the gate bypass; row 1884's command) and
  `tools/ranked_list.py --tail` appends its multi-exon chains that neither earlier tier holds, by reads, as `list_tier "3"`.
  Their transcript ids get the suffix `.tail` (a second assembly reuses ids). The relaxed run's files are removed after use.
- New switch `RUSTLE_POLISH_ISOFORM_FRACTION` (unset = 0.02, the value the command always had): exists so the relaxed run can
  turn the rule off; no other value is a judged default (rows 1837, 1876). A malformed value stops the run.
- Checks. Tool: on HG002 with row 1884's files the written list holds 1,365,473 multi-exon chains and, cut at FLAIR's 728,951,
  45,217 annotated in file order and 45,207 as the expected value over the tie group the cut splits: row 1884's figure exactly.
  Driver (`r1891_ranked_tail`, job 4494855, mouse heart): tiers 1 and 2 are line-identical with and without the switch
  (21,335 + 24,400), tier 3 adds 42,847 chains, the `oneread` step takes 489 s against 278 s. Driver tests 66 of 66.
- With the third tier the list is ahead of FLAIR on 15 of 15 datasets (row 1884). It stays off by default: about a million
  chains at under 2% annotated on deep libraries.
