# SQANTI3 on the polished `--assemble-only` output (§6q5, 2026-09-19)

Question: do the transcripts the polished assembler emits pass SQANTI3 as real transcripts, and how many
are full-splice matches?

SQANTI3 (conda env `sqanti3`, checkout `/mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3/`):
`sqanti3_qc.py --isoforms <gtf> --refGTF chrN_ref.gtf --refFasta chrN.fa --report skip -t 4`, then
`sqanti3_filter.py rules --sqanti_class <classification> --filter_gtf <corrected.gtf> --skip_report`
(the default rules filter, nothing customised). Same six chromosomes and the same StringTie 3.0.1 `-L -p 4`
runs as `bench/ASSEMBLY_POLISH.md`.

Our arm is the shipped setting:
`--assemble-only --assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7`

## What the polish does to the SQANTI3 profile (chr20, before vs after)

| | raw (`--assembly-polish none`) | **polished** | StringTie |
|---|---|---|---|
| transcripts | 976 | **658** | 699 |
| **full-splice_match** | 352 (36.1%) | **337 (51.2%)** | 332 (47.5%) |
| incomplete-splice_match | 217 (22.2%) | **78 (11.9%)** | 79 (11.3%) |
| novel_in_catalog | 77 (7.9%) | 64 (9.7%) | 84 (12.0%) |
| novel_not_in_catalog | 166 (17.0%) | 154 (23.4%) | 177 (25.3%) |
| antisense | 80 (8.2%) | **9 (1.4%)** | 12 (1.7%) |
| intergenic | 31 (3.2%) | **10 (1.5%)** | 9 (1.3%) |
| genic | 31 (3.2%) | **2 (0.3%)** | 1 (0.1%) |
| genic_intron | 16 (1.6%) | **1 (0.2%)** | 1 (0.1%) |
| fusion | 6 (0.6%) | 3 (0.5%) | 4 (0.6%) |
| **SQANTI3 rules filter: PASS** | 741 / 976 = **75.9%** | 624 / 658 = **94.8%** | 667 / 699 = 95.4% |

⭐The polish raises the FSM share **36.1% → 51.2%** and the rules-filter pass rate **75.9% → 94.8%**, and
it collapses exactly the categories it was designed to remove: **antisense 80 → 9, intergenic 31 → 10,
genic 31 → 2, genic_intron 16 → 1** — the shadow rule (§6q0) doing what it was built for, since those are
precisely the single-exon predictions in a spliced gene's shadow. FSM count falls only 352 → 337 while
318 transcripts are dropped.

## All six chromosomes

| chrom | arm | n | FSM | ISM | artifact cats¹ | rules filter PASS | FSM passing |
|---|---|---|---|---|---|---|---|
| chr20 | **ours** | 658 | **337 (51.2%)** | 78 (11.9%) | **25 (3.8%)** | 624 (94.8%) | **336** |
| | StringTie | 699 | 332 (47.5%) | 79 (11.3%) | 27 (3.9%) | 667 (**95.4%**) | 332 |
| chr11 | **ours** | 1,373 | **687** (50.0%) | 163 (11.9%) | **40 (2.9%)** | 1,309 (**95.3%**) | **683** |
| | StringTie | 1,307 | 657 (50.3%) | 157 (12.0%) | 52 (4.0%) | 1,235 (94.5%) | 653 |
| chr7 | ours | 1,172 | 518 (**44.2%**) | 155 (13.2%) | **112 (9.6%)** | 1,101 (93.9%) | 513 |
| | StringTie | 1,199 | **519** (43.3%) | 158 (13.2%) | 120 (10.0%) | 1,141 (**95.2%**) | **517** |
| chr14 | ours | 904 | 392 (**43.4%**) | **105 (11.6%)** | **50 (5.5%)** | 846 (**93.6%**) | 388 |
| | StringTie | 935 | **393** (42.0%) | 111 (11.9%) | 55 (5.9%) | 858 (91.8%) | **391** |
| chr5 | ours | 1,008 | 476 (**47.2%**) | **135 (13.4%)** | **68 (6.7%)** | 939 (**93.2%**) | 473 |
| | StringTie | 1,060 | **483** (45.6%) | 150 (14.2%) | 74 (7.0%) | 977 (92.2%) | **480** |
| chr9 | ours | 934 | 409 (**43.8%**) | **99 (10.6%)** | 59 (6.3%) | 886 (**94.9%**) | 408 |
| | StringTie | 953 | **410** (43.0%) | 115 (12.1%) | **56 (5.9%)** | 881 (92.4%) | **410** |

¹ genic + antisense + intergenic + genic_intron + fusion.

### Pooled

| | transcripts | FSM | ISM | artifact cats | rules filter PASS | FSM passing |
|---|---|---|---|---|---|---|
| **ours** | 6,049 | **2,819 (46.6%)** | **735 (12.2%)** | **354 (5.9%)** | **5,705 (94.3%)** | **2,801** |
| StringTie | 6,153 | 2,794 (45.4%) | 770 (12.5%) | 384 (6.2%) | 5,759 (93.6%) | 2,783 |

⭐**Pooled, we beat StringTie on every SQANTI3 quantity**: more full-splice matches (2,819 vs 2,794) at a
higher FSM share, fewer ISM, fewer artifact-category transcripts, a higher rules-filter pass rate, and
more FSM transcripts surviving the filter — from fewer emitted transcripts.

Per chromosome the pass rate is ours on chr11/chr14/chr5/chr9 and StringTie's on chr20/chr7, and the FSM
*share* is ours on five of six (chr11 is 50.0 vs 50.3).

## What is filtered out, and why it is not fixable into FSM

chr20, of our 34 rejects: novel_not_in_catalog 14, ISM 7, novel_in_catalog 4, antisense 4, fusion 2,
genic 2 — the same profile as StringTie's 32 rejects (NNC 15, ISM 7, NIC 5, antisense 2, fusion 2,
genic 1). **NIC/NNC transcripts cannot be "made FSM"**: they are novel junction combinations, i.e. either
real unannotated isoforms or mis-assemblies, and calling them FSM would require changing the reference,
not the assembler. The polish's contribution is to remove the categories that genuinely are artifacts
(ISM fragments and the single-exon shadow classes) — which is why the FSM *share* rises even though the
FSM count falls slightly.
