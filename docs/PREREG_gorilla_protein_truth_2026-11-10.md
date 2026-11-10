# PREREG r1611: a gorilla truth that names the LOC copies, then r1610 unchanged (2026-11-10)

Follows r1610 (`docs/PREREG_tie_groups_gorilla_2026-11-10.md`, STOP: M = 0). A gorilla gene takes its family from the human
gene of the same symbol; 5,304 of 22,650 protein-coding RefSeq gorilla genes are `LOC` genes with none. The existing
`~/runs/miniprot_gorilla` runs align gorilla proteins to the gorilla genome and do not give a human identity. Written and
committed before anything is run.

## Rule (`bench/famcmp/protein_truth.py`, `DENOVO_GGO_ENS_MAP`)
- Gorilla proteins: `gffread -y` on `DENOVO_GGO_GFF` + `GGO.fa`, the longest protein per gene.
- Human proteins: GENCODE v47 `pc_translations` (downloaded 2026-11-10, `annot/dl_gencode/`, 112,218 sequences).
- `mmseqs easy-search` gorilla -> human; per gorilla gene the hit with the highest bit score. It is accepted when identity
  >= 0.80 and the alignment covers >= 0.50 of the gorilla protein.
- A gene with **no symbol-inherited Ensembl id** takes the accepted hit's Ensembl gene. Genes that already have an id by
  symbol are not changed. The Compara truth then applies as before.

## Check before use (bar)
Among named genes that have a symbol-inherited id and an accepted hit, the hit's Ensembl gene is the symbol's gene in
>= 0.90 of cases. Below that the map is not used and the row STOPS.

## Then
r1610 is rerun unchanged with `DENOVO_GGO_ENS_MAP` set: same frozen recipe, same definition of M, same bars (B1g >= 0.80 of
M, B2s >= 0.90, B5r, B4 <= 1 h), same STOP rule at M < 20. Reported beside it: panel size before and after (176 expressed
members, 75 families).

## Outcome (2026-11-10, job 4456439, `/anvme/workspace/v138da13-rustle/r1611/`): PASS as registered
**Truth check: pass.** 20,741 of 22,650 gorilla proteins have an accepted human hit. Among 16,029 named genes with a hit, 15,918
(0.993, bar 0.90) hit the gene of their own symbol. 4,712 protein-coding genes without a symbol-inherited id get one.

**Panel with the map** (symbol-only in brackets): 171 judged families (75), 661 expressed members (176), 138 of them
reachable only by a tied secondary (1), 34 with MAPQ0 share >= 0.5 (0). Member hit rate under T: ours 0.716, StringTie 0.651,
FLAIR 0.579; on the 138 T-only members ours 0.601, StringTie 0.000, FLAIR 0.000.

**M = 128**: 214 members have no StringTie and no FLAIR model; 128 of them have no decisive molecule (73 have >= 2, 13 one).

| bar | gorilla OR6737 | human A119b (r1607) | verdict |
|---|---|---|---|
| B1g grouped reach >= 0.80 of M | 124 of 128 = 0.969 | 70 of 72 = 0.972 | pass |
| B2s purity, co-duplicated counted, >= 0.90 | 129 of 129 = 1.000 | 0.982 | pass |
| B5r largest group: one chromosome or >= 0.90 one family | 33 loci, 1 chromosome, 23 of 23 annotated loci in one family | 86 loci, chr5 (TAF11L) | pass |
| B4 wall <= 1 h | 528 s (truth 76 s, panel 141 s, records 97 s) | 1,346 s | pass |

- Strict purity (B2f) 117 of 129 = 0.907; 12 groups are impure only through co-duplicated loci.
- 123,558 molecules, 15,780 tied records, 843 loci in 471 groups; 366 loci set aside as interspersed repeat.
- Unreached: LOC129529770, LOC129530492.
- The recipe, frozen on human A119b, transfers to gorilla OR6737 with nothing retuned. Limits: the truth for LOC copies is a
  best-protein-hit family label; `or.gtf` is the base run's assembly; OR6737 has been used in earlier rows.
