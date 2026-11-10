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

## Outcome
(empty until the run)
