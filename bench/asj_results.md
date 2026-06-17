# Allele-specific junctions (genome-wide)

Splice junctions whose usage depends on the allele a molecule carries. Long HiFi reads link allele->junction PER MOLECULE (the read carries both its het-SNP allele and its junctions), so no statistical phasing is needed. The het loci that confounded copy-detection (task c) are the substrate: heterozygosity = the phasing signal (confound -> feature).

## Result
- genes with a phaseable het anchor + a tested junction: **7898**
- alternatively-spliced junctions tested (non-constitutive): **74674**
- **allele-specific junctions (BH-FDR q<0.05 AND |dPSI|>=0.3): 475** across **235** genes
- of those, **120** have a TRANSVERSION anchor (unambiguously genetic, not RNA-editing) — the high-confidence genetic set

## |dPSI| distribution of the ASJ calls
- median |dPSI|=0.64, max=1.00; >=0.5: 312, >=0.7: 213, ==1.0 (full switch): 56

## Top allele-specific junctions
| gene | chrom | anchor | alleles | junction | PSI_X | PSI_Y | dPSI | q | anchor type |
|---|---|---|---|---|---|---|---|---|---|
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192941955-192944457 | 1.0 | 0.0 | 1.00 | 1.9e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192944612-192944695 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192949190-192950494 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192947387-192949145 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192944823-192947259 | 1.0 | 0.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| LOC101141065 | NC_073224.2 | 192934097 | G/C | 192930566-192934097 | 0.0 | 1.0 | 1.00 | 6.6e-26 | transversion(genetic) |
| GLB1L3 | NC_073233.2 | 142241047 | C/A | 142240873-142240938 | 0.0 | 1.0 | 1.00 | 3.2e-14 | transversion(genetic) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16419347-16423955 | 0.0 | 1.0 | 1.00 | 4.7e-12 | transition(poss.edit) |
| TSPAN2 | NC_073224.2 | 119456078 | C/G | 119454970-119455031 | 0.0 | 1.0 | 1.00 | 7.0e-12 | transversion(genetic) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16426010-16429060 | 0.0 | 1.0 | 1.00 | 5.5e-11 | transition(poss.edit) |
| CNTNAP3C | NC_073236.2 | 33312391 | G/A | 33312509-33313164 | 0.0 | 1.0 | 1.00 | 9.2e-11 | transition(poss.edit) |
| CASP10 | NC_073235.2 | 103835969 | A/G | 103835036-103835169 | 1.0 | 0.0 | 1.00 | 8.1e-10 | transition(poss.edit) |
| LOC115933039 | NC_073244.2 | 20936571 | A/G | 20959729-20961004 | 1.0 | 0.0 | 1.00 | 1.0e-09 | transition(poss.edit) |
| LOC115934626 | NC_073228.2 | 45957287 | T/A | 45956648-45956726 | 1.0 | 0.0 | 1.00 | 3.3e-09 | transversion(genetic) |
| HIGD1A | NC_086017.1 | 52690273 | C/G | 52689198-52689242 | 1.0 | 0.0 | 1.00 | 2.5e-08 | transversion(genetic) |
| KCNAB2 | NC_073224.2 | 234526415 | C/G | 234525571-234526031 | 1.0 | 0.0 | 1.00 | 1.9e-07 | transversion(genetic) |
| PSMD2 | NC_086017.1 | 195406803 | T/G | 195406804-195406955 | 0.0 | 1.0 | 1.00 | 5.8e-07 | transversion(genetic) |
| CSNK1D | NC_073228.2 | 9740588 | G/C | 9740897-9740970 | 0.0 | 1.0 | 1.00 | 8.8e-07 | transversion(genetic) |
| LOC134758754 | NC_073224.2 | 19346753 | T/C | 19339269-19339542 | 1.0 | 0.0 | 1.00 | 1.3e-06 | transition(poss.edit) |
| LOC129527496 | NC_073242.2 | 4910233 | C/T | 4908324-4908418 | 0.0 | 1.0 | 1.00 | 3.1e-06 | transition(poss.edit) |
| C13H9orf43 | NC_073237.2 | 109779230 | T/C | 109779303-109779954 | 0.0 | 1.0 | 1.00 | 6.1e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16480542-16484928 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16471500-16475882 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| LOC115932992 | NC_073242.2 | 16399354 | C/T | 16476020-16480404 | 1.0 | 0.0 | 1.00 | 7.4e-06 | transition(poss.edit) |
| MYO18B | NC_086018.1 | 23923474 | G/T | 24012728-24023003 | 1.0 | 0.0 | 1.00 | 1.5e-05 | transversion(genetic) |

## Honest caveats
- **Editing confound:** a TRANSITION anchor (A/G or C/T) could be an RNA-edit site (phasing by edit-status, not allele) — that is edit-coupled splicing, still real but not *genetic* allele-specific. TRANSVERSION anchors are unambiguously genetic; treat the transversion subset as the high-confidence genetic ASJ. (Strand-aware editing detection would reclassify the transitions.)
- **Single-anchor phasing:** reads are phased by ONE balanced het SNP; a junction is tested only among reads spanning BOTH the anchor and the junction (long reads help). Multi-SNP haplotype phasing would add power/reach.
- **Collapsed-paralog masquerade:** at a collapsed locus the two 'alleles' could be paralog copies (=copy-specific splicing, also interesting, but not within-gene allele-specific). The het substrate from task c was uniquely-mapped (frac_mq0=0); genome-wide, low-MAPQ loci should be down-weighted.
- **Coverage-limited:** needs >=12x at the anchor and >=5 spanning reads per allele; deeper data finds more.

## Reproduce
- `MINIFORGE python bench/allele_specific_junctions.py --chrom <C>` (per-chrom) ; `--region` for one locus
- `python3 bench/asj_aggregate.py` (this: BH-FDR + calls)
