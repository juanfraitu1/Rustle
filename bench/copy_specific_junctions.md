# Copy-specific junctions — do paralog COPIES splice differently?

The unification: (1) discovered cross-chromosome copies → (2) reads assigned per copy (by locus here; PSV-based copy_split applies in the collapsed regime, data-limited in GGO) → (3) per-copy junction usage compared between copies (the ASJ machinery, allele→copy).

- cross-chrom recent-dup pairs (core_ident≥0.9, both ≥3 exons, BOTH copies expressed ≥15 reads): **248**
- homologous junctions compared between copies: **696**
- copy-specific junctions raw (q<0.05, |dPSI|≥0.3): **191**, of which:
  - **DIFFERENTIAL splicing (both copies splice): 146** across **66** copy pairs
  - one copy UNSPLICED/retrocopy (trivial dPSI): 45
- copy-private exon junctions (one copy includes an exon the other lacks): **1408**

## Top DIFFERENTIAL copy-specific junctions (both copies splice, differ at this junction)
| copy A | jA | PSI_A | copy B | jB | PSI_B | dPSI | q |
|---|---|---|---|---|---|---|---|
| LOC101147259 | 93534519-93538824 | 0.0 | LOC115933254 | 32884967-32889244 | 1.0 | 1.00 | 1.1e-19 |
| HERC2 | 29471916-29472174 | 1.0 | LOC101149126 | 68795812-68796072 | 0.0 | 1.00 | 1.8e-03 |
| HERC2 | 29473320-29476000 | 1.0 | LOC101149126 | 68791948-68794664 | 0.0 | 1.00 | 4.3e-03 |
| LOC115932779 | 28925295-28926477 | 0.0 | SORL1 | 129494627-129495811 | 0.995 | 0.99 | 1.1e-29 |
| LOC115933039 | 21007327-21008509 | 0.0 | SORL1 | 129494627-129495811 | 0.995 | 0.99 | 6.7e-28 |
| LOC101124243 | 6830621-6832115 | 0.0 | NOM1 | 165843104-165844598 | 0.974 | 0.97 | 2.0e-25 |
| LOC101147259 | 93526275-93527224 | 0.0 | LOC115933254 | 32876845-32877795 | 0.964 | 0.96 | 2.0e-19 |
| LOC101125415 | 17453004-17453086 | 0.961 | LOC101139163 | 141285207-141285366 | 0.0 | 0.96 | 2.3e-17 |
| LOC101125415 | 17453004-17453086 | 0.961 | LOC101152938 | 75725353-75725435 | 0.0 | 0.96 | 4.8e-12 |
| FRG1 | 210904802-210905583 | 0.955 | LOC115933219 | 42792484-42793264 | 0.0 | 0.95 | 2.1e-32 |
| LOC129529226 | 15687227-15687359 | 0.951 | PLA2G4B | 42658550-42658680 | 0.0 | 0.95 | 1.1e-17 |
| LOC129531205 | 4579716-4582007 | 0.0 | SAFB2 | 11378953-11381252 | 0.949 | 0.95 | 2.9e-17 |
| LOC101137900 | 32597303-32597305 | 0.0 | LOC129531862 | 31668721-31671325 | 0.941 | 0.94 | 3.5e-07 |
| LOC101142904 | 2232115-2232899 | 0.935 | LOC115933219 | 42792484-42793264 | 0.0 | 0.94 | 3.4e-21 |
| LOC129531982 | 29073750-29073882 | 0.931 | PLA2G4B | 42658550-42658680 | 0.0 | 0.93 | 6.3e-15 |
| LOC101125415 | 17453352-17454181 | 0.962 | LOC101139163 | 141284114-141284941 | 0.038 | 0.92 | 6.4e-16 |
| DPY19L2 | 37744041-37746734 | 0.922 | LOC115935345 | 37139557-37143806 | 0.0 | 0.92 | 5.3e-09 |
| DAZ1 | 42826894-42827637 | 0.074 | DAZL | 26430807-26431538 | 0.992 | 0.92 | 2.6e-151 |
| LOC101142904 | 2223093-2224079 | 1.0 | LOC129526395 | 24713570-24714553 | 0.087 | 0.91 | 2.3e-17 |
| LOC129531205 | 4578243-4579576 | 0.071 | SAFB2 | 11381383-11382713 | 0.974 | 0.90 | 8.7e-17 |

## Honest notes
- BOTH copies expressed is rare (~3% of recent-dup pairs) — most paralogs have one dominant copy (consistent with the headroom findings). So this is a small, real set, not a large catalog.
- read→copy assignment here is by genomic locus (cross-chrom copies map to distinct coords). The PSV-based copy_split assignment (interest #2) is exercised only in the COLLAPSED regime, which is essentially empty in GGO (task c) — it needs a deep co-located dataset (testis HiFi for DAZ/RBMY).
- homologous-junction mapping is via exon-length alignment; copy-private exons flag structure divergence (a copy gained/lost an exon).
- coverage-limited; deeper data finds more both-expressed pairs.

## Reproduce
- `MINIFORGE python bench/copy_specific_junctions.py`
