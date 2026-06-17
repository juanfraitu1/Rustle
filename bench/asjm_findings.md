# Allele-specific junctions — MULTI-SNP haplotype phasing

Each read is assigned to one of the two diploid haplotypes using ALL het SNPs in the gene (read-backed 2-means), not a single anchor SNP — so reads covering any subset of SNPs get phased (more reach/power), and the phasing quality (how cleanly reads bipartition) is a built-in confound check. With 1 SNP it degenerates to single-anchor, so this is a SUPERSET.

## Result
- junctions tested: **82179**
- **allele-specific junctions (q<0.05, |dPSI|>=0.3): 630** across **263** genes
- high-confidence (>=1 transversion SNP [genetic] AND phase_qual>=0.7): **453** / 184 genes
- from genes with >=2 het SNPs (multi-SNP phasing does real work): **559** / 241 genes

## Gain over single-anchor phasing
- single-anchor ASJ (committed): **475**; multi-SNP ASJ: **630**
- shared (gene+junction): **448**; **multi-SNP-only: 182** (the reach/power gain); single-anchor-only: 27
- net: multi-SNP ADDS +155 ASJ (263 genes)

## n_snp distribution of ASJ
```
  n_snp=1 : 71
  n_snp=2 : 64
  n_snp=3 : 54
  n_snp=4 : 35
  n_snp=5 : 35
  n_snp=6 : 36
  n_snp=7 : 38
  n_snp=8 : 18
  n_snp=9 : 49
  n_snp=10+: 230
```

## Top multi-SNP-only allele-specific junctions (reach/power gain)
| gene | chrom | n_snp | phase_qual | junction | psi0 | psi1 | dPSI | q |
|---|---|---|---|---|---|---|---|---|
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125054945-125055563 | 0.0 | 1.0 | 1.00 | 6.7e-12 |
| LOC109026500 | NC_073227.2 | 2 | 1.0 | 9976297-9976330 | 0.0 | 1.0 | 1.00 | 9.6e-08 |
| LOC129527636 | NC_073242.2 | 6 | 0.77 | 35228529-35228571 | 1.0 | 0.0 | 1.00 | 2.8e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21768385-21768637 | 1.0 | 0.0 | 1.00 | 9.7e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21767764-21767806 | 1.0 | 0.0 | 1.00 | 9.7e-04 |
| LOC115931062 | NC_073242.2 | 40 | 0.81 | 21766564-21766606 | 1.0 | 0.0 | 1.00 | 1.5e-03 |
| LOC109025073 | NC_073244.2 | 9 | 1.0 | 18096968-18106041 | 1.0 | 0.0 | 1.00 | 3.3e-03 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47843875-47845873 | 1.0 | 0.0 | 1.00 | 1.8e-02 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47843875-47844486 | 0.0 | 1.0 | 1.00 | 1.8e-02 |
| LOC101146691 | NC_073236.2 | 2 | 1.0 | 47844579-47845873 | 0.0 | 1.0 | 1.00 | 1.8e-02 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125054945-125065033 | 0.981 | 0.0 | 0.98 | 8.6e-10 |
| LOC109025073 | NC_073244.2 | 9 | 1.0 | 18106165-18107316 | 0.923 | 0.0 | 0.92 | 9.1e-03 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125057247-125057955 | 0.078 | 1.0 | 0.92 | 1.4e-16 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125056559-125057195 | 0.078 | 1.0 | 0.92 | 1.4e-16 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117517117-117517732 | 0.087 | 1.0 | 0.91 | 2.0e-19 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117513957-117514560 | 0.091 | 1.0 | 0.91 | 1.5e-17 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117513186-117513905 | 0.091 | 1.0 | 0.91 | 1.5e-17 |
| LOC115932779 | NC_073242.2 | 20 | 0.44 | 28928691-28933193 | 0.909 | 0.0 | 0.91 | 1.9e-02 |
| LOC134756861 | NC_073224.2 | 40 | 0.86 | 117512378-117513013 | 0.093 | 1.0 | 0.91 | 1.2e-16 |
| LOC129525331 | NC_073224.2 | 40 | 0.79 | 125058128-125058727 | 0.094 | 1.0 | 0.91 | 3.1e-16 |

## Caveats
- phase_qual is the diploid-bipartition cleanliness; low phase_qual = paralog/>2 haplotypes (down-weighted by the phase_qual>=0.7 high-confidence filter).
- transition-only-SNP loci (n_tv=0) could be edit-phased; the genetic core requires >=1 transversion het SNP.
- still coverage-limited (>=12x, >=5 reads per haplotype per junction).

## Reproduce
- `MINIFORGE python bench/allele_specific_junctions_multisnp.py --chrom <C>`
- `python3 bench/asjm_aggregate.py`
