# Truth frozen for PREREG_rna_allele_haplotype_count_2026-10-01 (before any RNA call)

Work dir `/mnt/linuxdisk/tmp/rna_allele/`; scripts `bench/rna_allele/` (chrmap.py, align_driver.sh, genes.py, truth_lift.py, sets.py,
truth_classes.py). sha1 of the frozen files:

```
e9182791acc5d6eff9ab96ca4ea7a69b976c3d1e  chrmap.tsv
dfa2ff0a27037aa5b01abd744da1f581dd17f1cc  genes.tsv
8a609fa10a050027645a28ea0060d5cb8a8bfcdd  lift.tsv
03699443a62b1d192eb8c167a1be7812a09ae81a  sets.tsv
2f86faa7a5943f88349ac465be1aec718ae183df  truth.tsv
bd483ffc0be7b6be50209f1073c1c46d5db18e93  truth_fam.tsv
64c76f878c3515fd9240a79076b3d90b4965994f  t1fam.mat.paf
297378dbcda6d794d6bcd1f1ac8ec3b0d44702b6  t1fam.pat.paf
bb66195a4208ecb66c8230438ac11c3b19853f9a  para.paf
```

- `_pri` = 16 paternal + 9 maternal chromosomes, each byte-identical (uppercase md5) to one haplotype's chromosome; 23 autosomes aligned
  to their other-haplotype partner (`minimap2 -x asm5 -c --cs`, query in 10 Mb pieces; chr22 and chr23 in 2 Mb pieces because their
  first piece took 430-530 s).
- Gene sets (paralog test `minimap2 -c -x splice -N 50`, default -p): S_fam 39, S_multi 4,181, S_single 34,091, S_X 1,427, excluded 1,438.
- Classes:

| set | T2d | T2i | T1 | T? |
|---|---|---|---|---|
| S_fam | 28 | 6 | 0 | 5 |
| S_multi | 2,351 | 1,549 | 169 | 112 |
| S_single | 22,901 | 10,976 | 11 | 203 |
| S_X | 0 | 0 | 1,427 | 0 |

- Families: NPIP 25 copies (T2d 14, T2i 6, T? 5) + 1 B-only locus (maternal chr18, CM054600.2:17,395,900-17,418,169) -> T = 41 haplotype
  copies; TBC1D3 14 copies (all T2d) -> T = 28.
- Deviation found and fixed before any RNA call: the haplotype splice indexes name chromosomes `chrN_<hap>_hsa*`; truth_classes.py now
  aliases them to the GenBank accession (checked by length). Before the fix no B-haplotype hit was recognised.
- Note: syntenic exonic differences at NPIP copies are high (e.g. NPIPA7 24 mismatches + 6 indels): in SD regions the chromosome
  alignment may pair non-allelic copies. Any difference makes a copy T2d, so this can only move copies from T2i to T2d.
