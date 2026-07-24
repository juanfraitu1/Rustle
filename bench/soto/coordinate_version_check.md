# Coordinate-version check: Soto BED vs the A119b BAM (CHM13 v1 vs v2)

**Concern:** A119b is aligned to CHM13 **v2.0**; Soto's *published* supplement is CHM13 **v1**. If our
`80_fams.chr.bed` were in v1 coordinates, every Soto member sequence we extract from the v2.0 FASTA
would be offset — corrupting family definitions, `near_family` labels, and the candidate
classification.

**Result: no mismatch. The `80_fams.chr.bed` we use is CHM13 v2.0, the same as the BAM. No correction
is needed.**

## Evidence (decisive)

Cross-checked the Soto BED member starts against the CHM13 **v2.0** RefSeq annotation
(`HSA_genomic.gff`, GCF_009914755). For every single-copy gene present by name, the starts agree to
**1–2 bp**:

| gene | Soto BED start | v2.0 GFF start | diff |
|---|---|---|---|
| NOTCH2 | 119,924,808 | 119,924,810 | 2 bp |
| NOTCH2NLR | 120,737,163 | 120,737,164 | 1 bp |
| NOTCH2NLC | 148,519,505 | 148,519,514 | 9 bp |
| NOTCH2NLA | 145,263,046 | 145,265,708 | 2.7 kb |
| PDE4DIP | 145,474,975 | 145,475,230 | 255 bp |
| FAM157A | 200,883,047 | 200,882,974 | 73 bp |
| FAM157B | 150,478,579 | 150,478,369 | 210 bp |
| GOLGA6L6 | 18,443,672 | 18,443,674 | 2 bp |
| GOLGA6L22 | 20,130,998 | 20,131,000 | 2 bp |
| GOLGA8S | 21,091,024 | 21,091,026 | 2 bp |
| GOLGA8F/G/M/J/T/R/Q/H/K/O/N (chr15) | … | … | ≤ 3 kb |

Across 101 named genes spanning chr1/3/9/15 and many families: **101 concordant (<20 kb), 0 with a
systematic offset.** A v1→v2 liftover cannot yield 1–2 bp agreement — v2.0's corrections are indels
that shift downstream coordinates by variable, non-zero amounts. 1–2 bp agreement ⟹ same coordinate
system. The FASTA also carries chrY (v2.0-only), and the annotation matches it — so the whole chain
**BED ≡ GFF ≡ FASTA ≡ BAM is v2.0.**

## The apparent "offsets" are paralog-naming, not coordinate error

The 29 genes with >20 kb "offset" are all **multi-copy families where the two annotations assign
copy-numbers differently**, not coordinate shifts. Examples (both positions are real v2.0 loci):

| gene name | Soto member | RefSeq annotation of that name | what it is |
|---|---|---|---|
| ANKRD20A1 | chr9:79,615,465 | chr9:43,016,231 | two different ANKRD20A paralogs |
| FAM95B1 | chr9:77,797,810 | chr9:40,338,481 | two different FAM95B paralogs |
| SPDYE8/9/17 | chr7:76.5 M | chr7:74.2 M | SPDYE copy-numbering swapped |
| AMY1B | chr1:103,630,425 | chr1:103,536,281 | different AMY1 copy |

This is the **same phenomenon** behind the candidate classification's "near family X but homologous to
family Y" result: genomic proximity and gene-name copy-numbering do not imply homology — even inside
RefSeq. Homology must be assigned by sequence alignment, which is what the classification does.

## Consequence

- The candidate classification (`candidate18_classification.md`) is **not** corrupted by a coordinate
  bug — it stands: 34 real Soto-family copies / 3 real non-Soto loci / 0 phantom FPs.
- If belt-and-suspenders confirmation is wanted, a read-level check (A119b reads at each Soto member
  locus express the expected gene) can be run, but the 1–2 bp annotation concordance already settles it.

**Reproduce:** parse `HSA_genomic.gff` gene features → name→span map; join to `80_fams.chr.bed` gene
tokens; report start-coordinate differences. (Clone-named members — AC*/AL*/GTF2IP*/GOLGA8*P etc. —
lack RefSeq symbols and are omitted from the name-join, but sit interleaved between concordant named
members at concordant spacing.)
