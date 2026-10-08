# PREREG — the RNA-only reference-absent-copy chain on the Y ampliconic genes (written before any run on these families)

Context: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendments 7-10 established, on gorilla autosomal families, a truth-free
chain — family read net -> IsoCon -> flag outputs not in the reference -> LINK outputs within delta of a reference locus (allele) -> MERGE
the rest into candidate copies at delta -> FLAG = candidate with >= 2 transcripts — with delta = the 99th percentile of allelic (haplotype)
divergence at single-copy genes. The Y ampliconic genes (YAGs) are where IsoCon (Sahlin et al. 2018) was built and validated, and they are
haploid: there are no alleles, so delta cannot be the allelic divergence. This prereg fixes how the chain is carried to the Y and what it
must predict there.

## Substrates

- **DEV: human A119b Iso-Seq on CHM13v2.0** (`winloci_data/A119b.t2t.bam`; chrY = HG002's Y; the donor is another male). It expresses the
  germ-cell YAGs (per copy: DAZ 201-354 primaries, RBMY 123-259, HSFY ~1,100, CDY 23-108, PRY 20-31, BPY2 11-30, VCY1B 26, TSPY <= 30).
  Annotation = CAT/Liftoff v2.0 (`gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz`).
- **HELD-OUT: gorilla OR6737 testis Iso-Seq on mGorGor1 `_pri`** (`winloci_data/GGO_mm.bam`; chrY NC_073248.2; another male than the
  reference animal), annotation `winloci_data/GGO_genomic.gff` (LOC genes by description). Run only after the dev run's rules are frozen;
  nothing is tuned on it.
- Not used: the ERR13885926 testis library expresses only TSPY and VCY (DAZ/CDY/HSFY/PRY 0-5 reads per copy).

## delta on a haploid chromosome: delta_Y

On autosomes delta is the individual-vs-reference divergence at single-copy genes, which there IS allelic divergence. On the Y the same
quantity is the donor-vs-reference divergence of a haploid chromosome plus the consensus error of the transcripts, and it is measured the
same way the chain measures everything else: **delta_Y = the 99th percentile of d over IsoCon outputs of the X-degenerate single-copy Y
genes** (SRY, RPS4Y1, ZFY, AMELY, TBL1Y, PRKY, USP9Y, DDX3Y, UTY, TMSB4Y, NLGN4Y, TXLNGY, KDM5D, EIF1AY, RPS4Y2; those with >= 20
primary reads; PAR and X-transposed genes excluded), each gene's net (reads with a record on its locus, <= 1,000) run through IsoCon, d
= 1 - matches / length of the output's best splice hit in the unmasked genome, outputs whose best hit is the gene itself. delta_Y is
computed and recorded BEFORE the YAG panel is run; it replaces 0.00958 in the link and merge steps. The autosomal delta is also run,
reported beside (expected: more copies linked away as "alleles", fewer flags).

## Deletion panel (the truth; no matched Y assembly exists for any donor)

- Families: CAT chrY gene records whose `gene_name` starts with TSPY, RBMY, DAZ, CDY, BPY2, HSFY, PRY or VCY (XKRY has no
  protein-coding record). Overlapping records merge into one copy (RBMY1A1/RBMY1B share a locus); a copy is protein-coding if any member
  is. Copies for the read net = every record of the family (pseudogenes included).
- **Masking rule (blind):** per family, among protein-coding copies with >= 20 primary reads in the dev BAM, mask the one with the MOST
  reads; if the family has >= 4 such copies also mask the one with the FEWEST. Masked copies = N over the merged gene interval of the
  unmasked genome. Families with < 2 protein-coding copies are dropped. All masked copies sit in ONE masked genome (as Amendment 7).
- **d_min of a masked copy** = d of its CAT transcript (the longest protein-coding transcript, spliced from the unmasked genome) against
  the MASKED genome (best splice hit; the link rule's own d). This is the covariate the prediction turns on.
- Reads scored: <= 500 baseline primaries per copy (masked = D reads, surviving = S reads; seed 1). Chain as Amendments 7-9 with delta_Y:
  arm R (masked genome), arm M (masked + the new-copy contigs, components as loci; flag = >= 2 transcripts). Labels as Amendment 7
  (a contig is D-derived when its best unmasked hit lies on the masked interval).

## Rules (fixed now)

- **Y1 (the chain's prediction on the Y):** a masked copy is FLAGGED (a >= 2-transcript candidate holding a D-derived contig) **iff d_min >
  delta_Y and it has >= 6 reads in IsoCon's input**. The prediction must hold for >= 80% of the masked copies (both directions count: a
  flag at d_min <= delta_Y is a miss, as is no flag at d_min > delta_Y with >= 6 input reads). Exact counts reported.
- **Y2 (cost, Amendment 5's rule):** S reads placed on a contig / candidate not derived from their own copy <= 5% of S reads.
- **Reported beside (not the verdict):** per masked copy — d_min, input reads, flagged, D right / wrong / unplaced in R and M; the same at
  the autosomal delta; the within-delta_Y masked copies' reads (they should sit on the nearest surviving copy, linked); IsoCon outputs per
  family and how many exceed delta_Y without any deletion (the family's other members' variants = the false-flag side).
- **Held-out (gorilla):** delta_Y,ggo from the gorilla X-degenerate single-copy genes (same rule); the same masking rule and Y1 / Y2 on
  OR6737 testis. Y1 must hold there as well for the claim "the chain transfers to the Y".

## What this does and does not claim

- It claims the chain detects a missing Y copy exactly when the copy's divergence from its nearest remaining copy exceeds the
  individual-vs-reference divergence of the chromosome, and not otherwise — the identifiability edge, stated as a prediction and tested.
  Byte-identical TSPY copies (register 572) and DAZ copies at 99.97% (647) are therefore PREDICTED undetectable; whether the chain agrees
  is the test.
- It does not estimate the donor's true copy numbers (no truth), and it does not re-run IsoCon's own benchmark.
