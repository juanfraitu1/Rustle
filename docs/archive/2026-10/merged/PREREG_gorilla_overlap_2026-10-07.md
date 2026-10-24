# Pre-registration: overlapping genes on real reads, gorilla OR6737 testis: our assembler against StringTie, FLAIR and isoseq collapse (2026-10-07)

**Written before the truth tables, the instrument or any arm of this study exists.** What was already seen: the input-only counts of the gorilla design sweep (2,202 annotated gene pairs share >= 100 exonic bp; 423 with both genes >= 3 exact-chain reads; 35 on the same strand; `docs/GORILLA_FEWCOPY_DESIGN_2026-10-07.md`), the 15-cell genome-wide bakeoffs of the earlier default against the lab's tools (register 1056-1057: aggregates, never per gene), and the ideal-read results of `docs/ENTANGLED_BASELINE_2026-10-06.md`. Nothing of any arm at these genes. User direction (2026-10-06/07): 'lets do 1' (overlapping genes on real gorilla reads).
**Label.** Descriptive real-read check on a partly spent library (OR6737 decided the seeding rule, F1/F1v2 and the strict-junction adoption; the BAM is cross-individual: reads OR6737, assembly and annotation KB3781). Not held-out, no bar except the one replication claim G1.

## 1. Inputs

BAM `winloci_data/GGO_mm.bam` (minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, 10,705,377 records); FASTA `_from_wsl/winloci_scratch/GGO.fasta`; annotation RefSeq `winloci_data/GGO_genomic.gff` (GCF_029281585.2). Per-read genome-wide best AS: `as_table` of the HEAD build on the BAM (name, best_as, ...). Lab arms (fixed products on this BAM, versions as in `docs/archive/2026-09/PREREG_copy_recovery_tools_2026-09-29.md` section 2): **S** `benchmark_collapse/stringtie_GGO/GGO.stringtie.gtf` (StringTie 3.0.1, `-L`, no `-G`), **F** `benchmark_collapse/flair_GGO/GGO.flair.isoforms.gtf` (FLAIR 3.0.1, bam2bed + `collapse --trust_ends`, no correct, no annotation), **I** `isoseq_upload/isoseq_GGO_OR6737/GGO_OR6737.collapsed.gff.gz` (isoseq collapse defaults). Our arms (HEAD release binaries, the driver's `assemble` stage genome-wide, default f1v2 / GOOD secondary seeding / strict junctions): **D** default, **P** `--no-seed-secondaries`, **PC** P plus `RUSTLE_POLISH_SUBCHAIN=drop`.

## 2. Truth and strata (inputs only; computed before any arm is scored)

- A VALID annotated chain: a RefSeq transcript with >= 2 exons whose every intron is >= 50 bp with a canonical motif (GT-AG, GC-AG, AT-AC on the transcript strand, read from the FASTA); chain = (chrom, strand, ordered intron list), distinct per gene.
- Read pools per alignment (E0 rule of `docs/PREREG_ideal_expression_2026-10-06.md`): an alignment carries a chain iff its junctions (N >= 50 bp) lying inside the gene's exon-union span equal the chain. **P1** = primary alignments (flag & 2308 = 0); **P2** = P1 plus secondary alignments (flag 256, not 2048) with AS >= 0.98 x the molecule's genome-wide best AS. Counted as distinct read names.
- A chain is EXPRESSED iff >= 3 reads carry it exactly. **Primary denominator: >= 3 reads of P1** (what StringTie, FLAIR and isoseq collapse see: they use primaries only); **secondary denominator: >= 3 reads of P2** (the pipeline's seeding pool, generous to our arms). Both are reported, never mixed.
- Genes with >= 1 expressed chain are evaluated. Stratum from the exon unions of ALL annotated genes (expressed or not): **E_both** = shares >= 100 exonic bp on the same strand with another gene that also has an expressed chain; **E_one** = same, partner without an expressed chain; **A** = no same-strand overlap but >= 100 bp with an opposite-strand gene; **N** = none. Inside E_both: **E_j** = an annotated valid chain of the gene shares an exact junction with an annotated valid chain of the partner (not separable by any junction rule), **E_x** = no shared junction.

## 3. Metrics (per arm and stratum; `bench/entangled/real_score.py`)

- **M1** expressed annotated chains recovered exactly by >= 1 transcript of the arm (same chromosome and strand); genes COMPLETE.
- **M2 annotation-match share** (a precision PROXY: unannotated isoforms count against every arm alike): among the arm's multi-exon transcripts whose exons overlap, on the same strand, an evaluated gene of the stratum, the share whose chain equals ANY valid annotated chain.
- **M3 resolution**: a gene is RESOLVED iff an arm gene_id carries an exact chain of it and an exact valid annotated chain of no other gene; reported for E_x and E_j. Caveat as in the ideal study: FLAIR names genes by a 1 kb start bin (over-splitting), StringTie fuses.
- Cross-check: `gffcompare` intron-chain sensitivity of D, S and F against the expressed-chain truth GTF of E_both.

## 4. Predictions (fixed now)

**G1** D recovers more expressed chains than S and than F at E_both and at N (primary denominator): the replication claim. **G2** P recovers fewer chains than D at E_both by >= 3 on the secondary denominator (the seeding gain is visible on real reads). **G3** M2 of D is below that of S and of F (the default over-emits); M2 of PC is >= the smaller of S and F. **G4** at E_x, D resolves at least as many genes as S. **G5** I recovers fewer chains than D at E_both.

## 5. Limits declared in advance

One library, cross-individual, RefSeq-only annotation (partly Gnomon, which may cite long reads: if OR6737 is among them the truth is partly circular, unverified); the secondary denominator favours arms that seed secondaries; the lab tools are in their default annotation-free recipes at versions that differ from the local ones; M2 is a proxy; 35 same-strand E_both pairs are few: E_both is reported with its size and counts, not percentages alone.

## Errata (2026-10-07, appended after the independent verification wf_994ece76-8d2; the text above is unchanged)

(1) The user's direction was one message, 'lets do 1 and 4'; 1 is this study. (2) Exposure: the seeding rule, F1 / F1v2 and the DAZ2 gate were decided with OR6737 (dev); strict junctions + retained-intron ratio 10, F1 (minus NC_073244.2) and R_J are held-out and spent (`docs/GORILLA_FEWCOPY_DESIGN_2026-10-07.md`). (3) The realised E_both is 54 same-strand pairs (108 genes), not the 35 of the design sweep (a stricter exact-chain-read rule); 202 gene-chains are 201 distinct chains. (4) Scoring choices left open and fixed by the implementation, effect on E_both none: M2 'exons overlap an evaluated gene' = overlap with the exon union (a span reading gives 48,720 / 432 instead of 48,670 / 431), M3 'an exact chain of it' = any valid chain of the gene. (5) G2 and G3 could not be met by any arm (the secondary-only headroom at E_both is 2 chains; F emits twice D's transcripts): they were mis-specified bars, see the results document. (6) The 'partly Gnomon' wording of section 5 is 99.9% Gnomon, with 'long SRA reads' cited by the models behind 98.9% of the expressed chains.
