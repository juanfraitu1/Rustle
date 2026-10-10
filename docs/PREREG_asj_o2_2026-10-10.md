# Allele-specific junctions as an O2 extension: haplotype assignment by the O2 certificate, allele vs copy by the haplotype assemblies, 2026-10-10

**Pre-registered before any read is phased or any junction tested.** User (2026-10-10): bring the dropped allele-specific-junction result (ASJ, June:
475 calls / 235 genes on `GGO.bam`, OR6737 testis, read-anchored, `project_allele_specific_junctions`) back as a sub-result of O2: re-run it on KB3781's own
fibroblast reads, phase reads with O2's rule instead of the June 0.8-agreement 2-means, and use mGorGor1's maternal / paternal assemblies (the same
animal) to tell an allele from a copy.

**Why O2.** O2 assigns a read to one of a set of sequences that differ at known columns (PSVs), or abstains. A haplotype is the same object one level down:
two sequences of one locus that differ at heterozygous SNVs. The ASJ test is then a downstream comparison of junction usage between the two assigned groups.
Junctions are NOT used to phase (they are what is tested).

## Inputs (fixed)

- Reads: `fibroblasts/GCA_029281585.2_flnc_mm.bam` (KB3781 fibroblast Iso-Seq on the primary, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`). Primary
  alignments only (no secondary / supplementary / unmapped); no MAPQ filter (paralogous reads are part of what is measured; per-call MAPQ-0 share reported).
- Genes: `GGO_genomic.gff` gene records on the 23 autosomes of `rna_allele/chrmap.tsv` (and chrX for the control), with >= 2 exons in some transcript and
  >= 20 primary reads overlapping the gene span.
- Haplotypes: each primary chromosome is byte-identical to one haplotype chromosome (`same_hap`); the other haplotype's chromosome is `B_name`
  (`gorilla_haps/{mat,pat}.fa`). Het sites of a gene = SNVs between the primary gene span (+ 1 kb flanks) and its alignment to `B_name`
  (`minimap2 -x asm5 -c --cs`, best alignment on `B_name` covering >= 90% of the span); SNVs within 5 bp of an indel in that alignment are dropped.
  chrX has no `B_name` (KB3781 is male): no het sites, no phasing.

## Arm A — assembly-phased, O2 rule (the new method)

Per read: at each het site it covers with a base (not deleted / not in an intron), the observation is the read base. **O2's certificate as in
`copy_pair_significance`** (`src/family.rs`): for the MLE haplotype `best` (log-likelihood, per column ln(1 - eps) if the read base equals that
haplotype's allele, else ln eps) vs the other, eps_j = 0.003 / 3 per distinguishing column (O2's `error_rate` default, no base qualities, editing filter
off), k = columns where the read matches `best`, p = P(X >= k), X ~ Poisson-binomial(eps). **Assign iff `best` is the strict MLE and p < alpha / (n - 1) =
1e-3** (n = 2); else abstain. One site never assigns (p = eps = alpha); two matching sites do.

## Arm R — read-anchored (June method, replicated on KB3781)

Per gene, candidate anchors from the reads (`--eqx`: X operations): a position with two alleles, each on >= 5 reads, minor fraction 0.35-0.65, depth >= 12;
the anchor = the candidate with the highest depth (ties: smallest position). Reads split by their base at the anchor.

## Junction test (both arms, June constants)

Per gene and group split: introns = the exact `N` operations of the reads, kept when >= 3 reads carry them. A read SPANS an intron when its alignment covers
both ends; PSI_g = reads of group g using it / reads of group g spanning it. Junctions with PSI = 1 in both groups are not tested. Fisher exact (two-sided)
on the 2 x 2 table; Benjamini-Hochberg over all tests of the arm (all chromosomes in the half scored); **ASJ = q < 0.05 and |dPSI| >= 0.3.**

## Allele vs copy (the haplotype assemblies as referee)

Each arm-R anchor is **ALLELIC** when it is an arm-A het site and the two read alleles are the primary and `B_name` bases; **NON-ALLELIC** when the `B_name`
alignment covers it (not within 5 bp of an indel) with the SAME base as the primary (homozygous in the animal: the second read allele comes from a paralog,
an edit or an error); **UNRESOLVED** otherwise (no qualifying `B_name` alignment, near an indel, or a third base).

## Tests

**Split.** Chromosomes by `chrmap.tsv` number: odd (1, 3, ..., 23) = DEV, scored and written first; even = HELD-OUT, computed after; at most ONE amendment
between them; the held-out half decides T1 and T2. Nothing is tuned: every constant is O2's default or the June ASJ value.

- **T1 (O2's bound holds on haplotypes).** Each gene's het sites split by site index parity; a read assigned by BOTH halves separately is concordant when
  the two agree. **Bar: discordance <= 0.005** over such reads (the certificate implies about 0.002 for two independent assignments at alpha 1e-3).
- **T2 (the assemblies separate allele from copy).** For each arm-R anchor with >= 10 reads assigned by arm A using the gene's OTHER het sites
  (anchor excluded), concordance = max over the two pairings of the share of assigned reads whose anchor allele matches the pairing. **Bar: median
  concordance >= 0.95 over ALLELIC anchors AND <= 0.80 over NON-ALLELIC anchors** (each with >= 10 anchors; fewer = UNDERPOWERED).
- **T3 (descriptive, all autosomes).** Arm-A ASJ calls (genes, calls, median |dPSI|); arm-R ASJ calls by referee class (the copy / edit-confounded share
  = NON-ALLELIC / resolved, beside June's 44 / 120); overlap of arm-R ALLELIC calls with arm-A calls (same gene and junction); per class the share of
  reads with MAPQ 0 and, among NON-ALLELIC anchors, transitions A/G vs others (editing vs paralog); abstention rate of arm A.
- **Control (descriptive).** chrX >= 5 Mb from either end (male, hemizygous, PAR excluded): arm-R candidate anchors per tested gene and arm-R ASJ calls,
  beside the autosomal rates; every such call is non-allelic by biology.

**Claim if T1 and T2 hold:** "O2's certificate assigns reads to haplotypes with its stated error bound (T1), and the animal's haplotype assemblies separate
allele-specific from copy-specific junction usage (T2); N junctions are allele-specific in KB3781, and X% of the read-anchored calls are copy- or
edit-confounded." If T1 fails, the bound does not transfer to haplotypes; if T2 fails, the referee is not established. Not claimed either way: a mechanism
(splice-site disruption; retracted for the June flagships), or anything about OR6737.

Code: `bench/asj_o2/` (test-first). Compute: per chromosome, foreground, `tools/rlock.sh heavy` for the alignment step.

## DEV result (2026-10-10; code ce1769f0 + two speed-ups shown byte-identical on chr21, scorer 23 tests; held-out NOT computed)

12 odd chromosomes, 6,807 expressed multi-exon genes (6,562 with an allelic alignment; median 27 het sites per gene).
**Arm A phasing:** 1,270,013 reads to the primary's haplotype, 1,230,894 to the other, 5,211,093 abstained (67.6%: reads covering < 2 het sites or tied).
- **T1 PASSES on dev: 872 discordant of 883,715 reads assigned by both halves = 0.00099** (bar 0.005; the certificate implies about 0.002).
- **T2 PASSES on dev:** median concordance **1.000 over 2,585 ALLELIC anchors** (bar >= 0.95) and **0.642 over 29 NON-ALLELIC anchors** (bar <= 0.80);
  UNRESOLVED 1.000 over 68.
- T3 (dev): **arm A 289 ASJ / 137 genes** (median |dPSI| 0.552; 69,633 tests). **Arm R 485 ASJ / 172 genes: ALLELIC 187, NON-ALLELIC 78, UNRESOLVED 220**;
  NON-ALLELIC = **29.4% of the resolved calls** (June, other animal and criteria: 44 / 120). 141 of the 187 ALLELIC arm-R calls are also arm-A calls.
  Anchors: ALLELIC 3,999, NON-ALLELIC 135, UNRESOLVED 822; transitions 0.75 of ALLELIC anchors vs **0.90 of NON-ALLELIC anchors** (an editing signature);
  median MAPQ-0 share 0 in every class (the non-allelic anchors are not low-MAPQ reads). Examples: ATXN3 two donors 7 bp apart switch with the allele
  (dPSI +0.89 / -0.97); AGTRAP two acceptors 21 bp apart (+0.68 / -0.68).
**No amendment.** The held-out half and the chrX control are run with the same code.
