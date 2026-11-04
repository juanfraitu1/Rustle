# PREREG PSV PILE AT ABSORBING LOCI (register 1538, 2026-11-05, before any number)

Question (user 2026-11-05 "lets measure the PSV-pile row", after r1537): a collapsed copy's reads are absorbed on a near-identical
relative at that relative's identity (no per-read hallmark). Do they leave a POPULATION signal at the absorbing locus: a second
haplotype of paralog-specific variants (PSVs) that the relative's own reads lack, above the heterozygosity every present locus shows?
Data: r1537 products (human only: the collapsed regime; gorilla's absent regime already separates on identity):
`/anvme/workspace/v138da13-rustle/nont2t/human/{nt.bam,hm.reads.tsv,read_origin.tsv}`, GRCh38 `~/annot/nont2t/grch38.fa`.
Loci (by read landing, after the r1537 correction): for every expressed panel copy (missing and present), cluster its reads'
primary landings on GRCh38 (same chromosome, gaps <= 10 kb merged); the copy's LOCUS is the cluster with most reads, and the
copy is scorable when that cluster holds >= 50% of its reads and >= 10 reads. A missing copy's locus is an ABSORBING locus; a
present copy's locus that overlaps no absorbing locus is a NULL locus (ordinary heterozygosity and present tandem copies).
Calling is blind (no origin labels): `bcftools mpileup -f grch38.fa -r LOCUS -Q 10 -q 0 --ff UNMAP,SECONDARY,SUPPLEMENTARY
-a AD,DP -d 100000 nt.bam` on all panel reads, then sites with DP >= 20 and a second allele with >= 5 carriers at fraction >= 0.10.
Per locus, blind statistics: S = such sites per kb of locus with DP >= 20; H = share of reads covering >= 2 sites that carry the
minor allele at >= 80% of the sites they cover (a second haplotype), among all reads covering >= 2 sites.
Labelled check (scoring only): purity = share of minor-allele carriers whose origin is the missing copy (reads of the locus's own
present copy give the heterozygosity part).
Hallmark: an absorbing locus is PSV-hallmarked when S or H exceeds the 95th percentile of the null loci.
Bars (as r1537): PSV-hallmarked >= 50% of scorable missing copies, split collapsed vs absent, AND null loci hallmarked <= 25%.
Stop rule: PASS -> a blind detector (site density + second-haplotype share at family loci) is the registered next row; FAIL -> the
collapsed regime closes on sequence data entirely (reads and population), and collapsed paralogs need a better assembly.
Script: `bench/nont2t/psv_pile.py` (loci, bcftools, per-locus table, summary); bcftools from the pixi env `bcftools-*`.

## Amendment 1 (2026-11-05, before any number)
1. SNV records only; bcftools indel records are ignored (HiFi homopolymer indels would set S at every locus).
2. Ploidy strata: A119b is male testis, chrX/chrY loci are haploid (no heterozygosity). Thresholds and bars are computed and
   reported separately for autosomal and X/Y loci; if the X/Y null holds fewer than 5 loci its threshold is the null maximum.
3. Unit: the bar is per copy as registered; the same count is also reported per unique absorbing locus (several copies share one:
   NPIPB15 x2, TSPY3 x several, NPIPB7 + NPIPB9).
4. Floor by construction: the minor allele needs fraction >= 0.10, so a copy whose reads are < ~10% of the pile cannot pass
   whatever its PSVs (NPIPB11, NPIPB7 at depth excess 0.02); the count of such copies is reported.
5. Predicted outcome, written now: r1537's de gap bounds the divergence between a missing copy and its absorbing locus at
   <= 0.002 for 19 of 21, i.e. <= ~2 PSV sites per kb, the human heterozygosity rate. Expected: FAIL on autosomal loci by that
   bound; the informative parts are the haploid stratum and the labelled purity. The collapsed/absent split is kept as registered
   but an "absent" copy with de gap ~0 is absorbed on a relative the truth alignment did not credit.
6. Calling: `bcftools mpileup -a AD,DP` parsed directly, no `bcftools call` (its genotype model is wrong for a pooled pile); every
   position is kept so the DP >= 20 denominator comes from the same pass. Per-read alleles for H and purity come from pysam
   pileup with the same filters. Purity is defined on the allele carried by fewer reads and reports which origin carries it
   (at CDY1 the missing copy is the majority, depth excess 6.0). The 1,751 ambiguous-origin reads enter the blind pileup and are
   excluded from the purity denominator.
