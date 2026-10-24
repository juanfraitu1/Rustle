# Phased-data toy — results (2026-10-12)

Prereg `docs/PREREG_phased_toy_2026-10-12.md` (c5431f41, committed before the run). Code `bench/phased_toy/phased_toy.py`
(seed 7, 1.1 s, 18 MB), outputs `/mnt/linuxdisk/tmp/phased_toy/` (`results.json`, `S*.calls.tsv`, haplotype FASTAs, reads).

Substrate verified by the script: 30-kb CHM13 chr20 slice, SNRPB unit 10,217 bp, V = 44 SNVs (5 exonic; 0.43%), mRNA 1,008 bp.

| scenario | truth | collapsed ref (hap1 only) | phased (hap1 + hap2, mapped separately) |
|---|---|---|---|
| S0 | homozygous | 0 variable columns | 100 SHARED |
| SA | heterozygous allele | 5 variable columns (= the exonic V), alt fraction 0.50 | 100 ALLELE (loci on both haplotypes, partner found) |
| SB | extra copy on hap2 | the SAME 5 columns, alt fraction 0.53 | 50 SHARED (copy A) + 50 EXTRA_COPY (hap2 locus, no partner on hap1) |

P0 0.00 non-SHARED (bar <= .02) PASS · P1 same columns, |diff| 0.03 (<= .10) PASS · P2 SA 1.00 ALLELE / 0 EXTRA, SB 1.00
EXTRA_COPY / 0 ALLELE PASS · P3 150/150 reads on their true haplotype PASS.

SB's 0.53 (not 0.50): three copy-A reads carry a simulated sequencing error at a V column (replayed: the injected
substitutions land on V positions 341/603/676); sampling noise, not copy-B reads.

**Reading.** The collapsed reference receives the same evidence in SA and SB — it can report "a second sequence at 5
columns" but not whether it is an allele or a copy. With both haplotype assemblies, a truth-free rule (which haplotype
explains the read better; does that locus sit between the same flanks on the other haplotype) separates them
perfectly. The toy is the easy case by design: dispersed duplication with unique flanks, 0.43% divergence, 0.1% read
error, perfect assemblies, equal expression, K = 2. Not tested: tandem duplications, gene conversion, K >= 3,
assembly errors, low divergence (PSVs < reads' error budget).
