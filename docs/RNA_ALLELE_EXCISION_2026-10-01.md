# Reference-absent copies from RNA: unmapped reads, the maternal-only NPIP copy, and the NPIPA2 deletion test (Amendment 4), 2026-10-01

Work dirs `/mnt/linuxdisk/tmp/rna_allele/isocon/`, `/mnt/linuxdisk/tmp/rna_allele/excise/`; scripts `bench/rna_allele/excise_readout.py`
(R1-R3); R4 and the read-assignment check were run inline (commands in the session log). KB3781 fibroblast Iso-Seq; minimap2 2.30.

## Unmapped reads (descriptive)
959 of 13.5 M fibroblast reads are unmapped; median length 69 bp. 109 touch an NPIP/TBC1D3 transcript, all at 0.85-0.95 identity over
100-280 bp: too short to place. No evidence of a missing family copy in the unmapped reads.

## Mapped reads that belong to a copy the reference lacks (descriptive)
Of the 2,951 NPIP-net reads, 3 match KB3781's maternal-only NPIP copy (chr18 mat 17.40 Mb) at 99.8-100% over their full 1.7-4.7 kb,
while their best reference placement is 93-96% (chr18). The signal of a missing copy is reads placed with a few percent of mismatches.

## NPIPA2 deleted from the reference (Amendment 4)
- **R1 fate:** 176 reads had their primary on NPIPA2; after masking 171 (97.2%) land on LOC124907808 (same chromosome, 280 kb away), 5
  scatter, 0 become unmapped. Median `de` 0.0009 -> 0.0020.
- **R2 fake allele:** none called at LOC124907808. The absorbed reads outnumber the copy's own: at 4 positions ~96% of reads disagree with
  the reference and the copy's own base drops to 2-6%, below the 0.20 minor-allele floor. A missing copy here HIDES the real one rather
  than posing as its allele. (In KB3781 a near-fixed mismatch against the assembly of the same animal is itself a flag; in another
  individual it reads as a homozygous SNP.)
- **R3 no-reference-match groups:** none possible: the registered paralog test (minimap2 default -p 0.8) lists no paralog for
  LOC124907808 (its hit to NPIPA2, 0.93 identity, is reported only from NPIPA2's side). Disclosed limitation of the paralog test.
- **R4 IsoCon:** 107 outputs on the re-selected net; 16 are NPIPA2 transcripts (>= 0.999 to NPIPA2 in the unmasked genome). Against the
  masked genome 35 outputs are flagged "not in the reference" (< 0.999): **15 of the 16 NPIPA2 outputs** (their best remaining match is
  LOC124907808 at 0.9982-0.9989) and 20 outputs that are flagged against the unmasked genome too. Masking adds 15 flags, all true.
  The 20 background flags: 5 match a KB3781 haplotype at >= 0.999 but not `_pri` (4 from maternal-only sequence at chr18 mat 6.72-6.76
  Mb, 79-93% to the reference: a second reference-absent expressed locus; 1 a maternal allele), 14 sit at 0.99-0.999 everywhere
  (consensus error or editing), 1 < 0.99. A "not in the reference" flag does not say whether it is a missing copy or an allele.

## IsoCon as an O2 assigner (simulation, read origin known; `cluster_info.tsv` read -> candidate)

| | right copy, right haplotype | right copy, haplotype ambiguous or wrong | wrong copy | tied between copies | not assigned |
|---|---|---|---|---|---|
| TBC1D3 (560 reads) | 30.2% | 24.6% | **24.8%** | 8.9% | 11.4% |
| NPIP 4 kb (880 reads) | 34.1% | 28.9% | 4.0% | 0% | 33.1% |

Among assigned reads the wrong-copy rate is 28% (TBC1D3) and 6% (NPIP): IsoCon's grouping is not an assign-or-abstain certificate.
