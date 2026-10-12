# PREREG — phased-data toy (2026-10-12)

**Question (collaborator, via the user):** would phased data help? Toy proof of principle only — one real gene on a
very small real chromosome; no claim about real gorilla data.

**Substrate.** CHM13 chr20:2,481,001-2,511,000 (30 kb, real sequence) containing SNRPB (CAT transcript
CHM13_T0135794, 7 exons, minus strand). Copy unit = transcript span +-500 bp (~10.2 kb). Variant set V = 5 exonic
SNVs + intronic/pad SNVs at rate 0.005 (~0.5% divergence, between gorilla nucleotide diversity and the KB3781
p99 exonic haplotype divergence .00958). Seed 7.

**Scenarios (two haplotypes each):**
- **S0 control** — hap2 = hap1 (homozygous).
- **SA allele** — hap2's only copy carries V (heterozygous SNVs). Reads: 50 hap1, 50 hap2.
- **SB copy** — hap2 carries an EXTRA copy B = unit + V inserted at chr20:2,506,000 (dispersed, unique spacer);
  hap2's own copy = hap1's. Reads: 25 hap1-A, 25 hap2-A, 50 hap2-B.

The same V is used in SA and SB and the per-sequence read counts match (50 reference-sequence, 50 V-sequence), so a
collapsed reference receives identically distributed reads in SA and SB **by construction** (allele fraction is not
used as a copy-number signal: allele-specific and copy-specific expression make it unreliable in RNA).
Reads = spliced mRNA, 5'/3' ends jittered 0-30 bp, substitutions 0.001, indels 0.0003 (famsim's settings).

**Arms.**
- **Collapsed** — map all reads to hap1 only (`minimap2 -ax splice:hq -uf`). Report the variable columns
  (alt >= 20% of >= 10 covering reads) and the alt-read fraction. This arm has no rule that can call allele vs copy.
- **Phased** — map reads to hap1 and hap2 separately (as with mat/pat). Truth-free rule per read:
  NM1 == NM2 -> SHARED; else home = haplotype with lower NM, home locus = merged primary spans there;
  ALLELE if the home locus +-2 kb aligns to the other haplotype in ONE alignment (`minimap2 -cx asm10`) covering the
  query end to end (qstart <= 100, qend >= qlen - 100) with a target span within 5% of the query span (the same locus
  between the same flanks, no large indel), else EXTRA_COPY. The 2-kb flank exceeds the unit's 500-bp pad, so it
  reaches unique sequence; on real data the flank must clear the whole duplicated unit (a toy limit).
- A read is **informative** if its window covers >= 1 exonic V position (truth, from the simulator).

**Predictions and bars (decided before any run):**
- **P0** S0: <= 2% of reads non-SHARED.
- **P1** collapsed SA vs SB: identical variable-column sets, alt fractions within 0.10 (the information boundary,
  `NEGATIVE_RESULTS_REGISTER.md:767`, made concrete — true by design, reported to show it, not as a finding).
- **P2** phased: SA >= 95% of informative reads ALLELE and <= 2% EXTRA_COPY; SB >= 95% of informative copy-B reads
  EXTRA_COPY and <= 5% of all reads ALLELE.
- **P3** non-SHARED reads in SA+SB assigned to their true haplotype >= 98%.

**Not claimed.** Real data has unequal expression, larger and tandem duplications (where "which copy is the extra
one" is undefined), gene conversion, K >= 3 copies (`NEGATIVE_RESULTS_REGISTER.md:561`), and needs a
haplotype-resolved assembly of the same animal. Code: `bench/phased_toy/phased_toy.py`; outputs
`/mnt/linuxdisk/tmp/phased_toy/`.
