# O3 port (Strainy's phasing + IsoCon's statistics): results on gorilla KB3781 (2026-10-23)

**Status: pre-registered** (`docs/PREREG_o3_port_2026-10-23.md`, ddc1b4b7). The code was frozen for P3 at 1b4535cc.
**Verdict:** P1 fails, P2 passes and P3 H1 passes, so the port is **not adopted** as the O3 stage: "detects but not sensitive"
(H3 holds). Design: `docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md`. Code: `bench/o3_graph/o3port.py`
(11 tests).

## What the frozen port does
- **Strainy side**, per reference copy:
  1. **Piles:** a read joins every copy it aligns to within 10 AS of its best.
  2. **Informative columns:** >= 2 alleles with >= 2 carriers each, the minor allele significant under the reads' own empirical
     error rates (IsoCon), Bonferroni, alpha 0.01.
  3. **Clusters:** a zero-disagreement read graph, label propagation and per-cluster refinement.
  4. **Haplotypes:** the largest set of haplotypes that differ pairwise at >= 2 pile-level columns.
  5. **Chimeras and duplicates** (DADA2's rule, before slots): a haplotype that is a left/right join of two parents takes no
     slot, and one equal to a more abundant haplotype shares that haplotype's slot. Parents are every family copy's reference
     and the more abundant haplotypes, mapped across copies by asm20 homology. A haplotype equal to a copy's reference is that
     copy's allele.
  6. **Certificate:** allele slots (2 per copy) filled reference-like first by b-matching. A haplotype is certified when it gets
     no slot, carries >= 2 private alleles and has >= 6 reads (the shipped stage's own flag floor).
- **IsoCon side:** the shipped stage's reference-free candidates, kept when most of their reads reached the net without
  aligning to any copy.

## Results (the same reads and the same scorer as the shipped stage K)
| rule | port | K (`candidate_copies`, Amendment 15) | bar | port verdict |
|---|---|---|---|---|
| P1 deletion panel (53 families, dev): D right | 3,995 | 7,882 | >= 5,842 | **FAILS** |
| P1: false moves | 7 / 41,727 = 0.02% | 0.77% | <= 5% | holds |
| P2 no-deletion control (53, dev): families with a false flag | 0 | 10 | <= 8 | **PASSES** |
| P3 held out (30 families, frozen code): H1 families with a false flag | 3 | 6 | <= 4 | **PASSES** |
| P3: H3 false moves | 178 / 27,757 = 0.64% | 0.26% | <= 5% | holds |

- **Held-out flags.** The port flagged 6 haplotypes in 4 families:
  - 1 is a real reference-absent copy: GWFAM390, 8 reads, 100% identical to a maternal-only locus (CM054598.2:30.25 Mb).
  - 3 are B-haplotype alleles of existing copies: GWFAM205 and GWFAM256, at 99.94-99.96%.
  - 2 are unmatched: GWFAM256 (6 reads) and GWFAM318 (IsoCon side).
- **On P1.** The port gives back half the stage's deleted-copy reads, with about 38x fewer false moves.
- **On H3.** 362 held-out reads land on their own candidate. Many of the 178 false moves are reads of the existing copies' B
  alleles, which the b_allele contigs attract.

## How v1 became v4 (development on P1 / P2 only; every change uses an existing constant)
| version | change | P2 families flagged | P1 families with a D candidate |
|---|---|---|---|
| v1 | base design | 35 (97 flags) | 16 |
| v1 + empirical error | IsoCon's per-read error rates, not PacBio QVs (capped at Q93) | | |
| v1.1 | conflicts proven only on pile-level columns | 23 | |
| v1.2 | >= 2 private alleles | 22 (C1' classified: 45 c, 5 b, 2 pri, 1 a) | |
| v2 | family-wide chimera rule before the slots | 19 | 13 |
| v3 | the stage's flag floor (6 reads); slot order by divergence (support order: 396 S reads certified on P1) | 0 | 8 |
| v4 | + IsoCon side (via-B reads) | 0 | 16 (P1 harness: D right 3,995) |

## What we learned
1. **Full-length reads were not the problem.** 68 of 72 surplus haplotypes were made of near-full-length reads. They were not
   A-to-I editing either: their alleles are transition-rich and none was all A>G.
2. **RNA piles are full of small artifact clusters** (2-7 reads): chimeric molecules, some joining paralogs, and reads that fit
   poorly. They match neither parent (median identity 0.991). Without a support floor, every allele-slot certificate drowns in
   them.
3. **The certificate is honest but conservative (Theorem A in practice).** A haplotype whose reads also tie to a copy with a free
   slot cannot be proven extra. That is how GWFAM175's real paternal-only copy (181 reads) was lost on development. On held-out
   data the port still found GWFAM390's maternal-only copy, with 1 extra false family fewer than the bar allows.
4. **Read support and divergence answer different questions.** Filling slots by support certifies the survivors' own minor
   alleles; filling them reference-like first certifies the collapsed copy.
5. **Sensitivity is lost in two places.**
   - Deleted copies whose reads tie to a survivor with a free slot.
   - The 31 of 53 deleted copies whose reads never touch a survivor. Those depend on the IsoCon side, and that side is only as
     sensitive as the stage's reference-free clusters.
