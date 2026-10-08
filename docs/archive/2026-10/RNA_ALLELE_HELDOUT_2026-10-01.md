# Held-out test: IsoCon's reference-absent transcripts as extra copies (Amendment 6), 2026-10-01

Prereg: `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 6 (commit 548f6767, before any run). Substrate never used to
develop Amendment 5: the 162 two-copy gorilla families of the 2026-08-14 whole-genome excision (one copy D hard-masked per family, the
other K kept; KB3781 fibroblast Iso-Seq). Scripts `bench/rna_allele/heldout_prep.py`, `heldout_score.py`; work dir
`/mnt/linuxdisk/tmp/rna_allele/heldout/` (`score.out`).

- Scored reads: 51,165 D + 52,871 K (<= 500 each per family, seed 1; 114 reads sampled under both copies of a family were counted once).
- IsoCon on each family's net (median 530, max 1,000 reads): 5,172 outputs; 3,644 not in the masked genome at identity x coverage
  0.999 -> contigs; 1,105 of them "is D" (>= 0.999 to D in the unmasked genome), from 114 families.

## Registered result

| | R (masked reference) | R+I (+ IsoCon contigs) |
|---|---|---|
| D reads wrong | 30,785 | 13,526 |
| D reads right (an "is D" contig) | 0 | 1,796 |
| D reads unplaced (unmapped or AS tie) | 20,380 | 35,843 |
| K reads staying on K | 52,225 | 15,172 |
| K reads unplaced | 645 | 35,182 |
| K reads moved to a contig that is not "is K" (false moves) | 0 | 2,516 (4.8%) |

**Rule: wrong D -56.1%, false moves 4.8% -> HELP.** 88 of 162 families individually lose >= 50% of their wrong D calls. By fate
(2026-08-14): absorbed families wrong 28,900 -> 12,616; orphaned 1,695 -> 705, right 0 -> 626; the 4 scattered families do not improve.

**What the rule does not count:** most of the drop in wrong D calls became abstentions, not right calls, and 67% of K reads became
abstentions. The cause is redundancy: a read ties between a reference locus and a contig made from the same locus.

## Post hoc (not the verdict): by primary placement, contigs labelled by the locus of their best unmasked hit

All flagged contigs but one come from a family's own copies: 2,249 from D's locus, 1,394 from K's locus. They are the copies' own
transcripts that miss the masked reference at 0.999 (the other haplotype's allele, ends, or the deleted copy itself).

| | R | R+I |
|---|---|---|
| D reads on D-derived contigs | 0 | **38,281 (74.8%)** |
| D reads on K's locus | 13,260 | 754 |
| D reads unmapped | 18,168 | 466 |
| D reads elsewhere (absorbed by loci outside the family; not in the net, as declared) | 19,737 | 11,479 |
| K reads on K's locus or K-derived contigs | 52,869 | 52,814 (99.9%) |
| K reads on D-derived contigs (true cross-copy moves) | 0 | **47 (0.09%)** |

## Reading

- **O3 (detect + flag a copy the reference lacks): supported on held-out data.** Family-scoped reference-free transcripts recover the
  deleted copy for 114 of 162 families, give 75% of its reads a home, and cut the reads forced onto the surviving copy by 94% and the
  orphaned (unmapped) reads by 97%, while moving 0.09% of the surviving copy's reads.
- **O2: only the copy set improves.** The tie-abstentions show the contigs must be linked to the locus they come from (a contig that
  matches a reference locus at >= 0.99 is an allele of that locus, not a new copy) before reads are assigned; otherwise every read of a
  locus with a contig becomes ambiguous. Reads absorbed by loci outside the family (22% of D reads) stay misplaced: the read net is
  family-scoped.
- Caveats: one individual (the reference animal), fibroblast expression only, IsoCon capped at 1,000 reads per family.
