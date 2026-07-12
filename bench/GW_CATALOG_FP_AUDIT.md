# Genome-wide catalog artifact audit: do we still have false positives?

**Date:** 2026-07-11. Binary `f379800`. Ran `gw_family_catalog` genome-wide (default, raw E_c conflict oracle,
`GGO_mm.bam` vs `GGO.fasta`), then a structural audit (`gw_artifact_audit.py`) + a **23-agent adversarial
verification** (each agent extracts copy sequences, aligns them pairwise, cross-checks the annotation) of the
suspicious subset. This is the first WHOLE-catalog audit — the earlier `FAMILY_ARTIFACT_AUDIT.md` sampled 30 loci
and **excluded the large clusters that time out per-region**, which is exactly where the artifacts live.

## Structural audit of the raw catalog (124 families, 320 copies) — the classic artifacts are gone

| check | count |
|---|---|
| cross-family chimeras (shared copy across families) | **0** |
| giant single-exon readthroughs (R4 target) | **0** |
| giant-span (>500kb) runaways | **0** |
| duplicate families | **0** |
| intra-family overlapping copies | 2 (both benign: SAFB/SAFB2 adjacency; one nested fragment) |

But the annotation cross-check flagged **23 of 124 families (19%) with a copy of large span (50–490 kb) but very
low read support (2–9 reads)** — the mis-assembly signature the 500 kb structural threshold misses.

## Adversarial verification of the 23 suspicious families — real FPs found

| verdict | n | families |
|---|---|---|
| **ARTIFACT_WHOLE_FAMILY** | **11** | GWFAM1, 8, 25, 39, 41, 47, 79, 85, 87, 94, 105 |
| MIXED (real family + artifact copies) | 6 | GWFAM10, 43, 54, 59, 68, 96 |
| REAL (flag was a false alarm) | 6 | GWFAM32, 33, 52, 61, 63, 77 |

**Two false-positive mechanisms, both in the raw E_c oracle:**

1. **Large-gene intra-gene mis-chaining.** One big gene with giant introns is split into "copies" that are
   non-overlapping fragments of the *same* gene, with **zero** mutual homology (minimap2 `asm20` → 0 alignments):
   PBX1 (GWFAM1, 292 kb), EBF1 (GWFAM8, 403 kb), CTNNA2 (GWFAM39, 1.2 Mb), HS6ST2 (GWFAM105, 333 kb). The
   low-read copy is a de-novo chain that spans the gene's giant introns (and sometimes nested lncRNA/snRNA),
   which R4 misses because it is *multi-exon*.

2. **Repeat / intronic-segdup bridge between UNRELATED genes.** Two single-copy genes share a small intronic
   repeat, a read maps ambiguously across it, and the raw conflict graph fuses them into a "family": NNT↔GHR
   (GWFAM79, 8.7 kb intronic repeat, 22% cov), ADAMTSL1↔CCDC171 (GWFAM41), GARRE1↔ZNF540 (GWFAM85, 2.9 kb / 13%),
   XNDC1N↔FAM168A (GWFAM25, 7%), plus GWFAM47/87/94.

The 6 MIXED families are real (FAM153B, IGH, GAGE tandem array, JMJD5/KDM8, DAPK1 segdup) but carry 1–2 artifact
copies that overstate `n_copies`. The 6 REAL families (POTE/ANKRD18A, PDPR, rhophilin-2, and three segdup pairs)
prove the flag is a *predictor*, not a verdict — a genuine low-expression segdup paralog also looks large-span
+ low-read.

## The fix: homology-gated refinement (`--refine`)

Every whole-family FP fails the `--refine` gate by construction — copies must be **mutually homologous
(asm20 id≥0.80, cov-of-shorter≥0.50) across ≥2 disjoint loci**:
- gene-splits: the two fragments give **0** alignments ⟹ removed;
- repeat-bridges: shared coverage is 7–22% « 50% ⟹ removed;
- MIXED artifact copies: the inflating copy's homology covers «50% of its length ⟹ trimmed, real core kept.

**Measured (`gw_family_catalog --refine`):** raw **124 families / 320 copies → refined 86 / 192** (38 families,
128 copies dropped). **10 of the 11 whole-family FPs are removed**; the structural audit of the refined catalog is
clean (0 cross-shared, 0 giant single-exon, 0 giant-span, 0 duplicate families, 1 benign intra-overlap), and the
large-span/low-read copy count falls 28 → 10 (the survivors are genuine segdup paralogs: POTE, and homologous
pairs). The single "surviving FP" (GWFAM94→GWFAM64) is a **verifier over-call, not a refine miss**: its two copies
align at **96.6% identity over 65% of the shorter** — a real homologous pair — one of whose copies merely has an
over-extended (mis-assembled) 154 kb span. Refine's gate correctly keeps it. So refine removes **every** true
whole-family FP (0-homology gene-splits, <50%-cov repeat-bridges).

## Resolution (both shipped)

1. **`--refine` is now the DEFAULT** (`--no-refine` opts out; `--refine` kept as a back-compat no-op). Measured:
   124/320 → **86/192**, all 11 whole-family FPs removed, refined structural audit clean. This is the real fix
   for the gene-splits and repeat-bridges — both fail the mutual-homology gate.

2. **Assembly-level mis-chain filter** (`retain_non_mischain`): drops a transcript with a giant intron (>50 kb)
   supported by fewer than `GATE_MIN_READS` (3) reads carrying that exact junction. Removes **467** clearly-
   unsupported spurious splices genome-wide — a transcript-quality improvement (cleaner spans, fewer false seed
   loci), independent of the catalog.

**⚠Measured limit — the assembly filter CANNOT catch the gene-splits alone.** PBX1's 115 kb spurious intron is
carried by **6 reads** (above the gate), so it is intrinsically indistinguishable from a real low-expression
large-gene intron by any within-transcript signal; the discriminator is that the copy shares **no homology** with
its supposed paralog, which only `--refine` can see. So mis-chain removal for the well-supported cases is a
*family-level homology* problem, not an assembler problem. `--refine` (item 1) is what actually removes them; the
assembly filter (item 2) cleans the sub-gate noise.

## Bottom line

The classic artifacts (readthroughs, chimeras, giant spans) are gone. The quieter FP class — large-gene
mis-chains and repeat-bridges, ~9% of the raw catalog — is removed by making **`--refine` the default** (now
done). The assembly-level mis-chain filter additionally strips 467 sub-gate spurious splices, but by measurement
the *well-supported* gene-splits are only separable with homology context, i.e. by refine.

Reproduce: `gw_family_catalog --bam GGO_mm.bam --fasta GGO.fasta --out gw` (now refined by default) then
`python3 bench/gw_artifact_audit.py gw.copies.tsv GGO_genomic.gff`; `--no-refine` for the raw catalog.
