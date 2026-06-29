# POA-core completion: reaching loosely-related paralogs the read-conflict graph misses (2026-06-29)

**The two-tier family definition.** The read-conflict graph answers *where is assignment needed?* — it links only
copies a read **confuses** (empirically down to ~87% identity; below that reads resolve the copies and raise no
edge). That deliberately excludes **loosely-related paralogs** (divergent copies at other loci). The POA-core
completion is the second tier: *given that read-conflict has determined a family*, complete its roster with the
divergent homologs by sequence — so the expensive homology step runs only **when needed**, not genome-wide.

**Mechanism (`family_detect::poa_core_completion_adds`).** For each read-conflict family, attach any genome rep
whose **contiguous POA core** to a family member clears `t_core` (0.13 — a ≥13%-length contiguous identical block,
i.e. a shared conserved exon even if the flanks diverge heavily). Bounded and seeded:
- the **minimizer-LSH prefilter** (`candidate_pairs`) restricts grading to homologous candidates;
- only pairs with **exactly one endpoint in a conflict family** (the other a free rep) are graded — free-vs-free
  pairs would be NEW families, not seeded by read-conflict, so they are skipped;
- the core is confirmed with the **linear-memory longest-common-substring** metric (the documented faithful
  equivalent of the poasta "longest ungapped equal run") rather than poasta itself — O(n) per pair, so it stays
  cheap. (Using poasta here was the first cut and took ~99 min genome-wide; the LCS swap fixes it.)
- a free rep matching several families attaches to the one with the strongest core.

**Wiring.** `DenovoConfig.complete_poa_core` (default **false** → no-op, byte-identical), exposed as
`gw_family_catalog --complete-core` (cross-chrom path). Added copies are appended after the conflict core and
logged. TDD: `poa_core_completion_attaches_a_divergent_paralog_at_a_new_locus` (attaches a divergent paralog that
shares the conserved core; rejects an unrelated rep and free-vs-free pairs). Full lib suite green.

## Real-data measurement (cross-chrom, GGO_mm.bam)

- **Read-conflict cross-chrom catalog:** 265 families / ~869 copies (reach ceiling ~87% identity).
- **+ POA-core completion (`--complete-core`, t_core=0.13):** **+55 copies attached to 35 families** → 265
  families / 924 copies (baseline 869 + 55, the diff is exactly the 55).
- **Reach of the attached copies (the honest result):** measured against their family core — asm20 54/55 aligned
  median 100% min **95%**; sensitive `-k11 -w5` 55/55 aligned median 100% min **93%**. So the 55 adds are
  **NEAR-IDENTICAL (93–100%), NOT loosely-related.**

## What this means (the important finding)

The mechanism is correct (the TDD test attaches a divergent-flank paralog that shares a conserved core), but on
**real GGO the nucleotide completion does not add loosely-related paralogs** — it recovers near-identical copies
the read-conflict graph **missed for COVERAGE reasons** (too few reads to clear the `min_reads≥3` conflict-edge
floor), which is a useful completeness gain but a different thing.

**Why no loosely-related adds, and it's illuminating:** read-conflict ALREADY links divergent-flank paralogs via
their **conserved-exon reads** (a read over the shared exon ties the two loci even if the flanks diverge) — that
is exactly why read-conflict reaches ~87% *overall* identity. So the conserved-exon-sharing divergent paralogs
are *already in the families*; the only copies left for a nucleotide completion are the near-identical ones missed
for coverage. A genuinely loosely-related paralog (divergent throughout, no conserved-exon read overlap) does NOT
clear `t_core=0.13` on the exact-LCS core either (30% scattered divergence breaks any long identical run) — it is
reachable only at the **protein level** (the tier deferred in this build).

**Net:** all NUCLEOTIDE family definitions on real GGO top out at ~81–87% identity:
- read-conflict ~87% (via conserved-exon ties);
- `--refine` validates within (asm20 ~81% floor; 422 copies, 2 below 87%);
- `--complete-core` adds 55 near-identical (93–100%) missed copies — completeness, not extended reach.

**Genuinely loosely-related paralogs need the PROTEIN tier** (6-frame ORF → mmseqs, fident≥0.50; reaches ~50%
protein / ~70% genomic) — already implemented as the opt-in `--protein-tail` for `--refine`, and the natural next
lever for the completion pass if loosely-related reach is the goal. (Or DNA/SEDEF for ancient paralogs.)

## How the tiers relate (corrected)
- **read-conflict** (tier 1): identifiability — which loci need assignment; ~87% via conserved-exon ties.
- **`--complete-core`** (nucleotide, to new loci): recovers near-identical copies missed for coverage; does NOT
  extend reach on real data because the divergent-but-core-sharing copies are already in tier 1.
- **`--refine` / `--protein-tail`** (homology within / protein): the protein tier is the only one that reaches
  genuinely loosely-related (70–87%+) paralogs; nucleotide cannot.
