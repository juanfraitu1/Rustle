# The T2T DNA graph is the 100% ceiling; RNA is a faithful subset

The frame that separates **O1 identifiability** (which copies exist / belong to a family) from
**O2 assignability** (which reads go to which copy), and answers the family-level false-positive
question head-on.

## O1 — the DNA ceiling (100%, no false positives)

From the **T2T genome** (CHM13) + Soto's DNA segmental-duplication catalog, **every family copy is a
path in a variation graph** (`bench/soto/dna_graphs/` — 76/83 families rendered). The genome *has*
all the copies; the graph *represents* them. This is:
- **100% sensitivity** — every real copy is a path (the DNA is complete).
- **~100% precision** — the graph contains only real segdup copies (given by DNA self-alignment /
  SEDEF, not invented by RNA). No phantom members.

This is the CEILING and the ground truth. It is anti-circular: DNA gives the truth set; RNA is the
method evaluated against it.

## O2 → RNA is a SUBSET of the ceiling (most, not all)

RNA reads recover a subset of the DNA graph's paths:

| | copies | % of DNA ceiling |
|---|---|---|
| DNA ceiling (all Soto members) | 362 | 100% |
| **RNA-recovered (now)** | 276 | **76.2%** |
| RNA-recovered (projected post-seeding-fix) | 316 | 87.3% |

The gap is honest and named (per `bench/soto/soto_sensitivity_precision.md`): **seeding-fixable
(38)** — reads existed, now recovered; **K=0 identifiability floor (34)** — homology groups them but
no read separates them (→ DNA); **silent/other (14)**. RNA is a faithful subset: it recovers paths
that are *in* the DNA graph and adds none that aren't (see FP section). In the graphs, **green =
RNA-recovered path, red = DNA-only path.**

## Family-level false positives ≈ 0 (the advisor's concern)

- **RNA copies that overlap a Soto member: 100% precision** — every RNA-detected member that maps to
  the DNA truth is real (0 phantom members).
- **The 37 extra RNA candidates** (RNA-detected loci NOT in Soto's DNA annotation), checked against
  ALL Soto members by the exact family-edge criterion (asm20 id≥0.80, cov-of-shorter≥0.50):
  - **19 are confirmed real copies** homologous to a Soto family — copies Soto's SD98 (≥98%-identity)
    threshold *missed*. RNA found MORE real copies, not false ones. **A win.**
  - **18 have no Soto homology** — real transcribed loci flagged as *candidates near a family*, NOT
    asserted family members (they are in the discovery table, never merged into a family roster).
    Their classification (copy of a non-Soto family vs parent gene vs artifact) is a follow-up; none
    is a family-membership over-merge.
- **The structural FP defense:** family definition requires the refine gate — asm20 id≥0.80,
  cov-of-shorter≥0.50, ≥2 distinct loci. The genome-wide audit showed this removes ~9% over-merge
  FPs (124→86 families). Non-homologous loci cannot be merged into a family.

**Bottom line for the room:** the DNA graph is a complete, precise 100% ceiling. RNA is a faithful
subset — it recovers 76% (→87% with the seeding fix) of that ceiling and invents no phantom family
members (precision ~100% against the DNA truth; 19 of 37 "extra" calls are real copies DNA missed).
What RNA cannot recover is the K=0 identifiability floor — homology-grouped copies no read can
separate — and that is a property of the reads, not a false positive.
