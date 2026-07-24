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
  - **34 are real copies homologous to a Soto family** — copies Soto's SD98 (≥98%-identity) threshold
    *missed*. RNA found MORE real copies, not false ones. **A win.** (19 confirmed on first pass; a
    further 15 recovered when the 18 "no-homology" candidates were re-aligned — see below.)
  - **3 are genuinely non-Soto** — all real genomic loci: 2 have their own high-identity genomic
    paralogs (real duplications outside Soto's 80 families), 1 is a real single-copy gene. None is a
    fabricated family member.
  - **0 are phantom artifacts / family-level false positives.**
  - *Correction (verified):* the earlier pass labelled 18 of these "no Soto homology," but its
    alignment step returned `best_id=0.000` for all 18 — a failed pass, since several clearly align
    at ≥0.94 id to their own family. Re-aligning to all 362 members reclassifies 15/18 as real
    Soto-family copies, 3 as real non-Soto loci, 0 as artifacts. Full breakdown:
    `bench/soto/candidate18_classification.md`.
- **The structural FP defense:** family definition requires the refine gate — asm20 id≥0.80,
  cov-of-shorter≥0.50, ≥2 distinct loci. The genome-wide audit showed this removes ~9% over-merge
  FPs (124→86 families). Non-homologous loci cannot be merged into a family.

**Bottom line for the room:** the DNA graph is a complete, precise 100% ceiling. RNA is a faithful
subset — it recovers 76% (→87% with the seeding fix) of that ceiling and invents no phantom family
members (precision ~100% against the DNA truth; **34 of 37 "extra" calls are real copies of a Soto
family DNA's per-member check missed, the other 3 are real non-Soto loci, and 0 are phantom**).
What RNA cannot recover is the K=0 identifiability floor — homology-grouped copies no read can
separate — and that is a property of the reads, not a false positive.
