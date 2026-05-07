# Per-copy HMM-EM for paralog read assignment — meeting notes (2026-05-06)

Seven slide-ready PNGs in this directory. Suggested order + a few sentences each.
The narrative arc directly answers two advisor concerns:

1. *"The VG just averages information across copies — paralogs collapse."* → Slides 2, 3, 4 show **explicitly** that we keep one profile per copy, never a single consensus.
2. *"This will only work for near-identical paralogs."* → Slide 7 shows recovery at jaccard 0.52 (squarely in the medium-similarity band).

Plus one slide (6) addresses the FLNC misunderstanding directly.

---

## Slide 1 — `01_problem_setup.png` *(motivation)*

**One-line:** an FLNC multi-mapper has the same alignment score against several paralogs — and the BAM aligner just splits weight evenly (1/NH).

**Key points:**
- A read maps to N paralogs → BAM hands us N placements with weight 1/N each.
- Treating that as ground truth has two failure modes:
  - All copies look equally expressed even when only one is on.
  - Real paralog-specific isoforms get diluted below the assembly threshold.
- We need a principled way to decide: which copy did this read actually come from?

---

## Slide 2 — `02_collapse_vs_per_copy.png` *(directly addresses advisor #1)*

**One-line:** a single shared profile averages SNP signal away — per-copy profiles preserve it.

**Key points:**
- Toy MSA: 3 paralogs A/B/C, 4 columns; column 2 is a paralog-distinguishing SNP (A=C, B=G, C=G).
- **Left, BAD:** consensus emission at column 2 is uniform-ish (≈2/3 G, 1/3 C) for every paralog. A read with G at col 2 scores the same against A, B, C — no discrimination.
- **Right, GOOD:** per-copy emission gives A=0.02 for G, B=0.98, C=0.98. The SNP signal is preserved, EM has something to work with.
- This is the answer to "the VG averages everything." We made it not.

---

## Slide 3 — `03_family_graph_per_copy.png` *(architecture)*

**One-line:** family-graph nodes share **structure** (which exons are orthologous, which junctions exist) but each node carries one profile **per** contributing copy — sequence is never collapsed.

**Key points:**
- Nodes = exon-classes E1…E5 (orthologous exons grouped). Some classes have all 3 paralogs, some are paralog-specific (E2 in A/B only, E4 in B/C only).
- Inside each node: a stack of per-copy profile lozenges, color-coded by copy.
- Edges = junctions; the dashed arcs show paralog-specific path topology (C skips E2, A skips E4).
- "Information sharing" lives at the level of structure (which class is which) and the POA-MSA prior used for smoothing — not at the per-base emission level.

---

## Slide 4 — `04_poa_tilt.png` *(the tunable trade-off)*

**One-line:** the per-copy profile is `α · singleton + (1−α) · POA-MSA`. α=1 is too harsh, α=0 is information-collapse, α=0.7 is the sweet spot.

**Key points:**
- α = 1.0 → pure singleton, sharp paralog signal, but per-base error costs ≈5 nats. Hostile to noisy reads.
- α = 0.0 → pure shared MSA, smooth, but identical for every paralog → INFO COLLAPSE.
- α = 0.7 (default) → keeps the dominant-base concentration at the copy's actual base while softening the per-error penalty.
- Knob: `RUSTLE_VG_HMM_PERCOPY_TILT` (env var); only triggers when n_copies ≥ 2 and mean pairwise identity ≥ 0.50 (else plain singleton).

---

## Slide 5 — `05_forward_dp_and_gap.png` *(the E-step + gap rule)*

**One-line:** score the read against each paralog path with the per-copy profile chain; redistribute weight only when the top score's lead is decisive.

**Key points:**
- Each paralog has its own path through the family graph (nodes it traverses).
- Forward DP along that path = log P(read | paralog), using the path's per-copy profiles.
- EM softmax gives a posterior over paralogs.
- **Score-gap rule (advisor's idea):** if (best − runner-up) < Δ, abstain — keep the BAM aligner's 1/NH instead. Default Δ = 10 nats (≈ 4 unambiguous SNPs).
- Mechanically, this means HMM-EM cannot regress on cases that StringTie already gets — the only way to lose a paralog is to mis-redistribute weight, and the gap rule prevents low-confidence redistribution.

---

## Slide 6 — `06_flnc_partial_coverage.png` *(advisor clarification — FLNC ≠ post-cluster2)*

**One-line:** FLNC reads have polyA + cap intact, but they are NOT post-cluster2 isoform consensus — many are partial, and multi-mappers often only span a conserved domain. Per-copy HMM-EM still works because scoring is local.

**Key points:**
- Row 1: post-cluster2 isoform = full-length consensus (downstream, ideal).
- Row 2: FLNC, full-length — possible, often the easy case.
- Row 3: FLNC, partial — common (internal degradation, truncation). Spans only a sub-region.
- Row 4: FLNC multi-mapper at a CONSERVED DOMAIN — the read sees only the part where paralogs look alike. This is the worst case.
- Why per-copy still works in row 4: forward DP integrates **only over the read's path slice**; per-copy SNPs anywhere within that slice contribute. No need to see the easy 5'/3' UTR differences.
- And if even the conserved-domain SNPs aren't enough: gap rule abstains. Worst case, we keep 1/NH — never worse than the BAM aligner.

---

## Slide 7 — `07_amy_result.png` *(empirical answer to advisor #2)*

**One-line:** on the medium-similarity demo, HMM-EM uniquely recovers AMY/LOC101133335 at jaccard 0.52 — squarely in the regime the advisor was worried about.

**Key points:**
- Per-paralog gffcompare exact-match (`=`) recovery on AMY + NBPF subset BAMs.
- LOC101133335 (jaccard 0.52, medium-similarity band) is missed by StringTie AND by heuristic junction-EM. Only the per-copy-profile HMM-EM gets it.
- The other three paralogs (AMY2B, LOC101133271, LOC115933275) are recovered by the gap rule + HMM-EM combination too — no regressions vs StringTie.
- Source: `bench/medium_similarity_demo/DEMO.md`, commit `0ac3475`.

---

## Closing pitch

> The advisor's worry was that VG-driven assembly would collapse paralog
> information. We've shown three independent things that make this not
> the case:
>
> 1. **Architecture** — the family graph keeps per-copy profiles inside each node (slides 2, 3).
> 2. **Mechanism** — POA-MSA tilt α=0.7 + path-local forward DP preserves per-copy SNP signal even at moderate divergence (slides 4, 5, 6).
> 3. **Empirical** — recovers a jaccard 0.52 paralog (LOC101133335) that no other method gets (slide 7).
>
> And the FLNC subtlety: per-copy scoring is local along the read's
> path slice, so partial reads / domain-only multi-mappers are still
> handled — and the gap rule guarantees we never lose to baseline.
