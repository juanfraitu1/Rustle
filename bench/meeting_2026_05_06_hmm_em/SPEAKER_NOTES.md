# Per-copy HMM-EM for paralog read assignment — meeting notes (2026-05-06)

Nine slide-ready PNGs in this directory. Suggested order + a few sentences each.
The narrative arc directly answers four advisor concerns:

1. *"The VG just averages information across copies — paralogs collapse."* → Slides 2, 3, 4 show explicitly that we keep one profile per copy, never a single consensus.
2. *"Does it still use a real HMM trellis or some shortcut?"* → Slide 5 shows the actual M/I/D forward DP, with the recurrence written out, and explains the boundary-threading graph extension that adds **no new equations**.
3. *"Where do the priors come from?"* → Slide 6 walks through every step of the EM math (initialization, M-step, E-step, score-gap rule, convergence) with the exact formulas and code line numbers.
4. *"This will only work for near-identical paralogs."* → Slide 9 shows recovery at jaccard 0.52 (squarely in the medium-similarity band).

Plus slide 8 addresses the FLNC misunderstanding directly.

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

**Reading the chart:** each "bar group" at a profile column is an emission distribution `P(base | column)` — one bar per nucleotide A/C/G/T (color-coded). Bar heights at a column sum to 1. Left panel has ONE bar group per column (the consensus profile, paralogs collapsed). Right panel has THREE bar groups per column (one per paralog — A, B, C — labelled above its mini-distribution).

**Key points:**
- Toy MSA: 3 paralogs A/B/C, 4 columns; column 2 (yellow band) is a paralog-distinguishing SNP — A has 'C', B and C both have 'G'.
- **Left, BAD:** consensus emission at column 2 is `P(C)=0.33, P(G)=0.67` — same distribution for every paralog. A read with G at col 2 scores `0.67` against A, B, AND C → no discrimination.
- **Right, GOOD:** per-copy singleton emission gives `P(G | A) = 0.02`, `P(G | B) = 0.98`, `P(G | C) = 0.98`. The SNP signal is preserved at the per-paralog level, EM has something to work with.
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

## Slide 5 — `05_hmm_trellis.png` *(directly addresses advisor #2 — yes, it's a real HMM)*

**One-line:** every per-copy profile is a textbook profile HMM with a full M/I/D forward trellis. The graph extension adds **no new equations** — it just chains per-node trellises by passing the exit-boundary distribution forward.

**Key points:**
- **Left panel** — the inner profile HMM trellis. Three states per profile column j: M_j (match), I_j (insert), D_j (delete). The forward recurrence is the standard one:

  ```
  log M_j(i) = e_M(read[i] | j)
              + logsumexp( M_{j-1}(i-1) + log a_MM,
                           I_{j-1}(i-1) + log a_IM,
                           D_{j-1}(i-1) + log a_DM )
  ```

  with analogous recurrences for I_j(i) and D_j(i). Banded so cost is O(M × W) per node, not O(M × L) — but the math is identical to Krogh '94 / HMMER.
- The **emission** `e_M(b | j)` is exactly the per-copy distribution from slide 4 (α-blend of singleton and POA-MSA column).
- **Right panel** — boundary threading. For a read scored against paralog c's path through K nodes:
  - Entry boundary at node 1 = `[0, −∞, …, −∞]` (no read consumed).
  - Each node fills its trellis using its entry boundary; the exit boundary is read off the last column.
  - That exit boundary becomes the entry boundary of the next node — the chain rule for log probabilities along a path.
  - After node K: `log P(read | path) = α^(K)_exit[L]`.
- Code: `profile_forward_with_boundary_banded` (per-node trellis) + `forward_against_path_for_copy` (chain along path c). Both in `src/rustle/vg_hmm/scorer.rs`.

---

## Slide 6 — `06_em_math.png` *(directly addresses advisor #3 — where do the priors come from)*

**One-line:** the priors come from textbook EM, exactly: each iteration's `π_c` is the share of total weight currently at copy c, with a small floor to keep logs finite.

**Key points (six numbered blocks):**

1. **Initialization (t = 0).** Each multi-mapped read r with placements P_r ⊂ {1…C} starts with `w_{r,c}^(0) = 1 / |P_r|` for every c ∈ P_r. That's the BAM aligner's 1/NH — pure agnostic prior.

2. **M-step (re-estimate copy priors).** For each copy:
   ```
   n_c     = Σ_r w_{r,c}^(t-1)
   π_c^(t) = n_c / Σ_{c'} n_{c'} + ε        with ε = 1e-3
   log π_c^(t) = ln( π_c^(t) )
   ```
   The floor ε guarantees `log π_c` is finite even when no read currently lands on copy c. Code: `vg.rs:2589–2598`.

3. **Per-paralog likelihood (constant across iterations).** For each (read r, placement c ∈ P_r), `log P(r | c) = forward_against_path_for_copy(graph, r, path_c, c)` — the chain of profile HMM trellises from slide 5. Computed once up-front in PHASE 2, parallelized via rayon.

4. **E-step (posterior).** For each (r, c):
   ```
   log Q_{r,c} = log P(r | c) + log π_c^(t)         ← un-normalized log-posterior
   w_{r,c}^(t) = exp(log Q_{r,c} − max_c log Q_{r,c}) / Σ_{c'} exp(·)
   ```
   The softmax is over c ∈ P_r only (paralog-restricted EM — we don't redistribute to copies the read never aligned to). Code: `vg.rs:2603–2640`.

5. **Score-gap rule (don't gamble).** Before applying the softmax:
   ```
   best   = max_c  log P(r | c)
   second = 2nd-max
   if (best − second) < Δ:    skip update — keep w_{r,c}^(t-1)
   ```
   Default Δ = 10 nats (≈ 4 unambiguous SNPs). When the HMM is uncertain we keep the BAM aligner's call. Code: `vg.rs:2614–2624`.

6. **Convergence.** Iterate while `max_{r,c} |w_{r,c}^(t) − w_{r,c}^(t-1)| ≥ 0.001` (and t < max_iter, default 5). Most families converge in 2–3 iterations on the AMY/NBPF demo.

The advisor's concern was *"is this just a heuristic with EM in the name?"* — no. Every step is the standard mixture-model EM, the only non-standard pieces are (a) using HMM forward as the data likelihood, (b) the score-gap abstention. Both are spelled out.

---

## Slide 7 — `07_forward_dp_and_gap.png` *(intuition, picture form)*

**One-line:** picture of slide 6 — score the read against three paralog paths, see the log-likelihoods, see the gap.

**Key points:**
- Three paralog paths through the family graph; the read is scored against each (using slides 5/6's forward DP).
- Numerical example: log P(R|A) = −118.3, log P(R|B) = −149.7, log P(R|C) = −120.9.
- After softmax with current priors: posterior ≈ {A: 0.93, B: 0.00, C: 0.07}.
- The score-gap (best − runner-up) here is 2.6 nats — too small. With Δ=10 default, we abstain.
- Mechanically prevents EM from regressing on cases StringTie already gets — the only way to lose a paralog is mis-redistribution, and the gap rule blocks low-confidence flips.

---

## Slide 8 — `08_flnc_partial_coverage.png` *(advisor clarification — FLNC ≠ post-cluster2)*

**One-line:** FLNC reads have polyA + cap intact, but they are NOT post-cluster2 isoform consensus — many are partial, and multi-mappers often only span a conserved domain. Per-copy HMM-EM still works because scoring is local.

**Key points:**
- Row 1: post-cluster2 isoform = full-length consensus (downstream, ideal).
- Row 2: FLNC, full-length — possible, often the easy case.
- Row 3: FLNC, partial — common (internal degradation, truncation). Spans only a sub-region.
- Row 4: FLNC multi-mapper at a CONSERVED DOMAIN — the read sees only the part where paralogs look alike. This is the worst case.
- Why per-copy still works in row 4: forward DP integrates **only over the read's path slice**; per-copy SNPs anywhere within that slice contribute. No need to see the easy 5'/3' UTR differences.
- And if even the conserved-domain SNPs aren't enough: gap rule abstains. Worst case, we keep 1/NH — never worse than the BAM aligner.

---

## Slide 9 — `09_amy_result.png` *(empirical answer to advisor #4)*

**One-line:** on the medium-similarity demo, HMM-EM uniquely recovers AMY/LOC101133335 at jaccard 0.52 — squarely in the regime the advisor was worried about.

**What is "the band"?** The jaccard-similarity scale at the top of the slide segments paralog pairs into four named regions, each with a different mechanistic story:

| Band | Jaccard | What it means | Difficulty |
|---|---|---|---|
| low-similarity | 0.00 – 0.30 | different gene families (incidentally co-multi-mapping) | should be filtered, not assigned |
| **medium-similarity** | **0.30 – 0.60** | paralogs with substantial divergence — many SNPs, indels, exon turnover | **per-copy profiles needed** |
| high-similarity | 0.60 – 0.90 | paralogs with mostly conserved sequence (SNP density much lower) | per-copy or heuristic both work |
| near-identical | 0.90 – 1.00 | recent duplicates — SNPs sparse, mostly identical | the VG primitive shines here |

The "k-mer Jaccard" axis is `|kmers(P) ∩ kmers(Q)| / |kmers(P) ∪ kmers(Q)|` between a paralog and its nearest sibling — a single number summarizing how much sequence is literally shared at the k-mer level.

**Key points:**
- All four paralogs in this demo land in the **medium-similarity** band (red ticks on the scale at the top). This is the band the advisor was specifically worried about.
- LOC101133335 (jaccard 0.52, dead-center medium-sim) is missed by StringTie AND by heuristic junction-EM. Only the per-copy-profile HMM-EM gets it.
- The other three paralogs (AMY2B 0.59, LOC101133271 0.50, LOC115933275 0.44) are recovered by the gap rule + HMM-EM combination too — no regressions vs StringTie.
- Source: `bench/medium_similarity_demo/DEMO.md`, commit `0ac3475`.

---

## Closing pitch

> The advisor's worries map cleanly onto specific slides:
>
> | Worry | Answered by |
> |---|---|
> | "VG collapses paralog info" | slides 2, 3, 4 (per-copy profiles) |
> | "is it really an HMM?" | slide 5 (M/I/D forward trellis) |
> | "where do the priors come from?" | slide 6 (mixture-model EM, code-tied) |
> | "FLNC = full transcript?" | slide 8 (no — partial reads work because scoring is local) |
> | "only works for near-identical?" | slide 9 (jaccard 0.52 recovery) |
>
> Every step of the math corresponds to a specific function and line range in
> `src/rustle/vg_hmm/scorer.rs` and `src/rustle/vg.rs`. There is no
> hand-wavy step.
