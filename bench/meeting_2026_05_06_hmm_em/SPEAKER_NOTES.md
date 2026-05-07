# Per-copy HMM-EM for paralog read assignment — meeting notes (2026-05-06)

Eleven slide-ready PNGs in this directory. Suggested order + a few sentences each.
The narrative arc directly answers four advisor concerns:

1. *"The VG just averages information across copies — paralogs collapse."* → Slides 2, 3, 4 show explicitly that we keep one profile per copy, never a single consensus.
2. *"Does it still use a real HMM trellis or some shortcut?"* → Slide 5 shows the actual M/I/D forward DP, with the recurrence written out, and explains the boundary-threading graph extension that adds **no new equations**.
3. *"Where do the priors come from?"* → Slide 6 walks through every step of the EM math (initialization, M-step, E-step, score-gap rule, convergence) with the exact formulas and code line numbers.
4. *"This will only work for near-identical paralogs."* → Slide 9 shows recovery at jaccard 0.52 (squarely in the medium-similarity band).

Plus slide 8 addresses the FLNC misunderstanding directly, slide 10 sketches an additive extension to the low-similarity band (jaccard < 0.30) that reuses the existing HMM trellis without overriding any of the medium/high-similarity machinery, and slide 11 is a didactic walk-through showing how each variant type (SNP, indel, alt-donor, alt-acceptor, UTR, exon-skipping) maps to a specific feature of the M/I/D trellis or graph topology.

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

**"What is a nat?" inset (bottom of slide).** Every log-space quantity on slides 5–7 is in nats (natural-log, not log₂) because the forward DP's logsumexp is cleanest in that base. Useful conversions to have ready when the advisor asks about Δ = 10:
- 1 nat ≈ 1.44 bits.
- `e^10 ≈ 22,000` → a 10-nat gap means one hypothesis is ~22,000× more likely than the runner-up.
- One unambiguous SNP at a profile-match column is worth ≈ 2.5 nats of evidence, so Δ = 10 ≈ 4 unambiguous SNPs of separation. That's the intuition for why the default isn't aggressive — it's exactly the sort of evidence threshold a careful biologist would want before flipping a copy assignment.

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

## Slide 10 — `10_low_similarity_extension.png` *(forward look — extending into the low-similarity band)*

**One-line:** for paralogs that are too diverged for the existing k-mer linker to even put them in the same family (jaccard < 0.30), an additive two-pronged extension reuses the same HMM trellis without overriding any of slides 5–7. Nothing existing changes.

**Why this slide exists:** the natural follow-up to slide 9 — "you've shown 0.52, what about 0.20?". The answer isn't "throw away the HMM" — it's "use the HMM differently, and only when the band warrants it."

**Two prongs (left/right boxes on the slide):**

1. **Family-discovery side.** The existing k-mer Jaccard linker (default ≥ 0.30) won't link two paralogs that share < 30% of their k-mers, so they never even meet inside a family graph. Add a complementary signal that survives sequence drift:
   - **splice-site Jaccard** — donor/acceptor positions stay conserved long after sequence does.
   - **conserved-domain anchors** — short Pfam / InterPro motifs identifiable directly from genomic k-mers.

   Two bundles enter the same family if EITHER signal fires. Gated by `--vg-low-sim-link` (off by default) so the production AMY/NBPF path is untouched.

2. **Scoring side.** The existing per-copy `forward_against_path_for_copy` is rigid — it forces the read to fit ONE paralog's exact path. At low similarity that's exactly what hurts: a normal-error read collects penalty against any one path. Swap it for `assign_via_graph_viterbi`:
   - Run `viterbi_path(fg, read)` — same M/I/D trellis from slide 5, no per-paralog path constraint. The read picks its own best path through the graph.
   - For each known paralog c, score by **node-overlap** with the Viterbi trace:
     `recall(c) = |viterbi ∩ path_c| / |path_c|`
     `precision(c) = |viterbi ∩ path_c| / |viterbi|`
   - Sort by recall descending. Top paralog = the assignment.

**Dispatch (bottom of slide):** auto-solver routes per family by similarity band. High/medium → existing path-forward. Low → graph-Viterbi. Near-identical → existing `--vg-snp` SNP-aware E-step. The M-step (priors) and gap rule are unchanged — only the per-(read, copy) score in the E-step differs.

**Generalized gap rule:** for graph-Viterbi the score-gap analogue uses recall (which is bounded in [0, 1] across families, unlike nats which scale with read length). Proposed default: abstain if `top.recall − second.recall < δ_recall` with `δ_recall = 0.20`.

**POC status (talking point — the slide footers say "POC, 3 tests passing"):**
- New function `assign_via_graph_viterbi` lives in `src/rustle/vg_hmm/scorer.rs`.
- Three unit tests on a synthetic 4-node, 2-paralog graph (paralog-specific middle nodes, deliberately low DNA-Jaccard between A and B) — all passing.
- Test 3 specifically shows: on a noisy read (~15% per-base error), forward log P gap = 61.42 nats but Viterbi recall gap = 0.333. Both methods agree on the assignment, but recall is bounded and portable across families; nats are not.
- **Honest caveat to lead with:** this is a mechanism POC, not yet a benchmark. The end-to-end test needs (1) above — the family-discovery side has to land first, otherwise the two paralogs never enter the same FamilyGraph in the first place. Source: `bench/low_similarity_poc/POC.md`.

---

## Slide 11 — `11_hmm_positional_info.png` *(didactic — variant-type to HMM-feature map)*

**One-line:** the HMM's M/I/D trellis is a position-aware ruler. State index = column index in the exon-class profile, regardless of read length. Every variant type maps to a specific feature of that ruler — none of them lose the position anchor.

**Why this slide exists:** the advisor has asked "but what about indels? alt splicing? UTR?" — and rather than handwave, this slide shows the literal HMM-state translation for each common variant. It's the most fundamental "how does the math handle X" slide in the deck.

**Reading the slide:**

- **Top reference panel** — one paralog's exon-class profile drawn as a chain of 10 match-states M_1…M_10. Each circle owns a `P(base | column)` distribution (slide 4's per-copy emission). The state index is the position anchor. Everything below the reference panel shows what each variant does to *this* picture.

- **Six variant panels (3 cols × 2 rows):**

  | # | Variant | What changes in the trellis / graph |
  |---|---|---|
  | 1 | **SNP at one column** | Same M_5 state for both paralogs; per-copy emission distribution differs: `e_M(b\|5)` for A favors A, for B favors G. Position preserved; SNP signal lives in emission. |
  | 2 | **Insertion / deletion** | Insertion → one or more I-states between M-cols (the `gat` example shows 3× I_4 emitting between M_4 and M_5). Deletion → silent D-states (D_5, D_6, D_7 advance the column without consuming a base). M-column anchors stay fixed. |
  | 3 | **Alt 5' splice donor** | Each paralog's exon-class node has its OWN length (per-copy profile from slide 3). Boundary threading exits at col 6 vs col 9. Different junction edges leave each copy. |
  | 4 | **Alt 3' splice acceptor** | Mirror of donor — per-copy profile starts at a different effective entry column. Equivalent: graph carries paralog-specific entry edges into the same exon-class node. |
  | 5 | **UTR variation** | UTR is its own exon-class node (not collapsed into the CDS exon). Each paralog gets its own per-copy profile of paralog-specific length. Reads are not required to span the UTR — partial coverage is allowed (slide 8). |
  | 6 | **Cassette exon / exon skipping** | Same M/I/D trellis inside every node. Different paralogs traverse different EDGES through the family graph: `forward_against_path_for_copy` chains the trellises along *each copy's* path. Skipping is a graph-topology change, not a trellis change. |

**Key didactic move:** the slide separates "what's preserved" (the M-state index = position anchor) from "what differs" (which feature absorbs the variant). Five of the six variants leave the trellis structure untouched — only the cassette case changes which nodes the path visits, and even then the per-node trellis is identical.

**Talking-track for the advisor:**
> "Each panel here is the literal answer to 'how does the HMM handle X?'. SNPs sit in the emission. Indels live in I and D states — that's textbook profile-HMM, exactly Krogh '94. Alt splicing changes the per-copy profile length and the graph edges — boundary threading from slide 5 already handles that with no new equations. UTR is just another node. Cassette exons = same trellis, different graph traversal. We didn't invent anything new for any of these; they all fall out of the standard profile-HMM-on-a-graph framework."

**Code anchors:** `forward_against_path_for_copy` and `profile_forward_with_boundary_banded` in `src/rustle/vg_hmm/scorer.rs` (the M/I/D recurrences) — same functions that handle every variant type above. The graph-topology cases are handled at the `path_c` level inside `family_graph.rs`.

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
> | "what about really diverged paralogs (jaccard < 0.30)?" | slide 10 (additive extension — same HMM, graph-Viterbi assignment, POC passing) |
> | "how does it handle indels / alt splicing / UTR?" | slide 11 (each variant maps to a specific M/I/D state or graph edge) |
>
> Every step of the math corresponds to a specific function and line range in
> `src/rustle/vg_hmm/scorer.rs` and `src/rustle/vg.rs`. There is no
> hand-wavy step.
