# Mechanism Transparency — Design (Deliverable A)

**Date:** 2026-07-20
**Status:** approved design, pre-implementation
**Motivation:** Advisor criticism — the method appears to "jump from approach to approach
(graph, then k-mers, then POA)" and to hide heuristics. Diagnosis (from a full repo sweep):
the *code* is one pipeline with sequential stages; the *artifacts* describe it in five
mutually inconsistent vocabularies, one artifact makes a claim that is false against the
code, and ~200 heuristics are disclosed unevenly across documents.

This deliverable (A) makes the single method legible and every heuristic explicit.
Code consolidation (collapsing duplicate implementations) is deliverable B, a later spec.

---

## 0. The one idea that does the heavy lifting (organizing principle)

The advisor wants what Eichler has — a method explainable end-to-end where one lever does the
heavy lifting (Eichler: "the read's copy is the one with AS > 10"). Rustle already has a better
one, and the artifact must **lead with it** so the method reads as simple *and* complete, not as
200 knobs.

**Two mechanisms, stated as two sentences:**

- **Define.** Two loci are one family iff they are homologous **by exon-sum** (asm20 identity
  ≥ 0.80 on the exon-union of each copy) and there are ≥ 2 of them.
- **Assign.** A read belongs to a copy iff it matches a base that **no other copy has, beyond
  chance** (Poisson-binomial upper tail, p < α/(n−1), α = 1e-3). If no such base exists, we
  **abstain** — we never guess, and never split 1/k.

This is the `AS > 10` equivalent, and stronger: it is a calibrated, Bonferroni-corrected
p-value rather than a magic number, and where `AS > 10` forces a wrong call on a tie, the
certificate **abstains** — the exact failure of Eichler's rule already documented in this repo.

**Two discoveries are CONSEQUENCES of these two mechanisms, not new methods:**

- **Copies missing from the *genome*** (collapsed / reference-absent, "O4") are a consequence
  of **Assign**. The certificate asks "do the copies I know about explain these reads?" When a
  co-segregating base pattern is carried by reads that **no** modeled copy explains at
  significance, an unmodeled copy must exist. Verified reuse: the same Poisson-binomial +
  Bonferroni α/(n−1) certificate appears in `copy_assign.rs` (assign), `absent_copy.rs`
  (absent-copy `min_p_distinct`), and `collapse_gate.rs` (collapse) — **one test, three uses.**
- **Genes missing from the *annotation*** (previously unknown, verified real because they are
  transcribed, "O1 de-novo / novel O4") are a consequence of **Define**. A locus that (a) is
  unannotated, (b) forms an exon-sum homology edge to a known family, and (c) has transcripts
  (reads assembled through it, ≥ GATE_MIN_READS) is a real family member regardless of
  annotation. Similarity says "it belongs"; exon-sum is *how* similarity is measured;
  transcripts prove it is a real gene.

**The honesty line (must be visible in the artifact):** detection is the consequence; some
downstream *characterization* is separate machinery (e.g. `hidden_copy.rs`'s six thresholds,
the linearize permutation test). The certificate *flags* the absent copy; that machinery
*characterizes* it. Machinery is disclosed in the table (§2.2–2.3), never smuggled into the
two-sentence story.

**What this does to the ~200 heuristics:** they are not decisions. They are the machinery that
computes "homologous by exon-sum?" and "distinguishable beyond chance?". Only ~4 numbers
(id≥0.80, ≥2 loci, α=1e-3, and the exon-sum edge itself) decide anything; the Soto inert-guard
proof (§3) is the evidence that the rest do not get a vote. The disclosure table thus reads as
"look how little decides," not "look, 200 knobs."

**Do NOT overclaim unification.** Define (homology) and Assign (significance) are genuinely two
ideas, not one — the artifact presents them as two. The "one test, three uses" claim is about
the *certificate* (Assign) only, and is verified in code. Family definition is not folded into
the certificate.

---

## 1. The finding this fixes

The live pipeline is **one chain of five stages, not five approaches**:

```
reads
  → [k-mer prefilter]        candidate generation      KMER=18            family_detect.rs:37
  → [asm20 identity edge]    THE family edge test      id≥0.80 cov≥0.50   denovo_pipeline.rs:2370-71
  → [γ-quasi-clique]         cluster the edges         γ=0.20             family_definition.rs:173
  → [de-tie conflict graph]  count copies (χ_H)        de tie, MIN_READS=3 read_conflict.rs:63
  → [IsoCon certificate]     assign each read          α/(n−1)            copy_assign.rs:213,453
assignment
```

k-mers **generate candidates**, alignment identity **is the edge**, the graph **clusters
edges**, `de` **counts copies**, the certificate **assigns reads**. Graph / k-mer / POA are
stages, not alternatives.

**Two facts established from the code during design (both must be reflected in artifacts):**

1. The live family edge is `refine_families_exon_sum` (asm20 id≥0.80, cov-of-shorter≥0.50) —
   called from `gw_family_catalog.rs:221` and `denovo_pipeline.rs:1516`.
2. `detect_families`/`detect_edges` (the `core_recip≥0.13` POA-core edge) has **no caller in
   any shipped binary** — that path is dead as an *edge test*. `T_CORE=0.13` survives live in
   exactly one narrow role: `rescue_thin_loci_iterative` admitting sub-read-floor thin loci as
   family members (`rescue_pipeline.rs`, `family_rescue.rs:41`). Therefore
   `DEFINITIONS_FORMAL.md:119` ("the definitional predicate … POA contiguous-core ≥ 0.13") is
   **incorrect** and will be corrected — 0.13 is thin-locus rescue admission, not the edge.

---

## 2. Deliverables

### 2.1 `bench/rustle_mechanism.html` — the single stage chain

One self-contained HTML page (no external assets; theme-aware; favicon 🧬). Structure,
**top to bottom, leading with the one idea (§0)**:

- **§0 first — "The method in two sentences."** Define + Assign, in full, as the opening. Then
  "The one test" panel: the certificate written once, with the three call sites where it is
  reused (`copy_assign.rs`, `absent_copy.rs`, `collapse_gate.rs`) — the `AS > 10` equivalent,
  shown to do the heavy lifting.
- **"Two consequences, not new methods."** The two consequence statements from §0 (copies
  missing from the genome ⇐ Assign; genes missing from the annotation ⇐ Define + exon-sum),
  each with its one-line "the honesty line" caveat naming the separate characterization
  machinery.
- **Header claim, bold:** "These are two mechanisms and their consequences — not four
  approaches."
- **The spine (now framed as "how the two sentences are computed"):** the five stages
  top-to-bottom. Each stage box carries exactly:
  1. one mechanism sentence (what it does),
  2. its one **primary** decision number (the decision-relevant heuristic for that stage),
  3. the `file:line` where that number lives.
  The spine is explicitly subordinate to §0: it is the machinery of Define + Assign, not a
  list of independent tricks.
- **Superseded-terminology table:** every vocabulary the advisor has seen mapped to its real
  stage or to "abandoned (date)":
  | term seen | reality |
  |---|---|
  | k-mer Jaccard / Union-Find (`figures/01-04`) | abandoned architecture, ~2026-04 |
  | profile-HMM emission (`advisor_response.md`) | abandoned proposal, 2026-05 |
  | POA contiguous-core 0.13 as "the edge" | thin-locus **rescue** admission only, not the edge |
  | "no ≥X% identity cutoff" (`quasiclique.html`) | false; edge IS id≥0.80 (corrected) |
- **Embedded heuristics appendix:** the full table from §2.2, tiered (§2.3), rendered inline.

### 2.2 `bench/heuristics.toml` + `bench/gen_heuristics.py` → `bench/heuristics.tsv`

**`heuristics.toml`** — hand-maintained registry, one entry per heuristic:

```toml
[[heuristic]]
name = "RefineParams::min_identity"
file = "src/rustle/vg_family/denovo_pipeline.rs"
line = 2370
value = "0.80"
stage = "O1-edge"           # O1-prefilter|O1-edge|O1-cluster|O1-count|O2|O3|O4|pre-O1
tier = "decision"           # decision|result-guard|inert-guard|derived
kind = "arbitrary"          # arbitrary|derived
rationale = "asm20 exon-sum homology edge identity floor; the family edge test"
```

**`gen_heuristics.py`** — reads the registry and, for every entry, **opens `file` at `line`
and asserts the literal `value` still appears there**. On mismatch it fails loudly (non-zero
exit, names the drifted entry). On success it emits `heuristics.tsv` (sorted by stage, tier)
which the HTML embeds. The verify-against-source step is the anti-drift device: the registry
cannot silently fall out of step with the code.

The generator is intentionally a *verifier*, not a scraper — no regex mining. A new constant
is disclosed by adding a registry line; the assertion guarantees the line is truthful.

### 2.3 The C-tiering (blast-radius classification)

Every heuristic gets a `tier`:

- **decision (~15):** numbers that pick scientific results. Named inline in the spine.
  (id≥0.80, cov≥0.50, γ=0.20, GATE_MIN_READS=3, de tie predicate, α=1e-3, junction_weight,
  error_rate, …)
- **result-guard (flagged loudly):** guards that *can* change the result. Known members:
  `POA_CAP=20000` (crosses → poasta swaps to minimap2 → different PSV set),
  `READ_CAP=6000` (truncates a family's read pool by BAM record order),
  the sensitive-tier mismatch (`0.70` literal at `denovo_pipeline.rs:2506` vs
  `sensitive_identity=0.60` at `:2645`), and the five inconsistent per-base error rates
  (0.003 / 0.01 / 0.005 / 0.05 / 0.001 across copy_assign / quant / mosaic / conflict).
- **inert-guard:** cost guards claimed not to move numbers — **admitted to this tier only by
  proof** (§3): the observed max of the guarded quantity, on the Soto RNA families, sits below
  the guard. If it fires, it is reclassified `result-guard`. No guard is called inert on
  assertion.
- **derived (~20):** Bonferroni α/(n−1), exact Poisson-binomial tail, BH-FDR, `tau_from_p`
  Bayes-optimal mapping, permutation p, log-gamma numerics.

### 2.4 Corrections to existing artifacts (in scope for A)

- `bench/quasiclique.html`: remove/replace the false "no arbitrary ≥X% identical cutoff"
  claim; state the edge is asm20 id≥0.80 cov≥0.50, and that γ=0.20 is the *clustering* density
  on top of those edges (a different quantity at a different stage).
- `bench/DEFINITIONS_FORMAL.md:119`: correct the definitional predicate to the asm20 edge;
  relabel 0.13 as thin-locus rescue admission.
- `figures/01_pipeline_overview.md` (+02/03/04) and
  `analysis/family_graphs/docs/advisor_response.md`: add a superseded banner at the top of each
  (date + pointer to `rustle_mechanism.html`). Do not delete — they are history.

### Out of scope for A (→ deliverable B)

Collapsing duplicate implementations (two refiners, two "≥2 loci" predicates, three EM paths),
unifying the five error-rate constants, moving buried constants behind CLI flags, and the six
hard-coded `/home/juanfra/winloci_scratch/` default paths. A *discloses* these (they appear in
the table, correctly tiered); B *fixes* them.

---

## 3. Verification (real data, Soto RNA families)

**Subset:** the Soto RNA-recovered families — the population the advisor evaluates the method
on. 66 detected families in `bench/soto/a119b_detected_families.tsv`; eval harness
`bench/soto/soto_detection_eval.py` (reads `/mnt/linuxdisk/home/juanfraitu/winloci_data`, now
mounted). Full BAM `GGO_mm.bam` + `GGO.fasta` in `winloci_scratch`; small chr19 BAM for fast
iteration.

**V1 — spine numbers are real.** Re-run the live reproduction commands and confirm each
primary number in the spine matches the artifact:
- `gw_family_catalog` on the Soto regions → family edges use id≥0.80/cov≥0.50 (confirm the
  values the binary actually applies).
- `copy_assign --region <GSTM> --homology-primary --min-copies 2` → n_copies=3, α=1e-3 gate.
- `asj --region <ABCC4>` → the O3 numbers.

**V2 — inert-guard proof.** For each guard proposed as `inert-guard`, instrument or post-hoc
measure the observed max of the guarded quantity across the 66 Soto families and record it in
the table's "why inert" column. Guard is inert iff observed-max < guard. Concretely at least:
- READ_CAP=6000: max reads pooled per Soto family.
- POA_CAP=20000 / LEN_CAP 9000/20000: max POA input length among Soto copies.
- MAX_MEMBERS=30 / MAX_LOCI=60: max family size in the Soto set (does any Soto family bypass
  the repeat-bridge / recombinant-split gate?).
Any guard that fires on a Soto family is reclassified `result-guard` and flagged in the HTML.

**V3 — the two known result-changers, toggled honestly.** For READ_CAP and POA_CAP, run the
affected Soto loci with the guard at its default vs relaxed, and report the delta in family
count / copy count / assignment md5. Report the truth, whatever it is.

**V4 — Tandem removal-recovery positive control (the headline evidence for §0's "consequence
of Assign" claim, and a figure in the artifact).** Ground-truth demonstration that a copy
deleted from the reference is put back by the *same* certificate — no copy-finder added. This
is the inverse of the existing `linearize` certificate (linearize folds a copy IN and asks "do
reads go unique?"; this pulls a copy OUT and asks "do the leftover reads point back to it?").

*Evaluation is permutation-invariant set matching:* a family is a **set** of copies (paths),
so each recovered copy is matched to its nearest ground-truth copy; recovering
{copy1, copy3, copy2} for {copy1, copy2, copy3} is a perfect recovery. Report copy-**number**
recovery and per-copy **sequence** identity to the deleted copy. "Linear combination of the
overhangs" is realized as: the column-wise pileup on the surviving copies is a superposition of
each true copy's contribution; the recovered copy is the residual addend the survivors do not
explain at significance (deconvolution), reconstructed as the consensus of the reads carrying
that residual.

- **V4a — Simulation (airtight feasibility proof; same individual by construction).** Plant a
  3-copy tandem array from a real seed, with **controlled divergence** (PSVs at known
  positions) and **shared splice structure** (repo pitfall: identical copies fabricate
  junctions → plant shared donors/acceptors + reverse-complement where relevant). Simulate
  HiFi/IsoSeq reads from all three; align with the project preset (`-N 50`, splice for RNA;
  **not** `--secondary=no`, which yields 0 families). Delete one copy, then two; recover;
  set-match against the planted copies. **Second sub-arm — the honest floor:** three
  *identical* copies (no PSVs) ⇒ recover copy-**number** correctly but **abstain** on separable
  sequence (the K=0 identifiability floor). Both sub-arms shipped; the floor is shown, not
  hidden — it is the strongest evidence the method claims only what the data support.
- **V4b — Real (gorilla A119b IsoSeq; the biological credibility the advisor requires).**
  Substrate: a real Soto tandem family (co-located, ≥3 copies) on gorilla. Re-align A119b to
  `mGorGor1` (indexed on disk). **Same-read-set ground truth (removes the cross-individual
  confound, since A119b ≠ Kamilah):** align A119b to the *full* mGorGor1, select loci where
  A119b's own reads support 3 divergent copies at significance; that set is the ground truth.
  Delete a copy from the assembly, re-align the same reads, recover, compare to the intact copy.
  Per the crash rule: per-locus, **foreground, serial, small batches**, outputs to
  `winloci_scratch`. Downsampled `A119b_ds.bam` + region extraction keeps it light.

**Success criteria:**
- `gen_heuristics.py` passes (every registry value verified against source) and emits a table
  covering all ~200 heuristics.
- Every `inert-guard` row has a measured observed-max justifying its tier.
- V4 runs on both arms: simulation recovers deleted divergent copies (number + sequence, set-
  matched) and correctly abstains on identical copies; the real gorilla arm recovers a deleted
  copy of a Soto tandem family against same-read-set ground truth. The recovery figure is
  embedded in the artifact's "consequence of Assign" section.
- `rustle_mechanism.html` **leads with §0** (the two-sentence method + the one reused
  certificate + the two consequences), then renders the five-stage spine (primary number +
  file:line per stage, subordinate to §0), the "two mechanisms and their consequences, not four
  approaches" claim, and the superseded-terminology table.
- The "one test, three uses" reuse claim is verified in code before it appears in the artifact
  (the Poisson-binomial + Bonferroni α/(n−1) form is confirmed present at all three call sites).
- The three corrections (§2.4) are applied.
- No artifact remains that contradicts the code (the quasiclique false claim and the
  DEFINITIONS_FORMAL 0.13 claim are gone).

---

## 4. Components & boundaries

| unit | purpose | input | output | depends on |
|---|---|---|---|---|
| `heuristics.toml` | source of truth for disclosure | (hand-maintained) | — | file:line in src/ |
| `gen_heuristics.py` | verify registry vs source, emit table | toml + src/ | `heuristics.tsv`, exit code | Python stdlib only |
| `rustle_mechanism.html` | the legible single method | `heuristics.tsv` (embedded at build) | the artifact | gen output |
| Soto verification | back the tiering with real data | Soto BAM/regions | observed-max per guard | mounted disk, built binaries |
| V4 tandem recovery | ground-truth proof of "absent copy = consequence of Assign" | simulated array + gorilla A119b/mGorGor1 | recovery figure (number + set-matched sequence id) | mounted disk, minimap2, certificate binaries |

Each is independently checkable: the generator has a pass/fail exit; the HTML is a static
render of a known TSV; the verification is a set of measured numbers.

---

## 5. Risks

- **Registry completeness.** ~200 entries is a large hand-write. Mitigation: seed the toml
  from the design-time sweep (already enumerated), then `gen_heuristics.py` guarantees each
  seeded line is truthful. Missing a constant is a disclosure gap, not a falsehood — acceptable
  to iterate.
- **A guard thought inert turns out to fire.** That is a *finding*, not a failure — it moves to
  result-guard and is flagged. The design explicitly permits this.
- **Line numbers drift as B lands.** Expected; `gen_heuristics.py` fails loudly and the
  registry line is updated. This is the mechanism working, not breaking.
