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

One self-contained HTML page (no external assets; theme-aware; favicon 🧬). Structure:

- **Header claim, bold, above everything:** "These are five stages of one pipeline, not five
  approaches." + the one-line role of each stage.
- **The spine:** the five stages top-to-bottom. Each stage box carries exactly:
  1. one mechanism sentence (what it does),
  2. its one **primary** decision number (the decision-relevant heuristic for that stage),
  3. the `file:line` where that number lives.
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

**Success criteria:**
- `gen_heuristics.py` passes (every registry value verified against source) and emits a table
  covering all ~200 heuristics.
- Every `inert-guard` row has a measured observed-max justifying its tier.
- `rustle_mechanism.html` renders the five-stage spine with primary number + file:line per
  stage, the "stages not approaches" claim, and the superseded-terminology table.
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
