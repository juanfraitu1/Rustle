# F4 — scope: the Canzar-shaped theory capstone (LP-rounding approximation for copy assignment)

**Goal.** Complete the identifiability theory's missing piece: an **approximation algorithm with a provable
ratio** for the NP-hard general case, cast in the advisor's own paradigm (Canzar 2016 = multimapping resolution
as **facility location**, solved by **LP relaxation + rounding** with approximation guarantees). Today we have
hardness (Thm 1), exact recovery under Strong Separation (Thm 2–3), and a per-read certificate (Thm 4) — but
**no guarantee for the hard general instance**, which is exactly the shape his 2016 paper fills.

This is a scope/plan, not a proof. It fixes the formulation, the target theorems, the algorithm, the machine-
check, and the honest boundaries, so the proof work is well-defined.

---

## 1. The formulation (faithful to Canzar's facility-location frame)

**MAX-WEIGHT COPY ASSIGNMENT (MWCA).**
- **Clients** = reads `R` (partial allele-functions over PSV columns).
- **Facilities** = candidate copies `C` (allele-vectors / haplotypes; from the family layer, or the maximal
  independent sets of the de-tie conflict graph `H̄`).
- **Compatibility** `N(r) ⊆ C` = the copies read `r` is consistent with (no contradicted PSV).
- **Weight** `w(r,c) = ` the read's log-likelihood for copy `c` (the production gate's per-column LLR), defined
  only for `c ∈ N(r)`. This is the principled replacement for `1/k` the advisor dislikes — a real evidence
  weight, not a uniform split.
- **Decision:** open a set `S ⊆ C` of copies and assign each read to one open compatible copy
  `σ(r) ∈ S ∩ N(r)` (or leave it unassigned), to **maximize `Σ_r w(r, σ(r))`** subject to a budget `|S| ≤ K`
  (the cover size, `= χ(H)` by Lemma 1) — or, in the facility-cost variant, minus `Σ_{c∈S} f_c`.

This is **uncapacitated facility location / max-coverage with a cardinality budget** — Canzar's exact paradigm.
Two objective variants, with very different approximability (the key honest fork, see §2).

**Bridge to the existing object.** `S` is a conflict-free cover (each copy = an independent set in the read-
conflict graph `H`; Lemma 1: a feasible `S` of size `K` ⟺ a proper `K`-coloring of `H`). So MWCA is the
**weighted, soft-assignment relaxation of the minimum cover** — the cover decides *which* copies, MWCA decides
the best evidence-weighted assignment given that choice.

---

## 2. The objective fork (decide first — it determines what's provable)

| Objective | What it is | Approximability | Verdict |
|---|---|---|---|
| **MIN-COVER** (minimize `K = χ(H)`) | fewest copies to explain all reads | graph coloring → **inapproximable** within `n^{1−ε}` for general `H` (unless P=NP) | **Do NOT target in general.** Only tractable on *structured* `H` (perfect/interval/bounded-degree). Keep as an honest hardness boundary + a structured-case result (see §5). |
| **MAX-ASSIGNMENT** (maximize weight, `≤K` copies) | best evidence-weighted assignment under a copy budget | submodular max-coverage → **(1 − 1/e)** via greedy/LP-pipage; with pairwise copy-conflicts → Canzar-style **constant-factor** LP-rounding (his 0.19–0.44 regime) | **TARGET THIS.** Clean constant-factor guarantee, Canzar's exact method, evidence-weighted (no 1/k). |

**Recommendation: target MAX-ASSIGNMENT.** It is the version that *has* a constant-factor guarantee, it is
literally Canzar's facility-location problem, and it uses the gate's LLR weights (his aesthetic). The MIN-COVER /
coloring side becomes the honest *inapproximability boundary* of the chapter (we already have its NP-hardness in
Thm 1), plus a tractable result on the structured conflict graphs our data actually produces.

---

## 3. Target theorems (numbered to extend the note)

- **Thm 5 (NP-hardness of MWCA).** MWCA is NP-hard (reduce from Thm 1 / max-coverage). *(Likely short.)*
- **Thm 6 (LP-rounding approximation — the capstone).** The natural LP relaxation of MWCA
  (`x_{r,c} ∈ [0,1]` assignment, `y_c ∈ [0,1]` opening; `Σ_c x_{r,c} ≤ 1`, `x_{r,c} ≤ y_c`, `Σ_c y_c ≤ K`,
  `x_{r,c}=0` for `c ∉ N(r)`) admits a rounding achieving a **constant factor** of the optimum:
  `(1 − 1/e)` for the cardinality-budget / no-conflict case (pipage / dependent rounding), degrading to a
  Canzar-style constant when pairwise copy-conflicts are added. *(The core proof effort.)*
- **Thm 7 (integrality bridge — ties F4 to Thm 2–4).** Under **Strong Separation**, the MWCA LP is **integral**
  and its optimum equals the true cover: the approximation collapses to **exact recovery**, recovering Thm 2/3
  as the integral special case, and the per-read certificate (Thm 4, `min_p`) is exactly the **LP dual /
  complementary-slackness witness** that certifies a read's assignment is tight. *(This is the elegant glue that
  makes the whole theory one object: hardness → approximation → exact-under-separation → per-read certificate.)*

---

## 4. The VG / flow-decomposition reading (the thesis framing, optional second lens)

Per `project_thesis_framing_family_vg`: copies = paths through one PSV-aware variation graph; reads = evidence on
edges. MWCA then reads as a **min-cost / max-weight flow decomposition**: open `≤K` paths (copies) carrying the
read-evidence flow, maximizing assigned weight. LP-rounding of path-flow = the same theorem in flow language, and
flow LPs are often integral (totally-unimodular structure) — a candidate *route to Thm 7* and to better-than-
constant ratios on the structured (interval/series-parallel) instances RNA produces. Worth exploring as the proof
vehicle for Thm 6/7; it also lands the "FLIP Canzar: multimapping = shared evidence" thesis message.

---

## 5. Structured-instance angle (where MIN-COVER becomes tractable — bonus)

Real read-conflict graphs are not arbitrary: PSV-column conflicts give them interval/comparability-like
structure. If `H` (or `H̄`) is shown to be **perfect** (or interval, or bounded-treewidth) on co-located
families, then `χ(H)` = MIN-COVER is **poly-time exact** there — a much stronger statement than the constant-
factor MAX-ASSIGNMENT bound, for the instances that matter. **Action:** characterize the conflict-graph class
empirically (are the GGO family conflict graphs perfect? chordal? interval?) and prove the structural lemma if so.
This is the highest-payoff side quest and may be easier than the general approximation.

---

## 6. Machine-check plan (mirror the existing exhaustive style)

`bench/mwca_approximation_check.py`, exhaustively over small universes (K≤3, m≤3 PSV columns, |R| small):
- enumerate instances; solve the **LP** (scipy.optimize.linprog or PuLP) and the **integer optimum** (brute force).
- verify the **rounded** solution achieves the claimed ratio vs the integer OPT on *every* instance (the bound is
  not just expected-case);
- verify **Thm 7**: on Strong-Separation instances the LP is integral and `= ` the true cover, and `min_p`
  matches the dual slack. (Reuse the Strong-Sep generator from `copy_assignment_theory_checks.py`.)
This makes Thm 6/7 machine-witnessed the way Thm 1–4 already are.

---

## 7. Honest boundaries (state up front, the way the note already does)

- MAX-ASSIGNMENT is approximable; **MIN-COVER (coloring) is not** in general — the chapter must lead with the
  objective choice, not blur them. (This is the single most important honesty point.)
- The constant factor (≈`1−1/e`, or Canzar's `~0.19` with conflicts) is a *worst-case* guarantee; real instances
  do far better (the structured-case §5 + integrality §3 explain why). Report both.
- `w(r,c)` inherits the gate's assumptions (error model, editing filter); the approximation is over the *given*
  weights, not a claim the weights are perfect.
- This is **theory**, not a production change — like Thm 1–4, `RECOVER`/MWCA need not run in the pipeline; the
  shipped gate is Thm 4. F4 guarantees the algorithm you *could* run, in his paradigm.

---

## 8. Work plan (phases, each independently checkable)

1. **Fix the formulation + objective** (this doc; recommend MAX-ASSIGNMENT). — *decision*
2. **Thm 5** NP-hardness of MWCA. — *short proof + note section*
3. **Thm 6** LP + rounding + ratio, MAX-ASSIGNMENT (start no-conflict `1−1/e`; then add conflicts → Canzar
   constant). — *core effort; try the flow-decomposition vehicle §4*
4. **Thm 7** integrality under Strong Separation + `min_p` = LP-dual witness. — *the glue to Thm 2–4*
5. **Structured-instance lemma** (§5): characterize + (if perfect) exact `χ(H)` on co-located families. — *bonus*
6. **`bench/mwca_approximation_check.py`** exhaustive machine-check (§6). — *witnesses 5–4 like the rest*
7. **Write §5c "MWCA / facility-location" into `bench/copy_assignment_theory.md`** + update the checks suite.

**Smallest valuable slice (if time-boxed):** steps 2 + 3-no-conflict + 6 give a clean, machine-checked
`(1−1/e)` LP-rounding theorem in Canzar's paradigm — already the missing capstone. Steps 4/5 are the elegance
multipliers.
