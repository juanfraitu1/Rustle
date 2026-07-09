# EM copy-assignment = soft SDA PSV-clustering: derivation and consistency theorem

*This note derives the EM copy-assignment engine (`src/rustle/vg_family/em_copy_assign.rs`) as the
maximum-likelihood **soft relaxation of SDA's PSV correlation-clustering** (Vollger et al., Nat Methods 2019 —
the paper the advisor sent as prior art), run on the thesis's PSV-aware variation graph, and proves it
**consistent** in the identifiable regime. The point is that the method is not "an EM we picked": every piece is
forced — the per-read likelihood is SDA's attraction/repulsion made continuous, the E-step is SDA's read
partition made soft, and the identifiable regime is exactly the Strong-Separation ⟹ unique-minimum-cover
(MCC = χ(H)) result already proved in `bench/THEORY.md`. The consistency theorem is the provable layer
**under** SDA's NP-complete heuristic, and it explains SDA's empirical 91–93% accuracy floor.*

**Cross-references.** `bench/THEORY.md` (the consolidated combinatorial theory — Lemma 1, Theorems 1–7, the
§5 Proposition/K-frontier, §5·SUN, §6b, IDENTIFIABILITY_LIMITS; formerly `bench/copy_assignment_theory.md`,
merged verbatim). `reference_sda_vollger` (Vollger et al. 2019). `reference_sudmant_2010_sun`,
`reference_isocon_sahlin`, `reference_clair3_rna` (per-column filters). Design spec:
`docs/superpowers/specs/2026-07-08-em-copy-assignment-design.md`.

---

## §1 Derivation from SDA — the EM is SDA's PSV graph made continuous

### 1.1 SDA's PSV graph (Vollger et al. 2019)

SDA (Segmental Duplication Assembler) resolves **collapsed** segmental duplications from long DNA reads
(`reference_sda_vollger`). After detecting a collapse by read-depth excess, it calls **PSVs** (paralogous
sequence variants — the second-most-frequent base at positions whose total depth is ≈ a multiple of single-copy
coverage, so the coverage gate separates PSVs from heterozygous alleles) and builds a **PSV graph**:

- **nodes** = PSVs;
- **reads** = edges between the PSVs they carry;
- a read carrying two PSVs **on one molecule** is an **attraction** edge — the two PSVs belong to the *same*
  paralog/copy;
- two PSVs that are **mutually exclusive** across reads (never co-carried) give a **repulsion** edge — they
  belong to *different* copies.

SDA then runs **correlation clustering** on this signed graph — *ab initio*, with **no preset number of
clusters**, which is **NP-complete**, solved heuristically with **15 random initializations** — to assign
PSVs → paralogs, and **WhatsHap** partitions the reads by the resulting PSV clusters. SDA validates at **91–93%**
of PSVs correctly assigned (SRGAP2 / NOTCH2NL, BAC ground truth) and states the identifiability limit outright:
*"virtually identical duplications cannot be distinguished and will require even longer reads (>100 kb)."*

### 1.2 The same object on the PSV-aware VG

On the RNA side the thesis already carries the same structure as a **PSV-aware variation graph**, one per family
(`project_psv_aware_vg`, `psv_linkage.rs`): a shared backbone with **parallel paths = copies** and **bubbles =
PSV columns**. Copy `k` is the path carrying a specific allele at each bubble — its allele-vector
`θ_k ∈ ∏_j A_j` over the PSV columns `[m]` (built from read-supported PSV bubbles, i.e. derived from the data,
SDA-style, not a handed-down catalog). A read `r` is a **partial path**: `obs(r) ⊆ [m]` is the set of bubbles it
spans on one molecule, and `r(j) ∈ A_j` the allele it carries at bubble `j`. Two PSVs co-carried by one read is
SDA's attraction edge; the read's `obs(r)` vector is precisely the "PSVs on the same molecule" that SDA reads off
its edges. This is the identical column/allele model of `bench/THEORY.md` §2 (reads as partial allele functions,
copies as consistent allele-vectors, conflict = disagreement at a co-observed column).

### 1.3 The per-read likelihood IS SDA's attraction/repulsion, made continuous

SDA's edges are **hard** (an attraction or a repulsion). Replace the hard signed edge by a per-column
substitution likelihood. Let `e_j` be the per-column sequencing-error probability at bubble `j` over an alphabet
`A_j` (for DNA `|A_j| = 4`, so a specific wrong base has probability `e_j/3`). Define the emission at one bubble

$$
q_j\bigl(a \mid b\bigr) \;=\;
\begin{cases}
1 - e_j & a = b \quad(\text{observed allele matches the copy's allele}),\\[2pt]
e_j/(|A_j|-1)\;(\approx e_j/3) & a \neq b \quad(\text{mismatch}).
\end{cases}
$$

The **per-read per-copy likelihood** is the product over the bubbles the read spans (reads independent given
origin; columns conditionally independent given the copy-path — the standard error-model factorization):

$$
L_{rk} \;=\; P\bigl(\mathrm{obs}(r) \,\big|\, z_r = k,\ \theta_k\bigr)
\;=\; \prod_{j \in \mathrm{obs}(r)} q_j\bigl(r(j) \,\big|\, (\theta_k)_j\bigr),
$$

so that

$$
\log L_{rk} \;=\; \sum_{j\in\mathrm{obs}(r)}
\Bigl[\, \mathbb{1}\!\left[r(j)=(\theta_k)_j\right]\log(1-e_j)
\;+\; \mathbb{1}\!\left[r(j)\neq(\theta_k)_j\right]\log\!\tfrac{e_j}{|A_j|-1}\Bigr].
$$

Read the two indicators as SDA's two edge types. At a bubble where read `r` carries **copy `k`'s private allele**
(a SUN column for `k`, `reference_sudmant_2010_sun`), `r(j) = (θ_k)_j` adds `log(1−e_j)` to copy `k` — a **soft
attraction to `k`** — while `r(j) ≠ (θ_{k'})_j` for every other copy `k'` adds `log(e_j/3)` — a **soft repulsion
from the rest**. This is exactly SDA's "attraction to the copy carrying this PSV, repulsion from the others,"
now a continuous log-likelihood contribution instead of a ±1 edge. A read spanning **no distinguishing** bubble
between two copies has `log L_{rk}` equal across those copies (its evidence factors are shared): SDA's
*un-attractable* read, the K-frontier of `bench/THEORY.md` §5. `L_{rk}` is the existing `ReadEvidence.logl`
(Task 1's `read_copy_evidence`); no new emission model is introduced.

### 1.4 The E-step is SDA's read-partition made soft; the M-step estimates path usage

Introduce copy **abundances** `π = (π_1,…,π_K)`, `π_k ≥ 0`, `∑_k π_k = 1` (path-usage fractions across the VG's
copy-paths). The mixture responsibilities are the **soft** version of SDA's WhatsHap read-partition:

- **E-step** (`e_step`): `γ_{rk} = π_k L_{rk} / ∑_{j} π_j L_{rj} = softmax_k(\log L_{rk} + \log π_k)`. Where SDA
  assigns each read **hard** to one PSV cluster, the E-step gives a **fractional** read-to-copy-path assignment —
  a *soft WhatsHap*. A read carrying a copy's private allele gets `γ` concentrated on that copy (hard attraction
  in the limit `e_j → 0`); an un-attractable read keeps its mass spread.
- **M-step** (`m_step`): `π_k = \tfrac{1}{N}∑_r γ_{rk}` — re-estimate **copy-path abundance** (mean
  responsibility = path-usage fraction). SDA has no abundance parameter (it assembles, it does not quantify);
  this is the quantity the VG layer was missing and the EM adds.
- **Convergence** (`loglik`): the observed-data log-likelihood
  `ℓ(π) = ∑_r \log ∑_k π_k L_{rk}` is **non-decreasing** each sweep (the EM ascent guarantee; the tested monotone
  invariant), and iteration stops at `Δℓ < ε(1+|ℓ|)`.

*(Copy-path refinement — re-estimating `θ_k` from γ-weighted reads, the direct EM analog of SDA re-clustering
its PSVs — is deferred; here `θ` is fixed from the VG and only `π` + soft assignments are estimated. This is
without loss of generality for the consistency statement below, which fixes the true paths `θ*` and asks whether
`π` and the read assignments are recovered.)*

**In one line:** the EM is SDA's PSV correlation-clustering with the ±1 signed edges replaced by the
substitution log-likelihood, the hard WhatsHap partition replaced by soft responsibilities, and a copy-abundance
parameter added — i.e. the **maximum-likelihood soft relaxation** of the SDA heuristic.

---

## §2 The model — a finite mixture over copy-paths

**Object.** One PSV-aware VG per family: bubbles `[m]` (PSV columns, alphabets `A_j`), `K` copy-paths
`θ_1,…,θ_K ∈ ∏_j A_j` (from the VG). A read `r` has a hidden **copy-of-origin** `z_r ∈ {1,…,K}` and observations
`obs(r) ⊆ [m]`, `r(j) ∈ A_j`.

**Generative model.**
1. draw origin `z_r ∼ Categorical(π)` (copy `k` with probability `π_k`);
2. draw the spanned bubble set `obs(r)` from a **coverage law** independent of `z_r` given the family backbone
   (which bubbles a molecule spans is a property of read length/position, not of copy identity — the standard
   long-read spanning model of `bench/THEORY.md` §6);
3. emit `r(j) ∼ q_j(\cdot \mid (θ_{z_r})_j)` independently across `j ∈ obs(r)`.

The marginal read law is the **finite mixture**
`P(r) = ∑_{k=1}^K π_k L_{rk}` with components `L_{rk} = ∏_{j∈obs(r)} q_j(r(j)\mid(θ_k)_j)`.

**Parameters.** `π` (estimated) and `θ` (fixed from the VG here). **Assumptions.** (A1) *error-free core /
well-specified copies*: each read originates from exactly one of the `K` copy-paths in the VG (completeness —
`origin(r) ∈ C`, the standing hypothesis of `bench/THEORY.md` Theorem 4; the reference-absent/O4 escape is
handled by abstain-and-re-thread, not here). (A2) *bounded, known per-column error* `e_j < (|A_j|-1)/|A_j|` so
the matched allele is the modal emission. (A3) *coverage law independent of origin* (step 2). These are the same
assumptions the combinatorial theory runs on, now given a probability measure.

---

## §3 Theorem (consistency of the EM/MLE)

Write `D_{ij} = {d : (θ_i)_d ≠ (θ_j)_d}` for the distinguishing bubbles (columns `d ∈ [m]`) of copies `i≠j`, and
`δ(r) = \min_{k≠b}|\mathrm{obs}(r)∩D_{bk}|` for the number of distinguishing bubbles read `r` spans against its
closest competitor (`b` = its MLE copy) — the same `δ` as `bench/THEORY.md` §5b. The per-read identifiability
certificate is `min_p(r) = ε^{δ(r)}` (Theorem 4(i)); a read is `Certified` iff `min_p(r) < α/(K−1)` and
`SoftZone` otherwise (`label_read`).[^bonf]

[^bonf]: **Denominator convention.** The shipped `label_read` uses `α/(K−1)` with `K` = number of copies —
a per-copy Bonferroni correction over the `K−1` competitor copies of a read's MLE copy. `bench/THEORY.md`
Theorem 4 writes the threshold as `α/(n−1)`. These are *not* silently the same symbol: the `K−1` (per-copy)
convention is what the code applies, and the `n−1` in THEORY.md Thm 4 is flagged here as the notation to
reconcile — do not equate them without checking which population (`K` copies vs `n` reads/comparisons) each
Bonferroni family is over.

> **Theorem (EM consistency in the identifiable regime).** Suppose the family satisfies **Strong Separation**
> (`bench/THEORY.md` §5): for every pair `i≠j`, every read of copy `i` conflicts with every read of copy `j` —
> equivalently, every copy pair is separated by a bubble that reads actually span, so `δ(r) ≥ 1` (i.e.
> `min_p(r) < 1`) is *achievable* for every cross-copy comparison, and `min_p < α/(K−1)` is achievable at
> sufficient coverage. Then:
>
> **(a) Identifiability.** The mixture `∑_k π_k L_{rk}` is identifiable — and, because `θ` is fixed/known from
> the VG, the components are **pre-labelled** (each `L_{k·}` is tied to a named copy-path `θ_k`), so there is
> **no relabelling ambiguity**: `π` and the per-copy assignments are uniquely determined by the read distribution.
>
> **(b) MLE consistency.** As the per-copy coverage `n → ∞` (reads per copy grow), the maximum-likelihood
> estimate is strongly consistent, `π̂_n → π*` almost surely, and the EM sequence from a generic start converges
> to it.
>
> **(c) Assignment consistency.** For every **identifiable** read (`δ(r) ≥ 1`, i.e. `min_p(r) < 1`), the MAP
> assignment `ẑ_r = \arg\max_k π̂_k L_{rk}` equals the true origin `z*_r` for all `n` large enough
> (`ẑ_r → z*_r`), and the soft posterior `γ_{r,z*_r} → 1`.
>
> **(d) The non-identifiable class is honest.** For reads with `δ(r)=0` (`min_p(r)=1`, `SoftZone` — the K=0 /
> K-frontier mass), the posterior stays at the prior `γ_{rk} = π_k` over the tied copies for all `n`; the EM
> makes **no hard call** and never imposes a `1/k` split.

This is exactly the regime `bench/THEORY.md` proves **Strong Separation ⟹ the true copy set `C*` is the unique
minimum copy cover, `MCC = χ(H) = K`** (Lemma 1 + Theorem 2); the EM converges to that cover (§5 below).

---

## §4 Proof sketch

The two ingredients are **finite-mixture identifiability** and **MLE consistency** (Redner & Walker, *SIAM
Review* 26:195–239, 1984), specialized to the discrete PSV emission; the identifiability partition is the
`min_p` per-read certificate, which is SDA's attraction/repulsion separability.

**(a) Identifiability.** A finite mixture `∑_k π_k f_k` is identifiable iff the component laws `{f_k}` are
**distinct and linearly independent** and no `π_k = 0` (Teicher 1963; Redner–Walker 1984 §2); with `θ` fixed the
components are pre-labelled, so identifiability here carries **no label-permutation caveat**. Here
`f_k = L_{k\cdot}` is the read law under copy-path `θ_k`. Under Strong Separation each pair `i≠j` has a
distinguishing bubble `d ∈ D_{ij}` that spanning reads observe with positive probability (assumption A3 + the
"reads actually span" clause); at that bubble `q_d(\cdot\mid(θ_i)_d) ≠ q_d(\cdot\mid(θ_j)_d)` because the matched
allele differs and `e_d < (|A_d|-1)/|A_d|` (A2) makes the modal base identify the copy. Hence `f_i ≠ f_j`.

For **linear independence** we do not appeal to "distinct product-multinomial laws are independent" — that is
false in general for `K ≥ 3`. Instead, evaluate each component law at the copy-paths themselves and consider the
`K×K` matrix `M = [\,f_k(θ_j)\,]_{k,j}`, where `f_k(θ_j) = ∏_d q_d((θ_j)_d\mid(θ_k)_d)` is the probability copy
`k` emits exactly the allele-vector `θ_j` over the columns. Each `f_k` **peaks at its own centre** `θ_k`: for the
full spanned column set, `f_k(θ_k) = ∏_d (1−e_d)` while any off-diagonal `f_k(θ_j)` (`j≠k`) loses a factor
`(1−e_d) → e_d/(|A_d|-1)` at every distinguishing column `d ∈ D_{kj}` (nonempty under Strong Separation). With
`e_d < (|A_d|-1)/|A_d|` (A2) each such factor strictly shrinks the product, so the diagonal entry strictly
dominates its row: `f_k(θ_k) > ∑_{j≠k} f_k(θ_j)` holds once separation is strong enough (equivalently, `M` is
**strictly diagonally dominant** in the well-separated regime). A strictly diagonally dominant matrix is
invertible (Levy–Desplanques), so the vectors `{f_k}` are linearly independent. The mixture is therefore
identifiable — the same content as Strong Separation making the copies distinct, spanned allele-vectors
(`bench/THEORY.md` §5, biological reading: "PSVs exist **plus** dense read coverage").

**(b) MLE consistency.** With identifiability (a), assumptions A1–A3, a compact simplex parameter space for `π`,
and the discrete (hence bounded, continuous-in-`π`) log-likelihood, the conditions of Redner–Walker 1984 (Wald
consistency of the MLE for identifiable mixtures) hold: `π̂_n → π*` a.s. as `n → ∞`.

*Global (not merely local) EM convergence — and why it is available here.* Finite-mixture EM is **not** globally
convergent in general: the observed-data likelihood of a mixture is typically multimodal and EM can stall at a
local maximum (this is the standard Redner–Walker caveat — EM converges to a *consistent root*, but only
**locally**, from a start in its basin). One must not argue "the global maximizer is unique, so a generic start
lands in its basin" — that is a non-sequitur for mixtures. The correct argument here rests on a special feature
of **this restricted problem: `θ` is fixed**, so only the mixing proportions `π` are estimated. Then

$$
\ell(\pi) \;=\; \sum_r \log \sum_k \pi_k L_{rk}
$$

is, for each read `r`, the logarithm of a function `\sum_k \pi_k L_{rk}` that is **affine in `π`** (the `L_{rk}`
are fixed constants when `θ` is fixed). A log of a nonnegative affine function is concave, and a sum of concave
functions is concave, so `ℓ(π)` is **concave on the probability simplex**. For fixed `θ` the mixture-proportion
EM is exactly coordinate ascent / an MM ascent on this concave `ℓ`; on a concave objective over a convex compact
set every ascent sequence with the EM fixed-point/KKT stationarity converges to the **global** maximizer (unique
whenever `ℓ` is strictly concave, i.e. under the identifiability of (a) with all `π*_k > 0`). Hence, **with `θ`
fixed there is no local-optimum problem**: the EM ascent (`loglik` non-decreasing, §1.4) reaches the global MLE
`π̂_n`, and by Redner–Walker `π̂_n → π*` a.s.

Two honest scope statements follow. **(i)** This global-convergence guarantee **relies on `θ` being fixed** — it
is the concavity of `ℓ` in `π` alone that removes the local maxima. **(ii)** If `θ` were **also** estimated (the
deferred copy-path refinement of §1.4 — the direct analog of SDA re-clustering its PSVs), the joint objective
`ℓ(π,θ)` is no longer concave and the full finite-mixture local-maximum problem returns; only the *local*
convergence-to-a-consistent-root guarantee of Redner–Walker 1984 would then hold, from a start in the true root's
basin. The shipped engine estimates `π` with `θ` fixed, so it lives in the concave, globally-convergent case.
(The `K=0`/flat directions where `ℓ` is constant in some `π`-direction — assumption A1's resolvable core excludes
them from `π` estimation; they surface as part (d).)

**(c) Assignment consistency.** Fix an identifiable read `r` (`δ(r) ≥ 1`). For any wrong copy `k ≠ z*_r`, the
likelihood ratio telescopes over the `|obs(r)∩D_{z*_r,k}| ≥ δ(r) ≥ 1` distinguishing bubbles the read spans:

$$
\frac{L_{r,z^*_r}}{L_{r,k}}
\;=\; \prod_{j\in\mathrm{obs}(r)\cap D_{z^*_r,k}} \frac{q_j(r(j)\mid(\theta_{z^*_r})_j)}{q_j(r(j)\mid(\theta_k)_j)}
\;=\; \left(\frac{1-e}{e/3}\right)^{\!\ge\,\delta(r)} \;>\; 1,
$$

since `r(j) = (θ_{z*_r})_j ≠ (θ_k)_j` at each such bubble and `(1−e)/(e/3) > 1` (A2). The likelihood **already**
favours the truth (this is the `min_p(r) = ε^{δ(r)} < 1` certificate, Theorem 4(i)). The only way a wrong prior
could flip the argmax is if `π̂_k / π̂_{z*_r}` exceeded this ratio; by (b) `π̂ → π*` with all `π*_k > 0` bounded
away from `0`, so for `n` large the prior factor is bounded and `π̂_{z*_r} L_{r,z*_r} > π̂_k L_{r,k}` for every
`k`. Hence `ẑ_r = z*_r` and `γ_{r,z*_r} → 1`. This is Theorem 4(ii) (soundness — under completeness the unique
consistent copy is the origin, so the argmax is correct), promoted to a limit statement via consistent `π̂`.

**Coverage attribution (what the limit is, and is not, buying).** The likelihood-ratio factor
`((1−e)/(e/3))^{≥δ(r)}` is a **per-read** quantity governed only by that read's distinguishing count `δ(r)` and
the error rate `e`: it does **not** grow with coverage. More coverage means more *reads*, not more distinguishing
columns *per read* — `δ(r)` is fixed by which bubbles the single molecule `r` spans. Consequently the per-read
identifiable correctness is a **`δ`/`e` property**, and under assumption A1 (the error-free core, `e → 0`) it is
**exact at any coverage**: for `e → 0` the ratio is `∞`, the prior cannot flip it, and the argmax is correct with
`n = 1` copy-read as readily as with `n = 100`. What coverage `n → ∞` actually buys is the **abundance /
family-level** consistency `π̂ → π*` of (b); it enters (c) only through the *prior factor* `π̂_k/π̂_{z*_r}`, which
matters solely in the noisy `e > 0` regime where a badly-estimated prior at tiny `n` could momentarily overpower
a small-`δ` read. So: assignment accuracy on identifiable reads is a `δ/e`(A1) effect (high at all coverages
under A1); the coverage effect is `π̂ → π*` and the abundance L1 error → 0.

**(d) Non-identifiable class.** If `δ(r) = 0` then `obs(r)∩D_{bk} = ∅` for `≥ 2` copies (Theorem 4(iii)): read
`r` spans no bubble distinguishing them, so `L_{rk}` is **equal** across those copies, and the E-step gives
`γ_{rk} = π_k` over the tied set — the posterior is the prior, for every `n`. This is the `min_p(r)=1` tied
certificate, and at the family level it is `bench/THEORY.md` §6b (Tier-3 unidentifiability): with copies
identical over the observed bubbles the mixture likelihood is **constant on the abundance simplex**, so RNA
contributes nothing to the per-copy direction — the EM correctly returns the prior, never a hard `1/k`. `∎`

**The identifiability partition = the `min_p` certificate = SDA's separability.** Parts (c)/(d) split the reads
by `δ(r) ≥ 1` vs `δ(r) = 0`, i.e. by `min_p < 1` vs `min_p = 1`, i.e. by whether the read carries a
distinguishing PSV (SDA-attractable) or not (SDA-un-attractable). The consistency holds **exactly on the
attractable set**; the un-attractable set is abstained. This is the whole content of the theorem.

---

## §5 Tie to `bench/THEORY.md` — EM converges to the unique minimum cover, and abstains at the K=0 floor

The consistency theorem lives directly on top of the combinatorial results in `bench/THEORY.md` (the
`copy_assignment_theory` section; formerly `bench/copy_assignment_theory.md`). Named ties:

- **Lemma 1 (MCC = χ(H)).** Copy covers of the reads are exactly proper colourings of the read-conflict graph
  `H`; the minimum copy cover equals the chromatic number. The EM's hard limit (part (c), `γ` → indicators)
  partitions reads into per-copy classes — a `χ(H) = K`-colouring of `H`.

- **Theorem 2 (recovery under Strong Separation, all K).** Under Strong Separation the true copy set `C*` is the
  **unique** minimum copy cover, `MCC = K`. The EM's identifiable regime is **this same hypothesis**: "every
  copy pair separated by a spanned bubble" is Strong Separation verbatim ("every cross-copy read pair conflicts";
  §5 biological reading = PSVs + dense coverage). By §4(c) the EM's MAP partition converges to `C*` — i.e. the
  soft SDA relaxation **converges to the provably-unique cover**.

- **Theorem 3 (`RECOVER`, poly-time).** Under Strong Separation the connected components of the compatibility
  graph `H̄` are exactly the true read-classes. The EM's hardened responsibilities reproduce that partition (the
  reads of one copy are mutually compatible → one component → one copy-path); the EM is the statistical route to
  the same object `RECOVER` computes combinatorially, without solving the NP-hard cover.

- **Theorem 4 (Bridge — `min_p` is a sound per-read certificate).** `min_p(r) = ε^{δ(r)}`; `δ ≥ 1 ⟺` unique
  consistent copy `⟹` correct assignment (soundness, given completeness); `δ = 0 ⟺` tied `⟹ min_p = 1`. This is
  the certificate the EM attaches per read via `label_read` (`Certified` = §4(c) consistency applies; `SoftZone`
  = §4(d) abstention). Theorem 4(iv) — `K_{ij}=0 ⟹` all reads of the pair have `δ=0`, `min_p=1`, tied — is
  precisely §4(d).

- **Theorem 7 (integrality bridge — §5c facility-location capstone).** Under Strong Separation every read is
  consistent with exactly one copy (`|N(r)|=1`), the MWCA LP is **integral**, greedy is exact, and `min_p` is
  the complementary-slackness dual witness. The EM's soft assignment, in the same regime, converges (§4(c)) to
  that integral optimum: the soft relaxation and the LP relaxation agree at their common integral point. So the
  EM is not a third, ad-hoc solver — it is the mixture-model face of the same `NP-hard (Thm 1/5) → (1−1/e)
  LP-rounding (Thm 6) → integral under Strong Separation (Thm 7) → per-read `min_p` (Thm 4)` object.

- **The K=0 floor = SDA's ">100 kb reads" = `SoftZone`.** When a copy pair shares every spanned bubble
  (`D_{ij}=∅` over `obs`, the K=0 vertex of the K-bound Corollary and the IDENTIFIABILITY_LIMITS §"K=0
  boundary"), no read is attractable: `δ=0`, `min_p=1`, the mixture likelihood is flat (§6b), and every read is
  `SoftZone`. This is **the same wall** SDA names — *"virtually identical duplications … will require even longer
  reads (>100 kb)"* — surfaced by the EM as a posterior-invariant class (`γ = π`), **never a hard 1/k**. The EM
  reaches the identifiability boundary and certifies it, exactly the Canzar-aesthetic result the thesis targets.

---

## §6 Consequence — explaining SDA's 91–93% floor and predicting the coverage sweep

**Why SDA floors at 91–93%.** SDA runs *hard* correlation clustering (NP-complete, 15 random inits) and
*forces* a cluster label on every PSV/read — including the un-attractable, K-frontier ones that carry no
distinguishing signal. Two things cost it accuracy: (i) the heuristic can get stuck in a non-global correlation
clustering, and (ii) hard-calling the genuinely unidentifiable mass **must** misassign a fraction of it. The
consistency theorem decomposes exactly this:

- On the **identifiable** set (Strong-Separated reads, `δ ≥ 1`), the ML soft relaxation is **consistent**: its
  per-read accuracy is high at every coverage (a `δ/e` property, §4(c) coverage-attribution note; → 100% under the
  error-free core A1) — it does not floor. SDA's shortfall *on this set* is the heuristic gap, which the
  likelihood ascent (monotone `ℓ`) closes. Note this is the EM's **identifiable-set** accuracy curve.
- On the **unidentifiable** set (`δ = 0`, the K-frontier / K=0 mass), no method resolves the reads (§4(d),
  Theorem 2 K-bound); SDA hard-calls it and eats a misassignment rate, whereas the EM abstains (`SoftZone`).

**What the theorem does — and does not — supply about the number.** The theorem explains why a floor **of this
KIND** exists: a method that hard-calls every read scores ≈100% on the identifiable set and pays a forced error
rate on the hard-called unidentifiable residue, so its *overall* accuracy is pinned below 100% by that residue.
The theorem does **not** derive the specific value **91–93%** — that value is set by the *instance-specific
unidentifiable fraction* (how much K-frontier/K=0 mass a given family carries, and how SDA splits it), a quantity
the theory does not supply. So the mechanism is derived; the number is not. Crucially, keep the two curves
**distinct**: the EM's **identifiable-set** accuracy → 100% is a *different curve* from SDA's **overall** 91–93%,
which mixes in the hard-called unidentifiable mass — they are not the same measurement and must not be equated.

**Prediction (Task 6 coverage sweep, `bench/em_coverage_sweep.py` → `bench/EM_COVERAGE_SWEEP.md`).** On a
planted sim genome with known `θ*, π*, z*`, as per-copy coverage grows `{1,2,5,10,20,50,100}×`, the sweep has
**two distinct axes** (the §4(c) coverage attribution):
- *(per-read axis, `δ/e`, not coverage-driven)* assignment accuracy on **identifiable** reads is **high at every
  coverage** and **→ 100%** under the error-free core A1 — this is the EM's identifiable-set curve, which
  explains (but is not numerically equal to) SDA's 91–93% overall floor;
- *(coverage axis)* abundance error `‖π̂ − π*‖₁ **→ 0`** as coverage grows (consistency of `π̂`, §4(b)) — this is
  the genuine coverage effect;
- K=0 families **stay `SoftZone`** at every coverage (the boundary is a boundary, not a coverage artifact).

These are the falsifiable content of the theorem; Task 6 is its demonstration.

**Which bubbles enter the graph (the per-column filters).** The theorem assumes the columns `[m]` are *true*
PSVs, not sequencing error — and `bench/THEORY.md` §5c warns that the **raw** allele-disagreement graph is
error-inflated (colouring `≈ 3×` the true `K`), so column selection is load-bearing. The bubbles are gated by
the same per-column filters already in the cascade, reused unchanged:
- **IsoCon** (`reference_isocon_sahlin`) — the per-position **real-variant-vs-error** significance test decides a
  column is a genuine PSV rather than error, keeping error edges out of the graph (the honest replacement for
  solving `χ` on the raw graph);
- **Sudmant 2010 SUN** (`reference_sudmant_2010_sun`) — single-position **private** markers, the strongest
  bubbles: a read over a SUN column is single-read pinned (`|N(r)|=1`, `bench/THEORY.md` §5·SUN Lemma), the
  concrete `δ`-contributing attraction the theorem's part (c) relies on;
- **Clair3-RNA** (`reference_clair3_rna`) — flags A→I editing so edited sites are not mistaken for PSVs.

The EM applies the IsoCon per-position and Clair3 A→I-editing column filters (`em_assign_family` now
computes `detect_editing_columns` and passes it into `read_copy_evidence`, so an editing column is
downweighted in the EM's likelihood exactly as it is in the hard gate). It does **not**, however, reuse
everything the hard gate does: the shipped `--em` path (`em_assign_family`) is PSV-only and threads no
copy-specific junctions and no per-base quality (`ReadFeatures::junctions`/`psv_qual` and
`CopyProfile::junctions` are left empty), so its per-read labels can differ from the hard
`.assignments.tsv` gate on reads whose call depends on junction or per-base-quality evidence. This does
not affect the abundance/consistency result above: junctions, quality, and editing change only `min_p`
and per-column weights in the likelihood, not the shape of the abundance fixed point the EM/M-step
converges to.

---

## References

- Vollger, Dishuck, Chaisson, Eichler, et al. **"Long-read sequence and assembly of segmental duplications."**
  *Nature Methods* 16:88–94 (2019). `reference_sda_vollger`. — SDA's PSV graph + correlation clustering; the
  91–93% floor; the ">100 kb reads" identifiability statement.
- Redner, R. A. & Walker, H. F. **"Mixture densities, maximum likelihood and the EM algorithm."** *SIAM Review*
  26(2):195–239 (1984). — finite-mixture identifiability + strong consistency of the MLE; EM convergence.
- Teicher, H. **"Identifiability of finite mixtures."** *Ann. Math. Statist.* 34:1265–1269 (1963). — linear
  independence ⟹ identifiability of finite mixtures.
- Sudmant, P. H., et al. **"Diversity of human copy number variation and multicopy genes."** *Science*
  330:641–646 (2010). `reference_sudmant_2010_sun`. — the SUN (single-position private allele).
- Sahlin & Medvedev. **IsoCon** (2018). `reference_isocon_sahlin`. — per-position real-vs-error PSV test.
- `bench/THEORY.md` — Lemma 1 (MCC = χ(H)), Theorem 1 (NP-hard), Theorem 2 (Strong Separation ⟹ unique minimum
  cover), Theorem 3 (`RECOVER`), Theorem 4 (`min_p` bridge), Theorems 5–7 (facility-location capstone,
  integrality bridge), §5 Proposition (K-frontier), §5·SUN, §6b (Tier-3 unidentifiability), IDENTIFIABILITY_LIMITS.
</content>
</invoke>
