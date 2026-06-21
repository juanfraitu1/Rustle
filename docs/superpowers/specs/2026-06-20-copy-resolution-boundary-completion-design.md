# Copy-resolution boundary completion — design

*Completes `bench/copy_resolution_boundary.md` into a working system with three focused, independently-testable
components: an exhaustive per-family K=0 census, a splice-divergence resolver (the Tier-2 lever), and the Tier-3
co-quantification proposition. Mostly Python bench tooling + one theory proposition; no Rust-pipeline
integration in this milestone.*

## Goal

Turn the measured copy-resolution boundary into reproducible, exhaustive artifacts: (1) a definitive per-family
K=0 label across all co-located multi-copy families, (2) an opt-in resolver that recovers the splice-divergent
subset of the K=0 residual, and (3) a formal statement of what is achievable for the irreducible core.

## Non-goals

- No Rust / `detect_and_assign` pipeline integration (the resolver stays a standalone bench tool).
- No full EM/Bayesian per-copy estimator for Tier-3 (the honest result is unidentifiability + aggregate
  identifiability + prior-conditioning — a proposition, not a solver).
- No new sequencing / external DNA data; the copy-number prior in Tier-3 is a parameter, not something we infer.

## Shared substrate

- BAM `/home/juanfra/winloci_scratch/GGO.bam` (FLNC, .bai), T2T reference `/home/juanfra/winloci_scratch/GGO.fasta`.
- Families: `bench/denovo_families.tsv` (1130 multi-copy families; member ids `DN_<chrom>_<start>_<n>`).
- Existing: `bench/copy_resolution_census.py` (per-pair classifier), `bench/copy_resolution_census.json`
  (258 assignment-relevant pairs), `bench/copy_resolution_boundary.md`, `bench/copy_assignment_theory.md` +
  `bench/copy_assignment_theory_checks.py`.
- K=0 definition: a co-located copy pair is K=0 iff its cross-mapping (MAPQ-0, same query_name at both loci)
  reads have `NM_A == NM_B` (the spliced read aligns equally to both copies' exons). The known K=0 panel
  (MAGEA, NC_073247.2): pair0 161251228-161257000~161458538-161464324; pair2 164381222-164384848~164442447-164446101;
  pair3 164397061-164401095~164426194-164430228 (inverted dups). Tier-2 splice-rescuable: pair0 (~33% reads),
  DNFAM451 (differential intron usage), DNFAM224/DMRTC (copy-private 489bp intron). Tier-3 irreducible: pair2/pair3.

---

## Component 1 — Exhaustive per-family census

**Files:** extend `bench/copy_resolution_census.py`; emit `bench/copy_resolution_census.tsv` (new, per-family) +
keep `bench/copy_resolution_census.json` (per-pair).

**What:** Extend the existing per-pair classifier from the 258 assignment-relevant pairs (>=3 cross-mapping
reads) to ALL co-located copy pairs across the 1130 families, lowering the read threshold to >=1 but recording
`n_xmap` so a confidence flag (e.g. `low_support` for `n_xmap < 3`) is retained. Aggregate pairs to a
per-FAMILY verdict and write one TSV row per co-located multi-copy family:

`family_id, n_copies, n_colocated_pairs, n_assignment_relevant, n_resolvable, n_k0, n_k0_strict, family_verdict, k0_pairs`

where `family_verdict` ∈ {`resolvable`, `k0`, `mixed`, `not_assignment_relevant`}: `not_assignment_relevant` if
no pair has cross-mapping MAPQ-0 reads; `k0` if all assignment-relevant pairs are K=0; `resolvable` if all are
resolvable; `mixed` otherwise. `k0_pairs` lists the K=0 copy-pair coordinates + gene names (from
`GGO_genomic.gff`). A summary block prints the family-level counts and the read-level resolvable fraction,
reproducing the established headline (≈77% resolvable pairs / ≈83% reads) on the assignment-relevant subset and
adding the exhaustive family-level numbers.

**Verification (check):** a deterministic check that (a) the known MAGEA family (DNFAM1) is labeled `k0`/`mixed`
with its pair2/pair3 K=0 pairs present, (b) the assignment-relevant-subset summary still reports
resolvable ≈ 77% / reads ≈ 83% (matching `copy_resolution_census.json`), (c) the per-family TSV is internally
consistent (`n_resolvable + n_k0 + ... = n_assignment_relevant`).

---

## Component 2 — Splice-divergence resolver

**Files:** create `bench/splice_divergence_resolver.py`.

**What:** For a K=0 copy pair `(locusA, locusB)` (handling the inverted-dup reflection — solve the axis so
A-donor maps to B-acceptor): (1) enumerate every read-supported splice junction at each locus (incl. alternative
donors/acceptors); (2) for each junction, extract the **splice-site dinucleotides** (donor `GT`, acceptor `AG`,
transcript-oriented) from `GGO.fasta` at copy A's coordinates and at copy B's homologous coordinates; (3) flag a
junction **copy-distinguishing** iff it is canonical (`GT-AG`/`GC-AG`) at one copy but degraded at the other's
homologous position (the intronic divergence at the boundary); (4) assign a read to the copy at which the
junction it uses is canonical. Also detect **differentially-used / copy-private introns** (a junction present in
one copy's reads with no homologous counterpart at the other, e.g. DMRTC's 489bp intron).

**Interfaces:** `resolve_pair(bam, fasta, locusA, locusB) -> {distinguishing_junctions: [...], n_reads: int,
n_resolved: int, resolved_fraction: float, per_read: {qname: copy}}`.

**Verification (check):** runs `resolve_pair` on the panel and asserts the established results: MAGEA pair0
resolved_fraction ≈ 0.33 (≥18 distinguishing junctions), DNFAM451 and DMRTC fire (>0 distinguishing or
copy-private junctions), MAGEA pair2 and pair3 resolved_fraction == 0 (splice-identical, GT-AG at both copies).
Deterministic.

---

## Component 3 — Tier-3 co-quantification proposition

**Files:** add a Proposition + proof to `bench/copy_assignment_theory.md` (Tier-3 / §6 area, after the Corollary);
add `check_tier3_coquant_unidentifiable` to `bench/copy_assignment_theory_checks.py`.

**What (the proposition):** Let a family have copies `C* = {c_1,…,c_K}` identical over the transcribed region
(K=0), and let the observed reads be `R` (each consistent with every copy). Define per-copy abundances
`a = (a_1,…,a_K)`, `a_k ≥ 0`, `Σ a_k = N` (the family aggregate). Under the standard mixture likelihood
`L(a) = Π_{r∈R} (Σ_k (a_k/N) · P(r | c_k))` with `P(r | c_k)` identical across `k` (copies identical), `L(a)`
is **constant** over the simplex `{a : Σ a_k = N}` — the per-copy split is **statistically unidentifiable** from
RNA. The aggregate `N = |R|` is identifiable (a sufficient statistic). Under a copy-number / dosage prior
`π(a)` (e.g. `a_k ∝ copy-number_k`), the MAP estimate is well-posed and equals the prior scaled to `N`; RNA
contributes nothing to the per-copy direction. Honest corollary: any reported per-copy split for a K=0 family is
prior-driven, with the full simplex as its (flat-likelihood) uncertainty set.

**Proof sketch (for the note):** since every `P(r | c_k)` is identical (call it `p_r`), each read's mixture term
is `Σ_k (a_k/N) p_r = p_r · (Σ_k a_k)/N = p_r`, independent of `a`. Hence `L(a) = Π_r p_r` is constant in `a`.
The likelihood factors through `N` only ⇒ `N` identifiable, `a` unidentifiable; with prior `π`, the posterior ∝
`π(a)` on the simplex ⇒ MAP = mode of `π` scaled to `N`. ∎

**Verification (check):** `check_tier3_coquant_unidentifiable` builds a small K=0 instance (≥2 identical copies,
reads consistent with both), evaluates the mixture log-likelihood at several distinct per-copy apportionments
(e.g. (N,0), (N/2,N/2), (0,N)), and asserts they are all equal (to floating tolerance) while the aggregate `N`
is fixed — the executable witness of unidentifiability. Reuses the existing reads/consistency helpers.

---

## Architecture & data flow

```
census (Comp 1)  → per-family K=0 labels
                      ├─ Tier-2 (splice-divergent) → splice resolver (Comp 2) recovers a fraction
                      └─ Tier-3 (splice-identical)  → co-quantification proposition (Comp 3): aggregate only
```

Each component is a self-contained unit with one responsibility and an executable check; they share the K=0
definition and panel but do not depend on each other's internals (the resolver and the proposition both consume
the census's K=0 labels but are independently testable on the fixed panel).

## Testing

Every component ships a deterministic check. The census check reproduces the headline numbers; the resolver
check pins the panel's per-family resolved fractions; the co-quant check is the flat-likelihood witness. The two
theory-checks live in `copy_assignment_theory_checks.py` (suite must still exit 0); the census/resolver checks
are runnable standalone (`python3 bench/<file>.py`).

## Risks / open points

- **Census threshold:** lowering to `n_xmap >= 1` admits noisier single-read pairs; mitigated by the
  `low_support` flag and by reporting the robust `>=3` summary alongside. The family verdict uses only
  assignment-relevant pairs.
- **Inverted-dup reflection axis** in the resolver must be solved per pair (A-donor ↔ B-acceptor); the workflow
  established the axes (pair0 322715544, pair2 328827293, pair3 328827289) — the resolver should compute the
  axis from the locus coordinates, not hardcode, and validate via the known panel.
- **Splice-site extraction strand:** the inverted copy is on the opposite genomic strand; donor/acceptor must be
  read transcript-oriented (reverse-complement for the inverted copy). The check on pair2/pair3 (must be 0)
  guards against a strand bug producing spurious distinguishing calls.
- **Tier-3 proposition** assumes the error-free identical-copy idealization (consistent with the rest of the
  note's error-free core); a footnote should note that read error perturbs `p_r` per read but not the
  `a`-independence, so the unidentifiability is robust to symmetric error.
