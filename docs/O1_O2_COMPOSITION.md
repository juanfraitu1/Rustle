# O1 and O2 as a composition — the statement of the method

> **Status: FRAMING DOCUMENT, not a result.** Every number here is cited to code or ledger.
> Written 2026-09-03. ⚠ Line numbers are from a DIRTY working tree on `dna-from-genome`
> (`M` on `copy_assign.rs`, `copy_assign_pipeline.rs`, `copy_graph.rs`, +20 more) — re-anchor
> before quoting in the thesis.

## 0. The sentence that must be said first

**O1's vertex set IS O2's path set.** That is a *composition*, not an identity, and the
distinction is the whole point:

| | O1's graph | O2's graph |
|---|---|---|
| vertices | **loci** (one rep transcript = one intron chain) | **allele-nodes** at PSV columns |
| edges | pairwise homology `E_r` | — (see §3: there is no edge set) |
| content | no sequence on the node | each copy is a path of alleles |
| built by | `family_definition.rs` / `mcl_families` | `BubbleGraph::from_copies`, `copy_assign.rs:53` |

⛔ **Do NOT say "O1 builds the variation graph and O2 reads paths through it."** It is false.
O1 emits node *sets* over a locus-level similarity graph with no sequence content. O2 builds the
variation graph *over* that set. Say: **O1 selects the node set; O2 builds the graph over it.**

⛔ **Do NOT say "a family is one variation graph."** `NEGATIVE_RESULTS_REGISTER.md:472` kills it —
*"a VG presupposes its members. Circular. Keep the homology/quasi-clique graph for the definition;
VG is representation."* And `:1086` kills the detector form. **Volunteer these rows before he finds
them.** The register's §ADVISOR table also records that "O1 and O2 can be presented as one decision
rule" was flagged, by us, as exactly what would irritate him.

## 1. What a family actually is

Not a connected component. A **certified block of a partition of a connected component** of
`(V, E_r)`: a vertex set `F` with

1. `F` inside one connected component of `E_r`;
2. induced density `ρ_in(F) = 2|E(F)|/(|F|(|F|−1)) ≥ γ = 0.20`, **or** `|F| ≤ 2`;
3. `|distinct_loci(F)| ≥ 2`.

`family_definition.rs:33-36` (γ), `:171-176` (density), `:52-84` (`distinct_loci`,
`LOCUS_OVERLAP = 0.50` reciprocal), `:368-400` (the recursion).

⚠ The exact max-γ-quasi-clique partition is NP-hard, so what ships is a **verifiable certificate
per block**, never an approximation ratio. The splitter (`family_split.rs:533-556`) is a heuristic
cascade — Louvain at resolutions [1,2,4,8], then components, then `members[..n/2] / members[n/2..]`,
an index-order halving. **Disclose the halving before he greps for it.**

## 2. What assign-or-abstain actually is

`copy_assign.rs:449-456`, the default branch (`use_margin_gate = false`):

```
thr        = α / max(n−1, 1)
resolvable = min_p < thr          where  min_p = ∏ εⱼ  over distinguishing observations
```

Defaults: α = 1e-3, `error_rate` 0.003, `junction_err` 1e-4, `edit_rate` 0.2 (`copy_assign.rs:210-222`).

⭐ **The true and defensible claim:** `min_p` does **not** depend on *which* alleles the read
carries — only on *which distinguishing sites it spans* and how confidently they were read
(`copy_pair_significance`, `copy_assign.rs:247-293`). Resolvability is therefore separated from
evidence: **`Tied` vs `resolvable` is structural + quality-weighted; `Assigned` vs `Ambiguous` is
the evidence-dependent layer.** Theorem 4 (`THEORY.md:678-733`) proves `min_p = ε^δ(r)`, machine-
checked exhaustively for K≤3, m≤3, |A|≤3 — 3,400 copy-sets / 67,320 reads
(`bench/bridge_theorem_check.py`).

⛔ **The threshold moved; it did not vanish.** Under uniform ε:
`assignable ⟺ δ(r) > log(α/(n−1)) / log ε`. Concretely — at flat error one PSV column is **Tied at
every n ≥ 2**; two columns resolve for n ≤ 1000; one junction for n ≤ 10; one column at Phred ≥ 40
for n ≤ 30. ⚠ The fraction of reads on the flat-error branch **has not been measured** on the
gorilla BAM. Measure it before quoting "one column" or "two".

⛔ **Recombinant reads are force-assigned.** A read carrying copy A's allele at one bubble and copy
B's at another belongs to no path — Theorem 4's precondition fails. `recombinant_abstain.rs` calls
itself "the DEFAULT-ON gate leg" (`:18`) but `MODULE_STATUS.md:112,157` records it **TEST-ONLY**:
its only non-test caller is `o2_materialize.rs:866`, which no binary imports, and
`RUSTLE_NO_RECOMBINANT_ABSTAIN` is a **vacuous opt-out**. State the claim as *"exactly one path is
consistent"* — a statement about the given path set, not about the molecule.

## 3. What the object is, honestly

`BubbleGraph { bubbles, n_copies }` (`copy_assign.rs:41-72`, **SHIPPED-DEFAULT**) is real and runs
by default; `copy_assign.rs:327` says verbatim *"thread the read as a WALK through the bubbles it
spans"*. The pipeline header (`copy_assign_pipeline.rs:5-9`) already states the Canzar flip. **So
this restatement is a renaming of running code — that is its strength and its limit.**

Four concessions that must travel with it:

- **It has no edge set.** One field, `bubbles: Vec<Bubble>`, plus `n_copies`. A referee will say:
  *you built a genotype matrix and called it a graph.* The junction term is computed **outside** it
  (`copy_assign.rs:351-363`) — so the "RNA family graph" omits exon/intron structure, the one signal
  that makes this transcriptomic.
- **Star alignment, not an MSA.** Every copy is aligned to `copies[0]` (`copy_assign_pipeline.rs:393`).
  A difference between copies 2 and 3 where copy 0 is gapped yields **no column**. **Substitution
  bubbles only** — indels require both bases ACGT (`:443-445`). ⚠ No test asserts column-set
  invariance to copy order; O1 treats seed-invariance as load-bearing, O2 has no analogue.
- **The bubbles are read-filtered by default.** `read_supported_columns` needs ≥2 alleles at ≥2 reads
  each above coverage 4 (`copy_assign_pipeline.rs:641-673`, wired `:1474-1488`). Near-inert —
  36,414 → 36,410, **4 dropped = 0.011%** — but it means the site set is not a pure sequence object.
- **The graph-shaped modules do not ship.** `o2_materialize`, `o2_columns`, `o2_margin_gate`,
  `positional` are TEST-ONLY; `copy_graph.rs` (GFA) needs `--phase`; `vg_realign` and `linearize`
  are opt-in. **14 of 52 modules reachable at defaults.**

⛔ **`family_graph.rs` is the trap.** Its header describes exactly the object a reader expects —
"union of per-copy splice graphs; nodes are exon-equivalence classes" — across 2,214 lines, and it
is **DEAD** (`MODULE_STATUS.md:143`, no production caller). It is the most seductive wrong
description available, because the header reads perfectly.

## 4. ⭐ The one consequence that is new

**Assignability is monotone-decreasing in family size, and that makes O1 error an O2 error channel.**

Because `thr = α/(n−1)`, adding a copy to a family converts previously-`Assigned` reads to `Tied`
**with no change to any read or any sequence**. So an O1 over-merge does not merely mislabel a
family — it *provably raises O2's abstention rate*, by a computable factor.

In the two-stage pipeline view this is invisible. In the composition view it is one line of
arithmetic, and it is **checkable**:

| family | copies | assignments | assigned | tied | ambiguous |
|---|---|---|---|---|---|
| MCLFAM2 | 3 | 1,242 | 467 | 768 | **7 = 0.56%** |
| MCLFAM1 | 38 | 4,036 | 272 | 3,407 | **357 = 8.8%** |

`o1_ledger.md:10243-10244`. The threshold's share of abstentions grows ~16× from 3 copies to 38,
and only **8 of 38** copies receive any decisive read. ⚠ **Never claim the threshold is inert
without naming the k it was measured at.**

⚠ This argument cuts against the sensitive 0.60 tier — the tier where the measured Alu false merges
live (§6cr: genuine paralogy 0.941–0.980, Alu artefact 0.710–0.852, 22/22; `-x asm20` finds zero).
Be willing to follow it.

⚠ **Discrepancy to reconcile before any slide:** `O1_TO_O2_BRIDGE.md §4` reports MCLFAM2 as
498/726/18; the ledger reports 467/768/7. Both sum to 1,242; neither names the arm (the bridge
runnable passes `--vg-realign-correct`, the ledger's run had it off).

## 5. Positioning against Canzar 2016 — a contrast, not a concession

His formulation is **maximum facility location** with LP-rounding, 0.19-approximation and a matching
(1−1/e) hardness. That machinery is **necessary when the assignment is coupled** — when opening a
facility for one read changes another read's cost.

Here `L(j) = Σ_r log P(obs(r) | copy_{j(r)})` has **no term coupling two reads and no term coupling
two copies**: no opening cost, no cardinality constraint, no abundance prior on the hard path. The
objective **decomposes**, so **per-read argmax is exactly optimal** — no LP, no rounding, no
approximation ratio needed. And 1/k is avoided by *abstaining*, not by splitting mass.

**The honest framing: we identified the regime where his machinery is not required — and it is
degenerate precisely because O1 fixed the copy set first.** Corroboration: an EM changed **0 of
3,081** evidenced decisions.

⛔ **`grep -rn facility --include=*.rs src/` returns 0.** Say so before he asks. `THEORY.md` §5c
Theorems 5-7 are checked only by a pure-Python script that never touches a binary — present them as
a *programme*, never as a measured Rustle result. `ONE_METHOD.md` already retracts the
facility-location framing.

⚠ Not fully decoupled: the Bonferroni divisor `α/(n−1)` makes each read's decision depend on the
family's cardinality — which is §4.

## 6. His stated prior is satisfied — a free win

He is "suspicious of methods with no similarity threshold for merging graphs". **This method is not
threshold-free**: `E_r` requires identity ≥ 0.80 (asm20) or ≥ 0.60 (sensitive `-k11 -w5`) **and**
coverage ≥ 0.50 of the shorter, on one alignment record, **plus** γ ≥ 0.20 cohesion. Say it early;
it converts his prior into agreement.

Two caveats in the same breath:
- The coverage floor is **one-sided**. A 2,037 bp NPIPB6 fragment reaches coverage 0.948 against a
  38,653 bp chimeric node while touching 5% of it, dragging EIF3CL into NPIP. The two-sided fix
  passed 5/5 and 9/9 pre-registered **and was still refused a default** (D4 null, D5 deleted 18
  families). That refusal is the strongest anti-overfitting artifact in the project.
- The 0.60 sensitive floor is where the Alu false merges live. Present it as an open precision
  problem **with a diagnosis**, not as solved.

## 7. Threshold inventory — correct §1.2 before he builds it himself

`ADVISOR_QUESTIONS.md §1.2` says "the shipped rule has **four** free numbers". A grep finds **~25**
default-reachable constants that decide outcomes, and **six more arrived 2026-09-03** as CLI
defaults in `mcl_families.rs:38-66` (inflation 2.8, min_identity 0.70, min_cov_longer 0.30,
min_bp 300, min_size 3, min_reads 3). Restate as **"four in the edge rule; here is the full
inventory"**.

⚠ Coverage is charged on the **shorter** side in shipped `E_r` and on the **longer** side in the MCL
path (`annotation_families.rs:47-49`, deliberately — every adjudicated NPIP false merge was
Alu-mediated). **Do not present the two conventions as one rule.**

## 8. Prior art — currently absent from the entire document set

⚠ `grep -rin 'giraffe|snarl|vg toolkit|Garrison|Paten|Sirén'` over `docs/` and
`src/rustle/vg_family/` returns **nothing**. "Variation graph" appears 15 times, never beside a
citation. ⚠ All of the following are **UNVERIFIED, from memory** — check before use: vg (Garrison
2018); snarls/ultrabubbles (Paten 2018) — a bubble *is* a snarl, so `BubbleGraph` is a degenerate
special case of a published decomposition; Giraffe (Sirén 2021); GraphTyper/Paragraph/PanGenie;
RSEM/kallisto/salmon as the EM baseline "never 1/k" is defined against.

⛔ **The closest prior art is in our own tree, uncited.** "Assignable iff it crosses a site private
to one copy" is the **SUN** (singly-unique nucleotide, Sudmant/Eichler ~2010 — UNVERIFIED). This
repo ships `SUN_TSV_BASENAME = "sun_identifiability.tsv"` and computes `n_sun_bubbles`
(`o2_materialize.rs:71, :575-584`) — and that module is **TEST-ONLY**. The idea is fifteen years
old, we know it, and our implementation of it does not run.

**What is genuinely new, narrow enough to survive:** (a) the setting — spliced RNA, where a
copy-specific *junction* is a separating feature and the graph is over spliced products, not genomic
paths; (b) abstention as the *default* outcome with a per-read certificate and a Bonferroni
correction over the family's copies, rather than a marker-count heuristic; (c) the composition —
the same γ-quasi-clique that defines the family supplies the exact byte sequences O2 scores against
(`--copies-fa`, `O1_TO_O2_BRIDGE.md §3`; unique-mapper agreement 68/68).

## 9. What this document does NOT claim

- ⛔ Not "the pipeline runs on variation graphs". Only the in-memory `BubbleGraph` is default.
- ⛔ Not "O1 decides membership by sequence alone" unqualified. True of the **edges**; the vertices
  are assembled from reads at `NODE_MIN_READS = 2` (`denovo_assemble.rs:1301`), and the O1⊥O2 guard
  discloses two live read-derived uses in the node set — `locus_unique_mapper_counts` and
  `reads_distinguish` (`denovo_pipeline.rs:8724-8735`), measured 0/109 pairs and 0/451 chr1 loci.
  ⚠ And it is false outright for a **bare `copy_assign` run**, whose default membership oracle is
  the read-conflict graph `E_c` (`denovo_pipeline.rs:2006-2010`; `--homology-primary` default
  false). **Always name the producer.**
- ⛔ Not "O2 assigns better than minimap2". Reassignment headroom is ~0.1% of reads. **Defend O2 on
  abstention only.**
- ⛔ Not a PSV VCF. `genome_pos` is the first copy carrying the column; others sit 33 kb – 81 Mb
  away, so every ALT would be megabases from the asserted POS.
- ⚠ **The one non-circular O2 abstention result is not a run of the shipped gate.** The excision
  numbers (TPR 0.5066 / FPR 0.0280 / AUC 0.7995 at <50% migrants; 0.2404 / 0.0239 / 0.6918 pooled)
  are the per-read **robust-z of `de`**, not the α-certificate at `copy_assign.rs:455`. The stratum
  conditions on an unobservable. **Either re-run the excision labels through the shipped
  `assign_read` — the highest-value experiment available before the meeting — or say out loud that
  O2's abstention is validated only by proxy.** Also: an AUC needs a sliding score; a structural
  predicate gives one ROC point.
