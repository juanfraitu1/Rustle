# E_r TE/repeat-exclusion gate: does soft-mask status resolve non-coding-paralog vs TE-junk?

**Scripts / artifacts** — `bench/family_er_te_gate.py` (deterministic, `PYTHONHASHSEED=0` pinned in-script),
`bench/family_er_te_gate.tsv` (one row per operating point), `bench/family_er_te_gate.json` (curve +
mask-discrimination + recommended point + off-curve proof + honest named confirmation + aligned cross-check
+ MAX-inflation), `bench/edge_bridge_mask.tsv` (per gene-pair edge: `core`, `mask_a`, `mask_b`,
`bridge_mask`). Reuses `bench/family_er_oracle_sweep.py` (gate/eval machinery, `compute_point`,
`prepare`) + `bench/family_er_pr.py` + `bench/genome_family_def.py` (γ-quasi-clique refine).
Reproduce: `PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_te_gate.py`. All three
artifacts are **byte-identical across independent launches** (md5 `family_er_te_gate.json`
`d344d710…`, `.tsv` `a66d7a71…`, `edge_bridge_mask.tsv` `7006feb0…`); the minimap2 aligned cross-check runs
`-t1` and is deterministic but **never feeds the gate**.

## Bottom line

The mask is a **real but weak, non-complementary** signal. It does **not** dissolve the open problem — it
**shifts the operating point** by a modest, genuine, off-the-no-mask-curve amount. The headline hypothesis
"genuine over-merge = TE/repeat-bridged (masked), real paralog = unique (unmasked)" is **only true for the
protein-*missed* TP subset**; it is **population-FALSE for the real protein-homologous paralogs**, which are
*as repeat-rich as the junk* and survive purely by the `in_ep` tautology. TEs at the RNA level are therefore
**partially handled**: the mask cleanly removes the flagrant zinc-finger / repeat-rich LOC over-merges
(mask ≥ 0.5) and keeps the low-repeat non-coding paralogs, but it **cannot** separate the hard residual —
near-identical unique-sequence LOC-family over-merges — and it wrongly cuts ~61% of the protein-missed TP it
was meant to rescue.

## The open problem (from `bench/FAMILY_ER_ORACLE_SWEEP.md`, commit 73a7120)

The pure-protein gate cuts genuine over-merges cleanly (genuine 4163→11) but drops **all 97 protein-missed
real TP** (non-coding / LOC\* paralogs, `in_ep=0`). No `core_recip` threshold separates real non-coding
paralogs (KEEP) from TE-bridged junk (CUT) because **both** are non-protein-homologous with overlapping
`core_recip`. **Hypothesis under test:** the soft-mask/TE status of the bridging exonic sequence
disambiguates them.

## The feature — whole-locus MAX exonic softmask (NOT the shared-region mask)

`GGO.fasta` is RepeatMasker soft-masked (lowercase = repeat, uppercase = unique; ~47% genome-wide). At the
**RNA level only exonized TEs matter** (spliced transcripts drop intronic repeats), so the mask is computed
on the **spliced exon-only** sequence of each de-novo locus (complement of the skeleton introns within
[start,end]) — verified intron-free: 0/94192 loci hit the whole-span fallback; exon/span ratios 0.14–0.48
confirm 52–86% of each genomic span (introns) is dropped before counting softmask.

For a gene-pair edge, `bridge_mask` = the **MAX** exonic softmask fraction over the two loci of the
argmax-`core` co-blocked DN–DN bridge edge (the same DN edge whose `core` the gate tests — no edge-mismatch).
It is a **whole-locus MAX proxy, not the aligned/shared-region mask**, for a concrete reason: divergent
paralogs (ENO2-ENO3, OGDH-OGDHL, GALNT14-GALNT16) **do not align at minimap2-asm20 at all**, so an
aligned-region mask is *undefined* for exactly the pairs the reframe exists to keep. The whole-locus proxy
covers all genes (incl. LOC\* with no annotated exon record) and needs no aligner → the gate is fully
deterministic. An **aligned-region shared-mask cross-check** (`minimap2 -cx asm20 -t1`) is computed as a
diagnostic to expose the proxy's MAX-inflation risk (below).

## The gate

> **`E_r edge valid  ⟺  in_ep  OR  (core_recip ≥ t  AND  bridge_mask < m)`**

`m = 1.0` reduces exactly to the existing `combined` gate (`in_ep OR core≥t`) — a built-in sanity anchor
(verified per `t`). Transitive-only pairs have no direct DN edge → `core=None` → `bridge_mask=None` → they
**fail the CORE branch regardless; the mask is never evaluated on them.**

## Part 1 — VERIFY: are genuine over-merges repeat-bridged? Honestly: only the protein-missed subset.

**Mask distribution by class** (pairs that carry a `bridge_mask`):

| class | n(with mask) | median mask | frac ≥ 0.10 | reading |
|---|---|---|---|---|
| **TP (all)** | 416 | **0.151** | 0.666 | real cDNA-truth paralogs — low |
| TP (protein-missed) | 90 | 0.172 | 0.611 | the recovery target — low-ish |
| **truthbar** (protein-homologous div. paralogs) | 3382 | **0.312** | 0.908 | **as repeat-rich as junk** |
| **genuine** (over-merges) | 1773 | **0.346** | 0.926 | junk — high |

The hypothesis is **POPULATION-FALSE for the real protein-homologous paralogs**: `truthbar` median mask
0.312 ≈ `genuine` 0.346 (frac ≥ 0.10: 0.908 vs 0.926). Real paralogs are as repeat-rich as the junk; they
survive **only** via the `in_ep` tautology, **not** via the mask. The mask separates **only** the
protein-missed TP subset (median 0.172) from genuine (0.346), and even there it is weak and **degrades at
the operating band**:

| population (in_ep=0, core present) | n genuine / n TP | **AUC(genuine > TP)** |
|---|---|---|
| whole re-admit region (core ≥ 0.13) | 1763 / 90 | **0.777** |
| core ≥ 0.30 | 521 / 64 | 0.705 |
| core ≥ 0.50 (recommended band) | 182 / 33 | **0.694** |
| core ≥ 0.70 (core-knee band) | 75 / 19 | **0.675** |

**Named "confirmation" — corrected for which branch actually decides each pair.** The showcase FPs are
**not** mask confirmations:

- **AMY2A–ZNF91, AMY2A–ZNF141, RPL14–ZNF669** have `core=None`/`mask=None` → cut by the **CORE branch**
  (transitive-only, no bridge). The mask is *never evaluated* on them. (The ZNF endpoints' high per-locus
  mask is incidental, not the gate feature.)
- **ENO2–ENO3, OGDH–OGDHL, GALNT14–GALNT16** are `in_ep=1` → kept by the **PROTEIN branch** regardless of
  mask (mask decorative).
- The **only named pairs the mask actually decides** are **SURF2–SURF4** (core 0.738, mask 0.068 → KEEP) and
  **CTSA–NEURL2** (core 0.808, mask 0.014 → KEEP).

**The honest mask-decided demonstration** (edges with `in_ep=0 AND core ≥ 0.5`, so the core branch admits
them and **only the mask can still cut**):

| what | count / rate | examples (core, mask) |
|---|---|---|
| **TE-junk the mask CUTS** (high core, high mask) | 138 / 182 genuine (75.8%) | LOC134757313–MTHFD2L (0.61, **0.95**); LOC109027212–LOC129528278 (1.00, **0.73**); LOC129528947–ZNF429 (0.71, **0.54**) |
| **real paralogs the mask KEEPS** (high core, low mask) | 21 / 33 TP kept | CENPP–IPPK (0.89, 0.04); SLC7A6–SLC7A6OS (0.80, 0.06); SURF2–SURF4 (0.74, 0.07); CTSA–NEURL2 (0.81, 0.01) |
| **junk the mask CANNOT catch** (high core, low mask genuine = residual) | 44 / 182 genuine survive | LOC129530205–LOC129530238 (**1.00, 0.035**) near-identical unique-seq LOC family |
| **real TP the mask WRONGLY cuts** (high core, high mask) | 12 / 33 TP lost | ASB6–NTMT1 (0.71, 0.34); LOC129523646–LOC129523647 (1.00, 0.33) |

Per-capita at the recommended band the mask cuts 75.8% of genuine vs 36.4% of real TP → **2.1× enrichment**
(a favorable but non-clean separator).

**Aligned-region cross-check (MAX-inflation).** 7 protein-missed TP have `|mask_a − mask_b| > 0.15` — the
whole-locus MAX may be driven by **one** repeat-rich partner, not a shared TE. The minimap2 shared-region
mask resolves direction:

| pair | whole-locus MAX | aligned shared-region | reading |
|---|---|---|---|
| SURF2–SURF4 | 0.068 | 0.074 | **faithful** (both low → correctly kept) |
| ENO2–ENO3, OGDH–OGDHL, GALNT14–GALNT16 | 0.076 / 0.052 / 0.046 | **DID NOT ALIGN (asm20)** | divergent → aligned mask *undefined* → whole-locus proxy **justified** |
| **FGF2–NUDT6** | **0.223** | **0.000** | MAX **over-inflates** (FGF2's repeat not in the shared bridge; NUDT6 is antisense) → this TP is **wrongly cut** at m=0.10 |
| DNAJC17–GCHFR | 0.181 | 0.039 | MAX over-inflates → wrongly cut |
| **ASB6–NTMT1** | 0.344 | **0.445** | shared bridge **genuinely repeat** → cut is correct |

So the whole-locus MAX is faithful for co-linear pairs and for the divergent paralogs (where it is the only
option), but over-inflates on ~2–3 of the 7 asymmetric TP; a shared-region or `min` mask would fix those at
the cost of being undefined for divergent paralogs. Documented as a known residual, not silently ignored.

## Part 2 — the combined oracle vs the reference operating points

| gate | rule | n | TPpm /97 | **ncTP /31** | truthbar | genuine | **blk prfam / edge-host** | recall .90/.95/.99 |
|---|---|---|---|---|---|---|---|---|
| none | keep all | 10755 | 97 | 31 | 6142 | 4163 | 15.89% / 48.73% | .951/.923/.972 |
| **pure protein** | `in_ep` | 6506 | **0** | **0** | 6142 | 11 | **0.25% / 2.01%** | .738/.731/.832 |
| no-mask core knee | `in_ep OR core≥0.70` | 6600 | 19 | 6 | 6142 | 86 | 2.75% / 11.21% | .803/.739/.881 |
| **RECOMMENDED** | **`in_ep OR (core≥0.50 AND mask<0.10)`** | **6571** | **21** | **7** | **6142** | **55** | **2.82% / 8.92%** | **.836/.746/.867** |
| strict-Pareto alt | `in_ep OR (core≥0.30 AND mask<0.05)` | — | 20 | 6 | 6142 | 50 | 1.69% / 6.27% | — |

**Recommended oracle:** `E_r edge = core_recip≥0.13 AND ( in_ep OR (core_recip≥0.50 AND bridge_mask<0.10) )`.
The mask threshold `m=0.10` is the smallest grid value that keeps SURF2–SURF4 (mask 0.068); `m=0.05` drops
it (the strict-Pareto alt) — an exposed, disclosed knob, not a hidden cherry-pick.

**Does it beat the no-mask sweep? Yes, off-curve, but modestly.** Exhaustive scan of the no-mask `combined`
curve finds **no** point with `TPpm ≥ 21` **and** `edge-host ≤ 8.92%` (t=0.6 gives TPpm 27 but edge-host
15.9%; t=0.8 gives edge-host 8.0% but TPpm 9; t=0.7 gives 19 / 11.2%). The recommended point therefore lies
**strictly off the no-mask curve**: vs the core knee it is **+2 protein-missed TP (+1 non-coding), −31
genuine over-merge edges, −2.29 pt edge-host block-rate** at ~equal prfam block-rate. Real, reproducible,
and not reachable by retuning `t` alone.

**But it is a precision knob, not a recovery mechanism.** At **fixed** `t=0.50`, turning the mask on
*reduces* recall (TPpm 33→21, genuine 193→55). The +2-TPpm net gain comes from being able to **lower `t`
0.70→0.50** — enabled by the mask's harder genuine-trim — not from the mask "recovering" any TP. And the
block-rate cut is **overwhelmingly the CORE branch, not the mask**: of the 4163 baseline genuine edges, the
core branch removes **3970**; the mask adds only **138** more. Do not credit the block-rate headline to the
mask.

## Honest residual (do not omit)

1. **TE-junk survives.** 44/182 genuine over-merges at the recommended band have low mask and cannot be cut
   — near-identical unique-sequence LOC-family homology (e.g. LOC129530205–LOC129530238, core 1.0 mask 0.035).
   The mask is blind to these by construction.
2. **Most real non-coding TP stay dropped.** 7/31 non-coding TP recovered (24 lost); 21/97 protein-missed TP
   overall (78% still cut). At m=0.10 the mask would cut **61%** of the protein-missed TP it aims to recover,
   because real LOC paralogs are themselves repeat-derived — the very premise "real non-coding paralog =
   unique/unmasked" fails for the majority.
3. **MAX-inflation.** ~2–3 of the 7 asymmetric protein-missed TP (FGF2–NUDT6, DNAJC17–GCHFR) are cut for the
   wrong reason — a repeat-rich partner, not a shared TE; the aligned cross-check shows the shared region is
   unique. A shared-region/`min` mask would fix these but is undefined for divergent paralogs.
4. **Intronic-TE caveat.** The mask counts only **spliced exon** sequence, so intronic/unexonized TEs are
   (correctly) invisible. A TE bridging two loci only through a shared intron would not be caught — but such
   a bridge also would not survive the exonic `core` homology, so it is not a live failure mode here.
5. **The truthbar 6142/6142 "win" is tautological** — it rides `in_ep`; the mask never touches it.

## Are TEs at the RNA level now handled?

**Partially, and honestly bounded.** The mask **is** a genuine, deterministic, exonic (intron-free), correctly
edge-attached repeat signal that (a) removes the flagrant repeat-rich over-merges (mask ≥ 0.5, e.g.
LOC…–ZNF429/675 zinc-finger bridges, MTHFD2L pseudogene) and (b) keeps the low-repeat non-coding paralogs
(SURF2–SURF4, CTSA–NEURL2), yielding an off-the-no-mask-curve operating point (**8.9% edge-host over-merge at
21/97 protein-missed TP**, vs the core knee's 11.2% at 19/97). It does **not** "handle TEs" in the strong
sense the open problem asked for: it cannot separate real protein-homologous paralogs from junk (they are
equally repeat-rich), it cannot cut the near-identical unique-sequence LOC residual, and it forfeits ~60% of
the protein-missed TP it targets. **The open problem is Pareto-shifted, not dissolved.**
