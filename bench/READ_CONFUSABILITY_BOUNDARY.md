# Read-confusability boundary vs the divergence floor

**Question (thesis scope).** The copy-assignment method targets **read-confusable** multi-copy
families: copies so similar the aligner maps a read equally well to ≥2 loci and gives up (MAPQ 0 /
tiny best-vs-2nd alignment-score margin). **Divergent** copies are aligner-resolvable (reads carry
enough mismatches to place uniquely, MAPQ>0) and never needed the method. If that is true there is a
**copy-identity boundary**: below it the aligner resolves copies itself; above it resolvability
collapses. The divergence **floor** wired into the family definition (`family_rna_refine`,
reciprocal whole-transcript identity `recip_id_best`, default 0.80) would be **principled** in
Canzar's sense (no arbitrary threshold) *iff* the P/R-optimal floor **coincides** with this
read-confusability boundary — i.e. the method models exactly the families the aligner cannot resolve.

**Verdict up front.** There **is** a read-confusability boundary (≈**0.98–0.99** copy identity,
confusable mass concentrated there), and the thesis direction holds (divergent copies are
aligner-resolved). But it **does not coincide** with the divergence floor: the floor is an
**O1 whole-transcript purity** cut, the boundary is an **O2 per-read shared-segment** property.
They are **two separate objects**. The "we only model the challenging cases the aligner cannot
resolve" claim is supported **directionally, not as a unified threshold**.

Deliverables (absolute paths):
- Script: `/mnt/c/Users/jfris/Desktop/Rustle/bench/read_confusability_boundary.py`
- Per-family TSV (393 fam): `/mnt/c/Users/jfris/Desktop/Rustle/bench/read_confusability_boundary.tsv`
- JSON: `/mnt/c/Users/jfris/Desktop/Rustle/bench/read_confusability_boundary.json`
- Figure: `/mnt/c/Users/jfris/Desktop/Rustle/bench/fig_read_confusability_boundary.png`

---

## Substrate and metrics

**Families.** 393 multi-copy families (`n_loci≥2`) that carry a copy-vs-copy identity signal, drawn
from the **floor-OFF (573-family) catalog** — the only catalog that still contains the families a
floor would cut, so the coincidence test is meaningful. The on-disk `family_rna_refine.tsv` is a
**floored default that a parallel sweep rewrites live** (it flipped between 573/floor-OFF and
307/floor-ON 0.80 during this analysis), so the substrate is a **race-safe snapshot** rebuilt
deterministically via `family_rna_refine.build_catalog(divergence_floor=0.0)` (recovers the
divergent-inclusive catalog md5 `dca64cbd` byte-identically). All 393 family_ids reproduce the
prior frozen substrate with identical `dominant_gene`; the load-bearing per-family counts
(BAM-MAPQ0, depth, AS-margin) are byte-identical to the prior run.

**Copy identity (x-axis), per family** over covered member pairs:
- `mean_matched_id` — whole-read proxy: mean identity of colinearly-matched exons `[gate]`
- `best_exon_id` — per-exon proxy: identity of the single most-similar exon `[gate]`
- `recip_best_mean` — **`recip_id_best` = matches_best/max(len)**, aggregated over member gene-pairs.
  **This is the ACTUAL divergence-floor metric** (`family_rna_refine.load_recip_id`; the floor cuts a
  cross-gene edge iff `recip_id_best < floor`). `[ri_sharedlen_recip_id.tsv]`
- `core_recip` — a gate-derived reciprocal identity, retained only as a **weaker proxy**
  (corr ≈ 0.49 to `recip_id_best`; the earlier draft mislabeled it "the floor metric"). `[gate]`

**Read resolvability (y-axis)**, three independent metrics that cross-check:
1. **BAM MAPQ0** from `GGO_mm.bam` (aligned `-N50 -p0.1 --secondary=yes`): fraction of primary
   reads at the family's loci with MAPQ==0 (multi-mapping / confusable). Standard signal, but forced
   secondary can depress MAPQ → controlled below.
2. **Re-align default MAPQ0**: a seeded 40-read/family sample re-aligned once with **default
   `minimap2 splice:hq`** (no forced secondary) — the un-depressed control.
3. **AS-margin < 10**: best AS − 2nd-best AS per read; small margin = confusable. Setting-robust and
   the direct connection to **Eichler's AS≥10 read-assignment rule** (TBC1D3: keep a read for a copy
   iff its best AS beats the runner-up by ≥10, else discard).

Deterministic: `PYTHONHASHSEED=0`, `seed=0`, `minimap2 -t4 --seed 42`. The coincidence stage runs
`--from-tsv` over the frozen per-family TSV (byte-reproducible; JSON md5 `580affa5…`, identical
across reruns modulo the figure-path key).

---

## 1. Resolvability-vs-copy-identity curve (read-weighted, whole-read axis)

| mean_matched_id | nFam | reads | BAM MAPQ0 | re-aln MAPQ0 | AS-margin<10 |
|---|---|---|---|---|---|
| <0.50       | 88 | 62.3k | 0.002 | 0.001 | **0.068** |
| 0.50–0.85   | 53 | 41.9k | 0.001 | 0.003 | 0.007 |
| 0.85–0.90   | 36 | 54.4k | 0.014 | 0.022 | 0.037 |
| 0.90–0.95   | 38 | 17.7k | 0.013 | 0.057 | 0.067 |
| 0.95–0.98   | 63 | 24.3k | 0.002 | 0.002 | 0.026 |
| 0.98–0.99   | 34 | 15.3k | 0.012 | 0.023 | 0.100 |
| **0.99–1.00** | 81 | 42.1k | **0.035** | **0.110** | **0.186** |

Confusion rises toward near-identity on all three metrics; the top bin (≥0.99) is the most
confusable. (The `<0.50` bin's elevated AS-margin is the segment-local retrocopy class — see §3.)

## 2. The located boundary (whole-read + per-exon + the actual floor axis)

Sharpest upward jump in confusable fraction, per copy-identity axis × read metric:

| copy-identity axis | BAM-MAPQ0 boundary | AS-margin boundary |
|---|---|---|
| `mean_matched_id` (whole-read) | 0.99 (+0.024) | 0.99 (+0.086) |
| `best_exon_id` (per-exon max)  | 0.99 (+0.025) | 0.99 (+0.122) |
| **`recip_id_best` (ACTUAL floor axis)** | **0.98 (+0.251)** | **0.98 (+0.294)** |
| `core_recip` (gate proxy)      | 0.90 (+0.235) | 0.99 (+0.430) |

**Boundary ≈ 0.98–0.99 copy identity.** 6 of 8 axis×metric estimates land at 0.99; the ACTUAL
floor axis (`recip_id_best`) puts it at **0.98** on both metrics; the lone outlier is the
weak `core_recip:BAM` proxy at 0.90. So the boundary is robust to axis choice within one bin.

**Mass concentration** (fraction of all confusable reads from families ≥ identity `t`):

| t | BAM-MAPQ0 reads | AS-margin reads |
|---|---|---|
| ≥0.85 | **94.8%** (2704/2853) | 78.4% (927/1182) |
| ≥0.95 | 60.0% | 65.7% |
| ≥0.98 | 58.5% | 60.4% |
| ≥0.99 | 52.2% | 49.1% |

Confusable mass concentrates at near-identity: 94.8% of BAM-MAPQ0 reads and 78% of AS-margin reads
are from families ≥0.85; ~half from ≥0.99. **The thesis direction holds** — divergent copies are
aligner-resolved.

**Caveat: it is a shallow ramp, not a cliff.** Even the ≥0.99 bin is only **3.5% BAM-MAPQ0 /
18.6% AS-margin** confusable, because full-length IsoSeq reads span the whole transcript and almost
always catch a distinguishing exon. The truly-unresolvable mass is the **per-exon/PSV-local +
collapsed** subset (consistent with the K=0 floor findings), not a whole-read identity cliff.

## 3. Clean identity line, or exon-dependent? → Decisively per-read / segment-local

- **No family-level identity metric predicts confusability.** corr(confusable, `mean_matched_id`)
  = **+0.14**; corr(confusable, `best_exon_id`) = **+0.12**; corr(confusable, `recip_id_best`)
  = **+0.30**. All weak — confusability is a **per-read** property, not a family property.
- **Segment-local proof (smoking gun).** 14 families with `mean_matched_id==0` **and**
  `best_exon_id==0` (no colinear exon match at all) are still highly confusable — ribosomal/hnRNP
  retrocopy families sharing one near-identical segment:

  | family | recip_id_best | core_recip | AS-margin confusable |
  |---|---|---|---|
  | RPL7     | **0.147** | 0.294 | 35/40 (87.5%) |
  | HNRNPA1  | **0.399** | 0.512 | 31/40 (77.5%) |
  | RPL9     | **0.402** | 0.606 | 29/40 (72.5%) |
  | RPL10L   | **0.244** | 0.773 | 28/40 (70.0%) |

  Confusability tracks the identity of the **shared segment a read covers**, decoupled from
  whole-transcript identity. All four sit at *low* whole-transcript `recip_id_best` (0.15–0.40).
- **Divergent-overall families still carry confusable reads.** The 141 families with
  `mean_matched_id<0.85` have only 0.14% BAM-MAPQ0 residual but a real AS-margin residual (the
  segment-local class), confirming exon-dependence rather than a clean line.

## 4. Confound checks (all pass)

- **Forced-secondary MAPQ depression = NEGLIGIBLE.** mean BAM-MAPQ0 conf **0.0365** vs re-aligned
  **default**-MAPQ0 conf **0.0361** (corr 0.994). The `-N50 -p0.1` does **not** inflate MAPQ0; the
  standard BAM signal is trustworthy. (At near-identity the default aligner produces *more* MAPQ0,
  not fewer.)
- **AS-margin is ~2.2× more sensitive** than MAPQ0 (mean 0.081 vs 0.037; BAM-vs-margin corr 0.77) —
  it catches near-ties minimap2 still assigns MAPQ>0. Same ordering / boundary. This is exactly the
  Eichler AS-margin regime.
- **Depth:** corr(conf, depth) = −0.083 (deeper families are *not* more confusable);
  partial(conf, identity | depth) = **+0.138** (identity survives controlling for depth).
  **Size:** corr = +0.019 (~0). Not a depth/size artifact.

## 5. The floor (O1) and the coincidence test

**The floor is a monotone precision/recall purity cut, not a resolvability boundary.** The
`divergence_floor.tsv` sweep (evaluated on the shipped family-level P/R) is monotone: E_p purity
0.892 → 0.979 (rising) and R_oracle 0.877 → 0.491 (falling) across floor 0 → 0.90. Both F1 series
peak at **floor = off** (F1(E_p,R)=0.884, F1(P_task,R)=0.886) — **no interior optimum**. The shipped
**0.80** is the **domain-share purity knee**: the smallest floor excluding all 4/4 named
domain-shares (0.70/0.75 exclude only 3/4, 0.80 excludes 4/4). That is an **O1-precision** rationale.

**Coincidence — two tests, both fail, on the ACTUAL floor metric `recip_id_best`:**

- **(A) Axis test.** Boundary on the floor's own axis (`recip_id_best`) = **0.98** vs operative
  floor **0.80** → gap **0.18** (> one identity bin). Across all axes the boundary–floor gap spans
  **0.10–0.19**. Not coincident.
- **(B) Mass test (decisive, axis-independent).** A `recip_id_best` floor at **0.80 discards
  67.4% of BAM-MAPQ0 and 55.5% of AS-margin confusable reads** — the opposite of a
  boundary-coincident floor (which would cut ~0):

  | floor (recip_id_best) | BAM-MAPQ0 reads cut | AS-margin reads cut |
  |---|---|---|
  | 0.70 | 52.9% | 50.0% |
  | **0.80 (shipped)** | **67.4%** | **55.5%** |
  | 0.85 | 74.6% | 65.5% |
  | 0.90 | 78.8% | 73.7% |

  Cross-check on the weaker `core_recip` proxy: 53.9% / 59.6% at 0.80. The **ACTUAL** metric cuts
  *more* BAM-MAPQ0 (67% vs 54%) because `recip_id_best` is lower on average — i.e. the earlier
  core_recip-proxy report, if anything, **understated** the anti-coincidence. The entire
  segment-local class (RPL7/HNRNPA1/RPL9/RPL10L, `recip_id_best` 0.15–0.40) is cut under both metrics.

Confusability's correlation with the floor axis (`recip_id_best` +0.30) is higher than with the
matched-exon axis (`mean_matched_id` +0.08), but **neither** family-level identity metric cleanly
predicts it — the driver is a per-read shared-segment identity no whole-transcript floor can gate.

---

## Boundary + floor + coincidence verdict

- **O2 read-confusability boundary** ≈ **0.98–0.99 copy identity** (0.98 on the actual `recip_id_best`
  axis, 0.99 whole-read / per-exon). A shallow ramp, not a cliff; confusable mass concentrated at
  near-identity (94.8% of BAM-MAPQ0 reads ≥0.85). Per-read / segment-local, not a family-level line.
- **O1 divergence floor** = **0.80 on `recip_id_best`** (whole-transcript reciprocal identity), a
  monotone-P/R domain-share **purity** knee (O1-precision), **not** a resolvability boundary.
- **Coincide?** **NO.** Axis gap 0.10–0.19; the 0.80 floor **cuts 67% BAM / 55% AS-margin of the
  read-confusable reads** and removes the entire segment-local confusable class. O1 ⟂ O2, exactly as
  the memory predicts (homology/precision axis ⟂ read-resolvability axis).

## Principled-scope statement

The divergence floor is **not** an arbitrary 85% and it is **not** the read-confusability boundary
either — it is a **principled O1 homology/precision floor** (`recip_id_best ≥ 0.80` = the smallest
whole-transcript reciprocal-identity cut that purges all named cross-gene domain-shares). The
read-confusability boundary is a **separate, principled O2 object** (≈0.98–0.99 copy identity, where
the aligner stops resolving copies) that principledly scopes **where the copy-assignment / PSV / SUN
machinery is needed** — the MAPQ-0 / AS-margin<10 substrate, the same regime as Eichler's AS≥10 rule.

The user's scope claim — *"we only model the challenging cases the aligner cannot resolve"* — is
supported **directionally, not as a unified threshold**, with two honest caveats:
1. The confusable mass is **small even at high identity** because full-length reads self-resolve via
   distinguishing exons; the real target is the **segment-local / collapsed** subset (the K=0 floor).
2. Some of the **most read-confusable families sit at low whole-transcript identity** (RPL7-type
   retrocopies) and are **removed by the 0.80 floor** — so the floor and the confusability boundary
   must be reported as **two objects**, not conflated: an O1-precision floor
   (0.80, whole-transcript `recip_id_best`) and an O2-read-confusability boundary
   (≈0.98–0.99, per-read shared-segment identity).

*Fix applied vs the prior draft:* the coincidence/mass test now uses the **actual floor metric
`recip_id_best`** (loaded from `ri_sharedlen_recip_id.tsv`, joined per family over member gene-pairs)
instead of the mislabeled `core_recip` proxy; `core_recip` is retained only as a documented
cross-check. Verdict unchanged and, on the BAM axis, strengthened.
