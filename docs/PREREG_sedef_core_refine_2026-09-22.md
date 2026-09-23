# Pre-registration — does consuming SEDEF improve the family definition in the de novo and guided modes?

**Written 2026-09-22, §6x2, before any arm is run.** User: *"would also consuming the SEDEF then improve
family definition in all modes?"* → *"yes lets actually run it."*

## Nothing is being built — the flag already exists

`mcl_families --core-refine --sedef <bed>` is §6eh's duplicon-first refinement: within each cluster a
member's CORE is the part linked by SEDEF pairs to ≥ half the other members; SD-poor clusters are
UNTOUCHED; otherwise core ≥ span/2 ⟹ kept, else core ≥ half the cluster's median core ⟹ kept and TRIMMED
to the core hull, else DROPPED. §6x1/r994 established that **de novo and guided consume no SEDEF today**
(`core_refine=false sedef=<unset>`). This turns it on and scores it.

## ⚠ What the operator can and cannot do — structural, not empirical

`CoreStatus` is `Untouched | KeptFull | KeptTrimmed | Dropped` (`annotation_families.rs:850-854`) and
`refined.clusters.tsv` emits only non-dropped members. **The refinement can never ADD a member to a
cluster, never merge two clusters, and never create an edge.** It is a PRECISION-ONLY operator by
construction, so **sensitivity can only fall or stay equal** and the entire question is whether the
precision gain outweighs it on F. Stated now so a sensitivity drop is not later reported as a discovery.

This also settles the recall half of the user's question in advance: consuming SEDEF **cannot** close
§6u4's edge-construction gap (48.8% of FNs) or move §6o8's 0.052 pairwise-recall ceiling, because both are
edge-ADMISSION losses and this operator runs post-MCL on edges that already exist.

## The decisive control — C1, and why it is the whole experiment

| arm | invocation |
|---|---|
| **A0** shipped | existing `clusters.tsv`, no refinement |
| **A1** SEDEF | `--core-refine --sedef /mnt/linuxdisk/home/juanfraitu/winloci_data/HSA_sedef_pairs.bed` |
| **C1** control | `--core-refine --core-from-paf` — derives the same SD-like pairs **from the input PAF itself** |

⚠⚠**If C1 reproduces A1, SEDEF contributed nothing an external caller was needed for**, and register 650
is confirmed on a second population: *"an SD-supported PAF cross-edge IS a SEDEF pair; 41/42 SD-rich
clusters are one SEDEF component alone."* A1 must beat C1 to be a SEDEF result at all rather than a
core-refinement result. Without C1 this experiment cannot distinguish the two and would be worthless.

## Substrates — and the truth may NOT be Soto

- **Development / visible: human chr16**, de novo (`dn16_fam3`) and guided (`chr16_guided`).
- **Held out: human chr2, chr8, chr10** (`/mnt/linuxdisk/tmp/heldout/chr*_fam.*`, zero exposure per §6s8).
  Reported with no re-tuning. ⚠Never pooled with gorilla.

⚠⚠**Soto is NOT the primary truth here and a Soto-scored gain would not count.** Soto's catalog is a
98%-identity SEGMENTAL-DUPLICATION-derived set (register 741), and register 1085 already records that 3 of
4 legs of our "Soto replication" consume Soto's own SEDEF track. Feeding SEDEF into the definition and
scoring against Soto measures **agreement with the comparator's own substrate**. Soto numbers are reported
as context and explicitly labelled circular.

**Primary truth = the protein-family referee** (§6u5/`bench/protein_families.py`): longest CDS per gene,
translated in-house, clustered independently of both our gate and Soto's SD track.

## ⚠ A measurement hazard found while reading the code, fixed before running

`bench/soto_vs_us_referee.py:131` joins a cluster member to a gene by **exact coordinates**
(`names.get(f'{chrom}:{start}-{end}')`). A KeptTrimmed member is emitted with its **core-hull** coordinates,
which match no gene's span, so it would silently vanish from `ours` — **making trimming indistinguishable
from dropping and understating every refined arm.** Scoring therefore uses `bench/mode_family_score.py`,
whose locus→gene resolver is max-OVERLAP, with the referee's labels dumped to Soto's two-column format.

## Metrics and the bar — committed now

Sensitivity · precision · one-to-one bipartite F · members dropped · members trimmed, per arm per chrom.
Judged on the **held-out chromosomes**, against the **protein referee**:

| outcome | verdict |
|---|---|
| A1 beats A0 on **F**, beats **C1**, and precision rises without F falling | ⭐ **ADOPT-WORTHY** — SEDEF earns its place |
| A1 beats A0 on precision only, F flat or down, **or A1 ≈ C1** | ⚠ **PARTIAL / NOT SEDEF** — a precision knob, not a definition change |
| A1 below A0 on F on the held-out set, or it drops real members | ⛔ **NO** |

**Predicted, before looking — ⚠ PARTIAL, and specifically `A1 ≈ C1`.** Reasons on record: register 650
(SD-supported PAF edge *is* a SEDEF pair); register 1060 (SEDEF **misses** real pairs — APOBEC3D/F is in
E_c at 88.4%, under the Bailey cutoff — so `E_c ⊄ E_a`); register 669 (SEDEF corroborates the DUPLICATION,
not the gene-ness: MCL32/MCL24 are fully SD-supported and 1.00/0.97 curated repeat); register 677 (the
DROP clause already removed real members of patchy-SD families — MCL5 −2, MCL25 −3, MCL13 −3). Register
318's ⭐ FP 0/14 admission result is the one reason to expect a precision gain at all.

I will not change the arms, the control, the substrates, the truth or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⛔ **NO. SEDEF changes nothing it should, and removes real genes where it changes anything.**

All ten arms reproduced their shipped `clusters.tsv` **byte-identically** (the built-in control the
`--core-refine` help promises), so the flag reconstruction is validated before any score.

## Held out — chr2 / chr8 / chr10, zero exposure

| truth | arm | chr2 | chr8 | chr10 |
|---|---|---|---|---|
| **protein referee** (primary) | A0 shipped | 0.236 | 0.436 | 0.197 |
| | **A1 SEDEF** | **0.236** | **0.436** | **0.197** |
| | C1 from-paf | 0.236 | 0.436 | 0.203 |
| Soto (⚠circular, context only) | A0 shipped | 0.778 | 0.554 | 0.750 |
| | **A1 SEDEF** | **0.778** | **0.554** | **0.750** |
| | C1 from-paf | 0.750 | 0.554 | 0.761 |

⛔**A1 is EXACTLY A0 on every held-out chromosome, on BOTH truths — to every digit.** Development moved
only on de novo chr16 (precision 0.949 → 0.974, **F +0.001**); guided chr16 was unchanged.

## ⛔ It is not inert, it is invisible — and what it removes is real

The refinement does real work (chr2/A1: 12 trimmed, 11 dropped). The score cannot see it. What it drops,
on the **held-out** chromosomes:

| arm | dropped locus | span | SD depth | SEDEF core | fraction of the gene |
|---|---|---|---|---|---|
| chr10/A1 | **AGAP11** | 39,816 bp | **9** | 9,262 bp | **1.00** |
| chr10/A1 | LINC00863 | 14,900 bp | 4 | 2,813 bp | 1.00 |
| chr10/A1 | CUBNP2 | 15,317 bp | 1 | 0 | 1.00 |
| chr2/A1 | FAM95A | 5,859 bp | 0 | 0 | 1.00 |

**`AGAP11` is a whole gene at SD depth 9 carrying 9.3 kb of SEDEF core, and it is DROPPED** — because
9,262 < span/2, and its core is under half its family's median. That is **register 677 reproducing on a
held-out chromosome**: the drop clause removes real members of families whose core coverage is uneven.
⭐⭐**And F does not move when it happens** — which is why the clause survived: the harm is real and the
metric is blind to it. A rule can only be as safe as the metric that is allowed to reject it.

⚠**Correction to my own first read of this table.** I initially reported "A1 drops CASP10, a real gene."
It does not. The chr2 drops are 517–5,859 bp fragments that merely *overlap* CASP10/IGK/LOC442028, and a
max-overlap resolver named them after the gene they sit inside. Printing the spans beside the names is
what caught it — **a dropped-member audit must show the locus span against the gene length, never the
gene name alone.**

## ⭐⭐ The control settles the actual question: SEDEF is a LOSSY COPY of the PAF we already compute

Over 1,013 members where both arms find a core:

- **median |A1 − C1| core = 0 bp**; **51.7% byte-identical**; 69.9% within 100 bp; 84.5% within 1 kb.
- `AGAP11`'s core is **9,262 bp under SEDEF and 9,257 bp under the PAF** — a 5 bp difference between two
  supposedly independent sources.
- **29.1% of members (425/1,461): SEDEF sees NO core where the PAF does.**

⭐**Register 650 confirmed on a second, held-out population** — *"an SD-supported PAF cross-edge IS a
SEDEF pair."* And register 1060 quantified: SEDEF is not a superset of our homology, it is a **subset with
a 29% blind spot**. So consuming SEDEF cannot add information to the definition; it can only re-derive
what the all-vs-all PAF already contains, minus what SEDEF's own thresholds discard.

## Verdict against the pre-registered bar

⛔ **NO** on the row "A1 below A0 on held-out F, **or it drops real members**" — it drops `AGAP11`,
`LINC00863`, `CUBNP2`, `FAM95A`. And ⚠ **NOT SEDEF** on the row "A1 ≈ C1": the two arms agree to a median
of 0 bp, so any gain would have been a *core-refinement* result, not a SEDEF result.

**Do not enable `--sedef` or `--core-refine` in the de novo or guided modes.** §6x1/r994's conclusion
stands unchanged and is now positively tested, not merely audited: SEDEF is not load-bearing, and adding
it is not an improvement — it is a lossy re-derivation that costs real members.

## What I got right and wrong in the prediction

Right: **`A1 ≈ C1`** (predicted explicitly, measured at median 0 bp), and no recall movement anywhere —
guaranteed a priori, since `CoreStatus` can only keep, trim or drop. Wrong: I predicted a **precision gain**
worth calling ⚠ PARTIAL. There is none on held-out — precision is bit-identical on all three chromosomes —
and the only precision movement anywhere (de novo chr16, +0.025) buys **+0.001 F**. Register 318's
⭐ FP 0/14 admission result does not transfer to this operator: 318 tested SD membership as an *admission*
certificate on edges, this tests SD core fraction as a *member-removal* rule post-MCL. **They are opposite
operators and the 318 result is not evidence for this one** — an error I made in the prediction and am
recording rather than quietly dropping.

> **Scorer (2026-09-22 port):** every `bench/mode_family_score.py` number above is reproduced byte-for-byte by `target/release/family_score` (same flags; 732/732 parity, scipy's assignment tie-breaking included — r1045/r1046). The Python was retired in §6z2.
