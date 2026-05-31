# Rustle ↔ StringTie Long-Read Parity — Consolidated Findings

Thesis-record synthesis of the parity investigation on GGO_19 (chr19 NC_073243.2, `rustle -L` de novo
vs `GGO_19.gtf` = StringTie's own output). Authoritative running log: `STRINGTIE_PARITY_FINDINGS.md`.
This document is the consolidated narrative + the bottom line.

## 0. Bottom line

**Rustle starts from StringTie's graph substrate and the remaining gap is structural, not a missing
filter.** Across five graph layers, Rustle is the same standpoint as StringTie. The residual
false-positive gap (the only material headline difference) is **over-enumeration: wrong *combinations*
of *real* introns** produced by Rustle's flow/path-selection on a prediction population that diverges
from StringTie's — not imperfect junctions, exons, colors, or graph connectivity. Building-block-level
fixes therefore find no realizable prize; the one realized win came from suppressing an **RU-only extra
mechanism** StringTie lacks.

Default operating point (unchanged this session): **Intron-chain 96.2 / 91.7, Transcript 95.6 / 90.5**.

## 1. Per-level gffcompare picture (where the gap actually is)

| Level | Rustle Sn/Pr | StringTie Sn/Pr |
|---|---|---|
| Base | 97.4 / 94.0 | 99.1 / 97.7 |
| Exon | 98.2 / 96.2 | 99.4 / 98.5 |
| **Intron** | **99.3 / 98.0** | 99.6 / 98.7 |
| **Intron chain** | 96.2 / **91.7** | 96.7 / 96.3 |
| Transcript | 95.6 / 90.5 | 96.7 / 96.4 |

Introns are near-perfect (99.3 Sn, only 7 missed of 6396, the *same* 54 novel as ST). The whole
precision gap is the **intron-Pr 98.0 → intron-chain-Pr 91.7 drop**: correct introns assembled into
extra chains. That is over-enumeration, not bad building blocks.

## 2. Foundation: same standpoint as StringTie (5 layers, measured)

| Layer | Verdict | Evidence |
|---|---|---|
| Bundle envelopes | same | 3351/3430 exact (98.4%); 25-gap = benign within-locus reshaping; 10 wholesale-missing all 0-tx/0-ref |
| Bundlenodes | same | 3226/3351 byte-identical node sets (96.3%); mismatches inert (one inert intron-interior node-policy diff) |
| **Graph edges** | same (measured) | wired `graph_edge` emit both tools; 3037/3226 identical edge sets on identical-node bundles; **junction-spanning backbone never misses a real edge** (1 difference = a correct alt-donor omission, RU more precise) |
| Colors | same | junction-clean spans identical 1609/1610; 124/125 partition diffs trace to only-ST good junctions |
| Junctions | same (raw) | RU raw observes **6456/6456 ST-final + 6396/6396 reference introns**; only-ST raw = 1003 inert off-by-1 ghosts (0 in final/ref). good-set is a strict subset (only-RU=0; only-ST=10170 = the support/mm_negative floor) |

Also audited byte-faithful: the three-strand mechanism, `compatible_long` len[] containment geometry,
the color/CGroup union-find. The only cosmetic gaps (no `XS`-tag read; `get_min_start` tie-break —
**shipped fix**, provably inert) don't affect minimap2 long-read data.

## 3. The FP floor — what it is NOT (experiments that ruled things out)

Every targeted experiment was oracle-first (force ST's decision/values into RU, bound the prize before
building) and default-OFF gated. Results:

| Experiment | Outcome | Finding |
|---|---|---|
| Color + mm_negative segmentation | **ABORT** | segmentation prize = 0; RU already has a spanning bundle at every divergent locus; colors structurally negligible for F1 |
| Source/sink terminal wiring | **ABORT** | net −4; RU's "extra" terminal edges are **read-backed real-transcript endpoints**, re-derived from reads at 3 stages — not the FP lever |
| Retained-intron predicate (lever #1) | **SHELVED** | +25 oracle *ceiling* but unrealizable — 18/27 of ST's killers don't survive in RU's final set; the geometry predicate would have killed 133 real isoforms (oracle-first averted a −114 catastrophe) |
| Donor-snap (lever #2A) | **SHELVED** | prize 0 — the "weak donors" are real reference donors; FPs are wrong combinations of real introns |

## 4. The one realized win — lever #2B (chimeric-suffix-rescue)

`RUSTLE_CSR_FOLD` (default OFF, opt-in): suppresses Rustle's `chimeric_suffix_rescue` when the
5'-truncated suffix is contained in a kept full-length flow path (the case StringTie folds into the
flow parent rather than re-extracting). **Realized: RU-only FP 158→155 (−3), Pr +0.1 (IC + Tx), Sn flat
95.6, 0 TP cost.** On GGO_19 the rescue recovers 0 reference TPs — it is purely net-negative here. This
worked because it targets an **RU-only mechanism StringTie does not have** — a clean behavioral
difference, not a population-matching problem. (Recommended for default-flip; left opt-in pending
decision.)

## 5. The meta-insight (the thesis of the FP analysis)

- Rustle's false positives are **wrong combinations of real introns** (chimeric path enumeration), not
  spurious junctions/exons/donors. The building blocks match StringTie's; the *paths* through them do
  not, because Rustle's flow produces a different **prediction population**.
- Therefore **building-block-level fixes find zero realizable prize** (retained-intron, donor-snap both
  measured net 0 or negative — the introns/donors are all real).
- **Realizable parity gains live where Rustle does something *extra* StringTie does not** (suppress it,
  e.g. 2B's chimeric-suffix rescue) — NOT where Rustle must reproduce StringTie's population-dependent
  path selection. The latter is the irreducible **flow / coverage-depletion floor**: StringTie kills
  over-enumerated paths using predictions that exist in *its* pre-filter population and are absent from
  Rustle's.

## 6. Methodology — oracle-first, and why it paid off

Every experiment bounded the prize by forcing StringTie's exact decision into Rustle *before* building
a fix, with a hard abort gate and adversarial verification. This:
- **Averted a catastrophe:** the retained-intron geometry predicate looked plausible but the oracle
  measured it would kill 133 real isoforms (net −114) — never shipped.
- **Corrected mechanism mis-attribution four times** (retained-intron: not flags, not pairing, the
  predicate; lever-2: not transfrag-collapse but a junction snap + a path rescue; 2A donors are real).
  Lesson distilled: **don't *infer* the residual mechanism after an oracle rules one out — *trace* it;
  and an oracle prize is a *ceiling*, not proof a downstream fix can realize it (verify the killers
  exist in RU's population / the read folds into a surviving TP).**

## 7. Shipped + opt-in inventory (this session)

- **Shipped, default:** `get_min_start` tie-break fix (matches ST, provably inert on GGO_19).
- **Opt-in, default OFF (realized but not flipped):** `RUSTLE_CSR_FOLD` (+3 FP / 0 TP / +0.1 Pr).
- **Opt-in, default OFF (instrumentation/experiment, reusable):** `RUSTLE_ST_BADJUNC`, `RUSTLE_TERMINAL_ORACLE`,
  `RUSTLE_RI_LOCAL`, `RUSTLE_HIGHERR_SNAP`; the `graph_edge` / `junction_raw` / `good_junction` parity
  emits + chain emit in ST `pred_kill`/`pred_intron_low`; bench oracles `edge_diff.py`,
  `segmentation_prize.py`, `terminal_oracle_report.py`, `retained_intron_prize.py`,
  `retained_intron_geometry_oracle.py`, `retained_intron_chainsubset_oracle.py`, `junction_set_diff.py`,
  `csr_classify.py`, `donor_snap_prize.py`, `capture_parity.sh`.

## 8. Open / future

- **Default-flip `RUSTLE_CSR_FOLD`** (the realized Pr+0.1 / Sn-flat win) — pending decision.
- **Lever #3 (checktrf `longrec_fail`-drop vs `zero_flux`-rescue):** the last trace lever and the
  meta-pattern's best remaining bet — another RU-only *extra* rescue (the shape that made 2B
  realizable). Not yet scoped.
- **The flow / coverage-depletion floor:** the true source every lever points back to (RU's prediction
  population differs from ST's). The documented hardest, most-F1-negative-to-date target; matching it
  is what would close the residual ~5pp precision gap, but prior attempts to mimic the coverage
  substrate are F1-negative.
