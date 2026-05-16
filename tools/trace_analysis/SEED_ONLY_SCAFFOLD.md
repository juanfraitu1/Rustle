# seed_only scaffold — refs whose exact chain IS seeded but flow extracts a different path

## Purpose

A `seed_only` missed ref has its EXACT intron chain registered as a
`transfrag_seed` (flow extraction WAS handed the right seed) yet no
`path_extracted` carries that chain — flow produced a different path from the
same seed. This is distinct from (and more tractable than) the `no_seed`
COMBINATORIAL class: here the seed exists, so the divergence is purely in path
extraction / extension, not transfrag construction or flow combinatorics.

chr19 GGO at the 1746/1948 baseline: **14 seed_only refs**.

## Component

`seed_only_analysis.py` — for each seed_only ref, aligns: ref chain vs the
seed (abund) vs checktrf_enter/result vs the best-overlapping `path_extracted`
chain (what flow built instead), reporting longest exact consecutive run and
the first divergence, then classifies.

```bash
RUSTLE_PARITY_LOG=/tmp/ns.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  RUSTLE_FLOW_RESIDUAL_SE=1 ./Rustle/target/dev-opt/rustle -L GGO_19.bam -o /tmp/ns.gtf
gffcompare -r GGO_19.gtf /tmp/ns.gtf -o /tmp/ns_cmp
python3 Rustle/tools/trace_analysis/seed_only_analysis.py \
  --parity /tmp/ns.jsonl --ref GGO_19.gtf --refmap /tmp/ns_cmp.ns.gtf.refmap
# --tid STRG.93.3 to deep-dive one
```

## Taxonomy (chr19 GGO, 1746/1948 baseline)

| class | n | meaning | tractability |
|---|---|---|---|
| `TERMINAL_5P` | 6 | flow path missing the ref's 5′-most 1–2 junction(s); a (often higher-cov) sibling won the dominant path. STRG.200.7, .22.2, .22.3, .555.5, .555.6, .566.9 | **most tractable** — 5′-boundary / first-junction / back-extension selection. Same family as Bug B alt-TSS work ([[project_bugB_terminal_alt_tss_2026_05_13]]) |
| `TERMINAL` | 4 | extracted = ref ± contiguous terminal jct. STRG.93.3 (ref ⊂ extracted nx=4 vs 3 — back-ext absorbs upstream exon, matches [[project_trace_strg93_2026_05_13]]), STRG.125.7 (checktrf_rescue superset), STRG.200.11 (lowcov superset), STRG.445.6 (extracted ⊂ ref, 3′ truncated) | tractable — terminal clip/extend gate |
| `READTHR_SKIP` | 3 | ref chain reached checktrf, killed by readthr gate. STRG.125.4, .125.6, STRG.409.4 | known F1-negative to loosen globally; per-case only |
| `NObest` | 1 | seed consumed, nothing emitted at that span. STRG.136.1 (only j-class via Bug B) | hard |

## Key observation

10 of 14 are **terminal divergence** (flow extracts ref ± terminal junctions),
NOT internal reroute. The seed carries the right chain; path extension either
(a) skips the ref's 5′-most junction in favour of a higher-cov sibling start
(TERMINAL_5P), or (b) over-/under-extends a terminal exon (TERMINAL). This is a
boundary-selection problem in `back_to_source` / `fwd_to_sink` / the seed
start-node choice — the same machinery as the landed Bug B + clip-gate
alt-TSS fixes, which used per-node read-start pileup as the gate.

## Coordinate convention

Intron-string `donor_exon_end+1 - acceptor_exon_start-1` (same as
no_seed_analysis.py / parity emits). `transfrag_seed` carries no tf index, so
seed/checktrf linkage is by exact-chain match; "what flow built instead" is the
best span-overlapping `path_extracted` (longest exact consecutive run).

## TERMINAL_5P deep-dive — DONE (2026-05-15): same long_max_flow wall

Deep-dived STRG.22.2/.3 + STRG.566.9. Result: **NOT the Bug-B back-extension
class; it is the same `long_max_flow` / path-completion node-selection wall as
the COMBINATORIAL no_seed class** (the 97% systemic fwd-extension divergence
from [[project_combinatorial_no_seed_2026_05_15]] option 2).

Evidence (STRG.22.2, representative):
- Its exact 8-jct chain IS a registered `transfrag_seed` (abund=5.0) AND
  appears in `transfrag_define` (×2) — keeptrf does NOT lose it.
- Yet its distinctive first junction `16864075-16864209` appears in
  `path_extracted` ZERO times at every stage (bundle_after_predcluster
  onward). checktrf=0 too.
- Read evidence: ~75 reads start at 16863900-16863999; first-intron
  distribution = `16864401-16864522` ×57 (dominant merged-first-exon,
  STRG.22.3-shaped), `16864075-16864209` ×11 + `16863941-16864209` ×9
  (STRG.22.2's distinctive split-first-exon = minority ~20 reads).
- A dominant merged-first-exon seed (abund=26) is processed first, depletes
  shared downstream flow; when STRG.22.2's registered abund-5 seed is
  flow-extracted, `long_max_flow` along its pattern drops the low-flow
  distinctive first junction and rejoins the dominant structure.
- STRG.566.9: additionally blocked by a REJECTED junction
  (`111149171-111150205`, `mm_negative`) — the known preemptive-mark bug the
  parity README says regresses F1 catastrophically if fixed in isolation.

Conclusion: TERMINAL_5P = registered-seed-but-flow-drops-its-own-junction
under depletion. Same root cause as COMBINATORIAL (long_max_flow node
selection), reached from the seed_only side. The Bug-B per-node read-start
pileup gate does NOT apply (loss is in flow path-selection, not terminal
boundary extension). Not a clean win. Do not pursue without solving the
underlying long_max_flow depletion/path-selection divergence (multi-session,
negative precedent).

## RELATED: killed-after-extract attribution (2026-05-15)

New scaffold `killed_after_extract_analysis.py` attributes the 42 missed
multi-exon refs that DID reach `path_extracted` (flow assembled them) to their
killing filter:

| killer | n |
|---|---|
| AFTER_isofrac | 17 |
| AFTER_pairwise_overlap_filter | 11 |
| AFTER_readthr_gate | 10 |
| LOST_after_last_stage (post-FINAL globals) | 4 |

**Key quantified finding (isofrac 17):** the isofrac filter is a faithful ST
port; it kills on too-low input cov. STRG.15.1 — identical exact 31-jct chain,
identical longcov=2.0 on BOTH sides — post-flow `cov`: **rustle 2.914 vs
StringTie 4.5223** (ST survives→FINAL→in ref; rustle <isofrac→killed). The
PATH is identical; only the cov-ESTIMATION diverges (~1.55×). This is a
tractable cov-formula divergence, NOT the path-selection wall. See
[[project_cov_estimation_divergence_2026_05_15]].

## Remaining (not deep-dived)

- **TERMINAL (4)** — STRG.93.3 (back-ext absorbs upstream exon, prior trace)
  / .445.6 (3′ trunc): plausibly per-node pileup-gateable but lower priority
  given TERMINAL_5P outcome.
- **READTHR_SKIP (3)** — per-case only; global loosening F1-negative
  ([[project_parity_ceiling_2026_05_08]]).
