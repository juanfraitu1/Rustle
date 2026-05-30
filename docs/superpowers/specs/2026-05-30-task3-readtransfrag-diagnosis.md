# Task 3 (Phase 2a) diagnosis: read→transfrag accumulation is already ST-faithful

Date: 2026-05-30. Branch: `parity/isofrac-chain-dedup`.

## Decision: STOP — do NOT build `update_abundance_st`.

Step-1 diagnosis confirms the read→transfrag abundance construction is already
equivalent to StringTie's `update_abundance`. The contested-minority coverage
divergence is NOT in accumulation; it is downstream in the flow allocation / cov
formula at emit. Building `update_abundance_st` (a re-port of accumulation) would
be a no-op at best and is harmful in practice (see ST_READ_PATTERN experiment).

## Evidence

Captured `path_extracted` (cov, entry_abund, chain) + `pred_kill` from both tools:

```
RUSTLE_PREDCLUSTER_ST=1 RUSTLE_PARITY_LOG=/tmp/ru_pe.jsonl \
  RUSTLE_PARITY_FILTER_STEPS=path_extracted,pred_kill rustle GGO_19.bam -L -o /tmp/ru.gtf
STRINGTIE_PARITY_LOG=/tmp/st_pe.jsonl STRINGTIE_PARITY_FILTER_STEPS=path_extracted \
  tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st.gtf
python3 bench/path_extracted_diff.py /tmp/ru_pe.jsonl /tmp/st_pe.jsonl ../GGO_19.gtf
```

### 1. entry_abund (the read→transfrag accumulation output) matches or exceeds ST

Of 3492 chains extracted by both tools, **80.7% have byte-identical raw
entry_abund**. The mean |entry_abund within-locus ratio gap| (0.0553) is no
smaller than the |cov within-locus ratio gap| (0.0397) — accumulation is not the
bottleneck.

For the **217 contested minorities** (ref-TP chains Rustle kills via
isofrac/retained_intron but ST keeps, ST within-locus cov ratio > Rustle's):
- 61% (132/217) have EXACTLY identical raw entry_abund.
- Rustle entry_abund is HIGHER than ST in 52 cases, LOWER in only 33.
- Mean Rustle entry_abund = 13.21 vs ST = 7.48 — Rustle accumulates **~1.8x more**
  abundance onto these chains, not less.

The contested minorities are not coverage-depressed at read→transfrag. Their
relative cov is lower only because the **locus-dominant denominator** and the
**flow-distributed per-bp cov** differ at emit — i.e. selection/depletion, not
accumulation.

### Traced examples (both tools, same chain)

- `+ 16862079-16863832,...` (nex=10): Rustle entry_abund=13, ST entry_abund=13
  (identical), yet ru cov-ratio=0.318 vs st cov-ratio=0.687. Divergence is purely
  in the within-locus dominant / cov formula.
- `- 58854984-58857035,...,58863933-58865395` (nex=6): Rustle flow entry_abund=12,
  ST flow entry_abund=8 (Rustle higher). Sibling chain ...58864773 carries 86 (ru)
  / 83 (st). Abundances comparable; ratio differs by how flow concentrates.

### 2. Forcing ST-faithful accumulation makes parity WORSE

`RUSTLE_ST_READ_PATTERN=1` is the existing implementation of exactly what
`update_abundance_st` would port (one transfrag per long-read pattern, no
coverage-valley splitting). Under PREDCLUSTER_ST it REGRESSES:

| config | Sn | Pr | matching chains |
|---|---|---|---|
| PREDCLUSTER_ST (split accum, current) | 95.7 | 92.1 | 1734 |
| PREDCLUSTER_ST + ST_READ_PATTERN (ST-faithful accum) | 91.6 | 90.8 | 1652 (−82) |

The contested minorities are not recovered; the count of contested chains shifts
217→161 but mean gap widens (0.044→0.052) and 82 matching chains are lost.

### 3. Combined-flag milestone (3 deterministic runs)

| config | Sn | Pr | F1 | matching |
|---|---|---|---|---|
| default OFF (regression guard) | 96.6 | 91.6 | 94.03 | 1750 |
| PREDCLUSTER_ST alone | 95.7 | 92.1 | 93.87 | 1734 |
| **FLOW_ST + PREDCLUSTER_ST** | 93.3 | 90.8 | **92.03** | 1691 (×3, deterministic) |

Combined F1 92.03 < baseline 93.78 and far from the 95.27 ceiling. The cap is not
read→transfrag.

## Redirected next target

The real divergence is the **flow allocation + cov formula at emit** (per-bp cov
distribution and the within-locus dominant the isofrac/RI thresholds compare
against), i.e. `long_max_flow` flux→cov and the predcluster dominant selection —
NOT `update_abundance`. The PATHPAT flow-parity scope
(`2026-05-28-pathpat-flow-parity-scope.md`) and predcluster selection parity
(`2026-05-30-predcluster-selection-parity-design.md`) are the correct follow-ups.

`bench/path_extracted_diff.py` is the cov-parity gate for that work.
