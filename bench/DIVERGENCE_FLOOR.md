# Absolute divergence floor — scoping the RNA family definition to near-identical paralogs

**What this is.** An **absolute reciprocal-identity divergence floor** on the cross-gene edge
rule of the RNA-only family definition (`bench/family_rna_refine.py`), **DEFAULT-ON** with an
opt-out. It scopes the catalog to gene families whose members are `>= FLOOR` **reciprocal
whole-transcript identity**, at HIGH precision (single-exon domain-shares and divergent paralogs
excluded). This is Soto's SD98 move (scope to the clean near-identical regime) — but with an
explicit opt-out that recovers the ambitious divergent-inclusive catalog **byte-identically**.

**Default FLOOR = 0.80** (`MIN_FAMILY_IDENTITY = 0.80`). Opt out with `--no-divergence-floor`
/ `--min-family-identity=0` / `RUSTLE_NO_DIVERGENCE_FLOOR=1`; set a different floor with
`--min-family-identity=X` (e.g. `0.85`).

---

## The metric (validated on data)

```
recip_id_best = matches_best / max(len_a, len_b)
              = min over the two members of (aligned-AND-matching bases / member length)
```

computed from the **SAME** spliced-exon minimap2 alignment (`-cx asm20`) that already defines
`aln_frac` / `core_recip` — the residue-match count and member lengths are the only piece added
(`bench/ri_build_recip_id.py` -> `bench/ri_sharedlen_recip_id.tsv`; the `aln_len`/`blocklen_best`
column reproduces the universal cache byte-for-byte, so it is the SAME RNA alignment, not a
re-derivation). RNA-only, **NOT DNA**.

It works because it is **RECIPROCAL** (min over both members): a domain-share where a small gene
sits inside one big exon of a large gene scores high on the small gene's side but LOW on the large
gene's side (`matches / len_large`), so the min cuts it. Shared-block identity or `aln_frac` alone
do **not** work — domain-shares share one exon at ~95%.

Class medians on the 1457 core+aln-passing edges (`bench/divergence_floor.tsv`):

| class | n | median `recip_id_best` |
|---|---|---|
| TP — real cDNA paralog | 375 | **0.615** |
| genuine over-merge | 619 | 0.295 |
| truthbar — divergent paralog | 463 | 0.320 |

The reciprocal length penalty separates near-identical real paralogs (0.615) from
divergent/domain-share over-merges (~0.30) by design.

---

## P/R-vs-floor curve (`bench/divergence_floor.tsv`)

Eval via the shipped `bench/family_level_pr_current.py` path (diploid gold oracle). Deterministic
(`PYTHONHASHSEED=0`, gamma 0.20, seed 0). `R_oracle` = diploid multi-copy oracle genes recovered
/ 57; `P_oracle(dedup)` = 1 - distinct-FP / oracle-mapped; `E_p` = protein-family block purity;
domShare = # of the 4 named domain-shares excluded from the multi-copy catalog.

| floor | nFam | R_oracle | P_task / P_dedup | E_p purity | distinctFP | FP a/o/m | domShare excl | md5 |
|-------|------|----------|------------------|-----------|-----------|----------|---------------|-----|
| **off** | 573 | **0.877 (50/57)** | 0.896 / 0.917 | 0.892 | 4 | 0/1/4 | 0/4 | **dca64cbd** (opt-out) |
| 0.70 | 335 | 0.596 (34/57) | 0.879 / 0.909 | 0.961 | 3 | 0/1/3 | 3/4 | 0c433468 |
| 0.75 | 320 | 0.579 (33/57) | 0.909 / 0.939 | 0.966 | 2 | 0/1/2 | 3/4 | 4c8c88ee |
| **0.80** ⭐ | **307** | **0.561 (32/57)** | **0.903 / 0.936** | **0.967** | **2** | 0/1/2 | **4/4** | **e84dc2bc** (default) |
| 0.85 | 293 | 0.509 (29/57) | 0.893 / 0.929 | 0.976 | 2 | 0/1/2 | 4/4 | b27512ce |
| 0.90 | 284 | 0.491 (28/57) | 0.926 / 0.963 | 0.979 | 1 | 0/1/1 | 4/4 | 0b9bf57b |

`floor=off` reproduces the current divergent-inclusive default catalog **dca64cbd BYTE-IDENTICAL**
— so the floored build path is validated and the opt-out is guaranteed exact.

---

## Domain-share exclusion + FP composition

The 4 named single-exon domain-shares and where they fall:

| pair | recip_id_best | excluded at floor |
|---|---|---|
| MOV10 + RHOC | 0.111 = min(590/1166=0.506, 590/5291=0.112) | every floor |
| RBBP4 + SYNC | 0.350 = min(0.605, 0.350) | every floor |
| DEDD + NIT1 | 0.533 = min(0.533, 0.688) | every floor |
| **RHD + SDHD** | **0.785** = min(0.785, 0.879) | **>= 0.80 only** (the tightest) |

RHD+SDHD (0.785) is the tightest domain-share; it is excluded only at **floor >= 0.80**, which is
exactly why 0.80 is the minimum knee (0.75 gives 3/4). Verified at the catalog level (each gene
absent from the multi-copy catalog), not merely a direct-edge cut: floor-OFF co-members RHD+SDHD,
MOV10+RHOC, RBBP4+SYNC, DEDD+NIT1; floored default co-members **none** of the four.

**FP residual is clean and STABLE from 0.80 on.** At 0.80 and 0.85 the distinct diploid-oracle FP
= **2**:
1. `LOC129529978 + LOC129529986` — the irreducible **MAGE-X** DNA-only array over-size floor (no
   RNA metric removes it; it is a genuine collapsed/inflated-copy DNA fact).
2. `GSTM2 + LOC101129940` — one **GST protein-domain hub**.

E_p-impure blocks drop **79 (off) -> 10 (0.80) -> 7 (0.85)**. At 0.90 the GST hub is also cut
(distinctFP -> 1, MAGE only) but at a cost of 4 more oracle genes.

---

## Chosen default: FLOOR = 0.80 (user's 0.85 also validated)

Per the criterion "clean FP residual — domain-shares excluded, high precision — at the best
retained recall", **0.80 is the empirically-best knee and DOMINATES 0.85**:
- identical distinct-FP residual (2 — the same two blocks),
- identical 4/4 domain-share exclusion,
- retains **3 MORE** diploid-oracle genes (32 vs 29, +10% relative recall),
- at essentially equal precision (E_p 0.967 vs 0.976, P_dedup 0.936 vs 0.929),
- and it is the **minimum** floor that excludes RHD+SDHD (0.785).

The user's suggested **0.85 is validated** — it meets every clean-FP bar (4/4 domain-shares out,
distinctFP 2, E_p 0.976) and is a defensible, marginally more conservative alternative — but it
buys no distinct-FP reduction for a real 3-gene recall cost. 0.90 is worth it only if killing the
GST-domain hub (distinctFP -> 1) justifies the extra recall loss. Shipped default = **0.80**;
change with `--min-family-identity=0.85`.

Precision improves vs floor-off on the load-bearing axes: E_p 0.892 -> **0.967**; P_dedup 0.917 ->
**0.936**; distinctFP 4 -> **2**. (The task-formula P is denominator-noisy — `oracle_mapped`
shrinks 48 -> 31 with the catalog — so E_p / distinct-FP are the load-bearing precision signals;
they improve unambiguously.)

---

## The scoped claim

> **Gene families >= 80% reciprocal whole-transcript identity: precision P_oracle(dedup) = 0.936**
> **(E_p protein-purity 0.967; distinct diploid-oracle FP = 2 = the irreducible MAGE-X DNA-only**
> **over-size floor + one GST domain hub; all 4 named domain-shares excluded), recall R_oracle =**
> **0.561 (32/57 diploid multi-copy oracle genes).**

User's more conservative operating point (validated): **"Gene families >= 85% reciprocal
whole-transcript identity: P_oracle(dedup) = 0.929 (E_p 0.976, distinct FP = 2, all domain-shares
excluded), R_oracle = 0.509 (29/57)."**

---

## Opt-out (byte-identical) and honesty

**Opt-out recovers the ambitious mode exactly.** `--no-divergence-floor`, `--min-family-identity=0`,
and `RUSTLE_NO_DIVERGENCE_FLOOR=1` all recover the divergent-inclusive catalog **dca64cbd**
byte-identical (573 families, R_oracle 0.877). Mechanism: when the floor is <= 0 the per-edge
`recip_ok()` returns `True` unconditionally and the reciprocal-identity cache is never loaded —
so the pre-floor code path is bit-for-bit unchanged. The floor **composes** with the other gates
(`--no-repeat-gate`, `--no-split-recombinants`, `--no-repeat-bridge-gate`, `--high-precision`
gamma 0.20 -> 0.40) — each ablation touches only its own gate.

**Honesty (this is the SCOPING, stated plainly).**
- The floor is **RNA-only**: `recip_id_best` comes from the RNA spliced-exon alignment, **NOT DNA**.
  The RNA-only guard hard-asserts the floor feature set is exactly `{recip_id_best}` and disjoint
  from any DNA / protein / genome / soft-mask column.
- Making the floor default **DROPS R_oracle from 0.877 (off) to 0.561 (0.80)** — 18 lost diploid
  oracle genes. This is **BY DESIGN**: it excludes divergent paralogs and partial/length-mismatched
  homologs (SD98-style scoping). The recall is materially reduced but remains useful (32/57
  multi-copy oracle genes at 0.80).
- The ambitious divergent-inclusive catalog is not lost — it is one flag away (`--no-divergence-floor`
  = dca64cbd exact).
- Deterministic: `PYTHONHASHSEED=0`, gamma 0.20, seed 0; the off / 0.80 / 0.85 catalog md5s
  (dca64cbd / e84dc2bc / b27512ce) are stable by construction across runs.

---

## Files (absolute)

- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_rna_refine.py` — the floored edge rule
  (`MIN_FAMILY_IDENTITY = 0.80`, default-ON; `--no-divergence-floor` / `--min-family-identity=X` /
  `RUSTLE_NO_DIVERGENCE_FLOOR=1`).
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/divergence_floor.py` — the measure-stage sweep.
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/divergence_floor.tsv` — the P/R-vs-floor curve.
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/ri_build_recip_id.py` — reciprocal-identity builder.
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/ri_sharedlen_recip_id.tsv` — the reciprocal-identity cache.
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.{json,tsv}` — the floored default
  metrics (307 fam, R_oracle 32/57, P 0.903/0.936, E_p 0.967, distinct-FP 2).
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/test_family_rna_refine.py` — 6 divergence-floor tests
  (opt-out byte-identity, domain-share exclusion, scoped precision improves, env==flag, composes,
  determinism).
