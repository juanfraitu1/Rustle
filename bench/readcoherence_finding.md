# Read-coherence + degradation-aware extraction — finding (2026-06-14)

**Status: DOCUMENTED, PARKED.** Strong empirical result, but parked in favor of the
multimapping/VG thesis direction (the advisor's interest). Resume from here.

## Result (genome-wide GGO, guided `-G stringtie.gtf`, SQANTI3-validated)

Read-coherence = `rustle -G st.gtf --read-chain` (per-molecule path extraction), compared
additively against the `-G` flow baseline (`gd`). The realness verdict (SQANTI3, paralog
universe + structural categories + canonical junctions + RT-switch):

| metric | value |
|--------|-------|
| additive multi-exon extras over flow baseline | **10,144** (drops only 85 flow chains) |
| canonical junctions | 9,650 / 10,144 (**95%**) |
| **STRICT real** (FSM/NIC + canonical + non-RT-switch) | **1,857 (18%)** |
| REAL (+ NNC canonical) | 4,685 (46%) |
| of real extras, cov≥5 | 3,607 / 4,685 (**77%** — well-supported) |

Of the 1,857 strict-real: **503 FSM** (exact RefSeq matches *both* StringTie and rustle-flow
missed — indisputable recall) + ~1,354 NIC (real novel combinations of known junctions).

**Contrast with VG-multimapper (same harness, guided): 0** FSM paralog-copy recoveries over
baseline. Read-coherence is the bigger recall lever by ~3 orders of magnitude — but it is
**not** the advisor's interest (multimapping is), hence parked.

## Noise anatomy (the ~5,459 non-real extras) — filter-vs-generate split

| artifact | count | fix path |
|---|---|---|
| ISM truncations (3′ 794 / 5′ 698 / internal 636 / IR 95) | 2,223 | **better generation** — degradation-aware collapse (fold fragments into full-length parent). rustle's existing 5′-degrade collapse (`RUSTLE_READCHAIN_DEGRADE_TOL`, default-on) is incomplete: leaves 3′ + internal + residual 5′. |
| RT-switching | 1,538 | **filter** (data artifact; the false junction is genuinely in the read — annotation-free detectable via genome sequence around the junction) |
| non-canonical junctions | 494 | **filter** (canonical GT-AG check, annotation-free) |
| intergenic / antisense (~55% single-read) | 1,331 | read-depth + locus-context gate |
| fusion (read-through) | 510 | read-through logic |

**Verdict:** hybrid. ~2,128 (ISM degradation) avoidable by better *generation* (the
"degraded-transcript modelling" idea); ~2,032 (RT-switch + non-canonical) are *data* artifacts
needing a *filter* (no extractor avoids them). Neither alone suffices.

## Proposed design (when resumed)

Shippable #1 = **degradation-aware read-coherent extraction** (extend the degrade-collapse to
3′ + internal fragments; the biggest, cleanest noise win, doubles as the 3′-degradation
feature) **+ an annotation-free realness gate** (canonical junctions + RT-switch detection +
read-depth). Additive over the byte-exact `-G` flow floor, opt-in/gated like the other
"better decisions", `RUSTLE_PRECISE`-exempt. SQANTI3 stays the validator each iteration.
Expected yield: 10,144 → ~4,000–4,700 real isoforms genome-wide over the StringTie-exact floor.

## Reproduce
- generate: `bench/rc_gen.sh` (per-chrom `rustle -G st --read-chain`, OOM-safe) → `/tmp/gw/rc_*.gtf`
- validate: `/tmp/cre_guided/rc_validate.sh` (merge + SQANTI3 + realness verdict)
- harness reuses the copy-recovery SQANTI3 stages + cached paralog universe.
