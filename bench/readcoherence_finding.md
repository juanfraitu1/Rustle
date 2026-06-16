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

## VALIDATED vs StringTie (2026-06-15, rigorous cut — survives the PSV-grade scrutiny)

gffcompare rc_$C vs RefSeq minus st_$C vs RefSeq, genome-wide (25 chroms):
- StringTie FSM(=)=23,371; read-coherence FSM(=)=23,951.
- **FSM (exact) read-coherence finds that StringTie MISSES: 580** (≈ the parked 503 estimate, confirmed).
- broad (=/c/j) StringTie misses: **2,784**.
- FSM = exact annotated-transcript match ⇒ canonical+real by construction (NOT the PSV non-canonical-
  artifact failure mode). This is the real beat-StringTie lever (~100–200× the PSV margin of 2–3 exact).
- COST: the 580 sit inside ~10,144 raw extras (~half noise) ⇒ need the annotation-free realness gate.

## SHIPPED — gated `--read-coherence` layer VALIDATED (2026-06-15, genome-wide 25 chroms)

The gated layer (degrade-collapse 3′+internal + annotation-free realness gate + holdout/union-back
additivity) is built, wired (default-off byte-identical, 422 lib tests), and validated. Rigorous
FSM(=) cut, `rcg`/`rc`/`gd`/`st` vs RefSeq (`bench/rcg_gen.sh` + `bench/rcg_validate.sh`):

| arm | FSM(=) total | beat-ST cut (\\ st_FSM) | total tx |
|---|---|---|---|
| StringTie (`st`) | 23,372 | — | — |
| flow baseline (`gd`) | 23,549 | 177 | — |
| ungated read-chain (`rc`) | 23,952 | **580** (reproduces the documented baseline ✓) | 78,559 |
| **gated read-coherence (`rcg`)** | **24,373** | **1001** | **76,680** |

- `st ⊆ rcg` exactly (guides preserved) ⇒ the 1001 is **pure additive recall over StringTie**.
- **Degrade-collapse is the hero:** gated finds +589 *new* FSM over ungated (folding degraded
  fragments back into full-length parents promotes ISM→FSM), at the cost of 168 ungated FSM the
  realness gate over-kills (71% retention). Net +421 over ungated, with ~1,900 FEWER total tx.
- **CRITICAL-1 cost (flow replacement):** `--read-coherence` REPLACES flow extraction per-bundle
  (pipeline.rs:7283 `RUSTLE_READCHAIN` else-if), so it gives up 96 flow-only FSM (all ST-misses).
  Truly-additive ceiling (`rcg ∪ flow`) = **1097** (+96). Option B = produce flow AND read-chain
  per bundle + merge (stronger "⊇ flow baseline" claim for the adversarial advisor).
- **Gate over-kill (168) partly fixable:** `is_canonical_junction` only accepts GT-AG/CT-AC, missing
  the minor-canonical GC-AG / AT-AC (U12) classes → some of the 168 are real GC-AG isoforms.

## OPTION B (truly-additive) + GC-AG fix — VALIDATED (2026-06-15, user chose "max")

Re-architected the layer from read-chain-REPLACES-flow to truly-additive `flow ∪ read-chain`:
per-bundle dual extraction (read-chain on a transfrag clone so flow sees pristine input) +
a per-bundle hold-out (`rc_layer_buf`) + a re-introduced global hold-out/union-back
(`read_coherence_holdout`, mirrors `union_baseline_holdout`) so read-chain NEVER displaces a
flow/guide transcript at any filter stage. Also fixed `is_canonical_junction` to accept the
minor canonical classes GC-AG / AT-AC (+ revcomps), recovering real isoforms the gate over-killed.

Genome-wide rigorous FSM(=) cut (25 chroms):

| arm | FSM(=) | beat-ST (\\ st) | total tx |
|---|---|---|---|
| StringTie | 23,372 | — | 68,157 |
| flow (gd) | 23,549 | 177 | 70,373 |
| ungated read-chain (rc) | 23,952 | 580 | 78,559 |
| replacement gated (prev) | — | 1,001 | 76,680 |
| **OPTION B (rcg)** | **25,620** | **2,248** | **94,145** |

- **+2,248 exact annotated isoforms StringTie misses** (~4× ungated, ~2.2× replacement). The jump:
  read-chain now bypasses the pairwise/isofrac/predcluster filters (gated only by canonical +
  RT-switch + read-depth), which were culling real FSM in the replacement path.
- **⊇ flow guarantee holds genome-wide: `gd \ rcg = 0` AND `st \ rcg = 0`** (never loses a flow
  find or a StringTie guide) — the provable additivity claim for the adversarial advisor.
- **Precision cost:** 94,145 tx vs StringTie 68,157 (+38%); = 70,373 flow floor + 23,791 gated
  read-chain extras. The gate `min_cov` (default 2, env RUSTLE_READ_COHERENCE_MIN_COV) is the
  tunable precision lever; raising it trades read-chain recall for fewer transcripts.
- Default-off BYTE-IDENTICAL (chr19 off==gd exactly); RUSTLE_PRECISE airtight (BOTH read-chain
  extraction branches precise-gated; precise+rc-env == precise-plain, diff 0); 424 lib tests.

### Precision of the 23,791 read-chain extras (SQANTI3, genome-wide — answering "is the +38% real?")
| structural_category | count | % |
|---|---|---|
| novel_not_in_catalog (NNC) | 10,059 | 42% |
| novel_in_catalog (NIC) | 6,275 | 26% |
| incomplete-splice_match (ISM) | 2,427 | 10% |
| full-splice_match (FSM) | 2,066 | 8% |
| intergenic | 1,134 | 4% |
| fusion | 949 | 4% |
| antisense | 709 | 3% |
| genic_intron / genic | 172 | 1% |

- **100% canonical** (all 23,791 — the gate guarantees it). RT-switch (SQANTI3 RTS): 3,971 (17%).
- **REAL (FSM/NIC/NNC + canonical + non-RT) = 14,999 (63%)**; STRICT (FSM/NIC + …) = 6,731 (28%).
- FSM+NIC+NNC = 18,400 (**77%** are matched/novel isoforms); likely-noise (intergenic+antisense+
  fusion+genic_intron) = 2,909 (**12%**); ISM/genic (degradation/partial) = 2,482 (10%).
- The gate IMPROVED composition vs the ungated 10,144: real 46%→63%, strict 18%→28%, canonical 95%→100%.
- Residual precision lever: rustle's narrow 8bp RT-switch heuristic misses the 3,971 SQANTI3 calls;
  a wider/two-orientation RTS detector (or higher gate min_cov) would shed those.

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
