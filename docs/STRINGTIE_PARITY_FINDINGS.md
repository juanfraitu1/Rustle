# StringTie Parity — Authoritative Findings & Ceilings

**Status:** current authority for the StringTie long-read parity effort (de novo, GGO_19 chr19).
**Last updated:** 2026-05-28. Supersedes the precision/parity-plan docs listed at the bottom.

This document is the single source of truth for: where Rustle stands vs StringTie, what the
remaining gaps are, what their *bounded ceilings* are, and which levers are exhausted vs open.
Detailed working trace: `docs/superpowers/specs/2026-05-28-pathpat-flow-parity-scope.md`.

## 1. Current operating point

De novo `rustle -L GGO_19.bam` vs `GGO_19.gtf` (StringTie's output), gffcompare intron-chain:

- **Sn 96.5 / Pr 90.7 / F1 93.5** (1974 transcripts).
- Rustle sits at a **higher-sensitivity / lower-precision** point than StringTie — and that is *partly a feature*: Rustle's flow recovers ~40 real reference isoforms StringTie does not extract.

**Shipped win:** `parity/isofrac-chain-dedup` branch (committed + pushed) — the isofrac multicov
**chain-dedup** fix: count each unique intron chain once in the isofrac denominator (StringTie
collapses dominant 3'-end variants before isofrac; Rustle now matches). +12 tx, Sn 96.1→96.5,
Pr 91.0→90.7, F1 +0.03; recovers STRG.267.3 / 343.2 / 15.1. Default ON, opt out via
`RUSTLE_ISOFRAC_CHAIN_DEDUP_OFF=1`.

## 2. The precision gap, fully characterized (~180 extra chains)

Every extra Rustle emits that StringTie does not was attributed (bench/pathpat_phase0*.py):

| bucket | count | mechanism |
|---|---:|---|
| filter divergence (co-extracted, ST kills) | ~97* | ST extracts the same chain then kills it; Rustle keeps it |
| alt-splice junction selection | ~48 | Rustle's flow generates a different/extra junction; near already-matched refs |
| subset/contained | ~32 | Rustle drops a junction (shorter chain) |
| genuine extra-junction over-generation | ~9 | true over-enumeration |

\* The ~97 "filter divergence" count is **inflated** by a 60bp coordinate-tolerance match in the
attribution; the genuinely filter-recoverable set is ~2 (see §3, bucket-A oracle).

## 3. Lever ceilings — what each fix could achieve (all oracle-bounded)

The reusable method: before building an expensive fix, feed StringTie's *exact* values in as an
**oracle** and measure the F1 ceiling. This killed two false starts cheaply.

| lever | ceiling (oracle) | verdict |
|---|---|---|
| **chain-dedup** | +0.03 (actual) | ✅ shipped — the one clean win |
| **bucket-A lowintron** (read-overhang bpcov) | **+0.1 F1** | ❌ falsified: oracle with ST's exact lowintron masks (93% match) gives only +0.1 Pr / 0 Sn / −2 j |
| **targeted FP-only suppression** | — | ❌ no discriminator: FP over-extractions (min_jct_mm med 4, longcov mostly 1) sit *inside* the real-isoform distribution; any gate that kills 31 FPs kills ~360 real isoforms |
| **flow path-enumeration rewrite** | **+1.1 F1** | ⚠️ open, but it's a **reshape**: Pr +4.2 / Sn −2.2, *loses ~40 real isoforms*; highest-risk, multi-session |

## 4. Abundance values (longcov / TPM / cov)

Measured on 1748 matched chains (bench, corrected GTF parsing — note ST writes `longcov` *before*
`cov`, so naive `cov "…"` regex mis-parses):

- **cov**: median ru/st 0.965, 95% within 2× — **well-aligned**. (The old "60× tlen-inflation" claim
  was largely a parsing artifact on final matched transcripts.)
- **TPM**: systematic **~1.26×**. Identical formula (`cov·1e6/Cov_Sum`); the gap is the denominator —
  Rustle sums final-transcript cov (clean); ST accumulates a per-bundle `sum_cov` over *pre-filter*
  predictions + gene-level sums (ST's own code: *"isn't this double counting?"*). Relative TPMs are
  preserved; ST's denominator is arguably *less* correct. Aligning = copying ST's messier accounting.
- **longcov**: 68% exact, 12% diverge (some extreme, e.g. ST STRG.26.1 longcov 1 vs Rustle 1211).
  **Tried** switching Rustle's longcov from `read_count` to `abundance` (ST's exact field): **no-op** —
  Rustle's transfrag *abundance itself* is 1211 vs ST's 1. The divergence is in read-to-transfrag
  assignment, upstream of longcov. Reverted.

## 5. The unifying conclusion

**Every remaining StringTie-parity gap — precision extras, longcov, TPM — traces to one foundation:
the read-to-isoform / read-to-transfrag assignment in the flow decomposition.**

- cov agrees on matched chains because the *path* is the same.
- longcov diverges where the flow attributes a different read mass to the chain's transfrag.
- the precision extras exist because the flow traverses alt-junction paths ST doesn't.
- the FP and real-minor isoforms are **indistinguishable** by any available feature (cov, longcov,
  min_jct_mm) — confirmed repeatedly.

There is **no cheap filter- or normalization-level lever left**. The only lever with a >1-F1 ceiling
is matching StringTie's read-assignment (the flow rewrite), and that is a sensitivity-for-precision
**reshape**, not a free win — a product decision about which operating point is desired, not a bug fix.

## 6. Instrumentation toolkit (kept, env-gated; default behavior unchanged)

- **Cross-tool parity diff**: `RUSTLE_PARITY_LOG` / `STRINGTIE_PARITY_LOG` (+ `_FILTER_CHROM/RANGE/STEPS`)
  emit JSONL decision events (path_extracted, pred_kill, checktrf_*, junction_accept, …) joined on
  `(step,chrom,start,end,strand)`.
- **Coverage/lowintron tracing**: `RUSTLE_TRACE_COV_NODES`, `RUSTLE_RI_TRACE`+`RUSTLE_TRACE_LOCUS`,
  `RUSTLE_BLF_TRACE`, `RUSTLE_JCTMM_DUMP`; StringTie side `ST_TRACE_COV_NODES`, `ST_ILMASK` (rlink.cpp).
- **Oracles**: `RUSTLE_LOWINTRON_ORACLE=<ST mask file>` (feed ST's lowintron in);
  `RUSTLE_LEFTOVER_REDIST`, `RUSTLE_RI_ZERO_INTRON_LOW`, `RUSTLE_LONGCOV_READCOUNT` (opt-in experiments).
- **Analysis scripts**: `bench/pathpat_phase0.py` / `phase0c.py` (extra attribution),
  `bench/pathpat_bucketB.py` (junction-shift characterization).
- **Build note (WSL):** StringTie submodule — `make clean release` (10s/run; DEBUG is ~9min). `/mnt/c`
  mtimes are flaky so `make` may skip recompiling edited `rlink.cpp`; `touch rlink.cpp` or `make clean`.
  Do NOT use `-o /dev/null` with StringTie (temp-dir error) — use a real path.

## 7. Superseded documents

The following are superseded by this file for the precision/parity-gap analysis (kept for history):
`EXACT_STRINGTIE_PARITY_PLAN.md` (its "100% parity via guided mode" premise is circular — `-G GGO_19.gtf`
scored against `GGO_19.gtf`), `PRECISION_DIAGNOSTIC_FINDINGS.md`, `PRECISION_IMPROVEMENT_PLAN.md`,
`PRECISION_WORK_CONTINUATION.md`, `docs/STRINGTIE_PARITY_SYSTEMATIC.md`,
`docs/STRINGTIE_LONGREAD_PARITY_CONTINUATION_GUIDE.md`, `docs/SUBBUNDLE_PARITY_HANDOFF.md`.
