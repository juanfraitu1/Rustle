# Oracle near-miss analysis — what it would take, and the empirical-rule verdict

**Date:** 2026-06-09. Baseline -L, 4 mid-size chroms. Tooling: oracle_nearmiss.py + oracle_support_sweep.py.

## Oracle half — what it would take (structural)
- 3,724 refs matched '='; **1,475 near-miss-only** (rustle shares a junction, no exact match).
- **529 (36%) are exactly 1 junction edit away** (both-axes-positive to convert: +TP/-FP):
  c_extend 191 (add missing junction) · altsplice 173 (snap boundary) · k_trim 165 (drop extra).
- ~3,000+ single-edit near-misses extrapolated genome-wide.

## Empirical-rule half — read-support sweep: NEGATIVE on all three
- **k_trim** (drop extra junction): 34% of extra junctions weak (reads/body <0.1), BUT 27% of
  MATCHED refs carry a real junction <0.1 too → a drop-rule breaks ~950 correct matches to fix 56.
  Indistinguishable.
- **altsplice** (snap to ref): read support points the WRONG way — rustle already picks the
  HIGHER-support junction; the annotation uses a LOWER-support one. "Prefer higher" wrong 85%.
  The annotated isoform is the minority in this sample's reads.
- **c_extend** (add missing junction): 45% of missing junctions have 0 reads (genuinely absent);
  median 1 read; only 32% (61) have >=2 reads = the flow-skipped-a-supported-junction cases =
  read-chain territory.

## Verdict
No read-support threshold rule converts near-misses without breaking correct calls — same
indistinguishability wall, confirmed from the both-axes near-miss angle. The ONLY orthogonal
escape is **read-chain** (per-molecule collapse, no support threshold), which maps to the
c_extend-supported 32% and the earlier flow_enumeration finding (project_readchain_flow_enum).

## CORRECTION via instrumented ST comparison (2026-06-09)
Traced 3 altsplice near-misses where ST='=' / rustle=near-miss (NC_073227.2): all are tiny
3-5bp alt-donor shifts. Instrumented junction_accept (both tools) on XM_055384841.2:
BOTH tools accept BOTH donors with identical support — 48967721 (annotated, 40 reads) and
48967718 (28 reads). **ST's path emits the higher-support donor (40 → annotated '='); rustle's
flow emits the lower one (28 → near-miss).** So the divergence is PATH SELECTION (parse_trflong),
NOT junction acceptance, and NOT the -E snap window (-E 0/5/10/30 all fail to recover).
**KEY:** my earlier "altsplice support points the wrong way (85%)" used RAW BAM CIGAR counts —
the WRONG metric. The tools' good-junction support (nreads_good/anchor) points toward the
annotation. The altsplice rule may be VIABLE if re-measured with the tool's good-junction
support. Re-measure: run both tools' junction_accept on the altsplice loci, compare ref-donor
vs rustle-donor by nreads_good, check if "prefer higher good-support" matches annotation.

## Corrected altsplice sweep (rustle good-junction support, not raw BAM)
ST-converts altsplice cases (ST='=', rustle near-miss) are RARE: only 6 across 4 chroms (most
altsplice near-misses, BOTH tools miss). With rustle's junction_accept `mm` support:
"prefer higher good-support at the shared boundary" matches annotation 3/6 (50%) — up from the
raw-BAM 12% (confirming that metric was wrong), but STILL only 50% = not a separator. ST's
parse_trflong path-selection is more nuanced than "prefer higher support." Lever is tiny (~40
genome-wide) and 50/50 => no net both-axes gain. **altsplice rule = NOT VIABLE even corrected.**
The instrumented comparison's value: corrected a wrong metric + localized the divergence to
parse_trflong path-selection (junctions rustle already has), reconfirming the architectural lever.
