#!/usr/bin/env python3
"""
Robustness of the diploid copy-number oracle to its two thresholds.

The old counter used a single merge_gap that (a) collapsed tandem arrays and
(b) an absolute 1 kb block floor that let a sub-gene domain count as a full copy.
The new counter (diploid_cn_oracle.count_cn) removes both: a copy = a distinct
target region a query covers by >= COV_FRAC * L (L = median locus length), with
query-OVERLAPPING alignments kept as separate tandem units. This script sweeps
COV_FRAC in {0.4, 0.5, 0.6} to show hap_CN is stable (not a knife-edge threshold),
and prints the counts for the families the adversarial review flagged.

Determinism: PYTHONHASHSEED=0; counts are set-based.
"""
import importlib.util
import os
from statistics import median

spec = importlib.util.spec_from_file_location("dco", "bench/diploid_cn_oracle.py")
dco = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dco)

SCRATCH = "/home/juanfra/winloci_scratch"
fam, loc2fam = dco.read_validated(dco.VALIDATED)
mat = dco.read_paf_hits(os.path.join(SCRATCH, "cn_oracle_paf/mat.paf"), loc2fam)
pat = dco.read_paf_hits(os.path.join(SCRATCH, "cn_oracle_paf/pat.paf"), loc2fam)

# families the review flagged: tandem-collapse, over-split, domain-inflation
FLAGGED = {
    "42": "LOC129529768 (tandem array, was mis-filed CN=1)",
    "1": "LOC109025447 (tandem)", "21": "", "22": "", "52": "", "57": "",
    "75": "", "87": "", "175": "",
    "76": "domain inflation (was 22)", "80": "domain inflation (was 37)",
    "14": "LOC115930164 domain inflation (was 37)",
    "2": "over-split (32 ref loci)", "5": "TUBGCP5 over-split",
    "20": "over-split", "65": "ZNF425 reference-collapsed",
    "94": "MAGEA9", "145": "LOC115935025", "58": "DHRSX het_risk",
    "91": "LOC129530050 het_risk",
}


def cn(family, cov):
    L = int(median([e - s for (_l, _c, s, e, _n) in fam[family]]))
    cm, fm = dco.count_cn(mat.get(family, []), L, cov)
    cp, fp = dco.count_cn(pat.get(family, []), L, cov)
    return cm, cp, max(fm, fp)


print("%-5s %-42s %-9s %-9s %-9s %-6s" %
      ("fam", "note", "cov0.4", "cov0.5", "cov0.6", "fragMax(0.5)"))
for f in sorted(FLAGGED, key=int):
    a = cn(f, 0.4)
    b = cn(f, 0.5)
    c = cn(f, 0.6)
    print("%-5s %-42s %-9s %-9s %-9s %-6d" %
          (f, FLAGGED[f][:42], "%d/%d" % a[:2], "%d/%d" % b[:2],
           "%d/%d" % c[:2], b[2]))
