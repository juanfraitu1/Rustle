#!/usr/bin/env python
"""DNA-truth help-vs-hurt test for the "suspect" (alignment/error-noisy) PSV columns.

The align_error_probe found that down-weighting low-coherence columns would change ~34% of assignments and
those columns are LOAD-BEARING — but "load-bearing" doesn't prove they're CORRECT (a load-bearing artifact also
flips reads). This test settles it with NON-CIRCULAR DNA truth (DNA-derived copy signatures, SDA/Vollger):

  - SUSPECT column = high OFF-SIGNATURE rate: a column where many reads show a base NO DNA copy has there
    (an error / mis-alignment signature, defined WITHOUT reference to any assignment).
  - For each read, split its CLEAN (non-suspect) observed columns into TRAIN_clean and a held-out TEST.
  - call_WITH    = decode on TRAIN_clean + the read's SUSPECT columns.
  - call_WITHOUT = decode on TRAIN_clean only.
  - For reads where the suspect columns CHANGE the call, ask which call the held-out CLEAN TEST confirms.
    TEST confirms WITH more  -> suspect columns HELP (keep them; the naive down-weight flag would hurt).
    TEST confirms WITHOUT more -> suspect columns HURT (model/drop them; the flag is worth building).

Run: /home/juanfra/miniforge3/bin/python bench/align_error_dna_test.py [max_families]
"""
import sys
from collections import defaultdict

sys.path.insert(0, ".")
sys.path.insert(0, "bench")
import dna_supervised_decode as d
import copy_assign as ca

OFF_SUSPECT = 0.05   # a column is suspect if > this fraction of reads show a base no DNA copy has
MIN_OBS = 5          # need >= this many reads observing a column to classify it


def best_on(cols, obs, copy_vecs):
    """decode a read using only `cols`; return (best_copy, n_decisive) or (None,0) if no decisive col."""
    sub = {c: obs[c] for c in cols if c in obs}
    if not sub:
        return None, 0
    best, _margin, _n, n_dec, _logL = ca.assign_read(sub, copy_vecs)
    return (best, n_dec)


def main():
    maxf = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    fams = d.colocated_families()
    changed = wins_with = wins_without = ties = 0          # over reads where suspect cols change the call
    n_reads = n_susp_cols = n_clean_cols = nfam = 0
    import json, os
    for fid, members in fams:
        if nfam >= maxf:
            break
        sig = d.build_signatures(fid, members)
        if not sig:
            continue
        chrom, ref_seq, psv_gpos, copy_vecs = sig
        rs = json.load(open(os.path.join(d.DNADIR, f"{fid}.json")))["iv"][0]
        obs_list = d.read_alleles_via_ref0(fid, chrom, members, ref_seq, psv_gpos, rs)
        if len(obs_list) < 10:
            continue
        nfam += 1
        ncol = len(psv_gpos)
        # off-signature per column (the artifact signature), over all reads
        off = [0] * ncol
        seen = [0] * ncol
        copy_alleles = [{copy_vecs[lab].get(j) for lab in copy_vecs} for j in range(ncol)]
        for obs in obs_list:
            for j, a in obs.items():
                seen[j] += 1
                if a not in copy_alleles[j]:
                    off[j] += 1
        suspect = {j for j in range(ncol) if seen[j] >= MIN_OBS and off[j] / seen[j] > OFF_SUSPECT}
        n_susp_cols += len(suspect)
        n_clean_cols += ncol - len(suspect)
        for obs in obs_list:
            n_reads += 1
            cols = sorted(obs)
            clean = [c for c in cols if c not in suspect]
            susp = [c for c in cols if c in suspect]
            if len(clean) < 3 or not susp:           # need clean train+test AND a suspect col to matter
                continue
            test = clean[0::2]                        # held-out CLEAN columns (the non-circular judge)
            train_clean = clean[1::2]
            cw, dw = best_on(train_clean + susp, obs, copy_vecs)
            co, do = best_on(train_clean, obs, copy_vecs)
            if cw is None or co is None or cw == co:
                continue                              # suspect cols didn't change the call
            changed += 1
            tb, td = best_on(test, obs, copy_vecs)    # held-out clean verdict
            if tb is None or td < 1:
                continue
            if tb == cw:
                wins_with += 1
            elif tb == co:
                wins_without += 1
            else:
                ties += 1
    print(f"== DNA-truth help-vs-hurt ({nfam} families, {n_reads} reads) ==")
    print(f"  suspect (off-signature>{OFF_SUSPECT}) cols: {n_susp_cols} / {n_susp_cols+n_clean_cols} "
          f"= {100*n_susp_cols/max(n_susp_cols+n_clean_cols,1):.1f}%")
    print(f"  reads where suspect cols CHANGE the call: {changed}")
    print(f"  held-out CLEAN test confirms the call made WITH suspect cols:    {wins_with}")
    print(f"  held-out CLEAN test confirms the call made WITHOUT suspect cols: {wins_without}  (ties/other {ties})")
    if wins_with + wins_without:
        frac = wins_with / (wins_with + wins_without)
        print(f"  -> WITH-suspect wins {100*frac:.1f}% of decisive cases  => suspect columns "
              f"{'HELP (keep them; naive down-weight would HURT)' if frac > 0.5 else 'HURT (model/drop them)'}")


if __name__ == "__main__":
    main()
