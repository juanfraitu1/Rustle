#!/usr/bin/env python3
"""Score per-copy read attribution for the synthetic paralog benchmark.

Given the fingerprint-EM per-read TSV (RUSTLE_VG_FP_ATTR_TSV) and the labeled
reads (names `c{src}_r{j}` from gen_synthetic.py), measure whether the EM assigns
each AMBIGUOUS (multimapping) read to its TRUE source copy.

Why multimappers only: a read minimap2 placed on a single copy is trivially and
correctly attributed (it carries enough copy-specific signal to map uniquely).
The EM's job — and the thesis question — is the AMBIGUOUS reads that map to ≥2
near-identical copies. Those are exactly the reads in the TSV.

The family copy index is an arbitrary internal label (bundle-discovery order), so
predicted-copy → source is resolved by the OPTIMAL BIJECTION over the confusion
matrix (Hungarian; greedy fallback). In the non-identifiable regime the EM
collapses every source onto one copy → the bijection can satisfy only one source
→ accuracy falls toward chance, which is the correct signal (collapse is NOT
hidden). Per read the PREDICTED copy is its max-`final_weight` placement.

Usage:
  score_attribution.py --tsv FP.tsv --reads reads.fa --meta meta.json [--json out.json]
"""
import argparse, json
from collections import defaultdict

def fnv1a64(s: bytes) -> int:
    h = 0xcbf29ce484222325
    for b in s:
        h ^= b
        h = (h * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def optimal_bijection(conf):
    """conf[src][copy] counts. Return (copy->src map, matched_total) maximizing the
    diagonal. Tries scipy Hungarian; greedy fallback."""
    n = len(conf)
    try:
        from scipy.optimize import linear_sum_assignment
        import numpy as np
        m = np.array(conf, dtype=float)
        # maximize matches → minimize negative
        rows, cols = linear_sum_assignment(-m)
        copy2src = {int(c): int(r) for r, c in zip(rows, cols)}
        matched = int(sum(conf[r][c] for r, c in zip(rows, cols)))
        return copy2src, matched
    except Exception:
        # greedy: repeatedly take the largest unused (src,copy) cell
        cells = sorted(((conf[r][c], r, c) for r in range(n) for c in range(n)),
                       reverse=True)
        used_r, used_c, copy2src, matched = set(), set(), {}, 0
        for cnt, r, c in cells:
            if r in used_r or c in used_c:
                continue
            used_r.add(r); used_c.add(c); copy2src[c] = r; matched += cnt
        return copy2src, matched

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--json", default=None)
    ap.add_argument("--min-decisive-sites", type=int, default=2,
                    help="a read must cover >= this many diagnostic sites to count as DECISIVE "
                         "(matches the EM's RUSTLE_VG_EM_MIN_DECISIVE_SITES gate). A read resting "
                         "on fewer sites is unreliable — one sequencing error flips a lone site — "
                         "so it ABSTAINS even if its prior-driven weight gap looks confident.")
    a = ap.parse_args()

    meta = json.load(open(a.meta))
    n_copies = meta["copies"]

    src_by_hash, n_reads_total = {}, 0
    with open(a.reads) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                src = int(name[1:].split("_")[0])
                src_by_hash[fnv1a64(name.encode())] = src
                n_reads_total += 1

    rows, gap_by_hash, sites_by_hash = defaultdict(list), {}, {}
    with open(a.tsv) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            h = int(c[1]); copy = int(c[2]); nsites = int(c[3]); w = float(c[4]); gap = float(c[5])
            rows[h].append((copy, w)); gap_by_hash[h] = gap
            # max diagnostic sites the read covers across its placements (the read,
            # not the placement, determines coverage).
            sites_by_hash[h] = max(sites_by_hash.get(h, 0), nsites)

    # confusion over multimapper reads in the TSV
    conf = [[0] * n_copies for _ in range(n_copies)]
    per_read = []  # (src, pred_copy, gap, max_sites)
    for h, src in src_by_hash.items():
        if h not in rows:
            continue
        pred, _ = max(rows[h], key=lambda t: t[1])
        if pred >= n_copies:
            continue
        conf[src][pred] += 1
        per_read.append((src, pred, gap_by_hash.get(h, 0.0), sites_by_hash.get(h, 0)))

    n_mm = len(per_read)
    copy2src, matched = optimal_bijection(conf)
    acc = matched / n_mm if n_mm else 0.0

    # decisive vs abstain among multimappers, using the resolved mapping. A read is
    # DECISIVE only if its weight gap is large AND it covers >= min_decisive_sites
    # diagnostic sites — the evidence gate. A confident gap on too-few sites is
    # prior-driven, not evidence-driven, so it abstains (matches the EM reporting fix).
    n_dec = n_dec_ok = n_abst = 0
    for src, pred, gap, nsites in per_read:
        ok = (copy2src.get(pred, -1) == src)
        if gap > 0.8 and nsites >= a.min_decisive_sites:
            n_dec += 1; n_dec_ok += ok
        elif gap < 0.5 or nsites < a.min_decisive_sites:
            n_abst += 1
    dec_acc = n_dec_ok / n_dec if n_dec else 0.0

    out = {
        "identity_target": meta["identity_target"],
        "identity_realized_exonic": meta["identity_realized_exonic"],
        "copies": n_copies,
        "min_decisive_sites": a.min_decisive_sites,
        "n_reads_total": n_reads_total,
        "n_multimapper_reads_scored": n_mm,
        "n_unique_reads": n_reads_total - n_mm,  # trivially-correct (single-copy) reads
        "attribution_accuracy_multimappers": round(acc, 4),
        "chance_level": round(1.0 / n_copies, 4),
        "decisive_accuracy": round(dec_acc, 4),
        "decisive_frac": round(n_dec / n_mm, 4) if n_mm else 0.0,
        "abstain_frac": round(n_abst / n_mm, 4) if n_mm else 0.0,
        "confusion_src_x_copy": conf,
        "resolved_copy_to_src": {str(k): v for k, v in copy2src.items()},
    }
    print(json.dumps(out, indent=2))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=2)

if __name__ == "__main__":
    main()
