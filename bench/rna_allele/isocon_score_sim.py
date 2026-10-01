#!/usr/bin/env python3
"""Score IsoCon on Arm S (Amendments 1-2): every final candidate against every truth transcript (edlib, candidate as an infix of the truth,
since reads are 5'/3'-shortened). A candidate's match = the truth with the smallest edit distance; tie = several truths at that distance.

  recovered (exact)  truth transcripts with a candidate at edit distance 0
  recovered (0.999)  ... at identity >= 0.999 (1 - ed / candidate length) covering >= 80% of the truth
  alleles separated  copies whose A and B transcripts differ and BOTH are recovered (0.999) by different candidates
  alleles merged     copies whose A and B differ and only one is recovered
  identical alleles  copies whose A and B are identical (one candidate can only ever stand for both: the blind spot)
  ties               candidates whose best match ties between two DIFFERENT copies (paralogs merged)
  extra              candidates whose best identity < 0.999 (consensus errors, chimeras)

    isocon_score_sim.py --truth NPIP.truth.fa --cands sim_NPIP/final_candidates.fa
"""
import argparse
import collections

import edlib


def fasta(p):
    s, cur = {}, None
    for ln in open(p):
        if ln[0] == ">":
            cur = ln[1:].split()[0]; s[cur] = []
        else:
            s[cur].append(ln.strip())
    return {k: "".join(v) for k, v in s.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--cands", required=True)
    a = ap.parse_args(argv)
    T, C = fasta(a.truth), fasta(a.cands)
    best = {}
    for cid, cs in C.items():
        res = []
        for tid, ts in T.items():
            r = edlib.align(cs, ts, mode="HW", task="distance")
            res.append((r["editDistance"], tid))
        res.sort()
        d0 = res[0][0]
        tied = [t for d, t in res if d == d0]
        best[cid] = (d0, tied, 1 - d0 / max(1, len(cs)), len(cs))
    rec_exact = {t for d, ts, i, L in best.values() if d == 0 for t in ts}
    rec = set()
    for cid, (d, ts, ident, L) in best.items():
        for t in ts:
            if ident >= 0.999 and L >= 0.8 * len(T[t]):
                rec.add(t)
    copy = lambda t: t.rsplit("_", 1)[0]
    pairs = sorted({copy(t) for t in T if t.endswith("_A") and copy(t) + "_B" in T})
    ident_pairs = [c for c in pairs if T[c + "_A"] == T[c + "_B"]]
    diff_pairs = [c for c in pairs if c not in ident_pairs]
    sep = [c for c in diff_pairs if c + "_A" in rec and c + "_B" in rec]
    merged = [c for c in diff_pairs if (c + "_A" in rec) != (c + "_B" in rec)]
    neither = [c for c in diff_pairs if c + "_A" not in rec and c + "_B" not in rec]
    para_ties = [cid for cid, (d, ts, i, L) in best.items() if len({copy(t) for t in ts}) > 1 and
                 not all(copy(t) in ident_pairs and len({copy(x) for x in ts}) == 1 for t in ts)]
    extra = [cid for cid, (d, ts, i, L) in best.items() if i < 0.999]
    print(f"truth transcripts {len(T)}; candidates {len(C)}")
    print(f"recovered exact {len(rec_exact)}/{len(T)}; recovered >=0.999 {len(rec)}/{len(T)}")
    print(f"copies with differing alleles {len(diff_pairs)}: separated {len(sep)}, merged {len(merged)}, neither recovered {len(neither)}")
    print(f"copies with identical alleles (blind spot) {len(ident_pairs)}")
    print(f"candidates tied between different copies (paralogs merged): {len(para_ties)}; extra candidates (< 0.999): {len(extra)}")
    bonly = [t for t in T if "Bonly" in t]
    if bonly:
        print(f"B-only (reference-absent) transcripts recovered: {sum(1 for t in bonly if t in rec)}/{len(bonly)}")
    per = collections.Counter(t for d, ts, i, L in best.values() if i >= 0.999 for t in ts)
    print("candidates per recovered truth:", dict(collections.Counter(per.values())))


if __name__ == "__main__":
    main()
