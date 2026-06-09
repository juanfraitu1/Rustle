"""Dissect the graph_missing category (ref junctions absent from rustle's graph).

For each graph_missing miss, find the ref intron(s) rustle's junction_accept set lacks, and
classify WHY by cross-referencing both tools' junction_accept logs:
  ru_rejected : rustle logged the junction but accepted=false -> junction FILTER lever (tunable)
  ru_unseen   : rustle never logged the junction at all -> upstream read/alignment (reads don't
                support it in rustle's view; multimapping/clipping/strand)
For ru_unseen, ST's nreads on that junction shows how much evidence ST had (big gap = alignment
divergence; small = marginal).
"""
import json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def junc_map(path):
    """junction (start,end) -> (accepted_bool, nreads) from a tool's junction_accept log."""
    m = {}
    for j in lib.read_parity(path, "junction_accept"):
        key = (j["start"], j["end"])
        p = j.get("payload", {})
        m[key] = (bool(p.get("accepted", False)), float(p.get("nreads", 0)))
    return m

def main():
    attr = [json.loads(l) for l in open(f"{lib.ROOT}/bench/flow_recall_phase0/attribution.jsonl")]
    gm = [a for a in attr if a["cause"] == "graph_missing"]
    sub = Counter()
    by_cat = defaultdict(Counter)
    st_nreads_unseen = []
    rows = []
    for a in gm:
        chrom, tid = a["chrom"], a["tid"]
        wd = f"{lib.CACHE}/{chrom}/{tid}"
        ch = lib.ref_chain(chrom, tid)
        if ch is None:
            continue
        ref_introns = ch[2]
        ru = junc_map(f"{wd}/r.jsonl")
        st = junc_map(f"{wd}/s.jsonl")
        ru_accepted = {k for k, (acc, _) in ru.items() if acc}
        missing = [iv for iv in ref_introns if iv not in ru_accepted]
        # classify each missing intron; assign the case its "worst" (most upstream) reason
        case_reason = None
        for iv in missing:
            if iv in ru:                       # logged but not accepted
                reason = "ru_rejected"
            else:
                reason = "ru_unseen"
                st_n = st.get(iv, (False, 0.0))[1]
                st_nreads_unseen.append(st_n)
            # ru_unseen dominates (more upstream) when both present
            if case_reason != "ru_unseen":
                case_reason = reason
        if case_reason is None:
            continue
        sub[case_reason] += 1
        by_cat[a["cat"]][case_reason] += 1
        rows.append({"chrom": chrom, "tid": tid, "cat": a["cat"], "reason": case_reason,
                     "n_missing": len(missing)})

    n = len(rows)
    print(f"=== graph_missing anatomy (n={n}) ===\n")
    for k, v in sub.most_common():
        print(f"  {k:12s} {v:4d}  ({100*v/n:.0f}%)")
    print("\nby chain-relationship category:")
    for cat in sorted(by_cat):
        cc = by_cat[cat]
        print(f"  {cat:22s} rejected={cc['ru_rejected']:3d}  unseen={cc['ru_unseen']:3d}")
    if st_nreads_unseen:
        s = sorted(st_nreads_unseen)
        med = s[len(s)//2]
        hi = sum(1 for x in s if x >= 10)
        print(f"\nru_unseen junctions — ST's read support (n={len(s)}):")
        print(f"  median ST nreads={med:.0f}; with ST nreads>=10: {hi} "
              f"({100*hi/len(s):.0f}%)  <- rustle saw 0, ST saw >=10 = alignment/multimap divergence")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/graph_missing_anatomy.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
