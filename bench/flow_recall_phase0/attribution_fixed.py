"""CORRECTED non-generation attribution (fixes the junction_accept coordinate convention).

junction_accept logs (start,end) = (intron_first-1, intron_last+1) (exon-flanking), NOT the
intron-interior coords ref_chain emits. The original attribute.py compared raw -> everything
fell to graph_missing. Convert (s,e)->(s+1,e-1) to intron-interior before matching.

For each NON-generated miss (gen_census verdict != generated), split:
  graph_missing    : a ref intron is not an accepted junction in rustle (after coord fix)
  flow_enumeration : all ref introns ARE accepted junctions, flow just didn't walk them
And for graph_missing, sub-split whether rustle SAW the junction (rejected) vs never saw it.
"""
import json, os, sys
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def accepted_introns(path):
    """Set of intron-interior (a,b) that rustle ACCEPTED, + set it merely SAW (any accepted flag)."""
    acc, seen = set(), set()
    for j in lib.read_parity(path, "junction_accept"):
        a, b = j["start"] + 1, j["end"] - 1          # exon-flank -> intron-interior
        seen.add((a, b))
        if bool(j.get("payload", {}).get("accepted", False)):
            acc.add((a, b))
    return acc, seen

def main():
    census = {}
    for line in open(f"{lib.ROOT}/bench/flow_recall_phase0/gen_census.jsonl"):
        r = json.loads(line)
        census[(r["chrom"], r["tid"])] = r["verdict"]
    nongen = [k for k, v in census.items() if v == "non_generated"]
    cause = Counter()
    gm_sub = Counter()
    rows = []
    for chrom, tid in nongen:
        wd = f"{lib.CACHE}/{chrom}/{tid}"
        ch = lib.ref_chain(chrom, tid)
        if ch is None:
            continue
        ref_introns = ch[2]
        acc, seen = accepted_introns(f"{wd}/r.jsonl")
        missing = [iv for iv in ref_introns if iv not in acc]
        if not missing:
            cause["flow_enumeration"] += 1
            rows.append({"chrom": chrom, "tid": tid, "cause": "flow_enumeration"})
        else:
            cause["graph_missing"] += 1
            # did rustle SEE any missing junction (rejected) or none (unseen)?
            sub = "ru_rejected" if any(iv in seen for iv in missing) else "ru_unseen"
            gm_sub[sub] += 1
            rows.append({"chrom": chrom, "tid": tid, "cause": "graph_missing", "gm_sub": sub})

    n = sum(cause.values())
    print(f"=== CORRECTED non-generation attribution (n={n}) ===")
    for k, v in cause.most_common():
        print(f"  {k:16s} {v:4d}  ({100*v/n:.0f}%)")
    print(f"\ngraph_missing sub-split (saw-but-rejected vs never-seen):")
    for k, v in gm_sub.most_common():
        print(f"  {k:12s} {v}")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/attribution_fixed.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
