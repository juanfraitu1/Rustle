"""Analysis (d): for non-generated misses, split graph_missing vs flow_enumeration."""
import json, sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def attribute_cause(ref_introns, rustle_junctions):
    for iv in ref_introns:
        if iv not in rustle_junctions:
            return "graph_missing"
    return "flow_enumeration"

def rustle_junctions(wd):
    # FIX (2026-06-09): junction_accept logs (start,end) in the exon-flanking convention
    # (intron_first-1, intron_last+1), NOT the intron-interior coords ref_chain emits. The
    # original raw comparison NEVER matched -> every multi-exon miss fell to graph_missing
    # (and only vacuous 0-intron single-exon cases matched -> false flow_enumeration). Emit
    # BOTH conventions per junction so the match is robust regardless of logging site.
    out = set()
    for j in lib.read_parity(f"{wd}/r.jsonl", "junction_accept"):
        s, e = j["start"], j["end"]
        out.add((s, e))
        out.add((s + 1, e - 1))
    return out

def main():
    census = {}
    gc = f"{lib.ROOT}/bench/flow_recall_phase0/gen_census.jsonl"
    if os.path.exists(gc):
        for line in open(gc):
            r = json.loads(line)
            census[(r["chrom"], r["tid"])] = r["verdict"]
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    cause = Counter()
    rows = []
    for m in misses:
        key = (m["chrom"], m["tid"])
        wd = f"{lib.CACHE}/{m['chrom']}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        if census.get(key) == "generated":
            continue   # only attribute the non-generated set
        ch = lib.ref_chain(m["chrom"], m["tid"])
        if ch is None:
            continue
        if not ch[2]:
            c = "single_exon"   # 0 introns: separate, not a vacuous flow_enumeration
        else:
            c = attribute_cause(ch[2], rustle_junctions(wd))
        cause[c] += 1
        rows.append({"chrom": m["chrom"], "tid": m["tid"], "cat": m["cat"], "cause": c})
    print("=== non-generation attribution ===")
    for k, v in cause.most_common():
        print(f"  {k:16s} {v}")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/attribution.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
