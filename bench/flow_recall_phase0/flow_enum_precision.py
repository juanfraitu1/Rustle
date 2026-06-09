"""Generation precision-cost of the parse_trflong lever, scoped to the flow_enumeration loci.

flow_enumeration = rustle's graph has the junctions, ST extracts the chain, rustle's flow
doesn't. Completing parse_trflong would extract more chains at these loci. Proxy for what it
would generate = ST's path_extracted at these loci (the prior work's oracle). At each
flow_enumeration locus, of the chains ST extracts that rustle does NOT:
  recovery (real) = chain is annotated  -> a true recovery
  FP (spurious)   = chain is NOT annotated -> collateral a port would likely add
Reports the recall:FP ratio -> net-F1 direction for the parse_trflong lever specifically.
Compare to the all-668 recall-oracle ratio of 0.39.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def main():
    # flow_enumeration loci from the FIXED attribution
    attr = [json.loads(l) for l in open(f"{lib.ROOT}/bench/flow_recall_phase0/attribution.jsonl")]
    fe = [(a["chrom"], a["tid"]) for a in attr if a["cause"] == "flow_enumeration"]
    # annotated intron-chain sets per chrom (reuse the same construction)
    ann_by_chrom = {}
    def ann(chrom):
        if chrom not in ann_by_chrom:
            ref = lib.parse_gtf(f"{lib.PERCHROM}/{chrom}/ref.gtf")
            ann_by_chrom[chrom] = {d["introns"] for d in ref.values() if d["introns"]}
        return ann_by_chrom[chrom]

    recovery = 0      # annotated chains ST extracts, rustle doesn't (real recoveries)
    fp = 0            # non-annotated chains ST extracts, rustle doesn't (collateral)
    target_hit = 0    # of the 190 targets, how many are themselves in ST's extraction
    seen = set()
    for chrom, tid in fe:
        wd = f"{lib.CACHE}/{chrom}/{tid}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        rustle = {lib.parse_chain_str(e["payload"]["introns"])
                  for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted")
                  if e["payload"].get("introns")}
        target = lib.ref_chain(chrom, tid)
        target_chain = tuple(target[2]) if target else None
        a = ann(chrom)
        for e in lib.read_parity(f"{wd}/s.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            key = (chrom, ch)
            if key in seen or ch in rustle:
                continue
            seen.add(key)
            if ch in a:
                recovery += 1
                if ch == target_chain:
                    target_hit += 1
            else:
                fp += 1

    print(f"=== parse_trflong generation precision-cost (flow_enumeration loci, n_targets={len(fe)}) ===")
    print(f"  real recoveries (annotated, ST-extracted, rustle-missed): {recovery}")
    print(f"  of which are the flow_enumeration targets themselves:     {target_hit}")
    print(f"  FP collateral (non-annotated, ST-extracted, rustle-missed): {fp}")
    ratio = recovery / fp if fp else float("inf")
    print(f"\n  recall:FP ratio = {recovery}:{fp}  ({ratio:.2f})")
    print(f"  (all-668 recall-oracle ratio was 0.39; >1.0 here would mean the lever is")
    print(f"   net-F1-favorable at these loci, <1.0 means collateral dominates)")

if __name__ == "__main__":
    main()
