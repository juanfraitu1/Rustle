"""Analysis (b): recall ceiling (FSM union) + precision cost of adopting ST's missing chains."""
import json, sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def fsm_union_per_chrom():
    agg = {}
    for line in open(f"{lib.RG}/results.jsonl"):
        r = json.loads(line)
        if not r.get("ok"):
            continue
        both, ro, so = r["both"], r["rustle_only"], r["st_only"]
        agg[r["chrom"]] = {
            "both": both, "rustle_only": ro, "st_only": so,
            "rustle_fsm": both + ro,
            "union_fsm": both + ro + so,
        }
    return agg

def annotated_chain_sets():
    """Per chrom: set of annotated intron chains (for precision-cost test)."""
    out = {}
    for line in open(f"{lib.RG}/results.jsonl"):
        r = json.loads(line)
        if not r.get("ok"):
            continue
        chrom = r["chrom"]
        ref = lib.parse_gtf(f"{lib.PERCHROM}/{chrom}/ref.gtf")
        out[chrom] = {d["introns"] for d in ref.values() if d["introns"]}
    return out

def precision_cost():
    """ST path_extracted chains at miss loci that are NOT annotated -> would-be FPs added."""
    misses = [json.loads(l) for l in open(f"{lib.RG}/st_only_classified.jsonl")]
    ann = annotated_chain_sets()
    extra = 0
    seen = set()
    for m in misses:
        wd = f"{lib.CACHE}/{m['chrom']}/{m['tid']}"
        if not os.path.exists(f"{wd}/.done"):
            continue
        rustle = {lib.parse_chain_str(e["payload"]["introns"])
                  for e in lib.read_parity(f"{wd}/r.jsonl", "path_extracted")
                  if e["payload"].get("introns")}
        for e in lib.read_parity(f"{wd}/s.jsonl", "path_extracted"):
            s = e["payload"].get("introns")
            if not s:
                continue
            ch = lib.parse_chain_str(s)
            key = (m["chrom"], ch)
            if key in seen:
                continue
            seen.add(key)
            if ch not in rustle and ch not in ann.get(m["chrom"], set()):
                extra += 1
    return extra

def main():
    agg = fsm_union_per_chrom()
    tot_rustle = sum(a["rustle_fsm"] for a in agg.values())
    tot_union = sum(a["union_fsm"] for a in agg.values())
    tot_so = sum(a["st_only"] for a in agg.values())
    print("=== recall-oracle ceiling (FSM union) ===")
    print(f"  current rustle FSM (genome-wide): {tot_rustle}")
    print(f"  ceiling if adopt ST misses:       {tot_union}  (+{tot_so})")
    fps = precision_cost()
    print(f"\n=== precision cost ===")
    print(f"  ST-extracted non-annotated chains at miss loci (would-be added FPs): {fps}")
    print(f"  crude recall:FP ratio = {tot_so}:{fps}"
          + (f"  ({tot_so/fps:.2f})" if fps else "  (no FPs)"))

if __name__ == "__main__":
    main()
