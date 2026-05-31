#!/usr/bin/env python3
"""Prize-bound for color/segmentation convergence: net reachable FP+FN at the
boundary-mismatch bundles. Abort gate: if net < 5, shelve Phase 1.
Usage: python3 bench/segmentation_prize.py /tmp/ru_gn.jsonl /tmp/st_gn.jsonl \
           /tmp/ru_final.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import json, sys, re
from collections import defaultdict

def load_gn(path):
    bundles = []
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get("step") != "graphnode_list":
                continue
            nodes = tuple(tuple(map(int, m.split("-")))
                          for m in re.findall(r"\d+-\d+", r.get("payload", {}).get("nodes", "")))
            bundles.append((r["start"], r["end"], r["strand"], nodes))
    return bundles

def load_gtf_chains(path):
    """transcript_id -> (strand, (start,end), intron_chain) ; intron_chain=() for single-exon."""
    tx = defaultdict(list)
    strand = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            tid = m.group(1)
            tx[tid].append((int(f[3]), int(f[4])))
            strand[tid] = f[6]
    chains = {}
    for tid, ex in tx.items():
        ex.sort()
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        chains[tid] = (strand[tid], (ex[0][0], ex[-1][1]), introns)
    return chains

def overlaps(a0, a1, b0, b1):
    return a0 <= b1 and b0 <= a1

def main():
    ru_gn = load_gn(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_gn.jsonl")
    st_gn = load_gn(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_gn.jsonl")
    ru_gtf = load_gtf_chains(sys.argv[3] if len(sys.argv) > 3 else "/tmp/ru_final.gtf")
    ref = load_gtf_chains(sys.argv[4] if len(sys.argv) > 4 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    ref_chains = {(s, ic) for (s, _sp, ic) in ref.values() if ic}

    # Find ST envelopes split by RU (>=2 same-strand RU bundles spanning one ST bundle).
    splits = []
    for (s0, s1, sstrand, snodes) in st_gn:
        ru_inside = [(r0, r1) for (r0, r1, rstrand, _rn) in ru_gn
                     if rstrand == sstrand and r0 >= s0 - 1 and r1 <= s1 + 1]
        # true split: >=2 RU bundles, none spanning the whole ST envelope
        if len(ru_inside) >= 2 and not any(r0 <= s0 + 1 and r1 >= s1 - 1 for (r0, r1) in ru_inside):
            splits.append((s0, s1, sstrand))

    # Attribute FP: RU final chains inside a split envelope that are NOT in the reference.
    fp = 0
    fp_examples = []
    for tid, (strand, (t0, t1), ic) in ru_gtf.items():
        if not ic:
            continue
        if (strand, ic) in ref_chains:
            continue  # a TP, not attributable
        for (s0, s1, sstrand) in splits:
            if sstrand == strand and overlaps(t0, t1, s0, s1):
                fp += 1
                fp_examples.append((tid, strand, t0, t1))
                break

    # FN: reference chains overlapping a split envelope that RU does NOT produce.
    ru_chains = {(s, ic) for (s, _sp, ic) in ru_gtf.values() if ic}
    fn = 0
    for (s, (r0, r1), ic) in ref.values():
        if not ic or (s, ic) in ru_chains:
            continue
        for (s0, s1, sstrand) in splits:
            if sstrand == s and overlaps(r0, r1, s0, s1):
                fn += 1
                break

    print(f"True-split ST envelopes (RU fragments): {len(splits)}")
    print(f"Attributable FP (RU non-ref chains in split envelopes): {fp}")
    print(f"Attributable FN (ref chains RU misses in split envelopes): {fn}")
    print(f"NET reachable prize (FP+FN): {fp + fn}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if fp + fn < 5 else 'PROCEED to Phase 1'}")
    for tid, strand, t0, t1 in fp_examples[:10]:
        print(f"    FP {tid} {strand} {t0}-{t1}")

if __name__ == "__main__":
    main()
