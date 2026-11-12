#!/usr/bin/env python3
"""r1657 follow-up: is a strict_truth "wrong" a truth artifact? For each decidable assigned wrong row of a dump, the
picked copy's state on the haplotypes: its sequence fails strict_truth's near-perfect gate on every haplotype (so truth
could never name it), and whether its own best hit (gate off, slack 0) overlaps the read's best hit.
Usage: truth_audit.py DUMP TRUTH_DIR [CORRECTED_DUMP: artifact rows rewritten as right]"""
import collections, sys

def main(dump, tdir, fixed=None):
    out = open(fixed, "w") if fixed else None
    hits = collections.defaultdict(list)
    fh = open(tdir + "/hits.tsv"); next(fh)
    for l in fh:
        q, hap, ch, s, e, nm, al = l.rstrip("\n").split("\t"); hits[q].append((hap, ch, int(s), int(e), int(nm), int(al)))
    def copy_state(fam, idx):
        hs = hits.get(f"C|{fam}|{idx}", [])
        if not hs: return "no_hit", []
        best = {}
        for h in hs: best[h[0]] = min(best.get(h[0], 10 ** 9), h[4])
        gated = all(best[h[0]] > max(2, 0.003 * h[5]) for h in hs)
        return ("gate_failed" if gated else "has_locus"), [h for h in hs if h[4] == best[h[0]]]
    c = collections.Counter(); fam = collections.defaultdict(collections.Counter); n = collections.Counter()
    for l in open(dump):
        q, f, st, v, al, pick, t = l.rstrip("\n").split("\t")
        if st != "assigned" or t == "OUTSIDE" or v != "wrong":
            if out: out.write(l)
            if st == "assigned" and t != "OUTSIDE": n[v] += 1
            continue
        n[v] += 1
        state, own = copy_state(f, pick)
        rh = hits["R|" + q]; b = min(h[4] for h in rh); near = [h for h in rh if h[4] == b]
        ov = any(o[0] == h[0] and o[1] == h[1] and o[2] < h[3] and h[2] < o[3] for o in own for h in near)
        k = (state, "pick's own best hit IS the read's best locus" if ov else "pick elsewhere")
        c[k] += 1; fam[k][f] += 1
        if out: out.write("\t".join([q, f, st, "right" if ov else v, al, pick, t]) + "\n")
    print(f"{dump}: right {n['right']} wrong {n['wrong']}")
    for k, x in c.most_common(): print(f"  {x}\t{k[0]}\t{k[1]}\t{fam[k].most_common(4)}")
    art = sum(x for k, x in c.items() if k[1].startswith("pick's"))
    print(f"  wrong rows that are truth artifacts (pick sits on the read's best haplotype locus): {art}; wrong rate {n['wrong'] / max(sum(n.values()), 1):.3f} -> {(n['wrong'] - art) / max(sum(n.values()), 1):.3f} if counted right")

if __name__ == "__main__": main(*sys.argv[1:4])
