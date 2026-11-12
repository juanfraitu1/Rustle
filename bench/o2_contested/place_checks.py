#!/usr/bin/env python3
"""r1647 post hoc checks on a picks file. Usage: place_checks.py PICKS COPIES.tsv MOLECULES.tsv HITS.tsv ASSIGN.tsv[,ASSIGN2] DUMP [DUMP2]
Prints: discordant molecule list (tool right, primary wrong) with primary position for an overlap recheck; split by
haplotype of the read's best hit; negatives by whether the pick's span holds the primary; tool vs O2 off the shared
set; the 10-AS rule (as_margin > 10 -> the primary)."""
import collections, csv, sys

def main(picks, copies, molecules, hits, assigns, *dumps):
    span = {}
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    ov = lambda a, b: span[a][0] == span[b][0] and span[a][1] < span[b][2] and span[b][1] < span[a][2]
    pick = {}
    for l in open(picks): q, f, i = l.rstrip("\n").split("\t"); pick[q] = (f, i)
    mol = {}
    for d in dumps:
        for l in open(d):
            q, fam, st, v, al, _, t = l.rstrip("\n").split("\t")
            m = mol.setdefault(q, dict(al=al, o2="abstain", T=None if t == "OUTSIDE" else {tuple(x.split(":")) for x in t.split(";")}, fams=set()))
            m["fams"].add(fam)
            if st == "assigned": m["o2"] = v if m["o2"] != "wrong" else "wrong"
    prim = {}
    for l in open(molecules):
        if l[0] != "#":
            f = l.rstrip("\n").split("\t")
            if f[0] in mol and len(f) >= 16 and f[4] != "NA": prim[f[0]] = (f[4], int(f[15]))
    hap = collections.defaultdict(dict)
    with open(hits) as fh:
        next(fh)
        for l in fh:
            q, h, ch, s, e, nm, al = l.rstrip("\n").split("\t")
            if q[:2] == "R|" and q[2:] in mol: hap[q[2:]][h] = min(hap[q[2:]].get(h, 10 ** 9), int(nm))
    def origin(q):
        h = hap.get(q, {}); b = min(h.values()); w = sorted(k for k in h if h[k] == b)
        return w[0] if len(w) == 1 else "both"
    right = lambda q: any(pick[q] == t or ov(pick[q], t) for t in mol[q]["T"])
    c = collections.Counter(); out = open(picks + ".discordant.tsv", "w")
    for q, m in mol.items():
        if not m["T"]:
            c["neg"] += 1
            if q in pick:
                c["neg assigned"] += 1; p = prim.get(q); s = span[pick[q]]
                c["neg assigned: pick span holds the primary"] += bool(p and p[0] == s[0] and s[1] <= p[1] <= s[2])
            continue
        o = origin(q); c[f"decidable origin {o}"] += 1
        if q in pick:
            r = right(q); al = m["al"] == "al_right"; c[f"assigned origin {o}"] += 1
            c[f"origin {o}: tool {'right' if r else 'wrong'}"] += 1; c[f"origin {o}: primary {'right' if al else 'wrong'}"] += 1
            if r and not al:
                c[f"origin {o}: tool+prim-"] += 1; t = sorted(m["T"])[0]; p = prim.get(q, ("NA", -1))
                out.write(f"{q}\t{p[0]}\t{p[1]}\t{span[t][0]}\t{span[t][1]}\t{span[t][2]}\n")
            if al and not r: c[f"origin {o}: tool-prim+"] += 1
            if m["o2"] == "abstain": c["tool assigned, O2 abstained: " + ("right" if r else "wrong")] += 1
        elif m["o2"] != "abstain": c["O2 assigned, tool abstained: " + m["o2"]] += 1
    asm = {}
    for a in assigns.split(","):
        for r in csv.DictReader(open(a), delimiter="\t"):
            if r["read_name"] in mol and r["as_margin"] not in ("", "NA"): asm[r["read_name"]] = max(asm.get(r["read_name"], -1), int(r["as_margin"]))
    for q, m in mol.items():
        if m["T"] and asm.get(q, -1) > 10: c["10-AS rule decides"] += 1; c["10-AS rule " + ("right" if m["al"] == "al_right" else "wrong")] += 1
        elif m["T"]: c["10-AS rule abstains"] += 1
    for k in sorted(c): print(k, c[k], sep="\t")

if __name__ == "__main__": main(*sys.argv[1:])
