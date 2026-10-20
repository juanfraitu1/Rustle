#!/usr/bin/env python3
"""PREREG_o2_second_round_2026-10-20 (user: "it needs to make loci where the corrected multimappers now have a say ... a
post-assembly second round"): per family, how each arm's assembled loci share the family's reads among its copies, against
the donor's truth (HG002's assembly) -- the total variation distance (TVD) between the arm's per-copy shares and the truth's.

    o2second.py score COPIES.tsv FAMILIES.json OUT.json TRUTH.tsv ... -- NAME=GTF:WEIGHT ...
WEIGHT: reads (the assembler's `reads` attribute), cov (StringTie's `cov`), flair (FLAIR's isoform.read.map.txt beside the GTF).
The truth share of a copy: the family's reads whose donor-truth class is unique / unique_elsewhere with that copy.
"""
import collections
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = ("LOOP-O2", "P", "GOOD", "ST", "FL")                     # LOOP-O2 is compared with each of the others


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_second", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tvd(a, b):
    """Total variation distance between two {copy: weight} shares; 1.0 when either is empty."""
    sa, sb = sum(a.values()), sum(b.values())
    if sa <= 0 or sb <= 0:
        return 1.0
    return 0.5 * sum(abs(a.get(k, 0) / sa - b.get(k, 0) / sb) for k in set(a) | set(b))


def copy_weights(tx, weight, cps, fam):
    """{copy: summed weight} of the transcripts lying on the family's copies (summary.attribute: strand, then exon overlap)."""
    S = _sibling("summary")
    out = collections.Counter()
    for t, v in tx.items():
        c = S.attribute(v["blocks"], v.get("strand", "+"), cps.get(v["chrom"], ()))
        if c in fam:
            out[c] += weight.get(t, 0.0)
    return dict(out)


def _weights(gtf_p, kind, tx):
    if kind == "flair":
        L = _sibling("o2loci")
        return L.read_map(gtf_p.replace(".isoforms.gtf", ".isoform.read.map.txt"))
    w = {}
    for ln in open(gtf_p):
        x = ln.split("\t")
        if len(x) > 8 and x[2] == "transcript":
            t = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
            m = re.search(rf'(?:^|;)\s*{"reads" if kind == "reads" else "cov"} "([^"]+)"', x[8])   # not matched_reads / longcov
            w[t] = float(m.group(1)) if m else 0.0
    return w


def cmd_score(copies_p, fam_p, out_p, truth_ps, arms):
    P, L = _sibling("o2psv"), _sibling("o2loci")
    cps = P.load_copies(copies_p)
    fams = set(json.load(open(fam_p)))
    members = collections.defaultdict(set)
    for v in cps.values():
        for c, *_ in v:
            if c.rsplit(".", 1)[0] in fams:
                members[c.rsplit(".", 1)[0]].add(c)
    truth = collections.defaultdict(collections.Counter)
    for p in truth_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["class"] in ("unique", "unique_elsewhere") and r["family"] in fams and r["truth"] in members[r["family"]]:
                truth[r["family"]][r["truth"]] += 1
    shares = {}
    for name, spec in arms.items():
        gtf_p, kind = spec.rsplit(":", 1)
        tx = L.read_gtf(gtf_p)
        strands = {}
        for ln in open(gtf_p):
            x = ln.split("\t")
            if len(x) > 8 and x[2] == "transcript":
                strands[re.search(r'transcript_id "([^"]+)"', x[8]).group(1)] = x[6]
        for t, v in tx.items():
            v["strand"] = strands.get(t, "+")
        w = _weights(gtf_p, kind, tx)
        shares[name] = {f: copy_weights(tx, w, cps, members[f]) for f in fams}
    res = {f: dict(truth_reads=sum(truth[f].values()), truth=dict(truth[f]),
                   **{a: dict(tvd=tvd(shares[a][f], truth[f]), share=shares[a][f]) for a in arms}) for f in sorted(fams)}
    scored = [f for f in res if res[f]["truth_reads"] > 0]
    wsum = sum(res[f]["truth_reads"] for f in scored)
    pooled = {a: sum(res[f][a]["tvd"] * res[f]["truth_reads"] for f in scored) / max(1, wsum) for a in arms}
    better = {a: (sum(res[f]["LOOP-O2"]["tvd"] < res[f][a]["tvd"] for f in scored),
                  sum(res[f]["LOOP-O2"]["tvd"] > res[f][a]["tvd"] for f in scored)) for a in arms if a != "LOOP-O2"}
    bars = dict(S1=all(pooled["LOOP-O2"] < pooled[a] for a in ("P", "GOOD", "ST", "FL") if a in pooled),
                S2=all(better[a][0] > better[a][1] for a in ("P", "GOOD") if a in better), S3=len(scored) >= 5)
    verdict = "UNDERPOWERED" if not bars["S3"] else ("PASS" if bars["S1"] and bars["S2"] else "FAIL")
    json.dump(dict(families=res, pooled=pooled, loop_better_worse=better, bars=bars, verdict=verdict), open(out_p, "w"), indent=1)
    print(f"[o2second] {len(scored)} families, {wsum} truth reads; pooled TVD " +
          ", ".join(f"{a} {pooled[a]:.3f}" for a in arms) + "; LOOP-O2 better / worse than " +
          ", ".join(f"{a} {b[0]}/{b[1]}" for a, b in better.items()) + f"; bars {bars}; verdict {verdict}", flush=True)
    for f in scored:
        print(f"    {f[3:]:12s} truth {res[f]['truth_reads']:5d}  " + "  ".join(f"{a} {res[f][a]['tvd']:.3f}" for a in arms), flush=True)
    return 0


def main(argv):
    if argv[:1] == ["score"] and "--" in argv:
        i = argv.index("--")
        head, arms = argv[1:i], dict(a.split("=", 1) for a in argv[i + 1:])
        if len(head) >= 4:
            return cmd_score(head[0], head[1], head[2], head[3:], arms)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
