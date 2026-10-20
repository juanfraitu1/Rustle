#!/usr/bin/env python3
"""PREREG_o2_second_round_2026-10-20 (user: "it needs to make loci where the corrected multimappers now have a say ... a
post-assembly second round"): per family, how each arm's assembled loci share the family's reads among its copies, against
the donor's truth (HG002's assembly) -- the total variation distance (TVD) between the arm's per-copy shares and the truth's.

    o2second.py score COPIES.tsv FAMILIES.json OUT.json TRUTH.tsv ... -- NAME=GTF:WEIGHT ...
    o2second.py score-reads BAM COPIES.tsv FAMILIES.json DOMAIN.tsv OUT.json TRUTH.tsv ... -- NAME=GTF:RULE[:HOME.tsv] ...
        register 1345's audit, the fair version: ONE read-crediting rule for every arm and only the reads O2 can move. Each
        corrector-domain read with a donor truth contributes the records its arm's read pool admits (RULE: primary, good =
        primaries + secondaries with AS >= 0.98 x the read's best, home = a homed read's records at home, else good); a record
        is credited to the copy of the arm's transcript covering >= half its aligned bases. A read is right (credited to its
        true copy alone), wrong (to another copy only), echo (to its true copy and another) or uncredited.
WEIGHT: reads (the assembler's `reads` attribute), cov (StringTie's `cov`), flair (FLAIR's isoform.read.map.txt beside the GTF).
The truth share of a copy: the family's reads whose donor-truth class is unique / unique_elsewhere with that copy.
"""
import bisect
import collections
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = ("LOOP-O2", "P", "GOOD", "ST", "FL")                     # LOOP-O2 is compared with each of the others
GOOD_RATIO = 0.98                                               # the shipped seeding pool (rustle_pipeline.sh --seed-pool good)
NOWHERE, OFFCOPY = "<nowhere>", "<off copy>"


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


def transcript_index(tx):
    """{chrom: (starts, [(start, end, blocks, copy)], longest)} of transcripts {t: {chrom, blocks (sorted), copy}}; a
    transcript on no copy is OFFCOPY."""
    by = collections.defaultdict(list)
    for v in tx.values():
        if v["blocks"]:
            by[v["chrom"]].append((v["blocks"][0][0], v["blocks"][-1][1], v["blocks"], v.get("copy") or OFFCOPY))
    out = {}
    for ch, lst in by.items():
        lst.sort(key=lambda x: (x[0], x[1]))
        out[ch] = ([x[0] for x in lst], lst, max(x[1] - x[0] for x in lst))
    return out


def credit(chrom, blocks, idx):
    """The copy of the transcript covering most of an aligned record, if it covers >= half its aligned bases; else NOWHERE."""
    if chrom not in idx or not blocks:
        return NOWHERE
    st, lst, longest = idx[chrom]
    lo, hi = blocks[0][0], blocks[-1][1]
    aligned = sum(e - s for s, e in blocks)
    best, who = 0, NOWHERE
    for s, e, tb, c in lst[bisect.bisect_left(st, lo - longest):bisect.bisect_left(st, hi)]:
        if e <= lo:
            continue
        ov = sum(max(0, min(b1, t1) - max(b0, t0)) for b0, b1 in blocks for t0, t1 in tb)
        if ov > best:
            best, who = ov, c
    return who if best >= 0.5 * aligned else NOWHERE


def admitted(recs, rule, home=None):
    """The records [{chrom, blocks, AS, prim}] an arm's read pool admits for one read."""
    if rule == "primary":
        return [r for r in recs if r["prim"]]
    best = max(r["AS"] for r in recs)
    good = [r for r in recs if r["prim"] or r["AS"] >= GOOD_RATIO * best]
    if rule == "good" or home is None:
        return good
    ch, s, e = home
    return [r for r in recs if r["chrom"] == ch and any(b0 < e and b1 > s for b0, b1 in r["blocks"])]


def verdict(credited, truth):
    c = {x for x in credited if x != NOWHERE}
    if not c:
        return "uncredited"
    if c == {truth}:
        return "right"
    return "echo" if truth in c else "wrong"


def cmd_score_reads(bam_p, copies_p, fam_p, domain_p, out_p, truth_ps, arms):
    import pysam
    P, L, S = _sibling("o2psv"), _sibling("o2loci"), _sibling("summary")
    cps = P.load_copies(copies_p)
    fams = set(json.load(open(fam_p)))
    truth = P._truth(list(truth_ps), True)
    reads = {}
    for r in csv.DictReader(open(domain_p), delimiter="\t"):
        t = truth.get(r["read"])
        if r["family"] in fams and t and t[0] in ("unique", "unique_elsewhere") and t[1] in r["candidates"].split(","):
            reads[r["read"]] = dict(family=r["family"], truth=t[1], moved=r["corrector"] != r["minimap2"])
    recs = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in reads and not (rec.is_unmapped or rec.is_supplementary) and rec.has_tag("AS"):
                recs[rec.query_name].append(dict(chrom=rec.reference_name, blocks=rec.get_blocks(), AS=rec.get_tag("AS"),
                                                 prim=not rec.is_secondary))
    res = {}
    for name, spec in arms.items():
        gtf_p, rule, *home_p = spec.split(":")
        tx = L.read_gtf(gtf_p)
        for ln in open(gtf_p):
            x = ln.split("\t")
            if len(x) > 8 and x[2] == "transcript":
                t = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
                if t in tx:
                    tx[t]["strand"] = x[6]
        for v in tx.values():
            v["copy"] = S.attribute(v["blocks"], v.get("strand", "+"), cps.get(v["chrom"], ())) if v["blocks"] else None
        idx = transcript_index(tx)
        homes = {}
        if home_p:
            for r in csv.DictReader(open(home_p[0]), delimiter="\t"):
                if r["read_name"] in reads:
                    homes[r["read_name"]] = (r["chrom"], int(r["start"]), int(r["end"]))
        n, fam_n, moved_n = collections.Counter(), collections.defaultdict(collections.Counter), collections.Counter()
        for q, info in reads.items():
            if q not in recs:
                continue
            got = {credit(r["chrom"], r["blocks"], idx) for r in admitted(recs[q], rule, homes.get(q) if rule == "home" else None)}
            v = verdict(got, info["truth"])
            n[v] += 1
            fam_n[info["family"]][v] += 1
            if info["moved"]:
                moved_n[v] += 1
        res[name] = dict(pooled=dict(n), families={f: dict(c) for f, c in fam_n.items()}, o2_moved=dict(moved_n))
    k = max(1, len([q for q in reads if q in recs]))
    rate = {a: {v: res[a]["pooled"].get(v, 0) / k for v in ("right", "wrong", "echo", "uncredited")} for a in arms}
    others = [a for a in ("P", "GOOD", "ST", "FL") if a in rate]
    gate = "ORACLE" in rate and rate["ORACLE"]["right"] >= max(rate[a]["right"] for a in others) + 0.05
    bars = dict(G=gate, B1=all(rate["LOOP-O2"]["right"] > rate[a]["right"] for a in others),
                B2=("P" in rate and rate["LOOP-O2"]["wrong"] + rate["LOOP-O2"]["echo"] <= rate["P"]["wrong"] + rate["P"]["echo"]),
                B3=k >= 100)
    verdict_ = ("UNDERPOWERED" if not bars["B3"] else "INFEASIBLE" if not bars["G"] else
                "PASS" if bars["B1"] and bars["B2"] else "FAIL")
    json.dump(dict(reads=k, rates=rate, arms=res, bars=bars, verdict=verdict_), open(out_p, "w"), indent=1)
    print(f"[o2second] {k} corrector-domain reads with a donor truth in {len(fams)} families; verdict {verdict_}; bars {bars}")
    for a in arms:
        r = rate[a]
        print(f"    {a:8s} right {r['right']:.3f}  wrong {r['wrong']:.3f}  echo {r['echo']:.3f}  uncredited {r['uncredited']:.3f}"
              f"   (reads O2 moved: {res[a]['o2_moved']})", flush=True)
    return 0


def main(argv):
    if argv[:1] == ["score-reads"] and "--" in argv:
        i = argv.index("--")
        head, arms = argv[1:i], dict(a.split("=", 1) for a in argv[i + 1:])
        if len(head) >= 6:
            return cmd_score_reads(head[0], head[1], head[2], head[3], head[4], head[5:], arms)
    if argv[:1] == ["score"] and "--" in argv:
        i = argv.index("--")
        head, arms = argv[1:i], dict(a.split("=", 1) for a in argv[i + 1:])
        if len(head) >= 4:
            return cmd_score(head[0], head[1], head[2], head[3:], arms)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
