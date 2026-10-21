#!/usr/bin/env python3
"""The vault's "statistics" for reads whose bases cannot pick a copy (goal 2026-10-20: "look into the ThesisVault for
pertinent methods that could be used for O2"), as assign-or-abstain rules over the copies O2's set keeps (o2window.py's
Pareto set, window 30):

- rescue (Deschamps 2020: weight a multiread by each candidate's uniquely mapped reads, "a single iteration of the
  expectation maximization algorithm"): the kept copy with the most uniquely placed reads; equal counts abstain.
- parsimony (Canzar 2016's facility location as a decision: a read may only go to an "open" facility): a kept copy is open
  when the sample has at least one uniquely placed read on it; the read is assigned when exactly one kept copy is open,
  else it abstains. No threshold.

Uniquely placed read: a primary alignment with MAPQ >= MIN_MAPQ on a copy (summary.attribute), not in the tied domain. A
fragment set's own source read is left out of the counts that judge its fragments (names SOURCE|f<i>|...).

    o2prior.py eval OUT.json --sample NAME CALLS TRUTH BAM COPIES [--fragments] [--sample ...]
        per sample (and pooled per substrate name prefix): the reads the set leaves with > 1 copy, determined, truth among
        the window's copies -> minimap2 right; each rule's decisions, right, and minimap2 right on the same reads.
"""
import collections
import csv
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_MAPQ = 20


def _sib(n):
    s = importlib.util.spec_from_file_location(f"o3_graph_{n}_prior", os.path.join(HERE, f"{n}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def rescue(kept, counts):
    """The kept copy with the most uniquely placed reads, or None when the most is shared or zero."""
    v = sorted(((counts.get(c, 0), c) for c in kept), reverse=True)
    if not v or v[0][0] == 0 or (len(v) > 1 and v[0][0] == v[1][0]):
        return None
    return v[0][1]


def parsimony(kept, counts):
    """The only kept copy with any uniquely placed read ("open"), or None when none or several are open."""
    open_ = [c for c in kept if counts.get(c, 0) > 0]
    return open_[0] if len(open_) == 1 else None


def equal_as(kept, best_as):
    """Whether every kept copy has an alignment and all share one best AS: the equal-score multimappers rescue and EM
    were designed for."""
    v = [best_as.get(c) for c in kept]
    return None not in v and len(set(v)) == 1


def own_truth(truth, tag):
    """A pooled-world truth table (keys '<world tag>/<read>', o2window.truth-sim) cut to one world and unkeyed; any other
    table (PacBio names carry '/') is returned as it is."""
    pre = tag + "/"
    if not any(k.startswith(pre) for k in truth):
        return truth
    return {k[len(pre):]: v for k, v in truth.items() if k.startswith(pre)}


def total(per_source):
    out = collections.Counter()
    for v in per_source.values():
        out.update(v)
    return dict(out)


def without(tot, own):
    return {c: n - own.get(c, 0) for c, n in tot.items()}


def unique_counts(bam_p, copies_p, domain, fragments=False):
    """{source: {copy: uniquely placed reads}} (source "" unless fragments): primaries with MAPQ >= MIN_MAPQ on a copy,
    read over the copies' spans, reads of the tied domain left out."""
    import pysam
    P, S, C = _sib("o2psv"), _sib("summary"), _sib("o2correct")
    cps = P.load_copies(copies_p)
    idx = C.span_index(cps)
    seen, out = set(), collections.defaultdict(collections.Counter)
    best_as = collections.defaultdict(dict)                     # domain reads: {read: {copy: best AS}}
    with pysam.AlignmentFile(bam_p) as bam:
        for ch, rows in cps.items():
            if ch not in bam.references:
                continue
            for _, _, _, (s, e) in rows:
                for rec in bam.fetch(ch, s, e):
                    q = rec.query_name
                    if rec.is_supplementary or rec.is_unmapped:
                        continue
                    if q in domain:
                        if rec.has_tag("AS"):
                            cp = S.attribute(rec.get_blocks(), P._strand(rec),
                                             C.overlapping(idx, ch, rec.reference_start, rec.reference_end))
                            if cp is not None:
                                best_as[q][cp] = max(best_as[q].get(cp, -10 ** 9), rec.get_tag("AS"))
                        continue
                    if q in seen or rec.is_secondary or rec.mapping_quality < MIN_MAPQ:
                        continue
                    cp = S.attribute(rec.get_blocks(), P._strand(rec),
                                     C.overlapping(idx, ch, rec.reference_start, rec.reference_end))
                    if cp is None:
                        continue
                    seen.add(q)
                    out[q.split("|")[0] if fragments else ""][cp] += 1
    return out, best_as


def evaluate(calls_p, truth, per_source, fragments=False, best_as=None):
    """Counts over the reads O2's set leaves with more than one copy (window 30), determined, truth among the window's
    copies; with best_as, only the equal-score ties among them (equal_as)."""
    tot = total(per_source)
    n = collections.Counter()
    for r in csv.DictReader(open(calls_p), delimiter="\t"):
        if r["window"] != "30":
            continue
        kept = set(r["o2_set"].split(","))
        st, t = truth.get(r["read"], ("missing", ""))
        if len(kept) < 2 or st != "determined" or t not in r["copies"].split(","):
            continue
        if best_as is not None and not equal_as(kept, best_as.get(r["read"], {})):
            continue
        counts = without(tot, per_source.get(r["read"].split("|")[0], {})) if fragments else tot
        mm = r["minimap2"] == t
        n["reads"] += 1
        n["truth kept"] += t in kept
        n["minimap2 right"] += mm
        for name, rule in (("rescue", rescue), ("parsimony", parsimony)):
            c = rule(kept, counts)
            if c is None:
                continue
            n[f"{name} decides"] += 1
            n[f"{name} right"] += c == t
            n[f"{name}: minimap2 right on those"] += mm
            if name == "rescue":                                # descriptive: precision by the top-to-second support ratio
                v = sorted((counts.get(x, 0) for x in kept), reverse=True)
                ratio = v[0] / max(1, v[1])
                b = ">= 8" if ratio >= 8 else ">= 4" if ratio >= 4 else ">= 2" if ratio >= 2 else "< 2"
                n[f"rescue ratio {b}: decides"] += 1
                n[f"rescue ratio {b}: right"] += c == t
                n[f"rescue ratio {b}: minimap2 right"] += mm
    return n


def cmd_eval(out_p, samples):
    res, pooled = {}, collections.defaultdict(collections.Counter)
    for name, calls_p, truth_p, bam_p, copies_p, frag in samples:
        truth = own_truth({r["read"]: (r["status"], r["truth"]) for r in csv.DictReader(open(truth_p), delimiter="\t")},
                          os.path.basename(calls_p)[:-len(".calls.tsv")])
        domain = {r["read"] for r in csv.DictReader(open(calls_p), delimiter="\t")}
        per_source, best_as = unique_counts(bam_p, copies_p, domain, frag)
        for label, ba in (("", None), (" [equal-AS ties]", best_as)):
            n = evaluate(calls_p, truth, per_source, frag, ba)
            res[name + label] = dict(n)
            pooled[name.split(":")[0] + label].update(n)
            print(f"[o2prior] {name}{label}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    for k, n in pooled.items():
        res[f"POOLED {k}"] = dict(n)
        print(f"[o2prior] POOLED {k}: " + ", ".join(f"{a} {b}" for a, b in sorted(n.items())), flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


def main(argv):
    if not argv or argv[0] != "eval":
        print(__doc__)
        return 2
    out_p, rest, samples = argv[1], argv[2:], []
    while rest:
        assert rest[0] == "--sample", rest[0]
        frag = len(rest) > 6 and rest[6] == "--fragments"
        samples.append((rest[1], rest[2], rest[3], rest[4], rest[5], frag))
        rest = rest[7 if frag else 6:]
    return cmd_eval(out_p, samples)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
