#!/usr/bin/env python3
"""PREREG_o2_per_family_2026-10-19: O2 (version 2) family by family in the near-identical regime.

    o2perfamily.py qualify AVA.paf FAMILIES.json SET OUT.json      the families whose members all align at [0.90, 1.00)
    o2perfamily.py table QUALIFY.json OUT.tsv                      every qualifying family's result in each data set
    o2perfamily.py fair O2.tsv TRUTH.tsv FRAG.bam QUALIFY.json OUT.json   the held-out verdict against minimap2 at MAPQ >= 1

A family works in a data set with >= MIN_PLACED placements, <= MAX_WRONG of them wrong, its sets holding the true copy for
>= MIN_SETS of its almost-tied reads, and O2's precision >= minimap2's primary's on those reads."""
import collections
import csv
import glob
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_ID, MIN_COVER, MIN_PLACED, MAX_WRONG, MIN_SETS = 0.90, 0.5, 5, 0.05, 0.95
T = "/mnt/linuxdisk/tmp/o3_graph"


def _o2psv():
    spec = importlib.util.spec_from_file_location("o3_graph_o2psv_pf", os.path.join(HERE, "o2psv.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pair_identities(paf_lines):
    """{(copy a, copy b): best identity} of same-family pairs whose alignment covers >= MIN_COVER of the shorter copy."""
    best = {}
    for ln in paf_lines:
        x = ln.split("\t")
        q, ql, qs, qe, t, tl, ts, te, m, bl = x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[6]), int(x[7]), int(x[8]), int(x[9]), int(x[10])
        if q == t or q.rsplit(".", 1)[0] != t.rsplit(".", 1)[0] or min(qe - qs, te - ts) < MIN_COVER * min(ql, tl):
            continue
        k = tuple(sorted((q, t)))
        best[k] = max(best.get(k, 0.0), m / bl)
    return best


def qualifies(n_copies, ids):
    """Every pair of the family's n_copies members aligned, at identity in [MIN_ID, 1)."""
    return len(ids) == n_copies * (n_copies - 1) // 2 and n_copies >= 2 and min(ids) >= MIN_ID and max(ids) < 1.0


def works(n):
    """A family's verdict in one data set from its score counts (o2psv.score)."""
    placed = n["right"] + n["wrong"]
    mm = n["minimap2 right"] + n["minimap2 wrong"]
    prec = n["right"] / placed if placed else 0.0
    return (placed >= MIN_PLACED and n["wrong"] <= MAX_WRONG * placed and n["judged sets"] > 0
            and n["set holds truth"] >= MIN_SETS * n["judged sets"] and prec >= (n["minimap2 right"] / mm if mm else 0.0))


def cmd_qualify(paf_p, fam_p, set_name, out_p):
    fams = json.load(open(fam_p))["sets"][set_name]
    ids = collections.defaultdict(list)
    for (a, b), v in pair_identities(open(paf_p)).items():
        ids[a.rsplit(".", 1)[0]].append(v)
    out = {f["family"]: dict(copies=len(f["copies"]), min_identity=min(ids[f["family"]]), max_identity=max(ids[f["family"]]))
           for f in fams if qualifies(len(f["copies"]), ids[f["family"]])}
    json.dump(out, open(out_p, "w"), indent=1)
    print(f"[o2perfamily] {len(out)} of {len(fams)} families qualify -> {out_p}", flush=True)
    return 0


def _by_family(rows):
    g = collections.defaultdict(list)
    for r in rows:
        g[r["family"]].append(r)
    return g


def cmd_table(qual_p, out_p):
    P = _o2psv()
    qual = json.load(open(qual_p))
    sim = collections.defaultdict(collections.Counter)                     # simulation: every CHM13-only world, rounds 23-25
    for t in sorted(glob.glob(f"{T}/fast/heldout2[345]/*_s1[345]0[123]/*/o2psv_v2.tsv")):
        w = t.rsplit("/", 1)[0]
        truth = {r["read"]: r["copy"] for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
        for f, rows in _by_family(csv.DictReader(open(t), delimiter="\t")).items():
            sim[f].update(P.score(rows, truth))
    hg = _by_family(csv.DictReader(open(f"{T}/o2psv_dev/hg002.v2.tsv"), delimiter="\t"))
    hg_truth = P._truth([f"{T}/hg002proof/truth.tsv", f"{T}/hg002proof/truth_untied.tsv"], True)
    a1_p = f"{T}/o2frag_a119b/o2psv_v2final.tsv"
    a1 = _by_family(csv.DictReader(open(a1_p), delimiter="\t")) if os.path.exists(a1_p) else {}
    a1_truth = P._truth([f"{T}/o2frag_a119b/truth.tsv"], False)
    keys = ("reads", "right", "wrong", "set holds truth", "judged sets", "minimap2 right", "minimap2 wrong")
    with open(out_p, "w") as fh:
        fh.write("family\tcopies\tmin_identity\tmax_identity\t" + "\t".join(f"{d}_{k.replace(' ', '_')}" for d in ("sim", "hg002", "a119b")
                                                                         for k in keys + ("works",)) + "\n")
        count = collections.Counter()
        for f in sorted(qual):
            cells = []
            for d, n in (("sim", sim.get(f)), ("hg002", P.score(hg.get(f, []), hg_truth, donor=True)),
                         ("a119b", P.score(a1.get(f, []), a1_truth, among_candidates=True))):
                n = n or collections.Counter()
                ok = works(n) if n.get("reads") else False
                count[(d, ok)] += 1
                cells += [str(n.get(k, 0)) for k in keys] + [str(ok)]
            fh.write(f"{f}\t{qual[f]['copies']}\t{qual[f]['min_identity']:.4f}\t{qual[f]['max_identity']:.4f}\t" + "\t".join(cells) + "\n")
    print(f"[o2perfamily] works in simulation: {count[('sim', True)]}; on HG002: {count[('hg002', True)]}; on A119b fragments: "
          f"{count[('a119b', True)]} (of {len(qual)} qualifying families) -> {out_p}", flush=True)
    return 0


FAIR_MIN_PLACED, FAIR_MAX_WRONG, FAIR_GAIN, FAIR_FAMILY_MIN = 100, 0.02, 0.05, 5


def fair_counts(rows, truth, mapq, families):
    """{family: {o2 right, o2 wrong, mm2 right, mm2 wrong}} counted by SOURCE read (a fragment's name is read|f<i>|...): a
    source is right for a method when every fragment of it the method places is on its truth copy, wrong when one is not.
    O2 places what it calls; minimap2's baseline places a fragment's primary only at MAPQ >= 1 (audit of register 1340: the
    fair, equally selective baseline). Fragments outside `families`, or whose truth is none of their tied copies, are left
    out."""
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        t = truth.get(r["read"])
        if t is None or r["family"] not in families or t not in r["candidates"].split(","):
            continue
        src = (r["family"], r["read"].split("|")[0])
        if r["call"] != "-":
            per[src]["o2"].append(r["call"] == t)
        if mapq.get(r["read"], 0) >= 1:
            per[src]["mm2"].append(r["primary"] == t)
    out = collections.defaultdict(lambda: dict.fromkeys(("o2 right", "o2 wrong", "mm2 right", "mm2 wrong"), 0))
    for (f, _), d in per.items():
        for m in ("o2", "mm2"):
            if d.get(m):
                out[f][f"{m} {'right' if all(d[m]) else 'wrong'}"] += 1
        out[f]                                                      # a family with fragments but no placement still shows
    return dict(out)


def cmd_fair(o2_p, truth_p, bam_p, qual_p, out_p):
    """PREREG_o2_fair_heldout_2026-10-19's verdict on one fragment set."""
    import pysam
    rows = list(csv.DictReader(open(o2_p), delimiter="\t"))
    truth = {r["read"]: r["copy"] for r in csv.DictReader(open(truth_p), delimiter="\t")}
    names = {r["read"] for r in rows}
    mapq = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in names and not rec.is_secondary and not rec.is_supplementary:
                mapq[rec.query_name] = rec.mapping_quality
    fam = fair_counts(rows, truth, mapq, set(json.load(open(qual_p))))
    tot = collections.Counter()
    for v in fam.values():
        tot.update(v)
    o2 = tot["o2 right"] + tot["o2 wrong"]
    mm = tot["mm2 right"] + tot["mm2 wrong"]
    p_o2, p_mm = tot["o2 right"] / max(1, o2), tot["mm2 right"] / max(1, mm)
    good = sorted(f for f, v in fam.items() if v["o2 right"] >= FAIR_FAMILY_MIN and v["o2 wrong"] == 0)
    bars = dict(F1=o2 >= FAIR_MIN_PLACED, F2=tot["o2 wrong"] <= FAIR_MAX_WRONG * o2, F3=p_o2 - p_mm >= FAIR_GAIN, F4=len(good) > 10)
    verdict = "UNDERPOWERED" if not bars["F1"] else ("PASS" if all(bars.values()) else "FAIL")
    json.dump(dict(families=fam, pooled=dict(tot), o2_precision=p_o2, mm2_precision=p_mm, families_working=good, bars=bars,
                   verdict=verdict), open(out_p, "w"), indent=1)
    print(f"[o2perfamily] fair: O2 placed {o2} source reads ({tot['o2 wrong']} wrong, precision {p_o2:.4f}); minimap2 MAPQ>=1 "
          f"placed {mm} ({tot['mm2 wrong']} wrong, precision {p_mm:.4f}); families with >= {FAIR_FAMILY_MIN} placed and none wrong: "
          f"{len(good)}; bars {bars}; verdict {verdict}", flush=True)
    return 0


def main(argv):
    if argv[:1] == ["qualify"] and len(argv) == 5:
        return cmd_qualify(*argv[1:])
    if argv[:1] == ["fair"] and len(argv) == 6:
        return cmd_fair(*argv[1:])
    if argv[:1] == ["table"] and len(argv) == 3:
        return cmd_table(*argv[1:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
