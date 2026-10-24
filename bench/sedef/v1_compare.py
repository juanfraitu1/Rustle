#!/usr/bin/env python3
"""V1 / V1b of docs/archive/2026-10/merged/PREREG_sedef_rs_layers_2026-10-13.md: sedef-rs SD calls (`utilities sd-scan` PREFIX.sd.pairs.tsv)
against SEDEF's (34-column BED), base-pair-wise.

    v1_compare.py --sedef REF.bed --sedefrs PREFIX.sd.pairs.tsv (--chain V1_TO_V2.chain | --fai GENOME.fai)
                  [--min-id 0.90] [--min-len 1000] --out OUT.json

Kept pairs: SEDEF identity = column 21 (fracMatch) >= --min-id; sedef-rs status `ok` and est_id / 100 >= --min-id; both
tools: min(side lengths) >= --min-len. Each tool's SD regions = both sides of its kept pairs, merged per chromosome and
intersected with the UNIVERSE: the v2.0 bases with a v1.0 ancestor in the chain (--chain: CHM13 v2.0 against lifted
v1.0 SEDEF) or every base of the genome (--fai). recall = SEDEF region bases covered by sedef-rs regions / SEDEF region
bases; precision = sedef-rs region bases covered by SEDEF regions / sedef-rs region bases. Descriptive: per chromosome,
and the recall of the SEDEF pairs binned by identity and by length (the bases of each pair's two sides inside the
universe that sedef-rs regions cover; overlapping pairs count their shared bases more than once)."""
import argparse
import bisect
import collections
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "o1_hier"))
from v1_to_v2 import read_chain  # noqa: E402

ID_BINS = [(0.90, 0.95), (0.95, 0.98), (0.98, 0.99), (0.99, 1.0000001)]
LEN_BINS = [(1_000, 5_000), (5_000, 20_000), (20_000, 100_000), (100_000, 10**12)]


def read_pairs(lines, fmt, min_id, min_len):
    """-> [((chrom, start, end), (chrom, start, end), identity)] kept by the filters; fmt 'sedef' or 'sedefrs'"""
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if fmt == "sedef":
            ident = float(f[20])
        else:
            if f[10] != "ok":
                continue
            ident = float(f[8]) / 100.0
        a, b = (f[0], int(f[1]), int(f[2])), (f[3], int(f[4]), int(f[5]))
        if ident >= min_id and min(a[2] - a[1], b[2] - b[1]) >= min_len:
            out.append((a, b, ident))
    return out


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def regions(pairs):
    by = collections.defaultdict(list)
    for a, b, _ in pairs:
        for c, s, e in (a, b):
            by[c].append((s, e))
    return {c: merge(v) for c, v in by.items()}


def intersect(x, y):
    """the bases in both merged, sorted interval lists, as a merged list"""
    out, i, j = [], 0, 0
    while i < len(x) and j < len(y):
        s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
        if s < e:
            out.append((s, e))
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return out


def total(iv):
    return sum(e - s for s, e in iv)


def covered_one(s, e, by, starts):
    """bases of [s, e) covered by the merged, sorted intervals `by` (`starts` = their starts)"""
    j = max(0, bisect.bisect_right(starts, s) - 1)
    tot = 0
    while j < len(by) and by[j][0] < e:
        tot += max(0, min(e, by[j][1]) - max(s, by[j][0]))
        j += 1
    return tot


def universe_from_chain(lines):
    """the v2 (query) side of every aligned block of the chain: {chrom: merged intervals}"""
    return {c: merge([(q, q + size) for _, q, size in blocks]) for c, blocks in read_chain(lines).items()}


def universe_from_fai(lines):
    return {f[0]: [(0, int(f[1]))] for f in (ln.split("\t") for ln in lines if ln.strip())}


def compare(ref_pairs, test_pairs, universe):
    """-> (result dict, the sedef-rs regions clipped to the universe)"""
    clip = lambda reg: {c: intersect(v, universe.get(c, [])) for c, v in reg.items()}  # noqa: E731
    r, t = clip(regions(ref_pairs)), clip(regions(test_pairs))
    per = {}
    for c in sorted(set(r) | set(t)):
        rv, tv = r.get(c, []), t.get(c, [])
        per[c] = dict(ref_bp=total(rv), test_bp=total(tv), ref_covered=total(intersect(rv, tv)),
                      test_covered=total(intersect(tv, rv)))
    rb, tb = sum(p["ref_bp"] for p in per.values()), sum(p["test_bp"] for p in per.values())
    rc, tc = sum(p["ref_covered"] for p in per.values()), sum(p["test_covered"] for p in per.values())
    return dict(recall=rc / rb if rb else None, precision=tc / tb if tb else None, ref_bp=rb, test_bp=tb,
                ref_covered=rc, test_covered=tc, per_chrom=per), t


def recall_bins(ref_pairs, test_regions, universe):
    """recall of the SEDEF pairs per identity bin and per length bin; `test_regions` are already inside the universe"""
    us = {c: [s for s, _ in v] for c, v in universe.items()}
    ts = {c: [s for s, _ in v] for c, v in test_regions.items()}
    acc = {"identity": collections.defaultdict(lambda: [0, 0]), "length": collections.defaultdict(lambda: [0, 0])}
    for a, b, ident in ref_pairs:
        length = min(a[2] - a[1], b[2] - b[1])
        for kind, value, bins in (("identity", ident, ID_BINS), ("length", length, LEN_BINS)):
            lab = next((f"{lo}-{hi}" for lo, hi in bins if lo <= value < hi), None)
            if lab is None:
                continue
            for c, s, e in (a, b):
                acc[kind][lab][0] += covered_one(s, e, universe.get(c, []), us.get(c, []))
                acc[kind][lab][1] += covered_one(s, e, test_regions.get(c, []), ts.get(c, []))
    return {k: {lab: dict(bp=v[0], recall=v[1] / v[0] if v[0] else None) for lab, v in sorted(d.items())}
            for k, d in acc.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sedef", required=True)
    ap.add_argument("--sedefrs", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--chain")
    g.add_argument("--fai")
    ap.add_argument("--min-id", type=float, default=0.90)
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    with open(a.sedef) as fh:
        ref = read_pairs(fh, "sedef", a.min_id, a.min_len)
    with open(a.sedefrs) as fh:
        test = read_pairs(fh, "sedefrs", a.min_id, a.min_len)
    with open(a.chain or a.fai) as fh:
        universe = universe_from_chain(fh) if a.chain else universe_from_fai(fh)
    res, test_regions = compare(ref, test, universe)
    res.update(inputs=vars(a), ref_pairs=len(ref), test_pairs=len(test),
               universe_bp=sum(total(v) for v in universe.values()),
               recall_bins=recall_bins(ref, test_regions, universe))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("recall", "precision", "ref_bp", "test_bp", "ref_pairs", "test_pairs")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
