#!/usr/bin/env python3
"""Gene ends per NPIP / TBC1D3 copy (human, A119b): the RefSeq gene, Soto's gene model (CAT v4, the genes Soto's families are made
of) and our best model at the copy, with each end judged by EXONS against the RefSeq gene (no bp threshold):

  longer   the model has k exon(s) beyond the RefSeq gene's first (5') or last (3') exon      -> +k
  shorter  the model starts (ends) past the RefSeq gene's first (last) exon, missing k exon(s) -> -k
  same     the model's end exon overlaps the RefSeq end exon; the bp offset is kept for display

Inputs are the 2026-09-29 copy-recovery instruments (scratch, see docs/COPY_RECOVERY_TOOLS_2026-09-29.md on main): the copy table and
truth GTF (RefSeq transcripts per copy, CHM13 v2.0), our scored models (models.hsa.ours.json, best = exact chain first, then the
best gffcompare class, then support) and their GTF; Soto's genes come from the meeting page's families.json (CAT v4, CHM13 v1.0),
moved to v2.0 by the per-chromosome offsets found by sequence (chr16 -5, chr17 -291, chr18 0).

    python3 bench/soto_m2/soto_m2_gene_ends.py --copies copies.hsa.tsv --truth truth.hsa.gtf --ours-json models.hsa.ours.json \
        --ours-gtf hsa.ours.gtf --cat families.json --out gene_ends.json
"""
import argparse
import collections
import csv
import json
import re

SHIFT = {"chr16": -5, "chr17": -291, "chr18": 0}   # CHM13 v1.0 -> v2.0, verified by SD98 region sequence (2026-09-30)
RANK = {c: i for i, c in enumerate("=ckjmneoixpu")}


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def gtf_exons(path, key):
    """key ('gene_id' or 'transcript_id') -> list of [start0, end] exons."""
    ex = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(key + r' "([^"]+)"', f[8])
        if m:
            ex[m.group(1)].append([int(f[3]) - 1, int(f[4])])
    return ex


def end_call(ref, lane, strand):
    """[(k, bp) for the 5' end, (k, bp) for the 3' end]; k > 0 longer by k exons, k < 0 shorter by k exons, 0 same exon (bp =
    how far the model end lies beyond the RefSeq end, positive = longer)."""
    ref, lane = merge(ref), merge(lane)
    def side(r_end, l_end, r_all, l_all, outward):
        # outward(a, b): True if interval a lies entirely beyond interval b on this side
        if outward(l_end, r_end):
            return sum(1 for x in l_all if outward(x, r_end)), None
        if outward(r_end, l_end):
            return -sum(1 for x in r_all if outward(x, l_end)), None
        return 0, None
    left = lambda a, b: a[1] <= b[0]     # a entirely left of b
    right = lambda a, b: a[0] >= b[1]    # a entirely right of b
    kl, _ = side(ref[0], lane[0], ref, lane, left)
    kr, _ = side(ref[-1], lane[-1], ref, lane, right)
    bl = ref[0][0] - lane[0][0]          # positive: lane extends further left
    br = lane[-1][1] - ref[-1][1]        # positive: lane extends further right
    L, R = (kl, bl if kl == 0 else None), (kr, br if kr == 0 else None)
    return [L, R] if strand == "+" else [R, L]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("copies", "truth", "ours_json", "ours_gtf", "cat", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)

    copies = list(csv.DictReader(open(a.copies), delimiter="\t"))
    truth = gtf_exons(a.truth, "gene_id")
    ours_ex = gtf_exons(a.ours_gtf, "transcript_id")
    scored = json.load(open(a.ours_json))["copies"]
    cat = [g for g in json.load(open(a.cat))["genes"] if g["c"] in SHIFT]
    out, skipped = [], []
    for c in copies:
        cid, chrom, strand = c["cid"], c["chrom"], c["strand"]
        if not truth.get(cid):
            skipped.append(c["name"])   # annotated without exon features (TBC1D3 pseudogenes): no gene ends to compare
            continue
        ref = merge(truth[cid])
        r_lo, r_hi = ref[0][0], ref[-1][1]
        row = dict(cid=cid, name=c["name"], family=c["family"], chrom=chrom, strand=strand, ref=ref,
                   readthrough=int(c.get("readthrough", 0) or 0))
        # Soto's gene: the CAT v4 gene on this strand with the most exonic overlap with the RefSeq gene (ties: closest extent)
        best, best_ov, best_ext = None, 0, None
        for g in cat:
            if g["c"] != chrom or g["st"] != strand:
                continue
            d = SHIFT[chrom]
            if g["e"] + d <= r_lo or g["s"] + d >= r_hi:
                continue
            gx = merge([[x + d, y + d] for x, y in g["x"]])
            ov = sum(max(0, min(y, q) - max(x, p)) for x, y in gx for p, q in ref)
            ext = abs(gx[0][0] - r_lo) + abs(gx[-1][1] - r_hi)   # ties on overlap go to the closest extent
            if ov > 0 and (ov > best_ov or (ov == best_ov and ext < best_ext)):
                best, best_ov, best_ext = (g, gx), ov, ext
        if best:
            g, gx = best
            row["soto"] = dict(name=g["n"], fams=[f for f in g["sf"] if not f.startswith("Unassigned")], exons=gx,
                               ends=end_call(ref, gx, strand))
        # our best model
        s = scored.get(cid)
        models = [m for m in (s["models"] if s else []) if m["tid"] in ours_ex]
        if models:
            m = min(models, key=lambda m: (not m["exact"], RANK.get(m["cls"], 99), -m["support"]))
            mx = merge(ours_ex[m["tid"]])
            row["ours"] = dict(tid=m["tid"], cls=m["cls"], exact=m["exact"], support=m["support"], exons=mx,
                               n_models=len(models), fused_with=s.get("fused_with", [])[:4], ends=end_call(ref, mx, strand))
        out.append(row)
    json.dump(out, open(a.out, "w"), separators=(",", ":"))
    n_s = sum("soto" in r for r in out)
    n_o = sum("ours" in r for r in out)
    print(f"{len(out)} copies ({len(skipped)} without annotated exons skipped: {', '.join(skipped)}); Soto gene found for {n_s}, "
          f"our model for {n_o}")
    for lane in ("soto", "ours"):
        cnt = collections.Counter()
        for r in out:
            if lane in r:
                for side, (k, _bp) in zip(("5'", "3'"), r[lane]["ends"]):
                    cnt[(side, "longer" if k > 0 else "shorter" if k < 0 else "same exon")] += 1
        print(lane, dict(sorted(cnt.items())))


if __name__ == "__main__":
    main()
