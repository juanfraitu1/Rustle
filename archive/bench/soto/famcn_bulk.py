#!/usr/bin/env python3
"""famCN for MANY intervals at once, by streaming whole WSSD tracks instead of per-interval queries.

`famcn_from_wssd.py` spawns one `bigBedToBed` per (interval, sample). That is fine for tens of intervals and
hopeless for the Soto replication, which needs famCN for every exon of every SD98 gene: 32k exons x N samples
is millions of subprocess calls. Here each sample's track is converted ONCE and swept against the sorted
interval list in memory, so the cost is (samples) conversions rather than (samples x intervals) queries.

famCN per gene = median over samples of the LENGTH-WEIGHTED mean CN across that gene's exons, matching how
Soto describe it ("read depth of multi-mapping reads with non-overlapping sliding windows").

COORDINATES: no liftover. The WSSD tracks and the CAT v4 annotation are both T2T-CHM13 **v1.0**, so intervals
are used as given. (`famcn_from_wssd.py` lifts v2.0->v1.0 because our own catalog is v2.0; applying that here
would shift every interval by a few hundred bp and silently corrupt the result.)
"""
import argparse, csv, os, statistics as st, subprocess, sys, tempfile
from collections import defaultdict


def load_intervals(path, id_col):
    """-> {chrom: sorted [(start, end, id)]}, plus the id order for output."""
    per = defaultdict(list)
    ids = []
    seen = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = r.get("chrom")
            if not c:
                continue
            gid = r[id_col]
            per[c].append((int(r["start"]), int(r["end"]), gid))
            if gid not in seen:
                seen.add(gid)
                ids.append(gid)
    for c in per:
        per[c].sort()
    return per, ids


def sweep(track_bed, per_chrom):
    """One pass over a sorted WSSD BED; returns {id: (weighted_sum, bases)}."""
    acc = defaultdict(lambda: [0.0, 0])
    cur_chrom, items, i = None, [], 0
    with open(track_bed) as fh:
        for ln in fh:
            f = ln.split("\t", 10)
            if len(f) < 10:
                continue
            c = f[0]
            if c != cur_chrom:
                cur_chrom, items, i = c, per_chrom.get(c, []), 0
                if not items:
                    continue
            if not items:
                continue
            try:
                ws, we, cn = int(f[1]), int(f[2]), float(f[9])
            except ValueError:
                continue
            # advance past intervals wholly left of this window
            while i < len(items) and items[i][1] <= ws:
                i += 1
            j = i
            while j < len(items) and items[j][0] < we:
                s, e, gid = items[j]
                ov = min(e, we) - max(s, ws)
                if ov > 0:
                    a = acc[gid]
                    a[0] += cn * ov
                    a[1] += ov
                j += 1
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervals", required=True, help="TSV with chrom/start/end + an id column")
    ap.add_argument("--wssd-dir", required=True)
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--id-col", default="gene_id")
    ap.add_argument("--tool", default="bigBedToBed")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    per, ids = load_intervals(a.intervals, a.id_col)
    n_iv = sum(len(v) for v in per.values())
    print(f"[famcn-bulk] {n_iv} intervals over {len(ids)} ids, {len(per)} chromosomes", file=sys.stderr)

    bbs = sorted(f for f in os.listdir(a.wssd_dir) if f.endswith("_wssd.bb"))[: a.samples]
    if not bbs:
        sys.exit(f"no *_wssd.bb in {a.wssd_dir}")
    print(f"[famcn-bulk] {len(bbs)} samples", file=sys.stderr)

    per_sample = defaultdict(list)   # id -> [famCN in each sample]
    with tempfile.TemporaryDirectory() as td:
        tmp = os.path.join(td, "track.bed")
        for k, bb in enumerate(bbs, 1):
            src = os.path.join(a.wssd_dir, bb)
            r = subprocess.run([a.tool, src, tmp], capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  [warn] {bb} failed to convert, skipping", file=sys.stderr)
                continue
            acc = sweep(tmp, per)
            for gid, (wsum, bases) in acc.items():
                if bases > 0:
                    per_sample[gid].append(wsum / bases)
            print(f"  [{k}/{len(bbs)}] {bb}: {len(acc)} ids covered", file=sys.stderr)
            os.remove(tmp)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "famCN", "famCN_mad", "n_samples"])
        n = 0
        for gid in ids:
            v = per_sample.get(gid)
            if not v:
                w.writerow([gid, "", "", 0])
                continue
            med = st.median(v)
            mad = st.median([abs(x - med) for x in v])
            w.writerow([gid, f"{med:.3f}", f"{mad:.3f}", len(v)])
            n += 1
    print(f"[done] famCN for {n}/{len(ids)} ids -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
