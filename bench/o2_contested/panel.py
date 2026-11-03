#!/usr/bin/env python3
"""r1509 panel: families ranked by their ms-near-tie reads (register 1508's O2 target set), from an extended
`as_table` (columns 7-17) and a `PREFIX.fam.copies.tsv`. A molecule is in the target set when it has a primary and a
rival, ms margin <= W, and is not the K = 0 wall (AS and ms margins 0, same splice shape, same chromosome). It counts
for every family with a copy whose window (copy +- FLANK) holds the primary's position.
Usage: panel.py MOLECULES.tsv COPIES.tsv OUT_DIR [--w 10] [--flank 50000] [--top 40]
Writes OUT_DIR/families.txt, OUT_DIR/regions.txt (per family and contig: hull of its copies +- 5 kb, merged per contig
as rustle_pipeline.sh fam_regions does), OUT_DIR/panel.tsv (family, n_target, n_wall, n_informative, n_copies)."""
import argparse, bisect, collections, csv, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("molecules"); ap.add_argument("copies"); ap.add_argument("out")
    ap.add_argument("--w", type=int, default=10); ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--top", type=int, default=40)
    a = ap.parse_args()
    copies = collections.defaultdict(list)  # chrom -> [(lo, hi, fam)]
    fam_copies = collections.Counter(); fam_hull = {}
    with open(a.copies) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, s, e, f = r["chrom"], int(r["start"]), int(r["end"]), r["family_id"]
            copies[c].append((max(0, s - a.flank), e + a.flank, f)); fam_copies[f] += 1
            k = (f, c); lo, hi = fam_hull.get(k, (s, e)); fam_hull[k] = (min(lo, s), max(hi, e))
    idx = {c: (sorted(v), ) for c, v in copies.items()}
    starts = {c: [x[0] for x in v[0]] for c, v in idx.items()}
    n_target = collections.Counter(); n_wall = collections.Counter()
    with open(a.molecules) as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 17 or f[4] == "NA" or f[11] == "NA": continue
            pas, pms, ras, rms = int(f[5]), int(f[8]), int(f[11]), int(f[12])
            mm = pms - rms
            if mm > a.w: continue
            wall = mm == 0 and pas - ras == 0 and f[10] == f[13] and f[4] == f[14]
            chrom, pos = f[4], int(f[15])
            if chrom not in idx: continue
            v = idx[chrom][0]; i = bisect.bisect_right(starts[chrom], pos)
            fams = set()
            for j in range(max(0, i - 50), i):  # windows are <= ~200 kb; 50 back is ample for nested copies
                lo, hi, fam = v[j]
                if lo <= pos <= hi: fams.add(fam)
            for fam in fams:
                (n_wall if wall else n_target)[fam] += 1
    rows = sorted(((n_target[f], n_wall[f], f) for f in set(n_target) | set(n_wall)), reverse=True)
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/panel.tsv", "w") as fh:
        fh.write("family\tn_informative\tn_wall\tn_copies\n")
        for t, w, f in rows: fh.write(f"{f}\t{t}\t{w}\t{fam_copies[f]}\n")
    pick = [f for _, _, f in rows[: a.top]]
    with open(f"{a.out}/families.txt", "w") as fh: fh.write("\n".join(pick) + "\n")
    per_chrom = collections.defaultdict(list)
    for (f, c), (lo, hi) in fam_hull.items():
        if f in set(pick): per_chrom[c].append((max(0, lo - 5000), hi + 5000))
    with open(f"{a.out}/regions.txt", "w") as fh:
        for c in sorted(per_chrom):
            cur = None
            for lo, hi in sorted(per_chrom[c]):
                if cur and lo <= cur[1]: cur[1] = max(cur[1], hi)
                else:
                    if cur: fh.write(f"{c}:{cur[0]}-{cur[1]}\n")
                    cur = [lo, hi]
            if cur: fh.write(f"{c}:{cur[0]}-{cur[1]}\n")
    tot = sum(n_target.values()); print(f"families with target reads: {len(rows)}; informative read-family pairs: {tot}; "
          f"top {len(pick)} hold {sum(n_target[f] for f in pick)} ({100*sum(n_target[f] for f in pick)/max(tot,1):.0f}%)")
    for t, w, f in rows[:15]: print(f"  {f}\t{t}\t{w}\t{fam_copies[f]} copies")

if __name__ == "__main__": main()
