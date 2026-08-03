#!/usr/bin/env python3
"""Compute famCN (WSSD read-depth copy number) at ARBITRARY coordinates, replicating Soto's CN leg.

WHY THIS EXISTS. Soto define a gene family as shared-exon genes whose family copy number agrees
(mean-absolute-deviation < 1), so famCN is their SPLITTER. Their supplementary table S1C publishes famCN,
but only at THEIR annotated genes -- and our over-merged families contain predicted copies with no Soto
gene, which is exactly where we need the number. This computes famCN ourselves from the per-sample WSSD
tracks, so it can be evaluated at any interval we predict.

VALIDATION (3 samples vs Soto's published famCN, n=269): SYCE1 2.19 vs 2.01, PMCHL1 3.93 vs 3.99,
OR4K7P 5.64 vs 5.98, USP32P3 8.37 vs 8.60, CBWD6 13.67 vs 13.70, AC133919.3 24.80 vs 22.76,
NPIPB15 46.67 vs 48.25 -- 7/8 within 9% over a 24x range. The 8th (BET1L, famCN 702, subtelomeric) needs
far more samples; treat CN > ~100 as unresolved at low sample counts.

COORDINATES. The WSSD tracks are T2T-CHM13 **v1.0**; we work in **v2.0**. Rather than a chain file, the
liftover is derived from Soto's own table S1E, which carries BOTH `SD98_v1.0` and `SD98_v2.0` for 1833
paralogs: per chromosome the difference is a CONSTANT (a few hundred bp), except the five acrocentrics,
which have exactly one breakpoint each in the short arm that v2.0 rebuilt. That gives a piecewise-constant
map fitted to real anchor pairs. Intervals landing in the uncertain zone between two regimes are reported
as unmapped rather than guessed -- a wrong offset would silently return another locus's copy number, which
is worse than an absent value.

CIRCULARITY WARNING. famCN is CONSTITUTIVE of the Soto ground truth: they built those families with it. Any
score that uses famCN to split and is then evaluated against Soto is partly circular and must say so.
The honest uses are (a) characterising what a CN signal would buy ("CN would resolve N of our M residual
over-merges"), and (b) splitting when evaluated on an INDEPENDENT truth set.
"""
import argparse, csv, os, re, statistics as st, subprocess, sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

BASE_URL = "http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/wssd"
GAP_GUARD = 50_000  # v2 bp around a regime switch treated as unmappable


def build_liftover(s1e_path):
    """Fit a piecewise-constant v2.0 -> v1.0 map from Soto's dual-coordinate anchors.

    Returns {chrom: [(v2_from, offset), ...]} sorted by v2_from, plus the anchor count per chromosome.
    """
    anchors = defaultdict(list)
    with open(s1e_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m1 = re.match(r"(chr[\w]+):(\d+)-(\d+)", r.get("SD98_v1.0", "") or "")
            m2 = re.match(r"(chr[\w]+):(\d+)-(\d+)", r.get("SD98_v2.0", "") or "")
            if not (m1 and m2) or m1.group(1) != m2.group(1):
                continue
            anchors[m1.group(1)].append((int(m2.group(2)), int(m1.group(2)) - int(m2.group(2))))
    table, spans = {}, {}
    for c, pts in anchors.items():
        pts.sort()
        regimes, prev = [], None
        for v2, off in pts:
            if off != prev:
                regimes.append((v2, off))
                prev = off
        table[c] = regimes
        # the v2 window where a regime switch happens is uncertain: between the last anchor of one
        # regime and the first of the next.
        bounds = []
        for i in range(1, len(regimes)):
            lo = max(v2 for v2, o in pts if o == regimes[i - 1][1])
            bounds.append((lo, regimes[i][0]))
        spans[c] = bounds
    return table, spans


def lift(table, spans, chrom, start, end):
    """v2.0 interval -> v1.0 interval, or None when the mapping is not trustworthy here."""
    regimes = table.get(chrom)
    if not regimes:
        return None
    for lo, hi in spans.get(chrom, ()):
        if start < hi + GAP_GUARD and end > lo - GAP_GUARD:
            return None  # straddles / sits inside a regime switch
    off = regimes[0][1]
    for v2_from, o in regimes:
        if start >= v2_from:
            off = o
    return (start + off, end + off)


def cn_for_interval(bb, chrom, start, end, tool):
    """Length-weighted mean CN over the WSSD windows covering [start, end)."""
    try:
        out = subprocess.run(
            [tool, f"-chrom={chrom}", f"-start={start}", f"-end={end}", bb, "/dev/stdout"],
            capture_output=True, text=True, timeout=300,
        ).stdout
    except subprocess.TimeoutExpired:
        return None
    tot = n = 0.0
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 10:
            continue
        try:
            s, e, cn = int(f[1]), int(f[2]), float(f[9])
        except ValueError:
            continue
        w = min(e, end) - max(s, start)
        if w > 0:
            tot += cn * w
            n += w
    return tot / n if n > 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--intervals", required=True,
                    help="TSV with chrom/start/end (+ any id columns); v2.0 coordinates")
    ap.add_argument("--wssd-dir", help="directory of *_wssd.bb; omit to stream from UCSC")
    ap.add_argument("--samples", type=int, default=8, help="how many samples to median over")
    ap.add_argument("--jobs", type=int, default=8, help="parallel bigBedToBed calls per interval")
    ap.add_argument("--s1e", default="bench/soto/soto_parCN_S1E.tsv")
    ap.add_argument("--tool", default="bigBedToBed")
    ap.add_argument("--id-col", default="family_id")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    table, spans = build_liftover(a.s1e)
    print(f"[liftover] fitted over {len(table)} chromosomes; "
          f"{sum(len(v) for v in spans.values())} regime switches guarded", file=sys.stderr)

    if a.wssd_dir:
        bbs = sorted(os.path.join(a.wssd_dir, f) for f in os.listdir(a.wssd_dir) if f.endswith("_wssd.bb"))
        if not bbs:
            sys.exit(f"no *_wssd.bb in {a.wssd_dir}")
    else:
        listing = subprocess.run(["curl", "-s", "--max-time", "60", BASE_URL + "/"],
                                 capture_output=True, text=True).stdout
        bbs = [f"{BASE_URL}/{m}" for m in re.findall(r'href="([^"]+_wssd\.bb)"', listing)]
    bbs = bbs[: a.samples]
    print(f"[wssd] {len(bbs)} sample track(s)", file=sys.stderr)

    rows = list(csv.DictReader(open(a.intervals), delimiter="\t"))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([a.id_col, "chrom", "start", "end", "v1_start", "v1_end",
                    "famCN_median", "famCN_mad", "n_samples", "status"])
        done = unmapped = 0
        for r in rows:
            if not r.get("chrom"):
                continue
            chrom, s, e = r["chrom"], int(r["start"]), int(r["end"])
            lifted = lift(table, spans, chrom, s, e)
            if lifted is None:
                w.writerow([r.get(a.id_col, ""), chrom, s, e, "", "", "", "", 0, "UNMAPPED_v1"])
                unmapped += 1
                continue
            v1s, v1e = lifted
            # One subprocess per (interval, sample); with the full 271-sample panel that is ~10.8k calls
            # for 40 intervals, so they run on a thread pool (each worker just waits on bigBedToBed).
            with ThreadPoolExecutor(max_workers=a.jobs) as ex:
                vals = [v for v in ex.map(lambda bb: cn_for_interval(bb, chrom, v1s, v1e, a.tool), bbs)
                        if v is not None]
            if not vals:
                w.writerow([r.get(a.id_col, ""), chrom, s, e, v1s, v1e, "", "", 0, "NO_CN"])
                continue
            med = st.median(vals)
            mad = st.median([abs(v - med) for v in vals])
            w.writerow([r.get(a.id_col, ""), chrom, s, e, v1s, v1e,
                        f"{med:.3f}", f"{mad:.3f}", len(vals), "OK"])
            done += 1
            if done % 50 == 0:
                print(f"  {done} intervals done", file=sys.stderr)
    print(f"[done] {done} with CN, {unmapped} unmappable -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
