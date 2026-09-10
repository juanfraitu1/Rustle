#!/usr/bin/env python3
"""PREREG_exon_sum_width_bias (64f80dd7): does the de novo (annotation-free-coordinate) read span
under-represent the true (annotation) copy width, by how much, which end, and does min_reads matter?
Also runs the identical rule over flair / StringTie / isoseq GTFs for comparison.

  python3 bench/width_deficit.py --copies copies16.tsv --bam hsa16.bam --gtf ours_final2_g2.gtf \
      --label ours flair=flair_family.gtf stringtie=stringtie_family.gtf isoseq=isoseq_family.gff
"""
import argparse, re, csv, subprocess, statistics
from collections import defaultdict, Counter


def load_copies(path):
    cop = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        cop[r["copy_idx"]] = dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                                   strand=r["strand"], n_reads=int(r["n_reads"]), family=r["family_id"])
    return cop


def copy_of(cop, chrom, s, e, min_contain=0.5):
    """Nearest copy by overlap, but ONLY if >= min_contain of THIS interval's own length lies inside the
    copy window (register row 879: raw max-overlap with no containment floor let one grazing outlier -- a
    read-through or chimeric alignment -- dominate a copy's span). None if no copy passes the floor."""
    length = e - s
    if length <= 0:
        return None
    best, best_ov = None, 0
    for ci, c in cop.items():
        if c["chrom"] != chrom:
            continue
        ov = min(e, c["end"]) - max(s, c["start"])
        if ov > best_ov and ov / length >= min_contain:
            best, best_ov = ci, ov
    return best


def trimmed_span(starts, ends, trim=0.05):
    """5th/95th percentile boundary (min_frac trim each side) -- robust to a single outlier alignment,
    the same failure class `snap_boundary` (Rust) guards the assembler's own boundary against."""
    if not starts:
        return None, None
    s = sorted(starts)
    e = sorted(ends)
    lo = int(len(s) * trim)
    hi = max(lo, len(e) - 1 - int(len(e) * trim))
    return s[lo], e[hi]


def span_from_gtf(path, cop, min_contain=0.5, trim=0.0):
    """per copy: robust (start, end) over every transcript's OWN genomic span, bucketed by copy_of
    (>= min_contain of the TRANSCRIPT's own length inside the copy -- excludes readthrough/chimeric
    transcripts that merely graze the window); `trim` percentile-trims each side (0 = raw min/max, safe
    once containment is enforced -- a transcript is a collapsed, already-deduplicated object, not a raw read)."""
    ex = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#") or f[2] != "exon":
            continue
        t = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
        if not t:
            continue
        ex[t.group(1)].append((int(f[3]) - 1, int(f[4]), f[0]))
    buckets = defaultdict(lambda: ([], []))
    for t, v in ex.items():
        v.sort()
        chrom = v[0][2]
        s0, e0 = v[0][0], v[-1][1]
        ci = copy_of(cop, chrom, s0, e0, min_contain)
        if ci is None:
            continue
        buckets[ci][0].append(s0)
        buckets[ci][1].append(e0)
    return {ci: trimmed_span(s, e, trim) for ci, (s, e) in buckets.items()}


def span_from_bam(bam, cop, min_contain=0.5, trim=0.05):
    """per copy: robust (start, end) over every primary bucketed by copy_of; no min_reads floor. `trim`
    percentile-trims each side -- individual reads (unlike collapsed transcripts) legitimately include
    chimeric/mis-mapped outliers, so a plain min/max is not safe here."""
    buckets = defaultdict(lambda: ([], []))
    out = subprocess.run(["samtools", "view", "-F", "2308", bam], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t", 6)
        chrom, pos, cig = f[2], int(f[3]) - 1, f[5]
        end = pos
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            if op in "M=XDN":
                end += int(n)
        ci = copy_of(cop, chrom, pos, end, min_contain)
        if ci is None:
            continue
        buckets[ci][0].append(pos)
        buckets[ci][1].append(end)
    return {ci: trimmed_span(s, e, trim) for ci, (s, e) in buckets.items()}


def deficits(cop, span, min_reads=0):
    """rows: (copy_idx, n_reads, true_width, denovo_width, deficit_total, deficit_5, deficit_3, rel)"""
    rows = []
    for ci, c in cop.items():
        if ci not in span or span[ci][0] is None:
            continue
        if c["n_reads"] < min_reads:
            continue
        s0, e0 = span[ci]
        true_w = c["end"] - c["start"]
        denovo_w = e0 - s0
        if c["strand"] == "+":
            d5 = s0 - c["start"]   # positive = read span starts AFTER the true TSS (under-shoot)
            d3 = c["end"] - e0     # positive = read span ends BEFORE the true TES (under-shoot)
        else:
            d5 = c["end"] - e0
            d3 = s0 - c["start"]
        rows.append((ci, c["n_reads"], true_w, denovo_w, true_w - denovo_w, d5, d3,
                      (true_w - denovo_w) / true_w if true_w else 0.0))
    return rows


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    rx = {v: i for i, v in enumerate(sorted(range(n), key=lambda k: xs[k]))}
    ry = {v: i for i, v in enumerate(sorted(range(n), key=lambda k: ys[k]))}
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n ** 2 - 1))


def report(label, rows, min_reads_thresh=10):
    if not rows:
        print(f"{label}: no copies with data")
        return
    sub = [r for r in rows if r[1] >= min_reads_thresh]
    n_under = sum(1 for r in sub if r[4] > 0)
    rel = [r[7] for r in sub]
    tot = [r[4] for r in rows]
    cv = (statistics.pstdev(tot) / abs(statistics.mean(tot))) if tot and statistics.mean(tot) else float("nan")
    nreads = [r[1] for r in rows]
    print(f"\n== {label}: {len(rows)} copies with data, {len(sub)} with n_reads>={min_reads_thresh}")
    if sub:
        print(f"   P1 under-represented (deficit>0): {n_under}/{len(sub)} = {100*n_under/len(sub):.0f}%  "
              f"median relative deficit {statistics.median(rel)*100:.1f}%  (IQR "
              f"{sorted(rel)[len(rel)//4]*100:.1f}..{sorted(rel)[3*len(rel)//4]*100:.1f}%)")
    print(f"   P2 deficit_total (bp): median {statistics.median(tot):.0f}  mean {statistics.mean(tot):.0f}  CV {cv:.2f}")
    rho = spearman(nreads, [r[7] for r in rows])
    print(f"   P3 spearman(n_reads, relative_deficit) = {rho:.2f}")
    d5 = [r[5] for r in rows]; d3 = [r[6] for r in rows]
    both_pos = sum(1 for a, b in zip(d5, d3) if a > 0 and b > 0)
    share5 = [a / (a + b) for a, b in zip(d5, d3) if (a + b) > 0]
    if share5:
        med5 = statistics.median(share5)
        print(f"   P5 median share of total deficit at the 5' end: {med5*100:.0f}%  "
              f"(3' end: {100-med5*100:.0f}%)  both-ends-positive: {both_pos}/{len(rows)}")
        print(f"      median |5' deficit| {statistics.median([abs(x) for x in d5]):.0f} bp, "
              f"median |3' deficit| {statistics.median([abs(x) for x in d3]):.0f} bp")
    return statistics.median(tot) if tot else None, statistics.median(rel) if rel else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--family")
    ap.add_argument("--min-nreads", type=int, default=10)
    ap.add_argument("--min-contain", type=float, default=0.8,
                     help="fraction of a transcript's/read's OWN length that must lie inside the copy "
                          "window to be bucketed there (register row 879 class: a low floor lets a single "
                          "readthrough/chimeric alignment dominate a copy's span)")
    ap.add_argument("tools", nargs="*", help="label=path.gtf, scored with the same rule")
    a = ap.parse_args()
    cop = load_copies(a.copies)
    print(f"copies loaded: {len(cop)}")

    span_gtf = span_from_gtf(a.gtf, cop, a.min_contain)
    rows_gtf = deficits(cop, span_gtf)
    med_gtf, medrel_gtf = report("ours (gtf, min_reads=3 collapse)", rows_gtf, a.min_nreads)

    span_raw = span_from_bam(a.bam, cop, a.min_contain)
    rows_raw = deficits(cop, span_raw)
    med_raw, medrel_raw = report("ours (raw BAM, no min_reads floor)", rows_raw, a.min_nreads)

    if med_gtf is not None and med_raw is not None and med_gtf:
        recovered = 100 * (med_gtf - med_raw) / med_gtf
        print(f"\nP4 min_reads contribution: gtf median deficit {med_gtf:.0f} bp -> raw {med_raw:.0f} bp "
              f"({recovered:.0f}% of the deficit recovered by dropping the min_reads floor)")

    results = {"ours": (med_gtf, medrel_gtf)}
    for spec in a.tools:
        label, path = spec.split("=", 1)
        span_t = span_from_gtf(path, cop, a.min_contain)
        rows_t = deficits(cop, span_t)
        m, mr = report(label, rows_t, a.min_nreads)
        results[label] = (m, mr)

    print("\nP6 median relative deficit by tool (report):")
    for l, (m, mr) in results.items():
        print(f"   {l:10s} {mr*100:.1f}%" if mr is not None else f"   {l:10s} n/a")


if __name__ == "__main__":
    main()
