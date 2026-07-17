#!/usr/bin/env python3
"""WSSD read-depth copy-number leg — the DEPTH axis (Bailey 2002 / Sudmant 2010 / QuicK-mer2 analog) that the
minimap2 + SEDEF + parcn pipeline does NOT cover. It is the orthogonal, DNA-read confirmation of the
assembly/sequence-based famCN (SEDEF copy count) and parcn (consensus projection).

Principle (WSSD): read depth over a locus, corrected for GC bias and normalized by the single-copy depth
lambda, equals that locus's copy NUMBER *in the sequenced individual*. Summed over a segdup family's reference
copies it gives famCN; evaluated at SUN positions (bases private to one copy) it partitions into per-copy
parCN. Unlike the assembly count (which is fixed at the reference's copies), depth detects copy-number
VARIATION in the individual (gain/loss vs the reference) and un-collapses copies the assembly merged.

    famCN_individual = n_ref_copies * (GC-corrected mean depth over the family) / lambda
    parCN_copy_i     = (depth of copy-i's SUN allele) / lambda        (sum_i parCN_i ~= famCN)

STATUS: no gorilla DNA WGS on hand yet ([[parcn_feasibility]] is data-blocked). This is the ready-to-run leg:
map WGS to GGO.fasta/mGorGor1, then
  python bench/wssd_depth_cn.py --dna-bam gorilla_wgs.bam --fasta GGO.fasta \\
        --segdup-bed /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed --out gorilla_wssd
  # add --sun sun.tsv (copy_id,chrom,pos,allele) from the parcn pipeline for per-copy parCN.
--dry-run: print the inputs it needs and exit.  --smoke BAM: exercise the mechanics on any BAM (NOT a CN result).
"""
import argparse
import csv
import statistics
from bisect import bisect_right
from collections import defaultdict

import pysam

WIN = 1000            # depth window (bp)
GC_BIN = 0.02         # GC-correction bin width
N_LAMBDA_WIN = 4000   # single-copy windows sampled for the lambda / GC-correction control set


def merge_intervals(ivs, gap=0):
    ivs = sorted(ivs); out = []
    for s, e in ivs:
        if out and s <= out[-1][1] + gap:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def load_segdup_families(bed):
    """SEDEF pairs (cols 1-3 = sideA, 4-6 = sideB) -> merged copy loci -> union-find -> families (>=2 loci)."""
    raw = defaultdict(list)   # chrom -> [(s,e), ...] (both sides of every pair)
    pairs = []
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or not f[0].startswith(("chr", "NC_", "NW_")):
            continue
        try:
            a = (f[0], int(f[1]), int(f[2])); b = (f[3], int(f[4]), int(f[5]))
        except ValueError:
            continue
        raw[a[0]].append((a[1], a[2])); raw[b[0]].append((b[1], b[2])); pairs.append((a, b))
    # merge each chrom's intervals into distinct copy loci; index for lookup
    loci = []                 # global list of (chrom, s, e)
    starts = defaultdict(list); lid = defaultdict(list)
    for ch, ivs in raw.items():
        for s, e in merge_intervals(ivs):
            loci.append((ch, s, e)); starts[ch].append(s); lid[ch].append(len(loci) - 1)

    def locus_of(ch, s, e):
        arr = starts.get(ch)
        if not arr:
            return None
        i = bisect_right(arr, s) - 1
        for j in (i, i + 1):
            if 0 <= j < len(arr):
                gi = lid[ch][j]; _, ls, le = loci[gi]
                if not (ls > e or le < s):
                    return gi
        return None

    parent = list(range(len(loci)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in pairs:
        ia, ib = locus_of(*a), locus_of(*b)
        if ia is not None and ib is not None and ia != ib:
            parent[find(ia)] = find(ib)
    fams = defaultdict(list)
    for i in range(len(loci)):
        fams[find(i)].append(loci[i])
    return [v for v in fams.values() if len(v) >= 2]


def gc_fraction(fa, chrom, s, e):
    seq = fa.fetch(chrom, s, e).upper()
    n = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / n if n else None


def mean_depth(bam, chrom, s, e):
    """Mean primary-read per-base coverage (count_coverage skips UNMAP/SECONDARY/QCFAIL/DUP)."""
    try:
        cov = bam.count_coverage(chrom, s, e, quality_threshold=15)
    except (ValueError, KeyError):
        return None
    L = e - s
    return sum(cov[0][i] + cov[1][i] + cov[2][i] + cov[3][i] for i in range(L)) / L if L else 0.0


def calibrate_lambda(bam, fa, seg_by_chrom):
    """lambda = GC-corrected median depth over single-copy (non-segdup) windows. Returns (lambda, gc_corr)."""
    samples = []  # (gc, depth)
    for chrom in fa.references:
        clen = fa.get_reference_length(chrom)
        if clen < 5 * WIN:
            continue
        step = max(WIN, clen // 200)
        for s in range(WIN, clen - 2 * WIN, step):
            e = s + WIN
            if any(not (a > e or b < s) for a, b in seg_by_chrom.get(chrom, ())):
                continue  # skip segdup windows
            gc = gc_fraction(fa, chrom, s, e)
            d = mean_depth(bam, chrom, s, e)
            if gc is not None and d is not None and d > 0:
                samples.append((gc, d))
            if len(samples) >= N_LAMBDA_WIN:
                break
        if len(samples) >= N_LAMBDA_WIN:
            break
    if not samples:
        return None, {}
    gmed = statistics.median(d for _, d in samples)
    by_bin = defaultdict(list)
    for gc, d in samples:
        by_bin[round(gc / GC_BIN)].append(d)
    gc_corr = {b: (gmed / statistics.median(ds)) for b, ds in by_bin.items() if statistics.median(ds) > 0}
    return gmed, gc_corr


def corrected_depth(bam, fa, chrom, s, e, gc_corr):
    d = mean_depth(bam, chrom, s, e); gc = gc_fraction(fa, chrom, s, e)
    if d is None or gc is None:
        return None
    return d * gc_corr.get(round(gc / GC_BIN), 1.0)


def sun_parcn(bam, sun_path, lam):
    """parCN per copy from SUN allele depth. --sun columns: copy_id, chrom, pos(0-based), allele."""
    obs = defaultdict(lambda: [0, 0])  # copy_id -> [allele_read_sum, n_sites]
    for r in csv.reader(open(sun_path), delimiter="\t"):
        if len(r) < 4 or r[0].startswith("copy"):
            continue
        cid, chrom, pos, allele = r[0], r[1], int(r[2]), r[3].upper()
        try:
            cov = bam.count_coverage(chrom, pos, pos + 1, quality_threshold=15)
        except (ValueError, KeyError):
            continue
        idx = "ACGT".find(allele)
        if idx < 0:
            continue
        obs[cid][0] += cov[idx][0]; obs[cid][1] += 1
    return {cid: (v[0] / v[1] / lam if v[1] and lam else 0.0, v[1]) for cid, v in obs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dna-bam"); ap.add_argument("--fasta"); ap.add_argument("--segdup-bed")
    ap.add_argument("--sun"); ap.add_argument("--out", default="wssd")
    ap.add_argument("--dry-run", action="store_true"); ap.add_argument("--smoke")
    a = ap.parse_args()

    if a.dry_run or not (a.dna_bam or a.smoke):
        print(__doc__)
        print("REQUIRED when data arrives: --dna-bam (WGS mapped to the reference), --fasta (+.fai),")
        print("  --segdup-bed (GGO_sedef_final.bed). OPTIONAL: --sun (from parcn) for per-copy parCN.")
        return

    bam_path = a.smoke or a.dna_bam
    bam = pysam.AlignmentFile(bam_path, "rb")
    fa = pysam.FastaFile(a.fasta)
    fams = load_segdup_families(a.segdup_bed) if a.segdup_bed else []
    seg_by_chrom = defaultdict(list)
    for fam in fams:
        for ch, s, e in fam:
            seg_by_chrom[ch].append((s, e))
    print(f"segdup families (>=2 copies): {len(fams)}")

    lam, gc_corr = calibrate_lambda(bam, fa, seg_by_chrom)
    if not lam:
        print("could not calibrate lambda (no usable single-copy windows)."); return
    tag = "  [SMOKE: mechanics only, depth != CN unless this is DNA WGS]" if a.smoke else ""
    print(f"lambda (single-copy depth) = {lam:.2f}   GC bins = {len(gc_corr)}{tag}")

    with open(f"{a.out}.wssd_famcn.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "n_ref_copies", "mean_corrected_depth", "lambda", "famCN_individual", "copies"])
        for i, fam in enumerate(sorted(fams, key=lambda f: -len(f))):
            ds = [corrected_depth(bam, fa, ch, s, e, gc_corr) for ch, s, e in fam]
            ds = [d for d in ds if d is not None]
            if not ds:
                continue
            md = statistics.mean(ds)
            famcn = len(fam) * md / lam
            copies = ";".join(f"{ch}:{s}-{e}" for ch, s, e in fam[:6])
            w.writerow([f"WSSD{i}", len(fam), f"{md:.2f}", f"{lam:.2f}", f"{famcn:.2f}", copies])
    print(f"wrote {a.out}.wssd_famcn.tsv")

    if a.sun:
        pc = sun_parcn(bam, a.sun, lam)
        with open(f"{a.out}.wssd_parcn.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["copy_id", "n_sun_sites", "parCN"])
            for cid, (cn, nsites) in sorted(pc.items()):
                w.writerow([cid, nsites, f"{cn:.2f}"])
        print(f"wrote {a.out}.wssd_parcn.tsv ({len(pc)} copies)")


if __name__ == "__main__":
    main()
