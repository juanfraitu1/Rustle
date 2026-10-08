#!/usr/bin/env python3
"""The RNA-only allele caller of docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md (reads `_pri`, the BAM, the annotation-derived
gene table, the paralog PAF, the genome-wide AS table and the O2 assignments; never the haplotype assemblies).

Per gene G (sets S_fam, S_multi, S_single, S_X):
  reads   single-copy sets: primary records (not unmapped/secondary/supplementary) with MAPQ >= 20;
          family sets: primary records of reads not AS-tied genome-wide (second AS < 0.98 x best) + the record at G of AS-tied reads the
          shipped O2 assigned to G (status "assigned"); everything else is dropped.
  pass 1  per-position A/C/G/T counts from those reads over G's exons (pysam count_coverage, C speed). Expressed = some exonic position
          with >= 10 reads.
  sites   >= 10 reads; second base >= 3 reads and >= 0.20 of the covering reads; not a PSV (every paralog hit of G's exon sum at
          identity >= 0.90 over >= 0.50 that covers the column carries G's base there, no gap); not A>G on the transcript strand (reference
          A, other allele G); not within 3 bp of a homopolymer >= 5 bp.
  pass 2  genes with >= 2 sites: reads spanning >= 2 sites are phased (best 2-haplotype split, exhaustive up to 12 sites, else greedy);
          consistent when >= 80% of spanning reads fit one haplotype. No spanning read: phasing not testable, the call stands (ruling).
  call    2 (>= 1 site, phasing consistent or untestable) / 1+ (expressed, no site) / NA (not expressed) / inconsistent (sites, phasing < 80%).
  novel   S_fam only: reads of the copy (any primary record overlapping its exons, tied or not) whose bases at G's PSV columns differ at
          >= 2 columns from G and from every paralog covering those columns; grouped by identical PSV pattern, >= 2 reads per group.

Resumable per chromosome: writes <out-dir>/<chrom>.calls.tsv (+ .novel.tsv) and skips chromosomes already done; stops at --deadline.

    python3 caller.py --genes genes.tsv --sets sets.tsv --paf para.paf --molecules fibro.molecules.tsv --assign 'o2_*.assignments.tsv' \
        --catalog-glob 'cat*.copies.tsv' --bam BAM --fasta GGO.fasta --out-dir calls --deadline EPOCH
"""
import argparse
import collections
import csv
import glob
import itertools
import os
import re
import sys
import time

import numpy as np
import pysam

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"


def blocks(s):
    return [tuple(int(x) for x in b.split("-")) for b in s.split(",") if b]


def cigar_ops(cg):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHX=])", cg)]


class Paralogs:
    """Maps an exon-sum index of G to the base each paralog hit carries there, one entry per hit in hit order (transcript orientation;
    None = gap at that column; "?" = the hit does not cover the column)."""

    def __init__(self, hits, fasta):
        self.hits = hits
        self.fa = fasta

    def bases_at(self, i, qlen):
        out = []
        for h in self.hits:
            t, ts, strand, qs, qe, ops = h["t"], h["ts"], h["strand"], h["qs"], h["qe"], h["ops"]
            qi = i if strand == "+" else qlen - 1 - i      # position on the query strand the CIGAR walks
            q = qs if strand == "+" else qlen - qe
            if not (q <= qi < (qe if strand == "+" else qlen - qs)):
                out.append("?"); continue
            tp, got = ts, "out"
            for n, op in ops:
                if op in "M=X":
                    if q <= qi < q + n:
                        got = tp + (qi - q); break
                    q += n; tp += n
                elif op == "I":
                    if q <= qi < q + n:
                        got = None; break
                    q += n
                elif op in "DN":
                    tp += n
            if got == "out":
                out.append("?"); continue
            if got is None:
                out.append(None)
            else:
                b = self.fa.fetch(t, got, got + 1).upper()
                out.append(b if strand == "+" else b.translate(COMP))
        return out


def homopolymer_near(seq, pos, lo):
    """True if a run of >= 5 identical bases lies within 3 bp of pos (seq = reference from lo)."""
    k = pos - lo
    s = max(0, k - 8)
    e = min(len(seq), k + 9)
    sub = seq[s:e]
    for m in re.finditer(r"(A{5,}|C{5,}|G{5,}|T{5,})", sub):
        a, b = s + m.start(), s + m.end()          # run [a, b) in seq coords
        if a - 3 <= k < b + 3:
            return True
    return False


def phase(reads_alleles, k):
    """reads_alleles: list of {site_index: 0/1}; returns (consistent_fraction, n_spanning)."""
    span = [r for r in reads_alleles if len(r) >= 2]
    if not span:
        return None, 0
    def score(h):
        return sum(1 for r in span if all(h[j] == v for j, v in r.items()) or all(1 - h[j] == v for j, v in r.items()))
    if k <= 12:
        best = max(score((0,) + rest) for rest in itertools.product((0, 1), repeat=k - 1))
    else:
        cov = collections.Counter(j for r in span for j in r)
        order = [j for j, _ in cov.most_common()] + [j for j in range(k) if j not in cov]
        h = {order[0]: 0}
        for j in order[1:]:
            v = collections.Counter()
            for r in span:
                if j in r:
                    for i2, a in r.items():
                        if i2 in h:
                            v[(a == r[j]) == (h[i2] == 0)] += 1     # same orientation vote
            h[j] = 0 if v[True] >= v[False] else 1
        best = score(tuple(h.get(j, 0) for j in range(k)))
    return best / len(span), len(span)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("genes", "sets", "paf", "molecules", "assign", "catalog_glob", "bam", "fasta", "out_dir"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    ap.add_argument("--deadline", type=float, default=0)
    ap.add_argument("--chroms", default="")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    genes = {r["gene_id"]: r for r in csv.DictReader(open(a.genes), delimiter="\t")}
    sets = {r["gene_id"]: r["set"] for r in csv.DictReader(open(a.sets), delimiter="\t")}
    todo = collections.defaultdict(list)
    for g, s in sets.items():
        if s in ("S_fam", "S_multi", "S_single", "S_X"):
            todo[genes[g]["chrom"]].append(g)
    chroms = [c for c in sorted(todo) if not os.path.exists(f"{a.out_dir}/{c}.calls.tsv")]
    if a.chroms:
        chroms = [c for c in chroms if c in a.chroms.split(",")]
    if not chroms:
        print("ALL_DONE"); return 0
    t0 = time.time()
    # paralog hits (the registered paralog definition) for family-set genes
    para = collections.defaultdict(list)
    qlen = {}
    for ln in open(a.paf):
        f = ln.rstrip("\n").split("\t")
        g = f[0]
        if sets.get(g) not in ("S_fam", "S_multi"):
            continue
        qlen[g] = int(f[1])
        ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        bl = blocks(genes[g]["exons"])
        own = f[5] == genes[g]["chrom"] and int(f[7]) < bl[-1][1] and bl[0][0] < int(f[8])
        if own or ident < 0.90 or cov < 0.50:
            continue
        cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        para[g].append(dict(t=f[5], ts=int(f[7]), strand=f[4], qs=int(f[2]), qe=int(f[3]), ops=cigar_ops(cg)))
    # genome-wide tied reads
    tied = set()
    with open(a.molecules) as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            n, b, s2 = ln.split("\t", 3)[:3]
            b, s2 = int(b), int(s2)
            if s2 > 0 and s2 >= 0.98 * b:
                tied.add(n)
    # O2 assignments -> gene id via the batch catalogs
    tid = {}
    for path in glob.glob(a.catalog_glob):
        for r in csv.DictReader(open(path), delimiter="\t"):
            tid[(r["family_id"], r["copy_idx"])] = r["tid"]
    assigned = collections.defaultdict(set)
    for path in glob.glob(a.assign):
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["status"] == "assigned" and (r["family_id"], r["catalog_copy_idx"]) in tid:
                assigned[tid[(r["family_id"], r["catalog_copy_idx"])]].add(r["read_name"])
    print(f"[caller] tied {len(tied):,}; O2-assigned reads {sum(len(v) for v in assigned.values()):,}; load {time.time() - t0:.0f}s",
          flush=True)
    bam = pysam.AlignmentFile(a.bam)
    fa = pysam.FastaFile(a.fasta)
    PART = 300
    for chrom in chroms:
      glist = todo[chrom]
      for pk in range(0, len(glist), PART):
        part = f"{a.out_dir}/{chrom}.part{pk // PART:04d}"
        if os.path.exists(part + ".calls.tsv"):
            continue
        if a.deadline and time.time() > a.deadline - 90:
            print(f"[caller] deadline before {chrom} part {pk // PART}"); sys.exit(75)
        rows, novel_rows = [], []
        for g in glist[pk:pk + PART]:
            r = genes[g]
            st = sets[g]
            bl = blocks(r["exons"])
            lo, hi = bl[0][0], bl[-1][1]
            fam = st in ("S_fam", "S_multi")
            mine = assigned.get(g, set())
            if fam:
                def keep(rd, mine=mine):
                    if rd.is_unmapped or rd.is_supplementary:
                        return False
                    if rd.is_secondary:
                        return rd.query_name in mine
                    return rd.query_name not in tied or rd.query_name in mine
            else:
                def keep(rd):
                    return not (rd.is_unmapped or rd.is_secondary or rd.is_supplementary) and rd.mapping_quality >= 20
            cnt = bam.count_coverage(chrom, lo, hi, quality_threshold=0, read_callback=keep)
            ref = fa.fetch(chrom, max(0, lo - 20), hi + 20).upper()
            roff = max(0, lo - 20)
            # exon-sum index of a genomic position
            exlen = sum(y - x for x, y in bl)
            def esum_index(p):
                acc = 0
                for x, y in bl:
                    if x <= p < y:
                        i = acc + (p - x)
                        return i if r["strand"] == "+" else exlen - 1 - i
                    acc += y - x
                return None
            P = Paralogs(para.get(g, []), fa) if fam else None
            maxcov, sites = 0, []
            C = np.array([np.frombuffer(cnt[k], dtype=np.dtype(f"u{cnt[k].itemsize}")) for k in range(4)],
                         dtype=np.int64)
            for x, y in bl:
                sub = C[:, x - lo:y - lo]
                if sub.shape[1] == 0:
                    continue
                tot = sub.sum(axis=0)
                maxcov = max(maxcov, int(tot.max()))
                srt = np.sort(sub, axis=0)
                sec = srt[2]
                cand = np.nonzero((tot >= 10) & (sec >= 3) & (sec >= 0.20 * tot))[0]
                for ci in cand:
                    p = x + int(ci)
                    c = [int(C[k][p - lo]) for k in range(4)]
                    order = sorted(range(4), key=lambda k: -c[k])
                    rb = ref[p - roff]
                    alleles = {BASES[order[0]], BASES[order[1]]}
                    alt = (alleles - {rb}).pop() if rb in alleles else BASES[order[1]]
                    # editing: reference A, other allele G on the transcript strand
                    tr_ref = rb if r["strand"] == "+" else rb.translate(COMP)
                    tr_alt = alt if r["strand"] == "+" else alt.translate(COMP)
                    if tr_ref == "A" and tr_alt == "G":
                        continue
                    if homopolymer_near(ref, p, roff):
                        continue
                    if fam:
                        i = esum_index(p)
                        gb = rb if r["strand"] == "+" else rb.translate(COMP)
                        pb = P.bases_at(i, qlen.get(g, exlen))
                        if any(b != "?" and b != gb for b in pb):
                            continue
                    sites.append((p, BASES[order[0]], BASES[order[1]]))
            expressed = maxcov >= 10
            call, frac, nspan = ("NA" if not expressed else "1+"), "", 0
            if sites:
                call = "2"
                if len(sites) >= 2:
                    pos = {p: j for j, (p, _, _) in enumerate(sites)}
                    ra = []
                    for rd in bam.fetch(chrom, lo, hi):
                        if not keep(rd):
                            continue
                        d = {}
                        for qp, rp in rd.get_aligned_pairs(matches_only=True):
                            j = pos.get(rp)
                            if j is not None:
                                b = rd.query_sequence[qp]
                                if b == sites[j][1]:
                                    d[j] = 0
                                elif b == sites[j][2]:
                                    d[j] = 1
                        if d:
                            ra.append(d)
                    fr, nspan = phase(ra, len(sites))
                    if fr is not None:
                        frac = f"{fr:.3f}"
                        if fr < 0.80:
                            call = "inconsistent"
            rows.append((g, st, int(expressed), maxcov, len(sites), call, frac, nspan,
                         ";".join(f"{p}:{a1}/{a2}" for p, a1, a2 in sites[:50])))
            # novel haplotypes (S_fam): PSV columns of G and per-read patterns
            if st == "S_fam" and expressed:
                psv = []
                for x, y in bl:
                    for p in range(x, y):
                        i = esum_index(p)
                        rb = ref[p - roff]
                        gb = rb if r["strand"] == "+" else rb.translate(COMP)
                        pb = P.bases_at(i, qlen.get(g, exlen))
                        if any(b != "?" and b != gb for b in pb):
                            psv.append((p, i, gb, pb))
                if psv:
                    pidx = {p: k for k, (p, _, _, _) in enumerate(psv)}
                    groups = collections.Counter()
                    for rd in bam.fetch(chrom, lo, hi):
                        if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                            continue
                        pat = {}
                        for qp, rp in rd.get_aligned_pairs(matches_only=True):
                            k = pidx.get(rp)
                            if k is not None:
                                b = rd.query_sequence[qp]
                                pat[k] = b if r["strand"] == "+" else b.translate(COMP)
                        if len(pat) < 2:
                            continue
                        dG = sum(1 for k, b in pat.items() if b != psv[k][2])
                        if dG < 2:
                            continue
                        npar = len(para.get(g, []))
                        ok = True
                        for hi_ in range(npar):
                            cols = [(k, b) for k, b in pat.items() if psv[k][3][hi_] not in ("?", None)]
                            d = sum(1 for k, b in cols if b != psv[k][3][hi_])
                            if d < 2:
                                ok = False; break
                        if ok:
                            groups[tuple(sorted(pat.items()))] += 1
                    for pat, n in groups.items():
                        if n >= 2:
                            novel_rows.append((g, n, len(pat), "".join(b for _, b in pat)))
        with open(part + ".calls.tsv.tmp", "w") as out:
            for row in rows:
                out.write("\t".join(str(x) for x in row) + "\n")
        with open(part + ".novel.tsv", "w") as out:
            for row in novel_rows:
                out.write("\t".join(str(x) for x in row) + "\n")
        os.replace(part + ".calls.tsv.tmp", part + ".calls.tsv")
        print(f"[caller] {chrom} part {pk // PART}: {len(rows)} genes, {sum(1 for x in rows if x[5] == '2')} called 2, "
              f"{time.time() - t0:.0f}s", flush=True)
      # all parts present: assemble the chromosome files
      parts = sorted(glob.glob(f"{a.out_dir}/{chrom}.part*.calls.tsv"))
      if len(parts) == (len(glist) + PART - 1) // PART:
          with open(f"{a.out_dir}/{chrom}.calls.tsv.tmp", "w") as out:
              out.write("gene_id\tset\texpressed\tmax_cov\tn_sites\tcall\tphase_consistent\tn_spanning\tsites\n")
              for pp in parts:
                  out.write(open(pp).read())
          with open(f"{a.out_dir}/{chrom}.novel.tsv", "w") as out:
              out.write("gene_id\tn_reads\tn_psv_columns\tpattern\n")
              for pp in sorted(glob.glob(f"{a.out_dir}/{chrom}.part*.novel.tsv")):
                  out.write(open(pp).read())
          os.replace(f"{a.out_dir}/{chrom}.calls.tsv.tmp", f"{a.out_dir}/{chrom}.calls.tsv")
    print("ALL_DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
