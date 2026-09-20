#!/usr/bin/env python
"""Build the READ MULTI-MAPPING degree cache for the read-intrinsic TE gate.

For each bridge de-novo locus (bridge_dn_a/b in edge_bridge_mask.tsv, core present), fetch reads whose
alignment overlaps its SPLICED EXON intervals in GGO_mm.bam (mapped minimap2 -N50 -> multimappers carry up
to ~50 alignment records), then count each read's genome-wide alignment MULTIPLICITY (number of alignment
records sharing its QNAME, primary+secondary+supplementary). A read over a TE region multi-maps to
hundreds-thousands of loci (caps at ~50 here); over a gene-family shared exon, to a few. Per-locus feature =
median / mean / max of the read multiplicities. DETERMINISTIC: reads collected in coordinate order (sorted
BAM), capped per locus; counts are integers from a single full-file pass.

Reproduce: PYTHONHASHSEED=0 python bench/ri_build_readmm_cache.py
Writes: bench/ri_readmm_cache.tsv  (dn, n_reads, mm_median, mm_mean, mm_max)
"""
import os, sys, statistics
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"; os.execv(sys.executable, [sys.executable] + sys.argv)
BENCH = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, BENCH)
import family_er_pr as F
import family_er_te_gate as TE
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
BAM = os.path.join(SCRATCH, "GGO_mm.bam")
CACHE = os.path.join(BENCH, "ri_readmm_cache.tsv")
CAP = 400   # max reads sampled per locus (coordinate order -> deterministic)


def target_loci():
    dn = set()
    with open(os.path.join(BENCH, "edge_bridge_mask.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[idx["core"]] == "":
                continue
            for c in ("bridge_dn_a", "bridge_dn_b"):
                v = f[idx[c]]
                if v:
                    dn.add(v)
    return sorted(dn)


def main():
    meta = F.load_meta()
    skel = TE.load_skeletons()
    loci = target_loci()
    print(f"[readmm] {len(loci)} target bridge loci", flush=True)
    bam = pysam.AlignmentFile(BAM, "rb")
    ref_set = set(bam.references)

    # 1) collect per-locus qnames (coordinate order, capped) + global target set
    dn_qn = {}
    target = set()
    for i, dn in enumerate(loci):
        chrom, exons = TE.dn_exons(dn, meta, skel)
        if chrom is None or chrom not in ref_set:
            dn_qn[dn] = []
            continue
        seen = []
        seen_set = set()
        for (a, b) in exons:
            if len(seen) >= CAP:
                break
            for r in bam.fetch(chrom, a, b):
                qn = r.query_name
                if qn not in seen_set:
                    seen_set.add(qn); seen.append(qn)
                    if len(seen) >= CAP:
                        break
        dn_qn[dn] = seen
        target.update(seen)
        if (i + 1) % 400 == 0:
            print(f"[readmm] collected {i+1}/{len(loci)} loci, |target qnames|={len(target)}", flush=True)
    print(f"[readmm] total target qnames={len(target)}; full-file pass to count multiplicity ...", flush=True)

    # 2) one full pass: alignment-record count per target qname. FRESH handle -- reusing the region-fetched
    # handle leaves the file pointer near the end so until_eof scans only the tail (undercount bug).
    bam.close()
    bam2 = pysam.AlignmentFile(BAM, "rb")
    cnt = {}
    n = 0
    for r in bam2.fetch(until_eof=True):
        qn = r.query_name
        if qn in target:
            cnt[qn] = cnt.get(qn, 0) + 1
        n += 1
        if n % 2000000 == 0:
            print(f"[readmm] scanned {n} records", flush=True)
    print(f"[readmm] scanned {n} records total", flush=True)

    # 3) per-locus summary
    with open(CACHE, "w") as out:
        out.write("dn\tn_reads\tmm_median\tmm_mean\tmm_max\n")
        for dn in loci:
            degs = [cnt.get(q, 1) for q in dn_qn[dn]]
            if degs:
                out.write(f"{dn}\t{len(degs)}\t{statistics.median(degs):.4f}\t"
                          f"{statistics.mean(degs):.4f}\t{max(degs)}\n")
            else:
                out.write(f"{dn}\t0\t\t\t\n")
    print(f"[readmm] wrote {CACHE}", flush=True)


if __name__ == "__main__":
    main()
