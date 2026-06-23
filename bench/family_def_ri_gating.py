#!/usr/bin/env python3
"""family_def_ri_gating.py — (a) does read-supported block gating clean S(v) and lift the
tau-edge families off the cliff?

Current S(v) = union of ALL reads' aligned blocks (any base covered by >=1 read), so a single
intron-read-through pads S(v) with a whole intron -> depresses ~B reciprocal coverage. Read-
supported S(v) keeps only bases covered by >=K reads (count_coverage, intron-aware: N-skips
contribute 0), so single-read padding is removed. We rebuild both copies of each family with
K=1 (current) and K=3, recompute reciprocal coverage, and ask:
  - do tau-edge families (recip in [0.30,0.55]) RISE off the edge?
  - do clean high-recip families stay high (control)? are models shrinking uniformly or
    specifically losing thin padding?
Run: /home/juanfra/miniforge3/bin/python bench/family_def_ri_gating.py
"""
import collections
import json
import os
import subprocess
import sys
import statistics as st

import pysam

K_SUP = int(sys.argv[1]) if len(sys.argv) > 1 else 3   # read-support threshold for gated model

HERE = os.path.dirname(os.path.abspath(__file__))
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
OLD = "/home/juanfra/winloci_scratch/GGO.bam"
TMP = "/home/juanfra/winloci_scratch/rigate"
PAD, MAX_DEPTH = 5000, 200000
OVERLAP_MERGE, COV_MIN = 0.5, 0.30


def model_at(bam, genome, chrom, start, end, k):
    """exon-union over [start,end]+-PAD keeping only bases with read coverage >= k."""
    s0, e0 = max(0, start - PAD), end + PAD
    try:
        cov = bam.count_coverage(chrom, s0, e0, quality_threshold=0,
                                 read_callback=lambda r: not (r.is_secondary or r.is_supplementary
                                                              or r.is_unmapped))
    except (ValueError, KeyError):
        return ""
    depth = [cov[0][i] + cov[1][i] + cov[2][i] + cov[3][i] for i in range(len(cov[0]))]
    intervals, run_s = [], None
    for i, d in enumerate(depth):
        if d >= k and run_s is None:
            run_s = s0 + i
        elif d < k and run_s is not None:
            intervals.append((run_s, s0 + i)); run_s = None
    if run_s is not None:
        intervals.append((run_s, e0))
    return "".join(genome.fetch(chrom, a, b).upper() for a, b in intervals)


def recip_cov(a, b):
    os.makedirs(TMP, exist_ok=True)
    if not a or not b:
        return 0.0
    open(f"{TMP}/a.fa", "w").write(f">a\n{a}\n")
    open(f"{TMP}/b.fa", "w").write(f">b\n{b}\n")
    out = subprocess.run(f"minimap2 -c -x asm20 {TMP}/a.fa {TMP}/b.fa 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    best = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        ql, qs, qe, tl, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        best = max(best, min((qe - qs) / max(ql, 1), (te - ts) / max(tl, 1)))
    return best


def dedup(loci):
    by = collections.defaultdict(list)
    for lid, c, s, e in loci:
        by[c].append((s, e))
    regs = []
    for c, it in by.items():
        it.sort(); cur = None
        for s, e in it:
            if cur and s < cur[1] and (min(cur[1], e) - s) > OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                cur = (cur[0], max(cur[1], e)); continue
            if cur:
                regs.append((c, cur[0], cur[1]))
            cur = (s, e)
        if cur:
            regs.append((c, cur[0], cur[1]))
    return regs


def main():
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((lid, c, int(s), int(e)))
    bam = pysam.AlignmentFile(OLD, "rb"); genome = pysam.FastaFile(GENOME)

    rows = []
    for fi in sorted(fams):
        regs = dedup(fams[fi])
        if len(regs) != 2:                  # size-2 families = clean single ~B pair
            continue
        (c1, s1, e1), (c2, s2, e2) = regs
        m1_1 = model_at(bam, genome, c1, s1, e1, 1); m2_1 = model_at(bam, genome, c2, s2, e2, 1)
        m1_3 = model_at(bam, genome, c1, s1, e1, K_SUP); m2_3 = model_at(bam, genome, c2, s2, e2, K_SUP)
        if not (m1_1 and m2_1 and m1_3 and m2_3):
            continue
        rc1 = recip_cov(m1_1, m2_1); rc3 = recip_cov(m1_3, m2_3)
        shrink = 1 - (len(m1_3) + len(m2_3)) / max(len(m1_1) + len(m2_1), 1)
        rows.append(dict(family=fi, recip_k1=round(rc1, 4), recip_k3=round(rc3, 4),
                         delta=round(rc3 - rc1, 4), shrink=round(shrink, 4)))
    bam.close(); genome.close()

    edge = [r for r in rows if 0.30 <= r["recip_k1"] <= 0.55]
    high = [r for r in rows if r["recip_k1"] > 0.55]
    json.dump(dict(rows=rows), open(os.path.join(HERE, "family_def_ri_gating.json"), "w"), indent=2)

    def report(g, name):
        if not g:
            print(f"  {name}: (none)"); return
        d = [r["delta"] for r in g]
        rose = sum(1 for r in g if r["delta"] > 0.05)
        offedge = sum(1 for r in g if r["recip_k1"] <= 0.55 < r["recip_k3"])
        print(f"  {name} (n={len(g)}): median recip K1={st.median(r['recip_k1'] for r in g):.3f}"
              f" -> K3={st.median(r['recip_k3'] for r in g):.3f}; "
              f"median delta={st.median(d):+.3f}; rose>0.05: {rose}; "
              f"crossed 0.55: {offedge}; median model shrink={st.median(r['shrink'] for r in g):.3f}")

    print("=== (a) read-supported block gating (K=1 all-blocks vs K=3 read-supported) ===")
    print(f"  size-2 families assessed: {len(rows)}")
    report(edge, "TAU-EDGE families (recip_k1 in [0.30,0.55])")
    report(high, "HIGH-recip families (recip_k1>0.55, control)")
    print("\n  edge families detail (family: K1 -> K3, shrink):")
    for r in sorted(edge, key=lambda x: x["recip_k1"])[:15]:
        print(f"    fam{r['family']}: {r['recip_k1']:.3f} -> {r['recip_k3']:.3f}  "
              f"(delta {r['delta']:+.3f}, shrink {r['shrink']:.3f})")
    print("\n[+] wrote family_def_ri_gating.json")


if __name__ == "__main__":
    main()
