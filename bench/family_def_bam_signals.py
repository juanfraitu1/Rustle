#!/usr/bin/env python3
"""family_def_bam_signals.py — is there discriminating signal in the BAM beyond `de`?

~R currently uses only `de` (gap-compressed divergence). The new BAM also carries NM (raw
edit distance = indel burden, which `de` gap-compresses away) and rl (length of query in
repetitive seeds = a repeat-mediated-bridge flag). Question: do real-copy cross-mappings
and bridge cross-mappings separate on NM or rl at the READ level — i.e. could a read-level
signal pre-empt some of what ~B prunes at the copy level?

Per chromosome (new BAM): scan cross-mapping reads, keep per-placement (de, NM, rl), build
~R edges, classify each edge ~B-copy vs ~B-bridge (reciprocal coverage of exon-union models),
and compare the cross-mapping reads' de / NM-rate / rl between the two classes.
Run: /home/juanfra/miniforge3/bin/python bench/family_def_bam_signals.py NC_073240.2
"""
import collections
import os
import subprocess
import sys
import statistics as st

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import de_of, best_gene, ref_span, DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import (load_denovo, build_model, recip_cov, load_models,
                                        GENOME, OLD, COV_MIN)

NEW = "/home/juanfra/winloci_scratch/GGO_mm.bam"
SAM = "/home/juanfra/miniforge3/bin/samtools"


def tag(fields, pre):
    for f in fields:
        if f.startswith(pre):
            return f.strip()
    return None


def scan_sig(by_chrom, region):
    """cross-mapping reads on `region`: qname -> {gene: (de, nm, rl, alnlen)}."""
    mm = collections.defaultdict(dict)
    for flag in ("-f", "0x100"), ("-F", "0x900"):
        p = subprocess.Popen([SAM, "view", flag[0], flag[1], NEW, region],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 9:
                continue
            tags = f[11:]
            de = de_of(tags)
            if de is None or de > DE_MAX:
                continue
            g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
            if g is None:
                continue
            nm = tag(tags, "NM:i:"); rl = tag(tags, "rl:i:")
            nm = int(nm[5:]) if nm else 0
            rl = int(rl[5:]) if rl else 0
            aln = ref_span(f[5])
            d = mm[f[0]]
            if g not in d or de < d[g][0]:
                d[g] = (de, nm, rl, aln)
        p.wait()
    return mm


def main():
    region = sys.argv[1] if len(sys.argv) > 1 else "NC_073240.2"
    by_chrom, info = load_denovo({region})
    models = load_models()
    bam_old = pysam.AlignmentFile(OLD, "rb"); genome = pysam.FastaFile(GENOME)
    mm = scan_sig(by_chrom, region)

    # ~R edges (de-tie quorum), same as the definition
    ev = collections.defaultdict(list)   # (ga,gb) -> list of reads' (placement_a, placement_b)
    for genes in mm.values():
        if len(genes) < 2:
            continue
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, pa = items[a]
            for b in range(a + 1, len(items)):
                gb, pb = items[b]
                if abs(pa[0] - pb[0]) <= DELTA and max(pa[0], pb[0]) <= DE_MAX:
                    ev[(ga, gb)].append((pa, pb))
    edges = [(ga, gb, recs) for (ga, gb), recs in ev.items() if len(recs) >= MIN_READS]

    # resolve models, classify edges ~B
    loci = {g for (ga, gb, _) in edges for g in (ga, gb)}
    for v in loci:
        if not models.get(v) and v in info:
            c, s, e, _ = info[v]
            models[v] = build_model(bam_old, genome, c, s, e)
    copy_recs, bridge_recs = [], []
    n_copy = n_bridge = 0
    for (ga, gb, recs) in edges:
        if not models.get(ga) or not models.get(gb):
            continue
        is_copy = recip_cov(models[ga], models[gb]) >= COV_MIN
        if is_copy:
            n_copy += 1
        else:
            n_bridge += 1
        for (pa, pb) in recs:
            row = ((pa[0] + pb[0]) / 2,                       # mean de
                   (pa[1] / max(pa[3], 1) + pb[1] / max(pb[3], 1)) / 2,   # mean NM-rate
                   (pa[2] + pb[2]) / 2,                       # mean rl (repeat-seed length)
                   (pa[2] / max(pa[3], 1) + pb[2] / max(pb[3], 1)) / 2)   # mean rl-fraction
            (copy_recs if is_copy else bridge_recs).append(row)
    bam_old.close(); genome.close()

    def summ(rows, name):
        if not rows:
            print(f"  {name}: (none)"); return
        de = [r[0] for r in rows]; nmr = [r[1] for r in rows]
        rl = [r[2] for r in rows]; rlf = [r[3] for r in rows]
        print(f"  {name} (n={len(rows)} read-placements):")
        print(f"    de:        median={st.median(de):.4f}")
        print(f"    NM-rate:   median={st.median(nmr):.4f}  (edit dist / aligned bp)")
        print(f"    rl(bp):    median={st.median(rl):.0f}  mean={st.mean(rl):.0f}")
        print(f"    rl-frac:   median={st.median(rlf):.3f}  (repeat-seed bp / aligned bp)")

    print(f"=== BAM-signal discrimination on {region}: ~B-copy vs ~B-bridge edges ===")
    print(f"  edges: {n_copy} copy, {n_bridge} bridge")
    summ(copy_recs, "REAL-COPY cross-mappings")
    summ(bridge_recs, "BRIDGE cross-mappings")
    # separation verdict
    if copy_recs and bridge_recs:
        for idx, lab in [(0, "de"), (1, "NM-rate"), (3, "rl-frac")]:
            c = st.median([r[idx] for r in copy_recs]); b = st.median([r[idx] for r in bridge_recs])
            print(f"  [{lab}] copy={c:.4f} vs bridge={b:.4f}  ratio={b/max(c,1e-9):.1f}x"
                  f"  -> {'DISCRIMINATES' if abs(b-c)/max(c,b,1e-9) > 0.5 else 'weak'}")


if __name__ == "__main__":
    main()
