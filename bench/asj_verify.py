#!/usr/bin/env python3
"""Deterministic confound control for the allele-specific-junction calls (task c style).

Two confounds the raw FDR cannot see:
  COLLAPSED-PARALOG MASQUERADE: at a multimapping locus the two 'alleles' are paralog copies, so
    'allele-specific' junctions are really copy-specific isoforms. Control: MAPQ-0 fraction at the
    anchor (uniquely mapped => genuine within-gene heterozygosity).
  MECHANISM PLAUSIBILITY: a SNP AT/near a splice site (<=~100bp) is a canonical splice-disrupting
    variant (cleanest ASS); a DISTAL anchor flipping a junction fully (dPSI=1.0) is more likely two
    distinct transcript populations than allele-modulated splicing of one gene.

Adds frac_mq0 + anchor-junction distance to each call; reports the confound-controlled high-confidence
set (transversion anchor AND uniquely mapped). Writes bench/asj_calls_verified.tsv + prints summary.
"""
import os
from collections import defaultdict

import pysam

CALLS = os.path.join(os.path.dirname(__file__), "asj_calls.tsv")
OUT = os.path.join(os.path.dirname(__file__), "asj_calls_verified.tsv")
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
MQ0_THR = 0.10
PROX = 100   # anchor within this many bp of the junction = splice-proximal (mechanism plausible)


def frac_mq0_at(bam, chrom, pos):
    n = z = 0
    for aln in bam.fetch(chrom, pos, pos + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        n += 1
        if aln.mapping_quality == 0:
            z += 1
        if n >= 600:
            break
    return (z / n) if n else 0.0


def main():
    rows = []
    with open(CALLS) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    bam = pysam.AlignmentFile(BAM, "rb")

    # MAPQ per unique anchor (dedup; many calls share an anchor)
    mq0 = {}
    for r in rows:
        key = (r["chrom"], int(r["anchor"]))
        if key not in mq0:
            mq0[key] = frac_mq0_at(bam, key[0], key[1])

    for r in rows:
        anc = int(r["anchor"]); d = int(r["donor"]); a = int(r["acceptor"])
        r["dist"] = min(abs(anc - d), abs(anc - a))
        r["frac_mq0"] = mq0[(r["chrom"], anc)]

    transversion = [r for r in rows if r["anchor_type"] == "transversion"]
    uniq = [r for r in rows if r["frac_mq0"] < MQ0_THR]
    hc = [r for r in rows if r["anchor_type"] == "transversion" and r["frac_mq0"] < MQ0_THR]
    prox = [r for r in hc if r["dist"] <= PROX]
    full = [r for r in rows if float(r["dPSI"]) >= 0.999]
    full_collapsed = [r for r in full if r["frac_mq0"] >= MQ0_THR]

    def genes(rs): return len({r["gene"] for r in rs})
    print(f"raw ASJ calls: {len(rows)} ({genes(rows)} genes)")
    print(f"uniquely-mapped (frac_mq0<{MQ0_THR}): {len(uniq)} ({genes(uniq)} genes) "
          f"| multimapping (collapsed-paralog suspect): {len(rows)-len(uniq)}")
    print(f"transversion (genetic): {len(transversion)} | uniquely-mapped: {sum(1 for r in transversion if r['frac_mq0']<MQ0_THR)}")
    print(f"** HIGH-CONFIDENCE (transversion AND uniquely-mapped): {len(hc)} junctions / {genes(hc)} genes **")
    print(f"   of HC, splice-proximal (anchor<= {PROX}bp from junction = candidate splice-disrupting): {len(prox)}")
    print(f"full-switch dPSI=1.0: {len(full)} total; of those multimapping (collapsed suspect): {len(full_collapsed)}")
    # the suspicious top hit
    loc = [r for r in rows if r["gene"] == "LOC101141065"]
    if loc:
        print(f"LOC101141065 (suspicious top hit): {len(loc)} junctions, frac_mq0={loc[0]['frac_mq0']:.2f}, "
              f"dists={sorted(set(r['dist'] for r in loc))[:6]} -> "
              f"{'COLLAPSED/multimap artifact' if loc[0]['frac_mq0']>=MQ0_THR else 'uniquely mapped'}")
    print("\nHIGH-CONFIDENCE examples (transversion + uniquely-mapped, by dPSI):")
    for r in sorted(hc, key=lambda r: (-float(r["dPSI"]), float(r["q"])))[:15]:
        print(f"  {r['gene']:14s} {r['chrom']} anc={r['anchor']} {r['alleleX']}/{r['alleleY']} "
              f"junc={r['donor']}-{r['acceptor']} PSI {r['psiX']}vs{r['psiY']} dPSI={r['dPSI']} "
              f"dist={r['dist']} mq0={r['frac_mq0']:.2f}")

    with open(OUT, "w") as fh:
        cols = ["gene", "chrom", "anchor", "alleleX", "alleleY", "donor", "acceptor",
                "psiX", "psiY", "dPSI", "q", "anchor_type", "dist", "frac_mq0"]
        fh.write("\t".join(cols) + "\thigh_confidence\n")
        for r in rows:
            h = int(r["anchor_type"] == "transversion" and r["frac_mq0"] < MQ0_THR)
            fh.write("\t".join(str(r[c]) for c in cols) + f"\t{h}\n")
    print(f"\n[wrote {OUT}]")


if __name__ == "__main__":
    main()
