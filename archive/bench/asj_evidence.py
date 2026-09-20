#!/usr/bin/env python3
"""Per-call evidence for an allele-specific junction: per-molecule 2x2 (anchor allele x junction
usage), the junction's splice dinucleotides (canonical GT-AG?), the anchor's reference context, and
the anchor->junction distance. Confirms the mechanism (splice-disrupting SNP vs long-range isoform
switch) from the molecules directly. Run with python3 (pysam). Use --row N into asj_calls_verified.tsv.
"""
import argparse
import os

import pysam

CALLS = os.path.join(os.path.dirname(__file__), "asj_calls_verified.tsv")
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"


def allele_at(aln, ref_pos):
    qpos, rpos = 0, aln.reference_start
    seq = aln.query_sequence
    for op, length in aln.cigartuples:
        if op in (0, 7, 8):
            if rpos <= ref_pos < rpos + length:
                return seq[qpos + (ref_pos - rpos)]
            qpos += length; rpos += length
        elif op in (2, 3):
            if rpos <= ref_pos < rpos + length:
                return None
            rpos += length
        elif op in (1, 4):
            qpos += length
    return None


def introns(aln):
    out = []; rpos = aln.reference_start
    for op, length in aln.cigartuples:
        if op == 3:
            out.append((rpos, rpos + length)); rpos += length
        elif op in (0, 2, 7, 8):
            rpos += length
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, required=True, help="1-based data row in asj_calls_verified.tsv")
    args = ap.parse_args()
    with open(CALLS) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        lines = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]
    r = lines[args.row - 1]
    chrom = r["chrom"]; anc = int(r["anchor"]); d = int(r["donor"]); a = int(r["acceptor"])
    ax, ay = r["alleleX"], r["alleleY"]
    bam = pysam.AlignmentFile(BAM, "rb"); fa = pysam.FastaFile(FASTA)

    print(f"# {r['gene']} {chrom} anchor={anc} {ax}/{ay} junction={d}-{a} "
          f"dPSI={r['dPSI']} q={r['q']} type={r['anchor_type']} dist={r['dist']} frac_mq0={r['frac_mq0']}")
    donor_di = fa.fetch(chrom, d, d + 2).upper()
    acc_di = fa.fetch(chrom, a - 2, a).upper()
    print(f"# splice dinucleotides: donor={donor_di} acceptor={acc_di} "
          f"({'canonical GT-AG (+)' if donor_di=='GT' and acc_di=='AG' else 'canonical CT-AC (-)' if donor_di=='CT' and acc_di=='AC' else 'non-canonical/other'})")
    print(f"# anchor ref context [{anc-5}..{anc+5}]: {fa.fetch(chrom, anc-5, anc+6).upper()} "
          f"(anchor base = pos {anc})")

    cnt = {ax: [0, 0], ay: [0, 0]}  # allele -> [used, spanning]
    for aln in bam.fetch(chrom, min(anc, d), max(anc, a) + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_sequence is None:
            continue
        if aln.reference_start > anc or aln.reference_end <= anc:
            continue
        b = allele_at(aln, anc)
        if b not in (ax, ay):
            continue
        if aln.reference_start <= d and aln.reference_end >= a:
            cnt[b][1] += 1
            if (d, a) in introns(aln):
                cnt[b][0] += 1
    print(f"# per-molecule 2x2 (allele x junction usage), reads spanning anchor AND junction:")
    for al in (ax, ay):
        u, s = cnt[al]
        print(f"#   allele {al}: uses junction {u}/{s}  PSI={u/s:.2f}" if s else f"#   allele {al}: 0 spanning")
    print("# => one allele uses the junction, the other does not"
          if all(cnt[al][1] for al in (ax, ay)) and
             abs((cnt[ax][0]/cnt[ax][1]) - (cnt[ay][0]/cnt[ay][1])) >= 0.3
          else "# => check power/effect")


if __name__ == "__main__":
    main()
