#!/usr/bin/env python3
"""Amendment 12 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): the `o3_candidates` input made from Amendment 7's held-out panel.

From `panel.json` (per family: `fam`, `mask` = the deleted copy, `keep` = the surviving copies' clean intervals, 0-based half-open) write:

  <out>.copies.tsv  one row per SURVIVING copy in the `P.fam.copies.tsv` layout (header = COPIES_HEADER of src/bin/mcl_families.rs):
                    tid / gene_id = the copy's name, start-end = the clean interval, one exon block, strand +, n_reads = distinct
                    primaries of the BAM overlapping the interval, source `panel`, core_hull NA, member_status `member`, locus = the
                    interval; columns no reader of this run consumes are NA
  <out>.copies.fa   the interval's sequence from the masked genome (upper case), header `>{family}|{idx}|{chrom}:{start}-{end}|+|nexon=1`
                    (the `parse_copies_fa` contract, src/rustle/vg_family/catalog_input.rs)
  <out>.regions     per family and chromosome, the merged span of its surviving copies +- 5 kb (clamped at 0): `{family}\t{chrom}:{lo}-{hi}`;
                    not consumed by the acceptance (an O2 run on the same inputs needs it)

The deleted copy is never written: the stage must find it from the reads alone.

    panel_to_copies.py --panel linktest/panel.json --bam linktest/R.bam --fasta linktest/masked.fa --out a12/A12
"""
import argparse
import collections
import json

import pysam

HEADER = ("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\tsource\tgene_id\tcore_hull\t"
          "sd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end")
PAD = 5000


def primaries(bam, chrom, start, end):
    return len({r.query_name for r in bam.fetch(chrom, start, end) if not (r.is_unmapped or r.is_secondary or r.is_supplementary)})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    panel = json.load(open(a.panel))
    bam, fa = pysam.AlignmentFile(a.bam), pysam.FastaFile(a.fasta)
    n_rows, n_zero = 0, 0
    with open(f"{a.out}.copies.tsv", "w") as t, open(f"{a.out}.copies.fa", "w") as f, open(f"{a.out}.regions", "w") as g:
        t.write(HEADER + "\n")
        for p in panel:
            fam = p["fam"]
            spans = collections.defaultdict(list)
            for idx, (chrom, start, end, gene) in enumerate(p["keep"]):
                seq = fa.fetch(chrom, start, end).upper()
                assert len(seq) == end - start, (fam, gene, len(seq), end - start)
                assert set(seq) != {"N"}, f"{fam} {gene}: a surviving copy is masked"
                n = primaries(bam, chrom, start, end)
                n_zero += n == 0
                row = [fam, idx, gene, chrom, start, end, 1, "+", n, f"{start}-{end}", "NA", "panel", gene, "NA", "NA", end - start, "NA",
                       "member", start, end]
                assert len(row) == len(HEADER.split("\t"))
                t.write("\t".join(map(str, row)) + "\n")
                f.write(f">{fam}|{idx}|{chrom}:{start}-{end}|+|nexon=1\n{seq}\n")
                spans[chrom].append((start, end))
                n_rows += 1
            for chrom, iv in sorted(spans.items()):
                iv.sort()
                lo, hi = iv[0]
                for s, e in iv[1:]:
                    if s <= hi + 2 * PAD:
                        hi = max(hi, e)
                    else:
                        g.write(f"{fam}\t{chrom}:{max(0, lo - PAD)}-{hi + PAD}\n"); lo, hi = s, e
                g.write(f"{fam}\t{chrom}:{max(0, lo - PAD)}-{hi + PAD}\n")
    print(f"families {len(panel)}, surviving copies {n_rows} ({n_zero} without a primary read); deleted copies (not written) {len(panel)}")


if __name__ == "__main__":
    main()
