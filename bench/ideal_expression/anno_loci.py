#!/usr/bin/env python3
"""G6 positive control of docs/PREREG_ideal_expression_2026-10-06.md: the canonicalized annotation written as an assembled GTF (every simulated transcript, equal read counts), so the
driver's `families` stage (mcl_families --from-gtf, the shipped flags) can be run on PERFECT loci and the same scorer applied: E2 and E3 must be 100% on R and its E4 is the family ceiling.

    anno_loci.py TRUTHPREFIX OUTPREFIX     writes OUTPREFIX.gtf and OUTPREFIX.families.gtf (identical; the second newer, as the driver's f1v2 guard requires)

A transcript's gene_id is its gene (several isoforms of a gene are one locus); transcript_id is the annotation transcript id; `reads "10"` on every transcript.
"""
import collections
import csv
import sys
import time

truth, out = sys.argv[1], sys.argv[2]
rows = [r for r in csv.DictReader(open(truth + ".transcripts.tsv"), delimiter="\t") if int(r["simulated"])]
lines = []
for r in sorted(rows, key=lambda r: (int(r["canon_blocks"].split(",")[0].split("-")[0]), r["transcript"])):
    bl = [tuple(map(int, b.split("-"))) for b in r["canon_blocks"].split(",")]
    gid, tid = r["gene"], r["transcript"]
    lines.append(f'{r["chrom"]}\trustle\ttranscript\t{bl[0][0] + 1}\t{bl[-1][1]}\t.\t{r["strand"]}\t.\tgene_id "{gid}"; transcript_id "{tid}"; reads "10"; cov "10.0"; TPM "1.0";')
    for i, (s, e) in enumerate(bl, 1):
        lines.append(f'{r["chrom"]}\trustle\texon\t{s + 1}\t{e}\t.\t{r["strand"]}\t.\tgene_id "{gid}"; transcript_id "{tid}"; exon_number "{i}";')
open(out + ".gtf", "w").write("\n".join(lines) + "\n")
time.sleep(1.1)
open(out + ".families.gtf", "w").write("\n".join(lines) + "\n")
print(f"{len(rows)} transcripts written to {out}.gtf and {out}.families.gtf")
