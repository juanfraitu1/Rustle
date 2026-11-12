#!/usr/bin/env python3
"""r1643 (docs/PREREG_exon_space_multi_isoform_2026-11-12.md): several sequences per catalog copy, named
FAMILY|COPY|KIND|ID so exon_place.py reads them unchanged. KIND: cat (catalog FASTA record), tx (assembly transcripts
on this copy's locus only), pre (genomic span). With a 5th argument `other` (r1647) every remaining assembly transcript is written as OTHER|0|tx|ID.
Usage: multi_iso_ref.py COPIES.tsv COPIES.fa ASSEMBLY.gtf GENOME.fa [other] > OUT.fa"""
import bisect, collections, csv, re, sys
import pysam

def main(copies, copies_fa, gtf, genome, other=None):
    fa = pysam.FastaFile(genome); by = collections.defaultdict(list); out = sys.stdout
    for r in csv.DictReader(open(copies), delimiter="\t"):  # spans are 0-based half-open (exon lengths sum to the FASTA record)
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"], r["copy_idx"]))
    for name, seq in ((x.name, x.sequence) for x in pysam.FastxFile(copies_fa)):
        f = name.split("|"); out.write(f">{f[0]}|{f[1]}|cat|0\n{seq}\n")
    for ch, v in by.items():
        v.sort()
        for s, e, fam, i in v: out.write(f">{fam}|{i}|pre|0\n{fa.fetch(ch, s, e)}\n")
    maxlen = {ch: max(e - s for s, e, _, _ in v) for ch, v in by.items()}; starts = {ch: [x[0] for x in v] for ch, v in by.items()}
    tx = collections.defaultdict(list)
    for line in open(gtf):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "exon" and (f[0] in by or other): tx[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    n = collections.Counter()
    for (ch, tid), ex in tx.items():
        ex.sort(); s, e = ex[0][0], ex[-1][1]; v = by.get(ch, []); hit = []
        if v:
            lo = bisect.bisect_left(starts[ch], s - maxlen[ch]); hi = bisect.bisect_left(starts[ch], e)
            hit = [c for c in v[lo:hi] if c[0] < e and s < c[1]]
        two = any(a[1] <= b[0] or b[1] <= a[0] for a in hit for b in hit)  # spans two loci: would tie them
        if two: n["two_loci_skipped"] += 1
        if not hit or two:
            if other: n["other"] += 1; out.write(f">OTHER|0|tx|{tid}\n" + "".join(fa.fetch(ch, a, b) for a, b in ex) + "\n")
            continue
        seq = "".join(fa.fetch(ch, a, b) for a, b in ex); n["tx"] += 1
        for _, _, fam, i in hit: out.write(f">{fam}|{i}|tx|{tid}\n{seq}\n")
    print(dict(n), file=sys.stderr)

if __name__ == "__main__": main(*sys.argv[1:6])
