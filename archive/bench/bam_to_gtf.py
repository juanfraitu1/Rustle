#!/usr/bin/env python3
"""Convert a transcript->genome spliced BAM (minimap2 -ax splice) to a GTF: one transcript per PRIMARY
alignment, exons from the CIGAR (N ops split exons). For recovering de-novo transcript intron chains so they
can be gffcompare'd against the reference annotation."""
import sys, pysam

bam_path, gtf_path = sys.argv[1], sys.argv[2]
bam = pysam.AlignmentFile(bam_path, "rb")
n_in = n_out = 0
with open(gtf_path, "w") as out:
    for rec in bam:
        n_in += 1
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.cigartuples is None:
            continue
        chrom = rec.reference_name
        strand = "-" if rec.is_reverse else "+"
        tid = rec.query_name
        exons = []
        ref = rec.reference_start
        cur = ref
        for op, ln in rec.cigartuples:
            if op in (0, 7, 8, 2):       # M / = / X / D  -> consume reference within the exon
                ref += ln
            elif op == 3:                # N -> intron: close current exon, skip the gap
                exons.append((cur, ref))
                ref += ln
                cur = ref
            # I/S/H/P do not consume reference
        exons.append((cur, ref))
        t0, t1 = exons[0][0] + 1, exons[-1][1]
        out.write('%s\tmm2\ttranscript\t%d\t%d\t.\t%s\t.\ttranscript_id "%s"; gene_id "%s";\n'
                  % (chrom, t0, t1, strand, tid, tid))
        for i, (s, e) in enumerate(exons, 1):
            out.write('%s\tmm2\texon\t%d\t%d\t.\t%s\t.\ttranscript_id "%s"; gene_id "%s"; exon_number "%d";\n'
                      % (chrom, s + 1, e, strand, tid, tid, i))
        n_out += 1
print("alignments=%d  transcripts_written=%d" % (n_in, n_out))
