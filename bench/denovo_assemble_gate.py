#!/usr/bin/env python3
"""Gate the de-novo read-coherence skeletons (general-purpose assembly) + extract their sequences.
This is the GENERAL assembler output — ALL transcripts (family or not), read-coherent, annotation-free.
The family/copy layer (next stage) is then applied only to the multi-copy SUBSET ('areas of focus');
single-copy transcripts are assembled and kept as-is.

Gate (read-coherence precision filter): >= MIN_READS reads AND all introns canonical (GT-AG/GC-AG/AT-AC
or the reverse-strand motif, consistent orientation). Run with /home/juanfra/miniforge3/bin/python.
"""
import pysam

SK = "/home/juanfra/winloci_scratch/denovo_skeletons.tsv"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
OUT_FA = "/home/juanfra/winloci_scratch/denovo_transcripts.fa"
OUT_META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
MIN_READS = 3
MAX_SPAN = 3_000_000   # reject skeletons spanning > this genomic distance (> largest real gene DMD ~2.3Mb;
                       #   read-coherence skeletons can chain across a whole chromosome -> 236Mb artifacts)
MAX_SPLICED = 300_000  # reject implausibly long spliced products (titin mRNA ~109kb) -- catches giant
                       #   terminal "exons" from artifactual long-span skeletons
PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}


def main():
    fa = pysam.FastaFile(FASTA)
    n_in = n_kept = n_multiexon = n_span_drop = 0
    fh = open(OUT_FA, "w"); mh = open(OUT_META, "w")
    mh.write("id\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\n")
    for line in open(SK):
        if line.startswith("chrom\t"):
            continue
        c, s, e, ne, nr, introns = line.rstrip("\n").split("\t")
        n_in += 1
        nr = int(nr); s = int(s); e = int(e)
        if nr < MIN_READS:
            continue
        if e - s > MAX_SPAN:           # drop chromosome-spanning artifact skeletons
            n_span_drop += 1
            continue
        ivs = [tuple(map(int, x.split("-"))) for x in introns.split(";")] if introns else []
        # canonical-junction gate (consistent strand)
        strand = None; ok = True
        for d, a in ivs:
            don = fa.fetch(c, d, d + 2).upper(); acc = fa.fetch(c, a - 2, a).upper()
            st = "+" if (don, acc) in PLUS else ("-" if (don, acc) in MINUS else None)
            if st is None or (strand and st != strand):
                ok = False; break
            strand = st
        if not ok:
            continue
        # exon coords from the intron chain + terminal extents
        exons = []
        prev = s
        for d, a in ivs:
            exons.append((prev, d)); prev = a
        exons.append((prev, e))
        seq = "".join(fa.fetch(c, xs, xe).upper() for xs, xe in exons)
        if len(seq) < 100 or len(seq) > MAX_SPLICED:   # drop too-short and implausibly-long spliced products
            continue
        if strand == "-":
            seq = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
        tid = f"DN_{c}_{s}_{len(ivs)+1}"
        fh.write(f">{tid}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")
        mh.write(f"{tid}\t{c}\t{s}\t{e}\t{strand or '+'}\t{len(ivs)+1}\t{nr}\n")
        n_kept += 1
        if len(ivs) >= 1:
            n_multiexon += 1
    fh.close(); mh.close()
    print(f"de-novo skeletons in: {n_in:,}  (dropped {n_span_drop:,} chromosome-spanning artifacts > {MAX_SPAN/1e6:.0f}Mb)")
    print(f"GATED assembled transcripts (>= {MIN_READS} reads + all-canonical junctions + <= {MAX_SPLICED/1e3:.0f}kb spliced): "
          f"{n_kept:,} (multi-exon {n_multiexon:,})")
    print(f"  = the general-purpose read-coherence assembly (family or not). -> {OUT_FA}")


if __name__ == "__main__":
    main()
