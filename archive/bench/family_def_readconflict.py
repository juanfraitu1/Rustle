#!/usr/bin/env python3
"""Read-conflict (mutual-mappability) family definition vs the sequence-similarity definition, on a
Compara-labelled set.

The current pipeline defines a family by SEQUENCE similarity (contiguous-core >= 0.13, or minimizer-Jaccard).
That conflates true recent paralogs with DOMAIN-SHARERS (genes sharing one exon/domain). The proposed
OPERATIONAL definition: a family is a connected component of the READ-CONFLICT graph -- two loci are linked
iff reads actually cross-map between them with tied alignment scores (the unit the copy-assignment problem
operates on; Canzar's multimapping-conflict graph). This script measures, per labelled gene pair:

  read-conflict  : # reads with a primary alignment in one gene and an AS-tied SECONDARY in the other (a
                   genuine alternative placement at a DISTINCT locus -- not a single read overlapping two
                   nested annotations), and the mutual coverage (fraction of each gene's span involved).
  similarity     : canonical minimizer-Jaccard of the two gene rep sequences (the proxy the current
                   definition thresholds).

Truth labels (Compara, external): CONFIRMED_PAIRS = recent paralogs; DOMAIN_SHARER_PAIRS = share a domain
but are not paralogs. Claim under test: read-conflict separates the two classes where similarity need not.

Run with /home/juanfra/miniforge3/bin/python."""
import collections
import hashlib
import pysam

BAM = "/mnt/c/Users/jfris/Desktop/GGO.bam"
GENE_FA = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_recovery_eval/results/gene_rep.fa"
AS_TIE = 0.90  # a secondary counts as a conflict iff its AS >= AS_TIE * the read's best AS at the other locus

GENES = {  # gene -> (chrom, start, end) from GGO_genomic.gff
    "RFPL1": ("NC_086018.1", 27445867, 27472024), "RFPL2": ("NC_086018.1", 30208643, 30221673),
    "RFPL3": ("NC_086018.1", 30368560, 30376055), "APOBEC3D": ("NC_086018.1", 37023493, 37035944),
    "APOBEC3F": ("NC_086018.1", 37043278, 37058621), "RABL2B": ("NC_086018.1", 48818440, 48832011),
    "RABL2A": ("NC_073235.2", 15131653, 15147533), "PPARA": ("NC_086018.1", 44159462, 44252981),
    "CDPF1": ("NC_086018.1", 44250569, 44262082), "CASP8": ("NC_073235.2", 103850223, 103901051),
    "FLACC1": ("NC_073235.2", 103898904, 103977183), "CREB1": ("NC_073235.2", 110133100, 110220981),
    "METTL21A": ("NC_073235.2", 110198842, 110241136), "GPR39": ("NC_073235.2", 34964263, 35195860),
    "LYPD1": ("NC_073235.2", 35194087, 35224540), "GCA": ("NC_073235.2", 64785386, 64835227),
    "KCNH7": ("NC_073235.2", 64832333, 65309180), "ASDURF": ("NC_073235.2", 92106380, 92111318),
    "ASNSD1": ("NC_073235.2", 92111059, 92115776), "ZDHHC8": ("NC_086018.1", 17564727, 17581526),
    "CCDC188": ("NC_086018.1", 17579900, 17589084),
}
CONFIRMED = [("APOBEC3D", "APOBEC3F"), ("RFPL1", "RFPL2"), ("RFPL1", "RFPL3"), ("RFPL2", "RFPL3"),
             ("RABL2A", "RABL2B")]
DOMAIN = [("ASDURF", "ASNSD1"), ("CASP8", "FLACC1"), ("CCDC188", "ZDHHC8"), ("CDPF1", "PPARA"),
          ("CREB1", "METTL21A"), ("GCA", "KCNH7"), ("GPR39", "LYPD1")]


def load_fa(path):
    seqs, name, buf = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = ln[1:].strip().split()[0]; buf = []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def minimizers(seq, k=15, w=10):
    seq = seq.upper().encode()
    rc = seq.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]
    n = len(seq) - k + 1
    if n <= 0:
        return set()

    def h(s, i):
        km = s[i:i + k]
        return int.from_bytes(hashlib.blake2b(km, digest_size=8).digest(), "little")
    hs = [min(h(seq, i), h(rc, len(seq) - k - i)) for i in range(n)]
    out = set()
    for i in range(n - w + 1):
        out.add(min(hs[i:i + w]))
    return out


def jaccard(a, b):
    ma, mb = minimizers(a), minimizers(b)
    return len(ma & mb) / len(ma | mb) if (ma or mb) else 0.0


def placements(bam, gene):
    """qname -> list of (ref_start, ref_end, AS, is_secondary) for alignments overlapping the gene span."""
    c, s, e = GENES[gene]
    out = collections.defaultdict(list)
    for r in bam.fetch(c, s, e):
        if r.is_supplementary or r.is_unmapped:
            continue
        try:
            as_ = r.get_tag("AS")
        except KeyError:
            as_ = 0
        out[r.query_name].append((r.reference_start, r.reference_end, as_, r.is_secondary))
    return out


def conflict(bam, ga, gb):
    """reads with a placement in ga AND a DISTINCT placement in gb, AS-tied -> count + mutual coverage."""
    pa, pb = placements(bam, ga), placements(bam, gb)
    shared = set(pa) & set(pb)
    ca, sa, ea = GENES[ga]; cb, sb, eb = GENES[gb]
    cover_a, cover_b, n = [], [], 0
    for q in shared:
        # the read must have two DISTINCT alignment records (different locus), one per gene, AS-tied.
        best = None
        for (rs_a, re_a, as_a, _) in pa[q]:
            for (rs_b, re_b, as_b, _) in pb[q]:
                if abs(rs_a - rs_b) < 200 and ca == cb:   # same alignment overlapping both nested genes: skip
                    continue
                hi, lo = max(as_a, as_b), min(as_a, as_b)
                if hi > 0 and lo >= AS_TIE * hi:
                    best = (rs_a, re_a, rs_b, re_b)
        if best:
            n += 1
            cover_a.append((best[0], best[1])); cover_b.append((best[2], best[3]))

    def cov(ivs, s, e):
        if not ivs:
            return 0.0
        ivs = sorted((max(a, s), min(b, e)) for a, b in ivs)
        merged = 0; cs, ce = ivs[0]
        for a, b in ivs[1:]:
            if a <= ce:
                ce = max(ce, b)
            else:
                merged += ce - cs; cs, ce = a, b
        merged += ce - cs
        return merged / max(1, e - s)
    return n, cov(cover_a, sa, ea), cov(cover_b, sb, eb), len(pa), len(pb)


def main():
    fa = load_fa(GENE_FA)
    bam = pysam.AlignmentFile(BAM, "rb")
    rows = []
    for label, pairs in [("paralog", CONFIRMED), ("domain-sharer", DOMAIN)]:
        for ga, gb in pairs:
            n, cva, cvb, na, nb = conflict(bam, ga, gb)
            jac = jaccard(fa.get(ga, ""), fa.get(gb, "")) if ga in fa and gb in fa else float("nan")
            rows.append(dict(label=label, a=ga, b=gb, n_conflict=n, mut_cov=min(cva, cvb),
                             cov_a=cva, cov_b=cvb, reads_a=na, reads_b=nb, jaccard=jac))
            print(f"{label:13} {ga:9}~{gb:9}  conflict_reads={n:3}  mut_cov={min(cva,cvb):.3f} "
                  f"(a={cva:.2f} b={cvb:.2f})  jaccard={jac:.3f}  reads={na}/{nb}")
    with open("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_readconflict.tsv", "w") as fh:
        fh.write("label\tgene_a\tgene_b\tn_conflict\tmut_cov\tcov_a\tcov_b\treads_a\treads_b\tjaccard\n")
        for r in rows:
            fh.write(f"{r['label']}\t{r['a']}\t{r['b']}\t{r['n_conflict']}\t{r['mut_cov']:.4f}\t{r['cov_a']:.4f}"
                     f"\t{r['cov_b']:.4f}\t{r['reads_a']}\t{r['reads_b']}\t{r['jaccard']:.4f}\n")
    print("[wrote bench/family_def_readconflict.tsv]")


if __name__ == "__main__":
    main()
