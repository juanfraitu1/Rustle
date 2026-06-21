"""Tier-2 resolver: assign K=0 cross-mapping reads via copy-specific splice sites. Intronic divergence at an
exon-intron boundary makes a junction canonical at one copy but degraded at the other; a read whose junctions
are fully canonical at copy X while its homologous alignment at copy Y carries a degraded junction came from X."""
import pysam

# canonical splice dinucleotide pairs, strand-agnostic (each major class + its reverse-complement):
#   GT-AG / CT-AC ;  GC-AG / CT-GC ;  AT-AC / GT-AT
CANONICAL = {("GT", "AG"), ("CT", "AC"), ("GC", "AG"), ("CT", "GC"), ("AT", "AC"), ("GT", "AT")}

def _introns(read):
    """Genomic introns as (donor, acceptor) 0-based half-open: donor=intron first base, acceptor=intron end."""
    out, pos = [], read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 2, 7, 8):       # M/D/=/X consume reference
            pos += ln
        elif op == 3:                 # N = intron
            out.append((pos, pos + ln)); pos += ln
    return out

def _splice(fasta, chrom, d, a):
    return (fasta.fetch(chrom, d, d + 2).upper(), fasta.fetch(chrom, a - 2, a).upper())

def _canon(introns, fasta, chrom):
    return sum(1 for (d, a) in introns if _splice(fasta, chrom, d, a) in CANONICAL)

def resolve_pair(bam_path, fasta_path, locusA, locusB):
    bam = pysam.AlignmentFile(bam_path, "rb"); fasta = pysam.FastaFile(fasta_path)
    cA, sA, eA = locusA; cB, sB, eB = locusB

    def reads_at(c, s, e):
        d = {}
        for r in bam.fetch(c, s, e):
            if r.is_unmapped or r.is_supplementary:
                continue
            if min(r.reference_end, e) - max(r.reference_start, s) < 200:
                continue
            d.setdefault(r.query_name, r)   # one alignment per query at this locus
        return d

    A, B = reads_at(cA, sA, eA), reads_at(cB, sB, eB)
    common = set(A) & set(B)
    per_read, distinguishing, n_junction, n_resolved = {}, set(), 0, 0
    for q in common:
        ia, ib = _introns(A[q]), _introns(B[q])
        if not ia and not ib:
            continue
        n_junction += 1
        ca, cb = _canon(ia, fasta, cA), _canon(ib, fasta, cB)
        if ia and ca == len(ia) and ib and cb < len(ib):        # A fully canonical, B degraded -> copy A
            per_read[q] = "A"; n_resolved += 1
            for (d, a) in ib:
                if _splice(fasta, cB, d, a) not in CANONICAL:
                    distinguishing.add((cB, d, a))
        elif ib and cb == len(ib) and ia and ca < len(ia):      # B fully canonical, A degraded -> copy B
            per_read[q] = "B"; n_resolved += 1
            for (d, a) in ia:
                if _splice(fasta, cA, d, a) not in CANONICAL:
                    distinguishing.add((cA, d, a))
    return dict(n_reads=len(common), n_junction_reads=n_junction, n_resolved=n_resolved,
                resolved_fraction=(n_resolved / n_junction if n_junction else 0.0),
                distinguishing_junctions=sorted(distinguishing), per_read=per_read)


PANEL = {
    "pair0": (("NC_073247.2", 161251228, 161257000), ("NC_073247.2", 161458538, 161464324)),
    "pair2": (("NC_073247.2", 164381222, 164384848), ("NC_073247.2", 164442447, 164446101)),
    "pair3": (("NC_073247.2", 164397061, 164401095), ("NC_073247.2", 164426194, 164430228)),
}

def check_resolver():
    BAM   = "/home/juanfra/winloci_scratch/GGO.bam"
    FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
    res = {k: resolve_pair(BAM, FASTA, a, b) for k, (a, b) in PANEL.items()}
    for k, r in res.items():
        print(f"    {k}: junction_reads={r['n_junction_reads']} resolved={r['n_resolved']} "
              f"({r['resolved_fraction']:.2f}) distinguishing={len(r['distinguishing_junctions'])}")
    # pair0 carries the copy-specific 5' splice site -> a meaningful fraction resolves;
    # pair2/pair3 are splice-identical (canonical at both copies) -> ~0 (also guards against a strand bug).
    assert res["pair0"]["resolved_fraction"] >= 0.15, f"pair0 should resolve a meaningful fraction, got {res['pair0']['resolved_fraction']:.2f}"
    assert res["pair2"]["resolved_fraction"] <= 0.05, f"pair2 is splice-identical, expected ~0, got {res['pair2']['resolved_fraction']:.2f}"
    assert res["pair3"]["resolved_fraction"] <= 0.05, f"pair3 is splice-identical, expected ~0, got {res['pair3']['resolved_fraction']:.2f}"
    print(f"OK  - splice resolver: pair0 {res['pair0']['resolved_fraction']:.2f} resolved "
          f"({len(res['pair0']['distinguishing_junctions'])} distinguishing junctions); "
          f"pair2 {res['pair2']['resolved_fraction']:.2f}, pair3 {res['pair3']['resolved_fraction']:.2f} (splice-identical)")
    return res

if __name__ == "__main__":
    check_resolver()
