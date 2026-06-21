"""Splice-divergence NEGATIVE-RESULT probe (Tier-2 collapses into Tier-3).

For K=0 co-located inverted-dup copies, intronic divergence near a splice site (e.g. the MAGEA pair0 3 bp
intronic indel) is REAL in the reference but is NOT usable for per-read copy assignment: minimap2's spliced
aligner independently snaps every junction to the nearest canonical GT-AG site at EACH copy, so a read is
fully canonical at both copies (intron lengths differing only by the indel, e.g. 3702 nt at A vs 3705 nt at B).
The indel sits inside the long intron -- absent from the spliced mature read, accommodated without an exonic
edit -- so NM stays tied (these are K=0 by construction) and there is no per-read direction signal.

This probe quantifies that on the MAGEA panel:
- resolved_fraction        (per-read canonicality direction)  -> ~0 for ALL pairs incl. pair0 (the negative
                                                                 result: aligner snapping masks the divergence)
- chain_divergent_fraction (read's A-intron-length multiset != B's) -> HIGH at pair0 (divergence IS detectable)
                                                                 but ~0 at splice-identical pair2/pair3
=> detectable yet NON-DIRECTIONAL: the K=0 residual has no per-read splice rescue; it collapses into Tier-3."""
import pysam

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

def probe_pair(bam_path, fasta_path, locusA, locusB):
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
    distinguishing, n_junction, n_resolved, n_chain_div = set(), 0, 0, 0
    for q in common:
        ia, ib = _introns(A[q]), _introns(B[q])
        if not ia and not ib:
            continue
        n_junction += 1
        # (1) DETECTABLE divergence: intron-length multiset differs between the two alignments.
        #     Reflection preserves intron lengths, so a length difference is a genuine splice divergence.
        if sorted(a - d for d, a in ia) != sorted(a - d for d, a in ib):
            n_chain_div += 1
        # (2) per-read DIRECTION via canonicality -> ~0: minimap2 snaps BOTH copies to canonical sites.
        ca, cb = _canon(ia, fasta, cA), _canon(ib, fasta, cB)
        if ia and ca == len(ia) and ib and cb < len(ib):        # A fully canonical, B degraded -> copy A
            n_resolved += 1
            for (d, a) in ib:
                if _splice(fasta, cB, d, a) not in CANONICAL:
                    distinguishing.add((cB, d, a))
        elif ib and cb == len(ib) and ia and ca < len(ia):      # B fully canonical, A degraded -> copy B
            n_resolved += 1
            for (d, a) in ia:
                if _splice(fasta, cA, d, a) not in CANONICAL:
                    distinguishing.add((cA, d, a))
    return dict(n_reads=len(common), n_junction_reads=n_junction, n_resolved=n_resolved,
                resolved_fraction=(n_resolved / n_junction if n_junction else 0.0),
                n_chain_divergent=n_chain_div,
                chain_divergent_fraction=(n_chain_div / n_junction if n_junction else 0.0),
                distinguishing_junctions=sorted(distinguishing))


PANEL = {
    "pair0": (("NC_073247.2", 161251228, 161257000), ("NC_073247.2", 161458538, 161464324)),
    "pair2": (("NC_073247.2", 164381222, 164384848), ("NC_073247.2", 164442447, 164446101)),
    "pair3": (("NC_073247.2", 164397061, 164401095), ("NC_073247.2", 164426194, 164430228)),
}

def check_probe():
    BAM   = "/home/juanfra/winloci_scratch/GGO.bam"
    FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
    res = {k: probe_pair(BAM, FASTA, a, b) for k, (a, b) in PANEL.items()}
    for k, r in res.items():
        print(f"    {k}: junction_reads={r['n_junction_reads']} "
              f"per-read-resolved={r['resolved_fraction']:.2f} "
              f"chain-divergent={r['chain_divergent_fraction']:.2f} "
              f"distinguishing_junc={len(r['distinguishing_junctions'])}")
    # NEGATIVE RESULT: per-read splice resolution is ~0 for ALL pairs incl. pair0 -- aligner junction-snapping
    # masks the intronic divergence (the prior "33%" was a reference-level junction count, not per-read).
    for k in PANEL:
        assert res[k]["resolved_fraction"] <= 0.05, f"{k} per-read resolved should be ~0 (aligner-masked), got {res[k]['resolved_fraction']:.2f}"
    # ...but the divergence IS real and DETECTABLE as intron-length divergence at pair0, and ABSENT at the
    # splice-identical pairs -- detectable yet non-directional (so it cannot rescue per-read assignment).
    assert res["pair0"]["chain_divergent_fraction"] >= 0.15, f"pair0 divergence should be detectable, got {res['pair0']['chain_divergent_fraction']:.2f}"
    assert res["pair2"]["chain_divergent_fraction"] <= 0.05, f"pair2 splice-identical, got {res['pair2']['chain_divergent_fraction']:.2f}"
    assert res["pair3"]["chain_divergent_fraction"] <= 0.05, f"pair3 splice-identical, got {res['pair3']['chain_divergent_fraction']:.2f}"
    print(f"OK  - splice probe (NEGATIVE per-read): pair0 per-read-resolved={res['pair0']['resolved_fraction']:.2f} "
          f"but chain-divergent={res['pair0']['chain_divergent_fraction']:.2f} (detectable, non-directional); "
          f"pair2/3 resolved~0 chain-divergent~0 (splice-identical). Tier-2 -> Tier-3.")
    return res

if __name__ == "__main__":
    check_probe()
