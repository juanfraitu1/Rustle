#!/usr/bin/env python3
"""Promote the 18 strong hidden-copy FLAGS to catalogued reference-absent COPIES (or reject).
Per candidate: (1) isolate the alt-haplotype reads (the second haplotype), (2) POA consensus
(pyabpoa) = the collapsed copy's transcript, (3) longest ORF, (4) DISTINCT-PARALOG homology
(minimap2 consensus vs the family's reference copy models: must map = homologous, but identity to
the nearest reference copy < ~98% = distinct/private), (5) depth-doubling (alt fraction + locus
depth). Emits a table + consensus FASTA + json. DETECT+characterise; promotion is evidence, not
fabrication."""
import sys, json, os, subprocess, tempfile
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam, pyabpoa

BAL_LO, BAL_HI, MIN_DEPTH, SHARE_HI = 0.20, 0.60, 8, 0.60
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/promoted"
os.makedirs(OUT, exist_ok=True)
COMP = str.maketrans("ACGTN", "TGCAN")


def alts_in_span(aln, s, e):
    out = []; rpos = aln.reference_start
    for op, ln in (aln.cigartuples or []):
        if op in (7, 0): rpos += ln
        elif op == 8:
            for i in range(ln):
                if s <= rpos + i < e: out.append(rpos + i)
            rpos += ln
        elif op in (2, 3): rpos += ln
    return out


def alt_haplotype_reads(bam, chrom, s, e):
    """Return (alt_read_seqs, n_primary, alt_frac): the reads forming the second haplotype."""
    recs = []
    for aln in bam.fetch(chrom, s, e):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped or aln.query_sequence is None:
            continue
        recs.append((aln.reference_start, aln.reference_end, alts_in_span(aln, s, e), aln.query_sequence))
    from collections import defaultdict
    altc = defaultdict(int)
    for st, en, alts, _ in recs:
        for p in alts: altc[p] += 1
    cand = []
    for pos, ac in altc.items():
        if ac < 2: continue
        cov = sum(1 for st, en, _, _ in recs if st <= pos < en)
        if cov >= MIN_DEPTH and BAL_LO <= ac / cov <= BAL_HI: cand.append(pos)
    cset = set(cand)
    alt_seqs = []
    for st, en, alts, seq in recs:
        covered = sum(1 for p in cand if st <= p < en)
        if covered == 0: continue
        if sum(1 for p in alts if p in cset) / covered >= SHARE_HI:
            alt_seqs.append(seq)
    return alt_seqs, len(recs), (len(alt_seqs) / len(recs) if recs else 0.0)


def longest_orf_aa(seq):
    """Longest ORF (aa) over both strands, 3 frames each. Simple ATG..stop scan."""
    from itertools import product
    codon = {}
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    for i, (a, b, c) in enumerate(product(bases, repeat=3)):
        codon[a + b + c] = aas[i]
    best = 0
    for strand in (seq, seq.translate(COMP)[::-1]):
        for f in range(3):
            prot = "".join(codon.get(strand[i:i+3], "X") for i in range(f, len(strand) - 2, 3))
            for piece in prot.split("*"):
                # require a start
                k = piece.find("M")
                if k >= 0: best = max(best, len(piece) - k)
    return best


def mm2_best_identity(consensus_fa, ref_fa):
    """minimap2 consensus vs reference copy models; return (best_ref, identity, coverage) of the
    top hit (de-gap-compressed identity via NM/aln-len)."""
    try:
        out = subprocess.run([MM2, "-cx", "asm20", "--eqx", ref_fa, consensus_fa],
                             capture_output=True, text=True, timeout=120).stdout
    except Exception:
        return (None, None, None)
    best = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12: continue
        qlen, nmatch, alnlen = int(f[1]), int(f[9]), int(f[10])
        ident = nmatch / alnlen if alnlen else 0
        cov = alnlen / qlen if qlen else 0
        score = ident * cov
        if best is None or score > best[3]:
            best = (f[5], ident, cov, score)
    return best[:3] if best else (None, None, None)


def main():
    cand = [r for r in json.load(open(f"{CAT}/hidden_copy_catalog_validated.json")) if r["verdict"] == "COLLAPSED-COPY"]
    meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
    resc = CA.load_rescued(meta, skel, spliced)
    fams = CA.colocated_families(meta, rescued_fam=resc)
    fam_copies = {}
    for fid, dom, copies in fams:
        fam_copies.setdefault(fid, []).extend(c[0] for c in copies)
    bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
    aligner = pyabpoa.msa_aligner()
    rows = []
    allcons = open(f"{OUT}/consensus.fa", "w")
    print(f"{'family':9s} {'chrom':12s} {'start':>9} {'altrd':>5} {'len':>5} {'ORFaa':>5} "
          f"{'nearestRefCopy':18s} {'ident':>6} {'cov':>5} {'verdict':16s}")
    for r in cand:
        chrom, s, e = r["chrom"], r["start"], r["end"]
        alt_seqs, nprim, altfrac = alt_haplotype_reads(bam, chrom, s, e)
        if len(alt_seqs) < 3:
            continue
        # POA consensus (cap reads for speed; POA over up to 60)
        res = aligner.msa(alt_seqs[:60], out_cons=True, out_msa=False)
        cons = res.cons_seq[0] if res.cons_seq else max(alt_seqs, key=len)
        cid = f"{r['family']}_{chrom}_{s}"
        allcons.write(f">{cid} altreads={len(alt_seqs)} altfrac={altfrac:.2f}\n{cons}\n")
        orf = longest_orf_aa(cons)
        # family reference-copy models for the distinct-paralog test
        ref_fa = f"{OUT}/{cid}.refcopies.fa"
        with open(ref_fa, "w") as fh:
            for tid in fam_copies.get(r["family"], []):
                if tid in spliced: fh.write(f">{tid}\n{spliced[tid]}\n")
        cons_fa = f"{OUT}/{cid}.cons.fa"
        open(cons_fa, "w").write(f">{cid}\n{cons}\n")
        nearest, ident, cov = mm2_best_identity(cons_fa, ref_fa) if os.path.getsize(ref_fa) else (None, None, None)
        # verdict: homologous (maps, cov>=0.5) but distinct (ident<0.98) + coherent ORF.
        # NOTE: homology + ORF establish a real, divergent, expressed SECOND HAPLOTYPE — they do NOT
        # prove it is a reference-absent COPY rather than a hyperdivergent ALLELE (the detector cannot
        # separate the two from RNA alone; needs DNA parCN / depth-doubling). So the strong verdict is a
        # CANDIDATE, not a confirmation.
        if ident is None:
            verdict = "no-homology?"
        elif cov is not None and cov >= 0.5 and ident < 0.98 and orf >= 80:
            verdict = "second-haplotype-candidate"
        elif ident >= 0.98:
            verdict = "~identical(het?)"
        elif orf < 80:
            verdict = "weak-ORF"
        else:
            verdict = "partial-homology"
        rows.append(dict(family=r["family"], chrom=chrom, start=s, alt_reads=len(alt_seqs),
                         alt_frac=round(altfrac, 2), cons_len=len(cons), orf_aa=orf,
                         nearest_ref=nearest, identity=ident, coverage=cov, verdict=verdict))
        print(f"{r['family']:9s} {chrom:12s} {s:>9} {len(alt_seqs):>5} {len(cons):>5} {orf:>5} "
              f"{str(nearest)[:18]:18s} {(ident or 0):>6.3f} {(cov or 0):>5.2f} {verdict:16s}")
    allcons.close()
    json.dump(rows, open(f"{OUT}/promoted.json", "w"), indent=1)
    from collections import Counter
    print("\nverdict tally:", dict(Counter(x["verdict"] for x in rows)))
    print(f"wrote {OUT}/promoted.json + consensus.fa")


if __name__ == "__main__":
    main()
