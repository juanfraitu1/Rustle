#!/usr/bin/env python3
"""REFERENCE-ABSENT COPY CATALOG (collapsed/CNV population) via the hidden-copy detector run
directly on the known gene-family copy loci. Faithful Python mirror of
src/rustle/vg_family/hidden_copy.rs::detect_hidden_copy.

Signal: a reference copy that is actually 2+ COLLAPSED copies (or has an extra CNV paralog the
reference lacks) shows, among its PRIMARY alignments (the paralog-bleed firewall), a coherent
SECOND HAPLOTYPE -- >=12 balanced-frequency (0.20-0.60) alt columns that co-segregate on >=5 reads.
A diploid het is 1-2 columns; sequencing error never reaches the balanced band. DETECT + FLAG only;
never place or fabricate the copy (DAZ3 discipline)."""
import sys, json
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam

# hidden_copy.rs defaults
BAL_LO, BAL_HI = 0.20, 0.60
MIN_DEPTH = 8
MIN_ALT_POS = 12
MIN_ALT_READS = 5
SHARE_HI = 0.60

GGO_BAM = CA.GGO_BAM


def alts_in_span(aln, s, e):
    """Genomic mismatch positions (eqx 'X' = op 8) within [s,e) for one alignment."""
    out = []
    rpos = aln.reference_start
    for op, ln in (aln.cigartuples or []):
        if op in (7, 0):       # = or M (match run)
            rpos += ln
        elif op == 8:          # X (mismatch run)
            for i in range(ln):
                p = rpos + i
                if s <= p < e:
                    out.append(p)
            rpos += ln
        elif op in (2, 3):     # D, N (ref advance, no query)
            rpos += ln
        # I(1), S(4), H(5): no ref advance
    return out


def detect_hidden(reads):
    """reads = list of (start, end, alts[]). Returns dict mirroring HiddenCopyEvidence."""
    n = len(reads)
    if n < MIN_DEPTH:
        return dict(n_primary_reads=n, n_alt_positions=0, n_alt_reads=0, alt_read_fraction=0.0, flagged=False)
    from collections import defaultdict
    alt_count = defaultdict(int)
    for _, _, alts in reads:
        for p in alts:
            alt_count[p] += 1
    candidates = []
    for pos, ac in alt_count.items():
        if ac < 2:
            continue
        cov = sum(1 for st, en, _ in reads if st <= pos < en)
        if cov < MIN_DEPTH:
            continue
        frac = ac / cov
        if BAL_LO <= frac <= BAL_HI:
            candidates.append(pos)
    cand = set(candidates)
    n_alt_positions = len(candidates)
    n_alt_reads = 0
    for st, en, alts in reads:
        covered = sum(1 for pos in candidates if st <= pos < en)
        if covered == 0:
            continue
        alt_at = sum(1 for p in alts if p in cand)
        if alt_at / covered >= SHARE_HI:
            n_alt_reads += 1
    flagged = n_alt_positions >= MIN_ALT_POS and n_alt_reads >= MIN_ALT_READS
    return dict(n_primary_reads=n, n_alt_positions=n_alt_positions, n_alt_reads=n_alt_reads,
                alt_read_fraction=(n_alt_reads / n if n else 0.0), flagged=flagged)


def main():
    meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
    resc = CA.load_rescued(meta, skel, spliced)
    fams = CA.colocated_families(meta, rescued_fam=resc)
    bam = pysam.AlignmentFile(GGO_BAM, "rb")
    catalog = []
    n_copy_loci = 0
    for fid, dom, copies in fams:
        for (tid, chrom, s1, e, strand) in copies:
            s = s1 - 1
            n_copy_loci += 1
            reads = []
            for aln in bam.fetch(chrom, s, e):
                if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                    continue  # primary-only firewall
                reads.append((aln.reference_start, aln.reference_end, alts_in_span(aln, s, e)))
            ev = detect_hidden(reads)
            if ev["flagged"]:
                catalog.append(dict(family=fid, dom=dom, chrom=chrom, start=s, end=e, strand=strand,
                                    copy_tid=tid, **ev))
    catalog.sort(key=lambda r: -r["n_alt_positions"])
    print(f"Scanned {n_copy_loci} reference-copy loci across {len(fams)} families.")
    print(f"FLAGGED (reads imply an unmodeled copy at the locus): {len(catalog)}\n")
    print(f"{'family':10s} {'dom':10s} {'chrom':13s} {'start':>10} {'altcols':>7} {'altreads':>8} {'altfrac':>7} {'nprim':>6}")
    for r in catalog:
        print(f"{r['family']:10s} {r['dom']:10s} {r['chrom']:13s} {r['start']:>10} "
              f"{r['n_alt_positions']:>7} {r['n_alt_reads']:>8} {r['alt_read_fraction']:>7.2f} {r['n_primary_reads']:>6}")
    json.dump(catalog, open("/home/juanfra/winloci_scratch/refabsent/hidden_copy_catalog.json", "w"), indent=1)
    with open("/home/juanfra/winloci_scratch/refabsent/hidden_copy_catalog.tsv", "w") as fh:
        fh.write("family\tdom\tchrom\tstart\tend\tstrand\tcopy_tid\tn_primary\talt_cols\talt_reads\talt_frac\n")
        for r in catalog:
            fh.write(f"{r['family']}\t{r['dom']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['strand']}\t"
                     f"{r['copy_tid']}\t{r['n_primary_reads']}\t{r['n_alt_positions']}\t{r['n_alt_reads']}\t{r['alt_read_fraction']:.3f}\n")
    print(f"\nwrote hidden_copy_catalog.{{json,tsv}} ({len(catalog)} flagged loci)")


if __name__ == "__main__":
    main()
