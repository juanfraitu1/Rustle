#!/usr/bin/env python3
"""Validate the hidden-copy (reference-absent) flags: discriminate a genuine COLLAPSED COPY (diverse
ref->alt substitution spectrum = real sequence divergence) from the main RNA confound, RNA EDITING
(A>G / T>C dominated), plus tier by alt-haplotype fraction. Adds the substitution spectrum each
flagged locus, computed at its balanced candidate columns from the reference base + the alt reads' base."""
import sys, json
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam

BAL_LO, BAL_HI, MIN_DEPTH = 0.20, 0.60, 8
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def alt_calls_in_span(aln, s, e):
    """list of (genomic_pos, query_base) at eqx 'X' mismatches within [s,e)."""
    out = []
    rpos = aln.reference_start
    qpos = 0
    q = aln.query_sequence or ""
    for op, ln in (aln.cigartuples or []):
        if op in (7, 0):
            rpos += ln; qpos += ln
        elif op == 8:
            for i in range(ln):
                p = rpos + i
                if s <= p < e and qpos + i < len(q):
                    out.append((p, q[qpos + i]))
            rpos += ln; qpos += ln
        elif op in (2, 3):
            rpos += ln
        elif op in (1, 4):
            qpos += ln
        # H(5): nothing
    return out


def main():
    cat = json.load(open(f"{CAT}/hidden_copy_catalog.json"))
    fa = pysam.FastaFile(FASTA)
    bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
    print(f"{'family':10s} {'dom':8s} {'chrom':13s} {'start':>10} {'altcols':>7} {'altfrac':>7} "
          f"{'%A>G+T>C':>8} {'n_subtypes':>10} {'verdict':14s}")
    out = []
    for r in cat:
        chrom, s, e = r["chrom"], r["start"], r["end"]
        reads = []
        for aln in bam.fetch(chrom, s, e):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            reads.append((aln.reference_start, aln.reference_end, alt_calls_in_span(aln, s, e)))
        # recompute balanced candidate columns + collect alt base per column
        from collections import defaultdict, Counter
        altpos = defaultdict(list)
        for st, en, calls in reads:
            for p, b in calls:
                altpos[p].append(b)
        spectrum = Counter()
        ncand = 0
        for pos, bases in altpos.items():
            ac = len(bases)
            if ac < 2:
                continue
            cov = sum(1 for st, en, _ in reads if st <= pos < en)
            if cov < MIN_DEPTH:
                continue
            frac = ac / cov
            if not (BAL_LO <= frac <= BAL_HI):
                continue
            ncand += 1
            refb = fa.fetch(chrom, pos, pos + 1).upper()
            altb = Counter(bases).most_common(1)[0][0]
            if refb in "ACGT" and altb in "ACGT" and refb != altb:
                spectrum[(refb, altb)] += 1
        tot = sum(spectrum.values())
        ag_tc = spectrum[("A", "G")] + spectrum[("T", "C")]
        ag_tc_frac = ag_tc / tot if tot else 0.0
        n_subtypes = len(spectrum)
        # verdict: editing if A>G/T>C dominate; weak if minority haplotype; else collapsed-copy candidate
        if ag_tc_frac >= 0.50:
            verdict = "RNA-editing?"
        elif r["alt_read_fraction"] < 0.20:
            verdict = "weak/minority"
        elif n_subtypes >= 6 and ag_tc_frac < 0.40:
            verdict = "COLLAPSED-COPY"
        else:
            verdict = "ambiguous"
        out.append(dict(**r, spectrum_total=tot, ag_tc_frac=round(ag_tc_frac, 2),
                        n_subtypes=n_subtypes, verdict=verdict))
        print(f"{r['family']:10s} {r['dom']:8s} {chrom:13s} {s:>10} {r['n_alt_positions']:>7} "
              f"{r['alt_read_fraction']:>7.2f} {ag_tc_frac:>8.2f} {n_subtypes:>10} {verdict:14s}")
    from collections import Counter as C2
    vc = C2(o["verdict"] for o in out)
    print("\nverdict tally:", dict(vc))
    json.dump(out, open(f"{CAT}/hidden_copy_catalog_validated.json", "w"), indent=1)
    print(f"wrote {CAT}/hidden_copy_catalog_validated.json")


if __name__ == "__main__":
    main()
