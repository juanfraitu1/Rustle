#!/usr/bin/env python3
"""Aggregate the genome-wide allele-specific-junction scan: BH-FDR over all tested junctions, call
ASJ (q<0.05 AND |dPSI|>=DPSI), characterize, and flag the editing confound (a transition anchor
A/G or C/T could be an RNA-edit site phasing by edit-status, not allele; a TRANSVERSION anchor is
unambiguously genetic). Writes bench/asj_results.md + asj_calls.tsv.

Strand-bias gate: every call is annotated with a GATK StrandOddsRatio (SOR) computed from the BAM
(see bench/asj_strand_bias.py for the formula + longcallR rationale, Nat Methods 2026: reject
allele-specific junctions whose usage is strand-confounded). SOR/sor_pass are ANNOTATED on every
call so existing outputs are not silently changed; calls are HARD-filtered only when --max-sor is
passed (default behaviour leaves the call set identical to before this gate existed). Annotation
needs --bam + pysam (miniforge python); without it SOR is left blank and nothing is dropped.
"""
import argparse
import glob
import math
import os
from collections import defaultdict

GW = "/tmp/gw"
DPSI = 0.30
Q = 0.05
OUT_MD = os.path.join(os.path.dirname(__file__), "asj_results.md")
OUT_TSV = os.path.join(os.path.dirname(__file__), "asj_calls.tsv")
TRANSITION = ({"A", "G"}, {"C", "T"})
DEFAULT_BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"


def _introns(aln):
    """CIGAR-derived intron chain (ref coords of each N block); mirrors allele_specific_junctions.py."""
    out = []
    rpos = aln.reference_start
    for op, length in aln.cigartuples:
        if op == 3:
            out.append((rpos, rpos + length)); rpos += length
        elif op in (0, 2, 7, 8):
            rpos += length
    return out


def _sor_for_call(bam, refs, chrom, donor, acceptor):
    """GATK StrandOddsRatio for one junction (used/notused x fwd/rev, +1 pseudocount).
    Reuses the logic of bench/asj_strand_bias.py — see that module for the formula + citation."""
    if donor > acceptor:
        donor, acceptor = acceptor, donor
    if chrom not in refs:
        return float("nan")
    rf = rr = af = ar = 0
    j = (donor, acceptor)
    for aln in bam.fetch(chrom, donor, acceptor + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.cigartuples is None:
            continue
        if aln.reference_start > donor or aln.reference_end < acceptor:
            continue
        used = j in _introns(aln)
        if used:
            rr += aln.is_reverse; rf += not aln.is_reverse
        else:
            ar += aln.is_reverse; af += not aln.is_reverse
    rp, rm, ap, am = rf + 1.0, rr + 1.0, af + 1.0, ar + 1.0
    ratio = (rp * am) / (ap * rm)
    sym = ratio + 1.0 / ratio
    return math.log(sym) + math.log(min(rp, rm) / max(rp, rm)) - math.log(min(ap, am) / max(ap, am))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default=DEFAULT_BAM,
                    help="BAM for SOR strand-bias annotation (needs pysam); '' to skip annotation")
    ap.add_argument("--max-sor", type=float, default=None,
                    help="if set, HARD-DROP calls with SOR > this (GATK SNV cutoff ~3.0, "
                         "longcallR-strict ~2.0); default None = annotate only, drop nothing")
    args = ap.parse_args()

    rows = []
    for p in sorted(glob.glob(os.path.join(GW, "asj_*.tsv"))):
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) != len(hdr):
                    continue
                d = dict(zip(hdr, c))
                d["dPSI"] = float(d["dPSI"]); d["pval"] = float(d["pval"])
                rows.append(d)
    m = len(rows)
    if m == 0:
        print("no asj rows yet"); return
    # BH-FDR
    rows.sort(key=lambda r: r["pval"])
    q = [0.0] * m
    prev = 1.0
    for i in range(m - 1, -1, -1):
        val = rows[i]["pval"] * m / (i + 1)
        prev = min(prev, val)
        q[i] = prev
    for i in range(m):
        rows[i]["q"] = q[i]

    asj = [r for r in rows if r["q"] < Q and r["dPSI"] >= DPSI]

    # Strand-bias (SOR) annotation + optional hard-filter (see module docstring / asj_strand_bias.py).
    n_sor_drop = 0
    if args.bam:
        try:
            import pysam
            bam = pysam.AlignmentFile(args.bam, "rb")
            refs = set(bam.references)
            for r in asj:
                s = _sor_for_call(bam, refs, r["chrom"], int(r["donor"]), int(r["acceptor"]))
                r["sor"] = s
                r["sor_pass"] = int(not math.isnan(s) and s < 3.0)
            if args.max_sor is not None:
                before = len(asj)
                asj = [r for r in asj if math.isnan(r["sor"]) or r["sor"] <= args.max_sor]
                n_sor_drop = before - len(asj)
        except ImportError:
            print("[warn] pysam unavailable; skipping SOR annotation (use miniforge python)")
            for r in asj:
                r["sor"] = float("nan"); r["sor_pass"] = ""
    else:
        for r in asj:
            r["sor"] = float("nan"); r["sor_pass"] = ""

    genes_tested = {r["gene"] for r in rows}
    genes_asj = {r["gene"] for r in asj}
    transversion_asj = [r for r in asj
                        if {r["alleleX"], r["alleleY"]} not in TRANSITION]

    L = []
    def P(s=""): L.append(s)
    P("# Allele-specific junctions (genome-wide)")
    P()
    P("Splice junctions whose usage depends on the allele a molecule carries. Long HiFi reads link "
      "allele->junction PER MOLECULE (the read carries both its het-SNP allele and its junctions), so "
      "no statistical phasing is needed. The het loci that confounded copy-detection (task c) are the "
      "substrate: heterozygosity = the phasing signal (confound -> feature).")
    P()
    P("## Result")
    P(f"- genes with a phaseable het anchor + a tested junction: **{len(genes_tested)}**")
    P(f"- alternatively-spliced junctions tested (non-constitutive): **{m}**")
    P(f"- **allele-specific junctions (BH-FDR q<{Q} AND |dPSI|>={DPSI}): {len(asj)}** "
      f"across **{len(genes_asj)}** genes")
    P(f"- of those, **{len(transversion_asj)}** have a TRANSVERSION anchor (unambiguously genetic, "
      f"not RNA-editing) — the high-confidence genetic set")
    P()
    P("## |dPSI| distribution of the ASJ calls")
    if asj:
        import numpy as np
        dd = np.array([r["dPSI"] for r in asj])
        P(f"- median |dPSI|={np.median(dd):.2f}, max={dd.max():.2f}; "
          f">=0.5: {(dd>=0.5).sum()}, >=0.7: {(dd>=0.7).sum()}, ==1.0 (full switch): {(dd>=0.999).sum()}")
    P()
    P("## Top allele-specific junctions")
    P("| gene | chrom | anchor | alleles | junction | PSI_X | PSI_Y | dPSI | q | anchor type |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(asj, key=lambda r: (-r["dPSI"], r["q"]))[:25]:
        at = "transition(poss.edit)" if {r["alleleX"], r["alleleY"]} in TRANSITION else "transversion(genetic)"
        P(f"| {r['gene']} | {r['chrom']} | {r['anchor']} | {r['alleleX']}/{r['alleleY']} | "
          f"{r['donor']}-{r['acceptor']} | {r['psiX']} | {r['psiY']} | {r['dPSI']:.2f} | "
          f"{r['q']:.1e} | {at} |")
    P()
    P("## Honest caveats")
    P("- **Editing confound:** a TRANSITION anchor (A/G or C/T) could be an RNA-edit site (phasing by "
      "edit-status, not allele) — that is edit-coupled splicing, still real but not *genetic* allele-"
      "specific. TRANSVERSION anchors are unambiguously genetic; treat the transversion subset as the "
      "high-confidence genetic ASJ. (Strand-aware editing detection would reclassify the transitions.)")
    P("- **Single-anchor phasing:** reads are phased by ONE balanced het SNP; a junction is tested only "
      "among reads spanning BOTH the anchor and the junction (long reads help). Multi-SNP haplotype "
      "phasing would add power/reach.")
    P("- **Collapsed-paralog masquerade:** at a collapsed locus the two 'alleles' could be paralog "
      "copies (=copy-specific splicing, also interesting, but not within-gene allele-specific). The "
      "het substrate from task c was uniquely-mapped (frac_mq0=0); genome-wide, low-MAPQ loci should "
      "be down-weighted.")
    P("- **Coverage-limited:** needs >=12x at the anchor and >=5 spanning reads per allele; deeper data "
      "finds more.")
    P()
    P("## Reproduce")
    P("- `MINIFORGE python bench/allele_specific_junctions.py --chrom <C>` (per-chrom) ; `--region` for one locus")
    P("- `python3 bench/asj_aggregate.py` (this: BH-FDR + calls)")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(OUT_TSV, "w") as fh:
        fh.write("gene\tchrom\tanchor\talleleX\talleleY\tdonor\tacceptor\tpsiX\tpsiY\tdPSI\tpval\tq\t"
                 "anchor_type\tsor\tsor_pass\n")
        for r in sorted(asj, key=lambda r: (-r["dPSI"], r["q"])):
            at = "transition" if {r["alleleX"], r["alleleY"]} in TRANSITION else "transversion"
            sor = r.get("sor", float("nan"))
            sor_s = "" if (isinstance(sor, float) and math.isnan(sor)) else f"{sor:.3f}"
            fh.write(f"{r['gene']}\t{r['chrom']}\t{r['anchor']}\t{r['alleleX']}\t{r['alleleY']}\t"
                     f"{r['donor']}\t{r['acceptor']}\t{r['psiX']}\t{r['psiY']}\t{r['dPSI']:.3f}\t"
                     f"{r['pval']:.3e}\t{r['q']:.3e}\t{at}\t{sor_s}\t{r.get('sor_pass', '')}\n")
    print("\n".join(L))
    if args.bam:
        n_pass = sum(1 for r in asj if r.get("sor_pass") == 1)
        print(f"[SOR strand-bias] annotated {len(asj)} calls; {n_pass} pass SOR<3.0"
              + (f"; --max-sor={args.max_sor} dropped {n_sor_drop}" if args.max_sor is not None else
                 " (annotate-only; nothing dropped)"))
    print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")


if __name__ == "__main__":
    main()
