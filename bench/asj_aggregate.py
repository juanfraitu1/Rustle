#!/usr/bin/env python3
"""Aggregate the genome-wide allele-specific-junction scan: BH-FDR over all tested junctions, call
ASJ (q<0.05 AND |dPSI|>=DPSI), characterize, and flag the editing confound (a transition anchor
A/G or C/T could be an RNA-edit site phasing by edit-status, not allele; a TRANSVERSION anchor is
unambiguously genetic). Writes bench/asj_results.md + asj_calls.tsv.
"""
import glob
import os
from collections import defaultdict

GW = "/tmp/gw"
DPSI = 0.30
Q = 0.05
OUT_MD = os.path.join(os.path.dirname(__file__), "asj_results.md")
OUT_TSV = os.path.join(os.path.dirname(__file__), "asj_calls.tsv")
TRANSITION = ({"A", "G"}, {"C", "T"})


def main():
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
        fh.write("gene\tchrom\tanchor\talleleX\talleleY\tdonor\tacceptor\tpsiX\tpsiY\tdPSI\tpval\tq\tanchor_type\n")
        for r in sorted(asj, key=lambda r: (-r["dPSI"], r["q"])):
            at = "transition" if {r["alleleX"], r["alleleY"]} in TRANSITION else "transversion"
            fh.write(f"{r['gene']}\t{r['chrom']}\t{r['anchor']}\t{r['alleleX']}\t{r['alleleY']}\t"
                     f"{r['donor']}\t{r['acceptor']}\t{r['psiX']}\t{r['psiY']}\t{r['dPSI']:.3f}\t"
                     f"{r['pval']:.3e}\t{r['q']:.3e}\t{at}\n")
    print("\n".join(L))
    print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")


if __name__ == "__main__":
    main()
