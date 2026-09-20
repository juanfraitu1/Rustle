#!/usr/bin/env python3
"""Aggregate the MULTI-SNP allele-specific-junction scan: BH-FDR, call ASJ, and COMPARE to the
single-anchor result (the reach/power gain from multi-SNP haplotype phasing). phase_qual is the
multi-SNP-native confound check (clean 2-way bipartition = diploid het; low = paralog/>2 haplotypes).
Writes bench/asjm_findings.md + asjm_calls.tsv.
"""
import csv
import glob
import os
from collections import defaultdict

GW = "/tmp/gw"
DPSI, Q = 0.30, 0.05
PHASE_QUAL_MIN = 0.70
OUT_MD = os.path.join(os.path.dirname(__file__), "asjm_findings.md")
OUT_TSV = os.path.join(os.path.dirname(__file__), "asjm_calls.tsv")
SINGLE = os.path.join(os.path.dirname(__file__), "asj_calls.tsv")


def bh(rows):
    rows.sort(key=lambda r: r["pval"])
    m = len(rows); prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, rows[i]["pval"] * m / (i + 1))
        rows[i]["q"] = prev
    return rows


def main():
    rows = []
    for p in sorted(glob.glob(os.path.join(GW, "asjm_NC_*.tsv"))):
        with open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) != len(hdr):
                    continue
                d = dict(zip(hdr, c))
                d["pval"] = float(d["pval"]); d["dPSI"] = float(d["dPSI"])
                d["n_snp"] = int(d["n_snp"]); d["n_tv"] = int(d["n_tv"]); d["phase_qual"] = float(d["phase_qual"])
                rows.append(d)
    if not rows:
        print("no asjm rows yet"); return
    bh(rows)
    asj = [r for r in rows if r["q"] < Q and r["dPSI"] >= DPSI]
    hc = [r for r in asj if r["n_tv"] >= 1 and r["phase_qual"] >= PHASE_QUAL_MIN]
    multi = [r for r in asj if r["n_snp"] >= 2]              # genes where multi-SNP phasing does real work
    multi_genes = {r["gene"] for r in multi}

    # compare to single-anchor
    single = set()
    if os.path.exists(SINGLE):
        with open(SINGLE) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                single.add((r["gene"], r["donor"], r["acceptor"]))
    mset = {(r["gene"], r["donor"], r["acceptor"]) for r in asj}
    both = mset & single
    multi_only = mset - single
    single_only = single - mset

    def G(rs): return len({r["gene"] for r in rs})
    L = []
    def P(s=""): L.append(s)
    P("# Allele-specific junctions — MULTI-SNP haplotype phasing")
    P()
    P("Each read is assigned to one of the two diploid haplotypes using ALL het SNPs in the gene "
      "(read-backed 2-means), not a single anchor SNP — so reads covering any subset of SNPs get phased "
      "(more reach/power), and the phasing quality (how cleanly reads bipartition) is a built-in "
      "confound check. With 1 SNP it degenerates to single-anchor, so this is a SUPERSET.")
    P()
    P("## Result")
    P(f"- junctions tested: **{len(rows)}**")
    P(f"- **allele-specific junctions (q<{Q}, |dPSI|>={DPSI}): {len(asj)}** across **{G(asj)}** genes")
    P(f"- high-confidence (>=1 transversion SNP [genetic] AND phase_qual>={PHASE_QUAL_MIN}): "
      f"**{len(hc)}** / {G(hc)} genes")
    P(f"- from genes with >=2 het SNPs (multi-SNP phasing does real work): **{len(multi)}** / "
      f"{len(multi_genes)} genes")
    P()
    P("## Gain over single-anchor phasing")
    P(f"- single-anchor ASJ (committed): **{len(single)}**; multi-SNP ASJ: **{len(asj)}**")
    P(f"- shared (gene+junction): **{len(both)}**; **multi-SNP-only: {len(multi_only)}** (the reach/power gain); "
      f"single-anchor-only: {len(single_only)}")
    P(f"- net: multi-SNP {'ADDS' if len(asj) > len(single) else 'changes'} "
      f"{len(asj) - len(single):+d} ASJ ({G(asj) - 0} genes)")
    P()
    P("## n_snp distribution of ASJ")
    nd = defaultdict(int)
    for r in asj:
        nd[min(r["n_snp"], 10)] += 1
    P("```")
    for k in sorted(nd):
        P(f"  n_snp={k}{'+' if k==10 else ' '}: {nd[k]}")
    P("```")
    P()
    P("## Top multi-SNP-only allele-specific junctions (reach/power gain)")
    P("| gene | chrom | n_snp | phase_qual | junction | psi0 | psi1 | dPSI | q |")
    P("|---|---|---|---|---|---|---|---|---|")
    mo_rows = [r for r in asj if (r["gene"], r["donor"], r["acceptor"]) in multi_only and r["n_snp"] >= 2]
    for r in sorted(mo_rows, key=lambda r: (-r["dPSI"], r["q"]))[:20]:
        P(f"| {r['gene']} | {r['chrom']} | {r['n_snp']} | {r['phase_qual']} | "
          f"{r['donor']}-{r['acceptor']} | {r['psi0']} | {r['psi1']} | {r['dPSI']:.2f} | {r['q']:.1e} |")
    P()
    P("## Caveats")
    P("- phase_qual is the diploid-bipartition cleanliness; low phase_qual = paralog/>2 haplotypes "
      "(down-weighted by the phase_qual>=0.7 high-confidence filter).")
    P("- transition-only-SNP loci (n_tv=0) could be edit-phased; the genetic core requires >=1 "
      "transversion het SNP.")
    P("- still coverage-limited (>=12x, >=5 reads per haplotype per junction).")
    P()
    P("## Reproduce")
    P("- `MINIFORGE python bench/allele_specific_junctions_multisnp.py --chrom <C>`")
    P("- `python3 bench/asjm_aggregate.py`")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(OUT_TSV, "w") as fh:
        fh.write("gene\tchrom\tn_snp\tn_tv\tphase_qual\tdonor\tacceptor\tpsi0\tpsi1\tdPSI\tq\tmulti_snp_only\n")
        for r in sorted(asj, key=lambda r: (-r["dPSI"], r["q"])):
            mo = int((r["gene"], r["donor"], r["acceptor"]) in multi_only)
            fh.write(f"{r['gene']}\t{r['chrom']}\t{r['n_snp']}\t{r['n_tv']}\t{r['phase_qual']}\t"
                     f"{r['donor']}\t{r['acceptor']}\t{r['psi0']}\t{r['psi1']}\t{r['dPSI']:.3f}\t"
                     f"{r['q']:.3e}\t{mo}\n")
    print("\n".join(L))
    print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")


if __name__ == "__main__":
    main()
