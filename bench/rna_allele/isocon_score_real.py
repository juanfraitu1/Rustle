#!/usr/bin/env python3
"""Score IsoCon on Arm I (real fibroblast reads; Amendments 1-3) against KB3781's own haplotypes.

Each S_fam copy has a locus on both haplotypes: on the haplotype `_pri` took its chromosome from, the `_pri` coordinates (identical
sequence); on the other (B), the lifted interval from lift.tsv (T2d/T2i copies only). B-only loci come from truth_fam.tsv. Every IsoCon
candidate is mapped to the maternal and paternal assemblies (`minimap2 -c -x splice -N 50`); per haplotype its best hit (identity =
matches / block length; ties broken by aligned length) is assigned to the copy locus it overlaps on that haplotype.
  haplotype of a candidate: the haplotype whose best on-copy hit has the higher identity; equal (within 0.0005) = both
  recovered: haplotype copies with a candidate at identity >= 0.999 on that haplotype's locus
  alleles separated / merged: expressed T2d copies with recovered candidates on both haplotypes / on one only
  paralogs merged: candidates whose best hit on one haplotype ties (within 0.0005) between two different copies
  other genes: candidates with no on-copy hit at >= 0.99 but a hit anywhere at identity x coverage >= 0.99 (the broad read net also
               catches reads of genes that share sequence with the family, e.g. PDXDC1 at NPIP)
  unmatched: candidates with no hit anywhere at identity x coverage >= 0.99
Expressed copy = >= 2 primary reads overlapping its exons (from the read sets used as IsoCon input, by their primary placement).

    isocon_score_real.py --work /mnt/linuxdisk/tmp/rna_allele --paf-mat real_cands.mat.paf --paf-pat real_cands.pat.paf
"""
import argparse
import collections
import csv
import re

import pysam


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("work", "paf_mat", "paf_pat"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    ap.add_argument("--bam", default="/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam")
    a = ap.parse_args(argv)
    W = a.work
    cm = {r["pri"]: r for r in csv.DictReader(open(f"{W}/chrmap.tsv"), delimiter="\t")}
    genes = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/genes.tsv"), delimiter="\t") if r["set"] == "S_fam"}
    truth = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t")}
    lift = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/lift.tsv"), delimiter="\t")}
    alias = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{W}/{h}.len.tsv"), delimiter="\t"):
            alias[(h, num)] = acc
    # copy loci per haplotype: (hap, accession, start, end, family, copy)
    loci = []
    for g, r in genes.items():
        fam = r["biotype"].split(":", 1)[1]
        c = cm[r["chrom"]]
        ex = [tuple(int(x) for x in b.split("-")) for b in r["exons"].split(",")]
        loci.append((c["same_hap"], c["same_name"], ex[0][0], ex[-1][1], fam, g))
        L = lift.get(g)
        if truth[g]["class"] in ("T2d", "T2i") and L and L["B_chrom"]:
            loci.append((c["B_hap"], L["B_chrom"], int(L["B_start"]), int(L["B_end"]), fam, g))
    acc_hap = {acc: h for (h, _), acc in alias.items()}
    for fr in csv.DictReader(open(f"{W}/truth_fam.tsv"), delimiter="\t"):
        for k, loc in enumerate(x for x in fr["B_only_detail"].split(";") if x):
            acc, se = loc.rsplit(":", 1)
            s0, e0 = (int(v) for v in se.split("-"))
            loci.append((acc_hap[acc], acc, s0, e0, fr["family"], f"{fr['family']}_Bonly{k}"))

    def hits(path, hap):
        out = collections.defaultdict(list)
        for ln in open(path):
            f = ln.rstrip("\n").split("\t")
            m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]+", f[5])
            if not m:
                continue
            acc = alias[(m.group(2), m.group(1))]
            ident = int(f[9]) / max(1, int(f[10]))
            fam = f[0].split("|")[0]
            on = [l for l in loci if l[0] == hap and l[1] == acc and l[4] == fam and int(f[7]) < l[3] and l[2] < int(f[8])]
            out[f[0]].append((ident, int(f[10]), [l[5] for l in on]))
            allbest[f[0]] = max(allbest.get(f[0], 0.0), ident * (int(f[3]) - int(f[2])) / max(1, int(f[1])))
        return out
    allbest = {}
    H = {"mat": hits(a.paf_mat, "mat"), "pat": hits(a.paf_pat, "pat")}
    cands = sorted(set(H["mat"]) | set(H["pat"]))
    # expressed copies: >= 2 primaries over exons, among the IsoCon input reads
    bam = pysam.AlignmentFile(a.bam)
    expressed = {}
    for g, r in genes.items():
        ex = [tuple(int(x) for x in b.split("-")) for b in r["exons"].split(",")]
        n = 0
        for rd in bam.fetch(r["chrom"], ex[0][0], ex[-1][1]):
            if rd.is_secondary or rd.is_supplementary or rd.is_unmapped:
                continue
            if any(rd.reference_start < y and x < rd.reference_end for x, y in ex):
                n += 1
        expressed[g] = n
    for fam in ("NPIP", "TBC1D3"):
        fc = [c for c in cands if c.startswith(fam + "|")]
        rec = collections.defaultdict(set)      # copy -> haplotypes recovered
        para_merged, unmatched, other_gene = 0, 0, 0
        for c in fc:
            best = {}
            for h in ("mat", "pat"):
                on = [(i, L, cp) for i, L, cp in H[h].get(c, []) if cp]
                if on:
                    top = max(on, key=lambda x: (x[0], x[1]))
                    tied = {cpy for i, L, cp in on if abs(i - top[0]) <= 0.0005 for cpy in cp}
                    best[h] = (top[0], tied)
            if not best or max(v[0] for v in best.values()) < 0.99:
                if allbest.get(c, 0.0) >= 0.99:
                    other_gene += 1          # a faithful transcript of a locus outside the family's copies (the read net is broad)
                else:
                    unmatched += 1
                continue
            if any(len(v[1]) > 1 for v in best.values()):
                para_merged += 1
            m = max(v[0] for v in best.values())
            haps = [h for h, v in best.items() if m - v[0] <= 0.0005]
            if m >= 0.999:
                for h in haps:
                    for cp in best[h][1]:
                        rec[cp].add(h)
        cps = [g for g, r in genes.items() if r["biotype"] == f"S_fam:{fam}"]
        expr = [g for g in cps if expressed[g] >= 2]
        t2d = [g for g in expr if truth[g]["class"] == "T2d"]
        sep = [g for g in t2d if len(rec.get(g, ())) == 2]
        mer = [g for g in t2d if len(rec.get(g, ())) == 1]
        hapcopies = sum(len(v) for k, v in rec.items() if not k.endswith(tuple(f"_Bonly{i}" for i in range(9))))
        bonly = [k for k in rec if "_Bonly" in k]
        n_bonly_truth = sum(1 for l in loci if l[4] == fam and "_Bonly" in l[5])
        T_expr = sum(2 if truth[g]["class"] in ("T2d", "T2i") else 1 for g in expr)
        print(f"\n== {fam}: {len(fc)} IsoCon candidates; copies {len(cps)}, expressed (>= 2 primaries) {len(expr)}")
        print(f"   haplotype copies recovered (>= 0.999): {hapcopies} of {T_expr} for the expressed copies "
              f"(T for the family: {sum(2 if truth[g]['class'] in ('T2d','T2i') else 1 for g in cps) + n_bonly_truth})")
        print(f"   expressed T2d copies {len(t2d)}: alleles separated {len(sep)}, merged {len(mer)}, not recovered {len(t2d) - len(sep) - len(mer)}")
        print(f"   candidates on family copies {len(fc) - other_gene - unmatched}; tied between copies (paralogs merged): {para_merged}")
        print(f"   off-copy candidates matching another haplotype locus at >= 0.99 x coverage (other genes in the read net): "
              f"{other_gene}; matching nothing at >= 0.99: {unmatched}")
        print(f"   B-only (reference-absent) copies recovered: {len(bonly)} of {n_bonly_truth}")
        print(f"   per-copy primary reads: " + ", ".join(f"{genes[g]['name']}:{expressed[g]}" for g in sorted(cps, key=lambda g: -expressed[g])[:12]))


if __name__ == "__main__":
    main()
