#!/usr/bin/env python3
"""Final truth classes for docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md (truth steps 3-5).

- lift.tsv (truth_lift.py) gives T2d / T2i / T? / T1-candidate per record.
- A T1-candidate becomes T1 unless a locus anywhere on the B haplotype is strictly closer to its exon sequence than its closest `_pri`
  paralog (sets.tsv best_paralog_identity), in which case it is T? (the copy moved). B-haplotype hits come from mapping the candidates'
  exon sums to the whole B assembly (`minimap2 -c -x splice -N 50 <B>.splice.mmi`); B = the other haplotype of the record's chromosome.
- chrX records outside the PAR are T1 by construction (KB3781 is male).
- B-only copies of the S_fam families: hits of the family's copy exon sums on the B haplotype of the hit's chromosome, identity >= 0.90 and
  coverage >= 0.80, that overlap the lift (B interval) of none of the family's `_pri` copies; overlapping hits merge into one locus.
- Haplotype copies of a family T = sum over its `_pri` copies (2 if T2d/T2i, 1 if T1) + its B-only loci (T? copies listed separately).

    python3 truth_classes.py --chrmap chrmap.tsv --genes genes.tsv --lift lift.tsv --sets sets.tsv \
        --t1-paf-mat t1.mat.paf --t1-paf-pat t1.pat.paf --fam-paf-mat fam.mat.paf --fam-paf-pat fam.pat.paf --hap-len-mat mat.len.tsv \
        --hap-len-pat pat.len.tsv --out truth.tsv --fam-out truth_fam.tsv
"""
import argparse
import collections
import csv
import re

ALIAS = {}   # haplotype splice-index names (chrN_mat_hsaX) -> GenBank accession of that chromosome (filled in main, length-checked)


def paf(path):
    out = collections.defaultdict(list)
    for ln in open(path):
        f = ln.split("\t")
        t = f[5]
        if t in ALIAS:
            acc, ln_ = ALIAS[t]
            assert ln_ == int(f[6]), (t, ln_, f[6])
            t = acc
        out[f[0]].append(dict(t=t, s=int(f[7]), e=int(f[8]), id=int(f[9]) / max(1, int(f[10])), cov=(int(f[3]) - int(f[2])) / max(1, int(f[1]))))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("chrmap", "genes", "lift", "sets", "t1_paf_mat", "t1_paf_pat", "fam_paf_mat", "fam_paf_pat", "out", "fam_out",
              "hap_len_mat", "hap_len_pat"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    for h, path in (("mat", a.hap_len_mat), ("pat", a.hap_len_pat)):
        for acc, num, ln_ in csv.reader(open(path), delimiter="\t"):
            ALIAS[f"__{h}_{num}"] = (acc, int(ln_))
    # the splice indexes name chromosomes chr<N>_<hap>_hsa<...> (no "_random"): alias them to the accession, checked by length
    def alias_name(t):
        m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]+", t)
        return ALIAS.get(f"__{m.group(2)}_{m.group(1)}") if m else None
    for path in (a.t1_paf_mat, a.t1_paf_pat, a.fam_paf_mat, a.fam_paf_pat):
        for ln in open(path):
            t = ln.split("\t", 6)[5]
            if t not in ALIAS and alias_name(t):
                ALIAS[t] = alias_name(t)
    cm = list(csv.DictReader(open(a.chrmap), delimiter="\t"))
    bhap_of_pri = {r["pri"]: r["B_hap"] for r in cm}
    # which haplotype chromosome names are "B" (absent from _pri): B_name of every row
    b_names = {r["B_name"]: r["B_hap"] for r in cm if r["B_name"]}
    genes = {r["gene_id"]: r for r in csv.DictReader(open(a.genes), delimiter="\t")}
    lift = {r["gene_id"]: r for r in csv.DictReader(open(a.lift), delimiter="\t")}
    sets = {r["gene_id"]: r for r in csv.DictReader(open(a.sets), delimiter="\t")}
    t1 = {"mat": paf(a.t1_paf_mat), "pat": paf(a.t1_paf_pat)}
    out = open(a.out, "w")
    out.write("gene_id\tname\tset\tclass\tlift_class\tlift_frac\tmismatches\tindels\tbest_paralog_identity\tbest_B_elsewhere_identity\n")
    final = {}
    for g, s in sets.items():
        if s["set"] == "excluded":
            continue
        L = lift.get(g)
        lc = L["class"] if L else "NA"
        bp = float(s["best_paralog_identity"])
        bb = ""
        if s["chrom_label"] == "X":
            cls = "T1"
        elif lc == "T1-candidate":
            h = bhap_of_pri.get(genes[g]["chrom"], "")
            hs = [x for x in t1.get(h, {}).get(g, []) if x["t"] in b_names and b_names[x["t"]] == h] if h else []
            best = max((x["id"] for x in hs if x["cov"] >= 0.50), default=0.0)
            bb = f"{best:.4f}"
            cls = "T?" if best > bp else "T1"
        else:
            cls = lc
        final[g] = cls
        out.write(f"{g}\t{s['name']}\t{s['set']}\t{cls}\t{lc}\t{L['lift_frac'] if L else ''}\t{L['mismatches'] if L else ''}\t"
                  f"{L['indels'] if L else ''}\t{bp:.4f}\t{bb}\n")
    out.close()
    # B-only copies of the S_fam families
    fam = collections.defaultdict(list)
    for g, r in genes.items():
        if r["set"] == "S_fam":
            fam[r["biotype"].split(":", 1)[1]].append(g)
    fpaf = {"mat": paf(a.fam_paf_mat), "pat": paf(a.fam_paf_pat)}
    with open(a.fam_out, "w") as fo:
        fo.write("family\tcopies\tT2d\tT2i\tT1\tTq\tB_only_loci\tT_haplotype_copies\tB_only_detail\n")
        for f_, cps in sorted(fam.items()):
            lifted = [(lift[c]["B_chrom"], int(lift[c]["B_start"]), int(lift[c]["B_end"])) for c in cps
                      if c in lift and lift[c]["B_chrom"] and lift[c]["B_start"]]
            hits = []
            for h in ("mat", "pat"):
                for c in cps:
                    for x in fpaf[h].get(c, []):
                        if x["t"] in b_names and b_names[x["t"]] == h and x["id"] >= 0.90 and x["cov"] >= 0.80:
                            if not any(x["t"] == bc and x["s"] < be and bs < x["e"] for bc, bs, be in lifted):
                                hits.append((x["t"], x["s"], x["e"]))
            loci = []
            for t, s, e in sorted(hits):
                if loci and loci[-1][0] == t and s < loci[-1][2]:
                    loci[-1][2] = max(loci[-1][2], e)
                else:
                    loci.append([t, s, e])
            k = collections.Counter(final.get(c, "NA") for c in cps)
            T = 2 * (k["T2d"] + k["T2i"]) + k["T1"] + len(loci)
            fo.write(f"{f_}\t{len(cps)}\t{k['T2d']}\t{k['T2i']}\t{k['T1']}\t{k['T?']}\t{len(loci)}\t{T}\t"
                     f"{';'.join(f'{t}:{s}-{e}' for t, s, e in loci)}\n")
    c = collections.Counter((sets[g]["set"], v) for g, v in final.items())
    for st in ("S_fam", "S_multi", "S_single", "S_X"):
        print(st, {k[1]: n for k, n in sorted(c.items()) if k[0] == st})


if __name__ == "__main__":
    main()
