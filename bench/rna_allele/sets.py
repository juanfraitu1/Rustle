#!/usr/bin/env python3
"""Gene sets for docs/PREREG_rna_allele_haplotype_count_2026-10-01.md from the paralog test:
`minimap2 -c -x splice -N 50 GGO.splice.mmi exonsum.fa` (default -p, as registered).

  own hit      : on the gene's own chromosome, overlapping its span
  paralog hit  : any other hit with identity (matches / block length) >= 0.90 over >= 0.50 of the exon sum
  S_fam        : the 39 NPIP/TBC1D3 copy records
  S_multi      : RefSeq records on autosomes with >= 1 paralog hit
  S_single     : RefSeq records on autosomes with none
  S_X          : RefSeq records on chrX with no paralog hit and no chrY hit at identity >= 0.90 (pseudo-autosomal test)
  excluded     : chrY records; records with no hit at all; chrX records failing the two tests
Also reports, per record, the closest paralog's identity (0 if none), used by truth_classes.py.

    python3 sets.py --genes genes.tsv --paf exonsum.pri.paf --out sets.tsv
"""
import argparse
import collections
import csv


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("genes", "paf", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    genes = {r["gene_id"]: r for r in csv.DictReader(open(a.genes), delimiter="\t")}
    span = {g: (int(r["exons"].split(",")[0].split("-")[0]), int(r["exons"].split(",")[-1].split("-")[1])) for g, r in genes.items()}
    ychrom = next((r["chrom"] for r in genes.values() if r["chrom_label"] == "Y"), None)
    hits = collections.defaultdict(list)
    for ln in open(a.paf):
        f = ln.split("\t")
        hits[f[0]].append((f[5], int(f[7]), int(f[8]), int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))))
    out = open(a.out, "w")
    out.write("gene_id\tname\tset\tchrom_label\tn_paralogs\tbest_paralog_identity\town_hit\tchrY_hit\n")
    cnt = collections.Counter()
    for g, r in genes.items():
        s0, e0 = span[g]
        hs = hits.get(g, [])
        own = [h for h in hs if h[0] == r["chrom"] and h[1] < e0 and s0 < h[2]]
        para = [h for h in hs if h not in own and h[3] >= 0.90 and h[4] >= 0.50]
        best = max((h[3] for h in hs if h not in own and h[4] >= 0.50), default=0.0)
        yhit = any(h[0] == ychrom and h[3] >= 0.90 for h in hs)
        lab = r["chrom_label"]
        if r["set"] == "S_fam":
            st = "S_fam"
        elif lab == "Y" or not hs:
            st = "excluded"
        elif lab == "X":
            st = "S_X" if not para and not yhit else "excluded"
        else:
            st = "S_multi" if para else "S_single"
        cnt[st] += 1
        out.write(f"{g}\t{r['name']}\t{st}\t{lab}\t{len(para)}\t{best:.4f}\t{int(bool(own))}\t{int(yhit)}\n")
    print(dict(cnt))


if __name__ == "__main__":
    main()
