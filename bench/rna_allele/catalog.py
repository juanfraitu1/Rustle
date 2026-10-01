#!/usr/bin/env python3
"""Catalog for the shipped O2 (`copy_assign --families`) in docs/PREREG_rna_allele_haplotype_count_2026-10-01.md, RNA step 3.

Families: NPIP and TBC1D3 (the S_fam copies), and S_multi paralog groups = connected components of "G has a paralog hit (identity >= 0.90,
>= 0.50 of its exon sum, not its own locus) that overlaps the exons of gene record H on the hit's genomic strand". RefSeq records that
overlap an S_fam copy on the same strand are left out of the S_multi groups (the S_fam family already holds that locus). A copy = the
record's exon union; its sequence = the exon sum on the transcript strand (copies.fa header `>{fam}|{idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`).

    python3 catalog.py --genes genes.tsv --sets sets.tsv --paf para.paf --exonsum exonsum.fa --out-prefix fibro_cat
"""
import argparse
import bisect
import collections
import csv


def blocks(s):
    return [tuple(int(x) for x in b.split("-")) for b in s.split(",") if b]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("genes", "sets", "paf", "exonsum", "out_prefix"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    genes = {r["gene_id"]: r for r in csv.DictReader(open(a.genes), delimiter="\t")}
    sets = {r["gene_id"]: r["set"] for r in csv.DictReader(open(a.sets), delimiter="\t")}
    ex = {g: blocks(r["exons"]) for g, r in genes.items()}
    fam_copies = [g for g, s in sets.items() if s == "S_fam"]
    # index gene records by (chrom, strand) for overlap lookups
    idx = collections.defaultdict(list)
    for g, s in sets.items():
        if s in ("S_multi", "S_single", "S_fam"):
            r = genes[g]
            idx[(r["chrom"], r["strand"])].append((ex[g][0][0], ex[g][-1][1], g))
    for v in idx.values():
        v.sort()
    starts = {k: [x[0] for x in v] for k, v in idx.items()}

    def overlapping(chrom, strand, s, e, exons=None):
        k = (chrom, strand)
        if k not in idx:
            return []
        out = []
        i = bisect.bisect_left(starts[k], e)
        for gs, ge, g in idx[k][:i]:
            if ge > s and (exons is None or any(min(y, q) > max(x, p) for x, y in ex[g] for p, q in exons)):
                out.append(g)
        return out
    shadowed = {g for c in fam_copies for g in overlapping(genes[c]["chrom"], genes[c]["strand"], ex[c][0][0], ex[c][-1][1], ex[c])
                if sets.get(g) != "S_fam"}
    adj = collections.defaultdict(set)
    for ln in open(a.paf):
        f = ln.split("\t")
        g = f[0]
        if sets.get(g) != "S_multi" or g in shadowed:
            continue
        ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if ident < 0.90 or cov < 0.50:
            continue
        r = genes[g]
        hs = r["strand"] if f[4] == "+" else ("-" if r["strand"] == "+" else "+")
        for h in overlapping(f[5], hs, int(f[7]), int(f[8]), [(int(f[7]), int(f[8]))]):
            if h != g and h not in shadowed and sets.get(h) in ("S_multi", "S_single"):
                adj[g].add(h); adj[h].add(g)
    # components
    seen, comps = set(), []
    for g in sorted(adj):
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop(); comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        comps.append(sorted(comp))
    fams = [("NPIP", sorted(c for c in fam_copies if genes[c]["biotype"] == "S_fam:NPIP")),
            ("TBC1D3", sorted(c for c in fam_copies if genes[c]["biotype"] == "S_fam:TBC1D3"))]
    fams += [(f"SM{i}", c) for i, c in enumerate(sorted(comps, key=lambda c: (-len(c), c[0])))]
    seqs, cur = {}, None
    for ln in open(a.exonsum):
        if ln[0] == ">":
            cur = ln[1:].strip(); seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    with open(a.out_prefix + ".copies.tsv", "w") as t, open(a.out_prefix + ".copies.fa", "w") as fa, \
            open(a.out_prefix + ".regions.txt", "w") as rg:
        t.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        spans = collections.defaultdict(list)
        for fid, cps in fams:
            for i, g in enumerate(cps):
                r, bl = genes[g], ex[g]
                s, e = bl[0][0], bl[-1][1]
                t.write(f"{fid}\t{i}\t{g}\t{r['chrom']}\t{s}\t{e}\t{len(bl)}\t{r['strand']}\t0\t{','.join(f'{x}-{y}' for x, y in bl)}\n")
                fa.write(f">{fid}|{i}|{r['chrom']}:{s}-{e}|{r['strand']}|nexon={len(bl)}\n{''.join(seqs[g])}\n")
                spans[r["chrom"]].append([max(0, s - 1000), e + 1000])
        for c, iv in spans.items():
            iv.sort(); m = []
            for x, y in iv:
                if m and x <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], y)
                else:
                    m.append([x, y])
            for x, y in m:
                rg.write(f"{c}:{x + 1}-{y}\n")
    n = sum(len(c) for _, c in fams)
    print(f"{len(fams)} families, {n} copies; S_multi genes in a family: {sum(len(c) for c in comps)} of "
          f"{sum(1 for s in sets.values() if s == 'S_multi')}; shadowed by S_fam: {len(shadowed)}")


if __name__ == "__main__":
    main()
