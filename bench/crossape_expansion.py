#!/usr/bin/env python3
"""Cross-ape multi-copy family expansion table. Annotates each species' de-novo families (gw_family_catalog
copies.tsv) with gene names from that species' RefSeq GFF, then links families ACROSS species by shared gene
names (union-find on gene overlap -- robust to the MAGEA1/SRGAP2B naming inconsistency that fuzzy stem
normalization breaks). Each cross-species group -> copy count per ape, so lineage-specific expansions surface.

Depth: chimp 4.0M, orang 5.1M, gorilla ~4M (downsampled to match) -> counts are comparable.
Run: /home/juanfra/miniforge3/bin/python bench/crossape_expansion.py
"""
import csv
import re
from bisect import bisect_right
from collections import defaultdict

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
SPECIES = [  # (name, copies_tsv, gff)  -- only those present are used
    ("chimp",   f"{D}/PTR_gwcat.copies.tsv", f"{D}/PTR_genomic.gff"),
    ("gorilla", f"{D}/GGO_gwcat.copies.tsv", f"{D}/GGO_genomic.gff"),
    ("orang",   f"{D}/PPY_gwcat.copies.tsv", f"{D}/PPY_genomic.gff"),
]
CH, SC, EC = 3, 4, 5  # copies.tsv: family_id=0, chrom=3, start=4, end=5


def gene_index(gff):
    """chrom -> (sorted starts, [(start,end,name)]) for `gene` features."""
    tmp = defaultdict(list)
    for line in open(gff):
        if "\tgene\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"gene=([^;]+)", f[8]) or re.search(r"Name=([^;]+)", f[8])
        if not m:
            continue
        tmp[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1)))
    idx = {}
    for c, ivs in tmp.items():
        ivs.sort()
        idx[c] = ([s for s, _, _ in ivs], ivs)
    return idx


def genes_at(idx, chrom, s, e):
    entry = idx.get(chrom)
    if not entry:
        return []
    starts, ivs = entry
    out = []
    i = bisect_right(starts, e)
    for j in range(i - 1, -1, -1):
        gs, ge, name = ivs[j]
        if ge < s:
            if gs < s - 3_000_000:  # genes are < a few Mb; stop scanning back
                break
            continue
        if not (gs > e or ge < s):
            out.append(name)
    return out


def annotate(name, copies_tsv, gff):
    """family_id -> {'genes': Counter-ish set with counts, 'n': n_copies}."""
    idx = gene_index(gff)
    fams = defaultdict(lambda: {"genes": defaultdict(int), "n": 0})
    r = csv.reader(open(copies_tsv), delimiter="\t"); next(r, None)
    for row in r:
        if len(row) <= EC:
            continue
        fid = row[0]
        try:
            c, s, e = row[CH], int(row[SC]), int(row[EC])
        except ValueError:
            continue
        fams[fid]["n"] += 1
        for g in genes_at(idx, c, s, e):
            if not g.upper().startswith("LOC"):     # skip uncharacterized LOC ids for the cross-species link
                fams[fid]["genes"][g] += 1
    return fams


def main():
    present = [(n, c, g) for (n, c, g) in SPECIES if __import__("os").path.exists(c)]
    print(f"species with catalogs: {[n for n,_,_ in present]}")
    per = {}   # species -> {fid: {...}}
    nodes = []  # (species, fid)
    gene_to_nodes = defaultdict(list)
    for name, cp, gff in present:
        fams = annotate(name, cp, gff)
        per[name] = fams
        for fid, d in fams.items():
            node = (name, fid); nodes.append(node)
            for g in d["genes"]:
                gene_to_nodes[g].append(node)
        print(f"  {name}: {len(fams)} families annotated ({sum(1 for d in fams.values() if d['genes'])} with a named gene)")

    # union-find: link families (across species) that share a gene name
    parent = {nd: nd for nd in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for g, nds in gene_to_nodes.items():
        for k in range(1, len(nds)):
            a, b = find(nds[0]), find(nds[k])
            if a != b:
                parent[a] = b

    comp = defaultdict(list)
    for nd in nodes:
        comp[find(nd)].append(nd)

    rows = []
    for members in comp.values():
        cnt = {n: 0 for n, _, _ in present}
        genes = defaultdict(int)
        for (sp, fid) in members:
            cnt[sp] += per[sp][fid]["n"]
            for g, c in per[sp][fid]["genes"].items():
                genes[g] += c
        if not genes:
            continue  # unnamed (LOC-only / intergenic) family group
        label = max(genes, key=genes.get)
        # family-family prefix for readability (strip trailing paralog digits/letters, keep known stems)
        stem = re.sub(r"\d+[A-Z]?$", "", label) or label
        rows.append((stem, label, cnt, sorted(genes, key=genes.get, reverse=True)[:6]))

    order = [n for n, _, _ in present]
    rows.sort(key=lambda r: -max(r[2].values()))
    print(f"\n=== cross-ape multi-copy family expansions (copies per species; {' / '.join(order)}) ===")
    print(f"{'family':14s} " + " ".join(f"{o:>7s}" for o in order) + "   modal_gene  sample_paralogs")
    for stem, label, cnt, sample in rows[:40]:
        counts = " ".join(f"{cnt[o]:>7d}" for o in order)
        print(f"  {stem:12s} {counts}   {label:11s} {','.join(sample)[:44]}")

    with open("bench/crossape_expansion.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_stem", "modal_gene"] + [f"copies_{o}" for o in order] + ["sample_paralogs"])
        for stem, label, cnt, sample in rows:
            w.writerow([stem, label] + [cnt[o] for o in order] + [",".join(sample)])
    print(f"\n  wrote bench/crossape_expansion.tsv ({len(rows)} named cross-ape family groups)")


if __name__ == "__main__":
    main()
