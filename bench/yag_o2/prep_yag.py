#!/usr/bin/env python3
"""YAG families of a catalog (user 2026-10-12: "we only need this for the YAGs, so DAZ, PRY, etc."):
catalog families with >= 1 copy overlapping a Y ampliconic gene of the annotation, named by those genes; their copies
table / FASTA; regions = each family's span per chromosome +- PAD, merged (not whole chromosomes).

    prep_yag.py CATALOG_PREFIX GFF YCONTIG_IN_GFF YCHROM OUTDIR
"""
import collections, json, os, re, sys

YAG = re.compile(r"^(DAZ|RBMY|TSPY|VCY|XKRY|CDY|HSFY|BPY2|PRY|PRORY)")
PAD = 50_000
MAX_COPIES = 100


def main(prefix, gff, ycontig, ychrom, out):
    os.makedirs(out, exist_ok=True)
    genes = []
    for line in open(gff):
        f = line.split("\t")
        if len(f) > 8 and f[0] == ycontig and f[2] == "gene":
            m = re.search(r"gene=([^;]+)", f[8])
            if m and YAG.match(m.group(1)):
                genes.append((int(f[3]) - 1, int(f[4]), m.group(1)))
    rows = [l.rstrip("\n").split("\t") for l in open(prefix + ".fam.copies.tsv")]
    head, rows = rows[0], rows[1:]
    ix = {h: i for i, h in enumerate(head)}
    names = collections.defaultdict(set)
    for r in rows:
        if r[ix["chrom"]] == ychrom:
            s, e = int(r[ix["start"]]), int(r[ix["end"]])
            for gs, ge, g in genes:
                if s < ge and gs < e:
                    names[r[ix["family_id"]]].add(g)
    size = collections.Counter(r[ix["family_id"]] for r in rows)
    keep = {f for f in names if 2 <= size[f] <= MAX_COPIES}
    with open(f"{out}/copies.tsv", "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            if r[ix["family_id"]] in keep:
                fh.write("\t".join(r) + "\n")
    with open(prefix + ".fam.copies.fa") as fi, open(f"{out}/copies.fa", "w") as fo:
        on = False
        for line in fi:
            if line.startswith(">"):
                on = line[1:].split("|")[0] in keep
            if on:
                fo.write(line)
    # one span per (family, chromosome) +- PAD: copy_assign binds a family to the ONE region containing its span
    fam_span = {}
    for r in rows:
        if r[ix["family_id"]] in keep:
            k = (r[ix["family_id"]], r[ix["chrom"]])
            s, e = int(r[ix["start"]]), int(r[ix["end"]])
            fam_span[k] = (min(s, fam_span[k][0]), max(e, fam_span[k][1])) if k in fam_span else (s, e)
    spans = collections.defaultdict(list)
    for (f, chrom), (s, e) in fam_span.items():
        spans[chrom].append((max(1, s - PAD), e + PAD))
    regions = []
    for chrom, iv in sorted(spans.items()):
        iv.sort()
        cur = list(iv[0])
        for s, e in iv[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                regions.append((chrom, *cur))
                cur = [s, e]
        regions.append((chrom, *cur))
    with open(f"{out}/regions.txt", "w") as fh:
        fh.writelines(f"{c}:{s}-{e}\n" for c, s, e in regions)
    fam_names = {f: sorted(names[f]) for f in sorted(keep, key=lambda x: -size[x])}
    json.dump({"families": fam_names, "sizes": {f: size[f] for f in keep},
               "dropped_arrays": {f: size[f] for f in names if size[f] > MAX_COPIES},
               "dropped_single": sorted(f for f in names if size[f] < 2)}, open(f"{out}/yag_families.json", "w"), indent=1)
    print(f"{len(keep)} YAG families, {sum(size[f] for f in keep)} copies, {len(regions)} regions")
    for f, g in fam_names.items():
        print(f"  {f} ({size[f]} copies): {','.join(g)}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
