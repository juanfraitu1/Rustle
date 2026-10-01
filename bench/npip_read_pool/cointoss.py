#!/usr/bin/env python3
"""Coin-toss check at the 25 CAT NPIP copies (human chr16, A119b): of the primary reads (-F 2308) whose aligned blocks overlap a copy's
exons, how many are AS-tied genome-wide (second-best AS >= 0.98 x best AS, from the as_table molecules table), and how many tie to
ANOTHER NPIP copy. A primaries-only locus is unsafe only where its reads are mostly coin tosses (the §6z7 r1061 audit, genome-wide:
0.3-0.5% of loci majority-tied). Descriptive; added after docs/PREREG_npip_read_pool_2026-10-01.md.

    python3 bench/npip_read_pool/cointoss.py --copies copies.hsa.tsv --cat-genes genes.tsv --bam A119b.t2t.bam --table molecules.tsv --out x.json
"""
import argparse
import csv
import json
import re
import subprocess


def blocks(pos1, cigar):
    out, p = [], pos1 - 1
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            out.append((p, p + n))
        if op in "MDN=X":
            p += n
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("copies", "cat_genes", "bam", "table", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    cat = {r["gene_id"]: r for r in csv.DictReader(open(a.cat_genes), delimiter="\t")}
    cps = [c for c in csv.DictReader(open(a.copies), delimiter="\t") if c["family"] == "NPIP" and c["chrom"] == "chr16"]
    reads = {}
    for c in cps:
        g = cat[c["isoform_gene"]]
        ex = [tuple(int(x) for x in b.split("-")) for b in g["exon_blocks"].split(",") if b]
        r = subprocess.run(["samtools", "view", "-F", "2308", a.bam, f"chr16:{int(g['start0']) + 1}-{g['end']}"], capture_output=True,
                           text=True, check=True)
        names = set()
        for ln in r.stdout.splitlines():
            f = ln.split("\t", 6)
            if any(min(e, y) > max(s, x) for s, e in blocks(int(f[3]), f[5]) for x, y in ex):
                names.add(f[0])
        reads[c["cid"]] = (c, names)
    want = set().union(*(n for _, n in reads.values()))
    tab = {}
    with open(a.table) as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            name, best, second = ln.split("\t", 3)[:3]
            if name in want:
                tab[name] = (int(best), int(second))
    out = []
    for cid, (c, names) in reads.items():
        tied = [n for n in names if n in tab and tab[n][1] >= 0.98 * tab[n][0] and tab[n][1] > 0]
        out.append(dict(cid=cid, name=c["refseq_name"], primaries=len(names), tied=len(tied),
                        tied_frac=round(len(tied) / len(names), 3) if names else None, untied=len(names) - len(tied)))
    json.dump(out, open(a.out, "w"), indent=0)
    for r in out:
        print(f"{r['name']:15} primaries {r['primaries']:6}  tied {r['tied']:6} ({r['tied_frac']})  untied {r['untied']}")
    tp = sum(r["primaries"] for r in out)
    tt = sum(r["tied"] for r in out)
    print(f"all copies: {tt} / {tp} primaries tied ({tt / tp:.1%}); copies with < 2 untied primaries: "
          f"{[r['name'] for r in out if r['untied'] < 2]}; copies majority-tied: {[r['name'] for r in out if r['tied_frac'] and r['tied_frac'] >= 0.5]}")


if __name__ == "__main__":
    main()
