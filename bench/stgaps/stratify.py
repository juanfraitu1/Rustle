#!/usr/bin/env python3
"""FAMILY_LAYER_SPEC Part A: per-stratum (single-copy / multi-copy reference genes) recall and precision of each tool.

usage: stratify.py REFTX.tsv REFTX.paf REF.gtf OUTDIR name=GTF [name=GTF ...]
A transcript is multi-copy when the PAF has a hit covering >= 80% of it at >= 95% identity (col10/col11) that does not
overlap its own locus; a gene is multi-copy when any of its transcripts is. GTF may be .gtf/.gff(.gz) or FLAIR .bed
(BED12, converted here). Recall = reference multi-exon transcripts with a gffcompare `=` query / all, per stratum.
Precision = `=` / multi-exon queries, a query taking the stratum of its nearest reference transcript (`none` if '-').
"""
import gzip
import os
import subprocess
import sys
from collections import Counter, defaultdict

reftx, paf, ref, out = sys.argv[1:5]
arms = [a.split("=", 1) for a in sys.argv[5:]]
os.makedirs(out, exist_ok=True)

gene, loc = {}, {}
for l in open(reftx):
    t, g, c, s, a, b = l.rstrip("\n").split("\t")
    gene[t], loc[t] = g, (c, int(a), int(b))
multi_tx = set()
for l in open(paf):
    f = l.split("\t")
    t, qlen, qs, qe, tn, ts, te, nm, al = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    # identity = 1 - dv: a spliced PAF's block length (col 11) counts intron bases (Amendment 1 of the spec)
    # minimap2 writes dv:f (approximate) without -c and de:f (gap-compressed, from the base alignment) with -c
    dv = next((float(x[5:]) for x in f[12:] if x.startswith(("dv:f:", "de:f:"))), 1.0)
    if t not in loc or 1 - dv < 0.95 or (qe - qs) < 0.8 * qlen:
        continue
    c, a, b = loc[t]
    if tn == c and ts < b and te > a:
        continue
    multi_tx.add(t)
multi_gene = {gene[t] for t in multi_tx}
stratum = {t: ("multi" if gene[t] in multi_gene else "single") for t in gene}
n_ref = Counter(stratum.values())
print(f"reference multi-exon transcripts: {dict(n_ref)}; genes multi {len(multi_gene)} / {len(set(gene.values()))}")
with open(os.path.join(out, "stratum.tsv"), "w") as fo:
    for t in sorted(stratum):
        fo.write(f"{t}\t{gene[t]}\t{stratum[t]}\n")


def to_gtf(src, dst):
    if os.path.exists(dst):  # converted / linked by an earlier run
        return dst
    op = gzip.open if src.endswith(".gz") else open
    if src.endswith(".bed"):
        with open(src) as fi, open(dst, "w") as fo:
            for l in fi:
                f = l.rstrip("\n").split("\t")
                c, s, n, st = f[0], int(f[1]), f[3], f[5]
                sizes = [int(x) for x in f[10].strip(",").split(",")]
                starts = [int(x) for x in f[11].strip(",").split(",")]
                a = f'gene_id "{n}"; transcript_id "{n}";'
                fo.write(f"{c}\tflair\ttranscript\t{s + 1}\t{f[2]}\t.\t{st}\t.\t{a}\n")
                for z, o in zip(sizes, starts):
                    fo.write(f"{c}\tflair\texon\t{s + o + 1}\t{s + o + z}\t.\t{st}\t.\t{a}\n")
        return dst
    if src.endswith(".gz"):
        with op(src, "rt") as fi, open(dst, "w") as fo:
            fo.writelines(fi)
        return dst
    # a link in OUTDIR, so gffcompare writes its .tmap/.refmap here and not next to the tool's file
    if not os.path.lexists(dst):
        os.symlink(os.path.abspath(src), dst)
    return dst


print("arm\tstratum\tref_tx\trecalled\trecall%\tqueries\t=\tprecision%")
for name, path in arms:
    g = to_gtf(path, os.path.join(out, f"{name}.gtf" if not path.endswith(".gff.gz") else f"{name}.gff"))
    pre = os.path.join(out, f"gc_{name}")
    if not os.path.exists(pre + ".stats") and not os.path.exists(pre):
        subprocess.run(["gffcompare", "-r", ref, "-o", pre, g], check=True, capture_output=True)
    tmap = os.path.join(os.path.dirname(g), f"gc_{name}.{os.path.basename(g)}.tmap")
    hit, q, qe = set(), Counter(), Counter()
    for l in open(tmap):
        f = l.rstrip("\n").split("\t")
        if f[0] == "ref_gene_id" or int(f[5]) < 2:
            continue
        s = stratum.get(f[1], "none")
        q[s] += 1
        if f[2] == "=":
            qe[s] += 1
            hit.add(f[1])
    for s in ("single", "multi", "none"):
        r = sum(1 for t in hit if stratum.get(t) == s)
        nr = n_ref.get(s, 0)
        print(f"{name}\t{s}\t{nr}\t{r}\t{100 * r / nr if nr else 0:.1f}\t{q[s]}\t{qe[s]}\t{100 * qe[s] / q[s] if q[s] else 0:.1f}")
