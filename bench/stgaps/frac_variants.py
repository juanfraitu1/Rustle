#!/usr/bin/env python3
"""Isoform-fraction diagnostic + variants on an `--assembly-polish none` GTF.

usage: frac_variants.py RAW.gtf REF.gtf RAW_TMAP OUTDIR
RAW_TMAP: gffcompare .tmap of RAW vs REF (class '=' marks a reference-matching intron chain).
Prints (1) how many shipped `frac` drops share no exonic bp with their gene's best transcript and how
many of those were '=' matches, (2) gffcompare chain/transcript Sn/Pr for each variant.
"""
import os
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polish_port as P  # noqa: E402

raw, ref, tmap, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
data = P.load(raw)
exons, key, reads, gene, chain, span = data
cls = {}
for line in open(tmap):
    f = line.split("\t")
    if f[0] != "ref_gene_id":
        cls[f[4]] = f[2]

base = P.polish(data)
frac = sorted(t for t, s in base.items() if s == "frac")
kept_pre = [t for t in gene if base.get(t) in (None, "frac", "ret")]
best = {}
for t in kept_pre:
    g = gene[t]
    if g not in best or reads.get(t, 0) > reads.get(best[g], 0) or \
            (reads.get(t, 0) == reads.get(best[g], 0) and t < best[g]):
        best[g] = t


def share(a, b):
    return key[a] == key[b] and any(x < d and c < y for x, y in exons[a] for c, d in exons[b])


disj = [t for t in frac if not share(t, best[gene[t]])]
print(f"frac drops {len(frac)}; '=' among them {sum(cls.get(t) == '=' for t in frac)}")
print(f"  no exonic overlap with gene best: {len(disj)}; '=' among them {sum(cls.get(t) == '=' for t in disj)}")
genes_multi = defaultdict(set)
for t in kept_pre:
    genes_multi[gene[t]].add(t)

VARIANTS = [
    ("shipped", dict()),
    ("local", dict(frac_mode="local")),
    ("credit", dict(credit_ism=True)),
    ("local_credit", dict(frac_mode="local", credit_ism=True)),
    ("nofrac", dict(frac=0.0)),
]
print("variant\tkept\tchain_match\tchainSn\tchainPr\ttxSn\ttxPr")
for name, kw in VARIANTS:
    s = P.polish(data, **kw)
    gtf = os.path.join(out, f"{name}.gtf")
    with open(gtf, "w") as fo:
        for line in open(raw):
            f = line.split("\t")
            if len(f) >= 9:
                t = P.attr(f[8], "transcript_id")
                if t is not None and t in s:
                    continue
            fo.write(line)
    pre = os.path.join(out, f"gc_{name}")
    subprocess.run(["gffcompare", "-r", ref, "-o", pre, gtf], check=True, capture_output=True)
    st = open(pre + ".stats").read().splitlines()
    get = lambda lab: [x.strip() for x in next(l for l in st if l.strip().startswith(lab)).split(":")[1].split("|")[:2]]
    m = next(l for l in st if "Matching intron chains" in l).split(":")[1].strip()
    print(name, len(exons) - len(s), m, *get("Intron chain level"), *get("Transcript level"), sep="\t")
