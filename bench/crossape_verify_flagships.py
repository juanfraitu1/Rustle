#!/usr/bin/env python3
"""Verify the flagship families against three independent signals per species, to catch the sequence-linking/
annotation errors (e.g. TBC1D3 coming out orangutan-high/gorilla-zero) before anything goes on a figure:

  ann   = # RefSeq gene features whose name matches the family (the assembly's annotated paralog count)
  cat   = # de-novo catalog copies overlapping one of those annotated genes (what our RNA catalog detected)
  link  = the sequence-linked cross-ape count (from crossape_seqlinked.tsv)

A family is VALIDATED when cat tracks ann and the branch matches known biology; flagged when they diverge
(cat<<ann => detection/expression gap; link disagrees with cat => a linking/labeling error).
Run: /home/juanfra/miniforge3/bin/python bench/crossape_verify_flagships.py
"""
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from crossape_expansion import annotate  # family -> {'genes': {gene:count}, 'n': ncopies}

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
SPECIES = [
    ("chimp",   f"{D}/PTR_gwcat.copies.tsv", f"{D}/PTR_genomic.gff"),
    ("gorilla", f"{D}/GGO_gwcat.copies.tsv", f"{D}/GGO_genomic.gff"),
    ("orang",   f"{D}/PPY_gwcat.copies.tsv", f"{D}/PPY_genomic.gff"),
]
ORDER = ["chimp", "gorilla", "orang"]
FLAGSHIPS = [  # (label, name-regex, known biology / expected pattern)
    ("MAGEA",    r"^MAGEA\d",    "cancer-testis; expanded across catarrhines (~ancestral)"),
    ("NPIP",     r"^NPIP",       "African-ape/human expansion; ~absent in orangutan"),
    ("NBPF",     r"^NBPF",       "Olduvai/DUF1220; African-ape+human expansion"),
    ("POTE",     r"^POTE",       "African-ape expansion"),
    ("TBC1D3",   r"^TBC1D3",     "African-ape/human (esp. human) expansion; NOT orang-high"),
    ("SRGAP2",   r"^SRGAP2",     "human-biased; few copies in apes"),
    ("ARHGAP11", r"^ARHGAP11",   "human-biased duplication"),
    ("FAM72",    r"^FAM72",      "human/African-ape expansion"),
    ("GSTM",     r"^GSTM\d",     "GSTM cluster; ~ancestral, ~3-5 copies"),
    ("TSPY",     r"^TSPY\d",     "Y, testis; expanded (needs assembled Y + testis expression)"),
    ("RBMY",     r"^RBMY",       "Y, testis; multi-copy (needs Y + testis)"),
    ("DAZ",      r"^DAZ\d",      "Y, testis; DAZ1-4 (needs Y + testis)"),
]


def gene_names(gff):
    out = []
    for line in open(gff):
        if "\tgene\t" not in line:
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"gene=([^;]+)", f[8]) or re.search(r"Name=([^;]+)", f[8])
        if m:
            out.append(m.group(1))
    return out


def load_seqlink():
    sl = defaultdict(lambda: {s: 0 for s in ORDER})
    try:
        for r in csv.DictReader(open("bench/crossape_seqlinked.tsv"), delimiter="\t"):
            sl[r["modal_gene"]] = {s: int(r[f"copies_{s}"]) for s in ORDER}
    except FileNotFoundError:
        pass
    return sl


def main():
    names = {}; cat = {}
    for sp, cp, gff in SPECIES:
        names[sp] = gene_names(gff)
        cat[sp] = annotate(sp, cp, gff)  # family -> genes/n
    sl = load_seqlink()

    print("\n=== flagship verification (ann = annotated paralogs / cat = catalog copies at those genes) ===")
    print(f"{'family':10s} {'chimp(ann/cat)':>16} {'gorilla(ann/cat)':>17} {'orang(ann/cat)':>16}  verdict")
    rows = []
    for label, rgx, bio in FLAGSHIPS:
        pat = re.compile(rgx, re.I)
        cells = {}
        for sp in ORDER:
            ann = sum(1 for n in names[sp] if pat.search(n))
            catn = 0
            for fam, d in cat[sp].items():
                catn += sum(c for g, c in d["genes"].items() if pat.search(g))
            cells[sp] = (ann, catn)
        # link count: seqlinked rows whose modal gene matches
        link = {s: 0 for s in ORDER}
        for g, cc in sl.items():
            if pat.search(g):
                for s in ORDER:
                    link[s] = max(link[s], cc[s])
        # verdict heuristic
        anns = [cells[s][0] for s in ORDER]; cats = [cells[s][1] for s in ORDER]
        verdict = "validated"
        if sum(anns) == 0:
            verdict = "not annotated (absent/LOC-named or no assembly)"
        elif sum(cats) == 0:
            verdict = "annotated but 0 catalog copies (not expressed / detection gap)"
        elif any(cells[s][0] >= 3 and cells[s][1] == 0 for s in ORDER):
            verdict = "DETECTION/LINK GAP (annotated multi-copy but 0 in catalog somewhere)"
        rows.append([label] + [f"{cells[s][0]}/{cells[s][1]}" for s in ORDER]
                    + [f"{link['chimp']}/{link['gorilla']}/{link['orang']}", verdict, bio])
        c = " ".join(f"{cells[s][0]:>7}/{cells[s][1]:<7}" for s in ORDER)
        print(f"  {label:9s} {c}  {verdict}")

    with open("bench/crossape_flagship_verify.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family"] + [f"{s}_ann/cat" for s in ORDER] + ["seqlink_c/g/o", "verdict", "known_biology"])
        w.writerows(rows)
    print("\n  ann/cat per species; seqlink = sequence-linked table's count (to spot linking errors)")
    print("  wrote bench/crossape_flagship_verify.tsv")


if __name__ == "__main__":
    main()
