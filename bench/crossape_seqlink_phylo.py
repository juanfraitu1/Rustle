#!/usr/bin/env python3
"""Sequence-based cross-ape family linking + light phylogenetic dating of expansions.

(1) LINK families across chimp/gorilla/orang by CONSENSUS SEQUENCE homology (minimap2 all-vs-all; cross-
    species copies aligning at id>=0.90, cov>=0.50) -> robust to the annotation-naming problem that made the
    gene-name table's zeros untrustworthy. Union-find over (species,family) nodes; a component = one ortholog
    group with a reliable copy count per ape.
(2) DATE each expansion on the great-ape species tree ((chimp,gorilla),orang) by Dollo parsimony on the
    multi-copy state (>=2 copies = expanded): the branch where multi-copy arose = WHEN the expansion started.
    Divergence times: African-ape ancestor ~9 Mya (chimp/gorilla split), great-ape ancestor ~17 Mya (orang).

Counts are RNA-expressed copies at matched ~4M depth (comparable across apes); a proxy for genomic famCN.
Run: /home/juanfra/miniforge3/bin/python bench/crossape_seqlink_phylo.py
"""
import csv
import os
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from crossape_expansion import annotate  # reuse GFF gene annotation

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
S = "/home/juanfra/winloci_scratch"
SPECIES = [
    ("human",   f"{D}/HSA_gwcat.copies.fa", f"{D}/HSA_gwcat.copies.tsv", f"{D}/HSA_genomic.gff"),
    ("chimp",   f"{D}/PTR_gwcat.copies.fa", f"{D}/PTR_gwcat.copies.tsv", f"{D}/PTR_genomic.gff"),
    ("gorilla", f"{D}/GGO_gwcat.copies.fa", f"{D}/GGO_gwcat.copies.tsv", f"{D}/GGO_genomic.gff"),
    ("orang",   f"{D}/PPY_gwcat.copies.fa", f"{D}/PPY_gwcat.copies.tsv", f"{D}/PPY_genomic.gff"),
]
SPECIES = [s for s in SPECIES if os.path.exists(s[1])]  # use whichever catalogs exist
ORDER = [s[0] for s in SPECIES]
ID_MIN, COV_MIN = 0.90, 0.50


def load_copies(fa):
    fams = defaultdict(list); fam = None; buf = []
    def flush():
        if fam is not None and buf:
            fams[fam].append("".join(buf))
    for line in open(fa):
        if line.startswith(">"):
            flush(); buf = []
            fam = line[1:].strip().split("|")[0]
        else:
            buf.append(line.strip())
    flush()
    return fams


BRANCHES = ["ancestral great-ape (>17 Mya, pre-orangutan)", "African great-ape ancestor (~9-17 Mya)",
            "human-chimp ancestor (~6.5-9 Mya)", "human lineage (<6.5 Mya)", "chimp lineage (<6.5 Mya)",
            "gorilla lineage (<9 Mya)", "orangutan lineage (<17 Mya)",
            "ancestral (mixed pattern; contraction or linking gap)"]

def origin_branch(cnt):
    """Dollo parsimony on (((human,chimp),gorilla),orang): where did the multi-copy (>=2) state arise?
    With human present, human-specific expansions (SRGAP2C etc.) become MEASURED, not inferred."""
    H = cnt.get("human", 0) >= 2; C = cnt.get("chimp", 0) >= 2
    G = cnt.get("gorilla", 0) >= 2; O = cnt.get("orang", 0) >= 2
    if H and C and G and O:         return "ancestral great-ape (>17 Mya, pre-orangutan)"
    if H and C and G and not O:     return "African great-ape ancestor (~9-17 Mya)"
    if H and C and not G and not O: return "human-chimp ancestor (~6.5-9 Mya)"
    if H and not C and not G and not O: return "human lineage (<6.5 Mya)"
    if C and not H and not G and not O: return "chimp lineage (<6.5 Mya)"
    if G and not H and not C and not O: return "gorilla lineage (<9 Mya)"
    if O and not H and not C and not G: return "orangutan lineage (<17 Mya)"
    if sum([H, C, G, O]) >= 2:      return "ancestral (mixed pattern; contraction or linking gap)"
    return "single-copy (no expansion)"


def main():
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    famcount = {}
    allfa = f"{S}/crossape_allcopies.fa"
    with open(allfa, "w") as out:
        for sp, fa, _, _ in SPECIES:
            for fam, seqs in load_copies(fa).items():
                famcount[(sp, fam)] = len(seqs)
                parent.setdefault(f"{sp}::{fam}", f"{sp}::{fam}")
                for i, s in enumerate(seqs):
                    if len(s) >= 100:
                        out.write(f">{sp}::{fam}::{i}\n{s}\n")
    print(f"copies: {sum(famcount.values())} across {len(famcount)} (species,family) nodes -> minimap2 all-vs-all...")
    paf = f"{S}/crossape_allcopies.paf"
    subprocess.run(f"{MM2} -cx asm20 -N20 -p0.3 -t6 {allfa} {allfa} > {paf} 2>/dev/null", shell=True, check=True)

    n_link = 0
    for line in open(paf):
        f = line.split("\t")
        if len(f) < 12:
            continue
        q = f[0].split("::"); t = f[5].split("::")
        if q[0] == t[0]:
            continue  # within-species: ignore (we only link ACROSS species)
        qlen, tlen, matches, alnlen = int(f[1]), int(f[6]), int(f[9]), int(f[10])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        idp = (1 - de) if de is not None else (matches / alnlen if alnlen else 0)
        cov = alnlen / min(qlen, tlen) if min(qlen, tlen) else 0
        if idp >= ID_MIN and cov >= COV_MIN:
            union(f"{q[0]}::{q[1]}", f"{t[0]}::{t[1]}"); n_link += 1
    print(f"  {n_link} cross-species copy alignments (id>={ID_MIN}, cov>={COV_MIN}) linking families")

    print("  annotating families with gene names (GFF overlap)...")
    ann = {sp: annotate(sp, cp, gff) for sp, _, cp, gff in SPECIES}
    def genes_of(sp, fam):
        return ann.get(sp, {}).get(fam, {}).get("genes", {})

    comp = defaultdict(list)
    for nd in list(parent):
        comp[find(nd)].append(nd)

    rows = []
    for members in comp.values():
        cnt = {s: 0 for s in ORDER}; genes = defaultdict(int)
        for nd in members:
            sp, fam = nd.split("::")
            cnt[sp] += famcount.get((sp, fam), 0)
            for g, c in genes_of(sp, fam).items():
                genes[g] += c
        if max(cnt.values()) < 2:
            continue
        label = max(genes, key=genes.get) if genes else "(unnamed)"
        rows.append((label, cnt, origin_branch(cnt), sorted(genes, key=genes.get, reverse=True)[:5]))

    rows.sort(key=lambda r: -max(r[1].values()))
    with open("bench/crossape_seqlinked.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["modal_gene"] + [f"copies_{s}" for s in ORDER] + ["origin_branch", "sample_paralogs"])
        for label, cnt, branch, sample in rows:
            w.writerow([label] + [cnt[s] for s in ORDER] + [branch, ",".join(sample)])

    print(f"\n=== SEQUENCE-LINKED cross-ape expansions ({' / '.join(ORDER)}) ===")
    print(f"{'family':13s} " + "".join(f"{s:>8}" for s in ORDER) + "   origin branch")
    for label, cnt, branch, _ in rows[:32]:
        print(f"  {label:11s} " + "".join(f"{cnt[s]:>8}" for s in ORDER) + f"   {branch}")

    print(f"\n=== WHEN did expansions arise? (families per branch of the ape tree) ===")
    by_branch = defaultdict(int)
    for _, _, branch, _ in rows:
        by_branch[branch] += 1
    for br in BRANCHES:
        if by_branch.get(br):
            print(f"  {by_branch[br]:3d}  {br}")
    print(f"\n  wrote bench/crossape_seqlinked.tsv ({len(rows)} expanded ortholog groups)")


if __name__ == "__main__":
    main()
