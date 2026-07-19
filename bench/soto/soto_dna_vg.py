#!/usr/bin/env python3
"""DNA variation-graph ceiling vs Soto SDs (visual demo). Per flagship Soto family, build a base-level
variation graph (abpoa MSA -> GFA) from the member GENOMIC sequences and colour each copy green
(RNA-recovered) / red (DNA-only). The graph REPRESENTS all copies (the DNA ceiling); it does not
"detect" families. See docs/superpowers/specs/2026-07-18-soto-dna-vg-ceiling-design.md."""
import os
import sys
from collections import defaultdict

# ---- paths ----
BED = "bench/soto/80_fams.chr.bed"
MEMFA = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_members.fa"
DETECT = "bench/soto/soto_member_detection.tsv"
OUT = "/home/juanfra/winloci_scratch/soto_vg"

# ---- colours (match copy_graph.rs) ----
GREEN = "#1e8e3e"   # RNA-recovered
RED = "#d93025"     # DNA-only (K=0 / silent / coverage)
GREY = "#9aa0a6"    # shared / conserved


def parse_family_members(bed_lines, family_id):
    """(gene, chrom, start, end) for every member of family_id, in file order. BED col4 = GENE|ID_k."""
    out = []
    for ln in bed_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        name = f[3]
        if "|" not in name:
            continue
        gene, fam = name.rsplit("|", 1)
        if fam == family_id:
            out.append((gene, f[0], int(f[1]), int(f[2])))
    return out


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def member_seq(fa, chrom, start, end):
    """Genomic sequence for a BED member. soto_members.fa headers are 1-based (start+1)."""
    return fa[f"{chrom}:{start + 1}-{end}"]


def abpoa_msa(seqs):
    """Aligned gapped rows for the members (str inputs; bytes silently yield all-N). Returns the first
    len(seqs) rows of abpoa's column-MSA (excludes any appended consensus row)."""
    import pyabpoa
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa([s.upper() for s in seqs], out_cons=False, out_msa=True)
    rows = list(res.msa_seq)[:len(seqs)]
    return [r.upper() for r in rows]


def msa_to_gfa(rows, names):
    """Column-MSA -> base-level variation-graph GFA. rows: equal-length uppercase gapped strings
    ('-' = gap), one per name. Returns (gfa_text, paths). Deterministic:
      - a maximal run of columns where all rows share the SAME non-gap base -> one shared node (all paths);
      - a maximal run of variant/gap columns -> one node per distinct gap-stripped allele (sorted); a
        member whose slice is all-gap skips the region."""
    assert rows and len({len(r) for r in rows}) == 1, "rows must be non-empty and equal length"
    assert len(rows) == len(names)
    m, L = len(rows), len(rows[0])

    def invariant(j):
        c0 = rows[0][j]
        return c0 != "-" and all(rows[i][j] == c0 for i in range(m))

    # segment columns into maximal same-class runs
    segments, j = [], 0
    while j < L:
        kind = "inv" if invariant(j) else "var"
        k = j + 1
        while k < L and ("inv" if invariant(k) else "var") == kind:
            k += 1
        segments.append((kind, j, k))
        j = k

    nodes, paths, nid = [], {n: [] for n in names}, 0
    for kind, a, b in segments:
        if kind == "inv":
            nid += 1
            node = str(nid)
            nodes.append((node, rows[0][a:b]))
            for n in names:
                paths[n].append(node)
        else:
            allele = {n: rows[i][a:b].replace("-", "") for i, n in enumerate(names)}
            node_of = {}
            for s in sorted(set(v for v in allele.values() if v)):
                nid += 1
                node_of[s] = str(nid)
                nodes.append((str(nid), s))
            for n in names:
                if allele[n]:
                    paths[n].append(node_of[allele[n]])

    links = set()
    for n in names:
        p = paths[n]
        for x, y in zip(p, p[1:]):
            links.add((x, y))

    out = ["H\tVN:Z:1.0"]
    for node, seq in nodes:
        out.append(f"S\t{node}\t{seq}")
    for x, y in sorted(links, key=lambda t: (int(t[0]), int(t[1]))):
        out.append(f"L\t{x}\t+\t{y}\t+\t0M")
    for n in names:
        out.append(f"P\t{n}\t{'+,'.join(paths[n])}+\t*")
    return "\n".join(out) + "\n", paths
