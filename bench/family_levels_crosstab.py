#!/usr/bin/env python
"""Authoritative cross-tab — each RNA (read-conflict) family scored at all three SEPARATE levels on ONE catalog:
  PROTEIN  : copies are a single protein-paralog component (longest ORF -> mmseqs, fident>=0.50, cov>=0.50)
  DNA      : copies are mutually homologous (Soto-style pairwise asm20, >=200 bp block at >=90% id)
  SEDEF    : two copies are the two ends of an INDEPENDENT segmental-duplication pair (final.bed) — the most
             orthogonal check (the segdup map is built from DNA alone, never sees the RNA catalog).

Answers: how many read-conflict families are genuine paralog families (coding / segdup) vs identifiability-only
(linked by shared non-coding / one exon / repeat)? The Venn over {protein, DNA, SEDEF} is the FP/FN picture.

Run: MINIFORGE python bench/family_levels_crosstab.py [catalog_prefix=gw_off] [sedef_bed=../final.bed]
"""
import bisect
import os
import subprocess
import sys
import tempfile

import pysam
from Bio.Seq import Seq

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = f"{SCRATCH}/GGO.fasta"
MMSEQS = os.environ.get("RUSTLE_MMSEQS", "/home/juanfra/miniforge3/bin/mmseqs")
MM2 = os.environ.get("RUSTLE_MINIMAP2", "/home/juanfra/miniforge3/bin/minimap2")
PREFIX = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/gw_off"
SEDEF = sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"
MIN_AA, MIN_FIDENT, MIN_COV = 40, 0.50, 0.50
DNA_BLOCK_BP, DNA_ID = 200, 0.90


def longest_orf_aa(nt):
    best = ""
    s = Seq(nt.upper().replace("N", "A"))
    for seq in (s, s.reverse_complement()):
        for f in range(3):
            for pep in str(seq[f:].translate()).split("*"):
                if len(pep) > len(best):
                    best = pep
    return best if len(best) >= MIN_AA else ""


def comps(nodes, edges):
    par = {n: n for n in nodes}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges:
        if a in par and b in par:
            par[find(a)] = find(b)
    seen = {find(n) for n in nodes}
    return len(seen)


def load_families():
    """famid -> list of (spliced_seq, (chrom,start,end)) from <prefix>.copies.fa headers."""
    fams = {}
    fa = pysam.FastaFile(f"{PREFIX}.copies.fa")
    for ref in fa.references:
        p = ref.split("|")                       # FAMID|idx|chrom:start-end|strand|nexon
        fid = p[0]
        chrom, span = p[2].split(":")
        s, e = span.split("-")
        fams.setdefault(fid, []).append((fa.fetch(ref), (chrom, int(s), int(e))))
    return {k: v for k, v in fams.items() if len(v) >= 2}


def protein_confirmed(fams):
    """family -> single protein-paralog component?"""
    node_fam, fam_nodes = {}, {f: [] for f in fams}
    with tempfile.TemporaryDirectory() as td:
        qf = os.path.join(td, "q.fa")
        with open(qf, "w") as out:
            for fi, fid in enumerate(fams):
                for pos, (seq, _loc) in enumerate(fams[fid]):
                    aa = longest_orf_aa(seq)
                    if aa:
                        nid = f"f{fi}_c{pos}"
                        out.write(f">{nid}\n{aa}\n"); node_fam[nid] = fid; fam_nodes[fid].append(nid)
        res = os.path.join(td, "r.m8")
        r = subprocess.run([MMSEQS, "easy-search", qf, qf, res, os.path.join(td, "tmp"),
                            "-s", "7.5", "--threads", "4", "--format-output", "query,target,fident,qcov,tcov"],
                           capture_output=True, text=True)
        edges = {f: [] for f in fams}
        if r.returncode == 0 and os.path.exists(res):
            for line in open(res):
                q, t, fi, qc, tc = line.rstrip("\n").split("\t")
                if q == t or float(fi) < MIN_FIDENT or min(float(qc), float(tc)) < MIN_COV:
                    continue
                if node_fam.get(q) == node_fam.get(t) and node_fam.get(q) is not None:
                    edges[node_fam[q]].append((q, t))
    out = {}
    for fid in fams:
        nodes = fam_nodes[fid]
        out[fid] = (len(nodes) >= 2 and comps(nodes, edges[fid]) == 1)
    return out


def dna_confirmed(fams):
    """family -> copies form one homology component via pairwise asm20 (>=200bp block, >=90% id)?"""
    fa = pysam.FastaFile(GENOME)
    out = {}
    for fid, members in fams.items():
        n = len(members); nodes = list(range(n)); edges = []
        with tempfile.TemporaryDirectory() as td:
            seqs = [fa.fetch(*m[1]) for m in members]
            import itertools
            for i, j in itertools.combinations(range(n), 2):
                rp = os.path.join(td, "r.fa"); qp = os.path.join(td, "q.fa")
                open(rp, "w").write(f">a\n{seqs[i]}\n"); open(qp, "w").write(f">b\n{seqs[j]}\n")
                o = subprocess.run([MM2, "-cx", "asm20", "--eqx", rp, qp], capture_output=True, text=True)
                ok = False
                for ln in o.stdout.splitlines():
                    ff = ln.split("\t")
                    if len(ff) < 12:
                        continue
                    aln_bp = int(ff[3]) - int(ff[2])
                    de = next((float(x[5:]) for x in ff[12:] if x.startswith("de:f:")), 1.0)
                    if aln_bp >= DNA_BLOCK_BP and (1 - de) >= DNA_ID:
                        ok = True; break
                if ok:
                    edges.append((i, j))
        out[fid] = (comps(nodes, edges) == 1)
    return out


def load_sedef(path):
    """per-chrom sorted intervals tagged (pair_id, side). Each final.bed line = one segdup pair (two regions)."""
    idx = {}
    pid = 0
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        try:
            ca, sa, ea, cb, sb, eb = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
        except (ValueError, IndexError):
            continue
        idx.setdefault(ca, []).append((sa, ea, pid, 0))
        idx.setdefault(cb, []).append((sb, eb, pid, 1))
        pid += 1
    for c in idx:
        idx[c].sort()
    return idx


def sedef_confirmed(fams, sedef):
    """family -> two copies are the two ends of the SAME segdup pair (copy_i side0, copy_j side1)?"""
    starts = {c: [iv[0] for iv in sedef[c]] for c in sedef}
    def overlaps(chrom, s, e):
        out = set()
        arr = sedef.get(chrom)
        if not arr:
            return out
        st = starts[chrom]
        i = bisect.bisect_right(st, e)            # intervals starting at/before e
        # scan back while intervals can still overlap (segdups are <~ a few Mb; cap the back-scan generously)
        j = i - 1
        while j >= 0 and st[j] >= s - 5_000_000:
            iv_s, iv_e, pid, side = arr[j]
            if iv_s < e and iv_e > s:
                out.add((pid, side))
            j -= 1
        return out
    out = {}
    for fid, members in fams.items():
        per_copy = [overlaps(*m[1]) for m in members]
        linked = False
        # any pair_id present on side0 in one copy AND side1 in another copy
        side0 = {}; side1 = {}
        for ci, ps in enumerate(per_copy):
            for pid, side in ps:
                (side0 if side == 0 else side1).setdefault(pid, set()).add(ci)
        for pid in set(side0) & set(side1):
            # the segdup pair `pid` has its two ends covered by >=2 DISTINCT family copies
            if len(side0[pid] | side1[pid]) >= 2:
                linked = True; break
        # also accept: >=2 copies fall inside segdup regions at all (weaker), tracked separately
        n_in_segdup = sum(1 for ps in per_copy if ps)
        out[fid] = (linked, n_in_segdup >= 2)
    return out


def main():
    fams = load_families()
    print(f"catalog {os.path.basename(PREFIX)}: {len(fams)} read-conflict families (>=2 copies)")
    print("scoring PROTEIN (mmseqs)..."); prot = protein_confirmed(fams)
    print("scoring DNA (pairwise asm20)..."); dna = dna_confirmed(fams)
    print(f"scoring SEDEF (final.bed)..."); sed_raw = sedef_confirmed(fams, load_sedef(SEDEF))
    sed = {f: v[0] for f, v in sed_raw.items()}            # strict: two copies = the two ends of a segdup pair
    sed_in = {f: v[1] for f, v in sed_raw.items()}         # weak: >=2 copies inside segdup regions

    from collections import Counter
    cells = Counter()
    for fid in fams:
        cells[(prot[fid], dna[fid], sed[fid])] += 1
    N = len(fams)
    print(f"\n=== CROSS-TAB over {N} families (protein, DNA, SEDEF-pair) ===")
    print(f"{'protein':>8} {'DNA':>5} {'SEDEF':>6} | {'count':>5}  {'%':>5}")
    for combo in sorted(cells, reverse=True):
        p, d, s = combo
        print(f"{str(p):>8} {str(d):>5} {str(s):>6} | {cells[combo]:>5}  {100*cells[combo]/N:>4.0f}%")
    pc = sum(1 for f in fams if prot[f]); dc = sum(1 for f in fams if dna[f])
    sc = sum(1 for f in fams if sed[f]); si = sum(1 for f in fams if sed_in[f])
    anyp = sum(1 for f in fams if prot[f] or dna[f] or sed[f])
    allp = sum(1 for f in fams if prot[f] and dna[f] and sed[f])
    print(f"\nmarginals: protein {pc}/{N} ({100*pc/N:.0f}%) | DNA {dc}/{N} ({100*dc/N:.0f}%) | "
          f"SEDEF-pair {sc}/{N} ({100*sc/N:.0f}%) | SEDEF-in-segdup {si}/{N} ({100*si/N:.0f}%)")
    print(f"confirmed by ANY level (a real paralog family): {anyp}/{N} ({100*anyp/N:.0f}%)")
    print(f"confirmed by ALL three: {allp}/{N} ({100*allp/N:.0f}%)")
    print(f"confirmed by NONE (identifiability-only — shared non-coding/exon/repeat): "
          f"{N-anyp}/{N} ({100*(N-anyp)/N:.0f}%)")


if __name__ == "__main__":
    main()
