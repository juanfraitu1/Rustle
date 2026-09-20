#!/usr/bin/env python3
"""family_def_vg_coherence.py — VG-style coherence + NON-LINEAR structure of each family.

Two jobs the pairwise-homology family definition cannot do:
  (A) COHERENCE: do a family's members thread ONE shared collinear backbone, or is the family held
      together non-coherently (a residual bridge)? Backbone-threading fraction is the intrinsic
      family-validity signal (no external protein needed) that should flag the ~10 low-density
      residual-bridge families.
  (B) NON-LINEAR structure the collinear homology pipeline is BLIND to, but a real variation graph
      models: INVERSION (a segment in reverse-complement in some copies), TANDEM REPEAT / copy-number
      (a segment matching the backbone more than once), EXON SHUFFLE (shared segments in a different
      order). We currently flatten all of these.

Method (annotation cDNA + minimap2; `vg` available for a fuller graph later): per family, reference =
longest-cDNA member; align every member to it with minimap2 -c -N 10 -p 0.1 (secondaries ON, so
repeats/shuffles surface as multiple chains). Parse ALL PAF blocks per member and compute:
  coh_fwd   : forward-strand collinear coverage of the reference (the threading fraction)
  inv_frac  : reverse-strand coverage fraction (INVERSION signal)
  rep_hits  : reference bp covered by >=2 distinct query blocks, /ref (TANDEM-REPEAT / CNV signal)
  shuffle   : forward blocks out of collinear order (Kendall-tau discordance > 0) (EXON-SHUFFLE signal)
  n_blocks  : number of merged alignment blocks (fragmentation)

Run: /home/juanfra/miniforge3/bin/python bench/family_def_vg_coherence.py
Out: bench/family_def_vg_coherence.{tsv,json}
"""
import collections
import json
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_build import (protein_whole_edges, load_coord_strand, SEG_BT, LOCUS_WIN)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import build_family_graph, load_intron_index
from family_def_strand_probe import edge_strand
from family_def_protein_validate import gene_biotype

import networkx as nx

GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
EXON_INDEX = "/home/juanfra/winloci_scratch/gene_exon_index.json"
EXON = None
CACHE = None
MIN_BLOCK = 50          # ignore alignment blocks shorter than this (bp)


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def build_cdna(name, fa):
    v = EXON.get(name)
    if v is None or not v["exons"]:
        return None
    blocks = [fa.fetch(v["chrom"], a - 1, b).upper() for a, b in v["exons"]]
    if v["strand"] == "-":
        blocks = [revcomp(b) for b in reversed(blocks)]
    return "".join(blocks)


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def analyze_member(ref_seq, mem_seq, tmp):
    """minimap2 mem->ref, parse all PAF blocks, return non-linear features vs ref length."""
    rp = os.path.join(tmp, "r.fa"); mp = os.path.join(tmp, "m.fa")
    open(rp, "w").write(f">r\n{ref_seq}\n")
    open(mp, "w").write(f">m\n{mem_seq}\n")
    res = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "10", "-p", "0.1", rp, mp],
                         capture_output=True, text=True)
    fwd, rev = [], []           # (tstart,tend) blocks on ref
    qfwd, qrev = [], []         # (qstart,qend) blocks on member (for member-side coverage)
    fwd_pairs = []              # (qstart, tstart) for collinearity/order
    rlen = len(ref_seq); mlen = len(mem_seq)
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qs, qe, strand, ts, te = int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
        if te - ts < MIN_BLOCK:
            continue
        if strand == "+":
            fwd.append((ts, te)); qfwd.append((qs, qe)); fwd_pairs.append((qs, ts))
        else:
            rev.append((ts, te)); qrev.append((qs, qe))
    # coherence = fraction of the MEMBER that threads the backbone (forward), not of the reference
    coh = sum(e - s for s, e in merge_intervals(qfwd)) / max(mlen, 1)
    inv = sum(e - s for s, e in merge_intervals(qrev)) / max(mlen, 1)
    # tandem-repeat: ref bp covered by >=2 distinct query blocks
    cov = collections.defaultdict(int)
    # count multiplicity by sweeping block starts/ends
    events = []
    for s, e in fwd + rev:
        events.append((s, 1)); events.append((e, -1))
    events.sort()
    depth = 0; last = None; multi_bp = 0
    for pos, d in events:
        if last is not None and depth >= 2:
            multi_bp += pos - last
        depth += d; last = pos
    rep = multi_bp / max(rlen, 1)
    # exon-shuffle: forward blocks out of collinear order (query order vs target order)
    shuffle = 0.0
    if len(fwd_pairs) >= 2:
        order = [t for q, t in sorted(fwd_pairs)]    # target coords sorted by query
        # count inversions in target order (out-of-order fraction)
        disc = sum(1 for i in range(len(order)) for j in range(i + 1, len(order)) if order[i] > order[j])
        tot = len(order) * (len(order) - 1) / 2
        shuffle = disc / max(tot, 1)
    n_blocks = len(merge_intervals(fwd + rev))
    return dict(coh_fwd=round(coh, 3), inv_frac=round(inv, 3), rep_hits=round(rep, 3),
                shuffle=round(shuffle, 3), n_blocks=n_blocks)


def work_family(args):
    fid, members = args
    recs = []
    tmp = tempfile.mkdtemp(prefix=f"vgc_{os.getpid()}_")
    seqs = {m: CACHE[m] for m in members if m in CACHE}
    if len(seqs) < 2:
        return recs
    ref = max(seqs, key=lambda m: len(seqs[m]))
    for m in seqs:
        if m == ref:
            continue
        feat = analyze_member(seqs[ref], seqs[m], tmp)
        recs.append(dict(family=fid, ref=ref, mem=m, **feat))
    for fn in ("r.fa", "m.fa"):
        try: os.remove(os.path.join(tmp, fn))
        except OSError: pass
    try: os.rmdir(tmp)
    except OSError: pass
    return recs


def build_families():
    Hd, _ = dna_homology(); idx = load_intron_index(); bt = gene_biotype()
    coord, strand = load_coord_strand()
    G, _, _ = build_family_graph(Hd, idx, filter_retrocopy=True)
    st = edge_strand()
    G.remove_edges_from([e for e in G.edges() if st.get(tuple(sorted(e))) == "-"])
    pe = protein_whole_edges()
    fams = []
    Ga = nx.Graph()
    for a, b in G.edges():
        if tuple(sorted((a, b))) in pe:
            Ga.add_edge(a, b)
    for i, c in enumerate(nx.connected_components(Ga)):
        if len(c) >= 2:
            fams.append((f"RCF_{i}", sorted(c)))
    # immune
    segset = {g for g, b in bt.items() if b in SEG_BT and g in G}
    Gi = nx.Graph(); Gi.add_nodes_from(segset)
    for a, b in G.subgraph(segset).edges():
        ca, cb = coord.get(a), coord.get(b)
        if ca and cb and ca[0] == cb[0] and abs(ca[1] - cb[1]) <= LOCUS_WIN:
            Gi.add_edge(a, b)
    for i, c in enumerate(nx.connected_components(Gi)):
        if len(c) >= 2:
            fams.append((f"IMM_{i}", sorted(c)))
    return fams


def main():
    global EXON, CACHE
    EXON = json.load(open(EXON_INDEX))
    fams = build_families()
    genes = {g for _, ms in fams for g in ms}
    fa = pysam.FastaFile(GENOME)
    CACHE = {}
    for g in genes:
        s = build_cdna(g, fa)
        if s:
            CACHE[g] = s
    fa.close()
    print(f"[setup] {len(fams)} families, cDNA for {len(CACHE)}/{len(genes)} genes; aligning members ...")
    recs = []
    with Pool(5) as pool:
        for batch in pool.imap_unordered(work_family, fams, chunksize=8):
            recs.extend(batch)
    cols = ["family", "ref", "mem", "coh_fwd", "inv_frac", "rep_hits", "shuffle", "n_blocks"]
    with open(os.path.join(HERE, "family_def_vg_coherence.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in recs:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # per-family aggregate + flags
    byfam = collections.defaultdict(list)
    for r in recs:
        byfam[r["family"]].append(r)
    low_coh = inv = rep = shuf = 0
    flagged = collections.defaultdict(list)
    for fid, rs in byfam.items():
        mc = min(r["coh_fwd"] for r in rs)            # worst-threading member
        if mc < 0.5:
            low_coh += 1; flagged["low_coherence"].append(fid)
        if any(r["inv_frac"] >= 0.10 for r in rs):
            inv += 1; flagged["inversion"].append(fid)
        if any(r["rep_hits"] >= 0.15 for r in rs):
            rep += 1; flagged["tandem_repeat"].append(fid)
        if any(r["shuffle"] >= 0.20 for r in rs):
            shuf += 1; flagged["exon_shuffle"].append(fid)
    n = len(byfam)
    print(f"\n=== VG coherence + non-linear structure over {n} families ({len(recs)} member alignments) ===")
    print(f"  (A) LOW-COHERENCE (worst member threads <50% of backbone) = residual-bridge: {low_coh} ({100*low_coh/n:.0f}%)")
    print(f"  (B) NON-LINEAR cases the collinear pipeline is blind to:")
    print(f"      INVERSION   (>=10% reverse-strand block) : {inv} families ({100*inv/n:.0f}%)")
    print(f"      TANDEM-REPEAT (ref region matched >=2x)   : {rep} families ({100*rep/n:.0f}%)")
    print(f"      EXON-SHUFFLE (>=20% out-of-order blocks)   : {shuf} families ({100*shuf/n:.0f}%)")
    for k in ("low_coherence", "inversion", "tandem_repeat", "exon_shuffle"):
        if flagged[k]:
            print(f"    {k}: {flagged[k][:8]}{' ...' if len(flagged[k])>8 else ''}")
    json.dump(dict(n_families=n, n_alignments=len(recs),
                   low_coherence=low_coh, inversion=inv, tandem_repeat=rep, exon_shuffle=shuf,
                   flagged={k: v for k, v in flagged.items()}),
              open(os.path.join(HERE, "family_def_vg_coherence.json"), "w"), indent=2)
    print(f"\n[+] wrote family_def_vg_coherence.{{tsv,json}}")


if __name__ == "__main__":
    main()
