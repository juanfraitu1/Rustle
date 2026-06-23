#!/usr/bin/env python3
"""
family_def_vg_reinforce.py — clean genome-wide multi-copy family definition:
de-novo loci (tight vertices) + per-family VG BACKBONE-REINFORCEMENT subgraph + guard.

  DISCOVER : read-conflict (de-tie) graph over de-novo loci -> candidate families.
  REINFORCE: for each family, build every member's all-isoform exon-union copy model,
             align them, and form the FAMILY BACKBONE GRAPH (edge = two members share a
             real backbone). Real copies mutually reinforce (a member backbone-linked to
             several mates is strongly real); a repeat-bridge member has NO backbone link
             and is rejected. Size-2 families survive (the pair reinforces itself iff it
             shares a backbone). A guard never rejects a member too sparse to validate.
  RESULT   : validated families = backbone-connected cores; bridge members rejected;
             over-merged families split into their real backbone cores.

OOM-safe: models built ONLY for family-member loci (~1.3k); reads capped; streamed.
Run (background): /home/juanfra/miniforge3/bin/python bench/family_def_vg_reinforce.py
"""
import collections
import json
import os
import subprocess
import sys

import pysam
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import (scan as conflict_scan, pair_evidence, components,
                                    best_gene, DELTA, DE_MAX, MIN_READS)

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
PAF = "/home/juanfra/winloci_scratch/vg_reinforce_copies.paf"
PAD = 5000
MAX_READS = 3000
GUARD_READS = 20      # below this, a model is "unvalidatable" -> never rejected (guard)
MIN_LEN = 300
ID_MIN = 0.80
RECIP_COV_MIN = 0.30  # backbone edge: id>=.80 AND BOTH copies aligned >= this


def load_denovo():
    by_chrom = collections.defaultdict(list)
    info = {}
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            vid, c, s, e = p[0], p[1], int(p[2]), int(p[3])
            by_chrom[c].append((s, e, vid))
            info[vid] = (c, s, e, int(p[6]))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, info


def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build_model(bam, genome, chrom, start, end):
    blocks = []; n = 0
    try:
        it = bam.fetch(chrom, max(0, start - PAD), end + PAD)
    except (ValueError, KeyError):
        return "", 0
    for r in it:
        if r.is_unmapped or r.is_supplementary:
            continue
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        b = r.get_blocks()
        if b:
            blocks.extend(b); n += 1
            if n >= MAX_READS:
                break
    if not blocks:
        return "", n
    union = merge([list(x) for x in blocks])
    return "".join(genome.fetch(chrom, s, e).upper() for s, e in union), n


def main():
    print("[1/4 discover] read-conflict over de-novo loci ...", flush=True)
    by_chrom, info = load_denovo()
    mm, _, _ = conflict_scan(by_chrom)
    ev = pair_evidence(mm)
    edges, fams = components(ev, DELTA, DE_MAX, MIN_READS)
    member_loci = sorted({g for c in fams for g in c})
    def chrom_of(v): return info[v][0] if v in info else None
    n_bridge0 = sum(1 for c in fams if len({chrom_of(g) for g in c}) >= 3)
    print(f"  candidate: {len(fams)} families, {len(edges)} edges, {len(member_loci)} loci; "
          f"cross-chrom bridges={n_bridge0}", flush=True)

    print(f"[2/4 build] exon-union copy models for {len(member_loci)} member loci ...", flush=True)
    bam = pysam.AlignmentFile(BAM, "rb"); genome = pysam.FastaFile(GENOME)
    nreads = {}; built = 0
    with open(SEQS, "w") as out:
        for i, v in enumerate(member_loci):
            c, s, e, _ = info[v]
            seq, nr = build_model(bam, genome, c, s, e)
            nreads[v] = nr
            if len(seq) >= MIN_LEN:
                out.write(f">{v}\n{seq}\n"); built += 1
            if (i + 1) % 300 == 0:
                print(f"    ... {i+1}/{len(member_loci)} ({built} models)", flush=True)
    bam.close(); genome.close()
    print(f"  built {built} models", flush=True)

    print("[3/4 reinforce] backbone alignment + per-family reinforcement subgraph ...", flush=True)
    with open(PAF, "w") as pf:
        subprocess.run(["minimap2", "-cx", "asm20", "-N", "30", "-p", "0.05", "-t", "6",
                        SEQS, SEQS], stdout=pf, stderr=subprocess.DEVNULL)
    H = {}
    with open(PAF) as f:
        for line in f:
            x = line.split("\t")
            if len(x) < 11:
                continue
            q, ql, qs, qe, t, m, bl = x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[9]), int(x[10])
            if q == t:
                continue
            ident = m / bl if bl else 0; cov = (qe - qs) / ql if ql else 0
            key = (q, t) if q < t else (t, q)
            r = H.setdefault(key, {"id": 0.0, "ca": 0.0, "cb": 0.0})
            r["id"] = max(r["id"], ident)
            if q == key[0]:
                r["ca"] = max(r["ca"], cov)
            else:
                r["cb"] = max(r["cb"], cov)

    def backbone(a, b):
        r = H.get((a, b) if a < b else (b, a))
        return bool(r and r["id"] >= ID_MIN and min(r["ca"], r["cb"]) >= RECIP_COV_MIN)

    validated = []          # backbone-reinforced families (>=2)
    rejected_members = []   # validatable loci with NO backbone link to family-mates
    for fam in fams:
        members = sorted(fam)
        BG = nx.Graph(); BG.add_nodes_from(members)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                if backbone(members[i], members[j]):
                    BG.add_edge(members[i], members[j])
        for comp in nx.connected_components(BG):
            if len(comp) >= 2:
                validated.append(set(comp))
        # isolated (backbone-degree 0) AND validatable -> confident bridge reject
        for v in members:
            if BG.degree[v] == 0 and nreads.get(v, 0) >= GUARD_READS:
                rejected_members.append(v)

    n_bridge1 = sum(1 for c in validated if len({chrom_of(g) for g in c}) >= 3)
    print("\n[4/4 result]")
    print(f"  candidate families   : {len(fams)}   ({n_bridge0} cross-chrom bridges)")
    print(f"  VALIDATED families   : {len(validated)} (backbone-reinforced cores)   "
          f"({n_bridge1} cross-chrom bridges)")
    print(f"  rejected bridge members (validatable, 0 backbone links): {len(rejected_members)}")
    sizes0 = collections.Counter(len(c) for c in fams)
    sizes1 = collections.Counter(len(c) for c in validated)
    print(f"  size-2 families: {sizes0[2]} -> {sizes1[2]} (kept iff the pair shares a backbone)")

    # DNA-truth audit: map each rejected/kept ORIGINAL edge to genes, classify
    gby = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            gby[c].append((int(s), int(e), name))
    for c in gby:
        gby[c].sort()
    def gene_of(v):
        c, s, e, _ = info[v]
        return best_gene(gby, c, s, e)
    try:
        from family_def_read_filters import dna_homology
        Hd, dna = dna_homology()
        def cls_edge(a, b):
            ga, gb = gene_of(a), gene_of(b)
            if not ga or not gb or ga == gb:
                return "no_gene"
            r = Hd.get((ga, gb) if ga < gb else (gb, ga))
            if r is None or r["id"] == 0:
                return "no_homology"
            if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
                return "DNA_paralog"
            return "sub_bar" if r["id"] >= 0.80 else "weak"
        rej_set = set(rejected_members)
        rej_edges = [(a, b) for a, b, n in edges if a in rej_set or b in rej_set]
        rc = collections.Counter(cls_edge(a, b) for a, b in rej_edges)
        print("\n  rejected-member edges by DNA class (want mostly no_homology):")
        for k, v in rc.most_common():
            print(f"     {v:5d}  {k}")
    except Exception as ex:
        print(f"  (DNA audit skipped: {ex})")

    json.dump(dict(candidate_families=len(fams), validated_families=len(validated),
                   rejected_members=len(rejected_members),
                   bridges_before=n_bridge0, bridges_after=n_bridge1,
                   size2_before=sizes0[2], size2_after=sizes1[2]),
              open(os.path.join(HERE, "family_def_vg_reinforce.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_vg_reinforce.json')}")


if __name__ == "__main__":
    main()
