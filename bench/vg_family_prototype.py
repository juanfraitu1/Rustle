#!/usr/bin/env python3
"""vg_family_prototype.py — VG-native, minimizer-free multi-copy family definition.

Build a pan-transcriptome exon-splice variation graph directly from de-novo
skeletons, merge exons by sequence similarity (vsearch, no minimizers), detect
repeat hubs from VG topology (node/edge multiplicity), and extract families as
cohesive components.  Compare to the shipped RNA catalog.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py
"""
import os
import sys

import argparse
import csv
import json
import re
import subprocess
import tempfile
from collections import defaultdict

import edlib
import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(__file__))
from genome_family_def import refine_families, distinct_loci  # noqa: E402

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"

SKEL_TSV = os.path.join(SCRATCH, "denovo_skeletons.tsv")
META_TSV = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
FASTA = os.path.join(SCRATCH, "GGO.fasta")
VSEARCH = "/home/juanfra/miniforge3/bin/vsearch"
CDHITEST = "/home/juanfra/miniforge3/bin/cd-hit-est"
MMSEQS = "/home/juanfra/miniforge3/bin/mmseqs"

OUT_TSV = os.path.join(BENCH, "vg_family_prototype.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype.json")
O1_EDGE_TSV = os.path.join(BENCH, "denovo_family_edges.tsv")

RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(RC)[::-1]


def canonical(seq):
    """Return (canonical_seq, is_rc) where canonical is lexicographically smaller
    of seq and its reverse-complement."""
    rc = revcomp(seq)
    if seq <= rc:
        return seq, False
    return rc, True


def load_meta():
    """DN id -> {chrom,start,end,strand,n_exon,n_reads}."""
    d = {}
    with open(META_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["id"]] = dict(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                n_exon=int(row["n_exon"]),
                n_reads=int(row["n_reads"]),
            )
    return d


def load_skeletons(meta):
    """Return list of locus dicts {lid, chrom, start, end, strand, introns, n_reads}
    for loci present in meta."""
    rows = []
    with open(SKEL_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            ne = int(row["n_exon"])
            lid = f"DN_{chrom}_{start}_{ne}"
            if lid not in meta:
                continue
            introns = []
            if row["introns"].strip():
                for tok in row["introns"].split(";"):
                    a, b = tok.split("-")
                    introns.append((int(a), int(b)))
            rows.append(dict(
                lid=lid,
                chrom=chrom,
                start=start,
                end=end,
                strand=meta[lid]["strand"],
                introns=introns,
                n_reads=int(row["n_reads"]),
            ))
    return rows


def extract_exon_paths(fa, loci):
    """For each locus, extract ordered exon sequences in transcription orientation.

    Returns:
        exon_occurrences: list of {oid, lid, eidx, seq, canon, is_rc, chrom, start, end}
        locus_paths: dict lid -> list of exon-occurrence ids (ordered 5'->3')
    """
    exons = []
    locus_paths = {}
    oid = 0
    for loc in loci:
        chrom, start, end, strand, introns = loc["chrom"], loc["start"], loc["end"], loc["strand"], loc["introns"]
        # exon genomic intervals
        exon_ivs = []
        prev = start
        for d, a in introns:
            exon_ivs.append((prev, d))
            prev = a
        exon_ivs.append((prev, end))
        # fetch sequences in genomic orientation
        seqs = []
        for s, e in exon_ivs:
            seq = fa.fetch(chrom, s, e).upper()
            seqs.append(seq)
        # orient 5'->3'
        if strand == "-":
            seqs = [revcomp(s) for s in seqs[::-1]]
            exon_ivs = exon_ivs[::-1]
        path_ids = []
        for eidx, (seq, (s, e)) in enumerate(zip(seqs, exon_ivs)):
            canon_seq, is_rc = canonical(seq)
            exons.append(dict(
                oid=oid,
                lid=loc["lid"],
                eidx=eidx,
                seq=seq,
                canon=canon_seq,
                is_rc=is_rc,
                chrom=chrom,
                start=s,
                end=e,
            ))
            path_ids.append(oid)
            oid += 1
        locus_paths[loc["lid"]] = path_ids
    return exons, locus_paths


def exact_dedup(exons):
    """Collapse identical canonical exon sequences.  Returns seq -> first oid,
    and a mapping from each oid to the representative oid."""
    seq_to_rep = {}
    oid_to_rep = {}
    for ex in exons:
        s = ex["canon"]
        if s not in seq_to_rep:
            seq_to_rep[s] = ex["oid"]
        oid_to_rep[ex["oid"]] = seq_to_rep[s]
    return seq_to_rep, oid_to_rep


def cluster_exons_exact(exons):
    """Exact-sequence graph-to-graph mode: every distinct canonical exon sequence is its own VG node.

    Returns (oid_to_rep, rep_to_cluster) compatible with build_graph.
    """
    seq_to_rep, oid_to_rep = exact_dedup(exons)
    reps = {oid: seq for seq, oid in seq_to_rep.items()}
    # node id = representative oid; unique per canonical sequence.
    rep_to_cluster = {oid: oid for oid in reps}
    print(f"    exact-dedup representatives: {len(reps)}", flush=True)
    return oid_to_rep, rep_to_cluster


_BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _kmer_iter(seq, k):
    """Yield integer encodings of canonical kmers in seq."""
    mask = (1 << (2 * k)) - 1
    kmer = 0
    valid = 0
    for ch in seq:
        b = _BASE.get(ch)
        if b is None:
            valid = 0
            kmer = 0
        else:
            kmer = ((kmer << 2) | b) & mask
            valid += 1
            if valid >= k:
                yield kmer


def _unique_kmers(seq, k):
    """Set of integer kmers in seq."""
    seen = set()
    out = []
    for kmer in _kmer_iter(seq, k):
        if kmer not in seen:
            seen.add(kmer)
            out.append(kmer)
    return out


def cluster_exons_seeded(exons, k=15, identity=0.90, max_occ=50, min_seed=2):
    """Seed-and-extend approximate exon clustering.

    - Index all non-over-represented canonical kmers (<= max_occ occurrences).
    - For each representative exon, find candidates that share at least one kmer.
    - Verify candidate pairs with edlib global alignment; keep if identity >= threshold.
    - Return connected components as exon clusters (nodes).

    This avoids global all-vs-all clustering and is fully VG-native in the sense
    that near-identical exons are merged into shared nodes via local alignment.
    """
    seq_to_rep, oid_to_rep = exact_dedup(exons)
    reps = {oid: seq for seq, oid in seq_to_rep.items()}
    print(f"    exact-dedup representatives: {len(reps)}", flush=True)
    print(f"    indexing kmers (k={k}, max_occ={max_occ}) ...", flush=True)

    # First pass: count kmer occurrences (one per sequence).
    counts = defaultdict(int)
    for seq in reps.values():
        for kmer in _unique_kmers(seq, k):
            counts[kmer] += 1

    # Second pass: build sorted kmer -> oid index, dropping repeat kmers.
    kmers_list = []
    oids_list = []
    for oid, seq in reps.items():
        for kmer in _unique_kmers(seq, k):
            if counts[kmer] <= max_occ:
                kmers_list.append(kmer)
                oids_list.append(oid)
    kmers_arr = np.array(kmers_list, dtype=np.uint64)
    oids_arr = np.array(oids_list, dtype=np.uint32)
    del kmers_list, oids_list, counts
    order = np.argsort(kmers_arr, kind="mergesort")
    kmers_arr = kmers_arr[order]
    oids_arr = oids_arr[order]
    # unique kmers + start positions for O(1) candidate lookup
    uniq, inv, ucounts = np.unique(kmers_arr, return_inverse=True, return_counts=True)
    starts = np.cumsum(ucounts) - ucounts
    print(f"    indexed kmers: {len(uniq)} unique", flush=True)

    # Union-find over representative oids.
    parent = {oid: oid for oid in reps}
    size = {oid: 1 for oid in reps}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]

    print(f"    verifying candidate pairs with edlib (identity>={identity}) ...", flush=True)
    min_len_ratio = identity  # identity floor; edlib distance <= (1-id)*maxlen
    sorted_oids = sorted(reps)
    for qi, oid in enumerate(sorted_oids):
        if qi % 10000 == 0 and qi:
            print(f"      processed {qi}/{len(sorted_oids)} reps", flush=True)
        seq = reps[oid]
        qkmers = _unique_kmers(seq, k)
        if not qkmers:
            continue
        qkmers_a = np.array(qkmers, dtype=np.uint64)
        idx = np.searchsorted(uniq, qkmers_a, side="left")
        valid_mask = idx < len(uniq)
        if valid_mask.any():
            vidx = idx[valid_mask]
            valid_mask[valid_mask] = (uniq[vidx] == qkmers_a[valid_mask])
        # Collect candidate oids from all seed hits and count occurrences in numpy.
        hits = []
        for i in np.where(valid_mask)[0]:
            s = starts[idx[i]]
            e = s + ucounts[idx[i]]
            if e > s:
                hits.append(oids_arr[s:e])
        if not hits:
            continue
        cand_arr = np.concatenate(hits).astype(np.int64)
        if cand_arr.size == 0:
            continue
        uniq_cands, seed_counts = np.unique(cand_arr, return_counts=True)
        keep = (seed_counts >= min_seed) & (uniq_cands > oid)
        for cand in uniq_cands[keep]:
            tseq = reps[cand]
            maxlen = max(len(seq), len(tseq))
            # quick length filter
            if min(len(seq), len(tseq)) / maxlen < identity:
                continue
            aln = edlib.align(seq, tseq, mode="NW", task="distance")
            dist = aln["editDistance"]
            if dist == -1:
                continue
            if 1 - dist / maxlen >= identity:
                union(int(cand), oid)

    rep_to_cluster = {oid: find(oid) for oid in reps}
    n_nodes = len(set(rep_to_cluster.values()))
    print(f"    seeded clusters (nodes): {n_nodes}", flush=True)
    return oid_to_rep, rep_to_cluster


def cluster_exons_vsearch(exons, identity=0.95, threads=4):
    """Cluster representative exon sequences with vsearch (--cluster_fast).

    Returns (oid_to_rep, rep_to_cluster).
    """
    seq_to_rep, oid_to_rep = exact_dedup(exons)
    reps = {oid: seq for seq, oid in seq_to_rep.items()}
    print(f"    exact-dedup representatives: {len(reps)}", flush=True)
    with tempfile.TemporaryDirectory(prefix="vgfp_vs_") as td:
        in_fa = os.path.join(td, "exons.fa")
        out_uc = os.path.join(td, "clusters.uc")
        # vsearch --cluster_fast requires input sorted by length descending.
        sorted_reps = sorted(reps.items(), key=lambda kv: (-len(kv[1]), kv[0]))
        with open(in_fa, "w") as fh:
            for oid, seq in sorted_reps:
                fh.write(f">{oid}\n{seq}\n")
        cmd = [
            VSEARCH,
            "--cluster_fast", in_fa,
            "--id", str(identity),
            "--strand", "both",
            "--minseqlength", "1",
            "--uc", out_uc,
            "--threads", str(threads),
            "--quiet",
        ]
        subprocess.run(cmd, check=True)
        rep_to_cluster = {}
        with open(out_uc) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                typ = f[0]
                cid = int(f[1])
                oid = int(f[8].split(";")[0])  # target/query id field
                if typ == "S":
                    rep_to_cluster[oid] = cid
                elif typ == "H":
                    rep_to_cluster[oid] = cid
                elif typ == "C":
                    pass
        # vsearch should keep everything with --minseqlength 1; guard against future changes.
        next_cid = max(rep_to_cluster.values(), default=-1) + 1
        for oid in reps:
            if oid not in rep_to_cluster:
                rep_to_cluster[oid] = next_cid
                next_cid += 1
        return oid_to_rep, rep_to_cluster


def cluster_exons_cdhit(exons, identity=0.95, threads=4):
    """Cluster representative exon sequences with cd-hit-est.

    Returns (oid_to_rep, rep_to_cluster).
    """
    seq_to_rep, oid_to_rep = exact_dedup(exons)
    reps = {oid: seq for seq, oid in seq_to_rep.items()}
    print(f"    exact-dedup representatives: {len(reps)}", flush=True)
    with tempfile.TemporaryDirectory(prefix="vgfp_cls_") as td:
        in_fa = os.path.join(td, "exons.fa")
        out_pref = os.path.join(td, "clusters")
        clstr = out_pref + ".clstr"
        with open(in_fa, "w") as fh:
            for oid in sorted(reps):
                fh.write(f">{oid}\n{reps[oid]}\n")
        # choose word size per cd-hit-est recommendations
        if identity >= 0.95:
            n = 8
        elif identity >= 0.90:
            n = 6
        else:
            n = 5
        cmd = [
            CDHITEST,
            "-i", in_fa,
            "-o", out_pref,
            "-c", str(identity),
            "-n", str(n),
            "-d", "0",
            "-T", str(threads),
            "-r", "1",          # both strands
            "-M", "0",          # unlimited memory
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # parse .clstr
        oid_re = re.compile(r">(\d+)\.\.\.")
        rep_to_cluster = {}
        cur_cid = None
        with open(clstr) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">Cluster"):
                    cur_cid = int(line.split()[1])
                elif line:
                    m = oid_re.search(line)
                    if m:
                        oid = int(m.group(1))
                        rep_to_cluster[oid] = cur_cid
        # cd-hit-est can drop very short / low-complexity reps; keep them as singleton clusters.
        next_cid = max(rep_to_cluster.values(), default=-1) + 1
        for oid in reps:
            if oid not in rep_to_cluster:
                rep_to_cluster[oid] = next_cid
                next_cid += 1
        return oid_to_rep, rep_to_cluster


def cluster_exons_mmseqs(exons, identity=0.90, threads=4):
    """Cluster representative exon sequences with mmseqs2 easy-cluster.

    This is a fast k-mer-seeded align-and-extend implementation of approximate
    exon clustering (the same biological operation as the pure-Python seeded
    mode, but executed by mmseqs2's C++ engine).

    Returns (oid_to_rep, rep_to_cluster).
    """
    seq_to_rep, oid_to_rep = exact_dedup(exons)
    reps = {oid: seq for seq, oid in seq_to_rep.items()}
    print(f"    exact-dedup representatives: {len(reps)}", flush=True)
    with tempfile.TemporaryDirectory(prefix="vgfp_mmseqs_") as td:
        in_fa = os.path.join(td, "exons.fa")
        out_pref = os.path.join(td, "out")
        tmp = os.path.join(td, "tmp")
        with open(in_fa, "w") as fh:
            for oid in sorted(reps):
                fh.write(f">{oid}\n{reps[oid]}\n")
        cmd = [
            MMSEQS,
            "easy-cluster",
            in_fa,
            out_pref,
            tmp,
            "--min-seq-id", str(identity),
            "-c", "0.80",
            "--cov-mode", "1",
            "--threads", str(threads),
            "-v", "0",
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cluster_tsv = out_pref + "_cluster.tsv"
        rep_to_cluster = {}
        with open(cluster_tsv) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                rep_oid = int(f[0])
                mem_oid = int(f[1])
                rep_to_cluster[mem_oid] = rep_oid
        # any missing reps become singleton clusters
        next_cid = max(rep_to_cluster.values(), default=-1) + 1
        for oid in reps:
            if oid not in rep_to_cluster:
                rep_to_cluster[oid] = next_cid
                next_cid += 1
        return oid_to_rep, rep_to_cluster


def build_graph(exons, locus_paths, oid_to_rep, rep_to_cluster):
    """Build exon-splice VG.  Returns:
        node_seq: dict node_id -> canonical sequence
        node_loci: dict node_id -> set of locus ids using the node
        edge_loci: dict (node_u, node_v) -> set of locus ids using the edge
        locus_path_nodes: dict lid -> ordered list of node ids
    """
    oid_to_node = {}
    node_seq = {}
    for ex in exons:
        rep = oid_to_rep[ex["oid"]]
        node = rep_to_cluster[rep]
        oid_to_node[ex["oid"]] = node
        node_seq[node] = ex["canon"]

    node_loci = defaultdict(set)
    edge_loci = defaultdict(set)
    locus_path_nodes = {}
    for lid, oids in locus_paths.items():
        nodes = [oid_to_node[oid] for oid in oids]
        locus_path_nodes[lid] = nodes
        for n in nodes:
            node_loci[n].add(lid)
        for i in range(len(nodes) - 1):
            edge_loci[(nodes[i], nodes[i + 1])].add(lid)
    return dict(node_seq), dict(node_loci), dict(edge_loci), locus_path_nodes


def apply_repeat_hub_gate(node_loci, edge_loci, locus_path_nodes, threshold):
    """Remove edges incident to a node whose path multiplicity >= threshold.

    Returns:
        kept_edges: set of (u,v) edges surviving the gate
        bad_nodes: set of hub nodes
    """
    bad_nodes = {n for n, loci in node_loci.items() if len(loci) >= threshold}
    kept_edges = set()
    for (u, v), loci in edge_loci.items():
        if u not in bad_nodes and v not in bad_nodes:
            kept_edges.add((u, v))
    return kept_edges, bad_nodes


def locus_graph_from_edges(locus_path_nodes, kept_edges):
    """Build locus-locus graph: edge if two loci share a surviving edge in the VG."""
    edge_to_loci = defaultdict(set)
    for lid, nodes in locus_path_nodes.items():
        for i in range(len(nodes) - 1):
            e = (nodes[i], nodes[i + 1])
            if e in kept_edges:
                edge_to_loci[e].add(lid)
    locus_adj = defaultdict(set)
    for loci in edge_to_loci.values():
        for a in loci:
            for b in loci:
                if a != b:
                    locus_adj[a].add(b)
    return locus_adj


def locus_graph_from_vg(node_loci, edge_loci, kept_edges, bad_nodes):
    """Graph-to-graph linkage: connect two loci if they share a non-repeat VG node (exon)
    OR a surviving splice edge.  This is the direct graph-of-graphs family edge."""
    locus_adj = defaultdict(set)
    # shared exon nodes
    for node, loci in node_loci.items():
        if node in bad_nodes:
            continue
        loci = list(loci)
        for i in range(len(loci)):
            for j in range(i + 1, len(loci)):
                a, b = loci[i], loci[j]
                locus_adj[a].add(b)
                locus_adj[b].add(a)
    # shared surviving splice edges
    edge_to_loci = defaultdict(set)
    for e, loci in edge_loci.items():
        if e in kept_edges:
            for lid in loci:
                edge_to_loci[e].add(lid)
    for loci in edge_to_loci.values():
        loci = list(loci)
        for i in range(len(loci)):
            for j in range(i + 1, len(loci)):
                a, b = loci[i], loci[j]
                locus_adj[a].add(b)
                locus_adj[b].add(a)
    return locus_adj


def load_o1_edges(threshold=0.13):
    """Load the O1 transcript-homology edges (core_recip >= threshold)."""
    edges = []
    with open(O1_EDGE_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if float(row["core_recip"]) >= threshold:
                edges.append((row["a"], row["b"]))
    return edges


def locus_graph_o1vg(edge_pairs, locus_path_nodes, node_loci, bad_nodes, kept_edges):
    """O1+VG integration: keep an O1 homology edge only if the two loci share a
    non-repeat VG exon node OR a surviving splice edge.  This replaces the
    minimizer-based repeat gate with a VG-topology support check."""
    good_nodes = set(node_loci) - bad_nodes
    kept_edges_set = set(kept_edges)
    locus_adj = defaultdict(set)
    for a, b in edge_pairs:
        if a not in locus_path_nodes or b not in locus_path_nodes:
            continue
        path_a = locus_path_nodes[a]
        path_b = locus_path_nodes[b]
        nodes_a = set(path_a)
        nodes_b = set(path_b)
        # shared unique exon
        shared_good = (nodes_a & nodes_b) & good_nodes
        if shared_good:
            locus_adj[a].add(b)
            locus_adj[b].add(a)
            continue
        # shared surviving splice edge
        edges_a = set(zip(path_a, path_a[1:]))
        edges_b = set(zip(path_b, path_b[1:]))
        if (edges_a & edges_b) & kept_edges_set:
            locus_adj[a].add(b)
            locus_adj[b].add(a)
    return locus_adj


def extract_components(locus_adj, loci_list):
    """Connected components of locus graph."""
    seen = set()
    comps = []
    for lid in loci_list:
        if lid in seen:
            continue
        stack = [lid]
        comp = []
        seen.add(lid)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in locus_adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        if len(comp) >= 2:
            comps.append(sorted(comp))
    return comps


def build_genes_array(loci_list, locus_info):
    """Build the `genes` array expected by genome_family_def.refine_families."""
    genes = []
    name_to_idx = {}
    for lid in loci_list:
        info = locus_info[lid]
        genes.append(dict(
            chrom=info["chrom"],
            start=info["start"],
            end=info["end"],
            name=lid,
            biotype="NA",
        ))
        name_to_idx[lid] = len(genes) - 1
    return genes, name_to_idx


def main():
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["exact", "seeded", "mmseqs", "vsearch", "cdhit", "o1vg"], default="o1vg",
                        help="family-linkage mode: O1+VG (default), kmer-seeded edlib/mmseqs, vsearch/cdhit, or exact shared-exon graph")
    parser.add_argument("--identity", type=float, default=0.95,
                        help="cd-hit-est/vsearch/seeded exon clustering identity (used with --mode vsearch/cdhit/seeded; default 0.95)")
    parser.add_argument("--k", type=int, default=15,
                        help="kmer length for --mode seeded (default 15)")
    parser.add_argument("--max-occ", type=int, default=50,
                        help="max kmer occurrences to use as seed in --mode seeded (default 50)")
    parser.add_argument("--min-seed", type=int, default=2,
                        help="minimum shared kmers before edlib verification in --mode seeded (default 2)")
    parser.add_argument("--o1-edge-thresh", type=float, default=0.13,
                        help="O1 edge creation core_recip threshold (used only with --mode o1vg; default 0.13)")
    parser.add_argument("--repeat-thresh", type=int, default=30,
                        help="node multiplicity threshold for repeat-hub gate (default 30)")
    parser.add_argument("--gamma", type=float, default=0.20,
                        help="gamma-quasi-clique cohesion threshold (default 0.20)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    print("[*] loading meta + skeletons ...", flush=True)
    meta = load_meta()
    loci = load_skeletons(meta)
    print(f"    loci from meta: {len(loci)}", flush=True)

    print("[*] extracting exon paths ...", flush=True)
    fa = pysam.FastaFile(FASTA)
    exons, locus_paths = extract_exon_paths(fa, loci)
    print(f"    exon occurrences: {len(exons)}", flush=True)

    if args.mode == "exact":
        print("[*] building exact pan-exon graph (graph-to-graph, no external clustering) ...", flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_exact(exons)
    elif args.mode == "seeded":
        print("[*] clustering exons with kmer-seeded edlib (k={}, identity={}) ...".format(args.k, args.identity), flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_seeded(
            exons, k=args.k, identity=args.identity, max_occ=args.max_occ, min_seed=args.min_seed)
    elif args.mode == "mmseqs":
        print("[*] clustering exons with mmseqs2 (identity={}) ...".format(args.identity), flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_mmseqs(exons, identity=args.identity, threads=args.threads)
    elif args.mode == "vsearch":
        print("[*] clustering exons with vsearch (identity={}) ...".format(args.identity), flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_vsearch(exons, identity=args.identity, threads=args.threads)
    elif args.mode == "cdhit":
        print("[*] clustering exons with cd-hit-est (identity={}) ...".format(args.identity), flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_cdhit(exons, identity=args.identity, threads=args.threads)
    else:
        # o1vg uses near-exact exon clustering to decide VG support for O1 edges
        print("[*] clustering exons with cd-hit-est for O1+VG support (identity={}) ...".format(args.identity), flush=True)
        oid_to_rep, rep_to_cluster = cluster_exons_cdhit(exons, identity=args.identity, threads=args.threads)
    n_nodes = len(set(rep_to_cluster.values()))
    print(f"    representative sequences: {len(set(oid_to_rep.values()))}", flush=True)
    print(f"    exon clusters (nodes): {n_nodes}", flush=True)

    print("[*] building VG ...", flush=True)
    node_seq, node_loci, edge_loci, locus_path_nodes = build_graph(
        exons, locus_paths, oid_to_rep, rep_to_cluster)
    print(f"    intron edges: {len(edge_loci)}", flush=True)

    print("[*] repeat-hub gate (threshold={}) ...".format(args.repeat_thresh), flush=True)
    kept_edges, bad_nodes = apply_repeat_hub_gate(
        node_loci, edge_loci, locus_path_nodes, args.repeat_thresh)
    print(f"    hub nodes: {len(bad_nodes)} / {len(node_loci)}", flush=True)
    print(f"    edges after gate: {len(kept_edges)} / {len(edge_loci)}", flush=True)

    print("[*] building locus graph ...", flush=True)
    if args.mode == "exact":
        locus_adj = locus_graph_from_vg(node_loci, edge_loci, kept_edges, bad_nodes)
    elif args.mode == "o1vg":
        print("[*] loading O1 homology edges ...", flush=True)
        o1_edges = load_o1_edges(args.o1_edge_thresh)
        print(f"    O1 edges (core_recip >= {args.o1_edge_thresh}): {len(o1_edges)}", flush=True)
        locus_adj = locus_graph_o1vg(o1_edges, locus_path_nodes, node_loci, bad_nodes, kept_edges)
    else:
        locus_adj = locus_graph_from_edges(locus_path_nodes, kept_edges)
    loci_list = sorted(locus_paths.keys())
    raw_comps = extract_components(locus_adj, loci_list)
    print(f"    raw components (>=2 loci): {len(raw_comps)}", flush=True)

    print("[*] refining families (gamma={}) ...".format(args.gamma), flush=True)
    genes, name_to_idx = build_genes_array(loci_list, meta)
    idx_comps = [[name_to_idx[lid] for lid in comp] for comp in raw_comps]
    all_edges = set()
    for lid, nbs in locus_adj.items():
        ui = name_to_idx[lid]
        for other in nbs:
            if other > lid:
                vi = name_to_idx[other]
                all_edges.add((ui, vi))
    refined = refine_families(idx_comps, all_edges, genes, args.gamma, seed=0)
    print(f"    refined families (>=2 loci): {len(refined)}", flush=True)

    # write catalog
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["fam_id", "member"])
        for fid, members in enumerate(refined):
            for idx in members:
                w.writerow([fid, genes[idx]["name"]])
    print(f"[+] wrote {OUT_TSV}", flush=True)

    # summary JSON
    metrics = dict(
        n_loci=len(loci),
        n_exon_occurrences=len(exons),
        n_nodes=n_nodes,
        n_edges=len(edge_loci),
        n_hub_nodes=len(bad_nodes),
        n_edges_after_gate=len(kept_edges),
        n_raw_components=len(raw_comps),
        n_refined_families=len(refined),
        params=dict(mode=args.mode, identity=args.identity, o1_edge_thresh=args.o1_edge_thresh,
                    repeat_thresh=args.repeat_thresh, gamma=args.gamma),
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(metrics, fh, indent=2, sort_keys=True)
    print(f"[+] wrote {OUT_JSON}", flush=True)
    print("\n=== SUMMARY ===")
    for k, v in metrics.items():
        if k != "params":
            print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
