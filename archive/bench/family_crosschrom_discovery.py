#!/usr/bin/env python3
"""Genome-wide cross-chromosome copy DISCOVERY at the RNA level (POA-only definition).

The production family grouping gathers copies per genomic region (position-overlap bundles), so
DISPERSED paralogs whose copies sit on DIFFERENT chromosomes (RABL2A vs RABL2B; the headroom
probe's 17 DISPERSED families) are never co-considered and never found. This harness removes that
restriction: it compares all 22,983 genome-wide gene representatives WITHOUT any chromosome
constraint, using the validated RNA-level family definition:

  minimizer-LSH prefilter (cheap, global)  ->  POA contiguous-core coverage >= T  (the definition)

A pair is a multi-copy transcript family iff their longest single ungapped co-aligned block covers
>= T of the shorter gene (bench/poa_family_definition.py, T=0.13: accepts recent tandem dups, rejects
domain-sharers). Minimizer-Jaccard is ONLY a candidate prefilter (loose, so it cannot miss real dups
like RFPL at Jaccard ~0.13); the POA contiguous-core is the actual gate.

Run with /home/juanfra/miniforge3/bin/python (needs Bio.Align). Deterministic (no RNG).

Stages (--stage): candidates (LSH+Jaccard -> /tmp/cc_candidates.tsv) | poa (confirm) | all.
"""
import argparse
import csv
import multiprocessing as mp
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P  # validated minimizers / make_aligner / poa_pair_stats

# worker globals (set by _winit) so large sequences are loaded ONCE per process
_W_SEQS = None
_W_ALN = None


def _winit(fa_path):
    global _W_SEQS, _W_ALN
    _W_SEQS = P.load_fasta(fa_path)
    _W_ALN = P.make_aligner()


def _wpoa(task):
    a, b, jac = task
    sa = _W_SEQS.get(a)
    sb = _W_SEQS.get(b)
    if not sa or not sb:
        return None
    st = P.poa_pair_stats(_W_ALN, sa, sb)
    if st["core_recip"] >= T_CORE:
        return (a, b, round(st["core_recip"], 3), jac)
    return None

FA = "/tmp/gene_reps_gw.fa"
META = "/tmp/gene_reps_gw.meta.tsv"
UNIVERSE = os.path.join(os.path.dirname(__file__), "copy_recovery_eval/results/universe.tsv")
CAND_TSV = "/tmp/cc_candidates.tsv"

# thresholds
REPETITIVE_CAP = 200    # skip minimizers shared by > this many genes (low-complexity/repeat)
MIN_SHARED = 4          # candidate iff >= this many shared (non-repetitive) minimizers
J_PRE = 0.03            # loose Jaccard prefilter (RFPL real dups ~0.13, so safe)
T_CORE = 0.13           # POA contiguous-core coverage gate (the validated definition)
LEN_CAP = 20000         # skip POA when min(len) > this (O(L^2)); reported as size-capped


def load_meta():
    g2chrom = {}
    g2len = {}
    with open(META) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            g2chrom[c[0]] = c[1]
            g2len[c[0]] = int(c[5])
    return g2chrom, g2len


def load_universe_families():
    g2f = {}
    with open(UNIVERSE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g2f.setdefault(r["gene_id"], r["family_id"])
    return g2f


def stage_candidates():
    seqs = P.load_fasta(FA)
    genes = sorted(seqs)
    gi = {g: i for i, g in enumerate(genes)}
    # minimizer sets (validated canonical/blake2b k=15 w=10)
    mins = [P.minimizers(seqs[g].encode()) for g in genes]
    # inverted index minimizer -> gene indices, skipping repetitive minimizers
    inv = defaultdict(list)
    for i, ms in enumerate(mins):
        for m in ms:
            inv[m].append(i)
    shared = defaultdict(int)   # (i,j) -> shared non-repetitive minimizer count
    for m, lst in inv.items():
        if len(lst) > REPETITIVE_CAP or len(lst) < 2:
            continue
        for a in range(len(lst)):
            for b in range(a + 1, len(lst)):
                shared[(lst[a], lst[b])] += 1
    rows = []
    for (i, j), s in shared.items():
        if s < MIN_SHARED:
            continue
        inter = len(mins[i] & mins[j])
        uni = len(mins[i] | mins[j])
        jac = inter / uni if uni else 0.0
        if jac >= J_PRE:
            rows.append((genes[i], genes[j], s, round(jac, 4)))
    rows.sort(key=lambda r: -r[2])
    with open(CAND_TSV, "w") as fh:
        fh.write("gene_a\tgene_b\tshared_min\tjaccard\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
    print(f"[candidates] {len(genes)} genes -> {len(rows)} candidate pairs "
          f"(>= {MIN_SHARED} shared minimizers, jaccard >= {J_PRE}) -> {CAND_TSV}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["candidates", "poa", "all"])
    ap.add_argument("--jobs", type=int, default=5)
    args = ap.parse_args()

    if args.stage in ("candidates", "all"):
        stage_candidates()
    if args.stage == "candidates":
        return

    # POA confirmation (parallel; deterministic — results sorted after)
    g2chrom, g2len = load_meta()
    g2fam = load_universe_families()
    cands = []
    with open(CAND_TSV) as fh:
        fh.readline()
        for line in fh:
            a, b, s, j = line.rstrip("\n").split("\t")
            cands.append((a, b, int(s), float(j)))

    tasks = []
    capped = 0
    for a, b, s, jac in cands:
        if min(g2len.get(a, 0), g2len.get(b, 0)) > LEN_CAP:
            capped += 1
            continue
        tasks.append((a, b, jac))

    jobs = args.jobs
    with mp.Pool(jobs, initializer=_winit, initargs=(FA,)) as pool:
        raw = pool.map(_wpoa, tasks, chunksize=64)
    confirmed = [(a, b, c, j, g2chrom.get(a) != g2chrom.get(b))
                 for r in raw if r for (a, b, c, j) in [r]]
    confirmed.sort()

    # families = connected components over confirmed pairs
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)
    for a, b, *_ in confirmed:
        union(a, b)
    comps = defaultdict(set)
    for a, b, *_ in confirmed:
        comps[find(a)].add(a); comps[find(a)].add(b)
    families = [sorted(m) for m in comps.values()]
    cross_fams = [m for m in families if len({g2chrom.get(g) for g in m}) > 1]

    cross_pairs = [(a, b, c, j) for a, b, c, j, x in confirmed if x]
    print(f"[poa] candidates={len(cands)} size-capped(min_len>{LEN_CAP})={capped} "
          f"confirmed_pairs={len(confirmed)} (cross-chrom={len(cross_pairs)})")
    print(f"[families] {len(families)} families; CROSS-CHROMOSOME families={len(cross_fams)}")

    # validation vs universe DISPERSED (cross-chrom) families
    # a universe cross-chrom family is "recovered" if >=2 of its genes land in one discovered family
    uni_fam_genes = defaultdict(set)
    for g, f in g2fam.items():
        uni_fam_genes[f].add(g)
    uni_cross = {f: gs for f, gs in uni_fam_genes.items()
                 if len(gs) >= 2 and len({g2chrom.get(g) for g in gs if g in g2chrom}) > 1}
    g2disc = {}
    for fi, m in enumerate(families):
        for g in m:
            g2disc[g] = fi
    recovered = []
    for f, gs in uni_cross.items():
        disc = defaultdict(int)
        for g in gs:
            if g in g2disc:
                disc[g2disc[g]] += 1
        best = max(disc.values()) if disc else 0
        recovered.append((f, len(gs), best))
    n_rec = sum(1 for _, _, b in recovered if b >= 2)
    print(f"[validation] universe cross-chrom families: {len(uni_cross)}; "
          f"recovered (>=2 genes co-grouped): {n_rec}")
    for f, n, b in sorted(recovered, key=lambda r: -r[2]):
        print(f"    {f}: {n} genes, {b} co-grouped {'OK' if b>=2 else 'MISS'}")

    # dump cross-chrom confirmed pairs (the headline output)
    out = os.path.join(os.path.dirname(__file__), "crosschrom_pairs.tsv")
    with open(out, "w") as fh:
        fh.write("gene_a\tchrom_a\tgene_b\tchrom_b\tcore_cov\tjaccard\tuniv_fam_a\tuniv_fam_b\n")
        for a, b, c, j in sorted(cross_pairs, key=lambda r: -r[2]):
            fh.write(f"{a}\t{g2chrom.get(a)}\t{b}\t{g2chrom.get(b)}\t{c}\t{j}\t"
                     f"{g2fam.get(a,'-')}\t{g2fam.get(b,'-')}\n")
    print(f"[wrote] {out} ({len(cross_pairs)} cross-chromosome confirmed pairs)")


if __name__ == "__main__":
    main()
