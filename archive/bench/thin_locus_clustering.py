#!/usr/bin/env python3
"""De novo clustering of thin loci (support 1-2) into multi-copy gene families.

Extracts genome-wide thin loci, builds their spliced sequences, finds homologous
pairs via canonical k-mer pre-filter + POA, and clusters into families. The goal
is to discover multi-copy families that are entirely missed by the conflict-graph
pipeline because every copy falls below the >=3 read assembly gate.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/thin_locus_clustering.py
"""
import csv
import os
import sys
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import multiprocessing as mp
import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P
from denovo_families import kmer_hashes, KMER

SCRATCH = "/home/juanfra/winloci_scratch"
HERE = os.path.dirname(__file__)
GGO_FA = os.path.join(SCRATCH, "GGO.fasta")
GGO_BAM = os.path.join(SCRATCH, "GGO.bam")
OUT_LOCI = os.path.join(SCRATCH, "thin_loci_genome_wide.tsv")
OUT_FAMILIES = os.path.join(SCRATCH, "thin_locus_families.tsv")
OUT_COPIES = os.path.join(SCRATCH, "thin_locus_copies.tsv")

MIN_SUPPORT = 1
MAX_SUPPORT = 2
MIN_LEN = 200
LEN_CAP = 9000
PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}
_RC = str.maketrans("ACGTN", "TGCAN")

# k-mer / clustering params
K_SHARE = 6          # propose a pair if two loci share >= this many informative k-mers
CNT_MAX = 40         # drop k-mers owned by > this many loci (repeats/TEs)
TOP_N = 50           # max candidate partners per locus for POA cost guard
T_CORE = 0.13        # POA contiguous-core threshold

# exon-structure guard: real paralog copies have similar exon count and spliced
# length; domain-sharers and read-throughs often do not.
MAX_EXON_DIFF = 2    # |n_exon_A - n_exon_B| must be <= this
MIN_LEN_RATIO = 0.5  # shorter / longer spliced length must be >= this
MAX_LEN_RATIO = 2.0  # shorter / longer spliced length must be <= this
MIN_SHARED_INTRONS = 1  # require at least this many identical intron coordinate pairs


def build_spliced_seq(fa, chrom, introns, start, end, strand):
    exons = []
    prev = start
    for d, a in introns:
        exons.append((prev, d))
        prev = a
    exons.append((prev, end))
    seq = "".join(fa.fetch(chrom, xs, xe).upper() for xs, xe in exons)
    if strand == "-":
        seq = seq.translate(_RC)[::-1]
    return seq


def canonical_kmer_set(seq):
    return set(kmer_hashes(seq).tolist()) if seq else set()


def extract_thin_loci(bam_path, fa_path, out_tsv, max_support=2, threads=4):
    """Genome-wide thin-locus extraction. Writes OUT_LOCI and returns list."""
    if os.path.exists(out_tsv):
        print(f"Reusing existing {out_tsv}")
        return load_thin_loci(out_tsv)

    bam = pysam.AlignmentFile(bam_path, "rb", threads=threads)
    fa = pysam.FastaFile(fa_path)

    loci = []
    for chrom in bam.references:
        clen = bam.lengths[bam.get_tid(chrom)]
        chains = defaultdict(lambda: [0, 10**12, 0])
        # fetch whole chromosome
        for aln in bam.fetch(chrom, 0, clen):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            rp = aln.reference_start
            intr = []
            for op, length in aln.cigartuples:
                if op == 3:
                    intr.append((rp, rp + length))
                    rp += length
                elif op in (0, 2, 7, 8):
                    rp += length
            if not intr:
                continue
            key = tuple(intr)
            chains[key][0] += 1
            chains[key][1] = min(chains[key][1], aln.reference_start)
            chains[key][2] = max(chains[key][2], aln.reference_end)

        for intr, (sup, s, e) in chains.items():
            if not (MIN_SUPPORT <= sup <= max_support):
                continue
            # strand inference
            strand = None
            ok = True
            for d, a in intr:
                don = fa.fetch(chrom, d, d + 2).upper()
                acc = fa.fetch(chrom, a - 2, a).upper()
                st = "+" if (don, acc) in PLUS else ("-" if (don, acc) in MINUS else None)
                if st is None or (strand and st != strand):
                    ok = False
                    break
                strand = st
            if not ok:
                continue
            seq = build_spliced_seq(fa, chrom, intr, s, e, strand)
            if len(seq) < MIN_LEN or len(seq) > LEN_CAP:
                continue
            loci.append({
                "lid": len(loci),
                "chrom": chrom,
                "start": s,
                "end": e,
                "strand": strand,
                "support": sup,
                "n_exon": len(intr) + 1,
                "seq": seq,
                "introns": intr,
            })
        print(f"  {chrom}: {len(loci)} total thin loci so far", flush=True)

    bam.close()
    fa.close()

    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["lid", "chrom", "start", "end", "strand", "support", "n_exon", "introns", "seq"])
        for L in loci:
            intr_str = ",".join(f"{d}-{a}" for d, a in L["introns"])
            w.writerow([L["lid"], L["chrom"], L["start"], L["end"], L["strand"],
                        L["support"], L["n_exon"], intr_str, L["seq"]])
    print(f"Wrote {out_tsv} ({len(loci)} loci)")
    return loci


def load_thin_loci(path):
    loci = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            introns = []
            if "introns" in row and row["introns"]:
                for part in row["introns"].split(","):
                    d, a = part.split("-")
                    introns.append((int(d), int(a)))
            loci.append({
                "lid": int(row["lid"]),
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "support": int(row["support"]),
                "n_exon": int(row["n_exon"]),
                "introns": introns,
                "seq": row["seq"],
            })
    return loci


def build_kmer_index(loci):
    """Return informative k-mers per locus + inverted index with repeat filter."""
    locus_kmers = []
    kmer_owners = defaultdict(set)
    for L in loci:
        km = canonical_kmer_set(L["seq"])
        locus_kmers.append(km)
        for k in km:
            kmer_owners[k].add(L["lid"])

    # filter pervasive k-mers
    info_kmers = {k for k, owners in kmer_owners.items() if 2 <= len(owners) <= CNT_MAX}
    print(f"Total k-mers: {len(kmer_owners)}; informative: {len(info_kmers)}", flush=True)

    # per-locus informative signature
    sig = [km & info_kmers for km in locus_kmers]
    return sig


def candidate_pairs(loci, sig):
    """For each locus, return candidate partner lids sharing >= K_SHARE informative k-mers
    and passing exon-structure similarity guards."""
    # inverted index over informative k-mers
    idx = defaultdict(list)
    for lid, km in enumerate(sig):
        for k in km:
            idx[k].append(lid)

    def exon_compatible(A, B):
        # same strand required for paralog copies
        if A["strand"] != B["strand"]:
            return False
        # similar exon count
        if abs(A["n_exon"] - B["n_exon"]) > MAX_EXON_DIFF:
            return False
        # similar spliced length
        la, lb = len(A["seq"]), len(B["seq"])
        if la == 0 or lb == 0:
            return False
        ratio = min(la, lb) / max(la, lb)
        if ratio < MIN_LEN_RATIO or ratio > MAX_LEN_RATIO:
            return False
        # at least one identical intron coordinate pair (real paralogs share splice sites)
        if A.get("introns") and B.get("introns"):
            if not (set(A["introns"]) & set(B["introns"])):
                return False
        return True

    pairs = defaultdict(int)
    for lid, km in enumerate(sig):
        if not km:
            continue
        A = loci[lid]
        # count co-occurrences
        partner_counts = defaultdict(int)
        for k in km:
            for pid in idx[k]:
                if pid != lid:
                    partner_counts[pid] += 1
        # keep partners above k-mer threshold AND exon-compatible
        cand = []
        for pid, cnt in partner_counts.items():
            if cnt < K_SHARE:
                continue
            if not exon_compatible(A, loci[pid]):
                continue
            cand.append((cnt, pid))
        cand.sort(reverse=True)
        for cnt, pid in cand[:TOP_N]:
            a, b = sorted((lid, pid))
            if cnt > pairs[(a, b)]:
                pairs[(a, b)] = cnt
    print(f"Candidate pairs after exon guard: {len(pairs)}", flush=True)
    return list(pairs.keys())


_A = None

def _init_poa():
    global _A
    _A = P.make_aligner()


def _poa_edge(task):
    i, j, loci = task
    s1 = loci[i]["seq"]
    s2 = loci[j]["seq"]
    if min(len(s1), len(s2)) > LEN_CAP:
        return None
    cr = P.poa_pair_stats(_A, s1, s2)["core_recip"]
    if cr < T_CORE:
        cr = max(cr, P.poa_pair_stats(_A, s1, s2.translate(_RC)[::-1])["core_recip"])
    return (i, j, round(cr, 4)) if cr >= T_CORE else None


def find_components(n, edges):
    parent = list(range(n))
    rank = [0] * n

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    for a, b in edges:
        union(a, b)

    comps = defaultdict(list)
    for i in range(n):
        comps[find(i)].append(i)
    return list(comps.values())


def same_locus(A, B):
    """Two loci are the same locus if they share a junction or one span contains the other."""
    if A["chrom"] != B["chrom"] or A["strand"] != B["strand"]:
        return False
    # share a junction
    if set(A.get("introns", [])) & set(B.get("introns", [])):
        return True
    # containment: one span is >=90% inside the other
    def containment(a, b):
        ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
        return ov >= 0.9 * (a["end"] - a["start"]) or ov >= 0.9 * (b["end"] - b["start"])
    return containment(A, B)


def distinct_loci(comp_loci):
    """Collapse same-locus loci; return representative per distinct locus (widest span)."""
    groups = []
    for L in comp_loci:
        merged = False
        for g in groups:
            if same_locus(L, g[0]):
                g.append(L)
                merged = True
                break
        if not merged:
            groups.append([L])
    reps = []
    for g in groups:
        # pick widest span
        rep = max(g, key=lambda x: x["end"] - x["start"])
        reps.append(rep)
    return reps


def cluster_loci(chrom_loci):
    """Cluster loci on one chromosome into families. Returns list of rep lists."""
    n = len(chrom_loci)
    if n < 2:
        return []
    sig = build_kmer_index(chrom_loci)
    pairs = candidate_pairs(chrom_loci, sig)
    if not pairs:
        return []

    tasks = [(i, j, chrom_loci) for i, j in pairs]
    with mp.Pool(5, initializer=_init_poa) as pool:
        edges = [e for e in pool.map(_poa_edge, tasks, chunksize=64) if e]
    if not edges:
        return []

    edge_list = [(i, j) for i, j, _ in edges]
    components = find_components(n, edge_list)

    families = []
    for comp in components:
        if len(comp) < 2:
            continue
        comp_loci = [chrom_loci[i] for i in comp]
        reps = distinct_loci(comp_loci)
        if len(reps) >= 2:
            families.append(reps)
    return families


def main():
    loci = extract_thin_loci(GGO_BAM, GGO_FA, OUT_LOCI, max_support=MAX_SUPPORT, threads=4)
    print(f"Genome-wide thin loci passing filters: {len(loci)}", flush=True)

    # Process per chromosome to keep memory bounded.
    by_chrom = defaultdict(list)
    for L in loci:
        by_chrom[L["chrom"]].append(L)

    all_families = []
    for chrom in sorted(by_chrom.keys()):
        chrom_loci = by_chrom[chrom]
        print(f"Clustering {chrom}: {len(chrom_loci)} loci ...", flush=True)
        families = cluster_loci(chrom_loci)
        print(f"  -> {len(families)} families", flush=True)
        all_families.extend(families)

    print(f"Total families with >=2 distinct loci: {len(all_families)}", flush=True)

    # write families/copies
    with open(OUT_FAMILIES, "w", newline="") as ff, open(OUT_COPIES, "w", newline="") as fc:
        wf = csv.writer(ff, delimiter="\t")
        wc = csv.writer(fc, delimiter="\t")
        wf.writerow(["family_id", "n_copies", "n_chroms", "chroms", "cross_chrom", "avg_reads"])
        wc.writerow(["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads"])
        for fi, reps in enumerate(sorted(all_families, key=lambda r: (r[0]["chrom"], r[0]["start"]))):
            fid = f"THINFAM{fi}"
            reps.sort(key=lambda x: (x["chrom"], x["start"]))
            chroms = sorted(set(r["chrom"] for r in reps))
            cross = len(chroms) > 1
            avg_reads = sum(r["support"] for r in reps) / len(reps)
            wf.writerow([fid, len(reps), len(chroms), ",".join(chroms), cross, f"{avg_reads:.1f}"])
            for ci, r in enumerate(reps):
                tid = f"TL_{r['chrom']}_{r['start']}_{r['n_exon']}"
                wc.writerow([fid, ci, tid, r["chrom"], r["start"], r["end"],
                             r["n_exon"], r["strand"], r["support"]])
    print(f"Wrote {OUT_FAMILIES} and {OUT_COPIES}")


if __name__ == "__main__":
    main()
