#!/usr/bin/env python
"""A1 read-heterogeneity PSV caller (SDA / Vollger-2019 correlation-clustering), NO asm20.

At a COLLAPSED reference locus (one locus, MAPQ-0 reads from K copies piling up), discover
copy-distinguishing columns FROM READ DISAGREEMENT, link them by cross-read correlation, and
partition reads into copies. This is the A1 read-derived analog of the (reference-seeded) asm20
PSV columns used by psv_graph_genomewide.py. Columns here are DISCOVERED from the pileup, not
seeded by any reference-copy alignment.

Smoke test only: 2-3 collapsed families, report method feasibility + read-derived K vs genome n_loci.

Run: /home/juanfra/miniforge3/bin/python bench/a1_read_sda_smoketest.py
"""
import collections
import sys

import numpy as np
import pysam

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"

# ---- method parameters (documented in the return) --------------------------------------------
MIN_BQ        = 20      # HiFi base-quality floor (Q20 = 1% err); drops low-qual basecalls
MIN_DEPTH     = 12      # a column must be covered by >= this many reads to be callable
MIN_MINOR_CNT = 4       # minor allele must appear in >= this many reads (kills singleton errors)
MIN_VAF       = 0.20    # minor allele fraction >= this (a real copy is ~1/K; error VAF ~ 1/depth)
PHI_LINK      = 0.60    # |phi| correlation to call two het columns LINKED (co-segregating)
MIN_COCOVER   = 10      # min reads covering BOTH columns to trust a phi
MIN_LINKED    = 2       # a PSV must link to >= this many other PSVs to be kept (error filter)
READ_MATCH    = 0.80    # read-read agreement over linked PSVs to join same copy-cluster
MIN_CLUST     = 4       # a read cluster must have >= this many reads to count as a copy
INCLUDE_SECONDARY = True  # include secondary alns: MAPQ-0 multi-mappers ARE the collapse mixture

# pileup flag filter: drop unmapped/qcfail/dup; keep (or drop) secondary
FLAG_KEEP_SEC  = 0x4 | 0x200 | 0x400          # UNMAP|QCFAIL|DUP  (secondary KEPT)
FLAG_DROP_SEC  = 0x4 | 0x100 | 0x200 | 0x400  # + SECONDARY dropped


def allele_matrix(bam, chrom, start, end, include_secondary):
    """dict[pos] -> {read_name: base}  from the read pileup. ACGT only, base-qual>=MIN_BQ.
    Secondary alignments optionally included (the collapsed MAPQ-0 mixture)."""
    flag = FLAG_KEEP_SEC if include_secondary else FLAG_DROP_SEC
    cols = collections.defaultdict(dict)
    depth = collections.Counter()
    for col in bam.pileup(chrom, start, end, min_base_quality=MIN_BQ, truncate=True,
                          max_depth=200000, flag_filter=flag, stepper="nofilter"):
        p = col.reference_pos
        for pr in col.pileups:
            if pr.query_position is None:
                continue
            aln = pr.alignment
            if aln.query_sequence is None:
                continue
            b = aln.query_sequence[pr.query_position]
            q = aln.query_qualities[pr.query_position] if aln.query_qualities is not None else 40
            if b in "ACGT" and q >= MIN_BQ:
                # unique key per physical read (secondary shares query_name w/ primary elsewhere,
                # but at THIS locus a read has one aln, so query_name is unique here)
                cols[p][aln.query_name] = b
                depth[p] += 1
    return cols, depth


def call_het_columns(cols, depth):
    """Candidate het columns: minor allele >= MIN_MINOR_CNT reads AND VAF >= MIN_VAF at depth>=MIN_DEPTH.
    Returns list of (pos, major_base, minor_base, minor_cnt, dp)."""
    cands = []
    for p, dp in depth.items():
        if dp < MIN_DEPTH:
            continue
        cnt = collections.Counter(cols[p].values())
        top = cnt.most_common()
        if len(top) < 2:
            continue
        (maj, majc), (minb, minc) = top[0], top[1]
        if minc >= MIN_MINOR_CNT and minc / dp >= MIN_VAF:
            cands.append((p, maj, minb, minc, dp))
    return sorted(cands)


def phi(a, b):
    """phi coefficient for two aligned binary vectors (2x2 contingency)."""
    a = np.asarray(a); b = np.asarray(b)
    n11 = int(np.sum((a == 1) & (b == 1))); n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1))); n00 = int(np.sum((a == 0) & (b == 0)))
    num = n11 * n00 - n10 * n01
    den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    return 0.0 if den == 0 else num / den


def link_columns(cols, cands):
    """SDA correlation step: binarize each het column (major=0, minor=1), compute pairwise phi over
    co-covering reads, keep columns that co-segregate (|phi|>=PHI_LINK) with >=MIN_LINKED others."""
    # binary read->{colidx: 0/1}
    binv = []
    for (p, maj, minb, minc, dp) in cands:
        d = {}
        for rn, base in cols[p].items():
            if base == maj:
                d[rn] = 0
            elif base == minb:
                d[rn] = 1
        binv.append(d)
    n = len(cands)
    link_deg = [0] * n
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            shared = set(binv[i]) & set(binv[j])
            if len(shared) < MIN_COCOVER:
                continue
            sl = sorted(shared)
            ai = [binv[i][r] for r in sl]; aj = [binv[j][r] for r in sl]
            ph = phi(ai, aj)
            if abs(ph) >= PHI_LINK:
                link_deg[i] += 1; link_deg[j] += 1
                edges.append((i, j, ph))
    keep = [k for k in range(n) if link_deg[k] >= MIN_LINKED]
    return keep, edges, binv


def cluster_reads(cands, binv, keep):
    """Partition reads into copies by agreement over the LINKED PSV columns (read-derived K)."""
    kset = set(keep)
    # read -> allele vector over kept columns (dict colidx->0/1)
    read_vec = collections.defaultdict(dict)
    for k in keep:
        for rn, v in binv[k].items():
            read_vec[rn][k] = v
    reads = [r for r in read_vec if len(read_vec[r]) >= max(2, int(0.4 * len(keep)))]
    # union-find on reads that agree over shared kept columns
    parent = {r: r for r in reads}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)
    for i in range(len(reads)):
        vi = read_vec[reads[i]]
        for j in range(i + 1, len(reads)):
            vj = read_vec[reads[j]]
            shared = set(vi) & set(vj)
            if len(shared) < max(2, int(0.4 * len(keep))):
                continue
            agree = sum(1 for c in shared if vi[c] == vj[c]) / len(shared)
            if agree >= READ_MATCH:
                union(reads[i], reads[j])
    groups = collections.defaultdict(list)
    for r in reads:
        groups[find(r)].append(r)
    sizes = sorted((len(v) for v in groups.values()), reverse=True)
    k_read = sum(1 for s in sizes if s >= MIN_CLUST)
    return k_read, sizes, len(reads)


def shuffle_cols(cols, cands, rng):
    """NULL control: independently permute each column's alleles ACROSS reads. Preserves per-column
    allele frequencies (so the same columns stay 'het') but DESTROYS cross-read co-segregation.
    Real SDA linkage -> collapses to ~0; a pileup/coverage artifact -> survives."""
    scols = {}
    for (p, maj, minb, minc, dp) in cands:
        items = list(cols[p].items())
        names = [k for k, _ in items]
        vals = [v for _, v in items]
        rng.shuffle(vals)
        scols[p] = dict(zip(names, vals))
    # merge non-candidate columns unchanged (not used downstream, but keep interface)
    return scols


def analyze_locus(bam, chrom, start, end, include_secondary, rng=None):
    cols, depth = allele_matrix(bam, chrom, start, end, include_secondary)
    cands = call_het_columns(cols, depth)
    if not cands:
        return dict(n_cand=0, n_linked=0, k_read=0, clust_sizes=[], n_reads=0,
                    med_depth=int(np.median(list(depth.values()))) if depth else 0,
                    null_linked=0, null_phi=0.0, null_k=0)
    keep, edges, binv = link_columns(cols, cands)
    k_read, sizes, n_reads = cluster_reads(cands, binv, keep)
    # NULL: same candidate columns, alleles shuffled across reads
    null_linked = null_phi = null_k = None
    if rng is not None:
        scols = shuffle_cols(cols, cands, rng)
        nkeep, nedges, nbinv = link_columns(scols, cands)
        nk, _, _ = cluster_reads(cands, nbinv, nkeep) if nkeep else (0, [], 0)
        null_linked = len(nkeep)
        null_phi = round(float(np.mean([abs(e[2]) for e in nedges])), 3) if nedges else 0.0
        null_k = nk
    return dict(n_cand=len(cands), n_linked=len(keep), n_edges=len(edges),
                k_read=k_read, clust_sizes=sizes[:8], n_reads=n_reads,
                med_depth=int(np.median([c[4] for c in cands])),
                mean_phi=round(float(np.mean([abs(e[2]) for e in edges])), 3) if edges else 0.0,
                null_linked=null_linked, null_phi=null_phi, null_k=null_k)


# collapsed families -> the loci where A1's independence lives. Windows target the EXONIC
# high-depth region of the backbone locus (RNA reads only pile up on exons).
# (family, gene, n_loci_genome, asm20_psv, (chrom, win_start, win_end))
TESTS = [
    # fam1: 7-copy segdup NC_073236.2; asm20 found ZERO PSVs (no_psv) -> read signal is PURE
    (1,  "LOC109025447", 7, 0, ("NC_073236.2", 51182668, 51205668)),
    # fam14: LOC115930164 4-copy tandem array; asm20 ZERO PSVs (no_psv). exons at 3' end
    (14, "LOC115930164", 4, 0, ("NC_073224.2", 223543000, 223553675)),
    # fam11: 4-copy segdup NC_073230.2; asm20 found 13 PSVs -> tests read/asm20 overlap
    (11, "LOC101130894", 4, 13, ("NC_073230.2", 86213469, 86236469)),
]


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    rng = np.random.default_rng(0)
    print("A1 READ-HETEROGENEITY PSV CALLER (SDA, no asm20) -- SMOKE TEST")
    print(f"params: MIN_BQ={MIN_BQ} MIN_DEPTH={MIN_DEPTH} MIN_MINOR_CNT={MIN_MINOR_CNT} "
          f"MIN_VAF={MIN_VAF} PHI_LINK={PHI_LINK} MIN_LINKED={MIN_LINKED} READ_MATCH={READ_MATCH}")
    print("null = alleles shuffled ACROSS reads within each column (destroys co-segregation)\n")
    print(f"{'fam':>4} {'gene':<15} {'nloci':>5} {'asm20':>5} {'meddp':>6} {'cand':>5} "
          f"{'link':>4} {'phi':>5} {'Kread':>5} | {'nullink':>7} {'nullphi':>7} {'nullK':>5} | clusters")
    for fam, gene, n_loci, asm20, (chrom, w1, w2) in TESTS:
        r = analyze_locus(bam, chrom, w1, w2, INCLUDE_SECONDARY, rng=rng)
        print(f"{fam:>4} {gene:<15} {n_loci:>5} {asm20:>5} {r['med_depth']:>6} "
              f"{r['n_cand']:>5} {r['n_linked']:>4} {r.get('mean_phi',0):>5} {r['k_read']:>5} | "
              f"{r['null_linked']:>7} {r['null_phi']:>7} {r['null_k']:>5} | {r['clust_sizes']}")
        sys.stdout.flush()
    bam.close()


if __name__ == "__main__":
    main()
