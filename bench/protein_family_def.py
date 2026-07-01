#!/usr/bin/env python
r"""PROTEIN-LEVEL multi-copy gene family definition on the haploid T2T genome (GCF_029281585.2).

This is the field-standard "gene family" notion (OrthoFinder / Pfam / mmseqs-cluster lineage):
two genes are in the same family iff their PROTEIN products are homologous. It is the NEW, fourth
oracle E_p in the family lattice, and it is built with the *same combinatorial skeleton and the same
cohesion operator R* as the DNA-level segdup family (bench/genome_family_def.py,
bench/DNA_FAMILY_DEFINITION_FORMAL.md) -- nodes = genes, family = refined connected component of a
homology graph -- differing ONLY in the edge oracle:

    E_a (DNA)     edge = a SEDEF segdup pair links the two gene loci        (read-independent, genomic)
    E_p (PROTEIN) edge = the two genes' PROTEINS are reciprocally homologous  (this file; coding homology)
    E_c (RNA)     edge = a read de-ties the two loci (read-conflict)          (operational copy unit)

Because /home/juanfra/winloci_scratch/proteins.fa headers ARE gene symbols (>A1BG, >A2M) and equal the
GFF `gene=` names, E_p, E_a and E_c live on ONE shared node set with zero mapping (22,614/22,614 map).

EDGE ORACLE E_p (two irreducible parameters, named honestly -- the analog of the segdup E-value + COVER):
    {u,v} in E_p  iff  the mmseqs all-vs-all reports hits in BOTH orientations (u->v AND v->u), each with
        E-value <= ALPHA_P               (significance floor; the genome-wide homology tail)
        qcov >= C_P AND tcov >= C_P       (whole-protein coverage of BOTH proteins -- NOT a shared domain)
Reciprocity (both orientations) + whole-protein reciprocal coverage are BOTH load-bearing:
  * The reciprocal whole-protein coverage clause demands coding homology over most of BOTH proteins, so the
    DGCR6<->DRD5<->SLC6A8 non-coding 22q11 segdup blob (1,547 genes in E_a's GFAM0) does NOT form in E_p
    (they are not mutually protein-homologous). E_p is CLEANER than the raw E_a superset by construction.
  * There is NO identity (fident) floor. Significance + reciprocal coverage alone define homology. This is
    exactly what lets E_p reach ANCIENT/divergent paralogs that E_a (>=90%, or even its 50% floor) misses:
    the globin superfamily (HBZ/HBM/HBQ1/HBE1/MB/CYGB/NGB, 5 chromosomes) is one E_p family down into the
    twilight zone (MB-HBE1 24% aa id, E~1e-7), so E_p \ E_a is non-empty and biologically real.

ALPHA_P is a statistically calibratable Type-I level, NOT an ad-hoc identity cliff: mmseqs's E-value is
already database-size corrected, so the expected number of chance (false-homology) DIRECTIONAL hits across
the whole all-vs-all is ~ ALPHA_P * n_queries = 1e-5 * 22,614 ~= 0.23 < 1 (Bonferroni-safe); requiring BOTH
orientations makes a chance edge far rarer still. C_P=0.50 is the one genuinely arbitrary knob (the
domain-vs-whole-protein boundary) -- the coding analog of the segdup COVER=0.50 gene-projection fraction.
mmseqs's BLOSUM62 + gap penalties are the scoring scheme S(.) that fixes the length<->identity exchange
rate -- the same irreducible-scoring caveat as the segdup Karlin-Altschul E-value. NO threshold-freedom is
claimed: ALPHA_P, C_P and gamma are named parameters, and their sensitivity is reported.

REFINEMENT: the SAME operator R as the DNA level (gamma-quasi-clique cohesion + >=2-distinct-loci multi-copy
predicate), imported verbatim from genome_family_def.refine_families -- so protein & DNA families are ONE
refinement operator over two oracles. gamma=0.20 TRANSFERS unchanged (it is near-inert on the protein graph:
reciprocal whole-protein homology yields raw components that are already near-cliques; gamma fires only on
the handful of domain-superfamily mega-components, the protein analog of the E_a repeat-array over-merge).

CROSS-VALIDATION: the E_p catalog is checked against a field-standard mmseqs *cluster* (greedy set-cover,
the standard protein-clustering baseline) by dominant-cluster purity, fully-contained fraction, and the
Adjusted Rand Index over the E_p member genes.

Run: /home/juanfra/miniforge3/bin/python bench/protein_family_def.py --refine
Inputs (in /home/juanfra/winloci_scratch): proteins.fa, aln.m8 (mmseqs easy-search, cols
query,target,evalue,pident,qcov,tcov), protclust_cluster.tsv (mmseqs easy-cluster, rep<TAB>member).
"""
import argparse
import itertools
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genome_family_def import (  # SAME operator + node loader as the DNA level
    load_genes, refine_families, distinct_loci, biotype_purity,
    _covered_genes, _init_worker, GFF, SEDEF, GAMMA, SEED, COVER, PURITY, UF,
)

SCRATCH = "/home/juanfra/winloci_scratch"
PROT_FA = os.path.join(SCRATCH, "proteins.fa")
ALN = os.path.join(SCRATCH, "aln.m8")                     # mmseqs easy-search all-vs-all (query,target,evalue,pident,qcov,tcov)
MMCLUST = os.path.join(SCRATCH, "protclust_cluster.tsv")  # mmseqs easy-cluster (field-standard set-cover clustering)

ALPHA_P = 1e-5  # SIGNIFICANCE floor (irreducible, honestly named): the homology-tail E-value cutoff, the
                # protein analog of the segdup Karlin-Altschul E-value. Calibratable Type-I level, not an
                # identity cliff: expected chance directional hits genome-wide ~ ALPHA_P*n_queries ~= 0.23 < 1
                # (Bonferroni-safe). Sensitivity reported (ALPHA_P in {1e-3,1e-5}) so the choice is not free.
C_P = 0.50      # RECIPROCAL WHOLE-PROTEIN coverage (irreducible): each retained alignment must cover >=C_P of
                # BOTH proteins. The domain-vs-whole-protein cut -- coding analog of the segdup COVER=0.50.
                # This is exactly what excludes single-domain sharers and the non-coding 22q11 blob from E_p.


# ----------------------------- inputs -----------------------------
def load_protein_symbols(fa=PROT_FA):
    """proteins.fa header = gene symbol (one protein per gene). Returns the list of symbols."""
    syms = []
    with open(fa) as fh:
        for ln in fh:
            if ln.startswith(">"):
                syms.append(ln[1:].strip())
    return syms


def build_name_index(genes):
    """gene symbol -> gene index (first occurrence; the duplicate names are tRNAs, never proteins)."""
    idx = {}
    for i, g in enumerate(genes):
        idx.setdefault(g["name"], i)
    return idx


def load_ep_edges(aln, name2idx, alpha_p=ALPHA_P, c_p=C_P):
    """mmseqs easy-search m8 (query,target,evalue,pident,qcov,tcov) -> E_p gene-index edge set.

    A directional hit (q,t) is retained iff q!=t, evalue<=alpha_p, qcov>=c_p AND tcov>=c_p.
    Edge {u,v} in E_p iff BOTH (q=u,t=v) AND (q=v,t=u) are retained (RECIPROCITY -- load-bearing:
    dropping it lets one-directional partial hits inflate the mega-components). Symmetric, irreflexive."""
    fwd = set()             # retained directional (gene_idx_q, gene_idx_t)
    n_rows = n_dir = 0
    with open(aln) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, t = f[0], f[1]
            if q == t:
                continue
            n_rows += 1
            try:
                ev, qc, tc = float(f[2]), float(f[4]), float(f[5])
            except ValueError:
                continue
            if ev > alpha_p or qc < c_p or tc < c_p:
                continue
            iu, iv = name2idx.get(q), name2idx.get(t)
            if iu is None or iv is None or iu == iv:
                continue
            n_dir += 1
            fwd.add((iu, iv))
    edges = set()
    for (u, v) in fwd:
        if (v, u) in fwd:
            edges.add((u, v) if u < v else (v, u))
    n_oneway = sum(1 for (u, v) in fwd if (v, u) not in fwd)
    return edges, n_rows, n_dir, n_oneway


def build_raw_families(edges):
    """connected components of an edge set with |C| >= 2 (same union-find skeleton as E_a / E_c)."""
    uf = UF()
    for u, v in edges:
        uf.union(u, v)
    fams = [sorted(set(m)) for m in uf.groups().values() if len(set(m)) >= 2]
    return fams


def fam_index(fams):
    return {i: k for k, m in enumerate(fams) for i in m}


def induced_adj(edges, fams):
    """adjacency restricted to INTRA-family E_p edges (for density / certificate)."""
    fof = fam_index(fams)
    adj = defaultdict(set)
    for u, v in edges:
        if fof.get(u) is not None and fof.get(u) == fof.get(v):
            adj[u].add(v); adj[v].add(u)
    return adj


def density(members, adj):
    m = set(members)
    n = len(m)
    if n <= 1:
        return 1.0
    e = sum(1 for u in m for v in adj.get(u, ()) if v in m and u < v)
    return 2 * e / (n * (n - 1))


# ----------------------------- mmseqs-cluster cross-validation (field standard) -----------------------------
def load_mmseqs_clusters(path=MMCLUST):
    """easy-cluster tsv: col1 = representative, col2 = member (both gene symbols) -> {symbol: cluster_rep}."""
    rep_of = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        rep_of[f[1]] = f[0]
    return rep_of


def adjusted_rand(labels_a, labels_b):
    """Adjusted Rand Index between two label dicts over their shared node set (manual, no sklearn)."""
    nodes = set(labels_a) & set(labels_b)
    cont = defaultdict(int)
    ca, cb = Counter(), Counter()
    for n in nodes:
        a, b = labels_a[n], labels_b[n]
        cont[(a, b)] += 1
        ca[a] += 1
        cb[b] += 1
    c2 = lambda x: x * (x - 1) // 2
    sum_ij = sum(c2(v) for v in cont.values())
    sum_a = sum(c2(v) for v in ca.values())
    sum_b = sum(c2(v) for v in cb.values())
    tot = c2(len(nodes))
    if tot == 0:
        return 1.0
    exp = sum_a * sum_b / tot
    maxi = (sum_a + sum_b) / 2
    return (sum_ij - exp) / (maxi - exp) if (maxi - exp) != 0 else 1.0


def crossval(fams, genes, rep_of, large=30):
    """Field-standard cross-validation of the E_p (gamma-community) families vs the mmseqs cluster partition.

    Reports THREE purity views (a per-family mean alone is inflated by the many size-2 families, so we also
    report the gene-weighted purity and the large-family mean -- the multi-copy families that carry the thesis):
      * mean_pur       -- UNWEIGHTED mean of per-family dominant-cluster purity (each family counts once).
      * gw_pur         -- GENE-WEIGHTED dominant-cluster purity = sum(dominant counts) / sum(reps); the
                          genome-wide fraction of E_p member genes that lie in their family's dominant cluster.
      * large_pur      -- mean per-family purity over families with >= `large` member genes (the big,
                          biologically-interesting superfamilies, where the granularity residual concentrates).
    Also returns #families fully inside one cluster and the ARI over E_p member genes."""
    pur = []
    n_in_one = 0
    large_pur, n_large = [], 0
    dom_sum = rep_sum = 0                       # for the gene-weighted purity
    lab_ep, lab_mm = {}, {}
    for k, m in enumerate(fams):
        reps = [rep_of.get(genes[i]["name"]) for i in m]
        reps = [r for r in reps if r is not None]
        if reps:
            dom = Counter(reps).most_common(1)[0][1]
            p = dom / len(reps)
            pur.append(p)
            dom_sum += dom
            rep_sum += len(reps)
            if len(set(reps)) == 1:
                n_in_one += 1
            if len(m) >= large:
                large_pur.append(p); n_large += 1
        else:
            pur.append(0.0)
        for i in m:
            r = rep_of.get(genes[i]["name"])
            if r is not None:
                lab_ep[genes[i]["name"]] = k
                lab_mm[genes[i]["name"]] = r
    mean_pur = sum(pur) / len(pur) if pur else 0.0
    gw_pur = dom_sum / rep_sum if rep_sum else 0.0
    large_mean = sum(large_pur) / len(large_pur) if large_pur else 0.0
    ari = adjusted_rand(lab_ep, lab_mm)
    return mean_pur, gw_pur, large_mean, n_large, n_in_one, ari


# ----------------------------- E_a (segdup) overlay + Bailey (>=90%) overlay -----------------------------
def _bailey_ok(f):
    """SEDEF Bailey SD(.) predicate: fracMatch>=0.90 (col21), aln_len>=1000 (col12), non-TE (col28/29>=0.5)."""
    try:
        aln_len, frac, um, am = int(f[11]), float(f[20]), int(f[27]), int(f[28])
    except (ValueError, IndexError):
        return False
    return frac >= 0.90 and aln_len >= 1000 and am > 0 and um / am >= 0.50


def load_ea_edges(genes_by, gstarts):
    """SINGLE pass over final.bed -> (E_a edge set, E_a-Bailey>=90% edge set) as gene-index pairs.
    A pair links gene u (>=COVER covered on one side) to gene v (>=COVER covered on the other)."""
    _init_worker(genes_by, gstarts)          # set the module globals _covered_genes reads
    ea, ea_b = set(), set()
    with open(SEDEF) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                cA, sA, eA, cB, sB, eB = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
            except ValueError:
                continue
            ga = _covered_genes(cA, sA, eA)
            if not ga:
                continue
            gb = _covered_genes(cB, sB, eB)
            if not gb:
                continue
            bok = _bailey_ok(f)
            for u in ga:
                for v in gb:
                    if u != v:
                        e = (u, v) if u < v else (v, u)
                        ea.add(e)
                        if bok:
                            ea_b.add(e)
    return ea, ea_b


def load_ec_edges(name2idx, path=None):
    """E_c read-conflict edges (gene-index pairs with n_conflict>0) from family_def_readconflict.tsv."""
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_readconflict.tsv")
    ec = set()
    if not os.path.exists(path):
        return ec
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ia, ib, ic = hdr.index("gene_a"), hdr.index("gene_b"), hdr.index("n_conflict")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= ic:
                continue
            try:
                if int(f[ic]) <= 0:
                    continue
            except ValueError:
                continue
            u, v = name2idx.get(f[ia]), name2idx.get(f[ib])
            if u is None or v is None or u == v:
                continue
            ec.add((u, v) if u < v else (v, u))
    return ec


def same_fam(idx, u, v):
    a, b = idx.get(u), idx.get(v)
    return a is not None and a == b


def grfam_containment(famP, genes, path=None):
    """Family-level operational-subset stat: how many canonical refined DNA families (GRFAM,
    genome_families_refined.tsv) with >=1 protein member sit FULLY inside one E_p family."""
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "genome_families_refined.tsv")
    if not os.path.exists(path):
        return None, None, None
    ep_of = {genes[i]["name"]: k for k, m in enumerate(famP) for i in m}
    tot = in_one = multi = 0
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi = hdr.index("genes")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= gi:
                continue
            ep_ids = {ep_of[g] for g in f[gi].split(",") if g in ep_of}
            if not ep_ids:
                continue
            tot += 1
            if len(ep_ids) == 1:
                in_one += 1
            else:
                multi += 1
    return tot, in_one, multi


# ----------------------------- catalog writer -----------------------------
def write_catalog(path, fams, genes, adj, rep_of, gamma, prefix):
    def named_dom(m):
        return Counter(genes[i]["biotype"] for i in m).most_common(1)[0][0]
    n_cert = 0
    with open(path, "w") as fh:
        fh.write("family_id\tn_genes\tn_loci\tn_contigs\tcross_chrom\tdensity\tgamma_cohesive\t"
                 "biotype_purity\trepeat_mosaic\tdom_biotype\tmmseqs_cluster_purity\tgenes\n")
        for k, m in enumerate(fams):
            d = density(m, adj)
            cert = int(len(m) <= 2 or d >= gamma)
            n_cert += cert
            pur = biotype_purity(m, genes)
            reps = [rep_of.get(genes[i]["name"]) for i in m]
            reps = [r for r in reps if r]
            mmpur = (Counter(reps).most_common(1)[0][1] / len(reps)) if reps else 0.0
            fh.write(f"{prefix}{k}\t{len(m)}\t{distinct_loci(m, genes)}\t"
                     f"{len({genes[i]['chrom'] for i in m})}\t{int(len({genes[i]['chrom'] for i in m})>=2)}\t"
                     f"{d:.3f}\t{cert}\t{pur:.3f}\t{int(pur < PURITY)}\t{named_dom(m)}\t{mmpur:.3f}\t"
                     f"{','.join(genes[i]['name'] for i in m)}\n")
    return n_cert


# ----------------------------- main -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evalue", type=float, default=ALPHA_P, help=f"E_p significance floor ALPHA_P (default {ALPHA_P})")
    ap.add_argument("--cov", type=float, default=C_P, help=f"E_p reciprocal whole-protein coverage C_P (default {C_P})")
    ap.add_argument("--gamma", type=float, default=GAMMA, help=f"cohesion parameter for R (default {GAMMA}, transferred from the DNA level)")
    ap.add_argument("--seed", type=int, default=SEED, help=f"splitter WITNESS seed (default {SEED})")
    ap.add_argument("--refine", action="store_true",
                    help="apply the SAME operator R (gamma-quasi-clique + >=2 loci) and ALSO write the canonical "
                         "refined catalog protein_families_refined.tsv (PRFAM). Raw protein_families.tsv (PFAM) "
                         "is always written; with --refine off it is the only catalog (baseline).")
    ap.add_argument("--no-lattice", action="store_true", help="skip the E_a/E_c lattice overlay (faster; witnesses off)")
    args = ap.parse_args()
    bench = os.path.dirname(os.path.abspath(__file__))

    print(f"[load] genes  {GFF}")
    genes, genes_by = load_genes(GFF)
    gstarts = {c: [x[0] for x in v] for c, v in genes_by.items()}
    name2idx = build_name_index(genes)
    syms = load_protein_symbols()
    mapped = sum(1 for s in syms if s in name2idx)
    print(f"       {len(genes)} gene nodes; {len(syms)} proteins; {mapped}/{len(syms)} map to a gene node by symbol")

    exp_false = args.evalue * len(syms)
    print(f"[E_p] edges from {ALN}  (ALPHA_P={args.evalue}, reciprocal cov>={args.cov}; "
          f"expected chance directional hits ~= {exp_false:.2f})")
    ep_edges, n_rows, n_dir, n_oneway = load_ep_edges(ALN, name2idx, args.evalue, args.cov)
    print(f"       {n_rows} non-self hits; {n_dir} pass filter; {n_oneway} one-directional-only (dropped by reciprocity) "
          f"-> {len(ep_edges)} E_p edges")

    famP_raw = build_raw_families(ep_edges)
    famP_raw.sort(key=lambda m: (-len(m), genes[m[0]]["name"]))
    adj_raw = induced_adj(ep_edges, famP_raw)
    print(f"       RAW: {len(famP_raw)} families ({sum(len(m) for m in famP_raw)} genes), max {len(famP_raw[0]) if famP_raw else 0}")

    famP = famP_raw
    if args.refine:
        print(f"[R] SAME operator as the DNA level (gamma={args.gamma}, seed={args.seed}, >=2 distinct loci) ...")
        famP = refine_families(famP_raw, list(ep_edges), genes, args.gamma, seed=args.seed)
        adj_ref = induced_adj(ep_edges, famP)
        cert = sum(1 for m in famP if len(m) <= 2 or density(m, adj_ref) >= args.gamma)
        print(f"       REFINED: {len(famP_raw)} raw -> {len(famP)} families ({sum(len(m) for m in famP)} genes), "
              f"max {len(famP[0]) if famP else 0}; certificate {cert}/{len(famP)} gamma-cohesive-or-size<=2")
    else:
        adj_ref = adj_raw

    # ---- cross-validation vs the field-standard mmseqs cluster ----
    rep_of = load_mmseqs_clusters()
    n_clusters = len(set(rep_of.values()))
    mean_pur, gw_pur, large_pur, n_large, n_in_one, ari = crossval(famP, genes, rep_of)
    print(f"[xval] mmseqs easy-cluster: {n_clusters} clusters over {len(rep_of)} proteins")
    print(f"       dominant purity: per-family-mean={mean_pur:.3f}  gene-weighted={gw_pur:.3f}  "
          f"large(>=30, n={n_large}) mean={large_pur:.3f}")
    print(f"       {n_in_one}/{len(famP)} families inside one cluster; ARI={ari:.3f}")

    # ---- lattice overlay: E_a (segdup) + E_c (read-conflict) on the SAME node set ----
    stats_lat = {}
    if not args.no_lattice:
        print(f"[E_a] one-pass segdup overlay {SEDEF} ...")
        ea_edges, ea_bailey = load_ea_edges(genes_by, gstarts)
        ea_fams = build_raw_families(ea_edges)
        ec_edges = load_ec_edges(name2idx)
        pc = lambda i: genes[i]["biotype"] == "protein_coding"
        ea_pc = {e for e in ea_edges if pc(e[0]) and pc(e[1])}
        print(f"       E_a: {len(ea_edges)} edges ({len(ea_pc)} protein-coding-both, {len(ea_bailey)} Bailey>=90%), "
              f"{len(ea_fams)} raw families; E_c: {len(ec_edges)} read-conflict edges")

        # EDGE-LEVEL lattice overlaps (the reviewer-facing |E_p only| / |E_a only| / |shared|)
        shared_pa = ep_edges & ea_edges
        ep_only = ep_edges - ea_edges
        ea_only = ea_edges - ep_edges
        shared_pa_pc = ep_edges & ea_pc
        ep_idx = fam_index(famP)
        ea_idx = fam_index(ea_fams)

        # E_c subset test: E_c edge co-protein-family (in an E_p family)?  and directly an E_p edge?
        ec_in_ep_fam = sum(1 for (u, v) in ec_edges if same_fam(ep_idx, u, v))
        ec_is_ep_edge = sum(1 for e in ec_edges if e in ep_edges)

        # FAMILY-LEVEL operational-subset: canonical refined DNA families inside one E_p family
        grf_tot, grf_in_one, grf_multi = grfam_containment(famP, genes)

        # --- WITNESS 1: ancient paralogs in E_p but NOT E_a  (globins across 5 chromosomes) ---
        globin_syms = ["HBZ", "HBM", "HBQ1", "HBE1", "MB", "CYGB", "NGB"]
        gidx = {s: name2idx[s] for s in globin_syms if s in name2idx}
        gfam = {s: ep_idx.get(i) for s, i in gidx.items()}
        glob_pairs = list(itertools.combinations(sorted(gidx), 2))
        glob_cofam = [(a, b) for a, b in glob_pairs if same_fam(ep_idx, gidx[a], gidx[b])]
        glob_ep_edge = [(a, b) for a, b in glob_pairs
                        if ((gidx[a], gidx[b]) if gidx[a] < gidx[b] else (gidx[b], gidx[a])) in ep_edges]
        glob_ea = [(a, b) for a, b in glob_pairs
                   if ((gidx[a], gidx[b]) if gidx[a] < gidx[b] else (gidx[b], gidx[a])) in ea_edges]
        glob_bailey = [(a, b) for a, b in glob_pairs
                       if ((gidx[a], gidx[b]) if gidx[a] < gidx[b] else (gidx[b], gidx[a])) in ea_bailey]
        glob_chroms = sorted({genes[i]["chrom"] for i in gidx.values()})

        # --- WITNESS 2: the 22q11 non-coding-segdup blob does NOT form in E_p ---
        blob_syms = ["DGCR6", "DRD5", "SLC6A8", "PRODH", "OTOP1", "BCAP31"]
        bidx = {s: name2idx[s] for s in blob_syms if s in name2idx}
        bfam = {s: ep_idx.get(i) for s, i in bidx.items()}
        blob_pairs = list(itertools.combinations(sorted(bidx), 2))
        blob_ep = [(a, b) for a, b in blob_pairs if same_fam(ep_idx, bidx[a], bidx[b])]
        blob_ea = [(a, b) for a, b in blob_pairs if same_fam(ea_idx, bidx[a], bidx[b])]
        n_distinct_ep_fams = len({v for v in bfam.values() if v is not None})

        print("\n== EDGE-LEVEL LATTICE OVERLAPS ==")
        print(f"  |E_p only| = {len(ep_only)}   |E_a only| = {len(ea_only)}   |E_p & E_a shared| = {len(shared_pa)} "
              f"(shared on protein-coding-both E_a: {len(shared_pa_pc)}/{len(ea_pc)})")
        # E_c and E_p are INCOMPARABLE in general (E_c\E_p via non-coding/UTR/retrocopy read-conflicts;
        # E_p\E_c via divergent non-read-confusable paralogs). The clean claim is CONDITIONAL: an E_c edge
        # whose shared stretch is coding at BOTH endpoints lies in E_p. Both available real E_c edges are
        # coding-both paralogs, so they confirm ONLY that coding-both conditional (not general containment).
        print(f"  E_c (coding-both) inside E_p: {ec_in_ep_fam}/{len(ec_edges)} read-conflict edges are "
              f"co-protein-family ({ec_is_ep_edge}/{len(ec_edges)} direct E_p edges); both are coding paralogs, "
              f"so this checks the coding-both CONDITIONAL, not unconditional E_c subset of E_p")
        if grf_tot:
            print(f"  family-level: {grf_in_one}/{grf_tot} refined DNA (GRFAM) families fully inside one E_p family "
                  f"({100*grf_in_one/grf_tot:.1f}%); {grf_multi} span >1 (segdup co-locating distinct protein families)")
        print("== WITNESS 1  E_p \\ E_a (globins) ==")
        print(f"  globins {sorted(gidx)} over chroms {glob_chroms}; E_p family ids {gfam}")
        print(f"  co-family pairs (transitive) {len(glob_cofam)}/{len(glob_pairs)}; direct E_p edges {len(glob_ep_edge)}; "
              f"E_a segdup pairs {len(glob_ea)}; Bailey>=90% segdup pairs {len(glob_bailey)}")
        print("== WITNESS 2  blob exclusion (22q11 DGCR6/DRD5/SLC6A8) ==")
        print(f"  E_p family ids {bfam}  -> {n_distinct_ep_fams} DISTINCT E_p families")
        print(f"  blob pairs in ONE E_p family: {len(blob_ep)}; in ONE E_a family (the GFAM0 blob): {len(blob_ea)}")

        stats_lat = dict(
            ea_edges=len(ea_edges), ea_pc_edges=len(ea_pc), ea_bailey_edges=len(ea_bailey), ea_families=len(ea_fams),
            ec_edges=len(ec_edges), ec_in_ep_family=ec_in_ep_fam, ec_is_ep_edge=ec_is_ep_edge,
            grfam_total=grf_tot, grfam_in_one_ep_family=grf_in_one, grfam_span_multi=grf_multi,
            lattice_ep_only=len(ep_only), lattice_ea_only=len(ea_only), lattice_shared=len(shared_pa),
            lattice_shared_pc=len(shared_pa_pc),
            witness_globin_chroms=glob_chroms, witness_globin_cofam_pairs=len(glob_cofam),
            witness_globin_ep_edges=len(glob_ep_edge), witness_globin_ea_pairs=len(glob_ea),
            witness_globin_bailey_pairs=len(glob_bailey), globin_ep_family_ids=gfam,
            witness_blob_distinct_ep_families=n_distinct_ep_fams, witness_blob_in_one_ep_family=len(blob_ep),
            witness_blob_in_one_ea_family=len(blob_ea), blob_ep_family_ids=bfam,
        )

    # ---- write catalogs (raw always; refined additionally under --refine) ----
    n_cert_raw = write_catalog(os.path.join(bench, "protein_families.tsv"), famP_raw, genes, adj_raw,
                               rep_of, args.gamma, "PFAM")
    outs = [os.path.join(bench, "protein_families.tsv")]
    n_cert_ref = None
    if args.refine:
        n_cert_ref = write_catalog(os.path.join(bench, "protein_families_refined.tsv"), famP, genes, adj_ref,
                                   rep_of, args.gamma, "PRFAM")
        outs.append(os.path.join(bench, "protein_families_refined.tsv"))

    stats = dict(
        gene_nodes=len(genes), proteins=len(syms), proteins_mapped=mapped,
        alpha_p=args.evalue, c_p=args.cov, gamma=args.gamma, seed=args.seed, refined=args.refine,
        expected_chance_directional_hits=round(exp_false, 4),
        ep_edges=len(ep_edges), ep_oneway_dropped=n_oneway,
        ep_raw_families=len(famP_raw), ep_raw_max=len(famP_raw[0]) if famP_raw else 0,
        ep_raw_cert=n_cert_raw,
        ep_families=len(famP), ep_max=len(famP[0]) if famP else 0,
        ep_member_genes=sum(len(m) for m in famP),
        ep_cross_chrom=sum(1 for m in famP if len({genes[i]["chrom"] for i in m}) >= 2),
        ep_refined_cert=n_cert_ref,
        ep_repeat_mosaic=sum(1 for m in famP if biotype_purity(m, genes) < PURITY),
        mmseqs_clusters=n_clusters, xval_mean_dominant_purity=round(mean_pur, 4),
        xval_gene_weighted_purity=round(gw_pur, 4),
        xval_large_family_purity=round(large_pur, 4), xval_large_family_count=n_large,
        xval_families_in_one_cluster=n_in_one, xval_adjusted_rand=round(ari, 4),
    )
    stats.update(stats_lat)
    json.dump(stats, open(os.path.join(bench, "protein_family_def.json"), "w"), indent=2, default=str)
    print(f"\n[wrote] " + "\n        ".join(outs) + f"\n        {os.path.join(bench, 'protein_family_def.json')}")


if __name__ == "__main__":
    main()
