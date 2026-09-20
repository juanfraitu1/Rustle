#!/usr/bin/env python
"""GENOME-grounded multi-copy gene family definition on the haploid T2T genome (GCF_029281585.2),
computed IN PARALLEL across the genome. Read-independent by construction: nodes come from the NCBI
RefSeq annotation and edges come from the SEDEF genome self-alignment (final.bed) — NEITHER input
ever sees an RNA read. This is the independent oracle the adversarial reviews kept asking for; the
RNA read-conflict family catalog can be overlaid on it without circularity.

It mirrors the RNA read-conflict definition (loci = nodes, a homology tie = an edge, family =
connected component with |C| >= 2), but with genome-only inputs, and builds BOTH:

  CATALOG A  (annotated-gene families)
    nodes = NCBI genes; edge (u,v) = a SINGLE SEDEF segdup pair covers gene u on one side and gene v
    on the other (each gene >=50% covered by its segdup side — the standard "duplicated gene" rule,
    not a tuned knob); family = connected component, |C| >= 2 genes.

  CATALOG B  (de-novo duplication-block families, annotation-free)
    nodes = segdup-footprint UNITS (overlapping segdup intervals merged per contig); edge = a SEDEF
    pair links two units; family = connected component of >= 2 units. Each block is then LABELED by
    the genes it covers (>=50%) — so a block is a gene set derived without using the annotation to
    build the graph.

PARALLEL: edge emission (A) and unit construction (B) are sharded by contig over a multiprocessing
pool (the genome analog of copy_assign --region-threads); a single global union-find merges the
per-contig results, so cross-chromosome paralog pairs are first-class.

Cross-checks (no extra edges): Compara paralogy (compara_paralog_relation.json) as an A-family
column; A-vs-B granularity agreement.

Run: /home/juanfra/miniforge3/bin/python bench/genome_family_def.py [--threads N]
"""
import argparse
import json
import os
import sys
from bisect import bisect_left
from collections import defaultdict
from multiprocessing import Pool

GFF = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)
COMPARA = os.path.join(os.path.dirname(__file__), "compara_paralog_relation.json")
RNA_FAM = os.path.join(os.path.dirname(__file__), "denovo_families.tsv")
COVER = 0.50   # a gene/locus is "in" a segdup side iff >=50% of it is covered (standard segdup-gene rule)
GAMMA = 0.20   # COHESION PARAMETER (irreducible, NOT threshold-free): a refined DNA family must be a
               # gamma-quasi-clique (internal edge density >= GAMMA). Chosen in the empirical GAP between the
               # verified-real-family density band [ZNF716 0.261 .. MAGEA/CSAG 0.733] and the transitive
               # over-merge blob (density 0.020). KNIFE-EDGE: the lowest real family (ZNF716, 0.261) survives
               # only for GAMMA < 0.261, so GAMMA=0.20 keeps a 0.061 margin and GAMMA>=0.27 splits the
               # KRAB-ZNF family -- a CALIBRATED, not free, choice. See DNA_FAMILY_DEFINITION_FORMAL.md sec 3/6.
SEED = 0       # SPLITTER SEED (disclosed constant, load-bearing for the WITNESS partition, not the object):
               # the exact-max-gamma-quasi-clique partition is NP-hard, so the splitter sigma yields ONE
               # certificate-passing witness among many. seed=0 fixes a deterministic, byte-identical witness;
               # the certificate PROPERTY (every family a gamma-quasi-clique with >=2 loci) is seed-invariant,
               # but the specific blob-derived family SET/COUNT is not (see doc sec 3.2/6).
LOCUS_OVERLAP = 0.50  # two member genes are the SAME physical locus iff same contig & reciprocal span
                      # overlap >= this (collapses redundant/nested annotations before the >=2-loci test)
PURITY = 0.50  # a genuine multi-copy gene family is copies of ONE gene => predominantly ONE biotype; a
               # certificate-passing family whose dominant-biotype fraction is < PURITY is DESCRIPTIVELY
               # flagged repeat_mosaic (necessary-not-sufficient: the gamma-quasi-clique certificate admits
               # repeat-mosaic hubs). Flag is a REPORTED column + sort-to-end; it does NOT filter families.

# ---- module globals shared with pool workers (set once per worker via initializer) ----
_GENES = None        # contig -> list[(start, end, gene_idx)] sorted by start
_GSTARTS = None      # contig -> [start, ...] parallel array for bisect


# ----------------------------- inputs -----------------------------
def load_genes(gff):
    """gene features -> per-contig sorted interval lists + flat metadata. 0-based half-open coords."""
    genes = []                       # idx -> dict(chrom,start,end,strand,name,biotype)
    by = defaultdict(list)
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]) - 1, int(f[4]), f[6], f[8]
            name, biotype = None, None
            for kv in attrs.split(";"):
                if kv.startswith("gene="):
                    name = kv[5:]
                elif kv.startswith("gene_biotype="):
                    biotype = kv[13:]
            idx = len(genes)
            genes.append(dict(chrom=chrom, start=start, end=end, strand=strand,
                              name=name or f"g{idx}", biotype=biotype or "?"))
            by[chrom].append((start, end, idx))
    for c in by:
        by[c].sort()
    return genes, by


def load_sedef_pairs(path, bailey=False):
    """final.bed -> (list of (cA,sA,eA,cB,sB,eB), n_total, n_kept). Each line is one SEDEF pair.

    DEFAULT (bailey=False) = the raw ~50%-floor SEDEF superset, chosen for the family-ORACLE purpose:
    HIGH RECALL of divergent paralog families (CEACAM5/6/7, KRAB-ZNF, PRSS1/2/3, IFITM, ULBP all fall
    BELOW the 90% genomic-segdup cliff and are LOST under SD(.)). The repeat-array over-merge it admits is
    cosmetic noise filtered per-component downstream; the divergent-family loss is a real recall hit.

    bailey=True applies the Bailey-2002 segmental-duplication predicate SD(.) -- a NARROWER 'recent segdup'
    oracle (see bench/SEGDUP_DEFINITION_FORMAL.md Check F): fracMatch (col 21) >= 0.90 AND aln_len (col 12)
    >= 1000 AND not a high-copy interspersed repeat (TE-exclusion: uppercaseMatches/aln_matches, col 28/29,
    >= 0.50). Keeps 27,623/253,029 (10.9%); removes the repeat-BRIDGED over-merge but NOT the transitive
    single-linkage chaining (max family still 317), and drops the divergent families above. Use for the
    'what is a recent high-identity segdup' view, not as the default family oracle."""
    pairs = []
    n_total = 0
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                rec = (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]))
            except ValueError:
                continue
            n_total += 1
            if bailey:
                try:                                  # SD(.): identity, length, TE-exclusion
                    aln_len, frac, um, am = int(f[11]), float(f[20]), int(f[27]), int(f[28])
                except (ValueError, IndexError):
                    continue
                if not (frac >= 0.90 and aln_len >= 1000 and am > 0 and um / am >= 0.50):
                    continue
            pairs.append(rec)
    return pairs, n_total, len(pairs)


# ----------------------------- overlap query -----------------------------
def _covered_genes(chrom, s, e):
    """gene indices on `chrom` whose span is >=COVER covered by [s,e]."""
    ivs = _GENES.get(chrom)
    if not ivs:
        return []
    starts = _GSTARTS[chrom]
    # candidates: genes starting before e; walk left while they can still overlap
    hi = bisect_left(starts, e)
    out = []
    for j in range(hi - 1, -1, -1):
        gs, ge, gi = ivs[j]
        if ge <= s:
            if s - gs > 5_000_000:      # start-sorted: stop after a long non-overlapping run
                break
            continue
        ov = min(ge, e) - max(gs, s)
        if ov > 0 and ov >= COVER * (ge - gs):
            out.append(gi)
    return out


def _init_worker(genes_by, gstarts):
    global _GENES, _GSTARTS
    _GENES, _GSTARTS = genes_by, gstarts


# ----------------------------- CATALOG A: parallel edge emission -----------------------------
def _edges_for_contig(args):
    """all SEDEF pairs anchored on side-A contig `c` -> gene-gene edges (side B may be any contig)."""
    c, pairs = args
    edges = []
    for (cA, sA, eA, cB, sB, eB) in pairs:
        ga = _covered_genes(cA, sA, eA)
        if not ga:
            continue
        gb = _covered_genes(cB, sB, eB)
        if not gb:
            continue
        for u in ga:
            for v in gb:
                if u != v:
                    edges.append((u, v) if u < v else (v, u))
    return edges


# ----------------------------- union-find -----------------------------
class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        r = x
        while self.p[r] != r:
            r = self.p[r]
        while self.p[x] != r:
            self.p[x], x = r, self.p[x]
        return r

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb

    def groups(self):
        g = defaultdict(list)
        for x in self.p:
            g[self.find(x)].append(x)
        return g


# --------------------- REFINEMENT operator R (cohesive-community, sec 3 of the formal note) ---------------------
def _induced_density(sub):
    """internal edge density rho_in(C) = 2|E(C)| / (|C|(|C|-1)); a singleton/empty = 1.0 (trivially cohesive)."""
    n = sub.number_of_nodes()
    return (2 * sub.number_of_edges() / (n * (n - 1))) if n > 1 else 1.0


def _split_once(sub, louvain, seed):
    """GUARANTEED-PROGRESS splitter sigma: adaptive-resolution modularity (Louvain), Kernighan-Lin bisection
    fallback, then a deterministic halving. Never returns a single block equal to the input on |C|>2,
    so the recursion in _refine_component provably terminates with every leaf gamma-cohesive or |C|<=2.
    `seed` is a DISCLOSED constant: exact max-gamma-quasi-clique partition is NP-hard, so sigma emits ONE
    certificate-passing witness; the certificate PROPERTY is seed-invariant, the blob family SET is not."""
    import networkx as nx
    n = sub.number_of_nodes()
    for res in (1.0, 2.0, 4.0, 8.0):                       # bump resolution until it splits
        parts = louvain(sub, resolution=res, seed=seed)
        if len(parts) >= 2:
            return parts
    if not nx.is_connected(sub):                           # disconnected -> split by components
        return [set(c) for c in nx.connected_components(sub)]
    try:                                                   # connected & modularity-stuck -> KL bisection
        a, b = nx.algorithms.community.kernighan_lin_bisection(sub, seed=seed)
        if a and b:
            return [a, b]
    except Exception:
        pass
    nodes = list(sub.nodes())                              # last resort: deterministic halving
    return [set(nodes[:n // 2]), set(nodes[n // 2:])]


def _refine_component(nodes, adj, gamma, louvain, seed):
    """R applied to ONE raw component's node set. Returns a list of node-set blocks (a PARTITION of `nodes`).
    Density-gated recursion: keep a block whole iff it is already a gamma-quasi-clique (or |C|<=2), else
    split and recurse. See formal note sec 3: idempotent, blocks subset the raw component, monotone in gamma."""
    import networkx as nx
    sub = nx.Graph()
    sub.add_nodes_from(nodes)
    for u in nodes:
        for v in adj.get(u, ()):                           # induced subgraph edges only
            if v in nodes and u < v:
                sub.add_edge(u, v)
    n = sub.number_of_nodes()
    if n <= 2 or _induced_density(sub) >= gamma:
        return [set(nodes)]
    out = []
    for p in _split_once(sub, louvain, seed):
        if len(p) == n:                                    # no-progress guard (cannot happen w/ guaranteed splitter)
            return [set(nodes)]
        out.extend(_refine_component(p, adj, gamma, louvain, seed))
    return out


def distinct_loci(members, genes):
    """Number of DISTINCT physical loci among member gene indices: genes on the same contig whose spans
    reciprocally overlap >= LOCUS_OVERLAP are collapsed to one locus (removes redundant/nested annotations,
    e.g. MAGEA9 == MAGd0b). The >=2-distinct-loci multi-copy predicate is |distinct_loci(C)| >= 2."""
    uf = UF()
    for i in members:
        uf.find(i)                                         # ensure singleton present
    by_c = defaultdict(list)
    for i in members:
        by_c[genes[i]["chrom"]].append(i)
    for c, idxs in by_c.items():
        idxs.sort(key=lambda i: genes[i]["start"])
        for a in range(len(idxs)):
            ia = idxs[a]; sa, ea = genes[ia]["start"], genes[ia]["end"]
            for b in range(a + 1, len(idxs)):
                ib = idxs[b]; sb, eb = genes[ib]["start"], genes[ib]["end"]
                if sb >= ea:                               # start-sorted: no further overlap on this contig
                    break
                ov = min(ea, eb) - max(sa, sb)
                if ov > 0 and ov >= LOCUS_OVERLAP * (ea - sa) and ov >= LOCUS_OVERLAP * (eb - sb):
                    uf.union(ia, ib)                       # same physical locus
    return len(uf.groups())


def biotype_purity(members, genes):
    """dominant-biotype fraction of a family: max_b |{i in C : biotype(i)=b}| / |C|. A genuine multi-copy
    family (copies of one gene) is biotype-pure; a repeat-mosaic hub mixes lncRNA/protein_coding/snRNA."""
    from collections import Counter
    if not members:
        return 1.0
    c = Counter(genes[i]["biotype"] for i in members)
    return c.most_common(1)[0][1] / len(members)


def refine_families(fams, all_edges, genes, gamma, seed=SEED):
    """R over a whole raw catalog: split each raw family (connected component) into cohesive communities,
    then apply the >=2-distinct-loci multi-copy predicate. Returns the refined catalog (list of gene-idx lists).
    `seed` fixes the splitter WITNESS (disclosed constant; the certificate is seed-invariant, the blob
    family set is not). Certificate-passing repeat_mosaic hubs (biotype purity < PURITY) are DESCRIPTIVELY
    sorted to the END of the catalog -- they are reported, not filtered (necessary-not-sufficient certificate)."""
    from networkx.algorithms.community import louvain_communities
    adj = defaultdict(set)
    fam_of = {}
    for k, m in enumerate(fams):
        for i in m:
            fam_of[i] = k
    for u, v in all_edges:                                 # adjacency restricted to nodes that are in some family
        if u in fam_of and v in fam_of and fam_of[u] == fam_of[v]:
            adj[u].add(v); adj[v].add(u)
    refined = []
    for m in fams:
        for block in _refine_component(set(m), adj, gamma, louvain_communities, seed):
            if distinct_loci(block, genes) >= 2:           # MULTI-COPY predicate: >=2 distinct loci
                refined.append(sorted(block))
    # repeat_mosaic hubs (purity < PURITY) sort to the END; genuine families first, largest-first within each.
    refined.sort(key=lambda m: (biotype_purity(m, genes) < PURITY, -len(m), genes[m[0]]["name"]))
    return refined


# ----------------------------- CATALOG B: de-novo duplication blocks -----------------------------
def _units_for_contig(args):
    """merge all segdup footprints on a contig into non-overlapping UNITS -> [(start,end), ...]."""
    c, ivals = args
    ivals.sort()
    units = []
    cs = ce = None
    for s, e in ivals:
        if cs is None:
            cs, ce = s, e
        elif s <= ce:
            ce = max(ce, e)
        else:
            units.append((cs, ce))
            cs, ce = s, e
    if cs is not None:
        units.append((cs, ce))
    return c, units


# ----------------------------- main -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=min(16, (os.cpu_count() or 4)))
    ap.add_argument("--bailey-sedef", action="store_true",
                    help="apply the Bailey SD(.) filter (>=90%% id, >=1kb, non-TE) instead of the default "
                         "raw ~50%%-floor superset. NARROWER 'recent segdup' oracle: cleaner repeat-wise but "
                         "DROPS divergent paralog families (CEACAM/KRAB-ZNF/PRSS fall below the 90%% cliff).")
    ap.add_argument("--refine", action="store_true",
                    help="apply the cohesive-community REFINEMENT R (density-gated recursive modularity split "
                         "with a gamma-quasi-clique stop) to each raw single-linkage component, and ADDITIONALLY "
                         "write the CANONICAL refined catalog genome_families_refined.tsv (prefix GRFAM). "
                         "OPT-IN: default (flag absent) is BYTE-IDENTICAL to the raw single-linkage catalog "
                         "(the default *.tsv/.json are untouched; refined output goes to its own file) -- mirrors "
                         "the RNA engine's --refine byte-identical-when-OFF discipline. "
                         "See bench/DNA_FAMILY_DEFINITION_FORMAL.md.")
    ap.add_argument("--gamma", type=float, default=GAMMA,
                    help=f"cohesion parameter for --refine: min internal edge density of a refined family "
                         f"(default {GAMMA}). IRREDUCIBLE (not threshold-free) -- honestly named per the formal note.")
    ap.add_argument("--seed", type=int, default=SEED,
                    help=f"splitter (sigma) seed for --refine (default {SEED}). DISCLOSED constant: the exact "
                         f"max-gamma-quasi-clique partition is NP-hard, so sigma emits ONE certificate-passing "
                         f"WITNESS; the certificate PROPERTY is seed-invariant but the blob-derived family "
                         f"SET/COUNT is not (~39-48 families from the 1547-gene blob across seeds).")
    args = ap.parse_args()
    assert SEDEF, "final.bed (SEDEF self-alignment) not found"
    bench = os.path.dirname(__file__)

    print(f"[load] genes  {GFF}")
    genes, genes_by = load_genes(GFF)
    gstarts = {c: [x[0] for x in v] for c, v in genes_by.items()}
    print(f"       {len(genes)} genes over {len(genes_by)} contigs")
    print(f"[load] sedef  {SEDEF}  (filter: {'Bailey SD(.)' if args.bailey_sedef else 'RAW superset (default, high-recall oracle)'})")
    pairs, n_total, n_kept = load_sedef_pairs(SEDEF, bailey=args.bailey_sedef)
    print(f"       {n_kept} segdup pairs kept of {n_total} "
          f"({100*n_kept/max(1,n_total):.1f}%{' pass >=90%/>=1kb/non-TE' if args.bailey_sedef else ' (raw superset)'})")

    # shard pairs by side-A contig for parallel edge emission
    pairs_by_cA = defaultdict(list)
    for p in pairs:
        pairs_by_cA[p[0]].append(p)
    shards = sorted(pairs_by_cA.items(), key=lambda kv: -len(kv[1]))
    print(f"[A] emitting gene-gene edges over {len(shards)} contig shards x {args.threads} threads ...")
    with Pool(args.threads, initializer=_init_worker, initargs=(genes_by, gstarts)) as pool:
        edge_lists = pool.map(_edges_for_contig, shards)
    all_edges = [(u, v) for el in edge_lists for (u, v) in el]
    n_edges = len(all_edges)

    def build(edges, keep=None):
        uf = UF()
        for u, v in edges:
            if keep is None or (keep(u) and keep(v)):
                uf.union(u, v)
        fams = [sorted(set(m)) for m in uf.groups().values() if len(set(m)) >= 2]
        fams.sort(key=lambda m: (-len(m), genes[m[0]]["name"]))
        return fams

    famA = build(all_edges)                                          # all biotypes (RAW single-linkage components)
    pc = lambda i: genes[i]["biotype"] == "protein_coding"
    famPC = build(all_edges, keep=pc)                                # protein-coding only (raw)
    print(f"    {n_edges} edges -> {len(famA)} RAW families all-biotype ({sum(len(m) for m in famA)} genes); "
          f"{len(famPC)} protein-coding raw families ({sum(len(m) for m in famPC)} genes)")

    # ---- REFINEMENT R (OPT-IN --refine): raw single-linkage components -> cohesive communities w/ >=2 loci ----
    # Default (flag absent) => famR stays None => refined.tsv NOT written and the json omits refined_* keys,
    # so the default catalog files are BYTE-IDENTICAL to the raw single-linkage output.
    famR = None
    if args.refine:
        print(f"[R] refining {len(famA)} raw components into cohesive communities "
              f"(gamma={args.gamma}, seed={args.seed}, >=2 distinct loci) ...")
        famR = refine_families(famA, all_edges, genes, args.gamma, seed=args.seed)
        # verify the certificate + idempotency the formal note claims (cheap, on the produced catalog)
        try:
            import networkx as _nx
            fam_of = {i: k for k, m in enumerate(famR) for i in m}
            adjR = defaultdict(set)
            for u, v in all_edges:
                if fam_of.get(u) is not None and fam_of.get(u) == fam_of.get(v):
                    adjR[u].add(v); adjR[v].add(u)
            def _dens(m):
                e = sum(1 for u in m for v in adjR.get(u, ()) if v in set(m) and u < v)
                n = len(m); return (2 * e / (n * (n - 1))) if n > 1 else 1.0
            cohesive = sum(1 for m in famR if len(m) <= 2 or _dens(m) >= args.gamma)
            print(f"    {len(famA)} raw -> {len(famR)} refined families "
                  f"({sum(len(m) for m in famR)} genes); certificate: {cohesive}/{len(famR)} "
                  f"gamma-cohesive-or-size<=2 (should be all)")
        except Exception as e:
            print(f"    refined to {len(famR)} families (certificate check skipped: {e})")

    def dom_biotype(m):
        from collections import Counter
        return Counter(genes[i]["biotype"] for i in m).most_common(1)[0][0]

    # ---- CATALOG B: duplication-block families (annotation-free graph) ----
    print(f"[B] building duplication-block families ...")
    foot_by_c = defaultdict(list)            # contig -> all segdup footprints (both sides)
    for (cA, sA, eA, cB, sB, eB) in pairs:
        foot_by_c[cA].append((sA, eA)); foot_by_c[cB].append((sB, eB))
    with Pool(args.threads) as pool:
        unit_res = pool.map(_units_for_contig, list(foot_by_c.items()))
    # global unit ids + per-contig sorted unit starts for lookup
    unit_id = {}                              # (contig, start, end) -> uid
    units_by_c = {}                           # contig -> ([start...], [(start,end,uid)...])
    uid = 0
    for c, units in unit_res:
        rows = []
        for (s, e) in units:
            unit_id[(c, s, e)] = uid
            rows.append((s, e, uid)); uid += 1
        units_by_c[c] = ([r[0] for r in rows], rows)

    def unit_of(c, s, e):
        info = units_by_c.get(c)
        if not info:
            return None
        starts, rows = info
        j = bisect_left(starts, e) - 1       # unit covering [s,e]; units are disjoint & sorted
        while j >= 0:
            us, ue, u = rows[j]
            if ue <= s:
                break
            if us <= s and e <= ue:
                return u
            j -= 1
        return None

    ufB = UF()
    for (cA, sA, eA, cB, sB, eB) in pairs:
        ua, ub = unit_of(cA, sA, eA), unit_of(cB, sB, eB)
        if ua is not None and ub is not None and ua != ub:
            ufB.union(ua, ub)
    # label each block-component by the genes covering its units; keep blocks with >=2 gene loci
    uid_to_iv = {u: (c, s, e) for (c, s, e), u in unit_id.items()}
    blocks = []
    for comp in ufB.groups().values():
        comp = set(comp)
        if len(comp) < 2:
            continue
        gset = set()
        contigs = set()
        for u in comp:
            c, s, e = uid_to_iv[u]; contigs.add(c)
            gset.update(_label_genes(genes_by, gstarts, c, s, e))
        if len(gset) >= 2:
            blocks.append(dict(n_units=len(comp), contigs=sorted(contigs),
                               genes=sorted(gset, key=lambda i: genes[i]["name"])))
    blocks.sort(key=lambda b: -len(b["genes"]))
    print(f"    {uid} units -> {len(blocks)} gene-bearing duplication blocks (>=2 units & >=2 genes)")

    # ---- cross-checks ----
    compara = json.load(open(COMPARA)) if os.path.exists(COMPARA) else {"families": {}}
    cpairs = set()
    for fam in compara.get("families", {}).values():
        mg = [g for g in fam.get("mapped_genes", [])]
        for i in range(len(mg)):
            for j in range(i + 1, len(mg)):
                cpairs.add(tuple(sorted((mg[i], mg[j]))))

    def named(idx):
        n = genes[idx]["name"]
        return None if n.startswith("LOC") or n.startswith("g") else n

    a_named = a_compara_hit = 0
    for m in famA:
        names = {named(i) for i in m} - {None}
        if len(names) >= 2:
            a_named += 1
            if any(tuple(sorted((x, y))) in cpairs
                   for x in names for y in names if x < y):
                a_compara_hit += 1

    # A-vs-B granularity: does each A-family's gene set sit inside ONE B-block?
    gene_block = {}
    for bi, b in enumerate(blocks):
        for gi in b["genes"]:
            gene_block.setdefault(gi, bi)
    a_in_one_block = sum(1 for m in famA
                         if len({gene_block.get(i) for i in m} - {None}) == 1
                         and any(i in gene_block for i in m))

    xchrom = sum(1 for m in famA if len({genes[i]["chrom"] for i in m}) >= 2)

    # ---- write catalogs ----
    def write_catalog(path, fams, prefix="GFAM"):
        with open(path, "w") as fh:
            fh.write("family_id\tn_genes\tn_contigs\tcross_chrom\tdom_biotype\tcompara_paralog\tgenes\n")
            for k, m in enumerate(fams):
                nn = {named(i) for i in m} - {None}
                comp_hit = any(tuple(sorted((x, y))) in cpairs for x in nn for y in nn if x < y)
                fh.write(f"{prefix}{k}\t{len(m)}\t{len({genes[i]['chrom'] for i in m})}\t"
                         f"{int(len({genes[i]['chrom'] for i in m})>=2)}\t{dom_biotype(m)}\t"
                         f"{int(comp_hit)}\t{','.join(genes[i]['name'] for i in m)}\n")

    outA = os.path.join(bench, "genome_families_annotated.tsv")
    write_catalog(outA, famA)
    outPC = os.path.join(bench, "genome_families_protein_coding.tsv")
    write_catalog(outPC, famPC)
    outR = None
    if famR is not None:                                            # the CANONICAL refined catalog
        adjR = defaultdict(set)
        fam_of = {i: k for k, m in enumerate(famR) for i in m}
        for u, v in all_edges:
            if fam_of.get(u) is not None and fam_of.get(u) == fam_of.get(v):
                adjR[u].add(v); adjR[v].add(u)
        def _dens(m):
            e = sum(1 for u in m for v in adjR.get(u, ()) if v in set(m) and u < v)
            n = len(m); return (2 * e / (n * (n - 1))) if n > 1 else 1.0
        outR = os.path.join(bench, "genome_families_refined.tsv")
        with open(outR, "w") as fh:
            fh.write("family_id\tn_genes\tn_loci\tn_contigs\tcross_chrom\tdensity\tbiotype_purity\t"
                     "repeat_mosaic\tdom_biotype\tcompara_paralog\tgenes\n")
            for k, m in enumerate(famR):
                nn = {named(i) for i in m} - {None}
                comp_hit = any(tuple(sorted((x, y))) in cpairs for x in nn for y in nn if x < y)
                pur = biotype_purity(m, genes)             # dominant-biotype fraction (coherence descriptor)
                mosaic = int(pur < PURITY)                 # necessary-not-sufficient: certificate-pass but mixed-biotype hub
                fh.write(f"GRFAM{k}\t{len(m)}\t{distinct_loci(m, genes)}\t{len({genes[i]['chrom'] for i in m})}\t"
                         f"{int(len({genes[i]['chrom'] for i in m})>=2)}\t{_dens(m):.3f}\t{pur:.3f}\t{mosaic}\t"
                         f"{dom_biotype(m)}\t{int(comp_hit)}\t{','.join(genes[i]['name'] for i in m)}\n")
    outB = os.path.join(bench, "genome_families_blocks.tsv")
    with open(outB, "w") as fh:
        fh.write("block_id\tn_units\tn_genes\tn_contigs\tgenes\n")
        for k, b in enumerate(blocks):
            fh.write(f"GBLK{k}\t{b['n_units']}\t{len(b['genes'])}\t{len(b['contigs'])}\t"
                     f"{','.join(genes[i]['name'] for i in b['genes'])}\n")

    xchrom_pc = sum(1 for m in famPC if len({genes[i]["chrom"] for i in m}) >= 2)
    stats = dict(genes=len(genes), segdup_pairs=len(pairs), edges=n_edges,
                 catalogA_families=len(famA), catalogA_member_genes=sum(len(m) for m in famA),
                 catalogA_cross_chrom=xchrom,
                 catalogPC_families=len(famPC), catalogPC_member_genes=sum(len(m) for m in famPC),
                 catalogPC_cross_chrom=xchrom_pc,
                 catalogB_blocks=len(blocks), units=uid,
                 A_named_ge2=a_named, A_compara_confirmed=a_compara_hit,
                 A_in_single_B_block=a_in_one_block, cover_frac=COVER)
    if famR is not None:
        stats.update(refined_families=len(famR), refined_member_genes=sum(len(m) for m in famR),
                     refined_cross_chrom=sum(1 for m in famR if len({genes[i]["chrom"] for i in m}) >= 2),
                     refined_repeat_mosaic=sum(1 for m in famR if biotype_purity(m, genes) < PURITY),
                     gamma=args.gamma, seed=args.seed, locus_overlap=LOCUS_OVERLAP, purity=PURITY)
    json.dump(stats, open(os.path.join(bench, "genome_family_def.json"), "w"), indent=2)

    def top(fams, n=6):
        return ", ".join(f"{genes[m[0]]['name']}({len(m)},{dom_biotype(m)})" for m in fams[:n])
    print("\n== GENOME-GROUNDED MULTI-COPY FAMILY DEFINITION (read-independent, parallel) ==")
    print(f"  RAW (all biotypes)   {len(famA)} families, {stats['catalogA_member_genes']} genes, "
          f"{xchrom} cross-chrom")
    print(f"      top: {top(famA)}")
    print(f"      -> dominated by tRNA/snRNA/lncRNA repeat arrays (RNA-gene multi-copy, repeat-bridged)")
    print(f"  RAW protein-coding   {len(famPC)} families, {stats['catalogPC_member_genes']} genes, "
          f"{xchrom_pc} cross-chrom")
    print(f"      top: {top(famPC)}")
    if famR is not None:
        print(f"  CANONICAL REFINED (--refine, seed={args.seed})  {len(famR)} cohesive families "
              f"(gamma={args.gamma}, >=2 loci), {stats['refined_member_genes']} genes, "
              f"{stats['refined_cross_chrom']} cross-chrom, {stats['refined_repeat_mosaic']} repeat_mosaic-flagged")
        print(f"      top: {top(famR)}")
    print(f"  CATALOG B (de-novo blocks)  {len(blocks)} duplication blocks (>=2 units & >=2 genes), {uid} units")
    print(f"  Compara cross-check  {a_compara_hit}/{a_named} multi-named families carry a Compara "
          f"paralog relation (sparse Compara coverage)")
    print(f"  A vs B granularity   {a_in_one_block}/{len(famA)} families inside one B-block "
          f"({100*a_in_one_block/max(1,len(famA)):.0f}%)")
    print(f"  wrote {outA}\n        {outPC}" + (f"\n        {outR}" if outR else "") +
          f"\n        {outB}\n        {os.path.join(bench,'genome_family_def.json')}")


def _label_genes(genes_by, gstarts, chrom, s, e):
    ivs = genes_by.get(chrom)
    if not ivs:
        return []
    starts = gstarts[chrom]
    hi = bisect_left(starts, e)
    out = []
    for j in range(hi - 1, -1, -1):
        gs, ge, gi = ivs[j]
        if ge <= s:
            if s - gs > 5_000_000:
                break
            continue
        ov = min(ge, e) - max(gs, s)
        if ov > 0 and ov >= COVER * (ge - gs):
            out.append(gi)
    return out


if __name__ == "__main__":
    main()
