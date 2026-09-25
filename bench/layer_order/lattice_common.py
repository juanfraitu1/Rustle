#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order + nested edge-test lattice (human CHM13, 2026-09-16): the shared library.

The CLI is `npip_tbc1d3.py` (one subcommand per stage). Wave 7 (2026-09-24) folded the 12 files of
bench/layer_order/ into this library and that CLI. The old files are at git tag `notebook-2026-09-24`
(`git show notebook-2026-09-24:bench/layer_order/<old>.py`). The module keeps the name `lattice_common`, and every
name it exported, because off-repo scripts import it by path: LAT/corrections_pass2/doc_numbers.py (OUT, components,
fnum) and strand_diag.py (OUT, tests, fnum), whose outputs docs/seeded_family_definition.md cites; and
LAT/audit_skeptic/*.py (tests, components, groups).

Old name -> new name (library; the stage commands are listed in npip_tbc1d3.py):
  lattice_common.*                                  -> same names (pairs_in dropped: no caller anywhere)
  soto_map.load / map_gene / FIELDS / parse_blocks / ov_bp
                                                    -> soto_load / soto_map_gene / SOTO_FIELDS / parse_blocks / ov_bp
  lo_analysis.c_tree                                -> c_tree(U, NAME)
  lo_analysis.bip_jaccard / f1 / score / score_member_anchored / fmt
                                                    -> bip_jaccard / f1 / score_lo / score_member_anchored / fmt
  lo_analysis.hgnc_all                              -> hgnc_all (the rule itself is hgnc_lookup; lo_corrected_tables
                                                       ran the same rule inline and also recorded how it matched)
  lo_analysis.soto_ok; the 'flag ok' rule inlined in lo_corrected_tables and lattice_edges
                                                    -> soto_ok / soto_flag
  lo_analysis.uf_components / _join_labels          -> uf_components / join_labels
  lo_corrected_tables.catalog_keys / membership / key_of / CATALOGS / LIT
                                                    -> same names (LIT now resolves to <repo>/docs/..., audit F1)
  lattice_edges.key_of / merge / exonic_bases_in / chain_groups / chain_stats
                                                    -> same names (merge is lib.merge; chain_eval dropped: no caller)
  lattice_edges head (genes, catalog keys, exon blocks; exec'd by lattice_check_c2)
                                                    -> catalog_context(say)
  lattice_truth.c2 / score                          -> c2 / score_counts
  lattice_levels.shortest_path = lattice_report_tables.bfs_path
                                                    -> shortest_path(adj, src, dst)
  lattice_expr / lo_expr_recount: GFF exon index, hits(), window BED, samtools loop
                                                    -> gff_exon_index / IntervalIndex.hits / write_windows_bed /
                                                       count_reads (one rule; the two scripts differed only in the
                                                       gene set, the 'gene-' prefix and unique_mr)
  denovo_shared_def.ExonIndex (archived 2026-09-19; lattice_check_c2 died on the import)
                                                    -> IntervalIndex (IntervalIndex.from_nodes = the old constructor)
  every per-script tsv / write / key_of / path block -> tsv / write (fmtv, full precision) / write_str (str(); the
                                                       lo_* scripts' writer) / key_of / paths below

Paths: every result directory hangs off ROOT = $LO_ROOT (default /mnt/linuxdisk/home/juanfraitu/layer_order/
npip_tbc1d3; `npip_tbc1d3.py --root DIR` calls set_root). The old scripts derived the repo root with one dirname too
few (`bench/`, not the repo; audit F1), so LIT and the bench imports never resolved from a clone; REPO below uses three.

Only the standard library is imported at module top (plus bench/lib.py, itself stdlib-only at top); numpy/scipy are
imported inside the scorers that use them.

------------------------------------------------------------------------------------------------ the lattice (verbatim)
Nested edge-test lattice (NPIP/TBC1D3, human CHM13) — shared paths, loaders and THE LEVEL TESTS.

User-approved design (2026-09-16 16:39): one graph (nodes = RefSeq gene records at annotated gene-body extent; every
candidate link carries all its evidence); level k = connected components of the edges passing t_k, with
t_{k+1} => t_k by construction (each test is the previous test AND one more clause):

  t0 superfamily            = not same_locus AND (protein §6ko edge OR t1)
  t1 family                 = not same_locus AND clause-2 DNA homology (APPROXIMATION from the E1 gene-body PAF; the hit
                              must overlap v's exons and pass the shipped strand check for spliced pairs, see
                              npip_tbc1d3.py lattice-edges)
  t2 shared-exon unit       = t1 AND shared-exon fraction >= 0.30  (mcl_families quantity: max over E1 records)
  t3 >= 0.98 identity unit  = t2 AND w_98 >= 0.98  (primary: the §0★★★.1 single-record w_98, gap-excluded identity)
(docs/seeded_family_definition.md §0★★★.1 names L2/L3 'shared-exon unit' and '>= 0.98 identity unit', not 'duplication
unit' / 'subfamily'.)

Correction pass 2026-09-16 (audit of the 17:12 run): every test thresholds UNROUNDED values (edges.tsv now stores floats
at full repr precision; the 17:03 build stored 4 dp, which admitted identities 0.97995-0.98 at L3). Pooled identities are
kept as variants and are NOT monotone (H3 fails for them).

A missing attribute makes a clause unsatisfiable (never imputed): no PAF record for the pair => t1..t3 false; no blastp
HSP pair => protein clause false; S2 variant: no S2 edge => false; no record witnessing sx >= 0.30 => w_98 is NA => t3 false.

same_locus: the two records' genomic intervals intersect on one chromosome (readthrough vs its component gene, nested
records). Their "alignment" is the identity map of shared bases, not homology between two copies, so the link is
excluded from every test (reported; variant 'with_same_locus' re-admits it).
"""
import bisect
import collections
import csv
import gzip
import itertools
import os
import re
import subprocess
import sys

# ------------------------------------------------------------------------------------------------ paths
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(REPO, "bench")
if BENCH not in sys.path:
    sys.path.insert(0, BENCH)
from lib import merge  # noqa: E402,F401  (the canonical interval union; lattice_edges and guided_pipeline had copies)

DEFAULT_ROOT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3"
H = "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge"
PAF = {"c15_17_22": f"{H}/human2/genes.asm20.paf", "c16_19_20": f"{H}/lit/aj_ho/refseq/all.paf"}
NODES_C16 = f"{H}/lit/aj_ho/refseq/nodes.tsv"
CAT_CHROMS = {"c15_17_22": {"chr15", "chr17", "chr22"}, "c16_19_20": {"chr16", "chr19", "chr20"}}
HGNC = "/mnt/linuxdisk/home/juanfraitu/winloci_data/hgnc/hgnc_complete_set.txt"
SOTO = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication"
S1C = os.path.join(REPO, "bench", "soto", "soto_famCN_S1C.tsv")
LIT = os.path.join(REPO, "docs", "lit_subclusters_npip_tbc1d3_truth.tsv")
BAM = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def set_root(root):
    """(Re)bind every ROOT-derived path: LIGHT (light/), HEAVY (heavy/), INT (integrate_slim/, 'IS/' in the reports),
    OUT (lattice/, 'LAT/'), the E1/S1 graph dumps and the two E1 catalogs."""
    global ROOT, LIGHT, HEAVY, INT, OUT, DUMP, DUMP_S1, CATALOGS
    ROOT = root
    LIGHT = f"{root}/light"
    HEAVY = f"{root}/heavy"
    INT = f"{root}/integrate_slim"
    OUT = f"{root}/lattice"
    DUMP = {"c15_17_22": f"{LIGHT}/work/D/c15_17_22.e1.graph.tsv", "c16_19_20": f"{LIGHT}/work/D/c16_19_20.e1.graph.tsv"}
    DUMP_S1 = {"c15_17_22": f"{LIGHT}/work/S1/c15_17_22.e1s.graph.tsv",
               "c16_19_20": f"{LIGHT}/work/S1/c16_19_20.e1s.graph.tsv"}
    CATALOGS = {
        "c15_17_22": dict(kind="regions", path=f"{H}/human2/genes.regions", D=f"{LIGHT}/work/D/c15_17_22.e1",
                          E0=f"{H}/human2/guided", E1=f"{H}/lit/aj_dev/refseq_e1"),
        "c16_19_20": dict(kind="nodes", path=f"{H}/lit/aj_ho/refseq/nodes.tsv", D=f"{LIGHT}/work/D/c16_19_20.e1",
                          E0=f"{H}/lit/aj_ho/refseq/e0", E1=f"{H}/lit/aj_ho/refseq/e1"),
    }


ROOT = LIGHT = HEAVY = INT = OUT = DUMP = DUMP_S1 = CATALOGS = None  # bound by set_root (next line and --root)
set_root(os.environ.get("LO_ROOT", DEFAULT_ROOT))

SEF_MIN = 0.30
ID_L3 = 0.98
LEVELS = ("L0", "L1", "L2", "L3")


# ------------------------------------------------------------------------------------------------ io
def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(path, rows, cols=None):
    """lattice writer: every value through fmtv (full precision; bool -> yes/no, None -> NA)."""
    cols = cols or (list(rows[0].keys()) if rows else ["empty"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmtv(r.get(c, "")) for c in cols) + "\n")


def write_str(path, rows, cols=None):
    """layer-order writer (lo_analysis.write; lo_corrected_tables.write with the argument order (path, cols, rows)):
    every value through str()."""
    cols = cols or (list(rows[0].keys()) if rows else ["empty"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def fmtv(x):
    """Full precision (repr round-trips exactly): tests and report tables must format from unrounded values."""
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        return "NA" if x != x else repr(x)
    if x is None:
        return "NA"
    return str(x)


def fnum(s):
    return None if s in ("", "NA", None) else float(s)


class Log:
    """The per-script say(): print a line and keep it for the stage's .out file."""

    def __init__(self, flush=True):
        self.lines = []
        self.flush = flush

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=self.flush)
        self.lines.append(s)

    def dump(self, path):
        with open(path, "w") as fh:
            fh.write("\n".join(self.lines) + "\n")


# ------------------------------------------------------------------------------------------------ the tests
def yes(s):
    return s == "yes"


def t_protein(r, p_aa_min=None):
    """§6ko edge as shipped (greedy HSP cover); p_aa_min strengthens the disjunct (T3(a): (t_P AND aa >= x) OR t_D)."""
    ok = yes(r["p_qualifies_6ko"])
    if ok and p_aa_min is not None:
        ok = (fnum(r["p_aa_identity"]) or 0.0) >= p_aa_min
    return ok


# L1 (clause-2) operationalisations: edges.tsv column per option
L1_COL = {
    "c2": "d_c2_approx",                # primary: hit overlaps v's exons + shipped strand check; shorter-body denominator
    "c2_nostrand": "d_c2nostrand_approx",  # hit overlaps v's exons, no strand check
    "c2_exontgt": "d_c2exontgt_approx",  # v-exon overlap required on the exon proxy only (gene-body chain unrestricted)
    "c2_loose": "d_c2loose_approx",     # no v-exon overlap, no strand check (the 17:03 build's primary)
    "c2x": "d_c2x_approx",              # guided finder's extrapolated-span denominator, v-exon overlap + strand check
}
# L3 identity fields
ID_COL = {
    "w98_gapexcl": "d_w98_gapexcl",      # primary: single-record w_98 (§0★★★.1), gap-excluded identity (SEDEF fracMatch analogue)
    "w98_gapincl": "d_w98_gapincl",      # single-record w_98, gap-inclusive identity (nm / block length)
    "pooled_gapexcl": "d_e1_identity_gapexcl",  # pooled over E1 records, gap-excluded (17:03 primary; NOT monotone)
    "pooled_gapincl": "d_e1_identity",   # pooled over E1 records, gap-inclusive (NOT monotone)
    "s2": "s2_max_identity",             # SEDEF SD98 map-back max identity
}


def t_dna_c2(r, l1="c2"):
    return yes(r[L1_COL[l1]])


def sef(r):
    return fnum(r["d_shared_exon_frac"])


def ident(r, which):
    return fnum(r[ID_COL[which]])


def tests(r, id_min=ID_L3, id_which="w98_gapexcl", sef_min=SEF_MIN, with_same_locus=False, l1="c2", p_aa_min=None):
    """(t0, t1, t2, t3) for one edge row of edges.tsv. l1: one of L1_COL (primary 'c2'), or 'e1' (E1 edge as built = the
    D layer's graph), 'e1c2' (E1 edge with identity >= 0.80 and cov_longer >= 0.50). id_which: one of ID_COL. All
    comparisons are on unrounded values."""
    if yes(r["same_locus"]) and not with_same_locus:
        return (False, False, False, False)
    if l1 in L1_COL:
        a1 = t_dna_c2(r, l1)
    elif l1 == "e1":
        a1 = yes(r["d_e1_edge"])
    elif l1 == "e1c2":
        a1 = yes(r["d_e1_edge"]) and (fnum(r["d_e1_identity"]) or 0) >= 0.80 and (fnum(r["d_e1_cov_longer"]) or 0) >= 0.50
    else:
        raise ValueError(l1)
    t1 = a1
    s = sef(r)
    t2 = t1 and s is not None and (sef_min is None or s >= sef_min)
    i = ident(r, id_which)
    t3 = t2 and i is not None and i >= id_min
    t0 = t1 or t_protein(r, p_aa_min)
    return (t0, t1, t2, t3)


# ------------------------------------------------------------------------------------------------ graph helpers
class UF:
    def __init__(self, nodes):
        self.p = {n: n for n in nodes}

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if ra < rb:
                ra, rb = rb, ra
            self.p[ra] = rb


def components(nodes, edges):
    """nodes: iterable; edges: iterable of (a, b). Returns {node: component representative (min node id)}."""
    uf = UF(nodes)
    for a, b in edges:
        if a in uf.p and b in uf.p:
            uf.union(a, b)
    lab = {n: uf.find(n) for n in uf.p}
    rep = {}
    for n, r in lab.items():
        rep[r] = min(rep.get(r, n), n)
    return {n: rep[r] for n, r in lab.items()}


def groups(lab):
    g = collections.defaultdict(set)
    for n, c in lab.items():
        g[c].add(n)
    return g


def truss3(edges):
    """3-truss edge set: drop every edge that lies in no triangle of the CURRENT edge set, until none is left to drop
    (fixed point; queue peeling with triangle-support counts). For k = 3 the fixed point is reached after one pass:
    deleting an edge that lies in no triangle cannot remove a triangle of any other edge, so no cascade occurs (the
    returned peel depth is therefore always 1 when anything is dropped; kept as a check).
    Returns (kept edges as sorted tuples, number of dropped edges, peel depth)."""
    E = {tuple(sorted(e)) for e in edges if e[0] != e[1]}
    adj = collections.defaultdict(set)
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    sup = {e: len(adj[e[0]] & adj[e[1]]) for e in E}
    depth = {e: 1 for e in E if sup[e] == 0}
    queue = collections.deque(sorted(depth))
    alive = set(E)
    dropped = 0
    while queue:
        e = queue.popleft()
        if e not in alive:
            continue
        a, b = e
        alive.discard(e)
        dropped += 1
        adj[a].discard(b)
        adj[b].discard(a)
        for w in adj[a] & adj[b]:
            for f in (tuple(sorted((a, w))), tuple(sorted((b, w)))):
                if f in alive:
                    sup[f] -= 1
                    if sup[f] == 0 and f not in depth:
                        depth[f] = depth[e] + 1
                        queue.append(f)
    return alive, dropped, (max(depth.values()) if depth else 0)


def refines(fine, coarse, nodes=None):
    """count fine groups (restricted to nodes) that span >1 coarse group."""
    nodes = set(fine) if nodes is None else set(nodes)
    viol = []
    for c, G in groups({n: fine[n] for n in nodes if n in fine}).items():
        if len({coarse[n] for n in G}) > 1:
            viol.append(G)
    return viol


def split_counts(fine, coarse, nodes, mem):
    """Non-vacuity of a refinement check: coarse blocks (restricted to nodes) with >= 2 genes, how many of them the fine
    partition actually splits, and how many of the split ones hold a member."""
    G = groups({n: coarse[n] for n in nodes})
    big = [S for S in G.values() if len(S) >= 2]
    split = [S for S in big if len({fine[n] for n in S}) > 1]
    return len(big), len(split), sum(1 for S in split if S & mem)


def shortest_path(adj, src, dst):
    """BFS with sorted neighbours (lattice_levels.shortest_path; lattice_report_tables.bfs_path was the same code with
    the arguments (src, dst, adj)). adj: node -> iterable of neighbours (a set, or a dict keyed by neighbour)."""
    prev = {src: None}
    q = collections.deque([src])
    while q:
        x = q.popleft()
        if x == dst:
            break
        for y in sorted(adj[x]):
            if y not in prev:
                prev[y] = x
                q.append(y)
    if dst not in prev:
        return None
    path = [dst]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return path[::-1]


def uf_components(nodes, links):
    """lo_analysis.uf_components (kept apart from UF/components: its component ORDER, i.e. the iteration order of
    `nodes`, decides tie order in IS/expr_groups.tsv and IS/expr_sweep.tsv)."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in links:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for n in nodes:
        comp[find(n)].add(n)
    return list(comp.values())


def join_labels(a, b):
    """partition join of two labelings on the same genes (union-find over shared labels). (lo_analysis._join_labels)"""
    parent = {g: g for g in a}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for lab in (a, b):
        first = {}
        for g in sorted(a):
            if lab[g] in first:
                parent[find(g)] = find(first[lab[g]])
            else:
                first[lab[g]] = g
    comp = collections.defaultdict(set)
    for g in a:
        comp[find(g)].add(g)
    out = {}
    for G in comp.values():
        names = sorted({lab[g] for lab in (a, b) for g in G if "singleton" not in lab[g] and "|" in lab[g]})
        if len(G) == 1:
            out[next(iter(G))] = a[next(iter(G))]
            continue
        for g in G:
            out[g] = "Cmid|" + ("+".join(x.split("|", 1)[1] for x in names) if names else "+".join(sorted(G)))
    return out


# ------------------------------------------------------------------------------------------------ intervals
def key_of(s):
    c, r = s.rsplit(":", 1)
    a, b = r.split("-")
    return (c, int(a), int(b))


def exonic_bases_in(blocks, gene_start1, s, e):
    off = gene_start1 - 1
    acc = 0
    for bs, be in blocks:
        lo, hi = max(bs - off, s), min(be - off, e)
        if hi > lo:
            acc += hi - lo
    return acc


class IntervalIndex:
    """Half-open intervals per chromosome; hits(chrom, s, e) = labels of the intervals overlapping [s, e) by >= 1 bp.
    The archived denovo_shared_def.ExonIndex (tag notebook-2026-09-19), generalised to (chrom, start, end, label) items;
    lattice_expr.py and lo_expr_recount.py built the same bisect index inline (sorted (start, end, label) per chromosome,
    back-off by the longest interval)."""

    def __init__(self, items):
        self.by = collections.defaultdict(list)
        for c, s, e, lab in items:
            self.by[c].append((s, e, lab))
        self.starts, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort()
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(e - s for s, e, _ in v)

    @classmethod
    def from_nodes(cls, nodes):
        """ExonIndex(nodes): nodes = [{"idx": label, "chrom": ..., "exons": [(s, e), ...]}]."""
        return cls((n["chrom"], s, e, n["idx"]) for n in nodes for s, e in n["exons"])

    def hits(self, chrom, s, e):
        if chrom not in self.by:
            return set()
        v = self.by[chrom]
        lo = bisect.bisect_left(self.starts[chrom], s - self.maxlen[chrom])
        hi = bisect.bisect_left(self.starts[chrom], e)
        return {i for x0, x1, i in v[lo:hi] if x1 > s and x0 < e}


# ------------------------------------------------------------------------------------------------ annotation
def load_genes():
    """light/work/refseq/genes.tsv -> (genes by gene_id, by_coord: (chrom, start1, end) -> gene rows)."""
    genes = {r["gene_id"]: r for r in tsv(f"{LIGHT}/work/refseq/genes.tsv")}
    by_coord = collections.defaultdict(list)
    for g in genes.values():
        by_coord[(g["chrom"], int(g["start0"]) + 1, int(g["end"]))].append(g)
    return genes, by_coord


def load_exons():
    """light/work/refseq/exons.tsv -> gene_id -> [(start0, end), ...]."""
    return {r["gene_id"]: parse_blocks(r["exons"]) for r in tsv(f"{LIGHT}/work/refseq/exons.tsv")}


def hgnc_tables():
    hg = tsv(HGNC)
    by_id = {r["hgnc_id"]: r for r in hg}
    by_sym = {r["symbol"]: r for r in hg}
    dbx = dict(line.rstrip("\n").split("\t") for line in open(f"{LIGHT}/work/refseq/gene_dbxref.tsv"))
    return by_id, by_sym, dbx


def hgnc_lookup(tables, gene_id, name):
    """HGNC row of a RefSeq gene: the first HGNC dbxref if HGNC knows it, else the symbol. Returns (row or None, how)."""
    by_id, by_sym, dbx = tables
    ids = [x[5:] for x in dbx.get(gene_id, "").split(",") if x.startswith("HGNC:")]
    if ids and ids[0] in by_id:
        return by_id[ids[0]], "dbxref"
    if name in by_sym:
        return by_sym[name], "symbol"
    return None, "none"


def hgnc_all(genes, tables=None):
    """gene_id -> HGNC gene_group_id for every gene with a group (lo_analysis.hgnc_all)."""
    tables = tables or hgnc_tables()
    out = {}
    for g, r in genes.items():
        h, _ = hgnc_lookup(tables, g, r["name"])
        if h and h["gene_group_id"]:
            out[g] = h["gene_group_id"]
    return out


# ---- Soto et al. 2025 gene ids (was bench/layer_order/soto_map.py, body verbatim)
# Map RefSeq CHM13 genes to the Soto et al. 2025 gene-ID convention (CAT CHM13_G* / Liftoff LOFF_G* ids).
# Sources (inventory §1, §7):
#   winloci_data/gencode_chm13/chm13v2.0_gencode.gff3   genome-wide CAT v2.0 annotation (source column 'CAT'); all 2,334
#                                                        Soto-universe ids are present in it (checked, work/cat/)
#   soto_replication/soto_gene_to_families.tsv          gene_id -> family ids (';'-separated) -> ambiguous (yes/no)
#   bench/soto/soto_famCN_S1C.tsv                       Table S1C: Gene Name <-> Gene ID
# Rule: exon unions (RefSeq: exon lines' gene=; CAT: exons of the gene's transcripts) on the same chrom and strand; a
# RefSeq gene maps to the CAT gene with the largest shared exonic bp (ties: larger Jaccard). Reported twice: best CAT gene
# overall (cat_*) and best CAT gene inside Soto's 2,334-gene universe (soto_*), with every Soto-universe gene overlapping
# it. Name agreement is reported, never used to decide. soto_match_quality: strong = shared exonic bp >= 0.5 of BOTH exon
# unions; partial = >= 0.5 of one; weak = neither.
def ov_bp(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            t += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


def parse_blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",")]


def soto_load():
    cat = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{LIGHT}/work/cat/cat_genes_exons.tsv"), delimiter="\t"):
        cat[r["chrom"]].append((int(r["start0"]), int(r["end"]), r["strand"], r["gene_id"], r["name"],
                                parse_blocks(r["exons"])))
    idx = {}
    for c, v in cat.items():
        v.sort()
        idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    fams = {}
    for r in csv.DictReader(open(f"{SOTO}/soto_gene_to_families.tsv"), delimiter="\t"):
        fams[r["gene_id"]] = (r["family_ids_semicolon_sep"], r["ambiguous"])
    names = collections.defaultdict(set)
    name_to_ids = collections.defaultdict(set)
    for r in csv.DictReader(open(S1C), delimiter="\t"):
        names[r["Gene ID"]].add(r["Gene Name"])
        name_to_ids[r["Gene Name"]].add(r["Gene ID"])
    return idx, fams, names, name_to_ids


def soto_map_gene(db, name, chrom, strand, exons):
    idx, fams, names, name_to_ids = db
    v, starts, ml = idx.get(chrom, ([], [], 0))
    s0, e0 = exons[0][0], exons[-1][1]
    lo, hi = bisect.bisect_left(starts, s0 - ml), bisect.bisect_left(starts, e0)
    L = sum(y - x for x, y in exons)
    hits = []
    for a0, a1, st, gid, cname, cex in v[lo:hi]:
        if st != strand or a1 <= s0:
            continue
        o = ov_bp(exons, cex)
        if o > 0:
            CL = sum(y - x for x, y in cex)
            hits.append((o, o / (L + CL - o), gid, cname, CL))
    hits.sort(key=lambda h: (h[0], h[1]), reverse=True)
    out = {k: "" for k in SOTO_FIELDS}
    by_name = sorted(name_to_ids.get(name, set()))
    out["soto_ids_by_name"] = ";".join(by_name)
    out["soto_families_by_name"] = ";".join(sorted({f for i in by_name for f in fams.get(i, ("", ""))[0].split(";") if f}))
    if hits:
        o, jac, gid, cname, CL = hits[0]
        out.update(cat_gene_id=gid, cat_name=cname, cat_in_soto="yes" if gid in fams else "no", cat_ov_bp=o,
                   cat_jaccard=f"{jac:.3f}")
    soto_hits = [h for h in hits if h[2] in fams]
    if soto_hits:
        o, jac, gid, cname, CL = soto_hits[0]
        fam, amb = fams[gid]
        out.update(soto_gene_id=gid, soto_name=";".join(sorted(names.get(gid, set()))) or cname, soto_ov_bp=o,
                   soto_frac_of_refseq=f"{o / L:.3f}", soto_frac_of_soto=f"{o / CL:.3f}",
                   name_match="yes" if (name in names.get(gid, set()) or name == cname) else "no",
                   soto_families=fam, soto_ambiguous=amb,
                   soto_match_quality=("strong" if min(o / L, o / CL) >= 0.5 else "partial" if max(o / L, o / CL) >= 0.5
                                       else "weak"),
                   all_soto_overlapping=";".join(f"{h[2]}({h[3]}:{h[0]}bp:{fams[h[2]][0]})" for h in soto_hits))
    return out


SOTO_FIELDS = ["cat_gene_id", "cat_name", "cat_in_soto", "cat_ov_bp", "cat_jaccard",
               "soto_gene_id", "soto_name", "soto_ov_bp", "soto_frac_of_refseq", "soto_frac_of_soto", "soto_match_quality",
               "name_match", "soto_families", "soto_ambiguous", "all_soto_overlapping", "soto_ids_by_name",
               "soto_families_by_name"]


def soto_flag(m):
    """The U table's Soto flag of one soto_map_gene row: ok iff matched, match quality not weak, not ambiguous."""
    return ("unmatched" if not m["soto_gene_id"] else "weak_match" if m["soto_match_quality"] == "weak"
            else "ambiguous_multi_family" if m["soto_ambiguous"] == "yes" else "ok")


def soto_ok(db, exons, genes, g):
    """Soto families of RefSeq gene g under the U table's flag rule (soto_flag == 'ok'). Returns the family string or
    None. (lo_analysis.soto_ok)"""
    r = genes[g]
    m = soto_map_gene(db, r["name"], r["chrom"], r["strand"], exons.get(g) or [(int(r["start0"]), int(r["end"]))])
    return m["soto_families"] if soto_flag(m) == "ok" and m["soto_families"] else None


# ------------------------------------------------------------------------------------------------ E1 catalogs
def catalog_keys(by_coord, kind, path):
    """record key -> RefSeq gene rows. Fix vs layer_dna.py/truths_universe.py: duplicate-coordinate node rows accumulate
    (setdefault/extend) instead of the last row overwriting the others."""
    key2genes = {}
    if kind == "nodes":
        names = {r["idx"]: r["name"] for r in tsv(path + ".names.tsv")}
        for r in tsv(path):
            k = (r["chrom"], int(r["start"]) + 1, int(r["end"]))
            lst = key2genes.setdefault(k, [])
            for g in by_coord.get(k, []):
                if g["name"] == names[r["idx"]] and g not in lst:
                    lst.append(g)
    else:
        for line in open(path):
            k = key_of(line.strip())
            key2genes[k] = list(by_coord.get(k, []))
    return key2genes


def membership(key2genes, prefix):
    rep = {key_of(r["annotation"]): key_of(r["representative"]) for r in tsv(prefix + ".loci.tsv")}
    cl = {(r["chrom"], int(r["start"]), int(r["end"])): r["cluster_id"] for r in tsv(prefix + ".clusters.tsv")}
    out = {}
    for k, gs in key2genes.items():
        r = rep.get(k, k)
        for g in gs:
            out[g["gene_id"]] = (cl.get(r, ""), f"{k[0]}:{k[1]}-{k[2]}", f"{r[0]}:{r[1]}-{r[2]}" if r != k else "")
    return out


class CatalogContext:
    """What the head of lattice_edges.py computed (lattice_check_c2.py exec'd that head to get key_blocks)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def catalog_context(say):
    """RefSeq genes, U, and the two E1 catalogs' record keys with their exon blocks: key2genes (key -> gene ids),
    key_blocks (key -> merged exon union: c16_19_20 from nodes.tsv, c15_17_22 from light/work/refseq/exons.tsv; a key
    without exons is one span block), key_cat, key_exlen, gene_key. Logs the two lines lattice_edges.py logged."""
    genes, by_coord = load_genes()
    exons_tsv = load_exons()
    U = {r["gene_id"]: r for r in tsv(f"{LIGHT}/universe.corrected.tsv")}
    say(f"[genes] RefSeq gene/pseudogene records {len(genes)}; with exon rows {len(exons_tsv)}; U {len(U)}")
    key2genes, key_blocks, key_cat = {}, {}, {}
    for cat, c in CATALOGS.items():
        k2g = catalog_keys(by_coord, c["kind"], c["path"])
        node_ex = {}
        if c["kind"] == "nodes":
            for r in tsv(c["path"]):
                k = (r["chrom"], int(r["start"]) + 1, int(r["end"]))
                bl = [tuple(map(int, x.split("-"))) for x in r["exons"].split(",")] if r["exons"] else []
                node_ex.setdefault(k, []).extend(bl)
        for k, gs in k2g.items():
            key2genes[k] = [g["gene_id"] for g in gs]
            key_cat[k] = cat
            if c["kind"] == "nodes":
                bl = node_ex.get(k, [])
            else:
                bl = [b for g in gs for b in exons_tsv.get(g["gene_id"], [])]
            if not bl:
                bl = [(k[1] - 1, k[2])]
            key_blocks[k] = merge(bl)
    key_exlen = {k: sum(e - s for s, e in b) for k, b in key_blocks.items()}
    gene_key = {}
    for k, gs in key2genes.items():
        for g in gs:
            gene_key[g] = k
    say(f"[catalogs] keys {len(key2genes)} (c15_17_22 {sum(1 for k in key_cat if key_cat[k] == 'c15_17_22')}, c16_19_20 "
        f"{sum(1 for k in key_cat if key_cat[k] == 'c16_19_20')}); keys with 0 genes {sum(1 for v in key2genes.values() if not v)}; "
        f"keys with >1 gene {sum(1 for v in key2genes.values() if len(v) > 1)}; genes with a key {len(gene_key)}")
    return CatalogContext(genes=genes, by_coord=by_coord, exons_tsv=exons_tsv, U=U, key2genes=key2genes,
                          key_blocks=key_blocks, key_cat=key_cat, key_exlen=key_exlen, gene_key=gene_key)


# ------------------------------------------------------------------------------------------------ DNA: gene-body chains
def chain_groups(recs, Lq):
    """The record groups (chains) of bench/guided_pipeline.gene_body_chains on one (query, target, strand) record list,
    in target order. recs: dicts qs qe ts te nm bl strand."""
    rs = sorted(recs, key=lambda r: (r["ts"], r["te"]))
    out, cur = [], []
    strand = rs[0]["strand"]
    for r in rs + [None]:
        ok = False
        if r is not None and cur:
            gap = r["ts"] - max(x["te"] for x in cur)
            span = r["te"] - min(x["ts"] for x in cur)
            order = r["qs"] >= cur[-1]["qs"] if strand == "+" else r["qs"] <= cur[-1]["qs"]
            ok = gap <= Lq and span <= 2 * Lq and order
        if r is None or (cur and not ok):
            out.append(cur)
            cur = []
        if r is not None:
            cur.append(r)
    return out


def chain_stats(cur, Lq, clen):
    """(passes with the finder's denominator, identity, aligned frac (finder), passes with the shorter-body denominator,
    aligned frac (shorter body), target span start, target span end) of one chain."""
    strand = cur[0]["strand"]
    nm, bl = sum(x["nm"] for x in cur), sum(x["bl"] for x in cur)
    qiv = merge([(x["qs"], x["qe"]) for x in cur])
    aligned = sum(e - s for s, e in qiv)
    ts, te = min(x["ts"] for x in cur), max(x["te"] for x in cur)
    if strand == "+":
        xs, xe = ts - qiv[0][0], te + (Lq - qiv[-1][1])
    else:
        xs, xe = ts - (Lq - qiv[-1][1]), te + qiv[0][0]
    xs, xe = max(0, xs), min(clen, xe)
    den = min(Lq, xe - xs)
    lit = min(Lq, clen)
    return (nm / bl >= 0.80 and aligned >= 0.50 * den, nm / bl, aligned / den if den > 0 else float("inf"),
            nm / bl >= 0.80 and aligned >= 0.50 * lit, aligned / lit, ts, te)


# ------------------------------------------------------------------------------------------------ expression (testis BAM)
def gff_exon_index(strip_prefix):
    """RefSeq GFF (chm13v2.0_RefSeq_full.gff.gz) -> (genes: id -> (chrom, start0, end), exons: id -> {(chrom, s, e)},
    IntervalIndex of every exon labelled by gene id). gene/pseudogene ids with any exon line below them (Parent chain up
    to 10 levels); an exon-less gene/pseudogene counts on its gene body. strip_prefix: label 'gene-X' as 'X'
    (lo_expr_recount) instead of 'gene-X' (lattice_expr)."""
    genes, parent, exraw = {}, {}, []
    with gzip.open(GFF, "rt") as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = f[2]
            if t not in ("gene", "pseudogene", "exon") and "Parent=" not in f[8]:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            fid, par = a.get("ID"), a.get("Parent")
            if t in ("gene", "pseudogene"):
                s, e = int(f[3]) - 1, int(f[4])
                if fid in genes:
                    c0, s0, e0 = genes[fid]
                    if c0 == f[0]:
                        genes[fid] = (c0, min(s0, s), max(e0, e))
                else:
                    genes[fid] = (f[0], s, e)
                continue
            if t == "exon" and par:
                exraw.append((par.split(",")[0], f[0], int(f[3]) - 1, int(f[4])))
            if fid and par:
                parent.setdefault(fid, par.split(",")[0])

    def up(p):
        k = 0
        while p not in genes:
            p = parent.get(p)
            k += 1
            if p is None or k > 10:
                return None
        return p

    exons = collections.defaultdict(set)
    for p, c, s, e in exraw:
        g = up(p)
        if g:
            exons[g].add((c, s, e))
    for g, (c, s, e) in genes.items():
        if not exons.get(g):
            exons[g] = {(c, s, e)}
    index = IntervalIndex((c, s, e, (g[5:] if strip_prefix and g.startswith("gene-") else g))
                          for g, lst in exons.items() for c, s, e in lst)
    return genes, exons, index


def write_windows_bed(path, spans):
    """spans: chrom -> [(start0, end)]; writes the merged windows (touching windows merge). Returns (windows, bp)."""
    nwin = nbp = 0
    with open(path, "w") as out:
        for c in sorted(spans):
            iv = sorted(spans[c])
            cs, ce = iv[0]
            for s, e in iv[1:] + [(None, None)]:
                if s is not None and s <= ce:
                    ce = max(ce, e)
                    continue
                out.write(f"{c}\t{cs}\t{ce}\n")
                nwin += 1
                nbp += ce - cs
                if s is not None:
                    cs, ce = s, e
    return nwin, nbp


def count_reads(bed, index, want, ignore=None):
    """The EXPR count rule (heavy/scripts/expr_counts.py): human_testis.t2t.bam records in the BED windows, primary only
    (samtools view -M -F 2308 -L), one count per (qname, flag, chrom, pos); read blocks split at N and D (M/=/X consume
    both); a read counts for G ('any') if >= 1 block overlaps >= 1 bp of an exon of G (strand ignored); 'unique' = the
    read's blocks hit exons of exactly one record genome-wide and it is G; 'unique_mr' (only with `ignore`) = the hits
    minus the ignored records (minus G itself) are empty.
    Returns (n_any, n_unique, n_unique_mr, primary records, records on >= 1 exon)."""
    p = subprocess.Popen(["samtools", "view", "-M", "-F", "2308", "-L", bed, BAM], stdout=subprocess.PIPE, text=True)
    n_any, n_uni, n_uni_mr = collections.defaultdict(int), collections.defaultdict(int), collections.defaultdict(int)
    seen = set()
    n_rec = n_hit = 0
    for ln in p.stdout:
        f = ln.split("\t", 6)
        key = (f[0], f[1], f[2], f[3])
        if key in seen:
            continue
        seen.add(key)
        n_rec += 1
        c, pos = f[2], int(f[3]) - 1
        hit = set()
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "M=X":
                hit |= index.hits(c, pos, pos + n)
                pos += n
            elif op in "DN":
                pos += n
        if not hit:
            continue
        n_hit += 1
        for g in hit & want:
            n_any[g] += 1
            if ignore and not (hit - ignore - {g}):
                n_uni_mr[g] += 1
        if len(hit) == 1:
            (g,) = hit
            if g in want:
                n_uni[g] += 1
    assert p.wait() == 0
    return n_any, n_uni, n_uni_mr, n_rec, n_hit


# ------------------------------------------------------------------------------------------------ clades (clause 5)
def c_tree(U, NAME):
    """§0★★ clause 5 estimator on the §6js reference trees: supported (SH-aLRT > 75, either tree), pairwise-compatible
    split system on the leaves common to both trees; cluster = smaller side of each split. U: the universe rows by
    gene_id (light/universe.corrected.tsv); NAME: gene_id -> name. Returns ({'Ctree_min'|'Ctree_top'|'Ctree_root':
    labels}, cluster rows, literature rows, {family: (leaves, compatible clusters)}). (lo_analysis.c_tree)"""
    rows = tsv(f"{LIGHT}/C.supported_clades.tsv")
    trees = collections.defaultdict(list)
    leaves = collections.defaultdict(set)
    for r in rows:
        s, c = set(r["side"].split(",")), set(r["complement"].split(","))
        trees[(r["family"], r["tree"])].append((s, c, r["sh_alrt"], r["ufboot"]))
        leaves[(r["family"], r["tree"])] |= s | c
    name2gid = {NAME[g]: g for g in U}
    lab_min, lab_top, lab_root, cl_rows, lit_rows, clusters = {}, {}, {}, [], [], {}
    lit = {r["name"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
    for fam in ("NPIP", "TBC1D3"):
        L = leaves[(fam, "exon")] & leaves[(fam, "intron")]
        splits = {}
        for tr in ("exon", "intron"):
            for s, c, sh, bs in trees[(fam, tr)]:
                s2, c2 = frozenset(s & L), frozenset(c & L)
                if min(len(s2), len(c2)) < 2:
                    continue
                small = min((s2, c2), key=lambda x: (len(x), sorted(x)))
                splits.setdefault(small, []).append(f"{tr} {sh}/{bs}")

        def compat(a, b):  # unrooted splits a|L-a and b|L-b: compatible iff one of the four intersections is empty
            A2, B2 = L - a, L - b
            return not (a & b) or not (a & B2) or not (A2 & b) or not (A2 & B2)
        K = [s for s in splits if all(compat(s, t) for t in splits if t != s)]
        conflicts = {s: sorted(",".join(sorted(t)) for t in splits if t != s and not compat(s, t)) for s in splits}
        minimal = [s for s in K if not any(t < s for t in K)]
        maximal = [s for s in K if not any(s < t for t in K)]
        clusters[fam] = (L, K)
        for s in minimal:
            for n in s:
                lab_min[name2gid[n]] = f"Ctree_min|{fam}|{'+'.join(sorted(s))}"
        for s in maximal:
            for n in s:
                lab_top[name2gid[n]] = f"Ctree_top|{fam}|{'+'.join(sorted(s))}"
        # rooted variant: root at the most balanced compatible split; its two sides are the top clusters
        root = max(K, key=lambda s: (min(len(s), len(L) - len(s)), sorted(s)))
        for side in (root, L - root):
            for n in side:
                lab_root[name2gid[n]] = f"Ctree_root|{fam}|{'+'.join(sorted(side))}"
        for n in L:
            lab_min.setdefault(name2gid[n], f"Ctree|singleton:{name2gid[n]}")
            lab_top.setdefault(name2gid[n], f"Ctree|singleton:{name2gid[n]}")
        for s in sorted(splits, key=lambda x: (len(x), sorted(x))):
            cl_rows.append({"family": fam, "cluster_smaller_side": ",".join(sorted(s)), "size": len(s),
                            "support": "; ".join(splits[s]), "compatible_with_all": "yes" if s in K else "no",
                            "minimal": "yes" if s in minimal else "no", "maximal": "yes" if s in maximal else "no",
                            "conflicts_with": " | ".join(conflicts[s])})
        for lvl in ("level1", "level2", "npipb_named_subfamily"):
            grp = collections.defaultdict(set)
            for n in L:
                if lit.get(n) and lit[n][lvl] and lit[n][lvl] != "no":
                    grp["named NPIPB {B3,B4,B5,B11,B12,B13}" if lvl == "npipb_named_subfamily" else lit[n][lvl]].add(n)
            for gname, G in sorted(grp.items()):
                if len(G) < 2:
                    continue
                Gf, comp = frozenset(G), frozenset(L - G)
                key = Gf if Gf in splits else comp if comp in splits else None
                lit_rows.append({"family": fam, "level": lvl, "literature_group": gname, "genes_on_common_leaves": len(G),
                                 "supported_split_any_tree": "yes" if key else "no",
                                 "support": "; ".join(splits.get(Gf, []) + splits.get(comp, [])),
                                 "in_compatible_system": "yes" if (Gf in K or comp in K) else "no",
                                 "conflicts_with": " | ".join(conflicts[key]) if key else ""})
    return {"Ctree_min": lab_min, "Ctree_top": lab_top, "Ctree_root": lab_root}, cl_rows, lit_rows, clusters


# ------------------------------------------------------------------------------------------------ truth scores
def c2(n):
    return n * (n - 1) // 2


def score_counts(pred, truth, genes):
    """lattice-truth scorer (lattice_truth.score): pairwise P/R from group counts; bipartite F with ONE-TO-ONE Jaccard
    matching (Hungarian on the Jaccard matrix; R = matched genes / genes, P = matched genes / genes of matched predicted
    groups); F is NA when either side has 0 same-group pairs. Values are raw floats (the writer formats them)."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    gs = [g for g in genes if g in pred and g in truth]
    out = {"n_genes": len(gs)}
    if not gs:
        return out
    cell = collections.Counter((pred[g], truth[g]) for g in gs)
    pc = collections.Counter(pred[g] for g in gs)
    tc = collections.Counter(truth[g] for g in gs)
    tp = sum(c2(n) for n in cell.values())
    npp = sum(c2(n) for n in pc.values())
    ntp = sum(c2(n) for n in tc.values())
    out.update({"truth_pairs": ntp, "pred_pairs": npp, "tp_pairs": tp,
                "pair_precision": tp / npp if npp else "NA", "pair_recall": tp / ntp if ntp else "NA"})
    if not npp or not ntp:
        out.update({"bip_R_jaccard": "NA", "bip_P_jaccard": "NA", "bip_F_jaccard": "NA (a side has 0 pairs)"})
        return out
    T, P = sorted(tc), sorted(pc)
    ti, pi = {t: i for i, t in enumerate(T)}, {p: j for j, p in enumerate(P)}
    J = np.zeros((len(T), len(P)))
    M = {}
    for (p, t), n in cell.items():
        J[ti[t], pi[p]] = n / (tc[t] + pc[p] - n)
        M[(ti[t], pi[p])] = n
    r, c = linear_sum_assignment(-J)
    matched = sum(M.get((i, j), 0) for i, j in zip(r, c))
    msize = sum(pc[P[j]] for i, j in zip(r, c) if M.get((i, j), 0) > 0)
    R = matched / len(gs)
    Pp = matched / msize if msize else float("nan")
    F = 2 * R * Pp / (R + Pp) if R + Pp and Pp == Pp else float("nan")
    out.update({"bip_R_jaccard": R, "bip_P_jaccard": Pp, "bip_F_jaccard": F})
    return out


def fmt(x):
    """layer-order number format (lo_analysis.fmt): 3 dp for floats and Fractions, NA for None/NaN."""
    from fractions import Fraction
    if x is None:
        return "NA"
    if isinstance(x, Fraction):
        return f"{float(x):.3f}"
    if isinstance(x, float):
        return "NA" if x != x else f"{x:.3f}"
    return str(x)


def bip_jaccard(pred, true):
    """One-to-one Jaccard bipartite matching (label vectors): returns (matched / all, matched / genes of matched
    predicted groups). (lo_analysis.bip_jaccard)"""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    P, T = sorted(set(pred), key=str), sorted(set(true), key=str)
    M = np.zeros((len(T), len(P)), dtype=int)
    for p, t in zip(pred, true):
        M[T.index(t), P.index(p)] += 1
    ts, ps = M.sum(axis=1), M.sum(axis=0)
    J = np.zeros(M.shape)
    for i in range(len(T)):
        for j in range(len(P)):
            if M[i, j]:
                J[i, j] = M[i, j] / (ts[i] + ps[j] - M[i, j])
    r, c = linear_sum_assignment(-J)
    matched = sum(M[i, j] for i, j in zip(r, c) if M[i, j] > 0)
    msize = sum(ps[j] for i, j in zip(r, c) if M[i, j] > 0)
    return matched / len(pred), (matched / msize if msize else float("nan"))


def f1(r, p):
    return 2 * r * p / (r + p) if r + p and r == r and p == p else float("nan")


def score_lo(lab, truth, genes):
    """layer-order scorer (lo_analysis.score): pair P/R, the §6ks count bipartite (lib.bipartite_items, which was
    guided_pipeline.bipartite) and the Jaccard bipartite, as 3-dp strings; no bipartite F when a side has 0 pairs."""
    from lib import bipartite_items
    gs = sorted(g for g in genes if g in lab and g in truth)
    if not gs:
        return {"n_genes": 0}
    pred = [lab[g] for g in gs]
    true = [truth[g] for g in gs]
    tp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if pred[i] == pred[j] and true[i] == true[j])
    npp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if pred[i] == pred[j])
    ntp = sum(1 for i, j in itertools.combinations(range(len(gs)), 2) if true[i] == true[j])
    out = {"n_genes": len(gs), "truth_pairs": ntp, "pred_pairs": npp, "tp_pairs": tp,
           "pair_precision": fmt(tp / npp) if npp else "NA", "pair_recall": fmt(tp / ntp) if ntp else "NA"}
    if ntp == 0 or npp == 0:  # audit: no bipartite F when either side has no pairs (singleton-to-singleton matches)
        out.update({"bip_F_count(§6ks)": "NA (a side has 0 pairs)", "bip_R": "NA", "bip_P": "NA", "bip_F_jaccard": "NA"})
        return out
    br, bp = bipartite_items(pred, true)
    jr, jp = bip_jaccard(pred, true)
    out.update({"bip_F_count(§6ks)": fmt(f1(br, bp)), "bip_R": fmt(br), "bip_P": fmt(bp), "bip_F_jaccard": fmt(f1(jr, jp))})
    return out


def score_member_anchored(lab, truth_parts, genes, members):
    """pairs (a, b), a != b, both in genes, >= 1 in members. pred: same lab (a gene absent from lab is its own group);
    truth: share >= 1 truth part (truth_parts: gene -> frozenset of group ids)."""
    gs = sorted(g for g in genes if truth_parts.get(g))
    tp = npp = ntp = 0
    for a, b in itertools.combinations(gs, 2):
        if a not in members and b not in members:
            continue
        pr = a in lab and b in lab and lab[a] == lab[b]
        tr = bool(truth_parts[a] & truth_parts[b])
        tp += pr and tr
        npp += pr
        ntp += tr
    return {"n_genes": len(gs), "n_members": len(set(gs) & members), "truth_pairs": ntp, "pred_pairs": npp,
            "tp_pairs": tp, "pair_precision": fmt(tp / npp) if npp else "NA", "pair_recall": fmt(tp / ntp) if ntp else "NA"}
