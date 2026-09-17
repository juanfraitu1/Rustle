#!/usr/bin/env python3
"""Nested edge-test lattice (NPIP/TBC1D3, human CHM13) — shared paths, loaders and THE LEVEL TESTS.

User-approved design (2026-09-16 16:39): one graph (nodes = RefSeq gene records at annotated gene-body extent; every
candidate link carries all its evidence); level k = connected components of the edges passing t_k, with
t_{k+1} => t_k by construction (each test is the previous test AND one more clause):

  t0 superfamily            = not same_locus AND (protein §6ko edge OR t1)
  t1 family                 = not same_locus AND clause-2 DNA homology (APPROXIMATION from the E1 gene-body PAF; the hit
                              must overlap v's exons and pass the shipped strand check for spliced pairs, see
                              lattice_edges.py)
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
import collections
import csv
import itertools

LIGHT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light"
HEAVY = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/heavy"
INT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/integrate_slim"
OUT = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/lattice"
H = "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge"
PAF = {"c15_17_22": f"{H}/human2/genes.asm20.paf", "c16_19_20": f"{H}/lit/aj_ho/refseq/all.paf"}
DUMP = {"c15_17_22": f"{LIGHT}/work/D/c15_17_22.e1.graph.tsv", "c16_19_20": f"{LIGHT}/work/D/c16_19_20.e1.graph.tsv"}
DUMP_S1 = {"c15_17_22": f"{LIGHT}/work/S1/c15_17_22.e1s.graph.tsv", "c16_19_20": f"{LIGHT}/work/S1/c16_19_20.e1s.graph.tsv"}
CAT_CHROMS = {"c15_17_22": {"chr15", "chr17", "chr22"}, "c16_19_20": {"chr16", "chr19", "chr20"}}

SEF_MIN = 0.30
ID_L3 = 0.98
LEVELS = ("L0", "L1", "L2", "L3")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(path, rows, cols=None):
    cols = cols or (list(rows[0].keys()) if rows else ["empty"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmtv(r.get(c, "")) for c in cols) + "\n")


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


def pairs_in(lab, genes):
    s = set()
    for G in groups({g: lab[g] for g in genes if g in lab}).values():
        for a, b in itertools.combinations(sorted(G), 2):
            s.add((a, b))
    return s
