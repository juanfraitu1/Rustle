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

    famA = build(all_edges)                                          # all biotypes (raw)
    pc = lambda i: genes[i]["biotype"] == "protein_coding"
    famPC = build(all_edges, keep=pc)                                # protein-coding only (refined)
    print(f"    {n_edges} edges -> {len(famA)} families all-biotype ({sum(len(m) for m in famA)} genes); "
          f"{len(famPC)} protein-coding families ({sum(len(m) for m in famPC)} genes)")

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
    def write_catalog(path, fams):
        with open(path, "w") as fh:
            fh.write("family_id\tn_genes\tn_contigs\tcross_chrom\tdom_biotype\tcompara_paralog\tgenes\n")
            for k, m in enumerate(fams):
                nn = {named(i) for i in m} - {None}
                comp_hit = any(tuple(sorted((x, y))) in cpairs for x in nn for y in nn if x < y)
                fh.write(f"GFAM{k}\t{len(m)}\t{len({genes[i]['chrom'] for i in m})}\t"
                         f"{int(len({genes[i]['chrom'] for i in m})>=2)}\t{dom_biotype(m)}\t"
                         f"{int(comp_hit)}\t{','.join(genes[i]['name'] for i in m)}\n")

    outA = os.path.join(bench, "genome_families_annotated.tsv")
    write_catalog(outA, famA)
    outPC = os.path.join(bench, "genome_families_protein_coding.tsv")
    write_catalog(outPC, famPC)
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
    json.dump(stats, open(os.path.join(bench, "genome_family_def.json"), "w"), indent=2)

    def top(fams, n=6):
        return ", ".join(f"{genes[m[0]]['name']}({len(m)},{dom_biotype(m)})" for m in fams[:n])
    print("\n== GENOME-GROUNDED MULTI-COPY FAMILY DEFINITION (read-independent, parallel) ==")
    print(f"  RAW (all biotypes)   {len(famA)} families, {stats['catalogA_member_genes']} genes, "
          f"{xchrom} cross-chrom")
    print(f"      top: {top(famA)}")
    print(f"      -> dominated by tRNA/snRNA/lncRNA repeat arrays (RNA-gene multi-copy, repeat-bridged)")
    print(f"  REFINED protein-coding  {len(famPC)} families, {stats['catalogPC_member_genes']} genes, "
          f"{xchrom_pc} cross-chrom")
    print(f"      top: {top(famPC)}")
    print(f"  CATALOG B (de-novo blocks)  {len(blocks)} duplication blocks (>=2 units & >=2 genes), {uid} units")
    print(f"  Compara cross-check  {a_compara_hit}/{a_named} multi-named families carry a Compara "
          f"paralog relation (sparse Compara coverage)")
    print(f"  A vs B granularity   {a_in_one_block}/{len(famA)} families inside one B-block "
          f"({100*a_in_one_block/max(1,len(famA)):.0f}%)")
    print(f"  wrote {outA}\n        {outPC}\n        {outB}\n        {os.path.join(bench,'genome_family_def.json')}")


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
