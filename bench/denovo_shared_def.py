#!/usr/bin/env python3
"""Prereg Addendum AB: the shared family definition (`docs/seeded_family_definition.md` §0★★) instantiated de novo, and
the (a) missing-node / (b) missing-edge / (c) unexpressed-bridge decomposition of the de novo <-> guided gap.

usage:
  denovo_shared_def.py nodes     --dump DIR/ggo.nodes.tsv --genome GGO.fasta --outdir OUT
  denovo_shared_def.py queries   --outdir OUT --genome GGO.fasta
  denovo_shared_def.py align     --outdir OUT --kind tx|body --batch K [--threads 4]
  denovo_shared_def.py families  --outdir OUT
  denovo_shared_def.py decompose --outdir OUT --clusters gw_units_v3.clusters.tsv --expr expr_c3.tsv
                                 --paf-subset guided3.paf name=copies.tsv:nodes.tsv ...
Arms: AB1 = dump nodes unchanged; AB2 = nodes cut at introns > 271,359 bp, pieces < 100 bp dropped, same-strand pieces
with exon overlap merged into gene-level loci. Edges (both arms) = the guided finders: exon edge from the representative
transcript's spliced hit (identity >= 0.80, query coverage >= 0.50) overlapping v's exons; gene-body edge from the body's
asm20 chain (`guided_pipeline.gene_body_chains`) overlapping v's exons. Families = connected components (>= 2 loci).
"""
import argparse
import bisect
import collections
import csv
import glob
import hashlib
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402

CONTIGS = ["NC_073241.2", "NC_073242.2", "NC_073244.2"]
MAX_INTRON, MIN_PIECE = 271_359, 100
BATCH_BP = 4_000_000
COLS = ["idx", "chrom", "start", "end", "strand", "n_exon", "n_reads", "exons"]


def blocks_of(s):
    return [tuple(map(int, b.split("-"))) for b in s.split(",") if b]


def fmt(blocks):
    return ",".join(f"{a}-{b}" for a, b in blocks)


def write_nodes(path, nodes):
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for i, n in enumerate(nodes):
            fh.write("\t".join(str(x) for x in (i, n["chrom"], n["exons"][0][0], n["exons"][-1][1], n["strand"],
                                               len(n["exons"]), n["n_reads"], fmt(n["exons"]))) + "\n")


def read_nodes(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        out.append({"idx": int(r["idx"]), "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                    "strand": r["strand"], "n_reads": int(r["n_reads"]), "exons": blocks_of(r["exons"])})
    return out


def key_of(n):
    return hashlib.md5(f"{n['chrom']}|{n['strand']}|{fmt(n['exons'])}".encode()).hexdigest()[:16]


def cmd_nodes(a):
    os.makedirs(a.outdir, exist_ok=True)
    dump = [r for r in csv.DictReader(open(a.dump), delimiter="\t") if r["chrom"] in CONTIGS]
    ab1 = [{"chrom": r["chrom"], "strand": r["strand"], "n_reads": int(r["n_reads"]), "exons": sorted(blocks_of(r["exons"]))}
           for r in dump]
    write_nodes(f"{a.outdir}/ab1.nodes.tsv", ab1)
    pieces = []
    for n in ab1:
        cur = [n["exons"][0]]
        for b in n["exons"][1:]:
            if b[0] - cur[-1][1] > MAX_INTRON:
                pieces.append(dict(n, exons=cur))
                cur = [b]
            else:
                cur.append(b)
        pieces.append(dict(n, exons=cur))
    n_cut = len(pieces) - len(ab1)
    pieces = [p for p in pieces if sum(e - s for s, e in p["exons"]) >= MIN_PIECE]
    parent = list(range(len(pieces)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    by = collections.defaultdict(list)
    for i, p in enumerate(pieces):
        for s, e in p["exons"]:
            by[(p["chrom"], p["strand"])].append((s, e, i))
    for k, iv in by.items():
        iv.sort()
        end, owner = -1, None
        for s, e, i in iv:
            if owner is not None and s < end:
                parent[find(i)] = find(owner)
            if e > end:
                end, owner = e, i
    groups = collections.defaultdict(list)
    for i in range(len(pieces)):
        groups[find(i)].append(pieces[i])
    ab2 = []
    for g in groups.values():
        rep = max(g, key=lambda p: (p["n_reads"], sum(e - s for s, e in p["exons"])))
        ab2.append({"chrom": rep["chrom"], "strand": rep["strand"], "n_reads": sum(p["n_reads"] for p in g),
                    "exons": gp.merge([b for p in g for b in p["exons"]]), "rep_exons": rep["exons"]})
    ab2.sort(key=lambda n: (n["chrom"], n["exons"][0][0], n["strand"]))
    write_nodes(f"{a.outdir}/ab2.nodes.tsv", ab2)
    with open(f"{a.outdir}/ab2.rep_exons.tsv", "w") as fh:
        fh.write("idx\trep_exons\n")
        for i, n in enumerate(ab2):
            fh.write(f"{i}\t{fmt(n['rep_exons'])}\n")
    print(f"AB1 nodes {len(ab1)}; AB2: chain cuts {n_cut}, pieces kept {len(pieces)}, gene-level loci {len(ab2)}")


def arm_queries(outdir, arm):
    nodes = read_nodes(f"{outdir}/{arm}.nodes.tsv")
    reps = {}
    if arm == "ab2":
        reps = {int(r["idx"]): blocks_of(r["rep_exons"]) for r in csv.DictReader(open(f"{outdir}/ab2.rep_exons.tsv"), delimiter="\t")}
    for n in nodes:
        n["rep"] = reps.get(n["idx"], n["exons"])
        n["tx_key"] = key_of({"chrom": n["chrom"], "strand": n["strand"], "exons": n["rep"]})
        n["body"] = (n["exons"][0][0], n["exons"][-1][1])
        n["body_key"] = hashlib.md5(f"{n['chrom']}|{n['body'][0]}|{n['body'][1]}".encode()).hexdigest()[:16]
    return nodes


def cmd_queries(a):
    genome = pysam.FastaFile(a.genome)
    tx, body = {}, {}
    for arm in ("ab1", "ab2"):
        for n in arm_queries(a.outdir, arm):
            if n["tx_key"] not in tx:
                s = "".join(genome.fetch(n["chrom"], x, y) for x, y in n["rep"]).upper()
                tx[n["tx_key"]] = s.translate(gp.COMP)[::-1] if n["strand"] == "-" else s
            if n["body_key"] not in body:
                body[n["body_key"]] = (n["chrom"], *n["body"])
    for kind, items in (("tx", tx), ("body", body)):
        batches, cur, bp = [], [], 0
        for k in sorted(items, key=lambda k: (len(items[k]) if kind == "tx" else items[k][2] - items[k][1], k)):
            L = len(items[k]) if kind == "tx" else items[k][2] - items[k][1]
            if cur and bp + L > BATCH_BP:
                batches.append(cur)
                cur, bp = [], 0
            cur.append(k)
            bp += L
        if cur:
            batches.append(cur)
        for i, b in enumerate(batches):
            with open(f"{a.outdir}/{kind}.{i:03d}.fa", "w") as fh:
                for k in b:
                    seq = items[k] if kind == "tx" else genome.fetch(*items[k]).upper()
                    fh.write(f">{k}\n{seq}\n")
        print(f"{kind}: {len(items)} unique queries in {len(batches)} batches")
    if not os.path.exists(f"{a.outdir}/ggo3.fa"):
        with open(f"{a.outdir}/ggo3.fa", "w") as fh:
            for c in CONTIGS:
                fh.write(f">{c}\n{genome.fetch(c)}\n")


def cmd_align(a):
    fa = f"{a.outdir}/{a.kind}.{a.batch:03d}.fa"
    out = f"{a.outdir}/{a.kind}.{a.batch:03d}.paf"
    if os.path.exists(out):
        print(f"{out} exists")
        return
    preset = ["-x", "splice", "-uf"] if a.kind == "tx" else ["-x", "asm20"]
    mmi = f"{a.outdir}/ggo3.{'splice' if a.kind == 'tx' else 'asm20'}.mmi"
    if not os.path.exists(mmi):
        subprocess.run(["minimap2", *preset, "-t", str(a.threads), "-d", mmi, f"{a.outdir}/ggo3.fa"], check=True,
                       stderr=subprocess.DEVNULL)
    with open(out + ".tmp", "w") as fh:
        subprocess.run(["minimap2", "-c", *preset, "-N", "50", "-p", "0.1", "-t", str(a.threads), mmi, fa], stdout=fh,
                       stderr=subprocess.DEVNULL, check=True)
    os.replace(out + ".tmp", out)
    print(f"wrote {out}")


class ExonIndex:
    def __init__(self, nodes):
        self.by = collections.defaultdict(list)
        for n in nodes:
            for s, e in n["exons"]:
                self.by[n["chrom"]].append((s, e, n["idx"]))
        self.starts, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort()
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(e - s for s, e, _ in v)

    def hits(self, chrom, s, e):
        if chrom not in self.by:
            return set()
        v = self.by[chrom]
        lo = bisect.bisect_left(self.starts[chrom], s - self.maxlen[chrom])
        hi = bisect.bisect_left(self.starts[chrom], e)
        return {i for x0, x1, i in v[lo:hi] if x1 > s and x0 < e}


def cmd_families(a):
    tx_hits = [h for f in sorted(glob.glob(f"{a.outdir}/tx.*.paf")) for h in gp.transcript_hits(f)]
    body_chains = [c for f in sorted(glob.glob(f"{a.outdir}/body.*.paf")) for c in gp.gene_body_chains(f)]
    by_tx, by_body = collections.defaultdict(list), collections.defaultdict(list)
    for h in tx_hits:
        by_tx[h["q"]].append(h)
    for c in body_chains:
        by_body[c["q"]].append(c)
    for arm in ("ab1", "ab2"):
        nodes = arm_queries(a.outdir, arm)
        idx = ExonIndex(nodes)
        spliced = {n["idx"]: len(n["exons"]) >= 2 for n in nodes}
        edges = collections.Counter()
        pairs = set()
        for u in nodes:
            us, ue = u["exons"][0][0], u["exons"][-1][1]
            cand = []
            for h in by_tx.get(u["tx_key"], []):
                orient = h["strand"]  # transcript-oriented query (-uf): '+' = transcript forward on the genome
                for s, e in gp.tx_exon_blocks(h):
                    cand.append(("exon", h["chrom"], s, e, h["s"], h["e"], orient))
            for c in by_body.get(u["body_key"], []):
                orient = u["strand"] if c["strand"] == "+" else {"+": "-", "-": "+"}.get(u["strand"], u["strand"])
                cand.append(("body", c["chrom"], c["s"], c["e"], c["s"], c["e"], orient))
            for kind, chrom, s, e, hs, he, orient in cand:
                if chrom == u["chrom"] and hs < ue and us < he:
                    continue  # u's own locus
                for v in idx.hits(chrom, s, e):
                    if v == u["idx"]:
                        continue
                    if spliced[u["idx"]] and spliced[v] and nodes[v]["strand"] != orient:
                        continue
                    p = (min(u["idx"], v), max(u["idx"], v))
                    if (p, kind) not in pairs:
                        pairs.add((p, kind))
                        edges[kind] += 1
        adj = collections.defaultdict(set)
        for (p, _k) in pairs:
            adj[p[0]].add(p[1])
            adj[p[1]].add(p[0])
        seen, fams = set(), []
        for n in nodes:
            if n["idx"] in seen or n["idx"] not in adj:
                continue
            comp, stack = [], [n["idx"]]
            seen.add(n["idx"])
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            fams.append(sorted(comp))
        fams = [f for f in fams if len(f) >= 2]
        with open(f"{a.outdir}/{arm}.copies.tsv", "w") as fh:
            fh.write("family_id\tchrom\tstart\tend\tstrand\tn_reads\tnode_idx\n")
            for k, f in enumerate(sorted(fams, key=lambda f: (-len(f), f[0]))):
                for i in f:
                    n = nodes[i]
                    fh.write(f"SDF{k}\t{n['chrom']}\t{n['exons'][0][0]}\t{n['exons'][-1][1]}\t{n['strand']}\t{n['n_reads']}\t{i}\n")
        both = len({p for p, _ in pairs})
        print(f"{arm}: nodes {len(nodes)}; exon edges {edges['exon']}, gene-body edges {edges['body']}, distinct pairs {both}; "
              f"families (>=2 loci) {len(fams)} holding {sum(len(f) for f in fams)} loci")


def cmd_decompose(a):
    expr = {(r["chrom"], int(r["start"]), int(r["end"])): int(r["u"]) for r in csv.DictReader(open(a.expr), delimiter="\t")}
    g = [(r["cluster_id"], r["chrom"], int(r["start"]) - 1, int(r["end"]))
         for r in csv.DictReader(open(a.clusters), delimiter="\t") if r["chrom"] in CONTIGS]
    allg = g
    ex = [x for x in g if expr.get((x[1], x[2], x[3]), 0) >= 3]
    cnt = collections.Counter(x[0] for x in ex)
    loci = [x for x in ex if cnt[x[0]] >= 2]
    by_cluster = collections.defaultdict(list)
    for li, x in enumerate(loci):
        by_cluster[x[0]].append(li)
    universe = [(p, q) for v in by_cluster.values() for i, p in enumerate(v) for q in v[i + 1:]]

    # guided edges among ALL guided loci on the contigs (re-derived from the all-genes PAF subset)
    L = collections.defaultdict(list)
    for gi, x in enumerate(allg):
        L[x[1]].append((x[2], x[3], gi))
    for c in L:
        L[c].sort()

    def locus_of(name):
        c, rng = name.rsplit(":", 1)
        s, e = rng.split("-")
        s, e = int(s) - 1, int(e)
        best = None
        for x0, x1, gi in L.get(c, []):
            if x0 >= e:
                break
            o = min(e, x1) - max(s, x0)
            if o > 0 and (best is None or o > best[0]):
                best = (o, gi)
        return None if best is None else best[1]
    gedge = set()
    memo = {}
    for line in open(a.paf_subset):
        f = line.split("\t")
        if f[0] not in memo:
            memo[f[0]] = locus_of(f[0])
        if f[5] not in memo:
            memo[f[5]] = locus_of(f[5])
        u, v = memo[f[0]], memo[f[5]]
        if u is None or v is None or u == v:
            continue
        nm, bl = int(f[9]), int(f[10])
        longer = max(allg[u][3] - allg[u][2], allg[v][3] - allg[v][2])
        aligned = max(int(f[3]) - int(f[2]), int(f[8]) - int(f[7]))
        if bl and nm / bl >= 0.70 and aligned >= 300 and aligned >= 0.30 * longer:
            gedge.add((min(u, v), max(u, v)))
    key_to_g = {(x[1], x[2], x[3]): gi for gi, x in enumerate(allg)}
    lg = [key_to_g[(x[1], x[2], x[3])] for x in loci]

    def comp_label(members):
        gi = {lg[m] for m in members}
        parent = {x: x for x in gi}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for u, v in gedge:
            if u in gi and v in gi:
                parent[find(u)] = find(v)
        return {m: find(lg[m]) for m in members}
    comp = {}
    for v in by_cluster.values():
        comp.update(comp_label(v))
    n_c = sum(1 for p, q in universe if comp[p] != comp[q])
    print(f"universe: {len(universe)} expressed co-clustered guided pairs; guided edges re-derived {len(gedge)}; "
          f"pairs whose loci connect only through unexpressed loci (class c, catalog-independent): {n_c}")

    def overlaps_any(rows, x):
        return any(r[0] == x[1] and r[1] < x[3] and x[2] < r[2] for r in rows)
    for spec in a.catalogs:
        name, rest = spec.split("=", 1)
        cpath, npath = rest.split(":", 1)
        copies = [r for r in csv.DictReader(open(cpath), delimiter="\t") if r["chrom"] in CONTIGS]
        fams = collections.defaultdict(set)
        cop = collections.defaultdict(list)
        for r in copies:
            cop[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
        for li, x in enumerate(loci):
            for s, e, f in cop.get(x[1], []):
                if s < x[3] and x[2] < e:
                    fams[li].add(f)
        nodes = collections.defaultdict(list)
        for r in csv.DictReader(open(npath), delimiter="\t"):
            nodes[r["chrom"]].append((r["chrom"], int(r["start"]), int(r["end"])))
        has_node = [overlaps_any(nodes.get(x[1], []), x) for x in loci]
        cls = collections.Counter()
        for p, q in universe:
            missing = not (has_node[p] and has_node[q])
            if fams[p] & fams[q]:
                cls["recovered"] += 1
            elif a.order == "cab" and comp[p] != comp[q]:
                cls["c_unexpressed_bridge"] += 1
            elif missing:
                cls["a_missing_node"] += 1
            elif comp[p] != comp[q]:
                cls["c_unexpressed_bridge"] += 1
            else:
                cls["b_missing_edge"] += 1
        rec_c = sum(1 for p, q in universe if comp[p] != comp[q] and fams[p] & fams[q])
        print(f"{name:10s} " + "  ".join(f"{k} {cls[k]} ({cls[k] / len(universe):.3f})" for k in
                                          ("recovered", "a_missing_node", "b_missing_edge", "c_unexpressed_bridge"))
              + f"  | recovered pairs that are class-c at guided level: {rec_c}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("nodes")
    p.add_argument("--dump", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True)
    p = sub.add_parser("queries")
    p.add_argument("--outdir", required=True)
    p.add_argument("--genome", required=True)
    p = sub.add_parser("align")
    p.add_argument("--outdir", required=True)
    p.add_argument("--kind", choices=("tx", "body"), required=True)
    p.add_argument("--batch", type=int, required=True)
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("families")
    p.add_argument("--outdir", required=True)
    p = sub.add_parser("decompose")
    p.add_argument("--outdir", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--expr", required=True)
    p.add_argument("--paf-subset", required=True)
    p.add_argument("--order", choices=("cab", "acb"), default="cab",
                   help="cab = prereg AB order; acb = missing node checked first (sensitivity)")
    p.add_argument("catalogs", nargs="+")
    a = ap.parse_args()
    {"nodes": cmd_nodes, "queries": cmd_queries, "align": cmd_align, "families": cmd_families, "decompose": cmd_decompose}[a.cmd](a)


if __name__ == "__main__":
    main()
