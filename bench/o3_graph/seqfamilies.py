#!/usr/bin/env python3
"""Name-free gene families for the integrated O1 -> O2 -> O3 simulation (2026-10-17; user: "we should go by families linked by
sequence identity, also if possible by intron/exon structure"). The name-stem pools (rounds 11-16) used the annotation's gene
names and missed every LOC / clone-named / hyphenated paralog: 1,486 of CAT's 20,067 protein-coding genes.

Genes: every CAT protein-coding gene, whatever its name, span <= families.MAX_SPAN, overlaps resolved as families.py does
(the gene with more transcripts kept); each gene's canonical transcript (families.canonical).
Sequence link (families.py's locus rule, genome-wide): gene A's canonical, spliced-aligned (minimap2 -cx splice:hq) to gene B's
locus (+- 1 kb), covers >= MIN_COVER of A at identity >= families.MIN_IDENT, in either direction.
Structure link: one gene's exon chain is carried by the other's GENOME -- in the spliced alignment of A's canonical to B's locus,
every intron sits within realworld.JUNCTION_TOL of one of A's junctions and none of A's junctions is left without one
(structure_match), in either direction. B's annotation is not used: comparing two annotated models failed identical genes whose
canonicals differ (OR4F21's 2-exon model against OR4F16's 1-exon one; ATAD3A's extra 80-bp exon), a choice of the annotation,
not of the genomes (amended before any draw). A processed retrocopy (one exon against a spliced locus) or a partial copy fails.
Families: connected components of the links that pass both; a component touching a gene of an earlier set is spent.

    seqfamilies.py build OUT.json SPENT.json ...   -> OUT.json (set seq17; every family), OUT.json.links.tsv (every sequence
                                                       link and its structure verdict)
    seqfamilies.py draw OUT_SET.json OUT.json      -> the round's set (draw: seeded order up to a copy budget)
"""
import collections
import csv
import importlib.util
import json
import os
import random
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_COVER, PAD, JUNCTION_TOL = 0.8, 1000, 10          # realworld.JUNCTION_TOL
SEED, BUDGET, MAX_TOTAL = 17, 90, 150                 # the draw, fixed before the build
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _families():
    spec = importlib.util.spec_from_file_location("o3_graph_families", os.path.join(HERE, "families.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def structure_match(junctions, qlen, strand, ts, cigar, tol=JUNCTION_TOL, q0=0):
    """Is gene A's model carried by gene B's genome? junctions: A's exon junctions in its transcript coordinates, qlen its
    spliced length; the alignment of that transcript to B's locus (target start ts, CIGAR [(n, op)], strand). Every intron (N)
    of the alignment sits within tol of one of A's junctions, and as many introns as junctions. q0: the alignment's first query
    base in the alignment's orientation (a - strand alignment runs A backwards)."""
    js = sorted(qlen - j for j in junctions) if strand == "-" else sorted(junctions)
    q, t, got = q0, ts, []
    for n, op in cigar:
        n = int(n)
        if op in "M=X":
            q, t = q + n, t + n
        elif op == "I":
            q += n
        elif op == "D":
            t += n
        elif op == "N":
            got.append(q)
            t += n
    return len(got) == len(js) and all(abs(qj - j) <= tol for qj, j in zip(got, js))


def components(links, spent=frozenset()):
    """Connected components (>= 2 genes) of the links (a, b, passes) that pass, less those touching a spent gene."""
    par = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b, ok in links:
        if ok:
            par[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for x in list(par):
        comp[find(x)].append(x)
    out = [sorted(c) for c in comp.values() if len(c) >= 2 and not spent & set(c)]
    return sorted(out, key=lambda c: c[0])


def draw(fams, seed=SEED, budget=BUDGET, max_total=MAX_TOTAL):
    """The round's families: in a seeded order of the families sorted by name, each taken unless it would take the copy total
    past max_total, until the total reaches budget."""
    order = sorted(fams, key=lambda f: f["family"])
    random.Random(seed).shuffle(order)
    got, total = [], 0
    for f in order:
        if total >= budget:
            break
        if total + len(f["copies"]) <= max_total:
            got.append(f)
            total += len(f["copies"])
    return sorted(got, key=lambda f: f["family"])


def junctions_of(exons):
    """A transcript's exon junctions in its own 5' -> 3' coordinates: cumulative exon lengths (exons in genome order; on the
    minus strand the transcript starts at the last exon)."""
    lens = [e - s for s, e in sorted(exons)]
    return [sum(lens[:i]) for i in range(1, len(lens))]


def spent_genes(paths):
    out = set()
    for p in paths:
        for fams in json.load(open(p)).get("sets", {}).values():
            for f in fams:
                out |= {c["gene_id"] for c in f["copies"]}
    return out


def cmd_build(out_p, spent_paths, threads=4):
    import pysam
    F = _families()
    gff = pysam.TabixFile(f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz")
    fa = pysam.FastaFile(F.GENOME)
    genes = []
    for r in csv.DictReader(open(f"{F.CAT}/chm13v2.0_CAT_Liftoff.genes.tsv"), delimiter="\t"):
        if r["gene_biotype"] == "protein_coding" and int(r["end"]) - int(r["start0"]) <= F.MAX_SPAN:
            genes.append(dict(gene_id=r["gene_id"], gene_name=r["gene_name"], chrom=r["chrom"], start=int(r["start0"]),
                              end=int(r["end"]), strand=r["strand"], n=int(r["n_transcripts"])))
    genes = F.drop_overlaps(genes)
    for g in genes:
        g["tx"] = F.gene_transcripts(gff, g)
    genes = [g for g in genes if g["tx"]]
    by = {g["gene_id"]: g for g in genes}
    canon = {g["gene_id"]: g["tx"][F.canonical(g["tx"])] for g in genes}
    print(f"[seqfamilies] {len(genes)} protein-coding genes", flush=True)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        with open(f"{d}/c.fa", "w") as fh:
            for gid, ex in canon.items():
                g = by[gid]
                fh.write(f">{gid}\n{F.spliced(fa, g['chrom'], g['strand'], ex)}\n")
        # candidates: an all-vs-all of the canonicals at the family identity over half the shorter
        paf = subprocess.run(["minimap2", "-cx", "asm20", "-X", "-t", str(threads), f"{d}/c.fa", f"{d}/c.fa"],
                             capture_output=True, text=True, check=True).stdout
        cand = set()
        for ln in paf.splitlines():
            f = ln.split("\t")
            if f[0] != f[5] and int(f[9]) >= F.MIN_IDENT * int(f[10]) and int(f[10]) >= 0.5 * min(int(f[1]), int(f[6])):
                cand |= {f[0], f[5]}
        print(f"[seqfamilies] {len(cand)} genes with a candidate partner", flush=True)
        # the locus rule on every candidate gene: its canonical onto every candidate locus
        with open(f"{d}/q.fa", "w") as fh:
            for gid in sorted(cand):
                g = by[gid]
                fh.write(f">{gid}\n{F.spliced(fa, g['chrom'], g['strand'], canon[gid])}\n")
        with open(f"{d}/t.fa", "w") as fh:
            for gid in sorted(cand):
                g = by[gid]
                fh.write(f">{gid}\n{fa.fetch(g['chrom'], max(0, g['start'] - PAD), g['end'] + PAD).upper()}\n")
        paf = subprocess.run(["minimap2", "-cx", "splice:hq", "-N", "1000", "-p", "0.01", "--secondary=yes", "-t", str(threads),
                              f"{d}/t.fa", f"{d}/q.fa"], capture_output=True, text=True, check=True).stdout
    blocks, best = collections.defaultdict(list), {}
    for ln in paf.splitlines():
        f = ln.split("\t")
        a, b = f[0], f[5]
        if a == b or int(f[9]) < F.MIN_IDENT * int(f[10]):
            continue
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        blocks[(a, b)].append((qs, qe))
        cov = (qe - qs) / qlen
        if cov > best.get((a, b), (0.0, None))[0]:
            best[(a, b)] = (cov, f)
    cover = {k: F.union_cover(v, int(best[k][1][1])) for k, v in blocks.items()}

    def struct(a, b):
        if (a, b) not in best:
            return False
        f = best[(a, b)][1]
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        q0 = qlen - qe if f[4] == "-" else qs
        # A's junctions in A's transcript coordinates (5' -> 3'): on the minus strand the transcript starts at the last exon
        js = junctions_of(canon[a])
        js = sorted(qlen - j for j in js) if by[a]["strand"] == "-" else js
        return structure_match(js, qlen, f[4], int(f[7]), [(int(n), op) for n, op in _CIGAR.findall(cg)], q0=q0)
    pairs = {tuple(sorted(k)) for k in cover}
    links, rows = [], []
    for a, b in sorted(pairs):
        cab, cba = cover.get((a, b), 0.0), cover.get((b, a), 0.0)
        if max(cab, cba) < MIN_COVER:
            continue
        sab, sba = struct(a, b), struct(b, a)
        links.append((a, b, sab or sba))
        rows.append([a, b, by[a]["gene_name"], by[b]["gene_name"], f"{cab:.3f}", f"{cba:.3f}", int(sab), int(sba), int(sab or sba)])
    spent = spent_genes(spent_paths)
    comps = components(links, spent=spent)
    n_seq = len(components([(a, b, True) for a, b, _ in links]))
    fams, used = [], collections.Counter()
    for comp in comps:
        gs = sorted((by[x] for x in comp), key=lambda g: (g["chrom"], g["start"]))
        base = re.sub(r"[^A-Za-z0-9]", "", gs[0]["gene_name"]) or "G"
        used[base] += 1
        name = f"S17{base}" + (f"x{used[base]}" if used[base] > 1 else "")
        copies = []
        for g in gs:
            main, iso = F.pick_isoforms(g["tx"])
            copies.append(dict(gene_id=g["gene_id"], gene_name=g["gene_name"], chrom=g["chrom"],
                               start=min(s for v in iso.values() for s, _ in v), end=max(e for v in iso.values() for _, e in v),
                               strand=g["strand"], main=main, transcripts={t: [list(x) for x in v] for t, v in iso.items()}))
        fams.append(dict(family=name, candidates=len(comp), component=len(comp), copies=copies))
    with open(f"{out_p}.links.tsv", "w") as fh:
        fh.write("a\tb\ta_name\tb_name\tcover_ab\tcover_ba\tstruct_ab\tstruct_ba\tstruct_link\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    json.dump(dict(genome=F.GENOME, annotation=f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz", rule="seqfamilies",
                   min_ident=F.MIN_IDENT, min_cover=MIN_COVER, junction_tol=JUNCTION_TOL, sets={"seq17": fams}),
              open(out_p, "w"), indent=1)
    unnamed = sum(any(re.match(r"^(LOC|ENSG|MSTRG|LINC|AC[0-9]|AL[0-9]|AP[0-9])", c["gene_name"]) or "-" in c["gene_name"]
                      for c in f["copies"]) for f in fams)
    print(f"[seqfamilies] sequence links {len(links)}, of them also structure {sum(ok for _, _, ok in links)}; "
          f"components by sequence {n_seq}, by sequence and structure {len(components(links))}, of them unspent {len(fams)} "
          f"({sum(len(f['copies']) for f in fams)} copies; {unnamed} with an unnamed copy) -> {out_p}", flush=True)
    return 0


def cmd_draw(out_p, all_p):
    d = json.load(open(all_p))
    got = draw(d["sets"]["seq17"])
    json.dump(dict({k: v for k, v in d.items() if k != "sets"}, draw=dict(seed=SEED, budget=BUDGET, max_total=MAX_TOTAL),
                   sets={"r17": got}), open(out_p, "w"), indent=1)
    print(f"[seqfamilies] drew {len(got)} families, {sum(len(f['copies']) for f in got)} copies: "
          + ", ".join(f"{f['family']} ({len(f['copies'])})" for f in got), flush=True)
    return 0


def main(argv):
    if argv[0] == "build":
        return cmd_build(argv[1], argv[2:])
    if argv[0] == "draw":
        return cmd_draw(argv[1], argv[2])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
