#!/usr/bin/env python3
"""Name-free gene families for the integrated O1 -> O2 -> O3 simulation (2026-10-17; user: "we should go by families linked by
sequence identity, also if possible by intron/exon structure"). The name-stem pools (rounds 11-16) used the annotation's gene
names and missed every LOC / clone-named / hyphenated paralog: 1,486 of CAT's 20,067 protein-coding genes.

Genes: every CAT protein-coding gene, whatever its name, span <= families.MAX_SPAN, overlaps resolved as families.py does
(the gene with more transcripts kept); each gene's canonical transcript (families.canonical).
Sequence link (families.py's locus rule, genome-wide): gene A's canonical, spliced-aligned (minimap2 -cx splice:hq) to gene B's
locus (+- 1 kb), covers >= MIN_COVER of A at identity >= families.MIN_IDENT, in either direction.
Structure link: one gene's intron chain is carried by the other's GENOME -- in the spliced alignment of A's canonical to B's locus,
every intron sits within realworld.JUNCTION_TOL of one of A's junctions and none of A's junctions is left without one
(structure_match, chain "exact"), in either direction. --intron-chain carried (2026-10-18): each of A's junctions has its
intron in B's genome, which may hold more (GOLGA6L; structure_match, chain "carried"). B's annotation is not used: comparing two annotated models failed identical genes whose
canonicals differ (OR4F21's 2-exon model against OR4F16's 1-exon one; ATAD3A's extra 80-bp exon), a choice of the annotation,
not of the genomes (amended before any draw). A processed retrocopy (one exon against a spliced locus) or a partial copy fails.
Families: connected components of the links that pass both; a component touching a gene of an earlier set is spent.
--paralogs (round 19; user 2026-10-18, after round 18's OR4F pseudogenes): genes of every other biotype join as copies like any
other (add_paralog_genes: none overlapping a coding gene; a pseudogene before an overlapping lncRNA). A paralog links to a coding
gene when its locus carries the coding gene's canonical by the same rule (sequence and structure, in that direction only:
pair_link); two paralogs do not link on their own. Families are named after their first coding gene (S19...; set seq19).

    seqfamilies.py build OUT.json SPENT.json ... [--paralogs] [--intron-chain exact|carried] [--retain-chains]
                                                   -> OUT.json (set seq17, or seq19: every unspent family), OUT.json.links.tsv
                                                       (every sequence link and its structure verdict)
    seqfamilies.py draw OUT_SET.json OUT.json [--seed S --set NAME]
                                                   -> the round's set (draw: seeded order up to a copy budget; round 17: seed 17,
                                                       set r17; round 18: seed 18, set r18)
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


def structure_match(junctions, qlen, strand, ts, cigar, tol=JUNCTION_TOL, q0=0, chain="exact"):
    """Is gene A's intron chain carried by gene B's genome? junctions: A's exon junctions in its transcript coordinates, qlen its
    spliced length; the alignment of that transcript to B's locus (target start ts, CIGAR [(n, op)], strand). Only junctions
    inside the aligned part of A count (one at the edge of an unaligned end is a terminal exon too short to align: LGALS7's
    6-bp last exon). q0: the alignment's first query base in the alignment's orientation (a - strand alignment runs A
    backwards).
    chain "exact" (rounds 17-19, and O1's genome search): every intron (N) of the alignment sits within tol of one of A's
    junctions, and there are as many introns as junctions -- B's genome holds A's intron chain and no other intron.
    chain "carried" (2026-10-18, the user: "intron chains"): each of A's junctions has an intron of the alignment within tol,
    and B's genome may hold more -- GOLGA6L's copies carry each other's every junction, but some split a long exon with one or
    two short introns of their own, and indels of 18-252 bp there come out of the spliced alignment as introns. A model with
    no junction in the aligned part is carried only by an alignment with no intron, so a retrocopy and its multi-exon parent
    still never carry each other."""
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
    js = [j for j in js if q0 + tol < j < q - tol]          # a junction at an unaligned end: a terminal exon too short to align
    if chain == "carried":
        return all(any(abs(g - j) <= tol for g in got) for j in js) if js else not got
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


def add_paralog_genes(coding, noncoding):
    """--paralogs (round 19, user 2026-10-18: paralogs of any biotype are copies like any other): the non-coding genes a family may
    take. None overlapping a kept coding gene (another annotation of that gene, as OR4F1P inside OR4F16); among overlapping ones
    a pseudogene first (the paralogous copy itself; a lncRNA over it is another transcript of the locus), then the one with more
    transcripts (families.drop_overlaps)."""
    free = [g for g in noncoding
            if not any(k["chrom"] == g["chrom"] and k["start"] < g["end"] and g["start"] < k["end"] for k in coding)]
    keep = []
    for g in sorted(free, key=lambda g: ("pseudogene" not in g["biotype"], -g["n"], g["chrom"], g["start"])):
        if not any(k["chrom"] == g["chrom"] and k["start"] < g["end"] and g["start"] < k["end"] for k in keep):
            keep.append(g)
    return sorted(keep, key=lambda g: (g["chrom"], g["start"]))


def pair_link(a, b, cover, struct, coding):
    """A gene pair's link: None (no sequence link), else whether it also carries structure. cover / struct: {(query gene, locus
    gene): value}. Two coding genes, as before: either direction (cover >= MIN_COVER; one gene's exon chain carried by the
    other's genome). A coding gene and a paralog (--paralogs): only the paralog's locus carrying the coding gene's canonical
    counts -- a paralog's own model inside a coding locus (a small RNA in an intron) is no copy of it. Two paralogs: no link of
    their own; they join a family through its coding genes."""
    ca, cb = a in coding, b in coding
    if ca and cb:
        if max(cover.get((a, b), 0.0), cover.get((b, a), 0.0)) < MIN_COVER:
            return None
        return bool(struct.get((a, b)) or struct.get((b, a)))
    if not (ca or cb):
        return None
    g, p = (a, b) if ca else (b, a)
    if cover.get((g, p), 0.0) < MIN_COVER:
        return None
    return bool(struct.get((g, p)))


def family_name(gs, used, prefix="S17"):
    """A family's name: prefix + its first coding gene's name (gs in genome order), alphanumerics only, and xN for the Nth family
    of that name (used counts them)."""
    first = next((g for g in gs if g.get("biotype", "protein_coding") == "protein_coding"), gs[0])
    base = re.sub(r"[^A-Za-z0-9]", "", first["gene_name"]) or "G"
    used[base] = used.get(base, 0) + 1
    return f"{prefix}{base}" + (f"x{used[base]}" if used[base] > 1 else "")


def spent_genes(paths):
    out = set()
    for p in paths:
        for fams in json.load(open(p)).get("sets", {}).values():
            for f in fams:
                out |= {c["gene_id"] for c in f["copies"]}
    return out


def splice_cmd(threads, retain=False):
    """The locus rule's spliced alignment (canonicals onto candidate loci). retain (--retain-chains, 2026-10-18): minimap2 -P,
    every chain kept -- among ~45 near-identical TSPY loci its default chain selection dropped the alignment of 62 of 63 TSPY
    transcripts to the most divergent copy's locus (TSPY2); -P kept 46."""
    return (["minimap2", "-cx", "splice:hq", "-N", "1000", "-p", "0.01", "--secondary=yes", "-t", str(threads)]
            + (["-P"] if retain else []))


def _paf_hits(paf, min_ident):
    """{(query, target): [(qs, qe)]} and {(query, target): (cover, PAF fields)} of a PAF at identity >= min_ident."""
    blocks, best = collections.defaultdict(list), {}
    for ln in paf.splitlines():
        f = ln.split("\t")
        a, b = f[0], f[5]
        if a == b or int(f[9]) < min_ident * int(f[10]):
            continue
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        blocks[(a, b)].append((qs, qe))
        cov = (qe - qs) / qlen
        if cov > best.get((a, b), (0.0, None))[0]:
            best[(a, b)] = (cov, f)
    return blocks, best


def cmd_build(out_p, spent_paths, threads=4, paralogs=False, chain="exact", retain=False):
    import pysam
    F = _families()
    gff = pysam.TabixFile(f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz")
    fa = pysam.FastaFile(F.GENOME)
    genes, others = [], []
    for r in csv.DictReader(open(f"{F.CAT}/chm13v2.0_CAT_Liftoff.genes.tsv"), delimiter="\t"):
        if int(r["end"]) - int(r["start0"]) > F.MAX_SPAN:
            continue
        g = dict(gene_id=r["gene_id"], gene_name=r["gene_name"], chrom=r["chrom"], start=int(r["start0"]), end=int(r["end"]),
                 strand=r["strand"], n=int(r["n_transcripts"]), biotype=r["gene_biotype"])
        (genes if r["gene_biotype"] == "protein_coding" else others).append(g)
    genes = F.drop_overlaps(genes)
    coding = {g["gene_id"] for g in genes}
    others = add_paralog_genes(genes, others) if paralogs else []
    for g in genes + others:
        g["tx"] = F.gene_transcripts(gff, g)
    genes, others = [g for g in genes if g["tx"]], [g for g in others if g["tx"]]
    by = {g["gene_id"]: g for g in genes + others}
    canon = {g["gene_id"]: g["tx"][F.canonical(g["tx"])] for g in genes + others}
    print(f"[seqfamilies] {len(genes)} protein-coding genes" + (f", {len(others)} paralog candidates of other biotypes"
                                                               if paralogs else ""), flush=True)

    def write_fa(path, ids, locus=False):
        with open(path, "w") as fh:
            for gid in sorted(ids):
                g = by[gid]
                sq = (fa.fetch(g["chrom"], max(0, g["start"] - PAD), g["end"] + PAD).upper() if locus
                      else F.spliced(fa, g["chrom"], g["strand"], canon[gid]))
                fh.write(f">{gid}\n{sq}\n")

    def candidates(paf):
        out = set()
        for ln in paf.splitlines():
            f = ln.split("\t")
            if f[0] != f[5] and int(f[9]) >= F.MIN_IDENT * int(f[10]) and int(f[10]) >= 0.5 * min(int(f[1]), int(f[6])):
                out.add((f[0], f[5]))
        return out
    splice = splice_cmd(threads, retain)
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as d:
        write_fa(f"{d}/c.fa", [g["gene_id"] for g in genes])
        # candidates: an all-vs-all of the canonicals at the family identity over half the shorter
        paf = subprocess.run(["minimap2", "-cx", "asm20", "-X", "-t", str(threads), f"{d}/c.fa", f"{d}/c.fa"],
                             capture_output=True, text=True, check=True).stdout
        cand = {x for pr in candidates(paf) for x in pr}
        print(f"[seqfamilies] {len(cand)} genes with a candidate partner", flush=True)
        # the locus rule on every candidate gene: its canonical onto every candidate locus
        write_fa(f"{d}/q.fa", cand)
        write_fa(f"{d}/t.fa", cand, locus=True)
        paf = subprocess.run(splice + [f"{d}/t.fa", f"{d}/q.fa"], capture_output=True, text=True, check=True).stdout
        blocks, best = _paf_hits(paf, F.MIN_IDENT)
        if paralogs:                                         # coding canonicals onto the loci of the paralog candidates
            write_fa(f"{d}/o.fa", [g["gene_id"] for g in others])
            pc = candidates(subprocess.run(["minimap2", "-cx", "asm20", "-t", str(threads), f"{d}/o.fa", f"{d}/c.fa"],
                                           capture_output=True, text=True, check=True).stdout)
            pq, pt = {a for a, _ in pc}, {b for _, b in pc}
            print(f"[seqfamilies] {len(pt)} paralogs of other biotypes with a candidate coding gene ({len(pq)})", flush=True)
            if pq:
                write_fa(f"{d}/pq.fa", pq)
                write_fa(f"{d}/pt.fa", pt, locus=True)
                b2, best2 = _paf_hits(subprocess.run(splice + [f"{d}/pt.fa", f"{d}/pq.fa"], capture_output=True, text=True,
                                                     check=True).stdout, F.MIN_IDENT)
                blocks.update(b2)
                best.update(best2)
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
        return structure_match(js, qlen, f[4], int(f[7]), [(int(n), op) for n, op in _CIGAR.findall(cg)], q0=q0, chain=chain)
    pairs = {tuple(sorted(k)) for k in cover}
    st = {k: struct(*k) for k in cover}
    links, rows = [], []
    for a, b in sorted(pairs):
        ok = pair_link(a, b, cover, st, coding)
        if ok is None:
            continue
        links.append((a, b, ok))
        rows.append([a, b, by[a]["gene_name"], by[b]["gene_name"], f"{cover.get((a, b), 0.0):.3f}", f"{cover.get((b, a), 0.0):.3f}",
                     int(st.get((a, b), False)), int(st.get((b, a), False)), int(ok)]
                    + ([by[a]["biotype"], by[b]["biotype"]] if paralogs else []))
    spent = spent_genes(spent_paths)
    comps = components(links, spent=spent)
    n_seq = len(components([(a, b, True) for a, b, _ in links]))
    fams, used = [], collections.Counter()
    for comp in comps:
        gs = sorted((by[x] for x in comp), key=lambda g: (g["chrom"], g["start"]))
        name = family_name(gs, used, "S19" if paralogs else "S17")
        copies = []
        for g in gs:
            main, iso = F.pick_isoforms(g["tx"])
            copies.append(dict(gene_id=g["gene_id"], gene_name=g["gene_name"], chrom=g["chrom"],
                               start=min(s for v in iso.values() for s, _ in v), end=max(e for v in iso.values() for _, e in v),
                               strand=g["strand"], main=main, transcripts={t: [list(x) for x in v] for t, v in iso.items()},
                               **({"biotype": g["biotype"]} if paralogs else {})))
        fams.append(dict(family=name, candidates=len(comp), component=len(comp), copies=copies))
    with open(f"{out_p}.links.tsv", "w") as fh:
        fh.write("a\tb\ta_name\tb_name\tcover_ab\tcover_ba\tstruct_ab\tstruct_ba\tstruct_link"
                 + ("\ta_biotype\tb_biotype" if paralogs else "") + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    json.dump(dict(genome=F.GENOME, annotation=f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz",
                   rule=("seqfamilies --paralogs" if paralogs else "seqfamilies") + (" --intron-chain carried" if chain == "carried" else "")
                   + (" --retain-chains" if retain else ""),
                   min_ident=F.MIN_IDENT, min_cover=MIN_COVER, junction_tol=JUNCTION_TOL,
                   **({"intron_chain": chain} if chain != "exact" else {}), **({"retain_chains": True} if retain else {}),
                   sets={"seq19" if paralogs else "seq17": fams}),
              open(out_p, "w"), indent=1)
    unnamed = sum(any(re.match(r"^(LOC|ENSG|MSTRG|LINC|AC[0-9]|AL[0-9]|AP[0-9])", c["gene_name"]) or "-" in c["gene_name"]
                      for c in f["copies"]) for f in fams)
    n_par = sum(1 for f in fams for c in f["copies"] if c.get("biotype", "protein_coding") != "protein_coding")
    print(f"[seqfamilies] sequence links {len(links)}, of them also structure {sum(ok for _, _, ok in links)}; "
          f"components by sequence {n_seq}, by sequence and structure {len(components(links))}, of them unspent {len(fams)} "
          f"({sum(len(f['copies']) for f in fams)} copies, {n_par} of them paralogs of other biotypes; {unnamed} with an unnamed "
          f"copy) -> {out_p}", flush=True)
    return 0


def cmd_draw(out_p, all_p, seed=SEED, set_name="r17"):
    d = json.load(open(all_p))
    (pool,) = d["sets"].values()                          # a pool holds one set (seq17, or seq19 under --paralogs)
    got = draw(pool, seed=seed)
    json.dump(dict({k: v for k, v in d.items() if k != "sets"}, draw=dict(seed=seed, budget=BUDGET, max_total=MAX_TOTAL),
                   sets={set_name: got}), open(out_p, "w"), indent=1)
    print(f"[seqfamilies] drew {len(got)} families, {sum(len(f['copies']) for f in got)} copies: "
          + ", ".join(f"{f['family']} ({len(f['copies'])})" for f in got), flush=True)
    return 0


def build_args(argv):
    """`build OUT SPENT... [--paralogs] [--intron-chain exact|carried] [--retain-chains]` (arguments after `build`) -> (OUT,
    SPENT paths, paralogs, intron-chain rule, retain chains)."""
    rest, paralogs, chain, retain = [], False, "exact", False
    it = iter(argv)
    for x in it:
        if x == "--paralogs":
            paralogs = True
        elif x == "--retain-chains":
            retain = True
        elif x == "--intron-chain":
            chain = next(it, "")
        else:
            rest.append(x)
    if chain not in ("exact", "carried") or not rest:
        raise SystemExit(__doc__)
    return rest[0], rest[1:], paralogs, chain, retain


def main(argv):
    if argv[0] == "build":
        out, spent, paralogs, chain, retain = build_args(argv[1:])
        return cmd_build(out, spent, paralogs=paralogs, chain=chain, retain=retain)
    if argv[0] == "draw":
        return cmd_draw(argv[1], argv[2], seed=int(argv[argv.index("--seed") + 1]) if "--seed" in argv else SEED,
                        set_name=argv[argv.index("--set") + 1] if "--set" in argv else "r17")
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
