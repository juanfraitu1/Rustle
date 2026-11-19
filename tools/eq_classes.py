#!/usr/bin/env python3
"""Equivalence classes of reads over the loci of an assembly, with per-locus count ranges (driver stage `eqclasses`).

    eq_classes.py GTF BAM AS_TABLE OUT_PREFIX [--score ms|AS] [--window W | --ratio R] [--threads N] [--exclude BED] [--contigs A,B] [--famwise FW_PREFIX]

A read's placements = the loci (gene_id of GTF, on the record's strand, most exon bases shared) of its records whose
score ties its best (primary or secondary alike). The score is minimap2's `ms` (max-scoring segment, the one it ranks
records on; default since 2026-11-16) or, with --score AS, the alignment score: AS charges for introns, so it ties a
spliced gene with its processed pseudogenes (KB3781: 80,611 reads FTH1 / LOC115931247), ms does not.
Two bands, side by side, since neither is right everywhere
(KB3781 2026-11-16: exact ties misplace reads of the other haplotype's allele, the wide band blurs pairs the aligner
separates): `exact` = score equal to the best, `wide` = score within W of the best (--window, default 10: the
project's standing ms near-tie window; a few mismatches whatever the read length) or, with --ratio R, score >= R x best
(0.98 was the first default: ~40 on a long read, it re-ties near-identical pseudogenes). Wide columns carry the suffix _w.
The primary record is a placement in both bands whatever its score (the aligner may prefer a
record that is not the best-scoring one; such a read is then tied between the two). A placement on no locus is the member `outside`. A read placed on one locus is `placed`; the others
form classes = sets of loci the band does not separate.
Per locus:  placed  reads whose best-scoring placement(s) are all here. A lower bound on what aligns best here, NOT a
                    proof of origin: a read from a copy the reference lacks is still placed on its nearest copy.
            max     placed + every read of every class holding the locus
            em      counts-only EM over the classes (abundance prior; quantification only, no read is assigned).
                    Where a tie group has no placed read (`group_placed` 0) em is the starting equal split.
Tie groups = components of loci sharing >= 2 class reads; loci mostly inside --exclude BED are `repeat` and link nothing.
Writes OUT_PREFIX.classes.tsv (column band) and OUT_PREFIX.loci.tsv; prints one summary line per band.
--famwise FW_PREFIX (driver stage `famwise`; the two tables above stay byte-identical) also writes, for tools/famwise.py,
FW_PREFIX.reads.tsv (every record within 20 of the best, with ms and AS, of each molecule tied that closely) and
FW_PREFIX.placed.tsv (the other molecules per locus and transcript of their primary record)."""
import argparse, collections, multiprocessing, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_names import overlap_genes, opn


def tx_strand(r):
    """The transcript's strand of a read: minimap2's ts:A tag says whether the transcript runs against the read (ts:A:-), so a
    reverse-complemented cDNA read counts on its gene's strand (audit 2026-11-19). Without the tag: the alignment strand."""
    rev = r.is_reverse
    if r.has_tag("ts") and r.get_tag("ts") == "-":
        rev = not rev
    return "-" if rev else "+"

OUTSIDE = "outside"
G_NODES = None  # {locus: {intron chain: node}} with --famwise

def load_loci(gtf):
    """-> meta {locus: [chrom, start, end, strand, name]}, gex {(chrom, strand): sorted [(s, e, locus)]} (1-based closed)."""
    ex, loc, names = collections.defaultdict(set), {}, collections.defaultdict(collections.Counter)
    for l in opn(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "transcript"):
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        m = re.search(r'ref_gene_name "([^"]+)"', f[8])
        if m and f[2] == "transcript":
            names[g][m.group(1)] += 1
        if f[2] == "exon":
            ex[g].add((int(f[3]), int(f[4]))); loc[g] = (f[0], f[6])
    meta, gex = {}, collections.defaultdict(list)
    for g, e in ex.items():
        c, st = loc[g]
        meta[g] = [c, min(s for s, _ in e), max(x for _, x in e), st, names[g].most_common(1)[0][0] if names[g] else "."]
        for s in ("+", "-") if st not in "+-" else (st,):
            gex[(c, s)].extend((a, b, g) for a, b in e)
    return meta, {k: sorted(v) for k, v in gex.items()}

COLS = {"AS": (1, 2, 5), "ms": (6, 7, 8)}  # as_table fields: best, second, primary's score

def floor_fn(window=10, ratio=None):
    """lowest score the wide band admits, given the read's best."""
    return (lambda b: ratio * b) if ratio is not None else (lambda b: b - window)

def load_ties(as_table, floor, score="ms"):
    """{read: best score} for the reads a second record can tie (second >= floor(best)) or whose primary is not their best record."""
    out, (ib, i2, ip) = {}, COLS[score]
    for l in open(as_table):
        if l[0] != "#":
            f = l.rstrip("\n").split("\t")
            if len(f) <= ip:
                sys.exit(f"eq_classes: {as_table} has no {score} columns ({len(f)} fields); rebuild it with as_table or pass --score AS")
            b, s = int(f[ib]), int(f[i2])
            if b > 0 and (s >= floor(b) or f[ip] == "NA" or int(f[ip]) < b):  # NA: the BAM holds only secondaries of the molecule
                out[f[0]] = b
    return out

def scan(contig):
    import pysam
    placed, pairs, members = collections.Counter(), [], collections.Counter()
    if G_NODES is not None:
        import famwise
    idx = {s: ([a for a, _, _ in G_GEX[(contig, s)]], max(b - a for a, b, _ in G_GEX[(contig, s)])) for s in "+-" if (contig, s) in G_GEX}
    for r in pysam.AlignmentFile(G_BAM).fetch(contig):
        if r.is_unmapped or r.is_supplementary:
            continue
        best = G_TIES.get(r.query_name)
        if best is None:
            if r.is_secondary:
                continue
        elif r.is_secondary and r.get_tag(G_TAG) < G_FLOOR(best):
            continue
        s = tx_strand(r)
        ov = overlap_genes([(a + 1, b) for a, b in r.get_blocks()], G_GEX[(contig, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else OUTSIDE
        m = None if G_NODES is None else famwise.member_of(g, famwise.record_chain(r.reference_start, r.cigartuples), G_NODES)
        if best is None:
            placed[g] += 1
            if m is not None:
                members[(g, m)] += 1
        elif m is None:
            pairs.append((r.query_name, g, not r.is_secondary or r.get_tag(G_TAG) == best))
        else:
            t = r.get_tag(G_TAG)
            pairs.append((r.query_name, g, not r.is_secondary or t == best, not r.is_secondary or t >= G_WIDE(best), m, best - t,
                          *(r.get_tag(x) if r.has_tag(x) else "NA" for x in ("ms", "AS")), int(not r.is_secondary)))
    return placed, pairs, members

def em(classes, placed, rounds=500, tol=1e-6):
    """counts-only EM: theta = placed + expected class reads; classes {frozenset: n}."""
    members = {l for c in classes for l in c}
    theta = {l: placed.get(l, 0) + sum(n / len(c) for c, n in classes.items() if l in c) for l in members}
    for _ in range(rounds):
        new = {l: float(placed.get(l, 0)) for l in members}
        for c, n in classes.items():
            z = sum(theta[l] for l in c)
            for l in c:
                new[l] += n * (theta[l] / z if z > 0 else 1 / len(c))
        d = max(abs(new[l] - theta[l]) for l in members) if members else 0
        theta = new
        if d < tol:
            break
    return theta

def summarise(placed, classes, meta, repeat=frozenset()):
    """-> loci rows, class rows. placed {locus: n}, classes {frozenset(loci): n}."""
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for c, n in classes.items():
        m = sorted(l for l in c if l != OUTSIDE and l not in repeat)
        if n >= 2:
            for l in m[1:]:
                parent[find(l)] = find(m[0])
    comp = collections.defaultdict(list)
    for l in list(parent):
        comp[find(l)].append(l)
    gid = {l: f"G{i + 1}" for i, (_, ls) in enumerate(sorted(comp.items(), key=lambda kv: (-len(kv[1]), min(kv[1])))) if len(ls) > 1 for l in ls}
    theta = em(classes, placed)
    extra = collections.Counter()
    for c, n in classes.items():
        for l in c:
            extra[l] += n
    total = sum(placed.values()) + sum(classes.values())
    assert abs(sum(theta.values()) + sum(n for l, n in placed.items() if l not in theta) - total) < 1e-3 * max(1, total), "EM lost reads"
    gsize, gplaced, gchrom = collections.Counter(gid.values()), collections.Counter(), collections.defaultdict(set)
    for l, g in gid.items():
        gplaced[g] += placed.get(l, 0); gchrom[g].add(meta[l][0])
    loci = []
    for l in sorted(meta, key=lambda l: (meta[l][0], meta[l][1], l)):   # locus name breaks start ties (audit 2026-11-19: run order)
        lo, hi = placed.get(l, 0), placed.get(l, 0) + extra[l]
        if hi == 0:
            continue
        e = theta.get(l, float(lo))
        assert lo - 1e-6 <= e <= hi + 1e-6, (l, lo, e, hi)
        g = gid.get(l, ".")
        label = "repeat" if l in repeat else "placed" if hi == lo else "tied_only" if lo == 0 else "mixed"
        loci.append([l, *meta[l][:4], meta[l][4], label, lo, hi, f"{e:.1f}", g, gsize.get(g, 0), len(gchrom.get(g, ())), gplaced.get(g, 0)])
    cls = [[n, len(c), ",".join(sorted(c)), next((gid[l] for l in sorted(c) if l in gid), ".")] for c, n in sorted(classes.items(), key=lambda kv: (-kv[1], sorted(kv[0])))]
    return loci, cls

def repeat_loci(meta, gex, bed):
    from gtf_exclude_bed import load_bed, covered  # merged intervals: nested BED elements used to count a base twice
    b, cov, size = load_bed(bed), collections.Counter(), collections.Counter()
    for (c, _), exs in gex.items():
        for s, e, g in exs:
            cov[g] += covered(b, c, s, e); size[g] += e - s + 1
    return frozenset(g for g in cov if cov[g] and cov[g] >= 0.5 * size[g])

def main():
    global G_GEX, G_BAM, G_TIES, G_FLOOR, G_TAG, G_NODES, G_WIDE
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gtf"); p.add_argument("bam"); p.add_argument("as_table"); p.add_argument("out")
    p.add_argument("--score", choices=("ms", "AS"), default="ms"); p.add_argument("--window", type=float, default=10); p.add_argument("--ratio", type=float); p.add_argument("--threads", type=int, default=4)
    p.add_argument("--exclude"); p.add_argument("--contigs"); p.add_argument("--famwise")
    a = p.parse_args()
    import pysam
    meta, G_GEX = load_loci(a.gtf)
    G_FLOOR = floor_fn(a.window, a.ratio)
    if a.famwise:
        if a.ratio is not None:
            p.error("--famwise needs the window band (no --ratio)")
        import famwise
        G_NODES, G_WIDE = famwise.load_nodes(a.gtf)[0], G_FLOOR
        G_FLOOR = floor_fn(max(a.window, famwise.FW_WINDOW))  # records down to here are written; the bands below still stop at --window
    G_BAM, G_TAG, G_TIES = a.bam, a.score, load_ties(a.as_table, G_FLOOR, a.score)
    contigs = a.contigs.split(",") if a.contigs else list(pysam.AlignmentFile(a.bam).references)
    with multiprocessing.get_context("fork").Pool(a.threads) as pool:
        res = pool.map(scan, contigs, chunksize=1)
    base, sets = collections.Counter(), {"exact": collections.defaultdict(set), "wide": collections.defaultdict(set)}
    recs, members = collections.defaultdict(list), collections.Counter()
    for pl, pairs, mem in res:
        base.update(pl); members.update(mem)
        for q, g, exact, *x in pairs:
            if not x or x[0]:
                sets["wide"][q].add(g)
            if exact:
                sets["exact"][q].add(g)
            if x:
                recs[q].append((g, *x[1:]))
    base.pop(OUTSIDE, None)
    repeat = repeat_loci(meta, G_GEX, a.exclude) if a.exclude else frozenset()
    rows, cls_all = {}, []
    for band in ("exact", "wide"):
        placed, classes = collections.Counter(base), collections.Counter()
        for s in sets[band].values():
            if len(s) == 1:
                placed.update(s - {OUTSIDE})  # its rival record lies on the same locus
            else:
                classes[frozenset(s)] += 1
        loci, cls = summarise(placed, classes, meta, repeat)
        rows[band] = {r[0]: r for r in loci}
        cls_all += [[band] + r for r in cls]
        lab = collections.Counter(r[6] for r in loci)
        big = max(loci, key=lambda r: r[11], default=None)
        print(f"{a.score} {band}{'' if band == 'exact' else (f' (>= {a.ratio} x best)' if a.ratio is not None else f' (within {a.window:g} of best)')}: {sum(placed.values())} reads placed, {sum(classes.values())} in {len(classes)} classes; loci " +
              ", ".join(f"{k} {v}" for k, v in sorted(lab.items())) + f"; {len({r[10] for r in loci} - {'.'})} tie groups" +
              (f", largest {big[11]} loci on {big[12]} contig(s)" if big and big[11] else ""))
    cols = "label\tplaced\tmax\tem\tgroup\tgroup_loci\tgroup_chroms\tgroup_placed".split("\t")
    with open(a.out + ".loci.tsv", "w") as fh:
        fh.write("locus\tchrom\tstart\tend\tstrand\tref_name\t" + "\t".join(cols) + "\t" + "\t".join(c + "_w" for c in cols) + "\n")
        for l in sorted(set(rows["exact"]) | set(rows["wide"]), key=lambda l: (meta[l][0], meta[l][1], l)):
            e, w = rows["exact"].get(l), rows["wide"].get(l)
            none = ["none", 0, 0, "0.0", ".", 0, 0, 0]  # a locus only weaker-than-best records reach has no exact row
            fh.write("\t".join(map(str, (e or w)[:6] + (e[6:] if e else none) + (w[6:] if w else none))) + "\n")
    with open(a.out + ".classes.tsv", "w") as fh:
        fh.write("band\treads\tsize\tloci\tgroup\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in cls_all)
    if a.famwise:
        _, i2, ip = COLS[a.score]
        second = {f[0]: (int(f[i2]), f[ip]) for f in (l.rstrip("\n").split("\t") for l in open(a.as_table) if l[0] != "#") if f[0] in recs}
        with open(a.famwise + ".reads.tsv", "w") as fh:
            fh.write("read\tlocus\tmember\tgap\tms\tAS\tprimary\n")
            for q in sorted(recs):
                rs, b, (s2, pr) = list(recs[q]), G_TIES[q], second[q]
                # records the scan never met (--contigs panel) get an `unseen` row: without it the molecule would read as placed.
                # The primary (a band member whatever its score), the best record, the second-best when it is within reach.
                # shortcut: the table holds no third score, so a panel of 3-copy ties should scan all three contigs
                unseen = lambda score, prim: rs.append((famwise.UNSEEN, ".", b - score, *((score, "NA") if a.score == "ms" else ("NA", score)), prim))
                if pr != "NA" and not any(r[-1] for r in rs):
                    unseen(int(pr), 1)
                if not any(r[2] == 0 for r in rs):
                    unseen(b, 0)
                if s2 >= G_FLOOR(b) and sum(b - r[2] >= s2 for r in rs) < 2:
                    unseen(s2, 0)
                fh.writelines("\t".join(map(str, (q, *r))) + "\n" for r in sorted(rs, key=lambda r: (r[2], r[0], r[1], str(r[3:]))))
        with open(a.famwise + ".placed.tsv", "w") as fh:
            fh.write("locus\tmember\tn\n")
            fh.writelines(f"{l}\t{m}\t{n}\n" for (l, m), n in sorted(members.items()))

if __name__ == "__main__":
    main()
