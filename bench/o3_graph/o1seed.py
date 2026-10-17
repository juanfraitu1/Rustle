#!/usr/bin/env python3
"""O1 seeding (round 7, docs/PREREG_o3_graph_round7_2026-10-16.md): a family's copy catalog from RNA support, before O2 and O3.

    o1seed.py WORLD    WORLD/graph/seeded.json

Per family (its curated copies, isoforms.json):
  family reads  reads with any alignment (primary or secondary) on a curated copy's locus (+- PAD) in WORLD/E0.bam
  support       a curated copy is KEPT when >= MIN_SUPPORT family reads have their PRIMARY alignment on its locus, or when one
                family read is placed on it uniquely (mapq >= MIN_MAPQ_NEW) and owned by it (round 9 review: realistic
                expression leaves expressed copies with 4-5 reads, and a coin-flip count alone dropped them); else it is
                DROPPED (nothing in the data says it is expressed: the universe O2 chooses from shrinks)
  new copies    family reads' primary alignments off every curated locus of every family, grouped by overlap (gaps <= JOIN),
                groups with >= MIN_SUPPORT reads placed there uniquely (mapq >= MIN_MAPQ_NEW): a locus the family's reads come from that the
                catalog lacks (a mapq-0 primary is a coin flip, no evidence for that locus; keeping a curated copy counts every
                primary, so an identical expressed copy is never dropped)
  structures    a kept copy: its curated isoforms plus the intron chains of near-uniquely placed primary alignments (mapq >=
                MIN_MAPQ, the user's suggestion: no fragments of ambiguous placement) that this copy overlaps most and that lie
                within its span (+- SPAN_SLACK), in >= MIN_CHAIN reads (ends: the reads' outermost ends); a new copy: its observed
                chains only
  one family    a new locus seeded in several families (adjacent genes: WASH.7 lies 1 kb from DDX11L.6, and WASH and DDX11L copies
                are neighbours in every subtelomere, so WASH.7's reads are DDX11L family reads too) stays only in the family whose
                curated copies its reads align to most (aligned bases, every alignment of every read that seeded it, unpadded copy
                spans); a tie seeds it in none
Output per family: copies [{idx, contig, start, end, strand, rep, isoforms, origin: curated | new, support}], dropped [curated
idx], family_reads [read names]; a new copy also lists the reads that seeded it. Seeded copies keep curated idx; new ones are
numbered after them.
"""
import collections
import json
import os
import sys

PAD, MIN_SUPPORT, MIN_CHAIN, JOIN, CHAIN_TOL, MIN_MAPQ, MIN_MAPQ_NEW, SPAN_SLACK = 1000, 5, 3, 1000, 10, 50, 1, 100
DELTA, TWIN_COVER = 0.00958, 0.5                                  # --twins: O3's delta; the alignment must cover half the locus


def twin_of(new_seq, curated):
    """--twins (round 10 development): (idx, de) of the curated copy whose genomic span the new copy's span aligns to (minimap2
    -x asm20, >= TWIN_COVER of it) with the lowest gap-compressed divergence `de`, when that is <= DELTA (a copy the reads cannot
    tell from it at O3's resolution), else (None, the lowest de or None). curated: {idx: sequence}."""
    import subprocess
    import tempfile
    if not curated:
        return None, None
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/t.fa", "w") as fh:
            fh.writelines(f">{i}\n{x}\n" for i, x in curated.items())
        open(f"{d}/q.fa", "w").write(f">new\n{new_seq}\n")
        out = subprocess.run(["minimap2", "-c", "-x", "asm20", "--secondary=yes", "-N", "200", "-p", "0.01", f"{d}/t.fa",
                              f"{d}/q.fa"], capture_output=True, text=True, check=True).stdout
    best = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if int(f[3]) - int(f[2]) < TWIN_COVER * int(f[1]):
            continue
        de = float(next(x for x in f[12:] if x.startswith("de:f:"))[5:])
        i = int(f[5])
        best[i] = min(de, best.get(i, 1.0))
    if not best:
        return None, None
    i, de = min(best.items(), key=lambda kv: (kv[1], kv[0]))
    return (i if de <= DELTA else None), de


def exons_from_cigartuples(ref_start, cigartuples):
    """Reference exons of an alignment: a new exon at every N; M/=/X/D extend the exon, I/S/H take no reference."""
    ex, pos, start = [], ref_start, ref_start
    for op, n in cigartuples:
        if op == 3:                                       # N: intron
            ex.append((start, pos))
            pos += n
            start = pos
        elif op in (0, 2, 7, 8):                          # M, D, =, X
            pos += n
    ex.append((start, pos))
    return [(a, b) for a, b in ex if b > a]


def _span(a):
    return a["blocks"][0][0], a["blocks"][-1][1]


def _overlaps(a, locus):
    s, e = _span(a)
    return a["chrom"] == locus[0] and s < locus[2] and e > locus[1]


def support(prim, loci):
    """{copy idx: number of primary alignments overlapping its locus (chrom, start, end)}."""
    return {i: sum(1 for a in prim if _overlaps(a, lo)) for i, lo in loci.items()}


def kept(sup, min_support=MIN_SUPPORT):
    return sorted(i for i, n in sup.items() if n >= min_support)


def new_loci(prim, loci, min_support=MIN_SUPPORT, join=JOIN):
    """Groups of primary alignments off every locus (overlapping or within `join`) with >= min_support reads placed there
    UNIQUELY (mapq >= MIN_MAPQ_NEW = 1): adding a copy needs positive evidence; a mapq-0 primary is a coin flip. Not mapq >= 50:
    a real uncatalogued copy near-identical to a paralog has unique placements at intermediate mapq (development: ANKRD20A.4,
    SPDYE.2, TRIM64.2, DEFA.1, FRG1.2, CSAG.3)."""
    off = sorted((a for a in prim if not any(_overlaps(a, lo) for lo in loci.values())), key=lambda a: (a["chrom"], _span(a)))
    groups, cur = [], None
    for a in off:
        s, e = _span(a)
        if cur and a["chrom"] == cur["chrom"] and s <= cur["end"] + join:
            cur["end"] = max(cur["end"], e)
            cur["alns"].append(a)
        else:
            cur = dict(chrom=a["chrom"], start=s, end=e, alns=[a])
            groups.append(cur)
    out = []
    for g in groups:
        reads = {a["read"] for a in g["alns"]}
        unique = {a["read"] for a in g["alns"] if a.get("mapq", 60) >= MIN_MAPQ_NEW}
        if len(unique) >= min_support:
            strands = collections.Counter(a["strand"] for a in g["alns"])
            out.append(dict(chrom=g["chrom"], start=g["start"], end=g["end"], strand=strands.most_common(1)[0][0],
                            n_reads=len(reads), alns=g["alns"]))
    return out


def observed_chains(prim, min_reads=MIN_CHAIN):
    """[(exons, n reads)] for each intron chain >= min_reads primary alignments share; ends are the reads' outermost ends."""
    by = collections.defaultdict(list)
    for a in prim:
        b = a["blocks"]
        chain = tuple((x[1], y[0]) for x, y in zip(b, b[1:]))
        if chain:
            by[chain].append(a)
    out = []
    for chain, alns in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        if len(alns) < min_reads:
            continue
        s, e = min(a["blocks"][0][0] for a in alns), max(a["blocks"][-1][1] for a in alns)
        bounds = [s] + [x for ab in chain for x in ab] + [e]
        out.append(([(bounds[2 * k], bounds[2 * k + 1]) for k in range(len(chain) + 1)], len(alns)))
    return out


def _chain(ex):
    ex = sorted(tuple(x) for x in ex)
    return tuple((a[1], b[0]) for a, b in zip(ex, ex[1:]))


def _known(chain, chains, tol=CHAIN_TOL):
    return any(len(k) == len(chain) and all(abs(x - y) <= tol for i, j in zip(chain, k) for x, y in zip(i, j)) for k in chains)


def owner(a, cps):
    """The one curated copy an alignment overlaps most (unpadded spans), or None: adjacent copies' padded windows overlap
    (CSAG copies 91 bp apart), so a read must count for one copy only."""
    s, e = _span(a)
    best, ov = None, 0
    for c in cps:
        if c["contig"] == a["chrom"]:
            o = min(e, c["end"]) - max(s, c["start"])
            if o > ov:
                best, ov = c["idx"], o
    return best


def seed_family(cps, prim, all_loci, min_support=MIN_SUPPORT, min_chain=MIN_CHAIN):
    """cps: the family's curated copies (isoforms.json entries); prim: its family reads' primary alignments; all_loci: every
    family's curated loci (padded), so another family's copy is never seeded here. -> (seeded copies, dropped idx)."""
    loci = {c["idx"]: (c["contig"], max(0, c["start"] - PAD), c["end"] + PAD) for c in cps}
    sup = support(prim, loci)
    unique = collections.Counter(owner(a, cps) for a in prim if a.get("mapq", 60) >= MIN_MAPQ_NEW)
    keep = set(kept(sup, min_support)) | {i for i in loci if unique[i] >= 1}   # (round 9 review) one read placed there
    # uniquely is evidence the copy is expressed: realistic expression leaves expressed copies under the count
    copies = []
    for c in cps:
        if c["idx"] not in keep:
            continue
        own = [a for a in prim if a.get("mapq", 60) >= MIN_MAPQ and owner(a, cps) == c["idx"]       # near-unique, this copy's,
               and _span(a)[0] >= c["start"] - SPAN_SLACK and _span(a)[1] <= c["end"] + SPAN_SLACK]  # and inside its span
        iso = {t: [list(x) for x in ex] for t, ex in c["isoforms"].items()}
        seen = [_chain(ex) for ex in iso.values()]
        for k, (ex, n) in enumerate(observed_chains(own, min_chain)):
            if not _known(_chain(ex), seen):
                seen.append(_chain(ex))
                iso[f"obs{k}"] = [list(x) for x in ex]
        copies.append(dict(idx=c["idx"], contig=c["contig"], start=c["start"], end=c["end"], strand=c["strand"], rep=c["rep"],
                           isoforms=iso, origin="curated", support=sup[c["idx"]]))
    nxt = max((c["idx"] for c in cps), default=-1) + 1
    for g in new_loci(prim, all_loci, min_support):
        chains = observed_chains([a for a in g["alns"] if a.get("mapq", 60) >= MIN_MAPQ_NEW], min_chain)
        if not chains:
            continue
        iso = {f"obs{k}": [list(x) for x in ex] for k, (ex, n) in enumerate(chains)}
        copies.append(dict(idx=nxt, contig=g["chrom"], start=min(x[0] for v in iso.values() for x in v),
                           end=max(x[1] for v in iso.values() for x in v), strand=g["strand"], rep="obs0", isoforms=iso,
                           origin="new", support=g["n_reads"], reads=sorted({a["read"] for a in g["alns"]})))
        nxt += 1
    return copies, sorted(set(loci) - keep)


def affinity(alns, cps):
    """{read: aligned bases of its alignments on the family's curated copies}: unpadded spans, merged so a base counts once."""
    spans = collections.defaultdict(list)
    for c in sorted(cps, key=lambda c: (c["contig"], c["start"])):
        sp = spans[c["contig"]]
        if sp and c["start"] <= sp[-1][1]:
            sp[-1][1] = max(sp[-1][1], c["end"])
        else:
            sp.append([c["start"], c["end"]])
    out = collections.Counter()
    for a in alns:
        for cs, ce in spans.get(a["chrom"], ()):
            out[a["read"]] += sum(max(0, min(e, ce) - max(s, cs)) for s, e in a["blocks"])
    return {q: n for q, n in out.items() if n > 0}


def one_family(seeded, aff):
    """A new locus seeded in several families (same contig, overlapping spans) stays only in the family whose curated copies
    the reads that seeded it align to most (aff: {family: affinity()}); a tie, also no alignment to any family's copies, seeds
    it in none. Curated copies and loci seeded in one family only are left alone."""
    new = sorted(((c["contig"], c["start"], c["end"], f, c["idx"], tuple(c["reads"])) for f, v in seeded.items()
                  for c in v["copies"] if c["origin"] == "new"))
    groups, end = [], None
    for x in new:
        if groups and x[0] == groups[-1][-1][0] and x[1] < end:
            groups[-1].append(x)
            end = max(end, x[2])
        else:
            groups.append([x])
            end = x[2]
    drop = set()
    for g in groups:
        fams = sorted({x[3] for x in g})
        if len(fams) < 2:
            continue
        reads = set().union(*(x[5] for x in g))
        score = {f: sum(aff.get(f, {}).get(q, 0) for q in reads) for f in fams}
        top = [f for f in fams if score[f] == max(score.values())]
        drop |= {(x[3], x[4]) for x in g if len(top) > 1 or x[3] != top[0]}
    return {f: dict(v, copies=[c for c in v["copies"] if (f, c["idx"]) not in drop]) for f, v in seeded.items()}


def seed_world(w, twins=False):
    import pysam
    iso = json.load(open(f"{w}/isoforms.json"))
    all_loci = {(f, c["idx"]): (c["contig"], max(0, c["start"] - PAD), c["end"] + PAD) for f, cps in iso.items() for c in cps}
    fam_reads = collections.defaultdict(set)
    seen, aff = collections.defaultdict(set), collections.defaultdict(collections.Counter)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for (f, _), (ch, s, e) in all_loci.items():
            for rec in bam.fetch(ch, s, e):
                if not rec.is_unmapped:
                    fam_reads[f].add(rec.query_name)
                    key = (rec.query_name, ch, rec.reference_start, rec.flag)
                    if key not in seen[f]:                            # an alignment in two padded windows counts once
                        seen[f].add(key)
                        aln = dict(read=rec.query_name, chrom=ch,
                                   blocks=exons_from_cigartuples(rec.reference_start, rec.cigartuples))
                        aff[f].update(affinity([aln], iso[f]))
    want = set().union(*fam_reads.values()) if fam_reads else set()
    prim = {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in want and not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
                prim[rec.query_name] = dict(read=rec.query_name, chrom=rec.reference_name, strand="-" if rec.is_reverse else "+",
                                            mapq=rec.mapping_quality,
                                            blocks=exons_from_cigartuples(rec.reference_start, rec.cigartuples))
    out = {}
    for f, cps in sorted(iso.items()):
        p = [prim[q] for q in sorted(fam_reads[f]) if q in prim]
        copies, dropped = seed_family(cps, p, all_loci)
        out[f] = dict(copies=copies, dropped=dropped, family_reads=sorted(fam_reads[f]))
    out = one_family(out, aff)
    if twins:                                                     # a new copy within delta of a kept curated copy is its twin
        fa = pysam.FastaFile(f"{w}/genome.ref.fa")
        for f, v in out.items():
            cur = {c["idx"]: fa.fetch(c["contig"], c["start"], c["end"]).upper() for c in v["copies"] if c["origin"] == "curated"}
            for c in v["copies"]:
                if c["origin"] == "new":
                    c["twin_of"], c["twin_de"] = twin_of(fa.fetch(c["contig"], c["start"], c["end"]).upper(), cur)
    os.makedirs(f"{w}/graph", exist_ok=True)                      # a fresh world: o1seed runs first (round 8)
    json.dump(out, open(f"{w}/graph/seeded.json", "w"))
    n_new = sum(1 for v in out.values() for c in v["copies"] if c["origin"] == "new")
    n_drop = sum(len(v["dropped"]) for v in out.values())
    print(f"[o1seed] {w}/graph/seeded.json: {sum(len(v['copies']) for v in out.values())} copies ({n_new} new), "
          f"{n_drop} curated copies dropped", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(seed_world(sys.argv[1], twins="--twins" in sys.argv))
