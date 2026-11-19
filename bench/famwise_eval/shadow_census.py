#!/usr/bin/env python3
"""Places reached only by secondary alignments: real missed loci, or shadows? (user / advisor, 2026-11-19: even an unexpressed
gene's paralogs or pseudogenes should show up through secondaries, and some could help seeding; are any loci covered only by
multi-mappers real?)

  shadow_census.py --truth READS.tsv --tie-split T.tsv --bam SEL.bam --out PREFIX [--window 10] [--min-molecules 2]
                   [--annot CAT.gff3.gz] [--gtf NAME=GTF ...]

Secondary-only place: records within --window ms of their molecule's best, merged on contig + strand + overlap, holding
>= --min-molecules molecules and NO primary record of any molecule. Each place is classed by the molecules' truth:
  REAL_DECIDED  >= 1 molecule whose single origin locus is this place (its primary went elsewhere: a missed real locus)
  REAL_SHARED   else >= 1 molecule whose origin set (SEVERAL COPIES tie) holds this place (real, but not separable)
  SHADOW        else >= 1 judged molecule, every origin elsewhere (the place only looks like the expressed copy)
  UNJUDGED      no judged molecule (outside / unjudged / one-place ties only)
With --annot each place gets the gene it overlaps (same strand) and the biotype group; with --gtf, whether an assembly has a
multi-exon transcript there. Writes PREFIX.places.tsv and PREFIX.summary.txt."""
import argparse, collections, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import records as R

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True); ap.add_argument("--tie-split", required=True); ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--window", type=float, default=10); ap.add_argument("--min-molecules", type=int, default=2)
    ap.add_argument("--annot"); ap.add_argument("--gtf", action="append", default=[])
    a = ap.parse_args(argv)
    T, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    prim_spans = collections.defaultdict(list)
    sec = []                                                     # (chrom, strand, start, end, read, record)
    for q, rs in recs.items():
        best = max(r.ms for r in rs)
        for r in rs:
            if r.primary:
                prim_spans[(r.chrom, r.strand)].append((r.start, r.end))
            elif best - r.ms <= a.window:
                sec.append((r.chrom, r.strand, r.start, r.end, q, r))
    for v in prim_spans.values():
        v.sort()
    import bisect
    def has_primary(c, st, s, e):
        v = prim_spans.get((c, st), [])
        i = bisect.bisect_left(v, (e,))
        return any(x0 < e and s < x1 for x0, x1 in v[max(0, i - 3000):i])
    sec.sort(key=lambda x: (x[0], x[1], x[2]))
    places = []
    for c, st, s, e, q, r in sec:
        if places and places[-1]["c"] == c and places[-1]["st"] == st and s < places[-1]["e"]:
            p = places[-1]; p["e"] = max(p["e"], e); p["m"].append((q, r))
        else:
            places.append(dict(c=c, st=st, s=s, e=e, m=[(q, r)]))
    out = []
    for p in places:
        mols = {q for q, _ in p["m"]}
        if len(mols) < a.min_molecules or has_primary(p["c"], p["st"], p["s"], p["e"]):
            continue
        kinds, cls = collections.Counter(), None
        dec = shr = judged = 0
        for q, r in p["m"]:
            t = T.get(q)
            if t is None:
                kinds["no_truth"] += 1; continue
            kinds[t.kind] += 1
            i = R.locus_index(r, t)
            if t.kind == "one":
                judged += 1; dec += i in t.origin
            elif t.kind == "several":
                judged += 1; shr += i in t.origin
        cls = "REAL_DECIDED" if dec else "REAL_SHARED" if shr else "SHADOW" if judged else "UNJUDGED"
        out.append(dict(p, n=len(mols), cls=cls, dec=dec, shr=shr, judged=judged, kinds=kinds))
    gene = {}
    if a.annot:
        from undecidable import load_annot, copy_annotation
        tx, genes = load_annot(a.annot)
        cps = [dict(id=str(i), chrom=p["c"], strand=p["st"], start=p["s"], end=p["e"]) for i, p in enumerate(out)]
        gene = copy_annotation(cps, tx, genes)
    asm = {}
    if a.gtf:
        from undecidable import read_gtf_tx, CopyIndex
        cps = [dict(id=str(i), chrom=p["c"], strand=p["st"], start=p["s"], end=p["e"]) for i, p in enumerate(out)]
        idx = CopyIndex(cps)
        for spec in a.gtf:
            n, g = spec.split("=", 1)
            hit = set()
            for v in read_gtf_tx(g).values():
                if len(v["exons"]) > 1:
                    hit |= {c["id"] for c in idx.hits(v["chrom"], v["strand"], v["exons"][0][0], v["exons"][-1][1])}
            asm[n] = hit
    with open(a.out + ".places.tsv", "w") as fh:
        fh.write("chrom\tstrand\tstart\tend\tmolecules\tclass\tdecided_here\tshared_here\tjudged\tgene\tgroup\t" + "\t".join(f"in_{n}" for n in asm) + "\n")
        for i, p in enumerate(out):
            g = gene.get(str(i), {})
            fh.write(f"{p['c']}\t{p['st']}\t{p['s']}\t{p['e']}\t{p['n']}\t{p['cls']}\t{p['dec']}\t{p['shr']}\t{p['judged']}\t{g.get('gene', '')}\t{g.get('group', '')}\t"
                     + "\t".join(str(int(str(i) in asm[n])) for n in asm) + "\n")
    L = [f"secondary-only places (>= {a.min_molecules} molecules within {a.window} ms, no primary): {len(out)}"]
    c = collections.Counter(p["cls"] for p in out)
    L.append("by class: " + ", ".join(f"{k} {c[k]}" for k in ("REAL_DECIDED", "REAL_SHARED", "SHADOW", "UNJUDGED")))
    for k in ("REAL_DECIDED", "REAL_SHARED", "SHADOW"):
        ps = [i for i, p in enumerate(out) if p["cls"] == k]
        mol = sorted(out[i]["n"] for i in ps)
        grp = collections.Counter(gene.get(str(i), {}).get("group", "-") for i in ps)
        L.append(f"  {k}: places {len(ps)}, molecules median {mol[len(mol)//2] if mol else 0} max {mol[-1] if mol else 0}; "
                 + "groups " + ", ".join(f"{g} {n}" for g, n in grp.most_common(6))
                 + "".join(f"; with a {n} multi-exon transcript {sum(str(i) in asm[n] for i in ps)}" for n in asm))
    big = sorted((p for p in out if p["cls"] == "REAL_DECIDED"), key=lambda p: -p["dec"])[:15]
    L.append("REAL_DECIDED, most molecules decided here: " + "; ".join(
        f"{p['c']}:{p['s']}-{p['e']}{p['st']} {gene.get(str(out.index(p)), {}).get('gene', '')} decided {p['dec']}/{p['n']}" for p in big))
    open(a.out + ".summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()
