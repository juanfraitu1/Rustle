#!/usr/bin/env python3
"""Are segmental duplications where recent multi-copy gene families live? (user, 2026-11-19)

  sd_families.py --sd SEDEF.bed --out PREFIX [--species human] [--min-cover 0.5]

Genes: protein-coding CAT genes with an Ensembl ID Compara knows (famtruth's judgeable universe). Each gene's AGE CLASS = the
youngest duplication level among its Compara within-species paralogue pairs (Homo sapiens < Homininae < Hominidae <
Hominoidea < Catarrhini < Simiiformes < Haplorrhini < Primates < ... ); genes with no paralogue pair = `no paralogue`.
SD (SEDEF, both copies of every pair; identity = fracMatch, column 21): a gene is IN AN SD when >= --min-cover of its exonic bp
lie in SD intervals; its SD identity = the highest identity among the SDs covering it.
Pair test: a paralogue pair is DUPLICATED BY AN SD when one SEDEF pair has one copy overlapping each gene's exons.
Writes PREFIX.genes.tsv, PREFIX.summary.txt."""
import argparse, bisect, collections, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import famtruth as F

ORDER = ["Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini", "Primates",
         "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia", "Amniota", "Tetrapoda", "Sarcopterygii",
         "Euteleostomi", "Vertebrata", "Gnathostomata", "Chordata", "Bilateria", "Opisthokonta"]
RANK = {x: i for i, x in enumerate(ORDER)}

def load_sd(path):
    """-> (intervals per chrom [(s, e, identity, pair id, side)], sorted), number of pairs"""
    iv = collections.defaultdict(list)
    n = 0
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 21 or f[0].startswith("#"):
            continue
        try:
            ident = float(f[20])
        except ValueError:
            continue
        iv[f[0]].append((int(f[1]), int(f[2]), ident, n, 0))
        iv[f[3]].append((int(f[4]), int(f[5]), ident, n, 1))
        n += 1
    for v in iv.values():
        v.sort()
    return iv, n

def overlaps(iv, chrom, s, e, maxlen):
    v = iv.get(chrom, [])
    i = bisect.bisect_left(v, (e,))
    j = bisect.bisect_left(v, (s - maxlen,))
    return [x for x in v[j:i] if x[1] > s]

def merged_cover(exons, hits):
    tot = sum(b - a for a, b in exons)
    if not tot:
        return 0.0
    cov = 0
    for a, b in exons:
        segs = sorted((max(a, x[0]), min(b, x[1])) for x in hits if x[0] < b and a < x[1])
        cur = None
        for s, e in segs:
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                if cur:
                    cov += cur[1] - cur[0]
                cur = [s, e]
        if cur:
            cov += cur[1] - cur[0]
    return cov / tot

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sd", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--species", default="human"); ap.add_argument("--min-cover", type=float, default=0.5)
    a = ap.parse_args(argv)
    G = F.load_genes(a.species)
    by_ens = collections.defaultdict(set)
    for g, d in G.items():
        if d["biotype"] == "protein_coding":
            for e in d["ens"]:
                by_ens[e].add(g)
    youngest, pairs, known = {}, collections.defaultdict(str), set()
    for line in F._open(F.PARALOGS):
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("Gene"):
            continue
        e0 = f[0].split(".")[0]; known.add(e0)
        if len(f) < 4 or not f[1] or f[3] not in RANK:
            continue
        e1 = f[1].split(".")[0]
        for x in by_ens.get(e0, ()):
            for y in by_ens.get(e1, ()):
                if x == y:
                    continue
                lv = f[3]
                for g in (x, y):
                    if g not in youngest or RANK[lv] < RANK[youngest[g]]:
                        youngest[g] = lv
                k = (min(x, y), max(x, y))
                if not pairs[k] or RANK[lv] < RANK[pairs[k]]:
                    pairs[k] = lv
    universe = {g for e, gs in by_ens.items() if e in known for g in gs}
    iv, npairs = load_sd(a.sd)
    maxlen = max((x[1] - x[0] for v in iv.values() for x in v), default=0)
    info = {}
    for g in universe:
        d = G[g]
        s, e = min(x for x, _ in d["exons"]), max(y for _, y in d["exons"])
        hits = overlaps(iv, d["chrom"], s, e, maxlen)
        hits = [h for h in hits if any(h[0] < y and x < h[1] for x, y in d["exons"])]
        cov = merged_cover(d["exons"], hits)
        info[g] = dict(cover=cov, ident=max((h[2] for h in hits), default=0.0), sd_pairs={h[3] for h in hits},
                       cls=youngest.get(g, "no paralogue"))
    L = [f"protein-coding genes Compara knows: {len(universe)}; SEDEF pairs {npairs}; in an SD = >= {a.min_cover:.0%} of exonic bp covered"]
    L.append("age class (youngest paralogue level) | genes | in an SD | in an SD >= 0.95 / 0.98 / 0.99 identity")
    groups = [("Homo sapiens",), ("Homininae", "Hominidae"), ("Hominoidea", "Catarrhini"), ("Simiiformes", "Haplorrhini", "Primates"),
              ("Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia"), tuple(ORDER[13:]), ("no paralogue",)]
    for grp in groups:
        gs = [g for g in universe if info[g]["cls"] in grp]
        if not gs:
            continue
        ins = [g for g in gs if info[g]["cover"] >= a.min_cover]
        hi = [sum(info[g]["ident"] >= t for g in ins) for t in (0.95, 0.98, 0.99)]
        L.append(f"  {'/'.join(grp) if len(grp) < 4 else grp[0] + ' .. ' + grp[-1]} | {len(gs)} | {len(ins)} ({len(ins) / len(gs):.1%}) | "
                 + " / ".join(f"{h} ({h / len(gs):.1%})" for h in hi))
    recent = {"Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini"}
    sd_genes = [g for g in universe if info[g]["cover"] >= a.min_cover and info[g]["ident"] >= 0.95]
    L.append(f"inverse: genes in an SD >= 0.95 identity {len(sd_genes)}; with a recent (<= Catarrhini) paralogue "
             f"{sum(info[g]['cls'] in recent for g in sd_genes)} ({sum(info[g]['cls'] in recent for g in sd_genes) / max(1, len(sd_genes)):.1%}); "
             f"no paralogue {sum(info[g]['cls'] == 'no paralogue' for g in sd_genes)}")
    rg = [g for g in universe if info[g]["cls"] in recent]
    L.append(f"recent-family genes outside any SD (< {a.min_cover:.0%} covered): {sum(info[g]['cover'] < a.min_cover for g in rg)} of {len(rg)}")
    L.append("pair test: paralogue pairs joined by ONE SEDEF pair (one copy on each gene), by the pair's level | pairs | joined | median SD identity")
    for grp in groups[:-1]:
        ps = [k for k, lv in pairs.items() if lv in grp and k[0] in info and k[1] in info]
        j = [k for k in ps if info[k[0]]["sd_pairs"] & info[k[1]]["sd_pairs"]]
        ids = sorted(max((h for h in [info[k[0]]["ident"]]), default=0) for k in j)
        L.append(f"  {'/'.join(grp) if len(grp) < 4 else grp[0] + ' .. ' + grp[-1]} | {len(ps)} | {len(j)} ({len(j) / max(1, len(ps)):.1%}) | "
                 f"{ids[len(ids) // 2] if ids else 0:.3f}")
    # mechanism of the joined pairs: does the joining SD carry genomic sequence beyond the exons (introns / flanks: a DNA-based,
    # segmental duplication) or only the exons (a processed retrocopy, RNA-based)? Per joined pair, for the gene with introns:
    # the share of its intronic bp inside the joining SD's copy on that gene; >= 0.5 = genomic, else exonic-only. Single-exon
    # genes: genomic when the SD copy extends >= 1 kb beyond the gene on either side.
    sd_by_id = {}
    for v in iv.values():
        for x in v:
            sd_by_id[(x[3], x[4])] = x
    L.append("mechanism of the joined pairs (genomic = the SD carries introns or >= 1 kb of flank; exonic-only = a retrocopy-like SD)")
    for grp in groups[:-1]:
        ps = [k for k, lv in pairs.items() if lv in grp and k[0] in info and k[1] in info and info[k[0]]["sd_pairs"] & info[k[1]]["sd_pairs"]]
        mech = collections.Counter()
        for x, y in ps:
            shared = info[x]["sd_pairs"] & info[y]["sd_pairs"]
            best = None
            for g in (x, y):
                d = G[g]; ex = sorted(d["exons"])
                ints = [(p_[1], q[0]) for p_, q in zip(ex, ex[1:]) if q[0] > p_[1]]
                segs = [sd_by_id[(pid, side)] for pid in shared for side in (0, 1) if (pid, side) in sd_by_id
                        and sd_by_id[(pid, side)][0] < ex[-1][1] and ex[0][0] < sd_by_id[(pid, side)][1]]
                if not segs:
                    continue
                if ints:
                    v = merged_cover(ints, segs)
                    best = max(best or 0, v)
                else:
                    ext = max(max(ex[0][0] - s_[0], s_[1] - ex[-1][1]) for s_ in segs)
                    best = max(best or 0, 1.0 if ext >= 1000 else 0.0)
            mech["genomic" if (best or 0) >= 0.5 else "exonic-only"] += 1
        tot = sum(mech.values())
        L.append(f"  {'/'.join(grp) if len(grp) < 4 else grp[0] + ' .. ' + grp[-1]} | joined {tot} | genomic {mech['genomic']} "
                 f"({mech['genomic'] / max(1, tot):.1%}) | exonic-only {mech['exonic-only']}")
    out_sd = [g for g in universe if info[g]["cls"] in recent and info[g]["cover"] < a.min_cover]
    L.append("recent-family genes outside SDs, most common name stems: " + ", ".join(
        f"{k} {n}" for k, n in collections.Counter(F.re.sub(r"\d+[A-Z]?\d*P?$", "", g) for g in out_sd).most_common(12)))
    with open(a.out + ".genes.tsv", "w") as fh:
        fh.write("gene\tchrom\tage_class\tsd_cover\tsd_identity\n")
        for g in sorted(universe):
            fh.write(f"{g}\t{G[g]['chrom']}\t{info[g]['cls']}\t{info[g]['cover']:.3f}\t{info[g]['ident']:.4f}\n")
    open(a.out + ".summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()
