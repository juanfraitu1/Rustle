#!/usr/bin/env python3
"""Soto 2025 family table (S1C) audits for the 2026-09-30 meeting (second-machine tasks B1-B3 of
docs/archive/2026-09/HANDOFF_SECOND_MACHINE_2026-09-30.md). Pure re-reads of S1C plus the frozen exon edges; no new data.

    python3 bench/soto_m2/soto_m2_audit.py biotype
    python3 bench/soto_m2/soto_m2_audit.py fragments --exons sd98_gene_exons.tsv
    python3 bench/soto_m2/soto_m2_audit.py narrower --geneset elig.tsv --full-geneset full.tsv \
        --famcn-ours famcn_ours_allwssd.tsv [--pairs attributed_pairs.tsv]

biotype    B1: per-family biotype composition from S1C's own `Biotype` column, checked against S1C's own
           per-family `No. Protein Coding` column.
fragments  B2: each member's merged exonic bp vs the family's largest member; families holding a member < 20% next
           to one >= 80% (the largest itself excluded, so the >= 80% member is a second full-length copy).
narrower   B3: our sequence-only clusters (exon map-back edges, no copy-number gate) that hold >= 2 Soto families:
           the Soto families, their sizes and copy numbers (S1C and ours), the sequence edges that cross Soto's
           family boundaries inside the cluster, and (with --pairs) the 83 attributed missing paralog pairs.
"""
import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
import soto_replication as sr  # noqa: E402

S1C = os.path.join(HERE, "..", "soto", "soto_famCN_S1C.tsv")
EDGES = os.path.join(HERE, "..", "soto", "shared_exons_5154_exon_mapback.tsv")


def read_s1c(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def families(rows):
    """Family ID -> list of S1C rows (a gene in k families appears in k rows); 'Unassigned*' ids skipped."""
    fam = defaultdict(list)
    for r in rows:
        f = r["Family ID"]
        if f and not f.startswith("Unassigned"):
            fam[f].append(r)
    return fam


def is_pseudo(bt):
    return "pseudogene" in bt


def pct(x, n):
    return f"{x}/{n} ({100.0 * x / n:.1f}%)" if n else "0/0"


def cmd_biotype(a):
    fam = families(read_s1c(a.truth))
    members = sum(len(v) for v in fam.values())
    pseudo = sum(is_pseudo(r["Biotype"]) for v in fam.values() for r in v)
    all_pseudo = [f for f, v in fam.items() if all(is_pseudo(r["Biotype"]) for r in v)]
    no_pc = [f for f, v in fam.items() if not any(r["Biotype"] == "protein_coding" for r in v)]
    majority = [f for f, v in fam.items() if 2 * sum(is_pseudo(r["Biotype"]) for r in v) > len(v)]
    # P0: our tabulation vs S1C's own per-family "No. Protein Coding" column
    mism = 0
    for f, v in fam.items():
        own = {r["No. Protein Coding"] for r in v}
        mine = sum(r["Biotype"] == "protein_coding" for r in v)
        if len(own) != 1 or float(next(iter(own))) != mine:
            mism += 1
    bt = defaultdict(int)
    for v in fam.values():
        for r in v:
            bt[r["Biotype"]] += 1
    print(f"## B1: biotype composition of Soto's S1C families\n")
    print(f"- members (rows with a Family ID): {members}; families: {len(fam)}")
    print(f"- pseudogene-biotype members: {pct(pseudo, members)}")
    print(f"- families entirely pseudogene-biotype: {pct(len(all_pseudo), len(fam))}")
    print(f"- families with no protein-coding member: {pct(len(no_pc), len(fam))}")
    print(f"- pseudogene-majority families: {pct(len(majority), len(fam))}")
    print(f"- P0 check, protein-coding count vs S1C's own `No. Protein Coding` column: {mism}/{len(fam)} mismatches")
    print("\n| biotype | members |\n|---|---|")
    for k, n in sorted(bt.items(), key=lambda t: -t[1]):
        print(f"| {k} | {n} |")
    sizes = defaultdict(int)
    for f in all_pseudo:
        sizes[len(fam[f])] += 1
    print("\nAll-pseudogene families by size: " + ", ".join(f"{k} members: {n}" for k, n in sorted(sizes.items())))


def exonic_bp(path):
    iv = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            iv[r["gene_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    out = {}
    for g, xs in iv.items():
        tot, cur = 0, None
        for c, s, e in sorted(xs):
            if cur and c == cur[0] and s <= cur[2]:
                cur = (c, cur[1], max(cur[2], e))
            else:
                if cur:
                    tot += cur[2] - cur[1]
                cur = (c, s, e)
        tot += cur[2] - cur[1]
        out[g] = tot
    return out


def cmd_fragments(a):
    fam = families(read_s1c(a.truth))
    bp = exonic_bp(a.exons)
    rows, n_cmp = [], 0
    for f, v in fam.items():
        mem = [(r["Gene Name"], r["Gene ID"], r["Biotype"], bp[r["Gene ID"]]) for r in v if r["Gene ID"] in bp]
        if len(mem) < 2:
            continue
        n_cmp += 1
        mem.sort(key=lambda t: -t[3])
        big = mem[0][3]
        frac = [m[3] / big for m in mem]
        small = [m for m, x in zip(mem, frac) if x < a.small]
        full_other = [m for m, x in zip(mem[1:], frac[1:]) if x >= a.full]
        full_any = [m for m, x in zip(mem, frac) if x >= a.full]
        rows.append((f, mem, small, full_other, full_any, big))
    both_excl = [t for t in rows if t[2] and t[3]]
    both_incl = [t for t in rows if t[2] and t[4]]
    print(f"## B2: fragment-sized members bundled with full-length members\n")
    print(f"- families with >= 2 members carrying exon coordinates: {n_cmp} (of {len(fam)})")
    print(f"- families with a member < {a.small:.0%} of the largest member's exonic bp (the largest counts as "
          f"full-length): {pct(len(both_incl), n_cmp)}")
    print(f"- ... and ALSO a second member >= {a.full:.0%} (two full-length copies next to a fragment): "
          f"{pct(len(both_excl), n_cmp)}")
    worst = sorted(both_incl, key=lambda t: -(t[5] / max(1, min(m[3] for m in t[2]))))
    print("\nWorst size ratios (largest member / smallest member, exonic bp):\n")
    print("| family | members | largest (bp, biotype) | smallest (bp, biotype) | ratio |\n|---|---|---|---|---|")
    for f, mem, small, _fo, _fa, big in worst[:a.top]:
        s = min(small, key=lambda m: m[3])
        print(f"| {f} | {len(mem)} | {mem[0][0]} ({big}, {mem[0][2]}) | {s[0]} ({s[3]}, {s[2]}) | "
              f"{big / max(1, s[3]):.0f}x |")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("family_id\tgene_name\tgene_id\tbiotype\texonic_bp\tfrac_of_largest\n")
            for f, mem, *_ in rows:
                for m in mem:
                    fh.write(f"{f}\t{m[0]}\t{m[1]}\t{m[2]}\t{m[3]}\t{m[3] / mem[0][3]:.4f}\n")


def cmd_narrower(a):
    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    edges = sr.read_edges(a.shared)
    clean, ambiguous = sr.load_truth(a.truth)
    s1c = read_s1c(a.truth)
    name = {r["Gene ID"]: r["Gene Name"] for r in s1c}
    cn_soto = sr.load_famcn(a.truth, "Median famCN")
    cn_ours = sr.load_famcn(a.famcn_ours, "famCN_sotoiv") if a.famcn_ours else {}
    members = defaultdict(set)
    for g, f in clean.items():
        if f and not f.startswith("Unassigned"):
            members[f].add(g)
    soto = {f: m for f, m in members.items() if len(m) >= 2}
    g2f = {g: f for f, m in soto.items() for g in m}
    cover, kept, leaf_of = sr.pair_families(edges, genes, full, {}, gate=False)
    pred = sr.collapse_cover(cover, leaf_of, genes)
    universe = (full | genes) - ambiguous
    clusters = defaultdict(set)
    for g in universe:
        if pred.get(g):
            clusters[pred[g]].add(g)
    multi = []
    for c, v in clusters.items():
        touched = sorted({g2f[g] for g in v if g in g2f}, key=lambda f: -len(soto[f]))
        if len(touched) >= 2:
            multi.append((c, v, touched))
    multi.sort(key=lambda t: (-len(t[2]), -len(t[1])))
    med = lambda xs: statistics.median(xs) if xs else float("nan")
    # sequence edges that cross a Soto family boundary inside one of our clusters (both ends in Soto families)
    cross = [(x, y) for x, y in kept if x in g2f and y in g2f and g2f[x] != g2f[y] and pred.get(x) and
             pred.get(x) == pred.get(y)]
    cross_pairs = {frozenset(p) for p in cross}
    gaps_s = [abs(cn_soto[x] - cn_soto[y]) for x, y in cross_pairs_list(cross_pairs) if x in cn_soto and y in cn_soto]
    gaps_o = [abs(cn_ours[x] - cn_ours[y]) for x, y in cross_pairs_list(cross_pairs) if x in cn_ours and y in cn_ours]
    in_multi = sum(len(t[2]) for t in multi)
    print("## B3: where Soto is finer than our sequence families\n")
    print(f"- our sequence-only clusters (>= 2 genes): {sum(1 for v in clusters.values() if len(v) >= 2)}; "
          f"holding >= 2 Soto families: {len(multi)} (holding {in_multi} of Soto's {len(soto)} multi-gene families)")
    print(f"- distinct >= 98% exon-sharing gene pairs that join two different Soto families inside one of our "
          f"clusters: {len(cross_pairs)}; copy-number gap across them, median |dCN| S1C {med(gaps_s):.1f} "
          f"(n={len(gaps_s)}), ours {med(gaps_o):.1f} (n={len(gaps_o)})")
    grown = [(f, len(m), len(clusters[pred[next(iter(m))]])) for f, m in soto.items()
             if all(pred.get(g) for g in m) and len({pred[g] for g in m}) == 1
             and len(clusters[pred[next(iter(m))]]) > len(m)]
    print(f"- Soto families strictly smaller than the sequence cluster that contains them: {pct(len(grown), len(soto))}"
          f"; genes added: median {med([c - s for _f, s, c in grown]):.0f}, total {sum(c - s for _f, s, c in grown)}")
    # Soto families with ONE exclusive member (the rest are genes Soto also puts in other families): is that lone
    # gene >= 98% exon-sharing with an exclusive member of another Soto family?
    all_members = defaultdict(set)
    for r in s1c:
        if r["Family ID"] and not r["Family ID"].startswith("Unassigned"):
            all_members[r["Family ID"]].add(r["Gene ID"])
    adj = sr.adjacency(edges)
    lone = [f for f, m in members.items() if len(m) == 1]
    linked = []
    for f in sorted(lone, key=lambda x: int(x.split("_")[1])):
        g = next(iter(members[f]))
        other = sorted({clean[h] for h in adj.get(g, ()) if h in clean and clean[h] != f
                        and not clean[h].startswith("Unassigned")})
        if other:
            linked.append((f, g, other))
    print(f"- Soto families with ONE exclusive member (all others shared with other Soto families): {len(lone)}; "
          f"that member shares >= 98% identical exons with an exclusive member of another Soto family: "
          f"{pct(len(linked), len(lone))}")
    for f, g, other in linked:
        if name.get(g, "").startswith(("NPIP", "GOLGA", "USP17", "TBC1D3", "PCMTD", "NIPA", "ANKRD20")):
            print(f"    {f}: {name.get(g, g)} (S1C famCN {cn_soto.get(g, float('nan')):.1f}; "
                  f"{len(all_members[f])} members) -> also linked to {', '.join(other[:6])}")
    if cn_ours:
        # how large a copy-number gap is, next to the disagreement between two measurements of the same gene
        both = [g for g in cn_soto if g in cn_ours]
        print(f"\nCopy-number re-measurement (ours: 268 SGDP samples, Soto's interval) vs S1C, {len(both)} genes:\n")
        print("| S1C famCN | genes | median abs difference | median relative difference | share differing by >= 2.6 |"
              "\n|---|---|---|---|---|")
        for lo, hi in ((0, 10), (10, 20), (20, 35), (35, 60), (60, None)):
            d = [abs(cn_soto[g] - cn_ours[g]) for g in both if cn_soto[g] >= lo and (hi is None or cn_soto[g] < hi)]
            r = [abs(cn_soto[g] - cn_ours[g]) / cn_soto[g] for g in both
                 if cn_soto[g] >= lo and (hi is None or cn_soto[g] < hi)]
            print(f"| {lo}-{hi if hi else ''} | {len(d)} | {med(d):.2f} | {med(r):.3f} | "
                  f"{sum(x >= 2.6 for x in d) / len(d):.2f} |")
        show = ("NPIPB3", "NPIPB4", "NPIPB5", "NPIPB12", "NPIPB13", "NPIPA1", "NPIPA7", "GOLGA8B")
        print("\nFlagship genes, famCN S1C / ours: " + "; ".join(
            f"{name[g]} {cn_soto[g]:.1f} / {cn_ours[g]:.1f}" for g in sorted(both, key=lambda g: name.get(g, g))
            if name.get(g) in show))
    print("\n| our cluster (genes) | Soto families (members; median famCN S1C / ours) | cross-family edges | "
          "example genes |\n|---|---|---|---|")
    for c, v, touched in multi:
        parts = []
        for f in touched:
            m = soto[f]
            cs = med([cn_soto[g] for g in m if g in cn_soto])
            co = med([cn_ours[g] for g in m if g in cn_ours])
            parts.append(f"{f} ({len(m)}; {cs:.1f} / {co:.1f})")
        ncross = sum(1 for p in cross_pairs if next(iter(p)) in v)
        ex = ", ".join(sorted({name.get(g, g) for g in v})[:5])
        print(f"| {len(v)} | {'; '.join(parts)} | {ncross} | {ex} |")
    if a.pairs:
        with open(a.pairs) as fh:
            prs = list(csv.DictReader(fh, delimiter="\t"))
        cause = defaultdict(list)
        for r in prs:
            cause[r["cause"]].append(float(r["identity"]))
        print(f"\n### The {len(prs)} attributed missing paralog pairs (`attributed_pairs.tsv`, 09-28)\n")
        print("| cause | pairs | median identity | chromosomes |\n|---|---|---|---|")
        for k, ids in sorted(cause.items()):
            chroms = sorted({r["chrom"] for r in prs if r["cause"] == k})
            print(f"| {k} | {len(ids)} | {med(ids):.3f} | {', '.join(chroms)} |")
        ids = [float(r["identity"]) for r in prs]
        print(f"\nAll: median identity {med(ids):.3f}; >= 0.90: {sum(x >= .90 for x in ids)}; "
              f">= 0.95: {sum(x >= .95 for x in ids)}; >= 0.98: {sum(x >= .98 for x in ids)}")
        top = sorted(prs, key=lambda r: -float(r["identity"]))[:a.top]
        print("\nMost identical pairs Soto keeps apart:\n\n| pair | identity | cause |\n|---|---|---|")
        for r in top:
            print(f"| {r['gene_a']} / {r['gene_b']} | {float(r['identity']):.3f} | {r['detail']} |")


def cross_pairs_list(ps):
    return [tuple(sorted(p)) for p in ps]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=S1C)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("biotype").set_defaults(fn=cmd_biotype)
    p = sub.add_parser("fragments")
    p.add_argument("--exons", required=True)
    p.add_argument("--small", type=float, default=0.20)
    p.add_argument("--full", type=float, default=0.80)
    p.add_argument("--top", type=int, default=10)
    p.add_argument("--out")
    p.set_defaults(fn=cmd_fragments)
    p = sub.add_parser("narrower")
    p.add_argument("--shared", default=EDGES)
    p.add_argument("--geneset", required=True)
    p.add_argument("--full-geneset", required=True)
    p.add_argument("--famcn-ours")
    p.add_argument("--pairs")
    p.add_argument("--top", type=int, default=10)
    p.set_defaults(fn=cmd_narrower)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
