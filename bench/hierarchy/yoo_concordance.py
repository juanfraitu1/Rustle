#!/usr/bin/env python3
"""Y block of docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md, Amendment 3: concordance of our loci with the Yoo et al. 2025 gorilla tables.

Match rule (registered): two loci on the same chromosome match iff their overlap is at least 50% of the shorter locus.
All coordinates are 1-based closed intervals on the RefSeq (NC_) names of the gorilla assembly; the Yoo tables use chrN / CM names and are mapped with GGO_chr_map.tsv.

  y1  our LRPAP1 loci against Table VIII.38 (block 'LRPAP1, DOK7, HGFAC') and the AMRP rows of Table VIII.35
  y2  the families of a BASE copies (registered) or clusters table against the 27 unit loci of Table VIII.38 (chr1 MAPKBP1 / JMJD7-PLA2G4B / SPTBN5 + ancestral chr16)

Tests: python3 -B -m unittest bench/hierarchy/test_yoo_concordance.py
"""
import argparse
import csv
import sys


def overlap(a, b):
    """bp overlap of two closed intervals (start, end)."""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def matches(a, b):
    """a, b = (chrom, start, end). Same chromosome and overlap >= 50% of the shorter locus."""
    if a[0] != b[0]:
        return False
    n = overlap((a[1], a[2]), (b[1], b[2]))
    shorter = min(a[2] - a[1] + 1, b[2] - b[1] + 1)
    return n > 0 and 2 * n >= shorter


def max_one_to_one(pairs):
    """Size of a maximum matching of the bipartite graph given as (left, right) pairs (Kuhn's augmenting paths)."""
    adj = {}
    for l, r in pairs:
        adj.setdefault(l, []).append(r)
    owner = {}

    def augment(l, seen):
        for r in adj[l]:
            if r in seen:
                continue
            seen.add(r)
            if r not in owner or augment(owner[r], seen):
                owner[r] = l
                return True
        return False

    return sum(1 for l in adj if augment(l, set()))


def read_chr_map(path):
    """GGO_chr_map.tsv (chr, genbank, refseq, length) -> {'chr1': 'NC_...'}."""
    out = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0] == "chr":
                continue
            out["chr" + row[0]] = row[2]
    return out


def read_genbank_map(path):
    """GGO_chr_map.tsv -> {'CM055446.2': 'NC_073224.2'} (the Yoo table 35 uses the CM names)."""
    out = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0] == "chr":
                continue
            out[row[1]] = row[2]
    return out


def read_unit_loci(path):
    """yoo_GGO_unit_loci.tsv (chr, refseq, start, end, gene, flag) -> [(gene, refseq, start, end)]."""
    out = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0] == "chr":
                continue
            out.append((row[4], row[1], int(row[2]), int(row[3])))
    return out


def read_clusters(path):
    """BASE .fam.clusters.tsv (cluster_id size density frac_in corroborated chrom start end) -> [(family, chrom, start, end)]."""
    out = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0] == "cluster_id":
                continue
            out.append((row[0], row[5], int(row[6]), int(row[7])))
    return out


def read_copies(path, use_locus=False):
    """BASE .fam.copies.tsv (family_id copy_idx tid chrom start end ... locus_start locus_end; start 0-based) -> [(family, chrom, start, end)] 1-based closed.
    The copy rows (transcript span) are the registered object of Amendment 3; use_locus=True gives the locus extent instead."""
    out = []
    with open(path) as f:
        rd = csv.reader(f, delimiter="\t")
        header = next(rd)
        col = {c: i for i, c in enumerate(header)}
        a, b = ("locus_start", "locus_end") if use_locus else ("start", "end")
        for row in rd:
            if row:
                out.append((row[col["family_id"]], row[col["chrom"]], int(row[col[a]]) + 1, int(row[col[b]])))
    return out


def per_gene(yoo, clusters):
    """yoo = [(gene, chrom, start, end)], clusters = [(family, chrom, start, end)] -> {gene: summary dict}."""
    hits = {}  # (yoo index) -> [member index]
    for i, (_g, c, s, e) in enumerate(yoo):
        hits[i] = [j for j, (_f, c2, s2, e2) in enumerate(clusters) if matches((c, s, e), (c2, s2, e2))]
    in_any = {j for js in hits.values() for j in js}
    by_family = {}
    for j, (f, _c, _s, _e) in enumerate(clusters):
        by_family.setdefault(f, []).append(j)
    out = {}
    for gene in dict.fromkeys(g for g, *_ in yoo):
        idx = [i for i, y in enumerate(yoo) if y[0] == gene]
        members = {j for i in idx for j in hits[i]}
        fams = sorted({clusters[j][0] for j in members})
        touched = [j for f in fams for j in by_family[f]]
        out[gene] = {
            "n_yoo": len(idx),
            "yoo_matched": sum(1 for i in idx if hits[i]),
            "one_to_one": max_one_to_one([(i, j) for i in idx for j in hits[i]]),
            "families": fams,
            "members_matching": len(members),
            "touched_members_without_any_yoo_locus": sum(1 for j in touched if j not in in_any),
        }
    return out


def shared_families(per):
    """{(gene_a, gene_b): [families touched by both]} for the gene pairs that share at least one family."""
    genes = sorted(per)
    out = {}
    for a in range(len(genes)):
        for b in range(a + 1, len(genes)):
            common = sorted(set(per[genes[a]]["families"]) & set(per[genes[b]]["families"]))
            if common:
                out[(genes[a], genes[b])] = common
    return out


def read_table38_lrpap1(path):
    """Table VIII.38, the block 'LRPAP1, DOK7, HGFAC in gorGor' (0-based columns 18 to 21): [(name, chr, start, end)] for the LRPAP1 records."""
    out = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) > 21 and row[18].startswith("chr") and row[21].startswith("LRPAP1"):
                out.append((row[21], row[18], int(row[19]), int(row[20])))
    return out


def read_table35_amrp(path, genbank_map, chr_names):
    """Table VIII.35 rows whose gene is AMRP (the UniProt mnemonic of LRPAP1): [(name, chr, start, end)] with the chr of the Yoo tables."""
    inv = {v: k for k, v in chr_names.items()}
    out = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) >= 7 and row[5] == "AMRP":
                out.append((f"AMRP_{row[6]}aa", inv[genbank_map[row[1]]], int(row[2]), int(row[3])))
    return out


def y1(ours, yoo_by_table):
    """ours = [(name, class, chr, start, end)]; yoo_by_table = {table: [(name, chr, start, end)]}. Rows: our locus x matching Yoo record."""
    rows = []
    for name, cls, c, s, e in ours:
        for table, recs in yoo_by_table.items():
            for yname, yc, ys, ye in recs:
                if yc != c:
                    continue
                n = overlap((s, e), (ys, ye))
                if n == 0:
                    continue
                rows.append({
                    "our_locus": name, "our_class": cls, "chr": c, "our_start": s, "our_end": e,
                    "yoo_locus": yname, "yoo_table": table, "yoo_start": ys, "yoo_end": ye,
                    "overlap_bp": n, "frac_of_ours": n / (e - s + 1), "frac_of_yoo": n / (ye - ys + 1),
                    "match": matches((c, s, e), (yc, ys, ye)),
                })
    return rows


def _write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("y1")
    p1.add_argument("--ours", required=True, help="TSV: our_locus our_class chr our_start our_end (further columns ignored)")
    p1.add_argument("--table38", required=True)
    p1.add_argument("--table35", required=True, help="yoo_GGO_table35_novel_paralogs.tsv")
    p1.add_argument("--chr-map", required=True)
    p1.add_argument("--out", required=True)
    p2 = sub.add_parser("y2")
    p2.add_argument("--unit-loci", required=True)
    p2.add_argument("--clusters", help="BASE .fam.clusters.tsv (locus extents, 1-based)")
    p2.add_argument("--copies", help="BASE .fam.copies.tsv (copy rows = transcript spans; the registered object)")
    p2.add_argument("--use-locus", action="store_true", help="with --copies: use locus_start/locus_end instead of the transcript span")
    p2.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args(argv)

    if a.cmd == "y1":
        names = read_chr_map(a.chr_map)
        ours = []
        with open(a.ours) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                ours.append((row["our_locus"], row["our_class"], row["chr"], int(row["our_start"]), int(row["our_end"])))
        t38 = read_table38_lrpap1(a.table38)
        t35 = read_table35_amrp(a.table35, read_genbank_map(a.chr_map), names)
        rows = y1(ours, {"VIII.38": t38, "VIII.35": t35})
        _write(a.out, ["our_locus", "our_class", "chr", "our_start", "our_end", "yoo_locus", "yoo_table", "yoo_start", "yoo_end",
                       "overlap_bp", "frac_of_ours", "frac_of_yoo", "match_rule_min50"],
               [[r["our_locus"], r["our_class"], r["chr"], r["our_start"], r["our_end"], r["yoo_locus"], r["yoo_table"], r["yoo_start"],
                 r["yoo_end"], r["overlap_bp"], f"{r['frac_of_ours']:.3f}", f"{r['frac_of_yoo']:.3f}", "yes" if r["match"] else "no"] for r in rows])
        m38 = {r["yoo_locus"] for r in rows if r["yoo_table"] == "VIII.38" and r["match"]}
        hit = {r["our_locus"] for r in rows if r["match"]}
        print(f"Table VIII.38 loci with a match: {len(m38)}/{len(t38)}; our loci with a match in VIII.38 or VIII.35: {len(hit)}/{len(ours)}")
        return 0

    yoo = read_unit_loci(a.unit_loci)
    if bool(a.clusters) == bool(a.copies):
        ap.error("y2 needs exactly one of --clusters and --copies")
    clusters = read_clusters(a.clusters) if a.clusters else read_copies(a.copies, a.use_locus)
    per = per_gene(yoo, clusters)
    _write(a.out + ".per_gene.tsv",
           ["gene", "n_yoo", "yoo_matched", "one_to_one", "families", "members_matching", "touched_members_without_any_yoo_locus"],
           [[g, r["n_yoo"], r["yoo_matched"], r["one_to_one"], ",".join(r["families"]), r["members_matching"],
             r["touched_members_without_any_yoo_locus"]] for g, r in per.items()])
    rows = []
    for g, c, s, e in yoo:
        js = [(f, cc, ss, ee) for (f, cc, ss, ee) in clusters if matches((c, s, e), (cc, ss, ee))]
        rows.append([g, c, s, e, ";".join(f"{f}:{ss}-{ee}" for f, _cc, ss, ee in js)])
    _write(a.out + ".per_locus.tsv", ["gene", "chrom", "start", "end", "matching_members"], rows)
    for g, r in per.items():
        print(g, {k: r[k] for k in ("n_yoo", "yoo_matched", "one_to_one", "families", "members_matching", "touched_members_without_any_yoo_locus")})
    print("shared families:", shared_families(per))
    return 0


if __name__ == "__main__":
    sys.exit(main())
