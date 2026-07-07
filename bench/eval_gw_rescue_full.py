#!/usr/bin/env python3
"""Evaluate gw_rescue vs gw_xcbase: oracle genes, known-family windows, annotations."""
import csv
import os
import sys
from bisect import bisect_right
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

SCRATCH = "/home/juanfra/winloci_scratch"
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ORACLE = os.path.join("bench", "diploid_cn_oracle.tsv")

KNOWN_WINDOWS = {
    "RABL2": {"seeds": {"RABL2A", "RABL2B"}, "prefix": "RABL2"},
    "APOBEC3": {"seeds": {"APOBEC3C", "APOBEC3D", "APOBEC3F"}, "prefix": "APOBEC3"},
    "MAGEA": {"seeds": {"MAGEA1", "MAGEA4", "MAGEA9", "MAGEA12"}, "prefix": "MAGEA"},
    "ANKRD18": {"seeds": {"ANKRD18A", "ANKRD18B"}, "prefix": "ANKRD18"},
    "RGPD8": {"seeds": {"RGPD8", "RANBP2"}, "prefix": "RGPD"},
    "ZNF92": {"seeds": {"ZNF92"}, "prefix": "ZNF92"},
    "GSTM": {"seeds": {"GSTM1", "GSTM2", "GSTM4", "GSTM5"}, "prefix": "GSTM"},
    "HERC2": {"seeds": {"HERC2"}, "prefix": "HERC2"},
}


def load_annot():
    by_c = defaultdict(list)
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            by_c[f[0]].append((int(f[1]), int(f[2]), f[3], f[4]))
    for c, ivals in by_c.items():
        ivals.sort()
    return by_c


def gene_of_factory(by_c):
    idx = {}
    for c, ivals in by_c.items():
        starts = [x[0] for x in ivals]
        maxlen = max((e - s) for s, e, _, _ in ivals) if ivals else 0
        idx[c] = (starts, ivals, maxlen)

    def gene_of(chrom, qs, qe):
        rec = idx.get(chrom)
        if rec is None:
            return None, None, 0.0
        starts, ivals, maxlen = rec
        hi = bisect_right(starts, qe)
        best_ov = 0
        best = (None, None)
        best_span = 0
        i = hi - 1
        lo_bound = qs - maxlen
        while i >= 0:
            s, e, bt, g = ivals[i]
            if s < lo_bound:
                break
            if e > qs:
                ov = min(e, qe) - max(s, qs)
                if ov > best_ov:
                    best_ov = ov
                    best = (g, bt)
                    best_span = e - s
            i -= 1
        if best[0] is None:
            return None, None, 0.0
        dn_span = qe - qs
        frac_dn = (best_ov / dn_span) if dn_span > 0 else 0.0
        frac_gene = (best_ov / best_span) if best_span > 0 else 0.0
        return best[0], best[1], max(frac_dn, frac_gene)
    return gene_of


def load_catalog(copies_tsv):
    copies = []
    with open(copies_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            copies.append(dict(
                family_id=row["family_id"],
                tid=row["tid"],
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                n_exon=int(row["n_exon"]),
                strand=row["strand"],
                n_reads=int(row["n_reads"]),
            ))
    return copies


def known_family_windows(annot):
    windows = {}
    for name, spec in KNOWN_WINDOWS.items():
        seeds = spec["seeds"]
        prefix = spec["prefix"]
        genes = set()
        for chrom, ivals in annot.items():
            for s, e, bt, g in ivals:
                if g in seeds or (prefix and g.startswith(prefix)):
                    genes.add((chrom, s, e, g))
        if not genes:
            windows[name] = None
            continue
        by_chrom = defaultdict(lambda: [1 << 60, -1, set()])
        for chrom, s, e, g in genes:
            by_chrom[chrom][0] = min(by_chrom[chrom][0], s)
            by_chrom[chrom][1] = max(by_chrom[chrom][1], e)
            by_chrom[chrom][2].add(g)
        windows[name] = [
            dict(chrom=c, start=lo, end=hi, genes=gs)
            for c, (lo, hi, gs) in by_chrom.items()
        ]
    return windows


def copies_in_window(copies, windows, pad=0):
    out = []
    for c in copies:
        for w in windows or []:
            if c["chrom"] == w["chrom"] and c["end"] >= w["start"] - pad and c["start"] <= w["end"] + pad:
                out.append(c)
                break
    return out


def load_oracle_genes():
    genes = []
    with open(ORACLE) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            gene, cls = f[1], f[12]
            if gene and gene != "NA" and cls == "multi_copy":
                genes.append(gene)
    return set(genes)


def covered_oracle_genes(copies, oracle_set, by_c):
    covered = set()
    for c in copies:
        for s, e, bt, g in by_c.get(c["chrom"], []):
            if g in oracle_set and c["end"] >= s and c["start"] <= e:
                covered.add(g)
    return covered


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-copies", default=os.path.join(SCRATCH, "gw_xcbase.copies.tsv"))
    ap.add_argument("--test-copies", default=os.path.join(SCRATCH, "gw_rescue.copies.tsv"))
    ap.add_argument("--rescued", default=os.path.join("bench", "gw_rescue_prototype.rescued.tsv"))
    args = ap.parse_args()

    annot = load_annot()
    gene_of = gene_of_factory(annot)
    oracle_genes = load_oracle_genes()
    windows = known_family_windows(annot)

    base_copies = load_catalog(args.baseline_copies)
    test_copies = load_catalog(args.test_copies)
    annotate(base_copies, gene_of)
    annotate(test_copies, gene_of)

    base_cov = covered_oracle_genes(base_copies, oracle_genes, annot)
    test_cov = covered_oracle_genes(test_copies, oracle_genes, annot)

    print(f"Oracle multi-copy genes: {len(oracle_genes)}")
    print(f"  baseline covered: {len(base_cov)} ({len(base_cov)/len(oracle_genes):.3f})")
    print(f"  test covered:     {len(test_cov)} ({len(test_cov)/len(oracle_genes):.3f})")
    print(f"  newly covered:    {sorted(test_cov - base_cov)}")
    print(f"  still missed:     {sorted(oracle_genes - test_cov)}")

    print("\nKnown-family window copies:")
    for name, wins in windows.items():
        base_n = len(copies_in_window(base_copies, wins, pad=100_000))
        test_n = len(copies_in_window(test_copies, wins, pad=100_000))
        print(f"  {name}: {base_n} -> {test_n} (+{test_n - base_n})")

    # Annotation purity of NEW copies
    base_tids = {c["tid"] for c in base_copies}
    new_copies = [c for c in test_copies if c["tid"] not in base_tids]
    print(f"\nNew copies added by rescue: {len(new_copies)}")
    biotypes = defaultdict(int)
    unannot = 0
    for c in new_copies:
        if c["gene"] is None:
            unannot += 1
        else:
            biotypes[c["biotype"]] += 1
    print(f"  unannotated: {unannot}")
    for bt, n in sorted(biotypes.items(), key=lambda kv: -kv[1]):
        print(f"  {bt}: {n}")


def annotate(copies, gene_of):
    for c in copies:
        g, bt, ov = gene_of(c["chrom"], c["start"], c["end"])
        c["gene"] = g
        c["biotype"] = bt
        c["gene_ov"] = ov


if __name__ == "__main__":
    main()
