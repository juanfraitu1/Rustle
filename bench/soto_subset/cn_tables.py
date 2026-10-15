#!/usr/bin/env python3
"""Copy number per Soto gene from three sources, keyed by CAT v4 gene id (an identifier, never a gene name).

    assembly    results_s1e.tsv (bench/soto/parcn_assembly.py, frozen 2026-09-29): one row per S1E region (gene ∩ SD98,
                CHM13 v2.0, BED half-open); `fam_HG002` = the HG002 diploid famCN analogue, empty when the region is
                unresolved (< 100 usable 30-mers). A gene's value = the median over its resolved rows.
    read depth  famcn_ours_allwssd.tsv: `gene_id`, `famCN_sotoiv` (268 SGDP genomes, Soto's interval).
    S1C         soto_famCN_S1C.tsv: `Gene ID`, `Median famCN`; a gene in several families repeats its one value.

check_coordinates() tests the id join against coordinates: every S1E region must overlap its gene's lifted body in
soto_truth.v2.tsv (GFF 1-based closed, CHM13 v2.0).

    cn_tables.py [--s1e P] [--read-depth P] [--s1c P] [--truth-v2 P]     -> JSON report; exit 1 on a disagreement
"""
import argparse
import collections
import csv
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
S1E = "/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv"
READ_DEPTH = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/famcn_ours_allwssd.tsv"
S1C = os.path.join(REPO, "bench", "soto", "soto_famCN_S1C.tsv")
TRUTH_V2 = "/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv"
MISSING = {"", "nan", "NaN", "NA", "N/A", "None"}


def num(v):
    """A table cell -> float, or None for a missing or non-numeric cell."""
    if v is None or v.strip() in MISSING:
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def assembly_cn(path, column="fam_HG002"):
    """-> (values {gene_id: median over resolved S1E rows}, unresolved {gene_id with S1E rows, none resolved})."""
    vals, seen = collections.defaultdict(list), set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "S1E":
                continue
            seen.add(r["gene_id"])
            v = num(r[column])
            if v is not None:
                vals[r["gene_id"]].append(v)
    out = {g: float(statistics.median(v)) for g, v in vals.items()}
    return out, seen - set(out)


def column_cn(path, id_col, value_col):
    """Any gene-keyed table -> {gene_id: value}. A gene repeated with a different value is an error."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = num(r[value_col])
            if v is None:
                continue
            g = r[id_col]
            if g in out and out[g] != v:
                raise ValueError(f"{path}: {g} has two values in {value_col}: {out[g]} and {v}")
            out[g] = v
    return out


def read_depth_cn(path, column="famCN_sotoiv"):
    return column_cn(path, "gene_id", column)


def s1c_cn(path):
    return column_cn(path, "Gene ID", "Median famCN")


def s1c_genes(path):
    """Every distinct `Gene ID` of S1C, sorted."""
    with open(path) as fh:
        return sorted({r["Gene ID"] for r in csv.DictReader(fh, delimiter="\t")})


def account(genes, vals, unresolved=frozenset()):
    """How each gene of `genes` stands in one table -> {valued, unresolved, absent}."""
    c = collections.Counter("valued" if g in vals else "unresolved" if g in unresolved else "absent" for g in genes)
    return {k: c.get(k, 0) for k in ("valued", "unresolved", "absent")}


def check_coordinates(s1e_path, truth_v2_path):
    """-> (n_ok, disagreements [(rid, gene_id, region, body)], n_gene_absent). A region (BED half-open) agrees when
    it shares >= 1 base with its gene's body (GFF closed) on the same chromosome."""
    body = {}
    with open(truth_v2_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            body[r["gene_id"]] = (r["chrom"], int(r["start"]), int(r["end"]))
    ok, bad, absent = 0, [], 0
    with open(s1e_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "S1E":
                continue
            b = body.get(r["gene_id"])
            if b is None:
                absent += 1
                continue
            s, e = int(r["start"]), int(r["end"])
            if r["chrom"] == b[0] and s < b[2] and e >= b[1]:
                ok += 1
            else:
                bad.append((r["rid"], r["gene_id"], f"{r['chrom']}:{s}-{e}", f"{b[0]}:{b[1]}-{b[2]}"))
    return ok, bad, absent


def write_table(path, vals, column="famCN"):
    """{gene_id: value} -> the 2-column TSV that soto_replication.py's load_famcn reads (gene_id, famCN)."""
    with open(path, "w") as fh:
        fh.write(f"gene_id\t{column}\n")
        fh.writelines(f"{g}\t{vals[g]!r}\n" for g in sorted(vals))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s1e", default=S1E)
    ap.add_argument("--read-depth", default=READ_DEPTH)
    ap.add_argument("--s1c", default=S1C)
    ap.add_argument("--truth-v2", default=TRUTH_V2)
    a = ap.parse_args(argv)
    genes = s1c_genes(a.s1c)
    asm, unresolved = assembly_cn(a.s1e)
    ok, bad, absent = check_coordinates(a.s1e, a.truth_v2)
    rep = dict(genes=len(genes), assembly=account(genes, asm, unresolved),
               read_depth=account(genes, read_depth_cn(a.read_depth)), s1c=account(genes, s1c_cn(a.s1c)),
               coordinate_check=dict(ok=ok, disagree=len(bad), gene_absent=absent, first=bad[:5]))
    print(json.dumps(rep, indent=1))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
