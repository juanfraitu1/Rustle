#!/usr/bin/env python3
"""Annotation-free targeted rescue for missed oracle genes.

For each missed oracle gene, take all thin loci overlapping its annotation span
as seeds. Align those seeds against the full genome-wide thin-locus database with
very permissive minimap2 parameters. Report any thin-locus hits (same or other
chromosome) that could represent additional copies.
"""
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

SCRATCH = "/home/juanfra/winloci_scratch"
HERE = os.path.dirname(__file__)
THIN_LOCI_TSV = os.path.join(SCRATCH, "thin_loci_genome_wide.tsv")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ORACLE = os.path.join(HERE, "diploid_cn_oracle.tsv")

ORACLE_GENES = ["LOC109029264", "UBE2Q2P16", "ZNF425"]


def load_annot():
    gene_to_coords = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            c, s, e, g = f[0], int(f[1]), int(f[2]), f[4]
            if g not in gene_to_coords or (e - s) > (gene_to_coords[g][2] - gene_to_coords[g][1]):
                gene_to_coords[g] = (c, s, e)
    return gene_to_coords


def load_thin_loci():
    loci = []
    with open(THIN_LOCI_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            introns = []
            if "introns" in row and row["introns"]:
                for part in row["introns"].split(","):
                    d, a = part.split("-")
                    introns.append((int(d), int(a)))
            loci.append({
                "lid": int(row["lid"]),
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "support": int(row["support"]),
                "n_exon": int(row["n_exon"]),
                "introns": introns,
                "seq": row["seq"],
            })
    return loci


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def parse_paf(paf_path):
    hits = []
    with open(paf_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qlen, qs, qe, strand, tname, tlen, ts, te, nmatch, aln_len, mapq = f[:12]
            qlen, qs, qe = int(qlen), int(qs), int(qe)
            tlen, ts, te = int(tlen), int(ts), int(te)
            nmatch, aln_len = int(nmatch), int(aln_len)
            identity = nmatch / aln_len if aln_len > 0 else 0.0
            qcov = (qe - qs) / qlen if qlen > 0 else 0.0
            tcov = (te - ts) / tlen if tlen > 0 else 0.0
            coverage = min(qcov, tcov)
            hits.append({
                "qname": qname, "tname": tname, "strand": strand,
                "identity": identity, "coverage": coverage,
                "aln_len": aln_len, "nmatch": nmatch,
                "qlen": qlen, "tlen": tlen,
            })
    return hits


def run_minimap2(query_fa, db_fa, out_paf, threads=8, asm_preset="asm20"):
    cmd = [
        "minimap2",
        "-x", asm_preset,
        "-t", str(threads),
        "-N", "50",
        "-p", "0.1",
        "--secondary", "yes",
        db_fa,
        query_fa,
    ]
    with open(out_paf, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)


def targeted_rescue(gene, coords, loci, min_identity=0.50, min_coverage=0.30,
                    max_exon_diff=2, min_len_ratio=0.5, max_len_ratio=2.0):
    chrom, gs, ge = coords
    gene_loci = [L for L in loci if L["chrom"] == chrom and L["start"] <= ge and L["end"] >= gs]
    print(f"\n=== Targeted rescue: {gene} {chrom}:{gs}-{ge} ===")
    print(f"Seed loci overlapping gene: {len(gene_loci)}")
    if not gene_loci:
        return []

    # Use all gene-overlapping loci as seeds, but avoid self-vs-self by excluding same locus id
    locus_by_name = {f"TL_{L['lid']}": L for L in loci}
    with tempfile.TemporaryDirectory() as tmpdir:
        seed_fa = os.path.join(tmpdir, "seeds.fa")
        db_fa = os.path.join(tmpdir, "thin_loci.fa")
        paf_out = os.path.join(tmpdir, "seeds_vs_thin.paf")
        write_fasta(seed_fa, [(f"SEED_{L['lid']}", L["seq"]) for L in gene_loci])
        write_fasta(db_fa, [(f"TL_{L['lid']}", L["seq"]) for L in loci])
        run_minimap2(seed_fa, db_fa, paf_out)
        hits = parse_paf(paf_out)

    print(f"Raw minimap2 alignments: {len(hits)}")

    # Filter
    good = []
    for h in hits:
        seed_lid = int(h["qname"].split("_")[1])
        target_lid = int(h["tname"].split("_")[1])
        if seed_lid == target_lid:
            continue
        if h["identity"] < min_identity or h["coverage"] < min_coverage:
            continue
        seed = next(L for L in gene_loci if L["lid"] == seed_lid)
        target = locus_by_name[h["tname"]]
        # exon compatibility
        if abs(target["n_exon"] - seed["n_exon"]) > max_exon_diff:
            continue
        ls, ll = len(seed["seq"]), len(target["seq"])
        ratio = min(ls, ll) / max(ls, ll)
        if ratio < min_len_ratio or ratio > max_len_ratio:
            continue
        good.append({
            "seed_lid": seed_lid,
            "target": target,
            "identity": h["identity"],
            "coverage": h["coverage"],
            "strand": h["strand"],
        })

    # Best hit per target locus
    best_by_lid = {}
    for g in good:
        lid = g["target"]["lid"]
        if lid not in best_by_lid or g["identity"] * g["coverage"] > best_by_lid[lid]["identity"] * best_by_lid[lid]["coverage"]:
            best_by_lid[lid] = g

    results = list(best_by_lid.values())
    results.sort(key=lambda x: -(x["identity"] * x["coverage"]))
    print(f"Rescued candidates after filtering: {len(results)}")
    for r in results[:20]:
        t = r["target"]
        print(f"  target lid={t['lid']} {t['chrom']}:{t['start']}-{t['end']} "
              f"strand={t['strand']} n_exon={t['n_exon']} support={t['support']} "
              f"identity={r['identity']:.3f} coverage={r['coverage']:.3f} strand={r['strand']}")
    if len(results) > 20:
        print(f"  ... and {len(results)-20} more")
    return results


def main():
    gene_coords = load_annot()
    loci = load_thin_loci()
    print(f"Loaded {len(loci)} thin loci")

    all_results = {}
    for gene in ORACLE_GENES:
        coords = gene_coords.get(gene)
        if coords is None:
            print(f"{gene}: no annotation coords")
            continue
        all_results[gene] = targeted_rescue(gene, coords, loci)

    # Summary of cross-chrom hits
    print("\n=== SUMMARY ===")
    for gene, results in all_results.items():
        chrom = gene_coords[gene][0]
        cross = [r for r in results if r["target"]["chrom"] != chrom]
        same = [r for r in results if r["target"]["chrom"] == chrom]
        print(f"{gene}: {len(results)} total candidates, {len(cross)} cross-chrom, {len(same)} same-chrom")


if __name__ == "__main__":
    main()
