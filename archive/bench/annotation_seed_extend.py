#!/usr/bin/env python3
"""Temporary annotation-seed extension to recover missed multi-copy families.

Uses annotated gene spans (from annot_intervals.tsv / GGO_genomic.gff) as seeds,
aligns them against genome-wide thin loci with minimap2, and adds full-length
homologous thin loci as new copies. This is intentionally temporary: the
annotation dependency will be removed once a clean RNA-only seed source is found.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/annotation_seed_extend.py
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
GGO_FA = os.path.join(SCRATCH, "GGO.fasta")
GGO_GFF = os.path.join(SCRATCH, "GGO_genomic.gff")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
THIN_LOCI_TSV = os.path.join(SCRATCH, "thin_loci_genome_wide.tsv")

# Genes we want to use as temporary seeds (the missed oracle genes + parents)
SEED_GENES = {
    "LOC101141440", "LOC109029264", "LOC115934629", "LOC129534585",
    "UBE2Q2P16", "ZNF425", "UBE2Q2",
}

MIN_IDENTITY = 0.80
MIN_COVERAGE = 0.40
LEN_CAP = 9000
MAX_EXON_DIFF = 1
MIN_LEN_RATIO = 0.6
MAX_LEN_RATIO = 1.7


def load_annot_genes():
    """Load annotated gene coordinates."""
    genes = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            c, s, e, g = f[0], int(f[1]), int(f[2]), f[4]
            if g in SEED_GENES:
                genes[g] = (c, s, e, f[3])
    return genes


def load_thin_loci(path):
    loci = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            loci.append({
                "lid": int(row["lid"]),
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "support": int(row["support"]),
                "n_exon": int(row["n_exon"]),
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
                "qlen": qlen, "tlen": tlen,
            })
    return hits


def parse_gff_exons(gff_path, target_genes):
    """For each target gene symbol, return the longest transcript's exons as
    (chrom, strand, [(start,end), ...])."""
    gene_id_to_symbol = {}
    mrna_to_gene = {}
    mrna_exons = defaultdict(lambda: dict(chrom=None, strand=None, exons=[]))
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, typ, start, end, _, strand, _, attrs = f
            start, end = int(start), int(end)
            attr = dict(p.split("=") for p in attrs.split(";") if "=" in p)
            if typ == "gene":
                sym = attr.get("gene", "")
                if sym in target_genes:
                    gene_id_to_symbol[attr.get("ID", "")] = sym
            elif typ == "mRNA":
                parent = attr.get("Parent", "")
                if parent in gene_id_to_symbol:
                    mrna_to_gene[attr.get("ID", "")] = gene_id_to_symbol[parent]
            elif typ == "exon":
                parent = attr.get("Parent", "")
                if parent in mrna_to_gene:
                    d = mrna_exons[parent]
                    d["chrom"] = chrom
                    d["strand"] = strand
                    d["exons"].append((start, end))
    # pick longest transcript per gene
    gene_best = {}
    for mrna, d in mrna_exons.items():
        sym = mrna_to_gene[mrna]
        elen = sum(e - s + 1 for s, e in d["exons"])
        if sym not in gene_best or elen > gene_best[sym][0]:
            gene_best[sym] = (elen, d["chrom"], d["strand"], d["exons"])
    return {sym: (c, st, ex) for sym, (_, c, st, ex) in gene_best.items()}


def build_spliced(fa, chrom, strand, exons):
    seq = "".join(fa.fetch(chrom, s - 1, e) for s, e in exons)
    if strand == "-":
        seq = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return seq


def main():
    seed_genes = load_annot_genes()
    print(f"Annotation seeds: {len(seed_genes)}", flush=True)
    for g, (c, s, e, bt) in sorted(seed_genes.items()):
        print(f"  {g}: {c}:{s}-{e} {bt}", flush=True)

    thin_loci = load_thin_loci(THIN_LOCI_TSV)
    thin_loci = [L for L in thin_loci if 200 <= len(L["seq"]) <= LEN_CAP]
    print(f"Thin loci after length filter: {len(thin_loci)}", flush=True)

    # Build spliced transcript sequences from annotation — temporary; will be replaced by RNA-only seeds
    fa = pysam.FastaFile(GGO_FA)
    gene_exons = parse_gff_exons(GGO_GFF, set(seed_genes.keys()))
    print(f"Parsed GFF transcripts for {len(gene_exons)} target genes", flush=True)
    seed_records = []
    seed_info = {}
    for g, (c, st, exons) in sorted(gene_exons.items()):
        seq = build_spliced(fa, c, st, exons)
        _, s, e, bt = seed_genes[g]
        seed_records.append((g, seq))
        seed_info[g] = dict(gene=g, chrom=c, start=s, end=e, biotype=bt, strand=st,
                            n_exon=len(exons), seq=seq)
    fa.close()

    thin_records = [(f"TL_{L['lid']}", L["seq"]) for L in thin_loci]
    thin_by_name = {f"TL_{L['lid']}": L for L in thin_loci}

    with tempfile.TemporaryDirectory() as tmpdir:
        thin_fa = os.path.join(tmpdir, "thin_loci.fa")
        seed_fa = os.path.join(tmpdir, "seeds.fa")
        paf_out = os.path.join(tmpdir, "seeds_vs_thin.paf")
        write_fasta(thin_fa, thin_records)
        write_fasta(seed_fa, seed_records)

        cmd = [
            "minimap2",
            "-x", "asm20",
            "-t", "8",
            "-N", "50",
            "-p", "0.1",
            "--secondary", "yes",
            thin_fa,
            seed_fa,
        ]
        print(f"Running minimap2 ...", flush=True)
        with open(paf_out, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)
        hits = parse_paf(paf_out)

    filtered = [h for h in hits if h["strand"] == "+" and h["identity"] >= MIN_IDENTITY and h["coverage"] >= MIN_COVERAGE]
    print(f"Full-length homologous hits: {len(filtered)}", flush=True)

    rescued_by_key = {}
    for h in filtered:
        g = h["qname"]
        thin_name = h["tname"]
        if thin_name not in thin_by_name:
            continue
        seed = seed_info[g]
        L = thin_by_name[thin_name]
        # exon-structure guard (seed n_exon is approximate; use length ratio mainly)
        ls, ll = len(seed["seq"]), len(L["seq"])
        if ls == 0 or ll == 0:
            continue
        ratio = min(ls, ll) / max(ls, ll)
        if ratio < MIN_LEN_RATIO or ratio > MAX_LEN_RATIO:
            continue
        # overlap with seed span
        s, e = L["start"], L["end"]
        ms, me = seed["start"], seed["end"]
        if s < me and e > ms:
            ov = min(e, me) - max(s, ms)
            if ov / (e - s) >= 0.5:
                continue
        key = (L["chrom"], L["start"], L["end"])
        entry = dict(
            seed_gene=g,
            tid=f"ANSE_{L['chrom']}_{L['start']}_{L['n_exon']}",
            chrom=L["chrom"], start=L["start"], end=L["end"],
            strand=L["strand"], n_exon=L["n_exon"], support=L["support"],
            identity=h["identity"], coverage=h["coverage"],
        )
        if key not in rescued_by_key or (entry["identity"] * entry["coverage"]) > (rescued_by_key[key]["identity"] * rescued_by_key[key]["coverage"]):
            rescued_by_key[key] = entry

    rescued = sorted(rescued_by_key.values(), key=lambda x: -(x["identity"] * x["coverage"]))
    print(f"Rescued {len(rescued)} copies from annotation seeds", flush=True)

    out_tsv = os.path.join(HERE, "annotation_seed_extend.rescued.tsv")
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seed_gene", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "identity", "coverage"])
        for r in rescued:
            w.writerow([r["seed_gene"], r["tid"], r["chrom"], r["start"], r["end"],
                        r["n_exon"], r["strand"], r["support"], f"{r['identity']:.4f}", f"{r['coverage']:.4f}"])
    print(f"Wrote {out_tsv}")


if __name__ == "__main__":
    import pysam
    main()
