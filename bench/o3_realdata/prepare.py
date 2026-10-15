#!/usr/bin/env python3
"""Inputs of the real-data O3 test (prereg docs/PREREG_o3_realdata_2026-10-15.md).

    prepare.py heldout --out DIR [--min-copies 3] [--cap 60] [--seed 1]   fresh erasure from the KB3781 catalog (A13 excluded)
                       [--min-expressed K] [--exclude-picked P.tsv ...]   K expressed copies suffice; skip families already used
    prepare.py control --out DIR [--n 300] [--n-high 100] [--seed 1]    expressed single-copy genes as one-copy families

heldout writes DIR/:
  picked.tsv          family, masked copy, survivors, masked span (the copy's locus hull, 0-based half-open)
  mask.bed            the masked spans (hard-masked to N in the genome)
  labels.tsv          read, family, role (D = the masked copy, S = a survivor), copy (family:idx)
  reads.fa            the labelled reads in their ORIGINAL orientation (the BAM's sequence reverse-complemented back on
                      reverse-strand records), for re-mapping with -uf
  removed.copies.tsv / .fa   the survivors (catalog rows and sequences): the catalog of the masked arm
  addedback.copies.tsv / .fa all copies of the picked families: the catalog of the added-back arm
  regions.txt         the picked families' catalog regions, merged (copy_assign refuses overlapping regions)
  manifest.json       counts at every step
control writes DIR/: genes.tsv (gene, heterozygosity, bin, reads), control.copies.tsv / .fa, regions.txt, manifest.json.
"""
import argparse
import collections
import csv
import json
import os
import random
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rd_lib as L  # noqa: E402

CATALOG = "/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/gorilla_KB3781"
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
GENOME = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
A13_ALL = "/mnt/linuxdisk/tmp/rna_allele/a14/A14.copies.tsv"
TRUTH = "/mnt/linuxdisk/tmp/rna_allele"
EXCLUDE = 0x4 | 0x100 | 0x800
M_OPS, D_OP, N_OP = (0, 7, 8), 2, 3
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def aligned_blocks(rec):
    """Reference blocks with an aligned base (M / = / X), 0-based half-open."""
    pos, out = rec.reference_start, []
    for op, n in rec.cigartuples:
        if op in M_OPS:
            out.append((pos, pos + n))
            pos += n
        elif op in (D_OP, N_OP):
            pos += n
    return out


def original_seq(rec):
    """The read as sequenced: the BAM stores reverse-strand records reverse-complemented."""
    s = rec.query_sequence
    return s[::-1].translate(COMP) if rec.is_reverse else s


def primaries(bam, spans):
    """Primary alignments overlapping any (chrom, start, end) span -> {read: (chrom, blocks)}."""
    out = {}
    for chrom, s, e in spans:
        for r in bam.fetch(chrom, s, e):
            if r.flag & EXCLUDE or r.cigartuples is None or r.query_name in out:
                continue
            out[r.query_name] = (r.reference_name, aligned_blocks(r))
    return out


def merge_regions(lines):
    """`family<TAB>chrom:start-end` lines -> merged `chrom:start-end` regions (sorted; overlapping or touching join)."""
    ivs = []
    for ln in lines:
        reg = ln.rstrip("\n").split("\t")[1]
        chrom, _, se = reg.rpartition(":")
        s, e = (int(x) for x in se.split("-"))
        ivs.append((chrom, s, e))
    out = []
    for c, s, e in sorted(ivs):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1] = (c, out[-1][1], max(out[-1][2], e))
        else:
            out.append((c, s, e))
    return [f"{c}:{s}-{e}" for c, s, e in out]


def read_catalog(prefix):
    with open(prefix + ".fam.copies.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    copies = []
    for r in rows:
        ex = [tuple(int(x) for x in b.split("-")) for b in r["exons"].split(",")]
        copies.append(dict(family=r["family_id"], idx=int(r["copy_idx"]), chrom=r["chrom"], start=int(r["locus_start"]),
                           end=int(r["locus_end"]), exons=sorted(ex), n=0, row=r))
    return copies


def write_catalog(path_tsv, path_fa, rows, header, fasta, keys):
    with open(path_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(path_fa, "w") as fh:
        for name in fasta.references:
            fam, idx = name.split("|")[:2]
            if (fam, int(idx)) in keys:
                fh.write(f">{name}\n{fasta.fetch(name)}\n")


def cmd_heldout(a):
    os.makedirs(a.out, exist_ok=True)
    copies = read_catalog(a.catalog)
    n_copies = collections.Counter(c["family"] for c in copies)
    cand = [c for c in copies if n_copies[c["family"]] >= a.min_copies]
    with open(a.exclude) as fh:
        excl = [(r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(fh, delimiter="\t")]
    bam = pysam.AlignmentFile(a.bam)
    alns = primaries(bam, [(c["chrom"], c["exons"][0][0], c["exons"][-1][1]) for c in cand])
    rb = L.reads_by_copy([(r, ch, bl) for r, (ch, bl) in sorted(alns.items())],
                         {(c["family"], c["idx"]): (c["chrom"], c["exons"]) for c in cand})
    n_reads = collections.Counter(rb.values())
    for c in cand:
        c["n"] = n_reads[(c["family"], c["idx"])]
    used = set()
    for path in a.exclude_picked:
        with open(path) as fh:
            used |= {r["family"] for r in csv.DictReader(fh, delimiter="\t")}
    picks = L.pick_families(cand, catalog=copies, exclude=excl, min_copies=a.min_copies, cap=a.cap, seed=a.seed,
                            min_expressed=a.min_expressed, exclude_families=used)
    fams = {p["family"] for p in picks}
    masked = {(p["family"], p["masked"]["idx"]) for p in picks}
    labels = L.sample_per_copy({r: k for r, k in rb.items() if k[0] in fams}, n=a.per_copy, seed=a.seed)
    with open(os.path.join(a.out, "picked.tsv"), "w") as fh:
        fh.write("family\tmasked_idx\tsurvivors\tchrom\tstart\tend\n")
        for p in picks:
            m = p["masked"]
            fh.write(f"{p['family']}\t{m['idx']}\t{','.join(str(s['idx']) for s in p['survivors'])}\t{m['chrom']}\t"
                     f"{m['start']}\t{m['end']}\n")
    with open(os.path.join(a.out, "mask.bed"), "w") as fh:
        for p in sorted(picks, key=lambda p: (p["masked"]["chrom"], p["masked"]["start"])):
            fh.write(f"{p['masked']['chrom']}\t{p['masked']['start']}\t{p['masked']['end']}\t{p['family']}\n")
    want = set(labels)
    with open(os.path.join(a.out, "labels.tsv"), "w") as fh:
        fh.write("read\tfamily\trole\tcopy\n")
        for r, (f, i) in sorted(labels.items()):
            fh.write(f"{r}\t{f}\t{'D' if (f, i) in masked else 'S'}\t{f}:{i}\n")
    seqs = {}
    for c in cand:
        if c["family"] in fams:
            for r in bam.fetch(c["chrom"], c["exons"][0][0], c["exons"][-1][1]):
                if r.query_name in want and not r.flag & EXCLUDE and r.query_name not in seqs:
                    seqs[r.query_name] = original_seq(r)
    missing = sorted(want - set(seqs))
    if missing:
        raise SystemExit(f"{len(missing)} labelled reads have no primary sequence, e.g. {missing[:3]}")
    with open(os.path.join(a.out, "reads.fa"), "w") as fh:
        for r in sorted(seqs):
            fh.write(f">{r}\n{seqs[r]}\n")
    header = list(copies[0]["row"].keys())
    fasta = pysam.FastaFile(a.catalog + ".fam.copies.fa")
    surv = {(p["family"], s["idx"]) for p in picks for s in p["survivors"]}
    allk = surv | masked
    write_catalog(os.path.join(a.out, "removed.copies.tsv"), os.path.join(a.out, "removed.copies.fa"),
                  [c["row"] for c in copies if (c["family"], c["idx"]) in surv], header, fasta, surv)
    write_catalog(os.path.join(a.out, "addedback.copies.tsv"), os.path.join(a.out, "addedback.copies.fa"),
                  [c["row"] for c in copies if (c["family"], c["idx"]) in allk], header, fasta, allk)
    with open(a.catalog + ".fam.copies.regions") as fh:
        regs = merge_regions([ln for ln in fh if ln.split("\t")[0] in fams])
    with open(os.path.join(a.out, "regions.txt"), "w") as fh:
        fh.writelines(r + "\n" for r in regs)
    role = collections.Counter("D" if k in masked else "S" for k in labels.values())
    man = dict(catalog=a.catalog, bam=a.bam, exclude=a.exclude, min_copies=a.min_copies,
               candidate_families=len({c["family"] for c in cand}),
               picked=len(picks), masked_copies=len(masked), survivors=len(surv), reads_D=role["D"], reads_S=role["S"],
               regions=len(regs), cap=a.cap, seed=a.seed, per_copy=a.per_copy, min_expressed=a.min_expressed,
               exclude_picked=a.exclude_picked, excluded_families=len(used))
    json.dump(man, open(os.path.join(a.out, "manifest.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(man, sort_keys=True))


def spliced(genome, chrom, exons, strand):
    s = "".join(genome.fetch(chrom, b, e) for b, e in exons).upper()
    return s[::-1].translate(COMP) if strand == "-" else s


def cmd_control(a):
    os.makedirs(a.out, exist_ok=True)
    sets = {r["gene_id"]: r for r in csv.DictReader(open(os.path.join(a.truth, "sets.tsv")), delimiter="\t")}
    lift = {r["gene_id"]: r for r in csv.DictReader(open(os.path.join(a.truth, "lift.tsv")), delimiter="\t")}
    genome = pysam.FastaFile(a.genome)
    contigs = set(genome.references)
    cand = []
    for g in csv.DictReader(open(os.path.join(a.truth, "genes.tsv")), delimiter="\t"):
        s, lf = sets.get(g["gene_id"]), lift.get(g["gene_id"])
        if not s or s["set"] != "S_single" or not lf or g["chrom"] not in contigs:
            continue
        if float(lf["lift_frac"]) < 0.9 or int(g["exonic_bp"]) < 500:
            continue
        d = (int(lf["mismatches"]) + int(lf["indels"])) / max(1, int(lf["lifted_bp"]))
        ex = sorted(tuple(int(x) for x in b.split("-")) for b in g["exons"].split(","))
        cand.append(dict(gene=g["gene_id"], chrom=g["chrom"], strand=g["strand"], exons=ex, het=d))
    bam = pysam.AlignmentFile(a.bam)

    def expressed(pool, n):
        out = []
        for c in pool:
            if len(out) >= n:
                break
            al = primaries(bam, [(c["chrom"], c["exons"][0][0], c["exons"][-1][1])])
            k = sum(1 for ch, bl in al.values() if L._overlap(sorted(bl), c["exons"]) > 0)
            if k >= a.min_reads:
                out.append(dict(c, reads=k))
        return out
    rng = random.Random(a.seed)
    pool = sorted(cand, key=lambda c: c["gene"])
    rng.shuffle(pool)
    random_set = expressed(pool, a.n)
    high = [c for c in pool if c["het"] > L.DELTA and c["gene"] not in {x["gene"] for x in random_set}]
    high_set = expressed(high, a.n_high)
    genes = [dict(g, stratum="random") for g in random_set] + [dict(g, stratum="high_het") for g in high_set]
    header = ["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "exons"]
    with open(os.path.join(a.out, "genes.tsv"), "w") as fh:
        fh.write("family\tgene\tstratum\thet\tbin\treads\n")
        for g in genes:
            fh.write(f"SC:{g['gene']}\t{g['gene']}\t{g['stratum']}\t{g['het']:.6f}\t{L.het_bin(g['het'])}\t{g['reads']}\n")
    with open(os.path.join(a.out, "control.copies.tsv"), "w") as fh, open(os.path.join(a.out, "control.copies.fa"), "w") as fa:
        fh.write("\t".join(header) + "\n")
        for g in genes:
            ex, s, e = g["exons"], g["exons"][0][0], g["exons"][-1][1]
            exs = ",".join(f"{x}-{y}" for x, y in ex)
            fh.write(f"SC:{g['gene']}\t0\tSC:{g['gene']}\t{g['chrom']}\t{s}\t{e}\t{len(ex)}\t{g['strand']}\t{g['reads']}\t{exs}\n")
            fa.write(f">SC:{g['gene']}|0|{g['chrom']}:{s}-{e}|{g['strand']}|nexon={len(ex)}\n"
                     f"{spliced(genome, g['chrom'], ex, g['strand'])}\n")
    regs = merge_regions([f"SC:{g['gene']}\t{g['chrom']}:{max(0, g['exons'][0][0] - 2000)}-{g['exons'][-1][1] + 2000}"
                          for g in genes])
    with open(os.path.join(a.out, "regions.txt"), "w") as fh:
        fh.writelines(r + "\n" for r in regs)
    man = dict(candidates=len(cand), random=len(random_set), high_het=len(high_set), regions=len(regs), seed=a.seed,
               bins=dict(collections.Counter(L.het_bin(g["het"]) for g in random_set)))
    json.dump(man, open(os.path.join(a.out, "manifest.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(man, sort_keys=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("heldout")
    h.add_argument("--out", required=True)
    h.add_argument("--catalog", default=CATALOG)
    h.add_argument("--bam", default=BAM)
    h.add_argument("--exclude", default=A13_ALL)
    h.add_argument("--cap", type=int, default=60)
    h.add_argument("--min-copies", type=int, default=3, help="Amendment 1 of the prereg: 2")
    h.add_argument("--per-copy", type=int, default=500)
    h.add_argument("--min-expressed", type=int, help="K: >= K copies with >= 20 reads (default: every copy)")
    h.add_argument("--exclude-picked", nargs="*", default=[], help="picked.tsv files of earlier builds: their families are skipped")
    h.add_argument("--seed", type=int, default=1)
    c = sub.add_parser("control")
    c.add_argument("--out", required=True)
    c.add_argument("--bam", default=BAM)
    c.add_argument("--genome", default=GENOME)
    c.add_argument("--truth", default=TRUTH)
    c.add_argument("--n", type=int, default=300)
    c.add_argument("--n-high", type=int, default=100)
    c.add_argument("--min-reads", type=int, default=20)
    c.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    {"heldout": cmd_heldout, "control": cmd_control}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
