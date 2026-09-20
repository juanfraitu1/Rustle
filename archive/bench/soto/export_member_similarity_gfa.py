#!/usr/bin/env python3
"""Export compact member-similarity graphs as Bandage-readable GFA1.

Unlike the base-level DNA variation graphs in ``dna_graphs/``, each segment in
these graphs is one annotated gene copy.  A link is retained only when a
pairwise whole-locus minimap2 alignment covers at least MIN_COVERAGE of both
copies, has at least MIN_IDENTITY identity, and reaches MIN_SCORE.  The GFA
link tags retain identity (ID), reciprocal coverage (CV), and SC = ID * CV.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import subprocess
from pathlib import Path

MIN_IDENTITY = 0.80
MIN_COVERAGE = 0.30
MIN_SCORE = 0.80
COLORS = {"NPIPA": "#4285f4", "NPIPB": "#ea4335", "NPIPP": "#fbbc04", "FAM72": "#34a853"}
OTHER = "#9aa0a6"


def family_members(bed: Path, family_id: str):
    ans = []
    for row in csv.reader(bed.open(), delimiter="\t"):
        if len(row) >= 4 and row[3].endswith("|" + family_id):
            gene, _ = row[3].rsplit("|", 1)
            ans.append((gene, row[0], int(row[1]), int(row[2])))
    if len(ans) < 2:
        raise ValueError(f"{family_id}: expected at least two members")
    return ans


def class_of(gene: str) -> str:
    if gene.startswith("NPIPA"): return "NPIPA"
    if gene.startswith("NPIPB"): return "NPIPB"
    if gene.startswith("NPIPP"): return "NPIPP"
    if gene.startswith("FAM72"): return "FAM72"
    return "other"


def extract_sequences(members, genome: Path, out_fa: Path):
    with out_fa.open("w") as fh:
        for gene, chrom, start, end in members:
            region = f"{chrom}:{start + 1}-{end}"  # BED is zero-based, faidx is one-based
            seq = subprocess.check_output(["samtools", "faidx", str(genome), region], text=True)
            fh.write(f">{gene}\n")
            fh.write("".join(x for x in seq.splitlines() if not x.startswith(">")) + "\n")


def align(a_name, a_seq, b_name, b_seq, tmp: Path):
    ref, qry = tmp / "ref.fa", tmp / "qry.fa"
    ref.write_text(f">{a_name}\n{a_seq}\n")
    qry.write_text(f">{b_name}\n{b_seq}\n")
    paf = subprocess.check_output(
        ["minimap2", "-x", "asm20", "--secondary=no", "-c", str(ref), str(qry)],
        text=True, stderr=subprocess.DEVNULL)
    best = None
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        nmatch, alen = int(f[9]), int(f[10])
        qcov = (int(f[3]) - int(f[2])) / int(f[1])
        tcov = (int(f[8]) - int(f[7])) / int(f[6])
        cand = (nmatch / alen, min(qcov, tcov), nmatch, alen)
        if best is None or cand[2] > best[2]:
            best = cand
    return best or (0.0, 0.0, 0, 0)


def write_graph(family_id, members, seqs, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_rows = []
    tmp = out_dir / ".similarity_tmp"
    tmp.mkdir(exist_ok=True)
    for (ga, _, _, _), (gb, _, _, _) in itertools.combinations(members, 2):
        identity, coverage, matches, aln_len = align(ga, seqs[ga], gb, seqs[gb], tmp)
        pair_rows.append({"member_1": ga, "member_2": gb, "identity": identity,
                          "reciprocal_coverage": coverage, "score": identity * coverage,
                          "matches": matches, "alignment_length": aln_len,
                          "retained": identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE and identity * coverage >= MIN_SCORE})
    # Do not leave temporary FASTAs in the deliverable directory.
    (tmp / "ref.fa").unlink(missing_ok=True); (tmp / "qry.fa").unlink(missing_ok=True); tmp.rmdir()

    stem = out_dir / family_id
    with (stem.with_suffix(".pairwise_similarity.tsv")).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_rows[0].keys(), delimiter="\t")
        w.writeheader(); w.writerows(pair_rows)
    with (stem.with_suffix(".gfa")).open("w") as fh:
        fh.write("H\tVN:Z:1.0\tTS:Z:member_similarity\n")
        for gene, chrom, start, end in members:
            fh.write(f"S\t{gene}\tN\tLN:i:{len(seqs[gene])}\tCL:Z:{class_of(gene)}\tLO:Z:{chrom}:{start}-{end}\n")
        for row in pair_rows:
            if row["retained"]:
                fh.write("L\t{member_1}\t+\t{member_2}\t+\t0M\tID:f:{identity:.4f}\tCV:f:{reciprocal_coverage:.4f}\tSC:f:{score:.4f}\n".format(**row))
        for gene, *_ in members:
            fh.write(f"P\t{gene}\t{gene}+\t*\n")
    with (stem.with_suffix(".colours.csv")).open("w") as fh:
        fh.write("Node,Colour\n")
        for gene, *_ in members:
            fh.write(f"{gene},{COLORS.get(class_of(gene), OTHER)}\n")
    with (stem.with_suffix(".legend.tsv")).open("w") as fh:
        fh.write("member\tclass\tlocus\tlength_bp\tcolour\n")
        for gene, chrom, start, end in members:
            fh.write(f"{gene}\t{class_of(gene)}\t{chrom}:{start}-{end}\t{len(seqs[gene])}\t{COLORS.get(class_of(gene), OTHER)}\n")
    retained = sum(r["retained"] for r in pair_rows)
    return len(members), len(pair_rows), retained


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bed", default="bench/soto/80_fams.chr.bed")
    p.add_argument("--genome", default="../Reference/chm13v2.0.fa")
    p.add_argument("--out-dir", default="bench/soto/member_similarity_graphs")
    args = p.parse_args()
    for fid in ("ID_154_NPIP", "ID_354_FAM72"):
        family_id = fid.split("_", 1)[0] + "_" + fid.split("_", 2)[1]
        members = family_members(Path(args.bed), family_id)
        fa = Path(args.out_dir) / f".{fid}.fa"
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)
        extract_sequences(members, Path(args.genome), fa)
        seqs, name, buf = {}, None, []
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:], []
            else: buf.append(line)
        seqs[name] = "".join(buf); fa.unlink()
        n, total, kept = write_graph(fid, members, seqs, Path(args.out_dir))
        print(f"{fid}: {n} members; retained {kept}/{total} pairwise links")


if __name__ == "__main__":
    main()
