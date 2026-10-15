#!/usr/bin/env python3
"""Template-guided copy assembly for one base GTF (spec docs/superpowers/specs/2026-10-14-template-copies-design.md).

    run_tg.py templates --base BASE.gtf --genome FA --out DIR [--batch-size 5000]
    run_tg.py project   --out DIR --index MMI_OR_FA [--batches all|0,3,...] [--threads 4]
    run_tg.py assemble  --base BASE.gtf --bam BAM --out DIR --label NAME

The stages are separate calls, so each stays under the 10-minute limit for one foreground job:
  templates  DIR/templates.fa + .tsv (one representative per locus) and DIR/batches/bNNNN.fa;
  project    minimap2 -x splice -c -N 50 -p 0.1 of each named batch -> DIR/batches/bNNNN.paf + .done
             (a batch with a .done marker is skipped; -c gives the CIGAR the projection needs; -p 0.1 is the guided
             finder's setting, plan ruling R2);
  assemble   every batch must be .done; hits -> candidates -> own reads -> loci; writes DIR/tg.gtf (the added loci),
             DIR/merged.gtf (base + added) and DIR/report.json (counts at every step).
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build as B  # noqa: E402
import candidates as C  # noqa: E402
import evidence as E  # noqa: E402
import templates as T  # noqa: E402

MINIMAP2 = ["minimap2", "-x", "splice", "-c", "-N", "50", "-p", "0.1"]


def stage_templates(base, genome, out, batch_size):
    tx = T.read_gtf(base)
    reps = T.representatives(tx)
    os.makedirs(os.path.join(out, "batches"), exist_ok=True)
    with pysam.FastaFile(genome) as fa:
        T.write_templates(reps, fa, os.path.join(out, "templates.fa"), os.path.join(out, "templates.tsv"))
    n, b, fh = 0, -1, None
    with open(os.path.join(out, "templates.fa")) as src:
        for line in src:
            if line.startswith(">"):
                if n % batch_size == 0:
                    if fh:
                        fh.close()
                    b += 1
                    fh = open(os.path.join(out, "batches", f"b{b:04d}.fa"), "w")
                n += 1
            fh.write(line)
    if fh:
        fh.close()
    rep = dict(transcripts=len(tx), loci=len(reps), templates=n, batches=b + 1, batch_size=batch_size)
    print(json.dumps(rep))
    return rep


def batch_files(out):
    return sorted(glob.glob(os.path.join(out, "batches", "b*.fa")))


def stage_project(out, index, which, threads):
    fas = batch_files(out)
    pick = range(len(fas)) if which == "all" else [int(x) for x in which.split(",")]
    for i in pick:
        fa = fas[i]
        paf, done = fa[:-3] + ".paf", fa[:-3] + ".done"
        if os.path.exists(done):
            continue
        with open(paf + ".tmp", "w") as fh:
            subprocess.run(MINIMAP2 + ["-t", str(threads), index, fa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.replace(paf + ".tmp", paf)
        open(done, "w").close()
        print(f"projected {os.path.basename(fa)}")


def gtf_lines(label, i, kept, others, gene_of, locus):
    gid = f"TG_{label}_{i:06d}"
    attrs = (f'template "{kept.template}"; template_gene "{gene_of.get(kept.template, "")}"; '
             f'identity "{kept.identity:.4f}"; reads "{locus["reads"]}"; also_templates "{len(others)}"; '
             f'source "template_copies";')
    out = []
    for k, ex in enumerate(locus["transcripts"]):
        tid = f"{gid}.{k + 1}"
        head = f"{kept.contig}\ttemplate_copies\t"
        out.append(f'{head}transcript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{kept.strand}\t.\tgene_id "{gid}"; '
                   f'transcript_id "{tid}"; {attrs}\n')
        out += [f'{head}exon\t{s + 1}\t{e}\t.\t{kept.strand}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n'
                for s, e in ex]
    return out


def stage_assemble(base, bam_path, out, label):
    fas = batch_files(out)
    missing = [f for f in fas if not os.path.exists(f[:-3] + ".done")]
    if not fas or missing:
        raise SystemExit(f"assemble needs every batch projected; missing: {[os.path.basename(m) for m in missing]}")
    tx = T.read_gtf(base)
    base_rows = T.loci(tx)
    gene_of = {t.tid: t.gene for t in T.read_templates(os.path.join(out, "templates.tsv")).values()}
    hs = [h for f in fas for h in C.hits(f[:-3] + ".paf")]
    cand = C.uncovered(hs, base_rows)
    groups = C.merge_candidates(cand)
    rep = dict(hits=len(hs), candidates=len(cand), merged_candidates=len(groups), contig_absent=0, few_reads=0,
               no_locus=0, with_loci=0, loci=0, transcripts=0)
    lines = []
    with pysam.AlignmentFile(bam_path) as bam:
        contigs = set(bam.references)
        for kept, others in groups:
            if kept.contig not in contigs:
                rep["contig_absent"] += 1
                continue
            reads = E.own_reads(bam, kept.contig, kept.exons)
            if len(reads) < 2:
                rep["few_reads"] += 1
                continue
            loci = B.build(reads, kept.exons)
            if not loci:
                rep["no_locus"] += 1
                continue
            rep["with_loci"] += 1
            for locus in loci:
                lines += gtf_lines(label, rep["loci"], kept, others, gene_of, locus)
                rep["loci"] += 1
                rep["transcripts"] += len(locus["transcripts"])
    if groups and rep["contig_absent"] > len(groups) / 2:
        raise SystemExit(f"{rep['contig_absent']} of {len(groups)} candidates lie on contigs absent from {bam_path}: "
                         "the index and the BAM name contigs differently")
    with open(os.path.join(out, "tg.gtf"), "w") as fh:
        fh.writelines(lines)
    with open(os.path.join(out, "merged.gtf"), "w") as fh, open(base) as src:
        for line in src:
            fh.write(line)
        fh.writelines(lines)
    with open(os.path.join(out, "report.json"), "w") as fh:
        json.dump(dict(base=base, bam=bam_path, label=label, **rep), fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps(rep))
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("templates")
    p.add_argument("--base", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-size", type=int, default=5000)
    p = sub.add_parser("project")
    p.add_argument("--out", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--batches", default="all")
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("assemble")
    p.add_argument("--base", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--label", required=True)
    a = ap.parse_args(argv)
    if a.stage == "templates":
        stage_templates(a.base, a.genome, a.out, a.batch_size)
    elif a.stage == "project":
        stage_project(a.out, a.index, a.batches, a.threads)
    else:
        stage_assemble(a.base, a.bam, a.out, a.label)
    return 0


if __name__ == "__main__":
    sys.exit(main())
