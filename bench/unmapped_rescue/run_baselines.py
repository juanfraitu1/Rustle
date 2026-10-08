#!/usr/bin/env python3
"""B0 / B1 against the pooled method on the SAME seeded sample of unmapped deleted-copy reads (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, Amendment 1).

    run_baselines.py <bed> [--n 500]      # writes W/<bed>/baselines.json"""
import argparse
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import attribute as A  # noqa: E402
import baselines as B  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
BLAST = "/home/juanfra/miniforge3/envs/blast/bin"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed")
    ap.add_argument("--n", type=int, default=500)
    a = ap.parse_args()
    d = f"{W}/{a.bed}"
    lab = {x["read"]: x for x in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    dreads = [n for n, x in lab.items() if x["role"] == "D"]
    samp = B.sample(dreads, a.n, seed=1)
    seqs = SD.read_fa(f"{d}/pool.fa")
    SD.write_fa(f"{d}/sample.fa", seqs, samp)
    if not os.path.exists(f"{d}/sample.paf"):
        subprocess.run(f"minimap2 -x map-ont -c -N 20 -t 4 {d}/targets.fa {d}/sample.fa > {d}/sample.paf 2> {d}/sample.paf.log", shell=True, check=True)
    if not os.path.exists(f"{d}/sample.blastn.tsv"):
        subprocess.run(f"PATH={BLAST}:$PATH blastn -task dc-megablast -query {d}/sample.fa -db {d}/targets_db -evalue 1e-5 -num_threads 4 -max_target_seqs 500 "
                       f"-outfmt '6 qseqid sseqid bitscore evalue length pident qstart qend' -out {d}/sample.blastn.tsv", shell=True, check=True)
    fam = {n: lab[n]["family"] for n in samp}
    b0 = B.single_read_nucleotide(open(f"{d}/sample.paf"))
    att1 = A.attribute_cover(A.cover_scores(list(A.read_blastn(f"{d}/sample.blastn.tsv"))), 1.10)
    b1 = {q: v[0] for q, v in att1.items()}
    # the pooled method on the same reads: the read's cluster, attributed by the frozen rule
    pooled_att = A.attribute_cover(A.cover_scores(list(A.read_blastn(f"{d}/registered/cons.blastn.tsv"))), 1.10)
    cl_of = {x["read"]: "cl" + x["cluster"] for x in csv.DictReader(open(f"{d}/registered/clusters.tsv"), delimiter="\t")}
    pooled = {n: pooled_att[cl_of[n]][0] for n in samp if n in cl_of and cl_of[n] in pooled_att}

    def tally(att):
        ok = sum(1 for n in samp if att.get(n) == fam[n])
        wrong = sum(1 for n in samp if att.get(n) not in (None, fam[n]))
        return dict(rescued_correct=ok, wrong=wrong, abstained_or_unplaced=len(samp) - ok - wrong, fraction_correct=ok / len(samp))
    out = dict(n=len(samp), B0_single_read_nucleotide=tally(b0), B1_single_read_dc_megablast=tally(b1), pooled_frozen=tally(pooled))
    json.dump(out, open(f"{d}/baselines.json", "w"), indent=1)
    for k, v in out.items():
        print(k, v)


if __name__ == "__main__":
    main()
