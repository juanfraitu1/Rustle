#!/usr/bin/env python3
"""Hybrid net on a bed (Amendment 6).   run_hybrid.py <bed>   resumable: exit 75 = run again (one blastn chunk per call). Writes W/<bed>/hybrid/result.json"""
import csv
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import attribute as A  # noqa: E402
import hybrid as H  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
BLAST = "/home/juanfra/miniforge3/envs/blast/bin"
CHUNK = 2000


def main():
    bed = sys.argv[1]
    d, out = f"{W}/{bed}", f"{W}/{bed}/hybrid"
    os.makedirs(out, exist_ok=True)
    lab = {x["read"]: x for x in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    dread = {n for n, x in lab.items() if x["role"] == "D"}
    cl = {}
    for x in csv.DictReader(open(f"{d}/registered/clusters.tsv"), delimiter="\t"):
        cl.setdefault("cl" + x["cluster"], []).append(x["read"])
    seqs = SD.read_fa(f"{d}/pool.fa")
    sc = A.cover_scores(list(A.read_blastn(f"{d}/registered/cons.blastn.tsv")))
    att = A.attribute_cover(sc, 1.10)
    cl_att = {c: att[c][0] if c in att else None for c in cl}
    reads = sorted(seqs)
    resid = H.residual(cl, cl_att, reads)
    chunks = [resid[i:i + CHUNK] for i in range(0, len(resid), CHUNK)]
    secs = 0.0
    for i, ch in enumerate(chunks):
        tsv = f"{out}/chunk{i}.tsv"
        if not os.path.exists(tsv + ".done"):
            SD.write_fa(f"{out}/chunk{i}.fa", seqs, ch)
            t0 = time.time()
            subprocess.run(f"PATH={BLAST}:$PATH blastn -task dc-megablast -query {out}/chunk{i}.fa -db {d}/targets_db -evalue 1e-5 -num_threads 4 "
                           f"-max_target_seqs 5000 -outfmt '6 qseqid sseqid bitscore evalue length pident qstart qend' -out {tsv}", shell=True, check=True)
            open(tsv + ".done", "w").write(str(time.time() - t0))
            print(f"chunk {i} of {len(chunks)} done in {time.time() - t0:.0f} s; run again")
            if i + 1 < len(chunks):
                sys.exit(75)
        secs += float(open(tsv + ".done").read())
    hs = []
    for i in range(len(chunks)):
        hs += list(A.read_blastn(f"{out}/chunk{i}.tsv"))
    assert not A.capped_queries(hs, 5000)
    rd = {q: v[0] for q, v in A.attribute_cover(A.cover_scores(hs), 1.10).items()}
    pooled = H.combine(cl, cl_att, {}, reads)
    hyb = H.combine(cl, cl_att, rd, reads)
    res = dict(bed=bed, pool_reads=len(reads), residual_reads=len(resid), residual_fraction=len(resid) / len(reads), extra_blast_seconds=secs,
               pooled=H.read_metrics(pooled, truth, dread), hybrid=H.read_metrics(hyb, truth, dread), n_d=len(dread))
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    main()
