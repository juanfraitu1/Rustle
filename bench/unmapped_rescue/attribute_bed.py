#!/usr/bin/env python3
"""Consensus, attribution and rescue accounting of a clustered bed (docs/PREREG_unmapped_rescue_2026-10-08.md section 1 steps 3-4 and Amendment 1).

    attribute_bed.py <bed> [--tag registered]   # needs W/<bed>/{pool.fa,labels.tsv,targets.fa} and W/<bed>/<tag>/clusters.tsv
Runs consensus (miniforge python, pyabpoa), dc-megablast and mmseqs2 translated search (foreground, ~3 min), then writes W/<bed>/<tag>/rescue.json:
the FROZEN rule (cover score, margin 1.10) as the headline, the sweep 1.00 / 1.10 / 1.50 and the translated comparator."""
import argparse
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import attribute as A  # noqa: E402
import score as S  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
PY = "/home/juanfra/miniforge3/bin/python"
BLAST = "/home/juanfra/miniforge3/envs/blast/bin"


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed")
    ap.add_argument("--tag", default="registered")
    a = ap.parse_args()
    d = f"{W}/{a.bed}"
    r = f"{d}/{a.tag}"
    tmp = os.environ.get("TMPDIR", "/tmp")
    if not os.path.exists(f"{r}/cons.fa"):
        sh(f"{PY} {HERE}/consensus.py {d}/pool.fa {r}/clusters.tsv {r}/cons.fa")
    if not os.path.exists(f"{r}/cons.blastn.tsv"):
        if not os.path.exists(f"{d}/targets_db.nsq"):
            sh(f"PATH={BLAST}:$PATH makeblastdb -in {d}/targets.fa -dbtype nucl -out {d}/targets_db > /dev/null")
        sh(f"PATH={BLAST}:$PATH blastn -task dc-megablast -query {r}/cons.fa -db {d}/targets_db -evalue 1e-5 -num_threads 4 -max_target_seqs 5000 "
           f"-outfmt '6 qseqid sseqid bitscore evalue length pident qstart qend' -out {r}/cons.blastn.tsv")
    if not os.path.exists(f"{r}/cons.m8"):
        sh(f"rm -rf {tmp}/mm_tmp_{a.bed}; mmseqs easy-search {r}/cons.fa {d}/targets.fa {r}/cons.m8 {tmp}/mm_tmp_{a.bed} --search-type 2 -e 1e-3 -s 7.5 "
           f"--threads 4 --format-output query,target,bits,evalue -v 1 > /dev/null")
    lab = {x["read"]: x for x in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    dread = {n for n, x in lab.items() if x["role"] == "D"}
    cl = {}
    for x in csv.DictReader(open(f"{r}/clusters.tsv"), delimiter="\t"):
        cl.setdefault("cl" + x["cluster"], []).append(x["read"])
    hsps = list(A.read_blastn(f"{r}/cons.blastn.tsv"))
    capped = A.capped_queries(hsps, 5000)
    assert not capped, f"BLAST target cap reached by {len(capped)} queries"
    sc = A.cover_scores(hsps)
    out = {"clusters": len(cl), "n_d": len(dread)}
    for name, margin in (("registered_1.10", 1.10), ("strict_1.00", 1.00), ("margin_1.50", 1.50)):
        att = A.attribute_cover(sc, margin)
        m = S.rescue_metrics(cl, {c: att[c][0] if c in att else None for c in cl}, truth, dread)
        m["wrong_fraction_of_joined"] = m["wrong_joins"] / max(1, m["rescued_correct"] + m["wrong_joins"])
        m["rescued_fraction_of_d"] = m["rescued_correct"] / len(dread)
        out[name] = m
    tr = A.attribute([(q.split("|")[0], t, b, e) for q, t, b, e in A.read_m8(f"{r}/cons.m8")], 1e-5)
    m = S.rescue_metrics(cl, {c: tr[c][0] if c in tr else None for c in cl}, truth, dread)
    m["wrong_fraction_of_joined"] = m["wrong_joins"] / max(1, m["rescued_correct"] + m["wrong_joins"])
    m["rescued_fraction_of_d"] = m["rescued_correct"] / len(dread)
    out["translated_mmseqs_strict"] = m
    att = A.attribute_cover(sc, 1.10)
    out["abstaining_clusters"] = sorted(((len(rs), S.majority(rs, truth) or "") for c, rs in cl.items() if not att.get(c, (None,))[0]), reverse=True)[:15]
    json.dump(out, open(f"{r}/rescue.json", "w"), indent=1)
    for k in ("registered_1.10", "strict_1.00", "margin_1.50", "translated_mmseqs_strict"):
        m = out[k]
        print(f"{k:26s} attributed {m['clusters_attributed']:3d} abstained {m['clusters_abstained']:3d} cluster_acc {m['cluster_accuracy']} rescued {m['rescued_correct']} "
              f"({m['rescued_fraction_of_d']:.3f} of D) wrong {m['wrong_joins']} ({m['wrong_fraction_of_joined']:.3f} of joined) copies {m['copies_reached']}/{m['copies_with_unmapped_reads']}")


if __name__ == "__main__":
    main()
