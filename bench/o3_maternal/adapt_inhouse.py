#!/usr/bin/env python3
"""In-house arm -> the common candidate tables.   adapt_inhouse.py <workdir>   (workdir = W/inhouse, after the o3_candidates batches)

Reads cand_g*.candidates.tsv (columns family candidate n_clusters n_reads flagged union_len nearest_locus d n_net n_used) and
cand_g*.contigs.fa; also writes cands_all.tsv. A candidate is a new-copy candidate iff flagged == 1 (distance to the nearest locus > delta); n_transcripts = n_clusters."""
import csv
import glob
import os
import sys


def read_fa(path):
    seq, cur = {}, None
    for ln in open(path):
        if ln[0] == ">":
            cur = ln[1:].strip()
            seq[cur] = []
        elif cur:
            seq[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seq.items()}


def main(w):
    rows, seqs = [], {}
    for f in sorted(glob.glob(f"{w}/cand_g*.candidates.tsv")):
        rows += list(csv.DictReader(open(f), delimiter="\t"))
        fa = f.replace(".candidates.tsv", ".contigs.fa")
        if os.path.exists(fa):
            seqs.update(read_fa(fa))
    new = [r for r in rows if r["flagged"] == "1"]
    with open(f"{w}/cands_all.tsv", "w") as o:        # every stage candidate, flagged or not (used by compare.py: "found in the reference")
        o.write("candidate\tfamily\tn_clusters\tflagged\tnearest_locus\td\n")
        for r in rows:
            o.write(f"{r['candidate']}\t{r['family']}\t{r['n_clusters']}\t{r['flagged']}\t{r['nearest_locus']}\t{r['d']}\n")
    with open(f"{w}/cands.tsv", "w") as o:
        o.write("candidate\tfamily\tn_transcripts\tcontigs\n")
        for r in new:
            o.write(f"{r['candidate']}\t{r['family']}\t{r['n_clusters']}\t{r['candidate']}\n")
    with open(f"{w}/cands.fa", "w") as o:
        for r in new:
            if r["candidate"] in seqs:
                o.write(f">{r['candidate']}\n{seqs[r['candidate']]}\n")
    print(f"stage candidates {len(rows)} in {len({r['family'] for r in rows})} families; new-copy (flagged==1) {len(new)}; "
          f">= 2 clusters {sum(1 for r in new if int(r['n_clusters']) >= 2)}")


if __name__ == "__main__":
    main(sys.argv[1])
