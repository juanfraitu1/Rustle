#!/usr/bin/env python3
"""Rule T (Amendment 13) on a real bed, descriptive: decisions, agreement with the genome-based finding, identity x coverage after the trim.
    apply_trim_real.py <bedA|bedH>      (run with /home/juanfra/miniforge3/bin/python or python3; needs polish/both.paf and polish/score.json from run_polish.py)"""
import collections
import csv
import json
import sys

import attribute as A
import run_polish as RP
import trim5g as T

W = "/mnt/linuxdisk/tmp/o3_rescue"


def main(bed):
    d = f"{W}/{bed}/registered"
    cons = RP.read_cons(f"{d}/cons.fa")
    hs = collections.defaultdict(list)
    allh = list(A.read_blastn(f"{d}/cons.blastn.tsv"))
    for q, t, a, b in allh:
        hs[q].append((q, t, a, b))
    att = A.attribute_cover(A.cover_scores(allh), 1.10)
    rec = RP.best_records(open(f"{d}/polish/both.paf").read().splitlines(), "orig.")
    sc = json.load(open(f"{d}/polish/score.json"))["clusters"]
    why = collections.Counter()
    agree = collections.Counter()
    ok_before = ok_after = n = 0
    for name, seq in cons.items():
        k = name.split("|")[0]
        if k not in sc or not sc[k]["orig"]["on"]:
            continue
        n += 1
        h = rec[k]
        qs, qe, qlen = h["qstart"], h["qend"], h["qlen"]
        genome_clip = qs if (0 < qs <= 3 and set(seq[:qs]) == {"G"}) else 0
        fam = att.get(k, (None,))[0]
        if fam is None:
            new, tl, reason = seq, 0, "abstained"
        else:
            new, tl, reason = T.decide(seq, T.best_copy_prefix(hs.get(k, []), fam))
        why[reason] += 1
        agree[("genome G clip" if genome_clip else "no genome G clip", "T trims" if tl else "T does not trim")] += 1
        if genome_clip and not tl:
            agree[("missed", reason)] += 1
        ident = h["ident"]
        before = ident * (qe - qs) / qlen
        aligned_after = (qe - qs) - max(0, tl - qs)          # trimmed bases that the aligner had covered are no longer counted
        after = ident * aligned_after / (qlen - tl)
        ok_before += before >= 0.999
        ok_after += after >= 0.999
    print(f"{bed}: {n} consensus sequences on the erased copy; T decisions {dict(why)}")
    for kk, v in sorted(agree.items(), key=str):
        print("  ", kk, v)
    print(f"  identity x coverage >= 0.999: before {ok_before}, after the trim {ok_after}")


if __name__ == "__main__":
    main(sys.argv[1])
