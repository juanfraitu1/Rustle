#!/usr/bin/env python3
"""Flag Q (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 21): is a cluster's consensus supported by its own reads?
    consensus_support.py <bedA|bedH>      (miniforge python; needs registered/{cons.fa,clusters.tsv} and polish/score.json)
A consensus is UNSUPPORTED if the median divergence of its reads (<= 100, splice:hq -uf) aligned to it exceeds 0.00958."""
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_augment as RA  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
DELTA = 0.00958


def median_read_divergence(cons, reads, tmp):
    os.makedirs(tmp, exist_ok=True)
    open(f"{tmp}/r.fa", "w").write("".join(f">{n}\n{s}\n" for n, s in reads.items()))
    open(f"{tmp}/c.fa", "w").write(f">c\n{cons}\n")
    sam = subprocess.run(f"minimap2 -ax splice:hq -uf --eqx -t 2 {tmp}/c.fa {tmp}/r.fa", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
    prim = RA.primaries(sam)
    des = [prim[n]["de"] for n in reads if n in prim]
    return statistics.median(des) if des else None


def main(bed):
    d = f"{W}/{bed}"
    reg = f"{d}/registered"
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{reg}/cons.fa").items()}
    cl = RP.cluster_reads(d)
    seqs = SD.read_fa(f"{d}/pool.fa")
    sc = json.load(open(f"{reg}/polish/score.json"))["clusters"]
    out = {}
    for k, v in sc.items():
        if not v["orig"]["on"] or len(cl[k]) < 3:
            continue
        names = sorted(cl[k])[:100]
        med = median_read_divergence(cons[k], {n: seqs[n] for n in names}, f"/tmp/cs_{bed}")
        out[k] = dict(reads=len(cl[k]), median_de=med, nm=v["orig"]["nm"] or 0, unsupported=med is None or med > DELTA)
    json.dump(out, open(f"{reg}/consensus_support.json", "w"), indent=1)
    big = {k for k, o in out.items() if o["nm"] > 200}
    flag = {k for k, o in out.items() if o["unsupported"]}
    rest = len(out) - len(big)
    print(f"{bed}: clusters {len(out)}; flagged {len(flag)}; clusters with > 200 edits {len(big)}, flagged {len(flag & big)}; "
          f"other clusters flagged {len(flag - big)} of {rest} ({len(flag - big) / max(1, rest):.1%}); flagged sizes {sorted(out[k]['reads'] for k in flag)}")


if __name__ == "__main__":
    main(sys.argv[1])
