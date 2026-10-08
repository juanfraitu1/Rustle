#!/usr/bin/env python3
"""Amendment 6 prep (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md): per family of the 2026-08-14 excision panel, the scored
read set S_f (<= 500 D reads + <= 500 K reads, seed 1, from baseline primaries), the IsoCon input net_f (S_f reads with a masked-arm
record overlapping K, plus S_f reads unmapped in the masked arm), one FASTA per family, the union FASTA to realign, and the labels.

    python3 heldout_prep.py --panel panel.json --baseline panel_primary.bam --masked masked.bam --fq panel_reads.fq --out-dir heldout
"""
import argparse
import collections
import json
import os
import random

import pysam


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("panel", "baseline", "masked", "fq", "out_dir"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    os.makedirs(f"{a.out_dir}/fam", exist_ok=True)
    panel = json.load(open(a.panel))
    base = pysam.AlignmentFile(a.baseline)
    rng = random.Random(1)
    S, role = {}, {}
    for p in panel:
        f = p["fam"]
        got = {}
        for lab in ("mask", "keep"):
            c, s, e = p[lab]
            names = sorted({rd.query_name for rd in base.fetch(c, s, e)
                            if not (rd.is_unmapped or rd.is_secondary or rd.is_supplementary)})
            rng.shuffle(names)
            got[lab] = names[:500]
        S[f] = got
        for lab, ns in got.items():
            for n in ns:
                role.setdefault(n, (f, "D" if lab == "mask" else "K"))
    want = set(role)
    # masked-arm records of the scored reads
    rec = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(a.masked).fetch(until_eof=True):
        if rd.query_name in want:
            rec[rd.query_name].append(("*", -1, -1) if rd.is_unmapped else (rd.reference_name, rd.reference_start, rd.reference_end))
    seq = {}
    with open(a.fq) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            n = h[1:].split()[0]
            if n in want:
                seq[n] = s
    pk = {p["fam"]: p["keep"] for p in panel}
    nets = {}
    with open(f"{a.out_dir}/labels.tsv", "w") as lab, open(f"{a.out_dir}/scored.fa", "w") as allfa:
        lab.write("read\tfamily\trole\tin_net\n")
        for f, got in S.items():
            c, s, e = pk[f]
            net = []
            for n in got["mask"] + got["keep"]:
                rs = rec.get(n, [])
                unm = bool(rs) and all(x[0] == "*" for x in rs)
                onk = any(x[0] == c and x[1] < e and s < x[2] for x in rs)
                if (onk or unm) and n in seq:
                    net.append(n)
            nets[f] = net
            with open(f"{a.out_dir}/fam/{f}.fa", "w") as o:
                for n in net:
                    o.write(f">{n}\n{seq[n]}\n")
            ns = set(net)
            for n in got["mask"] + got["keep"]:
                if role[n][0] != f or n not in seq:
                    continue
                lab.write(f"{n}\t{f}\t{role[n][1]}\t{int(n in ns)}\n")
                allfa.write(f">{n}\n{seq[n]}\n")
    nD = sum(len(g["mask"]) for g in S.values()); nK = sum(len(g["keep"]) for g in S.values())
    print(f"families {len(S)}; scored D {nD}, K {nK}; sequences found {len(seq)} of {len(want)}; IsoCon inputs: "
          f"total {sum(len(v) for v in nets.values())}, median {sorted(len(v) for v in nets.values())[len(nets) // 2]}, "
          f"max {max(len(v) for v in nets.values())}")


if __name__ == "__main__":
    main()
