#!/usr/bin/env python3
"""abPOA consensus per survivor copy. Run with /home/juanfra/miniforge3/bin/python.   survivor_consensus.py <surv_reads.fa (names copy|read)> <out.fa>"""
import sys


def main(inp, out):
    import pyabpoa
    by, cur = {}, None
    for ln in open(inp):
        if ln[0] == ">":
            cur = ln[1:].strip().split("|")[0]
            by.setdefault(cur, []).append([])
        else:
            by[cur][-1].append(ln.strip())
    aln = pyabpoa.msa_aligner(aln_mode="g", is_aa=False, cons_algrm="HB")
    n = 0
    with open(out, "w") as o:
        for name, reads in sorted(by.items()):
            seqs = sorted(("".join(r) for r in reads), key=lambda s: (-len(s), s))[:100]
            if len(seqs) < 3:
                continue
            res = aln.msa(seqs, out_cons=True, out_msa=False)
            if res.cons_seq:
                o.write(f">{name}\n{res.cons_seq[0]}\n")
                n += 1
    print("survivor consensus sequences:", n)


if __name__ == "__main__":
    main(*sys.argv[1:3])
