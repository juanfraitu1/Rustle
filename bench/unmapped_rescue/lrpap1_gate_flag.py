#!/usr/bin/env python3
"""O3 flag of the LRPAP1 consensus sequences, registered and gate-aware (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 17). Uses the existing alignments
cons.mat.sam / cons.pat.sam of run_lrpap1.py place; python3 lrpap1_gate_flag.py [full|half|fullpol]"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import apply_trim_real as AT  # noqa: E402
import libsig  # noqa: E402
import lrpap1 as L  # noqa: E402
import run_lrpap1 as RL  # noqa: E402
import run_polish as RP  # noqa: E402
from o3_maternal import common as C  # noqa: E402


def records(sam_path, al, gate_open):
    out = {}
    for ln in open(sam_path):
        if ln[0] == "@":
            continue
        f = ln.rstrip("\n").split("\t")
        if int(f[1]) & 4 or f[2] == "*":
            continue
        n, clip = libsig.five_prime_clip(int(f[1]), f[5], f[9])
        rlen = sum(int(x) for x, op in RP.R_CIG(f[5]) if op in "M=XIS")
        al_len = sum(int(x) for x, op in RP.R_CIG(f[5]) if op in "M=XI")
        span = sum(int(x) for x, op in RP.R_CIG(f[5]) if op in "MDN=X")
        de = float(next(t[5:] for t in f[11:] if t.startswith("de:f:")))
        cov = al_len / rlen
        cov_g = al_len / (rlen - n) if (gate_open and 0 < n <= 3 and set(clip) == {"G"}) else cov
        out.setdefault(f[0].split("|")[0], []).append(dict(ref=C.accession(f[2], al), start=int(f[3]) - 1, end=int(f[3]) - 1 + span, de=de, qcov=cov, qcov_gate=cov_g, clip5=clip))
    return out


def main(tag):
    d = f"{RL.OUT}/{tag}"
    al = C.alias()
    gate_open, p, sig = AT.library_gate()
    rows = RL.loci_rows()
    names = [ln[1:].strip().split()[0] for ln in open(f"{d}/cons.fa") if ln[0] == ">"]
    maj = RL.G_cluster_majority(d)
    best = {}
    for h in ("mat", "pat"):
        iv = L.intervals(rows, h)
        recs = records(f"{d}/cons.{h}.sam", al, gate_open)
        for k, rs in recs.items():
            for variant, key in (("registered", "qcov"), ("gate", "qcov_gate")):
                pl = {}
                for r in rs:
                    lab = L.locus_at(iv, r["ref"], r["start"], r["end"])
                    if lab:
                        pl[lab] = max(pl.get(lab, 0.0), (1 - r["de"]) * r[key])
                top = max(pl.items(), key=lambda kv: kv[1]) if pl else (None, None)
                best[(k, h, variant)] = top
    print(f"library gate {'OPEN' if gate_open else 'closed'}")
    print(f"{'cluster':30s} {'maj':4s} {'5p clip':8s} {'pat best (reg / gate)':32s} {'mat best (reg / gate)':32s} flag reg / gate")
    for n in names:
        k = n.split("|")[0]
        pr, pg, mr, mg = best[(k, "pat", "registered")], best[(k, "pat", "gate")], best[(k, "mat", "registered")], best[(k, "mat", "gate")]
        clip = next((r["clip5"] for r in records(f"{d}/cons.pat.sam", al, gate_open).get(k, [])), "")
        print(f"{n:30s} {str(maj.get(n)):4s} {clip[:6]:8s} {pr[0]!s:4s} {pr[1]:.4f} / {pg[1]:.4f}      {mr[0]!s:8s} {mr[1]:.4f} / {mg[1]:.4f}      "
              f"{L.o3_flag(mr[1], pr[1])} / {L.o3_flag(mg[1], pg[1])}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "full")
