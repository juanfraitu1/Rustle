#!/usr/bin/env python3
"""The gates of docs/PREREG_ideal_expression_2026-10-06.md for one arm, read from the files the runner wrote; prints VALID / INVALID with the reason for each gate.

    gates.py --dir W/FAMILY/repR --bam-reads N_FASTQ_READS_FROM_FQ
G1 g1.log; G2 strata.G2.json; G3 BAM record count == FASTQ reads, driver exit status 0, families.gtf newer than gtf; G4 two scorer runs under different PYTHONHASHSEED (done by the caller,
compared here from the two output prefixes given with --g4a/--g4b); G6 score_anno.summary.json (E2 and E3 equal |R| on R for the family); G7 reported against the control's value (no gate).
"""
import argparse
import json
import os
import re
import subprocess

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--g4a", required=True)
    ap.add_argument("--g4b", required=True)
    a = ap.parse_args()
    d = a.dir
    out = {}
    g1 = open(f"{d}/g1.log").read()
    out["G1"] = bool(re.search(r"0 violations", g1))
    g2 = json.load(open(f"{d}/strata.G2.json"))
    out["G2"] = bool(g2["ok"])
    n_fq = sum(1 for _ in open(f"{d}/reads.fq")) // 4
    bam = pysam.AlignmentFile(f"{d}/reads.bam")
    n_bam = sum(1 for r in bam.fetch(until_eof=True) if not (r.is_secondary or r.is_supplementary))
    exit0 = all("Exit status: 0" in open(f"{d}/asm.{s}.driver.stderr").read() for s in ("assemble", "families"))
    newer = os.path.getmtime(f"{d}/asm.families.gtf") >= os.path.getmtime(f"{d}/asm.gtf")
    out["G3"] = bool(n_fq == n_bam and exit0 and newer)
    same = all(open(f"{a.g4a}.{x}", "rb").read() == open(f"{a.g4b}.{x}", "rb").read() for x in (f"{a.family}.copies.tsv", "summary.json"))
    out["G4"] = bool(same)
    sa = json.load(open(f"{d}/score_anno.summary.json"))[a.family]
    out["G6"] = bool(sa["R"] > 0 and sa["R_counts"]["E2"] == sa["R"] and sa["R_counts"]["E3"] == sa["R"])
    arm = json.load(open(f"{d}/score.summary.json"))
    g7, g7c = arm.get("G7", {}), json.load(open(f"{d}/score_anno.summary.json")).get("G7", {})
    detail = dict(G1=g1.strip().splitlines()[-1], G2=f"{g2['primary_on_source']}/{g2['reads']} = {g2['share']}", G3=f"BAM {n_bam} vs FASTQ {n_fq}; exit0 {exit0}; families.gtf newer {newer}", G4="identical" if same else "DIFFERS",
                  G6=f"control E2 {sa['R_counts']['E2']}, E3 {sa['R_counts']['E3']}, E4 {sa['R_counts']['E4']} of {sa['R']} on R (E4* = {sa['R_counts']['E4'] / sa['R']:.2f})",
                  G7=f"arm {g7.get('exactly_one_locus')}/{g7.get('single_copy_genes')} = {g7.get('share_one')} against the annotation-as-loci control {g7c.get('exactly_one_locus')}/{g7c.get('single_copy_genes')} = {g7c.get('share_one')} (registered >= 0.95: not attainable even by the control)")
    valid = all(out[k] for k in ("G1", "G3", "G4", "G6")) 
    print(f"{a.family} {os.path.basename(d)}: " + ("VALID" if valid else "INVALID") + " | " + " | ".join(f"{k} {'ok' if out[k] else 'FAIL'}" for k in ("G1", "G2", "G3", "G4", "G6")))
    for k, v in detail.items():
        print(f"    {k}: {v}")
    json.dump(dict(valid=valid, gates=out, detail=detail), open(f"{d}/gates.json", "w"), indent=1)


if __name__ == "__main__":
    main()
