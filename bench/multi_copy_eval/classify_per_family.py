#!/usr/bin/env python3
"""Classify a query GTF (rustle/ST output) against per-family ref GFFs
(GOLGA6/GOLGA8/AMY/TBC1D3).

For each ref family, run gffcompare and parse the .stats summary; print:
  family | ref_genes | ref_mRNAs | query_loci | matching_loci |
         | sn_locus | pr_locus | matching_intron_chains | matching_transcripts

Also output a count of how many distinct copies (loci) the query represents
in the family's chromosomal region — the central question of objective 3/5.

Usage:  classify_per_family.py <query.gtf> <label>
"""
import os, re, sys, subprocess, tempfile

if len(sys.argv) != 3:
    sys.exit("usage: classify_per_family.py <query.gtf> <label>")

query, label = sys.argv[1], sys.argv[2]
families = ["GOLGA6", "GOLGA8", "AMY", "TBC1D3"]

GFFCOMPARE = "/storage/work/jxi21/.pixi_cache/pkgs/gffcompare-0.12.10-h9948957_0/bin/gffcompare"

def parse_stats(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(r"#\s*Query mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci", line)
            if m:
                out["q_mRNAs"] = int(m.group(1))
                out["q_loci"]  = int(m.group(2))
                continue
            m = re.match(r"#\s*Reference mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci", line)
            if m:
                out["r_mRNAs"] = int(m.group(1))
                out["r_loci"]  = int(m.group(2))
                continue
            for k, lbl in [("Locus level", "locus"),
                           ("Intron chain level", "intron_chain"),
                           ("Transcript level", "transcript"),
                           ("Exon level", "exon")]:
                m = re.match(rf"\s*{re.escape(k)}\s*:\s*([\d.]+)\s*\|\s*([\d.]+)", line)
                if m:
                    out[f"sn_{lbl}"] = float(m.group(1))
                    out[f"pr_{lbl}"] = float(m.group(2))
            m = re.match(r"\s*Matching intron chains:\s*(\d+)", line)
            if m: out["match_ichains"] = int(m.group(1))
            m = re.match(r"\s*Matching transcripts:\s*(\d+)", line)
            if m: out["match_tx"] = int(m.group(1))
            m = re.match(r"\s*Matching loci:\s*(\d+)", line)
            if m: out["match_loci"] = int(m.group(1))
    return out

print("="*100)
print(f"  {label}: {query}")
print("="*100)
header = f"{'family':10} {'r_loci':>6} {'r_mRNA':>6} {'q_loci':>6} {'q_mRNA':>6} {'match_loci':>10} {'sn_locus':>8} {'pr_locus':>8} {'match_ichain':>12} {'match_tx':>8}"
print(header)
print("-"*100)

with tempfile.TemporaryDirectory() as tmp:
    for fam in families:
        ref = os.path.join(os.path.dirname(query) or ".", f"ref_{fam}.gff")
        if not os.path.exists(ref):
            ref = f"ref_{fam}.gff"
        if not os.path.exists(ref):
            print(f"{fam:10} ref_{fam}.gff missing")
            continue
        prefix = os.path.join(tmp, f"cmp_{fam}")
        cmd = [GFFCOMPARE, "-r", ref, "-o", prefix, "--no-merge", query]
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        stats = prefix + ".stats"
        s = parse_stats(stats)
        print(f"{fam:10} "
              f"{s.get('r_loci',0):>6} {s.get('r_mRNAs',0):>6} "
              f"{s.get('q_loci',0):>6} {s.get('q_mRNAs',0):>6} "
              f"{s.get('match_loci',0):>10} "
              f"{s.get('sn_locus','?'):>8} {s.get('pr_locus','?'):>8} "
              f"{s.get('match_ichains','?'):>12} {s.get('match_tx','?'):>8}")

print()
