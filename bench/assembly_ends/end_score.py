#!/usr/bin/env python3
"""PREREG ASSEMBLY ENDS E3 (register 1487) scorer: 5' / 3' end accuracy of exactly matched multi-exon transcripts.
  end_score.py REF.gtf QUERY.gtf GFFCOMPARE_PREFIX
Over gffcompare `=` pairs (query transcript, reference transcript) with >= 2 exons: M = share whose 5' end lies within 50 bp of the
reference's 5' end (also the 3' share and medians); plus the stats file's base-level sensitivity / precision, matching intron chains
and matching transcripts. One TSV line."""
import re, sys, statistics as st

def ends(path):
    ex, strand = {}, {}
    for l in open(path):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search('transcript_id "([^"]+)"', f[8]).group(1)
        ex.setdefault(t, []).append((int(f[3]), int(f[4])))
        strand[t] = f[6]
    out = {}
    for t, e in ex.items():
        lo, hi = min(a for a, _ in e), max(b for _, b in e)
        out[t] = (len(e), (hi, lo) if strand[t] == "-" else (lo, hi))   # (n exons, (5' end, 3' end))
    return out

ref, qry, pre = sys.argv[1:4]
R, Q = ends(ref), ends(qry)
d5, d3 = [], []
for l in open(f"{pre}.tracking"):
    f = l.rstrip("\n").split("\t")
    if f[3] != "=" or f[2] == "-":
        continue
    rt = f[2].split("|")[1]
    for q in f[4:]:
        if q == "-":
            continue
        parts = q.split(":", 1)[1].split("|")       # gene|transcript|exons|FPKM|TPM|cov|len; a transcript id may hold '|'
        qt = "|".join(parts[1:-5])
        if qt in Q and rt in R and Q[qt][0] >= 2:
            d5.append(abs(Q[qt][1][0] - R[rt][1][0])); d3.append(abs(Q[qt][1][1] - R[rt][1][1]))
stats = open(pre).read()
num = lambda pat: re.search(pat, stats).group(1)
base_sn, base_pr = re.search(r"Base level:\s+([\d.]+)\s+\|\s+([\d.]+)", stats).groups()
print("\t".join(map(str, [qry, len(d5), f"{100 * sum(x <= 50 for x in d5) / len(d5):.2f}", st.median(d5),
                          f"{100 * sum(x <= 50 for x in d3) / len(d3):.2f}", st.median(d3), base_sn, base_pr,
                          num(r"Matching intron chains:\s+(\d+)"), num(r"Matching transcripts:\s+(\d+)")])))
