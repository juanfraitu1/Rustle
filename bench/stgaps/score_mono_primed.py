#!/usr/bin/env python3
"""Scorer for PREREG_mono_primed_2026-10-26.md (clauses IV, H1-H3). Predicate = polya_mono_genomic.py's, frozen.

usage: score_mono_primed.py BASE.gtf REF.gtf BAM FASTA OUTPREFIX
"""
import os
import re
import subprocess
import sys

import pysam

base, ref, bam_path, fa_path, pre = sys.argv[1:6]
RE_TID = re.compile(r'transcript_id "([^"]+)"')
fa, bam = pysam.FastaFile(fa_path), pysam.AlignmentFile(bam_path, "rb")


def rich(s, b):
    return bool(s) and (s.count(b) >= 0.6 * len(s) or b * 6 in s)


def gffc(gtf, out):
    subprocess.run(["gffcompare", "-r", ref, "-o", out, gtf], check=True, capture_output=True)
    # gffcompare writes OUT.stats, or OUT itself when the prefix already has a dot-extension
    st = open(out + ".stats" if os.path.exists(out + ".stats") else out).read().splitlines()
    tx = [x.strip() for x in next(l for l in st if l.strip().startswith("Transcript level")).split(":")[1].split("|")[:2]]
    ch = [x.strip() for x in next(l for l in st if l.strip().startswith("Intron chain level")).split(":")[1].split("|")[:2]]
    return tx, ch


exons, lines = {}, []
for line in open(base):
    f = line.split("\t")
    lines.append(line)
    if len(f) >= 9 and f[2] == "exon":
        exons.setdefault(RE_TID.search(f[8]).group(1), []).append((f[0], int(f[3]) - 1, int(f[4])))
mono = {t: e[0] for t, e in exons.items() if len(e) == 1}

tx_b, ch_b = gffc(base, pre + ".gc_base")
d = base.rsplit("/", 1)[0] if "/" in base else "."
tmap = f"{pre}.gc_base.{base.rsplit('/', 1)[-1]}.tmap"
if not os.path.exists(tmap):
    tmap = os.path.join(d, f"{os.path.basename(pre)}.gc_base.{os.path.basename(base)}.tmap")
cls = {}
for l in open(tmap):
    f = l.split("\t")
    if f[0] != "ref_gene_id":
        cls[f[4]] = f[2]


def flagged(shift):
    out = set()
    for t, (ctg, s, e) in mono.items():
        tot = pr = 0
        for r in bam.fetch(ctg, s, e):
            if r.flag & 2820 or r.cigartuples is None or any(op == 3 for op, _ in r.cigartuples):
                continue
            if min(e, r.reference_end) - max(s, r.reference_start) < 0.5 * (e - s):
                continue
            a, b = r.reference_start + shift, r.reference_end + shift
            tot += 1
            pr += rich(fa.fetch(ctg, b, b + 20).upper(), "A") or rich(fa.fetch(ctg, max(0, a - 20), a).upper(), "T")
        if tot and pr / tot >= 0.8:
            out.add(t)
    return out


rule, null = flagged(0), flagged(1000)
drop_gtf = pre + ".DROP.gtf"
with open(drop_gtf, "w") as fo:
    for line in lines:
        f = line.split("\t")
        m = RE_TID.search(f[8]) if len(f) >= 9 else None
        if m and m.group(1) in rule:
            continue
        fo.write(line)
multi = lambda p: [l for l in open(p) if (m := RE_TID.search(l)) and len(exons.get(m.group(1), [])) > 1]
iv = multi(base) == multi(drop_gtf)
tx_d, ch_d = gffc(drop_gtf, pre + ".gc_drop")
eq_mono = [t for t in mono if cls.get(t) == "="]
eq_lost = [t for t in rule if cls.get(t) == "="]
h1 = len(rule) >= 5 * len(null)
h2 = len(eq_lost) <= max(1, int(0.1 * len(eq_mono)))
h3 = float(tx_d[1]) > float(tx_b[1])
print(f"{pre}: single-exon {len(mono)} (= {len(eq_mono)}); flagged {len(rule)}, NULL {len(null)}; '=' lost {len(eq_lost)}")
print(f"  tx Sn/Pr BASE {tx_b[0]}/{tx_b[1]} DROP {tx_d[0]}/{tx_d[1]}; chain BASE {ch_b} DROP {ch_d}")
print(f"  IV {iv}  H1 {h1}  H2 {h2}  H3 {h3}  => {'PASS' if iv and h1 and h2 and h3 else 'FAIL'}")
