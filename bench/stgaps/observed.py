#!/usr/bin/env python3
"""Recall on OBSERVED reference transcripts (descriptive; FAMILY_LAYER_SPEC Part A follow-up).

  count  BAM OUT.tsv                                   -> (contig, chain, molecules, mapq0_molecules), primary spliced reads
  score  REF.gtf STRATUM.tsv COUNTS.tsv name=TMAP ...  -> recall per stratum x observation level

A molecule is a primary spliced read (-F 0x904) de-duplicated on (chain, start, end). A reference multi-exon transcript
is observed at level k when >= k molecules carry exactly its intron chain on its contig (k = 1: some read shows the
whole chain; k = 2: the assembler's pass-1 floor). `mapq0` = share of those molecules with MAPQ 0 (placement a coin
toss between copies). Recall = reference transcripts with a gffcompare `=` query (tmap) / reference transcripts.
"""
import re
import sys
from collections import Counter, defaultdict

RT = re.compile(r'transcript_id "([^"]+)"')


def count(bam_path, out):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb", threads=4)
    mol, q0 = defaultdict(set), defaultdict(set)
    n = 0
    for rec in bam.fetch(until_eof=True):
        if rec.flag & 0x904 or rec.is_unmapped:
            continue
        ct = rec.cigartuples
        if not ct or not any(op == 3 for op, _ in ct):
            continue
        r, ch = rec.reference_start, []
        for op, k in ct:
            if op == 3:
                ch.append(f"{r}-{r + k}")
            if op in (0, 2, 3, 7, 8):
                r += k
        key = (rec.reference_name, ",".join(ch))
        se = (rec.reference_start, r)
        mol[key].add(se)
        if rec.mapping_quality == 0:
            q0[key].add(se)
        n += 1
        if n % 2_000_000 == 0:
            print(f"[count] {n} spliced primaries, {len(mol)} chains", file=sys.stderr, flush=True)
    with open(out, "w") as fo:
        for (c, ch), m in mol.items():
            fo.write(f"{c}\t{ch}\t{len(m)}\t{len(q0.get((c, ch), ()))}\n")
    print(f"[count] {n} spliced primaries, {len(mol)} chains", file=sys.stderr)


def score(ref, stratum_path, counts, *arms):
    ex, key = defaultdict(list), {}
    for l in open(ref):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = RT.search(f[8])
        if not m:
            continue
        ex[m.group(1)].append((int(f[3]), int(f[4])))
        key[m.group(1)] = f[0]
    chain = {}
    for t, e in ex.items():
        if len(e) < 2:
            continue
        e.sort()
        chain[t] = (key[t], ",".join(f"{e[i][1]}-{e[i + 1][0] - 1}" for i in range(len(e) - 1)))
    want = defaultdict(list)
    for t, k in chain.items():
        want[k].append(t)
    obs, q0 = {}, {}
    for l in open(counts):
        c, ch, n, z = l.rstrip("\n").split("\t")
        for t in want.get((c, ch), ()):
            obs[t], q0[t] = int(n), int(z)
    stratum = {}
    for l in open(stratum_path):
        t, g, s = l.rstrip("\n").split("\t")
        stratum[t] = s
    levels = [("all", lambda t: True), ("obs>=1", lambda t: obs.get(t, 0) >= 1), ("obs>=2", lambda t: obs.get(t, 0) >= 2),
              ("obs>=2,mapq0<half", lambda t: obs.get(t, 0) >= 2 and q0[t] * 2 < obs[t]),
              ("obs>=2,mapq0>=half", lambda t: obs.get(t, 0) >= 2 and q0[t] * 2 >= obs[t])]
    hits = {}
    for a in arms:
        name, tm = a.split("=", 1)
        h = set()
        for l in open(tm):
            f = l.split("\t")
            if f[2] == "=":
                h.add(f[1])
        hits[name] = h
    print("stratum\tlevel\tref_tx\t" + "\t".join(f"{n} R%" for n in hits))
    for s in ("single", "multi"):
        for lv, fn in levels:
            ts = [t for t in chain if stratum.get(t) == s and fn(t)]
            row = [s, lv, str(len(ts))]
            for n, h in hits.items():
                row.append(f"{100 * sum(1 for t in ts if t in h) / len(ts):.1f}" if ts else "-")
            print("\t".join(row))


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    {"count": count, "score": score}[cmd](*a)
