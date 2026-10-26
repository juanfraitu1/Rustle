#!/usr/bin/env python3
"""Cross-copy dominance prune of pool-added transcripts (PREREG_dominance_prune_2026-10-26.md). Stages:

  support   POOL_ADDED.tsv BAM MOLECULES.tsv OUTDIR      -> OUTDIR/support.tsv (cand, molecule), OUTDIR/names.txt
  tx        POOLED.gtf OUTDIR GENOME.fa                  -> OUTDIR/tx.fa (multi-exon transcripts), OUTDIR/tx_span.tsv
  (then, outside: samtools view -F 0x900 -N names.txt BAM -> reads.fa; minimap2 -x map-hifi -c -N 100 -p 0.8 tx.fa reads.fa)
  decide    OUTDIR ALN.paf POOLED.gtf OUT.gtf            -> pruned GTF + OUTDIR/decisions.tsv

Support = molecules with a record at the candidate's contig whose exact intron chain is the candidate's: a primary, or a
secondary with AS >= 0.98 x the molecule's genome-wide best AS (as_table column 2), the pool's admission. A molecule is
won elsewhere when its best AS on a transcript NOT overlapping the candidate's span is strictly greater than its AS on the
candidate (no hit = -inf). A candidate with every support molecule won elsewhere is dropped; ties keep it.
"""
import os
import re
import sys
from collections import defaultdict

RT = re.compile(r'transcript_id "([^"]+)"')
RATIO = 0.98


def chain_of(rec):
    out, r = [], rec.reference_start
    for op, n in rec.cigartuples:
        if op == 3:
            out.append((r, r + n))
        if op in (0, 2, 3, 7, 8):
            r += n
    return tuple(out)


def support(added, bam_path, mol_path, out):
    import pysam
    os.makedirs(out, exist_ok=True)
    cands = []
    for l in open(added):
        t, c, s, ch = l.rstrip("\n").split("\t")
        chain = tuple(tuple(int(x) for x in p.split("-")) for p in ch.split(","))
        cands.append((t, c, chain))
    bam = pysam.AlignmentFile(bam_path, "rb")
    pairs, sec = [], defaultdict(list)  # sec: name -> [(cand, AS)] awaiting the genome-wide best
    # one pass per merged region (candidate spans within 100 kb merged, per contig): every record is read once and its
    # chain looked up among the contig's candidate chains (per-donor fetches rescanned the same deep piles thousands of times)
    by_chain = defaultdict(lambda: defaultdict(list))  # contig -> chain -> [cand]
    spans = defaultdict(list)
    for t, c, chain in cands:
        by_chain[c][chain].append(t)
        spans[c].append((chain[0][0], chain[-1][1]))
    for c in sorted(spans):
        iv = sorted(spans[c])
        merged = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= merged[-1][1] + 100_000:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        want = by_chain[c]
        seen = set()
        for a, b in merged:
            for rec in bam.fetch(c, a, b):
                if rec.is_unmapped or rec.is_supplementary or rec.cigartuples is None:
                    continue
                k = (rec.query_name, rec.reference_start, rec.flag)
                if k in seen:  # a long record overlapping two merged regions
                    continue
                seen.add(k)
                if not any(op == 3 for op, _ in rec.cigartuples):
                    continue
                ts = want.get(chain_of(rec))
                if not ts:
                    continue
                for t in ts:
                    if rec.is_secondary:
                        sec[rec.query_name].append((t, rec.get_tag("AS") if rec.has_tag("AS") else 0))
                    else:
                        pairs.append((t, rec.query_name))
        print(f"[support] {c}: {len(merged)} regions, {len(pairs)} pairs so far", file=sys.stderr, flush=True)
    for l in open(mol_path):
        if l.startswith("#"):
            continue
        f = l.split("\t", 2)
        if f[0] in sec:
            best = int(f[1])
            for t, a in sec.pop(f[0]):
                if a >= RATIO * best:
                    pairs.append((t, f[0]))
    pairs = sorted(set(pairs))
    with open(os.path.join(out, "support.tsv"), "w") as fo:
        for t, m in pairs:
            fo.write(f"{t}\t{m}\n")
    names = sorted({m for _, m in pairs})
    with open(os.path.join(out, "names.txt"), "w") as fo:
        fo.write("\n".join(names) + "\n")
    print(f"[support] {len(cands)} candidates, {len({t for t, _ in pairs})} with support, {len(pairs)} pairs, "
          f"{len(names)} molecules", file=sys.stderr)


def tx(gtf, out, fa_path):
    import pysam
    ex, key = defaultdict(list), {}
    for l in open(gtf):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = RT.search(f[8]).group(1)
        ex[t].append((int(f[3]) - 1, int(f[4])))
        key[t] = (f[0], f[6])
    fa = pysam.FastaFile(fa_path)
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    with open(os.path.join(out, "tx.fa"), "w") as fo, open(os.path.join(out, "tx_span.tsv"), "w") as fs:
        for t, e in ex.items():
            if len(e) < 2:
                continue
            e.sort()
            c, s = key[t]
            seq = "".join(fa.fetch(c, a, b) for a, b in e).upper()
            if s == "-":
                seq = seq[::-1].translate(comp)
            fo.write(f">{t}\n{seq}\n")
            fs.write(f"{t}\t{c}\t{e[0][0]}\t{e[-1][1]}\n")


def decide(out, paf, gtf, out_gtf):
    span = {}
    for l in open(os.path.join(out, "tx_span.tsv")):
        t, c, a, b = l.rstrip("\n").split("\t")
        span[t] = (c, int(a), int(b))
    sup = defaultdict(list)
    for l in open(os.path.join(out, "support.tsv")):
        t, m = l.rstrip("\n").split("\t")
        sup[t].append(m)
    want = {m for ms in sup.values() for m in ms}
    best = defaultdict(dict)  # molecule -> transcript -> best AS
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        if f[0] not in want:
            continue
        a = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), None)
        if a is None:
            continue
        d = best[f[0]]
        if a > d.get(f[5], -1):
            d[f[5]] = a
    drop, rows = set(), []
    for t in sorted(sup):
        c, a, b = span[t]
        won = 0
        for m in sup[t]:
            d = best.get(m, {})
            own = d.get(t, float("-inf"))
            other = max((v for u, v in d.items() if not (span[u][0] == c and span[u][1] < b and span[u][2] > a)),
                        default=float("-inf"))
            won += other > own
        dec = "drop" if won == len(sup[t]) else "keep"
        if dec == "drop":
            drop.add(t)
        rows.append((t, len(sup[t]), won, dec))
    with open(os.path.join(out, "decisions.tsv"), "w") as fo:
        fo.write("tid\tsupport\twon_elsewhere\tdecision\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    with open(out_gtf, "w") as fo:
        for l in open(gtf):
            m = RT.search(l)
            if m and m.group(1) in drop:
                continue
            fo.write(l)
    print(f"[decide] {len(rows)} candidates with support, dropped {len(drop)}", file=sys.stderr)


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    {"support": support, "tx": tx, "decide": decide}[cmd](*a)
