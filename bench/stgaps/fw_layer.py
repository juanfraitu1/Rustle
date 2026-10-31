#!/usr/bin/env python3
"""Family layer (FW): singleton confirmation across copies (PREREG_family_layer_2026-10-26.md). Stages:

  chains  COPIES.tsv BAM OUTDIR                    -> OUTDIR/chains.tsv (family, copy, chain, molecules)
  map     COPIES.tsv GENOME.fa OUTDIR              -> OUTDIR/reps.fa, copies.fa (then, outside:
                                                      minimap2 -x asm20 -c --cs -N 50 reps.fa copies.fa > OUTDIR/map.paf)
  admit   COPIES.tsv OUTDIR BASE.gtf MIN1.gtf FW.gtf -> FW.gtf + OUTDIR/candidates.tsv

Copy loci are `locus_start..locus_end` of fam.copies.tsv; a copy is (family_id, copy_idx). A chain is the tuple of
0-based half-open introns; a molecule is a primary read de-duplicated on (chain, start, end). A junction lifts when its
last exon base and next exon base both sit on aligned (match or mismatch) columns of the copy's best hit to its family
representative (the copy with most n_reads); the family junction is the sorted pair of representative-locus offsets.
"""
import bisect
import os
import re
import sys
from collections import defaultdict

RT = re.compile(r'transcript_id "([^"]+)"')


def read_copies(path):
    """-> {(fam, idx): (chrom, locus_start, locus_end, n_reads)}"""
    out, hdr = {}, None
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if hdr is None:
            hdr = {k: i for i, k in enumerate(f)}
            continue
        out[(f[hdr["family_id"]], f[hdr["copy_idx"]])] = (
            f[hdr["chrom"]], int(f[hdr["locus_start"]]), int(f[hdr["locus_end"]]), int(f[hdr["n_reads"]]))
    return out


def chain_of(rec):
    out, r = [], rec.reference_start
    for op, n in rec.cigartuples:
        if op == 3:
            out.append((r, r + n))
        if op in (0, 2, 3, 7, 8):
            r += n
    return tuple(out)


def fmt(chain):
    return ",".join(f"{a}-{b}" for a, b in chain)


def parse(s):
    return tuple(tuple(int(x) for x in p.split("-")) for p in s.split(","))


def chains(copies_path, bam_path, out):
    import pysam
    os.makedirs(out, exist_ok=True)
    copies = read_copies(copies_path)
    by_c = defaultdict(list)
    for k, (c, a, b, _) in copies.items():
        by_c[c].append((a, b, k))
    bam = pysam.AlignmentFile(bam_path, "rb")
    mol = defaultdict(set)  # (fam, idx, chain) -> {(start, end)}
    for c in sorted(by_c):
        loci = sorted(by_c[c])
        merged = []  # [start, end, [loci]]
        for a, b, k in loci:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
                merged[-1][2].append((a, b, k))
            else:
                merged.append([a, b, [(a, b, k)]])
        for a, b, ls in merged:
            starts = [x[0] for x in ls]
            for rec in bam.fetch(c, a, b):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.is_qcfail or rec.cigartuples is None:
                    continue
                if not any(op == 3 for op, _ in rec.cigartuples):
                    continue
                s, e = rec.reference_start, rec.reference_end
                # a read spanning two merged regions is fetched twice; the molecule set absorbs the repeat
                ch = chain_of(rec)
                for la, lb, k in ls[:bisect.bisect_left(starts, e)]:
                    if lb > s:
                        mol[(k[0], k[1], ch)].add((s, e))
        print(f"[chains] {c}: {len(merged)} regions", file=sys.stderr, flush=True)
    with open(os.path.join(out, "chains.tsv"), "w") as fo:
        for (fam, idx, ch), m in sorted(mol.items()):
            fo.write(f"{fam}\t{idx}\t{fmt(ch)}\t{len(m)}\n")
    print(f"[chains] {len(mol)} (copy, chain) keys", file=sys.stderr)


def reps_of(copies):
    fams = defaultdict(list)
    for (fam, idx), v in copies.items():
        fams[fam].append((idx, v))
    return {fam: max(xs, key=lambda x: (x[1][3], -int(x[0])))[0] for fam, xs in fams.items() if len(xs) >= 2}


def map_(copies_path, fa_path, out):
    import pysam
    os.makedirs(out, exist_ok=True)
    copies = read_copies(copies_path)
    reps = reps_of(copies)
    fa = pysam.FastaFile(fa_path)
    with open(os.path.join(out, "reps.fa"), "w") as fr, open(os.path.join(out, "copies.fa"), "w") as fc:
        for (fam, idx), (c, a, b, _) in sorted(copies.items()):
            if fam not in reps:
                continue
            seq = fa.fetch(c, a, b).upper()
            (fr if reps[fam] == idx else fc).write(f">{fam}|{idx}\n{seq}\n")


def load_lift(paf):
    """-> {(fam, idx): (q_starts, blocks)}; block = (q_lo, q_hi, t0, strand): q in [q_lo, q_hi) maps to
    t0 + (q - q_lo) for '+' and t0 + (q_hi - 1 - q) for '-'. Best hit (most matching bases) to the own representative."""
    best = {}
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        qf, tf = f[0].split("|")[0], f[5].split("|")[0]
        if qf != tf:
            continue
        key = tuple(f[0].split("|"))
        nm = int(f[9])
        if key in best and best[key][0] >= nm:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        best[key] = (nm, f[4], int(f[2]), int(f[3]), int(f[7]), cs)
    out = {}
    tok = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}[0-9]+[a-z]{2})")
    for key, (_, strand, qs, qe, ts, cs) in best.items():
        blocks, qc, tc = [], (qs if strand == "+" else qe), ts

        def aligned(n):
            nonlocal qc, tc
            if strand == "+":
                blocks.append((qc, qc + n, tc, "+"))
                qc += n
            else:
                blocks.append((qc - n, qc, tc, "-"))
                qc -= n
            tc += n

        for t in tok.findall(cs):
            if t[0] == ":":
                aligned(int(t[1:]))
            elif t[0] == "*":
                aligned(1)
            elif t[0] == "+":
                qc += len(t) - 1 if strand == "+" else -(len(t) - 1)
            elif t[0] == "-":
                tc += len(t) - 1
        blocks.sort()
        out[key] = ([b[0] for b in blocks], blocks)
    return out


def lift(lf, q):
    qst, blocks = lf
    i = bisect.bisect_right(qst, q) - 1
    if i < 0:
        return None
    lo, hi, t0, s = blocks[i]
    if not lo <= q < hi:
        return None
    return t0 + (q - lo) if s == "+" else t0 + (hi - 1 - q)


def fam_chain(copies, reps, lifts, fam, idx, chain):
    c, a, b, _ = copies[(fam, idx)]
    pts = []
    for x, y in chain:
        if reps.get(fam) == idx:
            p, q = x - 1 - a, y - a
            if p < 0 or q > b - a:        # junction outside the representative's locus: no family coordinate
                return None
        else:
            lf = lifts.get((fam, idx))
            if lf is None:
                return None
            p, q = lift(lf, x - 1 - a), lift(lf, y - a)
            if p is None or q is None:
                return None
        pts.append((min(p, q), max(p, q)))
    return tuple(sorted(pts))


def gtf_tx(path):
    ex, key = defaultdict(list), {}
    for l in open(path):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = RT.search(f[8]).group(1)
        ex[t].append((int(f[3]), int(f[4])))
        key[t] = f[0]
    out = {}
    for t, e in ex.items():
        if len(e) < 2:
            continue
        e.sort()
        out[t] = (key[t], e[0][0] - 1, e[-1][1], tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)))
    return out


def admit(copies_path, out, base_gtf, min1_gtf, fw_gtf):
    copies = read_copies(copies_path)
    reps = reps_of(copies)
    lifts = load_lift(os.path.join(out, "map.paf"))
    own = {}
    fam_sup = defaultdict(int)
    for l in open(os.path.join(out, "chains.tsv")):
        fam, idx, ch, n = l.rstrip("\n").split("\t")
        ch, n = parse(ch), int(n)
        own[(fam, idx, ch)] = n
        if fam in reps:
            fc = fam_chain(copies, reps, lifts, fam, idx, ch)
            if fc is not None:
                fam_sup[(fam, fc)] += n
    base = {(c, ch) for c, _, _, ch in gtf_tx(base_gtf).values()}
    by_c = defaultdict(list)
    for k, (c, a, b, _) in copies.items():
        by_c[c].append((a, b, k))
    n_copies = defaultdict(int)
    for fam, _ in copies:
        n_copies[fam] += 1
    rank = {"confirmed": 4, "unconfirmed": 3, "unlifted": 2, "own2plus": 1, "own0": 0}
    rows, confirmed = [], {}
    for t, (c, s, e, ch) in sorted(gtf_tx(min1_gtf).items()):
        if (c, ch) in base:
            continue
        hits = [k for a, b, k in by_c.get(c, ()) if a < e and b > s]
        if not hits:
            continue
        best = ("own0", 0, 0, "", 0)
        for fam, idx in hits:
            o = own.get((fam, idx, ch), 0)
            fc = fam_chain(copies, reps, lifts, fam, idx, ch) if fam in reps else None
            fs = fam_sup.get((fam, fc), 0) if fc is not None else 0
            if o == 1:
                st = "confirmed" if fs >= 2 else ("unlifted" if fc is None else "unconfirmed")
            else:
                st = "own0" if o == 0 else "own2plus"
            if rank[st] > rank[best[0]] or (st == best[0] and fs > best[2]):
                best = (st, o, fs, fam, n_copies[fam])
        rows.append((t, c, s, e) + best)
        if best[0] == "confirmed":
            confirmed[t] = best
    with open(os.path.join(out, "candidates.tsv"), "w") as fo:
        fo.write("tid\tchrom\tstart\tend\tstatus\town\tfamily_support\tfamily\tfamily_copies\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    with open(fw_gtf, "w") as fo:
        for l in open(base_gtf):
            fo.write(l)
        for l in open(min1_gtf):
            m = RT.search(l)
            if not m or m.group(1) not in confirmed:
                continue
            st, o, fs, fam, n_cp = confirmed[m.group(1)]
            l = l.rstrip("\n").replace(f'transcript_id "{m.group(1)}"', f'transcript_id "FW.{m.group(1)}"')
            fo.write(f'{l} family_support "{fs}"; family "{fam}";\n')
    from collections import Counter
    print(f"[admit] candidates {len(rows)}: {dict(Counter(r[4] for r in rows))}; admitted {len(confirmed)}",
          file=sys.stderr)


if __name__ == "__main__":
    cmd, *a = sys.argv[1:]
    {"chains": chains, "map": map_, "admit": admit}[cmd](*a)
