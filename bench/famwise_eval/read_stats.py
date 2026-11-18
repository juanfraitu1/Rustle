#!/usr/bin/env python3
"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md): read statistics of a selected-reads BAM (calibration, and
item E of the mimic bar: the share of reads with both ends complete against an annotation, strict rule).

  read_stats.py SEL.bam REF.gtf [ASSEMBLY.gtf] [--out STATS.tsv] [--pct 14]

Counts every record; the per-read statistics are of the primaries whose name hashes into the first --pct percent, so
the sample does not depend on the order of the file as long as it stays under the cap of 20,000 reads.
Strict rule: a read is judged on the transcript of its strand with the most exonic overlap (tie: the shortest); the 5'
end is complete when the read starts inside or before the end of the first exon, the 3' end when it reaches the last.
--out: stat, value (reads, read length, events per aligned base, secondaries, the three strict shares on REF.gtf)."""
import argparse, bisect, collections, re, sys, time, zlib
import pysam

CAP = 20000

def q(v, ps=(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), f="{:.4g}"):
    v = sorted(v)
    return " ".join(f"p{int(p * 100)}=" + f.format(v[min(len(v) - 1, int(p * len(v)))]) for p in ps) if v else "NA"

def blocks(start, cigar):
    """aligned blocks, merged over insertions and deletions, split at introns"""
    pos = s = start; b = []
    for op, l in cigar:
        if op == 3:
            b.append((s, pos)); pos += l; s = pos
        elif op in (0, 2, 7, 8):
            pos += l
    return b + [(s, pos)]

def sample(bam, pct):
    """-> (sampled primaries, record counts, secondaries per read name)"""
    reads, n, secper, over = [], collections.Counter(), collections.Counter(), 0
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        for r in fh:
            n["rec"] += 1
            if r.is_unmapped: n["unm"] += 1; continue
            if r.is_supplementary: n["sup"] += 1; continue
            if r.is_secondary: n["sec"] += 1; secper[r.query_name] += 1; continue
            n["prim"] += 1
            if zlib.crc32(r.query_name.encode()) % 100 >= pct: continue
            if len(reads) >= CAP: over += 1; continue
            c, ev, il, dl = (collections.Counter() for _ in range(4))
            for op, l in r.cigartuples:
                c[op] += l; ev[op] += 1
                if op == 1: il[min(l, 10)] += 1
                if op == 2: dl[min(l, 10)] += 1
            ends = (r.cigartuples[-1], r.cigartuples[0]) if r.is_reverse else (r.cigartuples[0], r.cigartuples[-1])
            reads.append(dict(chrom=r.reference_name, strand="-" if r.is_reverse else "+", blocks=blocks(r.reference_start, r.cigartuples),
                              qlen=r.infer_read_length(), c=c, ev=ev, il=il, dl=dl, mapq=r.mapping_quality,
                              de=r.get_tag("de") if r.has_tag("de") else None, clip5=ends[0], clip3=ends[1]))
    if over:
        print(f"read_stats: {over} hash-passing primaries beyond the cap of {CAP}: the sample depends on file order, lower --pct", file=sys.stderr)
    return reads, n, secper

def rate(r):
    return (r["c"][8] + r["ev"][1] + r["ev"][2]) / (r["c"][7] + r["c"][8])

def report(reads, n, secper):
    """the printed read report -> the values of --out"""
    sp = list(secper.values()); z = n["prim"] - len(secper)
    print("secondary records per primary read: mean", n["sec"] / n["prim"], "; share of reads with none", z / n["prim"])
    print("secondaries per primary: 0:", z, " among reads with any:", q(sp, (0.25, 0.5, 0.75, 0.9, 0.99)))
    L = [r["qlen"] for r in reads]; print("read length", q(L), "mean", sum(L) / len(L), "min", min(L), "max", max(L))
    print("aligned query bases", q([r["c"][7] + r["c"][8] + r["c"][1] for r in reads]))
    print("read exon count", sorted(collections.Counter(min(len(r["blocks"]), 15) for r in reads).items()))
    print("mapq", sorted(collections.Counter(0 if r["mapq"] == 0 else 1 if r["mapq"] < 10 else 10 if r["mapq"] < 60 else 60 for r in reads).items()))
    c5, c3 = ([r[k][1] if r[k][0] == 4 else 0 for r in reads] for k in ("clip5", "clip3"))
    clip = lambda c: f"share>0 {sum(1 for x in c if x) / len(c)} >=20 {sum(1 for x in c if x >= 20) / len(c)} {q(c, (0.5, 0.9, 0.99))}"
    print("softclip 5p:", clip(c5), "| 3p:", clip(c3))
    T, E, IL, DL = (collections.Counter() for _ in range(4))
    for r in reads:
        T.update(r["c"]); E.update(r["ev"]); IL.update(r["il"]); DL.update(r["dl"])
    ab = T[7] + T[8]; den = ab + T[1] + T[2]    # aligned bases; alignment columns
    print("CIGAR totals: = %d X %d I bases %d (events %d) D bases %d (events %d) N %d S %d; M ops %d" % (T[7], T[8], T[1], E[1], T[2], E[2], T[3], T[4], T[0]))
    print(f"per alignment column (=+X+I+D): mismatch {T[8] / den:.5%}, ins bases {T[1] / den:.5%}, del bases {T[2] / den:.5%}, total {(T[8] + T[1] + T[2]) / den:.5%}")
    print(f"per aligned base (=+X): mismatch {T[8] / ab:.5%}, ins events {E[1] / ab:.5%}, del events {E[2] / ab:.5%}, ins bases {T[1] / ab:.5%}, del bases {T[2] / ab:.5%}; events total {(T[8] + E[1] + E[2]) / ab:.5%}")
    print("insertion length hist (10=10+)", sorted(IL.items())); print("deletion length hist (10=10+)", sorted(DL.items()))
    pr = [rate(r) for r in reads]; share = lambda t: sum(1 for x in pr if x > t) / len(pr)
    print("per-read event rate (X + I events + D events)/(=+X):", q(pr, f="{:.5f}"), "share with 0 events", sum(1 for x in pr if x == 0) / len(pr), "share >1%", share(0.01), "share >3%", share(0.03))
    for lab, sel in (("mapq60", lambda r: r["mapq"] == 60), ("mapq<60", lambda r: r["mapq"] < 60)):
        v = [rate(r) for r in reads if sel(r)]
        if v: print("  ", lab, len(v), q(v, (0.25, 0.5, 0.75, 0.9, 0.99), "{:.5f}"), "zero", sum(1 for x in v if x == 0) / len(v))
    print("de tag", q([r["de"] for r in reads if r["de"] is not None], f="{:.5f}"))
    return [("reads", n["prim"]), ("reads_sampled", len(reads)), ("read_length_p50", sorted(L)[min(len(L) - 1, int(0.5 * len(L)))]),
            ("events_per_aligned_base", f"{(T[8] + E[1] + E[2]) / ab:.6f}"), ("secondaries_per_read", f"{n['sec'] / n['prim']:.4f}"),
            ("share_without_secondaries", f"{z / n['prim']:.4f}")]

def load(path):
    """-> ({transcript: sorted exons}, {(chrom, strand): sorted (start, end, transcript)}, {(chrom, strand): longest span})"""
    tx = collections.defaultdict(list); meta = {}
    with open(path) as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.split("\t", 9)
            if len(f) < 9 or f[2] != "exon": continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    idx = collections.defaultdict(list)
    for t, e in tx.items():
        e.sort(); idx[meta[t]].append((e[0][0], e[-1][1], t))
    mx = {}
    for k in idx:
        idx[k].sort(); mx[k] = max(b - a for a, b, _ in idx[k])
    return tx, idx, mx

def ov(a, b):
    return sum(min(y, v) - max(x, u) for x, y in a for u, v in b if x < v and u < y)

def judge(strand, s, e, ex):
    """a read spanning s-e against a transcript's exons -> (5' complete, 3' complete, transcript bases upstream of the
    read's 5' end, downstream of its 3' end, exons missing at 5', transcript length, exons)"""
    below = lambda p: sum(max(0, min(b, p) - a) for a, b in ex); above = lambda p: sum(max(0, b - max(a, p)) for a, b in ex)
    tl = sum(b - a for a, b in ex)
    if strand == "+":
        return s < ex[0][1], e > ex[-1][0], below(s), above(e), sum(1 for a, b in ex if b <= s), tl, len(ex)
    return e > ex[-1][0], s < ex[0][1], above(e), below(s), sum(1 for a, b in ex if a >= e), tl, len(ex)

def study(lab, path, reads):
    """the printed end-completeness report of the reads against one annotation -> its counts"""
    t1 = time.time(); tx, idx, mx = load(path); print(f"== {lab}: transcripts {len(tx)} ({time.time() - t1:.0f}s)")
    res, opt, me, bylen, d5, d3, frac = collections.Counter(), collections.Counter(), collections.Counter(), collections.defaultdict(collections.Counter), [], [], []
    for r in reads:
        k = (r["chrom"], r["strand"])
        if k not in idx: res["no transcript"] += 1; continue
        v = idx[k]; s, e = r["blocks"][0][0], r["blocks"][-1][1]; i = bisect.bisect_left(v, (s - mx[k],)); cand = []
        while i < len(v) and v[i][0] < e:
            if v[i][1] > s:
                o = ov(r["blocks"], tx[v[i][2]])
                if o: cand.append((o, v[i][2]))
            i += 1
        if not cand: res["no transcript"] += 1; continue
        best = max(o for o, _ in cand)
        t = min((tt for o, tt in cand if o == best), key=lambda tt: sum(b - a for a, b in tx[tt]))    # strict: most exonic overlap, tie: the shortest
        has5, has3, up, dn, nm, tl, nx = judge(r["strand"], s, e, tx[t])
        res["5p reaches first exon" if has5 else "5p truncated"] += 1; res["3p reaches last exon" if has3 else "3p short"] += 1
        res["both ends" if has5 and has3 else "not both"] += 1; res["n"] += 1
        if not has5: me[min(nm, 10)] += 1
        d5.append(up); d3.append(dn); frac.append((tl - up - dn) / tl)
        b = "<1kb" if tl < 1000 else "1-2kb" if tl < 2000 else "2-3kb" if tl < 3000 else "3-5kb" if tl < 5000 else ">=5kb"
        bylen[b]["n"] += 1; bylen[b]["has5"] += has5
        opt["any"] += any(judge(r["strand"], s, e, tx[tt])[0] for o, tt in cand if o >= 0.9 * best); opt["n"] += 1
    n = res["n"]; print(" ", dict(res))
    print(f"  strict (best-overlap transcript): 5' in first exon {res['5p reaches first exon'] / n:.3f}, 5' truncated {res['5p truncated'] / n:.3f}; 3' in last exon {res['3p reaches last exon'] / n:.3f}; both {res['both ends'] / n:.3f}")
    print(f"  optimistic (any transcript with >= 0.9 of the best overlap): 5' in first exon {opt['any'] / opt['n']:.3f}")
    print("  exons missing at 5p among truncated (10=10+)", sorted(me.items()))
    print("  transcript bases upstream of read 5p end", q(d5, (0.25, 0.5, 0.75, 0.9, 0.95, 0.99)), "| downstream of 3p end", q(d3, (0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))
    print("  fraction of transcript length spanned", q(frac, (0.05, 0.1, 0.25, 0.5, 0.75, 0.9), "{:.3f}"))
    print("  5p complete by transcript length", {b: (c["n"], round(c["has5"] / c["n"], 3)) for b, c in sorted(bylen.items())})
    return res

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam"); ap.add_argument("ref_gtf"); ap.add_argument("assembly_gtf", nargs="?"); ap.add_argument("--out")
    ap.add_argument("--pct", type=int, default=14, help="percent of the primaries sampled by name hash (default 14: 136k reads stay under the cap)")
    a = ap.parse_args()
    t0 = time.time(); reads, n, secper = sample(a.bam, a.pct)
    print(f"records {n['rec']} primary {n['prim']} secondary {n['sec']} supplementary {n['sup']} unmapped {n['unm']}; sampled primaries {len(reads)}; {time.time() - t0:.0f}s")
    rows = report(reads, n, secper)
    res = study("reference annotation " + a.ref_gtf, a.ref_gtf, reads)
    if a.assembly_gtf:
        study("assembly " + a.assembly_gtf, a.assembly_gtf, reads)
    if a.out:
        rows += [("reads_on_a_transcript", res["n"])] + [(name, f"{res[k] / res['n']:.4f}") for name, k in
                 (("first_exon_5p_strict", "5p reaches first exon"), ("last_exon_3p", "3p reaches last exon"), ("both_ends_strict", "both ends"))]
        with open(a.out, "w") as fh:
            fh.write("stat\tvalue\n"); fh.writelines(f"{k}\t{v}\n" for k, v in rows)
    print(f"total {time.time() - t0:.0f}s")

if __name__ == "__main__":
    main()
