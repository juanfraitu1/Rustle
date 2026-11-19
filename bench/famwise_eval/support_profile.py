#!/usr/bin/env python3
"""What do whole reads really support? A deterministic support profile per transcript (an assembly's GTF, or truth intervals),
from the reads themselves (user, 2026-11-19: make sure the reads really support every important structure, and that reads that
deviate even slightly from an interval do not introduce errors; the long-read counterpart of IRFinder-S's long-read mode, no
learned model).

  support_profile.py --gtf GTF --bam BAM --out PREFIX [--jitter 3] [--near 20] [--end-window 50] [--min-cover-reads 3]
                     [--all-records] [--multi-exon-only]

Reads: primary records (or every record with --all-records) on the transcript's strand overlapping its span.
A read is COMPATIBLE with the transcript when every intron of the read EQUALS one of the transcript's introns, its introns are
consecutive in the transcript, and none of its aligned blocks covers more than --jitter bp inside a transcript intron. Exact
matching is deliberate: a read whose junction lies 1 to --jitter bp off supports a neighbouring structure, and counting it here
let true reads vouch for 1-3 bp shifted false transcripts on SIRV (2026-11-19). Such reads are SHIFTED (counted per junction);
other incompatible reads get ONE reason, the first that applies:
  near_miss   a read intron whose ends lie within --near bp of a transcript intron but not within --jitter (a slightly shifted junction)
  skip        a read intron that joins non-consecutive transcript exons (a transcript intron's donor to a later acceptor)
  other_junction  a read intron unrelated to the transcript's introns
  retention   a read block covers more than --jitter bp inside a transcript intron
Per junction: reads with it exact, within --jitter (not exact), near-miss reads at it. Per transcript: compatible reads, reads with
the full chain, 5'-complete and 3'-complete compatible reads (their 5' / 3' end within --end-window bp of the transcript's), and
the share of exonic bp covered by >= 1 and >= --min-cover-reads compatible reads.
Flags (a transcript can carry several): WEAK_JUNCTION (a junction with 0 exact reads), SHIFTED_DOMINANT (a junction where reads
shifted 1 to --jitter bp outnumber exact ones), THIN_EXACT (the weakest junction has fewer exact reads than --thin-share of the
strongest junction, or than --min-exact; relative by default because an absolute floor measures library depth),
NEAR_MISS_DOMINANT (a junction where near-miss reads outnumber exact ones), UNCOVERED (< 90% of exonic bp covered by >= 1 compatible read),
NO_FULL_READ (no compatible read carries every junction). A near-miss read whose own junction is ESTABLISHED (>= 2 exact reads and
>= 10% of the nearby transcript junction's support) is counted as `alt_site` (another real splice site), not as a near-miss;
retention reads are counted, not flagged (real retention isoforms exist).
Writes PREFIX.transcripts.tsv, PREFIX.junctions.tsv, PREFIX.summary.txt."""
import argparse, bisect, collections, os, re, sys
import pysam

def read_gtf(path):
    ex, at = collections.defaultdict(list), {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[t].append((int(f[3]) - 1, int(f[4]))); at[t] = (f[0], f[6])
    return {t: dict(chrom=at[t][0], strand=at[t][1], exons=sorted(e)) for t, e in ex.items()}

def blocks_introns(rec):
    blocks, introns, pos, s = [], [], rec.reference_start, rec.reference_start
    for op, n in rec.cigartuples or ():
        if op in (0, 2, 7, 8):
            pos += n
        elif op == 3:
            blocks.append((s, pos)); introns.append((pos, pos + n)); pos += n; s = pos
    blocks.append((s, pos))
    return blocks, introns

def classify(blocks, introns, tint, jitter, near, established=frozenset()):
    """-> (verdict, list of (transcript intron index, 'exact'|'jitter') for the read's introns, near-miss junction indexes)
    verdict: 'compatible' or the reason of incompatibility"""
    hits, nm = [], []
    for d, a in introns:
        best = None
        for k, (td, ta) in enumerate(tint):
            dd, da = abs(d - td), abs(a - ta)
            if dd == 0 and da == 0:
                best = (k, "exact"); break
            if dd <= jitter and da <= jitter and best is None:
                best = (k, "jitter")
        if best:
            hits.append(best); continue
        near_k = [k for k, (td, ta) in enumerate(tint) if abs(d - td) <= near and abs(a - ta) <= near]
        if near_k:
            if (d, a) in established:                         # a well-supported junction of its own: another real splice site
                return "alt_site", hits, nm
            nm.append(near_k[0]); return "near_miss", hits, nm
        donor = [k for k, (td, _) in enumerate(tint) if abs(d - td) <= jitter]
        acc = [k for k, (_, ta) in enumerate(tint) if abs(a - ta) <= jitter]
        if donor and acc and acc[0] > donor[0]:
            return "skip", hits, nm
        return "other_junction", hits, nm
    ks = [k for k, _ in hits]
    if ks and ks != list(range(ks[0], ks[0] + len(ks))):
        return "skip", hits, nm
    for s, e in blocks:
        for td, ta in tint:
            inside = min(e, ta) - max(s, td)
            if inside > jitter:
                return "retention", hits, nm
    return "compatible", hits, nm

def profile(t, reads, jitter, near, end_window, min_cov, established=frozenset(), min_exact=0, thin_share=0.1):
    ex = t["exons"]; tint = [(a[1], b[0]) for a, b in zip(ex, ex[1:])]
    nj = len(tint)
    j_exact, j_jit, j_near = [0] * nj, [0] * nj, [0] * nj
    why = collections.Counter()
    covd = collections.Counter()
    full = five = three = comp = 0
    tlo, thi = ex[0][0], ex[-1][1]
    for blocks, introns in reads:
        v, hits, nm = classify(blocks, introns, tint, jitter, near, established)
        for k in nm:
            j_near[k] += 1
        if v != "compatible":
            why[v] += 1; continue
        if any(kind == "jitter" for _, kind in hits):
            # a junction 1-3 bp off the transcript's: the read supports a NEIGHBOURING structure, not this one; it vouches for
            # nothing here (SIRV, 2026-11-19: false transcripts are 1-3 bp shifts of true ones, carried by the true reads only
            # under a tolerance)
            why["shifted"] += 1
            for k, kind in hits:
                if kind == "jitter":
                    j_jit[k] += 1
            continue
        comp += 1
        for k, kind in hits:
            j_exact[k] += 1
        if nj and len(hits) == nj:
            full += 1
        rs, re_ = blocks[0][0], blocks[-1][1]
        left_ok, right_ok = abs(rs - tlo) <= end_window, abs(re_ - thi) <= end_window
        if t["strand"] == "-":
            left_ok, right_ok = right_ok, left_ok
        five += left_ok; three += right_ok
        for s, e in blocks:
            for a, b in ex:
                lo, hi = max(s, a), min(e, b)
                for p in range(lo, hi):
                    covd[p] += 1
    exlen = sum(b - a for a, b in ex)
    c1 = sum(1 for p, n in covd.items() if n >= 1) / exlen if exlen else 0
    ck = sum(1 for p, n in covd.items() if n >= min_cov) / exlen if exlen else 0
    flags = []
    if nj and min(j_exact) == 0:
        flags.append("WEAK_JUNCTION")
    if nj and any(j_jit[k] > j_exact[k] for k in range(nj)):
        flags.append("SHIFTED_DOMINANT")
    if nj and min(j_exact) > 0 and (min(j_exact) < thin_share * max(j_exact) or min(j_exact) < min_exact):
        flags.append("THIN_EXACT")
    if any(j_near[k] > j_exact[k] for k in range(nj)):
        flags.append("NEAR_MISS_DOMINANT")
    if c1 < 0.9:
        flags.append("UNCOVERED")
    if nj and full == 0:
        flags.append("NO_FULL_READ")
    return dict(reads=len(reads), compatible=comp, full_chain=full, five_complete=five, three_complete=three, cover1=c1, coverk=ck,
                why=why, j_exact=j_exact, j_jit=j_jit, j_near=j_near, flags=flags)

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True); ap.add_argument("--bam", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--jitter", type=int, default=3); ap.add_argument("--near", type=int, default=20)
    ap.add_argument("--end-window", type=int, default=50); ap.add_argument("--min-cover-reads", type=int, default=3)
    ap.add_argument("--thin-share", type=float, default=0.1, help="THIN_EXACT when the weakest junction has fewer exact reads than this share of the strongest")
    ap.add_argument("--min-exact", type=int, default=0, help="also THIN_EXACT below this absolute number of exact reads (off; 10 suits deep libraries such as SIRV)")
    ap.add_argument("--all-records", action="store_true"); ap.add_argument("--multi-exon-only", action="store_true")
    a = ap.parse_args(argv)
    T = read_gtf(a.gtf)
    if a.multi_exon_only:
        T = {k: v for k, v in T.items() if len(v["exons"]) > 1}
    bam = pysam.AlignmentFile(a.bam)
    rows, jrows, fl = [], [], collections.Counter()
    for tid, t in sorted(T.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["exons"][0][0])):
        if t["chrom"] not in bam.references:
            continue
        reads = []
        for r in bam.fetch(t["chrom"], t["exons"][0][0], t["exons"][-1][1]):
            if r.is_unmapped or r.is_supplementary or (r.is_secondary and not a.all_records):
                continue
            if ("-" if r.is_reverse else "+") != t["strand"]:
                continue
            reads.append(blocks_introns(r))
        jc = collections.Counter(j for _, ins in reads for j in ins)
        ex_ = t["exons"]; tint_ = [(x[1], y[0]) for x, y in zip(ex_, ex_[1:])]
        tsup = {j: jc[j] for j in tint_}
        # a read junction is ESTABLISHED (a real alternative site, not an artifact) when >= 2 reads carry it exactly and it has
        # >= 10% of the support of the transcript junction it lies near
        established = frozenset(j for j, n in jc.items() if n >= 2 and j not in tsup and
                                all(n >= 0.1 * tsup[k] for k in tint_ if abs(j[0] - k[0]) <= a.near and abs(j[1] - k[1]) <= a.near))
        p = profile(t, reads, a.jitter, a.near, a.end_window, a.min_cover_reads, established, a.min_exact, a.thin_share)
        for f in p["flags"]:
            fl[f] += 1
        fl["any flag"] += bool(p["flags"]); fl["transcripts"] += 1
        rows.append((tid, t, p))
        ex = t["exons"]
        for k, (d, ac) in enumerate((a_[1], b[0]) for a_, b in zip(ex, ex[1:])):
            jrows.append((tid, t["chrom"], d, ac, t["strand"], p["j_exact"][k], p["j_jit"][k], p["j_near"][k]))
    with open(a.out + ".transcripts.tsv", "w") as fh:
        fh.write("transcript\tchrom\tstrand\tstart\tend\texons\treads\tcompatible\tfull_chain\tfive_complete\tthree_complete\t"
                 "cover_ge1\tcover_ge_k\tmin_junction_exact\tshifted\tnear_miss\talt_site\tskip\tother_junction\tretention\tflags\n")
        for tid, t, p in rows:
            nj = len(t["exons"]) - 1
            mj = min((p["j_exact"][k] for k in range(nj)), default="")
            fh.write(f"{tid}\t{t['chrom']}\t{t['strand']}\t{t['exons'][0][0]}\t{t['exons'][-1][1]}\t{len(t['exons'])}\t{p['reads']}\t"
                     f"{p['compatible']}\t{p['full_chain']}\t{p['five_complete']}\t{p['three_complete']}\t{p['cover1']:.3f}\t{p['coverk']:.3f}\t"
                     f"{mj}\t{p['why']['shifted']}\t{p['why']['near_miss']}\t{p['why']['alt_site']}\t{p['why']['skip']}\t{p['why']['other_junction']}\t{p['why']['retention']}\t"
                     f"{','.join(p['flags']) or '-'}\n")
    with open(a.out + ".junctions.tsv", "w") as fh:
        fh.write("transcript\tchrom\tdonor\tacceptor\tstrand\texact\tshifted_1_to_jitter\tnear_miss\n")
        fh.writelines("\t".join(map(str, j)) + "\n" for j in jrows)
    n = fl["transcripts"]
    L = [f"transcripts {n}; with any flag {fl['any flag']} ({fl['any flag'] / max(1, n):.1%})"]
    L += [f"  {k}: {fl[k]} ({fl[k] / max(1, n):.1%})" for k in ("WEAK_JUNCTION", "SHIFTED_DOMINANT", "THIN_EXACT", "NEAR_MISS_DOMINANT", "UNCOVERED", "NO_FULL_READ")]
    why = collections.Counter()
    for _, _, p in rows:
        why.update(p["why"])
    tot_reads = sum(p["reads"] for _, _, p in rows); tot_comp = sum(p["compatible"] for _, _, p in rows)
    L.append(f"read-transcript pairs {tot_reads}: compatible {tot_comp} ({tot_comp / max(1, tot_reads):.1%}); incompatible by reason: "
             + ", ".join(f"{k} {why[k]}" for k in ("shifted", "near_miss", "alt_site", "skip", "other_junction", "retention")))
    open(a.out + ".summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()
