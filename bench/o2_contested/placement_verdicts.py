#!/usr/bin/env python3
"""r1701 (docs/PREREG_placement_verdicts_2026-11-13.md): one verdict per placement of a near-tied molecule.
Each secondary with AS >= RATIO x the molecule's best is judged against the primary with s1_o2.certify_pair (read
positions both records align as = / X where the two reference bases differ): excluded = the primary is certified,
over_primary = the secondary is, identical = no differing column, undecided = the rest and mosaic reads.
shortcut: prototype certificate (no junction term, base quality or origin test); port into copy_assign on a panel PASS.
  verdicts BAM FASTA MOLECULES.tsv RATIO > verdicts.tsv     (BAM read in file order; needs samtools on PATH)
  truth VERDICTS COPIES.tsv DUMP [DUMP...]                   step 0: verdicts against the strict_truth dumps
  filter VERDICTS IN.bam OUT.bam [RATIO]                     drop excluded secondaries (with RATIO: only those >= RATIO)"""
import collections, csv, os, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tool_difference"))
from s1_o2 import COMP, aligned_columns, certify_pair, parse_cigar  # noqa: E402

HDR = "read\tchrom\tpos\tend\tstrand\tprimary\tas\tas_ratio\tms\tms_ratio\tspliced\tka\tkb\tm\tp\tverdict"

def verdict(v, ka, kb, m):
    if m == 0: return "identical"
    if ka >= 2 and kb >= 2: return "undecided"  # mosaic: fits neither reference cleanly
    return {"A": "excluded", "B": "over_primary"}.get(v, "undecided")

def tag(ln, t):
    i = ln.find("\t" + t + ":i:")
    return int(ln[i + 6:].split("\t", 1)[0]) if i >= 0 else None

def verdicts(bam, fasta, molecules, ratio, out=sys.stdout):
    import pysam
    ratio = float(ratio); best = {}
    for ln in open(molecules):
        if ln[0] == "#": continue
        f = ln.split("\t", 8)
        if int(f[3]) >= 2 and int(f[2]) >= ratio * int(f[1]): best[f[0]] = (int(f[1]), int(f[6]))
    recs = collections.defaultdict(list)
    p = subprocess.Popen(["samtools", "view", "-F", "2052", bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for ln in p.stdout:
        f = ln.split("\t", 6)
        b = best.get(f[0])
        if b is None: continue
        a = tag(ln, "AS"); flag = int(f[1])
        if flag & 256 and (a is None or a < ratio * b[0]): continue
        recs[f[0]].append((flag, f[2], int(f[3]) - 1, f[5], a, tag(ln, "ms")))
    if p.wait(): sys.exit("samtools view failed")
    fa = pysam.FastaFile(fasta); print(HDR, file=out)

    def cols_of(r):
        cig = parse_cigar(r[3]); rev = bool(r[0] & 16); c = aligned_columns(r[2], cig, rev)
        end = r[2] + sum(l for op, l in cig if op in (0, 2, 3, 7, 8))
        if not c: return c, None, end, cig
        lo = min(x for x, _ in c.values()); hi = max(x for x, _ in c.values()) + 1; seq = fa.fetch(r[1], lo, hi).upper()
        return c, ((lambda o: seq[c[o][0] - lo].translate(COMP)) if rev else (lambda o: seq[c[o][0] - lo])), end, cig

    for name, rs in recs.items():
        prim = [r for r in rs if not r[0] & 256]; sec = [r for r in rs if r[0] & 256]
        if not prim or not sec: continue
        bas, bms = best[name]; ca, ba, pend, pcig = cols_of(prim[0])
        row = lambda r, end, cig, *v: print(name, r[1], r[2] + 1, end, "-" if r[0] & 16 else "+", int(not r[0] & 256), r[4], f"{r[4] / bas:.4f}" if r[4] is not None and bas else "NA",
                                           r[5] if r[5] is not None else "NA", f"{r[5] / bms:.4f}" if r[5] is not None and bms > 0 else "NA", int(any(op == 3 for op, _ in cig)), *v, sep="\t", file=out)
        row(prim[0], pend, pcig, "NA", "NA", "NA", "NA", "primary")
        for r in sec:
            cb, bb, end, cig = cols_of(r)
            if not ca or not cb: row(r, end, cig, 0, 0, 0, "1", "undecided"); continue
            v, ka, kb, m, pv = certify_pair(ca, cb, ba, bb); row(r, end, cig, ka, kb, m, f"{pv:.3g}", verdict(v, ka, kb, m))

def truth(vfile, copies, *dumps):
    span = {(r["family_id"], r["copy_idx"]): (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(copies), delimiter="\t")}
    T = {}
    for d in dumps:
        for l in open(d):
            q, _, _, _, _, _, t = l.rstrip("\n").split("\t")
            if t != "OUTSIDE": T.setdefault(q, set()).update(span[tuple(x.split(":"))] for x in t.split(";"))
    c = collections.Counter(); shares = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(vfile), delimiter="\t"):
        if r["primary"] == "1":
            if r["read"] in T: c["decidable molecules with a verdict row"] += 1
            continue
        a = float(r["as_ratio"]); band = ">= 0.98" if a >= 0.98 else "0.9-0.98" if a >= 0.9 else "0.8-0.9"
        shares[(band, "spliced" if r["spliced"] == "1" else "unspliced")][r["verdict"]] += 1
        if r["read"] not in T: continue
        on = any(ch == r["chrom"] and s < int(r["end"]) and int(r["pos"]) - 1 < e for ch, s, e in T[r["read"]])
        c[f"{r['verdict']}: {'ON the truth copy' if on else 'off the truth copy'}"] += 1
        c[f"{r['verdict']} ({band}): {'ON the truth copy' if on else 'off the truth copy'}"] += 1
    for k in sorted(c): print(k, c[k], sep="\t")
    ex = c["excluded: ON the truth copy"] + c["excluded: off the truth copy"]
    print(f"S1 excluded placements of decidable molecules {ex} (bar >= 100); S2 on the truth copy {c['excluded: ON the truth copy']} = {c['excluded: ON the truth copy'] / max(ex, 1):.4f} (bar <= 0.02)")
    print("STEP0", "PASS" if ex >= 100 and c["excluded: ON the truth copy"] <= 0.02 * ex else "STOP")
    print("all secondaries, band / shape\texcluded\tover_primary\tidentical\tundecided")
    for k in sorted(shares): print(*k, *(shares[k][v] for v in ("excluded", "over_primary", "identical", "undecided")), sep="\t")

def filt(vfile, inbam, outbam, ratio="0"):
    drop = {(r["read"], r["chrom"], r["pos"], r["strand"]) for r in csv.DictReader(open(vfile), delimiter="\t") if r["verdict"] == "excluded" and float(r["as_ratio"]) >= float(ratio)}
    i = subprocess.Popen(["samtools", "view", "-@", "8", "-h", inbam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    o = subprocess.Popen(["samtools", "view", "-@", "24", "-b", "-o", outbam, "-"], stdin=subprocess.PIPE, text=True, bufsize=1 << 20); n = 0
    for ln in i.stdout:
        if ln[0] != "@":
            f = ln.split("\t", 4); fl = int(f[1])
            if fl & 256 and (f[0], f[2], f[3], "-" if fl & 16 else "+") in drop: n += 1; continue
        o.stdin.write(ln)
    o.stdin.close()
    if i.wait() or o.wait(): sys.exit("samtools failed")
    print(f"excluded records listed {len(drop)}, removed {n}", file=sys.stderr)

def _selfcheck():
    assert verdict("A", 5, 0, 5) == "excluded" and verdict("B", 0, 5, 5) == "over_primary" and verdict(None, 0, 0, 0) == "identical"
    assert verdict("A", 6, 2, 8) == "undecided" and verdict(None, 1, 0, 1) == "undecided"
    assert tag("r\t0\tc\t1\tAS:i:7\tms:i:9\n", "ms") == 9 and tag("r\t0\n", "AS") is None

if __name__ == "__main__":
    _selfcheck(); {"verdicts": verdicts, "truth": truth, "filter": filt}[sys.argv[1]](*sys.argv[2:])
