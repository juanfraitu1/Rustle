#!/usr/bin/env python3
"""O1: is the passing block's divergence CONSISTENT with the rest of what the pair shares?

WHY A NEW KIND OF STATISTIC IS NEEDED. Everything tested so far is either an AMOUNT statistic
(coverage on any substrate/denominator, aligned bp, clipping — the whole 2x2 closed 2026-08-20) or a
MODEL-STRUCTURE statistic (exon count, junction crossing — all carry a 100x exon-count bias, because
structure correlates with exon count by construction).

THE IDEA. A duplication copies a unit at one moment, so every piece of sequence the two loci share
diverged together and should show ONE divergence level. A mobile element inserted into both loci
later is MORE SIMILAR than the pair's own ancestry. So the discriminating question is not "how much
aligns" but "is the passing block's identity consistent with the pair's OTHER shared sequence?"

  contrast = identity(passing block) - median identity(other shared blocks)

  true paralogue pair  -> contrast ~ 0   (one duplication, one divergence level)
  repeat-bridged pair  -> contrast > 0   (the bridge is younger than anything real they share)

Properties: pair-local (a, b only) so P1 is untouched; scale-free (a difference of rates, not a
fraction of a length); independent of exon count; and it needs no threshold on the amount of anything.

⚠ It ABSTAINS when a pair shares nothing but the passing block — which is exactly the FP shape, so
abstention here is informative and is reported, never silently dropped.

UNIT = PAIR. GGO. FP arm 14 (NOT held out), TP arm 150 (load-bearing). T8: offline.
"""
import csv, os, subprocess, sys, tempfile, statistics as st
from collections import defaultdict

D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
A   = "/home/juanfra/winloci_scratch/o1_antifp"
GEN = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
OUT = "/mnt/linuxdisk/home/juanfraitu/o1_subden"
MIN_BLOCK = 200          # a "shared block" must be at least this long to carry a divergence estimate
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o1_substrate_denominator import load_copies, fetch_spans


def align_all(seqs, names):
    """every record >= MIN_BLOCK, both tiers, no identity floor — we need the DISTRIBUTION"""
    idx = {n: i for i, n in enumerate(names)}
    recs = defaultdict(list)
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "s.fa")
        with open(fa, "w") as fh:
            for n in names:
                fh.write(f">{idx[n]}\n{seqs[n]}\n")
        for flags in (["-x", "asm20"], ["-k", "11", "-w", "5"]):
            r = subprocess.run(["minimap2", "-c", "-X", "--no-long-join", "-t", "4"] + flags + [fa, fa],
                               capture_output=True, text=True)
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 12:
                    continue
                q, t = int(f[0]), int(f[5])
                if q == t:
                    continue
                qs, qe = int(f[2]), int(f[3])
                if qe - qs < MIN_BLOCK:
                    continue
                de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
                ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
                recs[(min(q, t), max(q, t))].append((ident, qe - qs, qs, qe))
    return recs, idx


def auc(pos, neg):
    if not pos or not neg:
        return float("nan")
    return sum((1.0 if p > n else 0.5 if p == n else 0.0) for p in pos for n in neg) / (len(pos)*len(neg))


def main():
    cp = load_copies()
    fp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/fp14_detail.tsv"), delimiter="\t")]
    tp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/tp150_detail.tsv"), delimiter="\t")]
    pairs = [("FP", a, b) for a, b in fp if a in cp and b in cp] + \
            [("TP", a, b) for a, b in tp if a in cp and b in cp]
    names = sorted({n for _, a, b in pairs for n in (a, b)})
    span = fetch_spans(names, cp)
    recs, idx = align_all(span, names)
    print(f"pairs FP={sum(1 for k,_,_ in pairs if k=='FP')} TP={sum(1 for k,_,_ in pairs if k=='TP')}"
          f"   blocks >= {MIN_BLOCK} bp: {sum(len(v) for v in recs.values())}\n", flush=True)

    rows = []
    for klass, a, b in pairs:
        r = sorted(recs.get((min(idx[a], idx[b]), max(idx[a], idx[b])), []), key=lambda x: -x[1])
        if not r:
            rows.append({"klass": klass, "a": a, "b": b, "n_blocks": 0,
                         "top_ident": "NA", "other_med": "NA", "contrast": "NA"}); continue
        top = r[0]
        others = [x[0] for x in r[1:]]
        rows.append({"klass": klass, "a": a, "b": b, "n_blocks": len(r),
                     "top_ident": round(top[0], 4),
                     "other_med": round(st.median(others), 4) if others else "NA",
                     "contrast": round(top[0] - st.median(others), 4) if others else "NA"})
    with open(f"{OUT}/blockcons.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader()
        [w.writerow(x) for x in rows]

    for lab in ("FP", "TP"):
        R = [r for r in rows if r["klass"] == lab]
        nb = [r["n_blocks"] for r in R]
        ab = sum(1 for r in R if r["contrast"] == "NA")
        print(f"  {lab}: n={len(R)}  median blocks={st.median(nb):.0f}  "
              f"ABSTAINS (0 or 1 shared block) {ab}/{len(R)} = {ab/len(R):.4f}")
    print()
    P = [r["contrast"] for r in rows if r["klass"] == "FP" and r["contrast"] != "NA"]
    N = [r["contrast"] for r in rows if r["klass"] == "TP" and r["contrast"] != "NA"]
    print("=== CONTRAST = identity(top block) - median identity(other shared blocks) ===")
    print(f"  FP scored n={len(P)}  median {st.median(P):+.4f}")
    print(f"  TP scored n={len(N)}  median {st.median(N):+.4f}")
    print(f"  AUC(FP > TP) = {auc(P,N):.4f}   (0.5 = no signal)")
    print("\n=== operating point: reject when contrast >= x ===")
    print(f"  {'x':>7} {'FP rejected':>14} {'TP lost':>12} {'TP cost':>9}")
    for x in (0.00, 0.02, 0.05, 0.10, 0.15):
        fr = sum(1 for v in P if v >= x); tl = sum(1 for v in N if v >= x)
        print(f"  {x:>7.2f} {fr:>7}/{len(P):<5} {tl:>5}/{len(N):<5} {tl/max(len(N),1):>8.4f}")


if __name__ == "__main__":
    main()
