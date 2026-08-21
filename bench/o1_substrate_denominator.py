#!/usr/bin/env python3
"""O1: does the SUBSTRATE x DENOMINATOR combination capture the concept?

E_r's discriminating content is one number — coverage >= 0.50 of the SHORTER exon-sum — and it
provably does not order the classes: the accepted true pair GFPT1xGFPT2 scores 0.5353 while the
rejected false pair ATP1A1xATP4A scores 0.5689. A true pair below a false pair means NO threshold on
that statistic separates them. The census prescribed the shape of a fix: change the DENOMINATOR or
the SUBSTRATE, not the threshold.

THE 2x2 — both halves tested separately, never together:

                    coverage of SHORTER          coverage of LONGER
  exon-sum          the INCUMBENT                R2, refuted (makes E_r a function of which
                                                 REPRESENTATIVE you extract)
  genomic span      substrate fold, measured     ** NEVER TESTED **
                    better (P 0.916 vs 0.908)

Why the untested cell might work: R2 died of representative arbitrariness, and on exon-sums the
representative is a TRANSCRIPT MODEL (which isoform? how complete?). A genomic interval is set by the
locus, so that arbitrariness is largely gone — and it changes BOTH denominator and substrate.

PRE-STATED FAILURE MODE, checked before running: coverage-of-longer punishes pairs whose SPANS differ
rather than whose SEQUENCE does. On ANNOTATED spans that killed it (NPIP 10.6-49.4 kb around a ~16 kb
cassette; 134/171 true pairs could reach 0.50; NPIPB8-NPIPB2 capped at 0.215). On the READ-DERIVED
intervals used here the TP median span ratio is 1.35 vs FP 4.10, so it is 4x likelier to hit an FP —
but it genuinely caps 27/150 TP pairs and that is reported as a cost, not hidden.

UNIT = PAIR. GGO only. The FP arm is 14 pairs and is NOT held out; the TP side is load-bearing.
T8: offline, not the shipped binary.
"""
import csv, os, subprocess, sys, tempfile, statistics as st
from collections import defaultdict

D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
A   = "/home/juanfra/winloci_scratch/o1_antifp"
GEN = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
OUT = "/mnt/linuxdisk/home/juanfraitu/o1_subden"
TIERS = [(["-x", "asm20"], 0.80), (["-k", "11", "-w", "5"], 0.60)]
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(RC)[::-1]


def load_copies():
    return {f"{r['family_id']}~{r['copy_idx']}":
            (r["chrom"], int(r["start"]), int(r["end"]), r["strand"])
            for r in csv.DictReader(open(f"{D}/GGO_gwcat.copies.tsv"), delimiter="\t")}


def fetch_spans(names, cp):
    """genomic span per copy in TRANSCRIPTION orientation (matching the RNA-locus genomic-span path)"""
    os.makedirs(OUT, exist_ok=True)
    rf = os.path.join(OUT, "regions.txt")
    open(rf, "w").write("\n".join(f"{cp[n][0]}:{cp[n][1]+1}-{cp[n][2]}" for n in names) + "\n")
    p = subprocess.run(["samtools", "faidx", "-r", rf, GEN], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("samtools faidx failed: " + p.stderr[:300])
    seqs, key, buf = {}, None, []
    for line in p.stdout.splitlines():
        if line.startswith(">"):
            if key:
                seqs[key] = "".join(buf)
            key, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if key:
        seqs[key] = "".join(buf)
    out = {}
    for n in names:
        c, s, e, strand = cp[n]
        sq = seqs.get(f"{c}:{s+1}-{e}", "")
        out[n] = revcomp(sq) if strand == "-" else sq
    return out


def align(seqs, names):
    idx = {n: i for i, n in enumerate(names)}
    recs = defaultdict(list)
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "s.fa")
        with open(fa, "w") as fh:
            for n in names:
                fh.write(f">{idx[n]}\n{seqs[n]}\n")
        for flags, min_id in TIERS:
            r = subprocess.run(["minimap2", "-c", "-X", "--no-long-join", "-t", "4"] + flags + [fa, fa],
                               capture_output=True, text=True)
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 12:
                    continue
                q, t = int(f[0]), int(f[5])
                if q == t:
                    continue
                de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
                ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
                if ident < min_id:
                    continue
                recs[(min(q, t), max(q, t))].append(
                    (q, int(f[2]), int(f[3]), int(f[1]), t, int(f[7]), int(f[8]), int(f[6])))
    return recs, idx


def best_cov(records, longer):
    """M1-correct: the numerator's AXIS FOLLOWS THE DENOMINATOR."""
    best = 0.0
    for q, qs, qe, ql, t, ts, te, tl in records:
        pick_q = (ql >= tl) if longer else (ql <= tl)
        denom = max(ql, tl) if longer else min(ql, tl)
        span = (qe - qs) if pick_q else (te - ts)
        best = max(best, span / max(denom, 1))
    return best


def auc(pos, neg):
    if not pos or not neg:
        return float("nan")
    return sum((1.0 if p > n else 0.5 if p == n else 0.0) for p in pos for n in neg) / (len(pos) * len(neg))


def main():
    cp = load_copies()
    fp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/fp14_detail.tsv"), delimiter="\t")]
    tp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/tp150_detail.tsv"), delimiter="\t")]
    pairs = [("FP", a, b) for a, b in fp if a in cp and b in cp] + \
            [("TP", a, b) for a, b in tp if a in cp and b in cp]
    names = sorted({n for _, a, b in pairs for n in (a, b)})
    print(f"pairs FP={sum(1 for k,_,_ in pairs if k=='FP')} TP={sum(1 for k,_,_ in pairs if k=='TP')}"
          f"  distinct copies={len(names)}", flush=True)

    span = fetch_spans(names, cp)
    lens = sorted(len(v) for v in span.values())
    print(f"genomic spans: {sum(lens)/1e6:.1f} Mb, median {lens[len(lens)//2]} bp", flush=True)
    recs, idx = align(span, names)
    print(f"span all-vs-all: {sum(len(v) for v in recs.values())} qualifying records\n", flush=True)

    rows = []
    for klass, a, b in pairs:
        r = recs.get((min(idx[a], idx[b]), max(idx[a], idx[b])), [])
        rows.append({"klass": klass, "a": a, "b": b,
                     "span_short": round(best_cov(r, False), 4),
                     "span_long":  round(best_cov(r, True), 4),
                     "span_ratio": round(max(len(span[a]), len(span[b]))
                                         / max(min(len(span[a]), len(span[b])), 1), 3)})
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/subden.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        [w.writerow(r) for r in rows]

    print("=== SEPARATION (unit = pair; AUC = P(FP scores higher than TP); 0.5 = no signal) ===")
    print(f"  {'cell':<40} {'FP med':>8} {'TP med':>8} {'AUC':>8}")
    for lab, key in (("genomic span x cov-of-SHORTER", "span_short"),
                     ("genomic span x cov-of-LONGER  **untested**", "span_long")):
        P = [r[key] for r in rows if r["klass"] == "FP"]
        N = [r[key] for r in rows if r["klass"] == "TP"]
        print(f"  {lab:<40} {st.median(P):>8.4f} {st.median(N):>8.4f} {auc(P,N):>8.4f}")

    print("\n=== OPERATING POINT: reject when cov-of-longer on the SPAN falls below c ===")
    print(f"  {'c':>6} {'FP rejected':>14} {'TP lost':>12} {'TP cost':>9}")
    P = [r["span_long"] for r in rows if r["klass"] == "FP"]
    N = [r["span_long"] for r in rows if r["klass"] == "TP"]
    for c in (0.10, 0.15, 0.20, 0.30, 0.40, 0.50):
        fr = sum(1 for x in P if x < c); tl = sum(1 for x in N if x < c)
        print(f"  {c:>6.2f} {fr:>7}/{len(P):<5} {tl:>5}/{len(N):<5} {tl/len(N):>8.4f}")


if __name__ == "__main__":
    main()
