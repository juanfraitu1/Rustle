#!/usr/bin/env python3
"""Is the partner the locus's BEST match in the genome, or one of hundreds?

A NEW EXTERNAL REFERENCE. Everything internal to the pair is length in disguise; the external
references tried so far ask about the SHARED SEQUENCE (genome-anchored repeat multiplicity), about a
DUPLICATION CATALOG (SD containment) or about the NEIGHBOURHOOD (flank homology). This asks a
different question: **where does the partner rank among all the places in the genome this locus looks
like?**

  real paralogue pair -> the partner is at or near the top; there are few competitors
  repeat bridge       -> the partner is one of hundreds of equally good repeat hits

⚠ T19 — "argmax/best-hit counting is not counting", and best-hit has inflated 1.7x elsewhere here. So
this reports the RANK DISTRIBUTION and the competitor COUNT, never a bare best-hit call.

⚠ Known weakness of any reciprocity test: in a large family the top hit may be a DIFFERENT member than
the partner, so a strict best-hit rule punishes big families. Rank handles that; a binary would not.

UNIT = PAIR. GGO. FP 14 (not held out), TP 150 (load-bearing). T8: offline.
"""
import csv, os, subprocess, sys, collections, statistics as st

D   = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
A   = "/home/juanfra/winloci_scratch/o1_antifp"
GEN = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
OUT = "/mnt/linuxdisk/home/juanfraitu/o1_rank"
MIN_ID, MIN_BP = 0.80, 300


def main():
    os.makedirs(OUT, exist_ok=True)
    cp = {f"{r['family_id']}~{r['copy_idx']}": (r["chrom"], int(r["start"]), int(r["end"]))
          for r in csv.DictReader(open(f"{D}/GGO_gwcat.copies.tsv"), delimiter="\t")}
    fp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/fp14_detail.tsv"), delimiter="\t")]
    tp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/tp150_detail.tsv"), delimiter="\t")]
    pairs = [("FP", a, b) for a, b in fp if a in cp and b in cp] + \
            [("TP", a, b) for a, b in tp if a in cp and b in cp]
    names = sorted({n for _, a, b in pairs for n in (a, b)})

    # rep sequences for the arms only
    seqs, key, buf = {}, None, []
    for line in open(f"{D}/GGO_gwcat.copies.fa"):
        if line.startswith(">"):
            if key:
                seqs[key] = "".join(buf)
            p = line[1:].strip().split("|")
            key, buf = f"{p[0]}~{p[1]}", []
        else:
            buf.append(line.strip())
    seqs[key] = "".join(buf)

    fa = f"{OUT}/arms.fa"
    with open(fa, "w") as fh:
        for n in names:
            fh.write(f">{n}\n{seqs[n]}\n")
    paf = f"{OUT}/arms_vs_genome.paf"
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        print("aligning arm reps to the genome (one index build) ...", flush=True)
        r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--secondary=yes", "-N", "500",
                            "-p", "0.05", "-t", "4", GEN, fa],
                           capture_output=True, text=True)
        open(paf, "w").write(r.stdout)
    print(f"genomic hits: {sum(1 for _ in open(paf))} records", flush=True)

    hits = collections.defaultdict(list)
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q = f[0]
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
        span = int(f[3]) - int(f[2])
        if ident < MIN_ID or span < MIN_BP:
            continue
        c, ts, te = f[5], int(f[7]), int(f[8])
        cq, sq, eq = cp[q]
        if c == cq and ts < eq and sq < te:          # its own locus
            continue
        hits[q].append((span * ident, c, ts, te))    # score-ish
    for q in hits:
        hits[q].sort(key=lambda x: -x[0])

    def rank_of(a, b):
        cb, sb, eb = cp[b]
        for i, (_, c, ts, te) in enumerate(hits.get(a, []), 1):
            if c == cb and ts < eb and sb < te:
                return i
        return None

    rows = []
    for k, a, b in pairs:
        ra, rb = rank_of(a, b), rank_of(b, a)
        rows.append((k, a, b, ra, rb, len(hits.get(a, [])), len(hits.get(b, []))))

    print("\n=== where does the partner rank among the locus's genomic hits? ===")
    print(f"  {'arm':<5} {'partner found both ways':>24} {'median best rank':>18} {'median competitors':>20}")
    for lab in ("FP", "TP"):
        R = [r for r in rows if r[0] == lab]
        both = [r for r in R if r[3] and r[4]]
        best = [min(r[3], r[4]) for r in both]
        comp = [max(r[5], r[6]) for r in R]
        print(f"  {lab:<5} {len(both):>10}/{len(R):<6} {'':>6} "
              f"{(st.median(best) if best else float('nan')):>18.1f} {st.median(comp):>20.0f}")
    print("\n=== as a rule: accept only if the partner is within the top-k genomic hits ===")
    print(f"  {'k':>4} {'FP admitted':>14} {'TP kept':>12} {'TP rate':>9}")
    for k in (1, 2, 3, 5, 10):
        F = [r for r in rows if r[0] == "FP"]; T = [r for r in rows if r[0] == "TP"]
        fa_ = sum(1 for r in F if (r[3] and r[3] <= k) or (r[4] and r[4] <= k))
        tk = sum(1 for r in T if (r[3] and r[3] <= k) or (r[4] and r[4] <= k))
        print(f"  {k:>4} {fa_:>7}/{len(F):<6} {tk:>5}/{len(T):<6} {tk/len(T):>9.4f}")


if __name__ == "__main__":
    main()
