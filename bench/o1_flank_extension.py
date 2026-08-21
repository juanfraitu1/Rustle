#!/usr/bin/env python3
"""Does the homology EXTEND BEYOND the gene into the flanks?

THE HYPOTHESIS. A segmental duplication copies a genomic NEIGHBOURHOOD, so a true duplicated pair
should share sequence OUTSIDE the gene body as well as inside it. A mobile element inserted into two
unrelated genes shares only the element — its flanks are two different genomic neighbourhoods and do
not align.

WHY IT MATTERS EVEN THOUGH SD CONTAINMENT ALREADY WORKS. SD containment rejects 14/14 FPs but needs a
species-specific SD catalog (the gorilla one is a multi-day SEDEF run), which is exactly the
portability objection that stopped it becoming a membership condition. This is the SAME IDEA computed
LOCALLY from the genome alone: pair-local, no catalog, no third input, and it works for any species
with an assembly.

TEST. Take FLANK bp immediately outside each locus on both sides, exclude the gene bodies entirely,
and ask whether the flanks align to each other.

  true duplication  -> flanks align (the duplication is bigger than the gene)
  repeat bridge     -> flanks are unrelated neighbourhoods, no alignment

⚠ Directional caveat: a duplication whose boundaries fall INSIDE the gene, or an old one whose flanks
have diverged, will look like a repeat bridge here. So a negative is weak; a POSITIVE is the
informative direction. Reported both ways.

UNIT = PAIR. GGO. FP 14 (not held out), TP 150 (load-bearing). T8: offline.
"""
import csv, os, subprocess, sys, tempfile, collections, statistics as st

D    = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
A    = "/home/juanfra/winloci_scratch/o1_antifp"
GEN  = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
OUT  = "/mnt/linuxdisk/home/juanfraitu/o1_flank"
FLANK   = int(os.environ.get("O1_FLANK", "5000"))
MIN_ID  = 0.80
MIN_BP  = 300
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(RC)[::-1]


def main():
    os.makedirs(OUT, exist_ok=True)
    cp = {f"{r['family_id']}~{r['copy_idx']}": (r["chrom"], int(r["start"]), int(r["end"]), r["strand"])
          for r in csv.DictReader(open(f"{D}/GGO_gwcat.copies.tsv"), delimiter="\t")}
    fp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/fp14_detail.tsv"), delimiter="\t")]
    tp = [(r["a"], r["b"]) for r in csv.DictReader(open(f"{A}/tp150_detail.tsv"), delimiter="\t")]
    pairs = [("FP", a, b) for a, b in fp if a in cp and b in cp] + \
            [("TP", a, b) for a, b in tp if a in cp and b in cp]
    names = sorted({n for _, a, b in pairs for n in (a, b)})

    # flanks ONLY — the gene body is excluded, so any alignment is extra-genic homology
    regions, want = [], {}
    for n in names:
        c, s, e, strand = cp[n]
        L = (c, max(0, s - FLANK), s)
        R = (c, e, e + FLANK)
        want[n] = (f"{L[0]}:{L[1]+1}-{L[2]}", f"{R[0]}:{R[1]+1}-{R[2]}")
        regions += [want[n][0], want[n][1]]
    open(f"{OUT}/regions.txt", "w").write("\n".join(regions) + "\n")
    p = subprocess.run(["samtools", "faidx", "-r", f"{OUT}/regions.txt", GEN],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("faidx failed: " + p.stderr[:300])
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

    flank = {}
    for n in names:
        l, r = want[n]
        a, b = seqs.get(l, ""), seqs.get(r, "")
        if cp[n][3] == "-":
            a, b = revcomp(b), revcomp(a)
        flank[n] = (a, b)   # (upstream, downstream) in transcript orientation

    idx = {n: i for i, n in enumerate(names)}
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "f.fa")
        with open(fa, "w") as fh:
            for n in names:
                for side, s in zip(("U", "D"), flank[n]):
                    if len(s) >= MIN_BP:
                        fh.write(f">{idx[n]}_{side}\n{s}\n")
        r = subprocess.run(["minimap2", "-c", "-X", "--no-long-join", "-t", "4", "-k", "11", "-w", "5",
                            fa, fa], capture_output=True, text=True)
    hit = collections.defaultdict(int)
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qa, ta = f[0].split("_")[0], f[5].split("_")[0]
        if qa == ta:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
        span = int(f[3]) - int(f[2])
        if ident >= MIN_ID and span >= MIN_BP:
            k = (min(int(qa), int(ta)), max(int(qa), int(ta)))
            hit[k] = max(hit[k], span)

    print(f"flank {FLANK} bp per side, gene bodies EXCLUDED; alignment >= {MIN_BP} bp at id >= {MIN_ID}\n")
    print("=== does extra-genic sequence align between the two loci? ===")
    print(f"  {'arm':<6} {'flank-linked':>14} {'rate':>8}   median aligned bp")
    res = {}
    for lab in ("FP", "TP"):
        P = [(a, b) for k, a, b in pairs if k == lab]
        v = [hit.get((min(idx[a], idx[b]), max(idx[a], idx[b])), 0) for a, b in P]
        n = sum(1 for x in v if x > 0)
        res[lab] = (n, len(P))
        med = st.median([x for x in v if x > 0]) if n else float("nan")
        print(f"  {lab:<6} {n:>7}/{len(P):<6} {n/len(P):>8.4f}   {med:.0f}" if n
              else f"  {lab:<6} {n:>7}/{len(P):<6} {n/len(P):>8.4f}   n/a")
    (fn, ft), (tn, tt) = res["FP"], res["TP"]
    print(f"\n  As an ADMISSION certificate (require flank homology): admits {fn}/{ft} FP, keeps {tn}/{tt} TP")


if __name__ == "__main__":
    main()
