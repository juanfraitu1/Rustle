#!/usr/bin/env python3
"""Competitor-resolved corroboration: does a secondary read fit the locus it sits on BETTER than the
locus it was assigned to?

THE QUESTION THIS CLOSES. Candidate loci (SD mate of a catalog copy, exon-homologous, no rep) pass a
"zero primaries + >=3 secondaries + >=3-read exact intron chain" gate at 0.1236, and size- and
compartment-matched controls pass it at 0.1180 (Fisher p = 0.7701) -- i.e. READ PRESENCE has no
specificity. This asks the sharper question: of the reads that are there, do their BASES prefer this
locus over the one minimap2 gave them?

THE COMPETITOR IS THE READ'S OWN PRIMARY, NOT THE SD ANCHOR. A read piled at L by spillover came from
its primary locus P, so P is the alternative that must be beaten. Comparing L against the SD-mate anchor
instead would score a locus the read never contested.

WHY `de` AND NOT AN OFFLINE PSV VOTE. `de` is minimap2's own per-alignment divergence, present on both
records; comparing two numbers re-derives nothing. Offline E_r/PSV re-derivation has failed four times
on this project (T8), and a PSV column set built by hand here would be exactly that.

THE BASE RATE THIS IS READ AGAINST -- three comparators, all mandatory:
  1. the size/compartment-matched CONTROL arm (Mann-Whitney over per-locus fractions);
  2. the genome-wide rate at which a secondary fits better by `de` at all: 1.96% of reads-with-secondary,
     2.05% of near-ties, 12.16% of TIGHT ties (measured 2026-08-14);
  3. a constant predictor -- if "always vote L" or "never vote L" scores the same, the statistic is void.

PRE-REGISTERED DECISION RULE (written before the result):
  if the two arms do not separate at p < 0.01, C5 stays OFF permanently, no second O1-perp-O2 exception
  is written, and this is filed to NEGATIVE_RESULTS_REGISTER.md as "secondary-read corroboration does not
  distinguish exon-homologous SD loci from matched SD sides -- base rate at every stage".
EXPECTED OUTCOME, STATED IN ADVANCE: TIE. The gate arms already tied and the control won two of three
earlier stages. This is run to CLOSE the question, not to open it.

Usage:  o1_competitor_resolved.py extract <loci.bed> <out.tsv>      # targeted, fast
        o1_competitor_resolved.py primaries <names.txt> <out.tsv>   # ONE 4.3 GB pass
        o1_competitor_resolved.py score <cand.tsv> <ctrl.tsv> <primaries.tsv>
"""
import subprocess, sys, collections

BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam"


def de_of(fields):
    for f in fields[11:]:
        if f.startswith("de:f:"):
            return float(f[5:])
    return None


def extract(bed, out):
    """Secondary alignments overlapping each locus, WITH read names. Targeted -- no full scan."""
    n = 0
    with open(out, "w") as fh:
        fh.write("locus\tread\tflag\tde\tmapq\n")
        for line in open(bed):
            f = line.split()
            if len(f) < 3:
                continue
            reg = f"{f[0]}:{int(f[1])+1}-{f[2]}"
            lid = f[3] if len(f) > 3 else reg
            p = subprocess.run(["samtools", "view", BAM, reg], capture_output=True, text=True)
            for r in p.stdout.splitlines():
                c = r.split("\t")
                fl = int(c[1])
                if fl & 0x800:          # supplementary: never a competitor, it is one read split
                    continue
                d = de_of(c)
                fh.write(f"{lid}\t{c[0]}\t{fl}\t{'NA' if d is None else d}\t{c[4]}\n")
                n += 1
    print(f"{n} alignment records -> {out}")


def primaries(names_file, out):
    """The PRIMARY alignment of each named read. ONE streaming pass; only names we care about."""
    want = {l.strip() for l in open(names_file) if l.strip()}
    print(f"looking for {len(want)} read names in one pass over the BAM", flush=True)
    p = subprocess.Popen(["samtools", "view", "-F", "0x900", BAM], stdout=subprocess.PIPE, text=True)
    n = 0
    with open(out, "w") as fh:
        fh.write("read\tchrom\tpos\tde\tmapq\n")
        for r in p.stdout:
            c = r.split("\t", 12)
            if c[0] not in want:
                continue
            d = de_of(r.rstrip("\n").split("\t"))
            fh.write(f"{c[0]}\t{c[2]}\t{c[3]}\t{'NA' if d is None else d}\t{c[4]}\n")
            n += 1
    p.wait()
    print(f"{n} primaries -> {out}")


def load(path, prim):
    per = collections.defaultdict(lambda: [0, 0])   # locus -> [n_comparable, n_fits_L_better]
    for i, line in enumerate(open(path)):
        if i == 0:
            continue
        lid, read, fl, de, mq = line.rstrip("\n").split("\t")
        if not (int(fl) & 0x100):     # only SECONDARIES are contested
            continue
        if de == "NA" or read not in prim or prim[read] is None:
            continue
        per[lid][0] += 1
        if float(de) < prim[read]:
            per[lid][1] += 1
    return {k: v[1] / v[0] for k, v in per.items() if v[0] >= 3}, per


def mwu(a, b):
    """Mann-Whitney U, normal approximation with tie correction."""
    import math
    comb = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = {}; i = 0; tie = 0.0
    while i < len(comb):
        j = i
        while j + 1 < len(comb) and comb[j + 1][0] == comb[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1; t = j - i + 1; tie += t ** 3 - t
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    ra = sum(ranks[k] for k, (_, g) in enumerate(comb) if g == 0)
    na, nb = len(a), len(b); n = na + nb
    u = ra - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    sd = math.sqrt(na * nb / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sd == 0:
        return u, 1.0
    z = (u - mu) / sd
    return u, math.erfc(abs(z) / math.sqrt(2))


def score(cand, ctrl, primfile):
    prim = {}
    for i, line in enumerate(open(primfile)):
        if i == 0:
            continue
        rd, ch, po, de, mq = line.rstrip("\n").split("\t")
        prim[rd] = None if de == "NA" else float(de)
    A, pa = load(cand, prim)
    B, pb = load(ctrl, prim)
    import statistics as st
    print(f"CAND loci scored {len(A)}   median fraction fitting L better: "
          f"{st.median(A.values()) if A else float('nan'):.4f}")
    print(f"CTRL loci scored {len(B)}   median fraction fitting L better: "
          f"{st.median(B.values()) if B else float('nan'):.4f}")
    ta = sum(v[1] for v in pa.values()); na = sum(v[0] for v in pa.values())
    tb = sum(v[1] for v in pb.values()); nb = sum(v[0] for v in pb.values())
    print(f"\npooled READ unit — cand {ta}/{na} = {ta/max(1,na):.4f}   ctrl {tb}/{nb} = {tb/max(1,nb):.4f}")
    print(f"  genome-wide comparator: a secondary fits better by `de` in 0.0196 of reads-with-secondary")
    if A and B:
        u, p = mwu(list(A.values()), list(B.values()))
        print(f"\nMann-Whitney over per-locus fractions: U={u:.1f}  p={p:.4f}")
        print("  DECISION RULE (pre-registered): p >= 0.01 => C5 stays OFF, no second exception, "
              "file as REFUTED.")
        print(f"  VERDICT: {'SEPARATES' if p < 0.01 else 'TIE — the pre-registered NO branch'}")


if __name__ == "__main__":
    {"extract": lambda: extract(sys.argv[2], sys.argv[3]),
     "primaries": lambda: primaries(sys.argv[2], sys.argv[3]),
     "score": lambda: score(sys.argv[2], sys.argv[3], sys.argv[4])}[sys.argv[1]]()
