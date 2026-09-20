#!/usr/bin/env python3
"""Seeded transitive closure: find a gene family's members from ONE annotated copy.

WHAT THIS IS. Given a genome and one seed interval, align the seed to the genome, take every locus it
hits, use those as the next round's seeds, and repeat until a round adds nothing. The result is the
transitive closure of "is homologous to" starting from the seed.

WHY IT EXISTS, AND WHY IT IS NOT THE EXISTING SEED-EXTEND SCRIPTS. `bench/seed_extend_minimap2.py` and
`bench/annotation_seed_extend.py` both run ONE round against a pre-built database of RNA-derived thin
loci. This searches the GENOME and ITERATES. The iteration is the whole point: measured on gorilla NPIP
(`docs/o1_ledger.md` §5m), one round finds 16 of 31 true loci and iterating finds 25 — and the pipeline's
own RNA node construction finds 13.

VALIDATION (§5m gorilla, §5p human — ⚠ never pool the two):
  * gorilla NPIP, seeded from the single annotated copy NPIPB11: converges in 3 rounds at 25/31 loci.
    Non-member hits constant at 1 EVERY round (precision 0.941 -> 0.962 -> 0.962), so the closure does
    NOT chain through shared repeats. All 6 loci it never reaches have ZERO reads, so it recovers
    23/23 = 1.000 of the EXPRESSED members.
  * human, Soto's 65-family panel, one seed each: 65/65 CONVERGE, median recall 1.000 in every stratum.
    Stratified because segmental duplications are repeat-rich and were expected to chain: SD-like
    families score 0.885 against 0.895 for gene families — the concern does not manifest. Recall is not
    size-dependent (Pearson r = -0.043); families of 14, 16 and 17 members each reach 1.000.

WHAT IT DOES NOT DO. It finds LOCI, not nodes. Feeding these windows to RNA node construction was
measured and REJECTED (§5o): it loses to the unseeded footprint pass, 12/31 against 13/31. This narrows
the annotation dependency from many seeds to ONE; it does not remove it.

⚠ ALWAYS RECORD -p AND -N WITH ANY COPY COUNT. Counts are sensitive to both (a documented case gives 1/1
at -p 0.8 and 9/8 at -p 0.1). They are printed with every result.
"""
from __future__ import annotations
import argparse, collections, subprocess, sys, tempfile, os

PRESET, NSEC, PSEC = "asm20", "200", "0.1"


def _merge(hits, min_bp):
    by = collections.defaultdict(list)
    for chrom, s, e, bp in hits:
        if bp >= min_bp:
            by[chrom].append((s, e))
    out = []
    for chrom, iv in by.items():
        iv.sort()
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s < ce:
                ce = max(ce, e)
            else:
                out.append((chrom, cs, ce)); cs, ce = s, e
        out.append((chrom, cs, ce))
    return out


def _extract(genome, loci, out_fa, tag):
    with open(out_fa + ".regions", "w") as f:
        for c, s, e in loci:
            f.write(f"{c}:{s + 1}-{e}\n")
    r = subprocess.run(["samtools", "faidx", genome, "-r", out_fa + ".regions"],
                       capture_output=True, text=True)
    n = 0
    with open(out_fa, "w") as f:
        for line in r.stdout.split("\n"):
            if line.startswith(">"):
                f.write(f">{tag}_{n}\n"); n += 1
            elif line:
                f.write(line + "\n")
    return n


def _search(genome, seed_fa, threads, min_bp):
    p = subprocess.run(["minimap2", "-x", PRESET, "-c", "--eqx", "-N", NSEC, "-p", PSEC,
                        "-t", str(threads), genome, seed_fa],
                       capture_output=True, text=True)
    hits = []
    for ln in p.stdout.split("\n"):
        f = ln.split("\t")
        if len(f) < 11:
            continue
        hits.append((f[5], int(f[7]), int(f[8]), int(f[10])))
    return _merge(hits, min_bp)


def closure(genome, seed_locus, threads=4, min_bp=500, max_rounds=6, verbose=True):
    """Iterate to a fixed point. Returns (loci, per_round_counts)."""
    work = tempfile.mkdtemp(prefix="closure_")
    cur = os.path.join(work, "r0.fa")
    _extract(genome, [seed_locus], cur, "seed")
    rounds, prev = [], -1
    loci = []
    for rnd in range(1, max_rounds + 1):
        loci = _search(genome, cur, threads, min_bp)
        rounds.append(len(loci))
        if verbose:
            print(f"  round {rnd}: {len(loci)} loci", file=sys.stderr)
        if len(loci) == prev:                       # fixed point
            break
        prev = len(loci)
        cur = os.path.join(work, f"r{rnd}.fa")
        if _extract(genome, loci, cur, f"r{rnd}") == 0:
            break
    return loci, rounds


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="FASTA with a .fai")
    ap.add_argument("--seed", required=True, help="chrom:start-end, 0-based half-open")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-bp", type=int, default=500, help="minimum aligned bp for a hit to count")
    ap.add_argument("--max-rounds", type=int, default=6)
    a = ap.parse_args()
    c, rest = a.seed.split(":"); s, e = rest.split("-")
    loci, rounds = closure(a.genome, (c, int(s), int(e)), a.threads, a.min_bp, a.max_rounds)
    conv = len(rounds) > 1 and rounds[-1] == rounds[-2]
    print(f"# genome={a.genome} seed={a.seed}")
    print(f"# minimap2 -x {PRESET} -c --eqx -N {NSEC} -p {PSEC}   min_bp={a.min_bp}")
    print(f"# rounds={rounds} converged={conv}")
    if not conv:
        print(f"# ⚠ NOT CONVERGED after {a.max_rounds} rounds — the count is a lower bound, "
              f"and a non-converged closure may still be growing", file=sys.stderr)
    for c2, s2, e2 in sorted(loci):
        print(f"{c2}\t{s2}\t{e2}")


if __name__ == "__main__":
    main()
