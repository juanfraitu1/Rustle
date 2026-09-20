#!/usr/bin/env python3
"""Minimizer-FREE cross-chromosome copy discovery via INTRON-CHAIN ALIGNMENT (structural axis).

Replaces the minimizer-LSH candidate generator with a STRUCTURAL one: copies share their ordered
exon-length pattern (the intron chain) even on different chromosomes and even when point-diverged.
A domain-sharer shares only ONE exon's worth, a shared transposon has no chain at all -> structure
is naturally hostile to the two FP modes that bit the sequence pipeline. Blind spot (explicit):
single-exon genes / retrocopies have no intron chain.

Pipeline (NO minimizers, NO sequence k-mers):
  1. per-gene ordered exon-length vector (the intron chain), multi-exon only.
  2. candidate generation: length-window (total-exonic-length ratio) x exon-count +-1   (cheap, no all-pairs).
  3. full Needleman-Wunsch alignment of the two exon-length vectors:
       match(x,y) = +1 if |x-y| <= max(TOL_BP, TOL_FRAC*max(x,y)) else MISS; gap = intron gain/loss.
     structural score = matched aligned exons / min(n_exon_a, n_exon_b).
  4. structural copy iff score >= T_STRUCT; families = connected components; flag cross-chrom.
  5. validate recall (RABL2 + universe multi-exon cross-chrom families); compare to the sequence
     pipeline (crosschrom_graded.tsv): overlap / structure-only / sequence-only.

Pure-python DP, deterministic. Run with any python3.
"""
import argparse
import csv
import math
import os
from collections import defaultdict

CHAINS = "/tmp/gene_chains.tsv"
SEQ_PAIRS = os.path.join(os.path.dirname(__file__), "crosschrom_graded.tsv")
UNIVERSE = os.path.join(os.path.dirname(__file__), "copy_recovery_eval/results/universe.tsv")
OUT = os.path.join(os.path.dirname(__file__), "intronchain_crosschrom_pairs.tsv")

EX_BP, EX_FRAC = 6, 0.10       # exon-length match tolerance
IN_BP, IN_FRAC = 40, 0.20      # intron-length match tolerance (introns vary more; copies preserve)
T_STRUCT = 0.6                 # structural score gate (matched units / shorter chain)
MIN_MATCHED = 4                # >= this many matched (exon+intron) units (RABL2 matches 6; few/chance <4)
# Candidate generation tradeoff (a finding, both modes kept):
#   length  = total-exonic-length window x exon-count+-1. PRECISE + fast; the length constraint does
#             precision work but EXCLUDES partial / large-intron-gain-loss copies. (default)
#   shingle = 2-intron-SHINGLE inverted index (consecutive intron-bin pairs; single introns aren't
#             discriminative). Higher recall (catches partial copies) but broad -> needs a tighter NW
#             gate + ~10x slower at genome scale.
LEN_RATIO = 1.25               # length mode: candidate iff total-exonic-length ratio <= this
IBIN_RATIO = 1.2               # shingle mode: log-bin ratio (introns within ~20% collide)
K_SHINGLE = 2                  # shingle mode: candidate iff >= this many shared 2-intron shingles
CAP_BIN = 500                  # shingle mode: skip shingles shared by > this many genes
MISS = -2.0                    # NW mismatch penalty (structurally incompatible unit pair)
GAP = -1.0                     # NW gap (intron gain/loss / exon skip)


def _seg_match(sa, sb):
    """A chain UNIT = (exon_len, intron_after or None). Match iff exon lengths agree AND the
    flanking introns agree (both terminal=None, or both intron lengths within tolerance). Coupling
    exon+intron is what makes the signal discriminative: intron lengths span orders of magnitude."""
    ea, ia = sa
    eb, ib = sb
    if abs(ea - eb) > max(EX_BP, EX_FRAC * max(ea, eb)):
        return False
    if (ia is None) != (ib is None):
        return False
    if ia is not None and abs(ia - ib) > max(IN_BP, IN_FRAC * max(ia, ib)):
        return False
    return True


def nw_chain(a, b):
    """Needleman-Wunsch over the full intron chain (list of (exon_len, intron_len|None) units).
    Returns matched units along the best-score path. Gaps model intron gain/loss / exon skip."""
    na, nb = len(a), len(b)
    dp = [[0.0] * (nb + 1) for _ in range(na + 1)]
    mt = [[0] * (nb + 1) for _ in range(na + 1)]
    for i in range(1, na + 1):
        dp[i][0] = i * GAP
    for j in range(1, nb + 1):
        dp[0][j] = j * GAP
    for i in range(1, na + 1):
        ai = a[i - 1]
        for j in range(1, nb + 1):
            ok = _seg_match(ai, b[j - 1])
            diag = dp[i - 1][j - 1] + (1.0 if ok else MISS)
            up = dp[i - 1][j] + GAP
            left = dp[i][j - 1] + GAP
            best, who = diag, 0
            if up > best:
                best, who = up, 1
            if left > best:
                best, who = left, 2
            dp[i][j] = best
            if who == 0:
                mt[i][j] = mt[i - 1][j - 1] + (1 if ok else 0)
            elif who == 1:
                mt[i][j] = mt[i - 1][j]
            else:
                mt[i][j] = mt[i][j - 1]
    return mt[na][nb]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", default="length", choices=["length", "shingle"],
                    help="candidate generation: length-window (precise/fast) or 2-intron-shingle (recall)")
    args = ap.parse_args()
    genes = []
    chrom = {}
    chain = {}
    tot = {}
    with open(CHAINS) as fh:
        fh.readline()
        for line in fh:
            g, c, strand, ne, exl, inl = line.rstrip("\n").split("\t")
            if int(ne) < 2:
                continue
            ex = [int(x) for x in exl.split(",")]
            inn = [int(x) for x in inl.split(",")] if inl else []
            # full intron-chain units: (exon_len, intron_after) ; last exon has no intron
            units = [(ex[i], inn[i] if i < len(inn) else None) for i in range(len(ex))]
            genes.append(g)
            chrom[g] = c
            chain[g] = units
            tot[g] = sum(ex)

    # candidate generation (structural, NO sequence minimizers) — mode selectable.
    n = len(genes)
    if args.cand == "length":
        genes.sort(key=lambda g: tot[g])
        cand = []
        for i in range(n):
            gi = genes[i]; ti = tot[gi]; ni = len(chain[gi])
            k = i + 1
            while k < n and tot[genes[k]] <= ti * LEN_RATIO:
                if abs(len(chain[genes[k]]) - ni) <= 1:
                    cand.append((gi, genes[k]))
                k += 1
        print(f"[candidates] {n} multi-exon genes -> {len(cand)} candidate pairs "
              f"(length-window x exon-count+-1, NO minimizers)")
    else:
        def ibin(x):
            return int(round(math.log(max(x, 1)) / math.log(IBIN_RATIO)))
        inv = defaultdict(list)
        for g in genes:
            bins = [ibin(u[1]) for u in chain[g] if u[1] is not None]
            shingles = {(bins[i], bins[i + 1]) for i in range(len(bins) - 1)}
            for sh in shingles:
                inv[sh].append(g)
        shared = defaultdict(int)
        for sh, gl in inv.items():
            if len(gl) > CAP_BIN or len(gl) < 2:
                continue
            for a in range(len(gl)):
                ga = gl[a]
                for c in range(a + 1, len(gl)):
                    gc = gl[c]
                    shared[(ga, gc) if ga < gc else (gc, ga)] += 1
        cand = [p for p, s in shared.items() if s >= K_SHINGLE]
        print(f"[candidates] {n} multi-exon genes -> {len(cand)} candidate pairs "
              f"(2-intron-shingle index >= {K_SHINGLE} shingles, NO sequence minimizers)")

    confirmed = []
    for a, b in cand:
        ca, cb = chain[a], chain[b]
        m = nw_chain(ca, cb)
        score = m / min(len(ca), len(cb))
        if score >= T_STRUCT and m >= MIN_MATCHED:
            confirmed.append((a, b, round(score, 3), m, chrom[a] != chrom[b]))
    cross = [(a, b, s, m) for a, b, s, m, x in confirmed if x]
    print(f"[structural] confirmed copy pairs={len(confirmed)} (cross-chrom={len(cross)})")

    # families
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b, *_ in confirmed:
        parent[find(a)] = find(b)
    fg = defaultdict(set)
    for a, b, *_ in confirmed:
        fg[find(a)].add(a); fg[find(a)].add(b)
    cross_fams = [m for m in fg.values() if len({chrom[g] for g in m}) > 1]
    print(f"[families] {len(fg)} structural families; cross-chrom={len(cross_fams)}")

    # validation: universe cross-chrom families with >=2 MULTI-EXON members recovered?
    g2f = {}
    with open(UNIVERSE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g2f.setdefault(r["gene_id"], r["family_id"])
    g2disc = {}
    for fi, m in enumerate(fg.values()):
        for g in m:
            g2disc[g] = fi
    uni = defaultdict(set)
    for g, f in g2f.items():
        if g in chain:  # multi-exon only
            uni[f].add(g)
    uni_cross = {f: gs for f, gs in uni.items()
                 if len(gs) >= 2 and len({chrom[g] for g in gs}) > 1}
    rec = 0
    for f, gs in uni_cross.items():
        d = defaultdict(int)
        for g in gs:
            if g in g2disc:
                d[g2disc[g]] += 1
        if d and max(d.values()) >= 2:
            rec += 1
    print(f"[validation] universe cross-chrom families with >=2 multi-exon members: {len(uni_cross)}; "
          f"recovered: {rec}")
    # RABL2 anchor
    for a, b, s, m, x in confirmed:
        if "RABL2" in a and "RABL2" in b:
            print(f"    RABL2 anchor: {a}({len(chain[a])}ex) <-> {b}({len(chain[b])}ex) "
                  f"struct_score={s} matched={m} cross={x}")

    # compare to SEQUENCE pipeline (crosschrom_graded.tsv)
    seq_pairs = set()
    if os.path.exists(SEQ_PAIRS):
        with open(SEQ_PAIRS) as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                seq_pairs.add(tuple(sorted((c[0], c[2]))))
    struct_pairs = {tuple(sorted((a, b))) for a, b, *_ in cross}
    both = struct_pairs & seq_pairs
    struct_only = struct_pairs - seq_pairs
    seq_only = seq_pairs - struct_pairs
    print(f"[vs sequence] cross-chrom pairs: structural={len(struct_pairs)} sequence={len(seq_pairs)} "
          f"| both={len(both)} structure-only={len(struct_only)} sequence-only={len(seq_only)}")

    with open(OUT, "w") as fh:
        fh.write("gene_a\tchrom_a\tgene_b\tchrom_b\tstruct_score\tmatched_exons\tin_sequence_pipeline\n")
        for a, b, s, m in sorted(cross, key=lambda r: -r[2]):
            key = tuple(sorted((a, b)))
            fh.write(f"{a}\t{chrom[a]}\t{b}\t{chrom[b]}\t{s}\t{m}\t{int(key in seq_pairs)}\n")
    print(f"[wrote] {OUT}")


if __name__ == "__main__":
    main()
