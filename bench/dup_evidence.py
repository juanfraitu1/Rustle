#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §1a): duplication evidence D1-D4 in one internal format.

Pair table (TSV, header = PAIR_COLS): 0-based half-open intervals of side A and side B, strands, identity, CIGAR (M/I/D;
M both sides, D side A, I side B; side B reverse-complemented when strand_b is '-'), source (D1:sedef, D2:biser, D3:selfaln).
"""
import csv
import re

PAIR_COLS = ("chrom_a", "start_a", "end_a", "chrom_b", "start_b", "end_b", "strand_a", "strand_b", "identity", "cigar", "source")
MIN_SD_BP, MIN_SD_ID = 1000, 0.90


def normalize_cigar(cg):
    ops = [(int(n), "M" if o in "=X" else o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cg)]
    out = []
    for n, o in ops:
        if out and out[-1][1] == o:
            out[-1] = (out[-1][0] + n, o)
        else:
            out.append((n, o))
    return "".join(f"{n}{o}" for n, o in out)


def from_sedef(line, fmt):
    f = line.rstrip("\n").split("\t")
    if len(f) < 23 or not f[1].isdigit():
        return None
    if fmt == "gorilla":
        ident, cig = float(f[20]), (f[32] if len(f) > 32 else "")
    elif fmt == "human":
        ident, cig = float(f[22]), ""
    else:
        raise ValueError(fmt)
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], ident, cig, "D1:sedef")


def from_biser(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < 13 or not f[1].isdigit():
        return None
    m = re.search(r"X=([\d.]+)", f[13]) if len(f) > 13 else None
    err = float(m.group(1)) if m else float(f[7])
    return (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[8], f[9], round(1 - err / 100, 9),
            normalize_cigar(f[12]), "D2:biser")


def from_selfpaf(line):
    """minimap2 self-alignment record of a chunk (query `chrom@offset`) against the genome. Side A = target (forward),
    side B = query with the record's strand, so minimap2's CIGAR (D = target, I = query) already has SEDEF semantics.
    Keeps canonical non-self records >= MIN_SD_BP at identity >= MIN_SD_ID."""
    f = line.rstrip("\n").split("\t")
    qc, off = f[0].rsplit("@", 1)
    qs, qe = int(off) + int(f[2]), int(off) + int(f[3])
    tc, ts, te = f[5], int(f[7]), int(f[8])
    if qc == tc and qs < te and ts < qe:
        return None
    nm, bl = int(f[9]), int(f[10])
    if bl < MIN_SD_BP or nm / bl < MIN_SD_ID:
        return None
    if (tc, ts) > (qc, qs):
        return None
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    return (tc, ts, te, qc, qs, qe, "+", f[4], round(nm / bl, 9), normalize_cigar(cg), "D3:selfaln")


def to_sedef_gorilla(p):
    ca, a1, a2, cb, b1, b2, sa, sb, ident, cig, _src = p
    L = sum(int(n) for n, o in re.findall(r"(\d+)([MID])", cig) if o == "M") or max(a2 - a1, b2 - b1)
    m = round(ident * L)
    mm = L - m
    frac = m / (m + mm) if L else 0.0
    row = [""] * 34
    row[0:6] = [ca, str(a1), str(a2), cb, str(b1), str(b2)]
    row[6], row[7], row[8], row[9] = "S", f"{(1 - ident) * 100:.1f}", sa, sb
    row[16], row[17], row[20], row[32], row[33] = str(m), str(mm), f"{frac:.6f}", cig, f"{frac:.6f}"
    return "\t".join(row)


def write_pairs(pairs, path):
    with open(path, "w") as fh:
        fh.write("\t".join(PAIR_COLS) + "\n")
        for p in pairs:
            fh.write("\t".join(map(str, p)) + "\n")


def read_pairs(path):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        out.append((r["chrom_a"], int(r["start_a"]), int(r["end_a"]), r["chrom_b"], int(r["start_b"]), int(r["end_b"]),
                    r["strand_a"], r["strand_b"], float(r["identity"]), r["cigar"], r["source"]))
    return out
