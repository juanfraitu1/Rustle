#!/usr/bin/env python3
"""Per-family and per-copy ambiguity of a sample's reads: how "dark" each family is in long-read RNA.

A report over two files the pipeline already writes, no BAM pass:
  BAM.as_table.tsv     one row per molecule (as_table)
  PREFIX.fam.copies.tsv  the families' copy windows

A molecule belongs to the copy whose locus window holds its PRIMARY record's leftmost position. Per copy and per
family it counts three populations:
  dark_mapq   primary MAPQ < 10            (the Dark Region Finder read rule, Ebbert 2019 / Wadsworth 2026)
  as_tie      rival AS >= primary AS       (the O2 gate: exact alignment-score tie)
  as_near     rival AS >= primary AS - 10  (the O2 near-tie band)
A family or copy is flagged `dark` when it has >= MIN_READS molecules and >= 0.9 of them are dark_mapq (DRF's 90%).

Not comparable 1:1 with DRF: DRF is per reference base on DNA with secondaries off; this is per molecule on RNA.
shortcut: membership by leftmost position only, so a read starting upstream of the window is not counted;
upgrade to interval overlap if as_table ever writes the primary's end.

usage: family_ambiguity.py AS_TABLE COPIES_TSV OUT_PREFIX [EXCLUDE.txt]  -> OUT_PREFIX.ambiguity.{families,copies}.tsv
EXCLUDE.txt: family ids, one per line (e.g. PREFIX.assign.excluded.txt, the satellite / rDNA list); their rows are left out.
"""
import sys
from bisect import bisect_right
from collections import defaultdict

MIN_READS, DARK_FRAC, MAPQ_CUT, NEAR = 10, 0.9, 10, 10


def load_copies(path):
    by_chrom, meta = defaultdict(list), {}
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        c = {k: h.index(k) for k in ("family_id", "copy_idx", "tid", "chrom", "locus_start", "locus_end")}
        for line in f:
            p = line.rstrip("\n").split("\t")
            key = (p[c["family_id"]], p[c["copy_idx"]])
            s, e = int(p[c["locus_start"]]), int(p[c["locus_end"]])
            meta[key] = (p[c["tid"]], p[c["chrom"]], s, e)
            by_chrom[p[c["chrom"]]].append((s, e, key))
    idx = {}
    for ch, iv in by_chrom.items():
        iv.sort()
        idx[ch] = ([s for s, _, _ in iv], iv, max(e - s for s, e, _ in iv))
    return idx, meta


def hits(idx, chrom, pos):
    """Copies whose window holds pos."""
    if chrom not in idx:
        return []
    starts, iv, maxlen = idx[chrom]
    out, i = [], bisect_right(starts, pos) - 1
    while i >= 0 and iv[i][0] >= pos - maxlen:
        if iv[i][1] >= pos:
            out.append(iv[i][2])
        i -= 1
    return out


def classify(p):
    """(dark_mapq, as_tie, as_near) for one as_table row."""
    pas, mapq, rival = int(p[5]), int(p[9]), p[11]
    r = None if rival == "NA" else int(rival)
    return (mapq < MAPQ_CUT, r is not None and r >= pas, r is not None and r >= pas - NEAR)


def main(as_table, copies, out, exclude=None):
    idx, meta = load_copies(copies)
    skip = {l.split()[0] for l in open(exclude) if l.strip()} if exclude else set()
    cop = defaultdict(lambda: [0, 0, 0, 0])  # n, dark, tie, near
    fam = defaultdict(lambda: [0, 0, 0, 0])
    with open(as_table) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if p[15] == "NA":
                continue
            hs = hits(idx, p[4], int(p[15]))
            if not hs:
                continue
            add = (1,) + tuple(int(x) for x in classify(p))
            for k in hs:
                for j in range(4):
                    cop[k][j] += add[j]
            for fid in {k[0] for k in hs}:  # once per family
                for j in range(4):
                    fam[fid][j] += add[j]

    def row(v):
        n = v[0]
        fr = [f"{x / n:.3f}" if n else "NA" for x in v[1:]]
        dark = int(n >= MIN_READS and v[1] / n >= DARK_FRAC)
        return [str(n), str(v[1]), fr[0], str(v[2]), fr[1], str(v[3]), fr[2], str(dark)]

    cols = ["n_molecules", "dark_mapq", "frac_dark_mapq", "as_tie", "frac_as_tie", "as_near", "frac_as_near", "dark"]
    ncop = defaultdict(int)
    darkcop = defaultdict(int)
    with open(out + ".ambiguity.copies.tsv", "w") as o:
        o.write("\t".join(["family_id", "copy_idx", "tid", "chrom", "locus_start", "locus_end"] + cols) + "\n")
        for k in sorted(meta, key=lambda k: (k[0], int(k[1]))):
            if k[0] in skip:
                continue
            r = row(cop[k])
            ncop[k[0]] += 1
            darkcop[k[0]] += int(r[-1])
            o.write("\t".join([k[0], k[1]] + [str(x) for x in meta[k]] + r) + "\n")
    with open(out + ".ambiguity.families.tsv", "w") as o:
        o.write("\t".join(["family_id", "n_copies", "dark_copies"] + cols) + "\n")
        for fid in sorted(ncop, key=lambda x: -fam[x][0]):
            o.write("\t".join([fid, str(ncop[fid]), str(darkcop[fid])] + row(fam[fid])) + "\n")


def _selftest():
    idx = {"c": ([10, 20, 100], [(10, 50, "a"), (20, 30, "b"), (100, 110, "c")], 40)}
    assert sorted(hits(idx, "c", 25)) == ["a", "b"] and hits(idx, "c", 40) == ["a"]
    assert hits(idx, "c", 60) == [] and hits(idx, "c", 100) == ["c"] and hits(idx, "x", 1) == []
    row = "r 90 90 2 c 90 90 90 90 0 1 90 90 1 c 5 9".split()
    assert classify(row) == (True, True, True)
    row[9], row[11] = "60", "75"
    assert classify(row) == (False, False, False)
    row[11] = "NA"
    assert classify(row) == (False, False, False)


if __name__ == "__main__":
    _selftest()
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(*sys.argv[1:])
