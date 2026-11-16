#!/usr/bin/env python3
"""Drop from a GTF every gene (gene_id) whose exon bases lie mostly (>= FRAC, default 0.5) inside the elements of a BED.
The families stage's opt-in darklist (RUSTLE_FAMILIES_EXCLUDE_BED, register 1749): satellite / rDNA loci never become
family nodes. Same majority rule as assign's RUSTLE_EXCLUDE_BED (annotated_families.py --majority-in), per locus.

usage: gtf_exclude_bed.py IN.gtf ELEMENTS.bed OUT.gtf [FRAC]    (dropped gene ids -> OUT.gtf.excluded.txt)
"""
import re
import sys
from bisect import bisect_right
from collections import defaultdict

GID = re.compile(r'gene_id "([^"]+)"')


def load_bed(path):
    """chrom -> (starts, ends) of merged 1-based closed intervals."""
    raw = defaultdict(list)
    for line in open(path):
        f = line.split("\t")
        if len(f) >= 3 and f[1].isdigit() and f[2].strip().isdigit():
            raw[f[0]].append((int(f[1]) + 1, int(f[2])))
    out = {}
    for c, v in raw.items():
        v.sort()
        m = [list(v[0])]
        for s, e in v[1:]:
            if s <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = ([s for s, _ in m], [e for _, e in m])
    return out


def covered(bed, chrom, s, e):
    """Bases of [s, e] inside the BED."""
    if chrom not in bed:
        return 0
    starts, ends = bed[chrom]
    i, n = max(bisect_right(starts, s) - 1, 0), 0
    while i < len(starts) and starts[i] <= e:
        n += max(0, min(e, ends[i]) - max(s, starts[i]) + 1)
        i += 1
    return n


def main(gtf, bedpath, out, frac=0.5):
    bed = load_bed(bedpath)
    tot, ins = defaultdict(int), defaultdict(int)
    for line in open(gtf):
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "exon":
            continue
        g = GID.search(f[8]).group(1)
        s, e = int(f[3]), int(f[4])
        tot[g] += e - s + 1
        ins[g] += covered(bed, f[0], s, e)
    drop = {g for g in tot if ins[g] >= frac * tot[g]}
    with open(out, "w") as o:
        for line in open(gtf):
            m = GID.search(line)
            if not (m and m.group(1) in drop):
                o.write(line)
    with open(out + ".excluded.txt", "w") as o:
        o.write("".join(g + "\n" for g in sorted(drop)))
    print(f"{len(drop)} of {len(tot)} genes dropped (>= {frac} of exon bases in {bedpath})", file=sys.stderr)


def _selftest():
    import os, tempfile
    d = tempfile.mkdtemp()
    b, g, o = (os.path.join(d, x) for x in ("b.bed", "g.gtf", "o.gtf"))
    open(b, "w").write("c\t10\t20\nc\t15\t30\nc\t100\t110\n")
    bed = load_bed(b)
    assert bed["c"] == ([11, 101], [30, 110])
    assert covered(bed, "c", 1, 15) == 5 and covered(bed, "c", 25, 105) == 11 and covered(bed, "x", 1, 9) == 0
    ex = 'c\tr\texon\t{}\t{}\t.\t+\t.\tgene_id "{}"; transcript_id "t";\n'
    open(g, "w").write("#h\n" + ex.format(11, 20, "in") + ex.format(1, 20, "half") + ex.format(1, 19, "out") + ex.format(40, 50, "out"))
    main(g, b, o)
    kept = open(o).read()
    assert '"out"' in kept and '"in"' not in kept and '"half"' not in kept and kept.startswith("#h")
    assert open(o + ".excluded.txt").read() == "half\nin\n"


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        _selftest()
        sys.exit("ok")
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3], *(float(x) for x in sys.argv[4:]))
