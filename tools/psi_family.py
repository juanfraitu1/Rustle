#!/usr/bin/env python3
"""The same junction across the copies of a family: psi rows keyed by exon columns (driver stage `psi`, after `exonmap`).

    psi_family.py PSI.tsv EXONMAP.columns.tsv FAM.loci.gff3 OUT.tsv [--min-bunch N]

PSI.tsv = psi_bunches.py's table; columns.tsv = exonmap.py's exon columns (exons of a family's members that match each
other). An intron is keyed when both its flanking bases lie in column exons of one member on the intron's strand:
key = (family, anchor, col_anchor, col_other), the column of the bunch's shared site and the column the intron leads to,
so "from C5 to C6" at one copy is set beside "from C5 to C6" at another. Only the member's own exons (the representative
transcript `families` chose) are in columns: junctions into other exons are not keyed. Rows = the psi rows that key, plus
    bunch_placed_w  reads certain in this bunch (the sum of placed_w over all its introns, keyed or not)
    n_copies        members with a row under this key
    differs         1 when this copy's wide range [psi_lo_w, psi_hi_w] for the key lies wholly above or below another
                    copy's, both bunches holding >= --min-bunch (default 20) certain reads: a difference the placement
                    ties cannot explain, of any size (most ranges are points, so filter on the PSI columns for a size).
                    A copy's range = that of its deepest bunch under the key, summed over that bunch's introns under the
                    key (shifted sites inside one exon pair); rows of its other bunches carry 0. No test: counts are not
                    replicated, and a copy that never uses the junction has no row at all (it shows as a smaller
                    n_copies, not as PSI 0).
`placed` is where reads align best, not origin (see psi_bunches.py)."""
import argparse, bisect, collections, os, re, sys

def load_columns(columns, gff):
    """{(chrom, strand): sorted [(start, end, family, member, column, class)]} (1-based closed) from columns.tsv + member strands."""
    strand = {}
    for l in open(gff):
        f = l.split("\t")
        if len(f) > 8 and f[2] == "gene":
            strand[f"{f[0]}:{f[3]}-{f[4]}"] = f[6]
    out, lines = collections.defaultdict(list), open(columns).read().splitlines()
    for l in lines[1:]:
        f = l.split("\t")
        for x in f[8].split(","):
            m = re.fullmatch(r"(.+)#\d+:(\d+)-(\d+)", x)
            mem = m.group(1)
            out[(mem.rsplit(":", 1)[0], strand.get(mem, "."))].append((int(m.group(2)) + 1, int(m.group(3)), f[0], mem, f[2], f[3]))
    return {k: sorted(v) for k, v in out.items()}

def exons_at(ivs, starts, longest, pos):
    i, out = bisect.bisect_left(starts, pos - longest), []
    while i < len(ivs) and ivs[i][0] <= pos:
        if ivs[i][1] >= pos:
            out.append(ivs[i])
        i += 1
    return out

def key_rows(rows, cols, min_bunch=20):
    """rows: psi rows as dicts (strings). -> output rows (lists) in input order."""
    idx = {k: ([a for a, *_ in v], max(b - a for a, b, *_ in v)) for k, v in cols.items()}
    num = lambda v: None if v == "NA" else float(v)
    bsum = collections.Counter()
    for r in rows:
        bsum[r["bunch"]] += int(r["placed_w"])
    keyed, rng = [], collections.defaultdict(lambda: [0.0, 0.0])
    for r in rows:
        k = (r["chrom"], r["strand"])
        if k not in cols:
            continue
        s, e = int(r["start"]), int(r["end"])
        left, right = (exons_at(cols[k], *idx[k], p) for p in (s - 1, e + 1))
        don, acc = (left, right) if r["strand"] == "+" else (right, left)
        anc, oth = (don, acc) if r["anchor"] == "donor" else (acc, don)
        for a in anc:
            for o in oth:
                if a[3] == o[3]:
                    key = (a[2], r["anchor"], a[4], o[4])
                    keyed.append((key, a[3], a[5], o[5], r))
                    g = rng[(key, a[3], r["bunch"])]
                    lo, hi = num(r["psi_lo_w"]), num(r["psi_hi_w"])
                    g[0] += lo or 0.0; g[1] += 1.0 if hi is None else hi
    members = collections.defaultdict(dict)  # key -> member -> (lo, hi, certain reads, bunch) of the member's deepest bunch
    for (key, mem, b), (lo, hi) in sorted(rng.items()):
        if mem not in members[key] or bsum[b] > members[key][mem][2]:
            members[key][mem] = (lo, min(1.0, hi), bsum[b], b)
    out = []
    for key, mem, ca, co, r in keyed:
        lo, hi, n, b = members[key][mem]
        d = b == r["bunch"] and n >= min_bunch and any(m != mem and n2 >= min_bunch and (lo > hi2 or hi < lo2) for m, (lo2, hi2, n2, _) in members[key].items())
        out.append([key[0], mem, r["locus"], key[1], key[2], key[3], ca, co, r["chrom"], r["strand"], r["start"], r["end"], r["in_gtf"],
                    r["placed"], r["max"], r["psi_placed"], r["psi_lo"], r["psi_hi"], r["placed_w"], r["max_w"], r["psi_placed_w"], r["psi_lo_w"], r["psi_hi_w"],
                    bsum[r["bunch"]], len(members[key]), int(d)])
    return out

HEADER = ("family member locus anchor col_anchor col_other class_anchor class_other chrom strand start end in_gtf placed max psi_placed psi_lo psi_hi "
          "placed_w max_w psi_placed_w psi_lo_w psi_hi_w bunch_placed_w n_copies differs").split()

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("psi"); p.add_argument("columns"); p.add_argument("gff"); p.add_argument("out"); p.add_argument("--min-bunch", type=int, default=20)
    a = p.parse_args()
    lines = open(a.psi).read().splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    out = key_rows(rows, load_columns(a.columns, a.gff), a.min_bunch)
    with open(a.out + ".tmp", "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in out)
    os.replace(a.out + ".tmp", a.out)
    keys, multi, diff = set(), set(), set()
    for r in out:
        k = tuple(r[i] for i in (0, 3, 4, 5)); keys.add(k)
        if r[24] > 1:
            multi.add(k)
        if r[25]:
            diff.add(k)
    print(f"{len(out)} of {len(rows)} psi rows keyed in {len({r[0] for r in out})} families; {len(keys)} junction keys, {len(multi)} seen in >= 2 copies, "
          f"{len(diff)} with copies whose ranges do not overlap (>= {a.min_bunch} certain reads each)")

if __name__ == "__main__":
    main()
