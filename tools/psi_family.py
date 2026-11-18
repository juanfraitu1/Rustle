#!/usr/bin/env python3
"""The same junction across the copies of a family: psi rows keyed by exon columns (driver stage `psi`, after `exonmap`).

    psi_family.py PSI.tsv EXONMAP.columns.tsv FAM.loci.gff3 OUT.tsv [--min-bunch N] [--introns PSI.introns.tsv --paf FAM.loci.gated.paf]

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
                    replicated. n_copies and differs see only the rows the psi stage kept (an intron needs 3 reads).
Added rows (with --introns and --paf), after the keyed ones, for a member that owns both columns of a key and has no row
under it (the psi stage keeps an intron only with 3 reads and a second intron at one of its sites):
  - its own introns between its exons of the two columns, read from the dump with their real counts (a junction alone at
    its sites, or shifted a few bases against the other copies', is the copy's use of the key, as in the keyed rows): one
    row each, those of its deepest anchor site, the copy's range = their sum;
  - when it has none: the junction of the key's most-read copy carried over through the pair's alignment (the bases
    between the images of the two flanking exonic bases; first record, by matching bases and agreeing with the two
    strands, that maps both). No image: counted as not located. An image whose flanks are not in the copy's exons of the
    two columns: counted as off its exons. Neither gets a row.
A row is written only with >= --min-bunch certain reads on the introns leaving the row's anchor site. Counts and PSI
come from the dump over ALL introns at that site; locus `.`, in_gtf 0, n_copies and differs 0.
    silent          1 on a carried-over junction (no read carries it: range [0, 0]); 0 on every other row
    n_eligible      copies under the key with >= --min-bunch certain reads, added rows included
    range_lo, range_hi   the copy's wide range under the key (every row of the copy carries it)
    gap             range_lo minus the lowest range_hi among the other eligible copies (> 0: this copy lies wholly above
                    one of them); NA when the copy is not eligible or has no eligible partner.
`placed` is where reads align best, not origin (see psi_bunches.py)."""
import argparse, bisect, collections, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bench", "denovo_eval"))
from exonshare import blocks_from_cigar, seg_map, pos as project_point
from psi_bunches import psi

MIN_INTRON = 20   # a projected junction shorter than this is not a junction

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

def load_pair_records(paf, members):
    """{(src, dst): [(segs, src_start, dst_start), ...]} for the PAF records between two different `members` (keys
    CONTIG:START-END, START 1-based; PAF coordinates 0-based inside the span), both directions, most matching bases first."""
    out, start = collections.defaultdict(list), lambda k: int(k.rsplit(":", 1)[1].split("-")[0])
    strand = members if isinstance(members, dict) else {}      # {member: strand}: a record against the two strands is skipped
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        a, b = f[0], f[5]
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if a == b or a not in members or b not in members or not cg:
            continue
        if strand.get(a, ".") in "+-" and strand.get(b, ".") in "+-" and (f[4] == "+") != (strand[a] == strand[b]):
            continue
        bl = blocks_from_cigar(cg, int(f[2]), int(f[3]), int(f[7]), f[4])
        out[(a, b)].append((int(f[9]), seg_map(bl, f[4], True), start(a), start(b)))
        out[(b, a)].append((int(f[9]), seg_map(bl, f[4], False), start(b), start(a)))
    return {k: [r[1:] for r in sorted(v, key=lambda r: -r[0])] for k, v in out.items()}

def project_intron(records, src, dst, start, end):
    """the junction on `dst` that matches intron [start, end] (1-based closed) of `src`: the bases between the images of
    the two exonic bases flanking it, through the first record that maps both. None: no such record, or under MIN_INTRON."""
    for segs, s0, d0 in records.get((src, dst), ()):
        x, y = project_point(segs, start - 1 - s0), project_point(segs, end + 1 - s0)
        if x is not None and y is not None:
            lo, hi = sorted((x + d0, y + d0))
            return (lo + 1, hi - 1) if hi - lo - 1 >= MIN_INTRON else None
    return None

def exons_at(ivs, starts, longest, pos):
    i, out = bisect.bisect_left(starts, pos - longest), []
    while i < len(ivs) and ivs[i][0] <= pos:
        if ivs[i][1] >= pos:
            out.append(ivs[i])
        i += 1
    return out

def load_introns(path):
    """psi_bunches.py's --introns dump -> {(chrom, strand, start, end): [placed, max, placed_w, max_w]}"""
    out = {}
    for l in open(path).read().splitlines()[1:]:
        f = l.split("\t")
        out[(f[0], f[1], int(f[2]), int(f[3]))] = [int(x) for x in f[4:8]]
    return out

def site_of(intron, anchor):
    """the site an intron shares with its bunch: (chrom, strand, anchor, position of its donor or acceptor end)"""
    chrom, strand, s, e = intron
    return (chrom, strand, anchor, (s if strand == "+" else e) if anchor == "donor" else (e if strand == "+" else s))

def site_index(introns):
    """{site_of(intron, kind): [introns]} for both kinds"""
    site = collections.defaultdict(list)
    for k in introns:
        site[site_of(k, "donor")].append(k); site[site_of(k, "acceptor")].append(k)
    return site

def added_rows(members, top, cols, introns, records, min_bunch, stats):
    """rows for the members that own a key's two columns and have none under it (module docstring, "Added rows").
    members: key -> member -> (lo, hi, reads, bunch); top: (key, member) -> that member's most-read psi row in its deepest bunch.
    -> [(row, silent, (range_lo, range_hi, reads))]"""
    own, site, by = collections.defaultdict(dict), site_index(introns), collections.defaultdict(list)
    for (chrom, strand), ivs in cols.items():
        for a, b, fam, mem, col, cls in ivs:
            own[fam].setdefault(mem, (chrom, strand, {}, collections.defaultdict(list)))
            own[fam][mem][2][col] = cls; own[fam][mem][3][col].append((a, b))
    for k in introns:
        by[k[:2]].append(k)
    for v in by.values():
        v.sort()
    fmt = lambda v: "NA" if v is None else f"{v:.4f}"
    inside = lambda exons, p: any(a <= p <= b for a, b in exons)
    out = []
    for key, mm in members.items():
        fam, anchor, ca, co = key
        u = max(mm, key=lambda m: (mm[m][2], m))                 # the copy with the most reads gives the junction to carry over
        r = top[(key, u)]
        for s, (chrom, strand, cl, ex) in sorted(own[fam].items()):
            if s in mm or ca not in cl or co not in cl:
                continue
            left, right = (ex[ca], ex[co]) if (anchor == "donor") == (strand == "+") else (ex[co], ex[ca])
            lst = by.get((chrom, strand), [])
            mine = sorted({k for a, b in left for k in lst[bisect.bisect_left(lst, (chrom, strand, a + 1, 0)):bisect.bisect_left(lst, (chrom, strand, b + 2, 0))]
                           if inside(right, k[3] + 1)})          # the copy's own introns between the two columns, whatever their read count
            if mine:
                deep = max({site_of(k, anchor) for k in mine}, key=lambda t: (sum(introns[i][2] for i in site[t]), -t[3]))
                ks, silent = [k for k in mine if site_of(k, anchor) == deep], 0
            else:
                j = project_intron(records, u, s, int(r["start"]), int(r["end"]))
                if j is None:
                    stats["unlocated"] = stats.get("unlocated", 0) + 1
                    continue
                if not (inside(left, j[0] - 1) and inside(right, j[1] + 1)):
                    stats["off_exon"] = stats.get("off_exon", 0) + 1      # the image does not join the copy's exons of the two columns
                    continue
                ks, silent = [(chrom, strand, *j)], 1
            at = site[site_of(ks[0], anchor)]
            sums = [sum(introns[i][x] for i in at) for x in range(4)]
            if sums[2] < min_bunch:
                continue
            new, lo, hi = [], 0.0, 0.0
            for k in ks:
                c = introns.get(k, [0, 0, 0, 0])
                e, w = psi(c[0], c[1], sums[0], sums[1]), psi(c[2], c[3], sums[2], sums[3])
                lo += w[1] or 0.0; hi += 1.0 if w[2] is None else w[2]
                new.append([fam, s, ".", anchor, ca, co, cl[ca], cl[co], chrom, strand, str(k[2]), str(k[3]), "0", c[0], c[1], *map(fmt, e), c[2], c[3], *map(fmt, w), sums[2], 0, 0])
            out += [(row, silent, (lo, min(1.0, hi), sums[2])) for row in new]
    return out

def key_rows(rows, cols, min_bunch=20, introns=None, records=None, stats=None):
    """rows: psi rows as dicts (strings). -> output rows (lists): the keyed psi rows in input order, then the added rows
    (only with `introns`, psi_bunches.py's --introns dump as {intron: [placed, max, placed_w, max_w]}, and `records`)."""
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
    out, new, top = [], [], {}
    for key, mem, ca, co, r in keyed:
        if r["bunch"] == members[key][mem][3] and (top.get((key, mem)) is None or int(r["placed_w"]) > int(top[(key, mem)]["placed_w"])):
            top[(key, mem)] = r
    stats = {} if stats is None else stats
    add = added_rows(members, top, cols, introns, records, min_bunch, stats) if introns is not None and records is not None else []
    rng_of = {(key, mem): v[:3] for key, mm in members.items() for mem, v in mm.items()}      # (lo, hi, reads) per copy and key
    rng_of.update({((r[0], r[3], r[4], r[5]), r[1]): v for r, _, v in add})
    elig = collections.defaultdict(dict)
    for (key, mem), (lo, hi, n) in rng_of.items():
        if n >= min_bunch:
            elig[key][mem] = (lo, hi)
    def tail(key, mem, silent):
        lo, hi, _ = rng_of[(key, mem)]
        others = [h for m, (_, h) in elig[key].items() if m != mem]
        return [silent, len(elig[key]), round(lo, 4), round(hi, 4), round(lo - min(others), 4) if mem in elig[key] and others else "NA"]
    for key, mem, ca, co, r in keyed:
        lo, hi, n, b = members[key][mem]
        d = b == r["bunch"] and n >= min_bunch and any(m != mem and n2 >= min_bunch and (lo > hi2 or hi < lo2) for m, (lo2, hi2, n2, _) in members[key].items())
        out.append([key[0], mem, r["locus"], key[1], key[2], key[3], ca, co, r["chrom"], r["strand"], r["start"], r["end"], r["in_gtf"],
                    r["placed"], r["max"], r["psi_placed"], r["psi_lo"], r["psi_hi"], r["placed_w"], r["max_w"], r["psi_placed_w"], r["psi_lo_w"], r["psi_hi_w"],
                    bsum[r["bunch"]], len(members[key]), int(d)] + tail(key, mem, 0))
    return out + [r + tail((r[0], r[3], r[4], r[5]), r[1], silent) for r, silent, _ in add]

HEADER = ("family member locus anchor col_anchor col_other class_anchor class_other chrom strand start end in_gtf placed max psi_placed psi_lo psi_hi "
          "placed_w max_w psi_placed_w psi_lo_w psi_hi_w bunch_placed_w n_copies differs silent n_eligible range_lo range_hi gap").split()

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("psi"); p.add_argument("columns"); p.add_argument("gff"); p.add_argument("out"); p.add_argument("--min-bunch", type=int, default=20)
    p.add_argument("--introns", help="psi_bunches.py --introns dump: with --paf, adds the rows of copies under the stage's read floor")
    p.add_argument("--paf", help="FAM.loci.gated.paf: the pair alignments a junction is carried through")
    a = p.parse_args()
    lines = open(a.psi).read().splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    cols, introns, records, stats = load_columns(a.columns, a.gff), None, None, {}
    if a.introns and a.paf:
        introns = load_introns(a.introns)
        records = load_pair_records(a.paf, {iv[3]: k[1] for k, ivs in cols.items() for iv in ivs})
    out = key_rows(rows, cols, a.min_bunch, introns, records, stats)
    fmt = lambda v: f"{v:.4f}" if isinstance(v, float) else str(v)
    with open(a.out + ".tmp", "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.writelines("\t".join(map(fmt, r)) + "\n" for r in out)
    os.replace(a.out + ".tmp", a.out)
    keys, multi, diff = set(), set(), set()
    old = [r for r in out if r[24]]          # added rows carry n_copies 0: the first line counts what it always counted
    for r in old:
        k = tuple(r[i] for i in (0, 3, 4, 5)); keys.add(k)
        if r[24] > 1:
            multi.add(k)
        if r[25]:
            diff.add(k)
    print(f"{len(old)} of {len(rows)} psi rows keyed in {len({r[0] for r in old})} families; {len(keys)} junction keys, {len(multi)} seen in >= 2 copies, "
          f"{len(diff)} with copies whose ranges do not overlap (>= {a.min_bunch} certain reads each)")
    add = out[len(old):]
    if introns is None:
        print("added rows: none (needs --introns and --paf)")
    else:
        print(f"added rows {len(add)} (silent {sum(r[26] for r in add)}) for {len({tuple(r[i] for i in (0, 3, 4, 5)) for r in add})} keys; "
              f"{stats.get('unlocated', 0)} candidate copies not located, {stats.get('off_exon', 0)} with the image off their exons; eligible copies counted with added rows")

if __name__ == "__main__":
    main()
