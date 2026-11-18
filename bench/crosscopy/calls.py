#!/usr/bin/env python3
"""Calls of the cross-copy splicing test (docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md, section 1).

    calls.py FAM.psi.tsv [OUT.tsv]      class counts on stdout; the calls as a table when OUT is given

Input = PREFIX.fam.psi.tsv of stage `psi` (tools/psi_family.py, with its columns silent .. gap). A key is one junction read
across the copies of a family: (family, anchor, column of the anchor, column of the other end). Per key, among the copies with
>= --min-bunch certain reads: U = the copy with the highest range_lo, S = the other copy with the lowest range_hi,
gap = range_lo(U) - range_hi(S).
    LARGE  gap > 0.2      SMALL  0 < gap <= 0.2      SAME  the ranges overlap and the copies' psi_placed_w lie within 0.05
Keys of one eligible copy, and overlapping keys that are not SAME, give no call. An intron sits in a donor and an acceptor
bunch: per (family, U's intron) the key with the larger gap is kept; and one event read from both copies (the inclusion
key with U = a and the skip key with U = b at the same anchor column) is one call: per (family, anchor, anchor column,
the two copies) the orientation with the larger gap stays, the earlier key on a tie. COPY-SPECIFIC = LARGE with range_hi(S) <= 0.05.
A copy's numbers are read from its rows (never from the gap / n_eligible columns), so a table cut to some rows gives the
calls of those rows."""
import collections, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from psi_family import HEADER, site_of, project_intron

LARGE, SAME_POINTS, SPECIFIC = 0.2, 0.05, 0.05
Call = collections.namedtuple("Call", "family anchor col_anchor col_other cls gap u s u_intron s_intron reads_u reads_s lo_u hi_s silent copy_specific u_set s_set miskeyed")
Copy = collections.namedtuple("Copy", "lo hi reads point intron silent introns")
MIN_READS = 3   # an intron with fewer reads is not a used junction (the psi stage's floor)

def copies(path):
    """{key: {member: Copy}}; a copy's rows = those of its deepest bunch under the key."""
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    lines = open(path).read().splitlines()
    head = lines[0].split("\t")
    for l in lines[1:]:
        r = dict(zip(head, l.split("\t")))
        by[(r["family"], r["anchor"], r["col_anchor"], r["col_other"])][r["member"]].append(r)
    out = {}
    for key, mm in by.items():
        out[key] = {}
        for mem, rows in mm.items():
            n = max(int(r["bunch_placed_w"]) for r in rows)
            deep = [r for r in rows if int(r["bunch_placed_w"]) == n]
            top = max(deep, key=lambda r: int(r["placed_w"]))
            pts = [r["psi_placed_w"] for r in deep]
            out[key][mem] = Copy(float(top["range_lo"]), float(top["range_hi"]), n, None if "NA" in pts else sum(map(float, pts)),
                                 (top["chrom"], top["strand"], int(top["start"]), int(top["end"])), any(r["silent"] == "1" for r in rows),
                                 frozenset((r["chrom"], r["strand"], int(r["start"]), int(r["end"])) for r in deep))
    return out

def load_calls(path, min_bunch=20):
    best = {}
    for key, mm in sorted(copies(path).items()):
        el = {m: c for m, c in mm.items() if c.reads >= min_bunch}
        if len(el) < 2:
            continue
        u = max(el, key=lambda m: (el[m].lo, m))
        s = min((m for m in el if m != u), key=lambda m: (el[m].hi, m))
        gap = round(el[u].lo - el[s].hi, 4)
        pts = [c.point for c in el.values()]
        cls = "LARGE" if gap > LARGE else "SMALL" if gap > 0 else "SAME" if None not in pts and round(max(pts) - min(pts), 4) <= SAME_POINTS else None
        if cls is None:
            continue
        c = Call(*key, cls, gap, u, s, el[u].intron, el[s].intron, el[u].reads, el[s].reads, el[u].lo, el[s].hi, el[s].silent,
                 cls == "LARGE" and el[s].hi <= SPECIFIC, el[u].introns, el[s].introns, False)
        k = (key[0], c.u_intron)
        if k not in best or gap > best[k].gap:
            best[k] = c
    # one event read from both copies (the inclusion key with U = a, the skip key with U = b, at one anchor column): of the
    # two orientations of a pair of copies there, the one holding the larger gap stays (the earlier key on a tie)
    side = {}
    for c in sorted(best.values(), key=lambda c: c[:4]):
        k = (c.family, c.anchor, c.col_anchor, frozenset((c.u, c.s)))
        if k not in side or c.gap > side[k][0]:
            side[k] = (c.gap, c.u)
    return sorted((c for c in best.values() if side[(c.family, c.anchor, c.col_anchor, frozenset((c.u, c.s)))][1] == c.u), key=lambda c: c[:4])

def flag_miskeyed(calls, dump, records):
    """calls with `miskeyed` set: U's junction, carried to S through the pair's alignment, is an intron S uses (>= 3 reads)
    that is not among S's introns under the key. The exon columns then put the two copies' matching junctions under
    different keys, and the key's low range at S is not S's use of U's junction: such a call is counted and judged nowhere."""
    out = []
    for c in calls:
        j = project_intron(records, c.u, c.s, c.u_intron[2], c.u_intron[3])
        k = j and (c.s_intron[0], c.s_intron[1], *j)
        out.append(c._replace(miskeyed=bool(k) and k not in c.s_set and dump.get(k, (0, 0, 0, 0))[3] >= MIN_READS))
    return out

def unit(call, genome):
    """what is counted once in a pooled bar: the reference genome, U's junction, and S by contig, strand and anchor site
    (a member's name is its locus span in one sample's assembly and differs between samples)"""
    chrom, strand, _, pos = site_of(call.s_intron, call.anchor)
    return (genome, *call.u_intron, chrom, strand, pos)

def unit_name(call):
    """the unit without its genome, as one string (the `unit` column of site_cause.py's table)"""
    return "{}:{}:{}-{}|{}:{}:{}".format(*unit(call, "")[1:])

def main():
    C = load_calls(sys.argv[1])
    n = collections.Counter(c.cls for c in C)
    print(f"calls {len(C)}: LARGE {n['LARGE']} (copy-specific {sum(c.copy_specific for c in C)}, silent S {sum(c.silent for c in C if c.cls == 'LARGE')}), "
          f"SMALL {n['SMALL']}, SAME {n['SAME']}; families with a LARGE call {len({c.family for c in C if c.cls == 'LARGE'})}")
    if len(sys.argv) > 2:
        flat = lambda c: [*c[:8], *c.u_intron, *c.s_intron, *c[10:14], int(c.silent), int(c.copy_specific)]      # miskeyed needs the alignments: site_cause.py
        with open(sys.argv[2] + ".tmp", "w") as fh:
            fh.write("family anchor col_anchor col_other cls gap u s u_chrom u_strand u_start u_end s_chrom s_strand s_start s_end reads_u reads_s lo_u hi_s silent copy_specific".replace(" ", "\t") + "\n")
            fh.writelines("\t".join(map(str, flat(c))) + "\n" for c in C)
        os.replace(sys.argv[2] + ".tmp", sys.argv[2])

if __name__ == "__main__":
    main()
