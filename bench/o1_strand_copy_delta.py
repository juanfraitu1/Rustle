#!/usr/bin/env python3
"""Where did the 61 copies go? (2,019 OFF -> 1,958 ON under RUSTLE_READ_STRAND)

THE QUESTION THAT DECIDES THE DEFAULT. The read-strand fix passed every pre-registered acceptance
criterion, but the catalog lost 61 copies net. At 0.9650 strand accuracy some merges are wrong, and the
HUMAN negative panel structurally cannot see gorilla-side merges. So: were these loci CORRECTLY merged
(over-fragmentation fixed -- the change doing its job), or genuinely LOST?

⚠ TWO-SIDED. -61 is a NET. Copies can be gained too, and a one-sided count of losses would misread the
change exactly as a one-sided edge ledger would have reported "934 recovered" instead of +108.

⚠ MATCHING. A merged rep is LARGER than either of its parts, so reciprocal-overlap matching FAILS on
precisely the cases under test. "Is this locus still represented?" is an ASYMMETRIC question: does some
ON rep cover >=50% of the OFF copy's span? Using a reciprocal rule here would manufacture losses.

CLASSES for each OFF copy with no ON copy at its position:
  A ABSORBED   an ON *catalog copy* covers >=50% of it -> the locus is still in a family, just merged.
               This is the change working as designed.
  B DEMOTED    an ON *rep* covers it, but that rep is in no family -> the locus survived node
               construction and lost family membership.
  C VANISHED   no ON rep covers it -> the locus left the rep set entirely.
"""
import bisect, collections, csv

R = "/mnt/linuxdisk/home/juanfraitu/o1_strand"


def copies(p):
    return [(r["chrom"], int(r["start"]), int(r["end"]), r["family_id"], int(r["n_exon"]))
            for r in csv.DictReader(open(p), delimiter="\t")]


def nodes(p):
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(p), delimiter="\t"):
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), int(r["degree"]), int(r["n_exon"])))
    for c in by:
        by[c].sort()
    return by


def cover_index(items):
    """chrom -> (sorted spans, starts, max-end prefix) for exact asymmetric coverage queries."""
    idx = {}
    for c, v in items.items():
        v = sorted(v)
        mx, run = [], 0
        for x in v:
            run = max(run, x[1]); mx.append(run)
        idx[c] = (v, [x[0] for x in v], mx)
    return idx


def covers(idx, c, s, e, frac=0.50):
    """Any span covering >= frac of [s,e)? Exact backward sweep, no fixed window."""
    if c not in idx:
        return None
    v, st, mx = idx[c]
    hi = bisect.bisect_left(st, e)
    for j in range(hi - 1, -1, -1):
        if mx[j] <= s:
            break
        ov = min(e, v[j][1]) - max(s, v[j][0])
        if ov > 0 and ov >= frac * (e - s):
            return v[j]
    return None


off, on = copies(f"{R}/off/ggo_off.copies.tsv"), copies(f"{R}/on/ggo_on.copies.tsv")
print(f"OFF copies {len(off)}   ON copies {len(on)}   net {len(on)-len(off):+d}")

on_exact = {(c, s, e) for c, s, e, _, _ in on}
off_exact = {(c, s, e) for c, s, e, _, _ in off}
lost = [x for x in off if (x[0], x[1], x[2]) not in on_exact]
gained = [x for x in on if (x[0], x[1], x[2]) not in off_exact]
print(f"  positionally LOST (no ON copy at that exact span)   {len(lost)}")
print(f"  positionally GAINED (no OFF copy at that exact span) {len(gained)}")

on_cop = cover_index(collections.defaultdict(list, {c: [(s, e) for cc, s, e, _, _ in on if cc == c]
                                                    for c in {x[0] for x in on}}))
on_rep = cover_index({c: [(s, e) for s, e, _, _ in v] for c, v in nodes(f"{R}/on/dump/ggo.nodes.tsv").items()})

cls = collections.Counter(); nex_lost = collections.Counter()
for c, s, e, fam, nx in lost:
    nex_lost["single-exon" if nx == 1 else "spliced"] += 1
    if covers(on_cop, c, s, e):
        cls["A ABSORBED  — an ON CATALOG COPY covers it (merged, still in a family)"] += 1
    elif covers(on_rep, c, s, e):
        cls["B DEMOTED   — an ON REP covers it, but it is in no family"] += 1
    else:
        cls["C VANISHED  — no ON rep covers it"] += 1

n = len(lost)
print(f"\n=== the {n} lost copies ===")
for k, v in cls.most_common():
    print(f"  {k:<62} {v:>4}/{n} = {v/n:.4f}")

base = sum(1 for x in off if x[4] == 1) / len(off)
se = nex_lost["single-exon"] / n
print(f"\n  single-exon among LOST: {nex_lost['single-exon']}/{n} = {se:.4f}")
print(f"  COMPARATOR, single-exon among ALL OFF copies: {sum(1 for x in off if x[4]==1)}/{len(off)} = {base:.4f}"
      f"   ({se/base:.2f}x enrichment)")
