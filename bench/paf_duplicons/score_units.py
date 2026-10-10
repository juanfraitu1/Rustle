#!/usr/bin/env python3
"""Score the expansion level built on PAF units against the one built on DupMasker (docs/PREREG_paf_duplicon_units_2026-10-10.md).

    score_units.py <truth_dir> <ours_dir> <dev|heldout> <names.tsv> <loci.gff3> <dupmasker.bed> <units.bed>

Both dirs hold `ours.hierarchy.{members,external}.tsv`. Per family (>= 2 members) the cover = its expansions plus its `unique` class (matched only with
`unique`); family score = maximum-weight bipartite matching on Jaccard / max(#blocks). Only the requested half of the family-hash split is scored
(P1, P2, exact sensitivity / precision); P3 (Soto placement) and P4 (named fusions) are all-family checks. Pure functions are tested in test_score_units.py."""
import collections
import hashlib
import json
import math
import sys

NAMED_FUSIONS = {"CHM13_G0020726": "PKD1P6-NPIPP1", "CHM13_G0022124": "PDXDC2P-NPIPB14P"}
NPIP_FAMILY = "MCL9"


def covers(rows):
    """rows = [(family, member, expansions column)] -> {family: {block: set(members)}}; '<fam>.unique' -> block 'unique'"""
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    for fam, m, ex in rows:
        for x in ex.split(","):
            out[fam]["unique" if x.endswith(".unique") else x].add(m)
    return {f: dict(v) for f, v in out.items()}


def _jac(a, b):
    return len(a & b) / len(a | b) if a or b else 0.0


def family_score(d, p):
    from scipy.optimize import linear_sum_assignment
    dk, pk = sorted(d), sorted(p)
    w = [[_jac(d[x], p[y]) if (x == "unique") == (y == "unique") else 0.0 for y in pk] for x in dk]
    r, c = linear_sum_assignment(w, maximize=True)
    return sum(w[i][j] for i, j in zip(r, c)) / max(len(dk), len(pk))


def b1(d):
    return {"all": set().union(*d.values())}


def structured(d):
    return len(d) >= 2


def half(fam):
    return "dev" if int(hashlib.sha256(fam.encode()).hexdigest(), 16) % 2 == 0 else "heldout"


def exact_sens_prec(dc, pc, fams):
    """-> ((exact, DupMasker expansions), (exact, our expansions)) over the families, `unique` blocks excluded"""
    hit = nd = npp = 0
    for f in fams:
        ds = {frozenset(v) for k, v in dc[f].items() if k != "unique"}
        ps = {frozenset(v) for k, v in pc[f].items() if k != "unique"}
        hit += len(ds & ps)
        nd += len(ds)
        npp += len(ps)
    return (hit, nd), (hit, npp)


def sign_test(wins, losses):
    """one-sided P(X >= wins), X ~ Binomial(wins + losses, 1/2)"""
    n = wins + losses
    return sum(math.comb(n, k) for k in range(wins, n + 1)) / 2 ** n if n else 1.0


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered_bp(iv, bed):
    """bases of the merged intervals `iv` covered by the merged, sorted intervals `bed`"""
    import bisect
    starts = [s for s, _ in bed]
    tot = 0
    for s, e in iv:
        j = max(0, bisect.bisect_right(starts, s) - 1)
        while j < len(bed) and bed[j][0] < e:
            tot += max(0, min(e, bed[j][1]) - max(s, bed[j][0]))
            j += 1
    return tot


def _members(path):
    rows = [ln.rstrip("\n").split("\t") for ln in open(path)]
    h = {k: i for i, k in enumerate(rows[0])}
    return h, rows[1:]


def _bed(path):
    by = collections.defaultdict(list)
    for ln in open(path):
        f = ln.split("\t")
        if len(f) >= 4 and not ln.startswith(("#", "track")):
            by[f[0]].append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in by.items()}


def _region(r):
    c, x = r.rsplit(":", 1)
    s, e = x.split("-")
    return c, int(s) - 1, int(e)


def main(truth, ours, which, names, gff, dupbed, unitbed):
    hd, rd = _members(f"{truth}/ours.hierarchy.members.tsv")
    hp, rp = _members(f"{ours}/ours.hierarchy.members.tsv")
    key = lambda h, r: (r[h["family"]], r[h["member"]], r[h["expansions"]])  # noqa: E731
    dc, pc = covers([key(hd, r) for r in rd]), covers([key(hp, r) for r in rp])
    assert sorted(dc) == sorted(pc)
    fams = sorted(f for f in dc if half(f) == which)
    sc = {f: family_score(dc[f], pc[f]) for f in fams}
    st = [f for f in fams if structured(dc[f])]
    w = l = t = 0
    for f in st:
        b = family_score(dc[f], b1(dc[f]))
        w, l, t = w + (sc[f] > b), l + (sc[f] < b), t + (sc[f] == b)
    sens, prec = exact_sens_prec(dc, pc, fams)
    size = {f: len(set().union(*dc[f].values())) for f in fams}
    res = dict(half=which, families=len(fams), P1_mean_family_score=round(sum(sc.values()) / len(fams), 4),
               member_weighted_score=round(sum(sc[f] * size[f] for f in fams) / sum(size.values()), 4),
               families_identical=sum(v == 1.0 for v in sc.values()),
               P2_structured=len(st), P2_wins=w, P2_losses=l, P2_ties=t, P2_p=sign_test(w, l),
               b1_mean_on_structured=round(sum(family_score(dc[f], b1(dc[f])) for f in st) / max(1, len(st)), 4),
               ours_mean_on_structured=round(sum(sc[f] for f in st) / max(1, len(st)), 4),
               exact_sensitivity=f"{sens[0]}/{sens[1]}", exact_precision=f"{prec[0]}/{prec[1]}",
               blocks_per_family_dupmasker=dict(sorted(collections.Counter(len(dc[f]) for f in fams).items())),
               blocks_per_family_ours=dict(sorted(collections.Counter(len(pc[f]) for f in fams).items())),
               worst=sorted((round(v, 3), f) for f, v in sc.items())[:12])
    if NPIP_FAMILY in fams:
        res["npip"] = dict(dupmasker={k: len(v) for k, v in dc[NPIP_FAMILY].items()}, ours={k: len(v) for k, v in pc[NPIP_FAMILY].items()},
                           score=round(sc[NPIP_FAMILY], 4))
    # P3: Soto placement (all families)
    for tag, d in (("dupmasker", truth), ("ours", ours)):
        ex = _members(f"{d}/ours.hierarchy.external.tsv")
        h, rows = ex
        res[f"P3_soto_expansion_{tag}"] = sum(int(r[h["n_members"]]) >= 2 and r[h["placement"]] == "EXPANSION" for r in rows)
        res[f"P3_soto_placement_{tag}"] = dict(collections.Counter(r[h["placement"]] for r in rows if int(r[h["n_members"]]) >= 2))
    # P4: named fusions (all families)
    reg = dict(ln.rstrip("\n").split("\t")[:2] for ln in open(names))
    for gid, nm in NAMED_FUSIONS.items():
        c, s, e = _region(reg[gid])
        for tag, h, rows in (("dupmasker", hd, rd), ("ours", hp, rp)):
            best = max(rows, key=lambda r: (_region(r[h["member"]])[0] == c) * max(0, min(e, _region(r[h["member"]])[2]) - max(s, _region(r[h["member"]])[1])))
            res[f"P4_{nm}_{tag}"] = f"{best[h['family']]} {best[h['member']]} fusion_of={best[h['fusion_of']]} {best[h['fusion_kind']]}"
    # descriptive: members' exonic bases carrying a unit / a DupMasker duplicon (all members)
    ex = collections.defaultdict(list)
    for ln in open(gff):
        f = ln.split("\t")
        if len(f) > 4 and f[2] == "exon":
            ex[f[0]].append((int(f[3]) - 1, int(f[4])))
    ex = {c: merge(v) for c, v in ex.items()}
    mem = collections.defaultdict(list)
    for r in rd:
        c, s, e = _region(r[hd["member"]])
        mem[c].append((s, e))
    tot = cov_d = cov_p = 0
    dbed, pbed = _bed(dupbed), _bed(unitbed)
    for c, spans in mem.items():
        iv = [[max(s, a), min(e, b)] for a, b in merge(spans) for s, e in ex.get(c, []) if s < b and a < e]
        iv = merge(iv)
        tot += sum(e - s for s, e in iv)
        cov_d += covered_bp(iv, dbed.get(c, []))
        cov_p += covered_bp(iv, pbed.get(c, []))
    res["member_exonic_bp"] = tot
    res["share_exonic_with_dupmasker_duplicon"] = round(cov_d / tot, 4)
    res["share_exonic_with_paf_unit"] = round(cov_p / tot, 4)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:8])
