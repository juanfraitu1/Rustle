#!/usr/bin/env python3
"""Score `copy_assign --psv-family-best` against the old rule (prereg docs/archive/2026-10/merged/PREREG_o3_family_best_2026-10-15.md).

    score_fb.py LABELS OLD_RM NEW_RM OLD_AB NEW_AB OUT.json

Each PREFIX is a copy_assign --out (score.load: PREFIX.assignments.tsv, .psv_clusters.tsv, .psv_cluster_reads.tsv). The
families are those with a masked copy in LABELS (role D). S1, S1b and G read the added-back arms; S2 reads the removed arms
with register 1294's recovery rule over the masked copy's PLACED reads, which the filter does not change.
"""
import argparse
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rd_lib as L  # noqa: E402
import score as S  # noqa: E402

S2_MARGIN, S2_MIN_ELIGIBLE = 1, 5


def is_copy(c):
    return c["status"] == "unexplained" and float(c["est_div"]) > L.DELTA


def family_counts(families, clusters):
    """Families with >= 1 `copy` cluster and with >= 1 `explained` cluster, among `families`."""
    fams = set(families)
    copy = {c["family"] for c in clusters if c["family"] in fams and is_copy(c)}
    expl = {c["family"] for c in clusters if c["family"] in fams and c["status"] == "explained"}
    return dict(families=len(fams), copy=len(copy), explained=len(expl))


def removed_side(absent, placed_old, placed_new, old, new, min_placed=S.MIN_PLACED):
    """absent: {family: (masked copy, its reads)}; placed_*: {family: placed reads} of each rule's removed arm, which must be
    equal (the filter acts after placement); old, new: (clusters, members). Eligible = >= min_placed masked-copy reads placed."""
    rows = []
    for fam, (copy, d) in sorted(absent.items()):
        po, pn = placed_old.get(fam, set()), placed_new.get(fam, set())
        if po != pn:
            raise ValueError(f"{fam}: the placed reads differ between the rules ({len(po)} vs {len(pn)})")
        if len(d & po) < min_placed:
            continue
        rows.append(dict(family=fam, masked=copy, n_placed=len(d & po), old=L.score_absent(fam, d, po, *old),
                         new=L.score_absent(fam, d, pn, *new)))
    return dict(eligible=len(rows), recovered_old=sum(r["old"]["recovered"] for r in rows),
                recovered_new=sum(r["new"]["recovered"] for r in rows), rows=rows)


def bars(fams, copy_old, copy_new, expl_old, expl_new, eligible, rec_old, rec_new):
    """The prereg's bars, in integer arithmetic: S1 copy_new / fams <= 0.10; S1b copy_new <= 0.5 copy_old; S2 rec_new >=
    rec_old - 1 when >= 5 copies are eligible; G expl_new >= 0.9 expl_old."""
    s2m = eligible >= S2_MIN_ELIGIBLE
    return dict(S1=round(copy_new / fams, 4) if fams else None, S1_pass=bool(fams) and 10 * copy_new <= fams,
                S1b_pass=None if copy_old == 0 else 2 * copy_new <= copy_old, S2_measurable=s2m,
                S2_pass=(rec_new >= rec_old - S2_MARGIN) if s2m else None, G_pass=10 * expl_new >= 9 * expl_old)


def r1_r2(rm_rows, rule, absent, ab):
    """Register 1294's R1 (recovered / eligible) and R2 (cleared once added back / recovered) for one rule."""
    placed, clusters, members = ab
    rec = [r for r in rm_rows if r[rule]["recovered"]]
    clr = [r for r in rec if L.cleared(r["family"], absent[r["family"]][1], placed.get(r["family"], set()), clusters, members)]
    return dict(recovered=len(rec), cleared=len(clr), labels=dict(collections.Counter(r[rule]["label"] for r in rec)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("labels", "old_rm", "new_rm", "old_ab", "new_ab", "out"):
        ap.add_argument(k)
    a = ap.parse_args(argv)
    absent = S.absent_reads(a.labels)
    fams = sorted(absent)
    o_rm, n_rm, o_ab, n_ab = (S.load(p) for p in (a.old_rm, a.new_rm, a.old_ab, a.new_ab))
    fo, fn = family_counts(fams, o_ab[1]), family_counts(fams, n_ab[1])
    rs = removed_side(absent, o_rm[0], n_rm[0], o_rm[1:], n_rm[1:])
    res = dict(families=len(fams), added_back=dict(old=fo, new=fn),
               copy_families=dict(old=sorted({c["family"] for c in o_ab[1] if c["family"] in absent and is_copy(c)}),
                                  new=sorted({c["family"] for c in n_ab[1] if c["family"] in absent and is_copy(c)})),
               copy_clusters=dict(old=sum(map(is_copy, o_ab[1])), new=sum(map(is_copy, n_ab[1]))),
               removed=dict(eligible=rs["eligible"], recovered_old=rs["recovered_old"], recovered_new=rs["recovered_new"],
                            old=r1_r2(rs["rows"], "old", absent, o_ab), new=r1_r2(rs["rows"], "new", absent, n_ab)),
               bars=bars(len(fams), fo["copy"], fn["copy"], fo["explained"], fn["explained"], rs["eligible"],
                         rs["recovered_old"], rs["recovered_new"]),
               rows=[dict(family=r["family"], masked=r["masked"], n_placed=r["n_placed"],
                          old=dict((k, r["old"][k]) for k in ("recovered", "recall", "purity", "label")),
                          new=dict((k, r["new"][k]) for k in ("recovered", "recall", "purity", "label"))) for r in rs["rows"]])
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps({k: res[k] for k in ("families", "added_back", "copy_clusters", "removed", "bars")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
