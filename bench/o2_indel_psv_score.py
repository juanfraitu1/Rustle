#!/usr/bin/env python3
"""Score PREREG_indel_psv (021446fb): --indel-psv run vs its flag-off base, over the contested set
(origin_rejected == 0 and n_candidates >= 2, the stderr decomposition). Indel columns per molecule =
n_decisive(new) - n_decisive(base) (indel columns are decisive by construction and always observed).

  python3 bench/o2_indel_psv_score.py --base ours_odi.assignments.tsv --new ours_indelA.assignments.tsv \
      [--bam hsa16.bam --copies copies16.tsv [--region R]]   # optional: newly assigned vs PRIMARY copy
"""
import argparse, csv, statistics
from collections import Counter


def rows(path):
    with open(path) as fh:
        out = {r["read_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for r in out.values():
        r["copy"] = r.get("catalog_copy_idx") or r["assigned_copy"]
    return out


def contested(rs):
    return {n for n, r in rs.items() if r["origin_rejected"] == "0" and int(r["n_candidates"]) >= 2}


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--bam")
    ap.add_argument("--copies")
    ap.add_argument("--region")
    a = ap.parse_args()
    base, new = rows(a.base), rows(a.new)
    cb, cn = contested(base), contested(new)
    print(f"contested: base {len(cb)} / new {len(cn)} / symmetric diff {len(cb ^ cn)} (expect 0)")
    st = lambda rs, names: Counter(rs[n]["status"] for n in names)
    print(f"  base {dict(st(base, cb))}\n  new  {dict(st(new, cn))}")
    gained = [n for n in cb if int(new[n]["n_decisive"]) > int(base[n]["n_decisive"])]
    ncols = [int(new[n]["n_decisive"]) - int(base[n]["n_decisive"]) for n in gained]
    print(f"P1 contested molecules with >=1 indel column: {len(gained)}/{len(cb)} = {100*len(gained)/max(1,len(cb)):.1f}%  (pass <=30, refuted >50); columns/molecule med {med(ncols)} max {max(ncols) if ncols else 0}")
    asg = [n for n in cb if base[n]["status"] == "assigned"]
    same = [n for n in asg if new[n]["status"] == "assigned" and new[n]["copy"] == base[n]["copy"]]
    moved = [n for n in asg if new[n]["status"] == "assigned" and new[n]["copy"] != base[n]["copy"]]
    print(f"P2 base assigned {len(asg)}: same status+copy {len(same)} ({100*len(same)/max(1,len(asg)):.1f}%, pass >=99); changed copy {len(moved)} ({100*len(moved)/max(1,len(asg)):.1f}%, refuted >2); "
          f"now not assigned {Counter(new[n]['status'] for n in asg if new[n]['status'] != 'assigned')}")
    conv_pool = [n for n in cb if base[n]["status"] == "tied" and base[n].get("tie_outside_catalog", "0") == "0"]
    conv = [n for n in conv_pool if new[n]["status"] == "assigned"]
    print(f"P3 tied->assigned: {len(conv)} of {len(conv_pool)} convertible ({100*len(conv)/max(1,len(conv_pool)):.1f}%)  [human pass 8..80, refuted 0 or >120]")
    print(f"   other transitions from tied: {Counter(new[n]['status'] for n in conv_pool if new[n]['status'] != 'tied')}")
    amb = [n for n in cb if base[n]["status"] == "ambiguous"]
    print(f"   from ambiguous: {Counter(new[n]['status'] for n in amb)}")
    if conv:
        print(f"   converted: n_decisive med {med(int(new[n]['n_decisive']) for n in conv)}, margin med {med(float(new[n]['margin']) for n in conv):.1f}, "
              f"indel cols med {med(int(new[n]['n_decisive']) - int(base[n]['n_decisive']) for n in conv)}, copies {Counter(new[n]['copy'] for n in conv).most_common(6)}")
    if a.bam and a.copies and conv:
        import pysam
        units = [(r["chrom"], int(r["start"]), int(r["end"]), r["copy_idx"]) for r in csv.DictReader(open(a.copies), delimiter="\t")]
        want = set(conv)
        prim = {}
        with pysam.AlignmentFile(a.bam) as fh:
            it = fh.fetch(region=a.region) if a.region else fh
            for al in it:
                if al.query_name in want and not (al.is_secondary or al.is_supplementary or al.is_unmapped):
                    prim[al.query_name] = next((i for c, s, e, i in units if c == al.reference_name and al.reference_start < e and al.reference_end > s), None)
        k = sum(1 for n in conv if prim.get(n) == new[n]["copy"])
        print(f"   converted -> PRIMARY copy: {k}/{len(conv)} = {100*k/len(conv):.1f}% (report only; under a tie the primary is a coin toss)")
    # whole-file sanity: rows outside the contested set must be identical in status
    out = [n for n in base if n not in cb and (n not in new or new[n]["status"] != base[n]["status"])]
    print(f"sanity: non-contested rows with a changed status = {len(out)} (expect 0); rows base {len(base)} new {len(new)}")


if __name__ == "__main__":
    main()
