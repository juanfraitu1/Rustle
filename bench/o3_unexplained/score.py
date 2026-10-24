#!/usr/bin/env python3
"""Score the PSV clustering (psvclust.py) and the existing detector (`copy_assign --flag-missing-copies`) on an o3_unexplained world.

    score.py rows WORLD VARIANT ARM ROWS.tsv [--rust]    # one row per truth unit (absent source or control family); --rust scores
                                                         # copy_assign --psv-clusters' product (D/R.*) instead of psvclust.py's
    score.py summary ROWS.tsv                    # the bars of docs/archive/2026-10/merged/PREREG_o3_unexplained_2026-10-15.md and the reported-beside numbers

Truth units:
  source  an expressed sequence absent from the genome (copy E / F, or allele H) in a REMOVED_* / LOWEXPR / COMBO_EH / DOUBLE / ALLELE
          family. RECOVERED iff some `unexplained` cluster holds >= RECALL of the source's reads that O2 placed on the family, with purity
          (source reads / cluster reads) >= PURITY. Its label: `copy` if the cluster's est_div > DELTA, else `allele_or_copy`.
  control a NULL_* / ADDBACK / DECOY family: a FALSE FLAG is any `unexplained` cluster in it.
Baseline (the existing detector): the source's family has an `o3_flag == missing_copy` row at the copy hosting most of its reads.
"""
import argparse
import collections
import csv
import json
import os
import sys

RECALL, PURITY, DELTA = 0.8, 0.9, 0.00958
CONTROL = ("NULL_TWIN", "NULL_MULTI", "NULL_CLOSE", "ADDBACK", "DECOY",                      # base design
           "MANY_NULL", "NULL_CONV", "HIGHERR_NULL", "EDIT_NULL")                               # hard design
MULTI_GROUP = ("COMBO_EH", "DOUBLE", "HET_COPY")
# truth divergence of each absent source: the scenario's d where it carries one, else the design's value
FIXED_D = {("COMBO_EH", "E"): 0.02, ("COMBO_EH", "H"): 0.005, ("DOUBLE", "E"): 0.01, ("DOUBLE", "F"): 0.025, ("HET_COPY", "E"): 0.02,
           ("HET_COPY", "H"): 0.002, ("MANY_REMOVED", "E"): 0.01, ("CONV_ABSENT", "E"): 0.02, ("MOSAIC", "E"): 0.02,
           ("IMBALANCE", "E"): 0.02, ("ISOFORMS", "E"): 0.02, ("STRUCT_ABSENT", "E"): 0.02, ("EDIT_COPY", "E"): 0.02,
           ("LOW4", "E"): 0.02, ("SPLIT", "E"): 0.01}


def tdiv(scenario, source):
    kind = scenario.split("@")[0]
    if (kind, source[-1]) in FIXED_D:
        return FIXED_D[(kind, source[-1])]
    return float(scenario.split("@")[1])


def cmd_rows(a):
    w, d = a.world, f"{a.world}/{a.variant}.{a.arm}"
    scen = json.load(open(f"{w}/world.json"))
    kinds = json.load(open(f"{w}/kinds.json"))
    if a.rust:
        # the Rust product (`copy_assign --psv-clusters --out D/R`), its members turned into the prototype's `sources` column
        clusters = list(csv.DictReader(open(f"{d}/R.psv_clusters.tsv"), delimiter="\t"))
        srcs = collections.defaultdict(collections.Counter)
        for r in csv.DictReader(open(f"{d}/R.psv_cluster_reads.tsv"), delimiter="\t"):
            srcs[(r["family"], r["host"], r["cluster"])][r["read_name"].split("|")[0]] += 1
        for c in clusters:
            c["sources"] = ",".join(f"{k}:{v}" for k, v in srcs[(c["family"], c["host"], c["cluster"])].most_common())
    else:
        clusters = list(csv.DictReader(open(f"{d}/psvclust.tsv"), delimiter="\t"))
    placed = collections.Counter()
    host_of = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(f"{d}/A.assignments.tsv"), delimiter="\t"):
        src = r["read_name"].split("|")[0]
        placed[src] += 1
        host_of[src][r["catalog_copy_idx"]] += 1
    n_reads = collections.Counter(ln.split("\t")[1] for ln in list(open(f"{w}/reads.truth.{a.variant}.tsv"))[1:])
    flags = collections.defaultdict(dict)
    for r in csv.DictReader(open(f"{d}/A.family_join.tsv"), delimiter="\t"):
        flags[r["family_id"]][r["copy_index"]] = r["o3_flag"]
    out = []
    base = dict(world=os.path.basename(w.rstrip("/")), variant=a.variant, arm=a.arm)
    spurious = collections.Counter()     # unexplained clusters whose majority source is a PRESENT copy (any family)
    for c in clusters:
        if c["status"] == "unexplained":
            top = max((x.rsplit(":", 1) for x in c["sources"].split(",")), key=lambda x: int(x[1]))[0]
            if kinds.get(top) == "present":
                spurious[c["family"]] += 1
    for fam, sc in sorted(scen.items()):
        kind = sc.split("@")[0]
        fam_cl = [c for c in clusters if c["family"] == fam]
        if kind in CONTROL:
            fl = [c for c in fam_cl if c["status"] == "unexplained"]
            out.append(dict(base, unit="control", family=fam, scenario=sc, source="-", truth_kind="-", n_source=0, n_placed=0,
                            recovered=0, purity=0, recall=0, est_div=0, label="-", false_flag=int(bool(fl)), spurious=spurious[fam],
                            baseline=int(any(v == "missing_copy" for v in flags[fam].values()))))
            continue
        for src, k in kinds.items():
            if not src.startswith(fam) or k == "present":
                continue
            best = (0, None)
            for c in fam_cl:
                if c["status"] != "unexplained":
                    continue
                cnt = dict(x.rsplit(":", 1) for x in c["sources"].split(","))
                got = int(cnt.get(src, 0))
                if got > best[0]:
                    best = (got, c)
            got, c = best
            recall = got / placed[src] if placed[src] else 0.0
            purity = got / int(c["n_reads"]) if c else 0.0
            rec = int(c is not None and recall >= RECALL and purity >= PURITY)
            est = float(c["est_div"]) if c else 0.0
            host = host_of[src].most_common(1)[0][0] if host_of[src] else "-"
            out.append(dict(base, unit="source", family=fam, scenario=sc, source=src, truth_kind=k, n_source=n_reads[src],
                            n_placed=placed[src], recovered=rec, purity=round(purity, 3), recall=round(recall, 3), est_div=est,
                            label=("copy" if est > DELTA else "allele_or_copy") if rec else "-", false_flag=0,
                            spurious=spurious[fam] if src == min(x for x, kk in kinds.items() if x.startswith(fam) and kk != "present") else 0,
                            baseline=int(flags[fam].get(host) == "missing_copy")))
    new = not os.path.exists(a.rows)
    with open(a.rows, "a") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        if new:
            wr.writeheader()
        wr.writerows(out)
    src = [r for r in out if r["unit"] == "source"]
    ctl = [r for r in out if r["unit"] == "control"]
    print(f"{base['world']}\t{a.variant}\t{a.arm}\tsources recovered {sum(r['recovered'] for r in src)}/{len(src)}"
          f"\tbaseline {sum(r['baseline'] for r in src)}/{len(src)}\tcontrol false flags {sum(r['false_flag'] for r in ctl)}/{len(ctl)}"
          f" (baseline {sum(r['baseline'] for r in ctl)})")


def div_of(sc):
    return float(sc.split("@")[1]) if "@" in sc else {"LOWEXPR": 0.02}.get(sc, None)


def cmd_summary(a):
    rows = list(csv.DictReader(open(a.rows), delimiter="\t"))
    for r in rows:
        for k in ("recovered", "false_flag", "baseline", "n_placed"):
            r[k] = int(r[k])
    frac = lambda k, n: f"{k}/{n} = {k / n:.3f}" if n else "0/0"
    for arm in sorted({r["arm"] for r in rows}):
        R = [r for r in rows if r["arm"] == arm]
        S = [r for r in R if r["unit"] == "source"]
        C = [r for r in R if r["unit"] == "control"]
        print(f"===== arm {arm}")
        td = lambda r: tdiv(r["scenario"], r["source"])
        copies = [r for r in S if r["truth_kind"] == "absent"]
        d1 = [r for r in copies if td(r) >= 0.006]
        print(f"D1 absent copies >= 0.6%, recovered (pure unexplained cluster): {frac(sum(r['recovered'] for r in d1), len(d1))}"
              f"   baseline detector {frac(sum(r['baseline'] for r in d1), len(d1))}")
        print(f"D2 control families with an unexplained cluster: {frac(sum(r['false_flag'] for r in C), len(C))}"
              f"   baseline {frac(sum(r['baseline'] for r in C), len(C))}")
        multi = collections.defaultdict(list)
        for r in S:
            if r["scenario"] in MULTI_GROUP:
                multi[(r["world"], r["variant"], r["family"])].append(r["recovered"])
        both = sum(all(v) for v in multi.values())
        print(f"D3 two or more groups on one host ({', '.join(MULTI_GROUP)}): all recovered separately in {frac(both, len(multi))}")
        print(f"spurious unexplained clusters (majority a PRESENT copy's reads), all families: {sum(int(r['spurious']) for r in R)}")
        al = [r for r in S if r["truth_kind"] == "allele" and r["recovered"]]
        cp = [r for r in copies if td(r) >= 0.02 and r["recovered"]]
        print(f"D4 labels: recovered alleles labelled allele_or_copy {frac(sum(r['label'] == 'allele_or_copy' for r in al), len(al))};"
              f" recovered copies >= 2% labelled copy {frac(sum(r['label'] == 'copy' for r in cp), len(cp))}")
        print("-- per scenario (recovered / units; baseline)")
        by = collections.defaultdict(list)
        for r in S:
            by[(r["scenario"].split("@")[0], td(r), r["truth_kind"])].append(r)
        for k in sorted(by):
            x = by[k]
            print(f"  {k[0]:14s} d={k[1]:<6} {k[2]:7s} {sum(r['recovered'] for r in x):3d}/{len(x):<3d} baseline {sum(r['baseline'] for r in x):3d}"
                  f"  labels {dict(collections.Counter(r['label'] for r in x if r['recovered']))}")
        for k in sorted({r["scenario"] for r in C}):
            x = [r for r in C if r["scenario"] == k]
            print(f"  control {k:12s} false flags {sum(r['false_flag'] for r in x)}/{len(x)}  baseline {sum(r['baseline'] for r in x)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("rows"); p.add_argument("world"); p.add_argument("variant"); p.add_argument("arm"); p.add_argument("rows")
    p.add_argument("--rust", action="store_true", help="score the Rust product D/R.psv_clusters.tsv instead of D/psvclust.tsv")
    p.set_defaults(func=cmd_rows)
    p = sub.add_parser("summary"); p.add_argument("rows"); p.set_defaults(func=cmd_summary)
    a = ap.parse_args()
    sys.exit(a.func(a))
