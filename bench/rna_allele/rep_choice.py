#!/usr/bin/env python3
"""Representative per candidate component (design for the `candidates` stage): on Amendment 8's 53 families (linktest work dir), where
every read is aligned to every member contig of every component (RIL.bam), measure what each rule loses when only the representative
is kept as the O2 copy.

Rules: LONGEST member; SUPPORT = member with the most IsoCon read support; FIRST = IsoCon's first output (its own ordering); UNION =
ceiling (a read keeps its best AS over all members, i.e. an exon-union representative at best).
Per read with a record on a component: AS on the representative vs its best AS over the component's members. Loss = reads whose AS
on the representative falls below 0.98 x best (the tie rule would then prefer a reference locus or abstain) or that have no record on it.
D reads: reads of the deleted copy on a D-derived component; S reads: survivors' reads on components derived from their own copy.
Also: how often the longest / most-supported member looks chimeric (best unmasked-genome hit covers < 80% of it).

    rep_choice.py --w /mnt/linuxdisk/tmp/rna_allele/linktest
"""
import argparse
import collections
import csv
import re

import pysam


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    W = a.w
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{W}/contigs.tsv"), delimiter="\t")}
    comp = {r["contig"]: r["component"] for r in csv.DictReader(open(f"{W}/merge/components.tsv"), delimiter="\t")}
    members = collections.defaultdict(list)
    for c, cid in comp.items():
        members[cid].append(c)
    multi = {cid: cs for cid, cs in members.items() if len(cs) >= 2}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{W}/labels.tsv"), delimiter="\t")}
    sup = {c: int(re.search(r"support_(\d+)", ctg[c]["output"]).group(1)) for c in comp}
    order = {c: int(ctg[c]["contig"].rsplit("_", 1)[1]) for c in comp}       # iso_<fam>_<k>: k follows IsoCon's output order
    # unmasked-genome coverage of each contig (chimera check)
    cov = {}
    for ln in open(f"{W}/outputs.base.paf"):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        c = (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in cov or s > cov[f[0]][0]:
            cov[f[0]] = (s, c)
    out2ctg = {ctg[c]["output"]: c for c in comp}
    ccov = {out2ctg[o]: v[1] for o, v in cov.items() if o in out2ctg}
    reps = {
        "LONGEST": {cid: max(cs, key=lambda c: int(ctg[c]["length"])) for cid, cs in multi.items()},
        "SUPPORT": {cid: max(cs, key=lambda c: (sup[c], int(ctg[c]["length"]))) for cid, cs in multi.items()},
        "FIRST": {cid: min(cs, key=lambda c: order[c]) for cid, cs in multi.items()},
    }
    # per read: AS on each member contig of a multi-member component
    AS = collections.defaultdict(dict)        # read -> contig -> best AS
    for rd in pysam.AlignmentFile(f"{W}/RIL.bam").fetch(until_eof=True):
        if rd.is_unmapped or rd.is_supplementary or not rd.reference_name.startswith("iso_"):
            continue
        c = rd.reference_name
        if comp.get(c) not in multi:
            continue
        s = rd.get_tag("AS") if rd.has_tag("AS") else 0
        if s > AS[rd.query_name].get(c, -1):
            AS[rd.query_name][c] = s
    # evaluate
    print(f"components with >= 2 members: {len(multi)} (of {len(members)}); members total {sum(len(v) for v in multi.values())}")
    for name, R in reps.items():
        res = {"D": collections.Counter(), "S": collections.Counter()}
        for n, per in AS.items():
            r = lab.get(n)
            if not r:
                continue
            # the component this read belongs to by truth: D read -> D-derived component of its family; S read -> component holding its own copy's contig
            for cid, cs in multi.items():
                if ctg[cs[0]]["family"] != r["family"]:
                    continue
                srcs = {ctg[c]["source"] for c in cs}
                own = ("D" in srcs) if r["role"] == "D" else (("S:" + r["copy"]) in srcs)
                if not own:
                    continue
                best = max((per.get(c, -1) for c in cs), default=-1)
                if best < 0:
                    continue
                rep_as = per.get(R[cid], -1)
                k = "kept" if rep_as >= 0.98 * best else ("no_record" if rep_as < 0 else "lost")
                res[r["role"]][k] += 1
        chim = sum(1 for cid in multi if ccov.get(R[cid], 1.0) < 0.8)
        d, s = res["D"], res["S"]
        print(f"{name:8s} D reads: kept {d['kept']} lost {d['lost']} no_record {d['no_record']} (kept {d['kept'] / max(1, sum(d.values())):.1%}) | "
              f"S reads: kept {s['kept']} lost {s['lost']} no_record {s['no_record']} (kept {s['kept'] / max(1, sum(s.values())):.1%}) | "
              f"representatives with genome coverage < 0.8: {chim}/{len(multi)}")
    # the UNION ceiling is 100% kept by construction; show how much length a union would add over the longest
    extra = []
    for cid, cs in multi.items():
        L = max(int(ctg[c]["length"]) for c in cs)
        extra.append(sum(int(ctg[c]["length"]) for c in cs) / L)
    extra.sort()
    print(f"UNION    ceiling: every read kept; total member length / longest member: median {extra[len(extra) // 2]:.2f}, max {extra[-1]:.2f}")
    # agreement between rules
    same = sum(1 for cid in multi if reps["LONGEST"][cid] == reps["SUPPORT"][cid])
    print(f"LONGEST == SUPPORT in {same}/{len(multi)} components; member count per component median "
          f"{sorted(len(v) for v in multi.values())[len(multi) // 2]}, max {max(len(v) for v in multi.values())}")


if __name__ == "__main__":
    main()
