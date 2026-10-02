#!/usr/bin/env python3
"""Amendment 6 scoring (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md).

Steps (subcommands):
  contigs   collect every family's IsoCon outputs (iso/<fam>/final_candidates.fa) as `<fam>|<id>`; given their PAFs against the masked
            genome and the unmasked genome, flag outputs with no hit at identity x coverage >= 0.999 in the masked genome and write them as
            contigs `iso_<fam>_<k>` (+ a table with is_D / is_K from the unmasked hits).
  score     final call per scored read in arm R (o3_excise/masked.bam) and arm R+I (realigned BAM): primary locus, or "abstain" when the
            read's best and second-best AS over its records tie (second >= 0.98 x best); then the registered readouts and rule.

    heldout_score.py contigs --h heldout --panel panel.json --paf-masked outs.masked.paf --paf-base outs.base.paf
    heldout_score.py score --h heldout --panel panel.json --r-bam o3_excise/masked.bam --ri-bam heldout/ri.bam
"""
import argparse
import collections
import csv
import glob
import json
import os
import statistics

import pysam


def best_hits(path):
    b = {}
    for ln in open(path):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, f[5], int(f[7]), int(f[8]))
    return b


def contigs(a):
    panel = {p["fam"]: p for p in json.load(open(a.panel))}
    M, B = best_hits(a.paf_masked), best_hits(a.paf_base)
    n_out = 0
    with open(f"{a.h}/contigs.fa", "w") as fa, open(f"{a.h}/contigs.tsv", "w") as t:
        t.write("contig\tfamily\toutput\tlength\tbest_masked\tis_D\tis_K\n")
        for path in sorted(glob.glob(f"{a.h}/iso/*/final_candidates.fa")):
            fam = path.split("/")[-2]
            seqs, cur = {}, None
            for ln in open(path):
                if ln[0] == ">":
                    cur = f"{fam}|{ln[1:].strip()}"; seqs[cur] = []
                else:
                    seqs[cur].append(ln.strip())
            k = 0
            for o, sl in seqs.items():
                n_out += 1
                if M.get(o, (0,))[0] >= 0.999:
                    continue
                b = B.get(o)
                def on(lab):
                    c, s, e = panel[fam][lab]
                    return int(bool(b) and b[0] >= 0.999 and b[1] == c and b[2] < e and s < b[3])
                ctg = f"iso_{fam}_{k}"; k += 1
                s = "".join(sl)
                fa.write(f">{ctg}\n{s}\n")
                t.write(f"{ctg}\t{fam}\t{o}\t{len(s)}\t{M.get(o, (0,))[0]:.4f}\t{on('mask')}\t{on('keep')}\n")
    rows = list(csv.DictReader(open(f"{a.h}/contigs.tsv"), delimiter="\t"))
    print(f"IsoCon outputs {n_out}; flagged (not in the masked genome at 0.999) {len(rows)}; of which is_D "
          f"{sum(int(r['is_D']) for r in rows)}, is_K {sum(int(r['is_K']) for r in rows)}; families with >= 1 is_D contig "
          f"{len({r['family'] for r in rows if r['is_D'] == '1'})}")


def calls(bam, want):
    recs = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
        if rd.query_name not in want or rd.is_supplementary:
            continue
        if rd.is_unmapped:
            recs[rd.query_name].append(None)
            continue
        recs[rd.query_name].append((rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                    rd.get_tag("AS") if rd.has_tag("AS") else 0, rd.get_tag("de") if rd.has_tag("de") else None))
    out = {}
    for n, rs in recs.items():
        rs = [r for r in rs if r]
        if not rs:
            out[n] = ("unmapped", None, None)
            continue
        prim = next((r for r in rs if not r[0]), None)
        a = sorted((r[4] for r in rs), reverse=True)
        if len(a) > 1 and a[1] > 0 and a[1] >= 0.98 * a[0]:
            out[n] = ("abstain", prim, None)
        else:
            out[n] = ("placed", prim, prim[5] if prim else None)
    return out


def score(a):
    panel = {p["fam"]: p for p in json.load(open(a.panel))}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.h}/labels.tsv"), delimiter="\t")}
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.h}/contigs.tsv"), delimiter="\t")}
    want = set(lab)
    C = {"R": calls(a.r_bam, want), "RI": calls(a.ri_bam, want)}

    def classify(arm, n):
        f, role = lab[n]["family"], lab[n]["role"]
        st, prim, de = C[arm].get(n, ("unmapped", None, None))
        if st in ("unmapped", "abstain"):
            return "unplaced"
        chrom, s, e = prim[1], prim[2], prim[3]
        if chrom.startswith("iso_"):
            r = ctg[chrom]
            if role == "D":
                return "right" if (r["family"] == f and r["is_D"] == "1") else "wrong"
            return "stay" if (r["family"] == f and r["is_K"] == "1") else "false_move"
        kc, ks, ke = panel[f]["keep"]
        onk = chrom == kc and s < ke and ks < e
        if role == "D":
            return "wrong"
        return "stay" if onk else "elsewhere"
    res = {arm: collections.Counter() for arm in C}
    per = collections.defaultdict(lambda: {arm: collections.Counter() for arm in C})
    for n, r in lab.items():
        for arm in C:
            k = classify(arm, n)
            res[arm][(r["role"], k)] += 1
            per[r["family"]][arm][(r["role"], k)] += 1
    nD = sum(1 for r in lab.values() if r["role"] == "D"); nK = len(lab) - nD
    for arm in ("R", "RI"):
        d = {k[1]: v for k, v in res[arm].items() if k[0] == "D"}
        kk = {k[1]: v for k, v in res[arm].items() if k[0] == "K"}
        print(f"[{arm}] D reads {nD}: {d}\n[{arm}] K reads {nK}: {kk}")
    w0, w1 = res["R"][("D", "wrong")], res["RI"][("D", "wrong")]
    fm = res["RI"][("K", "false_move")] / nK
    drop = (w0 - w1) / w0 if w0 else 0.0
    v = "HELP" if drop >= 0.5 and fm <= 0.05 else "HURT" if fm > 0.10 else "MIXED"
    print(f"RULE: wrong D {w0} -> {w1} (drop {drop:.1%}); right D {res['R'][('D','right')]} -> {res['RI'][('D','right')]}; "
          f"false moves {res['RI'][('K','false_move')]}/{nK} = {fm:.1%} -> {v}")
    # strata: fate (2026-08-14) and divergence (median de of D reads on K in the masked arm)
    fate = {}
    for x in json.load(open(a.per_family)):
        fate[x["fam"]] = "orphaned" if x.get("unaln", 0) >= 0.5 else "absorbed" if x.get("conc", 0) >= 0.5 else "scattered"
    div = {}
    for f in panel:
        kc, ks, ke = panel[f]["keep"]
        des = [C["R"][n][2] for n, r in lab.items() if r["family"] == f and r["role"] == "D" and C["R"].get(n, ("",))[0] == "placed"
               and C["R"][n][1] and C["R"][n][1][1] == kc and C["R"][n][1][2] < ke and ks < C["R"][n][1][3] and C["R"][n][2] is not None]
        div[f] = "D-K >= 0.01" if des and statistics.median(des) >= 0.01 else ("D-K < 0.01" if des else "no D read on K")
    for name, strat in (("fate", fate), ("divergence", div)):
        groups = collections.defaultdict(lambda: collections.Counter())
        for f, pc in per.items():
            g = strat.get(f, "?")
            for arm in ("R", "RI"):
                for k, v in pc[arm].items():
                    groups[g][(arm,) + k] += v
        for g, c in sorted(groups.items()):
            print(f"  [{name}={g}] families {sum(1 for f in per if strat.get(f, '?') == g)}: wrong D {c[('R','D','wrong')]} -> "
                  f"{c[('RI','D','wrong')]}, right D {c[('R','D','right')]} -> {c[('RI','D','right')]}, unplaced D "
                  f"{c[('R','D','unplaced')]} -> {c[('RI','D','unplaced')]}, false moves {c[('RI','K','false_move')]}")
    helped = sum(1 for f, pc in per.items() if pc["R"][("D", "wrong")] and
                 (pc["R"][("D", "wrong")] - pc["RI"][("D", "wrong")]) / pc["R"][("D", "wrong")] >= 0.5)
    print(f"families where wrong D fell by >= 50%: {helped} of {sum(1 for f, pc in per.items() if pc['R'][('D','wrong')])} with wrong D in R")
    json.dump({f: {arm: {"|".join(k): v for k, v in pc[arm].items()} for arm in pc} for f, pc in per.items()},
              open(f"{a.h}/per_family_result.json", "w"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("contigs")
    for k in ("h", "panel", "paf_masked", "paf_base"):
        c.add_argument("--" + k.replace("_", "-"), required=True)
    s = sub.add_parser("score")
    for k in ("h", "panel", "r_bam", "ri_bam", "per_family"):
        s.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    contigs(a) if a.cmd == "contigs" else score(a)


if __name__ == "__main__":
    main()
