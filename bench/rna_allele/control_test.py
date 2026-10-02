#!/usr/bin/env python3
"""Amendment 9 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): the no-deletion control of the IsoCon chain on Amendment 7's 53
families against the unmasked `_pri`. Reuses the linktest work dir (--l: panel.json, scored.fa, labels.tsv) and merge_test's merge.

  net        IsoCon input per family from R0.bam (a record on ANY copy of the family, or unmapped), <= 1,000 reads (seed 1) -> fam/<fam>.fa
  outputs    iso/<fam>/final_candidates.fa -> outputs.fa (names <fam>|<id>)
  contigs    flag (identity x coverage < 0.999 vs `_pri`, outputs.pri.paf) -> link (d <= delta) -> contigs.tsv (+ source), contigs_I/L.fa;
             then merge_test pairs + components
  lift       every copy interval -> its B-haplotype interval through the frozen truth's asm5 alignments (truth_lift.py)
  classify   candidates (components) x haplotype hits (contigs_L.mat.paf / .pat.paf) -> a (haplotype-only locus) / b (allele) / c (unmatched)
             / pri (in `_pri` after all); rule C1; overlap with the deletion run's survivor-derived candidates (overlap.paf)
  score      arms R0 and C (components as loci): stay / own candidate / false move / elsewhere / unplaced; rule C2

    control_test.py net --w /mnt/linuxdisk/tmp/rna_allele/control --l /mnt/linuxdisk/tmp/rna_allele/linktest
"""
import argparse
import collections
import csv
import glob
import json
import os
import random
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_test  # noqa: E402
import truth_lift  # noqa: E402

DELTA = merge_test.DELTA
TRUTH = "/mnt/linuxdisk/tmp/rna_allele"


def panel(a):
    P = json.load(open(f"{a.l}/panel.json"))
    return {p["fam"]: [p["mask"]] + p["keep"] for p in P}          # every copy: [chrom, start, end, gene]


def labels(a):
    return {r["read"]: r for r in csv.DictReader(open(f"{a.l}/labels.tsv"), delimiter="\t")}


def net(a):
    P, lab = panel(a), labels(a)
    rec = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(f"{a.w}/R0.bam").fetch(until_eof=True):
        rec[rd.query_name].append(None if rd.is_unmapped else (rd.reference_name, rd.reference_start, rd.reference_end))
    seq, cur = {}, None
    for ln in open(f"{a.l}/scored.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip()
        else:
            seq[cur] = ln.strip()
    os.makedirs(f"{a.w}/fam", exist_ok=True)
    rng = random.Random(1)
    by = collections.defaultdict(list)
    for n, r in lab.items():
        rs = rec.get(n, [])
        unm = bool(rs) and all(x is None for x in rs)
        on = any(x and x[0] == c and x[1] < e and s < x[2] for x in rs for c, s, e, g in P[r["family"]])
        if unm or on:
            by[r["family"]].append(n)
    sizes = []
    for f, ns in by.items():
        ns = sorted(ns); rng.shuffle(ns); ns = ns[:1000]; sizes.append(len(ns))
        with open(f"{a.w}/fam/{f}.fa", "w") as o:
            for n in ns:
                o.write(f">{n}\n{seq[n]}\n")
    print(f"IsoCon inputs: {len(by)} families, total {sum(sizes)}, median {sorted(sizes)[len(sizes) // 2]}, max {max(sizes)}")


def outputs(a):
    n = 0
    with open(f"{a.w}/outputs.fa", "w") as o:
        for path in sorted(glob.glob(f"{a.w}/iso/*/final_candidates.fa")):
            fam = path.split("/")[-2]
            for ln in open(path):
                if ln[0] == ">":
                    o.write(f">{fam}|{ln[1:].strip()}\n"); n += 1
                else:
                    o.write(ln)
    print("outputs", n, "families", len(glob.glob(f"{a.w}/iso/*/final_candidates.fa")))


def best_hits(path):
    b = {}
    for ln in open(path):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, f[5], int(f[7]), int(f[8]), int(f[9]), int(f[1]))
    return b


def contigs(a):
    P = panel(a)
    B = best_hits(f"{a.w}/outputs.pri.paf")
    seqs, cur = {}, None
    for ln in open(f"{a.w}/outputs.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip(); seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    nI = nL = 0
    k = collections.Counter()
    with open(f"{a.w}/contigs_I.fa", "w") as fi, open(f"{a.w}/contigs_L.fa", "w") as fl, open(f"{a.w}/contigs.tsv", "w") as t:
        t.write("contig\tfamily\toutput\tlength\tbest_ref\td\tlinked\tsource\n")
        for o, sl in seqs.items():
            b = B.get(o)
            if b and b[0] >= 0.999:
                continue
            fam = o.split("|")[0]
            s = "".join(sl)
            d = 1 - (b[4] / len(s) if b else 0.0)
            linked = d <= DELTA
            src = "none"
            if b:
                src = next((f"S:{g}" for c, s0, e, g in P[fam] if b[1] == c and b[2] < e and s0 < b[3]), "elsewhere")
            ctg = f"iso_{fam}_{k[fam]}"; k[fam] += 1
            fi.write(f">{ctg}\n{s}\n"); nI += 1
            if not linked:
                fl.write(f">{ctg}\n{s}\n"); nL += 1
            t.write(f"{ctg}\t{fam}\t{o}\t{len(s)}\t{b[0] if b else 0:.4f}\t{d:.5f}\t{int(linked)}\t{src}\n")
    rows = list(csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t"))
    print(f"outputs {len(seqs)}; flagged {nI}; linked back {nI - nL}; kept as new copies {nL} in "
          f"{len({r['family'] for r in rows if r['linked'] == '0'})} families; new-copy sources",
          dict(collections.Counter(r["source"].split(":")[0] for r in rows if r["linked"] == "0")))
    if nL:
        merge_test.pairs(a); merge_test.components(a)


def lift(a):
    P = panel(a)
    with open(f"{a.w}/copies_genes.tsv", "w") as o:
        o.write("gene_id\tchrom\tstrand\texons\n")
        for fam, cps in P.items():
            for c, s, e, g in cps:
                o.write(f"{g}\t{c}\t+\t{s}-{e}\n")
    truth_lift.main(["--chrmap", f"{TRUTH}/chrmap.tsv", "--paf-dir", f"{TRUTH}/out", "--genes", f"{a.w}/copies_genes.tsv",
                     "--out", f"{a.w}/copies_lift.tsv"])
    rows = list(csv.DictReader(open(f"{a.w}/copies_lift.tsv"), delimiter="\t"))
    print(f"copies {len(rows)}; lifted >= 50%: {sum(1 for r in rows if float(r['lift_frac']) >= 0.5)}; "
          f"classes {dict(collections.Counter(r['class'] for r in rows))}")


def hap_alias():
    """chrN_<hap>_hsa* (index names) -> (hap, accession)"""
    alias = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            alias[(h, num)] = acc
    return alias


def comps_at(a, delta):
    rows = merge_test.rows_new(a.w)
    fam_of = {r["contig"]: r["family"] for r in rows}
    best = merge_test.best_pairs(a.w, sorted(set(fam_of.values())))
    par = {c: c for c in fam_of}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for (p, q), (m, cov, de) in best.items():
        if cov >= 0.5 and de <= delta:
            par[find(p)] = find(q)
    return {c: f"{fam_of[c]}:{find(c)}" for c in fam_of}, fam_of


def classify(a):
    P = panel(a)
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t") if r["linked"] == "0"}
    cm = {r["pri"]: r for r in csv.DictReader(open(f"{TRUTH}/chrmap.tsv"), delimiter="\t")}
    alias = hap_alias()
    bside = {(r["B_hap"], r["B_name"]) for r in cm.values() if r["B_name"]}
    liftd = {r["gene_id"]: r for r in csv.DictReader(open(f"{a.w}/copies_lift.tsv"), delimiter="\t")}
    # best haplotype hit per contig: (score, hap, acc, start, end)
    hb = {}
    for h in ("mat", "pat"):
        for ln in open(f"{a.w}/contigs_L.{h}.paf"):
            f = ln.rstrip("\n").split("\t")
            m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]*", f[5])
            if not m:
                continue
            acc = alias[(m.group(2), m.group(1))]
            s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
            if f[0] not in hb or s > hb[f[0]][0]:
                hb[f[0]] = (s, h, acc, int(f[7]), int(f[8]))

    def cls_contig(c):
        b = hb.get(c)
        if not b or b[0] < 0.999:
            return "c_unmatched"
        _, h, acc, s, e = b
        if (h, acc) not in bside:
            return "pri"
        fam = ctg[c]["family"]
        for cc, s0, e0, g in P[fam]:
            L = liftd.get(g)
            if L and float(L["lift_frac"]) >= 0.5 and L["B_chrom"] == acc and int(L["B_start"]) < e and s < int(L["B_end"]):
                return "b_allele"
        return "a_haplotype_only"
    # overlap with the deletion run's survivor-derived candidates (asm20 alignment, identity >= 0.999 over >= 90% of the shorter)
    ov = set()
    if os.path.exists(f"{a.w}/overlap.paf"):
        for ln in open(f"{a.w}/overlap.paf"):
            f = ln.split("\t")
            ident = int(f[9]) / max(1, int(f[10]))
            short = min(int(f[1]), int(f[6]))
            cov = (int(f[3]) - int(f[2])) / short if int(f[1]) <= int(f[6]) else (int(f[8]) - int(f[7])) / short
            if ident >= 0.999 and cov >= 0.9 and f[5].split("_")[1] == f[0].split("_")[1]:   # same family (iso_<fam>_<k>)
                ov.add(f[0])
    out = {}
    for name, delta in merge_test.DELTAS.items():
        comp, fam_of = comps_at(a, delta)
        members = collections.defaultdict(list)
        for c, cid in comp.items():
            members[cid].append(c)
        cls_comp, fam_false, fam_any = {}, set(), set()
        for cid, cs in members.items():
            best = max(cs, key=lambda c: hb.get(c, (0,))[0])
            k = cls_contig(best)
            cls_comp[cid] = k
            fam = fam_of[cs[0]]
            fam_any.add(fam)
            if k != "a_haplotype_only":
                fam_false.add(fam)
        out[name] = dict(delta=delta, candidates=len(members), classes=dict(collections.Counter(cls_comp.values())),
                         families_any=len(fam_any), families_false=len(fam_false),
                         candidates_matching_deletion_run=sum(1 for cid, cs in members.items() if any(c in ov for c in cs)))
        if name == "delta":
            with open(f"{a.w}/candidates.tsv", "w") as o:
                o.write("candidate\tfamily\tn_contigs\tclass\tbest_contig\tsource\thap_hit\tmatches_deletion_run\n")
                for cid, cs in sorted(members.items()):
                    best = max(cs, key=lambda c: hb.get(c, (0,))[0])
                    b = hb.get(best)
                    o.write(f"{cid}\t{fam_of[cs[0]]}\t{len(cs)}\t{cls_comp[cid]}\t{best}\t{ctg[best]['source']}\t"
                            f"{(b[1] + ':' + b[2] + ':' + str(b[3]) + '-' + str(b[4]) + ' ' + format(b[0], '.4f')) if b else ''}\t"
                            f"{int(any(c in ov for c in cs))}\n")
            per_fam = collections.Counter(fam_of[cs[0]] for cs in members.values())
            out["per_family_candidates"] = dict(per_fam)
        print(f"[{name} delta={delta:.5f}] candidates {len(members)} {dict(collections.Counter(cls_comp.values()))}; families with >= 1 "
              f"candidate {len(fam_any)}/53, with >= 1 FALSE candidate {len(fam_false)}/53 = {len(fam_false) / 53:.1%}; "
              f"matching the deletion run's survivor candidates {out[name]['candidates_matching_deletion_run']}")
    r = out["delta"]
    print(f"C1: families with >= 1 false candidate {r['families_false']}/53 = {r['families_false'] / 53:.1%} (bar <= 28%) -> "
          f"{'HOLDS' if r['families_false'] / 53 <= 0.28 else 'FAILS'}; any candidate {r['families_any']}/53 = {r['families_any'] / 53:.1%}")
    json.dump(out, open(f"{a.w}/classify.json", "w"), indent=0)


def score(a):
    P, lab = panel(a), labels(a)
    comp, fam_of = comps_at(a, DELTA)
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t")}
    holds = collections.defaultdict(set)
    for c, cid in comp.items():
        holds[cid].add(ctg[c]["source"])
    want = set(lab)

    def locus(chrom, s, e, fam):
        if chrom.startswith("iso_"):
            return ("ctg", comp[chrom])
        for c, s0, e0, g in P[fam]:
            if chrom == c and s < e0 and s0 < e:
                return ("copy", g)
        return ("other", f"{chrom}:{s // 100000}")

    def calls(bam):
        recs = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name not in want or rd.is_supplementary:
                continue
            recs[rd.query_name].append(None if rd.is_unmapped else
                                       (rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                        rd.get_tag("AS") if rd.has_tag("AS") else 0))
        out = {}
        for n, rs in recs.items():
            rs = [r for r in rs if r]
            if not rs:
                out[n] = ("unplaced", None); continue
            fam = lab[n]["family"]
            prim = next((r for r in rs if not r[0]), rs[0])
            srt = sorted(rs, key=lambda r: -r[4])
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= 0.98 * srt[0][4]:
                if len({locus(r[1], r[2], r[3], fam) for r in srt if r[4] >= 0.98 * srt[0][4]}) > 1:
                    out[n] = ("unplaced", None); continue
            out[n] = ("placed", locus(prim[1], prim[2], prim[3], fam))
        return out

    def cls(n, call):
        r = lab[n]
        st, L = call
        if st == "unplaced":
            return "unplaced"
        kind, key = L
        if kind == "copy":
            return "stay" if key == r["copy"] else "other_copy"
        if kind == "ctg":
            return "own_candidate" if ("S:" + r["copy"]) in holds[key] else "false_move"
        return "elsewhere"
    res = {}
    for arm in ("R0", "C"):
        if not os.path.exists(f"{a.w}/{arm}.bam"):
            continue
        Cc = calls(f"{a.w}/{arm}.bam")
        res[arm] = collections.Counter(cls(n, Cc.get(n, ("unplaced", None))) for n in lab)
        print(f"[{arm}]", dict(sorted(res[arm].items())))
    if "C" in res:
        fm = res["C"]["false_move"] / len(lab)
        print(f"C2: false moves {res['C']['false_move']}/{len(lab)} = {fm:.2%} (bar <= 5%) -> {'HOLDS' if fm <= 0.05 else 'FAILS'}; "
              f"own-candidate placements {res['C']['own_candidate']}; unplaced {res['R0']['unplaced']} -> {res['C']['unplaced']}")
    json.dump({k: dict(v) for k, v in res.items()}, open(f"{a.w}/score.json", "w"), indent=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["net", "outputs", "contigs", "lift", "classify", "score"])
    ap.add_argument("--w", required=True)
    ap.add_argument("--l", default="/mnt/linuxdisk/tmp/rna_allele/linktest")
    a = ap.parse_args(argv)
    os.makedirs(a.w, exist_ok=True)
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
