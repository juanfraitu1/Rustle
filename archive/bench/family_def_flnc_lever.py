#!/usr/bin/env python3
"""family_def_flnc_lever.py — leverage the FULL-LENGTH (FLNC) single-molecule property.

The reads are FLNC = complete transcript molecules, aligned with -Y (secondaries SOFT-clipped,
so the full read is retained at every placement). New levers the prior features missed:
  Q1 qaln_frac = (query-aligned bases)/(full read length) per placement. A real paralog's
     molecule aligns FULL-LENGTH at both copies (qaln_frac~1 at both); a bridge's molecule
     aligns full at its own locus but only the shared FRAGMENT at the other (low qaln_frac).
       qaln_min = median over crossmap reads of min(qaln_A, qaln_B)
       qaln_tie = median |qaln_A - qaln_B|
  Q2 cov_frac = fraction of the GENE covered by crossmap-read alignment blocks (real paralog:
     whole transcript covered at both copies; bridge: only the shared domain). cov_frac_min.
Tested on the 85 ~B-survivors (76 true DNA-paralog / 9 false), vs DNA label + permutation null,
to see if the FLNC signal separates the residual that recip-cov / mask / exonic could not.
Run: python bench/family_def_flnc_lever.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import re
import subprocess
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED, best_gene, ref_span, DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import build_model, recip_cov, GENOME, NEW, COV_MIN
from family_def_feature_discovery import auc, bimodality
from family_def_read_filters import dna_homology

HERE = os.path.dirname(os.path.abspath(__file__))
SAM = "/home/juanfra/miniforge3/bin/samtools"
CIG = re.compile(r"(\d+)([=XIDNMSH])")


def placement(f):
    """(de, qaln_frac, ref_start, blocks) from an --eqx -Y SAM line; None if unusable."""
    de = next((float(t[5:]) for t in f[11:] if t.startswith("de:f:")), None)
    if de is None or de > DE_MAX:
        return None
    qa = qt = 0
    pos = int(f[3]); blocks = []; bs = pos
    cur = pos
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op in "=X":
            qa += n; qt += n; cur += n
        elif op == "I":
            qa += n; qt += n
        elif op == "S":
            qt += n
        elif op in "DN":
            if cur > bs:
                blocks.append((bs, cur))
            cur += n; bs = cur
        # H: not in query
    if cur > bs:
        blocks.append((bs, cur))
    if qa < 100 or qt == 0:
        return None
    return dict(de=de, qaln=qa / qt, pos=pos, blocks=blocks)


def scan(by_chrom, region):
    mm = collections.defaultdict(dict)
    for flag in ("-f", "0x100"), ("-F", "0x900"):
        p = subprocess.Popen([SAM, "view", flag[0], flag[1], NEW, region],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 11:
                continue
            pf = placement(f)
            if pf is None:
                continue
            g = best_gene(by_chrom, f[2], pf["pos"] - 1, pf["pos"] - 1 + ref_span(f[5]))
            if g is None:
                continue
            d = mm[f[0]]
            if g not in d or pf["de"] < d[g]["de"]:
                d[g] = pf
        p.wait()
    return mm


def covered_frac(blocks_list, glen):
    iv = sorted(b for bl in blocks_list for b in bl)
    if not iv:
        return 0.0
    cov, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            cov += ce - cs; cs, ce = s, e
    cov += ce - cs
    return min(cov / max(glen, 1), 1.0)


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    by_chrom, coord = collections.defaultdict(list), {}
    with open(GENES_BED) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in set(chroms):
                by_chrom[p[0]].append((int(p[1]), int(p[2]), p[3])); coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    for c in by_chrom:
        by_chrom[c].sort()
    Hd, _ = dna_homology()
    bam = pysam.AlignmentFile(NEW, "rb"); genome = pysam.FastaFile(GENOME)
    mcache = {}

    def model(g):
        if g not in mcache:
            c, s, e = coord[g]; mcache[g] = build_model(bam, genome, c, s, e)
        return mcache[g]

    def dna(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1
        return 0 if (r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30) else None

    rows = []
    for region in chroms:
        mm = scan(by_chrom, region)
        ev = collections.defaultdict(list)
        for genes in mm.values():
            if len(genes) < 2:
                continue
            it = sorted(genes.items())
            for i in range(len(it)):
                for j in range(i + 1, len(it)):
                    ga, pa = it[i]; gb, pb = it[j]
                    if ga != gb and abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX:
                        ev[(ga, gb)].append((pa, pb))
        for (ga, gb), recs in ev.items():
            if len(recs) < MIN_READS or ga not in coord or gb not in coord:
                continue
            lab = dna(ga, gb)
            if lab is None:
                continue
            rc = recip_cov(model(ga), model(gb))
            lenA = coord[ga][2] - coord[ga][1]; lenB = coord[gb][2] - coord[gb][1]
            qmin = np.median([min(pa["qaln"], pb["qaln"]) for pa, pb in recs])
            qtie = np.median([abs(pa["qaln"] - pb["qaln"]) for pa, pb in recs])
            cfA = covered_frac([pa["blocks"] for pa, pb in recs], lenA)
            cfB = covered_frac([pb["blocks"] for pa, pb in recs], lenB)
            rows.append(dict(ga=ga, gb=gb, lab=lab, recip=rc,
                             qaln_min=float(qmin), qaln_tie=float(qtie),
                             cov_frac_min=min(cfA, cfB)))
    bam.close(); genome.close()

    surv = [r for r in rows if r["recip"] >= COV_MIN]
    sl = np.array([r["lab"] for r in surv])
    nt, nf = int((sl == 0).sum()), int((sl == 1).sum())
    print(f"=== FLNC full-length levers on {len(surv)} ~B-survivors ({nt} true / {nf} false) ===")
    feats = ["qaln_min", "qaln_tie", "cov_frac_min"]
    print(f"  {'feature':14} {'bimod':>6} {'AUC_vs_DNA':>11}  true-median / false-median")
    res = []
    for fk in feats:
        a = auc([r[fk] for r in surv], sl); sep, _ = bimodality([r[fk] for r in surv])
        tm = np.median([r[fk] for r in surv if r["lab"] == 0]); fm = np.median([r[fk] for r in surv if r["lab"] == 1])
        res.append((fk, a)); print(f"  {fk:14} {sep:>6.2f} {a:>11.3f}  {tm:.3f} / {fm:.3f}"
                                   f"{'  <== separates' if a >= 0.75 else ''}")
    # per false survivor: their qaln_min (do bridges map fragment at 2nd locus?)
    print(f"\n  the false survivors' qaln_min (1.0 = molecule maps full-length at BOTH copies):")
    for r in sorted([r for r in surv if r["lab"] == 1], key=lambda x: x["qaln_min"]):
        print(f"    {r['ga']}~{r['gb']:22} qaln_min={r['qaln_min']:.3f} cov_frac_min={r['cov_frac_min']:.3f} recip={r['recip']:.2f}")
    # operating point of cov_frac_min + the residual that survives EVEN coverage
    tcf = sorted(r["cov_frac_min"] for r in surv if r["lab"] == 0)
    print(f"\n  TRUE paralog cov_frac_min: min={tcf[0]:.3f} p10={tcf[len(tcf)//10]:.3f} median={tcf[len(tcf)//2]:.3f}")
    for thr in [0.20, 0.30, 0.40]:
        fc = sum(1 for r in surv if r["lab"] == 1 and r["cov_frac_min"] < thr)
        tcost = sum(1 for r in surv if r["lab"] == 0 and r["cov_frac_min"] < thr)
        print(f"    cov_frac_min<{thr}: drops {fc}/{nf} false survivors, costs {tcost}/{nt} true paralogs")
    resid = [r for r in surv if r["lab"] == 1 and r["cov_frac_min"] >= 0.30]
    print(f"  IRREDUCIBLE after ~B + coverage>=0.30: {len(resid)} false survivors remain "
          f"(full-length homologous, zero cDNA = segmental-dup-like):")
    for r in resid:
        print(f"    {r['ga']}~{r['gb']}  cov_frac={r['cov_frac_min']:.2f} qaln={r['qaln_min']:.2f}")
    rng = np.random.default_rng(1729)
    real_max = max(a for _, a in res)
    null = [max(auc([r[fk] for r in surv], rng.permutation(sl)) for fk in feats) for _ in range(2000)]
    p = float((np.array(null) >= real_max).mean())
    best = max(res, key=lambda t: t[1])
    print(f"\n  best FLNC lever: {best[0]} AUC={best[1]:.3f}; permutation p={p:.4f} -> "
          f"{'real lever' if p < 0.05 else 'not past chance'}")
    json.dump(dict(survivors=len(surv), n_true=nt, n_false=nf,
                   feature_auc=[dict(feature=f, auc=a) for f, a in res], perm_p=p,
                   false_detail=[dict(ga=r["ga"], gb=r["gb"], qaln_min=round(r["qaln_min"], 3),
                                      cov_frac_min=round(r["cov_frac_min"], 3)) for r in surv if r["lab"] == 1]),
              open(os.path.join(HERE, "family_def_flnc_lever.json"), "w"), indent=2)
    print("\n[+] wrote family_def_flnc_lever.json")


if __name__ == "__main__":
    main()
