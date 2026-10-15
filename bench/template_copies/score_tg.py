#!/usr/bin/env python3
"""Scores of the template-guided copy stage (spec section 2).

    score_tg.py liftoff --species S --sample ID --base BASE.gtf --merged MERGED.gtf --tg TG.gtf --templates T.tsv --out J
    score_tg.py chains  --ref REF.gtf[.gz] --work DIR --gtf NAME=PATH [--gtf ...] --out J
    score_tg.py tools   --species S --ref REF --exclude CHROM --strata PREFIX --work DIR --arm NAME=GTF ... --tg-arm NAME --out J

liftoff  recall of read-supported Liftoff extra copies (identity >= 0.95), base vs merged, and how many template-reachable
         copies (the copy unassembled, its source placement assembled) were reached; P1's located share, base vs
         merged; P2: the added loci's located share vs base loci whose representative has 2-5 reads, with a 95% interval
         over loci. Located = covered >= 50% by one placed Liftoff locus (figure 8's `located`).
chains   gffcompare 0.12.10 of each GTF against the reference: intron-chain sensitivity / precision (figure 1's rule).
tools    register 1268's scorer (bench/tool_difference/score_op.py) on the arms; B1 (chain sensitivity of the TG arm >
         StringTie's and > FLAIR's) and B2 (its chain precision > FLAIR's).
"""
import argparse
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "figures"))
import _liftoff as LF  # noqa: E402
import templates as T  # noqa: E402

LIFTOFF = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
SC = 0.95
BOOT_N, BOOT_SEED = 2000, 20260925
LEVEL = re.compile(r"^\s*(Base|Exon|Intron|Intron chain|Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)")
ROW = re.compile(r"^(\S+)\s+" + r"\s+".join([r"([\d.]+)\s*/\s*([\d.]+)"] * 5) + r"\s*$")


def placed_index(loci):
    return LF.Index([r for r in loci if r["cls"] in ("in_place", "moved", "extra_copy")])


def located(rows, gidx):
    """Per row (a GTF locus): covered >= 50% by one placed Liftoff locus."""
    return [gidx.best_cov(r["contig"], r["iv"])[0] >= LF.MATCH_COV for r in rows]


def gtf_rows(path, cache):
    return LF.denovo_loci(Path(path), Path(cache))


def interval(hits, n=BOOT_N, seed=BOOT_SEED):
    """95% percentile interval of the mean of 0/1 hits, resampling loci."""
    if not hits:
        return [float("nan"), float("nan")]
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(hits) for _ in hits) / len(hits) for _ in range(n))
    return [means[int(0.025 * n)], means[int(0.975 * n) - 1]]


def liftoff_scores(loci, sup, base_rows, merged_rows, tg_rows, base_reads):
    """loci, sup: figure 8's Liftoff table and the sample's read support; *_rows: GTF loci (`locus` = gene_id);
    base_reads: {gene_id: representative reads}."""
    copies = [r for r in loci if r["cls"] == "extra_copy" and not r["short"] and r["seqid_f"] >= SC - 1e-9
              and sup.get(LF.locus_key(r), 0) >= LF.SUPPORT_READS]
    placement = {r["source_id"]: r for r in loci if r["cls"] in ("in_place", "moved") and not r["short"]}
    bidx, midx = LF.Index(base_rows), LF.Index(merged_rows)
    cov = lambda idx, r: idx.best_cov(r["contig"], r["iv"])[0] >= LF.MATCH_COV  # noqa: E731
    in_base = [cov(bidx, c) for c in copies]
    in_merged = [cov(midx, c) for c in copies]
    reachable = [i for i, c in enumerate(copies) if not in_base[i] and c["source_id"] in placement
                 and cov(bidx, placement[c["source_id"]])]
    gidx = placed_index(loci)
    lb, lm, lt = located(base_rows, gidx), located(merged_rows, gidx), located(tg_rows, gidx)
    low = [h for r, h in zip(base_rows, lb) if 2 <= base_reads.get(r["locus"], 0) <= 5]
    share = lambda xs: sum(xs) / len(xs) if xs else float("nan")  # noqa: E731
    return dict(
        recall=dict(copies=len(copies), base=sum(in_base), merged=sum(in_merged),
                    reachable=len(reachable), reached=sum(in_merged[i] for i in reachable)),
        p1_located=dict(base=share(lb), merged=share(lm), passes=bool(share(lm) >= share(lb) - 0.002)),
        p2=dict(tg_loci=len(lt), tg_located=share(lt), tg_ci95=interval([int(x) for x in lt]),
                base_low_loci=len(low), base_low_located=share(low),
                passes=bool(lt) and share(lt) >= share(low)))


def parse_stats(path):
    out = {}
    for line in open(path):
        m = LEVEL.match(line)
        if m:
            k = m.group(1).lower().replace(" ", "_")
            out[f"{k}_sn"], out[f"{k}_pr"] = float(m.group(2)), float(m.group(3))
    return out


def plain(ref, work):
    """gffcompare 0.12.10 reads plain text only: a .gz reference is decompressed once into the work dir."""
    if not str(ref).endswith(".gz"):
        return ref
    dst = os.path.join(work, "ref.gtf")
    if not os.path.exists(dst):
        with gzip.open(ref, "rt") as src, open(dst + ".tmp", "w") as fh:
            shutil.copyfileobj(src, fh)
        os.replace(dst + ".tmp", dst)
    return dst


def chains(ref, work, gtfs):
    os.makedirs(work, exist_ok=True)
    ref = plain(ref, work)
    out = {}
    for name, path in gtfs:
        pre = os.path.join(work, f"gc_{name}")
        if not os.path.exists(pre + ".stats"):
            subprocess.run(["gffcompare", "-r", ref, "-o", pre, path], check=True, cwd=work,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        st = parse_stats(pre + ".stats")
        out[name] = dict(chain_sn=st["intron_chain_sn"], chain_pr=st["intron_chain_pr"])
    return out


def parse_score_op(text):
    """score_op.py's arm table -> {arm: (chain_sn, chain_pr)}."""
    out = {}
    for line in text.splitlines():
        m = ROW.match(line)
        if m and m.group(1) != "arm":
            v = [float(x) for x in m.groups()[1:]]
            out[m.group(1)] = (v[4], v[5])
    return out


def tools(args):
    cmd = [sys.executable, os.path.join(REPO, "bench", "tool_difference", "score_op.py"), "--species", args.species,
           "--ref", args.ref, "--exclude", args.exclude, "--strata", args.strata, "--work", args.work]
    for a in args.arm:
        cmd += ["--arm", a]
    os.makedirs(args.work, exist_ok=True)
    text = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    arms = parse_score_op(text)
    return dict(table=text, chain={k: list(v) for k, v in arms.items()}, **lead(arms, args.tg_arm))


def lead(arms, tg_arm):
    """B1 / B2 of register 1268 for the TG arm; None when an arm is missing (a cache-filling call with one arm)."""
    if not {tg_arm, "StringTie", "FLAIR"} <= set(arms):
        return dict(b1=None, b2=None)
    tg, st, fl = arms[tg_arm], arms["StringTie"], arms["FLAIR"]
    return dict(b1=bool(tg[0] > st[0] and tg[0] > fl[0]), b2=bool(tg[1] > fl[1]))


def write(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("liftoff")
    for k in ("--species", "--sample", "--base", "--merged", "--tg", "--templates", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--liftoff", default=LIFTOFF)
    p = sub.add_parser("chains")
    p.add_argument("--ref", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--gtf", action="append", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("tools")
    for k in ("--species", "--ref", "--exclude", "--strata", "--work", "--tg-arm", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--arm", action="append", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "liftoff":
        loci = LF.read_loci(Path(a.liftoff) / a.species / "liftoff_loci.tsv")
        sup = LF.read_support(Path(a.liftoff) / "support" / f"{a.sample}.tsv")
        cache = os.path.join(os.path.dirname(os.path.abspath(a.out)), "denovo_cache")
        base_reads = {t.gene: t.reads for t in T.read_templates(a.templates).values()}
        res = liftoff_scores(loci, sup, gtf_rows(a.base, cache), gtf_rows(a.merged, cache), gtf_rows(a.tg, cache),
                             base_reads)
        write(a.out, dict(sample=a.sample, **res))
    elif a.cmd == "chains":
        write(a.out, chains(a.ref, a.work, [g.split("=", 1) for g in a.gtf]))
    else:
        write(a.out, tools(a))
    print(open(a.out).read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
