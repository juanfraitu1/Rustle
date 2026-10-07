#!/usr/bin/env python3
"""Family scores of an ideal-expression run on the universe of simulated genes (docs/PREREG_ideal_expression_2026-10-06.md, endpoint E4 beside).

    fam_score.py --fs family_score --clusters RUN.fam.clusters.tsv --contig chr16 --windows WINDOWS.tsv|ALL --label L --out OUT.json

`family_score` is called as bench/default_rescore/score_families.py calls it (clusters and truth restricted to the contig, `--chrom ALL --pairwise --per-family`), with the truth
restricted to the UNIVERSE: only truth genes whose RefSeq record (families_gw genes_only.gff, matched by name; every record of a duplicated name counts) overlaps a simulated window
are kept (`--windows ALL` keeps the whole contig), and a truth family with fewer than 2 kept genes is dropped. Truths: compara, u2 (chr16 only), soto.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "default_rescore"))
import score_families as sf  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fs", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=None)
    a = ap.parse_args()
    work = a.work or os.path.dirname(os.path.abspath(a.out))
    os.makedirs(work, exist_ok=True)
    wins = []
    if a.windows != "ALL":
        wins = [(int(r["start0"]), int(r["end"])) for r in csv.DictReader(open(a.windows), delimiter="\t") if r["chrom"] == a.contig]
    recs = {}
    for ln in open(f"{sf.HUM}/genes_only.gff"):
        f = ln.split("\t")
        if len(f) > 8 and f[0] == a.contig:
            m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if m:
                recs.setdefault(m.group(1), []).append((int(f[3]) - 1, int(f[4])))
    inside = lambda nm: a.windows == "ALL" or any(s < w1 and e > w0 for s, e in recs.get(nm, []) for w0, w1 in wins)
    c_in = f"{work}/{a.label}.clusters.tsv"
    sf.keep(a.clusters, "chrom", a.contig, c_in)
    out = dict(label=a.label, clusters=a.clusters, contig=a.contig, windows=a.windows, truths={})
    for key, fn in sf.TRUTHS.items():
        rows = list(csv.reader(open(f"{sf.HUM}/{fn}"), delimiter="\t"))
        hdr = rows[0]; ci, gi, fi = hdr.index("Contig"), hdr.index("Gene Name"), hdr.index("Family ID")
        kept = [r for r in rows[1:] if r[ci] == a.contig and inside(r[gi])]
        cnt = {}
        for r in kept:
            cnt[r[fi]] = cnt.get(r[fi], 0) + 1
        kept = [r for r in kept if cnt[r[fi]] >= 2]
        if not kept:
            out["truths"][key] = dict(skipped="no truth family with >= 2 genes in the universe")
            print(f"{a.label} {key}: no truth family with >= 2 genes in the universe")
            continue
        t_in, pf = f"{work}/{a.label}.{key}.truth.tsv", f"{work}/{a.label}.{key}.pf.tsv"
        with open(t_in, "w") as fh:
            fh.write("\t".join(hdr) + "\n"); [fh.write("\t".join(r) + "\n") for r in kept]
        text = subprocess.run([a.fs, "--clusters", c_in, "--gff", f"{sf.HUM}/genes_only.gff", "--soto", t_in, "--chrom", "ALL", "--label", f"{a.label}_{key}",
                               "--pairwise", "--per-family", pf], capture_output=True, text=True, check=True).stdout
        m, q = sf.M_POOL.search(text), sf.M_PAIR.search(text)
        r = dict(kept_genes=len(kept), kept_families=len(cnt))
        if m:
            r.update(truth_fams=int(m.group(1)), truth_genes=int(m.group(2)), clusters=int(m.group(3)), sens=float(m.group(4)), prec=float(m.group(5)), f=float(m.group(6)))
        if q:
            r.update(pair_truth=int(q.group(1)), pair_pred=int(q.group(2)), pair_tp=int(q.group(3)), pair_f=float(q.group(6)))
        with open(pf) as fh:
            h2 = fh.readline().rstrip("\n").split("\t")
            r["per_family"] = [dict(zip(h2, ln.rstrip("\n").split("\t"))) for ln in fh]
        out["truths"][key] = r
        print(f"{a.label} {key}: universe {len(kept)} genes in {len(cnt)} families | sens {r.get('sens')} prec {r.get('prec')} F {r.get('f')} | pairs tp {r.get('pair_tp')}/{r.get('pair_truth')} truth, {r.get('pair_pred')} predicted")
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
