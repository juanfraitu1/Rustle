#!/usr/bin/env python3
"""Family scores of one clusters file on one contig against the human family truths, with `family_score` called as the container-headroom
run called it (docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md): clusters and truth restricted to the contig, `--chrom ALL --pairwise --per-family`.

    score_families.py --fs family_score --clusters X.fam.clusters.tsv --contig chr16 --label NAME --out NAME.json [--work DIR]

Truths (families_gw/species/human): compara = Ensembl Compara Primates, u2 = the NPIP union truth (Soto + RefSeq NPIP names + homology-admitted
genes, 3 families), soto = Soto 2025. Per truth: pooled bipartite sens / prec / F (family_score's own matching), pairwise tp / truth / predicted pairs,
and every truth family's row (matched cluster, hit, sens, prec, F). Writes the JSON and prints one line per truth.
"""
import argparse
import json
import os
import re
import subprocess

HUM = "/mnt/linuxdisk/tmp/rustle_figures/families_gw/species/human"
TRUTHS = {"compara": "compara.Primates.families.tsv", "u2": "npip_u2.families.tsv", "soto": "soto.families.tsv"}
M_POOL = re.compile(r"truth (\d+) fams / (\d+) genes \| clusters (\d+) \| sens ([\d.]+) prec ([\d.]+) F ([\d.]+)"
                    r"(?: \| collapsed (\d+))?(?: \| no-locus (\d+))?")
M_PAIR = re.compile(r"pairwise \| truth pairs (\d+) \| predicted pairs (\d+) \| tp (\d+) \| sens ([\d.]+) prec ([\d.]+) F ([\d.]+)")


def keep(path, col, contig, out):
    n = 0
    with open(path) as fh, open(out, "w") as o:
        hdr = fh.readline()
        o.write(hdr)
        ci = hdr.rstrip("\n").split("\t").index(col)
        for line in fh:
            if line.rstrip("\n").split("\t")[ci] == contig:
                o.write(line)
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fs", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--contig", default="chr16")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=None)
    ap.add_argument("--truths", default=",".join(TRUTHS), help="comma-separated subset of compara,u2,soto (default all three; the u2 truth is NPIP-only, so chr17 uses compara,soto)")
    a = ap.parse_args()
    work = a.work or os.path.dirname(os.path.abspath(a.out))
    os.makedirs(work, exist_ok=True)
    c_in = f"{work}/{a.label}.clusters.tsv"
    keep(a.clusters, "chrom", a.contig, c_in)
    out = dict(label=a.label, clusters=a.clusters, contig=a.contig, truths={})
    for key, fn in ((k, TRUTHS[k]) for k in a.truths.split(",")):
        t_in, pf = f"{work}/{a.label}.{key}.truth.tsv", f"{work}/{a.label}.{key}.pf.tsv"
        keep(f"{HUM}/{fn}", "Contig", a.contig, t_in)
        text = subprocess.run([a.fs, "--clusters", c_in, "--gff", f"{HUM}/genes_only.gff", "--soto", t_in, "--chrom", "ALL",
                               "--label", f"{a.label}_{key}", "--pairwise", "--per-family", pf],
                              capture_output=True, text=True, check=True).stdout
        m, q = M_POOL.search(text), M_PAIR.search(text)
        if not m:
            raise SystemExit(f"family_score printed no pooled line for {key}:\n{text}")
        r = dict(truth_fams=int(m.group(1)), truth_genes=int(m.group(2)), clusters=int(m.group(3)), sens=float(m.group(4)),
                 prec=float(m.group(5)), f=float(m.group(6)), collapsed=int(m.group(7) or 0), no_locus=int(m.group(8) or 0))
        if q:
            r.update(pair_truth=int(q.group(1)), pair_pred=int(q.group(2)), pair_tp=int(q.group(3)), pair_sens=float(q.group(4)),
                     pair_prec=float(q.group(5)), pair_f=float(q.group(6)))
        with open(pf) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            r["per_family"] = [dict(zip(hdr, ln.rstrip("\n").split("\t"))) for ln in fh]
        out["truths"][key] = r
        print(f"{a.label} {key}: sens {r['sens']:.3f} prec {r['prec']:.3f} F {r['f']:.3f} | pairs tp {r.get('pair_tp')}/{r.get('pair_truth')} truth, "
              f"{r.get('pair_pred')} predicted")
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
