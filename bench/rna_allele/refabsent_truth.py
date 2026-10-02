#!/usr/bin/env python3
"""Truth for item 2 (real reference-absent copies in KB3781): haplotype-only (B-only) loci of every family of the 2026-08-14 interval table,
and which of them are expressed in the fibroblast reads. B = the haplotype `_pri` did NOT take that chromosome from (chrmap.tsv).

  copies    copies.fa (clean interval of every copy, from `_pri`), copies_genes.tsv (for truth_lift), panel.json (control_test layout)
  lift      copies_lift.tsv: every copy interval -> its B interval (truth_lift.py over the frozen asm5 alignments)
  bonly     copies.{mat,pat}.paf (minimap2 -c -x asm20 <hap>.fa copies.fa) -> hits identity >= 0.90, coverage >= 0.80 on a B-side
            chromosome that overlap the lifted B interval of NO copy of the same family; overlapping hits merge -> bonly.tsv
  reads     net reads per family from the baseline BAM (a record on any copy), all of them up to --cap; scored.fa / labels.tsv (control_test
            layout) for the families listed in --fams (default: every family)
  express   reads aligned to mat and pat (splice:hq, -N 50): a read sits on a B-only locus when its best record (AS, untied at 0.98)
            overlaps it; expressed = >= 3 such reads -> bonly_expressed.tsv

    refabsent_truth.py copies --w /mnt/linuxdisk/tmp/rna_allele/refabsent
"""
import argparse
import collections
import csv
import json
import os
import random
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_lift  # noqa: E402

SRC = "/home/juanfra/winloci_scratch/o3_collapse/method/intervals/data/intervals.tsv"
PRI = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
TRUTH = "/mnt/linuxdisk/tmp/rna_allele"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def table():
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(SRC), delimiter="\t"):
        by[r["fam"]].append((r["chrom"], int(r["clean_start"]), int(r["clean_end"]), r["gene"], int(r["n_clean"])))
    return by


def chrmap():
    return {r["pri"]: r for r in csv.DictReader(open(f"{TRUTH}/chrmap.tsv"), delimiter="\t")}


def copies(a):
    by = table()
    fa = pysam.FastaFile(PRI)
    with open(f"{a.w}/copies.fa", "w") as o, open(f"{a.w}/copies_genes.tsv", "w") as g:
        g.write("gene_id\tchrom\tstrand\texons\n")
        for fam, cps in sorted(by.items()):
            for c, s, e, gene, n in cps:
                o.write(f">{gene}\n{fa.fetch(c, s, e)}\n")
                g.write(f"{gene}\t{c}\t+\t{s}-{e}\n")
    panel = [dict(fam=fam, mask=list(cps[0][:4]), keep=[list(x[:4]) for x in cps[1:]]) for fam, cps in sorted(by.items())]
    json.dump(panel, open(f"{a.w}/panel.json", "w"), indent=0)
    print(f"families {len(by)}, copies {sum(len(v) for v in by.values())}, bp {sum(e - s for v in by.values() for _, s, e, _, _ in v)}")


def lift(a):
    truth_lift.main(["--chrmap", f"{TRUTH}/chrmap.tsv", "--paf-dir", f"{TRUTH}/out", "--genes", f"{a.w}/copies_genes.tsv",
                     "--out", f"{a.w}/copies_lift.tsv"])
    rows = list(csv.DictReader(open(f"{a.w}/copies_lift.tsv"), delimiter="\t"))
    print(f"copies {len(rows)}; lifted >= 50%: {sum(1 for r in rows if float(r['lift_frac']) >= 0.5)}; "
          f"classes {dict(collections.Counter(r['class'] for r in rows))}")


def bonly(a):
    by, cm = table(), chrmap()
    fam_of = {gene: fam for fam, cps in by.items() for _, _, _, gene, _ in cps}
    bside = {(r["B_hap"], r["B_name"]): r["pri"] for r in cm.values() if r["B_name"]}
    liftd = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{a.w}/copies_lift.tsv"), delimiter="\t"):
        if float(r["lift_frac"]) >= 0.5:
            liftd[fam_of[r["gene_id"]]].append((r["B_chrom"], int(r["B_start"]), int(r["B_end"]), r["gene_id"]))
    hits = collections.defaultdict(list)
    for h in ("mat", "pat"):
        for ln in open(f"{a.w}/copies.{h}.paf"):
            f = ln.rstrip("\n").split("\t")
            ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
            if ident < 0.90 or cov < 0.80 or (h, f[5]) not in bside:
                continue
            fam = fam_of[f[0]]
            s, e = int(f[7]), int(f[8])
            if any(c == f[5] and s < e0 and s0 < e for c, s0, e0, _ in liftd[fam]):
                continue
            hits[fam].append((h, f[5], s, e, ident, cov, f[0]))
    n_loci = 0
    with open(f"{a.w}/bonly.tsv", "w") as o:
        o.write("family\tlocus\thap\tchrom\tstart\tend\tn_hits\tbest_identity\tfrom_copies\n")
        for fam, hs in sorted(hits.items()):
            hs.sort(key=lambda x: (x[0], x[1], x[2]))
            cur = None
            loci = []
            for h in hs:
                if cur and cur[0] == h[0] and cur[1] == h[1] and h[2] < cur[3]:
                    cur[3] = max(cur[3], h[3]); cur[4].append(h)
                else:
                    cur = [h[0], h[1], h[2], h[3], [h]]; loci.append(cur)
            for k, L in enumerate(loci):
                n_loci += 1
                o.write(f"{fam}\t{fam}_B{k}\t{L[0]}\t{L[1]}\t{L[2]}\t{L[3]}\t{len(L[4])}\t{max(x[4] for x in L[4]):.4f}\t"
                        f"{','.join(sorted({x[6] for x in L[4]}))}\n")
    fams = {r["family"] for r in csv.DictReader(open(f"{a.w}/bonly.tsv"), delimiter="\t")}
    print(f"B-only loci {n_loci} in {len(fams)} of {len(by)} families; by copy number of the family: "
          f"{dict(sorted(collections.Counter(len(by[f]) for f in fams).items()))}")
    # the locus sequences, for the back-alignment to `_pri` (bclass)
    hfa = {h: pysam.FastaFile(f"/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{h}.fa") for h in ("mat", "pat")}
    with open(f"{a.w}/bonly.fa", "w") as o:
        for r in csv.DictReader(open(f"{a.w}/bonly.tsv"), delimiter="\t"):
            o.write(f">{r['locus']}\n{hfa[r['hap']].fetch(r['chrom'], int(r['start']), int(r['end']))}\n")


def bclass(a):
    """bonly.pri.paf (minimap2 -c -x asm20 -p 0.1 -N 50 GGO.fasta bonly.fa): best `_pri` identity per locus (hits covering >= 50% of
    it) -> absent_beyond_delta (< 1 - 0.00958, detectable in principle) / absent_within_delta (an allele of a `_pri` locus the table does
    not list, or a near-identical recent duplicate: indistinguishable from an allele by construction)."""
    best = {}
    for ln in open(f"{a.w}/bonly.pri.paf"):
        f = ln.rstrip("\n").split("\t")
        ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if cov >= 0.5 and ident > best.get(f[0], (0.0, ""))[0]:
            best[f[0]] = (ident, f"{f[5]}:{f[7]}-{f[8]}")
    rows = list(csv.DictReader(open(f"{a.w}/bonly.tsv"), delimiter="\t"))
    with open(f"{a.w}/bonly.tsv", "w") as o:
        o.write("\t".join(list(rows[0].keys())[:9] + ["pri_best_identity", "pri_best_hit", "class"]) + "\n")
        for r in rows:
            b = best.get(r["locus"], (0.0, ""))
            cls = "absent_beyond_delta" if b[0] < 1 - 0.00958 else "absent_within_delta"
            o.write("\t".join([r[k] for k in list(r.keys())[:9]] + [f"{b[0]:.4f}", b[1], cls]) + "\n")
    rows = list(csv.DictReader(open(f"{a.w}/bonly.tsv"), delimiter="\t"))
    for cls in ("absent_beyond_delta", "absent_within_delta"):
        rs = [r for r in rows if r["class"] == cls]
        print(f"{cls}: {len(rs)} loci in {len({r['family'] for r in rs})} families; `_pri` identity median "
              f"{sorted(float(r['pri_best_identity']) for r in rs)[len(rs) // 2] if rs else 0:.4f}")


def reads(a):
    by = table()
    fams = sorted(by) if not a.fams else [l.strip() for l in open(a.fams) if l.strip()]
    bam = pysam.AlignmentFile(BAM)
    rng = random.Random(1)
    role = {}
    for fam in fams:
        names = collections.defaultdict(set)
        for c, s, e, gene, n in by[fam]:
            for rd in bam.fetch(c, s, e):
                if not (rd.is_unmapped or rd.is_supplementary):
                    names[rd.query_name].add(gene if not rd.is_secondary else "sec:" + gene)
        ns = sorted(names)
        rng.shuffle(ns)
        for n in ns[:a.cap]:
            prim = sorted(x for x in names[n] if not x.startswith("sec:"))
            role.setdefault(n, (fam, prim[0] if prim else "sec"))
    with open(f"{a.w}/names.txt", "w") as o:
        o.write("\n".join(sorted(role)) + "\n")
    subprocess.run(f"samtools view -F 2308 -N {a.w}/names.txt -@ 4 {BAM} -o {a.w}/prim.sam", shell=True, check=True)
    seq = {}
    for ln in open(f"{a.w}/prim.sam"):
        f = ln.split("\t", 11)
        s = f[9]
        if int(f[1]) & 16:
            s = s.translate(COMP)[::-1]
        seq[f[0]] = s
    with open(f"{a.w}/scored.fa", "w") as fa, open(f"{a.w}/labels.tsv", "w") as lab:
        lab.write("read\tfamily\trole\tcopy\n")
        for n, (fam, g) in sorted(role.items()):
            if n in seq:
                fa.write(f">{n}\n{seq[n]}\n"); lab.write(f"{n}\t{fam}\tS\t{g}\n")
    per = collections.Counter(f for f, _ in role.values())
    print(f"families {len(fams)}; reads {len(role)} (sequences {len(seq)}); per family min/median/max "
          f"{min(per.values())}/{sorted(per.values())[len(per) // 2]}/{max(per.values())}")


def express(a):
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.w}/labels.tsv"), delimiter="\t")}
    B = list(csv.DictReader(open(f"{a.w}/bonly.tsv"), delimiter="\t"))
    alias = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            alias[(h, num)] = acc
    recs = collections.defaultdict(list)
    for h in ("mat", "pat"):
        for rd in pysam.AlignmentFile(f"{a.w}/reads.{h}.bam").fetch(until_eof=True):
            if rd.is_unmapped or rd.is_supplementary:
                continue
            m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]*", rd.reference_name)
            acc = alias[(m.group(2), m.group(1))] if m else rd.reference_name
            recs[rd.query_name].append((rd.get_tag("AS") if rd.has_tag("AS") else 0, h, acc, rd.reference_start, rd.reference_end))
    on = collections.defaultdict(set)
    for n, rs in recs.items():
        rs.sort(key=lambda r: -r[0])
        if len(rs) > 1 and rs[1][0] >= 0.98 * rs[0][0]:
            continue
        _, h, acc, s, e = rs[0]
        fam = lab[n]["family"]
        for b in B:
            if b["family"] == fam and b["hap"] == h and b["chrom"] == acc and s < int(b["end"]) and int(b["start"]) < e:
                on[b["locus"]].add(n)
    with open(f"{a.w}/bonly_expressed.tsv", "w") as o:
        o.write("locus\tfamily\tn_reads\texpressed\n")
        for b in B:
            n = len(on.get(b["locus"], ()))
            o.write(f"{b['locus']}\t{b['family']}\t{n}\t{int(n >= 3)}\n")
    ex = [b for b in B if len(on.get(b["locus"], ())) >= 3]
    print(f"B-only loci {len(B)}; expressed (>= 3 reads best-placed there, untied) {len(ex)} in {len({b['family'] for b in ex})} families; "
          f"reads on expressed loci: median {sorted(len(on[b['locus']]) for b in ex)[len(ex) // 2] if ex else 0}, "
          f"max {max((len(on[b['locus']]) for b in ex), default=0)}; loci with 1-2 reads {sum(1 for b in B if 0 < len(on.get(b['locus'], ())) < 3)}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["copies", "lift", "bonly", "bclass", "reads", "express"])
    ap.add_argument("--w", required=True)
    ap.add_argument("--fams", default=None)
    ap.add_argument("--cap", type=int, default=2000)
    a = ap.parse_args(argv)
    os.makedirs(a.w, exist_ok=True)
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
