#!/usr/bin/env python3
"""Prereg Addendum AK: adjudicated two-annotation ground truth for multi-copy gene families, and its scorer.

build: adjudicated_truth.py build --out DIR --contigs c1,c2 --genome FA --sedef sd_v1.bed --hgnc hgnc.txt
         --ann NAME=NODES_TSV:TAG_PREFIX --ann NAME=NODES_TSV:TAG_PREFIX [--miniprot BIN] [--threads 4]
  NODES_TSV from `annotation_nodes.py` (with .names.tsv / .cds.tsv); TAG_PREFIX = the E1 construction's `<tag>`
  (reads `<tag>.clusters.tsv`, `<tag>.loci.tsv`). Exactly two --ann.
  Joint loci: records of both annotations whose exon unions share >= 1 bp. Opinion per annotation: SAME / DIFF / NONE.
  AGREED TRUE = SAME in both; DISPUTED = SAME in exactly one -> TRUE if evidence (P protein: miniprot Identity >= 0.70 and
  query coverage >= 0.30, either direction, protein = longest CDS of the locus; S: SEDEF pair with exon bases of each locus
  on opposite sides, linear projection within |lenA - lenB| + 1 kb), FALSE if both loci coding and no evidence, else
  UNSCORED. AK-0 gates are printed. Writes DIR/loci.tsv, DIR/pairs.tsv, DIR/clusters.tsv.
score: adjudicated_truth.py score --truth DIR --contigs c1,c2 [--expr loci_expr.tsv] name=copies.tsv ...
  Best-overlap assignment of loci to method families; TP/FP on TRUE/FALSE pairs, UNSCORED ignored; bipartite on the
  connected components of TRUE pairs (with --expr: loci with u >= 3 only, components recomputed).
"""
import argparse
import bisect
import collections
import csv
import itertools
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402

MIN_ID, MIN_COV, SD_SLACK = 0.70, 0.30, 1000
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, (_x, _y, _z) in enumerate(itertools.product(_b, _b, _b)):
    CODON[_x + _y + _z] = _aa[_i]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def translate(genome, chrom, strand, segs):
    segs = sorted(segs)
    seq = "".join(genome.fetch(chrom, s, e).upper() for s, e, _ in segs)
    phase = segs[0][2]
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
        phase = segs[-1][2]
    seq = seq[phase:]
    prot = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return prot.rstrip("*").replace("*", "X")


def load_records(name, spec, contigs):
    nodes_tsv, tag = spec.split(":", 1)
    names = {r["idx"]: (r["name"], r["biotype"]) for r in csv.DictReader(open(nodes_tsv + ".names.tsv"), delimiter="\t")}
    cds = {r["idx"]: (r["strand"], [(int(a), int(b.split(":")[0]), int(b.split(":")[1]))
                                    for a, b in (x.split("-") for x in r["cds"].split(","))])
           for r in csv.DictReader(open(nodes_tsv + ".cds.tsv"), delimiter="\t")}
    fam = {}
    for r in csv.DictReader(open(tag + ".clusters.tsv"), delimiter="\t"):
        fam[f"{r['chrom']}:{r['start']}-{r['end']}"] = r["cluster_id"]
    rep = {}
    if os.path.exists(tag + ".loci.tsv"):
        for r in csv.DictReader(open(tag + ".loci.tsv"), delimiter="\t"):
            rep[r["annotation"]] = r["representative"]
    recs = []
    for r in csv.DictReader(open(nodes_tsv), delimiter="\t"):
        if r["chrom"] not in contigs:
            continue
        key = f"{r['chrom']}:{int(r['start']) + 1}-{r['end']}"
        f = fam.get(key) or fam.get(rep.get(key, ""))
        ex = [tuple(map(int, b.split("-"))) for b in r["exons"].split(",")]
        recs.append({"ann": name, "idx": r["idx"], "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                     "exons": ex, "name": names[r["idx"]][0], "biotype": names[r["idx"]][1], "family": f,
                     "cds": cds.get(r["idx"])})
    return recs


def load_sedef(path, contigs):
    rows = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] in contigs and f[9] in contigs:
            rows[f[0]].append((int(f[1]), int(f[2]), f[9], int(f[10]), int(f[11]), f[13]))
    idx = {}
    for c, v in rows.items():
        v.sort()
        idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    return idx


def exon_hull_in(exons, s, e):
    parts = [(max(a, s), min(b, e)) for a, b in exons if a < e and s < b]
    return (min(p[0] for p in parts), max(p[1] for p in parts)) if parts else None


def sd_evidence(sd, u, v):
    if u["chrom"] not in sd:
        return False
    rows, starts, maxlen = sd[u["chrom"]]
    lo = bisect.bisect_left(starts, u["start"] - maxlen)
    hi = bisect.bisect_left(starts, u["end"])
    for a0, a1, oc, b0, b1, strand in rows[lo:hi]:
        if a1 <= u["start"] or oc != v["chrom"] or b1 <= v["start"] or b0 >= v["end"]:
            continue
        hu = exon_hull_in(u["exons"], a0, a1)
        hv = exon_hull_in(v["exons"], b0, b1)
        if not hu or not hv:
            continue
        ratio = (b1 - b0) / max(1, a1 - a0)
        if strand == "-":
            y = sorted((b1 - (hu[0] - a0) * ratio, b1 - (hu[1] - a0) * ratio))
        else:
            y = (b0 + (hu[0] - a0) * ratio, b0 + (hu[1] - a0) * ratio)
        tol = abs((b1 - b0) - (a1 - a0)) + SD_SLACK
        if y[0] - tol < hv[1] and hv[0] < y[1] + tol:
            return True
    return False


def cmd_build(a):
    os.makedirs(a.out, exist_ok=True)
    contigs = set(a.contigs.split(","))
    genome = pysam.FastaFile(a.genome)
    anns = [x.split("=", 1) for x in a.ann]
    assert len(anns) == 2, "exactly two --ann"
    (A, _), (Bn, _) = anns
    recs = [r for name, spec in anns for r in load_records(name, spec, contigs)]
    # joint loci: exon-union overlap across both annotations
    uf = UF()
    by_chrom = collections.defaultdict(list)
    for i, r in enumerate(recs):
        uf.find(i)
        for s, e in r["exons"]:
            by_chrom[r["chrom"]].append((s, e, i))
    for c, blocks in by_chrom.items():
        blocks.sort()
        cur_end, cur_i = -1, None
        for s, e, i in blocks:
            if cur_i is not None and s < cur_end:
                uf.union(cur_i, i)
            if e > cur_end:
                cur_end, cur_i = e, i
    groups = collections.defaultdict(list)
    for i in range(len(recs)):
        groups[uf.find(i)].append(i)
    loci = []
    for members in groups.values():
        rs = [recs[i] for i in members]
        ex = gp.merge([b for r in rs for b in r["exons"]])
        coding = [r for r in rs if r["cds"]]
        best = max(coding, key=lambda r: sum(e - s for s, e, _ in r["cds"][1])) if coding else None
        loci.append({"chrom": rs[0]["chrom"], "start": min(r["start"] for r in rs), "end": max(r["end"] for r in rs),
                     "exons": ex, "recs": rs, "coding": bool(coding),
                     "fam": {n: {r["family"] for r in rs if r["ann"] == n and r["family"]} for n in (A, Bn)},
                     "has": {n: any(r["ann"] == n for r in rs) for n in (A, Bn)},
                     "names": sorted({r["name"] for r in rs}), "best_cds": best})
    loci.sort(key=lambda l: (l["chrom"], l["start"], l["end"]))
    for k, l in enumerate(loci):
        l["id"] = k
    # SAME pairs per annotation
    same = {A: set(), Bn: set()}
    for n in (A, Bn):
        fam_loci = collections.defaultdict(set)
        for l in loci:
            for f in l["fam"][n]:
                fam_loci[f].add(l["id"])
        for ids in fam_loci.values():
            same[n].update(itertools.combinations(sorted(ids), 2))
    cand = same[A] | same[Bn]

    def opinion(n, u, v):
        if (u, v) in same[n]:
            return "SAME"
        return "DIFF" if loci[u]["has"][n] and loci[v]["has"][n] else "NONE"
    # HGNC hard negatives: coding loci sharing a gene group, DIFF in both annotations
    groups_of = collections.defaultdict(set)
    for r in csv.DictReader(open(a.hgnc), delimiter="\t"):
        for g in r["gene_group_id"].split("|"):
            if g:
                groups_of[r["symbol"]].add(g)
    grp_loci = collections.defaultdict(set)
    for l in loci:
        if l["coding"]:
            for nm in l["names"]:
                for g in groups_of.get(nm, ()):
                    grp_loci[g].add(l["id"])
    hard_neg = set()
    for ids in grp_loci.values():
        for u, v in itertools.combinations(sorted(ids), 2):
            if opinion(A, u, v) == "DIFF" and opinion(Bn, u, v) == "DIFF":
                hard_neg.add((u, v))
    agreed = {p for p in cand if p in same[A] and p in same[Bn]}
    disputed = cand - agreed
    need = disputed | {p for p in agreed if loci[p[0]]["coding"] and loci[p[1]]["coding"]} | hard_neg
    # protein evidence: every coding locus in a needed pair vs every locus in a needed pair
    involved = sorted({x for p in need for x in p})
    prot_fa, tgt_fa, mp_out = f"{a.out}/proteins.faa", f"{a.out}/targets.fa", f"{a.out}/miniprot.gff"
    with open(prot_fa, "w") as fp, open(tgt_fa, "w") as ft:
        for k in involved:
            l = loci[k]
            ft.write(f">L{k}\n{genome.fetch(l['chrom'], l['start'], l['end']).upper()}\n")
            if l["best_cds"]:
                r = l["best_cds"]
                p = translate(genome, r["chrom"], r["cds"][0], r["cds"][1])
                if len(p) >= 10:
                    fp.write(f">L{k}\n{p}\n")
    if not os.path.exists(mp_out):
        with open(mp_out + ".tmp", "w") as fh:
            subprocess.run([a.miniprot, "--gff", "-t", str(a.threads), "-N", "1000", "--outn=1000", "--outs=0", "-p", "0",
                            tgt_fa, prot_fa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.replace(mp_out + ".tmp", mp_out)
    prot_hit = set()
    paf = None
    for line in open(mp_out):
        if line.startswith("##PAF"):
            paf = line.rstrip("\n").split("\t")[1:]
        elif paf and "\tmRNA\t" in line:
            ident = next((float(x.split("=")[1]) for x in line.rstrip("\n").split("\t")[8].split(";") if x.startswith("Identity=")), 0.0)
            qlen, qs, qe = int(paf[1]), int(paf[2]), int(paf[3])
            if ident >= MIN_ID and qlen and (qe - qs) / qlen >= MIN_COV:
                u, v = int(paf[0][1:]), int(paf[5][1:])
                if u != v:
                    prot_hit.add((min(u, v), max(u, v)))
            paf = None
    sd = load_sedef(a.sedef, contigs)

    def evidence(p):
        u, v = loci[p[0]], loci[p[1]]
        P = p in prot_hit
        S = sd_evidence(sd, u, v) or sd_evidence(sd, v, u)
        return P, S
    status = {}
    ev = {}
    for p in sorted(need):
        ev[p] = evidence(p)
    for p in agreed:
        status[p] = "TRUE"
    for p in disputed:
        P, S = ev[p]
        if P or S:
            status[p] = "TRUE"
        elif loci[p[0]]["coding"] and loci[p[1]]["coding"]:
            status[p] = "FALSE"
        else:
            status[p] = "UNSCORED"
    with open(f"{a.out}/loci.tsv", "w") as fh:
        fh.write(f"locus\tchrom\tstart\tend\tcoding\tnames\t{A}_families\t{Bn}_families\texons\n")
        for l in loci:
            fh.write(f"L{l['id']}\t{l['chrom']}\t{l['start']}\t{l['end']}\t{int(l['coding'])}\t{','.join(l['names'])[:500]}\t"
                     f"{','.join(sorted(l['fam'][A])) or '.'}\t{','.join(sorted(l['fam'][Bn])) or '.'}\t"
                     f"{','.join(f'{x}-{y}' for x, y in l['exons'])}\n")
    with open(f"{a.out}/pairs.tsv", "w") as fh:
        fh.write(f"u\tv\t{A}\t{Bn}\tprotein\tsd\tstatus\n")
        for p in sorted(status):
            P, S = ev.get(p, ("NA", "NA"))
            fh.write(f"L{p[0]}\tL{p[1]}\t{opinion(A, *p)}\t{opinion(Bn, *p)}\t{P}\t{S}\t{status[p]}\n")
    uf2 = UF()
    for p, s in status.items():
        if s == "TRUE":
            uf2.union(p[0], p[1])
    comps = collections.defaultdict(list)
    for x in uf2.p:
        comps[uf2.find(x)].append(x)
    with open(f"{a.out}/clusters.tsv", "w") as fh:
        fh.write("cluster_id\tlocus\tchrom\tstart\tend\n")
        for k, (root, ms) in enumerate(sorted(comps.items())):
            for x in sorted(ms):
                l = loci[x]
                fh.write(f"T{k}\tL{x}\t{l['chrom']}\t{l['start'] + 1}\t{l['end']}\n")
    # AK-0
    cod_agreed = [p for p in agreed if loci[p[0]]["coding"] and loci[p[1]]["coding"]]
    sens = sum(1 for p in cod_agreed if any(ev[p])) / max(1, len(cod_agreed))
    fpr = sum(1 for p in hard_neg if any(ev[p])) / max(1, len(hard_neg))
    cnt = collections.Counter(status[p] for p in disputed)
    unsc = cnt["UNSCORED"] / max(1, len(disputed))
    by_src = collections.Counter((opinion(A, *p), opinion(Bn, *p), status[p]) for p in disputed)
    print(f"records {len(recs)}; joint loci {len(loci)}; SAME {A} {len(same[A])}, {Bn} {len(same[Bn])}; "
          f"agreed TRUE {len(agreed)}; disputed {len(disputed)} -> {dict(cnt)}")
    for k, v in sorted(by_src.items()):
        print(f"   disputed {k[0]:4s}/{k[1]:4s} -> {k[2]:8s} {v}")
    print(f"TRUE pairs {sum(1 for s in status.values() if s == 'TRUE')}; truth clusters {len(comps)} "
          f"(largest {max(map(len, comps.values())) if comps else 0}); proteins with a hit pair {len(prot_hit)}")
    ok = sens >= 0.80 and fpr <= 0.20 and unsc <= 0.50
    print(f"AK-0: evidence sensitivity on agreed coding pairs {sens:.3f} (n={len(cod_agreed)}, >= 0.80); "
          f"HGNC hard-negative evidence rate {fpr:.3f} (n={len(hard_neg)}, <= 0.20); "
          f"disputed unscored {unsc:.3f} (<= 0.50) -> {'VALID' if ok else 'NOT VALID'}")


def cmd_score(a):
    contigs = set(a.contigs.split(","))
    loci = {r["locus"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(f"{a.truth}/loci.tsv"), delimiter="\t")
            if r["chrom"] in contigs}
    expressed = None
    if a.expr:
        expressed = {r["locus"] for r in csv.DictReader(open(a.expr), delimiter="\t") if int(r["u"]) >= 3}
    status = {}
    for r in csv.DictReader(open(f"{a.truth}/pairs.tsv"), delimiter="\t"):
        if r["u"] in loci and r["v"] in loci and (expressed is None or (r["u"] in expressed and r["v"] in expressed)):
            status[(r["u"], r["v"])] = r["status"]
    true_pairs = [p for p, s in status.items() if s == "TRUE"]
    uf = UF()
    for u, v in true_pairs:
        uf.union(u, v)
    tl = sorted(uf.p, key=lambda x: int(x[1:]))
    tlabel = [uf.find(x) for x in tl]
    print(f"truth: {len(tl)} loci in {len(set(tlabel))} clusters, {len(true_pairs)} TRUE pairs, "
          f"{sum(1 for s in status.values() if s == 'FALSE')} explicit FALSE, {sum(1 for s in status.values() if s == 'UNSCORED')} UNSCORED")
    order = sorted(loci, key=lambda x: int(x[1:]))
    if expressed is not None:
        order = [x for x in order if x in expressed]
    print(f"{'catalog':16s} {'pair_sens':>9s} {'pair_prec':>9s} {'bip_R':>6s} {'bip_P':>6s} {'bip_F':>6s}  (TP / FP / ignored)")
    for spec in a.catalogs:
        name, path = spec.split("=", 1)
        by = collections.defaultdict(list)
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["chrom"] in contigs:
                by[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
        for c in by:
            by[c].sort()
        pred = {}
        for x in order:
            c, s, e = loci[x]
            h = [(min(e, b) - max(s, a0), f) for a0, b, f in by[c] if a0 < e and s < b]
            if h:
                pred[x] = max(h)[1]
        fam = collections.defaultdict(list)
        for x, f in pred.items():
            fam[f].append(x)
        tp = fp = ign = 0
        for ms in fam.values():
            for u, v in itertools.combinations(sorted(ms, key=lambda x: int(x[1:])), 2):
                s = status.get((u, v), "FALSE")
                if s == "TRUE":
                    tp += 1
                elif s == "FALSE":
                    fp += 1
                else:
                    ign += 1
        sens = tp / max(1, len(true_pairs))
        prec = tp / max(1, tp + fp)
        plabel = [pred.get(x, f"none:{x}") for x in tl]
        br, bp = gp.bipartite(plabel, tlabel)
        f1 = 2 * br * bp / (br + bp) if br + bp else float("nan")
        print(f"{name:16s} {sens:9.3f} {prec:9.3f} {br:6.3f} {bp:6.3f} {f1:6.3f}  ({tp} / {fp} / {ign})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    for k in ("--out", "--contigs", "--genome", "--sedef", "--hgnc"):
        p.add_argument(k, required=True)
    p.add_argument("--ann", action="append", required=True)
    p.add_argument("--miniprot", default="/home/juanfra/miniforge3/envs/prot/bin/miniprot")
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("score")
    p.add_argument("--truth", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--expr")
    p.add_argument("catalogs", nargs="+")
    a = ap.parse_args()
    {"build": cmd_build, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
