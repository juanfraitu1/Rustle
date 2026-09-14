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


def blast_evidence(a, loci, involved, genome):
    """X(u -> v): dc-megablast of u's exon-union sequence (both annotations; soft-masked repeats do not seed) against v's
    genomic span; non-overlapping HSPs (on the query, best nident first) sum to >= 300 aligned bp at identity >= 0.70 and
    cover >= 0.30 of u's exonic length. Returns the set of (u, v) with X."""
    q_fa, t_fa, db, out = (f"{a.out}/exons.fa", f"{a.out}/spans.fa", f"{a.out}/spans_db", f"{a.out}/blast.tsv")
    qlen = {}
    with open(q_fa, "w") as fq, open(t_fa, "w") as ft:
        for k in involved:
            l = loci[k]
            seq = "".join(genome.fetch(l["chrom"], s, e) for s, e in l["exons"])
            qlen[k] = len(seq)
            fq.write(f">L{k}\n{seq}\n")
            ft.write(f">L{k}\n{genome.fetch(l['chrom'], l['start'], l['end'])}\n")
    if not os.path.exists(out):
        subprocess.run([a.blast_bin + "/makeblastdb", "-dbtype", "nucl", "-in", t_fa, "-out", db],
                       stdout=subprocess.DEVNULL, check=True)
        with open(out + ".tmp", "w") as fh:
            subprocess.run([a.blast_bin + "/blastn", "-task", "dc-megablast", "-query", q_fa, "-db", db, "-lcase_masking",
                            "-evalue", "1e-5", "-max_target_seqs", "100000", "-num_threads", str(a.threads),
                            "-outfmt", "6 qseqid sseqid qstart qend nident length"], stdout=fh, check=True)
        os.replace(out + ".tmp", out)
    hsps = collections.defaultdict(list)
    for line in open(out):
        qs_, ss_, q0, q1, nid, ln = line.split("\t")
        if qs_ != ss_:
            q0, q1 = sorted((int(q0), int(q1)))
            hsps[(int(qs_[1:]), int(ss_[1:]))].append((int(nid), int(ln), q0 - 1, q1))
    hit = set()
    for (u, v), hs in hsps.items():
        taken, L, N = [], 0, 0
        for nid, ln, q0, q1 in sorted(hs, reverse=True):
            if any(q0 < y and x < q1 for x, y in taken):
                continue
            taken.append((q0, q1))
            L += ln
            N += nid
        cov = sum(y - x for x, y in gp.merge(taken)) / max(1, qlen[u])
        if L >= 300 and N / L >= MIN_ID and cov >= MIN_COV:
            hit.add((u, v))
    return hit


def overlap_groups(idx_list, recs):
    """Union-find over records whose exon unions share >= 1 bp (the `--merge-overlapping-loci` rule)."""
    uf = UF()
    by_chrom = collections.defaultdict(list)
    for i in idx_list:
        uf.find(i)
        for s, e in recs[i]["exons"]:
            by_chrom[recs[i]["chrom"]].append((s, e, i))
    for blocks in by_chrom.values():
        blocks.sort()
        cur_end, cur_i = -1, None
        for s, e, i in blocks:
            if cur_i is not None and s < cur_end:
                uf.union(cur_i, i)
            if e > cur_end:
                cur_end, cur_i = e, i
    groups = collections.defaultdict(list)
    for i in idx_list:
        groups[uf.find(i)].append(i)
    return list(groups.values())


def joint_groups(recs, A, Bn, mode):
    """union (AL): one union-find over both annotations — chains distinct genes through overlapping models of the other
    annotation (§6kn). matched (AM): each annotation's own loci, then a RefSeq locus and a GENCODE locus are one joint
    locus iff each is the other's best exonic-overlap partner; everything else stays alone."""
    if mode == "union":
        return overlap_groups(list(range(len(recs))), recs)
    per = {n: overlap_groups([i for i, r in enumerate(recs) if r["ann"] == n], recs) for n in (A, Bn)}
    blocks = collections.defaultdict(list)
    for n in (A, Bn):
        for k, members in enumerate(per[n]):
            for s, e in gp.merge([b for i in members for b in recs[i]["exons"]]):
                blocks[recs[members[0]]["chrom"]].append((s, e, n, k))
    ov = collections.Counter()
    for bl in blocks.values():
        bl.sort()
        active = []
        for s, e, n, k in bl:
            active = [x for x in active if x[1] > s]
            for s2, e2, n2, k2 in active:
                if n2 != n:
                    key = (k, k2) if n == A else (k2, k)
                    ov[key] += min(e, e2) - s
            active.append((s, e, n, k))
    best_a, best_b = {}, {}
    for (ka, kb), bp in ov.items():
        if bp > best_a.get(ka, (0, None))[0]:
            best_a[ka] = (bp, kb)
        if bp > best_b.get(kb, (0, None))[0]:
            best_b[kb] = (bp, ka)
    matched_b, groups = set(), []
    for ka, members in enumerate(per[A]):
        kb = best_a.get(ka, (0, None))[1]
        if kb is not None and best_b.get(kb, (0, None))[1] == ka:
            groups.append(members + per[Bn][kb])
            matched_b.add(kb)
        else:
            groups.append(members)
    groups += [m for kb, m in enumerate(per[Bn]) if kb not in matched_b]
    return groups


def cmd_build(a):
    import mcl_port
    os.makedirs(a.out, exist_ok=True)
    contigs = set(a.contigs.split(","))
    genome = pysam.FastaFile(a.genome)
    anns = [x.split("=", 1) for x in a.ann]
    assert len(anns) == 2, "exactly two --ann"
    (A, specA), (Bn, specB) = anns
    recs = [r for name, spec in anns for r in load_records(name, spec, contigs)]
    groups = joint_groups(recs, A, Bn, a.joint_loci)
    loci, rec_locus = [], {}
    for members in groups:
        rs = [recs[i] for i in members]
        loci.append({"chrom": rs[0]["chrom"], "start": min(r["start"] for r in rs), "end": max(r["end"] for r in rs),
                     "exons": gp.merge([b for r in rs for b in r["exons"]]), "members": members,
                     "coding": any(r["cds"] for r in rs), "has": {n: any(r["ann"] == n for r in rs) for n in (A, Bn)},
                     "names": sorted({r["name"] for r in rs})})
    loci.sort(key=lambda l: (l["chrom"], l["start"], l["end"]))
    for k, l in enumerate(loci):
        l["id"] = k
        for i in l["members"]:
            rec_locus[(recs[i]["ann"], f"{recs[i]['chrom']}:{recs[i]['start'] + 1}-{recs[i]['end']}")] = k
    # E1 graph edges of each annotation, on joint loci
    E = {}
    unmapped = collections.Counter()
    for name, spec in anns:
        tag = spec.split(":", 1)[1]
        E[name] = {}
        for line in open(tag + ".graph.tsv"):
            x, y, w = line.rstrip("\n").split("\t")
            if x.rsplit(":", 1)[0] not in contigs or y.rsplit(":", 1)[0] not in contigs:
                continue
            u, v = rec_locus.get((name, x)), rec_locus.get((name, y))
            if u is None or v is None:
                unmapped[name] += 1
                continue
            if u != v:
                k = (min(u, v), max(u, v))
                E[name][k] = max(E[name].get(k, 0.0), float(w))
    agreed = set(E[A]) & set(E[Bn])
    disputed = set(E[A]) ^ set(E[Bn])
    # HGNC hard negatives: coding loci sharing a gene group, both annotated in both, no edge in either graph
    groups_of = collections.defaultdict(set)
    for r in csv.DictReader(open(a.hgnc), delimiter="\t"):
        for g in r["gene_group_id"].split("|"):
            if g:
                groups_of[r["symbol"]].add(g)
    grp_loci = collections.defaultdict(set)
    for l in loci:
        if l["coding"] and l["has"][A] and l["has"][Bn]:
            for nm in l["names"]:
                for g in groups_of.get(nm, ()):
                    grp_loci[g].add(l["id"])
    hard_neg = set()
    for ids in grp_loci.values():
        for p in itertools.combinations(sorted(ids), 2):
            if p not in E[A] and p not in E[Bn]:
                hard_neg.add(p)
    involved = sorted({x for p in (agreed | disputed | hard_neg) for x in p})
    X = blast_evidence(a, loci, involved, genome)
    sd = load_sedef(a.sedef, contigs)

    def ev(p):
        u, v = p
        return ((u, v) in X or (v, u) in X), (sd_evidence(sd, loci[u], loci[v]) or sd_evidence(sd, loci[v], loci[u]))
    evd = {p: ev(p) for p in sorted(agreed | disputed | hard_neg)}
    w = lambda p: sum(E[n][p] for n in (A, Bn) if p in E[n]) / sum(1 for n in (A, Bn) if p in E[n])
    strict = {p: w(p) for p in agreed | {p for p in disputed if any(evd[p])}}
    permissive = {p: w(p) for p in agreed | disputed}
    label = {}
    for tag, graph in (("strict", strict), ("permissive", permissive)):
        for k, c in enumerate(mcl_port.mcl(graph, a.inflation, a.prune)):
            for x in c:
                label[(tag, x)] = (k, len(c))

    def co(tag, u, v):
        lu, lv = label.get((tag, u)), label.get((tag, v))
        return lu is not None and lv is not None and lu[0] == lv[0] and lu[1] >= 2
    cand = set()
    for tag in ("strict", "permissive"):
        cl = collections.defaultdict(list)
        for (t, x), (k, n) in label.items():
            if t == tag and n >= 2:
                cl[k].append(x)
        for ms in cl.values():
            cand.update(itertools.combinations(sorted(ms), 2))
    status = {p: ("TRUE" if co("strict", *p) and co("permissive", *p) else "UNSCORED") for p in cand}
    with open(f"{a.out}/loci.tsv", "w") as fh:
        fh.write(f"locus\tchrom\tstart\tend\tcoding\thas_{A}\thas_{Bn}\tnames\texons\n")
        for l in loci:
            fh.write(f"L{l['id']}\t{l['chrom']}\t{l['start']}\t{l['end']}\t{int(l['coding'])}\t{int(l['has'][A])}\t"
                     f"{int(l['has'][Bn])}\t{','.join(l['names'])[:500]}\t{','.join(f'{x}-{y}' for x, y in l['exons'])}\n")
    with open(f"{a.out}/edges.tsv", "w") as fh:
        fh.write(f"u\tv\tin_{A}\tin_{Bn}\tblast\tsd\tclass\n")
        for p in sorted(evd):
            cls = "agreed" if p in agreed else ("disputed" if p in disputed else "hard_negative")
            fh.write(f"L{p[0]}\tL{p[1]}\t{int(p in E[A])}\t{int(p in E[Bn])}\t{int(evd[p][0])}\t{int(evd[p][1])}\t{cls}\n")
    with open(f"{a.out}/pairs.tsv", "w") as fh:
        fh.write("u\tv\tstatus\n")
        for p in sorted(status):
            fh.write(f"L{p[0]}\tL{p[1]}\t{status[p]}\n")
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
    n_true = sum(1 for s in status.values() if s == "TRUE")
    n_uns = len(status) - n_true
    sens = sum(1 for p in agreed if any(evd[p])) / max(1, len(agreed))
    fpr = sum(1 for p in hard_neg if any(evd[p])) / max(1, len(hard_neg))
    dis_kept = sum(1 for p in disputed if any(evd[p]))
    print(f"records {len(recs)}; joint loci {len(loci)}; graph edges not mapped to a locus {dict(unmapped)}")
    print(f"edges {A} {len(E[A])}, {Bn} {len(E[Bn])}; agreed {len(agreed)}; disputed {len(disputed)} "
          f"(with evidence {dis_kept}); hard negatives {len(hard_neg)}")
    print(f"TRUE pairs {n_true}; UNSCORED {n_uns}; truth clusters {len(comps)} (largest {max(map(len, comps.values())) if comps else 0})")
    ok = sens >= 0.80 and fpr <= 0.20 and n_uns <= 0.50 * (n_true + n_uns)
    print(f"GATE: evidence sensitivity on agreed edges {sens:.3f} (>= 0.80); hard-negative evidence rate {fpr:.3f} "
          f"(<= 0.20); unscored share {n_uns / max(1, n_true + n_uns):.3f} (<= 0.50) -> {'VALID' if ok else 'NOT VALID'}")


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
    # pairs are scored among truth loci only (the AG / `rna_truth.py` convention): span overlap cannot place a nested gene
    # inside a family member's intron, so assigning every joint locus manufactures false pairs (§6kn)
    order = tl
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
    p.add_argument("--blast-bin", default="/home/juanfra/miniforge3/envs/blast/bin")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--inflation", type=float, default=2.8)
    p.add_argument("--joint-loci", choices=("union", "matched"), default="matched")
    p.add_argument("--prune", type=float, default=1e-9)
    p = sub.add_parser("score")
    p.add_argument("--truth", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--expr")
    p.add_argument("catalogs", nargs="+")
    a = ap.parse_args()
    {"build": cmd_build, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
