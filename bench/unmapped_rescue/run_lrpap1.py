#!/usr/bin/env python3
"""LRPAP1 worked example (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 8).

    run_lrpap1.py cluster <tag>    # tag = full | half: net = LRPAP1 reads with maternal de > delta; cluster (frozen edge rule), abPOA consensus
    run_lrpap1.py place <tag>      # align each consensus to the maternal and paternal genome (-N 50), per-locus identity x coverage
    run_lrpap1.py augment <tag>    # realign all LRPAP1 reads to the consensus sequences, combine with the maternal primaries, report moves
Writes W/lrpap1/<tag>/. 'half' builds the consensus from a seeded half (seed 1) of the net; the other half and all well-placed reads are the test."""
import csv
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import augment as G  # noqa: E402
import graph as GR  # noqa: E402
import lrpap1 as L  # noqa: E402
import run_augment as RA  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402
from o3_maternal import common as C  # noqa: E402

M = "/mnt/linuxdisk/tmp/o3_mat"
OUT = "/mnt/linuxdisk/tmp/o3_rescue/lrpap1"
DELTA = C.DELTA
MM2 = RA.MM2
PY = "/home/juanfra/miniforge3/bin/python"


def loci_rows():
    return list(csv.DictReader(open(f"{M}/mat/truth/lrpap1_loci.tsv"), delimiter="\t"))


def bam_primaries(h):
    txt = subprocess.run(f"samtools view -F 2308 {M}/map/R_LRP.{h}.bam", shell=True, capture_output=True, text=True, check=True).stdout
    prim = RA.primaries(txt.splitlines())
    al = C.alias()
    for r in prim.values():
        r["ref"] = C.accession(r["ref"], al)
    return prim


def setup(tag):
    d = f"{OUT}/{tag}"
    os.makedirs(d, exist_ok=True)
    rows = loci_rows()
    mat, pat = bam_primaries("mat"), bam_primaries("pat")
    ivp = L.intervals(rows, "pat")
    reads = SD.read_fa(f"{M}/reads/R_LRP.fa")
    truth = {}
    for n in reads:
        r = pat.get(n)
        truth[n] = L.group(L.locus_at(ivp, r["ref"], r["start"], r["end"]).split("|")[0]) if r and L.locus_at(ivp, r["ref"], r["start"], r["end"]) else None
    net = sorted(n for n in reads if mat[n]["de"] > DELTA)
    build = list(net)
    if tag == "half":
        random.Random(1).shuffle(build)
        build = sorted(build[:len(net) // 2])
    return d, reads, mat, pat, truth, net, build


def cluster(tag):
    d, reads, mat, pat, truth, net, build = setup(tag)
    SD.write_fa(f"{d}/pool.fa", reads, build)
    json.dump(dict(net=net, build=build), open(f"{d}/sets.json", "w"))
    lens = {n: len(reads[n]) for n in build}
    seqs = {n: reads[n] for n in build}
    os.makedirs(f"{d}/rounds", exist_ok=True)
    try:
        comp = SD.run_rounds(lens, SD.minimap_map_fn(seqs, f"{d}/rounds", DELTA, 0.5), n_seeds=2000, min_size=3, max_rounds=6)
    except SD.Pause as e:
        print(f"round {e} finished; run again")
        sys.exit(75)
    cl = GR.clusters(comp, 3)
    maj = {c: S.majority(rs, truth) for c, rs in cl.items()}
    with open(f"{d}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\tmajority\n")
        for c, rs in cl.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\t{maj[c] or ''}\n")
    subprocess.run([PY, f"{HERE}/consensus.py", f"{d}/pool.fa", f"{d}/clusters.tsv", f"{d}/cons.fa"], check=True)
    pur = {c: sum(truth[r] == maj[c] for r in rs) / len(rs) for c, rs in cl.items()}
    print(f"{tag}: net {len(net)} (maternal de > {DELTA}), built from {len(build)}, {len(cl)} clusters covering {sum(len(r) for r in cl.values())} reads")
    for c, rs in sorted(cl.items(), key=lambda x: -len(x[1])):
        print(f"  cluster {c}: n={len(rs)} majority {maj[c]} purity {pur[c]:.3f}")


def place(tag):
    d = f"{OUT}/{tag}"
    rows = loci_rows()
    al = C.alias()
    names = [ln[1:].strip().split()[0] for ln in open(f"{d}/cons.fa") if ln[0] == ">"]
    res = {n: {} for n in names}
    for h in ("mat", "pat"):
        sam = f"{d}/cons.{h}.sam"
        if not os.path.exists(sam + ".done"):
            subprocess.run(f"minimap2 {MM2} -t 4 {C.HAP_IDX.format(h)} {d}/cons.fa > {sam}", shell=True, check=True)
            open(sam + ".done", "w").write("ok")
        iv = L.intervals(rows, h)
        recs = {n: [] for n in names}
        for ln in open(sam):
            if ln[0] == "@":
                continue
            f = ln.rstrip("\n").split("\t")
            if int(f[1]) & 4 or f[2] == "*":
                continue
            rlen = sum(int(n) for n, op in RA.CIG.findall(f[5]) if op in "M=XIS")
            span = sum(int(n) for n, op in RA.CIG.findall(f[5]) if op in "MDN=X")
            recs[f[0]].append(dict(ref=C.accession(f[2], al), start=int(f[3]) - 1, end=int(f[3]) - 1 + span, de=float(RA.tag(f, "de", 1.0)), qcov=RA.qcov(f[5], rlen)))
        for n in names:
            res[n][h] = dict(sorted(L.placements(recs[n], iv).items(), key=lambda kv: -kv[1]))
    json.dump(res, open(f"{d}/placement.json", "w"), indent=1)
    cm = G_cluster_majority(d)
    print(f"{'cluster':32s} {'maj':8s} {'best pat locus (idxcov)':28s} {'best mat site (idxcov)':28s} O3")
    for n in names:
        bp = next(iter(res[n]["pat"].items()), (None, None))
        bm = next(iter(res[n]["mat"].items()), (None, None))
        print(f"{n:32s} {str(cm.get(n)):8s} {bp[0]!s:12s} {(bp[1] or 0):.4f}   {bm[0]!s:14s} {(bm[1] or 0):.4f}   {L.o3_flag(bm[1], bp[1])}")


def G_cluster_majority(d):
    return RA.cluster_family(f"{d}/clusters.tsv", f"{d}/cons.fa")


def augment(tag):
    d, reads, mat, pat, truth, net, build = setup(tag)
    fq = f"{d}/reads.fa"
    SD.write_fa(fq, reads, sorted(reads))
    sam = f"{d}/aug.sam"
    if not os.path.exists(sam + ".done"):
        subprocess.run(f"minimap2 {MM2} -t 4 {d}/cons.fa {fq} > {sam}", shell=True, check=True)
        open(sam + ".done", "w").write("ok")
    cons = RA.primaries(open(sam))
    cf = G_cluster_majority(d)
    netset, buildset = set(net), set(build)
    rows = []
    for n in sorted(reads):
        if n in buildset:
            cls = "net_build"
        elif n in netset:
            cls = "net_heldout"
        else:
            cls = "placed"
        rows.append((cls, truth[n], mat[n], cons.get(n), cons[n]["ref"] if n in cons else None))
    res = dict(tag=tag, rule="divergence", net=len(net), built_from=len(build), consensus=len(cf), reads=len(rows), moves=G.move_metrics(rows, cf),
               moves_score_rule=G.move_metrics(rows, cf, rule="score"))
    # placed reads by true copy: which copies attract false moves
    fm = {}
    for cls, fam, g, c, cname in rows:
        if cls == "placed":
            src, _ = G.new_primary(g, c)
            o = fm.setdefault(str(fam), [0, 0])
            o[0] += 1
            o[1] += src == "consensus"
    res["placed_by_true_copy"] = {k: dict(n=v[0], moved=v[1]) for k, v in sorted(fm.items())}
    json.dump(res, open(f"{d}/augment.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"cluster": cluster, "place": place, "augment": augment}[sys.argv[1]](sys.argv[2])
