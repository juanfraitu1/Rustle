#!/usr/bin/env python3
"""Discovery of reference-absent sequence from reads that do not align to the reference (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 28). Miniforge python.

    discover.py <name> <reads.fa>      resumable (exit 75 = run again); writes W/discover_<name>/{clusters.tsv, cons.fa, classes.json}
Cluster (adopted rule: --proper edge, star step), abPOA consensus, align each consensus to the primary (R) and to the mother's (Tm) and the father's (Tp) assembly,
classify CONFIRMED / NOVEL / ALLELE-LIKE / PRESENT (flagmetric.discovery_class)."""
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import chain as CH  # noqa: E402
import flagmetric as FM  # noqa: E402
import graph as GR  # noqa: E402
import consensus_support as CS  # noqa: E402
import partition as PT  # noqa: E402
import rescue as RS  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402
from o3_maternal import common as C  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
M = "/mnt/linuxdisk/tmp/o3_mat"
DELTA = 0.00958
RCT = str.maketrans("ACGTacgt", "TGCAtgca")
PRIMARY_FA = "/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/genome.fa"
HAP_FA = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{}.fa"
WINDOW = 5000


def cluster(name, reads_fa):
    d = f"{W}/discover_{name}"
    os.makedirs(f"{d}/rounds", exist_ok=True)
    open(f"{d}/reads_path.txt", "w").write(reads_fa)
    seqs = SD.read_fa(reads_fa)
    lens = {n: len(s) for n, s in seqs.items()}
    try:
        comp = SD.run_rounds(lens, SD.minimap_map_fn(seqs, f"{d}/rounds", DELTA, 0.5, proper=True), n_seeds=2000, min_size=3, max_rounds=6)
    except SD.Pause as e:
        print(f"round {e} finished; run again")
        sys.exit(75)
    cl = {c: rs for c, rs in GR.clusters(comp, 3).items()}
    new = CH.refine(cl, seqs, CH.minimap_allvsall(f"{d}/chain_tmp"), DELTA)
    with open(f"{d}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\n")
        for c, rs in new.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\n")
    cons = {f"cl{c}|n={len(rs)}": PT.abpoa_consensus([seqs[r] for r in rs]) for c, rs in new.items()}
    SD.write_fa(f"{d}/cons.fa", cons, list(cons))
    print(f"{name}: {len(seqs)} reads, {len(cl)} components, {len(new)} clusters after the star step, {sum(len(v) for v in new.values())} reads clustered")
    return d, seqs, new


def best_records(paf):
    b = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        k = f[0].split("|")[0]
        if k not in b or int(f[9]) > int(b[k][9]):
            b[k] = f
    return b


def idcov_with_rescue(cons, f, fasta, to_fasta_name, gate, tmp):
    """gate-aware identity x coverage of a consensus on one assembly, with the unaligned ends of the best record rescued inside the locus window (Amendment 29)"""
    qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
    matches, blk = int(f[9]), int(f[10])
    segs = RS.end_segments(qs, qe, qlen)
    iv, pm, pb = [], 0, 0
    if segs:
        chrom, ts, te = to_fasta_name(f[5]), int(f[7]), int(f[8])
        lo, hi = max(0, ts - WINDOW), min(fasta.get_reference_length(chrom), te + WINDOW)
        iv, pm, pb = RS.rescued(segs, RS.map_pieces(fasta.fetch(chrom, lo, hi).upper(), segs, cons, tmp))
    aligned, ident = RS.combine(qlen, (qs, qe, matches, blk), (iv, pm, pb))
    lead = min([qs] + [a for a, _ in iv])
    return FM.gate_aware_total(ident, aligned, qlen, lead, cons[:lead], gate)


def classify(name, confirmable=()):
    import pysam
    d = f"{W}/discover_{name}"
    gate = json.load(open(f"{W}/o3hap_mat/gate.json"))["gate"]
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{d}/cons.fa").items()}
    al = C.alias()
    assemblies = {"R": (RP.GENOME_IDX, PRIMARY_FA, lambda n: n), "Tm": (C.HAP_IDX.format("mat"), HAP_FA.format("mat"), lambda n: C.accession(n, al)),
                  "Tp": (C.HAP_IDX.format("pat"), HAP_FA.format("pat"), lambda n: C.accession(n, al))}
    recs = {}
    for which, (idx, _fa, _nm) in assemblies.items():
        paf = f"{d}/cons.{which}.paf"
        if not os.path.exists(paf + ".done"):
            subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {idx} {d}/cons.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
            open(paf + ".done", "w").write("ok")
        recs[which] = best_records(paf)
    fas = {w: pysam.FastaFile(a[1]) for w, a in assemblies.items()}
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
        cl["cl" + r["cluster"]].append(r["read"])
    reads_fa = SD.read_fa(open(f"{d}/reads_path.txt").read().strip())
    rows = []
    for k, rs in cl.items():
        sc = {w: (idcov_with_rescue(cons[k], recs[w][k], fas[w], assemblies[w][2], gate, f"{d}/tmp_{w}") if k in recs[w] else None) for w in assemblies}
        med = CS.median_read_divergence(cons[k], {n: reads_fa[n] for n in sorted(rs)[:100]}, f"{d}/tmp_support")
        unsupported = med is None or med > DELTA
        cls = "UNSUPPORTED" if unsupported else FM.discovery_class(sc["R"], sc["Tm"], sc["Tp"])
        hm = recs["Tm"].get(k)
        rows.append(dict(k=k, reads=len(rs), length=len(cons[k]), R=sc["R"], Tm=sc["Tm"], Tp=sc["Tp"], median_read_divergence=med, cls=cls, mat_hit=(hm[5], int(hm[7]), int(hm[8])) if hm else None))
    json.dump(rows, open(f"{d}/classes.json", "w"), indent=1)
    by = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        by[r["cls"]][0] += 1
        by[r["cls"]][1] += r["reads"]
    print(f"{name}: gate {'OPEN' if gate else 'closed'}")
    for c in ("CONFIRMED", "NOVEL", "DIVERGED", "UNSUPPORTED", "ALLELE-LIKE", "PRESENT"):
        print(f"  {c:12s} clusters {by[c][0]:4d}  reads {by[c][1]:5d}")
    if confirmable:
        cls_of = {r_: next(x["cls"] for x in rows if x["k"] == k) for k, rs in cl.items() for r_ in rs}
        inn = sum(cls_of.get(r) == "CONFIRMED" for r in confirmable)
        print(f"  confirmable reads (map on the mother's or the father's assembly): {len(confirmable)}; in CONFIRMED clusters {inn} ({inn / len(confirmable):.0%}); in any cluster "
              f"{sum(r in cls_of for r in confirmable)}")
    big = sorted(rows, key=lambda x: -x["reads"])[:5]
    f3 = lambda v: None if v is None else round(v, 4)
    print("  largest clusters (reads, length, class, R, Tm, Tp, median read divergence):", [(x["reads"], x["length"], x["cls"], f3(x["R"]), f3(x["Tm"]), f3(x["Tp"]), f3(x["median_read_divergence"])) for x in big])
    for c in ("NOVEL", "DIVERGED", "UNSUPPORTED"):
        sel = sorted((x for x in rows if x["cls"] == c), key=lambda x: -x["reads"])
        print(f"  {c} clusters: {len(sel)}; reads {[x['reads'] for x in sel][:15]}; lengths {[x['length'] for x in sel][:15]}")


if __name__ == "__main__":
    name, reads_fa = sys.argv[1], sys.argv[2]
    conf = ()
    if name == "unm":
        import pysam
        s = set()
        for h in ("mat", "pat"):
            for rd in pysam.AlignmentFile(f"{M}/map/R_unm.{h}.bam", "rb").fetch(until_eof=True):
                if not rd.is_unmapped and not rd.is_secondary and not rd.is_supplementary:
                    s.add(rd.query_name)
        conf = sorted(s)
    os.makedirs(f"{W}/discover_{name}", exist_ok=True)
    cluster(name, reads_fa)
    classify(name, conf)
