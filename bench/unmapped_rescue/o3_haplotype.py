#!/usr/bin/env python3
"""The three-way O3 call on the maternal-reference study (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 26). Miniforge python.

    o3_haplotype.py <mat|pat> prep     # net reads of the 34 families + LRPAP1 on the reference haplotype -> W/o3hap_<ref>/{pool.fa, labels.tsv}
    (then:  run_bed.py o3hap_<ref> --proper --tag chain ; run_chain.py o3hap_<ref> chain chain3 ; consensus.py pool chain3/clusters.tsv chain3/cons.fa)
    o3_haplotype.py <mat|pat> eval     # consensus aligned to both haplotypes, COPY / ALLELE / PRESENT / None per cluster, table by majority label"""
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import flagmetric as FM  # noqa: E402
import libsig  # noqa: E402
import reciprocal as RCP  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402
from o3_maternal import common as C  # noqa: E402

M = "/mnt/linuxdisk/tmp/o3_mat"
W = "/mnt/linuxdisk/tmp/o3_rescue"
DELTA = 0.00958
RCT = str.maketrans("ACGTacgt", "TGCAtgca")


def view(ref):
    return subprocess.Popen(f"samtools view -F 2304 {M}/map/reads.{ref}.all.bam", shell=True, stdout=subprocess.PIPE, text=True)


def prep(ref):
    d = f"{W}/o3hap_{ref}"
    os.makedirs(d, exist_ok=True)
    lab = {r["read"]: r["label"] for r in csv.DictReader(open(f"{M}/{ref}/truth/labels.tsv"), delimiter="\t")}
    seqs, n_map, n_net = {}, 0, 0
    for ln in view(ref).stdout:
        f = ln.rstrip("\n").split("\t")
        if f[0] not in lab:
            continue
        n_map += 1
        flag = int(f[1])
        de = next((float(t[5:]) for t in f[11:] if t.startswith("de:f:")), 1.0)
        if flag & 4 or de > DELTA:
            s = f[9]
            seqs[f[0]] = s.translate(RCT)[::-1] if flag & 16 else s
            n_net += 1
    SD.write_fa(f"{d}/pool.fa", seqs, sorted(seqs))
    with open(f"{d}/labels.tsv", "w") as o:
        o.write("read\tfamily\tcopy\trole\n")
        for n in sorted(seqs):
            l = lab[n]
            keep = l not in ("shared", "ambiguous")
            o.write(f"{n}\t{l if keep else ''}\t{l if keep else ''}\t{'D' if keep else 'bg'}\n")
    shared = {n for n, l in lab.items() if l == "shared"}
    sig = libsig.signature(view(ref).stdout, keep=lambda n: n in shared)
    ok, p = libsig.gate(sig)
    json.dump(dict(gate=ok, p=p, signature=sig), open(f"{d}/gate.json", "w"), indent=1)
    c = collections.Counter(lab[n] for n in seqs)
    print(f"{ref}: labelled reads with a primary {n_map}; net (unmapped or de > {DELTA}) {n_net}; by label: {dict(c.most_common(12))}; gate {'OPEN' if ok else 'closed'} (pure clips {sig['pure']}, clean reads {sig['reads']})")


def best_hits(paf):
    b = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        m = int(f[9])
        if f[0] not in b or m > b[f[0]][0]:
            b[f[0]] = (m, m / max(1, int(f[10])), int(f[2]), int(f[3]), int(f[1]))
    return b


def evaluate(ref):
    truth_h = "pat" if ref == "mat" else "mat"
    d = f"{W}/o3hap_{ref}"
    out = f"{d}/chain3"
    gate = json.load(open(f"{d}/gate.json"))["gate"]
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{out}/cons.fa").items()}
    fa = f"{out}/cons.o3.fa"
    SD.write_fa(fa, cons, list(cons))
    hits = {}
    for which, h in (("R", ref), ("T", truth_h)):
        paf = f"{out}/cons.o3.{which}.paf"
        if not os.path.exists(paf + ".done"):
            subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {C.HAP_IDX.format(h)} {fa} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
            open(paf + ".done", "w").write("ok")
        hits[which] = best_hits(paf)

    def sc(which, k):
        h = hits[which].get(k)
        return FM.gate_aware(h[1], h[2], h[3], h[4], cons[k][:h[2]], gate) if h else None
    lab = {r["read"]: r["label"] for r in csv.DictReader(open(f"{M}/{ref}/truth/labels.tsv"), delimiter="\t")}
    lrp = {ln.split("\t")[0] for ln in open(f"{M}/reads/R_LRP.names.tsv")}
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{out}/clusters.tsv"), delimiter="\t"):
        cl["cl" + r["cluster"]].append(r["read"])
    rows = []
    for k, rs in cl.items():
        labs = collections.Counter(lab[r] for r in rs)
        maj, nmaj = labs.most_common(1)[0]
        is_lrp = sum(r in lrp for r in rs) > len(rs) / 2
        R, T = sc("R", k), sc("T", k)
        rows.append(dict(k=k, reads=len(rs), label=("LRPAP1" if is_lrp else maj), label_share=round(nmaj / len(rs), 2), R=R, T=T, cls=FM.o3_class(R, T), registered=FM.registered_flag(R, T)))
    json.dump(rows, open(f"{out}/o3_classes.json", "w"), indent=1)
    tab = collections.defaultdict(collections.Counter)
    for r in rows:
        g = r["label"] if r["label"] not in ("shared", "ambiguous") else r["label"]
        tab[g][str(r["cls"])] += 1
        tab[g]["clusters"] += 1
        tab[g]["registered_flag"] += r["registered"]
    print(f"reference {ref}, truth {truth_h}, gate {'OPEN' if gate else 'closed'}; clusters {len(rows)} (>= 3 reads)")
    print(f"{'majority label':16s} {'clusters':>8s} {'COPY':>5s} {'ALLELE':>6s} {'PRESENT':>7s} {'none':>5s} {'registered flag':>15s}")
    for g, c in sorted(tab.items(), key=lambda kv: (kv[0] in ("shared", "ambiguous"), kv[0])):
        print(f"{g:16s} {c['clusters']:>8d} {c['COPY']:>5d} {c['ALLELE']:>6d} {c['PRESENT']:>7d} {c['None']:>5d} {c['registered_flag']:>15d}")
    print("COPY calls:")
    for r in rows:
        if r["cls"] == "COPY":
            print(f"  {r['k'][-14:]:14s} reads {r['reads']:4d} label {r['label']} ({r['label_share']}) R={r['R'] if r['R'] is None else round(r['R'], 4)} T={round(r['T'], 4)}")


def best_lines(paf):
    b = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        if f[0] not in b or int(f[9]) > int(b[f[0]][9]):
            b[f[0]] = f
    return b


def reciprocal(ref):
    """Amendment 27: for every COPY cluster, map the reference locus' transcript back to the truth haplotype; reciprocal -> diverged ortholog (ALLELE), else COPY"""
    import pysam
    truth_h = "pat" if ref == "mat" else "mat"
    out = f"{W}/o3hap_{ref}/chain3"
    rows = json.load(open(f"{out}/o3_classes.json"))
    Rl, Tl = best_lines(f"{out}/cons.o3.R.paf"), best_lines(f"{out}/cons.o3.T.paf")
    fa = pysam.FastaFile(f"/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{ref}.fa")
    al = C.alias()
    xt = {}
    for r in rows:
        if r["cls"] != "COPY" or r["k"] not in Rl:
            continue
        f = Rl[r["k"]]
        cg = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), None)
        segs = RCP.target_segments(cg, int(f[7]))
        acc = C.accession(f[5], al)
        seq = "".join(fa.fetch(acc, s, e) for s, e in segs).upper()
        xt[r["k"]] = RCP.orient(seq, f[4])
    SD.write_fa(f"{out}/xtranscripts.fa", xt, list(xt))
    paf = f"{out}/xtranscripts.truth.paf"
    subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {C.HAP_IDX.format(truth_h)} {out}/xtranscripts.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
    Y1 = best_lines(paf)
    res = []
    for r in rows:
        if r["k"] not in xt:
            continue
        t0 = Tl[r["k"]]
        y0 = (t0[5], int(t0[7]), int(t0[8]))
        y1l = Y1.get(r["k"])
        y1 = (y1l[5], int(y1l[7]), int(y1l[8])) if y1l else None
        rec = RCP.same_locus(y1, y0)
        res.append(dict(k=r["k"], reads=r["reads"], label=r["label"], R=r["R"], T=r["T"], xlen=len(xt[r["k"]]), Y0=y0, Y1=y1, Y1_identity=(int(y1l[9]) / max(1, int(y1l[10])) if y1l else None),
                        reciprocal=rec, final="ALLELE (diverged ortholog)" if rec else "COPY"))
    json.dump(res, open(f"{out}/reciprocal.json", "w"), indent=1)
    print(f"reference {ref}, truth {truth_h}: {len(res)} COPY clusters tested")
    for x in res:
        print(f"  {x['final']:26s} {x['label']:11s} reads {x['reads']:4d} R={x['R']:.4f}  consensus' truth hit {x['Y0'][0]}:{x['Y0'][1]}  X-transcript ({x['xlen']} bp) back-maps to "
              f"{x['Y1'][0] + ':' + str(x['Y1'][1]) if x['Y1'] else 'nothing'} (identity {x['Y1_identity']:.4f})" if x['Y1'] else f"  {x['final']:26s} {x['label']:11s} reads {x['reads']:4d} X-transcript back-maps to nothing")


if __name__ == "__main__":
    {"prep": prep, "eval": evaluate, "reciprocal": reciprocal}[sys.argv[2]](sys.argv[1])
