#!/usr/bin/env python3
"""rcf611_graph_vs_linear.py — does GRAPH alignment assign reads to copies better than LINEAR?

Decisive test of the user's hypothesis on RCF_611 (4 co-located paralogs w/ tandem-repeat CNV).
Simulate KNOWN-ORIGIN reads from the 4 copy cDNAs (full-length + fragments + HiFi error), then assign:
  LINEAR-primary    : minimap2 reads -> 4 copy cDNAs (-N50 -p0.1); copy = PRIMARY target (the pipeline).
  LINEAR-ASdecisive : best alignment-score copy, with a decisiveness gate (AS gap >= 20), else ambiguous.
  GRAPH             : vg giraffe reads -> vg msga graph of the 4 copies; copy = which copy's UNIQUE
                      graph nodes the read traverses (decisive if one copy), else ambiguous.
Compare per-copy assignment ACCURACY against the simulated truth.

Run: /home/juanfra/miniforge3/bin/python bench/rcf611_graph_vs_linear.py
"""
import collections
import json
import os
import random
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
import family_def_vg_coherence as V
from family_def_vg_coherence import build_cdna, build_families

VG = "/home/juanfra/miniforge3/bin/vg"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
FAM = "RCF_611"
N_PER_COPY = 200
ERR = 0.012           # ~1.2% HiFi-ish error
random.seed(7)


def mutate(seq):
    out = []
    for b in seq:
        r = random.random()
        if r < ERR * 0.6:
            out.append(random.choice("ACGT"))           # substitution
        elif r < ERR * 0.8:
            continue                                      # deletion
        elif r < ERR:
            out.append(b); out.append(random.choice("ACGT"))  # insertion
        else:
            out.append(b)
    return "".join(out)


def simulate(copies):
    reads = []
    for cp, seq in copies.items():
        for i in range(N_PER_COPY):
            if random.random() < 0.4 or len(seq) < 250:   # full-length
                sub = seq
            else:                                          # fragment
                L = random.randint(200, len(seq))
                st = random.randint(0, len(seq) - L)
                sub = seq[st:st + L]
            reads.append((f"{cp}|{i}", cp, mutate(sub)))
    random.shuffle(reads)
    return reads


def linear_assign(reads, copies, tmp):
    cref = os.path.join(tmp, "copies.fa")
    with open(cref, "w") as o:
        for g, s in copies.items():
            o.write(f">{g}\n{s}\n")
    rq = os.path.join(tmp, "reads.fa")
    with open(rq, "w") as o:
        for n, _, s in reads:
            o.write(f">{n}\n{s}\n")
    res = subprocess.run([MM2, "-ax", "map-hifi", "-N", "50", "-p", "0.1", cref, rq],
                         capture_output=True, text=True)
    primary = {}; mapq = {}; asmap = collections.defaultdict(dict)
    for line in res.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x4:
            continue
        AS = next((int(t[5:]) for t in f[11:] if t.startswith("AS:i:")), -9)
        asmap[f[0]][f[2]] = max(asmap[f[0]].get(f[2], -9), AS)
        if not (flag & 0x100) and not (flag & 0x800):
            primary[f[0]] = f[2]; mapq[f[0]] = int(f[4])
    prim = {n: (primary.get(n), mapq.get(n, -1)) for n, _, _ in reads}
    # AS-decisive
    asdec = {}
    for n, _, _ in reads:
        d = asmap.get(n, {})
        if not d:
            asdec[n] = None; continue
        v = sorted(d.values(), reverse=True)
        top = max(d, key=d.get)
        asdec[n] = top if (len(v) == 1 or v[0] - v[1] >= 20) else None
    return prim, asdec


def graph_assign(reads, copies, tmp):
    cref = os.path.join(tmp, "g_copies.fa")
    with open(cref, "w") as o:
        for g, s in copies.items():
            o.write(f">{g}\n{s}\n")
    g = os.path.join(tmp, "g.vg")
    with open(g, "wb") as out:
        subprocess.run([VG, "msga", "-f", cref, "-t", "4"], stdout=out, stderr=subprocess.DEVNULL, check=True)
    gfa = subprocess.run([VG, "convert", g, "-f"], capture_output=True, text=True).stdout
    # copy -> node set ; copy-unique nodes
    cpaths = {}
    for line in gfa.splitlines():
        f = line.split("\t")
        if f and f[0] == "P":
            cpaths[f[1]] = {x[:-1] for x in f[2].split(",")}
    node_owner = collections.Counter()
    for ns in cpaths.values():
        node_owner.update(ns)
    unique = {cp: {n for n in ns if node_owner[n] == 1} for cp, ns in cpaths.items()}
    # build giraffe indexes
    gfap = os.path.join(tmp, "g.gfa"); open(gfap, "w").write(gfa)
    gbz = os.path.join(tmp, "g.gbz"); dist = os.path.join(tmp, "g.dist"); mn = os.path.join(tmp, "g.min")
    subprocess.run([VG, "gbwt", "-G", gfap, "--gbz-format", "-g", gbz], stderr=subprocess.DEVNULL, check=True)
    subprocess.run([VG, "index", "-j", dist, gbz], stderr=subprocess.DEVNULL, check=True)
    subprocess.run([VG, "minimizer", "-d", dist, "-o", mn, gbz], stderr=subprocess.DEVNULL, check=True)
    rq = os.path.join(tmp, "g_reads.fq")
    with open(rq, "w") as o:
        for n, _, s in reads:
            o.write(f"@{n}\n{s}\n+\n{'I'*len(s)}\n")
    gam = os.path.join(tmp, "aln.gam")
    with open(gam, "wb") as out:
        subprocess.run([VG, "giraffe", "-Z", gbz, "-m", mn, "-d", dist, "-f", rq],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    js = subprocess.run([VG, "view", "-aj", gam], capture_output=True, text=True).stdout
    assign = {}
    for line in js.splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rn = d.get("name")
        nodes = {str(m["position"]["node_id"]) for m in d.get("path", {}).get("mapping", [])}
        hits = {cp: len(nodes & u) for cp, u in unique.items()}
        v = sorted(hits.values(), reverse=True)
        top = max(hits, key=hits.get) if hits else None
        assign[rn] = top if (top and hits[top] >= 1 and (len(v) < 2 or v[0] > v[1])) else None
    return {n: assign.get(n) for n, _, _ in reads}


def main():
    V.EXON = json.load(open("/home/juanfra/winloci_scratch/gene_exon_index.json"))
    fa = pysam.FastaFile("/home/juanfra/winloci_scratch/GGO.fasta")
    members = dict(build_families())[FAM]
    copies = {g: build_cdna(g, fa) for g in members if build_cdna(g, fa)}
    fa.close()
    tmp = tempfile.mkdtemp(prefix="rcf611_")
    reads = simulate(copies)
    truth = {n: cp for n, cp, _ in reads}
    print(f"[{FAM}] {len(copies)} copies, {len(reads)} simulated reads ({N_PER_COPY}/copy, err~{ERR})")

    prim, asdec = linear_assign(reads, copies, tmp)
    graph = graph_assign(reads, copies, tmp)

    def score(pred, name):
        correct = wrong = amb = 0
        for n in truth:
            p = pred[n][0] if isinstance(pred[n], tuple) else pred[n]
            if p is None:
                amb += 1
            elif p == truth[n]:
                correct += 1
            else:
                wrong += 1
        tot = len(truth)
        print(f"  {name:20} correct {correct:4d} ({100*correct/tot:.1f}%)  "
              f"WRONG {wrong:4d} ({100*wrong/tot:.1f}%)  ambiguous {amb:4d} ({100*amb/tot:.1f}%)")
        return correct, wrong, amb

    print("\n=== copy-assignment accuracy vs simulated truth ===")
    score(prim, "LINEAR-primary")
    score(asdec, "LINEAR-ASdecisive")
    score(graph, "GRAPH (vg giraffe)")

    # head-to-head: linear-primary WRONG but graph RIGHT (the rescue), and reverse
    lp_wrong_g_right = sum(1 for n in truth if prim[n][0] != truth[n] and prim[n][0] is not None and graph[n] == truth[n])
    g_wrong_lp_right = sum(1 for n in truth if graph[n] != truth[n] and graph[n] is not None and prim[n][0] == truth[n])
    print(f"\n  head-to-head: LINEAR-primary WRONG but GRAPH correct (graph rescues): {lp_wrong_g_right}")
    print(f"                GRAPH wrong but LINEAR-primary correct (graph regresses): {g_wrong_lp_right}")
    # on the linear MAPQ0 subset
    m0 = [n for n in truth if prim[n][1] == 0]
    if m0:
        gc = sum(1 for n in m0 if graph[n] == truth[n])
        pc = sum(1 for n in m0 if prim[n][0] == truth[n])
        print(f"  on LINEAR MAPQ0 reads ({len(m0)}): linear-primary correct {pc}, GRAPH correct {gc}")


if __name__ == "__main__":
    main()
