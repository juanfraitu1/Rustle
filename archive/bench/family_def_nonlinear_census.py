#!/usr/bin/env python3
"""family_def_nonlinear_census.py — build a real `vg` graph per non-linear-flagged family and
classify the structure (tandem-repeat/CNV vs rearrangement vs artifact) + measure whether the
STRUCTURE distinguishes copies that SNP/indel PSVs cannot.

For each flagged family (tandem-repeat ∪ exon-shuffle from family_def_vg_coherence.json):
  1. build member cDNA FASTA (cap members for msga speed).
  2. `vg msga -f fam.fa` -> graph; `vg stats -A` (cyclic?) ; `vg view` -> GFA P-lines = per-member
     node-traversal order.
  3. per member: node-path, max node MULTIPLICITY (>=2 => the copy's sequence revisits a node =
     internal tandem repeat), distinct-node order.
  4. classify:
       CYCLE_CNV     : cyclic AND members differ in max-multiplicity  -> repeat copy-number varies
                       BETWEEN copies = a STRUCTURAL distinguishing feature SNP-PSVs miss.
       REPEAT_SHARED : cyclic but all members same multiplicity (repeat present, not distinguishing)
       REARRANGEMENT : acyclic but a member-pair shares >=3 nodes in a DIFFERENT order (exon shuffle)
       LINEAR        : acyclic + collinear (the coherence flag was a minimap2 multi-map artifact)
  5. artifact guard: complexity of the repeated nodes (low-complexity simple repeat vs real exon unit).
  6. distinguishing = does the structural feature take >=2 distinct values across members?

Run: /home/juanfra/miniforge3/bin/python bench/family_def_nonlinear_census.py
Out: bench/family_def_nonlinear_census.{tsv,json}
"""
import collections
import json
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import family_def_vg_coherence as V
from family_def_vg_coherence import build_cdna, build_families

GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
EXON_INDEX = "/home/juanfra/winloci_scratch/gene_exon_index.json"
VG = "/home/juanfra/miniforge3/bin/vg"
MAX_MEMBERS = 10        # msga is slow on big families; sample to this many (longest first)
MSGA_TIMEOUT = 180


def complexity(seq):
    """fraction of 2-mers that are NOT the single most-common 2-mer (1.0=complex, ~0=simple repeat)."""
    if len(seq) < 4:
        return 1.0
    c = collections.Counter(seq[i:i+2] for i in range(len(seq) - 1))
    return 1.0 - c.most_common(1)[0][1] / max(sum(c.values()), 1)


def census_family(args):
    fid, members = args
    tmp = tempfile.mkdtemp(prefix=f"nlc_{os.getpid()}_")
    fa = pysam.FastaFile(GENOME)
    seqs = {}
    for g in members:
        s = build_cdna(g, fa)
        if s:
            seqs[g] = s
    fa.close()
    if len(seqs) < 2:
        return None
    # cap to longest MAX_MEMBERS for msga tractability
    keep = sorted(seqs, key=lambda g: -len(seqs[g]))[:MAX_MEMBERS]
    fap = os.path.join(tmp, "f.fa")
    with open(fap, "w") as o:
        for g in keep:
            o.write(f">{g}\n{seqs[g]}\n")
    vgp = os.path.join(tmp, "g.vg")
    try:
        with open(vgp, "wb") as out:
            subprocess.run([VG, "msga", "-f", fap, "-t", "2"], stdout=out,
                           stderr=subprocess.DEVNULL, timeout=MSGA_TIMEOUT, check=True)
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
        return dict(family=fid, n=len(keep), status="msga_failed")
    stat = subprocess.run([VG, "stats", "-z", "-A", vgp], capture_output=True, text=True).stdout
    cyclic = "cyclic" in stat and "acyclic" not in stat.replace("cyclic", "", 1)
    nodes = edges = 0
    for line in stat.splitlines():
        if line.startswith("nodes"):
            nodes = int(line.split()[1])
        elif line.startswith("edges"):
            edges = int(line.split()[1])
    # GFA P-lines: per-member node order
    gfa = subprocess.run([VG, "view", vgp], capture_output=True, text=True).stdout
    seg = {}
    paths = {}
    for line in gfa.splitlines():
        f = line.split("\t")
        if f[0] == "S":
            seg[f[1]] = f[2]
        elif f[0] == "P":
            nodes_o = [x[:-1] for x in f[2].split(",")]   # strip orientation
            paths[f[1]] = nodes_o
    # per member: max node multiplicity (internal repeat copy number)
    mult = {}
    for m, ns in paths.items():
        c = collections.Counter(ns)
        mult[m] = max(c.values()) if c else 1
    multiset = set(mult.values())
    # rearrangement: a member-pair sharing >=3 nodes in different order
    rearr = False
    pl = list(paths.items())
    for i in range(len(pl)):
        for j in range(i + 1, len(pl)):
            a = [n for n in pl[i][1]]; b = [n for n in pl[j][1]]
            common = set(a) & set(b)
            if len(common) < 3:
                continue
            sa = [n for n in a if n in common]
            sb = [n for n in b if n in common]
            # order discordance: are the common nodes in the same relative order?
            if sa != sb and sa != sb[::-1]:
                # count inversions between the two orderings
                rank = {n: k for k, n in enumerate(sa)}
                seq = [rank[n] for n in sb if n in rank]
                disc = sum(1 for x in range(len(seq)) for y in range(x + 1, len(seq)) if seq[x] > seq[y])
                if disc > 0 and disc < len(seq) * (len(seq) - 1) / 2:   # not pure reverse
                    rearr = True
    # repeated-node complexity (artifact guard): nodes appearing >=2x in any path
    rep_nodes = set()
    for ns in paths.values():
        c = collections.Counter(ns)
        rep_nodes |= {n for n, k in c.items() if k >= 2}
    rep_len = sum(len(seg.get(n, "")) for n in rep_nodes)
    rep_cplx = (sum(complexity(seg[n]) * len(seg[n]) for n in rep_nodes if n in seg)
                / max(rep_len, 1)) if rep_nodes else 1.0
    max_mult = max(mult.values()) if mult else 1
    # classify
    if cyclic and len(multiset) >= 2:
        cls = "CYCLE_CNV"
    elif cyclic and max_mult >= 2:
        cls = "REPEAT_SHARED"
    elif rearr:
        cls = "REARRANGEMENT"
    else:
        cls = "LINEAR"
    distinguishing = (cls == "CYCLE_CNV") or (cls == "REARRANGEMENT")
    for fn in ("f.fa", "g.vg"):
        try: os.remove(os.path.join(tmp, fn))
        except OSError: pass
    try: os.rmdir(tmp)
    except OSError: pass
    return dict(family=fid, n=len(keep), status="ok", cyclic=cyclic, nodes=nodes, edges=edges,
                max_mult=max_mult, mult_values=sorted(multiset), rearrangement=rearr,
                rep_len=rep_len, rep_complexity=round(rep_cplx, 2), cls=cls,
                distinguishing=distinguishing)


def main():
    V.EXON = json.load(open(EXON_INDEX))
    fams = dict(build_families())
    flag = json.load(open(os.path.join(HERE, "family_def_vg_coherence.json")))["flagged"]
    flagged = sorted(set(flag.get("tandem_repeat", []) + flag.get("exon_shuffle", [])))
    jobs = [(fid, fams[fid]) for fid in flagged if fid in fams]
    print(f"[census] {len(jobs)} non-linear-flagged families -> building vg graphs ...")
    recs = []
    with Pool(3) as pool:
        for r in pool.imap_unordered(census_family, jobs):
            if r:
                recs.append(r)
                print(f"  {r['family']:8} n={r['n']:2} {r.get('cls','?'):14} "
                      f"{'cyclic' if r.get('cyclic') else 'acyclic':8} mult={r.get('mult_values','?')} "
                      f"rep_cplx={r.get('rep_complexity','?')} {r['status']}")
    ok = [r for r in recs if r["status"] == "ok"]
    byc = collections.Counter(r["cls"] for r in ok)
    dist = sum(1 for r in ok if r["distinguishing"])
    # artifact: low-complexity repeats among CYCLE classes
    lowc = [r for r in ok if r["cls"] in ("CYCLE_CNV", "REPEAT_SHARED") and r["rep_complexity"] < 0.5]
    print(f"\n=== NON-LINEAR CENSUS ({len(ok)} graphs built, {len(recs)-len(ok)} failed) ===")
    for c, n in byc.most_common():
        print(f"  {c:14}: {n}")
    print(f"  copy-DISTINGUISHING structure (CYCLE_CNV or REARRANGEMENT): {dist}")
    print(f"  low-complexity-repeat (artifact-suspect, rep_complexity<0.5): {len(lowc)}")
    print(f"  => genuine structural distinguishers: {dist - sum(1 for r in lowc if r['distinguishing'])}")
    json.dump(dict(records=recs, by_class=dict(byc), distinguishing=dist,
                   low_complexity=len(lowc)),
              open(os.path.join(HERE, "family_def_nonlinear_census.json"), "w"), indent=2)
    with open(os.path.join(HERE, "family_def_nonlinear_census.tsv"), "w") as out:
        cols = ["family", "n", "status", "cls", "cyclic", "nodes", "edges", "max_mult",
                "mult_values", "rearrangement", "rep_len", "rep_complexity", "distinguishing"]
        out.write("\t".join(cols) + "\n")
        for r in recs:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"\n[+] wrote family_def_nonlinear_census.{{tsv,json}}")


if __name__ == "__main__":
    main()
