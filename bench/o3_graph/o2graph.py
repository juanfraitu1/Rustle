#!/usr/bin/env python3
"""Graph-based O2 (prereg docs/PREREG_o3_graph_real_2026-10-15.md): assign a family's reads to its copies on the family's
variation graph, where every curated isoform's exons carry their PSVs (a representative lacks some of them).

    o2graph.py run WORLD [paths|nodes] [MARGIN]   WORLD/graph/o2.tsv; `paths` (default): reads aligned to every curated isoform
                                                  path of every catalogued copy; `nodes`: allele nodes of a free GraphAligner
                                                  alignment (needs detect.py's graph/<family>/s.gfa + aln.gaf)
    o2graph.py score OUT.json O2_NAME WORLD ...   O2 (graph/O2_NAME) against the truth, the current method beside it
    (run ... reps MARGIN: the ablation with each copy's representative path only; outputs graph/o2.<method>.m<margin>.tsv)

A read's informative nodes are the graph nodes it traverses (its best GraphAligner alignments, detect.best_records) that lie
on some copies' genomic paths but not all. It is assigned to the copy whose path misses the fewest of them when the next copy
misses >= MARGIN more; otherwise it abstains (`tied`, `thin`, `no_informative`). A read the graph does not explain (identity
< 1 - DELTA or cover < MIN_COVER) abstains as `unexplained`: it belongs to no catalogued copy (O3's business).
Scoring, per read whose source copy is catalogued: assigned (a copy is named) and correct (it is the source); per read of a
copy NOT in the catalog (absent, uncatalogued, novel-exon): assigned at all (wrong by construction). The current method's call
is copy_assign's `assigned` row (catalog_copy_idx).
"""
import collections
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MARGIN, MIN_COVER, DELTA = 2, 0.8, 0.00958


def _detect():
    spec = importlib.util.spec_from_file_location("o3_graph_detect", os.path.join(HERE, "detect.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def assign(nodes, paths, margin=MARGIN):
    """-> (copy or None, {copy: informative nodes off its path})."""
    inf = [n for n in set(nodes) if 0 < sum(n in p for p in paths.values()) < len(paths)]
    miss = {c: sum(1 for n in inf if n not in p) for c, p in paths.items()}
    if not inf or len(paths) < 2:
        return None, miss
    order = sorted(miss.items(), key=lambda kv: (kv[1], kv[0]))
    return (order[0][0] if order[1][1] - order[0][1] >= margin else None), miss


def decide(dist, qlen, margin=MARGIN, max_div=DELTA):
    """dist: {copy: edits + unaligned bases of the read against the copy's best-fitting isoform path}. -> (copy or None, reason):
    `assigned` when the best copy beats the next by >= margin, else `tied` / `thin`; `unexplained` when even the best copy is
    more than max_div away (a read of no catalogued copy)."""
    if not dist:
        return None, "unexplained"
    srt = sorted(dist.items(), key=lambda kv: (kv[1], kv[0]))
    if srt[0][1] > max_div * qlen:
        return None, "unexplained"
    if len(srt) == 1:
        return srt[0][0], "assigned"
    gap = srt[1][1] - srt[0][1]
    if gap == 0:
        return None, "tied"
    return (srt[0][0], "assigned") if gap >= margin else (None, "thin")


def run_paths(w, margin=MARGIN, reps_only=False, out_name="o2.tsv"):
    """O2 on the family graph's ISOFORM PATHS: each hosted read aligned (minimap2 -x map-hifi, all secondaries) to every
    curated isoform of every catalogued copy as its own sequence -- the paths `vg rna` threads, without free graph alignment."""
    import subprocess
    import pysam
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    fa = pysam.FastaFile(f"{w}/genome.ref.fa")
    iso = json.load(open(f"{w}/isoforms.json"))
    reads = {}
    name = None
    with open(f"{w}/reads.E0.fq") as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 0:
                name = ln[1:].strip().split()[0]
            elif i % 4 == 1:
                reads[name] = ln.strip()
    hosted = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        hosted[r["family_id"]].add(r["read_name"])
    d = f"{w}/graph/o2paths"
    os.makedirs(d, exist_ok=True)
    out = []
    for fam in sorted(iso):
        names = sorted(n for n in hosted.get(fam, ()) if n in reads)
        if not names:
            continue
        with open(f"{d}/{fam}.tx.fa", "w") as fh:
            for c in iso[fam]:
                for t, ex in c["isoforms"].items():
                    if reps_only and t != c["rep"]:         # the ablation: each copy's representative only
                        continue
                    seq = "".join(fa.fetch(c["contig"], a, b).upper() for a, b in ex)
                    fh.write(f">copy{c['idx']}|{t}\n{seq[::-1].translate(comp) if c['strand'] == '-' else seq}\n")
        with open(f"{d}/{fam}.reads.fa", "w") as fh:
            for n in names:
                fh.write(f">{n}\n{reads[n]}\n")
        res = subprocess.run(["minimap2", "-c", "-x", "map-hifi", "-N", "200", "-p", "0.3", "--secondary=yes",
                              f"{d}/{fam}.tx.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True, check=True).stdout
        dist = collections.defaultdict(dict)
        for ln in res.splitlines():
            f = ln.split("\t")
            nm = int(next(x for x in f[12:] if x.startswith("NM:i:"))[5:])
            dd, c = nm + int(f[1]) - (int(f[3]) - int(f[2])), f[5].split("|")[0]
            dist[f[0]][c] = min(dd, dist[f[0]].get(c, 10 ** 9))
        for n in names:
            c, reason = decide(dist.get(n, {}), len(reads[n]), margin)
            out.append(dict(read=n, family=fam, call=c[4:] if c else "-", reason=reason))
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w} ({'representative' if reps_only else 'isoform'} paths, margin {margin}) -> {out_name}: {len(out)} reads, " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())),
          flush=True)
    return 0


def run_world(w):
    D = _detect()
    out = []
    for fam in sorted(os.listdir(f"{w}/graph")):
        d = f"{w}/graph/{fam}"
        if not os.path.isfile(f"{d}/aln.gaf"):
            continue
        paths = {}
        for ln in open(f"{d}/s.gfa"):
            f = ln.rstrip("\n").split("\t")
            if f[0] == "P" and re.fullmatch(r"copy\d+", f[1]):
                paths[f[1]] = {x[:-1] for x in f[2].split(",")}
        recs = collections.defaultdict(list)
        for ln in open(f"{d}/aln.gaf"):
            f = ln.rstrip("\n").split("\t")
            tags = {x[:5]: x[5:] for x in f[12:]}
            recs[f[0]].append(dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), m=int(f[9]), b=int(f[10]), path=f[5],
                                   score=float(tags.get("AS:f:", 0)), cigar=tags.get("cg:Z:", "")))
        for q, rs in sorted(recs.items()):
            best = D.best_records(rs)
            ident = sum(r["m"] for r in best) / max(1, sum(r["b"] for r in best))
            iv = sorted((r["qs"], r["qe"]) for r in best)
            cov, end = 0, -1
            for s, e in iv:
                if e > end:
                    cov += e - max(s, end)
                    end = e
            cov /= best[0]["qlen"]
            if ident < 1 - DELTA or cov < MIN_COVER:
                out.append(dict(read=q, family=fam, call="-", reason="unexplained"))
                continue
            nodes = [n for r in best for n in re.findall(r"[<>]([^<>]+)", r["path"])]
            c, miss = assign(nodes, paths)
            reason = "assigned" if c else ("no_informative" if not any(miss.values()) else
                                           "tied" if sorted(miss.values())[:2][0] == sorted(miss.values())[:2][-1] else "thin")
            out.append(dict(read=q, family=fam, call=c[4:] if c else "-", reason=reason))
    with open(f"{w}/graph/o2.tsv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w}: {len(out)} reads, " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return 0


def score_world(w, o2_name="o2.tsv"):
    """Per family: reads whose source copy is catalogued; assigned and correct, by the current method and the graph."""
    idx_of = {}
    for r in csv.DictReader(open(f"{w}/P.fam.copies.tsv"), delimiter="\t"):
        idx_of[r["tid"]] = (r["family_id"], r["copy_idx"])
    truth, noncat = {}, {}
    for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t"):
        if r["copy"] in idx_of:
            truth[r["read"]] = idx_of[r["copy"]]
        else:                                             # absent / uncatalogued / novel-exon copy: any assignment is wrong
            noncat[r["read"]] = r["family"]
    cur = {}
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        if r["status"] == "assigned" and r["catalog_copy_idx"].isdigit():
            cur[(r["read_name"], r["family_id"])] = r["catalog_copy_idx"]
    gr = {(r["read"], r["family"]): r["call"] for r in csv.DictReader(open(f"{w}/graph/{o2_name}"), delimiter="\t")
          if r["call"] != "-"}
    per = collections.defaultdict(lambda: collections.Counter())
    for q, (fam, idx) in truth.items():
        t = per[fam]
        t["reads"] += 1
        cc, gc = cur.get((q, fam)), gr.get((q, fam))
        t["cur_assigned"] += cc is not None
        t["cur_correct"] += cc == idx
        t["graph_assigned"] += gc is not None
        t["graph_correct"] += gc == idx
        if cc is None:
            t["cur_unassigned"] += 1
            t["graph_correct_where_cur_unassigned"] += gc == idx
    for q, fam in noncat.items():
        t = per[fam]
        t["noncat_reads"] += 1
        t["cur_assigned_noncat"] += (q, fam) in cur
        t["graph_assigned_noncat"] += (q, fam) in gr
    return {f: dict(v) for f, v in sorted(per.items())}


def main(argv):
    if argv[0] == "run":
        # run WORLD [paths|nodes] [MARGIN]: `nodes` = free GraphAligner alignment, kept for the record (development: it makes
        # artifactual allele calls through near-identical copies, GAGE 74 wrong of 289)
        method = argv[2] if len(argv) > 2 else "paths"
        margin = int(argv[3]) if len(argv) > 3 else MARGIN
        if method == "nodes":
            return run_world(argv[1])
        return run_paths(argv[1], margin, reps_only=method == "reps", out_name=f"o2.{method}.m{margin}.tsv")
    # score OUT.json O2_NAME WORLD ...
    out_p, o2_name, worlds = argv[1], argv[2], argv[3:]
    res = {w: score_world(w, o2_name) for w in worlds}
    tot = collections.Counter()
    for w, fams in res.items():
        for f, v in fams.items():
            tot.update(v)
            print(f"{os.path.basename(w)} {f}: reads {v['reads']} | current assigned {v['cur_assigned']} correct {v['cur_correct']} | "
                  f"graph assigned {v['graph_assigned']} correct {v['graph_correct']} | where current abstains: "
                  f"{v.get('cur_unassigned', 0)} reads, graph correct {v.get('graph_correct_where_cur_unassigned', 0)}")
    print("TOTAL", dict(tot))
    json.dump(dict(worlds=res, total=dict(tot)), open(out_p, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
