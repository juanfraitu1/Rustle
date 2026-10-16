#!/usr/bin/env python3
"""Graph-based O2 (prereg docs/PREREG_o3_graph_real_2026-10-15.md): assign a family's reads to its copies on the family's
variation graph, where every curated isoform's exons carry their PSVs (a representative lacks some of them).

    o2graph.py run WORLD [paths|nodes] [MARGIN]   WORLD/graph/o2.tsv; `paths` (default): reads aligned to every curated isoform
                                                  path of every catalogued copy; `nodes`: allele nodes of a free GraphAligner
                                                  alignment (needs detect.py's graph/<family>/s.gfa + aln.gaf)
    o2graph.py score OUT.json O2_NAME WORLD ...   O2 (graph/O2_NAME) against the truth, the current method beside it
    (run ... reps MARGIN: the ablation with each copy's representative path only; run ... union MARGIN: every copy's own
    isoforms plus every paralog isoform structure projected onto its locus (union_structures); run ... unione MARGIN: union,
    read bases past a path's end not counted (_paf_dist ends=True); run ... unionp MARGIN: union, every path's two ends
    extended by PATH_PAD bases of the copy's genome (padded); run ... locus MARGIN: the
    annotation-free baseline on seqsupport's genomic-locus distances; outputs graph/o2.<method>.m<margin>.tsv)

A read's informative nodes are the graph nodes it traverses (its best GraphAligner alignments, detect.best_records) that lie
on some copies' genomic paths but not all. It is assigned to the copy whose path misses the fewest of them when the next copy
misses >= MARGIN more; otherwise it abstains (`tied`, `thin`, `no_informative`). A read the graph does not explain (identity
< 1 - DELTA or cover < MIN_COVER) abstains as `unexplained`: it belongs to no catalogued copy (O3's business).
Scoring, per read whose source copy is catalogued: assigned (a copy is named) and correct (it is the source); per read of a
copy NOT in the catalog (absent, uncatalogued, novel-exon): assigned at all (wrong by construction); any assignment in a family
other than the read's own (`*_wrong_family`, counted in the family that made it) is wrong. The current method's call is
copy_assign's `assigned` row (catalog_copy_idx). Reads of a structure annotated only on a paralog (CROSS, truth isoform
`proj:...`) are also counted on their own (`*_proj`).
    o2graph.py seqsupport WORLD    WORLD/graph/locus.dist.json: every hosted read spliced-aligned to each catalogued copy's
                                   GENOMIC locus (+-1 kb); `score` then also counts the correct calls the genome itself
                                   separates (`*_correct_seq`: the read's own locus beats every other copy's by >= MARGIN),
                                   i.e. without the calls only the annotated isoform structure decides (review C1)
"""
import collections
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MARGIN, MIN_COVER, DELTA, END_SLACK, PATH_PAD = 2, 0.8, 0.00958, 5, 50


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


def seq_supported(gd, own, margin=MARGIN):
    """gd: {copy: edits + unaligned bases of the read on that copy's genomic locus}. Does the genome itself put the read on
    `own` -- its own locus explains it, every other copy's locus by >= margin worse (or not at all)?"""
    if own not in gd:
        return False
    others = [v for c, v in gd.items() if c != own]
    return not others or min(others) - gd[own] >= margin


def _paf_dist(paf, ends=False, slack=END_SLACK):
    """{read: {target: min(NM + unaligned read bases)}} from minimap2 -c output (target `copyN|...` -> copyN). ends: a read end
    left unaligned where the alignment reaches the path's end (within `slack` bases) is not counted -- bases past an annotated
    transcript end are annotation, not a difference between copies (register 1298 / union prereg U4: AMY1's 16-bp end)."""
    dist = collections.defaultdict(dict)
    for ln in paf.splitlines():
        f = ln.split("\t")
        nm = int(next(x for x in f[12:] if x.startswith("NM:i:"))[5:])
        qlen, qs, qe, tlen, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        clip_start, clip_end = qs, qlen - qe                  # the read's own start and end
        if ends:
            at_t_start, at_t_end = ts <= slack, te >= tlen - slack
            if f[4] == "+":
                clip_start, clip_end = (0 if at_t_start else clip_start), (0 if at_t_end else clip_end)
            else:                                             # minus: the read's start faces the path's end
                clip_start, clip_end = (0 if at_t_end else clip_start), (0 if at_t_start else clip_end)
        dd, c = nm + clip_start + clip_end, f[5].split("|")[0]
        dist[f[0]][c] = min(dd, dist[f[0]].get(c, 10 ** 9))
    return dist


def _read_fq(p):
    reads, name = {}, None
    with open(p) as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 0:
                name = ln[1:].strip().split()[0]
            elif i % 4 == 1:
                reads[name] = ln.strip()
    return reads


def run_seqsupport(w, pad=1000):
    """WORLD/graph/locus.dist.json: {family: {read: {copyN: distance}}} for every hosted read on each catalogued copy's genomic
    locus (minimap2 -c -x splice:hq -uf, all secondaries): what the sequence alone says, whatever isoforms are annotated."""
    import subprocess
    import pysam
    fa = pysam.FastaFile(f"{w}/genome.ref.fa")
    iso = json.load(open(f"{w}/isoforms.json"))
    reads = _read_fq(f"{w}/reads.E0.fq")
    hosted = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        hosted[r["family_id"]].add(r["read_name"])
    d = f"{w}/graph/o2loci"
    os.makedirs(d, exist_ok=True)
    out = {}
    for fam in sorted(iso):
        names = sorted(n for n in hosted.get(fam, ()) if n in reads)
        if not names:
            continue
        with open(f"{d}/{fam}.loci.fa", "w") as fh:
            for c in iso[fam]:
                fh.write(f">copy{c['idx']}\n{fa.fetch(c['contig'], max(0, c['start'] - pad), c['end'] + pad).upper()}\n")
        with open(f"{d}/{fam}.reads.fa", "w") as fh:
            for n in names:
                fh.write(f">{n}\n{reads[n]}\n")
        res = subprocess.run(["minimap2", "-c", "-x", "splice:hq", "-uf", "-N", "50", "-p", "0.3", "--secondary=yes", "-t", "2",
                              f"{d}/{fam}.loci.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True, check=True).stdout
        out[fam] = _paf_dist(res)
    json.dump(out, open(f"{w}/graph/locus.dist.json", "w"))
    print(f"[o2graph] {w}/graph/locus.dist.json: {sum(len(v) for v in out.values())} reads", flush=True)
    return 0


def _realworld():
    spec = importlib.util.spec_from_file_location("o3_graph_realworld", os.path.join(HERE, "realworld.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def padded(exons, pad):
    """A path's exons with both transcript ends extended by `pad` genomic bases (clamped at the contig start): a read running
    past an annotated end is then compared with each copy's own sequence there, not penalised (or excused) for the annotation."""
    ex = sorted(tuple(x) for x in exons)
    if not pad:
        return ex
    ex[0] = (max(0, ex[0][0] - pad), ex[0][1])
    ex[-1] = (ex[-1][0], ex[-1][1] + pad)
    return ex


def union_structures(cps, projector):
    """{copy idx: {name: exons}}: each catalogued copy's own curated isoforms plus every paralog isoform structure it can
    carry -- projector(copy, {name: (paralog copy, isoform exons)}) -> {name: exons on this copy or None} -- whose intron
    chain is new to it (realworld.cross_structures, uncapped). The family graph threads every structure through every copy,
    so only the sequence separates copies, not which copy the annotation lists a structure on."""
    R = _realworld()
    out = {}
    for c in cps:
        own = {t: [tuple(x) for x in ex] for t, ex in c["isoforms"].items()}
        queries = {f"proj:copy{p['idx']}|{t}": (p, ex) for p in cps if p is not c for t, ex in p["isoforms"].items()}
        proj = projector(c, queries) if queries else {}
        out[c["idx"]] = dict(own, **R.cross_structures(own, proj, max_proj=None))
    return out


def run_paths(w, margin=MARGIN, reps_only=False, out_name="o2.tsv", union=False, ends=False, pad=0):
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
    d = f"{w}/graph/o2paths" + (".reps" if reps_only else ".union" if union else "") + (".ends" if ends else "") + (".pad" if pad else "")
    os.makedirs(d, exist_ok=True)
    R = _realworld() if union else None

    def projector(c, queries):                       # union: paralog isoforms spliced-aligned onto the copy's locus
        lo = max(0, c["start"] - R.PROJ_PAD)
        locus = fa.fetch(c["contig"], lo, c["end"] + R.PROJ_PAD).upper()
        q = {}
        for name, (p, ex) in queries.items():
            seq = "".join(fa.fetch(p["contig"], a, b).upper() for a, b in ex)
            q[name] = (seq[::-1].translate(comp) if p["strand"] == "-" else seq, len(ex))
        return R.project_many(q, locus, lo, c["strand"])
    out, n_struct = [], collections.Counter()
    for fam in sorted(iso):
        names = sorted(n for n in hosted.get(fam, ()) if n in reads)
        if not names:
            continue
        paths = (union_structures(iso[fam], projector) if union else
                 {c["idx"]: {t: ex for t, ex in c["isoforms"].items() if not reps_only or t == c["rep"]} for c in iso[fam]})
        json.dump({str(k): {t: [list(x) for x in ex] for t, ex in v.items()} for k, v in paths.items()},
                  open(f"{d}/{fam}.paths.json", "w"))                  # the structures each copy's paths follow
        n_tx = 0
        with open(f"{d}/{fam}.tx.fa", "w") as fh:
            for c in iso[fam]:
                for t, ex in paths[c["idx"]].items():
                    seq = "".join(fa.fetch(c["contig"], a, b).upper() for a, b in padded(ex, pad))
                    fh.write(f">copy{c['idx']}|{t}\n{seq[::-1].translate(comp) if c['strand'] == '-' else seq}\n")
                    n_tx += 1
                    n_struct["projected" if t.startswith("proj:") else "own"] += 1
        with open(f"{d}/{fam}.reads.fa", "w") as fh:
            for n in names:
                fh.write(f">{n}\n{reads[n]}\n")
        res = subprocess.run(["minimap2", "-c", "-x", "map-hifi", "-N", str(max(200, 2 * n_tx)), "-p", "0.3", "--secondary=yes",
                              f"{d}/{fam}.tx.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True, check=True).stdout
        dist = _paf_dist(res, ends=ends)
        for n in names:
            c, reason = decide(dist.get(n, {}), len(reads[n]), margin)
            out.append(dict(read=n, family=fam, call=c[4:] if c else "-", reason=reason))
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    kind = (("representative" if reps_only else "union" if union else "isoform") + (" (path ends not counted)" if ends else "")
            + (f" (ends padded {pad} bp)" if pad else ""))
    print(f"[o2graph] {w} ({kind} paths {dict(n_struct)}, margin {margin}) -> {out_name}: {len(out)} reads, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return 0


def run_locus(w, margin=MARGIN, out_name=None):
    """Locus O2, the annotation-free baseline: decide() on the distances of every hosted read to each catalogued copy's
    genomic locus (graph/locus.dist.json from run_seqsupport; any canonical splicing allowed)."""
    lj = json.load(open(f"{w}/graph/locus.dist.json"))
    reads = _read_fq(f"{w}/reads.E0.fq")
    out = []
    for fam in sorted(lj):
        for n in sorted(lj[fam]):
            c, reason = decide(lj[fam][n], len(reads[n]), margin)
            out.append(dict(read=n, family=fam, call=c[4:] if c else "-", reason=reason))
    out_name = out_name or f"o2.locus.m{margin}.tsv"
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w} (genomic loci, margin {margin}) -> {out_name}: {len(out)} reads, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
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


def score_world(w, o2_name="o2.tsv", locus=None):
    """Per family: reads whose source copy is catalogued; assigned and correct, by the current method and the graph; every
    assignment made in a family that is not the read's own (wrong family). `locus` (run_seqsupport's {family: {read: {copyN:
    distance}}}, the read's own family used): also the correct calls the genome separates (`*_correct_seq`)."""
    idx_of = {}
    for r in csv.DictReader(open(f"{w}/P.fam.copies.tsv"), delimiter="\t"):
        idx_of[r["tid"]] = (r["family_id"], r["copy_idx"])
    truth, noncat, fam_of, proj = {}, {}, {}, set()
    for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t"):
        fam_of[r["read"]] = idx_of[r["copy"]][0] if r["copy"] in idx_of else r["family"]
        if r.get("isoform", "").startswith("proj:"):      # CROSS: a structure annotated only on a paralog of its copy
            proj.add(r["read"])
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
        if q in proj:
            t["proj_reads"] += 1
            for name, call in (("cur", cc), ("graph", gc)):
                t[f"{name}_assigned_proj"] += call is not None
                t[f"{name}_correct_proj"] += call == idx
        if locus is not None and seq_supported(locus.get(fam, {}).get(q, {}), f"copy{idx}"):
            t["seq_supported_reads"] += 1
            t["cur_correct_seq"] += cc == idx
            t["graph_correct_seq"] += gc == idx
        if cc is None:
            t["cur_unassigned"] += 1
            t["graph_correct_where_cur_unassigned"] += gc == idx
    for q, fam in noncat.items():
        t = per[fam]
        t["noncat_reads"] += 1
        t["cur_assigned_noncat"] += (q, fam) in cur
        t["graph_assigned_noncat"] += (q, fam) in gr
    for name, calls in (("cur", cur), ("graph", gr)):
        for (q, fam) in calls:
            if q in fam_of and fam_of[q] != fam:
                per[fam][f"{name}_wrong_family"] += 1
    return {f: dict(v) for f, v in sorted(per.items())}


def gain_kind(true_exons, runner_up_structures, tol=10):
    """A correct call's separation: `psv` when the runner-up copy carries the read's own structure (an intron chain within tol
    bp at every junction), so only sequence separated them; `structure` when it cannot carry it."""
    def chain(ex):
        ex = sorted(tuple(x) for x in ex)
        return [(a[1], b[0]) for a, b in zip(ex, ex[1:])]
    c = chain(true_exons)
    for ex in runner_up_structures:
        k = chain(ex)
        if len(k) == len(c) and all(abs(x - y) <= tol for i, j in zip(c, k) for x, y in zip(i, j)):
            return "psv"
    return "structure"


def run_gains(w, method="union", margin=MARGIN):
    """Split `method`'s correct calls where the current O2 abstains into gain_kind's `psv` / `structure` (re-aligning the
    run's saved reads to its saved paths for each copy's distance; the read's structure from isoforms.json or cross.json,
    projected onto the runner-up copy's locus as union O2 projects -- structures on different copies share no coordinates)."""
    import subprocess
    import pysam
    R = _realworld()
    fa = pysam.FastaFile(f"{w}/genome.ref.fa")
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    d = f"{w}/graph/o2paths." + {"union": "union", "unionp": "union.pad", "unione": "union.ends"}[method]
    iso = json.load(open(f"{w}/isoforms.json"))
    cross = json.load(open(f"{w}/cross.json"))["expressed"] if os.path.exists(f"{w}/cross.json") else {}
    idx_of = {r["tid"]: (r["family_id"], r["copy_idx"]) for r in csv.DictReader(open(f"{w}/P.fam.copies.tsv"), delimiter="\t")}
    truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    cur = {(r["read_name"], r["family_id"]) for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t")
           if r["status"] == "assigned" and r["catalog_copy_idx"].isdigit()}
    calls = {(r["read"], r["family"]): r["call"] for r in csv.DictReader(open(f"{w}/graph/o2.{method}.m{margin}.tsv"), delimiter="\t")}
    out = collections.Counter()
    for fam in sorted(iso):
        if not os.path.exists(f"{d}/{fam}.tx.fa"):
            continue
        paths = json.load(open(f"{d}/{fam}.paths.json"))
        res = subprocess.run(["minimap2", "-c", "-x", "map-hifi", "-N", "1000", "-p", "0.3", "--secondary=yes",
                              f"{d}/{fam}.tx.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True, check=True).stdout
        dist = _paf_dist(res, ends=method == "unione")
        cp = {str(c["idx"]): c for c in iso[fam]}
        todo = collections.defaultdict(dict)              # runner-up copy -> {read: (its structure's sequence, exon count)}
        for q, t in truth.items():
            if t["copy"] not in idx_of or idx_of[t["copy"]][0] != fam:
                continue
            i = idx_of[t["copy"]][1]
            if calls.get((q, fam)) != i or (q, fam) in cur:
                continue
            ex = cross[t["copy"]][t["isoform"]] if t["isoform"].startswith("proj:") else cp[i]["isoforms"][t["isoform"]]
            ranked = sorted((v, c) for c, v in dist.get(q, {}).items() if c != f"copy{i}")
            if not ranked:
                out["single copy aligned"] += 1
                continue
            seq = "".join(fa.fetch(cp[i]["contig"], a, b).upper() for a, b in sorted(tuple(x) for x in ex))
            todo[ranked[0][1][4:]][q] = (seq[::-1].translate(comp) if cp[i]["strand"] == "-" else seq, len(ex))
        for r, queries in todo.items():
            c = cp[r]
            lo = max(0, c["start"] - R.PROJ_PAD)
            proj = R.project_many(queries, fa.fetch(c["contig"], lo, c["end"] + R.PROJ_PAD).upper(), lo, c["strand"])
            for q, ex in proj.items():
                out[gain_kind(ex, list(paths[r].values())) if ex else "structure"] += 1
    return dict(out)


def main(argv):
    if argv[0] == "run":
        # run WORLD [paths|nodes] [MARGIN]: `nodes` = free GraphAligner alignment, kept for the record (development: it makes
        # artifactual allele calls through near-identical copies, GAGE 74 wrong of 289)
        method = argv[2] if len(argv) > 2 else "paths"
        margin = int(argv[3]) if len(argv) > 3 else MARGIN
        if method == "nodes":
            return run_world(argv[1])
        if method == "locus":
            return run_locus(argv[1], margin)
        return run_paths(argv[1], margin, reps_only=method == "reps", union=method in ("union", "unione", "unionp"),
                         ends=method == "unione", pad=PATH_PAD if method == "unionp" else 0,
                         out_name=f"o2.{method}.m{margin}.tsv")
    if argv[0] == "seqsupport":
        return run_seqsupport(argv[1])
    if argv[0] == "gains":                               # gains WORLD METHOD: the psv / structure split of METHOD's gains
        print(os.path.basename(argv[1]), argv[2], run_gains(argv[1], argv[2]))
        return 0
    # score OUT.json O2_NAME WORLD ...   (+ the sequence-supported counts where WORLD/graph/locus.dist.json exists)
    out_p, o2_name, worlds = argv[1], argv[2], argv[3:]
    res = {}
    for w in worlds:
        lp = f"{w}/graph/locus.dist.json"
        res[w] = score_world(w, o2_name, json.load(open(lp)) if os.path.exists(lp) else None)
    tot = collections.Counter()
    for w, fams in res.items():
        for f, v in fams.items():
            tot.update(v)
            print(f"{os.path.basename(w)} {f}: reads {v.get('reads', 0)} | current assigned {v.get('cur_assigned', 0)} correct "
                  f"{v.get('cur_correct', 0)} wrong-family {v.get('cur_wrong_family', 0)} | graph assigned "
                  f"{v.get('graph_assigned', 0)} correct {v.get('graph_correct', 0)} wrong-family {v.get('graph_wrong_family', 0)} | "
                  f"where current abstains: {v.get('cur_unassigned', 0)} reads, graph correct "
                  f"{v.get('graph_correct_where_cur_unassigned', 0)} | genome-separated correct: current "
                  f"{v.get('cur_correct_seq', '-')} graph {v.get('graph_correct_seq', '-')}")
    print("TOTAL", dict(tot))
    json.dump(dict(worlds=res, total=dict(tot)), open(out_p, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
