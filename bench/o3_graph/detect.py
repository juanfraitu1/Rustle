#!/usr/bin/env python3
"""The graph-informed O3 detector on a world (prereg docs/PREREG_o3_graph_sim_2026-10-15.md): the current method's per-host PSV
clustering (bench/o3_unexplained/psvclust.py, = `copy_assign --psv-clusters`) plus three steps, as separate arms.

    detect.py WORLD             reads WORLD/cur/A.assignments.tsv (the hosts), writes WORLD/graph/calls.tsv (+ graph/<family>/)

Per host class (copies with one representative sequence; its reads = O2's hosts, as the current method):
  the reads and the copies are aligned to the host locus, clustered (cluster_host), and each cluster's status follows the
  prototype's rule against a set of copy sequences:
    R0    each copy's REPRESENTATIVE (the catalog): the current method, reproduced
    S1    step 1: each copy's EVERY curated isoform (the union of their bases on the host) -- the isoform paths of the family's
          variation graph (`vg rna`)
  A PSV call = an `unexplained` cluster with est_div > DELTA (register 1294's `copy` label).
Step 2 (S12 = S1 + branch calls): the family's variation graph (vg msga of its copies' loci + every curated isoform as a path,
  `vg rna`); GraphAligner maps the family's hosted reads; a read is a BRANCH read when >= MIN_NOVEL of its bases lie on no
  graph path (large insertions, dense error windows where the aligner forced novel sequence onto the graph, gaps between
  aligned pieces, unaligned ends; an unaligned read counts whole). Branch reads
  sharing a novel segment (all-vs-all identity >= 0.97 over >= half the shorter segment) are linked; a group of >= MIN_BRANCH
  reads is a branch call.
Step 3 (S123 = S12 kept by the genome; S3 = R0 kept by the genome): a call is kept only when the genome does not explain its
  reads: PSV call -> the median `de` of their primary alignments > DELTA; branch call -> their median unaligned bases on the
  genome (insertions + soft clips >= MIN_NOVEL_PART) >= MIN_NOVEL, or the median `de` > DELTA.
"""
import collections
import csv
import json
import os
import re
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "o3_unexplained"))
import psvclust as PC  # noqa: E402

DELTA = 0.00958
MIN_NOVEL, MIN_NOVEL_PART, MIN_BRANCH, LOCUS_PAD = 50, 30, 5, 1000
DENSE_W, DENSE_FRAC = 50, 0.3
ARMS = ("R0", "S1", "S12", "S123", "S3")


def explains(cons, sites, ref, cc):
    """The prototype's (and `host_rows`') rule: the copy carries every consensus alt where it has a base, has a base at every
    alt, and no other base where the consensus is the reference."""
    alts = [p for p, a in cons.items() if a == "alt"]
    return (all(cc.get(p, sites[p]) == sites[p] for p in alts)
            and not any(a == "ref" and p in cc and cc[p] != ref[p] for p, a in cons.items())
            and all(p in cc for p in alts))


def union_map(maps):
    out = {}
    for m in maps:
        for p, b in m.items():
            out.setdefault(p, b)
    return out


def dense_intervals(rec, w=DENSE_W, frac=DENSE_FRAC):
    """Query intervals where the aligner FORCED sequence the graph lacks: >= frac * w errors (mismatched or inserted bases,
    deletion events) within w query bases, trimmed to the first and last error. Sequencing errors (~0.1%) and paralog PSVs
    (a few %) never reach it; a novel exon aligned base by base does (~50%)."""
    import bisect
    errs, q = [], rec["qs"]
    for k, op in re.findall(r"(\d+)([=XIDM])", rec["cigar"]):
        k = int(k)
        if op in "=M":
            q += k
        elif op in "XI":
            errs.extend(range(q, q + k))
            q += k
        elif op == "D":
            errs.append(q)
    errs.sort()
    out = []
    for p in errs:
        hi = bisect.bisect_left(errs, p + w)
        if hi - bisect.bisect_left(errs, p) >= frac * w:
            s_, e_ = p, errs[hi - 1] + 1
            if out and s_ <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e_))
            else:
                out.append((s_, e_))
    return out


def novel_intervals(recs, min_len=MIN_NOVEL_PART):
    """Query intervals of a read on no graph path: insertions >= min_len inside its alignments, dense error windows (sequence
    the aligner forced onto the graph), gaps >= min_len between alignments, and unaligned read ends >= min_len (the aligner
    often stops at a novel exon near a short last exon); overlapping intervals merged."""
    raw = []
    for r in recs:
        q = r["qs"]
        for k, op in re.findall(r"(\d+)([=XIDM])", r["cigar"]):
            k = int(k)
            if op in "=XM":
                q += k
            elif op == "I":
                if k >= min_len:
                    raw.append((q, q + k))
                q += k
        raw += [x for x in dense_intervals(r) if x[1] - x[0] >= min_len]
    iv = sorted((r["qs"], r["qe"]) for r in recs)
    for (_, e1), (s2, _) in zip(iv, iv[1:]):
        if s2 - e1 >= min_len:
            raw.append((e1, s2))
    end = max(e for _, e in iv)
    if iv[0][0] >= min_len:
        raw.append((0, iv[0][0]))
    if recs[0]["qlen"] - end >= min_len:
        raw.append((end, recs[0]["qlen"]))
    out = []
    for s_, e_ in sorted(raw):
        if out and s_ <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e_))
        else:
            out.append((s_, e_))
    return out


def novel_bp(recs, min_len=MIN_NOVEL_PART):
    return sum(e - s for s, e in novel_intervals(recs, min_len))


def link_segments(names, hits, min_ident=0.97):
    """Reads linked when any of their novel segments align (identity >= min_ident over >= half the shorter segment): the reads
    sharing one missing branch. hits: (query segment, target segment, identity, block, query length, target length);
    segment ids are `read#k`. -> groups of read names."""
    par = {n: n for n in names}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for qa, tb, ident, blk, ql, tl in hits:
        a, b = qa.rsplit("#", 1)[0], tb.rsplit("#", 1)[0]
        if a != b and a in par and b in par and ident >= min_ident and blk >= 0.5 * min(ql, tl):
            par[find(a)] = find(b)
    groups = collections.defaultdict(list)
    for n in names:
        groups[find(n)].append(n)
    return list(groups.values())


def cigar_novel(cigartuples, min_len=MIN_NOVEL_PART):
    """Bases of a read the genome alignment leaves out: insertions and soft clips >= min_len (pysam op codes 1, 4)."""
    return sum(n for op, n in cigartuples if op in (1, 4) and n >= min_len)


def genome_keeps(kind, median_de, median_novel):
    if kind == "branch":
        return median_novel >= MIN_NOVEL or median_de > DELTA
    return median_de > DELTA


def status_of(members, cons, sites, ref, maps, cov, vec):
    alts = [p for p, a in cons.items() if a == "alt"]
    by = [n for n, cc in maps.items() if explains(cons, sites, ref, cc)]
    unexpl = [p for p in alts if not any(cc.get(p) == sites[p] for cc in maps.values())]
    covered = sum(1 for p in range(len(ref)) if any(p in cov[m] for m in members[:50]))
    link = PC.linkage(unexpl, vec) if len(unexpl) >= PC.MIN_ALT else None
    if by:
        st = "explained"
    elif len(unexpl) < PC.MIN_ALT or link is None:
        st = "unresolved"
    else:
        st = "unexplained" if link >= PC.MIN_LINK else "unlinked"
    return st, (len(unexpl) / covered if covered else 0.0)


def best_records(recs):
    out = []
    for r in sorted(recs, key=lambda r: -r["score"]):
        if all(min(r["qe"], o["qe"]) - max(r["qs"], o["qs"]) <= 0.5 * (r["qe"] - r["qs"]) for o in out):
            out.append(r)
    return out


def sh(cmd, binary_out=None):
    if binary_out:
        with open(binary_out, "wb") as fh:
            r = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, executable="/bin/bash")
    else:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
    if r.returncode:
        err = r.stderr if isinstance(r.stderr, str) else r.stderr.decode(errors="replace")
        raise SystemExit(f"failed: {cmd}\n{err[-600:]}")
    return r


def family_graph(d, genome, copies):
    """vg graph of the family's loci (genomic orientation; every copy of a world family shares its strand) with every curated
    isoform as a path -> GFA path."""
    os.makedirs(d, exist_ok=True)
    with open(f"{d}/loci.fa", "w") as fa:
        for c in copies:
            lo = max(0, c["start"] - LOCUS_PAD)
            hi = min(len(genome[c["contig"]]), c["end"] + LOCUS_PAD)
            c["lo"] = lo
            fa.write(f">copy{c['idx']}\n{genome[c['contig']][lo:hi]}\n")
            with open(f"{d}/copy{c['idx']}.gtf", "w") as g:
                for iso, ex in c["isoforms"].items():
                    for s, e in ex:
                        g.write(f"copy{c['idx']}\tsim\texon\t{s - lo + 1}\t{e - lo}\t.\t{c['strand']}\t.\t"
                                f"transcript_id \"copy{c['idx']}|{iso}\"; gene_id \"copy{c['idx']}\";\n")
    sh(f"vg msga -f {d}/loci.fa -t 2", binary_out=f"{d}/g.vg")
    gtfs = " ".join(f"-n {d}/copy{c['idx']}.gtf" for c in copies)
    sh(f"vg rna {gtfs} -r -t 2 {d}/g.vg", binary_out=f"{d}/s.vg")
    sh(f"vg view {d}/s.vg > {d}/s.gfa")
    return f"{d}/s.gfa"


def branch_calls(d, gfa, names, reads):
    """Branch reads of a family (>= MIN_NOVEL bases on no graph path) and their single-linkage clusters (>= MIN_BRANCH)."""
    with open(f"{d}/reads.fa", "w") as fh:
        for n in names:
            fh.write(f">{n}\n{reads[n]}\n")
    sh(f"GraphAligner -g {gfa} -f {d}/reads.fa -a {d}/aln.gaf -x vg -t 2")
    recs = collections.defaultdict(list)
    for ln in open(f"{d}/aln.gaf"):
        f = ln.rstrip("\n").split("\t")
        tags = {x[:5]: x[5:] for x in f[12:]}
        recs[f[0]].append(dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), score=float(tags.get("AS:f:", 0)),
                               cigar=tags.get("cg:Z:", "")))
    segs = {}
    for q in names:
        iv = novel_intervals(best_records(recs[q])) if q in recs else [(0, len(reads[q]))]   # unaligned: all novel
        if sum(e - s for s, e in iv) >= MIN_NOVEL:
            segs[q] = [reads[q][s:e] for s, e in iv]
    if len(segs) < MIN_BRANCH:
        return [], segs
    with open(f"{d}/segs.fa", "w") as fh:
        for q, ss in segs.items():
            for k, x in enumerate(ss):
                fh.write(f">{q}#{k}\n{x}\n")
    hits = []
    # short segments: the asm presets need an alignment score >= 200
    for ln in sh(f"minimap2 -c -X -k 13 -w 5 -s 50 -m 20 -t 2 {d}/segs.fa {d}/segs.fa").stdout.splitlines():
        f = ln.split("\t")
        hits.append((f[0], f[5], int(f[9]) / int(f[10]), int(f[10]), int(f[1]), int(f[6])))
    groups = link_segments(sorted(segs), hits)
    return [sorted(m) for m in groups if len(m) >= MIN_BRANCH], segs


def genome_stats(bam_p, want):
    import pysam
    out = {}
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name in want and not rec.flag & 0x904 and rec.cigartuples:
                out[rec.query_name] = (rec.get_tag("de") if rec.has_tag("de") else 0.0, cigar_novel(rec.cigartuples))
    return out


def main(argv):
    w = argv[0]
    genome = PC.read_fa(f"{w}/genome.ref.fa")
    reads = PC.read_fq(f"{w}/reads.E0.fq")
    rows = list(csv.DictReader(open(f"{w}/P.fam.copies.tsv"), delimiter="\t"))
    reps = list(PC.read_fa(f"{w}/P.fam.copies.fa").values())
    isoj = json.load(open(f"{w}/isoforms.json"))
    fam = collections.defaultdict(list)
    for r, seq in zip(rows, reps):
        f, i = r["family_id"], int(r["copy_idx"])
        c = next(x for x in isoj[f] if x["idx"] == i)
        fam[f].append(dict(idx=i, contig=r["chrom"], rs=int(r["start"]), re=int(r["end"]), seq=seq, start=c["start"],
                           end=c["end"], strand=c["strand"], isoforms={k: [tuple(x) for x in v] for k, v in c["isoforms"].items()}))
    by_host = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        if r["catalog_copy_idx"].isdigit():
            by_host[(r["family_id"], int(r["catalog_copy_idx"]))].append(r["read_name"])
    out_dir = f"{w}/graph"
    os.makedirs(out_dir, exist_ok=True)
    tmp = f"{out_dir}/tmp"
    calls = []                                     # (arm, family, kind, host, cluster, members, est_div)
    for f in sorted(fam):
        cps = fam[f]
        classes = collections.defaultdict(list)
        for c in cps:
            classes[c["seq"]].append(c)
        for cls in classes.values():
            rep = min(cls, key=lambda c: c["idx"])
            names = sorted({n for c in cls for n in by_host.get((f, c["idx"]), [])})
            if len(names) < PC.MIN_SUPPORT:
                continue
            g = genome[rep["contig"]]
            lo, hi = max(0, rep["rs"] - PC.PAD), min(len(g), rep["re"] + PC.PAD)
            ref = g[lo:hi].upper()
            cov = PC.align("host", ref, {n: reads[n] for n in names if n in reads}, tmp)
            m_rep = PC.align("host", ref, {f"copy{c['idx']}": c["seq"] for c in cps}, tmp)
            iso_seqs = {f"copy{c['idx']}|{k}": "".join(g2[a:b] for a, b in ex) for c in cps
                        for g2 in [genome[c["contig"]]] for k, ex in c["isoforms"].items()}
            m_iso_each = PC.align("host", ref, iso_seqs, tmp)
            m_iso = {f"copy{c['idx']}": union_map([m for k, m in m_iso_each.items() if k.startswith(f"copy{c['idx']}|")])
                     for c in cps}
            sites, clusters, vec = PC.cluster_host(ref, cov)
            host = ",".join(str(c["idx"]) for c in sorted(cls, key=lambda c: c["idx"]))
            for k, (members, cons) in enumerate(sorted(clusters, key=lambda x: -len(x[0]))):
                for arm, maps in (("R0", m_rep), ("S1", m_iso)):
                    st, div = status_of(members, cons, sites, ref, maps, cov, vec)
                    if st == "unexplained" and div > DELTA:
                        calls.append([arm, f, "psv", host, k, sorted(members), div])
        names = sorted({n for c in cps for n in by_host.get((f, c["idx"]), []) if n in reads})
        if names:
            gfa = family_graph(f"{out_dir}/{f}", genome, cps)
            br, _ = branch_calls(f"{out_dir}/{f}", gfa, names, reads)
            for k, m in enumerate(br):
                calls.append(["S2", f, "branch", "-", k, m, 0.0])
    want = {m for c in calls for m in c[5]}
    gst = genome_stats(f"{w}/E0.bam", want)
    rows_out = []
    for arm, f, kind, host, k, m, div in calls:
        des = [gst[x][0] for x in m if x in gst]
        nov = [gst[x][1] for x in m if x in gst]
        mde = statistics.median(des) if des else 0.0
        mnov = statistics.median(nov) if nov else 0
        keep = genome_keeps(kind, mde, mnov)
        arms = {"R0": ["R0"] + (["S3"] if keep else []),
                "S1": ["S1", "S12"] + (["S123"] if keep else []),
                "S2": ["S12"] + (["S123"] if keep else [])}[arm]
        for a in arms:
            rows_out.append(dict(arm=a, family=f, kind=kind, host=host, cluster=k, n_reads=len(m), est_div=round(div, 5),
                                 genome_median_de=round(mde, 5), genome_median_novel=mnov, members=",".join(m)))
    with open(f"{out_dir}/calls.tsv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["arm", "family", "kind", "host", "cluster", "n_reads", "est_div", "genome_median_de",
                                            "genome_median_novel", "members"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows_out)
    for x in (".t.fa", ".q.fa"):
        if os.path.exists(tmp + x):
            os.remove(tmp + x)
    print(f"[o3_graph] {w}: " + ", ".join(f"{a} {sum(1 for r in rows_out if r['arm'] == a)}" for a in ARMS) + " calls", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
