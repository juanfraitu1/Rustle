#!/usr/bin/env python3
"""The graph-informed O3 detector on a world (prereg docs/PREREG_o3_graph_sim_2026-10-15.md): the current method's per-host PSV
clustering (bench/o3_unexplained/psvclust.py, = `copy_assign --psv-clusters`) plus three steps, as separate arms.

    detect.py WORLD             reads WORLD/cur/A.assignments.tsv (the hosts), writes WORLD/graph/calls.tsv (+ graph/<family>/)

Per host class (copies with one representative sequence; its reads = O2's hosts, as the current method):
  the reads and the copies are aligned to the host locus, clustered (cluster_host), and each cluster's status follows the
  prototype's rule against a set of copy sequences:
    R0    each copy's REPRESENTATIVE (the catalog): the current method, reproduced
    S1    step 1: each copy's EVERY curated isoform (the union of their bases on the host), aligned to the host window with
          minimap2 -- the same isoforms the family's vg graph carries as paths
  A PSV call = an `unexplained` cluster with est_div > DELTA (register 1294's `copy` label).
Step 2 (S12 = S1 + branch calls): the family's variation graph (vg msga of its copies' loci + every curated isoform as a path,
  `vg rna`; the graph holds the whole loci, introns included); GraphAligner maps the family's hosted reads; a read is a BRANCH
  read when >= MIN_NOVEL of its bases do not align to the loci graph (large insertions, dense error windows where the aligner forced novel sequence onto the graph, gaps between
  aligned pieces, unaligned ends; an unaligned read counts whole). Branch reads
  sharing a novel segment (all-vs-all identity >= 0.97 over >= half the shorter segment) are linked; a group of >= MIN_BRANCH
  reads is a branch call.
S4 (graph O2 into O3): the S1 clustering run on the reads path-O2 (o2graph.py) cannot explain with any catalogued copy's
  isoform paths (`unexplained`), plus the branch calls (S4n), kept by the genome (S4): reads a catalogued copy explains
  can no longer dilute an absent copy's cluster.
G (union-O2-informed O3; prereg docs/PREREG_o3_graph_union_2026-10-15.md): the same on the reads UNION O2 cannot explain
  (every family structure threaded through every copy, o2graph.py run WORLD union 2), each cluster's status with the MOSAIC
  rule (status_of mosaic=True: no single copy within DELTA), plus the branch calls (Gn), kept by the genome (G). S123m = S123
  with the mosaic rule, without the O2 restriction.
H (prereg docs/PREREG_o3_graph_union2_2026-10-15.md): as G, but on the reads union O2 with padded path ends (o2graph.py run
  WORLD unionp 2) cannot explain, POOLED over every host of the family and aligned to each host (family_pool): an absent
  copy copy_assign spreads over several hosts clusters once (Hn without the genome check).
Step 3 (S123 = S12 kept by the genome; S3 = R0 kept by the genome; S23 = S3 + kept branch calls, i.e. without step 1):
  a call is kept only when the genome does not explain its
  reads: PSV call -> the median `de` of their primary alignments > DELTA; branch call -> their median unaligned bases on the
  genome (clips and dense error windows >= MIN_NOVEL_PART, minus what the read's supplementary alignments cover) >= MIN_NOVEL,
  or the median `de` > DELTA. A member without a primary alignment counts as wholly unexplained.
"""
import collections
import csv
import hashlib
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
MIN_OVERLAP = 0.5
ARMS = ("R0", "S1", "S12", "S123", "S3", "S23", "S4n", "S4", "S123m", "Gn", "G", "Hn", "H", "Kn", "K", "In", "I")


class _Contig:
    def __init__(self, fa, name):
        self.fa, self.name, self.n = fa, name, fa.get_reference_length(name)

    def __len__(self):
        return self.n

    def __getitem__(self, sl):
        lo, hi = sl.start or 0, self.n if sl.stop is None else min(self.n, sl.stop)
        return self.fa.fetch(self.name, lo, hi).upper()


class Genome:
    """The world's genome read lazily by slice (a whole CHM13 does not fit a dict of strings)."""

    def __init__(self, path):
        import pysam
        self.fa, self._c = pysam.FastaFile(path), {}

    def __getitem__(self, name):
        if name not in self._c:
            self._c[name] = _Contig(self.fa, name)
        return self._c[name]


def orient(seq, exons, strand, lo):
    """A locus and an isoform's exons (genomic, 0-based half-open) -> (the locus in its transcript's orientation, the exons
    relative to it): real families put copies on both strands, and vg msga aligns sequences as given."""
    rel = sorted((s_ - lo, e_ - lo) for s_, e_ in exons)
    if strand != "-":
        return seq, rel
    n = len(seq)
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq[::-1].translate(comp), sorted((n - e_, n - s_) for s_, e_ in rel)


def restrict(names, family, keep_set):
    """S4: only the reads no catalogued copy's isoform paths explain IN THIS FAMILY (path-O2 `unexplained`, keyed by (read,
    family): a read unexplained only in a sister family it was hosted in is explained here); none without the table."""
    return [n for n in names if keep_set is not None and (n, family) in keep_set]


def dedupe(calls):
    """One call per (arm, family, kind, member set): H aligns the same pooled reads to every host, and a host class that
    clusters them identically would flag the absent copy once per host (round 3 review I5)."""
    seen, out = set(), []
    for c in calls:
        key = (c[0], c[1], c[2], frozenset(c[5]))
        if key not in seen:
            seen.add(key)
            out.append(c)
    return out


def family_pool(by_host, family, keep_set):
    """H: every read hosted anywhere in the family that union O2 (padded path ends) cannot explain there -- pooled across
    hosts, so an absent copy whose reads copy_assign spreads over several hosts clusters once (PSG, union prereg G3)."""
    if keep_set is None:
        return []
    return sorted({n for (f, _), names in by_host.items() if f == family for n in names if (n, family) in keep_set})


def arms_of(arm, keep):
    """The arms a call belongs to: R0 (+ S3, S23 when the genome keeps it); S1 (+ S12, S123); S2 branch calls (S12, S4n, Gn; +
    S123, S23, S4, S123m, G); S4p, a PSV call among path-O2-unexplained reads (S4n; + S4); S1m, a PSV or MOSAIC call of S1's
    clusters (+ S123m); Gp, a PSV or mosaic call among the reads UNION O2 cannot explain (Gn; + G). G = union-O2-informed O3:
    clustering only reads no copy explains on any family structure, the mosaic rule, branch calls, the genome check."""
    return {"R0": ["R0"] + (["S3", "S23"] if keep else []),
            "S1": ["S1", "S12"] + (["S123"] if keep else []),
            "S2": ["S12"] + (["S123", "S23"] if keep else []) + ["S4n"] + (["S4"] if keep else []) + (["S123m"] if keep else [])
                  + ["Gn"] + (["G"] if keep else []) + ["Hn"] + (["H"] if keep else []) + ["Kn"] + (["K"] if keep else [])
                  + ["In"] + (["I"] if keep else []),
            "S4p": ["S4n"] + (["S4"] if keep else []),
            "S1m": ["S123m"] if keep else [],
            "Gp": ["Gn"] + (["G"] if keep else []),
            "Hp": ["Hn"] + (["H"] if keep else []),
            "Kp": ["Kn"] + (["K"] if keep else []),
            "Ip": ["In"] + (["I"] if keep else [])}[arm]


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
    """Bases of a read the genome alignment does not explain: soft clips >= min_len and dense error windows >= min_len (the
    aligner forcing novel sequence onto the genome, often as many small insertions; needs an =/X CIGAR, `--eqx`). pysam op
    codes: 0 M, 1 I, 2 D, 3 N, 4 S, 7 =, 8 X."""
    code = {0: "M", 1: "I", 2: "D", 7: "=", 8: "X"}
    clips = sum(n for op, n in cigartuples if op == 4 and n >= min_len)
    lead = cigartuples[0][1] if cigartuples and cigartuples[0][0] == 4 else 0
    cig = "".join(f"{n}{code[op]}" for op, n in cigartuples if op in code)
    dense = sum(e - s for s, e in dense_intervals(dict(qs=lead, cigar=cig)) if e - s >= min_len)
    return clips + dense


def genome_keeps(kind, median_de, median_novel):
    if kind == "branch":
        return median_novel >= MIN_NOVEL or median_de > DELTA
    return median_de > DELTA


def member_filter(kind, members, des, novs):
    """The genome check per READ (2026-10-16, `--member-genome-check`): an absent-copy call keeps only the members the genome
    does not explain themselves (genome_keeps on the read's own `de` / unexplained bases). A present copy's read swept into a
    cluster sits within delta of its own locus; it is no evidence of an absent copy."""
    return [m for m in members if genome_keeps(kind, des.get(m, 1.0), novs.get(m, 0))]


def nearest_copy(cons, sites, ref, maps, covered, min_overlap=MIN_OVERLAP):
    """(copy, distance) of the single copy closest to a cluster: the cluster's consensus base (its alt at an `alt` site, else
    the host's) against the copy's, over the cluster's covered positions the copy has a base at (a consensus alt the copy has
    no base at counts as a difference). A copy over < min_overlap of the covered positions is no candidate. (None, None)
    when no copy is."""
    best = (None, None)
    for n, cc in sorted(maps.items()):
        diff = cmp = 0
        for p in covered:
            alt = cons.get(p) == "alt"
            if p in cc:
                cmp += 1
                diff += cc[p] != (sites[p] if alt else ref[p])
            elif alt:
                cmp += 1
                diff += 1
        if cmp >= min_overlap * len(covered) and cmp and (best[1] is None or diff / cmp < best[1]):
            best = (n, diff / cmp)
    return best


def status_of(members, cons, sites, ref, maps, cov, vec, mosaic=False):
    """The prototype's status of a cluster against copy maps; est_div = its unexplained alts per covered base. mosaic: a cluster
    no copy explains, with >= MIN_ALT alts, whose NEAREST single copy is > DELTA away is `mosaic` (est_div = that distance):
    an absent copy whose alleles each occur in some other copy (GOLGA8, register 1298) is no single copy's."""
    alts = [p for p, a in cons.items() if a == "alt"]
    by = [n for n, cc in maps.items() if explains(cons, sites, ref, cc)]
    unexpl = [p for p in alts if not any(cc.get(p) == sites[p] for cc in maps.values())]
    cov_pos = {p for m in members[:50] for p in cov[m]}
    covered = len(cov_pos)
    link = PC.linkage(unexpl, vec) if len(unexpl) >= PC.MIN_ALT else None
    if by:
        st = "explained"
    elif len(unexpl) < PC.MIN_ALT or link is None:
        st = "unresolved"
    else:
        st = "unexplained" if link >= PC.MIN_LINK else "unlinked"
    div = len(unexpl) / covered if covered else 0.0
    if mosaic and st != "explained" and not (st == "unexplained" and div > DELTA) and len(alts) >= PC.MIN_ALT:
        _, d = nearest_copy(cons, sites, ref, maps, cov_pos)     # every cluster no copy explains, not only unresolved ones
        if d is not None and d > DELTA:                           # (round 3 review C2: KIR, unexplained at est_div 0.0032)
            return "mosaic", d
    return st, div


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


def graph_key(d, idxs):
    """A family graph's cache key: the content of its inputs (loci.fa and each copy's isoform GTF)."""
    h = hashlib.sha1()
    for name in ["loci.fa"] + [f"copy{i}.gtf" for i in idxs]:
        h.update(name.encode())
        h.update(open(f"{d}/{name}", "rb").read())
    return h.hexdigest()


def msga_orders(names):
    """Loci orders to retry `vg msga -a` in when FASTA order fails (round 9 amendment 1): reversed, then each rotation."""
    if len(names) < 2:
        return []
    return [names[::-1]] + [names[k:] + names[:k] for k in range(1, len(names))]


def family_graph(d, genome, copies):
    """vg graph of the family's loci (each in its transcript's orientation) with every curated isoform as a path -> GFA path."""
    os.makedirs(d, exist_ok=True)
    with open(f"{d}/loci.fa", "w") as fa:
        for c in copies:
            lo = max(0, c["start"] - LOCUS_PAD)
            hi = min(len(genome[c["contig"]]), c["end"] + LOCUS_PAD)
            raw = genome[c["contig"]][lo:hi]
            seq, _ = orient(raw, [], c["strand"], lo)
            fa.write(f">copy{c['idx']}\n{seq}\n")
            with open(f"{d}/copy{c['idx']}.gtf", "w") as g:
                for iso, ex in c["isoforms"].items():
                    for s, e in orient(raw, ex, c["strand"], lo)[1]:
                        g.write(f"copy{c['idx']}\tsim\texon\t{s + 1}\t{e}\t.\t+\t.\t"
                                f"transcript_id \"copy{c['idx']}|{iso}\"; gene_id \"copy{c['idx']}\";\n")
    cache = os.environ.get("O3_GRAPH_CACHE")              # (fastset.py) worlds with one catalog share the family graph
    key = graph_key(d, [c["idx"] for c in copies]) if cache else None
    if cache and os.path.exists(f"{cache}/{key}.gfa"):
        with open(f"{cache}/{key}.gfa") as src, open(f"{d}/s.gfa", "w") as dst:
            dst.write(src.read())
        return f"{d}/s.gfa"
    try:
        sh(f"vg msga -f {d}/loci.fa -t 2", binary_out=f"{d}/g.vg")
    except SystemExit:                    # vg msga (v1.73) aborts on some loci sets (UGT2B, round 3): build in FASTA order
        try:
            sh(f"vg msga -f {d}/loci.fa -t 2 -a", binary_out=f"{d}/g.vg")
        except SystemExit:                # round 9 amendment 1: KRTAP4 asserts in FASTA order too; other orders build
            seqs = PC.read_fa(f"{d}/loci.fa")
            for order in msga_orders(list(seqs)):
                with open(f"{d}/loci.order.fa", "w") as fh:
                    fh.writelines(f">{n}\n{seqs[n]}\n" for n in order)
                try:
                    sh(f"vg msga -f {d}/loci.order.fa -t 2 -a", binary_out=f"{d}/g.vg")
                    break
                except SystemExit:
                    continue
            else:
                return None               # no graph: the family gets no branch calls
    gtfs = " ".join(f"-n {d}/copy{c['idx']}.gtf" for c in copies)
    sh(f"vg rna {gtfs} -r -t 2 {d}/g.vg", binary_out=f"{d}/s.vg")
    sh(f"vg view {d}/s.vg > {d}/s.gfa")
    if cache:
        os.makedirs(cache, exist_ok=True)
        with open(f"{d}/s.gfa") as src, open(f"{cache}/{key}.gfa.{os.getpid()}", "w") as dst:
            dst.write(src.read())
        os.replace(f"{cache}/{key}.gfa.{os.getpid()}", f"{cache}/{key}.gfa")    # atomic: concurrent worlds never read half
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


def _clips(cigar):
    lead = trail = 0
    for op, n in cigar:
        if op not in (4, 5):
            break
        lead += n
    for op, n in reversed(cigar):
        if op not in (4, 5):
            break
        trail += n
    return lead, trail


def aligned_interval(rec, read_len):
    """The read interval an alignment record covers, in the read's ORIGINAL orientation. rec: {cigar (pysam tuples), reverse}."""
    lead, _ = _clips(rec["cigar"])
    e = lead + sum(n for op, n in rec["cigar"] if op in (0, 1, 7, 8))
    return (read_len - e, read_len - lead) if rec["reverse"] else (lead, e)


def primary_novel_intervals(rec, read_len, min_len=MIN_NOVEL_PART):
    """Read intervals (original orientation) the primary alignment does not explain: end clips >= min_len and dense error
    windows >= min_len (`--eqx` CIGAR)."""
    lead, trail = _clips(rec["cigar"])
    out = [(0, lead)] if lead >= min_len else []
    if trail >= min_len:
        out.append((read_len - trail, read_len))
    code = {0: "M", 1: "I", 2: "D", 7: "=", 8: "X"}
    cg = "".join(f"{n}{code[op]}" for op, n in rec["cigar"] if op in code)
    out += [x for x in dense_intervals(dict(qs=lead, cigar=cg)) if x[1] - x[0] >= min_len]
    if rec["reverse"]:
        out = [(read_len - e, read_len - s_) for s_, e in out]
    return sorted(out)


def genome_novel(prim, supps, read_len, min_len=MIN_NOVEL_PART):
    """Bases the genome does not explain: the primary's unexplained intervals minus what the read's supplementary alignments
    cover (a readthrough's clipped part aligned elsewhere IS explained), pieces >= min_len."""
    cover = sorted(aligned_interval(r, read_len) for r in supps)
    n = 0
    for s_, e in primary_novel_intervals(prim, read_len, min_len):
        pieces, cur = [], s_
        for cs, ce in cover:
            if ce <= cur or cs >= e:
                continue
            if cs > cur:
                pieces.append((cur, cs))
            cur = max(cur, ce)
        if cur < e:
            pieces.append((cur, e))
        n += sum(b - a for a, b in pieces if b - a >= min_len)
    return n


def genome_stats(bam_p, want):
    """{read: (de of its primary, genome_novel)} for the wanted reads that have a primary alignment."""
    import pysam
    prim, supp = {}, collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.query_name not in want or rec.is_unmapped or rec.is_secondary or not rec.cigartuples:
                continue
            r = dict(cigar=rec.cigartuples, reverse=rec.is_reverse)
            if rec.is_supplementary:
                supp[rec.query_name].append(r)
            else:
                prim[rec.query_name] = (rec.get_tag("de") if rec.has_tag("de") else 0.0, r, rec.infer_read_length())
    return {q: (de, genome_novel(r, supp[q], L)) for q, (de, r, L) in prim.items()}


def member_genome(members, gst, read_len):
    """Per member: (de, genome-unexplained bases); a member without a primary alignment is wholly unexplained (de 1.0, its
    whole length)."""
    de = [gst[m][0] if m in gst else 1.0 for m in members]
    nov = [gst[m][1] if m in gst else read_len[m] for m in members]
    return de, nov


_ONE_FAMILY = None


def _run_family(f):                                   # a fork worker's entry: the family closure main() set before forking
    return _ONE_FAMILY(f)


def main(argv):
    w = argv[0]
    only_i = "--arms" in argv and argv[argv.index("--arms") + 1] == "I"   # fastset.py: the integrated pipeline's arm only
    genome = Genome(f"{w}/genome.ref.fa")
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
    def unexplained(name):                             # S4 / G: o2graph.py run WORLD paths|union 2, first
        p = f"{w}/graph/{name}"
        return ({(r["read"], r["family"]) for r in csv.DictReader(open(p), delimiter="\t") if r["reason"] == "unexplained"}
                if os.path.exists(p) else None)
    o2_unexpl, o2u_unexpl = unexplained("o2.paths.m2.tsv"), unexplained("o2.union.m2.tsv")
    o2e_unexpl = unexplained("o2.unionp.m2.tsv")      # H: union O2 on paths padded with each copy's genome at both ends
    o2k_unexpl = unexplained("o2.unionf.m2.tsv")      # K: what no family's graph explains, nor the family's own loci
    o2i_unexpl = unexplained("o2.integ.m2.tsv")       # I (round 7): what O2 on O1's RNA-seeded catalog cannot explain
    by_host = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        if r["catalog_copy_idx"].isdigit():
            by_host[(r["family_id"], int(r["catalog_copy_idx"]))].append(r["read_name"])
    out_dir = f"{w}/graph"
    os.makedirs(out_dir, exist_ok=True)
    tmp = f"{out_dir}/tmp"
    calls = []                                     # (arm, family, kind, host, cluster, members, est_div)

    def one_family(f):                             # families are independent: O3_PROCS of them at a time, in family order
        calls = []
        tmp = f"{out_dir}/{f}.tmp"
        cps = fam[f]
        pool = [n for n in family_pool(by_host, f, o2e_unexpl) if n in reads]
        pool_k = [n for n in family_pool(by_host, f, o2k_unexpl) if n in reads]
        pool_i = [n for n in family_pool(by_host, f, o2i_unexpl) if n in reads]
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
            if not only_i:                                    # R0/S1/S4/G: every host read; arm I needs only its own pool
                cov = PC.align("host", ref, {n: reads[n] for n in names if n in reads}, tmp)
                m_rep = PC.align("host", ref, {f"copy{c['idx']}": c["seq"] for c in cps}, tmp)
            iso_seqs = {f"copy{c['idx']}|{k}": "".join(g2[a:b] for a, b in ex) for c in cps
                        for g2 in [genome[c["contig"]]] for k, ex in c["isoforms"].items()}
            m_iso_each = PC.align("host", ref, iso_seqs, tmp)
            m_iso = {f"copy{c['idx']}": union_map([m for k, m in m_iso_each.items() if k.startswith(f"copy{c['idx']}|")])
                     for c in cps}
            host = ",".join(str(c["idx"]) for c in sorted(cls, key=lambda c: c["idx"]))
            sites, clusters, vec = PC.cluster_host(ref, cov) if not only_i else (None, [], None)
            for k, (members, cons) in enumerate(sorted(clusters, key=lambda x: -len(x[0]))):
                for arm, maps in (("R0", m_rep), ("S1", m_iso)):
                    st, div = status_of(members, cons, sites, ref, maps, cov, vec)
                    if st == "unexplained" and div > DELTA:
                        calls.append([arm, f, "psv", host, k, sorted(members), div])
                st, div = status_of(members, cons, sites, ref, m_iso, cov, vec, mosaic=True)
                if st in ("unexplained", "mosaic") and div > DELTA:
                    calls.append(["S1m", f, "psv", host, k, sorted(members), div])
            un = restrict([n for n in names if n in cov], f, o2_unexpl) if not only_i else []
            if len(un) >= PC.MIN_SUPPORT:
                cov_u = {n: cov[n] for n in un}
                sites_u, clusters_u, vec_u = PC.cluster_host(ref, cov_u)
                for k, (members, cons) in enumerate(sorted(clusters_u, key=lambda x: -len(x[0]))):
                    st, div = status_of(members, cons, sites_u, ref, m_iso, cov_u, vec_u)
                    if st == "unexplained" and div > DELTA:
                        calls.append(["S4p", f, "psv", host, k, sorted(members), div])
            un = restrict([n for n in names if n in cov], f, o2u_unexpl) if not only_i else []
            if len(un) >= PC.MIN_SUPPORT:
                cov_u = {n: cov[n] for n in un}
                sites_u, clusters_u, vec_u = PC.cluster_host(ref, cov_u)
                for k, (members, cons) in enumerate(sorted(clusters_u, key=lambda x: -len(x[0]))):
                    st, div = status_of(members, cons, sites_u, ref, m_iso, cov_u, vec_u, mosaic=True)
                    if st in ("unexplained", "mosaic") and div > DELTA:
                        calls.append(["Gp", f, "psv", host, k, sorted(members), div])
            if len(pool) >= PC.MIN_SUPPORT and not only_i:    # H: the family's pooled unexplained reads on this host
                cov_h = PC.align("host", ref, {n: reads[n] for n in pool}, tmp)
                if len(cov_h) >= PC.MIN_SUPPORT:
                    sites_h, clusters_h, vec_h = PC.cluster_host(ref, cov_h)
                    for k, (members, cons) in enumerate(sorted(clusters_h, key=lambda x: -len(x[0]))):
                        st, div = status_of(members, cons, sites_h, ref, m_iso, cov_h, vec_h, mosaic=True)
                        if st in ("unexplained", "mosaic") and div > DELTA:
                            calls.append(["Hp", f, "psv", host, k, sorted(members), div])
            if len(pool_k) >= PC.MIN_SUPPORT and not only_i:  # K: as H, on the pool no family's graph nor the loci explain
                cov_k = PC.align("host", ref, {n: reads[n] for n in pool_k}, tmp)
                if len(cov_k) >= PC.MIN_SUPPORT:
                    sites_k, clusters_k, vec_k = PC.cluster_host(ref, cov_k)
                    for k, (members, cons) in enumerate(sorted(clusters_k, key=lambda x: -len(x[0]))):
                        st, div = status_of(members, cons, sites_k, ref, m_iso, cov_k, vec_k, mosaic=True)
                        if st in ("unexplained", "mosaic") and div > DELTA:
                            calls.append(["Kp", f, "psv", host, k, sorted(members), div])
            if len(pool_i) >= PC.MIN_SUPPORT:                 # I: as H, on what the integrated O1 -> O2 leaves unexplained
                cov_i = PC.align("host", ref, {n: reads[n] for n in pool_i}, tmp)
                if len(cov_i) >= PC.MIN_SUPPORT:
                    sites_i, clusters_i, vec_i = PC.cluster_host(ref, cov_i)
                    for k, (members, cons) in enumerate(sorted(clusters_i, key=lambda x: -len(x[0]))):
                        st, div = status_of(members, cons, sites_i, ref, m_iso, cov_i, vec_i, mosaic=True)
                        if st in ("unexplained", "mosaic") and div > DELTA:
                            calls.append(["Ip", f, "psv", host, k, sorted(members), div])
        names = sorted({n for c in cps for n in by_host.get((f, c["idx"]), []) if n in reads})
        if names:
            gfa = family_graph(f"{out_dir}/{f}", genome, cps)
            if gfa is None:
                print(f"[o3_graph] {w}: no vg graph for {f} (vg msga failed in every order): no branch calls", flush=True)
                return calls
            br, _ = branch_calls(f"{out_dir}/{f}", gfa, names, reads)
            for k, m in enumerate(br):
                calls.append(["S2", f, "branch", "-", k, m, 0.0])
        return calls
    global _ONE_FAMILY
    _ONE_FAMILY = one_family
    procs = int(os.environ.get("O3_PROCS", "1"))
    if procs > 1:
        import multiprocessing
        with multiprocessing.get_context("fork").Pool(procs) as pool:
            for part in pool.map(_run_family, sorted(fam)):
                calls += part
    else:
        for f in sorted(fam):
            calls += one_family(f)
    calls = dedupe(calls)
    want = {m for c in calls for m in c[5]}
    gst = genome_stats(f"{w}/E0.bam", want)
    rows_out = []
    member_check = "--member-genome-check" in argv            # arm I only: the genome check per read (member_filter)
    for arm, f, kind, host, k, m, div in calls:
        des, nov = member_genome(m, gst, {x: len(reads[x]) for x in m})
        mde = statistics.median(des) if des else 0.0
        mnov = statistics.median(nov) if nov else 0
        keep = genome_keeps(kind, mde, mnov)
        for a in arms_of(arm, keep):
            mm = m
            if member_check and a == "I":
                mm = member_filter(kind, m, dict(zip(m, des)), dict(zip(m, nov)))
                if len(mm) < PC.MIN_SUPPORT:
                    continue
            rows_out.append(dict(arm=a, family=f, kind=kind, host=host, cluster=k, n_reads=len(mm), est_div=round(div, 5),
                                 genome_median_de=round(mde, 5), genome_median_novel=mnov, members=",".join(mm)))
    if only_i:
        rows_out = [r for r in rows_out if r["arm"] in ("In", "I")]
    with open(f"{out_dir}/calls.tsv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["arm", "family", "kind", "host", "cluster", "n_reads", "est_div", "genome_median_de",
                                            "genome_median_novel", "members"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows_out)
    for x in (".t.fa", ".q.fa"):
        for t in [tmp] + [f"{out_dir}/{f}.tmp" for f in fam]:
            if os.path.exists(t + x):
                os.remove(t + x)
    print(f"[o3_graph] {w}: " + ", ".join(f"{a} {sum(1 for r in rows_out if r['arm'] == a)}" for a in ARMS) + " calls", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
