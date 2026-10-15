#!/usr/bin/env python3
"""O3 on real reads (prereg docs/PREREG_o3_realdata_2026-10-15.md): pick families, attribute reads to copies, score.

The scoring rules are the simulation's (`bench/o3_unexplained/score.py`): a copy absent from the genome is RECOVERED when the
`unexplained` PSV cluster holding the most of its reads that O2 placed on the family holds >= RECALL of them with purity
(its reads / cluster reads) >= PURITY; the cluster's label is `copy` when est_div > DELTA, else `allele_or_copy`.
"""
import bisect
import collections
import random

RECALL, PURITY, DELTA, CLEAR = 0.8, 0.9, 0.00958, 0.5


def _overlap(a, b):
    """Bases shared by two sorted, disjoint interval lists (0-based half-open)."""
    i = j = n = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            n += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def pick_families(copies, catalog, exclude, min_copies=3, min_reads=20, max_span=200_000, cap=60, seed=1,
                  min_expressed=None, exclude_families=()):
    """Families of `copies` (dicts: family, idx, chrom, start, end, exons, n = reads on the copy's exons) with >= min_copies
    copies, every copy >= min_reads reads and span <= max_span, no copy overlapping an `exclude` interval (chrom, s, e); the
    copy last by (chrom, start) is masked, and a family whose masked span overlaps any other `catalog` copy is dropped.
    `min_expressed` = k relaxes the read rule: >= k copies with >= min_reads reads, the last of THOSE is masked, and the
    silent copies stay survivors. Families in `exclude_families` are skipped.
    More than `cap` families: a seeded sample. -> [{family, masked, survivors}] sorted by family."""
    by_fam = collections.defaultdict(list)
    for c in copies:
        by_fam[c["family"]].append(c)
    keep = []
    for fam, cs in sorted(by_fam.items()):
        if fam in exclude_families or len(cs) < min_copies or any(c["end"] - c["start"] > max_span for c in cs):
            continue
        expressed = [c for c in cs if c["n"] >= min_reads]
        if not expressed or ((len(expressed) < len(cs)) if min_expressed is None else (len(expressed) < min_expressed)):
            continue
        if any(c["chrom"] == x and c["start"] < e and s < c["end"] for c in cs for x, s, e in exclude):
            continue
        order = lambda c: (c["chrom"], c["start"], c["idx"])  # noqa: E731
        m = sorted(expressed, key=order)[-1]
        if any(o["chrom"] == m["chrom"] and o["start"] < m["end"] and m["start"] < o["end"]
               and not (o["family"] == fam and o["idx"] == m["idx"]) for o in catalog):
            continue
        keep.append(dict(family=fam, masked=m, survivors=[c for c in sorted(cs, key=order) if c is not m]))
    if len(keep) > cap:
        keep = sorted(random.Random(seed).sample(keep, cap), key=lambda p: p["family"])
    return keep


def reads_by_copy(alns, copies):
    """alns: [(read, chrom, aligned blocks)] (primary alignments); copies: {(family, idx): (chrom, exons)}. -> {read: copy}
    for reads with an aligned base on exactly one copy's exons."""
    by_chrom = collections.defaultdict(list)
    for k, (chrom, ex) in copies.items():
        by_chrom[chrom].append((ex[0][0], ex[-1][1], k, sorted(ex)))
    for v in by_chrom.values():
        v.sort()
    starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
    out = {}
    for read, chrom, blocks in alns:
        v = by_chrom.get(chrom, [])
        hi = bisect.bisect_right(starts.get(chrom, []), blocks[-1][1])
        hit = {k for s, e, k, ex in v[:hi] if e > blocks[0][0] and _overlap(sorted(blocks), ex) > 0}
        if len(hit) == 1:
            out[read] = hit.pop()
    return out


def sample_per_copy(read_copy, n=500, seed=1):
    """At most n reads per copy, a seeded sample over the reads sorted by name. -> {read: copy}."""
    per = collections.defaultdict(list)
    for r, k in sorted(read_copy.items()):
        per[k].append(r)
    rng, out = random.Random(seed), {}
    for k in sorted(per):
        rs = per[k] if len(per[k]) <= n else sorted(rng.sample(per[k], n))
        out.update((r, k) for r in rs)
    return out


def _fam_unexplained(fam, clusters):
    return [c for c in clusters if c["family"] == fam and c["status"] == "unexplained"]


def score_absent(fam, absent, placed, clusters, members):
    """absent: the reads of the copy missing from the genome; placed: the reads O2 placed on `fam`; members:
    {(family, host, cluster): read names}. -> n_absent, n_placed, recovered, recall, purity, est_div, label."""
    ap = absent & placed
    best, bc = 0, None
    for c in _fam_unexplained(fam, clusters):
        got = len(members.get((c["family"], c["host"], c["cluster"]), set()) & ap)
        if got > best:
            best, bc = got, c
    recall = best / len(ap) if ap else 0.0
    purity = best / int(bc["n_reads"]) if bc else 0.0
    rec = int(bc is not None and recall >= RECALL and purity >= PURITY)
    est = float(bc["est_div"]) if bc else 0.0
    return dict(n_absent=len(absent), n_placed=len(ap), recovered=rec, recall=round(recall, 3), purity=round(purity, 3),
                est_div=est, label=("copy" if est > DELTA else "allele_or_copy") if rec else "-")


def cleared(fam, absent, placed, clusters, members):
    """After the copy is added back: no `unexplained` cluster of `fam` holds >= CLEAR of the absent copy's placed reads."""
    ap = absent & placed
    if not ap:
        return True
    return all(len(members.get((c["family"], c["host"], c["cluster"]), set()) & ap) < CLEAR * len(ap)
               for c in _fam_unexplained(fam, clusters))


def control_flags(families, clusters):
    """Per one-copy family: `copy` if an unexplained cluster has est_div > DELTA, else `allele_or_copy` if any unexplained
    cluster, else `-`."""
    out = {}
    for f in families:
        u = _fam_unexplained(f, clusters)
        out[f] = "copy" if any(float(c["est_div"]) > DELTA for c in u) else "allele_or_copy" if u else "-"
    return out


def het_bin(d):
    return "0" if d == 0 else "<=0.1%" if d <= 0.001 else "<=0.5%" if d <= 0.005 else "<=delta" if d <= DELTA else ">delta"
