#!/usr/bin/env python3
"""Graph-based O2 (prereg docs/archive/2026-10/merged/PREREG_o3_graph_real_2026-10-15.md): assign a family's reads to its copies on the family's
variation graph, where every curated isoform's exons carry their PSVs (a representative lacks some of them).

O2's domain (user, 2026-10-19): AS-tied multimapping reads only -- a rival alignment within 10 AS of the best (Eichler's rule).
A read with a clear best placement is minimap2's; applying O2 to it is tautological. `run W intega` (fastset --o2-scope tied
--as-slack 9) confines the calls to that domain; the other modes predate it and decide every hosted read.

    o2graph.py run WORLD [paths|nodes] [MARGIN]   WORLD/graph/o2.tsv; `paths` (default): reads aligned to every curated isoform
                                                  path of every catalogued copy; `nodes`: allele nodes of a free GraphAligner
                                                  alignment (needs detect.py's graph/<family>/s.gfa + aln.gaf)
    o2graph.py score OUT.json O2_NAME WORLD ...   O2 (graph/O2_NAME) against the truth, the current method beside it
    (run ... reps MARGIN: the ablation with each copy's representative path only; run ... union MARGIN: every copy's own
    isoforms plus every paralog isoform structure projected onto its locus (union_structures); run ... unione MARGIN: union,
    read bases past a path's end not counted (_paf_dist ends=True); run ... unionp MARGIN: union, every path's two ends
    extended by PATH_PAD bases of the copy's genome (padded); run ... locus MARGIN: the
    annotation-free baseline on seqsupport's genomic-locus distances; run ... unionf MARGIN: padded union with the locus
    fallback and the cross-family guard (fallback_call); run ... flag BASE: o2.<BASE> with the reads of flagged PSV clusters
    reported as unannotated-copy candidates (flag_reads, apply_flag); run ... integ MARGIN: padded union on O1's seeded catalog
    (o1seed.py; run_integrated); run ... integl MARGIN: integ with the CIGAR-aware locus fallback for the reads it leaves
    unexplained (run_locus_fallback); run ... integc MARGIN: every read on min(path, genome) distances
    (run_locus_fallback mode "all"); run ... integs MARGIN: mode "seeded"; run ... intega MARGIN [--base integv]: O2 only on
    the AS-tied reads (run_scope, as_scope); outputs graph/o2.<method>.m<margin>.tsv)

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
MARGIN, MIN_COVER, DELTA, END_SLACK, PATH_PAD, FLAG_MIN_READS = 2, 0.8, 0.00958, 5, 50, 20
EXT_CAP, EXT_BAND = 50, 30                               # --distance gapcomp: clipped ends aligned up to EXT_CAP bases
_INDEL = re.compile(r"(\d+)[ID]")
NOVEL_INS = 30                     # locus fallback: an insertion or unaligned end this long is sequence no copy holds
UNION_EXACT_MICRO = 30             # --union-exact: a terminal exon this short may be missing from a projection
ELSEWHERE = -1                     # locus fallback: the read's best alignment off every seeded copy of the family
DROPPED = -2                       # --twin-guard: an alignment on a catalogued copy O1 dropped (silent), not off the catalogue
ALPHA = 1e-3                       # --o2-significance: a path call's misassignment bound (copy_assign's certificate alpha)
HP_RUN = 3                         # a homopolymer run, as the v2 profiles' error model places homopolymer indels (calibrate.py)
EDIT_RATE_READS = 3000             # --o2-significance: uniquely placed reads the world's edit rates are counted on
_CIGAR_OPS = re.compile(r"(\d+)([MIDNSHP=X])")


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


def gate(floor=0.0):
    """The unexplained gate's divergence (round 11, `--error-floor`): DELTA plus the world's sequencing-error floor -- DELTA is
    an error-free divergence (p99 between KB3781's haplotype assemblies), a read's distance adds its errors."""
    return DELTA + floor


_TWIN_ABSTAIN = False                                      # --twin-abstain (any run mode): see decide


def decide(dist, qlen, margin=MARGIN, max_div=DELTA, twins=None, twin_abstain=None, gate=None):
    """dist: {copy: edits + unaligned bases of the read against the copy's best-fitting isoform path}. -> (copy or None, reason):
    `assigned` when the best copy beats the next by >= margin, else `tied` / `thin`; `unexplained` when even the best copy is
    more than max_div away (a read of no catalogued copy). twins ({new copy: curated copy}, O1 --twins: a new copy within
    delta of a curated one): a twin takes a read from its curated copy only by the margin; otherwise the curated copy competes
    with the better of the two distances (the annotation is the prior a within-delta twin must overcome). twin_abstain
    (--twin-abstain, 2026-10-17: in the O2 proof's OUTSIDE worlds VCY's uncatalogued twin is real, and the prior gave its reads
    to the curated copy): when neither twin beats the other by the margin, the read is `tied_twin` and abstained. gate
    (--length-gate): gate(edits, qlen) -> within reach, in place of edits <= max_div x qlen (within_reach)."""
    twin_abstain = _TWIN_ABSTAIN if twin_abstain is None else twin_abstain
    if twins:
        dist = dict(dist)
        for t, c in twins.items():
            if t in dist and c in dist and dist[c] - dist[t] < margin:
                if twin_abstain and dist[t] - dist[c] < margin:
                    return None, "tied_twin"
                dist[c] = min(dist[c], dist.pop(t))
    if not dist:
        return None, "unexplained"
    srt = sorted(dist.items(), key=lambda kv: (kv[1], kv[0]))
    if (srt[0][1] > max_div * qlen) if gate is None else not gate(srt[0][1], qlen):
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


def end_cost(p, t, anchored_at_end):
    """Edit distance of a whole clipped read end `p` against path bases `t`, anchored where the end meets the aligned core (a
    prefix at t's end, a suffix at t's start) and free on the far side; read bases beyond the path cost one each. An end longer
    than EXT_CAP costs its length (as clipping counts it: a different structure, not errors)."""
    m = len(p)
    if m > EXT_CAP:
        return m
    if not m:
        return 0
    if anchored_at_end:
        p, t = p[::-1], t[::-1]
    # the DP D[i][j] (i read bases, j path bases; D[0][j] = j, D[i][0] = i), min over j of D[m][j], one path base at a time
    # with Myers' bit vectors (Hyyro's global variant: a +1 horizontal delta enters at row 0)
    full, top = (1 << m) - 1, 1 << (m - 1)
    peq = {}
    for i, c in enumerate(p):
        peq[c] = peq.get(c, 0) | (1 << i)
    pv, mv, score = full, 0, m
    best = m
    for c in t:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = ((((eq & pv) + pv) & full) ^ pv) | eq
        ph = (mv | ~(xh | pv)) & full
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        if score < best:
            best = score
        ph = ((ph << 1) | 1) & full
        mh = (mh << 1) & full
        pv = (mh | ~(xv | ph)) & full
        mv = ph & xv
    return best


def gap_compressed(line, read, target, struct=False):
    """A read's distance to a path from one minimap2 -c PAF line, CIGAR-aware (2026-10-16): mismatches plus indel EVENTS (an
    indel counts once whatever its length, as minimap2's `de` counts it: real reads' errors are mostly homopolymer indels, and a
    retained intron is one event, not its length), plus each clipped read end aligned against the path bases next to the core
    (end_cost: a clipped base is not an edit unless it disagrees with the path)."""
    f = line.rstrip("\n").split("\t")
    qlen, qs, qe, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8])
    if f[4] == "-":
        read, qs, qe = read[::-1].translate(str.maketrans("ACGTN", "TGCAN")), qlen - qe, qlen - qs
    pre, suf = read[:qs], read[qe:]
    return (_gc_core(f, struct) + end_cost(pre, target[max(0, ts - len(pre) - EXT_BAND):ts], True)
            + end_cost(suf, target[te:te + len(suf) + EXT_BAND], False))


def _gc_core(f, struct=False):
    """The aligned core's share of gap_compressed: mismatches plus indel events (NM less the indel bases, plus the events).
    struct (--distance gapstruct, round 12 review C1): only indels shorter than NOVEL_INS are compressed; a longer one (a
    retained intron, a novel exon) counts its length, as the same sequence costs when the aligner clips it instead."""
    nm = int(next(x for x in f[12:] if x.startswith("NM:i:"))[5:])
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    indels = [int(n) for n in _INDEL.findall(cg)]
    if struct:
        indels = [n for n in indels if n < NOVEL_INS]
    return nm - sum(indels) + len(indels)


def _paf_dist_gc(paf, reads, tx, struct=False):
    """{read: {copy: min gap_compressed distance}} (--distance gapcomp); tx: {path name: sequence as aligned}. Exactly the
    minimum over every line, scored lazily: a copy's lines in ascending core distance, stopping once the core alone reaches
    the best full distance (the clipped ends only add)."""
    by = collections.defaultdict(list)
    for ln in paf.splitlines():
        f = ln.split("\t")
        by[(f[0], f[5].split("|")[0])].append((_gc_core(f, struct), ln, f[5]))
    dist = collections.defaultdict(dict)
    for (q, c), cand in by.items():
        best = 10 ** 9
        for core, ln, t in sorted(cand, key=lambda x: x[0]):
            if core >= best:
                break
            best = min(best, gap_compressed(ln, reads[q], tx[t], struct))
        dist[q][c] = best
    return dist


def hp_indel(target, pos, n, ins=None):
    """Is an indel inside a homopolymer run of the target -- the deletion of target[pos:pos + n], or the insertion of `ins`
    before target[pos]: all its bases one base, whose run in the target at the indel is at least HP_RUN long?"""
    seq = target[pos:pos + n] if ins is None else ins
    if not seq or seq.count(seq[0]) != len(seq):
        return False
    b, lo, hi = seq[0], pos, pos + (n if ins is None else 0)
    while lo > 0 and target[lo - 1] == b:
        lo -= 1
    while hi < len(target) and target[hi] == b:
        hi += 1
    return hi - lo >= HP_RUN


def gap_weighted(line, read, target, wts):
    """gap_compressed with each edit priced by its type (--o2-significance, 2026-10-17): -ln of the world's own rate of it
    (edit_weights) -- a mismatch, or an edit of a clipped end, wts["X"]; an indel event inside a homopolymer run (hp_indel)
    wts["hp"]; any other indel event wts["indel"]. A difference of two copies' distances is then the read's log-likelihood
    ratio between them (each edit an independent error of its type)."""
    f = line.rstrip("\n").split("\t")
    qlen, qs, qe, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8])
    if f[4] == "-":
        read, qs, qe = read[::-1].translate(str.maketrans("ACGTN", "TGCAN")), qlen - qe, qlen - qs
    nm = int(next(x for x in f[12:] if x.startswith("NM:i:"))[5:])
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    ti, qi, indel_bases, hp, ind = ts, qs, 0, 0, 0
    for n, op in _CIGAR_OPS.findall(cg):
        n = int(n)
        if op in "M=X":
            ti, qi = ti + n, qi + n
        elif op in "ID":
            if hp_indel(target, ti, n, read[qi:qi + n] if op == "I" else None):
                hp += 1
            else:
                ind += 1
            indel_bases += n
            ti, qi = (ti, qi + n) if op == "I" else (ti + n, qi)
    pre, suf = read[:qs], read[qe:]
    ends = (end_cost(pre, target[max(0, ts - len(pre) - EXT_BAND):ts], True)
            + end_cost(suf, target[te:te + len(suf) + EXT_BAND], False))
    return (nm - indel_bases + ends) * wts["X"] + hp * wts["hp"] + ind * wts["indel"]


def _paf_dist_w(paf, reads, tx, wts):
    """{read: {copy: min gap_weighted distance}}, exactly the minimum over every line: a copy's lines in ascending gap-compressed
    core, stopping once the cheapest edit's weight times the core reaches the best weighted distance."""
    wmin = min(wts.values())
    by = collections.defaultdict(list)
    for ln in paf.splitlines():
        f = ln.split("\t")
        by[(f[0], f[5].split("|")[0])].append((_gc_core(f), ln, f[5]))
    dist = collections.defaultdict(dict)
    for (q, c), cand in by.items():
        best = float("inf")
        for core, ln, t in sorted(cand, key=lambda x: x[0]):
            if wmin * core >= best:
                break
            best = min(best, gap_weighted(ln, reads[q], tx[t], wts))
        dist[q][c] = best
    return dist


def edit_counts(w, n_reads=EDIT_RATE_READS):
    """The edits of WORLD/E0.bam's first n_reads uniquely placed reads (primary, MAPQ 60) against WORLD/genome.ref.fa:
    aligned bases, mismatches (X, or NM less the indel bases without --eqx), indel events under NOVEL_INS inside a homopolymer
    run (hp) and elsewhere (indel), and the homopolymer runs on the reference the alignments cover (hp_runs: each a chance
    of a homopolymer indel). Sequencing errors and the donor's alleles both: what a read of its own copy differs by."""
    import pysam
    fa = pysam.FastaFile(f"{w}/genome.ref.fa")
    c = collections.Counter(aligned=0, X=0, hp=0, indel=0, hp_runs=0)
    seen, pad = 0, HP_RUN + 10
    runs = re.compile("|".join(f"{b}{{{HP_RUN},}}" for b in "ACGT"))
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if (rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.mapping_quality < 60
                    or not rec.query_sequence):
                continue
            seen += 1
            if seen > n_reads:
                break
            lo = max(0, rec.reference_start - pad)
            ref = fa.fetch(rec.reference_name, lo, rec.reference_end + pad).upper()
            q = rec.query_sequence.upper()
            ti, qi, indel_bases, eqx, segs = rec.reference_start - lo, 0, 0, False, []
            for op, n in rec.cigartuples:
                if op in (0, 7, 8, 2):
                    if segs and segs[-1][1] == ti:
                        segs[-1][1] = ti + n
                    else:
                        segs.append([ti, ti + n])
                if op in (0, 7, 8):
                    c["aligned"] += n
                    eqx |= op in (7, 8)
                    c["X"] += n if op == 8 else 0
                    ti, qi = ti + n, qi + n
                elif op in (1, 2):
                    if n < NOVEL_INS:
                        c["hp" if hp_indel(ref, ti, n, q[qi:qi + n] if op == 1 else None) else "indel"] += 1
                    indel_bases += n
                    ti, qi = (ti, qi + n) if op == 1 else (ti + n, qi)
                elif op == 3:
                    ti += n
                elif op == 4:
                    qi += n
            if not eqx and rec.has_tag("NM"):
                c["X"] += rec.get_tag("NM") - indel_bases
            c["hp_runs"] += sum(1 for a, b in segs for _ in runs.finditer(ref, a, b))
    return c


def weights_from_counts(c):
    """{edit type: -ln of its rate} from edit_counts, each rate by the rule of succession: a mismatch to one given base (a third
    of the mismatch rate per aligned base), a homopolymer indel per homopolymer run, another indel per aligned base."""
    import math
    return {"X": max(0.0, math.log(3 * (c["aligned"] + 2) / (c["X"] + 1))),
            "hp": max(0.0, math.log((c["hp_runs"] + 2) / (c["hp"] + 1))),
            "indel": max(0.0, math.log((c["aligned"] + 2) / (c["indel"] + 1)))}


def poisson_quantile(lam, p):
    """The smallest n with P(Poisson(lam) <= n) >= p."""
    import math
    n, term = 0, math.exp(-lam)
    cdf = term
    while cdf < p:
        n += 1
        term *= lam / n
        cdf += term
    return n


def within_reach(edits, qlen, floor, alpha=ALPHA):
    """--length-gate (2026-10-17, round 16's B1 shortfall): can a read of qlen bases carry `edits` from its own copy? Allelic
    divergence up to delta (the p99 between haplotypes) plus sequencing errors up to the 1 - alpha quantile of
    Poisson(floor x qlen). The fixed gate, (delta + floor) x qlen, puts the errors at their mean: a 281-bp read was allowed 3.6
    edits, and round 16's short reads with 4-9 edits (SLC35G.1, OR13C.0, TAS2R.0) were left unexplained where the current
    method placed them. For long reads the two gates converge."""
    return edits <= DELTA * qlen + poisson_quantile(floor * qlen, 1 - alpha)


def rescue(reason, pdist, wdist, qlen, floor, twins=None, margin=MARGIN, length_gate=False, significance_margin=False,
           realigned=False, sdist=None):
    """The genome step's second look at a read path O2 did not assign (2026-10-17, round 16): -> (copy idx, "length_gate" |
    "significant") or None. length_gate (--length-gate): an `unexplained` read is decided again on its path distances with
    within_reach for the gate. significance_margin (--significance-margin): a `thin` or `tied` read goes to the copy of least
    weighted distance (gap_weighted) when that copy beats every rival by ln(R / alpha) -- one substitution (7.6-7.8) does so
    against one or two rivals, one homopolymer indel (3.9-4.6) never; round 16's SPAG.1 reads were 3 vs 4 edits, 14.8 vs
    22.5 weighted, with 2 rivals. realigned (--short-path-align): an `unexplained` read the path step left without any
    alignment, now aligned (short_path_distances), is decided on its path distances as path O2 would have ("short_path").
    sdist (--struct-rescue): a `thin` or `tied` read's structure-aware path distances (struct_distances); a copy that wins on
    them by the margin takes it ("struct_aware"). The call then goes through the genome step's checks like any path call."""
    pd = {int(k): v for k, v in pdist.items()}
    if not pd:
        return None
    reach = (lambda k, n: within_reach(k, n, floor)) if length_gate else None
    if reason == "unexplained":
        if not (length_gate or realigned):
            return None
        c, why = decide(pd, qlen, margin, max_div=gate(floor), twins=twins, gate=reach)
        if c is not None:
            return c, "short_path" if realigned else "length_gate"
        if why not in ("tied", "thin"):
            return None
        reason = why
    if sdist and reason in ("tied", "thin"):
        c, why = decide({int(k): v for k, v in sdist.items()}, qlen, margin, max_div=gate(floor), twins=twins, gate=reach)
        if c is not None:
            return c, "struct_aware"
    if significance_margin and reason in ("tied", "thin") and wdist:
        wd = {str(k): v for k, v in wdist.items()}
        c = min(wd, key=lambda k: (wd[k], int(k)))
        if (int(c) in pd and (reach(pd[int(c)], qlen) if reach else pd[int(c)] <= gate(floor) * qlen)
                and significant(c, {str(k): v for k, v in pd.items()}, wd, twins, margin)):
            return int(c), "significant"
    return None


SHORT_SCORE = 40                   # --short-path-align: minimap2's minimum DP score (its map-ont default; map-hifi's is 200)


def short_path_distances(d, fam, reads, wts=None):
    """--short-path-align (2026-10-17, round 16): {read: {copy idx as str: path distance}} and the weighted ones (wts given)
    of `reads` ({name: sequence}) against the family's paths (d/fam.tx.fa), aligned as the path step aligns them but with
    the minimum DP score at SHORT_SCORE: map-hifi's -s 200 is beyond a read under ~200 bp, which then has no path alignment at
    all (round 16: OR4A.5's and GGTLC.1's 150-200-bp reads; GGTLC.1's tie on the genome while its own isoform's junction
    separates them on the paths)."""
    import subprocess
    tx = _fasta(f"{d}/{fam}.tx.fa")
    fa = f"{d}/{fam}.short.{os.getpid()}.fa"
    with open(fa, "w") as fh:
        fh.writelines(f">{n}\n{q}\n" for n, q in sorted(reads.items()))
    args = ["minimap2", "-c", "-x", "map-hifi", "-s", str(SHORT_SCORE), "-N", str(max(200, 2 * len(tx))), "-p", "0.3",
            "--secondary=yes", "-t", "1", f"{d}/{fam}.tx.fa", fa]
    paf = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    os.remove(fa)
    dist = path_dists(_paf_dist_gc(paf, reads, tx))
    return dist, (path_dists(_paf_dist_w(paf, reads, tx, wts)) if wts else {})


def _fasta(p):
    out, name = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            name = ln[1:].strip()
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def struct_distances(d, fam, reads, struct=True):
    """--struct-rescue (2026-10-17, round 16's thin calls): {read: {copy idx as str: path distance}} of `reads` ({name:
    sequence}) against the family's paths (d/fam.tx.fa), aligned exactly as the path step aligns them (through O2_ALN_CACHE
    when it is set), scored by gap_compressed with struct (an indel >= NOVEL_INS counts its length). Gap compression counts a
    whole missing exon as one insertion, so a read's best path could be an isoform lacking one of its exons, and the
    paralogous bases in that exon went uncompared: SPAG.1's reads were 3 vs 4 on copy 1's exon-skipping isoform (a 76-bp
    insertion), 4 vs 9 on their own isoform."""
    import subprocess
    tx = _fasta(f"{d}/{fam}.tx.fa")
    mm2 = ["minimap2", "-c", "-x", "map-hifi", "-N", str(max(200, 2 * len(tx))), "-p", "0.3", "--secondary=yes"]
    if os.environ.get("O2_ALN_CACHE"):
        res = cached_minimap2(os.environ["O2_ALN_CACHE"], mm2, f"{d}/{fam}.tx.fa", reads, d)
    else:
        fa = f"{d}/{fam}.struct.{os.getpid()}.fa"
        with open(fa, "w") as fh:
            fh.writelines(f">{n}\n{q}\n" for n, q in sorted(reads.items()))
        res = subprocess.run(mm2 + ["-t", os.environ.get("O2_THREADS", "1"), f"{d}/{fam}.tx.fa", fa], capture_output=True,
                             text=True, check=True).stdout
        os.remove(fa)
    return path_dists(_paf_dist_gc(res, reads, tx, struct=struct))


def significant(call, pdist, wdist, twins=None, margin=MARGIN, alpha=ALPHA):
    """--o2-significance (2026-10-17, the O2 margin): a path call stands when its weighted distance (gap_weighted) beats every
    rival's by ln(R / alpha), R the rivals with a path alignment -- the read's likelihood ratio against each rival at least
    R / alpha, so a misassignment has probability at most alpha (Bonferroni over the rivals, as copy_assign's certificate
    divides alpha by them). The fixed margin counted two homopolymer indels as two substitutions: r12 A119b's DUX4 errors,
    2 edits of 1-bp indels against 63 rivals. pdist / wdist: {copy idx as str: path / weighted distance}; a twin folded into
    the call (decide's rule) is no rival; no weighted distance for the call: no gate."""
    import math
    c = str(call)
    if c not in wdist:
        return True
    folded = {str(t) for t, k in (twins or {}).items()
              if str(k) == c and str(t) in pdist and c in pdist and pdist[c] - pdist[str(t)] < margin}
    rivals = [v for k, v in wdist.items() if k != c and k not in folded]
    return not rivals or min(rivals) - wdist[c] >= math.log(len(rivals) / alpha)


def edit_weights(w):
    """edit_weights of WORLD (edit_counts -> weights_from_counts), recorded in WORLD/graph/edit_weights.json."""
    c = edit_counts(w)
    wts = weights_from_counts(c)
    json.dump(dict(counts=dict(c), weights=wts), open(f"{w}/graph/edit_weights.json", "w"), indent=1)
    return wts


_last_aligned = 0                                          # reads the last cached_minimap2 call aligned (tests)


def cached_minimap2(cache_dir, args, target, reads, work):
    """minimap2 `args` of `reads` {name: sequence} against the FASTA `target`, through a cache shared by the worlds of a family
    set (O2_ALN_CACHE): keyed by the target's bytes and the arguments, one entry per read SEQUENCE (the worlds share their reads
    under other names). Only the sequences the cache lacks are aligned, named by their hash; each line is returned under every
    name holding that sequence. Exact for O2's distances: minimap2 aligns each read on its own (the set of lines is the same),
    and the distances read NM, the CIGAR and the coordinates; minimap2 breaks equal-score ties by a hash of the read's name,
    so a tied line's tp / s1 / s2 tags may differ."""
    import fcntl
    import hashlib
    import pickle
    import subprocess
    global _last_aligned
    os.makedirs(cache_dir, exist_ok=True)
    os.makedirs(work, exist_ok=True)
    h = hashlib.sha1(open(target, "rb").read())
    h.update("\0".join(args).encode())
    path = f"{cache_dir}/{h.hexdigest()}.pkl"
    sha = {n: hashlib.sha1(q.encode()).hexdigest() for n, q in reads.items()}
    entries = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
    need = {k: reads[n] for n, k in sha.items() if k not in entries}
    _last_aligned = len(need)
    if need:
        fa = f"{work}/aln_cache.{os.getpid()}.fa"
        with open(fa, "w") as fh:
            fh.writelines(f">{k}\n{q}\n" for k, q in sorted(need.items()))
        out = subprocess.run(args + ["-t", os.environ.get("O2_THREADS", "1"), target, fa], capture_output=True, text=True,
                             check=True).stdout
        os.remove(fa)
        new = {k: [] for k in need}
        for ln in out.splitlines():
            k, rest = ln.split("\t", 1)
            new[k].append(rest)
        with open(f"{path}.lock", "w") as lk:            # concurrent worlds merge their entries, never overwrite each other's
            fcntl.flock(lk, fcntl.LOCK_EX)
            cur = pickle.load(open(path, "rb")) if os.path.exists(path) else {}
            cur.update(new)
            tmp = f"{path}.{os.getpid()}"
            pickle.dump(cur, open(tmp, "wb"))
            os.replace(tmp, path)
        entries = {**entries, **new}
    return "".join(f"{n}\t{rest}\n" for n in reads for rest in entries[sha[n]])


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


def projection_spec(exons, micro=0):
    """What a projection query carries about its source isoform: its exon count, or its exon lengths for realworld's micro-exon
    rule (micro > 0: the projection may lack a terminal exon shorter than micro bases)."""
    return [b - a for a, b in exons] if micro else len(exons)


def union_structures(cps, projector, tol=None):
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
        out[c["idx"]] = dict(own, **R.cross_structures(own, proj, max_proj=None,
                                                       **({} if tol is None else {"tol": tol})))
    return out


_PATHS_FAMILY = None


def _paths_family(fam):
    return _PATHS_FAMILY(fam)


def run_paths(w, margin=MARGIN, reps_only=False, out_name="o2.tsv", union=False, ends=False, pad=0, catalog=None, hosted=None,
              distance="clip", proj_micro=0, proj_tol=None, max_div=DELTA, weights=None):
    """O2 on the family graph's ISOFORM PATHS: each hosted read aligned (minimap2 -x map-hifi, all secondaries) to every
    curated isoform of every catalogued copy as its own sequence -- the paths `vg rna` threads, without free graph alignment."""
    import subprocess
    import pysam
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    handles = {}

    def fa():                                        # one FASTA handle per process: forked workers share no file offset
        if os.getpid() not in handles:
            handles[os.getpid()] = pysam.FastaFile(f"{w}/genome.ref.fa")
        return handles[os.getpid()]
    iso = catalog if catalog is not None else json.load(open(f"{w}/isoforms.json"))
    reads = {}
    name = None
    with open(f"{w}/reads.E0.fq") as fh:
        for i, ln in enumerate(fh):
            if i % 4 == 0:
                name = ln[1:].strip().split()[0]
            elif i % 4 == 1:
                reads[name] = ln.strip()
    if hosted is None:
        hosted = collections.defaultdict(set)
        for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
            hosted[r["family_id"]].add(r["read_name"])
    d = (f"{w}/graph/o2paths" + (".reps" if reps_only else ".union" if union else "") + (".ends" if ends else "")
         + (".pad" if pad else "") + (".integ" if catalog is not None else ""))
    os.makedirs(d, exist_ok=True)
    R = _realworld() if union else None

    def projector(c, queries):                       # union: paralog isoforms spliced-aligned onto the copy's locus
        lo = max(0, c["start"] - R.PROJ_PAD)
        locus = fa().fetch(c["contig"], lo, c["end"] + R.PROJ_PAD).upper()
        q = {}
        for name, (p, ex) in queries.items():
            seq = "".join(fa().fetch(p["contig"], a, b).upper() for a, b in ex)
            q[name] = (seq[::-1].translate(comp) if p["strand"] == "-" else seq, projection_spec(ex, proj_micro))
        return R.project_many(q, locus, lo, c["strand"], micro=proj_micro)
    def one(fam):                                    # families are independent: O2_PROCS of them at a time, in family order
        names = sorted(n for n in hosted.get(fam, ()) if n in reads)
        rows, cnt = [], collections.Counter()
        if not names:
            return rows, cnt
        import time
        t0 = time.time()
        paths = (union_structures(iso[fam], projector, proj_tol) if union else
                 {c["idx"]: {t: ex for t, ex in c["isoforms"].items() if not reps_only or t == c["rep"]} for c in iso[fam]})
        t1 = time.time()
        json.dump({str(k): {t: [list(x) for x in ex] for t, ex in v.items()} for k, v in paths.items()},
                  open(f"{d}/{fam}.paths.json", "w"))                  # the structures each copy's paths follow
        n_tx, tx = 0, {}
        with open(f"{d}/{fam}.tx.fa", "w") as fh:
            for c in iso[fam]:
                for t, ex in paths[c["idx"]].items():
                    seq = "".join(fa().fetch(c["contig"], a, b).upper() for a, b in padded(ex, pad))
                    seq = seq[::-1].translate(comp) if c["strand"] == "-" else seq
                    tx[f"copy{c['idx']}|{t}"] = seq
                    fh.write(f">copy{c['idx']}|{t}\n{seq}\n")
                    n_tx += 1
                    cnt["projected" if t.startswith("proj:") else "own"] += 1
        with open(f"{d}/{fam}.reads.fa", "w") as fh:
            for n in names:
                fh.write(f">{n}\n{reads[n]}\n")
        t2 = time.time()
        mm2 = ["minimap2", "-c", "-x", "map-hifi", "-N", str(max(200, 2 * n_tx)), "-p", "0.3", "--secondary=yes"]
        if os.environ.get("O2_ALN_CACHE"):            # the family set's worlds share their reads: align each sequence once
            res = cached_minimap2(os.environ["O2_ALN_CACHE"], mm2, f"{d}/{fam}.tx.fa", {n: reads[n] for n in names}, d)
        else:
            res = subprocess.run(mm2 + ["-t", os.environ.get("O2_THREADS", "1"),  # output does not depend on the thread count
                                        f"{d}/{fam}.tx.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True,
                                 check=True).stdout
        t3 = time.time()
        dist = (_paf_dist_gc(res, reads, tx, struct=distance == "gapstruct") if distance in ("gapcomp", "gapstruct")
                else _paf_dist(res, ends=ends))
        if os.environ.get("O2_TIMING"):
            print(f"[o2graph] timing {fam}: {len(names)} reads, {n_tx} paths, {res.count(chr(10))} alignments; union paths "
                  f"{t1 - t0:.1f} s, write {t2 - t1:.1f} s, minimap2 {t3 - t2:.1f} s, distances {time.time() - t3:.1f} s", flush=True)
        json.dump(path_dists(dist), open(f"{d}/{fam}.dist.json", "w"))     # run_locus_fallback mode "all" combines them
        if weights is not None:                        # --o2-significance: the veto step's gate reads them
            json.dump(path_dists(_paf_dist_w(res, reads, tx, weights)), open(f"{d}/{fam}.wdist.json", "w"))
        tw = {f"copy{c['idx']}": f"copy{c['twin_of']}" for c in iso[fam] if c.get("twin_of") is not None}
        for n in names:
            c, reason = decide(dist.get(n, {}), len(reads[n]), margin, max_div=max_div, twins=tw)
            rows.append(dict(read=n, family=fam, call=c[4:] if c else "-", reason=reason))
        return rows, cnt
    global _PATHS_FAMILY
    _PATHS_FAMILY = one
    fams, procs = sorted(iso), int(os.environ.get("O2_PROCS", "1"))
    if procs > 1 and len(fams) > 1:
        import multiprocessing
        with multiprocessing.get_context("fork").Pool(procs) as pool:
            parts = pool.map(_paths_family, fams, chunksize=1)
    else:
        parts = [one(f) for f in fams]
    out, n_struct = [], collections.Counter()
    for rows, cnt in parts:
        out += rows
        n_struct.update(cnt)
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


def fallback_call(union_call, union_reason, locus_d, qlen, other_explained, margin=MARGIN):
    """Graph O2 with a structure-free fallback (round 5): union O2's decision stands unless union cannot explain the read in this
    family. Then the read is another family's (`other_family`) when that family's union explains it (within delta), else it is
    decided on its genomic-locus distances (`locus`: any canonical junction, so a structure no path threads is still placed)."""
    if union_reason != "unexplained":
        return union_call, union_reason
    if other_explained:
        return "-", "other_family"
    c, reason = decide(locus_d or {}, qlen, margin)
    return (c[4:], "locus") if c else ("-", reason)


def run_fallback(w, margin=MARGIN, union_name="o2.unionp.m2.tsv", out_name=None):
    """WORLD/graph/o2.unionf.m<margin>.tsv: padded union O2 with the locus fallback and the cross-family guard (fallback_call),
    from the world's union table and graph/locus.dist.json."""
    rows = list(csv.DictReader(open(f"{w}/graph/{union_name}"), delimiter="\t"))
    explained_in = collections.defaultdict(set)
    for r in rows:
        if r["reason"] != "unexplained":
            explained_in[r["read"]].add(r["family"])
    lj = json.load(open(f"{w}/graph/locus.dist.json"))
    qlen = {n: len(s) for n, s in _read_fq(f"{w}/reads.E0.fq").items()}
    out = []
    for r in rows:
        c, reason = fallback_call(r["call"], r["reason"], lj.get(r["family"], {}).get(r["read"]), qlen[r["read"]],
                                  bool(explained_in[r["read"]] - {r["family"]}), margin)
        out.append(dict(read=r["read"], family=r["family"], call=c, reason=reason))
    out_name = out_name or f"o2.unionf.m{margin}.tsv"
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w} (padded union + locus fallback + cross-family guard, margin {margin}) -> {out_name}: {len(out)} reads, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return 0


def flag_reads(clusters, members, min_reads=FLAG_MIN_READS):
    """Round 6: {(family, read)} in a PSV cluster (copy_assign --psv-clusters) whose copy-private alleles are linked
    (`unexplained`) with >= min_reads reads, WHATEVER its divergence -- below delta that is an unannotated copy or an allele."""
    out = set()
    for c in clusters:
        if c["status"] == "unexplained" and int(c["n_reads"]) >= min_reads:
            out |= {(c["family"], q) for q in members.get((c["family"], c["host"], c["cluster"]), ())}
    return out


def apply_flag(rows, flagged):
    """O2 + flag: a read in a flagged cluster is reported as an unannotated-copy candidate (`flagged`), not assigned."""
    return [dict(r, call="-", reason="flagged") if (r["family"], r["read"]) in flagged and r["call"] != "-" else dict(r)
            for r in rows]


def run_flagged(w, base, margin=MARGIN):
    """WORLD/graph/o2.flag_<base>.m<margin>.tsv: o2.<base>.m<margin>.tsv with the reads of flagged clusters left as candidates."""
    members = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.psv_cluster_reads.tsv"), delimiter="\t"):
        members[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
    flagged = flag_reads(list(csv.DictReader(open(f"{w}/cur/A.psv_clusters.tsv"), delimiter="\t")), members)
    rows = apply_flag(list(csv.DictReader(open(f"{w}/graph/o2.{base}.m{margin}.tsv"), delimiter="\t")), flagged)
    out_name = f"o2.flag_{base}.m{margin}.tsv"
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    c = collections.Counter(r["reason"] for r in rows)
    print(f"[o2graph] {w} ({base} + flag, margin {margin}) -> {out_name}: {len(rows)} reads, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return 0


def map_new_copies(new_copies, truth_copies):
    """Scoring only (round 7): {new seeded copy idx: the truth copy whose locus it overlaps most (copies.tsv rows), or None}."""
    out = {}
    for c in new_copies:
        best, ov = None, 0
        for t in truth_copies:
            if t["chrom"] != c["contig"]:
                continue
            o = min(c["end"], int(t["end"])) - max(c["start"], int(t["start"]))
            if o > ov:
                best, ov = t["copy"], o
        out[c["idx"]] = best
    return out


def score_assign(calls, idx2copy, truth, unknowable=()):
    """calls {(read, family): copy idx}; idx2copy {(family, idx): truth copy or None}. A call is correct when it names the read's
    own copy, wrong otherwise (another copy, a seeded locus no copy explains, or a read of no catalogued copy). outside_own: the
    reads of an uncatalogued in-genome copy (kind `outside`) called to that very copy (only an O1-seeded catalog can).
    unknowable: absent copies the genome explains within delta (undetectable by O3's definition); a placement of their reads
    cannot be told from that locus's own reads, so it is counted apart (`unknowable_placed`), not as an error (round 7)."""
    out = collections.Counter(correct=0, wrong=0, outside_own=0)
    if unknowable:
        out["unknowable_placed"] = 0
    for (q, f), i in calls.items():
        t = truth.get(q)
        if t is None:
            continue
        own = idx2copy.get((f, i)) == t["copy"] and t["family"] == f
        if not own and t["copy"] in unknowable:
            out["unknowable_placed"] += 1
            continue
        out["correct" if own else "wrong"] += 1
        out["outside_own"] += own and t["kind"] == "outside"
    return dict(out)


def run_integrated(w, margin=MARGIN, distance="clip", proj_micro=0, proj_tol=None, floor=0.0):
    """Round 7: padded union O2 on O1's seeded catalog (o1seed.py: graph/seeded.json) over the family reads it found ->
    graph/o2.integ.m<margin>.tsv. distance: "clip" (NM + clipped bases, round 7) or "gapcomp" (gap_compressed). proj_micro,
    proj_tol (--union-exact: 30, 0; round 10 development): a paralog isoform projects onto a copy though the aligner drops a
    terminal exon shorter than proj_micro bases, and a projected structure stays apart from the copy's own unless its
    junctions are within proj_tol bp (default realworld.JUNCTION_TOL)."""
    sj = json.load(open(f"{w}/graph/seeded.json"))
    catalog = {f: v["copies"] for f, v in sj.items() if v["copies"]}
    hosted = {f: set(v["family_reads"]) for f, v in sj.items()}
    return run_paths(w, margin, union=True, pad=PATH_PAD, out_name=f"o2.integ.m{margin}.tsv", catalog=catalog, hosted=hosted,
                     distance=distance, proj_micro=proj_micro, proj_tol=proj_tol, max_div=gate(floor),
                     weights=edit_weights(w) if _SIGNIFICANCE else None)


def locus_distance(cigartuples):
    """Gap-compressed edits of one genome alignment (an --eqx CIGAR): mismatches + insertion and deletion EVENTS (as minimap2's
    `de`); None when an insertion or a soft-clipped end of >= NOVEL_INS bases leaves read sequence the locus does not hold
    (a novel exon, an unannotated copy's own sequence: O3's evidence, never an O2 placement)."""
    d = 0
    for op, n in cigartuples:
        if op == 8:
            d += n
        elif op in (1, 2):
            if op == 1 and n >= NOVEL_INS:
                return None
            d += 1
        elif op == 4 and n >= NOVEL_INS:
            return None
    return d


def veto_distance(cigartuples):
    """The genome veto's distance of one alignment: mismatches, indel events, and every base the locus does not hold counted in
    full -- soft-clipped ends and insertions >= NOVEL_INS. locus_distance leaves a clip under NOVEL_INS free, so a paralog's
    alignment that clips a read's differing end beat the read's own locus (round 12 development: KRTAP2.1, 27S... at 2 vs 5)."""
    d = 0
    for op, n in cigartuples:
        if op == 8:
            d += n
        elif op == 1:
            d += n if n >= NOVEL_INS else 1
        elif op == 2:
            d += 1
        elif op == 4:
            d += n
    return d


def pairwise_cost(h, rival):
    """A copy's genome cost against one rival (--locus-pairwise): its alignment's edits plus the read bases the rival's alignment
    explains and its own leaves out. h, rival: (edits, read start, read end) as genome_hits gives them."""
    e, s, t = h
    _, rs, rt = rival
    return e + max(0, min(rt, s) - rs) + max(0, rt - max(rs, t))


def locus_decide(dist, qlen, margin=MARGIN, twins=None, max_div=DELTA, vdist=None, hits=None, gate=None):
    """The locus fallback's call from {seeded copy idx or ELSEWHERE: best locus_distance}: decide() over the family's seeded copies
    and the read's best alignment off them. -> (copy idx or None, `locus` | `locus_tied` | `locus_thin` | `unexplained` |
    `locus_clipped`); a read whose best genome alignment lies off the seeded copies (an unseeded or dropped copy, another family)
    stays `unexplained`.
    vdist (--locus-clipped, 2026-10-17, NOT ADOPTED): the same alignments' veto_distance. A copy's genome cost lies between its
    locus_distance (a clip under NOVEL_INS free) and its veto_distance (every clipped base paid); a call made on the genome
    alone stands only when its worst case beats every rival's best case by the margin, else `locus_clipped` (round 15 A119b
    502 CROSS: FAM90A.45's alignment clipped the 11 bases of FAM90A.51's read that told the copies apart, 3 edits on the rest
    against 5 for nine siblings). A twin folded into the call (decide's rule) is no rival. Development: it abstained correct
    calls where every copy's alignment clips the read's ends (HNRNPCL, TAF11L, FAM236) and failed r14 A119b s205's B1.
    hits (--locus-pairwise, 2026-10-17): genome_hits of the read. The call stands only when every rival's pairwise_cost against
    it exceeds the call's against that rival by the margin: a read base counts against a copy only where the rival explains
    it, so bases every copy clips are nobody's evidence, and FAM90A.45's 11 clipped bases still are."""
    c, reason = decide(dist, qlen, margin, max_div=max_div, twins=twins, gate=gate)
    if reason == "unexplained" or c == ELSEWHERE:
        return None, "unexplained"
    folded = {t for t, k in (twins or {}).items() if k == c and t in dist and dist[c] - dist[t] < margin}
    if c is not None and vdist is not None and c in vdist:
        rivals = [g for i, g in dist.items() if i != c and i not in folded]
        if rivals and vdist[c] + margin > min(rivals):
            return None, "locus_clipped"
    if c is not None and hits is not None and c in hits:
        if any(pairwise_cost(h, hits[c]) - pairwise_cost(hits[c], h) < margin for i, h in hits.items()
               if i != c and i not in folded):
            return None, "locus_clipped"
    return (c, "locus") if c is not None else (None, f"locus_{reason}")


VETO_MAPQ = 30                                           # the genome veto's confident placement (Phred 30: 1 in 1,000)


_VETO_ABSTAIN = False                                      # --veto-abstain (modes veto, guard): veto_decide
_LOCUS_CLIPPED = False                                     # --locus-clipped (modes veto, guard; not adopted): locus_decide
_LOCUS_PAIRWISE = False                                    # --locus-pairwise (modes veto, guard): locus_decide's hits
_STRUCTURE_CHECK = False                                   # --structure-check (modes veto, guard): structure_only
_LENGTH_GATE = False                                       # --length-gate (modes veto, guard): rescue, locus_decide
_SIGNIFICANCE_MARGIN = False                               # --significance-margin (modes veto, guard): rescue
_SHORT_PATH_ALIGN = False                                  # --short-path-align (modes veto, guard): short_path_distances
_STRUCT_RESCUE = False                                     # --struct-rescue (modes veto, guard): struct_distances
_TWIN_GUARD = False                                        # --twin-guard (modes veto, guard): tie_guard off the catalogue
_SIGNIFICANCE = False                                      # --o2-significance (integ: wdist.json; veto, guard: the gate)


def veto_decide(call, gdist, qlen, margin=MARGIN, twins=None, max_div=DELTA, primary=None, contradict_abstain=None,
                ldist=None):
    """mode veto (round 12 review C1): a path call stands unless the read's own genome alignments put it elsewhere by the
    margin -- another seeded copy, or a locus off them (ELSEWHERE), at least `margin` closer than the called copy's locus (no
    alignment there: infinitely far) and within max_div -- AND the genome places the read there confidently: its primary
    alignment (primary = (seeded copy idx or ELSEWHERE, MAPQ)) lies on the genome's winner at MAPQ >= VETO_MAPQ. A retained
    intron a paralog's path absorbs as one gap-compressed insertion is contiguous on the own locus (CDRT15: MAPQ 53-60); a short
    terminal exon the aligner misses at the true locus is not evidence (DDX11L, development: MAPQ 1). A genome tie keeps the
    path call. -> (call, "veto_locus" | "veto_unexplained") or None.
    contradict_abstain (--veto-abstain, 2026-10-17, round 15's FAM90A errors): when the genome contradicts the called copy but
    no confident single winner takes the read (the better copies tie among themselves, or the primary is not confident), the
    path call is abstained (`genome_contradicts`): the paths decided where the sequence contradicts them (FAM90A.45: 16 edits on
    the genome against 5 for four tied siblings, the siblings' paths lacking the read's exon skip, 64-378 away). A copy's
    genome cost lies between its locus_distance (ldist: a clip under NOVEL_INS free) and its veto_distance (gdist: every
    clipped base paid); the genome contradicts the call when a rival's worst case beats the call's best case by the margin. A
    clip alone is no contradiction: DDX11L.6's alignment clips a 15-bp terminal exon the aligner missed (1 edit on its
    aligned bases, 16 with the clip), DAZ.2's a read end. With no alignment on the called copy it is no evidence (minimap2's
    secondary cap); ldist None stands for gdist."""
    contradict_abstain = _VETO_ABSTAIN if contradict_abstain is None else contradict_abstain
    own = gdist.get(int(call), 1 << 30)
    rivals = [g for i, g in gdist.items() if i != int(call) and g + margin <= own]
    if not rivals or min(rivals) > max_div * qlen:
        return None
    c, why = decide(gdist, qlen, margin, max_div=max_div, twins=twins)
    if c is None or not primary or primary[0] != c or primary[1] < VETO_MAPQ:
        lo = (gdist if ldist is None else ldist).get(int(call)) if int(call) in gdist else None
        return ("-", "genome_contradicts") if contradict_abstain and lo is not None and min(rivals) + margin <= lo else None
    return ("-", "veto_unexplained") if c == ELSEWHERE else (str(c), "veto_locus")


def structure_only(call, pdist, cvdist, cldist, twins=None, margin=MARGIN):
    """--structure-check (2026-10-17): is a path call's margin over some rival only the paths' structure? A rival whose paths
    explain the read worse, by the margin, than its own locus does (an alignment inside its padded locus: cvdist its
    veto_distance, cldist its locus_distance) lacks the read's structure among its paths -- union projection adds a paralog's
    isoform only for an intron chain new to the copy, so a longer first exon that absorbs a retained intron, or a read-derived
    structure, can stay on one copy -- and the call is compared with that rival on the genome: it stands only when its best full
    explanation (its paths or its own locus) beats the rival's best case by the margin. r12 A119b KRTAP2.0's read retains an
    intron: copy 1's long first exon absorbs it (2), copy 0's paths clip it (85), both loci hold it exactly (0). r15 public
    503 FAM90A.12's 201-bp read: only copy 11's paths align it, eight loci hold it exactly. pdist: {copy: path distance};
    a twin folded into the call is no rival."""
    c, inf = int(call), float("inf")
    pd = {int(k): v for k, v in pdist.items()}
    own = min(pd.get(c, inf), cvdist.get(c, inf))
    folded = {t for t, k in (twins or {}).items() if k == c and t in pd and c in pd and pd[c] - pd[t] < margin}
    return any(pd.get(r, inf) - g >= margin and own + margin > cldist[r] for r, g in cvdist.items()
               if r not in (c, ELSEWHERE) and r not in folded and r in cldist)


def tie_guard(call, gdist):
    """mode guard (2026-10-17, the O2 proof's end-to-end OUTSIDE worlds: O1 did not seed C4's / HSFY's uncatalogued
    near-identical twin, and path O2 then placed the twin's reads on the catalogued copy): a path call is abstained when a
    locus off every seeded copy (ELSEWHERE) explains the read's genome alignments at least as well as the called copy's own
    locus. The catalog then lacks a copy the read may come from. With no alignment on the called copy it is the veto's case
    (a missed short exon is no evidence). --twin-guard (2026-10-18, round 17: every cell failed OUTSIDE on a pair's unseeded
    twin) applies it in mode veto with a catalogued copy O1 dropped as silent no longer off the catalogue (DROPPED): mode
    guard was rejected (register 1315) because it abstained reads of an expressed copy tied with such a silent copy."""
    own = gdist.get(int(call))
    return own is not None and ELSEWHERE in gdist and gdist[ELSEWHERE] <= own


def path_dists(dist):
    """run_paths' {read: {"copyN": distance}} -> {read: {"N": distance}}, as graph/o2paths*/<family>.dist.json holds it."""
    return {q: {c[4:]: v for c, v in cs.items()} for q, cs in dist.items()}


def genome_distances(w, todo, spans, dist_fn=None, contained=None):
    """{(read, family): {seeded copy idx or ELSEWHERE: best locus_distance}} over every primary and secondary alignment in
    WORLD/E0.bam of the reads in todo ({read: families}), and {read: its length from the alignments}. dist_fn: another
    per-alignment distance (veto_distance). contained (bases): an alignment counts for a copy only inside its locus padded by
    that much -- one that overlaps copies but lies inside none counts for none (--structure-check)."""
    dist_fn = dist_fn or locus_distance
    import pysam
    dist, qlen = collections.defaultdict(dict), {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in todo or rec.is_unmapped or rec.is_supplementary:
                continue
            qlen[q] = max(qlen.get(q, 0), rec.infer_read_length() or 0)
            d = dist_fn(rec.cigartuples)
            if d is None:
                continue
            for f in todo[q]:
                on = [i for ch, s, e, i in spans.get(f, ()) if ch == rec.reference_name and rec.reference_start < e
                      and rec.reference_end > s]
                if contained is not None:
                    on = [i for ch, s, e, i in spans.get(f, ()) if ch == rec.reference_name
                          and s - contained <= rec.reference_start and rec.reference_end <= e + contained] or (
                        [] if on else [ELSEWHERE])
                for i in on or ([] if contained is not None else [ELSEWHERE]):
                    if d < dist[(q, f)].get(i, 1 << 30):
                        dist[(q, f)][i] = d
    return dist, qlen


def genome_hits(w, todo, spans):
    """{(read, family): {seeded copy idx or ELSEWHERE: (edits, read start, read end)}}: of the primary and secondary alignments
    in WORLD/E0.bam of the reads in todo, each copy's of least veto_distance (then the widest) -- its edits (veto_distance less
    its soft-clipped bases) and the read interval it aligns, in the read's own orientation (--locus-pairwise)."""
    import pysam
    best = collections.defaultdict(dict)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in todo or rec.is_unmapped or rec.is_supplementary:
                continue
            cig = rec.cigartuples
            vd, lead = veto_distance(cig), 0
            for op, n in cig:
                if op not in (4, 5):
                    break
                lead += n
            qs = lead
            qe = qs + sum(n for op, n in cig if op in (0, 1, 7, 8))
            if rec.is_reverse:
                ln = rec.infer_read_length() or 0
                qs, qe = ln - qe, ln - qs
            key, hit = (vd, qs - qe), (vd - sum(n for op, n in cig if op == 4), qs, qe)
            for f in todo[q]:
                on = [i for ch, s, e, i in spans.get(f, ()) if ch == rec.reference_name and rec.reference_start < e
                      and rec.reference_end > s]
                for i in on or [ELSEWHERE]:
                    if i not in best[(q, f)] or key < best[(q, f)][i][0]:
                        best[(q, f)][i] = (key, hit)
    return {k: {i: h for i, (_, h) in v.items()} for k, v in best.items()}


def primary_loci(w, todo, spans):
    """{(read, family): (the seeded copy idx its primary alignment overlaps most, or ELSEWHERE; MAPQ)} for the reads in todo."""
    import pysam
    out = {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in todo or rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            for f in todo[q]:
                ov = [(min(rec.reference_end, e) - max(rec.reference_start, s), -i) for ch, s, e, i in spans.get(f, ())
                      if ch == rec.reference_name and rec.reference_start < e and rec.reference_end > s]
                out[(q, f)] = (-max(ov)[1] if ov else ELSEWHERE, rec.mapping_quality)
    return out


def run_locus_fallback(w, margin=MARGIN, out_name=None, mode="unexplained", floor=0.0):
    """Integrated O2 (o2.integ.m<margin>.tsv) with the read's own genome alignments (genome_distances: every primary and secondary
    alignment in E0.bam scored by locus_distance, the best per seeded copy of the family and the best off them, ELSEWHERE).
    mode "unexplained" -> WORLD/graph/o2.integl.m<margin>.tsv: only the reads O2 leaves `unexplained` in a family and explains in
    no other (else `other_family`), placed by locus_decide. mode "all" -> o2.integc.m<margin>.tsv: every hosted read decided on
    each copy's min(path distance, genome distance) and ELSEWHERE (graph/o2paths.union.pad.integ/<family>.dist.json), so a copy
    whose paths lack the read's structure is judged on its sequence; a call that differs from O2's is labelled genome_<reason>.
    mode "seeded" -> o2.integs.m<margin>.tsv: the reads O2 explains decided on min(path, genome) over the seeded copies alone
    (no ELSEWHERE: an array's uncatalogued identical units would tie every read), the unexplained ones as in "unexplained".
    mode "resolve" -> o2.integr.m<margin>.tsv: O2's assignments stand; its unexplained reads as in "unexplained", and its
    tied / thin ones decided on min(path, genome) with ELSEWHERE (`genome_assigned` when a seeded copy wins by the margin).
    mode "veto" -> o2.integv.m<margin>.tsv: as "unexplained", and an assigned read's path call is rechecked on its genome
    alignments (veto_decide); with --veto-abstain a contradicted path call is abstained (--locus-clipped: a clipped locus
    call too; not adopted), with --o2-significance an insignificant one (significant).
    mode "guard" -> o2.integt.m<margin>.tsv: as "veto", and a call the veto keeps is abstained (`guard_unseeded`) when a locus
    off every seeded copy explains the read as well (tie_guard). O3 still reads o2.integ.m<margin>.tsv, so no mode changes O3's evidence; an O3 call on a read takes precedence over its
    placement (summary.read_fate)."""
    rows = list(csv.DictReader(open(f"{w}/graph/o2.integ.m{margin}.tsv"), delimiter="\t"))
    by_read = collections.defaultdict(list)
    for r in rows:
        by_read[r["read"]].append(r)
    explained_in = collections.defaultdict(set)
    for r in rows:
        if r["reason"] != "unexplained":
            explained_in[r["read"]].add(r["family"])
    sj = json.load(open(f"{w}/graph/seeded.json"))
    spans = {f: [(c["contig"], c["start"], c["end"], c["idx"]) for c in v["copies"]] for f, v in sj.items()}
    twins = {f: {c["idx"]: c["twin_of"] for c in v["copies"] if c.get("twin_of") is not None} for f, v in sj.items()}
    todo = collections.defaultdict(set)                   # read -> the families decided here
    for r in rows:
        unexplained = r["reason"] == "unexplained"
        if (unexplained and not explained_in[r["read"]] - {r["family"]} or mode in ("all", "seeded") and not unexplained
                or mode == "resolve" and r["reason"] in ("tied", "thin")
                or mode in ("veto", "guard") and r["reason"] == "assigned"
                or mode in ("veto", "guard") and (_SIGNIFICANCE_MARGIN or _STRUCT_RESCUE) and r["reason"] in ("tied", "thin")):
            todo[r["read"]].add(r["family"])
    dist, qlen = genome_distances(w, todo, spans)
    resc = (_LENGTH_GATE or _SIGNIFICANCE_MARGIN or _SHORT_PATH_ALIGN or _STRUCT_RESCUE) and mode in ("veto", "guard")
    if mode in ("veto", "guard"):                         # the assigned reads' alignments, every unaligned base counted
        vtodo = todo if _LOCUS_CLIPPED or resc else {q: fs for q, fs in todo.items()
                                                     if any(r["reason"] == "assigned" for r in by_read[q])}
        vdist, _ = genome_distances(w, vtodo, spans, veto_distance)
        prim = primary_loci(w, vtodo, spans)
        hits = genome_hits(w, todo, spans) if _LOCUS_PAIRWISE else {}
        if _TWIN_GUARD:                                   # off the catalogue: neither a seeded copy nor one O1 dropped
            cat = json.load(open(f"{w}/isoforms.json")) if os.path.exists(f"{w}/isoforms.json") else {}
            tspans = {f: spans.get(f, []) + [(c["contig"], c["start"], c["end"], DROPPED) for c in cat.get(f, [])
                                             if c["idx"] in set(v.get("dropped", []))] for f, v in sj.items()}
            tdist, _ = genome_distances(w, vtodo, tspans, veto_distance)
        if _STRUCTURE_CHECK:                              # each copy's own locus: alignments inside its padded span
            cvd, _ = genome_distances(w, vtodo, spans, veto_distance, contained=PATH_PAD)
            cld, _ = genome_distances(w, vtodo, spans, contained=PATH_PAD)
    sig = _SIGNIFICANCE and mode in ("veto", "guard")
    if (sig or _STRUCTURE_CHECK or resc) and mode in ("veto", "guard"):   # the path distances, plain and weighted
        pg = {}
        for f in {r["family"] for r in rows}:
            p0 = f"{w}/graph/o2paths.union.pad.integ/{f}"
            pg[f] = tuple(json.load(open(f"{p0}.{k}.json")) if os.path.exists(f"{p0}.{k}.json") else {} for k in ("dist", "wdist"))
    if mode not in ("unexplained", "veto", "guard"):
        lens = {n: len(x) for n, x in _read_fq(f"{w}/reads.E0.fq").items()}
        pd = {}
        for f in {r["family"] for r in rows}:
            p = f"{w}/graph/o2paths.union.pad.integ/{f}.dist.json"
            pd[f] = json.load(open(p)) if os.path.exists(p) else {}
    realigned = collections.defaultdict(set)
    if resc:
        rseq = _read_fq(f"{w}/reads.E0.fq")
        rlen = {n: len(x) for n, x in rseq.items()}
    if resc and _SHORT_PATH_ALIGN:                        # the reads the path step left without any alignment
        wp = f"{w}/graph/edit_weights.json"
        wts = json.load(open(wp))["weights"] if os.path.exists(wp) and _SIGNIFICANCE else None
        for f in sorted({r["family"] for r in rows}):
            d0 = f"{w}/graph/o2paths.union.pad.integ"
            need = {r["read"] for r in rows if r["family"] == f and r["reason"] == "unexplained" and f in todo.get(r["read"], ())
                    and not pg[f][0].get(r["read"]) and r["read"] in rseq}
            if need and os.path.exists(f"{d0}/{f}.tx.fa"):
                sd, sw = short_path_distances(d0, f, {q: rseq[q] for q in need}, wts)
                for q in sd:
                    pg[f][0][q] = sd[q]
                    if q in sw:
                        pg[f][1][q] = sw[q]
                    realigned[f].add(q)
    sdist = {}
    if resc and _STRUCT_RESCUE:                           # the thin and tied reads, on structure-aware path distances
        d0 = f"{w}/graph/o2paths.union.pad.integ"
        for f in sorted({r["family"] for r in rows}):
            need = {r["read"] for r in rows if r["family"] == f and r["reason"] in ("thin", "tied") and f in todo.get(r["read"], ())
                    and r["read"] in rseq}
            if need and os.path.exists(f"{d0}/{f}.tx.fa"):
                sdist[f] = struct_distances(d0, f, {q: rseq[q] for q in need})
    out = []
    for r in rows:
        call, reason, q, f = r["call"], r["reason"], r["read"], r["family"]
        rescued = None
        if resc and f in todo.get(q, ()) and reason in ("unexplained", "thin", "tied"):
            rc = rescue(reason, pg[f][0].get(q, {}), pg[f][1].get(q, {}), qlen.get(q) or rlen.get(q, 0), floor, twins.get(f),
                        margin, length_gate=_LENGTH_GATE, significance_margin=_SIGNIFICANCE_MARGIN, realigned=q in realigned[f],
                        sdist=sdist.get(f, {}).get(q))
            if rc:
                call, reason, rescued = str(rc[0]), "assigned", rc[1]
        if f not in todo.get(q, ()):
            if reason == "unexplained":
                call, reason = "-", "other_family"
        elif mode == "all" or mode == "seeded" and reason != "unexplained":
            comb = {int(c): v for c, v in pd[f].get(q, {}).items()}
            for i, g in dist.get((q, f), {}).items():
                if mode == "all" or i != ELSEWHERE:
                    comb[i] = min(comb.get(i, 1 << 30), g)
            c, why = decide(comb, lens[q], margin, max_div=gate(floor), twins=twins.get(f))
            c, why = (None, "unexplained") if c == ELSEWHERE else (c, why)
            new = ("-" if c is None else str(c), why)
            call, reason = (call, reason) if new == (call, reason) else (new[0], f"genome_{why}")
        elif mode == "resolve" and reason != "unexplained":
            comb = {int(c): v for c, v in pd[f].get(q, {}).items()}
            for i, g in dist.get((q, f), {}).items():
                comb[i] = min(comb.get(i, 1 << 30), g)
            c, why = decide(comb, lens[q], margin, max_div=gate(floor), twins=twins.get(f))
            if c is not None and c != ELSEWHERE:
                call, reason = str(c), "genome_assigned"
        elif mode in ("veto", "guard") and reason == "assigned":
            v = veto_decide(call, vdist.get((q, f), {}), qlen.get(q, 0), margin, twins=twins.get(f), max_div=gate(floor),
                            primary=prim.get((q, f)), ldist=dist.get((q, f), {}))
            call, reason = v if v else (call, reason)
            if sig and reason == "assigned" and not significant(call, pg[f][0].get(q, {}), pg[f][1].get(q, {}), twins.get(f),
                                                                margin):
                call, reason = "-", "insignificant"
            if _STRUCTURE_CHECK and reason == "assigned" and structure_only(call, pg[f][0].get(q, {}), cvd.get((q, f), {}),
                                                                            cld.get((q, f), {}), twins.get(f), margin):
                call, reason = "-", "structure_only"
            if _TWIN_GUARD and reason == "assigned" and tie_guard(call, tdist.get((q, f), {})):
                call, reason = "-", "uncatalogued_tie"
            if mode == "guard" and reason == "assigned" and tie_guard(call, vdist.get((q, f), {})):
                call, reason = "-", "guard_unseeded"
        elif reason == "unexplained":
            vd = vdist.get((q, f), {}) if _LOCUS_CLIPPED and mode in ("veto", "guard") else None
            hs = hits.get((q, f), {}) if _LOCUS_PAIRWISE and mode in ("veto", "guard") else None
            c, reason = locus_decide(dist.get((q, f), {}), qlen.get(q, 0), margin, twins=twins.get(f), max_div=gate(floor),
                                     vdist=vd, hits=hs,
                                     gate=(lambda k, n: within_reach(k, n, floor)) if _LENGTH_GATE and resc else None)
            call = "-" if c is None else str(c)
        if rescued and reason == "assigned":
            reason = rescued
        out.append(dict(read=q, family=f, call=call, reason=reason))
    out_name = out_name or f"o2.integ{dict(all='c', seeded='s', unexplained='l', resolve='r', veto='v', guard='t')[mode]}.m{margin}.tsv"
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w} (integrated + genome alignments, mode {mode}, margin {margin}) -> {out_name}: {len(out)} rows, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    return 0


def as_scope(loci, call, reason, qlen, floor=0.0):
    """The user's split (2026-10-19: "go ahead and use only for tied AS reads"; 2026-10-18: "O2 can only work for tied AS
    reads, everything else is already resolved by minimap, then O3 is just an extension in which if we notice the PSVs nothing
    in the genome explains them"). loci: {seeded copy idx or ELSEWHERE: the least locus_distance of the read's alignments
    there at its best AS (None: the locus lacks the read's sequence)}; call / reason: the base O2 table's. -> (call, reason):
    - no alignment: `as_unaligned`;
    - a best-AS alignment off the family's seeded copies (alone, or tied with a copy: a locus O1 did not seed): `as_elsewhere`;
    - one best copy: minimap2 resolved the read -- that copy (`as_resolved`) when its alignment explains the read within
      gate(floor), else `as_unexplained` (nothing in the genome explains it: O3's);
    - an AS tie between copies: O2's call when it is one of them, else `as_tie_outside`; O2's abstention stands.
    A catalogued copy O1 dropped as silent (DROPPED, scope_spans) is a copy, not a locus off the catalogue: a read tied with it
    is O2's (register 1315: abstaining there lost an expressed copy's reads), and one placed best there alone is `as_dropped`."""
    if not loci:
        return "-", "as_unaligned"
    if ELSEWHERE in loci:
        return "-", "as_elsewhere"
    if set(loci) == {DROPPED}:
        return "-", "as_dropped"
    if len(loci) == 1:
        (i, d), = loci.items()
        return ("-", "as_unexplained") if d is None or d > gate(floor) * qlen else (str(i), "as_resolved")
    if call == "-":
        return "-", reason
    return (call, reason) if int(call) in loci else ("-", "as_tie_outside")


FILL = ("thin", "tied", "locus_thin", "locus_tied")      # O2 abstained for want of its margin, not for a conflict


def fewest_edit_call(call, reason, edits):
    """The single-genome rule (user 2026-10-19: "until we make O2 work we stay with only 1 genome, in this case CHM13"; the
    design approved the same day). edits: {seeded copy idx or DROPPED: the fewest edits of the read's tied alignments there,
    every unaligned base counted (veto_distance)}. -> (call, reason, set). O2's own call stands, its structural checks
    included. Where O2 abstained for want of its two-edit margin (FILL), the copy with the fewest edits takes the read when it
    is the only one: in one genome a single differing base places a HiFi read (CHM13-only development on round 22's
    families: 154 of 155). Otherwise the read keeps the set of its fewest-edit copies; the other placements explain it worse
    and are dropped. Any other abstention (a structural conflict, a call off the tied copies) keeps every tied copy, and a
    silent copy (DROPPED) is no seeded copy to place a read on."""
    if call != "-":
        return call, reason, [int(call)]
    if reason not in FILL:
        return "-", reason, sorted(edits)
    m = min(edits.values())
    best = sorted(c for c, e in edits.items() if e == m)
    if len(best) == 1 and best[0] != DROPPED:
        return str(best[0]), "fewest_edits", best
    return "-", "fewest_edit_set", best


def _idx(copies):
    return ",".join("dropped" if c == DROPPED else str(c) for c in copies)


def best_as_loci(w, todo, spans, slack=0, seq_tie=False, edits=False):
    """{(read, family): {seeded copy idx or ELSEWHERE: least locus_distance}} over the alignments in WORLD/E0.bam (primary and
    secondary) at the read's best AS -- within `slack` of it (--as-slack; 0: equal AS only) -- for the reads in todo ({read:
    families}); and {read: length}. Under splice:hq one substitution moves AS by 5 (+1 match, -4 mismatch), so a slack of 4
    ties alignments no single sequence difference separates. seq_tie (--seq-tie, user 2026-10-19 after round 21): an
    alignment that explains the read with no more edits than the best-AS one, counting every base it leaves unaligned
    (veto_distance), ties too, whatever its AS -- AS charges each intron and nothing for a soft clip, so it can favour the
    copy that explains less of the read (round 21's NCOR1P4 read: 278 with 38 bases clipped against its own copy's 262 with
    all 340 aligned, both NM 0). edits (--fewest-edits): a third value, {(read, family): {locus: the fewest edits of the
    read's tied alignments there, every unaligned base counted (veto_distance)}}."""
    import pysam
    alns, qlen = collections.defaultdict(list), {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in todo or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            qlen[q] = max(qlen.get(q, 0), rec.infer_read_length() or 0)
            alns[q].append((rec.get_tag("AS"), rec.reference_name, rec.reference_start, rec.reference_end,
                            locus_distance(rec.cigartuples), veto_distance(rec.cigartuples)))
    out, eout = {}, {}
    for q, al in alns.items():
        best = max(a[0] for a in al)
        ref = min(a[5] for a in al if a[0] == best)
        for f in todo[q]:
            loci, ed = {}, {}
            for sc, ch, s0, e0, d, vd in al:
                if sc < best - slack and not (seq_tie and vd <= ref):
                    continue
                for i in [i for c2, s, e, i in spans.get(f, ()) if c2 == ch and s0 < e and e0 > s] or [ELSEWHERE]:
                    old = loci.get(i, "unset")
                    loci[i] = d if old == "unset" or old is None else (old if d is None else min(old, d))
                    ed[i] = min(ed.get(i, vd), vd)
            out[(q, f)] = loci
            eout[(q, f)] = ed
    return (out, qlen, eout) if edits else (out, qlen)


def scope_spans(sj, cat):
    """{family: [(contig, start, end, seeded copy idx or DROPPED)]}: graph/seeded.json's copies, and the catalogued copies
    (isoforms.json) O1 dropped as silent."""
    out = {}
    for f, v in sj.items():
        drop = set(v.get("dropped", []))
        out[f] = ([(c["contig"], c["start"], c["end"], c["idx"]) for c in v["copies"]]
                  + [(c["contig"], c["start"], c["end"], DROPPED) for c in cat.get(f, []) if c["idx"] in drop])
    return out


def run_scope(w, margin=MARGIN, base="integv", floor=0.0, slack=0, seq_tie=False, fewest=False):
    """WORLD/graph/o2.<base>.m<margin>.tsv scoped by as_scope -> o2.intega.m<margin>.tsv: O2 only on the AS-tied reads
    (best_as_loci's slack). fewest (--fewest-edits): also o2.intege.m<margin>.tsv, the single-genome rule (fewest_edit_call)
    on the tied reads, with the o2.intega call beside it (base_call), the tied copies (candidates) and the fewest-edit set
    (set); o2.intega is unchanged."""
    rows = list(csv.DictReader(open(f"{w}/graph/o2.{base}.m{margin}.tsv"), delimiter="\t"))
    sj = json.load(open(f"{w}/graph/seeded.json"))
    cat = json.load(open(f"{w}/isoforms.json")) if os.path.exists(f"{w}/isoforms.json") else {}
    spans = scope_spans(sj, cat)
    todo = collections.defaultdict(set)
    for r in rows:
        todo[r["read"]].add(r["family"])
    loci, qlen, *ed = best_as_loci(w, todo, spans, slack, seq_tie, edits=fewest)
    out, ext = [], []
    for r in rows:
        q, f = r["read"], r["family"]
        lo = loci.get((q, f), {})
        call, reason = as_scope(lo, r["call"], r["reason"], qlen.get(q, 0), floor)
        out.append(dict(read=q, family=f, call=call, reason=reason))
        if fewest and len(lo) >= 2 and ELSEWHERE not in lo:
            c2, r2, st = fewest_edit_call(call, reason, ed[0][(q, f)])
            ext.append(dict(read=q, family=f, call=c2, reason=r2, base_call=call, candidates=_idx(sorted(lo)), set=_idx(st)))
        elif fewest:
            ext.append(dict(read=q, family=f, call=call, reason=reason, base_call=call, candidates="-", set="-"))
    out_name = f"o2.intega.m{margin}.tsv"
    with open(f"{w}/graph/{out_name}", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason"], delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(out)
    c = collections.Counter(r["reason"] for r in out)
    print(f"[o2graph] {w} (O2 on AS-tied reads only, base {base}, AS slack {slack}{', sequence ties' if seq_tie else ''}, margin {margin}) -> {out_name}: {len(out)} rows, "
          + ", ".join(f"{k} {v}" for k, v in sorted(c.items())), flush=True)
    if fewest:
        with open(f"{w}/graph/o2.intege.m{margin}.tsv", "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=["read", "family", "call", "reason", "base_call", "candidates", "set"],
                                delimiter="\t", lineterminator="\n")
            wr.writeheader()
            wr.writerows(ext)
        c = collections.Counter(r["reason"] for r in ext if r["candidates"] != "-")
        print(f"[o2graph] {w} (fewest edits on the tied reads) -> o2.intege.m{margin}.tsv: "
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
    global _TWIN_ABSTAIN, _VETO_ABSTAIN, _LOCUS_CLIPPED, _LOCUS_PAIRWISE, _STRUCTURE_CHECK, _SIGNIFICANCE, _LENGTH_GATE
    global _SIGNIFICANCE_MARGIN, _SHORT_PATH_ALIGN, _STRUCT_RESCUE, _TWIN_GUARD
    _TWIN_ABSTAIN = "--twin-abstain" in argv
    _VETO_ABSTAIN = "--veto-abstain" in argv
    _LOCUS_CLIPPED = "--locus-clipped" in argv
    _LOCUS_PAIRWISE = "--locus-pairwise" in argv
    _STRUCTURE_CHECK = "--structure-check" in argv
    _LENGTH_GATE = "--length-gate" in argv
    _SIGNIFICANCE_MARGIN = "--significance-margin" in argv
    _SHORT_PATH_ALIGN = "--short-path-align" in argv
    _STRUCT_RESCUE = "--struct-rescue" in argv
    _TWIN_GUARD = "--twin-guard" in argv
    _SIGNIFICANCE = "--o2-significance" in argv
    if argv[0] == "run":
        # run WORLD [paths|nodes] [MARGIN]: `nodes` = free GraphAligner alignment, kept for the record (development: it makes
        # artifactual allele calls through near-identical copies, GAGE 74 wrong of 289)
        method = argv[2] if len(argv) > 2 else "paths"
        margin = int(argv[3]) if len(argv) > 3 and argv[3].isdigit() else MARGIN
        if method == "nodes":
            return run_world(argv[1])
        if method == "locus":
            return run_locus(argv[1], margin)
        if method == "integ":                             # round 7: needs o1seed.py WORLD first
            return run_integrated(argv[1], margin, argv[argv.index("--distance") + 1] if "--distance" in argv else "clip",
                                  *((UNION_EXACT_MICRO, 0) if "--union-exact" in argv else (0, None)),
                                  floor=_detect().error_floor(f"{argv[1]}/E0.bam") if "--error-floor" in argv else 0.0)
        if method in ("integl", "integc", "integs", "integr", "integv", "integt"):   # + the genome alignments; needs o2.integ
            return run_locus_fallback(argv[1], margin, mode=dict(integl="unexplained", integc="all", integs="seeded",
                                                                 integr="resolve", integv="veto", integt="guard")[method],
                                      floor=_detect().error_floor(f"{argv[1]}/E0.bam") if "--error-floor" in argv else 0.0)
        if method == "intega":                            # O2 only on AS-tied reads: o2.<--base, integv>.m<margin>.tsv scoped
            return run_scope(argv[1], margin, argv[argv.index("--base") + 1] if "--base" in argv else "integv",
                             floor=_detect().error_floor(f"{argv[1]}/E0.bam") if "--error-floor" in argv else 0.0,
                             slack=int(argv[argv.index("--as-slack") + 1]) if "--as-slack" in argv else 0,
                             seq_tie="--seq-tie" in argv, fewest="--fewest-edits" in argv)
        if method == "flag":                              # run WORLD flag BASE [MARGIN]: o2.<BASE> + the round-6 flag
            return run_flagged(argv[1], argv[3], int(argv[4]) if len(argv) > 4 else MARGIN)
        if method == "unionf":                            # needs o2.unionp.m<margin>.tsv and graph/locus.dist.json first
            return run_fallback(argv[1], margin, union_name=f"o2.unionp.m{margin}.tsv")
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
