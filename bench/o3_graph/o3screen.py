#!/usr/bin/env python3
"""EXPLORATORY (spike, 2026-10-22, post hoc): an RNA-only O3 "look harder" screen, tiers 1 and 2 of
docs/O3_FINGERPRINT_SURVEY_2026-10-22.md (definitions: docs/o3survey_2026-10-22/C_certificate.md §1-3; results and every choice
made here: docs/O3_SCREEN_SPIKE_2026-10-22.md). Read-only over a frozen O2-VC run directory (o2vc.py run); reads no truth.

    o3screen.py OUTDIR COPIES GENOME OUT.tsv [--sex male|female|unknown] [--exclude-editing] [--pair-certified]
                [--merge-single] [--families A,B]

- Informative positions of copy c: every position of c in the family's columns.tsv (any kind, any tied read; with
  --exclude-editing the rows of kind "editing" are left out).
- A pattern: a read's bases at the informative positions it covers, base quality >= MIN_BASEQ (10). Teacher reads (OUTDIR/teach,
  each record attributed to its copy exactly as o2vc_data.scan does) are read from the BAM in the genome's strand; tied reads
  from columns.tsv, put in each tied copy's genome strand (flip found from the genome base). A read covering no position has no
  pattern.
- Two patterns CONFLICT when they differ at >= 2 positions both cover. Classes: greedy first fit over the distinct patterns
  (most positions first, then most reads, then the pattern), a pattern joining the first class none of whose member patterns it
  conflicts with. A class of >= 2 reads is SUPPORTED. A class's pattern is its certified consensus: the positions where >= 2
  members carry the same base and it is the strict majority (every fact observed twice). Two classes conflict when their
  consensuses do (v1, default); --pair-certified (v2) also needs two such positions carried together by >= 2 members of each
  class; --merge-single (v3) first merges classes whose consensuses differ at exactly one position, by single linkage
  (C_certificate §1), until none is left; --leak-closure (v4) flags tier 1 only on a Hall deficiency > 0 when each clique
  class may also take the slots of any other family copy having a teacher class it does not conflict with, compared in that
  copy's coordinates through the columns.tsv homology (C §1, S4); a copy no tied read links to c shares no position and
  counts as compatible (the conservative choice for a certificate) unless --leak-needs-overlap.
- Tier 1 (per copy): a greedy multi-start clique of pairwise-conflicting supported teacher classes (a lower bound on the clique
  number), flagged when larger than the copy's capacity: 2, or 1 on a non-PAR chrX/chrY copy when --sex male.
- Tier 2 (per family): the tied reads of calls.tsv, classed per anchor copy; a supported tied class is NOVEL when its consensus,
  put on each candidate copy of its reads (the union), conflicts with every supported teacher class of every candidate copy
  (literal: a copy with no supported class adds no constraint, vacuously). The "taught" variant also needs every candidate copy
  to have >= 1 supported teacher class. Classes found from several anchors are counted once (by their read set).
- Today's flag: OUTDIR/<family>.flags.tsv (o2vc_rule.flags), counted as is.
"""
import argparse
import collections
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUPPORT = 2                                                     # a fact (a class, a consensus base) is observed by >= 2 reads
CONFLICT = 2                                                    # two patterns conflict at >= 2 co-covered positions


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_o3screen", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _sibling("o2vc_rule")
MIN_BASEQ, COMP = R.MIN_BASEQ, R.COMP
PAR = {"chrX": ((0, 2394410), (153925834, 154259566)), "chrY": ((0, 2458320), (62122809, 62460029))}   # hg002truth.PAR


def conflict(a, b):
    """Whether patterns a and b ({pos: base}) differ at >= CONFLICT positions both cover."""
    if len(a) > len(b):
        a, b = b, a
    n = 0
    for p, x in a.items():
        y = b.get(p)
        if y is not None and y != x:
            n += 1
            if n >= CONFLICT:
                return True
    return False


def classes(patterns, merge=False):
    """Greedy first-fit classes of {read: {pos: base}} (empty patterns left out): [klass()] of the SUPPORTED classes (>= SUPPORT
    reads), in the order they were opened; merge=True then merges them by single linkage (merge_single)."""
    distinct = collections.defaultdict(list)
    for r, pat in patterns.items():
        if pat:
            distinct[tuple(sorted(pat.items()))].append(r)
    order = sorted(distinct, key=lambda k: (-len(k), -len(distinct[k]), k))
    opened = []                                                 # [[pattern keys]]
    for k in order:
        pk = dict(k)
        for cl in opened:
            if not any(conflict(pk, dict(m)) for m in cl):
                cl.append(k)
                break
        else:
            opened.append([k])
    out = []
    for cl in opened:
        reads = sorted(r for k in cl for r in distinct[k])
        if len(reads) < SUPPORT:
            continue
        out.append(klass(reads, patterns))
    return merge_single(out, patterns) if merge else out


def n_diff(a, b):
    """The number of co-covered positions where patterns a and b differ, counted up to CONFLICT."""
    if len(a) > len(b):
        a, b = b, a
    n = 0
    for p, x in a.items():
        y = b.get(p)
        if y is not None and y != x:
            n += 1
            if n >= CONFLICT:
                break
    return n


def merge_single(cls, patterns):
    """C_certificate §1: supported classes whose consensuses differ at exactly one co-covered position are merged by single
    linkage, repeated until no such pair is left (a class one column from another is not a separate haplotype)."""
    while True:
        n = len(cls)
        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        joined = False
        for i in range(n):
            for j in range(i + 1, n):
                if n_diff(cls[i][1], cls[j][1]) == 1 and find(i) != find(j):
                    parent[max(find(i), find(j))] = min(find(i), find(j))
                    joined = True
        if not joined:
            return cls
        groups = collections.defaultdict(list)
        for i in range(n):
            groups[find(i)].extend(cls[i][0])
        cls = [klass(sorted(groups[g]), patterns) for g in sorted(groups)]


def klass(reads, patterns):
    """(reads, certified consensus, carriers {pos: the members carrying the consensus base there}) of a class."""
    cons = consensus(patterns[r] for r in reads)
    car = collections.defaultdict(set)
    for r in reads:
        for p, b in patterns[r].items():
            if cons.get(p) == b:
                car[p].add(r)
    return reads, cons, car


def class_conflict(a, b, pair):
    """Whether classes a and b conflict: their consensuses differ at >= CONFLICT co-covered positions (pair=False), or
    (pair=True, C_certificate §1 S1: a pattern is certified by >= 2 reads carrying ALL of it) some two such positions are each
    carried together by >= SUPPORT members of a and by >= SUPPORT members of b."""
    if not pair:
        return conflict(a[1], b[1])
    d = sorted(p for p, x in a[1].items() if b[1].get(p, x) != x)
    for x in range(len(d)):
        for y in range(x + 1, len(d)):
            i, j = d[x], d[y]
            if len(a[2][i] & a[2][j]) >= SUPPORT and len(b[2][i] & b[2][j]) >= SUPPORT:
                return True
    return False


def consensus(pats):
    """The certified consensus of patterns: {pos: base} where >= SUPPORT carry the base and it is the strict majority."""
    count = collections.defaultdict(collections.Counter)
    for pat in pats:
        for p, b in pat.items():
            count[p][b] += 1
    out = {}
    for p, c in count.items():
        top = c.most_common(2)
        if top[0][1] >= SUPPORT and (len(top) == 1 or top[0][1] > top[1][1]):
            out[p] = top[0][0]
    return out


def greedy_clique(cls, pair=False):
    """The largest set of pairwise-conflicting classes found by greedy multi-start (start from each, add in the given order
    every class conflicting with all chosen): a lower bound on the clique number. Returns the chosen indices."""
    n = len(cls)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = adj[j][i] = class_conflict(cls[i], cls[j], pair)
    best = []
    for s in range(n):
        chosen = [s]
        for j in range(n):
            if j != s and all(adj[j][i] for i in chosen):
                chosen.append(j)
        if len(chosen) > len(best):
            best = chosen
    return sorted(best)


def deficiency(lists):
    """Hall deficiency of classes into slots: len(lists) minus a maximum matching, lists[i] being the slots class i may use
    (Kuhn's augmenting paths; the instances hold a few classes)."""
    owner = {}

    def augment(i, seen):
        for s in lists[i]:
            if s not in seen:
                seen.add(s)
                if s not in owner or augment(owner[s], seen):
                    owner[s] = i
                    return True
        return False

    return len(lists) - sum(1 for i in range(len(lists)) if augment(i, set()))


def capacity(chrom, start, end, sex):
    """A copy's allele slots: 1 on a non-PAR chrX/chrY copy of a male, else 2."""
    if sex != "male" or chrom not in PAR:
        return 2
    return 2 if any(start < e and end > s for s, e in PAR[chrom]) else 1


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def family_rows(out_dir, f, exclude_editing):
    """columns.tsv of family f: (informative {copy: {pos}}, rows [(read, base, qual, [(copy, chrom, pos, ref)])])."""
    inf, rows = collections.defaultdict(set), []
    p = os.path.join(out_dir, f"{f}.columns.tsv")
    if not os.path.exists(p):
        return inf, rows
    for r in _tsv(p):
        if exclude_editing and r["kind"] == "editing":
            continue
        cps = []
        for x in r["copies"].split(";"):
            c, chrom, pos, ref, _, _ = x.split("|")
            cps.append((c, chrom, int(pos), ref))
            inf[c].add(int(pos))
        rows.append((r["read"], r["base"], int(r["qual"]), cps))
    return inf, rows


def teacher_patterns(bam_p, fam_copies, inf, cps_all, idx):
    """{copy: {teacher read: {pos: base}}} from the family's teacher BAM, genome strand, base quality >= MIN_BASEQ."""
    import bisect
    import pysam
    P, S, C = _sibling("o2psv"), _sibling("summary"), _sibling("o2correct")
    want = {c: sorted(v) for c, v in inf.items()}
    out = collections.defaultdict(dict)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.query_sequence is None:
                continue
            cands = C.overlapping(idx, rec.reference_name, rec.reference_start, rec.reference_end)
            c = S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None
            if c not in fam_copies:
                continue
            pos = want.get(c, [])
            seq, qual = rec.query_sequence.upper(), rec.query_qualities
            pat, q, r = {}, 0, rec.reference_start
            for op, n in rec.cigartuples:
                if op in (0, 7, 8):
                    i = bisect.bisect_left(pos, r)
                    while i < len(pos) and pos[i] < r + n:
                        k = q + pos[i] - r
                        if (qual is None or qual[k] >= MIN_BASEQ) and seq[k] in "ACGT":
                            pat[pos[i]] = seq[k]
                        i += 1
                    q, r = q + n, r + n
                elif op in (1, 4):
                    q += n
                elif op in (2, 3):
                    r += n
            out[c][rec.query_name] = pat
    return out


def tied_patterns(rows, genome):
    """{read: {copy: {pos: base in the copy's genome strand}}} of the tied reads' columns at quality >= MIN_BASEQ; a
    (read, copy) whose strand cannot be told from the genome base is left out."""
    flip, out = {}, collections.defaultdict(lambda: collections.defaultdict(dict))
    for q, base, qual, cps in rows:
        for c, chrom, pos, ref in cps:
            if (q, c) not in flip:
                g = genome.fetch(chrom, pos, pos + 1).upper()
                if g == ref:
                    flip[q, c] = False
                elif g.translate(COMP) == ref:
                    flip[q, c] = True
                else:
                    continue
            if qual >= MIN_BASEQ and base in "ACGT":
                out[q][c][pos] = base.translate(COMP) if flip[q, c] else base
    return out, flip


def homologs(rows, flip):
    """{(c, d): {position on c: (position on d, complemented)}}: the columns.tsv rows put one read base on every tied copy, so a
    row maps c's position to d's, complemented when the read lies on the two copies in opposite orientations; the most frequent
    mapping wins."""
    seen = collections.defaultdict(collections.Counter)
    for q, _, _, cps in rows:
        for c, _, pc, _ in cps:
            for d, _, pd, _ in cps:
                if c != d and (q, c) in flip and (q, d) in flip:
                    seen[c, d, pc][pd, flip[q, c] != flip[q, d]] += 1
    out = collections.defaultdict(dict)
    for (c, d, pc), cnt in seen.items():
        out[c, d][pc] = max(cnt.items(), key=lambda kv: (kv[1], kv[0]))[0]
    return out


def translate(cons, hom):
    """A consensus {pos on c: base} put on d's positions and strand through hom = homologs()[c, d]; unmapped positions drop."""
    return {hom[p][0]: (b.translate(COMP) if hom[p][1] else b) for p, b in cons.items() if p in hom}


def screen_family(f, out_dir, cand, fam_copies, span, cps_all, idx, genome, sex, exclude_editing, pair=False, merge=False,
                  leak=False, leak_overlap=False):
    """One output row for family f (cand: {tied read: [candidate copies]} from calls.tsv)."""
    inf, rows = family_rows(out_dir, f, exclude_editing)
    teach = teacher_patterns(os.path.join(out_dir, "teach", f"{f}.bam"), set(fam_copies), inf, cps_all, idx)
    tclass = {c: classes(teach.get(c, {}), merge) for c in fam_copies}
    tp, flip = tied_patterns(rows, genome)
    hom = homologs(rows, flip) if leak else {}
    t1 = []
    for c in fam_copies:
        cl = tclass[c]
        q = greedy_clique(cl, pair) if cl else []
        cap = capacity(*span[c], sex)
        if len(q) <= cap:
            continue
        if not leak:
            t1.append(f"{c}:{len(q)}/{cap}")
            continue
        lists = []                                              # leak closure (C §1): c's slots, and those of every
        for i in q:                                             # other copy with a teacher class compatible with class i
            sl = [(c, x) for x in range(cap)]
            for d in fam_copies:
                moved = translate(cl[i][1], hom.get((c, d), {}))     # class i in d's coordinates (consensus rule)
                if d != c and any(not conflict(moved, td[1]) and (not leak_overlap or any(p in td[1] for p in moved))
                                  for td in tclass[d]):
                    sl += [(d, x) for x in range(capacity(*span[d], sex))]
            lists.append(sl)
        dfc = deficiency(lists)
        if dfc > 0:
            t1.append(f"{c}:{len(q)}/{cap}:d{dfc}")
    novel, novel_t, seen = {}, {}, set()
    for a in fam_copies:                                        # anchor copy
        pats = {q: tp[q].get(a, {}) for q, cs in cand.items() if a in cs}
        for reads, _, _ in classes(pats, merge):
            key = tuple(reads)
            if key in seen:
                continue
            seen.add(key)
            union = sorted({c for q in reads for c in cand[q]})
            is_novel = True
            for c in union:
                kc = klass([q for q in reads if c in cand[q]], {q: tp[q].get(c, {}) for q in reads})
                if not all(class_conflict(kc, tc, pair) for tc in tclass.get(c, [])):
                    is_novel = False
                    break
            if is_novel:
                novel[key] = union
                if all(tclass.get(c) for c in union):
                    novel_t[key] = union
    fl = _tsv(os.path.join(out_dir, f"{f}.flags.tsv")) if os.path.exists(os.path.join(out_dir, f"{f}.flags.tsv")) else []

    def summ(d):
        sizes = sorted((len(k) for k in d), reverse=True)
        return [len(d), sizes[0] if sizes else 0, len({q for k in d for q in k}), ",".join(str(s) for s in sizes[:10])]

    return [f, len(cand), len(fam_copies), sum(1 for c in fam_copies if tclass[c]), sum(len(teach.get(c, {})) for c in fam_copies),
            int(bool(t1)), ";".join(t1), *summ(novel), *summ(novel_t), len(fl), len({q for r in fl for q in r["reads"].split(",")})]


HEADER = ("family tied_reads copies taught_copies teachers tier1 tier1_copies tier2_patterns tier2_max_reads tier2_reads "
          "tier2_sizes tier2t_patterns tier2t_max_reads tier2t_reads tier2t_sizes flags flag_reads").split()


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out_dir")
    ap.add_argument("copies")
    ap.add_argument("genome")
    ap.add_argument("out_tsv")
    ap.add_argument("--sex", choices=("male", "female", "unknown"), default="unknown")
    ap.add_argument("--exclude-editing", action="store_true")
    ap.add_argument("--families", default=None)
    ap.add_argument("--merge-single", action="store_true", help="merge supported classes one position apart by single "
                    "linkage (C_certificate §1) before the clique and the novelty test")
    ap.add_argument("--leak-closure", action="store_true", help="tier 1 needs a Hall deficiency > 0 once each clique class may "
                    "also use the slots of any other copy with a compatible teacher class (C_certificate §1, S4)")
    ap.add_argument("--leak-needs-overlap", action="store_true", help="with --leak-closure: a class is compatible with a "
                    "copy's class only when they share >= 1 position (default: no shared position counts as compatible)")
    ap.add_argument("--pair-certified", action="store_true", help="class conflicts need two positions carried together by "
                    ">= 2 members of each class (C_certificate §1); default: the consensuses differ at >= 2 positions")
    a = ap.parse_args(argv)
    import pysam
    P, C = _sibling("o2psv"), _sibling("o2correct")
    cps_all = P.load_copies(a.copies)
    idx = C.span_index(cps_all)
    span = {c: (ch, s, e) for ch, v in cps_all.items() for c, _, _, (s, e) in v}
    fam = collections.defaultdict(list)
    for c in sorted(span):
        fam[c.rsplit(".", 1)[0]].append(c)
    cand = collections.defaultdict(dict)
    for r in _tsv(os.path.join(a.out_dir, "calls.tsv")):
        cand[r["family"]][r["read"]] = r["candidates"].split(",")
    todo = sorted(cand) if a.families is None else sorted(set(a.families.split(",")) & set(cand))
    genome = pysam.FastaFile(a.genome)
    with open(a.out_tsv, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for f in todo:
            row = screen_family(f, a.out_dir, cand[f], fam[f], span, cps_all, idx, genome, a.sex, a.exclude_editing,
                                a.pair_certified, a.merge_single, a.leak_closure, a.leak_needs_overlap)
            fh.write("\t".join(str(x) for x in row) + "\n")
            fh.flush()
            print(f"[o3screen] {f}: tier1 {row[5]} {row[6]}; tier2 {row[7]} ({row[9]} reads); tier2t {row[11]}; flags {row[15]}",
                  file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
