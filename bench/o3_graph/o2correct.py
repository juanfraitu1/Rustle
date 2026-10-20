#!/usr/bin/env python3
"""O2 as a corrector of minimap2 (2026-10-20; user: "the entire point is to use the tied and nearly tied AS alignments and find
where they truly belong as an addendum to minimap2 since it gave everything a coin toss ... there must be a way with PSVs,
de, nm, md tags to find the true placement of the multimappers").

For a read whose alignments within 10 AS of its best (Eichler) lie on two or more copies of one family and on no locus off
it, the corrector takes the copy whose best tied alignment has the FEWEST SUBSTITUTIONS (mismatches at aligned bases: X ops,
or NM less its indel bases), and keeps minimap2's primary when two copies share that minimum. minimap2's AS charges an indel
as it charges a substitution, and a HiFi read's errors are mostly homopolymer indels; a PSV is a substitution. So AS can
prefer a copy that differs from the read at more PSVs (HG002 development, the donor's assembly as truth: on 8,732 near-tied
reads the corrector is right 81.5% vs minimap2 74.7% -- 604 of minimap2's errors fixed, 18 right picks broken; NM and de do
not separate the copies). Exact ties (copies identical over the read) stay minimap2's coin toss.

But 588 of those 605 fixes are two families (SMN, THOC3) where the donor's copy carries a 1-bp insertion CHM13's lacks; in
the other families the rule fixed 17 and broke 18. One PSV mismatch against one indel is no evidence either way until the
sample says which copy carries the indel. The sample-aware corrector (donor) asks it. Per copy, the sample's confident reads
(primary on it, every alignment on another copy of its family > 9 AS lower) give its events -- substitutions and
left-aligned indels. Per tied copy, a near-tied read scores the events of its alignment the copy's reads do not explain (an
event is explained when >= SHARED = 2 of them carry it and >= HET = 10% of their depth, depth >= K = 3), plus the events
>= HOM = 80% of them carry that the read lacks. The lowest score takes the read, but overrules minimap2's primary only when
it is lower even granting the primary every event its reads do not cover (dominance: no copy wins for want of coverage);
a shared lowest keeps minimap2's pick. In one genome a copy's reads explain only their shared errors, so the rule is
minimap2's ordering by events. HG002 development (HET 0.2/0.1/0.05 x dominance; donor truth, PAR out): near ties 87.3%
right vs minimap2's 66.4% (1,161 fixed, 15 broken), exact ties 52.6% vs 39.6% (97 / 22); 15 families gain, 4 lose.
Held out (PREREG_o2_corrector_2026-10-20, register 1343): FAIL overall -- fresh A119b fragments 81.7% vs 71.5% (130 / 3,
40 families gain, 1 loses) and no harm in one genome, but a copy the reference lacks (BRANCH, REMOVED) teaches its nearest
copy foreign alleles and the rule then moves that copy's neighbours' reads (heldout22: 5 near-tie breaks, ANKRD30BL).

    o2correct.py run BAM COPIES.tsv OUT.tsv           read, family, corrector, minimap2, reason, candidates
    o2correct.py donor BAM COPIES.tsv GENOME.fa OUT.tsv [--teach LIBRARY.bam] [--cache DIR --budget SECONDS]
                                   the sample-aware corrector (also the AS margin and per-copy scores); --teach: the library
                                   a fragment set was cut from, each fragment's source read left out of its own judgement;
                                   --cache/--budget: stop after SECONDS, keeping each family done, and resume (exit 3)
    o2correct.py score OUT.tsv TRUTH.tsv ... [--donor] [--among-candidates] [--families Q.json] [--exclude PAR.json]
    o2correct.py score-worlds OUT.json WORLD_DIR ... [--exclude PAR.json]     simulated worlds pooled, truth among the tied
"""
import bisect
import collections
import csv
import importlib.util
import json
import os
import pickle
import sys
import time

SLACK = 9
K, SHARED, HET, HOM = 3, 2, 0.1, 0.8
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
HERE = os.path.dirname(os.path.abspath(__file__))


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_corr", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def tied(alns):
    """({copy: its best-AS alignment within SLACK of the read's best}, minimap2's primary copy), the first None when an
    alignment within SLACK lies off every copy, on another family, or the read is tied on fewer than two copies."""
    prim = next((a["copy"] for a in alns if a["prim"]), None)
    best = max(a["AS"] for a in alns)
    per = {}
    for a in alns:
        if a["AS"] < best - SLACK:
            continue
        if a["copy"] is None:
            return None, prim
        if a["copy"] not in per or a["AS"] > per[a["copy"]]["AS"]:
            per[a["copy"]] = a
    if len(per) < 2 or len({c.rsplit(".", 1)[0] for c in per}) > 1 or prim not in per:
        return None, prim
    return per, prim


def pick(per, prim):
    """(copy, reason): the tied copy with the fewest substitutions, or minimap2's primary when that minimum is shared."""
    m = min(a["sub"] for a in per.values())
    top = [c for c, a in per.items() if a["sub"] == m]
    return (top[0], "fewest_substitutions") if len(top) == 1 else (prim, "minimap2")


def tally(rows):
    """Counts over (corrector, minimap2, truth) triples: right for each, and where they disagree, the corrector's fixes (it
    right, minimap2 wrong), breaks (it wrong, minimap2 right) and both wrong."""
    n = collections.Counter()
    for c, m, t in rows:
        n["reads"] += 1
        n["corrector right"] += c == t
        n["minimap2 right"] += m == t
        if c != m:
            n["fixes" if c == t else "breaks" if m == t else "both wrong"] += 1
    return n


def span_index(cps):
    """{chrom: (span starts sorted, [(start, end, table index, copy row)], longest span)} for overlapping()."""
    out = {}
    for ch, rows in cps.items():
        v = sorted((r[3][0], r[3][1], i, r) for i, r in enumerate(rows))
        out[ch] = ([x[0] for x in v], v, max((x[1] - x[0] for x in v), default=0))
    return out


def overlapping(idx, chrom, lo, hi):
    """The copy rows whose span overlaps [lo, hi), in their table order (summary.attribute breaks ties by it)."""
    if chrom not in idx:
        return []
    st, v, longest = idx[chrom]
    i, j = bisect.bisect_left(st, lo - longest), bisect.bisect_left(st, hi)
    return [x[3] for x in sorted((x for x in v[i:j] if x[1] > lo), key=lambda x: x[2])]


def read_span(cigar):
    """(first, last + 1) read position an alignment aligns, in its own orientation, clipped bases counted."""
    lead = 0
    for op, n in cigar:
        if op not in (4, 5):
            break
        lead += n
    return lead, lead + sum(n for op, n in cigar if op in (0, 1, 7, 8))


def common_window(spans, length):
    """{copy: (lo, hi)} in each alignment's own read coordinates: the read bases every tied alignment aligns. spans: {copy:
    ((first, last + 1) in its own orientation, whether that orientation is the primary's reverse)}."""
    prim = {c: ((length - e, length - s) if rev else (s, e)) for c, ((s, e), rev) in spans.items()}
    lo, hi = max(s for s, _ in prim.values()), min(e for _, e in prim.values())
    hi = max(lo, hi)
    return {c: ((length - hi, length - lo) if rev else (lo, hi)) for c, (_, rev) in spans.items()}


def alignment_events(start, cigar, seq, fetch, window=None):
    """(events, blocks) of one alignment against the reference: substitutions ("S", position, read base), insertions ("I",
    position, inserted bases) and deletions ("D", position, length), each indel shifted left as far as the reference repeats
    it within its block (so every read writes an indel one way), and its aligned blocks [(start, end)] between introns. seq is
    the read in the alignment's orientation, clipped bases included; fetch(lo, hi) -> the reference; window = (lo, hi): only
    the read positions lo <= q < hi count (common_window: tied alignments compared over the same read bases)."""
    lo, hi = window if window is not None else (0, len(seq) + 1)

    def at(x):
        return fetch(x, x + 1).upper()

    ev, seg = set(), []

    def cover(a, b):
        if a >= b:
            return
        if seg and seg[-1][1] == a:
            seg[-1] = (seg[-1][0], b)
        else:
            seg.append((a, b))

    q, r, b0 = 0, start, start
    for op, n in cigar:
        if op in (0, 7, 8):
            j0, j1 = max(0, lo - q), min(n, hi - q)
            if j0 < j1:
                a, b = seq[q + j0:q + j1].upper(), fetch(r + j0, r + j1).upper()
                if a != b:
                    ev.update(("S", r + j0 + j, x) for j, (x, y) in enumerate(zip(a, b)) if x != y and "N" not in (x, y))
                cover(r + j0, r + j1)
            q, r = q + n, r + n
        elif op == 1:
            if lo <= q < hi:
                ins, p = seq[q:q + n].upper(), r
                while p > b0 and at(p - 1) == ins[-1]:
                    ins, p = ins[-1] + ins[:-1], p - 1
                ev.add(("I", p, ins))
            q += n
        elif op == 2:
            if lo <= q < hi:
                p = r
                while p > b0 and at(p - 1) == at(p + n - 1):
                    p -= 1
                ev.add(("D", p, n))
                cover(r, r + n)
            r += n
        elif op == 3:
            r += n
            b0 = r
        elif op in (4, 5):
            q += n
    return ev, seg


class Donor:
    """The sample's own copies: per copy, how many of its confident reads carry each event, and their aligned blocks. own =
    (copy, events, blocks) of one read leaves that read out (a fragment's source)."""

    def __init__(self):
        self.ev = collections.defaultdict(collections.Counter)
        self.st, self.en = collections.defaultdict(list), collections.defaultdict(list)
        self._ready, self._near, self._hom = True, {}, {}

    def add(self, copy, events, blocks):
        self.ev[copy].update(events)
        for s, e in blocks:
            self.st[copy].append(s)
            self.en[copy].append(e)
        self._ready, self._near, self._hom = False, {}, {}

    def count(self, copy, e, own=None):
        n = self.ev[copy][e] if copy in self.ev else 0
        return n - (own is not None and own[0] == copy and e in own[1])

    def depth(self, copy, pos, own=None):
        if not self._ready:
            for c in self.st:
                self.st[c].sort()
                self.en[c].sort()
            self._ready = True
        d = bisect.bisect_right(self.st.get(copy, []), pos) - bisect.bisect_right(self.en.get(copy, []), pos)
        if own is not None and own[0] == copy:
            d -= sum(1 for s, t in own[2] if s <= pos < t)
        return d

    def explained(self, copy, e, own=None):
        n = self.count(copy, e, own)
        if n < SHARED:
            return False
        d = self.depth(copy, e[1], own)
        return d >= K and n >= HET * d

    def _carried(self, copy, e, own):
        n, d = self.count(copy, e, own), self.depth(copy, e[1], own)
        return n >= SHARED and d >= K and n >= HOM * d

    def hom(self, copy, own=None):
        """The events nearly all of the copy's reads carry (>= HOM of a depth >= K)."""
        if copy not in self._near:                              # those that are, or without one read could be
            self._near[copy] = sorted(e for e, n in self.ev.get(copy, {}).items()
                                      if n >= SHARED and n >= HOM * (self.depth(copy, e[1]) - 1))
        if own is None or own[0] != copy:
            if copy not in self._hom:
                self._hom[copy] = [e for e in self._near[copy] if self._carried(copy, e, None)]
            return self._hom[copy]
        return [e for e in self._near[copy] if self._carried(copy, e, own)]


def donor_parts(copy, events, blocks, donor, own=None):
    """(refuted, lacking, unknown): the read's events the copy's own reads cover >= K times and do not explain; the events
    nearly all of them carry that the read lacks over its aligned blocks; the read's events they do not cover."""
    refuted = unknown = 0
    for e in events:
        if donor.explained(copy, e, own):
            continue
        if donor.depth(copy, e[1], own) >= K:
            refuted += 1
        else:
            unknown += 1
    lack = sum(1 for e in donor.hom(copy, own) if e not in events and any(s <= e[1] < t for s, t in blocks))
    return refuted, lack, unknown


def donor_score(copy, events, blocks, donor, own=None):
    return sum(donor_parts(copy, events, blocks, donor, own))


def donor_pick(per, prim, donor, own=None):
    """(copy, reason) from {copy: (events, blocks)}: the copy with the lowest score (donor_parts summed), kept against
    minimap2's primary only when it wins even if every event of the primary's its reads do not cover were the donor's
    ("sample"); a shared lowest ("tie") or an alternative that does not dominate ("undominated") keeps minimap2's pick."""
    pt = {c: donor_parts(c, ev, bl, donor, own) for c, (ev, bl) in per.items()}
    sc = {c: sum(v) for c, v in pt.items()}
    m = min(sc.values())
    top = [c for c, v in sc.items() if v == m]
    if len(top) > 1:
        return prim, "tie"
    if top[0] != prim and not m < pt[prim][0] + pt[prim][1]:
        return prim, "undominated"
    return top[0], "sample"


def cmd_run(bam_p, copies_p, out_p):
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    cps = P.load_copies(copies_p)
    alns = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            sub, _ = P.substitutions(rec.cigartuples, rec.get_tag("NM") if rec.has_tag("NM") else 0)
            alns[rec.query_name].append(dict(AS=rec.get_tag("AS"), sub=sub, prim=not rec.is_secondary,
                                             copy=S.attribute(rec.get_blocks(), P._strand(rec), cps.get(rec.reference_name, ()))))
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\tcorrector\tminimap2\treason\tcandidates\n")
        for q, al in sorted(alns.items()):
            per, prim = tied(al)
            if per is None:
                continue
            c, why = pick(per, prim)
            n[why] += 1
            n["overrules"] += c != prim
            fh.write(f"{q}\t{prim.rsplit('.', 1)[0]}\t{c}\t{prim}\t{why}\t{','.join(sorted(per))}\n")
    print(f"[o2correct] {bam_p} -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def collect(bam_p, copies_p, genome_p, teach_p=None, cache=None, budget=None):
    """(domain {read: (tied copies, minimap2's copy, AS margin between its two best copies)}, the sample's Donor, {read:
    {copy: (events, blocks)}}, {source read: (copy, events, blocks)}) from the reads' BAM and the sample's (teach_p: the
    library a fragment set was cut from; default the reads' own). The sample's confident reads on a copy: primary on it,
    every alignment on another copy of its family > SLACK AS lower (o2psv.teacher_bases' rule), found family by family over
    the copies' spans so a large library is never read whole; a fragment's source read is kept apart to leave it out. cache: a
    directory keeping each family's confident reads as it is done, so a run stopped at `budget` seconds resumes there (None
    is returned for an unfinished run)."""
    import pysam
    t0 = time.time()

    def log(msg):
        print(f"[o2correct] {time.time() - t0:7.1f}s {msg}", file=sys.stderr, flush=True)

    P, S = _sibling("o2psv"), _sibling("summary")
    cps = P.load_copies(copies_p)
    span = {c: (ch, s, e) for ch, v in cps.items() for c, _, _, (s, e) in v}
    fam = collections.defaultdict(list)
    for c in span:
        fam[c.rsplit(".", 1)[0]].append(c)
    fa = pysam.FastaFile(genome_p)

    idx = span_index(cps)

    def copy_of(rec):
        cands = overlapping(idx, rec.reference_name, rec.reference_start, rec.reference_end)
        return S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None

    def fetcher(chrom):
        return lambda lo, hi: fa.fetch(chrom, lo, hi)

    def records(bam):
        for rec in bam.fetch(until_eof=True):
            if not (rec.is_unmapped or rec.is_supplementary) and rec.has_tag("AS"):
                yield rec

    alns = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:                      # pass 1: the near-tied reads
        for rec in records(bam):
            alns[rec.query_name].append(dict(AS=rec.get_tag("AS"), prim=not rec.is_secondary, copy=copy_of(rec)))
    domain = {}
    for q, al in alns.items():
        per, prim = tied(al)
        if per is not None:
            a = sorted((x["AS"] for x in per.values()), reverse=True)
            domain[q] = (set(per), prim, a[0] - a[1])
    del alns
    log(f"pass 1: {len(domain)} near-tied reads")
    recs, seqs = collections.defaultdict(dict), {}
    with pysam.AlignmentFile(bam_p) as bam:                      # pass 2: their alignments
        for rec in records(bam):
            q = rec.query_name
            if q not in domain:
                continue
            c, a = copy_of(rec), rec.get_tag("AS")
            if c in domain[q][0] and (c not in recs[q] or a > recs[q][c][0]):
                recs[q][c] = (a, rec.reference_start, tuple(rec.cigartuples), rec.is_reverse, rec.reference_name)
            if not rec.is_secondary:
                seqs[q] = (rec.query_sequence, rec.is_reverse)
    events = {}
    for q, (cset, prim, _) in domain.items():
        if q in seqs and set(recs[q]) == cset:
            seq, prev = seqs[q]
            win = common_window({c: (read_span(cig), rev != prev) for c, (_, _, cig, rev, _) in recs[q].items()}, len(seq))
            events[q] = {c: alignment_events(st, cig, seq if rev == prev else seq.translate(COMP)[::-1], fetcher(ch), win[c])
                         for c, (_, st, cig, rev, ch) in recs[q].items()}
    log(f"pass 2: events of {len(events)} reads")
    need_iv = collections.defaultdict(list)                     # per copy, the reference the near-tied reads align
    for per in events.values():
        for c, (_, bl) in per.items():
            need_iv[c].extend(bl)
    for c, iv in need_iv.items():
        iv.sort()
        merged = [list(iv[0])]
        for x, y in iv[1:]:
            if x <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], y)
            else:
                merged.append([x, y])
        need_iv[c] = ([x for x, _ in merged], [y for _, y in merged])

    def needed(c, lo, hi):
        st, en = need_iv[c]
        i = bisect.bisect_right(st, hi - 1) - 1
        return i >= 0 and en[i] > lo

    def region(chrom, lo, hi):                                  # one reference fetch per copy span
        lo = max(0, lo)
        seq = fa.fetch(chrom, lo, hi).upper()
        return lambda a, b: seq[a - lo:b - lo] if lo <= a and b <= lo + len(seq) else fa.fetch(chrom, a, b)

    sources = {P.own_name(q) for q in domain}
    donor, own, n = Donor(), {}, 0
    if cache:
        os.makedirs(cache, exist_ok=True)
    bam = pysam.AlignmentFile(teach_p or bam_p)

    def family_sample(f):
        """[(copy, events, blocks, the read's name when it is a source, else None)] of the family's confident reads."""
        cand = {}
        for c in fam[f]:                                        # primaries over what the near-tied reads align
            if c not in need_iv:
                continue
            ch, s, e = span[c]
            fetch = region(ch, s - 50000, e + 50000)
            for rec in bam.fetch(ch, s, e):
                q = rec.query_name
                if (rec.is_secondary or rec.is_supplementary or rec.is_unmapped or q in cand or q in domain
                        or not needed(c, rec.reference_start, rec.reference_end) or not rec.has_tag("AS")
                        or copy_of(rec) != c or rec.query_sequence is None):
                    continue
                cand[q] = (c, rec.get_tag("AS"), alignment_events(rec.reference_start, rec.cigartuples,
                                                                  rec.query_sequence, fetch))
        best = collections.defaultdict(dict)                    # their alignments on the family's other copies
        for c in fam[f]:
            ch, s, e = span[c]
            for rec in bam.fetch(ch, s, e):
                q = rec.query_name
                if q not in cand or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                    continue
                cp = copy_of(rec)
                if cp is not None and cp != cand[q][0] and cp.rsplit(".", 1)[0] == f:
                    best[q][cp] = max(best[q].get(cp, -10 ** 9), rec.get_tag("AS"))
        return [(cp, ev, bl, q if q in sources else None) for q, (cp, a, (ev, bl)) in cand.items()
                if a - max(best[q].values(), default=-10 ** 9) > SLACK]

    with bam:                                                   # the sample's confident reads, family by family
        for i, f in enumerate(sorted({c.rsplit(".", 1)[0] for c in need_iv})):
            if i % 25 == 0:
                log(f"sample: family {i} ({f}), {n} confident reads so far")
            done = os.path.join(cache, f"{f}.pkl") if cache else None
            if done and os.path.exists(done):
                part = pickle.load(open(done, "rb"))
            elif budget is not None and time.time() - t0 > budget:
                log(f"stopped at the {budget} s budget before family {i} ({f}); rerun to resume from {cache}")
                return None
            else:
                part = family_sample(f)
                if done:
                    pickle.dump(part, open(done + ".tmp", "wb"))
                    os.replace(done + ".tmp", done)
            for cp, ev, bl, q in part:
                donor.add(cp, ev, bl)
                n += 1
                if q is not None:
                    own[q] = (cp, ev, bl)
    print(f"[o2correct] {bam_p}: near-tied {len(domain)} ({len(events)} with every tied alignment); the sample "
          f"({teach_p or bam_p}): {n} confident reads on {len(need_iv)} copies, {len(own)} of them sources", flush=True)
    return domain, donor, events, own


def cmd_donor(bam_p, copies_p, genome_p, out_p, teach_p=None, cache=None, budget=None):
    P = _sibling("o2psv")
    got = collect(bam_p, copies_p, genome_p, teach_p, cache, budget)
    if got is None:
        return 3
    domain, donor, events, own = got
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\tcorrector\tminimap2\treason\tcandidates\tmargin\tscores\n")
        for q, per in sorted(events.items()):
            cset, prim, margin = domain[q]
            o = own.get(P.own_name(q))
            c, why = donor_pick(per, prim, donor, o)
            n[why] += 1
            n["overrules"] += c != prim
            sc = ",".join(f"{x}:{donor_score(x, *per[x], donor, o)}" for x in sorted(per))
            fh.write(f"{q}\t{prim.rsplit('.', 1)[0]}\t{c}\t{prim}\t{why}\t{','.join(sorted(per))}\t{margin}\t{sc}\n")
    print(f"[o2correct] donor -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_score(out_p, truth_ps, donor=False, among=False, fam_p=None, excl_p=None):
    """Fixes and breaks against the truth, pooled and per family, for every read and by AS band (exact tie: margin 0; near
    tie: 1-9) when the output has a margin; fam_p keeps only those families, excl_p drops these (both JSON lists or keys)."""
    P = _sibling("o2psv")
    truth = P._truth(truth_ps, donor)
    fams = set(json.load(open(fam_p))) if fam_p else None
    excl = set(json.load(open(excl_p))) if excl_p else set()
    per_fam, band, srcs = collections.defaultdict(list), collections.defaultdict(list), set()
    for r in csv.DictReader(open(out_p), delimiter="\t"):
        t = truth.get(r["read"])
        if t is None:
            continue
        cls, tc = t if donor else (None, t)
        if donor and cls not in ("unique", "unique_elsewhere"):
            continue
        if among and tc not in r["candidates"].split(","):
            continue
        if (fams is not None and r["family"] not in fams) or r["family"] in excl:
            continue
        per_fam[r["family"]].append((r["corrector"], r["minimap2"], tc))
        if r["corrector"] != r["minimap2"]:
            srcs.add(r["read"].split("|")[0])
        if r.get("margin") not in (None, ""):
            band["exact tie" if r["margin"] == "0" else "near tie"].append((r["corrector"], r["minimap2"], tc))
    tot = tally([x for v in per_fam.values() for x in v])
    fam = {f: dict(tally(v)) for f, v in per_fam.items()}
    better = sorted(f for f, v in fam.items() if v.get("fixes", 0) > v.get("breaks", 0))
    worse = sorted(f for f, v in fam.items() if v.get("breaks", 0) > v.get("fixes", 0))

    def line(name, t):
        k = max(1, t["reads"])
        return (f"{name}: reads {t['reads']}; corrector right {t['corrector right']} ({t['corrector right'] / k:.4f}), "
                f"minimap2 right {t['minimap2 right']} ({t['minimap2 right'] / k:.4f}); fixes {t['fixes']}, breaks {t['breaks']}, "
                f"both wrong {t['both wrong']}")
    print(f"[o2correct] {out_p}: " + line("all", tot) + f"; overruled reads from {len(srcs)} source reads; families fixed "
          f"more than broken {len(better)}, broken more than fixed {len(worse)}", flush=True)
    for b, rows in sorted(band.items()):
        print("    " + line(b, tally(rows)), flush=True)
    json.dump(dict(pooled=dict(tot), bands={b: dict(tally(v)) for b, v in band.items()}, overruled_sources=len(srcs),
                   families=fam, families_better=better, families_worse=worse), open(out_p + ".score.json", "w"), indent=1)
    return 0


def score_worlds(dirs, excl=frozenset()):
    """{"exact tie" / "near tie": tally} pooled over simulated worlds (directories with o2correct.tsv and reads.truth.E0.tsv),
    a read scored when its truth is one of its tied copies; a world's read names stay its own."""
    band = collections.defaultdict(list)
    for w in dirs:
        truth = {r["read"]: r["copy"] for r in csv.DictReader(open(os.path.join(w, "reads.truth.E0.tsv")), delimiter="\t")}
        for r in csv.DictReader(open(os.path.join(w, "o2correct.tsv")), delimiter="\t"):
            t = truth.get(r["read"])
            if t is None or t not in r["candidates"].split(",") or r["family"] in excl:
                continue
            band["exact tie" if r["margin"] == "0" else "near tie"].append((r["corrector"], r["minimap2"], t))
    return {b: tally(v) for b, v in band.items()}


def cmd_score_worlds(out_p, dirs, excl_p=None):
    excl = set(json.load(open(excl_p))) if excl_p else set()
    skipped = [d for d in dirs if not os.path.exists(os.path.join(d, "o2correct.tsv"))]
    dirs = [d for d in dirs if d not in skipped]
    if skipped:
        print(f"[o2correct] {len(skipped)} worlds without o2correct.tsv, not scored: {', '.join(skipped)}", flush=True)
    res = score_worlds(dirs, excl)
    for b, t in sorted(res.items()):
        k = max(1, t["reads"])
        print(f"[o2correct] {len(dirs)} worlds, {b}: reads {t['reads']}; corrector right {t['corrector right'] / k:.4f}, "
              f"minimap2 right {t['minimap2 right'] / k:.4f}; fixes {t['fixes']}, breaks {t['breaks']}, both wrong "
              f"{t['both wrong']}", flush=True)
    json.dump(dict(worlds=dirs, skipped=skipped, bands={b: dict(t) for b, t in res.items()}), open(out_p, "w"), indent=1)
    return 0


def main(argv):
    donor, argv = "--donor" in argv, [a for a in argv if a != "--donor"]
    among, argv = "--among-candidates" in argv, [a for a in argv if a != "--among-candidates"]
    fam_p = excl_p = None
    if "--families" in argv:
        i = argv.index("--families")
        fam_p, argv = argv[i + 1], argv[:i] + argv[i + 2:]
    if "--exclude" in argv:
        i = argv.index("--exclude")
        excl_p, argv = argv[i + 1], argv[:i] + argv[i + 2:]
    if argv[:1] == ["run"] and len(argv) == 4:
        return cmd_run(*argv[1:])
    opt = {}
    for k in ("--teach", "--cache", "--budget"):
        if k in argv:
            i = argv.index(k)
            opt[k], argv = argv[i + 1], argv[:i] + argv[i + 2:]
    if argv[:1] == ["donor"] and len(argv) == 5:
        return cmd_donor(*argv[1:], teach_p=opt.get("--teach"), cache=opt.get("--cache"),
                         budget=float(opt["--budget"]) if "--budget" in opt else None)
    if argv[:1] == ["score-worlds"] and len(argv) >= 3:
        return cmd_score_worlds(argv[1], argv[2:], excl_p)
    if argv[:1] == ["score"] and len(argv) >= 3:
        return cmd_score(argv[1], argv[2:], donor, among, fam_p, excl_p)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
