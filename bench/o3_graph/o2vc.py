#!/usr/bin/env python3
"""O2 from the sample's own variant calls (design docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md, plan
docs/superpowers/plans/2026-10-22-o2-variant-calls.md): the CLI over o2vc_rule.py (the decision rule) and o2vc_data.py (the
BAM, VCF and caller side).

    o2vc.py run BAM COPIES.tsv GENOME.fa OUTDIR [--families A,B] [--exclude-chrM] [--threads N]
    o2vc.py score OUTDIR [TRUTH...] --mode donor|none [--n4 O2CORRECT.tsv] [--panel A|C|E] [--focus FAMILY] [--exclude PAR.json]

run never reads truth. Per family, in sorted order -- a family with no AS-tied read is skipped, and under --exclude-chrM a
family with a copy on chrM:
1. its copies' teachers (o2vc_data.scan) go to OUTDIR/teach/<family>.bam, base qualities of 0 written as 1 (longcallR
   2.0.0 panics on a Q0 at a variant site);
2. longcallR runs once per copy with >= 1 teacher (caller_region: its span and PAD bases each side, START >= 1 and END
   capped at the contig's length + 1, past which the caller panics) against the copy's chromosome alone,
   OUTDIR/ref/<chrom>.fa, into OUTDIR/vcf/<copy>.vcf. A VCF already there is reused, never re-called:
   the caller is not deterministic, so the saved VCFs are the record. A copy whose call fails is unknown everywhere
   (run.json failed_copies) and the run goes on;
3. every record's REF is checked against the genome, upper-cased: a mismatch stops the run;
4. each tied read is held against every tied copy's genotype at its columns (o2vc_data.columns), in the read's frame. VC
   takes the sample's calls (genotype(): the record's GT, else the reference base); N3, the reference base wherever VC's
   genotype is known, isolates what the calls add. Foreign alleles that recur and co-occur across the family's reads flag a
   possible missing copy (flags()), and a flagged read is not placed (vc empty);
5. the teacher side is reported, never counted: FILTER counts per copy, and flags over the teachers' own-copy foreign alleles.

A copy is unknown at a position (None in every column's g, never left out of it) when it has no teacher, the position lies
outside the region the caller was given, the depth -- the record's DP, else o2vc_data.teacher_depth, the caller's own count
-- is under MIN_DEPTH or over the caller's MAX_DEPTH, or the record is RnaEdit or Dense or calls no allele.

Files in OUTDIR (positions 0-based; a copy's alleles in columns.tsv in the read's frame, as the rule compares them, and in
flags.tsv in the genome's strand):
- calls.tsv: read family candidates (the tied copies) minimap2 (its primary's copy, "" off every copy) mapq margin (AS,
  the two best tied copies) outside (1: the primary is not a tied copy) vc vc_set vc_cols n3 n3_set flag -- vc / n3 the
  placed copy or "", *_set the copies no copy beats (empty on a beat cycle), vc_cols the separating columns voting for vc;
  flag "<family>#<id>" (comma-joined when several), and then vc is "" (vc_set and vc_cols stay the rule's).
- <family>.reads.tsv: read copy primary AS NM de ms subs indels ref_for ref_against call_for call_against -- per tied read
  and tied copy: the alignment's tags, its substitutions (o2psv.substitutions) and indel bases, and the columns voting for
  and against the copy under N3 (ref_*) and VC (call_*).
- <family>.columns.tsv: read p base qual kind copies -- per tied read, the columns that separate two copies under VC or
  N3 or hold a foreign allele. kind is the first that holds: separating_call (VC separates a pair of copies the
  reference does not), separating_ref (every pair VC separates, the reference does too), erased_ref (the reference
  separates a pair and VC none: the sample's calls share an allele there, a reference PSV this sample does not carry;
  N3 knows only what VC knows), foreign or editing (o2vc_rule.foreign(); a foreign allele feeds flags() whatever the
  kind). copies: copy|chrom|pos|ref|genotype|depth per tied copy, ";"-joined (genotype "." when unknown, alleles
  "/"-joined; depth "." when not measured).
- <family>.flags.tsv: flag n_reads reads alleles (copy|chrom|pos|base, ";"-joined).
- <family>.teachers.tsv: copy teachers records pass lowqual homref multiallelic rnaedit dense teacher_flag_reads -- the
  FILTERs of the copy's VCF, and its teachers in a flag of the teachers' own-copy foreign alleles.
- vcf/<copy>.vcf, teach/<family>.bam, ref/<chrom>.fa, run.json.

score is the one step that reads truth. It reads OUTDIR's calls.tsv, <family>.reads.tsv and <family>.columns.tsv; --mode donor also
reads TRUTH... (o2psv._truth: read, class, truth -- HG002's truth.tsv and truth_untied.tsv together, the later file winning a
read both hold), --mode none (DAZ: no truth exists) takes none and says nothing of right or wrong. It writes
OUTDIR/score[.PANEL[.FOCUS]].json and, for every family scored (a --exclude family too, apart), <family>.case.tsv and
<family>.summary.md.
- Arms, all on the same reads. N1 = calls.tsv's minimap2: always a placement, so a primary on no copy or outside the tied copies
  is a wrong placement unless it is the truth copy; N2 = N1 at mapq >= 1, else it abstains; N3 = n3; N4 = the corrector's pick
  (--n4; N1 for a read it does not hold, and how many tied reads it holds is reported: n4_held); VC = vc. Only a non-empty vc /
  n3 is a placement: an empty set (a beat cycle) and a flagged read abstain.
- Scored reads (donor; the amended spec, 2026-10-22): the tied reads of class unique, unique_elsewhere, elsewhere or
  no_counterpart. A placement is right only when the placed copy is the read's truth copy, which only unique and
  unique_elsewhere reads have; every other placement on a scored read is wrong -- a truth copy outside the tied copies, an
  elsewhere read, a no_counterpart read. Per arm over ALL scored reads, the same for every arm: placed, right, wrong,
  abstained; fixes (placed right, N1 wrong), breaks (placed wrong, N1 right: minimap2's primary can be the truth copy outside
  the tied set, and VC moves it), both_wrong (placed wrong, N1 wrong, on the same copy or not; so placed = both right + fixes +
  breaks + both_wrong, and N1's own three are against itself). Precision is wrong / placed, always with both counts. The
  decidable subset (the truth copy is one of the tied copies, the earlier registers' reads) is reported beside it (decidable).
  Not scored: own_tie (K1 only), poor_match, unmapped and reads with no truth row, counted only; the VC placements among them
  are listed (vc.not_scored) so that the arm table cannot hide them.
- Calibration (donor): K1 = own_tie reads VC places, over the own_tie reads; K2 = flagged reads over the reads whose class has a
  copy in the reference (unique, unique_elsewhere, own_tie); C1 = flags holding >= 2 no_counterpart reads; D2 = scored reads N1
  places right that VC moves or flags. The same numbers are kept per family and pooled over the families not in --exclude.
- Bars: A (--focus family: A1 VC wrong/placed <= 2%, A2 fixes >= 10, A3 VC fixes > N3 fixes; pooled: K1 <= 2%, K2 <= 2%),
  C (C1 >= 1), E (pooled: E1 wrong/placed <= 2%, E2 fixes > breaks, K2 <= 2%), all over the scored reads. Ratios are compared
  exactly. A bar with nothing to measure is not evaluable (pass null, with a reason): a zero denominator, no no_counterpart
  tied read for C1, no scored read on the focus family for A2 and A3, no VC placement of a scored read for E2. The verdict of a
  panel is FAIL when a bar failed, else UNDERPOWERED when a bar was not evaluable, else PASS. --focus belongs to --panel A and
  must be a family of the run.
- Not built (deferred with panels D and C-gorilla): --mode world and units, --panel D and Cg, score-worlds.
- The summary and the case file record the score options (mode, truth files, --n4, --panel, --focus, --exclude). The summary:
  the arm table; the VC placements resting on one column (vc_cols == 1), whose read has an erased_ref column, and that depend
  on erased columns (VC places, N3 abstains or places elsewhere, and the read has one) -- a column counts only when its base
  can vote (quality >= MIN_BASEQ), as it does for a fix's evidence(); each fix with its evidence() and whether AS (max), NM
  (min) or subs (min) alone picks the truth (a tie is no pick, as o2tags.pick); three worked reads.
"""
import bisect
import collections
import csv
import fractions
import importlib.util
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PAD = 1000                                                      # the caller's region: a copy's span and PAD bases each side
CHRM = "chrM"
FILTERS = (("PASS", "pass"), ("LowQual", "lowqual"), ("HomRef", "homref"), ("Multiallelic", "multiallelic"),
           ("RnaEdit", "rnaedit"), ("Dense", "dense"))


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_o2vc", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D, R, P, T = _sibling("o2vc_data"), _sibling("o2vc_rule"), _sibling("o2psv"), _sibling("o2tags")


def caller_region(chrom, start, end, length):
    """(longcallR's -r, the 0-based [lo, hi) where its calls are read) for a copy spanning [start, end) on a contig of
    `length`: PAD bases each side, START 1-based and >= 1 (0 panics), END exclusive and <= length + 1 (beyond it the pileup
    indexes past the contig and panics, util.rs:664, observed). Its first base is left out of [lo, hi): the caller fetches
    reads from START as if 0-based (o2vc_data.teacher_depth)."""
    s, e = max(1, start + 1 - PAD), min(end + PAD, length + 1)
    return f"{chrom}:{s}-{e}", (s, e - 1)


class Genome:
    """The genome, upper-cased (the FASTA is soft-masked; longcallR upper-cases it). fetch(chrom, lo, hi) serves [lo, hi)
    from windows read once -- hold() one per copy -- and reads a window PAD wider each side for a request outside them, so
    a read's columns across its introns are a slice, never a genome read (o2correct.collect's region())."""

    def __init__(self, path):
        import pysam
        self.fa = pysam.FastaFile(path)
        self.length = dict(zip(self.fa.references, self.fa.lengths))
        self.win = collections.defaultdict(list)

    def hold(self, chrom, lo, hi):
        lo = max(0, lo)
        self.win[chrom].append((lo, self.fa.fetch(chrom, lo, hi).upper()))

    def fetch(self, chrom, lo, hi):
        for s, seq in self.win[chrom]:
            if s <= lo and hi <= s + len(seq):
                return seq[lo - s:hi - s]
        self.hold(chrom, lo - PAD, hi + PAD)
        s, seq = self.win[chrom][-1]
        return seq[lo - s:hi - s]


def chrom_fasta(genome, chrom, ref_dir):
    """ref_dir/<chrom>.fa with its .fai: the whole contig as the genome holds it (longcallR upper-cases it), written once.
    longcallR reads its whole FASTA: 16 s and 3.4 GB per call on CHM13, 0.7 s and 0.5 GB on one contig, the same calls."""
    import pysam
    path = os.path.join(ref_dir, f"{chrom}.fa")
    if not os.path.exists(path + ".fai"):
        os.makedirs(ref_dir, exist_ok=True)
        seq = genome.fa.fetch(chrom)
        with open(path + ".tmp", "w") as fh:
            fh.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
        os.replace(path + ".tmp", path)
        pysam.faidx(path)
    return path


def check_ref(calls, genome):
    """Stop (ValueError) at the first record, in position order, whose REF is not the genome's base there, upper-cased."""
    for (chrom, p), call in sorted(calls.items()):
        base = genome.fetch(chrom, p, p + len(call.ref))
        if call.ref.upper() != base:
            raise ValueError(f"{chrom}:{p + 1} REF {call.ref} != genome {base}")


def separated(g):
    """The pairs (c, d), c < d, of known copies in g whose genotypes share no allele."""
    known = [(c, x) for c, x in sorted(g.items()) if x is not None]
    return {(c, d) for i, (c, x) in enumerate(known) for d, y in known[i + 1:] if not x & y}


def votes(copies, cols):
    """{copy: (columns voting for it, columns voting against it)} in o2vc_rule.assign()'s sense: a column separating c
    from d votes, at quality >= MIN_BASEQ, for the one whose genotype holds the read's base and against the other."""
    cs = sorted(copies)
    vf, va = {c: set() for c in cs}, {c: set() for c in cs}
    for n, col in enumerate(cols):
        if col["qual"] < R.MIN_BASEQ:
            continue
        for c, d in separated(col["g"]):
            hit = [x for x in (c, d) if col["base"] in col["g"][x]]
            if hit:
                vf[hit[0]].add(n)
                va[d if hit[0] == c else c].add(n)
    return {c: (len(vf[c]), len(va[c])) for c in cs}


def mismatches(a, b, off=0):
    """The offsets where the equal-length strings a and b differ, in order, found by halving (C-speed comparisons)."""
    if a == b:
        return []
    if len(a) <= 32:
        return [off + i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    h = len(a) // 2
    return mismatches(a[:h], b[:h], off) + mismatches(a[h:], b[h:], off + h)


def _alleles(g):
    return "." if g is None else "/".join(sorted(g))


def _key(k):
    return "|".join(str(x) for x in k)


class Copies:
    """A family's copies as the rule sees them: per copy its calls ({(chrom, pos): Call}, absent when it has no teacher),
    the [lo, hi) its calls cover and the teacher depth at the positions asked for (measure())."""

    def __init__(self, span, calls, win, teach_bam, genome):
        self.span, self.calls, self.win, self.teach_bam, self.genome = span, calls, win, teach_bam, genome
        self.depth = collections.defaultdict(dict)
        self.records = {c: sorted(p for _, p in v) for c, v in calls.items()}

    def called(self, c, p):
        return c in self.calls and self.win[c][0] <= p < self.win[c][1]

    def measure(self, c, positions):
        """Teacher depth at the positions of copy c (inside its region, without a record) not measured yet."""
        chrom = self.span[c][0]
        want = sorted(p for p in set(positions) if p not in self.depth[c] and self.called(c, p)
                      and (chrom, p) not in self.calls[c])
        if want:
            self.depth[c].update(D.teacher_depth(self.teach_bam, chrom, want, self.genome.fetch))

    def genotype(self, c, chrom, p, ref):
        """(copy c's alleles at p in the genome's strand, None when unknown; the depth used, None when not measured). ref: the
        genome's base at p, upper-cased. A position without a record must have been measured."""
        if not self.called(c, p) or ref not in "ACGT":
            return None, None
        call = self.calls[c].get((chrom, p))
        if call is not None:
            return R.genotype(call, call.dp, ref), call.dp
        d = self.depth[c][p]
        return (R.genotype(None, d, ref) if d <= D.MAX_DEPTH else None), d


def _interesting(col, cps):
    """Whether any copy's genotype can matter at the column: a record there, two copies' reference bases differing, or the
    read's base not every copy's reference base. Elsewhere every known genotype is the one reference base the read holds,
    so the column neither separates nor holds a foreign allele, known or not."""
    refs = set(col["ref"].values())
    return len(refs) > 1 or col["base"] not in refs or any(col["at"][c] in cps.calls.get(c, ()) for c in col["at"])


def hold_reads(f, reads, tied, cps):
    """Per tied read of family f: its calls.tsv fields, reads.tsv and columns.tsv rows, and its foreign columns for flags()."""
    picked = {}
    for q in reads:
        t = tied[q]
        picked[q] = [col for col in D.columns(t, cps.genome.fetch) if _interesting(col, cps)]
    need = collections.defaultdict(set)
    for q in reads:
        for col in picked[q]:
            for c, (_, p) in col["at"].items():
                need[c].add(p)
    for c in sorted(need):
        if c in cps.calls:
            cps.measure(c, need[c])
    out = {}
    for q in reads:
        t, cs = tied[q], sorted(tied[q].copies)
        vc_cols, n3_cols, rows, foreign = [], [], [], []
        for col in picked[q]:
            gv, gn, dep = {}, {}, {}
            for c in cs:
                flip, (chrom, p) = t.copies[c].flip, col["at"][c]
                ref = col["ref"][c]
                x, dep[c] = cps.genotype(c, chrom, p, ref.translate(R.COMP) if flip else ref)
                gv[c] = R.to_read_frame(x, flip)
                gn[c] = frozenset({ref}) if gv[c] is not None else None
            vc_cols.append(dict(base=col["base"], qual=col["qual"], g=gv))
            n3_cols.append(dict(base=col["base"], qual=col["qual"], g=gn))
            fk, sv, sn = R.foreign(vc_cols[-1]), separated(gv), separated(gn)
            if fk == "foreign":                                 # flags() takes every foreign allele, whatever the kind
                foreign.append((col["p"], frozenset(
                    (c, *col["at"][c], col["base"].translate(R.COMP) if t.copies[c].flip else col["base"]) for c in cs)))
            # the first that holds; N3 knows only what VC knows, so an erased pair's calls share an allele
            kind = "separating_call" if sv - sn else "separating_ref" if sv else "erased_ref" if sn else fk
            if not kind:
                continue
            rows.append([q, col["p"], col["base"], col["qual"], kind, ";".join(
                f"{c}|{col['at'][c][0]}|{col['at'][c][1]}|{col['ref'][c]}|{_alleles(gv[c])}|{'.' if dep[c] is None else dep[c]}"
                for c in cs)])
        vc, n3 = R.assign(cs, vc_cols), R.assign(cs, n3_cols)
        vv, vn = votes(cs, vc_cols), votes(cs, n3_cols)
        tags = []
        for c in cs:
            a = t.copies[c]
            subs = "" if a.NM is None else P.substitutions(a.cigar, a.NM)[0]
            indels = sum(n for op, n in a.cigar if op in (1, 2))
            tags.append([q, c, int(t.prim == c), a.AS, "" if a.NM is None else a.NM, "" if a.de is None else f"{a.de:.6g}",
                         "" if a.ms is None else a.ms, subs, indels, *vn[c], *vv[c]])
        out[q] = dict(vc=vc, n3=n3, reads=tags, columns=rows, foreign=foreign)
    return out


def teacher_side(owner, cps):
    """{copy: its teachers in a flag}: o2vc_rule.flags() over the teachers' own-copy foreign alleles -- a base of a
    teacher's primary that its copy's genotype lacks where the copy is known (foreign(): quality >= MIN_BASEQ, A>G in the
    transcript frame left out as editing). Reported, never counted: teachers cover every untied read. Depth is measured only
    at the no-record alleles of a flag made over every candidate: keeping the known ones leaves a subset, whose flags lie
    within those."""
    import pysam
    cand, pend = {}, collections.defaultdict(set)
    with pysam.AlignmentFile(cps.teach_bam) as bam:
        for rec in bam.fetch(until_eof=True):
            c = owner.get(rec.query_name)
            if c not in cps.calls or rec.query_sequence is None:
                continue
            seq, qual, flip = rec.query_sequence.upper(), rec.query_qualities, rec.is_reverse
            chrom, (lo, hi), recs = rec.reference_name, cps.win[c], cps.records[c]
            found, q, r = [], 0, rec.reference_start
            for op, n in rec.cigartuples:
                if op in D.ALIGN_OPS:
                    a, b = max(r, lo), min(r + n, hi)
                    if a < b:
                        ref = cps.genome.fetch(chrom, a, b)
                        hits = {a + i for i in mismatches(seq[q + a - r:q + b - r], ref)}
                        hits.update(recs[bisect.bisect_left(recs, a):bisect.bisect_left(recs, b)])
                        for p in sorted(hits):
                            base, call = seq[q + p - r], cps.calls[c].get((chrom, p))
                            g = R.genotype(call, call.dp, ref[p - a]) if call is not None else frozenset({ref[p - a]})
                            col = dict(base=base.translate(R.COMP) if flip else base,
                                       qual=D.MISSING_QUAL if qual is None else qual[q + p - r], g={c: R.to_read_frame(g, flip)})
                            if g is not None and R.foreign(col) == "foreign":
                                found.append((p, frozenset({(c, chrom, p, base)})))
                                if call is None:
                                    pend[c].add(p)
                    q, r = q + n, r + n
                elif op in (1, 4):
                    q += n
                elif op in (2, 3):
                    r += n
            if found:
                cand[rec.query_name] = found
    need = collections.defaultdict(set)
    for fl in R.flags(cand):
        for c, _, p, _ in fl["alleles"]:
            if p in pend[c]:
                need[c].add(p)
    for c in sorted(need):
        cps.measure(c, need[c])

    def known(k):
        c, _, p, _ = k
        return p not in pend[c] or R.MIN_DEPTH <= cps.depth[c].get(p, 0) <= D.MAX_DEPTH

    keep = {q: [(p, ks) for p, ks in v if all(known(k) for k in ks)] for q, v in cand.items()}
    out = collections.Counter()
    for fl in R.flags(keep):
        for q in fl["reads"]:
            out[owner[q]] += 1
    return out


def run_family(f, cs, reads, sc, ctx):
    """One family: teacher BAM, calls, tied reads held, flags, teacher side, the family's files. Returns its calls.tsv
    rows."""
    out, g, log = ctx["out_dir"], ctx["genome"], ctx["log"]
    span = {c: ctx["span"][c] for c in cs}
    teachers = {c: sc.teachers.get(c, set()) for c in cs}
    teach_bam = os.path.join(out, "teach", f"{f}.bam")
    n_teach = D.write_teacher_bam(ctx["bam_p"], span, teachers, teach_bam, counts=ctx["q0_raised"])
    g.win.clear()
    calls, win = {}, {}
    for c in cs:
        chrom, s, e = span[c]
        region, win[c] = caller_region(chrom, s, e, g.length[chrom])
        g.hold(chrom, win[c][0] - 1, win[c][1])
        if not teachers[c]:
            continue                                            # no teacher: no call, unknown everywhere
        vcf = os.path.join(out, "vcf", f"{c}.vcf")
        if os.path.exists(vcf):                                 # never re-called: the caller is not deterministic
            ctx["reused"].append(c)
            calls[c] = D.parse_vcf(vcf)
            check_ref(calls[c], g)
            log(f"{f}: {c}: reusing {vcf} ({len(calls[c])} records), the saved calls being the record")
            continue
        try:                                                    # called into vcf/.tmp, saved once checked
            got = D.run_longcallr(teach_bam, chrom_fasta(g, chrom, os.path.join(out, "ref")),
                                  os.path.join(out, "vcf", ".tmp", c), binary=ctx["longcallr"], threads=ctx["threads"],
                                  region=region)
        except RuntimeError as err:
            ctx["failed"][c] = str(err)
            log(f"{f}: {c} {region}: longcallR failed, the copy is unknown everywhere: {str(err).splitlines()[0]}")
            continue
        calls[c] = D.parse_vcf(got)
        check_ref(calls[c], g)
        os.replace(got, vcf)
        log(f"{f}: {c} {region}: {len(teachers[c])} teachers, {len(calls[c])} records")
    cps = Copies(span, calls, win, teach_bam, g)
    for c in sorted(calls):
        td = D.teacher_depth(teach_bam, span[c][0], [p for _, p in calls[c]], g.fetch)
        ctx["depth_agreement"][0] += sum(abs(call.dp - td[p]) <= 1 for (_, p), call in calls[c].items())
        ctx["depth_agreement"][1] += len(calls[c])
    held = hold_reads(f, reads, sc.tied, cps)
    flags = R.flags({q: held[q]["foreign"] for q in reads if held[q]["foreign"]})
    flag_of = collections.defaultdict(list)
    for fl in flags:
        for q in fl["reads"]:
            flag_of[q].append(f"{f}#{fl['id']}")
    owner = {q: c for c in cs for q in teachers[c]}
    teacher_flags = teacher_side(owner, cps)
    rows = []
    for q in reads:
        t, h = sc.tied[q], held[q]
        (vc, vc_set, vc_n), (n3, n3_set, _) = h["vc"], h["n3"]
        rows.append([q, f, ",".join(sorted(t.copies)), t.prim or "", t.mapq, t.margin, int(t.prim not in t.copies),
                     "" if flag_of[q] or vc is None else vc, ",".join(sorted(vc_set)), vc_n, n3 or "", ",".join(sorted(n3_set)),
                     ",".join(flag_of[q])])
    _write(os.path.join(out, f"{f}.reads.tsv"), "read copy primary AS NM de ms subs indels ref_for ref_against call_for "
           "call_against", [r for q in reads for r in held[q]["reads"]])
    _write(os.path.join(out, f"{f}.columns.tsv"), "read p base qual kind copies", [r for q in reads for r in held[q]["columns"]])
    _write(os.path.join(out, f"{f}.flags.tsv"), "flag n_reads reads alleles",
           [[f"{f}#{fl['id']}", len(fl["reads"]), ",".join(fl["reads"]), ";".join(_key(k) for k in fl["alleles"])] for fl in flags])
    filt = {c: collections.Counter(call.filter for call in calls.get(c, {}).values()) for c in cs}
    _write(os.path.join(out, f"{f}.teachers.tsv"), "copy teachers " + "records " + " ".join(k for _, k in FILTERS)
           + " teacher_flag_reads", [[c, len(teachers[c]), len(calls.get(c, {})), *(filt[c][v] for v, _ in FILTERS),
                                      teacher_flags[c]] for c in cs])
    placed = sum(1 for r in rows if r[7])
    log(f"{f}: {len(reads)} tied reads ({n_teach} teachers written); VC placed {placed}, N3 placed "
        f"{sum(1 for r in rows if r[10])}, {len(flags)} flags on {sum(1 for r in rows if r[12])} reads")
    return rows


def _write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header.split()) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def cmd_run(bam_p, copies_p, genome_p, out_dir, families=None, exclude_chrM=False, longcallr=D.LONGCALLR, threads=2):
    t0 = time.time()

    def log(msg):
        print(f"[o2vc] {time.time() - t0:7.1f}s {msg}", file=sys.stderr, flush=True)

    os.makedirs(out_dir, exist_ok=True)
    cps = P.load_copies(copies_p)
    span = {c: (ch, s, e) for ch, v in cps.items() for c, _, _, (s, e) in v}
    fam = collections.defaultdict(list)
    for c in sorted(span):
        fam[c.rsplit(".", 1)[0]].append(c)
    asked = sorted(fam) if families is None else sorted(set(families) & set(fam))
    if families is not None and set(families) - set(fam):
        log(f"not in {copies_p}: {', '.join(sorted(set(families) - set(fam)))}")
    chrm = [f for f in asked if exclude_chrM and any(span[c][0] == CHRM for c in fam[f])]
    if chrm:
        log(f"left out (a copy on {CHRM}): {', '.join(chrm)}")
    run_f = [f for f in asked if f not in chrm]
    sc = D.scan(bam_p, copies_p, None if families is None and not chrm else set(run_f))
    by_fam = collections.defaultdict(list)
    for q in sorted(sc.tied):
        by_fam[sc.tied[q].family].append(q)
    todo = [f for f in run_f if by_fam[f]]
    no_qual = sum(1 for q in sorted(sc.tied) if all(x == D.MISSING_QUAL for x in sc.tied[q].qual))
    log(f"{len(todo)} of {len(run_f)} families have tied reads; {no_qual} tied reads without QUAL (each base votes)")
    try:
        version = D.longcallr_version(longcallr)
    except (OSError, subprocess.CalledProcessError) as err:
        version = f"unavailable: {err}"
    ctx = dict(out_dir=out_dir, genome=Genome(genome_p), log=log, span=span, bam_p=bam_p, longcallr=longcallr,
               threads=threads, depth_agreement=[0, 0], failed={}, reused=[], q0_raised=collections.Counter())
    rows = []
    for f in todo:
        rows += run_family(f, fam[f], by_fam[f], sc, ctx)
    _write(os.path.join(out_dir, "calls.tsv"), "read family candidates minimap2 mapq margin outside vc vc_set vc_cols n3 "
           "n3_set flag", rows)
    run = dict(longcallr=version, binary=longcallr, preset=D.PRESET, min_mapq=D.MIN_MAPQ, threads=threads, families=todo,
               excluded_chrM=chrm, depth_agreement=ctx["depth_agreement"], failed_copies=ctx["failed"],
               reused_vcfs=ctx["reused"], q0_raised=dict(reads=ctx["q0_raised"]["reads"], bases=ctx["q0_raised"]["bases"]),
               tied_without_qual=no_qual)
    with open(os.path.join(out_dir, "run.json"), "w") as fh:
        json.dump(run, fh, indent=1)
    try:
        os.rmdir(os.path.join(out_dir, "vcf", ".tmp"))          # empty unless a failed call left a partial VCF there
    except OSError:
        pass
    log(f"-> {out_dir}: {len(rows)} tied reads in {len(todo)} families; VC placed {sum(1 for r in rows if r[7])}, "
        f"flagged {sum(1 for r in rows if r[12])}; depth agreement {ctx['depth_agreement'][0]} of {ctx['depth_agreement'][1]} "
        f"records; {len(ctx['failed'])} copies failed, {len(ctx['reused'])} VCFs reused, {ctx['q0_raised']['bases']} teacher "
        f"bases of quality 0 written as 1")
    return 0


# ---- score: the arms, the bars and the case files with truth ----------------------------------------------------------------

ARMS = ("N1", "N2", "N3", "N4", "VC")
ARM_NAME = dict(N1="minimap2", N2="minimap2 at MAPQ >= 1", N3="reference only", N4="corrector", VC="variant calls")
MODES, PANELS = ("donor", "none"), ("A", "C", "E")              # world / units modes and the D, Cg panels are deferred
HAS_COPY = ("unique", "unique_elsewhere")                       # classes whose truth names one copy of the reference
SCORED = HAS_COPY + ("elsewhere", "no_counterpart")             # a definite truth: a copy, or none (every placement is wrong)
PRESENT = HAS_COPY + ("own_tie",)                               # classes whose copy the reference holds (K2)
NO_TRUTH = "no_truth_row"                                       # a tied read the truth tables do not hold
OFF = "(off every copy)"                                        # N1's placement when minimap2's primary lies on no copy
LIMIT = fractions.Fraction(2, 100)                              # the 2% of A1, E1, K1 and K2, compared exactly
MIN_FIXES = 10                                                  # A2
EVIDENCE = ("sample call", "reference PSV corrected by the sample", "reference PSV")
TAGS = ("AS", "NM", "subs")                                     # o2tags.pick: AS the highest, NM and subs the lowest
BAR_NAME = dict(A1="wrong/placed", A2="VC fixes", A3="VC fixes vs N3 fixes", K1="VC placed/own_tie reads",
                K2="flagged/reads with a present copy", C1="flags with >= 2 no_counterpart reads", E1="wrong/placed",
                E2="fixes vs breaks")


def evidence(kinds):
    """What decided a placement, from the kinds of its read's columns in columns.tsv: "sample call" when any is separating_call
    (VC separates a pair of copies the reference does not); else "reference PSV corrected by the sample" when any is erased_ref
    (a reference PSV the sample's calls do not carry was taken out of the vote); else "reference PSV"."""
    if "separating_call" in kinds:
        return EVIDENCE[0]
    return EVIDENCE[1] if "erased_ref" in kinds else EVIDENCE[2]


def _tsv(path):
    """(header, rows as dicts) of a tab-separated file."""
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return rd.fieldnames, list(rd)


def need_columns(path, *cols):
    """Stop (ValueError) when the tab-separated file `path` lacks a named column (a truth table of another mode, say)."""
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
    missing = [c for c in cols if c not in head]
    if missing:
        raise ValueError(f"{path}: no {' / '.join(missing)} column (has {', '.join(head)})")


def tag_picks(rows, truth):
    """{AS, NM, subs: whether that statistic alone picks `truth`} over a read's tied copies (reads.tsv rows {copy: row}): AS the
    highest, NM and subs the lowest, a tie or a statistic not measured on some copy no pick (o2tags.pick)."""
    out = {}
    for key in TAGS:
        v = {c: row[key] for c, row in rows.items()}
        out[key] = "" not in v.values() and T.pick({c: {key: float(x)} for c, x in v.items()}, key) == truth
    return out


def outcome(r):
    """VC's fate on a read, against minimap2's and the truth: flagged, abstain, placed (a read that is not scored: no truth to
    judge it against), or for a scored read fix, break, both right, both wrong. A read with no truth copy (elsewhere,
    no_counterpart) has `truth` None, so no placement is right and N1 is wrong: a placed one is both wrong."""
    vc, n1, t = r["place"]["VC"], r["place"]["N1"], r["truth"]
    if r["flag"]:
        return "flagged"
    if not vc:
        return "abstain"
    if not r["scored"]:
        return "placed"
    return ("both right" if n1 == t else "fix") if vc == t else ("break" if n1 == t else "both wrong")


def load_reads(out_dir, truth, n4):
    """calls.tsv as records, sorted by (family, read). place: each arm's copy, "" when it abstains (a flagged read is not
    placed, so a row holding both is inconsistent); cls / truth from the truth {read: (class, copy)}, None in none mode.
    scored: the class has a definite truth (unique, unique_elsewhere, elsewhere, no_counterpart); decidable: it is a class with
    a truth copy and that copy is one of the tied copies (the subset earlier registers scored)."""
    recs, seen = [], set()
    for r in _tsv(os.path.join(out_dir, "calls.tsv"))[1]:
        q, cands, mm = r["read"], r["candidates"].split(","), r["minimap2"]
        if q in seen:
            raise ValueError(f"{q}: twice in calls.tsv")
        seen.add(q)
        if (r["outside"] == "1") != (mm not in cands):
            raise ValueError(f"{q}: outside = {r['outside']}, but minimap2 {mm or '(no copy)'} and candidates {r['candidates']}")
        if r["flag"] and r["vc"]:
            raise ValueError(f"{q}: flagged ({r['flag']}) and placed ({r['vc']}); a flagged read is not placed")
        n1 = mm or OFF
        place = dict(N1=n1, N2=n1 if int(r["mapq"]) >= 1 else "", N3=r["n3"])
        if n4 is not None:
            place["N4"] = n4.get(q, n1)
        place["VC"] = r["vc"]
        cls, tc = (None, None) if truth is None else truth.get(q, (NO_TRUTH, None))
        if cls in HAS_COPY and tc is None:
            raise ValueError(f"{q}: class {cls}, but the truth table names no truth copy")
        recs.append(dict(read=q, family=r["family"], cands=cands, mm=mm, mapq=int(r["mapq"]), outside=r["outside"] == "1",
                         place=place, flag=r["flag"], vc_cols=int(r["vc_cols"]), cls=cls, truth=tc, scored=cls in SCORED,
                         decidable=cls in HAS_COPY and tc in cands, n4_held=n4 is not None and q in n4))
    recs.sort(key=lambda r: (r["family"], r["read"]))
    for r in recs:
        r["outcome"] = outcome(r)
    return recs


# a family's reads.tsv (header, rows, rows by read) and columns.tsv (header, rows by read)
Files = collections.namedtuple("Files", "rh rrows reads ch cols")


def by_family(recs):
    """{family: its records}, families in sorted order."""
    out = collections.defaultdict(list)
    for r in recs:
        out[r["family"]].append(r)
    return dict(sorted(out.items()))


def load_files(out_dir, fams):
    """Each record of `fams` ({family: records}) gets its reads.tsv rows (tags: {copy: row}) and the kinds of its columns.tsv
    rows whose base can vote (quality >= R.MIN_BASEQ: the rest decide nothing, so they set neither the evidence of a fix nor
    the erased counts). Returns {family: Files}; a read in those files that calls.tsv lacks is an error."""
    files = {}
    for f in fams:
        (rh, rrows), (ch, crows) = (_tsv(os.path.join(out_dir, f"{f}.{x}.tsv")) for x in ("reads", "columns"))
        reads, cols = collections.defaultdict(list), collections.defaultdict(list)
        for x in rrows:
            reads[x["read"]].append(x)
        for x in crows:
            cols[x["read"]].append(x)
        extra = sorted((set(reads) | set(cols)) - {r["read"] for r in fams[f]})
        if extra:
            raise ValueError(f"{f}.reads.tsv / {f}.columns.tsv: {extra[0]} is not in calls.tsv")
        bare = [r["read"] for r in fams[f] if not reads[r["read"]]]
        if bare:
            raise ValueError(f"{f}.reads.tsv: {bare[0]} (a tied read of calls.tsv) has no row")
        files[f] = Files(rh, rrows, reads, ch, cols)
        for r in fams[f]:
            r["tags"] = {x["copy"]: x for x in reads[r["read"]]}
            r["kinds"] = frozenset(x["kind"] for x in cols[r["read"]] if int(x["qual"]) >= R.MIN_BASEQ)
    return files


def arm_counts(recs, arms):
    """{arm: placed, right, wrong, abstained, fixes, breaks, both_wrong} over `recs` (scored reads, the same for every arm).
    A placement is right only when it is the read's truth copy -- which a read with a truth copy outside the tied copies, an
    elsewhere or a no_counterpart read never gets from a tied copy. fixes: placed right where N1 is wrong; breaks: placed wrong
    where N1 is right (minimap2's primary may be the truth copy outside the tied set); both_wrong: placed wrong where N1 is
    wrong too (the same copy or not)."""
    out = {}
    for a in arms:
        c = dict(placed=0, right=0, wrong=0, abstained=0, fixes=0, breaks=0, both_wrong=0)
        for r in recs:
            t, x = r["truth"], r["place"][a]
            if not x:
                c["abstained"] += 1
                continue
            c["placed"] += 1
            if x == t:
                c["right"] += 1
                c["fixes"] += r["place"]["N1"] != t
            else:
                c["wrong"] += 1
                c["breaks" if r["place"]["N1"] == t else "both_wrong"] += 1
        out[a] = c
    return out


def flag_order(f):
    fam, _, n = f.rpartition("#")
    return fam, len(n), n


def flag_table(recs, donor):
    """[{id, family, reads, classes}] of the flags the reads carry, in flag order; classes: the truth classes of a flag's reads."""
    by = collections.defaultdict(list)
    for r in recs:
        for f in r["flag"].split(","):
            if f:
                by[f].append(r)
    out = []
    for f in sorted(by, key=flag_order):
        e = dict(id=f, family=f.rpartition("#")[0], reads=len(by[f]))
        if donor:
            e["classes"] = dict(sorted(collections.Counter(r["cls"] for r in by[f]).items()))
        out.append(e)
    return out


def vc_groups(recs, donor):
    """VC's placements, and those resting on a single column (vc_cols == 1), whose read has an erased_ref column, and that depend
    on erased columns: VC places the read, N3 abstains or places it elsewhere, and the read has one. Each with, for a donor
    truth, how many are on scored reads and right / wrong there (a placement is right only when it is the truth copy).
    not_scored: the placements on reads that are not scored (own_tie, poor_match, unmapped, no truth row), which the arm table
    cannot see, by class."""
    placed = [r for r in recs if r["place"]["VC"]]

    def group(pick):
        rows = [r for r in placed if pick(r)]
        g = dict(placements=len(rows))
        if donor:
            sc = [r for r in rows if r["scored"]]
            right = sum(1 for r in sc if r["place"]["VC"] == r["truth"])
            g.update(scored=len(sc), right=right, wrong=len(sc) - right)
        return g

    def erased(r):
        return "erased_ref" in r["kinds"]
    out = dict(placements=len(placed), single_column=group(lambda r: r["vc_cols"] == 1), erased_ref=group(erased),
               erased_dependent=group(lambda r: erased(r) and r["place"]["N3"] != r["place"]["VC"]))
    if donor:
        ns = collections.Counter(r["cls"] for r in placed if not r["scored"])
        out["not_scored"] = dict(placements=sum(ns.values()), classes=dict(sorted(ns.items())))
    return out


def fix_rows(recs):
    """[(record, evidence, tag picks)] of the fixes (scored reads VC places right where minimap2 is wrong), by read name."""
    return [(r, evidence(r["kinds"]), tag_picks(r["tags"], r["truth"])) for r in recs if r["outcome"] == "fix"]


def precision(arms):
    """{arm: wrong, placed, rate}: wrong / placed always with both counts (rate None when nothing is placed)."""
    return {a: dict(wrong=c["wrong"], placed=c["placed"], rate=c["wrong"] / c["placed"] if c["placed"] else None)
            for a, c in arms.items()}


def block(recs, arms, donor):
    """The score of one set of tied reads (a family, or the pooled families). Donor: the arms are over the SCORED reads, for
    every arm the same (class unique, unique_elsewhere, elsewhere or no_counterpart: a placement is right only when it is the
    truth copy, wrong otherwise); `decidable` repeats them over the subset whose truth copy is one of the tied copies."""
    flags = flag_table(recs, donor)
    out = dict(tied=len(recs))
    if not donor:
        n_flag = sum(1 for r in recs if r["flag"])
        placed = {a: sum(1 for r in recs if r["place"][a]) for a in arms}
        out["arms"] = {a: dict(placed=placed[a], abstained=len(recs) - placed[a], flagged=n_flag) for a in arms}
        out.update(outside=dict(tied=sum(1 for r in recs if r["outside"])), vc=vc_groups(recs, False), flags=flags)
    else:
        scored = [r for r in recs if r["scored"]]
        dec = [r for r in scored if r["decidable"]]
        outs = [r for r in scored if r["outside"]]
        out.update(scored=len(scored),
                   truth_outside_candidates=sum(1 for r in scored if r["cls"] in HAS_COPY and not r["decidable"]))
        out["arms"] = arm_counts(scored, arms)
        out["precision"] = precision(out["arms"])
        dec_arms = arm_counts(dec, arms)
        out["decidable"] = dict(reads=len(dec), arms=dec_arms, precision=precision(dec_arms))
        out["outside"] = dict(tied=sum(1 for r in recs if r["outside"]), scored=len(outs), arms=arm_counts(outs, arms))
        cls = {}
        for r in recs:
            c = cls.setdefault(r["cls"], dict(reads=0, placed=0, flagged=0))
            c["reads"] += 1
            c["placed"] += bool(r["place"]["VC"])
            c["flagged"] += bool(r["flag"])
        out["classes"] = dict(sorted(cls.items()))
        tie, here = cls.get("own_tie", {}), [cls[c] for c in PRESENT if c in cls]
        out["K1"] = dict(reads=tie.get("reads", 0), placed=tie.get("placed", 0))
        out["K2"] = dict(reads=sum(c["reads"] for c in here), flagged=sum(c["flagged"] for c in here))
        right = [r for r in scored if r["place"]["N1"] == r["truth"]]
        out["D2"] = dict(reads=len(right), moved=sum(1 for r in right if r["place"]["VC"] and r["place"]["VC"] != r["truth"]),
                         flagged=sum(1 for r in right if r["flag"]))
        out["C1"] = dict(flags=len(flags), counted=sum(1 for f in flags if f["classes"].get("no_counterpart", 0) >= 2))
        out["vc"] = vc_groups(recs, True)
        fixes = fix_rows(scored)
        ev = collections.Counter(e for _, e, _ in fixes)
        out["fixes"] = dict(reads=len(fixes), evidence={e: ev[e] for e in EVIDENCE},
                            tags_alone=dict({k: sum(p[k] for _, _, p in fixes) for k in TAGS},
                                            none=sum(1 for _, _, p in fixes if not any(p.values()))))
        out["flags"] = flags
    if "N4" in arms:                                            # how many tied reads the corrector holds; N1 stands in for the rest
        out["n4_held"] = sum(1 for r in recs if r["n4_held"])
    return out


def ratio_bar(num, den, why):
    """A bar of the form num / den <= 2%. With nothing to measure (den 0) it is not evaluable: pass None, and `why` says why."""
    b = {"value": num / den if den else None, "num": num, "den": den, "limit": float(LIMIT), "op": "<=",
         "pass": fractions.Fraction(num, den) <= LIMIT if den else None}
    return b if den else dict(b, reason=why)


def count_bar(value, op, limit, why=None):
    """A bar of the form value >= limit (or value > limit). `why`, when given, says there is nothing to measure: the bar is not
    evaluable and passes None."""
    b = {"value": value, "limit": limit, "op": op, "pass": None if why else value >= limit if op == ">=" else value > limit}
    return dict(b, reason=why) if why else b


def panel_bars(panel, pooled, focus):
    """The bars of a panel from the pooled score and the focus family's (A1-A3 are the focus family's alone), all over the
    scored reads. A bar with nothing to measure is not evaluable: a zero denominator, no no_counterpart tied read for C1, no
    scored read on the focus family for A2 and A3, no VC placement of a scored read for E2."""
    if panel == "A":
        a, n3 = focus["arms"]["VC"], focus["arms"]["N3"]
        bare = None if focus["scored"] else "the focus family has no scored read"
        return {"A1": ratio_bar(a["wrong"], a["placed"], "VC placed no scored read of the focus family"),
                "A2": count_bar(a["fixes"], ">=", MIN_FIXES, bare), "A3": count_bar(a["fixes"], ">", n3["fixes"], bare),
                "K1": ratio_bar(pooled["K1"]["placed"], pooled["K1"]["reads"], "no own_tie read"),
                "K2": ratio_bar(pooled["K2"]["flagged"], pooled["K2"]["reads"], "no read with a present copy")}
    if panel == "C":
        bare = None if pooled["classes"].get("no_counterpart", {}).get("reads") else "no no_counterpart tied read in the pool"
        return {"C1": count_bar(pooled["C1"]["counted"], ">=", 1, bare)}
    if panel == "E":
        a = pooled["arms"]["VC"]
        none_placed = "VC placed no scored read"
        return {"E1": ratio_bar(a["wrong"], a["placed"], none_placed),
                "E2": count_bar(a["fixes"], ">", a["breaks"], None if a["placed"] else none_placed),
                "K2": ratio_bar(pooled["K2"]["flagged"], pooled["K2"]["reads"], "no read with a present copy")}
    return {}


def verdict_of(bars):
    """FAIL when any bar failed on evidence; else UNDERPOWERED when any bar had nothing to measure; else PASS."""
    passes = [b["pass"] for b in bars.values()]
    return "FAIL" if any(p is False for p in passes) else "UNDERPOWERED" if any(p is None for p in passes) else "PASS"


def _ratio(num, den):
    return f"{num}/{den} = {num / den:.4f}" if den else f"{num}/{den} = n/a"


def _group(g, of):
    """"k of n" for a vc_groups() group, with the right / wrong of its placements on scored reads when there is a truth."""
    return f"{g['placements']} of {of}" + (f" (scored {g['scored']}: {g['right']} right, {g['wrong']} wrong)"
                                           if "scored" in g else "")


def arm_table(arms, donor):
    if not donor:
        return ["| arm | placed | abstained | flagged |", "|---|---|---|---|"] + [
            f"| {a} {ARM_NAME[a]} | {c['placed']} | {c['abstained']} | {c['flagged']} |" for a, c in arms.items()]
    lines = ["| arm | placed | right | wrong | abstained | wrong / placed | fixes | breaks | both wrong |",
             "|---|---|---|---|---|---|---|---|---|"]
    for a, c in arms.items():
        cmp = ["-"] * 3 if a == "N1" else [c["fixes"], c["breaks"], c["both_wrong"]]
        lines.append("| " + " | ".join(str(x) for x in [f"{a} {ARM_NAME[a]}", c["placed"], c["right"], c["wrong"], c["abstained"],
                                                        _ratio(c["wrong"], c["placed"]), *cmp]) + " |")
    return lines


def _fenced(head, rows):
    """A fenced block of tab-separated rows under their header."""
    return ["```", "\t".join(head)] + ["\t".join(x[k] for k in head) for x in rows] + ["```", ""]


def worked(label, r, files, donor):
    """A worked read: the decision, the truth, its reads.tsv rows and its columns.tsv rows, verbatim."""
    if r is None:
        return [f"### {label}: none", "", "No such read in this family.", ""]
    p = r["place"]
    vc = f"VC {p['VC']} (vc_cols {r['vc_cols']})" if p["VC"] else "VC abstains"
    if r["flag"]:
        vc += f" (flagged {r['flag']})"
    lines = [f"### {label}: {r['read']}", "", f"minimap2 {r['mm'] or OFF} (MAPQ {r['mapq']}); N3 {p['N3'] or 'abstains'}; {vc}"]
    if donor:
        truth = f"truth {r['cls']}" + (f" {r['truth']}" if r["truth"] else "")
        lines.append(f"{truth}; outcome {r['outcome']}" + (f"; evidence {evidence(r['kinds'])}" if r["outcome"] == "fix" else ""))
    lines += ["", "reads.tsv rows:", ""] + _fenced(files.rh, files.reads[r["read"]])
    cols = files.cols[r["read"]]
    return lines + (["columns.tsv rows:", ""] + _fenced(files.ch, cols) if cols else ["columns.tsv rows: none.", ""])


LEGEND = ("Every arm is counted on the same scored reads. A placement is right only when it is the read's truth copy; every other "
          "placement on a scored read is wrong -- a truth copy outside the tied copies, an elsewhere read, a no_counterpart read. "
          "Fixes: placed right where N1 is wrong. Breaks: placed wrong where N1 is right (N1's primary can be the truth copy "
          "outside the tied copies). Both wrong: placed wrong where N1 is wrong too, on the same copy or not (N1's own cells are "
          "against itself). N1 always places: a primary on no copy is a wrong placement. Wrong / placed is each arm's precision on "
          "its own placements, so it is not comparable between arms that place different reads; compare arms by fixes and breaks.")


def head_lines(b, donor, opts):
    """The opening of a summary: what was scored and with which options, the arm table, the primaries outside the tied copies in
    a row of their own, and for a donor truth the decidable subset (earlier registers' definition)."""
    o = b["outside"]
    n4 = []
    if "n4_held" in b:
        n4 = [f"N4: the corrector holds {b['n4_held']} of {b['tied']} tied reads of this family; the others are N1's.", ""]
    if not donor:
        return [f"No truth was read. {b['tied']} tied reads; {o['tied']} have a minimap2 primary outside the tied copies.",
                f"Score options: {opts}", "", "## Arms on the tied reads", ""] + arm_table(b["arms"], False) + [
                "", "Flagged is a property of the read (a possible missing copy), the same for every arm; VC does not place "
                "a flagged read.", ""] + n4
    cls = b["classes"]
    every = ", ".join(f"{k} {v['reads']}" for k, v in cls.items())
    scored = ", ".join(f"{k} {v['reads']}" for k, v in cls.items() if k in SCORED) or "none"
    rest = ", ".join(f"{k} {v['reads']}" + (" (K1 only)" if k == "own_tie" else "") for k, v in cls.items() if k not in SCORED)
    L = [f"Scored against the donor's truth. {b['tied']} tied reads ({every}).",
         f"Scored: {b['scored']} ({scored}). Not scored: {rest or 'none'}.",
         f"Scored reads whose truth copy is outside the tied copies: {b['truth_outside_candidates']} (VC cannot be right on them).",
         f"Score options: {opts}", "", "## Arms on the scored reads", ""] + arm_table(b["arms"], True) + ["", LEGEND, ""] + n4
    L.append(f"Reads whose minimap2 primary lies outside the tied copies: {o['tied']} of {b['tied']} tied reads, "
             f"{o['scored']} scored.")
    L += ([""] + arm_table(o["arms"], True) if o["scored"] else [])
    d = b["decidable"]
    return L + ["", f"Decidable subset (the truth copy is one of the tied copies; the definition of earlier registers): "
                f"{d['reads']} reads.", ""] + arm_table(d["arms"], True) + [""]


def vc_lines(b, donor):
    v = b["vc"]
    L = ["## VC placements", "", f"VC placed {v['placements']} of {b['tied']} tied reads."]
    if donor:
        ns = v["not_scored"]
        by_class = ", ".join(f"{k} {n}" for k, n in ns["classes"].items()) or "none"
        L.append("VC placements on reads that are not scored, which the arm table cannot see: "
                 f"{ns['placements']} ({by_class}).")
    return L + [f"- resting on a single column (`vc_cols == 1`): {_group(v['single_column'], v['placements'])}",
                f"- whose read has an `erased_ref` column: {_group(v['erased_ref'], v['placements'])}",
                "- depending on erased columns (VC places, N3 abstains or places elsewhere, and the read has an `erased_ref` "
                f"column): {_group(v['erased_dependent'], v['placements'])}", ""]


def fix_lines(b, recs):
    fx = b["fixes"]
    t = fx["tags_alone"]
    ev = ", ".join(f"{e} {fx['evidence'][e]}" for e in EVIDENCE)
    L = ["## Fixes", "", f"Fixes (VC right, minimap2 wrong): {fx['reads']}. Evidence: {ev}. Alone, AS (max) picks the truth in "
         f"{t['AS']}, NM (min) in {t['NM']}, subs (min) in {t['subs']}; none of them in {t['none']} (a tie is no pick).", ""]
    rows = fix_rows(recs)
    if rows:
        L += ["| read | truth | minimap2 | evidence | AS alone | NM alone | subs alone |", "|---|---|---|---|---|---|---|"]
        for r, e, pk in rows:
            picks = " | ".join("yes" if pk[k] else "no" for k in TAGS)
            L.append(f"| {r['read']} | {r['truth']} | {r['mm'] or OFF} | {e} | {picks} |")
        L.append("")
    return L


def calibration_lines(b):
    k1, k2, d2, c1 = b["K1"], b["K2"], b["D2"], b["C1"]
    L = ["## Calibration of this family", "",
         f"- K1: VC places {k1['placed']} of {k1['reads']} own_tie reads.",
         f"- K2: {k2['flagged']} of {k2['reads']} reads with a present copy are flagged.",
         f"- D2: of {d2['reads']} scored reads N1 places right, VC moves {d2['moved']} and flags {d2['flagged']}.",
         f"- C1: {c1['counted']} of {c1['flags']} flags hold >= 2 no_counterpart reads."]
    for x in b["flags"]:
        L.append(f"- flag {x['id']}: reads {x['reads']} (" + ", ".join(f"{k} {n}" for k, n in x["classes"].items()) + ")")
    return L + [""]


def summary_text(f, recs, b, files, donor, opts):
    """<family>.summary.md for family f from its records and its block; opts: the score options, as options_text() writes them."""
    L = [f"# {f}: O2 from the sample's variant calls", ""] + head_lines(b, donor, opts) + vc_lines(b, donor)
    if donor:
        L += fix_lines(b, recs) + calibration_lines(b)

    def first(outcome_):                                                    # recs are sorted by read name
        return next((r for r in recs if r["outcome"] == outcome_), None)
    L += ["## Worked reads", ""]
    L += worked("First fix", first("fix"), files, donor) if donor else worked("First VC placement", first("placed"), files, donor)
    L += worked("First abstention", first("abstain"), files, donor)
    L += worked("First flagged read", first("flagged"), files, donor)
    return "\n".join(L)


def options_text(mode, truth_ps, n4_p, panel, focus, exclude_p):
    """The score options on one line (a path as it was given, "-" for an option not set), for the summary and the case file."""
    return (f"mode={mode} truth={','.join(truth_ps) or '-'} n4={n4_p or '-'} panel={panel or '-'} focus={focus or '-'} "
            f"exclude={exclude_p or '-'}")


def write_family(out_dir, f, recs, b, files, donor, opts):
    """<family>.case.tsv (reads.tsv's rows, in its order, plus the decision, the truth, the outcome and the score options -- a
    column, so that the file stays a plain table) and <family>.summary.md."""
    by = {r["read"]: r for r in recs}
    rows = []
    for x in files.rrows:
        r = by[x["read"]]
        rows.append([x[k] for k in files.rh] + [r["place"]["VC"], r["place"]["N3"], r["flag"], (r["cls"] or "") if donor else "",
                                                r["truth"] or "", r["outcome"], opts])
    head = files.rh + ["vc", "n3", "flag", "truth_class", "truth", "vc_outcome", "score_options"]
    _write(os.path.join(out_dir, f"{f}.case.tsv"), " ".join(head), rows)
    with open(os.path.join(out_dir, f"{f}.summary.md"), "w") as fh:
        fh.write(summary_text(f, recs, b, files, donor, opts) + "\n")


def say(msg):
    print(msg, flush=True)


def report(s):
    """The printed lines of a score, as o2correct.cmd_score prints its own."""
    donor = s["mode"] == "donor"
    what = s["mode"] + (f", panel {s['panel']}" if s["panel"] else "") + (f", focus {s['focus']}" if s["focus"] else "")
    excluded = f", excluded {len(s['excluded'])} ({', '.join(s['excluded'])})" if s["excluded"] else ""
    say(f"[o2vc] score {s['out_dir']} ({what}): families scored {len(s['families'])}{excluded}; {s['tied']} tied reads"
        + (f", {s['scored']} scored" if donor else ""))
    for a, c in s["arms"].items():
        if not donor:
            say(f"    {a}: placed {c['placed']}, abstained {c['abstained']}, flagged {c['flagged']}")
            continue
        line = (f"    {a}: placed {c['placed']}, right {c['right']}, wrong {c['wrong']}, abstained {c['abstained']}; "
                f"wrong/placed {_ratio(c['wrong'], c['placed'])}")
        say(line + ("" if a == "N1" else f"; fixes {c['fixes']}, breaks {c['breaks']}, both wrong {c['both_wrong']}"))
    if "n4_held" in s:
        say(f"    N4: the corrector holds {s['n4_held']} of {s['tied']} tied reads; the others are N1's")
    if donor:
        d = s["decidable"]
        dv, dn = d["arms"]["VC"], d["arms"]["N1"]
        say(f"    decidable subset ({d['reads']} reads): VC placed {dv['placed']}, right {dv['right']}, wrong {dv['wrong']}; "
            f"N1 right {dn['right']}, wrong {dn['wrong']}")
    v = s["vc"]
    say(f"    VC placements {v['placements']}: single column {v['single_column']['placements']}, with an erased_ref column "
        f"{v['erased_ref']['placements']}, depending on erased columns {v['erased_dependent']['placements']}"
        + (f"; on reads that are not scored {v['not_scored']['placements']}" if donor else ""))
    if donor:
        o, vo = s["outside"], s["outside"]["arms"]["VC"]
        say(f"    primary outside the tied copies: {o['tied']} tied reads, {o['scored']} scored; N1 right "
            f"{o['arms']['N1']['right']}, VC places {vo['placed']} and is right on {vo['right']}")
        say(f"    K1: VC places {s['K1']['placed']} of {s['K1']['reads']} own_tie reads")
        say(f"    K2: {s['K2']['flagged']} of {s['K2']['reads']} reads with a present copy are flagged")
        say(f"    D2: {s['D2']['reads']} scored reads N1 places right; VC moves {s['D2']['moved']}, flags {s['D2']['flagged']}")
        say(f"    C1: {s['C1']['counted']} of {s['C1']['flags']} flags hold >= 2 no_counterpart reads")
    for name, b in s["bars"].items():
        verdict = "PASS" if b["pass"] else "FAIL" if b["pass"] is False else f"not evaluable ({b['reason']})"
        measured = _ratio(b["num"], b["den"]) if "num" in b else str(b["value"])
        say(f"    bar {name}: {BAR_NAME[name]} {measured} {b['op']} {b['limit']}: {verdict}")
    if s["panel"]:
        say(f"    verdict {s['verdict']}")


def cmd_score(out_dir, truth_ps, mode, n4_p=None, panel=None, focus=None, exclude_p=None):
    if mode not in MODES:
        raise ValueError(f"--mode {mode}: only donor and none are built (world and units are deferred)")
    if mode == "none" and truth_ps:
        raise ValueError("--mode none takes no truth tables")
    if mode == "donor" and not truth_ps:
        raise ValueError("--mode donor needs a truth table")
    if panel is not None and panel not in PANELS:
        raise ValueError(f"--panel {panel}: only A, C and E are built (D and Cg are deferred)")
    if panel is not None and mode == "none":
        raise ValueError(f"--panel {panel} needs --mode donor: its bars are judged against truth")
    if panel == "A" and not focus:
        raise ValueError("--panel A needs --focus FAMILY")
    if focus is not None and panel != "A":
        raise ValueError(f"--focus {focus} belongs to --panel A (this is panel {panel})")
    donor = mode == "donor"
    excl = set()
    if exclude_p:
        with open(exclude_p) as fh:
            excl = set(json.load(fh))                                       # a JSON list of families, or a dict keyed by them
    if panel == "A" and focus in excl:
        raise ValueError(f"--focus {focus} is in --exclude: an excluded family is scored apart and has no bars")
    for t in truth_ps:
        need_columns(t, "read", "class", "truth")
    if n4_p:
        need_columns(n4_p, "read", "corrector")
    truth = P._truth(truth_ps, True) if donor else None
    n4 = {r["read"]: r["corrector"] for r in _tsv(n4_p)[1]} if n4_p else None
    recs = load_reads(out_dir, truth, n4)
    fams = by_family(recs)
    if panel == "A" and focus not in fams:
        raise ValueError(f"--focus {focus}: no tied read of that family in {out_dir}")
    files = load_files(out_dir, fams)
    arms = [a for a in ARMS if a != "N4" or n4 is not None]
    opts = options_text(mode, truth_ps, n4_p, panel, focus, exclude_p)
    pooled_recs = [r for r in recs if r["family"] not in excl]
    blocks = {f: block(v, arms, donor) for f, v in fams.items()}
    s = dict(mode=mode, panel=panel, focus=focus, out_dir=out_dir, truth_files=list(truth_ps), n4_file=n4_p,
             exclude_file=exclude_p, options=opts)
    s.update(block(pooled_recs, arms, donor))
    s["families"] = {f: b for f, b in blocks.items() if f not in excl}
    s["excluded"] = {f: b for f, b in blocks.items() if f in excl}
    if n4 is not None and not any(r["n4_held"] for r in recs):
        print(f"[o2vc] --n4 {n4_p} holds none of the tied reads: N4 is N1's on every read", file=sys.stderr, flush=True)
    s["bars"] = panel_bars(panel, s, blocks[focus] if panel == "A" else None) if panel else {}
    s["verdict"] = verdict_of(s["bars"]) if panel else None
    for f in fams:
        write_family(out_dir, f, fams[f], blocks[f], files[f], donor, opts)
    name = ".".join(["score"] + [x for x in (panel, focus) if x]) + ".json"
    with open(os.path.join(out_dir, name), "w") as fh:
        json.dump(s, fh, indent=1)
    report(s)
    return 0


RUN_OPTS = ("--families", "--threads", "--exclude-chrM")
SCORE_OPTS = ("--mode", "--n4", "--panel", "--focus", "--exclude")


def main(argv):
    opt, pos, i = {}, [], 0
    while i < len(argv):
        if argv[i] in RUN_OPTS[:2] + SCORE_OPTS:
            if i + 1 >= len(argv):
                raise SystemExit(__doc__)
            opt[argv[i]], i = argv[i + 1], i + 2
        elif argv[i] == "--exclude-chrM":
            opt[argv[i]], i = True, i + 1
        else:
            pos.append(argv[i])
            i += 1
    if pos[:1] == ["run"] and len(pos) == 5 and set(opt) <= set(RUN_OPTS):
        fams = set(opt["--families"].split(",")) if "--families" in opt else None
        return cmd_run(*pos[1:], families=fams, exclude_chrM=opt.get("--exclude-chrM", False),
                       threads=int(opt.get("--threads", 2)))
    if pos[:1] == ["score"] and len(pos) >= 2 and "--mode" in opt and set(opt) <= set(SCORE_OPTS):
        return cmd_score(pos[1], pos[2:], opt["--mode"], n4_p=opt.get("--n4"), panel=opt.get("--panel"), focus=opt.get("--focus"),
                         exclude_p=opt.get("--exclude"))
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
