#!/usr/bin/env python3
"""O2 from the sample's own variant calls (design docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md, plan
docs/superpowers/plans/2026-10-22-o2-variant-calls.md): the CLI over o2vc_rule.py (the decision rule) and o2vc_data.py (the
BAM, VCF and caller side).

    o2vc.py run BAM COPIES.tsv GENOME.fa OUTDIR [--families A,B] [--exclude-chrM] [--threads N]

run never reads truth. Per family, in sorted order -- a family with no AS-tied read is skipped, and under --exclude-chrM a
family with a copy on chrM:
1. its copies' teachers (o2vc_data.scan) go to OUTDIR/teach/<family>.bam, base qualities of 0 written as 1 (longcallR
   2.0.0 panics on a Q0 at a variant site);
2. longcallR runs once per copy with >= 1 teacher (caller_region: its span and PAD bases each side) against the copy's
   chromosome alone, OUTDIR/ref/<chrom>.fa, into OUTDIR/vcf/<copy>.vcf. A VCF already there is reused, never re-called:
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
  hold a foreign allele: kind separating_call (VC separates a pair of copies the reference does not), separating_ref
  (every pair VC separates, the reference does too), foreign or editing; copies: copy|chrom|pos|ref|genotype|depth per
  tied copy, ";"-joined (genotype "." when unknown, alleles "/"-joined; depth "." when not measured).
- <family>.flags.tsv: flag n_reads reads alleles (copy|chrom|pos|base, ";"-joined).
- <family>.teachers.tsv: copy teachers records pass lowqual homref multiallelic rnaedit dense teacher_flag_reads -- the
  FILTERs of the copy's VCF, and its teachers in a flag of the teachers' own-copy foreign alleles.
- vcf/<copy>.vcf, teach/<family>.bam, ref/<chrom>.fa, run.json.
"""
import bisect
import collections
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


D, R, P = _sibling("o2vc_data"), _sibling("o2vc_rule"), _sibling("o2psv")


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
            kind = R.foreign(vc_cols[-1])
            if not kind:
                sv = separated(gv)
                kind = "" if not sv else "separating_ref" if sv <= separated(gn) else "separating_call"
            if not kind:
                continue
            if kind == "foreign":
                foreign.append((col["p"], frozenset(
                    (c, *col["at"][c], col["base"].translate(R.COMP) if t.copies[c].flip else col["base"]) for c in cs)))
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


def main(argv):
    opt, pos, i = {}, [], 0
    while i < len(argv):
        if argv[i] in ("--families", "--threads"):
            opt[argv[i]], i = argv[i + 1], i + 2
        elif argv[i] == "--exclude-chrM":
            opt[argv[i]], i = True, i + 1
        else:
            pos.append(argv[i])
            i += 1
    if pos[:1] == ["run"] and len(pos) == 5:
        fams = set(opt["--families"].split(",")) if "--families" in opt else None
        return cmd_run(*pos[1:], families=fams, exclude_chrM=opt.get("--exclude-chrM", False),
                       threads=int(opt.get("--threads", 2)))
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
