#!/usr/bin/env python3
"""The data side of O2 from the sample's own variant calls (design 2026-10-22). The caller side runs longcallR on the teacher
reads of a copy and reads its VCF into the Calls that o2vc_rule.genotype() turns into a copy's alleles at a position. The BAM
side finds the AS-tied reads (O2's domain; minimap2's primary may lie outside the tied set) and each copy's teachers (scan),
puts a tied read's bases in its own frame, as sequenced, beside each tied copy's reference base there (columns), writes the
teachers' primaries for the caller (write_teacher_bam) and counts teacher depth with the caller's read filters as the design
names them (teacher_depth). Nothing here decides a placement (o2vc_rule.py does) and nothing reads truth.

What longcallR 2.0.0 does, observed on SMN2 (HG002, chr5:70808816-70837687, primary records, `-p hifi-masseq --min-mapq 0
--no-bam-output`), and what the code below rests on:
- `-o PREFIX` writes PREFIX.vcf and nothing else: VCF 4.3, plain text, one sample ("Sample"). It does not make the directory of
  PREFIX (exit 101). A region without reads still exits 0 and leaves a header-only VCF, so a missing file is an error and an
  empty one is not.
- FILTER is one of PASS, LowQual, HomRef, RnaEdit, Multiallelic, Dense (the header's six; SMN2 shows the first four). FORMAT is
  GT:GQ:PS:DP:AF:PQ on a PASS or LowQual record (PS "." when unphased) and GT:GQ:DP:AF on the others, so FORMAT keys are read by
  name, never by position. REF is upper-cased though the CHM13 FASTA is soft-masked (46 of 47 SMN2 records sit on a lowercase
  base): compare it with the genome base upper-cased.
- `-r CHROM:START-END`: START is 1-based and must be >= 1 (0 panics: exit 101, no VCF), END is exclusive, so a record is in the
  region when START <= POS < END.
- It needs the BAM sorted and indexed and a `.fai` beside the FASTA (it does not make one: exit 101), and it reads the whole FASTA:
  16 s and 3.4 GB RSS against the 3.1 GB CHM13 file, 0.7 s and 0.5 GB against a one-contig FASTA of chr5, with the same calls.
- It is not deterministic: on identical inputs the phase labels of phased records flip (0|1 / 1|0) and in about 1 run in 10 two weak
  T>C records move from HomRef / LowQual to RnaEdit (SMN2, 41 runs at -t 1 and -t 2). The rule reads a GT as a set of alleles, so
  the labels never matter, but a re-call can differ at a weak RnaEdit-like site: keep the VCF a run used, do not re-call to
  reproduce it.
"""
import bisect
import collections
import dataclasses
import gzip
import importlib.util
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LONGCALLR = "/home/juanfra/miniforge3/envs/longcallr/bin/longcallR"
PRESET = "hifi-masseq"                                          # transcript-strand HiFi reads: no strand-bias filter
MIN_MAPQ = 0                                                    # the teachers are chosen by AS margin, not by MAPQ
SLACK = 9                                                       # tied: within 10 AS of the read's best (Eichler)
MIN_READ_LENGTH = 500                                           # the caller's --min-read-length
DIVERGENCE = 0.05                                               # --divergence as its --help states it (see teacher_depth)
DISTANCE_TO_READ_END = 20                                       # --distance-to-read-end as its --help states it
MISSING_QUAL = 255                                              # htslib's byte for a missing QUAL (SAM "*"); it passes MIN_BASEQ
READ_OPS, ALIGN_OPS = (0, 1, 4, 5, 7, 8), (0, 7, 8)             # CIGAR ops holding read bases (both clip kinds) / aligning one


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_vc", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_RULE = _sibling("o2vc_rule")
Call, COMP, MIN_BASEQ = _RULE.Call, _RULE.COMP, _RULE.MIN_BASEQ


def longcallr_version(binary=LONGCALLR):
    """The caller's own version line ("longcallR 2.0.0"), kept in run.json so a result names the caller that made it."""
    return subprocess.run([binary, "--version"], capture_output=True, text=True, check=True).stdout.strip()


def run_longcallr(teacher_bam, genome, out_prefix, binary=LONGCALLR, threads=2, region=None):
    """Run longcallR (PRESET, MIN_MAPQ, no BAM output, every other option its default) on teacher_bam against genome and return
    the VCF it wrote, out_prefix + ".vcf". teacher_bam must be sorted and indexed, genome a FASTA with its .fai. region is
    "CHROM:START-END" with START 1-based and >= 1 and END exclusive (the records with START <= POS < END), or None for the whole
    BAM. The directory of out_prefix is made first (longcallR does not). A non-zero exit raises RuntimeError with the last 20
    lines of stderr, as does an exit 0 that wrote no VCF."""
    cmd = [binary, "-b", teacher_bam, "-f", genome, "-p", PRESET, "--min-mapq", str(MIN_MAPQ), "--no-bam-output",
           "-t", str(threads), "-o", out_prefix] + (["-r", region] if region else [])
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = "\n".join(r.stderr.splitlines()[-20:])
    if r.returncode:
        raise RuntimeError(f"longcallR exited {r.returncode} ({' '.join(cmd)}):\n{tail}")
    vcf = out_prefix + ".vcf"
    if not os.path.exists(vcf):
        raise RuntimeError(f"longcallR wrote no VCF at {vcf} ({' '.join(cmd)}):\n{tail}")
    return vcf


def parse_vcf(path):
    """{(chrom, POS - 1): Call} for every record of a longcallR VCF, .vcf or .vcf.gz (POS is 1-based in the file, 0-based in the
    key). FORMAT keys are read by name. gt is the tuple of the GT's allele indices (0 = REF), () when it calls none ("./.",
    ".|.", "."): a GT missing only one allele counts as calling none too, for half a genotype would make a copy look homozygous.
    The record is kept either way, and o2vc_rule.genotype() leaves such a copy unknown. dp is the record's own DP (0 when absent):
    the caller's count, not the teacher depth the rule's MIN_DEPTH uses."""
    out = {}
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\r\n").split("\t")
            f = dict(zip(c[8].split(":"), c[9].split(":")))
            idx = re.split(r"[/|]", f["GT"])
            gt = tuple(int(i) for i in idx) if all(i.isdigit() for i in idx) else ()
            alts = () if c[4] == "." else tuple(c[4].split(","))
            out[c[0], int(c[1]) - 1] = Call(c[3], alts, gt, c[6], int(f["DP"]) if f.get("DP", ".").isdigit() else 0)
    return out


@dataclasses.dataclass
class Aln:
    """A read's best alignment on one tied copy. cigar is in reference order, as the BAM holds it; flip is the record's
    is_reverse, so the CIGAR walks the reverse complement of the read as sequenced. NM, de and ms are None when absent."""
    copy: str
    chrom: str
    start: int
    cigar: tuple
    flip: bool
    AS: int
    NM: int
    de: float
    ms: int


@dataclasses.dataclass
class TiedRead:
    """An AS-tied read: its family; its best alignment on each tied copy; minimap2's primary copy (None off every copy; it may
    lie outside the tied set); the primary's MAPQ; the AS margin between its two best tied copies; and its bases and base
    qualities in the read frame (as sequenced), full length including soft clips."""
    family: str
    copies: dict[str, Aln]
    prim: str | None
    mapq: int
    margin: int
    seq: str
    qual: list[int]


@dataclasses.dataclass
class Scan:
    """scan()'s result: {read: TiedRead}, and {copy: the names of its teachers} listing every copy of the families scanned."""
    tied: dict[str, TiedRead]
    teachers: dict[str, set[str]]


def tied_wide(alns):
    """({copy: its best alignment within SLACK of the read's best}, minimap2's primary copy) for a read's alignments
    [dict(AS, prim, copy)], the first None when an alignment within SLACK lies off every copy or on another family, or the
    read ties on fewer than two copies. o2correct.tied() but for the design's one change (§1): minimap2's primary may lie
    outside the tied set, more than SLACK below the copies that share the best AS."""
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
    if len(per) < 2 or len({c.rsplit(".", 1)[0] for c in per}) > 1:
        return None, prim
    return per, prim


def is_teacher(primary_as, other_copy_as):
    """True when the primary's AS is more than SLACK (so >= 10) above every AS in other_copy_as: the read's alignments on the
    other copies of its primary's family (its own copy and other families are left out by the caller)."""
    return all(primary_as - a > SLACK for a in other_copy_as)


def _read_bases(cigar):
    """The read bases a CIGAR holds, both clip kinds counted."""
    return sum(n for op, n in cigar if op in READ_OPS)


def read_frame_positions(start, cigar, length, flip):
    """{read-frame position: reference position} of the read bases the alignment at `start` aligns (CIGAR ops 0, 7, 8). Ops 1,
    4 and 5 advance the read only (both clip kinds count: a secondary without SEQ may be hard-clipped), 2 and 3 the reference only.
    The CIGAR walks the read in reference order, so on a reverse record (flip) its read base q is position length - 1 - q of
    the read as sequenced. A CIGAR that does not hold `length` read bases raises ValueError: every column would shift."""
    out, q, r = {}, 0, start
    for op, n in cigar:
        if op in ALIGN_OPS:
            for i in range(n):
                out[length - 1 - q - i if flip else q + i] = r + i
            q, r = q + n, r + n
        elif op in (1, 4, 5):
            q += n
        elif op in (2, 3):
            r += n
    if q != length:
        raise ValueError(f"the CIGAR holds {q} read bases, the read {length}")
    return out


def columns(read, fetch):
    """The read's columns: the read-frame positions every tied alignment aligns, in order, each dict(p, base, qual, at={copy:
    (chrom, reference position)}, ref={copy: its reference base in the read frame}). fetch(chrom, lo, hi) returns the genome's
    bases [lo, hi); it is called once per copy and its bases upper-cased (the genome is soft-masked). A copy aligned opposite
    to the read (flip) has its base complemented, as o2vc_rule.to_read_frame() does a genotype."""
    cs = sorted(read.copies)
    pos = {c: read_frame_positions(read.copies[c].start, read.copies[c].cigar, len(read.seq), read.copies[c].flip) for c in cs}
    common = sorted(set.intersection(*(set(pos[c]) for c in cs))) if cs else []
    if not common:
        return []
    ref = {}
    for c in cs:
        a, rs = read.copies[c], [pos[c][p] for p in common]
        lo = min(rs)
        seq = fetch(a.chrom, lo, max(rs) + 1).upper()
        ref[c] = (lo, seq.translate(COMP) if a.flip else seq)
    return [dict(p=p, base=read.seq[p], qual=read.qual[p], at={c: (read.copies[c].chrom, pos[c][p]) for c in cs},
                 ref={c: ref[c][1][pos[c][p] - ref[c][0]] for c in cs}) for p in common]


def _read_frame(rec):
    """(SEQ, QUAL as a list, MAPQ, the read bases its CIGAR holds) of a primary record, SEQ and QUAL in the read frame
    (reverse-complemented when the record is reverse; a missing QUAL is MISSING_QUAL), or None when it has no SEQ."""
    seq, qual = rec.query_sequence, rec.query_qualities
    if seq is None:
        return None
    qual = [MISSING_QUAL] * len(seq) if qual is None else list(qual)
    if rec.is_reverse:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    return seq, qual, rec.mapping_quality, _read_bases(rec.cigartuples)


def scan(bam_p, copies_p, families=None):
    """Scan(tied, teachers) of a BAM against a copy table (o2psv.load_copies' format), both restricted to `families` when
    given. Two passes, as o2correct.collect(). Pass 1 keeps (AS, primary, copy) per alignment -- the copy is summary.attribute's
    among the copies whose span the alignment overlaps -- and decides the domain (tied_wide) and the teachers: the primary on
    copy c, the read not tied, is_teacher against its alignments on the other copies of c's family. Pass 2 reads only the tied
    reads' records: each tied copy's best (the first in the file on equal AS) and the primary's SEQ, QUAL and MAPQ. A tied
    read is skipped, and counted on stderr, when it has no primary SEQ, or when the primary's CIGAR or a tied alignment's holds
    a different number of read bases (both clip kinds) than that SEQ: its read-frame positions would shift."""
    import pysam
    P, S, C = _sibling("o2psv"), _sibling("summary"), _sibling("o2correct")
    cps = P.load_copies(copies_p)
    idx = C.span_index(cps)

    def family(c):
        return c.rsplit(".", 1)[0]

    def copy_of(rec):
        cands = C.overlapping(idx, rec.reference_name, rec.reference_start, rec.reference_end)
        return S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None

    def records(bam):
        for rec in bam.fetch(until_eof=True):
            if not (rec.is_unmapped or rec.is_supplementary) and rec.has_tag("AS"):
                yield rec

    def tag(rec, t):
        return rec.get_tag(t) if rec.has_tag(t) else None

    alns = collections.defaultdict(list)                        # pass 1: read -> [(AS, primary, copy)]
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in records(bam):
            alns[rec.query_name].append((rec.get_tag("AS"), not rec.is_secondary, copy_of(rec)))
    teachers = {c: set() for v in cps.values() for c, _, _, _ in v if families is None or family(c) in families}
    domain = {}                                                 # read -> (its tied copies sorted, minimap2's copy, AS margin)
    for q, al in alns.items():
        per, prim = tied_wide([dict(AS=a, prim=p, copy=c) for a, p, c in al])
        if per is not None:                                     # a tied read never teaches
            if families is None or family(min(per)) in families:
                top = sorted((x["AS"] for x in per.values()), reverse=True)
                domain[q] = (tuple(sorted(per)), prim, top[0] - top[1])
            continue
        own = next(((a, c) for a, p, c in al if p), None)
        if own is not None and own[1] in teachers:
            a0, c0 = own
            if is_teacher(a0, [a for a, _, c in al if c is not None and c != c0 and family(c) == family(c0)]):
                teachers[c0].add(q)
    del alns
    best, prims = collections.defaultdict(dict), {}
    with pysam.AlignmentFile(bam_p) as bam:                     # pass 2: the tied reads' records only
        for rec in records(bam):
            q = rec.query_name
            if q not in domain:
                continue
            if not rec.is_secondary:
                prims[q] = _read_frame(rec)
            c, a = copy_of(rec), rec.get_tag("AS")
            if c in domain[q][0] and (c not in best[q] or a > best[q][c].AS):
                best[q][c] = Aln(c, rec.reference_name, rec.reference_start, tuple(rec.cigartuples), rec.is_reverse, a,
                                 tag(rec, "NM"), tag(rec, "de"), tag(rec, "ms"))
    tied, no_seq, bad_len = {}, 0, 0
    for q in sorted(domain):
        cset, prim, margin = domain[q]
        got = prims.get(q)
        if got is None:
            no_seq += 1
            continue
        seq, qual, mapq, n = got
        if n != len(seq) or any(_read_bases(best[q][c].cigar) != len(seq) for c in cset):
            bad_len += 1
            continue
        tied[q] = TiedRead(family(cset[0]), {c: best[q][c] for c in cset}, prim, mapq, margin, seq, qual)
    print(f"[o2vc] scan {bam_p}: {len(domain)} tied reads, {len(tied)} kept; {bad_len} skipped (a CIGAR's read bases, both "
          f"clip kinds counted, differ from the primary's SEQ length), {no_seq} skipped (no primary SEQ); "
          f"{sum(len(v) for v in teachers.values())} teachers on {sum(1 for v in teachers.values() if v)} of "
          f"{len(teachers)} copies", file=sys.stderr, flush=True)
    return Scan(tied, teachers)


def write_teacher_bam(bam_p, spans, teachers, out_bam):
    """Write the primary records of each copy's teachers ({copy: names}) found over its span ({copy: (chrom, start, end)}) to
    out_bam, sorted and indexed (the caller's input), and return how many there are. A teacher's primary lies on its copy, so
    it overlaps the span; a read is written once. The directory of out_bam is made when missing."""
    import pysam
    os.makedirs(os.path.dirname(out_bam) or ".", exist_ok=True)
    tmp, n, seen = out_bam + ".unsorted.bam", 0, set()
    with pysam.AlignmentFile(bam_p) as bam, pysam.AlignmentFile(tmp, "wb", template=bam) as out:
        for c in sorted(spans):
            names = teachers.get(c)
            if not names:
                continue
            ch, s, e = spans[c]
            for rec in bam.fetch(ch, s, e):
                q = rec.query_name
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or q not in names or q in seen:
                    continue
                seen.add(q)
                out.write(rec)
                n += 1
    pysam.sort("-o", out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    return n


def _divergence(rec):
    """The record's de, rounded to 6 decimals (a BAM keeps it as a float32: 0.05 reads 0.0500000007), or NM / aligned bases
    when it has none (NM 0 when that is absent too)."""
    if rec.has_tag("de"):
        return round(rec.get_tag("de"), 6)
    aligned = sum(n for op, n in rec.cigartuples if op in ALIGN_OPS)
    return (rec.get_tag("NM") if rec.has_tag("NM") else 0) / max(aligned, 1)


def teacher_depth(teacher_bam, chrom, positions):
    """{position: the teacher reads counted there} on chrom, for every position given (0 where none is). A read counts when it
    is a primary with SEQ, its full length (the read bases its CIGAR holds, both clip kinds) is >= MIN_READ_LENGTH, its
    divergence (_divergence) is <= DIVERGENCE, and it aligns a base at the position (op 0, 7 or 8) of quality >= MIN_BASEQ (a
    missing QUAL passes as MISSING_QUAL) that is >= DISTANCE_TO_READ_END read bases from either end of the read: read base q
    of a read of length L counts when q >= 20 and L - 1 - q >= 20, clipped bases included.

    These are the filters the design names, longcallR's --help defaults. The HiFi presets of its 2.0.0 source differ (main.rs
    sets divergence 0.5 and distance-to-read-end 40; util.rs fill_data_into_freq_vec skips de >= 0.5 with no NM fallback,
    trims read ends only for ONT, and on HiFi drops only bases in a 5-base poly-A/T/C/G run within 40 bp of an aligned end;
    the depth it calls with, DP, is A+C+G+T at any base quality). So this count is at most the caller's, but for such
    poly-run bases, N bases (the caller counts neither) and a hard-clipped primary (its length floor is on SEQ; minimap2
    writes none). On Task 3's fixture world with an ALT on every teacher, longcallR wrote DP 9 at chrA:110, 125 and 400
    where this counts 0, 8 and 7."""
    import pysam
    want = sorted(set(positions))
    out = dict.fromkeys(want, 0)
    if not want:
        return out
    with pysam.AlignmentFile(teacher_bam) as bam:
        for rec in bam.fetch(chrom, want[0], want[-1] + 1):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.query_length:    # 0: no SEQ
                continue
            length = _read_bases(rec.cigartuples)
            i, j = bisect.bisect_left(want, rec.reference_start), bisect.bisect_left(want, rec.reference_end)
            if i == j or length < MIN_READ_LENGTH or _divergence(rec) > DIVERGENCE:
                continue
            qual = rec.query_qualities
            q, s, r = 0, 0, rec.reference_start                # read base (both clip kinds), SEQ index, reference position
            for op, n in rec.cigartuples:
                if op in ALIGN_OPS:
                    k = bisect.bisect_left(want, r, i, j)
                    while k < j and want[k] < r + n:
                        x = want[k] - r
                        if (DISTANCE_TO_READ_END <= q + x <= length - 1 - DISTANCE_TO_READ_END
                                and (MISSING_QUAL if qual is None else qual[s + x]) >= MIN_BASEQ):
                            out[want[k]] += 1
                        k += 1
                    q, s, r = q + n, s + n, r + n
                elif op in (1, 4):
                    q, s = q + n, s + n
                elif op == 5:
                    q += n
                elif op in (2, 3):
                    r += n
    return out
