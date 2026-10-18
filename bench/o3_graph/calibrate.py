#!/usr/bin/env python3
"""Measure real long-read RNA data so the simulator can match it (2026-10-16 goal: more realistic simulations).

    calibrate.py profile BAM OUT.json [--chroms chr20,chr21,chr22] [--max-reads 30000] [--min-mapq 60]
    calibrate.py expression BAM FAMILIES.json OUT.json
    calibrate.py errors BAM OUT.json [--chroms ...] [--max-reads 30000]   errors and donor alleles over every uniquely placed read
    calibrate.py tags BAM OUT.json [--regions chr20,...|--all] [--min-mapq 60] [--max-reads 20000]
    calibrate.py paralogs BAM FAMILIES.json OUT.json [--min-mapq 1] [--min-coding 10]

`profile`, from primary MAPQ >= min-mapq reads lying inside single-isoform genes (CAT/Liftoff, CHM13 v2.0; unique placement,
so every difference is a sequencing error or the donor's own allele, never a paralog's):
  errors      substitutions, insertions and deletions per aligned base, after removing the donor's alleles (a mismatch shared
              by >= MIN_VAR_READS reads at >= HET_FRAC of the reads covering it: het below HOM_FRAC, hom at or above); the
              fraction of indels in homopolymer context (a run of >= 3 of the base, reference or read side) and their lengths;
              soft-clipped bases
  donor       het and hom allele sites seen per read base (the divergence a read inherits from the donor, not from errors)
  ends        transcript bases missing at the 5' and 3' end of the read (the transcript's strand), quantiles
  retention   reads with >= RETAIN_MIN bases aligned inside an annotated intron
  length      read length quantiles
`tags`: per uniquely placed primary read (MAPQ >= min-mapq): minimap2's de and NM, mismatches (X), insertion and deletion
events, the homopolymer share of indels and soft clips; medians and per-kb rates (real data and simulated worlds alike).
`expression`: per family of FAMILIES.json, every copy's primary alignments (all, and MAPQ >= 20), and the spread of log counts.
`paralogs` (round 19): in families of seqfamilies.py build --paralogs, how the library expresses a copy of another biotype than
protein-coding (a pseudogene, a lncRNA) against its family's coding copies: uniquely placed reads (MAPQ >= 1) per copy, the
fraction of such copies silent and the median ratio of the expressed (paralog_expression).
"""
import argparse
import bisect
import collections
import importlib.util
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MIN_VAR_READS, HET_FRAC, HOM_FRAC, RETAIN_MIN, HP_RUN = 3, 0.2, 0.8, 50, 3
QUANTS = (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(ref_base, p, b, step, cap=50):
    n = 0
    while n < cap and ref_base(p + step * n) == b:
        n += 1
    return n


def cigar_events(cig, ref_start, query, ref_base):
    """An alignment's events from its =/X CIGAR: aligned bases, mismatches [(ref pos, read base)], deletions and insertions
    [(ref pos, length, homopolymer?)], soft-clipped bases."""
    out = dict(aligned=0, mismatch=[], ins=[], clip=0, blocks=[])
    out["del"] = []
    r, q, bs = ref_start, 0, ref_start
    for op, n in cig:
        if op in (7, 0):                                   # = (or M, counted as aligned)
            out["aligned"] += n
            r += n
            q += n
        elif op == 8:                                      # X
            for k in range(n):
                out["mismatch"].append((r + k, query[q + k]))
            out["aligned"] += n
            r += n
            q += n
        elif op == 2:                                      # D
            bases = {ref_base(r + k) for k in range(n)}
            hp = False
            if len(bases) == 1:
                b = bases.pop()
                hp = _run(ref_base, r - 1, b, -1) + n + _run(ref_base, r + n, b, 1) >= HP_RUN
            out["del"].append((r, n, hp))
            r += n
        elif op == 1:                                      # I
            ins = set(query[q:q + n])
            hp = False
            if len(ins) == 1:
                b = ins.pop()
                hp = _run(ref_base, r - 1, b, -1) + _run(ref_base, r, b, 1) + n >= HP_RUN
            out["ins"].append((r, n, hp))
            q += n
        elif op == 3:                                      # N: an intron closes the block
            out["blocks"].append((bs, r))
            r += n
            bs = r
        elif op == 4:
            out["clip"] += n
            q += n
    out["blocks"].append((bs, r))
    return out


def classify_mismatches(counts, cov, min_reads=MIN_VAR_READS, het=HET_FRAC, hom=HOM_FRAC):
    """{(chrom, pos, alt): reads} with coverage {(chrom, pos): reads} -> het / hom (the donor's) or error."""
    out = {}
    for (ch, p, alt), n in counts.items():
        f = n / max(1, cov.get((ch, p), 0))
        out[(ch, p, alt)] = ("hom" if f >= hom else "het") if n >= min_reads and f >= het else "error"
    return out


def transcript_ends(blocks, exons, strand):
    """(5' missing, 3' missing) transcript bases of a read's aligned blocks on a transcript (exons, genome order), or None."""
    exons = sorted(exons)
    tx_len = sum(e - s for s, e in exons)
    cov = []
    off = 0
    for s, e in exons:
        for bs, be in blocks:
            lo, hi = max(s, bs), min(e, be)
            if lo < hi:
                cov += [off + lo - s, off + hi - s - 1]
        off += e - s
    if not cov:
        return None
    first, last = min(cov), max(cov)
    left, right = first, tx_len - last - 1
    return (left, right) if strand == "+" else (right, left)


def intron_bases(blocks, exons):
    exons = sorted(exons)
    introns = [(a[1], b[0]) for a, b in zip(exons, exons[1:])]
    return sum(max(0, min(be, ie) - max(bs, is_)) for bs, be in blocks for is_, ie in introns)


def retained_introns(blocks, exons):
    """Annotated introns that one aligned block spans end to end: the read kept the intron (not an exon inside it)."""
    exons = sorted(exons)
    return sum(1 for a, b in zip(exons, exons[1:]) if any(bs <= a[1] and be >= b[0] for bs, be in blocks))


def tag_summary(rows):
    """Medians of de and NM per kb, per-kb rates of mismatches and indel events, the homopolymer share of indels, clip fraction."""
    kb = sum(r["aligned"] for r in rows) / 1000
    return dict(reads=len(rows), de_median=statistics.median(r["de"] for r in rows),
                nm_per_kb_median=statistics.median(1000 * r["nm"] / max(1, r["aligned"]) for r in rows),
                x_per_kb=sum(r["x"] for r in rows) / kb, ins_per_kb=sum(r["ins"] for r in rows) / kb,
                del_per_kb=sum(r["dele"] for r in rows) / kb,
                hp_indel_frac=sum(r["hp"] for r in rows) / max(1, sum(r["indels"] for r in rows)),
                clip_frac=sum(r["clip"] for r in rows) / max(1, sum(r["qlen"] for r in rows)))


def expression_log_sd(counts):
    """The between-copy spread of log expression (round 9 review I2): over the copies with reads (an unread copy is silent,
    not spread), the variance of log(n) less its Poisson counting noise (~1/n per copy, delta method), floored at 0."""
    n = [c for c in counts if c > 0]
    if len(n) < 2:
        return 0.0
    logs = [math.log(c) for c in n]
    return math.sqrt(max(0.0, statistics.pvariance(logs) - statistics.mean(1 / c for c in n)))


def adapter_consensus(clips):
    """The majority base per position of the clipped sequences at their most common length."""
    L = collections.Counter(len(c) for c in clips).most_common(1)[0][0]
    same = [c for c in clips if len(c) == L]
    return "".join(collections.Counter(c[i] for c in same).most_common(1)[0][0] for i in range(L))


def _quantiles(v):
    v = sorted(v)
    return {str(q): v[min(len(v) - 1, int(q * len(v)))] for q in QUANTS} if v else {}


def error_stats(reads):
    """Errors and the donor's alleles over reads [{ch, qlen, ev: cigar_events}]: a mismatch shared by >= MIN_VAR_READS reads at
    >= HET_FRAC of its coverage is the donor's (classify_mismatches), the rest are errors; indel rates, homopolymer shares and
    lengths; soft-clipped bases."""
    counts = collections.Counter((rd["ch"], p, b) for rd in reads for p, b in rd["ev"]["mismatch"])
    cand = collections.defaultdict(set)
    for (ch, p, b), n in counts.items():
        if n >= 2:
            cand[ch].add(p)
    cand = {ch: sorted(v) for ch, v in cand.items()}
    cov = collections.Counter()
    for rd in reads:
        pos = cand.get(rd["ch"], [])
        for bs, be in rd["ev"]["blocks"]:
            for p in pos[bisect.bisect_left(pos, bs):bisect.bisect_left(pos, be)]:
                cov[(rd["ch"], p)] += 1
    kind = classify_mismatches(counts, cov)
    agg, ins_len, del_len = collections.Counter(), collections.Counter(), collections.Counter()
    for rd in reads:
        ev = rd["ev"]
        agg["aligned"] += ev["aligned"]
        agg["clip"] += ev["clip"]
        agg["query"] += rd["qlen"]
        for p, b in ev["mismatch"]:
            agg["mm_" + kind[(rd["ch"], p, b)]] += 1
        for key, lens in (("ins", ins_len), ("del", del_len)):
            for _, n, hp in ev[key]:
                agg[key + "_events"] += 1
                agg[key + "_hp"] += hp
                lens[min(n, 4)] += 1
    A = max(1, agg["aligned"])
    return dict(reads=len(reads), aligned_bases=agg["aligned"], sub_error=agg["mm_error"] / A, ins_event=agg["ins_events"] / A,
                del_event=agg["del_events"] / A, ins_hp_frac=agg["ins_hp"] / max(1, agg["ins_events"]),
                del_hp_frac=agg["del_hp"] / max(1, agg["del_events"]),
                ins_len={str(k): v / max(1, agg["ins_events"]) for k, v in sorted(ins_len.items())},
                del_len={str(k): v / max(1, agg["del_events"]) for k, v in sorted(del_len.items())},
                clip_frac=agg["clip"] / max(1, agg["query"]), het_per_base=agg["mm_het"] / A, hom_per_base=agg["mm_hom"] / A)


def cmd_errors(a):
    """Errors and donor alleles over EVERY primary read with MAPQ >= min-mapq on the chromosomes (not only single-isoform genes;
    round 9 review I5: the single-isoform subset under-counts mismatches, 0.90 vs 1.42 per kb)."""
    import pysam
    F = _load("families")
    fa = pysam.FastaFile(F.GENOME)
    reads = []
    with pysam.AlignmentFile(a.bam) as bam:
        for ch in a.chroms.split(","):
            cache = {}

            def ref_base(p, ch=ch, cache=cache):
                if p not in cache:
                    cache[p] = fa.fetch(ch, max(0, p), p + 1).upper() if p >= 0 else "N"
                return cache[p]
            for r in bam.fetch(ch):
                if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < a.min_mapq or r.query_sequence is None:
                    continue
                reads.append(dict(ch=ch, qlen=r.query_length, ev=cigar_events(r.cigartuples, r.reference_start, r.query_sequence, ref_base)))
                if len(reads) >= a.max_reads:
                    break
            if len(reads) >= a.max_reads:
                break
    out = dict(source=a.bam, chroms=a.chroms, **error_stats(reads))
    json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps(out, indent=1))
    return 0


def cmd_profile(a):
    import pysam
    F = _load("families")
    fa = pysam.FastaFile(F.GENOME)
    gff = pysam.TabixFile(f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz")
    import csv
    genes = [r for r in csv.DictReader(open(f"{F.CAT}/chm13v2.0_CAT_Liftoff.genes.tsv"), delimiter="\t")
             if r["chrom"] in a.chroms.split(",") and r["n_transcripts"] == "1"]
    reads = []
    with pysam.AlignmentFile(a.bam) as bam:
        for g in genes:
            if len(reads) >= a.max_reads:
                break
            gs, ge = int(g["start0"]), int(g["end"])
            tx = F.gene_transcripts(gff, dict(gene_id=g["gene_id"], chrom=g["chrom"], start=gs, end=ge))
            if len(tx) != 1:
                continue
            exons = next(iter(tx.values()))
            if len(exons) < 2:
                continue
            ch = g["chrom"]
            cache = {}

            def ref_base(p, ch=ch, cache=cache):
                if p not in cache:
                    cache[p] = fa.fetch(ch, max(0, p), p + 1).upper() if p >= 0 else "N"
                return cache[p]
            for r in bam.fetch(ch, gs, ge):
                if (r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < a.min_mapq
                        or r.reference_start < gs - 1000 or r.reference_end > ge + 1000
                        or ("-" if r.is_reverse else "+") != g["strand"]):
                    continue
                ev = cigar_events(r.cigartuples, r.reference_start, r.query_sequence, ref_base)
                ends = transcript_ends(ev["blocks"], exons, g["strand"])
                if ends is None:
                    continue
                reads.append(dict(ch=ch, ev=ev, ends=ends, tx_len=sum(e - s for s, e in exons), qlen=r.query_length,
                                  intron=intron_bases(ev["blocks"], exons), kept=retained_introns(ev["blocks"], exons),
                                  n_introns=len(exons) - 1))
                if len(reads) >= a.max_reads:
                    break
    agg = collections.Counter()
    for rd in reads:
        agg["retained"] += rd["intron"] >= RETAIN_MIN
        agg["kept_reads"] += rd["kept"] > 0
        agg["kept_introns"] += rd["kept"]
        agg["introns"] += rd["n_introns"]
    prof = dict(source=a.bam, chroms=a.chroms, **error_stats(reads),
                missing5_frac=_quantiles([rd["ends"][0] / rd["tx_len"] for rd in reads]),
                missing3_frac=_quantiles([rd["ends"][1] / rd["tx_len"] for rd in reads]),
                full5_frac=sum(1 for rd in reads if rd["ends"][0] <= 50) / max(1, len(reads)),
                retained_frac=agg["retained"] / max(1, len(reads)), read_length=_quantiles([rd["qlen"] for rd in reads]),
                kept_intron_read_frac=agg["kept_reads"] / max(1, len(reads)),
                kept_intron_per_intron=agg["kept_introns"] / max(1, agg["introns"]),
                end_pairs=[[round(rd["ends"][0] / rd["tx_len"], 4), round(rd["ends"][1] / rd["tx_len"], 4)] for rd in reads],
                missing5_frac_samples=sorted(round(rd["ends"][0] / rd["tx_len"], 4) for rd in reads),
                missing3_frac_samples=sorted(round(rd["ends"][1] / rd["tx_len"], 4) for rd in reads))
    json.dump(prof, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in prof.items() if not k.endswith("_samples") and k != "end_pairs"}, indent=1))
    return 0


def cmd_tags(a):
    import pysam
    F = _load("families")
    fa = pysam.FastaFile(F.GENOME)
    rows = []
    with pysam.AlignmentFile(a.bam) as bam:
        it = bam.fetch(until_eof=True) if a.all else (r for reg in a.regions.split(",") for r in bam.fetch(region=reg))
        for r in it:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
                continue
            ch, cache = r.reference_name, {}

            def ref_base(p, ch=ch, cache=cache):
                if p not in cache:
                    cache[p] = fa.fetch(ch, max(0, p), p + 1).upper() if p >= 0 else "N"
                return cache[p]
            ev = cigar_events(r.cigartuples, r.reference_start, r.query_sequence, ref_base)
            rows.append(dict(de=r.get_tag("de"), nm=r.get_tag("NM"), aligned=ev["aligned"], x=len(ev["mismatch"]),
                             ins=len(ev["ins"]), dele=len(ev["del"]), hp=sum(h for _, _, h in ev["ins"] + ev["del"]),
                             indels=len(ev["ins"]) + len(ev["del"]), clip=ev["clip"], qlen=r.query_length))
            if len(rows) >= a.max_reads:
                break
    s = tag_summary(rows)
    json.dump(dict(source=a.bam, **s), open(a.out, "w"), indent=1)
    print(json.dumps(s))
    return 0


def cmd_adapters(a):
    """Library adapters left on the reads: the soft clips at each end of primary MAPQ >= min-mapq reads (read orientation), their
    consensus (adapter_consensus over clips >= 40 bases) and every read's (5' clip, 3' clip) pair (sampled)."""
    import pysam
    import random
    comp = str.maketrans("ACGTN", "TGCAN")
    five, three, pairs = [], [], []
    with pysam.AlignmentFile(a.bam) as bam:
        for ch in a.chroms.split(","):
            for r in bam.fetch(ch):
                if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < a.min_mapq or r.query_sequence is None:
                    continue
                c, q = r.cigartuples, r.query_sequence
                if r.is_reverse:
                    q, c = q[::-1].translate(comp), c[::-1]
                l5 = c[0][1] if c[0][0] == 4 else 0
                l3 = c[-1][1] if c[-1][0] == 4 else 0
                pairs.append([l5, l3])
                if l5 >= 40:
                    five.append(q[:l5])
                if l3 >= 40:
                    three.append(q[len(q) - l3:])
                if len(pairs) >= a.max_reads:
                    break
            if len(pairs) >= a.max_reads:
                break
    out = dict(source=a.bam, reads=len(pairs), both_clipped_frac=sum(1 for x, y in pairs if x >= 20 and y >= 20) / max(1, len(pairs)),
               adapter5=adapter_consensus(five) if five else "", adapter3=adapter_consensus(three) if three else "",
               adapter_len_pairs=random.Random(1).sample(pairs, min(5000, len(pairs))))
    json.dump(out, open(a.out, "w"), indent=0)
    print(json.dumps({k: v for k, v in out.items() if k != "adapter_len_pairs"}, indent=1))
    return 0


def cmd_expression(a):
    import pysam
    fams = json.load(open(a.families))["sets"]
    out = {}
    with pysam.AlignmentFile(a.bam) as bam:
        for st, fs in fams.items():
            for f in fs:
                rows = []
                for c in f["copies"]:
                    n = n20 = 0
                    for r in bam.fetch(c["chrom"], c["start"], c["end"]):
                        if not (r.is_secondary or r.is_supplementary or r.is_unmapped):
                            n += 1
                            n20 += r.mapping_quality >= 20
                    rows.append((c["gene_name"], n, n20))
                logs = [math.log(n + 1) for _, n, _ in rows]
                out[f"{st}/{f['family']}"] = dict(copies=rows, log_sd=statistics.pstdev(logs) if len(logs) > 1 else 0.0,
                                                  log_sd_corrected=expression_log_sd([n for _, n, _ in rows]),
                                                  zero=sum(1 for _, n, _ in rows if n < 5))
    sds = [v["log_sd"] for v in out.values()]
    summ = dict(families=len(out), log_sd_median=statistics.median(sds), copies=sum(len(v["copies"]) for v in out.values()),
                copies_lt5=sum(v["zero"] for v in out.values()))
    json.dump(dict(summary=summ, families=out), open(a.out, "w"), indent=1)
    print(json.dumps(summ, indent=1))
    return 0


MIN_CODING_READS = 10                                   # paralogs: a family's coding copies must have this many reads to compare


def overlaps_exons(blocks, exons):
    """Does an aligned block of a read overlap one of a copy's exons (half-open)? A spliced read whose intron spans the copy's
    exons does not."""
    return any(s < e2 and s2 < e for s, e in blocks for s2, e2 in exons)


def paralog_expression(rows, min_coding=MIN_CODING_READS):
    """The simulator's paralog expression (round 19; user 2026-10-18: paralogs of any biotype are copies like any other, expressed
    as the real library expresses them). rows: (family, biotype, uniquely placed reads) per copy. Against its family's coding
    copies (their mean read count; a family whose coding copies average under min_coding reads is left out), each copy of
    another biotype is silent (no read) or expressed at its ratio: the fraction silent and the median ratio of the expressed,
    pooled and by biotype."""
    by = collections.defaultdict(list)
    for f, bt, n in rows:
        by[f].append((bt, n))
    ratios, fams = [], set()
    for f, v in by.items():
        cod = [n for bt, n in v if bt == "protein_coding"]
        if not cod or statistics.mean(cod) < min_coding:
            continue
        c = statistics.mean(cod)
        for bt, n in v:
            if bt != "protein_coding":
                ratios.append((bt, n / c))
                fams.add(f)

    def summ(xs):
        nz = [r for r in xs if r > 0]
        return dict(members=len(xs), silent=len(xs) - len(nz), ratio_median=statistics.median(nz) if nz else 0.0)
    tot = summ([r for _, r in ratios])
    return dict(members=tot["members"], families=len(fams),
                silent_frac=tot["silent"] / tot["members"] if tot["members"] else 0.0, ratio_median=tot["ratio_median"],
                by_class={bt: summ([r for b, r in ratios if b == bt]) for bt in sorted({b for b, _ in ratios})})


def cmd_paralogs(a):
    """Per copy of every family with a paralog of another biotype: the primary reads at MAPQ >= min_mapq (uniquely placed) with
    an aligned block on one of its exons; then paralog_expression."""
    import pysam
    d = json.load(open(a.families))
    rows = []
    with pysam.AlignmentFile(a.bam) as bam:
        for fs in d["sets"].values():
            for f in fs:
                if all(c.get("biotype", "protein_coding") == "protein_coding" for c in f["copies"]):
                    continue
                for c in f["copies"]:
                    ex = sorted({tuple(x) for v in c["transcripts"].values() for x in v})
                    names = set()
                    for s, e in ex:
                        for r in bam.fetch(c["chrom"], s, e):
                            if (r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < a.min_mapq
                                    or r.query_name in names):
                                continue
                            if overlaps_exons(r.get_blocks(), ex):
                                names.add(r.query_name)
                    rows.append((f["family"], c.get("biotype", "protein_coding"), len(names)))
    summ = paralog_expression(rows, a.min_coding)
    json.dump(dict(summary=summ, bam=a.bam, families=a.families, min_mapq=a.min_mapq, min_coding=a.min_coding, rows=rows),
              open(a.out, "w"), indent=1)
    print(json.dumps(summ, indent=1))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("bam")
    p.add_argument("out")
    p.add_argument("--chroms", default="chr20,chr21,chr22")
    p.add_argument("--max-reads", type=int, default=30000)
    p.add_argument("--min-mapq", type=int, default=60)
    p.set_defaults(func=cmd_profile)
    p = sub.add_parser("errors")
    p.add_argument("bam")
    p.add_argument("out")
    p.add_argument("--chroms", default="chr20,chr21,chr22")
    p.add_argument("--max-reads", type=int, default=30000)
    p.add_argument("--min-mapq", type=int, default=60)
    p.set_defaults(func=cmd_errors)
    p = sub.add_parser("adapters")
    p.add_argument("bam")
    p.add_argument("out")
    p.add_argument("--chroms", default="chr21")
    p.add_argument("--max-reads", type=int, default=20000)
    p.add_argument("--min-mapq", type=int, default=60)
    p.set_defaults(func=cmd_adapters)
    p = sub.add_parser("tags")
    p.add_argument("bam")
    p.add_argument("out")
    p.add_argument("--regions", default="chr20,chr21,chr22")
    p.add_argument("--all", action="store_true")
    p.add_argument("--min-mapq", type=int, default=60)
    p.add_argument("--max-reads", type=int, default=20000)
    p.set_defaults(func=cmd_tags)
    p = sub.add_parser("expression")
    p.add_argument("bam")
    p.add_argument("families")
    p.add_argument("out")
    p.set_defaults(func=cmd_expression)
    p = sub.add_parser("paralogs")
    p.add_argument("bam")
    p.add_argument("families")
    p.add_argument("out")
    p.add_argument("--min-mapq", type=int, default=1)
    p.add_argument("--min-coding", type=int, default=MIN_CODING_READS)
    p.set_defaults(func=cmd_paralogs)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
