#!/usr/bin/env python3
"""Tandem arrays (user, 2026-11-19): how assemblers build transcripts in tandem arrays; can Iso-Seq reads that map to an array
be told from mapping artifacts; a rule for the size of each member.

  tandem.py --units UNITS.tsv --bam BAM --fasta GENOME.fa --gtf NAME=GTF [...] --out PREFIX [--min-units 4] [--max-cv 0.35]

Arrays: rows of `undecidable.py arrays`'s units.tsv (CAT genes sharing a name stem, <= 200 kb apart); kept as TANDEM when
>= --min-units units, >= 80% on one strand, and the spacing of consecutive unit starts has a coefficient of variation <= --max-cv.
Period P = median spacing.
Per transcript (multi- or single-exon) overlapping a tandem array: units its exons touch (1 = one member; >= 2 = a multi-unit
model), its genomic span / the median unit length.
Per primary spliced record in an array, each intron is a PERIOD JUMP when its length is within 2% of k x P (k = 1..5) AND its
donor and acceptor sit at the same offset inside their units (within 2% of P): the read continued at the homologous position of
another unit, the signature of a read from one unit aligned across two. Its splice motif (GT-AG and the other canonical pairs,
either strand) is printed beside. Per primary unspliced record: span / median unit length.
Writes PREFIX.arrays.tsv (per array), PREFIX.models.tsv, PREFIX.jumps.tsv and PREFIX.summary.txt."""
import argparse, collections, csv, os, re, statistics, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from undecidable import read_gtf_tx

CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC"), ("CT", "AC"), ("CT", "GC"), ("GT", "AT")}

def tandem_arrays(units_tsv, min_units, max_cv):
    by = collections.defaultdict(list)
    with open(units_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["group"] in ("protein_coding", "unprocessed_pseudogene", "processed_pseudogene", "pseudogene_unsplit", "lncRNA"):
                by[r["array"]].append((int(r["start"]), int(r["end"]), r["strand"], r["gene"], r["chrom"]))
    out = {}
    for a, us in by.items():
        us.sort()
        if len(us) < min_units:
            continue
        st = collections.Counter(u[2] for u in us).most_common(1)[0]
        if st[1] < 0.8 * len(us):
            continue
        sp = [b[0] - a_[0] for a_, b in zip(us, us[1:])]
        cv = statistics.pstdev(sp) / statistics.mean(sp) if statistics.mean(sp) else 9
        if cv > max_cv:
            continue
        out[a] = dict(chrom=us[0][4], strand=st[0], units=[(s, e) for s, e, *_ in us], period=statistics.median(sp),
                      unit_len=statistics.median(e - s for s, e, *_ in us), cv=cv, start=us[0][0], end=max(e for _, e, *_ in us),
                      genes=collections.Counter(re.sub(r"\d+[A-Z]?\d*P?$", "", u[3]) for u in us).most_common(1)[0][0])
    return out

def units_touched(exons, units):
    return {i for i, (s, e) in enumerate(units) for a, b in exons if a < e and s < b}

def unit_offset(pos, units):
    for i, (s, e) in enumerate(units):
        if s - 1000 <= pos < e + 1000:
            return i, pos - s
    return None, None

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--units", required=True); ap.add_argument("--bam", required=True); ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", action="append", default=[]); ap.add_argument("--out", required=True)
    ap.add_argument("--min-units", type=int, default=4); ap.add_argument("--max-cv", type=float, default=0.35)
    a = ap.parse_args(argv)
    A = tandem_arrays(a.units, a.min_units, a.max_cv)
    fa = pysam.FastaFile(a.fasta)
    L = [f"tandem arrays: {len(A)} (>= {a.min_units} units, >= 80% one strand, spacing CV <= {a.max_cv})"]
    # models
    mrows, msum = [], collections.defaultdict(collections.Counter)
    for spec in a.gtf:
        name, g = spec.split("=", 1)
        for t, v in read_gtf_tx(g).items():
            for an, x in A.items():
                if v["chrom"] != x["chrom"] or v["exons"][-1][1] <= x["start"] or v["exons"][0][0] >= x["end"]:
                    continue
                u = units_touched(v["exons"], x["units"])
                if not u:
                    continue
                span = v["exons"][-1][1] - v["exons"][0][0]
                k = len(u)
                msum[(name, an)]["models"] += 1; msum[(name, an)]["multi_unit"] += k >= 2
                mrows.append((name, an, t, len(v["exons"]), k, span, span / x["unit_len"]))
    with open(a.out + ".models.tsv", "w") as fh:
        fh.write("arm\tarray\ttranscript\texons\tunits_touched\tspan\tspan_over_unit\n")
        fh.writelines("\t".join(map(str, r[:6])) + f"\t{r[6]:.2f}\n" for r in mrows)
    # reads
    jumps, rsum = [], collections.defaultdict(collections.Counter)
    with pysam.AlignmentFile(a.bam) as bam:
        for an, x in A.items():
            P = x["period"]
            for r in bam.fetch(x["chrom"], x["start"], x["end"]):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rs = rsum[an]; rs["primaries"] += 1
                pos, introns = r.reference_start, []
                for op, n in r.cigartuples:
                    if op in (0, 2, 7, 8):
                        pos += n
                    elif op == 3:
                        introns.append((pos, pos + n)); pos += n
                if not introns:
                    rs["unspliced"] += 1
                    rs["unspliced_gt_unit"] += (r.reference_end - r.reference_start) > 1.2 * x["unit_len"]
                    continue
                rs["spliced"] += 1
                has_jump = False
                for d, ac in introns:
                    ln = ac - d
                    k = round(ln / P)
                    ui, od = unit_offset(d, x["units"]); uj, oa = unit_offset(ac, x["units"])
                    jump = 1 <= k <= 5 and abs(ln - k * P) <= 0.02 * k * P and od is not None and oa is not None and abs(od - oa) <= 0.02 * P
                    if jump:
                        m5, m3 = fa.fetch(x["chrom"], d, d + 2).upper(), fa.fetch(x["chrom"], ac - 2, ac).upper()
                        canon = (m5, m3) in CANON
                        rs["jump_introns"] += 1; rs["jump_canonical"] += canon; has_jump = True
                        jumps.append((an, r.query_name, d, ac, ln, k, f"{m5}-{m3}", int(canon)))
                    rs["introns"] += 1
                rs["reads_with_jump"] += has_jump
    with open(a.out + ".jumps.tsv", "w") as fh:
        fh.write("array\tread\tdonor\tacceptor\tlength\tk\tmotif\tcanonical\n")
        fh.writelines("\t".join(map(str, j)) + "\n" for j in jumps)
    with open(a.out + ".arrays.tsv", "w") as fh:
        arms = [s.split("=", 1)[0] for s in a.gtf]
        fh.write("array\tchrom\tunits\tperiod\tunit_len\tcv\tprimaries\tspliced\treads_with_jump\tjump_introns\tjump_canonical\tunspliced\tunspliced_gt_unit\t"
                 + "\t".join(f"{n}_models\t{n}_multi_unit" for n in arms) + "\n")
        for an, x in sorted(A.items(), key=lambda kv: -rsum[kv[0]]["primaries"]):
            rs = rsum[an]
            fh.write(f"{an}\t{x['chrom']}\t{len(x['units'])}\t{x['period']:.0f}\t{x['unit_len']:.0f}\t{x['cv']:.2f}\t{rs['primaries']}\t{rs['spliced']}\t"
                     f"{rs['reads_with_jump']}\t{rs['jump_introns']}\t{rs['jump_canonical']}\t{rs['unspliced']}\t{rs['unspliced_gt_unit']}\t"
                     + "\t".join(f"{msum[(n, an)]['models']}\t{msum[(n, an)]['multi_unit']}" for n in arms) + "\n")
    tot = collections.Counter()
    for rs in rsum.values():
        tot.update(rs)
    L.append(f"primaries in tandem arrays {tot['primaries']}: spliced {tot['spliced']}, with a period jump {tot['reads_with_jump']} "
             f"({tot['jump_introns']} jump introns, canonical motif {tot['jump_canonical']}); unspliced {tot['unspliced']}, longer than 1.2 units {tot['unspliced_gt_unit']}")
    for n in [s.split("=", 1)[0] for s in a.gtf]:
        rows = [r for r in mrows if r[0] == n]
        L.append(f"{n}: models in tandem arrays {len(rows)}, touching >= 2 units {sum(r[4] >= 2 for r in rows)}, "
                 f"span > 1.5 units {sum(r[6] > 1.5 for r in rows)}")
    open(a.out + ".summary.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))

if __name__ == "__main__":
    main()


def fold_member(reads, units, period, min_share=0.2, bin_bp=25):
    """The member as the reads see it. Each read (start, end, strand) is folded onto the unit whose start is nearest at or
    before its start (+ half a period of slack): offsets relative to that unit's start. Per strand, the member = the contiguous
    run of folded positions covered by >= min_share of the strand's peak folded coverage around the peak; reads whose folded
    end exceeds one period are MULTI-UNIT (spanning into the next unit). -> {strand: dict(start, end, size, reads, multi_unit,
    tss_mode, tes_mode)}"""
    starts = [u[0] for u in units]
    out = {}
    for st in "+-":
        cov, ss, es, n, multi = collections.Counter(), collections.Counter(), collections.Counter(), 0, 0
        for s, e, rs in reads:
            if rs != st:
                continue
            i = max((k for k, u0 in enumerate(starts) if u0 <= s + period / 2), default=None)
            if i is None:
                continue
            fs, fe = s - starts[i], e - starts[i]
            n += 1
            if fe - fs > period or fe > period + period / 2:
                multi += 1
                continue
            for b in range(fs // bin_bp, fe // bin_bp + 1):
                cov[b] += 1
            ss[fs // bin_bp] += 1; es[fe // bin_bp] += 1
        if not cov:
            continue
        peak = max(cov, key=cov.get); thr = min_share * cov[peak]
        lo = peak
        while cov.get(lo - 1, 0) >= thr:
            lo -= 1
        hi = peak
        while cov.get(hi + 1, 0) >= thr:
            hi += 1
        tss = max((b for b in ss if lo <= b <= hi), key=ss.get, default=lo)
        tes = max((b for b in es if lo <= b <= hi), key=es.get, default=hi)
        out[st] = dict(start=lo * bin_bp, end=(hi + 1) * bin_bp, size=(hi + 1 - lo) * bin_bp, reads=n, multi_unit=multi,
                       tss_mode=tss * bin_bp, tes_mode=(tes + 1) * bin_bp)
    return out


def fold_pos(pos, units, period):
    starts = [u[0] for u in units]
    i = max((k for k, u0 in enumerate(starts) if u0 <= pos + period / 2), default=None)
    return (None, None) if i is None else (i, pos - starts[i])

def classify_model(exons, units, period, motif_of, tol=0.02, same_bp=10):
    """Artifact signatures of a model in a tandem array (each returns True when present):
    jump    an intron whose length is within tol of k x period (k = 1..9), donor and acceptor at the same offset inside their
            units (within tol x period), and a non-canonical motif: one unit's read aligned across to another
    repeat  the same intron (folded donor and acceptor within same_bp) in >= 2 different units: members concatenated
    bridge  an exon that covers a whole inter-unit gap (from the end of a unit to the start of the next)"""
    ints = [(a[1], b[0]) for a, b in zip(exons, exons[1:])]
    jump = False
    for d, ac in ints:
        ln = ac - d; k = round(ln / period)
        if 1 <= k <= 9 and abs(ln - k * period) <= tol * k * period:
            (_, od), (_, oa) = fold_pos(d, units, period), fold_pos(ac, units, period)
            if od is not None and oa is not None and abs(od - oa) <= tol * period and not motif_of(d, ac):
                jump = True
    folded = collections.defaultdict(set)
    for d, ac in ints:
        (i, od), (_, oa) = fold_pos(d, units, period), fold_pos(ac, units, period)
        if i is not None:
            folded[(od // same_bp, oa // same_bp)].add(i)
    repeat = any(len(v) >= 2 for v in folded.values())
    gaps = [(units[i][1], units[i + 1][0]) for i in range(len(units) - 1) if units[i][1] < units[i + 1][0]]
    bridge = any(s <= g0 and g1 <= e for s, e in exons for g0, g1 in gaps)
    return dict(jump=jump, repeat=repeat, bridge=bridge)
