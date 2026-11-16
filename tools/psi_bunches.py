#!/usr/bin/env python3
"""PSI of introns within bunches, with the range placement ties allow (driver stage `psi`; a report, nothing reads it).

    psi_bunches.py GTF BAM AS_TABLE OUT.tsv [--score ms|AS] [--window W | --ratio R] [--min-reads N] [--threads N] [--contigs A,B]

Bunch = the introns that share a donor, or that share an acceptor, on one strand (MntJULiP's unit, Yang 2022; an intron
can sit in one bunch of each kind). PSI of an intron = its reads over the reads of its bunch. Introns are read from the
alignments (CIGAR N), not from the GTF, so reads the assembler's filters dropped count too; the strand is the read's.
A read's records = its primary plus the secondaries whose score ties the primary's band, the bands of eq_classes.py:
`exact` = score equal to the read's best, `wide` = within W of it (--window, default 10) or >= R x best (--ratio);
wide columns carry the suffix _w. Per intron and band:
    placed  reads whose every record in the band carries the intron (a read with one record: that record)
    max     reads with any record in the band carrying it
    psi_placed  placed / the bunch's placed (NA when the bunch has none)
    psi_lo      placed / (placed + the other introns' max)      the lowest PSI the tied reads allow
    psi_hi      max / (max + the other introns' placed)         the highest
A wide [psi_lo, psi_hi] means the PSI is not measurable from these alignments. `placed` is where reads align best, NOT
origin: reads of a copy the reference lacks sit on their nearest copy, inside `placed`. Intron retention is not seen
(no read spans a retained intron). Kept: introns with max_w >= --min-reads (default 3) in bunches of >= 2 kept introns.
Columns locus / in_gtf: the GTF gene_id holding the intron (else the commonest among the bunch's introns, else `.`) and
whether an emitted transcript has the intron. With --contigs, tied records on other contigs are not seen (a read whose
primary is there is skipped; one whose primary is here counts as if its other records did not exist)."""
import argparse, collections, multiprocessing, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eq_classes import load_ties, floor_fn
from gtf_names import opn

def gtf_introns(gtf):
    """{(chrom, strand, start, end): gene_id} for the introns of the GTF's transcripts (1-based, first / last intron base)."""
    ex, loc = collections.defaultdict(list), {}
    for l in opn(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            ex[t].append((int(f[3]), int(f[4]))); loc[t] = (f[0], f[6], re.search(r'gene_id "([^"]+)"', f[8]).group(1))
    out = {}
    for t, e in ex.items():
        e.sort(); c, st, g = loc[t]
        for (_, a), (b, _) in zip(e, e[1:]):
            out[(c, st, a + 1, b - 1)] = g
    return out

def introns_of(pos, cigartuples):
    """introns (start, end), 1-based closed, of a record starting at 0-based `pos`."""
    out = []
    for op, n in cigartuples:
        if op == 3:
            out.append((pos + 1, pos + n))
        if op in (0, 2, 3, 7, 8):
            pos += n
    return out

def scan(contig):
    """-> Counter {(contig, strand, s, e): reads without a tie}, [(read, exact, (intron keys))] for the tied reads' records."""
    import pysam
    solo, tied = collections.Counter(), []
    for r in pysam.AlignmentFile(G_BAM).fetch(contig):
        if r.is_unmapped or r.is_supplementary:
            continue
        best = G_TIES.get(r.query_name)
        if best is None:
            if r.is_secondary:
                continue
        elif r.is_secondary and r.get_tag(G_TAG) < G_FLOOR(best):
            continue
        st = "-" if r.is_reverse else "+"
        keys = tuple((contig, st, s, e) for s, e in introns_of(r.reference_start, r.cigartuples))
        if best is None:
            solo.update(keys)
        else:
            tied.append((r.query_name, not r.is_secondary or r.get_tag(G_TAG) == best, keys))
    return solo, tied

def tally(solo, tied):
    """{intron: [placed, max, placed_w, max_w]}; tied = [(read, exact, introns)], one entry per record in the wide band."""
    cnt = collections.defaultdict(lambda: [0, 0, 0, 0])
    for k, n in solo.items():
        c = cnt[k]
        for i in range(4):
            c[i] += n
    recs = collections.defaultdict(list)
    for q, exact, keys in tied:
        recs[q].append((exact, frozenset(keys)))
    for rs in recs.values():
        if not any(x for x, _ in rs):
            continue  # the primary is on a contig --contigs left out: its secondaries alone say nothing
        for off, sets in ((0, [k for x, k in rs if x]), (2, [k for _, k in rs])):
            for k in frozenset.intersection(*sets):
                cnt[k][off] += 1
            for k in frozenset.union(*sets):
                cnt[k][off + 1] += 1
    return cnt

def psi(p, m, sp, sm):
    """(psi_placed, psi_lo, psi_hi) of an intron with placed p / max m in a bunch whose sums are sp / sm; None = undefined."""
    lo, hi = p + (sm - m), m + (sp - p)
    return (p / sp if sp else None, p / lo if lo else None, m / hi if hi else None)

def bunch_rows(cnt, known, min_reads=3):
    """rows [bunch, anchor, chrom, strand, start, end, locus, in_gtf, n, (placed, max, psi_placed, psi_lo, psi_hi) x 2 bands]."""
    groups = collections.defaultdict(list)
    for k, c in cnt.items():
        if c[3] >= min_reads:
            chrom, st, s, e = k
            don, acc = (s, e) if st == "+" else (e, s)
            groups[(chrom, st, "donor", don)].append(k); groups[(chrom, st, "acceptor", acc)].append(k)
    rows = []
    for (chrom, st, kind, pos), ks in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][3], kv[0][1], kv[0][2])):
        if len(ks) < 2:
            continue
        loci = collections.Counter(known[k] for k in ks if k in known)
        fallback = min(loci, key=lambda g: (-loci[g], g)) if loci else "."
        sums = [sum(cnt[k][i] for k in ks) for i in range(4)]
        for k in sorted(ks):
            c, row = cnt[k], [f"{chrom}:{st}:{kind[0].upper()}:{pos}", kind, chrom, st, k[2], k[3], known.get(k, fallback), int(k in known), len(ks)]
            for o in (0, 2):
                row += [c[o], c[o + 1], *psi(c[o], c[o + 1], sums[o], sums[o + 1])]
            rows.append(row)
    return rows

HEADER = "bunch anchor chrom strand start end locus in_gtf n_introns".split() + [c + s for s in ("", "_w") for c in "placed max psi_placed psi_lo psi_hi".split()]

def summary(rows, slack=0.1):
    """one line per band: bunches, and those where some intron's [psi_lo, psi_hi] is wider than `slack`."""
    out = []
    for name, o in (("exact", 9), ("wide", 14)):
        b, loose = set(), set()
        for r in rows:
            b.add(r[0])
            if r[o + 3] is None or r[o + 4] is None or r[o + 4] - r[o + 3] > slack:
                loose.add(r[0])
        out.append(f"{name}: {len(b)} bunches, {len(loose)} with a PSI range wider than {slack:g} ({len(loose) / max(1, len(b)):.4f})")
    return out

def main():
    global G_BAM, G_TIES, G_FLOOR, G_TAG
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gtf"); p.add_argument("bam"); p.add_argument("as_table"); p.add_argument("out")
    p.add_argument("--score", choices=("ms", "AS"), default="ms"); p.add_argument("--window", type=float, default=10); p.add_argument("--ratio", type=float)
    p.add_argument("--min-reads", type=int, default=3); p.add_argument("--threads", type=int, default=4); p.add_argument("--contigs")
    a = p.parse_args()
    import pysam
    known = gtf_introns(a.gtf)
    G_FLOOR = floor_fn(a.window, a.ratio)
    G_BAM, G_TAG, G_TIES = a.bam, a.score, load_ties(a.as_table, G_FLOOR, a.score)
    contigs = a.contigs.split(",") if a.contigs else list(pysam.AlignmentFile(a.bam).references)
    with multiprocessing.get_context("fork").Pool(a.threads) as pool:
        res = pool.map(scan, contigs, chunksize=1)
    solo, tied = collections.Counter(), []
    for s, t in res:
        solo.update(s); tied += t
    cnt = tally(solo, tied)
    rows = bunch_rows(cnt, known, a.min_reads)
    fmt = lambda v: "NA" if v is None else f"{v:.4f}" if isinstance(v, float) else str(v)
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.writelines("\t".join(map(fmt, r)) + "\n" for r in rows)
    os.replace(tmp, a.out)
    kept = {tuple(r[2:6]) for r in rows}
    print(f"{len(cnt)} introns in the reads, {sum(1 for c in cnt.values() if c[3] >= a.min_reads)} with >= {a.min_reads} reads, {len(kept)} in bunches "
          f"({sum(1 for k in kept if k in known)} in the GTF); tied reads {len({q for q, _, _ in tied})}")
    print("\n".join(summary(rows)))

if __name__ == "__main__":
    main()
