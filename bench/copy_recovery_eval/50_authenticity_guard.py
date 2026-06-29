#!/usr/bin/env python3
"""U5: tool-independent authenticity guard.

For each recovered copy we measure DECISIVE own-copy read evidence rather than a
primary-alignment *count*. A read is decisive for this copy iff its best alignment
here scores strictly better (higher AS) than its best alignment to any sister copy
in the same paralog family, or it maps uniquely to this copy. Reads that tie (equal
best AS on this copy and a sister) are genuinely undecidable. The copy is then
bucketed authentic / unresolvable / phantom (see lib_eval.classify_authenticity).

Why not count primaries: minimap2's primary flag among MAPQ-0 multimappers is
arbitrary, so a count both over-credits phantoms (a copy that wins the coin-flip)
and starves real copies whose reads all tie. `primary_support` is kept only as a
reference column.
"""
import argparse, csv, os, sys
import pysam

sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def primary_support(bam_path, chrom, exons):
    """Count distinct reads with a PRIMARY alignment overlapping any exon (reference column)."""
    seen = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for (start, end) in exons:
            for read in bam.fetch(chrom, start, end):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                seen.add(read.query_name)
    return len(seen)


def _best_as(bam, chrom, exons):
    """read_name -> (best_AS, best_mapq) over the given exon windows. Reads lacking an
    AS tag are ignored (HiFi minimap2 always emits AS)."""
    best = {}
    for (start, end) in exons:
        for r in bam.fetch(chrom, max(0, start), end):
            if r.is_unmapped:
                continue
            try:
                a = r.get_tag("AS")
            except KeyError:
                continue
            cur = best.get(r.query_name)
            if cur is None or a > cur[0]:
                best[r.query_name] = (a, r.mapping_quality)
    return best


def decisive_support(bam, chrom, copy_exons, sisters):
    """sisters: list of (chrom, exons). Returns dict with n_decisive / n_tied /
    n_foreign / n_total for reads overlapping this copy's exons.
      decisive -- best AS on this copy > best AS on every sister, or unique to copy.
      tied     -- best AS equal to the best sister AS.
      foreign  -- best AS lower than some sister (sister-owned spillover).
    `bam` is an open pysam.AlignmentFile (so all copies share one handle)."""
    copy_as = _best_as(bam, chrom, copy_exons)
    sister_as = {}
    for (sch, sex) in sisters:
        for rd, (a, _mq) in _best_as(bam, sch, sex).items():
            if rd not in sister_as or a > sister_as[rd]:
                sister_as[rd] = a
    n_dec = n_tied = n_foreign = 0
    for rd, (cas, mq) in copy_as.items():
        sas = sister_as.get(rd)
        if sas is None:
            n_dec += 1            # maps only to this copy within the family
        elif cas > sas:
            n_dec += 1
        elif cas == sas:
            n_tied += 1
        else:
            n_foreign += 1
    return {"n_decisive": n_dec, "n_tied": n_tied, "n_foreign": n_foreign,
            "n_total": len(copy_as)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--recovery", required=True, help="rustle-VG recovery.tsv (ref_transcript,family_id,fsm,ism)")
    ap.add_argument("--exons", required=True, help="TSV: transcript_id,chrom,exon_starts(csv),exon_ends(csv)")
    ap.add_argument("--universe", required=True, help="universe.tsv: transcript_id,gene_id,family_id,...")
    ap.add_argument("--k-decisive", type=int, default=2)
    ap.add_argument("--k-tied", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    exon_map = {}
    with open(args.exons) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            starts = [int(x) for x in row["exon_starts"].split(",") if x]
            ends = [int(x) for x in row["exon_ends"].split(",") if x]
            exon_map[row["transcript_id"]] = (row["chrom"], list(zip(starts, ends)))

    # family_id -> [transcript_id,...] (for sister lookup)
    fam_members = {}
    with open(args.universe) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            fam_members.setdefault(row["family_id"], []).append(row["transcript_id"])

    bam = pysam.AlignmentFile(args.bam, "rb")
    with open(args.recovery) as f, open(args.out, "w", newline="") as o:
        rows = list(csv.DictReader(f, delimiter="\t"))
        w = csv.writer(o, delimiter="\t")
        w.writerow(["ref_transcript", "family_id", "n_decisive", "n_tied",
                    "n_foreign", "primary_support", "status", "authentic"])
        out_rows = []
        for r in rows:
            tx = r["ref_transcript"]
            fam = r["family_id"]
            if tx not in exon_map:
                sys.stderr.write(f"[guard] no exons for {tx}, skipping\n"); continue
            chrom, exons = exon_map[tx]
            sisters = [exon_map[s] for s in fam_members.get(fam, [])
                       if s != tx and s in exon_map]
            d = decisive_support(bam, chrom, exons, sisters)
            status = lib_eval.classify_authenticity(
                d["n_decisive"], d["n_tied"], args.k_decisive, args.k_tied)
            prim = primary_support(args.bam, chrom, exons)
            out_rows.append((tx, fam, d["n_decisive"], d["n_tied"], d["n_foreign"],
                             prim, status, "true" if status == lib_eval.AUTHENTIC else "false"))
        out_rows.sort(key=lambda x: (x[1], x[0]))
        for row in out_rows:
            w.writerow(row)


if __name__ == "__main__":
    main()
