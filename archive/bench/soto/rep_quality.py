#!/usr/bin/env python3
"""PERSISTED per-locus record: how well does this representative represent its locus?

    rep_quality.py --copies <copies.tsv> --out <prefix> [--bam <bam>] [--chrom chr16]
                   [--genes NPIPB6,NPIPA2 | --gene-prefix NPIP] [--universe annotation,reads]

WHY THIS EXISTS AS A FILE AND NOT AS A NUMBER IN A NOTEBOOK. The same representative at the same locus
(NPIPB6) has been scored 0.333, 0.583, 0.625, 1.000 and 0.367 depending only on which denominator was
picked, and each of those numbers was computed once in a scratchpad and lost. This writes the COUNTS, with
their parameters, so the ratio is derived by whoever reads the file and the choice of denominator is
visible rather than assumed.

THE OBJECT. A representative is ONE LINEAR PATH. It can carry only a set of pairwise NON-OVERLAPPING
junctions, so the largest set it could ever carry is the maximum chain of the locus's junction set — the
SINGLE-PATH CEILING (rustlib.single_path_ceiling). That is an intrinsic property of the locus: it is not an
argmax over models, it cannot be inflated by lengthening the representative, and it needs no threshold.
The dual number, the minimum number of linear paths needed to carry ALL of J (rustlib.min_path_cover, =
max antichain by Dilworth), is emitted alongside because it is arguably the sharper statement: NPIPB6
needs 5, so one representative is at best a fifth of it.

⚠⚠ THIS FILE IS A MEASUREMENT. NOTHING IN THE PIPELINE MAY CONSUME IT TO DECIDE WHAT A LOCUS IS. The
  matching inside min_path_cover measures an already-built representative; O1 defines families by sequence
  homology alone. If a bipartite matching ever influences what a locus or a family IS, the method is
  broken. This script only READS copies.tsv; it writes nothing the engine reads.

⚠ THE SCORE IS A **SATURATION**, NEVER A PRECISION / ACCURACY / RECALL. "Junction recall against the
  locus's own >=3-read set" was killed as tautological: the rep is BUILT from those reads, so its
  junctions are a subset of the truth by construction and no false junction can ever register. Dividing by
  the ceiling instead of |J| RESCALES that quantity, it does not repair it. What makes the record
  non-vacuous is that the CEILING is a property of the locus and moves independently of the model — hence
  `on_ceiling / ceiling` is "of what one linear path could carry here, how much did it carry".

⚠ THE COUNTS ARE THE RECORD; the two ratios are derived columns and are written as explicit NA whenever
  their denominator is 0. A locus with no junction clearing the floor has ceiling 0, and 0/0 must never be
  written out as a 0.0 that reads like a real score. Because the counts are kept, any OTHER denominator
  (a curated transcript, a co-observation-constrained ceiling) can be re-derived without re-running.

TWO JUNCTION UNIVERSES, NEVER MIXED IN ONE ROW (a `universe` column says which):
  annotation  the union of the overlapping RefSeq gene's transcript introns
              (rustlib.introns_of_transcripts, min_len=20). This is the PAPER number and the one that
              reproduces the recorded NPIPB6 12/21. ⚠ Annotation is a HYPOTHESIS: a locus with ONE
              annotated transcript necessarily has ceiling == |J|, which says nothing about the locus and
              everything about RefSeq. `n_tx` is emitted so that is auditable.
  reads       the locus's read-supported junctions at >= --min-reads, from the BAM (rustlib.chains_in,
              primary only). Annotation-free, denser, and a DIFFERENT OBJECT — the recorded six pairs are
              annotation numbers and must not be restated from these.

⚠ THE `read_exons` GAP ROUTE IS BANNED. `read_exons` returns MERGED blocks, so its gaps are pairwise
  disjoint by construction, ceiling == |J| at every locus, and every score would be 1.000. There is a
  regression test asserting that degeneracy (test_rustlib_parity.py, section (d)).

⚠ THE CEILING IS ONLY AS GOOD AS THE LOCUS BOUNDARY. Boundaries are the known open problem here (fixing
  locus sizes moved recall 0.698 -> 0.788 but precision 0.709 -> 0.269). `rep_start/rep_end`, `n_copies`
  and `rep_off_universe` are emitted so an over-extended representative importing a neighbour's junctions
  is visible rather than absorbed into the score.

⚠ THE CEILING IS PURELY GEOMETRIC AND THEREFORE GENEROUS. It ignores read co-observation: two junctions
  can be geometrically compatible yet never co-observed on any read. The achievable single-path maximum is
  <= ceiling, sometimes strictly. Calling the gap "algorithmic shortfall" would overstate it.

⚠⚠ DO NOT POINT THIS AT A `--from-genome` CATALOG. DNA-substrate reps are genomic intervals with an
  EMPTY intron chain BY DESIGN (from_genome.rs asserts it), so every junction field would be a structural
  NA reported as a hard 0 — "the model explains none of its locus" when the truth is that the substrate
  has no junctions to explain. The record is meaningful only on a spliced (RNA) catalog.

⚠⚠ CHECK HOW THE CATALOG'S BAM WAS BUILT BEFORE QUOTING ANY ROW. A `samtools view -M -L` subset BAM
  yields loci only where its BED had reads, and that has caused three retractions in this project. The
  header of the emitted file records `--copies`, `--bam` and every parameter; it cannot record what BAM the
  catalog was built from, so state it yourself.

Stdlib only, via rustlib. Reads copies.tsv; never writes anything the engine reads.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rustlib as R  # noqa: E402

COLUMNS = [
    "locus",             # gene symbol (the SEED naming the locus; annotation is a hypothesis, not truth)
    "node_key",          # L~{chrom}_{start}_{end} of the REP — joins to copies.tsv and the E_r edge dump
    "universe",          # annotation | reads  — which junction set the denominators came from
    "chrom", "gene_start", "gene_end",
    "n_junctions_union",  # |J|            UNFAIR alone: mixes mutually exclusive isoforms
    "ceiling",            # max chain      THE fair structural denominator
    "min_paths",          # max antichain  reps needed to carry ALL of J (Dilworth dual)
    "n_junctions_rep",    # |rep chain|    the rep's own junctions, ALL of them
    "rep_on_universe",    # |rep ∩ J|      the numerator: same universe as the denominator
    "rep_off_universe",   # |rep \ J|      rep junctions the universe does not contain (boundary/novel)
    "raw",                # rep_on_universe / n_junctions_union  — the UNFAIR ratio, kept for audit
    "normalised",         # rep_on_universe / ceiling            — the SATURATION. NA when ceiling == 0.
    "n_tx",               # annotated transcripts at the locus (0 => ceiling==|J| is an artefact)
    "locus_reads",        # primary reads on chains CONTAINED in the gene span — the same population the
                          # `reads` universe is built from, so the two are consistent. A thin locus is not
                          # a badly-represented one, and this is how you tell them apart. ⚠ It is NOT the
                          # total read pileup: chains merely overlapping the span are excluded on purpose
                          # (the neighbour leak that once made an adjacent transcript NPIPB9's rep).
    "rep_reads",          # the representative's own read support
    "n_copies",           # catalog copies overlapping the locus (>1 => the locus was split)
    "rep_start", "rep_end", "rep_strand", "rep_tid",
    "status",             # ok | no_rep | no_junctions
]

# Parameters that move every number in the file. Emitted as `#param` lines, following the .params.tsv
# discipline of the E_r edge dump — a record whose knobs are not recorded is not reproducible.
PARAM_KEYS = ["copies", "bam", "chrom", "genes", "gene_prefix", "universe", "min_reads",
              "min_intron_len", "touch_ok", "contained", "gff", "note"]

NA = "NA"


def rep_junctions(row):
    """The representative's intron chain, from copies.tsv's `exons` column (0-based, half-open).

    ⚠ TAUTOLOGY GUARD: this is the NUMERATOR only. It must never become the denominator — "fraction of the
      rep's own junctions that are read-supported" is identically 1 because the rep is built from those
      reads. Two metrics have already been killed for exactly that shape.
    """
    ex = row.get("exons")
    if not ex:
        return None
    blocks = []
    for b in ex.split(","):
        if "-" not in b:
            return None
        s, e = b.split("-", 1)
        blocks.append((int(s), int(e)))
    blocks.sort()
    return [(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1)]


def read_universe(bam, chrom, start, end, min_reads, min_len, contained=True):
    """Read-supported junctions at >= min_reads, and the locus's primary read total.

    ⚠ A region query returns reads OVERLAPPING the region, so a chain can carry junctions far outside it
      (measured: 43-44% do). `contained` keeps only chains lying inside the gene span — the same rule
      `read_exons(contained=True)` applies, and the reason it exists is that the neighbouring locus's
      transcripts once became NPIPB9's "representative".
    """
    chains = R.chains_in(bam, chrom, start + 1, end)
    if contained:
        chains = [c for c in chains
                  if c["start"] is not None and c["start"] >= start and c["end"] <= end]
    depth = {}
    total = 0
    for c in chains:
        total += c["n_reads"]
        for j in c["introns"]:
            if j[1] - j[0] < min_len:
                continue
            depth[j] = depth.get(j, 0) + c["n_reads"]
    return sorted(j for j, n in depth.items() if n >= min_reads), total


def build_rows(args):
    genes = [(s, e, n) for s, e, n in R.genes_on(args.chrom)]
    if args.genes:
        want = set(args.genes.split(","))
        genes = [g for g in genes if g[2] in want]
    elif args.gene_prefix:
        genes = [g for g in genes if g[2].startswith(args.gene_prefix)]
    copies = [c for c in R.load_copies(args.copies, require_exons=True) if c["chrom"] == args.chrom]

    rows = []
    for gs, ge, name in sorted(genes, key=lambda g: g[0]):
        ov = [c for c in copies if c["start"] < ge and c["end"] > gs]
        # ⚠ LARGEST overlap, never first — first-overlap has produced a plausible wrong number three
        #   times in this project (NOTCH2NLB, AMY1, NPIPB9). rustlib.best_overlap is the one rule.
        bi = R.best_overlap(gs, ge, [(c["start"], c["end"]) for c in ov]) if ov else None
        rep = ov[bi] if bi is not None else None
        rj = rep_junctions(rep) if rep is not None else None

        tx = R.introns_of_transcripts(args.chrom, name, min_len=args.min_intron_len)
        ann = sorted({iv for t in tx.values() for iv in t})
        # The locus's read total is a property of the LOCUS, not of the universe, and it is what tells a
        # thin locus apart from a badly-represented one. Compute it once and put it on every row.
        rj_reads, locus_reads = ((None, NA) if not args.bam else
                                 read_universe(args.bam, args.chrom, gs, ge,
                                               args.min_reads, args.min_intron_len, args.contained))

        for uni in args.universe.split(","):
            if uni == "annotation":
                J = ann
            elif uni == "reads":
                if rj_reads is None:
                    continue
                J = rj_reads
            else:
                raise ValueError(f"unknown universe {uni!r}")

            cl = R.locus_ceiling(J, touch_ok=args.touch_ok) if J else None
            js = set(J)
            row = {c: NA for c in COLUMNS}
            row.update(locus=name, universe=uni, chrom=args.chrom,
                       gene_start=gs, gene_end=ge, n_tx=len(tx),
                       locus_reads=locus_reads, n_copies=len(ov))
            if cl is not None:
                row.update(n_junctions_union=cl["n_junctions"], ceiling=cl["ceiling"],
                           min_paths=cl["min_paths"])
            if rep is not None:
                row.update(node_key=R.canon_node(f"{args.chrom}:{rep['start'] + 1}-{rep['end']}"),
                           rep_start=rep["start"], rep_end=rep["end"],
                           rep_strand=rep.get("strand", NA), rep_tid=rep.get("tid", NA),
                           rep_reads=rep.get("n_reads", NA))
            if rj is not None:
                row.update(n_junctions_rep=len(rj),
                           rep_on_universe=sum(1 for j in rj if j in js),
                           rep_off_universe=sum(1 for j in rj if j not in js))
            # The two ratios the brief asks for, DERIVED from the counts above and written only when the
            # denominator is a real number. ⚠ 0/0 must never appear as 0.0: a locus with no junction
            # clearing the floor is NA, not "scored zero". The counts stay in the file so any other
            # denominator can be re-derived without re-running anything.
            if row["rep_on_universe"] != NA:
                if row["n_junctions_union"] not in (NA, 0):
                    row["raw"] = f"{row['rep_on_universe'] / row['n_junctions_union']:.4f}"
                if row["ceiling"] not in (NA, 0):
                    row["normalised"] = f"{row['rep_on_universe'] / row['ceiling']:.4f}"
            # An explicit status, so a 0 is never read as a score. `no_rep` means the catalog holds no
            # copy overlapping the locus — which is NOT the same as "no representative was built": copies
            # that never joined a family, or fell below --min-copies, never reach copies.tsv at all.
            row["status"] = ("no_rep" if rep is None else
                             "no_junctions" if cl is None else "ok")
            rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True, help="a shipped copies.tsv WITH the `exons` column")
    ap.add_argument("--out", required=True, help="output prefix; writes <out>.rep_quality.tsv")
    ap.add_argument("--bam", default=None, help="BAM for the `reads` universe and locus read totals")
    ap.add_argument("--chrom", default="chr16")
    ap.add_argument("--genes", default=None, help="comma-separated gene symbols")
    ap.add_argument("--gene-prefix", default=None, help="e.g. NPIP")
    ap.add_argument("--universe", default="annotation",
                    help="comma-separated: annotation,reads (emitted as separate ROWS, never mixed)")
    ap.add_argument("--min-reads", type=int, default=3,
                    help="junction read floor for the `reads` universe")
    ap.add_argument("--min-intron-len", type=int, default=20,
                    help="drops sub-20bp RefSeq transcript-to-genome ALIGNMENT GAPS, which are not "
                         "splices (21 of 229 across the NPIP loci). This parameter MOVES the numbers: "
                         "NPIPB2 is 9/14 at 20 and 10/15 at 0.")
    ap.add_argument("--strict-exon", dest="touch_ok", action="store_false", default=True,
                    help="require >= 1 exonic base between two introns (a zero-length exon is not an "
                         "object). MEASURED: changes no ceiling at any of the 19 chr16 NPIP loci.")
    ap.add_argument("--overlapping", dest="contained", action="store_false", default=True,
                    help="`reads` universe: keep chains merely OVERLAPPING the locus (the neighbour leak)")
    ap.add_argument("--note", default="",
                    help="free text recorded as a #param: state the REP CONSTRUCTION and the BAM the "
                         "catalog was built from. This script cannot infer either from copies.tsv, and "
                         "both move every number: at NPIPB6 the SAME rep DN_chr16_28623009_8 carries 7 "
                         "junctions under `pick_locus_rep` and 8 under RUSTLE_COTHREAD_REP.")
    args = ap.parse_args()

    rows = build_rows(args)
    path = f"{args.out}.rep_quality.tsv"
    vals = vars(args)
    with open(path, "w") as fh:
        fh.write("#record\tper-locus single-path saturation of the representative — A MEASUREMENT.\n")
        fh.write("#warn\traw = rep_on_universe/n_junctions_union (UNFAIR: the union mixes mutually "
                 "exclusive isoforms). normalised = rep_on_universe/ceiling. Both are NA when the "
                 "denominator is 0 — a 0/0 locus is never scored 0.0.\n")
        fh.write("#warn\tthe normalised score is a SATURATION, not a precision/accuracy/recall.\n")
        for k in PARAM_KEYS:
            v = vals.get(k, R.GFF if k == "gff" else None)
            fh.write(f"#param\t{k}\t{v}\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    print(f"wrote {path}  ({len(rows)} rows)")

    # Console table. The ratios are computed HERE, from the counts, and never stored.
    for uni in dict.fromkeys(args.universe.split(",")):
        sel = [r for r in rows if r["universe"] == uni]
        if not sel:
            continue
        print(f"\nuniverse = {uni}")
        print(f"{'locus':10s} {'ntx':>4s} {'|J|':>4s} {'ceil':>5s} {'nrep':>5s} {'on':>4s} {'off':>4s} "
              f"{'raw':>6s} {'norm':>6s} {'paths':>6s} {'reads':>7s} {'rep_rd':>7s} status")
        num = []
        for r in sel:
            ok = r["status"] == "ok" and r["rep_on_universe"] != NA and r["ceiling"] not in (NA, 0)
            raw = (f"{r['rep_on_universe'] / r['n_junctions_union']:.3f}"
                   if ok and r["n_junctions_union"] else NA)
            nrm = f"{r['rep_on_universe'] / r['ceiling']:.3f}" if ok else NA
            if ok:
                num.append(r["rep_on_universe"] / r["ceiling"])
            print(f"{r['locus']:10s} {str(r['n_tx']):>4s} {str(r['n_junctions_union']):>4s} "
                  f"{str(r['ceiling']):>5s} {str(r['n_junctions_rep']):>5s} "
                  f"{str(r['rep_on_universe']):>4s} {str(r['rep_off_universe']):>4s} "
                  f"{raw:>6s} {nrm:>6s} {str(r['min_paths']):>6s} "
                  f"{str(r['locus_reads']):>7s} {str(r['rep_reads']):>7s} {r['status']}")
        if num:
            num.sort()
            med = num[len(num) // 2] if len(num) % 2 else (num[len(num) // 2 - 1] + num[len(num) // 2]) / 2
            print(f"  scored {len(num)}/{len(sel)} loci; median normalised {med:.3f}; "
                  f"range {min(num):.3f}-{max(num):.3f}")


if __name__ == "__main__":
    main()
