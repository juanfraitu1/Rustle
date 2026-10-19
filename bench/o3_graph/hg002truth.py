#!/usr/bin/env python3
"""A truth for AS-tied reads from the donor's own genome (2026-10-18; user: "we need to prove it in that way ... it just needs to
assign reads to copies based on PSVs for tied AS reads"). HG002: PacBio Kinnex reads and the donor's diploid T2T assembly
(hg002v1.1, maternal and paternal haplotypes), with chains to CHM13 v2.0.

A read is TIED in CHM13 when its best alignment score is shared by alignments on two or more copies (copy_assign's exact AS
tie); copies by exon overlap (realreads.copy_exons). Its truth comes from its alignments to the donor's own haplotypes: the
best-scoring ones, lifted to CHM13 through the chain, name the copy --
  unique          every best placement (either haplotype: the copy's two alleles) lifts onto one CHM13 copy -> the truth
  own_tie         they lift onto two or more copies: the donor's copies tie too -> undecidable from the read
  elsewhere       they lift onto CHM13, but onto no copy of the family: the read comes from another CHM13 locus (whose CHM13
                  allele, unlike the donor's, loses to the tied copies) -> never a copy to assign
  no_counterpart  none lifts onto CHM13 at all: a copy of the donor's that CHM13 lacks -> O3's case, never a copy to assign
  unmapped        no alignment to the donor's genome

    hg002truth.py tied BAM FAMILIES.json SET OUT.tsv           -> the reads tied in BAM (CHM13 coordinates): read, AS, copies,
                                                                  and minimap2's primary (its coin toss among them) + MAPQ
    hg002truth.py truth TIED.tsv MAT.bam PAT.bam CHAIN FAMILIES.json SET OUT.tsv -> each tied read's class and truth copy
"""
import collections
import csv
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("o3_graph_realreads_h", os.path.join(HERE, "realreads.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)


def _copy_of(chrom, blocks, idx):
    return R.truth_copy(chrom, blocks, 0, idx, min_mapq=0)


def tied_copies(alns, idx):
    """(best AS, the copies its best-scoring alignments lie on) of one read's alignments [(chrom, blocks, AS)]; tied when two
    or more."""
    if not alns:
        return None, set()
    best = max(a[2] for a in alns)
    return best, {c for ch, bl, sc in alns if sc == best for c in [_copy_of(ch, bl, idx)] if c}


OFF = "off"                                                   # a best placement lifting onto CHM13 outside every copy


def own_truth(own, lift, idx):
    """(class, truth copy) of a read from its alignments to the donor's haplotypes [(contig, blocks, AS)]; lift(contig, s,
    e) -> [(chrom, s, e)] on CHM13."""
    if not own:
        return "unmapped", None
    best = max(a[2] for a in own)
    found = set()
    for contig, bl, sc in own:
        if sc != best:
            continue
        pieces = [(ch, s2, e2) for s, e in bl for ch, s2, e2 in lift(contig, s, e)]
        cps = collections.Counter()
        for ch, s2, e2 in pieces:
            c = _copy_of(ch, [(s2, e2)], idx)
            if c:
                cps[c] += e2 - s2
        found.add(max(sorted(cps), key=lambda c: cps[c]) if cps else (OFF if pieces else None))
    real = found - {None, OFF}
    if not real:
        return ("elsewhere" if OFF in found else "no_counterpart"), None
    if len(found) > 1:
        return "own_tie", None
    return "unique", next(iter(real))


def _alns(bam_path, primaries=None):
    """{read: [(contig, blocks, AS)]} of every primary and secondary alignment; primaries, a dict, gets {read: (contig,
    blocks, MAPQ)} of each primary."""
    import pysam
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            out[r.query_name].append((r.reference_name, r.get_blocks(), r.get_tag("AS")))
            if primaries is not None and not r.is_secondary:
                primaries[r.query_name] = (r.reference_name, r.get_blocks(), r.mapping_quality)
    return out


def cmd_tied(bam_p, fam_p, set_name, out_p):
    (fams,) = [json.load(open(fam_p))["sets"][set_name]]
    idx = R.copy_exons(fams)
    fam_of = {f"{f['family']}.{k}": f["family"] for f in fams for k in range(len(f["copies"]))}
    n, prim = collections.Counter(), {}
    with open(out_p, "w") as fh:
        fh.write("read\tas\tfamily\tcopies\tprimary\tmapq\n")
        for q, al in _alns(bam_p, prim).items():
            best, cps = tied_copies(al, idx)
            n["reads"] += 1
            fs = {fam_of[c] for c in cps}
            if len(cps) >= 2 and len(fs) == 1:
                ch, bl, mq = prim[q]
                fh.write(f"{q}\t{best}\t{next(iter(fs))}\t{','.join(sorted(cps))}\t{_copy_of(ch, bl, idx) or ''}\t{mq}\n")
                n["tied"] += 1
            elif len(cps) >= 2:
                n["tied_across_families"] += 1
    print(f"[hg002truth] {dict(n)} -> {out_p}", flush=True)
    return 0


def cmd_truth(tied_p, mat_p, pat_p, chain_p, fam_p, set_name, out_p):
    sys.path.insert(0, "/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/src")
    from chainlift import Chain
    ch = Chain(chain_p)
    fams = json.load(open(fam_p))["sets"][set_name]
    idx = R.copy_exons(fams)
    lift = lambda contig, s, e: [(q, a, b) for q, a, b, _ in ch.lift(contig, s, e)]  # noqa: E731
    own = _alns(mat_p)
    for q, v in _alns(pat_p).items():
        own[q] += v
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\ttied_copies\tclass\ttruth\n")
        for r in csv.DictReader(open(tied_p), delimiter="\t"):
            cls, cp = own_truth(own.get(r["read"], []), lift, idx)
            if cls == "unique" and cp not in r["copies"].split(","):
                cls = "unique_elsewhere"                    # the donor places it on a copy CHM13 did not tie it with
            n[cls] += 1
            fh.write(f"{r['read']}\t{r['family']}\t{r['copies']}\t{cls}\t{cp or ''}\n")
    print(f"[hg002truth] {dict(n)} -> {out_p}", flush=True)
    return 0


def main(argv):
    if len(argv) == 5 and argv[0] == "tied":
        return cmd_tied(*argv[1:])
    if len(argv) == 8 and argv[0] == "truth":
        return cmd_truth(*argv[1:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
