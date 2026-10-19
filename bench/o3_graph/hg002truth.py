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
    hg002truth.py score TRUTH.tsv TIED.tsv WORLD OUT.json [--o2 NAME]  -> per method (minimap2's primary, the primary at MAPQ
                                                                  >= 1, current, integrated: summary.world_calls with the O2
                                                                  table NAME, o2.integv.m2.tsv by default) and truth class,
                                                                  the fates (tied_fate) of the world's tied reads
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


METHODS = ("mm2", "mm2_q1", "current", "integrated")


def _summary():
    spec = importlib.util.spec_from_file_location("o3_graph_summary_h", os.path.join(HERE, "summary.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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


def tied_fate(cls, truth, assigned, in_o3, tied):
    """A tied read's fate under one method (an absent-copy call takes precedence over a placement, as summary.read_fate):
    unique / unique_elsewhere: right only on the donor's copy, and an O3 call is an error (its copy is in CHM13); own_tie:
    placing it on one of the copies it ties on is a guess, anywhere else or an O3 call an error, abstaining right;
    no_counterpart (a donor copy CHM13 lacks): any placement is an error, an O3 call right; elsewhere (another CHM13 locus):
    a placement on the family is an error, an O3 call is reported apart."""
    if in_o3:
        return {"no_counterpart": "correct", "elsewhere": "o3"}.get(cls, "error")
    if assigned is None:
        return "abstain"
    if cls in ("unique", "unique_elsewhere"):
        return "correct" if assigned == truth else "error"
    return "guess" if cls == "own_tie" and assigned in tied else "error"


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


def _fam_call(calls, family):
    """A system's placement of a read of `family`: None, its copy, `no-truth-copy` or `wrong-family` (as o2proof)."""
    if not calls:
        return None
    if family in calls:
        return calls[family] or "no-truth-copy"
    return "wrong-family"


def _world_reads(fq):
    with open(fq) as fh:
        return {h[1:].split()[0] for i, h in enumerate(fh) if i % 4 == 0}


def cmd_score(truth_p, tied_p, world, out_p, o2_name="o2.integv.m2.tsv"):
    wc = _summary().world_calls(world, o2_name)
    tied = {r["read"]: r for r in csv.DictReader(open(tied_p), delimiter="\t")}
    reads = _world_reads(f"{world}/reads.E0.fq")
    card = {m: collections.defaultdict(collections.Counter) for m in METHODS}
    skipped = collections.Counter()
    for r in csv.DictReader(open(truth_p), delimiter="\t"):
        q, f, cls = r["read"], r["family"], r["class"]
        if cls == "unmapped" or q not in reads:
            skipped["unmapped" if cls == "unmapped" else "not_in_world"] += 1
            continue
        t, ts = tied[q], set(r["tied_copies"].split(","))
        prim = t["primary"] or "wrong-locus"
        calls = dict(mm2=(prim, False), mm2_q1=(prim if int(t["mapq"]) >= 1 else None, False),
                     current=(_fam_call(wc["current"].get(q), f), q in wc["current_o3"]),
                     integrated=(_fam_call(wc["integrated"].get(q), f), q in wc["integrated_o3"]))
        for m, (a, o3) in calls.items():
            card[m][cls][tied_fate(cls, r["truth"] or None, a, o3, ts)] += 1
    res = dict(card={m: {c: dict(v) for c, v in card[m].items()} for m in METHODS}, skipped=dict(skipped))
    json.dump(res, open(out_p, "w"), indent=1)
    for m in METHODS:
        print(f"  {m:10s} " + " | ".join(f"{c}: " + " ".join(f"{k} {v}" for k, v in sorted(res['card'][m][c].items()))
                                        for c in sorted(res["card"][m])), flush=True)
    print(f"  skipped {dict(skipped)}", flush=True)
    return 0


def main(argv):
    if argv[:1] == ["score"] and (len(argv) == 5 or len(argv) == 7 and argv[5] == "--o2"):
        return cmd_score(*argv[1:5], o2_name=argv[6] if len(argv) == 7 else "o2.integv.m2.tsv")
    if len(argv) == 5 and argv[0] == "tied":
        return cmd_tied(*argv[1:])
    if len(argv) == 8 and argv[0] == "truth":
        return cmd_truth(*argv[1:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
