#!/usr/bin/env python3
"""Compare copy-ASSIGNABILITY on the SAME reads aligned by minimap2 vs winnowmap to the same
NC_073247.2 sub-reference. Isolates the aligner's effect (repeat-aware placement) on PSV-based
copy assignment. Uses the production assignment engine (copy_assign.py::assign_family).

Metrics per family: reads placed, MAPQ=0 fraction, soft-clip, PSV columns, resolvable%, assigned%
(PSV+junction), and unique-mapper agreement (PSV-assigned copy == best-overlap copy, MQ>0)."""
import sys, json, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import copy_assign as CA
import pysam

WT = "/home/juanfra/winloci_scratch/win_test"
fam_info = json.load(open(f"{WT}/fam_info.json"))
spliced = CA.load_spliced()
meta = CA.load_meta(); skel = CA.load_skel()
# rescued copies were merged into meta/skel/spliced when fam_info was built; reload to match
resc = CA.load_rescued(meta, skel, spliced)


def placement_stats(bam, contig, lo, hi, copy_spans):
    """raw placement over the family region: distinct reads, primaries, MAPQ0 frac, mean soft-clip frac."""
    b = pysam.AlignmentFile(bam, "rb")
    reads = set(); prim = 0; mq0 = 0; clip_fracs = []
    placed_on_copy = 0
    for aln in b.fetch(contig, lo, hi):
        if aln.is_unmapped: continue
        reads.add(aln.query_name)
        if not aln.is_secondary and not aln.is_supplementary:
            prim += 1
            if aln.mapping_quality == 0: mq0 += 1
            # soft-clip fraction
            ql = aln.infer_read_length() or (aln.query_length or 1)
            sc = sum(l for op, l in (aln.cigartuples or []) if op in (4, 5))
            clip_fracs.append(sc / ql if ql else 0)
        if any(min(aln.reference_end, e) - max(aln.reference_start, s) > 0 for s, e in copy_spans):
            placed_on_copy += 1
    mean_clip = sum(clip_fracs) / len(clip_fracs) if clip_fracs else float("nan")
    return dict(distinct_reads=len(reads), primaries=prim, mq0=mq0,
                mq0_frac=(mq0 / prim if prim else float("nan")),
                mean_softclip_frac=mean_clip, alns_on_copy=placed_on_copy)


def run():
    print(f"{'family':10s} {'aligner':9s} | {'reads':>5} {'prim':>5} {'MQ0%':>5} {'clip%':>6} "
          f"| {'PSVc':>4} {'resolv%':>7} {'asgn%':>6} {'uniq_agree':>12} {'uniq':>5}")
    summary = {}
    for fid, info in fam_info.items():
        contig = info["contig"]; lo = info["lo"]; hi = info["hi"]
        copies = [(c["tid"], contig, c["start"] + 1, c["end"], c["strand"]) for c in info["copies"]]
        copy_spans = [(c["start"], c["end"]) for c in info["copies"]]
        summary[fid] = {}
        for aligner, bam in (("minimap2", f"{WT}/mm2.bam"), ("winnowmap", f"{WT}/wm.bam")):
            ps = placement_stats(bam, contig, lo, hi, copy_spans)
            b = pysam.AlignmentFile(bam, "rb")
            st = CA.assign_family(spliced, b, copies, meta, skel)
            n = st["n"]; pc = st["psv_cols"]
            resolv = 100 * st["resolvable_j"] / n if n else 0
            asgn = 100 * st["assigned_j"] / n if n else 0
            agree = (st["uniq_agree_j"] / st["uniq_j"]) if st["uniq_j"] else float("nan")
            summary[fid][aligner] = dict(placement=ps, assign=st,
                                         resolv_pct=resolv, asgn_pct=asgn, uniq_agreement=agree)
            print(f"{fid:10s} {aligner:9s} | {ps['distinct_reads']:>5} {ps['primaries']:>5} "
                  f"{100*ps['mq0_frac']:>4.0f}% {100*ps['mean_softclip_frac']:>5.1f}% "
                  f"| {pc:>4} {resolv:>6.1f}% {asgn:>5.1f}% {agree:>12.4f} {st['uniq_j']:>5}")
        print()
    json.dump(summary, open(f"{WT}/winnowmap_vs_minimap2_summary.json", "w"), indent=1, default=str)
    print(f"wrote {WT}/winnowmap_vs_minimap2_summary.json")


if __name__ == "__main__":
    run()
