#!/usr/bin/env python
"""Score the pipeline against the PLANTED ground truth on the fully-simulated genome (bench/sim_genome.py).

Non-circular: the read name encodes the TRUE family/copy (SIMGW|fam|copy|iso|i) and the copy loci are planted
(simgw_truth.tsv). We map each detected copy to a planted copy by genomic locus and check:
  O1 family detection — are the planted multi-copy families recovered? do single-copy/domain-sharer NOT merge?
  O2 assignment       — is each read assigned to its TRUE copy? (the real non-circular accuracy)
  K=0 floor           — are the identical-copy (K0tandem) reads certified TIED, not guessed?
  posterior/zone      — does a Tied read's consistent zone contain its true copy?

Run: /home/juanfra/miniforge3/bin/python bench/sim_eval.py
"""
import os
from collections import defaultdict

import pysam

OUT = "/home/juanfra/winloci_scratch"
W = 3000  # locus-match tolerance (bp)


def load_readinfo():
    """read_name -> (primary MAPQ, multimaps?). A read MULTIMAPS (the gate's regime) when minimap2 keeps a
    secondary alignment for it OR gives it MAPQ 0 — i.e. it places ambiguously across ≥2 copies. A read with a
    secondary but a low positive MAPQ (e.g. a collapsed segdup: few PSVs over a long copy) is STILL the gate's
    job — minimap2 cannot confidently resolve it, even though the PSV columns can."""
    mq, multi = {}, defaultdict(bool)
    b = pysam.AlignmentFile(f"{OUT}/simgw.bam", "rb")
    for a in b.fetch():
        if a.is_unmapped:
            continue
        if a.is_secondary or a.is_supplementary:
            multi[a.query_name] = True
            continue
        mq[a.query_name] = a.mapping_quality
    return {n: (mq.get(n, 0), multi[n] or mq.get(n, 0) == 0) for n in mq}


def load_truth():
    t = {}  # (family, copy) -> (chrom, start, end, divergence)
    for line in open(f"{OUT}/simgw_truth.tsv").readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        t[(f[0], int(f[1]))] = (f[2], int(f[3]), int(f[4]), float(f[6]))
    return t


def parse_tid(tid):
    """tid like DN_/RC_/AC_<chrom>_<start>[_<n>]; sim chroms are single-token (simA/simB)."""
    rest = tid
    for pre in ("DN_", "RC_", "AC_"):
        if rest.startswith(pre):
            rest = rest[len(pre):]
            break
    parts = rest.split("_")
    chrom = parts[0]
    for x in parts[1:]:
        if x.isdigit():
            return (chrom, int(x))
    return None


def main():
    truth = load_truth()
    rinfo = load_readinfo()
    # detected copy loci per (family_id, copy_index) from quant.tsv
    copy_loc = {}
    for line in open(f"{OUT}/simgw_as.quant.tsv").readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        loc = parse_tid(f[2])
        if loc:
            copy_loc[(f[0], int(f[1]))] = loc

    def match_planted(chrom, start):
        # nearest planted copy by start on the same chrom, within W (copies can be <W apart, so a window test
        # alone would assign several detected copies to the same planted copy — nearest-start disambiguates).
        best, bestd = None, W + 1
        for (fam, ci), (c, s, e, _) in truth.items():
            if c == chrom and abs(start - s) < bestd:
                best, bestd = (fam, ci), abs(start - s)
        return best

    # which planted copies were RECOVERED as a distinct detected copy (have a detected locus near them)?
    # A read whose true copy was NOT recovered CANNOT be assigned to it (the target is absent from the model) —
    # that is an O1 detection miss, not an O2 misassignment, so we score O2 conditional on recovery.
    recovered = set()
    for (fam_det, ci), loc in copy_loc.items():
        mp = match_planted(*loc)
        if mp:
            recovered.add(mp)

    # ---- O2: per-read assignment accuracy vs TRUE copy ----
    # Two regimes (the whole point of the gate): a read is genuinely AMBIGUOUS when minimap2 gives it MAPQ==0
    # (it multimaps to >1 copy) — that is the gate's job. A read with MAPQ>0 is already resolved by mapping
    # (a copy diverged enough to be uniquely mappable), so it never needs the gate. We score them separately.
    n = correct = assigned = tied = 0
    by_div = defaultdict(lambda: [0, 0, 0])  # divergence -> [correct, assigned, ambiguous reads] (recovered copies only)
    k0_total = k0_tied = 0
    uniq_resolved = uniq_total = 0           # MAPQ>0 reads correctly resolved by mapping
    o1_miss_reads = 0                        # ambiguous reads whose TRUE copy was not recovered (O1, not O2)
    dom_single_assigned_to_family = 0
    seen = set()  # a read maps to every copy locus -> appears once per placement; the call is identical, dedupe
    for line in open(f"{OUT}/simgw_as.assignments.tsv").readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        name, fam_det, copy_det, status = f[0], f[1], int(f[2]), f[3]
        parts = name.split("|")
        if parts[0] != "SIMGW":
            continue
        if name in seen:
            continue
        seen.add(name)
        true_fam, true_copy = parts[1], int(parts[2])
        ambiguous = rinfo.get(name, (0, True))[1]   # multimaps (has secondary or MAPQ0) = the gate's regime
        n += 1
        if true_fam == "K0tandem":
            k0_total += 1
            k0_tied += (status == "tied")
        if true_fam.startswith(("single", "domshare")):
            # these should NOT be in a multi-copy family at all; if they were assigned in one, flag it
            if status == "assigned":
                dom_single_assigned_to_family += 1
        det = copy_loc.get((fam_det, copy_det))
        mp = match_planted(*det) if det else None
        ok = (mp is not None and mp[0] == true_fam and mp[1] == true_copy)
        d = truth.get((true_fam, true_copy), (None, None, None, -1))[3]
        if not ambiguous:
            # uniquely-mappable copy: correct outcome is being placed on its own locus (by minimap2 or the gate)
            uniq_total += 1
            uniq_resolved += (status == "assigned" and ok)
            continue
        # multimapping ambiguous read — the assignment regime
        if true_fam != "K0tandem" and (true_fam, true_copy) not in recovered:
            o1_miss_reads += 1               # true copy absent from the model — O1 miss, not an O2 error
            continue
        if status == "assigned":
            assigned += 1
            correct += ok
            by_div[d][0] += ok
            by_div[d][1] += 1
        elif status == "tied":
            tied += 1
        by_div[d][2] += 1

    print(f"== O2: per-read assignment vs TRUE copy (non-circular; copies in the model) ==")
    print(f"  MULTIMAPPING reads (secondary/MAPQ0 — the gate's regime), true copy RECOVERED: "
          f"assigned={assigned}  tied={tied}")
    if assigned:
        print(f"  ACCURACY on assigned: {correct}/{assigned} = {100*correct/assigned:.1f}%  (= TRUE copy)")
    print(f"  by planted divergence:")
    for d in sorted(by_div):
        c, a, amb = by_div[d]
        print(f"    div={d:.3f}: assigned {c}/{a} correct" +
              (f" = {100*c/a:.0f}%" if a else "") + f"  ({amb} reads)")
    print(f"  reads whose TRUE copy was NOT recovered (O1 miss, not an O2 error): {o1_miss_reads}")
    print(f"  uniquely-mapped reads (MAPQ>0, no secondary — resolved by mapping): {uniq_resolved}/{uniq_total}")
    print(f"\n== K=0 floor (K0tandem, identical copies) ==")
    print(f"  TIED: {k0_tied}/{k0_total} = {100*k0_tied/max(k0_total,1):.0f}%  (should be ~100% — certified unresolvable, not guessed)")
    print(f"\n== FP controls ==")
    print(f"  single/domain-sharer reads wrongly assigned within a multi-copy family: {dom_single_assigned_to_family}")

    # ---- O1: family recovery (from the gw_family_catalog, which is the cross-chrom-aware family detector;
    #         copy_assign above is co-located/per-contig only, so it cannot see cross-chrom families) ----
    print(f"\n== O1: planted multi-copy families recovered? (gw_family_catalog) ==")
    planted_fams = defaultdict(set)
    for (fam, ci) in truth:
        planted_fams[fam].add(ci)
    multi = {f: cs for f, cs in planted_fams.items() if len(cs) >= 2}
    # detected copies straight from the catalog (explicit chrom/start), mapped to planted by locus
    det_planted_copies = defaultdict(lambda: defaultdict(set))   # det_family -> planted_family -> {copy idx}
    det_fam_planted = defaultdict(set)                           # det_family -> {planted families}
    for line in open(f"{OUT}/simgw_fam.copies.tsv").readlines()[1:]:
        f = line.rstrip("\n").split("\t")
        fam_det, chrom, start = f[0], f[3], int(f[4])
        mp = match_planted(chrom, start)
        if mp:
            det_planted_copies[fam_det][mp[0]].add(mp[1])
            det_fam_planted[fam_det].add(mp[0])
    for fam, cs in multi.items():
        ncov = max((len(det_planted_copies[d].get(fam, set())) for d in det_planted_copies), default=0)
        rec = ncov >= 2
        print(f"  {fam} ({len(cs)} copies): {'RECOVERED' if rec else 'MISSED'}  ({ncov}/{len(cs)} copies in one detected family)")
    # over-merge: a detected family spanning >1 planted family
    overmerge = {d: ps for d, ps in det_fam_planted.items() if len(ps) > 1}
    print(f"  over-merged detected families (span >1 planted family): {len(overmerge)}" +
          (f" -> {dict(overmerge)}" if overmerge else " (none — no domain-sharer/unrelated merge)"))


if __name__ == "__main__":
    main()
