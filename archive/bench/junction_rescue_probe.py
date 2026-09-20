#!/usr/bin/env python3
"""junction_rescue_probe.py — measure the INFO-THEORETIC CEILING of the junction-axis
copy-resolution lever: of the ~31% no-PSV-coverage reads (full-length reads spanning
only shared exons, currently `unassignable_no_psv`), how many span a COPY-DISTINGUISHING
JUNCTION (an intron boundary where >=2 copies of the family differ — copy-private exon
boundary or shifted splice site).

This is the smallest experiment that validates or kills the lever:
  - For each resolvable family (>=2 copies, has PSVs so the family IS multi-haplotype),
    build per-copy intron-boundary sets from the COPY MODELS (donor/acceptor coords in
    the backbone frame), exactly like the PSV columns are built from copy alleles.
  - A junction is COPY-DISTINGUISHING iff some copy has it and >=1 other copy lacks it
    (within tolerance) -> the binary junction-as-column of copy_assignment_theory.md §6.
  - Re-thread each read. Partition the no-PSV reads into:
      rescuable_decisive : spans >=1 distinguishing junction AND argmin-unique copy
                           (the read's junction set matches exactly one copy best)
      rescuable_ambiguous: spans a distinguishing junction but ties >=2 copies
      lost_no_distinguishing_junction : spans only SHARED junctions / single-exon
                           -> info-theoretically lost on this axis (the honest ceiling)

Output: bench/junction_rescue_probe.json  (+ stdout summary)
Run    : /home/juanfra/miniforge3/bin/python bench/junction_rescue_probe.py
"""
import collections
import json
import os
import subprocess
import sys
import traceback

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
TMP = "/home/juanfra/winloci_scratch/jrescue"
MIN_READS_PSV = 3
PAD = 2000
TOL_FRAC = 0.15
OVERLAP_MERGE = 0.5
READ_CAP = 1500         # bound per-family compute (junction ceiling is read-fraction, cap is unbiased)
JTOL = 8                # junction-boundary match tolerance (bp) for splice-site jitter
MIN_J_READS = 3         # a distinguishing junction must recur across reads (not a 1-read artifact)


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def copy_column_alleles(bampath, ref, positions=None):
    """allele per (pos,name). If `positions` given, restrict the pileup to those columns
    (cheap: only the PSV-candidate columns of the copy alignment, not the whole backbone)."""
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    for col in bam.pileup(ref, min_base_quality=0, truncate=False, max_depth=300000):
        for pr in col.pileups:
            if pr.query_position is not None:
                out[col.reference_pos][pr.alignment.query_name] = \
                    pr.alignment.query_sequence[pr.query_position]
    bam.close()
    return out


def read_base_at(bampath, ref, positions):
    """For the READS bam, return {pos: {name: base}} restricted to `positions` (the PSV
    columns). Pileup is restricted by intersecting columns -> avoids the full-backbone sweep."""
    posset = set(positions)
    out = collections.defaultdict(dict)
    if not posset:
        return out
    bam = pysam.AlignmentFile(bampath, "rb")
    lo, hi = min(posset), max(posset) + 1
    for col in bam.pileup(ref, lo, hi, min_base_quality=0, truncate=True, max_depth=300000):
        if col.reference_pos not in posset:
            continue
        for pr in col.pileups:
            if pr.query_position is not None:
                out[col.reference_pos][pr.alignment.query_name] = \
                    pr.alignment.query_sequence[pr.query_position]
    bam.close()
    return out


def introns_of(bampath, ref):
    """name -> sorted set of (donor,acceptor) intron boundaries in backbone frame."""
    out = {}
    bam = pysam.AlignmentFile(bampath, "rb")
    for aln in bam.fetch(ref):
        if aln.is_unmapped or aln.is_supplementary:
            continue
        js = []
        rpos = aln.reference_start
        for op, length in aln.cigartuples:
            if op == 3:
                js.append((rpos, rpos + length)); rpos += length
            elif op in (0, 2, 7, 8):
                rpos += length
        out.setdefault(aln.query_name, set()).update(js)
    bam.close()
    return out


def dedup_copies(loci):
    by_chrom = collections.defaultdict(list)
    for lid, c, s, e, nr in loci:
        by_chrom[c].append((s, e, nr))
    regions = []
    for c, items in by_chrom.items():
        items.sort()
        cur = None
        for s, e, nr in items:
            if cur is not None and s < cur[1]:
                ov = min(cur[1], e) - s
                if ov > OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                    cur = (cur[0], max(cur[1], e), cur[2] + nr); continue
            if cur is not None:
                regions.append((c, cur[0], cur[1], cur[2]))
            cur = (s, e, nr)
        if cur is not None:
            regions.append((c, cur[0], cur[1], cur[2]))
    return regions


def jmatch(j, jset):
    return any(abs(j[0] - x[0]) <= JTOL and abs(j[1] - x[1]) <= JTOL for x in jset)


def run_family(fam_id, regions, genome, contig_len):
    os.makedirs(TMP, exist_ok=True)
    regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
    names = [f"copy{i}" for i in range(len(regions))]
    bc, bs, be, _ = regions[0]
    bs, be = max(0, bs), min(contig_len.get(bc, be), be)
    refseq = genome.fetch(bc, bs, be).upper()
    if len(refseq) < 200:
        return None
    open(f"{TMP}/ref.fa", "w").write(f">{names[0]}\n{refseq}\n")
    sh(f"samtools faidx {TMP}/ref.fa")
    with open(f"{TMP}/copies.fa", "w") as o:
        for nm, (c, s, e, _) in list(zip(names, regions))[1:]:
            s2, e2 = max(0, s), min(contig_len.get(c, e), e)
            o.write(f">{nm}\n{genome.fetch(c, s2, e2).upper()}\n")
    # genomic-span copy alignment carries INTRONS -> gives us copy intron-boundary sets
    sh(f"minimap2 -ax splice {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")

    bam = pysam.AlignmentFile(BAM, "rb")
    best = {}
    seqs = {}
    for i, (c, s, e, _) in enumerate(regions):
        s2, e2 = max(0, s - PAD), e + PAD
        for r in bam.fetch(c, s2, min(contig_len.get(c, e2), e2)):
            if r.is_unmapped or r.is_supplementary:
                continue
            if r.reference_end is None or r.reference_end < s or r.reference_start > e:
                continue
            de = dict(r.get_tags()).get("de")
            if de is None:
                continue
            if r.query_name not in best or de < best[r.query_name][1]:
                best[r.query_name] = (i, de)
            if not r.is_secondary and r.query_sequence and r.query_name not in seqs:
                seqs[r.query_name] = r.query_sequence
    bam.close()
    if len(seqs) > READ_CAP:
        seqs = dict(list(seqs.items())[:READ_CAP])
    with open(f"{TMP}/reads.fa", "w") as o:
        for n, sq in seqs.items():
            o.write(f">{n}\n{sq}\n")
    # splice:hq keeps the reads' introns -> we can read each read's junctions in backbone frame
    sh(f"minimap2 -ax splice:hq {TMP}/ref.fa {TMP}/reads.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/reads.bam - && samtools index {TMP}/reads.bam")

    # copy alleles: only a handful of copies -> cheap full pileup of copies.bam
    copy_cols = copy_column_alleles(f"{TMP}/copies.bam", names[0])

    # candidate PSV positions = backbone columns where copies carry >=2 alleles (no read
    # cost yet). Then read the reads' bases ONLY at those columns (restricted pileup).
    cand = []
    cand_hap = {}
    for p, cc in copy_cols.items():
        hap = {names[0]: refseq[p]} if p < len(refseq) else {}
        for nm in names[1:]:
            if nm in cc:
                hap[nm] = cc[nm]
        if len(hap) == len(names) and len(set(hap.values())) >= 2:
            cand.append(p); cand_hap[p] = hap
    read_cols = read_base_at(f"{TMP}/reads.bam", names[0], cand)

    # PSV bubbles = candidate columns that are read-supported (recurrence test)
    psvs = []
    for p in cand:
        if len(read_cols.get(p, {})) >= MIN_READS_PSV:
            psvs.append((p, cand_hap[p]))
    if len(psvs) == 0:        # not a resolvable family on PSV axis; skip (frontier)
        return dict(family=fam_id, resolvable=False)

    # per-copy intron-boundary sets (backbone frame), from the copy models' splice alignment
    copy_introns = introns_of(f"{TMP}/copies.bam", names[0])
    # backbone copy0 introns: from its own self ... copy0 IS the ref, so 0 introns in ref frame.
    # We instead define each copy's junctions as the union over its READS at that locus, NOT the
    # copy-model alignment, to be faithful to expressed structure. Build per-copy read-junction
    # consensus: a junction belongs to copy i iff >=MIN_J_READS reads whose best=i carry it.
    read_introns = introns_of(f"{TMP}/reads.bam", names[0])
    by_copy_j = collections.defaultdict(collections.Counter)
    for n, (ci, _de) in best.items():
        for j in read_introns.get(n, ()):
            by_copy_j[ci][j] += 1
    copy_jset = {}
    for ci in range(len(names)):
        copy_jset[ci] = {j for j, c in by_copy_j[ci].items() if c >= MIN_J_READS}

    # DISTINGUISHING junctions: present in some copy, absent (within tol) in >=1 other copy.
    all_j = set()
    for ci in copy_jset:
        all_j |= copy_jset[ci]
    distinguishing = []
    for j in all_j:
        present = [ci for ci in range(len(names)) if jmatch(j, copy_jset[ci])]
        if 0 < len(present) < len(names):
            distinguishing.append((j, set(present)))

    # Now classify the NO-PSV reads (the target population)
    cats = collections.Counter()
    # validation: among junction-rescued reads, does the junction-chosen copy == best-mapping copy?
    rescue_agree = rescue_tot = 0
    for n in best:
        # does this read cover any PSV column?
        covers_psv = any(read_cols.get(p, {}).get(n) is not None for (p, _h) in psvs)
        if covers_psv:
            cats["has_psv"] += 1
            continue
        # NO-PSV read: the target. What junctions does it carry?
        rjs = read_introns.get(n, set())
        if not rjs:
            cats["noPSV_single_exon"] += 1     # no splice signal at all -> lost on this axis
            continue
        # which distinguishing junctions does it span, and which copies do they point to?
        copy_support = collections.Counter()
        copy_against = collections.Counter()
        spanned_distinguishing = 0
        for (j, present) in distinguishing:
            if jmatch(j, rjs):
                spanned_distinguishing += 1
                for ci in present:
                    copy_support[ci] += 1
                for ci in range(len(names)):
                    if ci not in present:
                        copy_against[ci] += 1
        if spanned_distinguishing == 0:
            cats["noPSV_shared_junctions_only"] += 1   # spans junctions but none distinguishing -> LOST
            continue
        # score = (#distinguishing junctions the read shares with copy) - (# it contradicts)
        score = {ci: copy_support[ci] - copy_against[ci] for ci in range(len(names))}
        mx = max(score.values())
        winners = [ci for ci in range(len(names)) if score[ci] == mx]
        if len(winners) == 1:
            cats["noPSV_rescuable_decisive"] += 1
            rescue_tot += 1
            rescue_agree += int(best[n][0] == winners[0])
        else:
            cats["noPSV_rescuable_ambiguous"] += 1

    return dict(family=fam_id, resolvable=True, n_copies=len(names), n_psvs=len(psvs),
                n_distinguishing_junctions=len(distinguishing),
                cats=dict(cats), rescue_agree=rescue_agree, rescue_tot=rescue_tot)


def main():
    fams = collections.defaultdict(list)
    with open(FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))

    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))

    results = []
    only = set(int(x) for x in sys.argv[1:]) if len(sys.argv) > 1 else None
    fam_ids = sorted(fams)
    for k, fi in enumerate(fam_ids):
        if only is not None and fi not in only:
            continue
        regions = dedup_copies(fams[fi])
        if len(regions) < 2:
            continue
        try:
            r = run_family(fi, regions, genome, contig_len)
            if r is not None:
                results.append(r)
        except Exception as ex:
            sys.stderr.write(f"[fam {fi}] {ex}\n")
        if only is None and (k + 1) % 20 == 0:
            print(f"  ... {k+1}/{len(fam_ids)} ({len(results)} done)", flush=True)
    genome.close()

    res = [r for r in results if r.get("resolvable")]
    agg = collections.Counter()
    for r in res:
        for kk, vv in r["cats"].items():
            agg[kk] += vv
    ra = sum(r.get("rescue_agree", 0) for r in res)
    rt = sum(r.get("rescue_tot", 0) for r in res)

    noPSV_total = (agg["noPSV_single_exon"] + agg["noPSV_shared_junctions_only"]
                   + agg["noPSV_rescuable_decisive"] + agg["noPSV_rescuable_ambiguous"])
    summary = dict(
        resolvable_families=len(res),
        noPSV_reads_in_resolvable_families=noPSV_total,
        noPSV_single_exon=agg["noPSV_single_exon"],
        noPSV_shared_junctions_only=agg["noPSV_shared_junctions_only"],
        noPSV_rescuable_decisive=agg["noPSV_rescuable_decisive"],
        noPSV_rescuable_ambiguous=agg["noPSV_rescuable_ambiguous"],
        pct_decisive_of_noPSV=round(100 * agg["noPSV_rescuable_decisive"] / max(noPSV_total, 1), 1),
        pct_lost_of_noPSV=round(100 * (agg["noPSV_single_exon"] + agg["noPSV_shared_junctions_only"])
                                / max(noPSV_total, 1), 1),
        rescue_validation=f"{ra}/{rt}",
        rescue_validation_pct=round(100 * ra / max(rt, 1), 1),
    )
    out = dict(summary=summary, families=res)
    json.dump(out, open(os.path.join(HERE, "junction_rescue_probe.json"), "w"), indent=2)
    print("\n=== JUNCTION-AXIS RESCUE CEILING (no-PSV reads in resolvable families) ===")
    for kk, vv in summary.items():
        print(f"  {kk:38s}: {vv}")
    print("[+] wrote junction_rescue_probe.json")


if __name__ == "__main__":
    main()
