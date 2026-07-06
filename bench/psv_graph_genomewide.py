#!/usr/bin/env python3
"""psv_graph_genomewide.py — scale the PSV-aware variation graph (psv_graph_demo.py) over
EVERY validated multi-copy gene family, to measure the genome-wide RESOLVABLE vs
K-FRONTIER split.

For each family: pick the longest copy as the backbone, align the other copies (minimap2
asm20) and the reads (splice:hq) onto it, surface the read-supported PSV bubbles (the
columns of copy_assignment_theory.md), read off each copy's path, and thread the reads.

The headline quantity is, per family, K = the number of RESOLVABLE haplotypes (distinct
copy paths through the bubbles):
  - K == #copies            -> fully resolvable  (every copy read-distinguishable)
  - 2 <= K <  #copies       -> partially resolvable (some copies collapse)
  - K == 1 (or 0 PSVs)      -> K-FRONTIER (copies exonically identical -> unresolvable)

Input : /home/juanfra/winloci_scratch/validated_families.tsv  (from family_def_vg_reinforce.py)
Output: bench/psv_graph_genomewide.json  +  bench/psv_graph_genomewide_perfamily.tsv
Run   : /home/juanfra/miniforge3/bin/python bench/psv_graph_genomewide.py
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
FAM_TSV = os.environ.get("FAM_TSV", "/home/juanfra/winloci_scratch/validated_families.tsv")
TMP = "/home/juanfra/winloci_scratch/psvgw"
MIN_READS_PSV = 3        # a PSV bubble must recur across reads (not a one-read error)
PAD = 2000
TOL_FRAC = 0.15          # HiFi-tolerant best-PSV-match
OVERLAP_MERGE = 0.5      # de-novo loci with reciprocal overlap > this = one copy (isoforms)
READ_CAP = 6000          # bound per-family compute; families above this are subsampled


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def column_alleles(bampath, ref):
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    # flag_filter 0xF04 = UNMAP|SECONDARY|QCFAIL|DUP (pysam default 0x704) PLUS SUPPLEMENTARY (0x800):
    # a re-aligned read's supplementary alignment must not contribute a paralogous base to a column
    # (the primary already carries the read; a supplementary is a chimeric/split fragment). The pysam
    # default did not exclude 0x800 -> a handful of split re-alignments could leak a base per family.
    for col in bam.pileup(ref, min_base_quality=0, truncate=False, max_depth=300000,
                          flag_filter=0xF04):
        for pr in col.pileups:
            if pr.query_position is not None:
                out[col.reference_pos][pr.alignment.query_name] = \
                    pr.alignment.query_sequence[pr.query_position]
    bam.close()
    return out


def dedup_copies(loci):
    """Collapse de-novo loci that are the same genomic copy (isoforms) into one copy-region:
    same chrom + reciprocal overlap > OVERLAP_MERGE. Distinct tandem copies stay separate."""
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
                    cur = (cur[0], max(cur[1], e), cur[2] + nr)
                    continue
            if cur is not None:
                regions.append((c, cur[0], cur[1], cur[2]))
            cur = (s, e, nr)
        if cur is not None:
            regions.append((c, cur[0], cur[1], cur[2]))
    return regions


def run_family(fam_id, regions, genome, contig_len):
    os.makedirs(TMP, exist_ok=True)
    # backbone = longest copy-region (most stable coordinate frame)
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
    sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")

    # reads from all copy loci -> best-mapping copy (lowest de) + cross-copy ties
    bam = pysam.AlignmentFile(BAM, "rb")
    best = {}
    ties = collections.defaultdict(set)
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
            ties[r.query_name].add(i)
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
    sh(f"minimap2 -ax splice:hq {TMP}/ref.fa {TMP}/reads.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/reads.bam - && samtools index {TMP}/reads.bam")

    copy_cols = column_alleles(f"{TMP}/copies.bam", names[0])
    read_cols = column_alleles(f"{TMP}/reads.bam", names[0])

    # PSV bubbles: all copies defined, >=2 distinct alleles, read-supported (recurrence)
    psvs = []
    for p in range(len(refseq)):
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            if nm in copy_cols.get(p, {}):
                hap[nm] = copy_cols[p][nm]
        if (len(hap) == len(names) and len(set(hap.values())) >= 2
                and len(read_cols.get(p, {})) >= MIN_READS_PSV):
            psvs.append((p, hap))
    hap_rows = {nm: tuple(hap[nm] for _, hap in psvs) for nm in names}
    groups = collections.defaultdict(list)
    for nm in names:
        groups[hap_rows[nm]].append(nm)
    n_groups = len(groups)

    # thread reads: best PSV-allele match (HiFi tolerant)
    cats = collections.Counter()
    per_copy = collections.Counter()
    agree = utot = 0
    for n in best:
        cov = [(hap, read_cols.get(p, {}).get(n))
               for (p, hap) in psvs if read_cols.get(p, {}).get(n) is not None]
        if not cov:
            cats["unassignable_no_psv"] += 1
            continue
        mism = {nm: sum(1 for hap, ra in cov if hap[nm] != ra) for nm in names}
        mn = min(mism.values())
        tol = max(1, int(TOL_FRAC * len(cov)))
        if mn > tol:
            cats["unexplained"] += 1
            continue
        consistent = [nm for nm in names if mism[nm] == mn]
        if len(consistent) == 1:
            cats["one_copy"] += 1
            per_copy[consistent[0]] += 1
            utot += 1
            agree += int(best[n][0] == names.index(consistent[0]))
        elif len({hap_rows[nm] for nm in consistent}) == 1:
            cats["collapsed_group"] += 1
        else:
            cats["ambiguous"] += 1

    n_copies = len(names)
    if len(psvs) == 0:
        cls = "no_psv"
    elif n_groups == 1:
        cls = "collapsed"
    elif n_groups == n_copies:
        cls = "fully_resolvable"
    else:
        cls = "partial"
    return dict(family=fam_id, n_copies=n_copies, chroms=sorted({r[0] for r in regions}),
                cross_chrom=len({r[0] for r in regions}) >= 2,
                psvs=len(psvs), K=n_groups, group_sizes=sorted((len(v) for v in groups.values()),
                reverse=True), cls=cls, reads=len(best),
                multimappers=sum(1 for n in ties if len(ties[n]) >= 2),
                assign=dict(cats), per_copy=dict(per_copy), agree=agree, utot=utot)


def main():
    if not os.path.exists(FAM_TSV):
        sys.exit(f"missing {FAM_TSV} — run family_def_vg_reinforce.py first")
    fams = collections.defaultdict(list)
    with open(FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))

    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))

    results = []
    skipped = collections.Counter()
    fam_ids = sorted(fams)
    for k, fi in enumerate(fam_ids):
        regions = dedup_copies(fams[fi])
        if len(regions) < 2:
            skipped["single_region_after_dedup"] += 1
            continue
        try:
            r = run_family(fi, regions, genome, contig_len)
            if r is None:
                skipped["backbone_too_short"] += 1
                continue
            results.append(r)
        except Exception as ex:
            skipped["error"] += 1
            sys.stderr.write(f"[fam {fi}] {ex}\n{traceback.format_exc()}\n")
        if (k + 1) % 20 == 0:
            print(f"  ... {k+1}/{len(fam_ids)} families ({len(results)} processed)", flush=True)
    genome.close()

    # aggregate
    by_cls = collections.Counter(r["cls"] for r in results)
    resolvable = [r for r in results if r["cls"] in ("fully_resolvable", "partial")]
    frontier = [r for r in results if r["cls"] in ("collapsed", "no_psv")]
    tot_reads = sum(r["reads"] for r in results)
    tot_one = sum(r["assign"].get("one_copy", 0) for r in results)
    tot_collapsed = sum(r["assign"].get("collapsed_group", 0) for r in results)
    tot_unassign = sum(r["assign"].get("unassignable_no_psv", 0) for r in results)
    agree = sum(r["agree"] for r in results)
    utot = sum(r["utot"] for r in results)

    summary = dict(
        validated_families=len(fams), processed=len(results), skipped=dict(skipped),
        by_class=dict(by_cls),
        resolvable_families=len(resolvable), frontier_families=len(frontier),
        pct_resolvable=round(100 * len(resolvable) / max(len(results), 1), 1),
        pct_frontier=round(100 * len(frontier) / max(len(results), 1), 1),
        total_reads=tot_reads, reads_to_one_copy=tot_one,
        reads_to_collapsed_group=tot_collapsed, reads_unassignable=tot_unassign,
        single_copy_validation=f"{agree}/{utot}",
        single_copy_validation_pct=round(100 * agree / max(utot, 1), 1),
        copies_per_family=round(sum(r["n_copies"] for r in results) / max(len(results), 1), 2),
    )
    out = dict(summary=summary, families=results)
    json.dump(out, open(os.path.join(HERE, "psv_graph_genomewide.json"), "w"), indent=2)
    with open(os.path.join(HERE, "psv_graph_genomewide_perfamily.tsv"), "w") as o:
        o.write("family\tn_copies\tK\tpsvs\tclass\treads\treads_one_copy\tvalidation\tchroms\n")
        for r in sorted(results, key=lambda x: (x["cls"], -x["n_copies"])):
            o.write(f"{r['family']}\t{r['n_copies']}\t{r['K']}\t{r['psvs']}\t{r['cls']}\t"
                    f"{r['reads']}\t{r['assign'].get('one_copy',0)}\t{r['agree']}/{r['utot']}\t"
                    f"{','.join(r['chroms'])}\n")

    print("\n=== GENOME-WIDE PSV RESOLUTION (validated multi-copy families) ===")
    print(f"  validated families      : {summary['validated_families']}")
    print(f"  processed (>=2 copies)   : {summary['processed']}   skipped: {dict(skipped)}")
    print(f"  mean copies/family       : {summary['copies_per_family']}")
    print(f"  by class                 : {dict(by_cls)}")
    print(f"  RESOLVABLE families      : {len(resolvable)} ({summary['pct_resolvable']}%) "
          f"[fully + partial]")
    print(f"  K-FRONTIER families      : {len(frontier)} ({summary['pct_frontier']}%) "
          f"[collapsed + no-PSV]")
    print(f"  reads: {tot_reads} total | {tot_one} -> ONE copy | {tot_collapsed} -> collapsed "
          f"group | {tot_unassign} unassignable (no PSV covered)")
    print(f"  single-copy assignment agrees with best-mapping copy: "
          f"{agree}/{utot} ({summary['single_copy_validation_pct']}%)")
    print(f"\n[+] wrote psv_graph_genomewide.json + _perfamily.tsv")


if __name__ == "__main__":
    main()
