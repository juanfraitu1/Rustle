#!/usr/bin/env python3
"""
psv_graph_demo.py — PSV-aware variation graph for real great-ape families, end to end:
collapse the copies onto one backbone, surface the PSV bubbles (the columns of
copy-assignment theory), read off each copy's path, then THREAD the real reads and
assign each to a copy by its PSV alleles — bounded by identifiability.

Runs two families to show BOTH regimes:
  - RABL2 (2 copies): divergent enough -> copies resolve -> reads assigned to a copy.
  - RFPL4A array (5 copies): a founder + 4 near-identical duplicates -> the graph
    exposes only 2 resolvable haplotypes (RFPL4A | {copy2..5 collapsed}) = the K-frontier.

assignment = best PSV-allele match (HiFi-error tolerant); a read resolves to ONE copy,
to a COLLAPSED group of identical copies (identifiability-limited), or is unassignable
(covers no distinguishing PSV). Validates uniquely-assigned reads against best-mapping copy.

Run: /home/juanfra/miniforge3/bin/python bench/psv_graph_demo.py
"""
import collections
import json
import os
import subprocess

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
TMP = "/home/juanfra/winloci_scratch/psvdemo"
MIN_READS_PSV = 3   # a PSV must be read-supported (recurrence, not a one-read error)
PAD = 3000

FAMILIES = {
    "RABL2 (2 copies)": [
        ("RABL2A", "NC_073235.2", 15131652, 15147533),
        ("RABL2B", "NC_086018.1", 48818439, 48832011)],
    "RFPL4A array (5 copies)": [
        ("RFPL4A", "NC_073244.2", 67776054, 67779514),
        ("copy2", "NC_073244.2", 67785173, 67788657),
        ("copy3", "NC_073244.2", 67794315, 67797800),
        ("copy4", "NC_073244.2", 67803458, 67806943),
        ("copy5", "NC_073244.2", 67812601, 67816085)],
}


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def column_alleles(bampath, ref):
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    for col in bam.pileup(ref, min_base_quality=0, truncate=False, max_depth=200000):
        for pr in col.pileups:
            if pr.query_position is not None:
                out[col.reference_pos][pr.alignment.query_name] = pr.alignment.query_sequence[pr.query_position]
    bam.close()
    return out


def run_family(label, copies):
    os.makedirs(TMP, exist_ok=True)
    refname = copies[0][0]
    g = pysam.FastaFile(GENOME)
    refseq = g.fetch(copies[0][1], copies[0][2], copies[0][3]).upper()
    open(f"{TMP}/ref.fa", "w").write(f">{refname}\n{refseq}\n"); sh(f"samtools faidx {TMP}/ref.fa")
    with open(f"{TMP}/copies.fa", "w") as o:
        for nm, c, s, e in copies[1:]:
            o.write(f">{nm}\n{g.fetch(c, s, e).upper()}\n")
    g.close()
    sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null | samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")

    # reads from ALL copy loci -> best-mapping copy (lowest de) + which copies a read ties
    bam = pysam.AlignmentFile(BAM, "rb")
    best = {}; ties = collections.defaultdict(set); seqs = {}
    for i, (nm, c, s, e) in enumerate(copies):
        for r in bam.fetch(c, max(0, s - PAD), e + PAD):
            if r.is_unmapped or r.is_supplementary:
                continue
            if r.reference_end < s or r.reference_start > e:
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
    with open(f"{TMP}/reads.fa", "w") as o:
        for n, sq in seqs.items():
            o.write(f">{n}\n{sq}\n")
    sh(f"minimap2 -ax splice:hq {TMP}/ref.fa {TMP}/reads.fa 2>/dev/null | samtools sort -o {TMP}/reads.bam - && samtools index {TMP}/reads.bam")

    copy_cols = column_alleles(f"{TMP}/copies.bam", refname)
    read_cols = column_alleles(f"{TMP}/reads.bam", refname)
    names = [nm for nm, *_ in copies]
    # PSV columns: all copies defined, >=2 alleles, read-supported
    psvs = []
    for p in range(len(refseq)):
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            if nm in copy_cols.get(p, {}):
                hap[nm] = copy_cols[p][nm]
        if len(hap) == len(names) and len(set(hap.values())) >= 2 and len(read_cols.get(p, {})) >= MIN_READS_PSV:
            psvs.append((p, hap))
    hap_rows = {nm: tuple(hap[nm] for _, hap in psvs) for nm in names}
    # resolvable groups = copies with identical PSV haplotypes
    groups = collections.defaultdict(list)
    for nm in names:
        groups[hap_rows[nm]].append(nm)
    n_groups = len(groups)

    # THREAD reads: best PSV-allele match (HiFi-tolerant)
    cats = collections.Counter(); per_copy = collections.Counter(); agree = utot = 0
    for n in best:
        cov = [(hap, read_cols.get(p, {}).get(n)) for (p, hap) in psvs if read_cols.get(p, {}).get(n) is not None]
        if not cov:
            cats["unassignable (no PSV covered)"] += 1; continue
        mism = {nm: sum(1 for hap, ra in cov if hap[nm] != ra) for nm in names}
        mn = min(mism.values()); tol = max(1, int(0.15 * len(cov)))
        if mn > tol:
            cats["unexplained (novel/error)"] += 1; continue
        consistent = [nm for nm in names if mism[nm] == mn]
        # is the consistent set a single resolvable group (collapsed copies) or genuinely 1 copy?
        if len(consistent) == 1:
            cats["resolved to ONE copy"] += 1; per_copy[consistent[0]] += 1
            utot += 1; agree += int(best[n][0] == names.index(consistent[0]))
        elif len(set(hap_rows[nm] for nm in consistent)) == 1:
            cats["resolved to a COLLAPSED group"] += 1
        else:
            cats["ambiguous (insufficient PSV coverage)"] += 1
    mm = [n for n in ties if len(ties[n]) >= 2]
    mm_res = 0
    for n in mm:
        cov = [(hap, read_cols.get(p, {}).get(n)) for (p, hap) in psvs if read_cols.get(p, {}).get(n) is not None]
        if not cov:
            continue
        mism = {nm: sum(1 for hap, ra in cov if hap[nm] != ra) for nm in names}
        mn = min(mism.values())
        if mn <= max(1, int(0.15 * len(cov))) and sum(1 for nm in names if mism[nm] == mn) == 1:
            mm_res += 1

    print(f"\n=== {label} ===")
    print(f"  reads={len(best)}  cross-copy multimappers={len(mm)}  PSV columns={len(psvs)}")
    print(f"  copies={len(names)}  ->  RESOLVABLE haplotype-groups={n_groups}  "
          f"({'; '.join('{'+','.join(v)+'}' for v in groups.values())})")
    for k, v in cats.most_common():
        print(f"    {k:38}: {v}")
    print(f"  validation: {agree}/{utot} single-copy assignments agree with best-mapping copy "
          f"({100*agree/max(utot,1):.0f}%)")
    print(f"  multimappers uniquely resolved by PSVs: {mm_res}/{len(mm)}")
    return dict(label=label, reads=len(best), multimappers=len(mm), psvs=len(psvs),
                copies=len(names), resolvable_groups=n_groups,
                groups=[v for v in groups.values()], assignment=dict(cats),
                per_copy=dict(per_copy), validation=f"{agree}/{utot}", agree=agree, utot=utot,
                mm_resolved=mm_res, copy_names=names,
                haplotypes={nm: "".join(hap_rows[nm]) for nm in names})


def main():
    res = [run_family(lab, cp) for lab, cp in FAMILIES.items()]
    json.dump(res, open(os.path.join(HERE, "psv_graph_demo.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE, 'psv_graph_demo.json')}")


if __name__ == "__main__":
    main()
