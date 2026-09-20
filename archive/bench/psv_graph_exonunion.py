#!/usr/bin/env python3
"""psv_graph_exonunion.py — re-run the PSV-aware graph on EXON-UNION copy models
(introns removed) instead of genomic spans, to recover the `align_fail` families
whose copies could not be aligned by genomic span (divergent introns break asm20).

Copy sequence = the cached exon-union model S(v) from family_def_vg_reinforce.py
(`vg_reinforce_copies.fa`: each read's aligned exon blocks, unioned, genome-sequence
concatenated -> intron-free). Backbone = longest copy model; copies aligned asm20;
reads (intron-free mature mRNA query seqs) aligned map-hifi. Everything else (PSV
bubbles, copy paths, threading, classification) is identical to psv_graph_genomewide.py
so the two methods are directly comparable.

Usage: python bench/psv_graph_exonunion.py 1 14 21 34 75 175            # the align_fail set
       python bench/psv_graph_exonunion.py --controls                   # method-agreement check
Output: bench/psv_graph_exonunion.json
"""
import collections
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
TMP = "/home/juanfra/winloci_scratch/psvexon"
MIN_READS_PSV = 3
PAD = 2000
TOL_FRAC = 0.15
OVERLAP_MERGE = 0.5
READ_CAP = 6000
ALIGN_FAIL = [1, 14, 21, 34, 75, 175]
CONTROLS = [18, 0, 4, 62]   # 3 fully_resolvable + 1 genuine_K0 (genomic-span succeeded)


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def load_models():
    models = {}
    name = None
    buf = []
    with open(SEQS) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    models[name] = "".join(buf)
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        models[name] = "".join(buf)
    return models


def column_alleles(bampath, ref):
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    for col in bam.pileup(ref, min_base_quality=0, truncate=False, max_depth=300000):
        for pr in col.pileups:
            if pr.query_position is not None:
                out[col.reference_pos][pr.alignment.query_name] = \
                    pr.alignment.query_sequence[pr.query_position]
    bam.close()
    return out


def dedup_with_ids(loci):
    byc = collections.defaultdict(list)
    for lid, c, s, e, nr in loci:
        byc[c].append((s, e, nr, lid))
    regs = []
    for c, it in byc.items():
        it.sort()
        cur = None
        cids = []
        for s, e, nr, lid in it:
            if cur is not None and s < cur[1] and (min(cur[1], e) - s) > \
                    OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                cur = (cur[0], max(cur[1], e), cur[2] + nr)
                cids.append(lid)
                continue
            if cur is not None:
                regs.append((c, cur[0], cur[1], cur[2], cids))
            cur = (s, e, nr)
            cids = [lid]
        if cur is not None:
            regs.append((c, cur[0], cur[1], cur[2], cids))
    return regs


def run_family(fam_id, regions, models, contig_len):
    os.makedirs(TMP, exist_ok=True)
    # each copy-region -> its longest exon-union model
    copies = []
    for (c, s, e, nr, cids) in regions:
        cand = [(len(models[l]), models[l]) for l in cids if models.get(l)]
        if cand:
            copies.append((c, s, e, nr, max(cand)[1]))
    if len(copies) < 2:
        return None
    copies.sort(key=lambda r: -len(r[4]))
    names = [f"copy{i}" for i in range(len(copies))]
    refseq = copies[0][4]
    if len(refseq) < 200:
        return None
    open(f"{TMP}/ref.fa", "w").write(f">{names[0]}\n{refseq}\n")
    sh(f"samtools faidx {TMP}/ref.fa")
    with open(f"{TMP}/copies.fa", "w") as o:
        for nm, cp in list(zip(names, copies))[1:]:
            o.write(f">{nm}\n{cp[4]}\n")
    sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")

    # reads from genome loci -> best-mapping copy (de) + mature query seqs
    bam = pysam.AlignmentFile(BAM, "rb")
    best = {}
    seqs = {}
    for i, (c, s, e, nr, _) in enumerate(copies):
        for r in bam.fetch(c, max(0, s - PAD), min(contig_len.get(c, e + PAD), e + PAD)):
            if r.is_unmapped or r.is_supplementary or r.reference_end is None:
                continue
            if r.reference_end < s or r.reference_start > e:
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
    # reads are intron-free mature mRNA; backbone is intron-free -> map-hifi (colinear)
    sh(f"minimap2 -ax map-hifi {TMP}/ref.fa {TMP}/reads.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/reads.bam - && samtools index {TMP}/reads.bam")

    copy_cols = column_alleles(f"{TMP}/copies.bam", names[0])
    read_cols = column_alleles(f"{TMP}/reads.bam", names[0])
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

    cats = collections.Counter()
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
            utot += 1
            agree += int(best[n][0] == names.index(consistent[0]))
        elif len({hap_rows[nm] for nm in consistent}) == 1:
            cats["collapsed_group"] += 1
        else:
            cats["ambiguous"] += 1

    n_copies = len(names)
    n_aln = 1 + sum(1 for nm in names[1:] if any(nm in copy_cols.get(p, {}) for p in copy_cols))
    if len(psvs) == 0:
        cls = "no_psv"
    elif n_groups == 1:
        cls = "collapsed"
    elif n_groups == n_copies:
        cls = "fully_resolvable"
    else:
        cls = "partial"
    return dict(family=fam_id, n_copies=n_copies, copies_aligned=n_aln, psvs=len(psvs),
                K=n_groups, cls=cls, reads=len(best), assign=dict(cats),
                agree=agree, utot=utot)


def regions_match(ra, rb):
    if len(ra) != len(rb):
        return False
    a = sorted((r[0], r[1], r[2]) for r in ra)
    b = sorted((r[0], r[1], r[2]) for r in rb)
    for (ca, sa, ea), (cb, sb, eb) in zip(a, b):
        if ca != cb or min(ea, eb) - max(sa, sb) <= 0:
            return False
    return True


def run_all(fams, models, contig_len):
    """Full genome-wide exon-union scan with cross-family dedup (the authoritative pass)."""
    regs = {fi: dedup_with_ids(fams[fi]) for fi in fams}
    multi = {fi: r for fi, r in regs.items() if len(r) >= 2}
    ids = sorted(multi)
    parent = {fi: fi for fi in ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if regions_match(multi[ids[i]], multi[ids[j]]):
                parent[find(ids[i])] = find(ids[j])
    groups = collections.defaultdict(list)
    for fi in ids:
        groups[find(fi)].append(fi)
    reps = [sorted(m)[0] for m in groups.values()]
    redundant = len(ids) - len(reps)

    results = []
    for k, fi in enumerate(sorted(reps)):
        r = run_family(fi, multi[fi], models, contig_len)
        if r is None:
            results.append(dict(family=fi, cls="single_after_models", n_copies=0))
        else:
            results.append(r)
        if (k + 1) % 20 == 0:
            print(f"  ... {k+1}/{len(reps)} unique families", flush=True)

    proc = [r for r in results if r.get("n_copies", 0) >= 2]
    refined = []
    for r in proc:
        cls = r["cls"]
        if cls == "no_psv":
            cls = "genuine_K0" if r["copies_aligned"] == r["n_copies"] else "align_fail"
        refined.append({**r, "cls_refined": cls})
    cnt = collections.Counter(r["cls_refined"] for r in refined)
    resolvable = cnt["fully_resolvable"] + cnt["partial"]
    frontier = cnt["genuine_K0"]
    N = len(refined)
    tot_reads = sum(r["reads"] for r in proc)
    tot_one = sum(r["assign"].get("one_copy", 0) for r in proc)
    tot_collapsed = sum(r["assign"].get("collapsed_group", 0) for r in proc)
    tot_unassign = sum(r["assign"].get("unassignable_no_psv", 0) for r in proc)
    agree = sum(r["agree"] for r in proc)
    utot = sum(r["utot"] for r in proc)
    out = dict(
        method="exon_union", redundant_families_collapsed=redundant, unique_families=N,
        by_class=dict(cnt), resolvable=resolvable, frontier=frontier,
        pct_resolvable=round(100 * resolvable / max(N, 1), 1),
        pct_frontier=round(100 * frontier / max(N, 1), 1),
        no_psv_triage={"genuine_K0": cnt["genuine_K0"], "align_fail": cnt["align_fail"]},
        total_reads=tot_reads, reads_to_one_copy=tot_one,
        reads_to_collapsed_group=tot_collapsed, reads_unassignable=tot_unassign,
        single_copy_validation=f"{agree}/{utot}",
        single_copy_validation_pct=round(100 * agree / max(utot, 1), 1),
        families_refined=[dict(family=r["family"], n_copies=r["n_copies"], K=r["K"],
                              psvs=r["psvs"], reads=r["reads"], cls_refined=r["cls_refined"])
                          for r in refined],
    )
    json.dump(out, open(os.path.join(HERE, "psv_graph_exonunion_all.json"), "w"), indent=2)
    print("\n=== EXON-UNION genome-wide (authoritative) ===")
    print(f"  unique families : {N}  (redundant collapsed {redundant})")
    print(f"  by class        : {dict(cnt)}")
    print(f"  RESOLVABLE : {resolvable} ({out['pct_resolvable']}%)   "
          f"genuine K0 : {frontier} ({out['pct_frontier']}%)   align_fail : {cnt['align_fail']}")
    print(f"  reads {tot_reads}: {tot_one} one-copy, {tot_collapsed} collapsed-group, "
          f"{tot_unassign} no-PSV ; single-copy val {agree}/{utot} ({out['single_copy_validation_pct']}%)")
    print("[+] wrote psv_graph_exonunion_all.json")


def main():
    argv = sys.argv[1:]
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))
    models = load_models()
    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))
    genome.close()

    if argv == ["--all"]:
        run_all(fams, models, contig_len)
        return

    if argv == ["--controls"]:
        targets, mode = CONTROLS, "controls"
    elif argv:
        targets, mode = [int(x) for x in argv], "custom"
    else:
        targets, mode = ALIGN_FAIL, "align_fail"
    results = []
    for fi in targets:
        r = run_family(fi, dedup_with_ids(fams[fi]), models, contig_len)
        results.append(r if r else dict(family=fi, cls="still_failed"))
    json.dump({"mode": mode, "results": results},
              open(os.path.join(HERE, f"psv_graph_exonunion_{mode}.json"), "w"), indent=2)
    print(f"=== EXON-UNION re-run ({mode}) ===")
    for r in results:
        if r.get("cls") == "still_failed":
            print(f"  fam{r['family']:>3}: still failed")
            continue
        print(f"  fam{r['family']:>3}: {r['n_copies']}cp aln={r['copies_aligned']} "
              f"PSV={r['psvs']:>3} K={r['K']} -> {r['cls']:16} "
              f"reads={r['reads']:>5} one_copy={r['assign'].get('one_copy',0):>4} "
              f"val={r['agree']}/{r['utot']}")
    print(f"\n[+] wrote psv_graph_exonunion_{mode}.json")


if __name__ == "__main__":
    main()
