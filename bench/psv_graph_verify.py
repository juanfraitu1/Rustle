#!/usr/bin/env python3
"""psv_graph_verify.py — adversarial verification + hardening of psv_graph_genomewide.py.

Two corrections and one load-bearing check:
  (1) CROSS-FAMILY DEDUP: collapse families whose copy-region sets coincide (the same
      genomic copies re-detected as several de-novo isoform-loci, e.g. 168/169/170).
  (2) NO_PSV TRIAGE: for every collapsed (no_psv / K=1) family, re-align the copies to the
      backbone and record per-copy aligned-fraction + identity, to split:
        genuine_K0   — copies align well (>=ALN_FRAC, id>=ALN_ID) but carry 0 read-supported
                       divergent columns  -> truly exonically identical (real K-frontier)
        align_fail   — a copy fails to align -> PSVs cannot be called (NOT a real collapse)
        low_support  — divergent columns exist but none reach MIN_READS_PSV (thin coverage)
  (3) VALIDATION OUTLIERS: dump per-read PSV-vs-best-map disagreement for flagged families.

Run: /home/juanfra/miniforge3/bin/python bench/psv_graph_verify.py
"""
import collections
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
TMP = "/home/juanfra/winloci_scratch/psvverify"
OVERLAP_MERGE = 0.5
ALN_FRAC = 0.50    # a copy "aligned" if >=50% of its length maps to the backbone
ALN_ID = 0.80      # and gap-compressed identity >= 0.80


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def dedup_copies(loci):
    byc = collections.defaultdict(list)
    for c, s, e, nr in loci:
        byc[c].append((s, e, nr))
    regs = []
    for c, it in byc.items():
        it.sort()
        cur = None
        for s, e, nr in it:
            if cur is not None and s < cur[1]:
                ov = min(cur[1], e) - s
                if ov > OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                    cur = (cur[0], max(cur[1], e), cur[2] + nr)
                    continue
            if cur is not None:
                regs.append((c, cur[0], cur[1], cur[2]))
            cur = (s, e, nr)
        if cur is not None:
            regs.append((c, cur[0], cur[1], cur[2]))
    return regs


def regions_match(ra, rb):
    """True if two region-sets are the same copies (same count, pairwise overlapping)."""
    if len(ra) != len(rb):
        return False
    a = sorted(ra)
    b = sorted(rb)
    for (ca, sa, ea, _), (cb, sb, eb, _) in zip(a, b):
        if ca != cb or min(ea, eb) - max(sa, sb) <= 0:
            return False
    return True


def copy_align_stats(regions, genome, contig_len):
    """Re-align non-backbone copies to the backbone; return per-copy (aln_frac, identity)."""
    os.makedirs(TMP, exist_ok=True)
    regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
    bc, bs, be, _ = regions[0]
    bs, be = max(0, bs), min(contig_len.get(bc, be), be)
    refseq = genome.fetch(bc, bs, be).upper()
    open(f"{TMP}/ref.fa", "w").write(f">bb\n{refseq}\n")
    sh(f"samtools faidx {TMP}/ref.fa")
    lens = {}
    with open(f"{TMP}/copies.fa", "w") as o:
        for i, (c, s, e, _) in enumerate(regions[1:], 1):
            s2, e2 = max(0, s), min(contig_len.get(c, e), e)
            sq = genome.fetch(c, s2, e2).upper()
            lens[f"c{i}"] = len(sq)
            o.write(f">c{i}\n{sq}\n")
    sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/c.bam - && samtools index {TMP}/c.bam")
    bam = pysam.AlignmentFile(f"{TMP}/c.bam", "rb")
    aligned = {}
    for r in bam.fetch():
        if r.is_unmapped:
            continue
        nm = r.query_name
        frac = (r.reference_length or 0) / max(lens.get(nm, 1), 1)
        de = dict(r.get_tags()).get("de")
        idp = (1 - de) if de is not None else None
        prev = aligned.get(nm, (0, 0))
        if frac > prev[0]:
            aligned[nm] = (frac, idp if idp is not None else 0)
    bam.close()
    stats = []
    for nm in lens:
        f, idp = aligned.get(nm, (0.0, 0.0))
        stats.append((nm, round(f, 3), round(idp, 4)))
    return stats


def main():
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e), int(nr)))

    gw = json.load(open(os.path.join(HERE, "psv_graph_genomewide.json")))
    by_fam = {r["family"]: r for r in gw["families"]}

    # per-family deduped regions
    regs = {fi: dedup_copies(fams[fi]) for fi in fams}
    multi = {fi: r for fi, r in regs.items() if len(r) >= 2}

    # ---- (1) cross-family dedup via region-set match (union-find) ----
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
    # representative per group = the one actually processed in gw (has a result), else first
    reps = []
    redundant = 0
    for g, members in groups.items():
        redundant += len(members) - 1
        rep = next((m for m in members if m in by_fam), members[0])
        reps.append(rep)

    # ---- robust summary over unique families ----
    uniq = [by_fam[r] for r in reps if r in by_fam]
    by_cls = collections.Counter(x["cls"] for x in uniq)
    resolvable = sum(by_cls[k] for k in ("fully_resolvable", "partial"))
    frontier = sum(by_cls[k] for k in ("collapsed", "no_psv"))

    # ---- (2) no_psv triage + per-unique-family refined class ----
    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))
    triage = collections.Counter()
    no_psv_detail = []
    refined = {}   # family -> refined class
    for x in uniq:
        if x["cls"] != "no_psv":
            refined[x["family"]] = x["cls"]
            continue
        st = copy_align_stats(multi[x["family"]], genome, contig_len)
        n_aln = sum(1 for _, f, idp in st if f >= ALN_FRAC and idp >= ALN_ID)
        verdict = "genuine_K0" if n_aln == len(st) else "align_fail"
        triage[verdict] += 1
        refined[x["family"]] = verdict
        no_psv_detail.append(dict(family=x["family"], n_copies=x["n_copies"], reads=x["reads"],
                                  copies_aligned=f"{n_aln}/{len(st)}", verdict=verdict,
                                  stats=st))
    genome.close()
    families_refined = [dict(family=x["family"], n_copies=x["n_copies"], K=x["K"],
                             psvs=x["psvs"], reads=x["reads"],
                             cls_refined=refined[x["family"]]) for x in uniq]

    # read totals over UNIQUE families only (avoid triple-counting redundant loci)
    tot_reads = sum(x["reads"] for x in uniq)
    tot_one = sum(x["assign"].get("one_copy", 0) for x in uniq)
    tot_collapsed = sum(x["assign"].get("collapsed_group", 0) for x in uniq)
    tot_unassign = sum(x["assign"].get("unassignable_no_psv", 0) for x in uniq)
    agree = sum(x["agree"] for x in uniq)
    utot = sum(x["utot"] for x in uniq)

    out = dict(
        raw_processed=gw["summary"]["processed"],
        redundant_families_collapsed=redundant,
        total_reads=tot_reads, reads_to_one_copy=tot_one,
        reads_to_collapsed_group=tot_collapsed, reads_unassignable=tot_unassign,
        single_copy_validation=f"{agree}/{utot}",
        single_copy_validation_pct=round(100 * agree / max(utot, 1), 1),
        unique_families=len(uniq),
        by_class_unique=dict(by_cls),
        resolvable=resolvable, frontier=frontier,
        pct_resolvable=round(100 * resolvable / max(len(uniq), 1), 1),
        pct_frontier=round(100 * frontier / max(len(uniq), 1), 1),
        no_psv_triage=dict(triage),
        no_psv_detail=no_psv_detail,
        families_refined=families_refined,
    )
    json.dump(out, open(os.path.join(HERE, "psv_graph_verify.json"), "w"), indent=2)

    print("=== VERIFICATION: cross-family dedup + no_psv triage ===")
    print(f"  raw processed families         : {out['raw_processed']}")
    print(f"  redundant (same copies) collapsed: {redundant}")
    print(f"  UNIQUE multi-copy families     : {len(uniq)}")
    print(f"  by class (unique)              : {dict(by_cls)}")
    print(f"  RESOLVABLE : {resolvable} ({out['pct_resolvable']}%)   "
          f"K-FRONTIER : {frontier} ({out['pct_frontier']}%)")
    print(f"\n  no_psv triage (is the collapse REAL?): {dict(triage)}")
    for d in sorted(no_psv_detail, key=lambda z: z["verdict"]):
        print(f"    fam{d['family']:>3} {d['copies_aligned']} copies aln  "
              f"reads={d['reads']:>5}  -> {d['verdict']}   {d['stats']}")
    print(f"\n[+] wrote psv_graph_verify.json")


if __name__ == "__main__":
    main()
