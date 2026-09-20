#!/usr/bin/env python3
"""vg_family_prototype_retrocopy_mark.py — mark retrocopy / processed-pseudogene members in the VG-native catalog.

Retrocopies are processed pseudogenes: intronless, often annotated as pseudogene, and
highly similar to a spliced parent gene.  At the RNA level the de-novo pipeline here is
spliced-only, so classic intronless retrocopies are largely dropped upstream.  The ones
that survive are typically mapped to an annotated pseudogene or are intronless members
of an otherwise spliced family.

This script marks two classes:

  1. annotated_pseudogene  — gene_biotype is "pseudogene" or "transcribed_pseudogene".
  2. intronless_in_spliced_family — annotated gene has 0 introns but sits in a family
     whose other mapped genes are spliced (n_intron > 0).  These are retrocopy-like
     candidates even if the annotation calls them protein_coding.

Outputs:
  bench/vg_family_prototype_retrocopy_member_flags.tsv
  bench/vg_family_prototype_retrocopy_family_flags.tsv
  bench/vg_family_prototype_retrocopy_mark.json

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_retrocopy_mark.py
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict, Counter

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

from family_level_pr_current import build_ctx

CATALOG_TSV = os.path.join(BENCH, "vg_family_prototype.tsv")
GFF = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
INTRON_INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"
META_TSV = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"

MEMBER_OUT = os.path.join(BENCH, "vg_family_prototype_retrocopy_member_flags.tsv")
FAMILY_OUT = os.path.join(BENCH, "vg_family_prototype_retrocopy_family_flags.tsv")
SUMMARY_OUT = os.path.join(BENCH, "vg_family_prototype_retrocopy_mark.json")

PSEUDO_BIOTYPES = {"pseudogene", "transcribed_pseudogene"}
GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")
BIOTYPE_RE = re.compile(r"(?:^|;)gene_biotype=([^;]+)")


def load_gene_biotypes(gff_path):
    """Return dict gene Name -> biotype from gene/pseudogene features."""
    biotypes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            attrs = f[8]
            m = GENE_RE.search(attrs)
            b = BIOTYPE_RE.search(attrs)
            if m:
                biotypes[m.group(1)] = b.group(1) if b else "unknown"
    return biotypes


def load_meta(meta_path):
    d = {}
    with open(meta_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            lid = f"DN_{row['chrom']}_{row['start']}_{row['n_exon']}"
            d[lid] = dict(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                n_exon=int(row["n_exon"]),
                n_reads=int(row["n_reads"]),
            )
    return d


def main():
    print("[*] building annotation context ...", flush=True)
    ctx, *_ = build_ctx()
    gene_of_dn = ctx["gene_of_dn"]

    print("[*] loading gene biotypes ...", flush=True)
    gene_biotype = load_gene_biotypes(GFF)

    print("[*] loading gene intron index ...", flush=True)
    with open(INTRON_INDEX) as fh:
        intron_index = json.load(fh)

    print("[*] loading de-novo meta ...", flush=True)
    meta = load_meta(META_TSV)

    print("[*] loading prototype catalog ...", flush=True)
    fam_members = defaultdict(list)
    with open(CATALOG_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fam_members[int(row["fam_id"])].append(row["member"])

    # per-family spliced status (using mapped genes)
    family_gene_introns = defaultdict(list)
    for fid, members in fam_members.items():
        for m in members:
            g = gene_of_dn.get(m)
            if g and g in intron_index:
                family_gene_introns[fid].append(intron_index[g]["n_intron"])

    # flag each member
    member_rows = []
    family_flags = defaultdict(set)
    for fid, members in fam_members.items():
        fam_n_introns = family_gene_introns.get(fid, [])
        fam_has_spliced = any(n > 0 for n in fam_n_introns)
        for m in members:
            g = gene_of_dn.get(m)
            biotype = gene_biotype.get(g, "unknown") if g else "unknown"
            n_introns = intron_index.get(g, {}).get("n_intron", -1) if g else -1
            dn_n_exon = meta.get(m, {}).get("n_exon", -1)

            flags = []
            reasons = []
            if biotype in PSEUDO_BIOTYPES:
                flags.append("annotated_pseudogene")
                reasons.append("annotated_pseudogene")
            if n_introns == 0 and biotype == "pseudogene":
                flags.append("intronless_pseudogene")
                reasons.append("intronless_pseudogene")
            if n_introns == 0 and biotype == "protein_coding" and fam_has_spliced:
                flags.append("intronless_in_spliced_family")
                reasons.append("intronless_in_spliced_family")
            if dn_n_exon == 1:
                flags.append("denovo_single_exon")
                reasons.append("denovo_single_exon")

            flag = ";".join(flags) if flags else "clean"
            reason = ";".join(reasons) if reasons else "clean"
            member_rows.append(dict(
                member=m,
                fam_id=fid,
                gene=g or "NA",
                gene_biotype=biotype,
                gene_n_introns=n_introns,
                dn_n_exon=dn_n_exon,
                flag=flag,
                reason=reason,
            ))
            if flags:
                family_flags[fid].update(flags)

    # write member flags
    cols = ["member", "fam_id", "gene", "gene_biotype", "gene_n_introns", "dn_n_exon", "flag", "reason"]
    with open(MEMBER_OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in member_rows:
            w.writerow(row)
    print(f"[write] {MEMBER_OUT}", flush=True)

    # write family flags
    family_rows = []
    for fid in sorted(fam_members):
        members = fam_members[fid]
        flagged_members = sum(1 for r in member_rows if r["fam_id"] == fid and r["flag"] != "clean")
        flags = ";".join(sorted(family_flags.get(fid, []))) or "clean"
        family_rows.append(dict(
            fam_id=fid,
            n_members=len(members),
            n_flagged_members=flagged_members,
            flag=flags,
        ))
    with open(FAMILY_OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fam_id", "n_members", "n_flagged_members", "flag"], delimiter="\t")
        w.writeheader()
        for row in family_rows:
            w.writerow(row)
    print(f"[write] {FAMILY_OUT}", flush=True)

    # summary
    flag_counter = Counter()
    for row in member_rows:
        if row["flag"] != "clean":
            for f in row["flag"].split(";"):
                flag_counter[f] += 1

    summary = dict(
        total_families=len(fam_members),
        families_with_flag=len(family_flags),
        total_members=len(member_rows),
        flagged_members=sum(1 for r in member_rows if r["flag"] != "clean"),
        flag_counts=dict(flag_counter),
    )
    with open(SUMMARY_OUT, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"[write] {SUMMARY_OUT}", flush=True)

    print("\n=== SUMMARY ===")
    print(f"Families with any retrocopy/pseudogene flag: {summary['families_with_flag']} / {summary['total_families']}")
    print(f"Flagged members: {summary['flagged_members']} / {summary['total_members']}")
    print("Flag counts:")
    for f, c in flag_counter.most_common():
        print(f"  {f:30s} {c:5d}")


if __name__ == "__main__":
    main()
