#!/usr/bin/env python3
"""Inputs of the recovery chain with C.REF (env O3_REF, default mat) as the reference (docs/PREREG_o3_maternal_reference_2026-10-08.md section 6).

    O3_REF=mat chain_inputs.py panel     # WR/isoc/{panel.json,labels.tsv,scored.fa,R0.bam}; WR/inhouse/{panel.json,labels.tsv,R0.bam}
    O3_REF=mat chain_inputs.py fasta     # W/<REF>.idx.fa: gorilla_haps/<REF>.fa with the haplotype-index sequence names (chrN_<REF>_hsaX), so that
                                         # BAM, FASTA and .mmi agree for panel_to_copies.py and o3_candidates
"""
import collections
import csv
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import extract_reads as E  # noqa: E402

REF = f"{C.TRUTH}/refabsent"


def merge_loci(hits, gap=5000):
    """hits: [(chrom, start, end)] -> merged [(chrom, start, end)] (hits within `gap` bp on one chromosome join)"""
    out = []
    for c, s, e in sorted(hits):
        if out and out[-1][0] == c and s <= out[-1][2] + gap:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return [tuple(x) for x in out]


def build_panel(copy_loci, fams):
    """copy_loci: {family: [(chrom, start, end)]} (merged or not); control_test.panel layout, mask = first copy, keep = the rest;
    families without a reference locus are dropped"""
    out = []
    for f in fams:
        cps = [(c, s, e, f"{f}:{k}") for k, (c, s, e) in enumerate(merge_loci(copy_loci.get(f, [])))]
        if cps:
            out.append(dict(fam=f, mask=list(cps[0]), keep=[list(x) for x in cps[1:]]))
    return out


def name_map(fai_rows, sq):
    """fai_rows: [(name, length)] of the FASTA; sq: [(name, length)] of the BAM header (index order, which differs from the FASTA order).
    -> {fasta name: index name}, paired by sequence length; a length shared by several sequences pairs them in relative order.
    Refuses unless both sides hold the same multiset of lengths."""
    assert sorted(l for _, l in fai_rows) == sorted(l for _, l in sq), "FASTA and index hold different sequence lengths"
    by_len = {}
    for name, ln in sq:
        by_len.setdefault(ln, []).append(name)
    return {a: by_len[ln].pop(0) for a, ln in fai_rows}


def rename_hits(hits, mapping):
    """[(chrom, s, e)] with FASTA accessions replaced by index names (names already in index form stay)"""
    return [(mapping.get(c, c), s, e) for c, s, e in hits]


def index_names():
    """{FASTA accession: index name} of the reference haplotype: lengths of gorilla_haps/<REF>.fa.fai against the header of its alignments"""
    import pysam
    fai = [(r[0], int(r[1])) for r in csv.reader(open(C.HAP_FA.format(C.REF) + ".fai"), delimiter="\t")]
    with pysam.AlignmentFile(f"{C.W}/map/reads.{C.REF}.all.bam") as b:
        sq = [(s["SN"], s["LN"]) for s in b.header.to_dict()["SQ"]]
    return name_map(fai, sq)


def families():
    return [ln.strip() for ln in open(f"{REF}/fams_bonly.txt") if ln.strip()] + ["LRPAP1"]


def ref_copy_loci():
    """family -> raw hits (index names) on the reference haplotype of the family's copies at identity >= .90, coverage >= .80 (-p/-N of the
    source PAFs: copies.<hap>.paf = Amendment 10; LRPAP1 = asm20 -N 50 -p 0.5)"""
    out = collections.defaultdict(list)
    for q, hs in C.paf_hits(f"{REF}/copies.{C.REF}.paf", {}).items():
        out[q.split(":")[0]].extend((h[0], h[1], h[2]) for h in hs)
    for f in (f"copies8.{C.REF}.paf", f"partial3.{C.REF}.paf"):
        for _q, hs in C.paf_hits(f"{E.LRP}/{f}", {}).items():
            out["LRPAP1"].extend((h[0], h[1], h[2]) for h in hs)
    return out


def cmd_panel():
    fams = families()
    m = index_names()
    panel = build_panel({f: rename_hits(h, m) for f, h in ref_copy_loci().items()}, fams)
    missing = sorted(set(fams) - {p["fam"] for p in panel})
    for d in ("isoc", "inhouse"):
        os.makedirs(f"{C.WR}/{d}", exist_ok=True)
        json.dump(panel, open(f"{C.WR}/{d}/panel.json", "w"), indent=0)
        bam = f"{C.W}/map/reads.{C.REF}.all.bam"
        for ext in ("", ".bai"):
            link = f"{C.WR}/{d}/R0.bam{ext}"
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(bam + ext, link)
    r34 = list(csv.DictReader(open(E.R34), delimiter="\t"))
    seen = {r["read"] for r in r34}
    with open(f"{C.WR}/isoc/labels.tsv", "w") as o:
        o.write("read\tfamily\trole\tcopy\n")
        for r in r34:
            o.write(f"{r['read']}\t{r['family']}\t{r['role']}\t{r['copy']}\n")
        new = 0
        for r in csv.DictReader(open(f"{C.W}/reads/R_LRP.names.tsv"), delimiter="\t"):
            if r["read"] not in seen:
                o.write(f"{r['read']}\tLRPAP1\tS\t{r['cids']}\n")
                new += 1
    with open(f"{C.WR}/isoc/scored.fa", "w") as o:
        for src in (f"{REF}/scored.fa", f"{C.W}/reads/R_LRP.fa"):
            with open(src) as f:
                shutil.copyfileobj(f, o)
    shutil.copy(f"{C.WR}/isoc/labels.tsv", f"{C.WR}/inhouse/labels.tsv")
    print(f"reference {C.REF}: families in the panel {len(panel)} of {len(fams)}; dropped (no {C.REF} locus): {missing}; copies per family: "
          f"{collections.Counter(1 + len(p['keep']) for p in panel)}; labels {len(r34)} + {new} LRPAP1-only reads")


def cmd_fasta():
    m = index_names()
    src, out = C.HAP_FA.format(C.REF), f"{C.W}/{C.REF}.idx.fa"
    with open(src) as f, open(out, "w") as o:
        for ln in f:
            o.write(">" + m[ln[1:].split()[0]] + "\n" if ln[0] == ">" else ln)
    subprocess.run(["samtools", "faidx", out], check=True)
    print(f"{out}: {len(m)} sequences renamed")


if __name__ == "__main__":
    {"panel": cmd_panel, "fasta": cmd_fasta}[sys.argv[1]]()
