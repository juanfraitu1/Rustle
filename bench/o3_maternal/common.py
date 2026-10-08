#!/usr/bin/env python3
"""Shared helpers of the maternal-reference study (docs/PREREG_o3_maternal_reference_2026-10-08.md). Pure functions are unit-tested in test_common.py."""
import collections
import csv
import os
import re

import pysam

TRUTH = "/mnt/linuxdisk/tmp/rna_allele"
REF = os.environ.get("O3_REF", "mat")      # the reference haplotype of this run: mat (copies only the father has are missing) or pat (the reverse)
assert REF in ("mat", "pat"), REF
OTHER = "pat" if REF == "mat" else "mat"   # the truth haplotype: where the missing copies are present
W = "/mnt/linuxdisk/tmp/o3_mat"            # shared by both runs: reads/, map/, artifact/
WR = f"{W}/{REF}"                          # per-run outputs: truth/, fate/, isoc/, inhouse/, score/
HAP_FA = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{}.fa"
HAP_IDX = "/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.{}.splice.mmi"
DELTA = 0.00958        # merge_test.DELTA: the registered allele cutoff (Amendments 7-10)
COV_MIN = 0.80         # registered hit-coverage rule
TIE = 0.98             # registered tie rule (refabsent_truth.express)
IDX = re.compile(r"chr(\w+?)_(mat|pat)_hsa[^_]*")

Rec = collections.namedtuple("Rec", "primary mapq score qcov ref start end de")


def alias():
    """(hap, chromosome number) -> GenBank accession, from TRUTH/{mat,pat}.len.tsv"""
    out = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            out[(h, num)] = acc
    return out


def accession(name, al):
    """haplotype index name chrN_<hap>_hsaX -> accession; any other name (an unplaced scaffold, an unknown number) passes through"""
    m = IDX.fullmatch(name)
    if m and (m.group(2), m.group(1)) in al:
        return al[(m.group(2), m.group(1))]
    return name


def read_records(path, al=None):
    """{read: [Rec, ...]} from a BAM: supplementary records dropped, primary first; an unmapped read maps to []"""
    al = al or {}
    out = {}
    with pysam.AlignmentFile(path) as bam:
        for rd in bam.fetch(until_eof=True):
            if rd.is_supplementary:
                continue
            recs = out.setdefault(rd.query_name, [])
            if rd.is_unmapped:
                continue
            n = rd.infer_read_length() or 1
            recs.append(Rec(not rd.is_secondary, rd.mapping_quality, rd.get_tag("AS") if rd.has_tag("AS") else 0,
                            rd.query_alignment_length / n, accession(rd.reference_name, al), rd.reference_start,
                            rd.reference_end, rd.get_tag("de") if rd.has_tag("de") else None))
    for recs in out.values():
        recs.sort(key=lambda r: not r.primary)
    return out


def place(recs):
    """the read's single best placement over all its records, or None when there is none or the runner-up scores within TIE of it"""
    rs = sorted(recs, key=lambda r: -r.score)
    if not rs or (len(rs) > 1 and rs[1].score >= TIE * rs[0].score):
        return None
    return rs[0]


def classify_fate(recs, paralog):
    """Fate of one read on the reference (prereg S5). recs: Rec list (primary first); paralog: (acc, start, end) of the
    locus' nearest reference paralog, or None."""
    prim = [r for r in recs if r.primary]
    if not prim:
        return "UNMAPPED"
    p = prim[0]
    if p.qcov < COV_MIN:
        return "PARTIAL"
    if p.mapq == 0 or any((not r.primary) and r.score >= TIE * p.score for r in recs):
        return "TIED"
    if paralog and p.ref == paralog[0] and p.start < paralog[2] and paralog[1] < p.end:
        return "ABSORBED_NEAREST"
    return "ABSORBED_OTHER"


def paf_hits(path, al, idmin=0.90, covmin=0.80):
    """{query: [(acc, start, end, identity, query coverage)]} of the PAF records at identity >= idmin and coverage >= covmin"""
    out = collections.defaultdict(list)
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        ident = int(f[9]) / max(1, int(f[10]))
        cov = (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if ident >= idmin and cov >= covmin:
            out[f[0]].append((accession(f[5], al), int(f[7]), int(f[8]), ident, cov))
    return dict(out)


def best_hits(path, al):
    """{query: (score, acc, start, end, identity)}: best PAF record per query by identity x query coverage (control_test.best_hits)"""
    b = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        ident = int(f[9]) / max(1, int(f[10]))
        s = ident * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, accession(f[5], al), int(f[7]), int(f[8]), ident)
    return b
