#!/usr/bin/env python3
"""Opt-in family HMM layer (docs/superpowers/specs/2026-11-05-family-hmm-design.md): one DNA and one protein
profile per multi-copy family, Pfam labels, and a pre-registered concordance table against the MCL families.
Not a member-discovery or assignment emitter (unmapped-rescue Amendment 5 registered that arm)."""
import collections
import os
import re
import subprocess
import sys


def read_fasta(path):
    """[(header without '>', seq)]"""
    out, hdr, seq = [], None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(seq)))
                hdr, seq = ln[1:].strip(), []
            else:
                seq.append(ln.strip())
    if hdr is not None:
        out.append((hdr, "".join(seq)))
    return out


def copies_by_family_records(records):
    """{family: [(header, seq)]} keyed by the header's first '|' field"""
    fams = collections.defaultdict(list)
    for h, s in records:
        fams[h.split("|")[0]].append((h, s))
    return dict(fams)


def copies_by_family(copies_fa):
    return copies_by_family_records(read_fasta(copies_fa))


def longest_orf(records):
    """esl-translate records of ONE copy -> (header, frame, 0-based nt offset, aa) of the longest, or None"""
    if not records:
        return None
    hdr, aa = max(records, key=lambda r: len(r[1]))
    m = re.search(r"coords=(\d+)\.\.(\d+) .*frame=(\d+)", hdr)
    return hdr, int(m.group(3)), int(m.group(1)) - 1, aa


def esl_translate(copy_fa, min_aa):
    """raw FASTA text of the top-strand ORFs of at least min_aa residues"""
    return subprocess.run(["esl-translate", "-l", str(min_aa), "--watson", copy_fa],
                          capture_output=True, text=True, check=True).stdout
