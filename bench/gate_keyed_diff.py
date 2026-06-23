#!/usr/bin/env python3
"""gate_keyed_diff.py OFF.gtf ON.gtf — coordinate-keyed diff of two rustle GTFs.

GTF line ORDER is rayon-nondeterministic, so a textual diff is meaningless. Key each
transcript by (chrom, strand, exon-chain) and compare the SETS, plus per-transcript
attributes. Reports: transcripts lost (only OFF), gained (only ON), shared, and of the
shared how many have a REAL attribute change (copy_status / rescue_class / copy_id /
reference_id) vs only a cosmetic family_id renumber. Prints one TSV summary line:
    <chrom>\t<tx_off>\t<tx_on>\t<lost>\t<gained>\t<real_attr_changed>\t<cosmetic_only>
"""
import re
import sys

ATTR = re.compile(r'(\w+) "([^"]*)"')
REAL_KEYS = ("copy_status", "rescue_class", "copy_id", "reference_id", "class_code")


def load(path):
    tx = {}            # key -> dict(attrs)
    exons = {}         # tid -> list of (start,end)
    meta = {}          # tid -> (chrom, strand, attrs)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom, feat, start, end, strand, attrs = f[0], f[2], f[3], f[4], f[6], f[8]
        d = dict(ATTR.findall(attrs))
        tid = d.get("transcript_id", "")
        if feat == "transcript":
            meta[tid] = (chrom, strand, d)
        elif feat == "exon":
            exons.setdefault(tid, []).append((int(start), int(end)))
    for tid, (chrom, strand, d) in meta.items():
        ch = tuple(sorted(exons.get(tid, [])))
        key = (chrom, strand, ch)
        tx[key] = d
    return tx


def main():
    off, on = load(sys.argv[1]), load(sys.argv[2])
    ko, kn = set(off), set(on)
    lost = ko - kn
    gained = kn - ko
    shared = ko & kn
    real = cosmetic = 0
    real_examples = []
    for k in shared:
        a, b = off[k], on[k]
        rc = any(a.get(x) != b.get(x) for x in REAL_KEYS)
        fc = a.get("family_id") != b.get("family_id")
        if rc:
            real += 1
            if len(real_examples) < 8:
                diffs = {x: (a.get(x), b.get(x)) for x in REAL_KEYS if a.get(x) != b.get(x)}
                real_examples.append((k[0], k[1], diffs))
        elif fc:
            cosmetic += 1
    chrom = next(iter(off))[0] if off else (next(iter(on))[0] if on else "NA")
    sys.stdout.write(f"{chrom}\t{len(off)}\t{len(on)}\t{len(lost)}\t{len(gained)}\t{real}\t{cosmetic}\n")
    if lost or gained or real:
        sys.stderr.write(f"# {chrom}: lost={len(lost)} gained={len(gained)} real_attr={real}\n")
        for k in list(lost)[:5]:
            sys.stderr.write(f"#   LOST  {k[0]} {k[1]} {len(k[2])}ex {k[2][:2]}...\n")
        for k in list(gained)[:5]:
            sys.stderr.write(f"#   GAIN  {k[0]} {k[1]} {len(k[2])}ex {k[2][:2]}...\n")
        for c, s, d in real_examples:
            sys.stderr.write(f"#   ATTR  {c} {s} {d}\n")


if __name__ == "__main__":
    main()
