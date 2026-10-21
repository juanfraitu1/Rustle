#!/usr/bin/env python3
"""PREREG_o2_known_family_smn_2026-10-20 (user: "can we outperform minimap2 in picking which is the real primary or can we
reduce the pool ... to make new loci and thus outperform stringtie and flair?"): each read of a family gets ONE primary --
minimap2's (mm2), O2's corrector's (o2), the donor's truth (truth) -- and each arm's BAM is assembled by StringTie and FLAIR;
each copy's exon-7 inclusion (SMN2 skips it, SMN1 keeps it) is read from the reads and from the assembled transcripts.

    o2loci.py arms BAM COPIES.tsv FAMILY O2.tsv OUT_DIR [--truth T.tsv ...]   OUT_DIR/{mm2,o2[,truth]}.bam + reads.fq
    o2loci.py score OUT_DIR SPEC.json                                         per copy, per arm: reads, StringTie, FLAIR
SPEC.json: {copy: {"chrom", "span": [s, e], "e6": [s, e], "e7": [s, e], "e8": [s, e]}} (0-based, half-open).
"""
import collections
import csv
import importlib.util
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MIN_E7 = 20


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_loci", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _ov(blocks, iv):
    return sum(max(0, min(e, iv[1]) - max(s, iv[0])) for s, e in blocks)


def classify(blocks, spec):
    """"FL" (spans exons 6-8 with >= MIN_E7 bp of exon 7), "D7" (spans them without), or None (does not span)."""
    if not (_ov(blocks, spec["e6"]) and _ov(blocks, spec["e8"])):
        return None
    return "FL" if _ov(blocks, spec["e7"]) >= MIN_E7 else "D7"


def inclusion(n):
    k = n.get("FL", 0) + n.get("D7", 0)
    return n.get("FL", 0) / k if k else None


def choose(recs, copy):
    """The best-AS record of a read on `copy` ([{copy, AS, rec}]), or None."""
    on = [r for r in recs if r["copy"] == copy]
    return max(on, key=lambda r: r["AS"]) if on else None


def as_primary(rec, prim):
    """A primary copy of `rec` carrying the read's sequence and qualities from its primary record `prim`."""
    import pysam
    out = pysam.AlignedSegment.fromstring(rec.to_string(), rec.header)
    seq, qual = prim.query_sequence, prim.query_qualities
    if rec.is_reverse != prim.is_reverse:
        seq, qual = seq.translate(COMP)[::-1], (qual[::-1] if qual is not None else None)
    out.flag = rec.flag & ~0x100 & ~0x800
    out.query_sequence = seq
    if qual is not None:
        out.query_qualities = qual
    return out


def read_gtf(path):
    """{transcript: {chrom, blocks (0-based, half-open, sorted), cov}} of a StringTie or FLAIR GTF."""
    tx = collections.defaultdict(lambda: dict(chrom=None, blocks=[], cov=0.0, reads=0.0))
    for ln in open(path):
        if ln.startswith("#"):
            continue
        x = ln.rstrip("\n").split("\t")
        if len(x) < 9:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', x[8])
        if not tid:
            continue
        t = tx[tid.group(1)]
        t["chrom"] = x[0]
        if x[2] == "exon":
            t["blocks"].append((int(x[3]) - 1, int(x[4])))
        elif x[2] == "transcript":
            cov = re.search(r'(?:^|;)\s*cov "([^"]+)"', x[8])            # by name: StringTie also writes longcov
            t["cov"] = float(cov.group(1)) if cov else 0.0
            reads = re.search(r'(?:^|;)\s*reads "([^"]+)"', x[8])        # Rustle: the transcript's reads
            t["reads"] = float(reads.group(1)) if reads else 0.0
    for t in tx.values():
        t["blocks"].sort()
    return dict(tx)


def read_map(path):
    """{isoform: reads} of FLAIR's isoform.read.map.txt, whose first column is isoform_gene (the GTF's transcript_id, then
    "_" and its gene_id)."""
    out = {}
    for ln in open(path):
        x = ln.rstrip("\n").split("\t")
        if len(x) >= 2:
            out[x[0].rsplit("_", 1)[0]] = len([r for r in x[1].split(",") if r])
    return out


def tx_inclusion(tx, weight, copies):
    """{copy: {"FL": weight, "D7": weight}} of the transcripts on each copy's span."""
    out = collections.defaultdict(collections.Counter)
    for t, v in tx.items():
        for c, sp in copies.items():
            if v["chrom"] == sp["chrom"] and _ov(v["blocks"], sp["span"]):
                k = classify(v["blocks"], sp)
                if k:
                    out[c][k] += weight.get(t, 0.0)
    return {c: dict(n) for c, n in out.items()}


def arm_copy(arm, calls, read, prim_copy):
    """The copy a read's one record goes to in an arm, or None for no record: mm2 = minimap2's primary; o2 = O2's call
    (calls: {read: copy, or "" when O2 abstained -- PREREG_o2_smn_showcase: an abstained read is left out of the family;
    a read without a row was not tied and keeps the primary); o2coin = as o2 but an abstained read keeps the primary."""
    if arm == "mm2":
        return prim_copy
    c = calls.get(read)
    if c is None:
        return prim_copy
    if c == "":
        return None if arm == "o2" else prim_copy
    return c


def cmd_arms(bam_p, copies_p, family, o2_p, out_dir, truth_ps=()):
    import pysam
    P, S = _sibling("o2psv"), _sibling("summary")
    cps = P.load_copies(copies_p)
    fam = {c for v in cps.values() for c, *_ in v if c.rsplit(".", 1)[0] == family}
    spans = [(ch, s, e) for ch, v in cps.items() for c, _, _, (s, e) in v if c in fam]
    calls = {r["read"]: r["corrector"] for r in csv.DictReader(open(o2_p), delimiter="\t") if r["family"] == family}
    truth = P._truth(list(truth_ps), True) if truth_ps else {}
    os.makedirs(out_dir, exist_ok=True)
    recs, seen = collections.defaultdict(list), set()
    with pysam.AlignmentFile(bam_p) as bam:
        header = bam.header
        for ch, s, e in spans:
            for rec in bam.fetch(ch, s, e):
                key = (rec.query_name, rec.reference_name, rec.reference_start, rec.cigarstring, rec.flag)
                if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS") or key in seen:
                    continue
                seen.add(key)
                c = S.attribute(rec.get_blocks(), P._strand(rec), cps.get(rec.reference_name, ()))
                if c in fam:
                    recs[rec.query_name].append(dict(copy=c, AS=rec.get_tag("AS"), rec=rec))
    arms = ["mm2", "o2", "o2coin"] + (["truth"] if truth else [])
    outs = {a: pysam.AlignmentFile(os.path.join(out_dir, f"{a}.unsorted.bam"), "wb", header=header) for a in arms}
    n = collections.Counter()
    with open(os.path.join(out_dir, "reads.fq"), "w") as fq:
        for q, rs in sorted(recs.items()):
            prim = next((r for r in rs if not r["rec"].is_secondary), None)
            if prim is None:
                continue                                        # the read's primary is off the family
            n["reads"] += 1
            seq, qual = prim["rec"].query_sequence, prim["rec"].query_qualities
            if prim["rec"].is_reverse:
                seq, qual = seq.translate(COMP)[::-1], (qual[::-1] if qual is not None else None)
            fq.write(f"@{q}\n{seq}\n+\n{pysam.qualities_to_qualitystring(qual) if qual is not None else 'I' * len(seq)}\n")
            want = {a: arm_copy(a, calls, q, prim["copy"]) for a in ("mm2", "o2", "o2coin")}
            t = truth.get(q)
            if truth:
                want["truth"] = t[1] if t and t[0] in ("unique", "unique_elsewhere") else prim["copy"]
            for a in arms:
                if want[a] is None:
                    n[f"{a}: abstained, left out"] += 1
                    continue
                r = choose(rs, want[a])
                if r is None:
                    n[f"{a}: chosen copy without an alignment"] += 1
                    r = prim
                n[f"{a}: moved"] += r["copy"] != prim["copy"]
                outs[a].write(as_primary(r["rec"], prim["rec"]))
    for a, fh in outs.items():
        fh.close()
        pysam.sort("-o", os.path.join(out_dir, f"{a}.bam"), os.path.join(out_dir, f"{a}.unsorted.bam"))
        pysam.index(os.path.join(out_dir, f"{a}.bam"))
        os.remove(os.path.join(out_dir, f"{a}.unsorted.bam"))
    json.dump(dict(n), open(os.path.join(out_dir, "arms.json"), "w"), indent=1)
    print(f"[o2loci] {family}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())) + f" -> {out_dir}", flush=True)
    return 0


def cmd_score(out_dir, spec_p):
    import pysam
    copies = {c: {k: tuple(v) if isinstance(v, list) else v for k, v in d.items()} for c, d in json.load(open(spec_p)).items()}
    res = {}
    for a in ("mm2", "o2", "o2coin", "truth"):
        bam_p = os.path.join(out_dir, f"{a}.bam")
        if not os.path.exists(bam_p):
            continue
        reads = collections.defaultdict(collections.Counter)
        with pysam.AlignmentFile(bam_p) as bam:
            for rec in bam.fetch(until_eof=True):
                for c, sp in copies.items():
                    if rec.reference_name == sp["chrom"] and _ov(rec.get_blocks(), sp["span"]):
                        k = classify(rec.get_blocks(), sp)
                        if k:
                            reads[c][k] += 1
        res[a] = dict(reads={c: dict(v) for c, v in reads.items()})
        st = os.path.join(out_dir, f"{a}.st.gtf")
        if os.path.exists(st):
            tx = read_gtf(st)
            res[a]["stringtie"] = tx_inclusion(tx, {t: v["cov"] for t, v in tx.items()}, copies)
        fl = os.path.join(out_dir, f"{a}.fl.isoforms.gtf")
        if os.path.exists(fl):
            res[a]["flair"] = tx_inclusion(read_gtf(fl), read_map(os.path.join(out_dir, f"{a}.fl.isoform.read.map.txt")), copies)
        ru = os.path.join(out_dir, f"{a}.rustle.gtf")
        if os.path.exists(ru):
            tx = read_gtf(ru)
            res[a]["rustle"] = tx_inclusion(tx, {t: v["reads"] for t, v in tx.items()}, copies)
    json.dump(res, open(os.path.join(out_dir, "score.json"), "w"), indent=1)
    for kind in ("reads", "stringtie", "flair", "rustle"):
        for c in sorted(copies):
            cells = []
            for a in ("mm2", "o2", "o2coin", "truth"):
                n = res.get(a, {}).get(kind, {}).get(c, {})
                i = inclusion(n)
                cells.append(f"{a} {('%.3f' % i) if i is not None else '-'} (FL {n.get('FL', 0):.0f}, D7 {n.get('D7', 0):.0f})")
            print(f"[o2loci] {kind:9s} {c}: " + "; ".join(cells), flush=True)
    return 0


def main(argv):
    truth = []
    if "--truth" in argv:
        i = argv.index("--truth")
        truth, argv = argv[i + 1:], argv[:i]
    if argv[:1] == ["arms"] and len(argv) == 6:
        return cmd_arms(*argv[1:], truth_ps=truth)
    if argv[:1] == ["score"] and len(argv) == 3:
        return cmd_score(*argv[1:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
