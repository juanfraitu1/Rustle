#!/usr/bin/env python3
"""Q1 of docs/PREREG_o3_maternal_reference_2026-10-08.md: the fate on the reference haplotype C.REF (env O3_REF, default mat) of the reads of
the loci absent from it.

    O3_REF=mat fate.py fasta    # WR/truth/loci.fa: the C.OTHER sequence of every absent locus
    O3_REF=mat fate.py run      # WR/fate/fate.{tsv,json} (needs WR/truth/loci.ref.paf, see Task 5 step 5)
    fate.py unm                 # descriptive: where the 959 unmapped primaries go on mat / pat (reference-independent) -> W/unm.txt
"""
import collections
import csv
import json
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

FATES = ("UNMAPPED", "PARTIAL", "TIED", "ABSORBED_NEAREST", "ABSORBED_OTHER")


def fractions(fates, n):
    """PARTIAL is summed into unmapped (prereg S5). None when the locus has no reads."""
    if n == 0:
        return None
    return {"absorbed": (fates["ABSORBED_NEAREST"] + fates["ABSORBED_OTHER"]) / n,
            "unmapped": (fates["UNMAPPED"] + fates["PARTIAL"]) / n, "tied": fates["TIED"] / n}


def verdict(fr):
    """The bar of the 09-23 registration, per LARGE locus. 'unmapped or tied' is read as the two lost classes together:
    >= 50% lost -> SUPPORTED, 20-50% -> PARTLY, >= 80% absorbed with <= 10% unmapped and <= 10% tied -> REFUTED, else BETWEEN_BARS."""
    if fr is None:
        return "NO_READS"
    lost = fr["unmapped"] + fr["tied"]
    if lost >= 0.5:
        return "SUPPORTED"
    if fr["absorbed"] >= 0.8 and fr["unmapped"] <= 0.1 and fr["tied"] <= 0.1:
        return "REFUTED"
    if 0.2 <= lost < 0.5:
        return "PARTLY"
    return "BETWEEN_BARS"


def cluster_regions(hits, gap=100000, largest_first=False):
    """[(ref, start, end)] -> [(ref, start, end, n)]: hits on one reference within `gap` bp of the region so far join it"""
    out = []
    for ref, a, b in sorted(hits):
        if out and out[-1][0] == ref and a <= out[-1][2] + gap:
            out[-1][2] = max(out[-1][2], b)
            out[-1][3] += 1
        else:
            out.append([ref, a, b, 1])
    out = [tuple(x) for x in out]
    return sorted(out, key=lambda r: -r[3]) if largest_first else out


def bar_applies(kind, n):
    """the registered bar is judged per LARGE absent locus: not for the sex control, not for descriptive loci (Amendment 2)"""
    return kind not in ("sex", "lrpap1_desc") and n >= 20


def fate_rows(recs, labels, paralog):
    """recs: {read: [Rec]} on mat (accessions); labels: {read: locus|'shared'|'ambiguous'}; paralog: {locus: (acc, s, e)}.
    -> {group: {'n', 'fates': Counter, 'de': [floats of absorbed reads], 'reads': [[read, fate, ref, start, de]]}}"""
    out = {}
    for n, g in labels.items():
        if g == "ambiguous" or n not in recs:
            continue
        r = recs[n]
        f = C.classify_fate(r, paralog.get(g))
        d = out.setdefault(g, {"n": 0, "fates": collections.Counter(), "de": [], "reads": []})
        d["n"] += 1
        d["fates"][f] += 1
        p = r[0] if r else None
        if f.startswith("ABSORBED") and p.de is not None:
            d["de"].append(p.de)
        d["reads"].append([n, f, p.ref if p else None, p.start if p else None, p.de if p else None])
    return out


def cmd_fasta():
    fa = pysam.FastaFile(C.HAP_FA.format(C.OTHER))
    with open(f"{C.WR}/truth/loci.fa", "w") as o:
        for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t"):
            o.write(f">{r['locus']}\n{fa.fetch(r['chrom'], int(r['start']), int(r['end'])).upper()}\n")
    print("loci.fa:", sum(1 for ln in open(f"{C.WR}/truth/loci.fa") if ln[0] == ">"), "sequences")


def cmd_run():
    al = C.alias()
    loci_ = {r["locus"]: r for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t")}
    labels = {r["read"]: r["label"] for r in csv.DictReader(open(f"{C.WR}/truth/labels.tsv"), delimiter="\t")}
    best = C.best_hits(f"{C.WR}/truth/loci.ref.paf", al)
    paralog = {k: tuple(v[1:4]) for k, v in best.items()}
    ident = {k: v[4] for k, v in best.items()}
    recs = C.read_records(f"{C.W}/map/reads.{C.REF}.all.bam", al)
    res = fate_rows(recs, labels, paralog)
    out = {"loci": {}, "shared": None}
    os.makedirs(f"{C.WR}/fate", exist_ok=True)
    with open(f"{C.WR}/fate/fate.tsv", "w") as o:
        o.write("group\tkind\tn\t" + "\t".join(FATES) + "\tde_median\tparalog\tparalog_identity\tverdict\n")
        for g in list(loci_) + ["shared"]:
            d = res.get(g, {"n": 0, "fates": collections.Counter(), "de": [], "reads": []})
            fr = fractions(d["fates"], d["n"])
            kind = loci_[g]["kind"] if g in loci_ else "control"
            large = g in loci_ and bar_applies(kind, d["n"])
            v = verdict(fr) if large else ""
            dm = statistics.median(d["de"]) if d["de"] else None
            row = dict(kind=kind, n=d["n"], fates={f: d["fates"][f] for f in FATES}, fractions=fr, verdict=v, de_median=dm,
                       paralog=list(paralog[g]) if g in paralog else None, paralog_identity=ident.get(g), reads=d["reads"])
            if g in loci_:
                out["loci"][g] = row
            else:
                out["shared"] = row
            o.write("\t".join(str(x) for x in [g, kind, d["n"], *[d["fates"][f] for f in FATES],
                                               "" if dm is None else f"{dm:.4f}", ":".join(str(x) for x in paralog.get(g, ())),
                                               "" if g not in ident else f"{ident[g]:.4f}", v]) + "\n")
    json.dump(out, open(f"{C.WR}/fate/fate.json", "w"))
    print(open(f"{C.WR}/fate/fate.tsv").read())
    print("selection: every read of the 34-family and LRPAP1 sets was mapped on `_pri` first, so UNMAPPED/PARTIAL there means "
          f"'mapped on _pri, lost on {C.REF}'; reads unmapped on _pri are only in R_unm (fate.py unm).")


def cmd_unm():
    al = C.alias()
    lines = []
    for hap in ("mat", "pat"):
        recs = C.read_records(f"{C.W}/map/R_unm.{hap}.bam", al)
        ok = [(n, r[0]) for n, r in recs.items() if r and r[0].primary and r[0].qcov >= C.COV_MIN]
        lines.append(f"R_unm on {hap}: {len(recs)} reads, {len(ok)} with a primary at query coverage >= {C.COV_MIN}")
        reg = cluster_regions([(p.ref, p.start, p.end) for _, p in ok], largest_first=True)[:3]
        for ref, a, b, n in reg:
            mq = sorted(p.mapq for _, p in ok if p.ref == ref and a <= p.start < b)
            lines.append(f"  region {ref}:{a}-{b}: {n} reads, median MAPQ {mq[len(mq) // 2]}")
    open(f"{C.W}/unm.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"fasta": cmd_fasta, "run": cmd_run, "unm": cmd_unm}[sys.argv[1]]()
