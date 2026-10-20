#!/usr/bin/env python3
"""The 16q22 PDPR-like / NPIPB15 / CLEC18 repeat units in CHM13 against four phased haplotypes (HG002 Q100 v1.1 maternal and
paternal, HG00733 HPRC year-1 maternal and paternal): the data and the SVG for the NPIP Read Pools page (user 2026-10-20: "my
advisor keeps saying that it is impossible and that there must be an error, can we add a visualization that proves that?").

    npip_cnv_fig.py OUT_DIR      OUT_DIR/{data.json, section.html}

CHM13 chr16:80,100,000-80,550,000 holds three tandem units (~114.4 kb apart): AC009053.1 / MSTRG PDPR-like (-), an NPIPB15 copy
(+: NPIPB15, LOC124907808, LOC124907807) and a CLEC18 gene (-). Each haplotype's syntenic region is aligned to it (minimap2
-x asm10 -c, secondaries kept); each CHM13 gene body is aligned to each haplotype region; exact 32-mer matches make the dot plots.
"""
import collections
import json
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
CHM13 = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa"
HG002 = "/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/hg002/hg002v1.1.fa"
HG733 = "/mnt/linuxdisk/tmp/o3_graph/donorsim/HG00733/asm"
REGION = ("chr16", 80100000, 80550000)
HAPS = [  # name, fasta, contig, start, end, reverse-complement to CHM13's orientation
    ("HG002 mother", HG002, "chr16_MATERNAL", 74370000, 74630000, False),
    ("HG002 father", HG002, "chr16_PATERNAL", 77600000, 77975000, False),
    ("HG00733 mother", f"{HG733}/maternal.fa", "HG00733#2#JAHEPP010000066.1", 31960000, 32225000, False),
    ("HG00733 father", f"{HG733}/paternal.fa", "HG00733#1#JAHEPQ010000060.1", 15740000, 16120000, True),
]
GENES = [  # CHM13 gene bodies (pool copy table / CAT names), unit 1-3, class
    ("AC009053.1", 80150189, 80186019, "-", 1, "PDPR-like"), ("NPIPB15", 80195644, 80209850, "+", 1, "NPIP"),
    ("CLEC18C", 80226127, 80239904, "-", 1, "CLEC18"),
    ("MSTRG.2242", 80252177, 80302014, "-", 2, "PDPR-like"), ("LOC124907808", 80310037, 80324225, "+", 2, "NPIP"),
    ("CLEC18B", 80340622, 80354307, "-", 2, "CLEC18"),
    ("MSTRG.2243", 80366579, 80415239, "-", 3, "PDPR-like"), ("LOC124907807", 80424387, 80438556, "+", 3, "NPIP"),
    ("CLEC18B'", 80454955, 80468628, "-", 3, "CLEC18"),
]
SEGMENTS = [(80100000, 80146000, 0), (80146000, 80246040, 1), (80246040, 80360443, 2), (80360443, 80474000, 3),
            (80474000, 80550000, 0)]                         # unit 0 = unique flank
K, BIN, STRIDE, CAP = 32, 2000, 4, 6
COMP = str.maketrans("ACGTN", "TGCAN")


def _sibling(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_cnv", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rc(s):
    return s.translate(COMP)[::-1]


def fetch(path, contig, lo, hi, rev=False):
    s = pysam.FastaFile(path).fetch(contig, lo, hi).upper()
    return rc(s) if rev else s


def paf(target, query, extra=()):
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/t.fa", "w").write(f">t\n{target}\n")
        open(f"{d}/q.fa", "w").write(f">q\n{query}\n")
        out = subprocess.run(["minimap2", "-x", "asm10", "-c", "--secondary=yes", "-N", "10", *extra, f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout
    rows = []
    for ln in out.splitlines():
        x = ln.split("\t")
        tags = {t[:5]: t[5:] for t in x[12:]}
        rows.append(dict(qlen=int(x[1]), qs=int(x[2]), qe=int(x[3]), strand=x[4], ts=int(x[7]), te=int(x[8]),
                         nm=int(x[9]), alen=int(x[10]), mapq=int(x[11]), tp=tags.get("tp:A:", "P"), cg=tags.get("cg:Z:", "")))
    return rows


def ribbons(seq, ref, lo):
    """Each CHM13 segment (the unique flanks and units 1-3) aligned on its own to the haplotype region: every hit covering
    >= half the segment -> (unit, q0, q1, t0, t1) in region coordinates. A haplotype with one unit takes all three units'
    ribbons onto it (one alignment of the whole region chains only one of them)."""
    out = []
    for s, e, u in SEGMENTS:
        for r in paf(seq, ref[s - lo:e - lo]):
            if r["qe"] - r["qs"] >= 0.5 * r["qlen"]:
                out.append(dict(unit=u, q0=s - lo + r["qs"], q1=s - lo + r["qe"], t0=r["ts"], t1=r["te"], strand=r["strand"],
                                mapq=r["mapq"], primary=r["tp"] == "P", identity=round(r["nm"] / r["alen"], 4)))
    return out


def dots(xs, ys):
    """{(xbin, ybin): [forward, reverse]} of exact K-mer matches (x sampled every STRIDE; k-mers seen > CAP times dropped)."""
    idx = collections.defaultdict(list)
    for i in range(len(ys) - K + 1):
        idx[ys[i:i + K]].append(i)
    idx = {k: v for k, v in idx.items() if len(v) <= CAP and "N" not in k}
    xcount = collections.Counter(xs[i:i + K] for i in range(0, len(xs) - K + 1, STRIDE))
    out = collections.defaultdict(lambda: [0, 0])
    for i in range(0, len(xs) - K + 1, STRIDE):
        k = xs[i:i + K]
        if xcount[k] > CAP:
            continue
        for j in idx.get(k, ()):
            out[(i // BIN, j // BIN)][0] += 1
        for j in idx.get(rc(k), ()):
            out[(i // BIN, j // BIN)][1] += 1
    return {f"{a},{b}": v for (a, b), v in out.items()}


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    ch, lo, hi = REGION
    ref = fetch(CHM13, ch, lo, hi)
    data = dict(region=dict(chrom=ch, start=lo, end=hi), genes=[dict(name=g[0], start=g[1] - lo, end=g[2] - lo, strand=g[3],
                unit=g[4], cls=g[5]) for g in GENES], unit_period=[GENES[4][1] - GENES[1][1], GENES[7][1] - GENES[4][1]],
                self_dots=dots(ref, ref), haps=[])
    for name, fa, contig, s, e, rev in HAPS:
        seq = fetch(fa, contig, s, e, rev)
        rows = paf(seq, ref)
        prim = [r for r in rows if r["tp"] == "P"]
        span = (min(r["ts"] for r in prim), max(r["te"] for r in prim))
        hits = []
        for g in GENES:
            body = ref[g[1] - lo:g[2] - lo]
            for r in paf(seq, body):
                if r["qe"] - r["qs"] >= 0.9 * r["qlen"]:
                    hits.append(dict(gene=g[0], unit=g[4], cls=g[5], ts=r["ts"], te=r["te"], strand=r["strand"],
                                     identity=round(r["nm"] / r["alen"], 4), mapq=r["mapq"]))
        copies = []                                          # one donor copy per cluster of overlapping same-class hits
        for h in sorted(hits, key=lambda h: (h["cls"], h["ts"])):
            for c in copies:
                if c["cls"] == h["cls"] and h["ts"] < c["te"] and h["te"] > c["ts"]:
                    c["hits"].append(h)
                    c["ts"], c["te"] = min(c["ts"], h["ts"]), max(c["te"], h["te"])
                    break
            else:
                copies.append(dict(cls=h["cls"], ts=h["ts"], te=h["te"], hits=[h]))
        for c in copies:
            best = max(c["hits"], key=lambda h: h["identity"])
            c.update(best_gene=best["gene"], unit=best["unit"], identity=best["identity"], strand=best["strand"],
                     others={h["gene"]: h["identity"] for h in c["hits"]})
            del c["hits"]
        data["haps"].append(dict(name=name, contig=contig, start=s, end=e, reverse=rev, length=len(seq),
                                 syntenic=span, syntenic_len=span[1] - span[0], ribbons=ribbons(seq, ref, lo),
                                 copies=sorted(copies, key=lambda c: c["ts"]), dots=dots(ref, seq)))
        print(f"[npip_cnv] {name}: syntenic {span[1] - span[0]:,} bp (CHM13 {hi - lo:,}); copies "
              + ", ".join(f"{c['cls']}~{c['best_gene']} {c['identity']:.4f}" for c in sorted(copies, key=lambda c: c["ts"])), flush=True)
    json.dump(data, open(f"{out_dir}/data.json", "w"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
