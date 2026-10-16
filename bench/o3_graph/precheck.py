#!/usr/bin/env python3
"""Detectability pre-check (round 7; register 1303 review M6): choose the world seed BEFORE any read exists so that REMOVED has
>= NEED detectable removed copies.

    precheck.py FAMILIES.json SET SEED_FROM SEED_TO [--index MMI]

For each seed, each family's special copy (realworld.special_copy) is the copy REMOVED masks. Its main isoform's sequence
is spliced-aligned to the whole genome; its divergence is the `de` of its best alignment away from its own locus (+- 1 kb). That
is what its reads meet once the copy is masked, without sequencing errors. It is detectable when that divergence > DELTA
(O3's definition). The chosen seed is the smallest with >= NEED detectable families.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DELTA, NEED, PAD = 0.00958, 3, 1000
INDEX = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi"


def pick_seed(detectable, need=NEED):
    """{seed: number of detectable families} -> the smallest seed with >= need, or None."""
    ok = [s for s, n in sorted(detectable.items()) if n >= need]
    return ok[0] if ok else None


def away_divergence(alns, own):
    """The `de` of the best-scoring alignment not overlapping `own` (chrom, start, end) +- PAD; 1.0 when there is none."""
    away = [a for a in alns if not (a["chrom"] == own[0] and a["start"] < own[2] + PAD and a["end"] > own[1] - PAD)]
    return max(away, key=lambda a: a["score"])["de"] if away else 1.0


def parse_sam(sam):
    """{query: [{chrom, start (0-based), end, score (AS), de}]} for every mapped record of a SAM text."""
    import re
    alns = {}
    for ln in sam.splitlines():
        if ln.startswith("@"):
            continue
        f = ln.split("\t")
        if int(f[1]) & 4:
            continue
        tags = {x[:5]: x[5:] for x in f[11:]}
        start = int(f[3]) - 1
        span = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", f[5]))
        alns.setdefault(f[0], []).append(dict(chrom=f[2], start=start, end=start + span, score=int(tags.get("AS:i:", 0)),
                                              de=float(tags.get("de:f:", 1.0))))
    return alns


def main(argv):
    import pysam
    spec = importlib.util.spec_from_file_location("o3_graph_realworld", os.path.join(HERE, "realworld.py"))
    R = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(R)
    fams = json.load(open(argv[0]))["sets"][argv[1]]
    seeds = range(int(argv[2]), int(argv[3]))
    index = argv[argv.index("--index") + 1] if "--index" in argv else INDEX
    fa = pysam.FastaFile(R.GENOME)
    special = {s: {f["family"]: R.special_copy(f, s) for f in fams} for s in seeds}
    todo = {(f, c["gene_id"], c["chrom"], c["start"]): c for v in special.values() for f, c in v.items()}
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/q.fa", "w") as fh:
            for (f, g, ch, st), c in todo.items():
                fh.write(f">{f}|{g}|{ch}|{st}\n{R.splice(lambda a, b, ch=ch: fa.fetch(ch, a, b), [tuple(x) for x in c['transcripts'][c['main']]], c['strand'], {})}\n")
        sam = subprocess.run(f"minimap2 -ax splice:hq -uf -N 50 -p 0.1 --secondary=yes -t 4 {index} {d}/q.fa 2>/dev/null",
                             shell=True, capture_output=True, text=True, check=True, executable="/bin/bash").stdout
    alns = parse_sam(sam)
    div = {}
    for (f, g, ch, st), c in todo.items():
        div[(f, g)] = away_divergence(alns.get(f"{f}|{g}|{ch}|{st}", []), (ch, c["start"], c["end"]))
    detectable = {s: sum(1 for f, c in v.items() if div[(f, c["gene_id"])] > DELTA) for s, v in special.items()}
    for s, v in special.items():
        print(s, detectable[s], {f: round(div[(f, c["gene_id"])], 4) for f, c in v.items()})
    chosen = pick_seed(detectable)
    print("chosen seed:", chosen)
    return 0 if chosen is not None else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
