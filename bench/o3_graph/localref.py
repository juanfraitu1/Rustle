#!/usr/bin/env python3
"""A homology-closed LOCAL reference for a family set: map a world's reads in seconds, not the whole genome's 79 s and
14.5 GB index load per world (and REMOVED's 146-s, 19-GB masked index), with the reads' alignments in genome coordinates.

    localref.py build FAMILIES.json SET OUT_DIR [--pad 10000] [--index MMI]
        OUT_DIR/closure.tsv   every locus any isoform of any copy of the set aligns to (minimap2 splice:hq, -N 500 -p 0.01:
                              more secondaries than the reads' -N 50 -p 0.1, so every locus a read could reach is in), plus each
                              copy's own span, padded by --pad and merged
        OUT_DIR/local.fa(.fai), local.mmi   those regions as contigs `chrom:start-end` (0-based, half-open)
    localref.py mask OUT_DIR MASK_BED MASKED_DIR   the same contigs with MASK_BED's intervals hard-masked (a REMOVED world)
    localref.py map WORLD LOCAL_DIR [--threads 4]  WORLD/E0.bam: the reads mapped to LOCAL_DIR/local.mmi with the baseline's
                                                   exact minimap2 command (run.MM2), lifted to genome coordinates (header: the
                                                   world genome's chromosomes; SA tags lifted), sorted and indexed

Mapping is per read, so a read lands where whole-genome mapping puts it as long as every locus it could reach is in the
closure; `localref.py` is validated against whole-genome mapping before use (docs/PREREG_* / the run's validation table).
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
INDEX = "/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi"
CLOSURE_OPTS = "-x splice:hq -uf -N 500 -p 0.01 --secondary=yes"
PAD = 10000


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MM2 = _load("run").MM2                                    # the reads' exact minimap2 command


def merge_regions(hits, pad, lengths):
    """[(chrom, start, end)] padded by `pad`, clamped to the chromosome, merged where they overlap or touch."""
    out = []
    for ch, s, e in sorted((ch, max(0, s - pad), min(lengths[ch], e + pad)) for ch, s, e in hits):
        if out and out[-1][0] == ch and s <= out[-1][2]:
            out[-1] = (ch, out[-1][1], max(out[-1][2], e))
        else:
            out.append((ch, s, e))
    return out


def contig_name(ch, s, e):
    return f"{ch}:{s}-{e}"


def parse_contig(name):
    ch, span = name.rsplit(":", 1)
    s, e = span.split("-")
    return ch, int(s), int(e)


def _lift_ref(name, pos):
    """A local contig and 1-based position -> the genome chromosome and position."""
    ch, s, _ = parse_contig(name)
    return ch, int(pos) + s


def lift_line(line):
    """One SAM record from the local reference -> genome coordinates (RNAME/POS, RNEXT/PNEXT, SA:Z)."""
    f = line.rstrip("\n").split("\t")
    local_name = f[2]
    if f[2] != "*":
        f[2], p = _lift_ref(f[2], f[3])
        f[3] = str(p)
    if f[6] == "=" and local_name != "*":                  # the mate on the same contig: shift by the same offset
        _, p = _lift_ref(local_name, f[7])
        f[7] = str(p)
    elif f[6] != "*":
        f[6], p = _lift_ref(f[6], f[7])
        f[7] = str(p)
    for i in range(11, len(f)):
        if f[i].startswith("SA:Z:"):
            parts = []
            for sa in f[i][5:].split(";"):
                if sa:
                    x = sa.split(",")
                    x[0], p = _lift_ref(x[0], x[1])
                    x[1] = str(p)
                    parts.append(",".join(x))
            f[i] = "SA:Z:" + "".join(p + ";" for p in parts)
    return "\t".join(f) + "\n"


def closure_regions(queries, target, pad, lengths, extra=(), opts=CLOSURE_OPTS):
    """Every locus the query sequences ({name: seq}) align to in `target` (a FASTA or a prebuilt .mmi), plus `extra`
    [(chrom, start, end)], padded and merged."""
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/q.fa", "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in queries.items())
        paf = subprocess.run(f"minimap2 {opts} -t 4 {target} {d}/q.fa 2>/dev/null", shell=True, check=True,
                             capture_output=True, text=True, executable="/bin/bash").stdout
    hits = [(f[5], int(f[7]), int(f[8])) for f in (ln.split("\t") for ln in paf.splitlines()) if len(f) > 8]
    return merge_regions(hits + list(extra), pad, lengths)


def write_local(genome_fa, regions, out_dir, mask=()):
    """out_dir/local.fa (+ .fai) and local.mmi: the regions as contigs; `mask` [(chrom, start, end)] hard-masked."""
    import pysam
    os.makedirs(out_dir, exist_ok=True)
    fa = pysam.FastaFile(genome_fa)
    with open(f"{out_dir}/local.fa", "w") as fh:
        for ch, s, e in regions:
            seq = list(fa.fetch(ch, s, e))
            for mc, ms, me in mask:
                if mc == ch and ms < e and me > s:
                    for p in range(max(ms, s), min(me, e)):
                        seq[p - s] = "N"
            fh.write(f">{contig_name(ch, s, e)}\n{''.join(seq)}\n")
    with open(f"{out_dir}/regions.tsv", "w") as fh:
        fh.writelines(f"{ch}\t{s}\t{e}\n" for ch, s, e in regions)
    subprocess.run(f"samtools faidx {out_dir}/local.fa && minimap2 -x splice:hq -d {out_dir}/local.mmi {out_dir}/local.fa "
                   f"2> {out_dir}/local.mmi.log", shell=True, check=True, executable="/bin/bash")


def map_reads(fq, local_dir, genome_fai, out_bam, threads=4):
    """The reads mapped to local_dir/local.mmi with MM2, lifted to genome coordinates, sorted and indexed into out_bam."""
    header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in
                                                     (ln.split("\t")[:2] for ln in open(genome_fai)))
    log = open(os.path.splitext(out_bam)[0] + ".mm2.log", "w")
    mm = subprocess.Popen(f"{MM2} -t {threads} {local_dir}/local.mmi {fq}", shell=True, stdout=subprocess.PIPE,
                          stderr=log, text=True, executable="/bin/bash")
    so = subprocess.Popen(f"samtools sort -@2 -o {out_bam} -", shell=True, stdin=subprocess.PIPE, text=True,
                          executable="/bin/bash")
    so.stdin.write(header)
    for line in mm.stdout:
        if not line.startswith("@"):
            so.stdin.write(lift_line(line))
    so.stdin.close()
    if mm.wait() or so.wait():
        raise RuntimeError(f"local mapping failed for {fq}")
    log.close()
    subprocess.run(["samtools", "index", out_bam], check=True)


def set_queries(fams, genome_fa):
    """{name: sense sequence} of every isoform of every copy of a family set, and each copy's span [(chrom, start, end)]."""
    import pysam
    R = _load("realworld")
    fa = pysam.FastaFile(genome_fa)
    q, spans = {}, []
    for f in fams:
        for c in f["copies"]:
            fetch = lambda s, e, ch=c["chrom"]: fa.fetch(ch, s, e)  # noqa: E731
            for t, ex in c["transcripts"].items():
                q[f"{f['family']}|{c['gene_id']}|{t}"] = R.splice(fetch, [tuple(x) for x in ex], c["strand"], {})
            spans.append((c["chrom"], c["start"], c["end"]))
    return q, spans


def cmd_build(a):
    import json
    R = _load("realworld")
    fams = json.load(open(a.families))["sets"][a.set]
    lengths = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in open(R.GENOME + ".fai")}
    q, spans = set_queries(fams, R.GENOME)
    regions = closure_regions(q, a.index, a.pad, lengths, extra=spans)
    write_local(R.GENOME, regions, a.out_dir)
    with open(f"{a.out_dir}/closure.tsv", "w") as fh:
        fh.write(f"# {len(q)} isoforms of {sum(len(f['copies']) for f in fams)} copies, {len(fams)} families; pad {a.pad}\n")
        fh.writelines(f"{ch}\t{s}\t{e}\n" for ch, s, e in regions)
    print(f"[localref] {a.out_dir}: {len(regions)} regions, {sum(e - s for _, s, e in regions):,} bp", flush=True)
    return 0


def cmd_mask(a):
    R = _load("realworld")
    regions = [tuple(ln.split("\t")[:3]) for ln in open(f"{a.out_dir}/regions.tsv")]
    regions = [(ch, int(s), int(e)) for ch, s, e in regions]
    mask = [(f[0], int(f[1]), int(f[2])) for f in (ln.split("\t") for ln in open(a.mask_bed)) if len(f) >= 3]
    write_local(R.GENOME, regions, a.masked_dir, mask=mask)
    print(f"[localref] {a.masked_dir}: {len(mask)} intervals masked", flush=True)
    return 0


def cmd_map(a):
    map_reads(f"{a.world}/reads.E0.fq", a.local_dir, f"{a.world}/genome.ref.fa.fai", f"{a.world}/E0.bam", a.threads)
    print(f"[localref] {a.world}/E0.bam (local reference {a.local_dir})", flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    p.add_argument("families")
    p.add_argument("set")
    p.add_argument("out_dir")
    p.add_argument("--pad", type=int, default=PAD)
    p.add_argument("--index", default=INDEX)
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("mask")
    p.add_argument("out_dir")
    p.add_argument("mask_bed")
    p.add_argument("masked_dir")
    p.set_defaults(func=cmd_mask)
    p = sub.add_parser("map")
    p.add_argument("world")
    p.add_argument("local_dir")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_map)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
