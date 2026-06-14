#!/usr/bin/env python3
"""4-way attribution aggregator for the genome-wide scorecard.

Splits the VG-only-vs-StringTie NCBI recoveries into:
  - "assembler" : also recovered by rustle's GUIDED BASELINE (rustle -G st, NO --vg)
                  -> rustle's assembler beats StringTie; the VG layer was not needed.
  - "VG-layer"  : recovered ONLY with --vg --vg-layer2 (absent from the guided baseline)
                  -> the secondary-alignment mechanism's own contribution.

Reuses per-chrom gffcompare .tmap files produced by gw_threeway.sh (gc_st_, gc_vg_)
and gw_attribution.sh (gc_gd_). Run after gw_attribution.sh.
"""
import glob
import os
from collections import defaultdict

OUT = "/tmp/gw"
STRICT = {"=", "c"}
BROAD = {"=", "c", "j"}


def tmap_refset(pattern_prefix, chrom, codes):
    """Return the set of (chrom, ref_id) matched by this query's tmap at the given class codes."""
    paths = glob.glob(os.path.join(OUT, f"{pattern_prefix}_{chrom}.*.tmap"))
    s = set()
    if not paths:
        return None  # signal "missing" (distinct from empty)
    with open(paths[0]) as fh:
        fh.readline()  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            ref_id, cls = cols[1], cols[2]
            if ref_id == "-" or ref_id == "":
                continue
            if cls in codes:
                s.add((chrom, ref_id))
    return s


def ref_gene_of(chrom, ref_id):
    """Best-effort ref gene name for a ref_id, from the gc_vg tmap (col0 ref_gene_id)."""
    return _GENE_CACHE.get((chrom, ref_id), "?")


_GENE_CACHE = {}


def load_gene_cache(chrom):
    for pat in (f"gc_vg_{chrom}", f"gc_st_{chrom}", f"gc_gd_{chrom}"):
        for p in glob.glob(os.path.join(OUT, f"{pat}.*.tmap")):
            with open(p) as fh:
                fh.readline()
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 3:
                        continue
                    if c[1] not in ("-", ""):
                        _GENE_CACHE.setdefault((chrom, c[1]), c[0] if c[0] != "-" else "?")


def done_chroms():
    cs = []
    for f in glob.glob(os.path.join(OUT, "gd_*.done")):
        cs.append(os.path.basename(f)[len("gd_"):-len(".done")])
    return sorted(cs)


def main():
    chroms = done_chroms()
    if not chroms:
        print("no gd_*.done chroms yet — run bench/gw_attribution.sh first")
        return

    tot = defaultdict(int)
    assembler_genes = defaultdict(int)
    vglayer_genes = defaultdict(int)
    per_chrom = []

    for c in chroms:
        load_gene_cache(c)
        st = tmap_refset("gc_st", c, STRICT)
        vg = tmap_refset("gc_vg", c, STRICT)
        gd = tmap_refset("gc_gd", c, STRICT)
        if st is None or vg is None or gd is None:
            per_chrom.append((c, "MISSING tmap", 0, 0, 0))
            continue
        vg_only = vg - st                 # VG recovers, StringTie misses (the net-new)
        assembler = vg_only & gd          # guided baseline also has it -> assembler win
        vglayer = vg_only - gd            # only with --vg --vg-layer2 -> VG-layer win
        gd_only = gd - st                 # guided baseline beats ST (assembler, regardless of VG)
        tot["st"] += len(st)
        tot["vg"] += len(vg)
        tot["gd"] += len(gd)
        tot["vg_only"] += len(vg_only)
        tot["assembler"] += len(assembler)
        tot["vglayer"] += len(vglayer)
        tot["gd_only"] += len(gd_only)
        for (ch, rid) in assembler:
            assembler_genes[ref_gene_of(ch, rid)] += 1
        for (ch, rid) in vglayer:
            vglayer_genes[ref_gene_of(ch, rid)] += 1
        per_chrom.append((c, "", len(vg_only), len(assembler), len(vglayer)))

    print("=" * 72)
    print("4-WAY ATTRIBUTION (strict =/c, NCBI RefSeq) — VG-only-vs-StringTie split")
    print("=" * 72)
    print(f"{'chrom':<16}{'VG_only':>9}{'assembler':>11}{'VG_layer':>10}")
    for c, note, vo, asm, vgl in per_chrom:
        if note:
            print(f"{c:<16}  {note}")
        else:
            print(f"{c:<16}{vo:>9}{asm:>11}{vgl:>10}")
    print("-" * 72)
    print("GENOME-WIDE TOTALS (unique NCBI ref transcripts, strict =/c)")
    print(f"  StringTie-matched              : {tot['st']}")
    print(f"  guided-baseline-matched (no vg): {tot['gd']}")
    print(f"  rustle-VG-matched             : {tot['vg']}")
    print(f"  VG-only vs StringTie (net-new): {tot['vg_only']}")
    print(f"    -> attributable to ASSEMBLER (also in guided baseline): {tot['assembler']}")
    print(f"    -> attributable to VG-LAYER  (only with --vg-layer2)  : {tot['vglayer']}")
    print(f"  guided-baseline-only vs ST (assembler net-new): {tot['gd_only']}")
    print()
    print(f"VG-LAYER net-new genes (top 25 of {len(vglayer_genes)}):")
    for g, n in sorted(vglayer_genes.items(), key=lambda kv: (-kv[1], kv[0]))[:25]:
        print(f"  {g:<28}{n}")
    # write full lists
    with open(os.path.join(OUT, "vglayer_genes.tsv"), "w") as fh:
        for g, n in sorted(vglayer_genes.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"{g}\t{n}\n")
    with open(os.path.join(OUT, "assembler_genes.tsv"), "w") as fh:
        for g, n in sorted(assembler_genes.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"{g}\t{n}\n")
    print(f"\n(full lists: {OUT}/vglayer_genes.tsv, {OUT}/assembler_genes.tsv)")


if __name__ == "__main__":
    main()
