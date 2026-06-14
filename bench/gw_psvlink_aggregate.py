#!/usr/bin/env python3
"""Genome-wide PSV-linkage attribution aggregator.

Compares the PSV-linkage-ON output (pl_$C, gc_pl_$C) against the PSV-linkage-OFF
VG baseline (vg_$C, gc_vg_$C) from the three-way scorecard, per chromosome, to
isolate the PSV-linkage channel's marginal NCBI-corroborated contribution:
  PSV net-new = ref transcripts matched by pl but NOT by vg (=/c strict, =/c/j broad).
Also checks pl >= vg (additive; vg-only-not-pl should be 0). Run after gw_psvlink.sh.
"""
import glob
import os
from collections import defaultdict

OUT = "/tmp/gw"
STRICT = {"=", "c"}
BROAD = {"=", "c", "j"}


def tmap_refset(prefix, chrom, codes):
    paths = glob.glob(os.path.join(OUT, f"{prefix}_{chrom}.*.tmap"))
    if not paths:
        return None
    s = set()
    with open(paths[0]) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            ref_id, cls = c[1], c[2]
            if ref_id not in ("-", "") and cls in codes:
                s.add((chrom, ref_id))
    return s


_GENE = {}
def load_genes(chrom):
    for pat in (f"gc_pl_{chrom}", f"gc_vg_{chrom}"):
        for p in glob.glob(os.path.join(OUT, f"{pat}.*.tmap")):
            with open(p) as fh:
                fh.readline()
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) >= 3 and c[1] not in ("-", ""):
                        _GENE.setdefault((chrom, c[1]), c[0] if c[0] != "-" else "?")


def txcount(path):
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if "\ttranscript\t" in line:
                n += 1
    return n


def done_chroms():
    return sorted(os.path.basename(f)[len("pl_"):-len(".done")]
                  for f in glob.glob(os.path.join(OUT, "pl_*.done")))


def main():
    chroms = done_chroms()
    if not chroms:
        print("no pl_*.done chroms yet — run bench/gw_psvlink.sh first")
        return
    tot = defaultdict(int)
    netnew_genes = defaultdict(int)
    per_chrom = []
    for c in chroms:
        load_genes(c)
        for tier_name, codes in (("strict", STRICT), ("broad", BROAD)):
            pl = tmap_refset("gc_pl", c, codes)
            vg = tmap_refset("gc_vg", c, codes)
            if pl is None or vg is None:
                continue
            netnew = pl - vg
            dropped = vg - pl
            tot[f"pl_{tier_name}"] += len(pl)
            tot[f"vg_{tier_name}"] += len(vg)
            tot[f"netnew_{tier_name}"] += len(netnew)
            tot[f"dropped_{tier_name}"] += len(dropped)
            if tier_name == "strict":
                for (ch, rid) in netnew:
                    netnew_genes[_GENE.get((ch, rid), "?")] += 1
                per_chrom.append((c, len(pl), len(vg), len(netnew), len(dropped),
                                  txcount(os.path.join(OUT, f"pl_{c}.gtf")),
                                  txcount(os.path.join(OUT, f"vg_{c}.gtf"))))
    print("=" * 78)
    print("GENOME-WIDE PSV-LINKAGE ATTRIBUTION (pl = PSV-linkage ON vs vg = OFF; NCBI strict =/c)")
    print("=" * 78)
    print(f"{'chrom':<16}{'pl_match':>9}{'vg_match':>9}{'PSV_netnew':>11}{'dropped':>8}{'pl_tx':>8}{'vg_tx':>8}")
    for c, pl, vg, nn, dr, pltx, vgtx in per_chrom:
        print(f"{c:<16}{pl:>9}{vg:>9}{nn:>11}{dr:>8}{pltx:>8}{vgtx:>8}")
    print("-" * 78)
    print("GENOME-WIDE TOTALS")
    print(f"  STRICT (=/c):  vg-matched={tot['vg_strict']}  pl-matched={tot['pl_strict']}  "
          f"PSV-net-new={tot['netnew_strict']}  dropped(pl<vg)={tot['dropped_strict']}")
    print(f"  BROAD (=/c/j): vg-matched={tot['vg_broad']}  pl-matched={tot['pl_broad']}  "
          f"PSV-net-new={tot['netnew_broad']}  dropped(pl<vg)={tot['dropped_broad']}")
    print(f"  (dropped should be 0 — PSV-linkage is additive over the VG layer)")
    print()
    print(f"PSV-linkage net-new genes (strict =/c), top 30 of {len(netnew_genes)}:")
    for g, n in sorted(netnew_genes.items(), key=lambda kv: (-kv[1], kv[0]))[:30]:
        print(f"  {g:<30}{n}")
    with open(os.path.join(OUT, "psvlink_netnew_genes.tsv"), "w") as fh:
        for g, n in sorted(netnew_genes.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"{g}\t{n}\n")
    print(f"\n(full list: {OUT}/psvlink_netnew_genes.tsv)")


if __name__ == "__main__":
    main()
