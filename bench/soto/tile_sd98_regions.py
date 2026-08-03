#!/usr/bin/env python3
"""Tile sd98_regions.fa into overlapping windows before the Soto map-back. THIS IS NOT A TUNING KNOB --
it repairs a defect that silently deleted a third of the duplication signal.

THE BUG. Soto's step 3 command ends in `-p 0.5`, which discards any secondary alignment scoring below half
the primary. That is safe when the query is ONE SD unit. Our sd98_v1.bed had MERGED Soto's ~11k SD units
into 817 non-overlapping blocks (median 12.8 kb but mean 119.7 kb, max 4.25 Mb). minimap2 aligns each whole
block, so the block's own perfect self-alignment becomes the primary with an enormous score (AS 2,747,198
for a 1.37 Mb block); every paralogous chain overlaps it in query space, is therefore marked secondary, and
`-p 0.5` deletes it. Measured best-non-self/self score ratios in offending regions: 0.426, 0.264, 0.205,
0.053 -- all below the 0.5 cut.
Consequence: 106 of 817 regions emitted ONLY their self-alignment, and those regions hold 31.4 Mb = 32.1%
of all SD98 sequence. 69 entire Soto families (AMY1A/B/C+AMY2A, PPIAL4A-F, TRIM64, TRIM49D/53, NAIP,
C4A/C4B/STK19, MBD3L2/4/5, TP53TG3B-F, PGA3/4/5, HERC2P4/5/8, NSUN5P1/2+TRIM73/74, the ID_226 histone
cluster) were STRUCTURALLY INVISIBLE -- not missed, but incapable of being found.

CAUSALLY PROVEN, not inferred: re-mapping four offending regions with `-p 0` took their alignment counts
from 1 to 5,316 / 7,707 / 24,516 / 131,185 and gave 94/94 of their previously edgeless genes both an edge
and a partner in their own Soto family. It is `-p`, not `-N`: on the amylase region (p 0.5, N 200) still
yields 1 alignment while (p 0.0, N 50) yields 7,707. It is not the aligner preset: asm10 and asm20 rescue
only 33 and 37 of the 321 affected genes.

THE FIX. Tiling makes the self-alignment score comparable to a paralog hit, so `-p 0.5` stops eating the
duplications, and Soto's parameters can be used verbatim. Window/step names are kept as chr:start-end so
soto_replicate_clustering.py's parse_region() needs no change.
MEASURED, Soto's exact minimap2 params, MAD 1.0, universe U = their 1,793 Yes genes:
    ARI            0.5140 -> 0.7571
    pair P/R/F1    0.854/0.369/0.516 -> 0.820/0.706/0.758
    exact families 104/491 -> 155/491
    U placed       67.9% -> 91.3%
    edgeless genes 334 -> 8
Tiling is a PROXY for the right input: if Soto's unmerged SD98 unit BED can be obtained, map the units.
"""
import sys

WIN, STEP = 20000, 10000
IN = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/sd98_regions.fa"
OUT = sys.argv[1]


def emit(name, seq, fh):
    c, rng = name.rsplit(":", 1)
    a, b = rng.split("-")
    off = int(a) - 1                     # 0-based genome start of the region
    n = len(seq)
    i = 0
    while True:
        j = min(i + WIN, n)
        if j - i >= 200:
            fh.write(f">{c}:{off+i+1}-{off+j}\n")
            for k in range(i, j, 60):
                fh.write(seq[k:min(k+60, j)] + "\n")
        if j >= n:
            break
        i += STEP


with open(IN) as fin, open(OUT, "w") as fh:
    name, buf = None, []
    for ln in fin:
        if ln[0] == ">":
            if name:
                emit(name, "".join(buf), fh)
            name, buf = ln[1:].strip().split()[0], []
        else:
            buf.append(ln.strip())
    if name:
        emit(name, "".join(buf), fh)
print("done", file=sys.stderr)
