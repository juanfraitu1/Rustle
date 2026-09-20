#!/usr/bin/env python3
"""Junction accounting at hand-picked control loci, against BOTH the reads and the RefSeq annotation.

Usage: control_loci.py <copies.tsv> <bam> <gff.gz> [out.tsv]

WHY CONTROLS. The genome-wide number -- a representative carries ~2.8 of a locus's ~24 read-supported
junctions -- is only alarming if those 24 junctions are real. Two controls decide it:

  SINGLE-COPY genes. If an ordinary non-duplicated gene also shows ~24 junctions and ~12% model recall,
  the deficit is in our modelling generally and has nothing to do with segmental duplications. If it shows
  few junctions and high recall, the deficit is specific to duplicated loci, where mis-mapping between
  near-identical copies manufactures junctions that no single copy actually splices.

  THE ANNOTATION. For every locus, how many of the >= 3-read junctions are annotated introns of that gene?
  Junctions that are read-supported AND annotated are real structure the model is missing. Junctions that
  are read-supported but unannotated are candidates for mis-mapping or noise, and missing them is not
  obviously a defect.

Columns: read junctions at >= 3 reads, how many are annotated, annotated introns of the gene, and what
fraction of each the model's own chain recovers. Sorted so the two classes can be compared directly.
"""
import gzip
import re
import subprocess
import sys
from collections import defaultdict

COP, BAM, GFF = sys.argv[1], sys.argv[2], sys.argv[3]
OUT = sys.argv[4] if len(sys.argv) > 4 else None
MIN_J = 3

# Controls: known-good family members, known-bad ones, and ordinary single-copy genes on the same
# chromosomes so that library depth and alignment settings are held constant.
CONTROLS = {
    "AMY1A": "family-good", "AMY2B": "family-good", "GOLGA6L10": "family-good",
    "FAM72A": "family-good", "SEC22B": "family-good",
    "SRGAP2": "family-bad", "SRGAP2C": "family-bad", "NOTCH2NLC": "family-bad",
    "NBPF26": "family-bad",
    # single-copy, well-expressed, chr1/chr15, not in the Soto set
    "RPL22": "single-copy", "SDHB": "single-copy", "CAPZB": "single-copy",
    "PSMB4": "single-copy", "B2M": "single-copy", "IGF1R": "single-copy",
    "THBS1": "single-copy", "CHD2": "single-copy",
}

genes, introns = {}, defaultdict(set)
tx2gene = {}
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in ("chr1", "chr15"):
            continue
        if f[2] == "gene":
            m = re.search(r"Name=([^;]+)", f[8])
            if m and m.group(1) in CONTROLS:
                genes.setdefault(m.group(1), (f[0], int(f[3]) - 1, int(f[4])))
        elif f[2] in ("mRNA", "transcript"):
            mi = re.search(r"ID=([^;]+)", f[8])
            mp = re.search(r"gene=([^;]+)", f[8])
            if mi and mp and mp.group(1) in CONTROLS:
                tx2gene[mi.group(1)] = mp.group(1)
        elif f[2] == "exon":
            mp = re.search(r"Parent=([^;]+)", f[8])
            if mp and mp.group(1) in tx2gene:
                introns[(tx2gene[mp.group(1)], mp.group(1))].add((int(f[3]) - 1, int(f[4])))

# exon sets -> intron sets, per transcript, unioned per gene
annot = defaultdict(set)
for (g, _tx), exs in introns.items():
    e = sorted(exs)
    for a, b in zip(e, e[1:]):
        annot[g].add((a[1], b[0]))

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 10:
        continue
    ex = [tuple(int(x) for x in b.split("-")) for b in f[9].split(",")]
    copies.append((f[3], int(f[4]), int(f[5]), ex))

rows = []
for g, kind in CONTROLS.items():
    if g not in genes:
        continue
    c, s, e = genes[g]
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{s}-{e}"],
                         capture_output=True, text=True).stdout
    jr = defaultdict(int)
    nreads = 0
    for ln in out.splitlines():
        nreads += 1
        f = ln.split("\t")
        p, n = int(f[3]) - 1, 0          # SAM POS is 1-based; everything else is 0-based
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch == "N":
                    jr[(p, p + n)] += 1
                    p += n
                elif ch in "MDX=":
                    p += n
                n = 0
    read_j = {j for j, k in jr.items() if k >= MIN_J}
    ann = annot.get(g, set())
    ov = [(min(e, pe) - max(s, ps), ex) for (pc, ps, pe, ex) in copies if pc == c and s < pe and e > ps]
    model = set()
    if ov:
        ex = sorted(max(ov)[1])
        model = {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
    rows.append((kind, g, nreads, len(read_j), len(read_j & ann), len(ann),
                 len(model & read_j), len(model & ann), len(model)))

hdr = f"{'class':<12}{'gene':<11}{'reads':>7}{'readJ':>7}{'annot':>7}{'annJ':>6}{'mdl':>5}{'mdl∩rd':>8}{'mdl∩ann':>9}"
print(hdr)
print("-" * len(hdr))
for kind in ("single-copy", "family-good", "family-bad"):
    sel = [r for r in rows if r[0] == kind]
    for (k, g, nr, nj, nja, na, mr, ma, nm) in sorted(sel, key=lambda r: -r[3]):
        print(f"{k:<12}{g:<11}{nr:>7}{nj:>7}{nja:>7}{na:>6}{nm:>5}{mr:>8}{ma:>9}")
    if sel:
        aj = sum(r[3] for r in sel) / len(sel)
        ann_frac = sum(r[4] for r in sel) / max(sum(r[3] for r in sel), 1)
        rec = sum(r[7] for r in sel) / max(sum(r[5] for r in sel), 1)
        print(f"{'':<12}{'MEAN/POOL':<11}{'':>7}{aj:>7.1f}"
              f"{'':>7}{'':>6}{'':>5}{'':>8}{'':>9}"
              f"   read-J annotated {100*ann_frac:.0f}%   model recall of ANNOTATED introns {100*rec:.0f}%\n")

if OUT:
    with open(OUT, "w") as fh:
        fh.write("class\tgene\treads\tread_junctions\tread_and_annotated\tannotated_introns"
                 "\tmodel_junctions\tmodel_and_read\tmodel_and_annotated\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in (r[0], r[1], r[2], r[3], r[4], r[5], r[8], r[6], r[7])) + "\n")
    print(f"wrote {OUT}")
