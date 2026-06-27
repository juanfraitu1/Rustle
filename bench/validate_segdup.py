#!/usr/bin/env python3
"""validate_segdup.py — DNA-orthogonal validation of the refined cross-chrom paralog catalog against a
BISER segmental-duplication map.

This is the ONE truly orthogonal check: BISER computes segmental duplications from the GENOME ASSEMBLY
alone, independent of the RNA reads, annotation, gene names, and the exon-sum homology that BUILT the
catalog. A refined family is DNA-CONFIRMED if a segdup pair links >=2 of its copies at the genome level.

Two outputs matter:
 1. DNA-confirmed fraction = an orthogonal precision (not depressed by annotation gaps).
 2. The segdup-SILENT families — exon-sum-homologous (so real RNA-level copies) but NOT DNA-duplicated.
    These are the candidate RETROCOPIES / retrogenes (intronless, RNA-mediated, invisible to DNA segdup
    callers) — i.e. exactly what an RNA-level paralog catalog FINDS that a DNA segdup map MISSES. The
    n_exon signature (retrocopies are exon-poor) is reported to support the interpretation.

Run: python bench/validate_segdup.py <refined_prefix> <biser_segdup_file>
"""
import sys, csv, bisect, collections

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/gw_xchrom_refined"
SEGDUP = sys.argv[2] if len(sys.argv) > 2 else "/home/juanfra/winloci_scratch/GGO_segdup"

# --- 1. refined catalog: family -> [(chrom, start, end, n_exon, strand)] ---
copies = collections.defaultdict(list)
for c in csv.DictReader(open(f"{PREFIX}.copies.tsv"), delimiter="\t"):
    copies[c["family_id"]].append(
        (c["chrom"], int(c["start"]), int(c["end"]), int(c["n_exon"]), c["strand"]))
fams = {f["family_id"]: f for f in csv.DictReader(open(f"{PREFIX}.families.tsv"), delimiter="\t")}

# --- 2. BISER segdups: pairs of genomic intervals (chr1,s1,e1)<->(chr2,s2,e2). Index per chrom with a
#        pair_id+side so two loci can be tested for "linked by the same segdup". ---
anchors = collections.defaultdict(list)  # chrom -> list of (start, end, pair_id, side)
npairs = 0
for line in open(SEGDUP):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 6:
        continue
    try:
        c1, s1, e1, c2, s2, e2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])
    except ValueError:
        continue
    pid = npairs
    anchors[c1].append((s1, e1, pid, 0))
    anchors[c2].append((s2, e2, pid, 1))
    npairs += 1
for c in anchors:
    anchors[c].sort()
starts = {c: [a[0] for a in v] for c, v in anchors.items()}

def hits(chrom, s, e):
    """set of (pair_id, side) for segdup anchors overlapping [s,e) on chrom."""
    out = set()
    iv = anchors.get(chrom)
    if not iv:
        return out
    # anchors with start < e; scan back while end > s (segdups can be long, so scan a generous window).
    hi = bisect.bisect_right(starts[chrom], e)
    j = hi - 1
    # walk left until clearly past (anchor starts far below s by more than any plausible segdup length).
    while j >= 0 and iv[j][0] > s - 5_000_000:
        as_, ae, pid, side = iv[j]
        if ae > s and as_ < e:
            out.add((pid, side))
        j -= 1
    return out

def dna_linked_component(cps):
    """largest set of copies pairwise/transitively linked by a shared segdup pair. Edge(i,j) iff some
    pair_id has copy_i on one side and copy_j on the other."""
    n = len(cps)
    hitsets = [hits(c[0], c[1], c[2]) for c in cps]
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i in range(n):
        pi = {pid: side for pid, side in hitsets[i]}
        for j in range(i + 1, n):
            linked = any(pid in pi and pi[pid] != side for pid, side in hitsets[j])
            if linked:
                parent[find(i)] = find(j)
    comp = collections.Counter(find(i) for i in range(n))
    return max(comp.values(), default=0)

# --- 3. classify each refined family ---
n = conf = conf_xc = n_xc = 0
silent = []
confirmed_examples = []
exons_conf, exons_silent = [], []
for fid, cps in copies.items():
    cross = fams[fid]["cross_chrom"] == "true"
    n += 1; n_xc += cross
    lc = dna_linked_component(cps)
    is_conf = lc >= 2
    conf += is_conf; conf_xc += is_conf and cross
    mean_exon = sum(c[3] for c in cps) / len(cps)
    if is_conf:
        exons_conf.append(mean_exon)
        if len(confirmed_examples) < 8:
            confirmed_examples.append((fid, cross, len(cps), lc, round(mean_exon, 1)))
    else:
        exons_silent.append(mean_exon)
        silent.append((fid, cross, len(cps), round(mean_exon, 1),
                       sorted(set(c[0] for c in cps))))

def pct(a, b):
    return f"{a}/{b} = {100*a/b:.1f}%" if b else "n/a"
def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0

print(f"=== DNA-orthogonal (BISER segdup) validation of {PREFIX} ===")
print(f"segdup pairs loaded: {npairs}\n")
print(f"refined families: {n}  (cross-chrom: {n_xc})")
print(f"DNA-CONFIRMED (>=2 copies linked by a genome segdup): all {pct(conf, n)} | "
      f"same-chrom {pct(conf - conf_xc, n - n_xc)} | cross-chrom {pct(conf_xc, n_xc)}")
print(f"\nDNA-confirmed is ORTHOGONAL to annotation (genome-only). The segdup-SILENT families are the")
print(f"candidate RETROCOPIES — RNA-level copies a DNA segdup caller misses.")
print(f"  mean n_exon: DNA-confirmed={mean(exons_conf):.2f}  segdup-silent={mean(exons_silent):.2f}  "
      f"(retrocopies are exon-poor → expect silent << confirmed)")
print(f"\nDNA-confirmed examples (fid, cross, n_copies, linked, mean_exon):")
for e in confirmed_examples:
    print(f"  {e[0]} xc={e[1]} n={e[2]} linked={e[3]} mean_exon={e[4]}")
print(f"\nsegdup-silent families (candidate retrocopies/sub-threshold), first 15:")
for s in silent[:15]:
    print(f"  {s[0]} xc={s[1]} n_copies={s[2]} mean_exon={s[3]} chroms={s[4][:4]}")

# RABL2 (find the 5-copy family on RABL2's chroms)
print("\n=== RABL2 flagship ===")
RABL2_CHROMS = {"NC_073235.2", "NC_086018.1"}
for fid, cps in copies.items():
    chset = set(c[0] for c in cps)
    if len(cps) == 5 and RABL2_CHROMS <= chset:
        lc = dna_linked_component(cps)
        print(f"  {fid}: 5 copies, DNA-linked component = {lc}/5, mean_exon "
              f"{sum(c[3] for c in cps)/5:.1f}")
        print(f"    -> RABL2A(NC_073235.2)/RABL2B(NC_086018.1) DNA-segdup linkage vs retrocopy silence")
        break
