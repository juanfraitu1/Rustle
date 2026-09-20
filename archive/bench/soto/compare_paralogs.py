#!/usr/bin/env python3
"""Everything the pipeline knows about two paralogous loci, side by side.

Usage: compare_paralogs.py <copies.tsv> <bam> <genome.fa> <gff.gz> <GENE_A> <GENE_B>

This is the artifact a reviewer asks for: not a benchmark score, but the evidence for one specific pair.
It answers four questions in order, and each is a different object.

  WHAT IS THERE          read-supported junctions (>= 3 reads) at each locus, and how many are annotated
                         introns of that gene. This is the sashimi view.
  WHAT WE MODELLED       the representative's own chain, and which of the above it captures. The gap is
                         what a browser would show as unaccounted.
  WHY THEY ARE A FAMILY  the E_r edge itself -- identity and coverage of the pairwise alignment, which is
                         the ONLY thing that decides membership under --homology-primary. This is O1.
  WHAT DISTINGUISHES     junctions private to each locus, and whether the shared ones CO-THREAD: each
                         junction of A projected through the alignment onto B. Shared structure argues
                         they are copies; private structure is what could tell a read which copy it came
                         from. That second use is O2 and is NOT part of the membership decision.

Nothing here consults the family assignment, so the comparison is independent of whether the pipeline
happened to group them.
"""
import gzip
import re
import subprocess
import sys
from collections import defaultdict

COP, BAM, GENOME, GFF, GA, GB = sys.argv[1:7]
MIN_J = 3
TOL = 10

genes, tx2gene, exons = {}, {}, defaultdict(set)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "gene":
            m = re.search(r"Name=([^;]+)", f[8])
            if m and m.group(1) in (GA, GB):
                genes.setdefault(m.group(1), (f[0], int(f[3]) - 1, int(f[4]), f[6]))
        elif f[2] in ("mRNA", "transcript"):
            mi, mp = re.search(r"ID=([^;]+)", f[8]), re.search(r"gene=([^;]+)", f[8])
            if mi and mp and mp.group(1) in (GA, GB):
                tx2gene[mi.group(1)] = mp.group(1)
        elif f[2] == "exon":
            mp = re.search(r"Parent=([^;]+)", f[8])
            if mp and mp.group(1) in tx2gene:
                exons[(tx2gene[mp.group(1)], mp.group(1))].add((int(f[3]) - 1, int(f[4])))
annot = defaultdict(set)
for (g, _t), ex in exons.items():
    e = sorted(ex)
    for a, b in zip(e, e[1:]):
        annot[g].add((a[1], b[0]))

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 10:
        copies.append((f[0], f[3], int(f[4]), int(f[5]), f[7],
                       [tuple(int(x) for x in b.split("-")) for b in f[9].split(",")]))


def locus(g):
    c, s, e, strand = genes[g]
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{s}-{e}"],
                         capture_output=True, text=True).stdout
    jr = defaultdict(int)
    n = 0
    for ln in out.splitlines():
        n += 1
        f = ln.split("\t")
        p, k = int(f[3]) - 1, 0
        for ch in f[5]:
            if ch.isdigit():
                k = k * 10 + ord(ch) - 48
            else:
                if ch == "N":
                    jr[(p, p + k)] += 1
                    p += k
                elif ch in "MDX=":
                    p += k
                k = 0
    ov = [(min(e, pe) - max(s, ps), cp) for cp in copies
          for (fam, pc, ps, pe, st, ex) in [cp] if pc == c and s < pe and e > ps]
    best = max(ov)[1] if ov else None
    model = set()
    if best:
        ex = sorted(best[5])
        model = {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
    return {"gene": g, "chrom": c, "s": s, "e": e, "strand": strand, "reads": n,
            "readj": {j for j, v in jr.items() if v >= MIN_J}, "depth": jr,
            "annot": annot.get(g, set()), "model": model, "copy": best}


A, B = locus(GA), locus(GB)
print(f"{'':<24}{GA:>22}{GB:>22}")
print("-" * 68)
for lab, key in (("locus", None), ("annotated span (bp)", "span"), ("primary reads", "reads"),
                 ("junctions >= 3 reads", "nj"), ("  ...annotated", "nja"),
                 ("annotated introns", "na"), ("model junctions", "nm"),
                 ("  ...of read junctions", "nmr"), ("  ...of annotated introns", "nma")):
    def val(L):
        if key is None:
            return f"{L['chrom']}:{L['s']}{L['strand']}"
        if key == "span":
            return f"{L['e']-L['s']:,}"
        if key == "reads":
            return f"{L['reads']:,}"
        if key == "nj":
            return len(L["readj"])
        if key == "nja":
            return len(L["readj"] & L["annot"])
        if key == "na":
            return len(L["annot"])
        if key == "nm":
            return len(L["model"])
        if key == "nmr":
            return f"{len(L['model'] & L['readj'])}/{len(L['readj'])}"
        return f"{len(L['model'] & L['annot'])}/{len(L['annot'])}"
    print(f"{lab:<24}{str(val(A)):>22}{str(val(B)):>22}")

if not (A["copy"] and B["copy"]):
    print("\nAt least one locus has no emitted copy; the membership question cannot be asked.")
    sys.exit(0)
print(f"\nemitted copy            {A['copy'][0]:>22}{B['copy'][0]:>22}")
print(f"  same predicted family: {'YES' if A['copy'][0]==B['copy'][0] else 'NO'}")

# ---- the E_r edge: pairwise alignment of the two emitted copies ------------------------------------
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
fa = f"{CACHE}/pair_{GA}_{GB}.fa"
with open(f"{CACHE}/pair_{GA}_{GB}.txt", "w") as fh:
    for L in (A, B):
        fh.write(f"{L['copy'][1]}:{L['copy'][2]+1}-{L['copy'][3]}\n")
subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/pair_{GA}_{GB}.txt > {fa}", shell=True, check=True)
paf = subprocess.run(["minimap2", "-x", "asm20", "-c", "--eqx", "-X", "-N", "50", "-p", "0.1", fa, fa],
                     capture_output=True, text=True).stdout
print("\nTHE E_r EDGE (genomic spans of the two emitted copies; this alone decides membership)")
if not paf.strip():
    print("  NO ALIGNMENT RECORD -- no edge is possible, whatever the reads say")
else:
    rec = max(paf.splitlines(), key=lambda l: int(l.split("\t")[3]) - int(l.split("\t")[2]))
    f = rec.split("\t")
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    cov = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
    print(f"  identity {ident:.4f}   coverage-of-shorter {cov:.3f}   aligned {int(f[3])-int(f[2]):,} bp")
    print(f"  passes the shipped rule (id >= 0.70, cov >= 0.50): "
          f"{'YES' if ident>=0.70 and cov>=0.50 else 'NO'}")

# ---- shared vs private junction structure ----------------------------------------------------------
print("\nJUNCTION STRUCTURE (>= 3 reads)")
# Junctions MUST be projected through the alignment, not compared by offset from each locus start. The two
# copies have different internal offsets and may be on opposite strands, so a naive offset comparison
# reports 0 shared junctions even for a 99.94%-identical pair. Same trap as the LRRC37 strand bug.
blocks = []
if paf.strip():
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
    q0, t0, qs, ts, rev = int(f[2]), int(f[7]), int(f[2]), int(f[7]), f[4] == "-"
    qlen, tlen = int(f[1]), int(f[6])
    if cg:
        n = 0
        # For a '-' record the CIGAR walks the REVERSE-COMPLEMENTED query, so the aligned block starts at
        # qlen-qend in that frame -- not at qstart, which is given on the original strand.
        qp, tp = ((qlen - int(f[3])) if rev else qs), ts
        for ch in cg:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch in "=XM":
                    blocks.append((qp, qp + n, tp))
                    qp += n
                    tp += n
                elif ch in "DN":
                    tp += n
                elif ch == "I":
                    qp += n
                n = 0
    # which of A, B is the PAF query
    qname_start = int(f[0].rsplit(":", 1)[1].split("-")[0]) - 1
    qL, tL = (A, B) if qname_start == A["copy"][2] else (B, A)

    def project(x):
        # PAF gives query coords on the ORIGINAL strand, but for a '-' record the CIGAR walks the
        # reverse-complemented query. So flip the QUERY axis; target coords are already forward.
        # (Flipping the target instead is what made a 99.94%-identical pair report 0 shared junctions.)
        xq = (qlen - 1 - x) if rev else x
        for a, b, t in blocks:
            if a <= xq < b:
                return t + (xq - a)
        return None

    # Project BOTH endpoints and compare the projected INTERVAL. Under an inverted duplication (PAF
    # strand '-') a junction's donor lands on the other copy's ACCEPTOR, so matching donor-to-donor
    # reports 0 shared junctions for a 99.94%-identical pair. AMY1A/AMY1B is exactly that case.
    tgt_iv = [(k[0] - tL["copy"][2], k[1] - tL["copy"][2]) for k in tL["readj"]]
    shared = 0
    projected = 0
    for j in qL["readj"]:
        p0, p1 = project(j[0] - qL["copy"][2]), project(j[1] - qL["copy"][2])
        if p0 is None or p1 is None:
            continue
        projected += 1
        lo, hi = min(p0, p1), max(p0, p1)
        if any(abs(lo - a) <= TOL and abs(hi - b) <= TOL for (a, b) in tgt_iv):
            shared += 1
    print(f"  {qL['gene']} junctions inside the aligned block: {projected}/{len(qL['readj'])}")
    print(f"  of those, landing within {TOL} bp of a {tL['gene']} junction: {shared} "
          f"({100*shared/max(projected,1):.0f}%)")
    print(f"  private to {qL['gene']}: {projected-shared}   private to {tL['gene']}: "
          f"{len(tL['readj'])-shared}")
else:
    print("  no alignment, so junctions cannot be compared in a common frame")
print(f"  => shared structure argues they are copies; private structure is what could assign a read")
print(f"     to one copy rather than the other (that use is O2, not membership).")
