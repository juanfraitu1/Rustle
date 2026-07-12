#!/usr/bin/env python3
"""Genome-wide artifact audit of a gw_family_catalog `.copies.tsv`.

Non-circular structural checks (no truth needed):
  [1] intra-family overlapping copies  -> DuplicateLocus (recip>0.9) / Containment (readthrough enclosing a copy)
  [2] copies SHARED across families    -> a copy of family A overlaps a copy of family B (chimera/readthrough)
  [3] duplicate families               -> identical copy sets under >1 family_id
  [4] giant single-exon copies         -> n_exon==1 and span > SINGLE_EXON_MAX (readthrough signature)
  [5] giant multi-exon spans           -> span > SPAN_MAX (assembly runaway)
Plus an optional annotation cross-check (corroborating): copies landing on 0 annotated genes, or on the SAME
gene as another copy of the family.
"""
import sys, csv, collections

SINGLE_EXON_MAX = 20_000   # a real intronless mRNA is short; a single-exon 20kb+ span is a readthrough
SPAN_MAX = 500_000         # a transcript longer than this is almost certainly an assembly runaway

def load(path):
    fam = collections.defaultdict(list)
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            fam[r["family_id"]].append((r["tid"], r["chrom"], int(r["start"]), int(r["end"]),
                                        int(r["n_exon"]), r["strand"], int(r["n_reads"])))
    return fam

def overlap(a, b):
    _, ca, sa, ea, *_ = a; _, cb, sb, eb, *_ = b
    if ca != cb: return 0, 0.0
    ov = max(0, min(ea, eb) - max(sa, sb))
    if ov <= 0: return 0, 0.0
    return ov, ov / max(ea - sa, eb - sb)

def load_genes(gff):
    genes = collections.defaultdict(list)
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene": continue
            name = ""
            for f in c[8].split(";"):
                if f.startswith("gene="): name = f[5:]
            genes[c[0]].append((int(c[3]), int(c[4]), name))
    return genes

def genes_at(genes, chrom, s, e):
    return sorted({n for gs, ge, n in genes.get(chrom, []) if gs <= e and ge >= s and n})

def main():
    path = sys.argv[1]
    gff = sys.argv[2] if len(sys.argv) > 2 else None
    fam = load(path)
    ncop = sum(len(v) for v in fam.values())
    print(f"catalog: {len(fam)} families, {ncop} copies\n")

    print("[1] intra-family overlapping copies (copies of one family must be disjoint loci)")
    n1 = 0
    for fid, cs in fam.items():
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                ov, recip = overlap(cs[i], cs[j])
                if ov > 0:
                    n1 += 1
                    kind = "DUPLICATE_LOCUS" if recip > 0.9 else "CONTAINMENT"
                    print(f"    {fid:12s} {kind:15s} recip={recip:.2f} {cs[i][0]} <-> {cs[j][0]} ({ov}bp)")
    print(f"    => {n1} pair(s)\n")

    print("[2] copies SHARED across families (a copy overlaps a copy of a DIFFERENT family = chimera/readthrough)")
    flat = [(fid, c) for fid, cs in fam.items() for c in cs]
    n2 = 0
    for i in range(len(flat)):
        for j in range(i + 1, len(flat)):
            if flat[i][0] == flat[j][0]: continue
            ov, recip = overlap(flat[i][1], flat[j][1])
            if ov > 5000:  # only material overlaps
                n2 += 1
                print(f"    {flat[i][0]}/{flat[i][1][0]} <-> {flat[j][0]}/{flat[j][1][0]} recip={recip:.2f} ({ov}bp)")
    print(f"    => {n2} cross-family shared pair(s)\n")

    print("[3] duplicate families (identical copy interval set under >1 id)")
    sig = collections.defaultdict(list)
    for fid, cs in fam.items():
        sig[tuple(sorted((c[1], c[2], c[3]) for c in cs))].append(fid)
    dups = {k: v for k, v in sig.items() if len(v) > 1}
    for k, v in dups.items():
        print(f"    {v}  ({len(k)} copies)")
    print(f"    => {len(dups)} duplicated set(s)\n")

    print(f"[4] giant SINGLE-EXON copies (n_exon==1, span > {SINGLE_EXON_MAX} = readthrough signature)")
    n4 = 0
    for fid, cs in fam.items():
        for c in cs:
            span = c[3] - c[2]
            if c[4] == 1 and span > SINGLE_EXON_MAX:
                n4 += 1; print(f"    {fid:12s} {c[0]} span={span}bp n_exon=1 n_reads={c[6]}")
    print(f"    => {n4} giant single-exon copy(ies)\n")

    print(f"[5] giant-span copies (span > {SPAN_MAX})")
    n5 = 0
    for fid, cs in fam.items():
        for c in cs:
            span = c[3] - c[2]
            if span > SPAN_MAX:
                n5 += 1; print(f"    {fid:12s} {c[0]} span={span}bp n_exon={c[4]} n_reads={c[6]}")
    print(f"    => {n5} giant-span copy(ies)\n")

    if gff:
        genes = load_genes(gff)
        print("[6] annotation cross-check (corroborating): copies on 0 genes / same gene as a sibling copy")
        n6a = n6b = 0
        for fid, cs in fam.items():
            hit = {}
            for c in cs:
                gs = genes_at(genes, c[1], c[2], c[3])
                if not gs:
                    n6a += 1; print(f"    {fid:12s} {c[0]} overlaps NO annotated gene")
                for g in gs:
                    hit.setdefault(g, []).append(c[0])
            for g, tids in hit.items():
                if len(tids) > 1:
                    n6b += 1; print(f"    {fid:12s} gene {g} claimed by {len(tids)} copies {tids}")
        print(f"    => {n6a} copy(ies) on no gene, {n6b} multi-copy-per-gene collision(s)\n")

    print(f"SUMMARY: intra-overlap={n1} cross-shared={n2} dup-fam={len(dups)} giant-1exon={n4} giant-span={n5}")

if __name__ == "__main__":
    main()
