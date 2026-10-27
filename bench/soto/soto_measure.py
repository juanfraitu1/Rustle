#!/usr/bin/env python3
"""What Soto's families are made of, how much of each member they share, and what long RNA reads see (descriptive,
2026-10-26).

usage: soto_measure.py MEMBERS.tsv CAT_V4.bed CHAIN_v1.0_to_v1.1 BAM OUTDIR

MEMBERS.tsv = exon_sharing.py's per-member table (S1C families, single-label genes). Blocks = merged exons of all CAT
v4 transcripts of the gene (CHM13 v1.0), lifted to v2.0 (= v1.1 on autosomes) through the T2T chain block by block;
a block that does not lift whole is left out of the RNA counts.
RNA: primary records (-F 0x904) overlapping >= 1 lifted block (any strand). Per member: reads, reads at MAPQ >= 1
(placed with some confidence) and >= 20, and the share of its blocks touched by a read. A member is expressed at
>= 2 reads; distinguishable when >= 2 of its reads have MAPQ >= 20.
"""
import bisect
import csv
import os
import sys
from collections import Counter, defaultdict

members_path, cat_bed, chain_path, bam_path, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
PSEUDO = {"unprocessed_pseudogene", "processed_pseudogene", "transcribed_unprocessed_pseudogene",
          "transcribed_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene",
          "transcribed_unitary_pseudogene", "polymorphic_pseudogene", "pseudogene"}


def group(bt):
    if bt == "protein_coding":
        return "protein_coding"
    if bt in PSEUDO or "pseudogene" in bt:
        return "pseudogene"
    if bt == "lncRNA":
        return "lncRNA"
    return "other"


mem = list(csv.DictReader(open(members_path), delimiter="\t"))
fam = defaultdict(list)
for r in mem:
    fam[r["family"]].append(r)

# --- composition
comp = Counter()
for f, rs in fam.items():
    gs = Counter(group(r["biotype"]) for r in rs)
    comp["families"] += 1
    comp["with >= 1 protein-coding"] += gs["protein_coding"] > 0
    comp["with >= 2 protein-coding"] += gs["protein_coding"] >= 2
    comp["pseudogenes only"] += gs["pseudogene"] == len(rs)
    comp["no protein-coding member"] += gs["protein_coding"] == 0
    comp["with >= 1 lncRNA"] += gs["lncRNA"] > 0
mcomp = Counter(group(r["biotype"]) for r in mem)

# --- completeness: share of each member's exonic bp shared with its family
def frac(r):
    return int(r["shared_bp"]) / int(r["exonic_bp"]) if int(r["exonic_bp"]) else 0.0


bins = [(0.999, "all (100%)"), (0.8, "80-99%"), (0.5, "50-79%"), (0.2, "20-49%"), (0.0, "<20%")]
cbin = defaultdict(Counter)
for r in mem:
    x = frac(r)
    lab = next(l for t, l in bins if x >= t) if x > 0 else "<20%"
    cbin[group(r["biotype"])][lab] += 1
    cbin["all"][lab] += 1
fam_min = Counter()
for f, rs in fam.items():
    m = min(frac(r) for r in rs)
    fam_min["a member shares < 50% of its exonic bp"] += m < 0.5
    fam_min["a member shares < 20%"] += m < 0.2

# --- lift blocks v1.0 -> v1.1 through the chain
chains = defaultdict(list)  # chrom -> [(t_start, t_end, q_chrom, q_start)] ungapped blocks, + strand only
with open(chain_path) as fh:
    cur = None
    for ln in fh:
        f = ln.split()
        if not f:
            continue
        if f[0] == "chain":
            tname, tstart, qname, qstrand, qstart = f[2], int(f[5]), f[7], f[9], int(f[10])
            cur = [tname, tstart, qname, qstart] if qstrand == "+" else None
            continue
        if cur is None:
            continue
        size = int(f[0])
        chains[cur[0]].append((cur[1], cur[1] + size, cur[2], cur[3]))
        if len(f) == 3:
            cur[1] += size + int(f[1])
            cur[3] += size + int(f[2])
for c in chains:
    chains[c].sort()
cstart = {c: [b[0] for b in v] for c, v in chains.items()}


def lift(c, s, e):
    v = chains.get(c)
    if not v:
        return None
    i = bisect.bisect_right(cstart[c], s) - 1
    if i < 0:
        return None
    ts, te, qc, qs = v[i]
    if s >= ts and e <= te:
        return qc, qs + (s - ts), qs + (e - ts)
    return None


blocks = defaultdict(list)
with open(cat_bed) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 21:
            continue
        g, ts = f[18], int(f[1])
        sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
        starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
        blocks[g].extend((f[0], ts + o, ts + o + z) for z, o in zip(sizes, starts))


def merged(v):
    out_ = []
    for c, a, b in sorted(v):
        if out_ and out_[-1][0] == c and a <= out_[-1][2]:
            out_[-1][2] = max(out_[-1][2], b)
        else:
            out_.append([c, a, b])
    return out_


import pysam

bam = pysam.AlignmentFile(bam_path, "rb")
rna = {}
lifted_ok = Counter()
for r in mem:
    g = r["gene_id"]
    bl = [lift(c, a, b) for c, a, b in merged(blocks[g])]
    lifted_ok["blocks"] += len(bl)
    lifted_ok["lifted"] += sum(1 for x in bl if x)
    bl = [x for x in bl if x]
    names, q1, q20, touched = set(), set(), set(), set()
    for i, (c, a, b) in enumerate(bl):
        for rec in bam.fetch(c, a, b):
            if rec.flag & 0x904 or rec.is_unmapped:
                continue
            if not any(min(b, be) > max(a, bs) for bs, be in rec.get_blocks()):
                continue
            names.add(rec.query_name)
            touched.add(i)
            if rec.mapping_quality >= 1:
                q1.add(rec.query_name)
            if rec.mapping_quality >= 20:
                q20.add(rec.query_name)
    rna[g] = (len(names), len(q1), len(q20), len(touched) / len(bl) if bl else 0.0, len(bl))

with open(os.path.join(out, "members_rna.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene_id", "name", "biotype", "family", "shared_frac_bp", "reads", "reads_mapq1", "reads_mapq20",
                "blocks_touched_frac", "blocks_lifted"])
    for r in mem:
        x = rna[r["gene_id"]]
        w.writerow([r["gene_id"], r["name"], r["biotype"], r["family"], f"{frac(r):.3f}", *x[:3], f"{x[3]:.3f}", x[4]])

expr = {g for g, x in rna.items() if x[0] >= 2}
dist = {g for g, x in rna.items() if x[2] >= 2}
fv = Counter()
for f, rs in fam.items():
    n = len(rs)
    ne = sum(r["gene_id"] in expr for r in rs)
    nd = sum(r["gene_id"] in dist for r in rs)
    fv["families"] += 1
    fv["no member expressed"] += ne == 0
    fv["exactly 1 member expressed"] += ne == 1
    fv[">= 2 members expressed (visible as multi-copy)"] += ne >= 2
    fv["all members expressed"] += ne == n
    fv[">= 2 members with >= 2 MAPQ>=20 reads (copies told apart)"] += nd >= 2
mv = Counter()
for r in mem:
    x = rna[r["gene_id"]]
    gr = group(r["biotype"])
    mv[(gr, "members")] += 1
    mv[(gr, "expressed")] += x[0] >= 2
    mv[(gr, "distinguishable")] += x[2] >= 2
reads_all = sum(x[0] for x in rna.values())
reads_q0 = sum(x[0] - x[1] for x in rna.values())

with open(os.path.join(out, "summary.md"), "w") as fo:
    def p(s=""):
        print(s)
        fo.write(s + "\n")
    p(f"## Composition ({comp['families']} S1C families, {len(mem)} single-label members)")
    for k in ["with >= 1 protein-coding", "with >= 2 protein-coding", "no protein-coding member", "pseudogenes only",
              "with >= 1 lncRNA"]:
        p(f"- {k}: {comp[k]} ({100 * comp[k] / comp['families']:.1f}%)")
    p("- members: " + ", ".join(f"{k} {v} ({100 * v / len(mem):.1f}%)" for k, v in mcomp.most_common()))
    p("\n## Completeness: share of each member's exonic bp shared with its own family")
    p("| members | " + " | ".join(l for _, l in bins) + " |")
    p("|---|" + "---|" * len(bins))
    for gr in ["all", "protein_coding", "pseudogene", "lncRNA", "other"]:
        tot = sum(cbin[gr].values())
        if tot:
            p(f"| {gr} ({tot}) | " + " | ".join(f"{cbin[gr][l]} ({100 * cbin[gr][l] / tot:.0f}%)" for _, l in bins) + " |")
    for k, v in fam_min.items():
        p(f"- families where {k}: {v} ({100 * v / comp['families']:.1f}%)")
    p(f"\n## Long RNA reads (A119b testis Iso-Seq; blocks lifted {lifted_ok['lifted']}/{lifted_ok['blocks']})")
    for k in ["no member expressed", "exactly 1 member expressed", ">= 2 members expressed (visible as multi-copy)",
              "all members expressed", ">= 2 members with >= 2 MAPQ>=20 reads (copies told apart)"]:
        p(f"- families with {k}: {fv[k]} ({100 * fv[k] / fv['families']:.1f}%)")
    for gr in ["protein_coding", "pseudogene", "lncRNA", "other"]:
        n = mv[(gr, "members")]
        if n:
            p(f"- {gr} members: {n}; expressed (>= 2 reads) {mv[(gr, 'expressed')]} "
              f"({100 * mv[(gr, 'expressed')] / n:.1f}%); distinguishable (>= 2 MAPQ>=20 reads) "
              f"{mv[(gr, 'distinguishable')]} ({100 * mv[(gr, 'distinguishable')] / n:.1f}%)")
    p(f"- primary reads on family members: {reads_all:,}; MAPQ 0 {reads_q0:,} ({100 * reads_q0 / max(reads_all, 1):.1f}%)")
