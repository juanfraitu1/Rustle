#!/usr/bin/env python3
"""How much of each gene do Soto's shared-exon links share? (descriptive, 2026-10-26)

usage: exon_sharing.py CAT_V4.bed SD98.bed SD98_EXONS.paf S1C.tsv OUTDIR

Re-traces the exon map-back of `soto_replication.py edges --exon-mapback` (same exon set, same PAF, same
`-f 0.99` same-strand hit test) but keeps WHICH exons each gene pair shares. A gene's exons are the merged exonic
blocks of all its CAT v4 transcripts (all blocks, also those outside SD98). For a linked pair (a, b), a's shared
blocks = blocks of a holding an exon that a b-exon maps onto, or an a-exon that maps onto a b-exon (and the same for
b); frac = shared blocks / all blocks. Pair class, multi-exon genes on both sides:
  full        both genes share >= 80% of their blocks (whole-structure copies)
  nested      one gene shares >= 80%, the other < 80% (a partial copy, or exons gained / lost on one side)
  partial     both < 80%, >= 2 shared blocks on some side (a duplicated piece of both genes)
  one_exon    one shared block on each side (the link rests on a single exon)
  single_exon_gene  either gene has one block (no structure to compare)
Then per Soto family (S1C, >= 2 members): which link classes are needed to connect it.
"""
import csv
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import soto_replication as S  # noqa: E402

cat_bed, sd98_bed, paf, s1c, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)

blocks, biotype, name = defaultdict(list), {}, {}
with open(cat_bed) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 21:
            continue
        g, ts = f[18], int(f[1])
        biotype[g] = f[19]
        name[g] = f[3]
        sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
        starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
        blocks[g].extend((ts + o, ts + o + z) for z, o in zip(sizes, starts))
blocks = {g: S.merge_intervals(v) for g, v in blocks.items()}


def block_of(g, s, e):
    for i, (a, b) in enumerate(blocks[g]):
        if a <= s and e <= b:
            return i
    return None


regions = S.load_regions_bed(sd98_bed)
exons, _bt = S.sd98_exons_from_cat(cat_bed, regions)
by_cs = defaultdict(list)
for c, s, e, st, g in exons:
    by_cs[(c, st)].append((s, e, g))
for v in by_cs.values():
    v.sort()
starts = {k: [x[0] for x in v] for k, v in by_cs.items()}
maxlen = max(e - s for c, s, e, st, g in exons)
import bisect

shared = defaultdict(set)  # (gene, partner) -> shared block indices of gene
with open(paf) as fh:
    for ln in fh:
        f = ln.split("\t")
        qname = f[0].split("|")
        gq = qname[0]
        qc, qse = qname[1].split(":")
        qs, qe = (int(x) for x in qse.split("-"))
        c, ts, te, st = f[5], int(f[7]), int(f[8]), f[4]
        v = by_cs.get((c, st))
        if not v:
            continue
        lo = bisect.bisect_left(starts[(c, st)], ts - maxlen)
        hi = bisect.bisect_left(starts[(c, st)], te)
        for s, e, ga in v[lo:hi]:
            ov = min(e, te) - max(s, ts)
            if ov <= 0 or 100 * ov < 99 * (e - s) or ga == gq:
                continue
            bi = block_of(ga, s, e)
            bq = block_of(gq, qs, qe)
            if bi is not None:
                shared[(ga, gq)].add(bi)
            if bq is not None:
                shared[(gq, ga)].add(bq)

edges = S.read_edges(os.path.join(os.path.dirname(os.path.abspath(__file__)), "shared_exons_5154_exon_mapback.tsv"))


def klass(a, b):
    na, nb = len(blocks[a]), len(blocks[b])
    sa, sb = len(shared[(a, b)]), len(shared[(b, a)])
    if na == 1 or nb == 1:
        return "single_exon_gene", sa / na, sb / nb
    fa, fb = sa / na, sb / nb
    if fa >= 0.8 and fb >= 0.8:
        k = "full"
    elif max(fa, fb) >= 0.8:
        k = "nested"
    elif max(sa, sb) >= 2:
        k = "partial"
    else:
        k = "one_exon"
    return k, fa, fb


fam = defaultdict(set)
truth, _amb = S.load_truth(s1c)
for g, f in truth.items():
    fam[f].add(g)
ELIG = {"protein_coding", "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene",
        "translated_unprocessed_pseudogene"}
rows, within = [], Counter()
within_elig = Counter()
for a, b in edges:
    k, fa, fb = klass(a, b)
    same = a in truth and b in truth and truth[a] == truth[b]
    rows.append((a, b, name.get(a, ""), name.get(b, ""), biotype.get(a, ""), biotype.get(b, ""), len(blocks[a]),
                 len(blocks[b]), len(shared[(a, b)]), len(shared[(b, a)]), f"{fa:.3f}", f"{fb:.3f}", k,
                 truth.get(a, "") if same else ""))
    if same:
        within[k] += 1
        if biotype.get(a) in ELIG and biotype.get(b) in ELIG:
            within_elig[k] += 1
with open(os.path.join(out, "pairs.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene_a", "gene_b", "name_a", "name_b", "biotype_a", "biotype_b", "blocks_a", "blocks_b",
                "shared_a", "shared_b", "frac_a", "frac_b", "class", "same_soto_family"])
    w.writerows(rows)

order = ["full", "nested", "partial", "one_exon", "single_exon_gene"]
print("links inside one Soto family, by class (all genes | both coding/unprocessed):")
for k in order:
    print(f"  {k:17s} {within[k]:6d} | {within_elig[k]:6d}")


def connected(members, allowed):
    adj = defaultdict(set)
    for a, b in edges:
        if a in members and b in members and klass(a, b)[0] in allowed:
            adj[a].add(b)
            adj[b].add(a)
    start = next(iter(members))
    seen, st = {start}, [start]
    while st:
        x = st.pop()
        for y in adj[x] - seen:
            seen.add(y)
            st.append(y)
    return len(seen) == len(members)


ladder = [("full", {"full", "single_exon_gene"}), ("+nested", {"full", "single_exon_gene", "nested"}),
          ("+partial", {"full", "single_exon_gene", "nested", "partial"}),
          ("+one_exon", set(order))]
need = Counter()
per = []
for f, mem in sorted(fam.items()):
    if len(mem) < 2 or f.startswith("Unassigned"):
        continue
    lvl = "not connected by edges"
    for lab, allowed in ladder:
        if connected(mem, allowed):
            lvl = lab
            break
    need[lvl] += 1
    per.append((f, len(mem), lvl, ",".join(sorted(name.get(g, g) for g in mem))[:200]))
with open(os.path.join(out, "families.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["family", "members", "weakest_link_class_needed", "member_names"])
    w.writerows(per)
print(f"Soto families (>= 2 members, {sum(need.values())}): weakest link class needed to connect all members")
for lab in ["full", "+nested", "+partial", "+one_exon", "not connected by edges"]:
    print(f"  {lab:24s} {need[lab]}")

# per member: exonic blocks (and bp) shared with ANY other member of its own S1C family
with open(os.path.join(out, "members.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene_id", "name", "biotype", "family", "blocks", "shared_blocks", "exonic_bp", "shared_bp"])
    for f, mem in sorted(fam.items()):
        if f.startswith("Unassigned"):
            continue
        for g in sorted(mem):
            if g not in blocks:
                continue
            sb = set()
            for h in mem:
                if h != g:
                    sb |= shared.get((g, h), set())
            bp = sum(b - a for a, b in blocks[g])
            sbp = sum(blocks[g][i][1] - blocks[g][i][0] for i in sb)
            w.writerow([g, name.get(g, ""), biotype.get(g, ""), f, len(blocks[g]), len(sb), bp, sbp])
