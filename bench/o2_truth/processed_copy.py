#!/usr/bin/env python3
"""Label processed copies from high-score multimappers and score the label (register 1962).

    processed_copy.py BAM RELATED.gtf RELATIONS.tsv OUT_PREFIX [--asm ASSEMBLY.gtf] [--planted COUNTS.tsv]

Rule: a read with one primary P; every secondary S with ms >= 0.9 x P's ms; a pair counts when one of the two alignments
has no intron, the other has one or more, and they do not overlap. A place = overlapping intronless alignments of such
pairs from >= 2 reads; `processed_copy_of` = the gene under the spliced alignments of most of its reads; `own_reads` =
reads whose primary is the intronless alignment. Writes OUT_PREFIX.places.tsv and prints the claims L1-L4 of
docs/PREREG_processed_copy_label_2026-11-25.md. With --planted, reads are named SCENARIO|gene|... (row 1909's set)."""
import sys, re, collections, argparse, pysam
ap = argparse.ArgumentParser(); ap.add_argument('bam'); ap.add_argument('gtf'); ap.add_argument('rel'); ap.add_argument('out')
ap.add_argument('--asm'); ap.add_argument('--planted'); a = ap.parse_args()
BIN = 20000; PP = {'processed_pseudogene', 'transcribed_processed_pseudogene', 'translated_processed_pseudogene'}

class Index:                                                    # intervals with a label, looked up by overlap
    def __init__(self): self.b = collections.defaultdict(list)
    def add(self, c, s, e, v):
        for k in range(s // BIN, e // BIN + 1): self.b[(c, k)].append((s, e, v))
    def hits(self, c, blocks):                                  # label -> overlapping bases
        n = collections.Counter(); seen = set()
        for s, e in blocks:
            for k in range(s // BIN, e // BIN + 1):
                for x, y, v in self.b.get((c, k), ()):
                    if x < e and s < y and (x, y, v, s) not in seen: seen.add((x, y, v, s)); n[v] += min(y, e) - max(x, s)
        return n

at = lambda s, k: (re.search(k + r' "([^"]*)"', s) or [None, None])[1]
G = {}; nex = collections.Counter(); EX = Index()               # gene -> [name, biotype, parent, max exons of a transcript]
for l in open(a.gtf):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t'); g = at(f[8], 'gene_id')
    if f[2] == 'gene': G[g] = [at(f[8], 'gene_name') or g, at(f[8], 'gene_biotype') or '', at(f[8], 'pseudogene_parent'), 0]
    elif f[2] == 'exon': EX.add(f[0], int(f[3]) - 1, int(f[4]), g); nex[(g, at(f[8], 'transcript_id'))] += 1
for (g, t), n in nex.items():
    if g in G: G[g][3] = max(G[g][3], n)
related = set()
for l in open(a.rel):
    f = l.rstrip('\n').split('\t'); related.add((f[1], f[3])); related.add((f[3], f[1]))
genes = lambda c, blocks: [g for g, n in EX.hits(c, blocks).most_common() if n >= 50 and g in G]
ASM = Index()
if a.asm:
    for l in open(a.asm):
        f = l.split('\t', 6)
        if len(f) > 5 and f[2] == 'exon': ASM.add(f[0], int(f[3]) - 1, int(f[4]), 1)

B = pysam.AlignmentFile(a.bam); R = collections.defaultdict(list)
for r in B:
    if r.is_unmapped or r.is_supplementary: continue
    p = r.reference_start; intr = []
    for op, k in r.cigartuples:
        if op == 3: intr.append((p, p + k))
        if op in (0, 2, 3, 7, 8): p += k
    R[r.query_name].append((r.is_secondary, r.reference_id, r.get_tag('ms'), r.reference_start, p, tuple(intr)))
name = B.references
def blocks(r):
    o = []; s = r[3]
    for x, y in r[5]: o.append((s, x)); s = y
    return o + [(s, r[4])]

rows = []; expr = collections.Counter()                         # intronless side: (contig, start, end, read, own, spliced-side gene)
for q, rs in R.items():
    P = [r for r in rs if not r[0]]
    if len(P) != 1: continue
    P = P[0]
    if P[5]:
        for g in genes(name[P[1]], blocks(P))[:1]: expr[g] += 1  # spliced primary reads per gene, for L3
    for S in rs:
        if not S[0] or S[2] < 0.9 * P[2] or bool(S[5]) == bool(P[5]): continue
        if S[1] == P[1] and S[3] < P[4] and P[3] < S[4]: continue
        U, Sp = (P, S) if not P[5] else (S, P)
        rows.append((name[U[1]], U[3], U[4], q, U is P, (genes(name[Sp[1]], blocks(Sp)) or [None])[0]))
rows.sort(key=lambda r: (r[0], r[1])); places = []
for r in rows:
    if places and places[-1][0] == r[0] and r[1] < places[-1][2]: places[-1][2] = max(places[-1][2], r[2]); places[-1][3].append(r)
    else: places.append([r[0], r[1], r[2], [r]])

planted = {}
if a.planted:
    for l in open(a.planted):
        f = l.rstrip('\n').split('\t')
        if f[4] == 'pseudogene': planted[f[2]] = f[1]            # gene -> unit
    unit_gene = {f.split('\t')[1]: f.split('\t')[2] for f in open(a.planted) if f.split('\t')[4] == 'gene'}
n = collections.Counter(); out = open(a.out + '.places.tsv', 'w'); labelled = set()
out.write('contig\tstart\tend\treads\town_reads\tprocessed_copy_of\tgenes_here\tclass\tlisted_parent\tassembled\n')
for c, s, e, rr in places:
    reads = {r[3] for r in rr}
    if len(reads) < 2: n['places with one read'] += 1; continue
    own = len({r[3] for r in rr if r[4]}); par = collections.Counter(r[5] for r in rr if r[5]).most_common(1)
    par = par[0][0] if par else None; here = genes(c, [(s, e)]); labelled.update(here)
    pp = [g for g in here if G[g][1] in PP]
    k = 'processed pseudogene' if pp else 'single-exon gene' if here and all(G[g][3] == 1 for g in here) else 'exon of a spliced gene' if here else 'nothing annotated'
    listed = sorted({G[g][2] for g in pp if G[g][2]}); asm = bool(a.asm) and sum(ASM.hits(c, [(s, e)]).values()) >= 50
    n['places'] += 1; n['class: ' + k] += 1; n['own reads: ' + ('none' if not own else 'under 5%' if own < 0.05 * len(reads) else '5-50%' if own < 0.5 * len(reads) else '50% or more')] += 1
    if a.asm: n[f'assembled today, {"with" if own else "without"} own reads'] += asm
    pn = G[par][0] if par else None
    if listed:
        n['L2 places with a listed parent'] += 1; n['L2 label = the listed parent'] += pn in listed
        n['L2 label = the parent or a listed relative'] += pn in listed or any((pn, x) in related for x in listed)
    if planted:
        src = collections.Counter(r[3].split('|')[1] for r in rr).most_common(1)[0][0]
        for g in here:
            if g in planted: n['L4 labelled planted pseudogenes'] += 1; n['L4 label = the planted source gene'] += par == unit_gene.get(planted[g]) == src
    out.write('\t'.join(map(str, (c, s, e, len(reads), own, pn, ','.join(G[g][0] for g in here), k, ','.join(listed), int(asm)))) + '\n')
by_name = collections.Counter()
for g, v in expr.items(): by_name[G[g][0]] += v
elig = [g for g, v in G.items() if v[1] in PP and v[2] and by_name[v[2]] >= 50]
n['L3 processed pseudogenes whose parent has >= 50 spliced primary reads'] = len(elig); n['L3 of them labelled'] = sum(g in labelled for g in elig)
if planted: n['planted pseudogenes'] = len(planted); n['planted pseudogenes labelled'] = sum(g in labelled for g in planted)
assert n['places'] == sum(v for k, v in n.items() if k.startswith('class: '))
for k, v in n.items(): print(f'{v:>9,}  {k}')
P = n['places'] or 1; ok = n['class: processed pseudogene'] + n['class: single-exon gene']
print(f'L1 {ok / P:.3f}  ({ok} of {n["places"]})')
d = n['L2 places with a listed parent'] or 1; print(f'L2 {n["L2 label = the listed parent"] / d:.3f} parent, {n["L2 label = the parent or a listed relative"] / d:.3f} parent or relative')
print(f'L3 {n["L3 of them labelled"] / (len(elig) or 1):.3f}')
if planted: print(f'L4 {n["L4 label = the planted source gene"] / (n["L4 labelled planted pseudogenes"] or 1):.3f}')
