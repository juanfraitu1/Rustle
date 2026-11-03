#!/usr/bin/env python3
"""PREREG r1521 (docs/PREREG_family_panel_tools_2026-11-03.md): what each single-sample assembler emits at the truth multi-copy
families (famtruth primates-level components). Descriptive; definitions fixed before the run.

  family_panel.py BAM OUT_PREFIX TOOL=GTF[.gz] ...

Per truth gene: expressed = >= 2 PRIMARY records overlapping its exons on its strand; mapq0 share = share of those at MAPQ 0.
Per tool per gene: hit = >= 1 spliced model (>= 2 exons) on the gene's strand overlapping its exons.
Per tool per family (families with >= 2 expressed members, judged on the expressed members only): complete (every expressed
member hit) / partial / none; fused = >= 1 model of the tool overlapping exons of >= 2 members of the family on their strand.
OUT_PREFIX.genes.tsv (one row per expressed truth gene x tool), OUT_PREFIX.families.tsv, summary on stdout."""
import sys, os, re, gzip, collections
import pysam
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F

def read_models(path):
    ex = collections.defaultdict(list); strand = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[(f[0], t)].append((int(f[3]) - 1, int(f[4]))); strand[(f[0], t)] = f[6]
    by = collections.defaultdict(list)  # (chrom, strand) -> [(start, end, exons)]
    for (c, t), e in ex.items():
        if len(e) < 2: continue
        e.sort(); by[(c, strand[(c, t)])].append((e[0][0], e[-1][1], e))
    for v in by.values(): v.sort()
    return by

def overlaps_exons(models, gexons, lo, hi):
    out = []
    for s, e, ex in models:
        if e <= lo: continue
        if s >= hi: break
        if any(a < y and x < b for a, b in ex for x, y in gexons): out.append((s, e, ex))
    return out

bam_path, out = sys.argv[1:3]; tools = [a.split('=', 1) for a in sys.argv[3:]]
G = F.load_cat(); pairs, universe = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['primates'])); fams = F.components(pairs)
bam = pysam.AlignmentFile(bam_path)
expr = {}
for f in fams:
    for g in f:
        v = G[g]; n = m0 = 0; seen = set()
        for r in bam.fetch(v['chrom'], v['exons'][0][0], v['exons'][-1][1]):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_name in seen: continue
            if (r.is_reverse and v['strand'] == '+') or (not r.is_reverse and v['strand'] == '-'): continue
            if not any(a < y and x < b for a, b in r.get_blocks() for x, y in v['exons']): continue
            seen.add(r.query_name); n += 1; m0 += r.mapping_quality == 0
        expr[g] = (n, m0 / n if n else 0.0)
models = {t: read_models(p) for t, p in tools}
hit = {t: {} for t, _ in tools}; fused = {t: set() for t, _ in tools}
for f in fams:
    for t, _ in tools:
        touched = collections.defaultdict(set)  # model id -> members
        for g in f:
            v = G[g]
            ms = overlaps_exons(models[t].get((v['chrom'], v['strand']), []), v['exons'], v['exons'][0][0], v['exons'][-1][1])
            hit[t][g] = len(ms)
            for s, e, ex in ms: touched[(v['chrom'], v['strand'], s, e)].add(g)
        if any(len(gs) >= 2 for gs in touched.values()): fused[t].add(id(f))
with open(f'{out}.genes.tsv', 'w') as w:
    w.write('family\tgene\tname\tchrom\tprimary_reads\tmapq0_share\t' + '\t'.join(t for t, _ in tools) + '\n')
    for i, f in enumerate(fams):
        for g in f:
            if expr[g][0] >= 2: w.write(f"F{i}\t{g}\t{G[g]['name']}\t{G[g]['chrom']}\t{expr[g][0]}\t{expr[g][1]:.2f}\t" + '\t'.join(str(hit[t][g]) for t, _ in tools) + '\n')
judged = [f for f in fams if sum(expr[g][0] >= 2 for g in f) >= 2]
print(f"truth families {len(fams)} (genes {sum(len(f) for f in fams)}); judged = >= 2 expressed members: {len(judged)} families, {sum(sum(expr[g][0] >= 2 for g in f) for f in judged)} expressed members")
print("tool\tcomplete\tpartial\tnone\tfused_families\tmember_hit_rate\thit_rate_mapq0_share>=0.5\thit_rate_mapq0_share<0.5")
with open(f'{out}.families.tsv', 'w') as w:
    w.write('family\tgenes\texpressed\t' + '\t'.join(f'{t}_class\t{t}_fused' for t, _ in tools) + '\n')
    rows = {}
    for t, _ in tools:
        c = collections.Counter(); hm = [0, 0]; hl = [0, 0]
        for f in judged:
            ex = [g for g in f if expr[g][0] >= 2]; h = [hit[t][g] > 0 for g in ex]
            cls = 'complete' if all(h) else 'none' if not any(h) else 'partial'; c[cls] += 1
            rows.setdefault(id(f), []).append((cls, id(f) in fused[t]))
            for g, hh in zip(ex, h):
                (hm if expr[g][1] >= 0.5 else hl)[0] += hh; (hm if expr[g][1] >= 0.5 else hl)[1] += 1
        n = hm[1] + hl[1]
        print(f"{t}\t{c['complete']}\t{c['partial']}\t{c['none']}\t{sum(id(f) in fused[t] for f in judged)}\t{(hm[0] + hl[0]) / max(n, 1):.3f}\t{hm[0] / max(hm[1], 1):.3f} (n={hm[1]})\t{hl[0] / max(hl[1], 1):.3f} (n={hl[1]})")
    for i, f in enumerate(fams):
        if f in judged: w.write(f"F{i}\t{','.join(G[g]['name'] for g in f)}\t{sum(expr[g][0] >= 2 for g in f)}\t" + '\t'.join(f'{c}\t{int(u)}' for c, u in rows[id(f)]) + '\n')
