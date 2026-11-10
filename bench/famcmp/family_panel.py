#!/usr/bin/env python3
"""PREREG r1521 (docs/PREREG_family_panel_tools_2026-11-03.md): what each single-sample assembler emits at the truth multi-copy
families (famtruth primates-level components). Descriptive; definitions fixed before the run.

  family_panel.py BAM AS_TABLE OUT_PREFIX TOOL=GTF[.gz] ...

Per truth gene, two expression definitions (Amendment 1: a member whose reads' primaries all landed on a paralog is invisible to a
primary-only definition, which would judge primary-only tools on exactly the copies their input can see):
  P = >= 2 molecules with a PRIMARY record overlapping its exons on its strand (mapq0 share = share of those at MAPQ 0);
  T = >= 2 molecules with ANY placement there (primary, or secondary at AS >= 0.98 x the molecule's genome-wide best, the good
      pool's tie set from AS_TABLE). T is a superset of P; T-only members are the copies no primary-only tool can reach.
Per tool per gene: hit = >= 1 model on the gene's strand overlapping its exons (Amendment 2, 260ac8c9: any model; FAMPANEL_SPLICED_ONLY=1
restores the original spliced-only rule).
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
        if len(e) < 2 and os.environ.get('FAMPANEL_SPLICED_ONLY'): continue  # the registered definition; Amendment 2 = any model
        e.sort(); by[(c, strand[(c, t)])].append((e[0][0], e[-1][1], e))
    for v in by.values(): v.sort()
    return by

MONO_SPAN = float(os.environ.get('PANEL_MONO_SPAN', '0'))  # r1535 span guard; 0 = off (the r1521 definition)

def overlaps_exons(models, gexons, lo, hi):
    out = []
    for s, e, ex in models:
        if e <= lo: continue
        if s >= hi: break
        if len(ex) < 2 and MONO_SPAN and e - s > MONO_SPAN * (hi - lo): continue  # r1535: a single-exon model hits only within MONO_SPAN x the gene span
        if any(a < y and x < b for a, b in ex for x, y in gexons): out.append((s, e, ex))
    return out

bam_path, as_table, out = sys.argv[1:4]; tools = [a.split('=', 1) for a in sys.argv[4:]]
best = {}
with open(as_table) as fh:
    for l in fh:
        f = l.rstrip('\n').split('\t')
        if len(f) > 3 and f[0][0] != '#' and f[3].isdigit() and int(f[3]) > 1 and f[1].isdigit(): best[f[0]] = int(f[1])
G = F.load_genes(os.environ.get('PANEL_SPECIES', 'human')); pairs, universe = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['primates'])); fams = F.components(pairs)
if os.environ.get('PANEL_CONTIGS'):  # r1531: judge only families whose every gene lies on these contigs (sliced runs)
    keep = set(os.environ['PANEL_CONTIGS'].split(',')); fams = [f for f in fams if all(G[g]['chrom'] in keep for g in f)]
bam = pysam.AlignmentFile(bam_path)
expr = {}
for f in fams:
    for g in f:
        v = G[g]; m0 = 0; seen_p = set(); seen_t = set()
        for r in bam.fetch(v['chrom'], v['exons'][0][0], v['exons'][-1][1]):
            if r.is_supplementary or r.is_unmapped: continue
            if (r.is_reverse and v['strand'] == '+') or (not r.is_reverse and v['strand'] == '-'): continue
            if not any(a < y and x < b for a, b in r.get_blocks() for x, y in v['exons']): continue
            if r.is_secondary:
                b = best.get(r.query_name)
                if b and r.has_tag('AS') and r.get_tag('AS') >= 0.98 * b: seen_t.add(r.query_name)
            elif r.query_name not in seen_p:
                seen_p.add(r.query_name); seen_t.add(r.query_name); m0 += r.mapping_quality == 0
        n = len(seen_p); expr[g] = (n, m0 / n if n else 0.0, len(seen_t))
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
    w.write('family\tgene\tname\tchrom\tprimary_molecules\tmapq0_share\ttied_molecules\t' + '\t'.join(t for t, _ in tools) + '\n')
    for i, f in enumerate(fams):
        for g in f:
            if expr[g][2] >= 2: w.write(f"F{i}\t{g}\t{G[g]['name']}\t{G[g]['chrom']}\t{expr[g][0]}\t{expr[g][1]:.2f}\t{expr[g][2]}\t" + '\t'.join(str(hit[t][g]) for t, _ in tools) + '\n')
print(f"truth families {len(fams)} (genes {sum(len(f) for f in fams)})")
with open(f'{out}.families.tsv', 'w') as w:
    w.write('definition\tfamily\tgenes\texpressed\t' + '\t'.join(f'{t}_class\t{t}_fused' for t, _ in tools) + '\n')
    for dname, col in (('P (primary)', 0), ('T (primary or tied secondary)', 2)):
        judged = [f for f in fams if sum(expr[g][col] >= 2 for g in f) >= 2]
        nmem = sum(sum(expr[g][col] >= 2 for g in f) for f in judged)
        tonly = sum(sum(expr[g][2] >= 2 and expr[g][0] < 2 for g in f) for f in judged)
        print(f"\n== definition {dname}: judged {len(judged)} families, {nmem} expressed members" + (f" ({tonly} reachable only by a tied secondary)" if col == 2 else ""))
        print("tool\tcomplete\tpartial\tnone\tfused_families\tmember_hit_rate\thit_rate_mapq0_share>=0.5\thit_rate_mapq0_share<0.5" + ("\thit_rate_T-only_members" if col == 2 else ""))
        rows = {}
        for t, _ in tools:
            c = collections.Counter(); hm = [0, 0]; hl = [0, 0]; ht = [0, 0]
            for f in judged:
                ex = [g for g in f if expr[g][col] >= 2]; h = [hit[t][g] > 0 for g in ex]
                cls = 'complete' if all(h) else 'none' if not any(h) else 'partial'; c[cls] += 1
                rows.setdefault(id(f), []).append((cls, id(f) in fused[t]))
                for g, hh in zip(ex, h):
                    if expr[g][0] < 2: ht[0] += hh; ht[1] += 1
                    else: (hm if expr[g][1] >= 0.5 else hl)[0] += hh; (hm if expr[g][1] >= 0.5 else hl)[1] += 1
            n = hm[1] + hl[1] + ht[1]
            print(f"{t}\t{c['complete']}\t{c['partial']}\t{c['none']}\t{sum(id(f) in fused[t] for f in judged)}\t{(hm[0] + hl[0] + ht[0]) / max(n, 1):.3f}\t{hm[0] / max(hm[1], 1):.3f} (n={hm[1]})\t{hl[0] / max(hl[1], 1):.3f} (n={hl[1]})" + (f"\t{ht[0] / max(ht[1], 1):.3f} (n={ht[1]})" if col == 2 else ""))
        for i, f in enumerate(fams):
            if f in judged: w.write(f"{dname[0]}\tF{i}\t{','.join(G[g]['name'] for g in f)}\t{sum(expr[g][col] >= 2 for g in f)}\t" + '\t'.join(f'{c}\t{int(u)}' for c, u in rows[id(f)]) + '\n')
