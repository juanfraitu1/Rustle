#!/usr/bin/env python3
"""PREREG r1772 (docs/PREREG_sibling_reach_2026-11-17.md): are the siblings of an expressed truth-family member reached by
its multimapping reads? Descriptive; definitions fixed before the run.

  sibling_reach.py BAM FASTA OUT_PREFIX [OURS.gtf]      (env as famtruth; PANEL_SPECIES=human|gorilla;
                                                         SIB_EXCLUDE=chrom:lo-hi drops genes inside; SIB_MAX_FAMS=N smoke test;
                                                         SIB_REVERSE=1 adds the r1773 tables)

Per truth family (simiiformes components) and member g: source molecules = primaries on the OTHER members; reached(K, band) =
>= K source molecules with a secondary record on g's exons and strand at AS >= band x that molecule's primary AS.
Identity of g = best spliced-exon identity to a member with >= 2 primaries (seqscore's minimap2 rule, shorter covered >= 0.5).
OUT_PREFIX.siblings.tsv (one row per member of a family with >= 1 source molecule), summary tables on stdout."""
import sys, os, re, gzip, collections, subprocess, tempfile
import pysam
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F

bam_path, fasta, out = sys.argv[1:4]; gtf = sys.argv[4] if len(sys.argv) > 4 else None
MIN_SRC = 10; BANDS = (0.0, 0.5, 0.98); KS = (1, 2, 5)
G = F.load_genes(os.environ.get('PANEL_SPECIES', 'human'))
pairs, _ = F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['simiiformes'])); fams = F.components(pairs)
if os.environ.get('SIB_EXCLUDE'):
    c, lo, hi = re.match(r'(.+):(\d+)-(\d+)$', os.environ['SIB_EXCLUDE']).groups()
    inside = lambda g: G[g]['chrom'] == c and G[g]['exons'][0][0] < int(hi) and int(lo) < G[g]['exons'][-1][1]
    fams = [f - {g for g in f if inside(g)} for f in fams]
fams = [f for f in fams if len(f) >= 2 and all(G[g]['exons'] for g in f)]
if os.environ.get('SIB_MAX_FAMS'): fams = sorted(fams, key=lambda f: min(f))[:int(os.environ['SIB_MAX_FAMS'])]

bam = pysam.AlignmentFile(bam_path)
prim = {}; sec = {}  # gene -> {molecule: primary AS} / {molecule: best secondary AS here}
for f in fams:
    for g in f:
        v = G[g]; p = prim[g] = {}; s = sec[g] = {}
        if v['chrom'] not in bam.references: continue
        for r in bam.fetch(v['chrom'], v['exons'][0][0], v['exons'][-1][1]):
            if r.is_supplementary or r.is_unmapped or r.is_reverse == (v['strand'] == '+'): continue
            if not any(a < y and x < b for a, b in r.get_blocks() for x, y in v['exons']): continue
            a_s = r.get_tag('AS') if r.has_tag('AS') else 0
            if r.is_secondary: s[r.query_name] = max(a_s, s.get(r.query_name, 0))
            else: p[r.query_name] = a_s

# identity of each member to the family's expressed members (seqscore's alignment, within the truth families only)
members = sorted({g for f in fams for g in f}); fam_of = {g: i for i, f in enumerate(fams) for g in f}
ident = {}
with tempfile.TemporaryDirectory() as d:
    fa = os.path.join(d, 'g.fa'); ref = pysam.FastaFile(fasta)
    with open(fa, 'w') as w:
        for g in members:
            v = G[g]
            if v['chrom'] not in ref.references: continue
            q = ''.join(ref.fetch(v['chrom'], a, b) for a, b in sorted(v['exons'])).upper()
            if v['strand'] == '-': q = q.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            if len(q) >= 100: w.write(f">{g}\n{q}\n")
    paf = subprocess.run(['minimap2', '-c', '-x', 'asm20', '-N', '500', '-p', '0.05', '--secondary=yes', '-t', str(os.cpu_count()), fa, fa],
                         capture_output=True, text=True, check=True).stdout
for line in paf.split('\n'):
    f = line.split('\t')
    if len(f) < 12 or f[0] == f[5] or fam_of[f[0]] != fam_of[f[5]] or f[4] != '+': continue
    ql, tl = int(f[1]), int(f[6]); bl = int(f[10])
    cov = (int(f[3]) - int(f[2]) if ql <= tl else int(f[8]) - int(f[7])) / min(ql, tl)
    if cov < 0.5 or not bl: continue
    for a, b in ((f[0], f[5]), (f[5], f[0])):  # a's identity to an EXPRESSED relative b
        if len(prim[b]) >= 2: ident[a] = max(ident.get(a, 0), int(f[9]) / bl)

models = collections.defaultdict(list)
if gtf:
    ex = collections.defaultdict(list)
    for l in (gzip.open if gtf.endswith('.gz') else open)(gtf, 'rt'):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    for (c, st, _), e in ex.items(): models[(c, st)].extend(e)
    for k, v in models.items():  # merged, so one bisect answers "does any model exon overlap"
        v.sort(); m = [list(v[0])]
        for x, y in v[1:]:
            if x <= m[-1][1]: m[-1][1] = max(m[-1][1], y)
            else: m.append([x, y])
        models[k] = m
import bisect
def has_model(v):
    m = models.get((v['chrom'], v['strand']), [])
    for x, y in v['exons']:
        i = bisect.bisect_left(m, [y, 0])
        if i and m[i - 1][1] > x: return True
    return False

BIN = lambda i: 'no alignment' if i is None else '>=0.98' if i >= 0.98 else '0.95-0.98' if i >= 0.95 else '0.90-0.95' if i >= 0.90 else '0.85-0.90' if i >= 0.85 else '<0.85'
rows = []
with open(f'{out}.siblings.tsv', 'w') as w:
    w.write('family\tgene\tname\tchrom\tfamily_size\town_primaries\tsource_molecules\tidentity\tbin\t' + '\t'.join(f'sec_band{b}' for b in BANDS) + '\tsec_any_molecule\tmodel\n')
    for i, f in enumerate(fams):
        for g in f:
            src = {}
            for h in f:
                if h != g: src.update(prim[h])
            if not src: continue
            n = [sum(1 for m, a in sec[g].items() if m in src and a >= b * src[m]) for b in BANDS]
            r = dict(own=len(prim[g]), src=len(src), bin=BIN(ident.get(g)), n=n, model=has_model(G[g]) if gtf else None)
            rows.append(r)
            w.write(f"F{i}\t{g}\t{G[g]['name']}\t{G[g]['chrom']}\t{len(f)}\t{r['own']}\t{r['src']}\t{ident.get(g, 0):.4f}\t{r['bin']}\t" + '\t'.join(map(str, n)) + f"\t{len(sec[g])}\t{'' if r['model'] is None else int(r['model'])}\n")

print(f"truth families {len(fams)}, members {len(members)}, members with >= 1 source molecule {len(rows)}, judged (>= {MIN_SRC} source molecules) {sum(r['src'] >= MIN_SRC for r in rows)}")
order = ['>=0.98', '0.95-0.98', '0.90-0.95', '0.85-0.90', '<0.85', 'no alignment', 'all']
sh = lambda k, n: f"{k / n:.3f}" if n else '-'
for title, keep in (('siblings with < 2 own primaries', lambda r: r['own'] < 2), ('siblings with >= 2 own primaries', lambda r: r['own'] >= 2)):
    print(f"\n== {title}, family has >= {MIN_SRC} source molecules; reached = >= K source molecules with a secondary here")
    print('identity_bin\tn\t' + '\t'.join(f'K{k}_any' for k in KS) + '\tK2_band0.5\tK2_band0.98\tmodel_ours\tmodel_ours_if_K2_any\tmodel_ours_if_not')
    for b in order:
        R = [r for r in rows if r['src'] >= MIN_SRC and keep(r) and b in (r['bin'], 'all')]; n = len(R)
        hit = [r for r in R if r['n'][0] >= 2]; miss = [r for r in R if r['n'][0] < 2]
        print(f"{b}\t{n}\t" + '\t'.join(sh(sum(r['n'][0] >= k for r in R), n) for k in KS) + f"\t{sh(sum(r['n'][1] >= 2 for r in R), n)}\t{sh(sum(r['n'][2] >= 2 for r in R), n)}"
              + (f"\t{sh(sum(r['model'] for r in R), n)}\t{sh(sum(r['model'] for r in hit), len(hit))} (n={len(hit)})\t{sh(sum(r['model'] for r in miss), len(miss))} (n={len(miss)})" if gtf else '\t-\t-\t-'))
print(f"\n== by source depth (all judged or not), siblings with < 2 own primaries: reach K2_any")
for lo, hi in ((1, 9), (10, 99), (100, 999), (1000, 10 ** 9)):
    R = [r for r in rows if lo <= r['src'] <= hi and r['own'] < 2]
    print(f"source molecules {lo}-{hi if hi < 10 ** 9 else ''}\tn={len(R)}\t{sh(sum(r['n'][0] >= 2 for r in R), len(R))}")

if os.environ.get('SIB_REVERSE'):
    # PREREG r1773 (docs/PREREG_sibling_reverse_2026-11-17.md): which loci do a family's multimapping reads reach?
    judged = {i for i, f in enumerate(fams) if len({m for g in f for m in prim[g]}) >= MIN_SRC}
    mol = {}  # molecule -> (family, primary AS, members holding its primary)
    for i in judged:
        for g in fams[i]:
            for m, a in prim[g].items():
                if m in mol and mol[m][0] != i: continue  # a primary over members of two families: the first family keeps it
                mol[m] = (i, a, mol[m][2] | {g} if m in mol else {g})
    idx = F.exon_index(G)
    older = collections.defaultdict(set)
    for a, b in F.truth_pairs(G, F.PARALOGS, set(F.LEVELS['eutheria']))[0]:
        for x, y in ((a, b), (b, a)):
            if x in fam_of and fam_of[x] != fam_of.get(y): older[fam_of[x]].add(y)
    def klass(i, g):
        bt = G[g].get('biotype', '')
        return 'member' if fam_of.get(g) == i else 'older paralog' if g in older[i] else 'pseudogene' if 'pseudogene' in bt else 'other coding' if bt == 'protein_coding' else 'other non-coding'
    RANK = {'member': 0, 'older paralog': 1, 'pseudogene': 2, 'other coding': 3, 'other non-coding': 4}
    hits = collections.defaultdict(dict)  # (family, class, target) -> {molecule: best AS ratio}
    proc = subprocess.Popen(['samtools', 'view', '-@', '8', '-f', '0x100', '-F', '0x800', bam_path], stdout=subprocess.PIPE, text=True)
    for l in proc.stdout:
        name = l[:l.index('\t')]
        if name not in mol: continue
        f = l.split('\t', 12); i, pas, pg = mol[name]
        strand = '-' if int(f[1]) & 16 else '+'; pos = int(f[3]) - 1; blocks = []; s0 = pos
        for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
            n = int(n)
            if op in 'M=XD': pos += n
            elif op == 'N': blocks.append((s0, pos)); pos += n; s0 = pos
        blocks.append((s0, pos))
        a_s = re.search(r'\tAS:i:(-?\d+)', l); ratio = int(a_s.group(1)) / pas if a_s and pas else 0
        genes = set()
        if (f[2], strand) in idx:
            v, st, mx = idx[(f[2], strand)]
            for x, y in blocks:
                j = bisect.bisect_left(st, x - mx)
                while j < len(v) and v[j][0] < y:
                    if v[j][1] > x: genes.add(v[j][2])
                    j += 1
        if genes & pg: continue  # on the member that holds its own primary
        if genes:
            g = min(genes, key=lambda g: (RANK[klass(i, g)], g)); key = (i, klass(i, g), g)
        else: key = (i, 'unannotated', (f[2], strand, blocks[0][0] // 1000))
        d = hits[key]; d[name] = max(ratio, d.get(name, 0))
    proc.wait()
    # unannotated: adjacent 1 kb bins of one family, chromosome and strand are one locus
    un = sorted(k for k in hits if k[1] == 'unannotated'); last = None
    for k in un:
        if last and k[0] == last[0] and k[2][:2] == last[2][:2] and k[2][2] - prev <= 1:
            for m, r in hits.pop(k).items(): hits[last][m] = max(r, hits[last].get(m, 0))
        else: last = k
        prev = k[2][2]
    nsrc = collections.Counter(v[0] for v in mol.values())
    with open(f'{out}.reverse.tsv', 'w') as w:
        w.write('family\tclass\ttarget\tfamily_source_molecules\t' + '\t'.join(f'mol_band{b}' for b in BANDS) + '\n')
        T = []
        for (i, c, t), d in hits.items():
            n = [sum(r >= b for r in d.values()) for b in BANDS]; T.append((i, c, n))
            w.write(f"F{i}\t{c}\t{G[t]['name'] if c != 'unannotated' else '%s:%s:%d' % (t[0], t[1], t[2] * 1000)}\t{nsrc[i]}\t" + '\t'.join(map(str, n)) + '\n')
    classes = ['member', 'older paralog', 'pseudogene', 'other coding', 'other non-coding', 'unannotated']
    quiet = sum(1 for i in judged for g in fams[i] if len(prim[g]) < 2)  # members without own expression in judged families
    print(f"\n== r1773 reverse: targets reached by the secondaries of a family's molecules ({len(judged)} families with >= {MIN_SRC} source molecules, {len(mol)} molecules)")
    print('band\tK\t' + '\t'.join(classes) + '\tall\tmember_share\trelative_share\tmembers_without_own_primaries_reached')
    for bi, b in enumerate(BANDS):
        for kname, need in (('2', lambda i: 2), ('5', lambda i: 5), ('5%', lambda i: max(2, 0.05 * nsrc[i]))):
            c = collections.Counter(cl for i, cl, n in T if n[bi] >= need(i)); tot = sum(c.values())
            q = sum(1 for (i, cl, t), d in hits.items() if cl == 'member' and len(prim[t]) < 2 and sum(r >= b for r in d.values()) >= need(i))
            print(f"{b or 'any'}\t{kname}\t" + '\t'.join(str(c[x]) for x in classes) + f"\t{tot}\t{sh(c['member'], tot)}\t{sh(c['member'] + c['older paralog'] + c['pseudogene'], tot)}\t{q} of {quiet} ({sh(q, quiet)})")
