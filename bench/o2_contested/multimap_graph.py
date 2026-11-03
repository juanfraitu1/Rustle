"""r1511 + r1512 (docs/PREREG_secondary_usefulness_2026-11-03.md): where secondaries land by ms-ratio bin, and the gene graph
their multimapping draws.  Usage: multimap_graph.py RECORDS_SORTED_BY_NAME COPIES_TSV OUT_PREFIX  (truth env as famtruth)."""
import sys, os, bisect, collections, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F

records, copies_tsv, out = sys.argv[1:4]
K = 3
BINS = [0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.01]
LAB = ['<0.5', '0.5-0.8', '0.8-0.9', '0.9-0.95', '0.95-0.99', '>=0.99']
G = F.load_genes('gorilla')
idx = F.exon_index(G)
truth, universe = F.truth_pairs(G, F.need(F.PARALOGS, 'DENOVO_PARALOGS'), F.levels('primates'))
truth = set(truth)

def genes_at(chrom, s, e):
    """genes with an exon overlapping [s, e) on either strand"""
    out = set()
    for strand in '+-':
        if (chrom, strand) not in idx: continue
        v, st, m = idx[(chrom, strand)]
        i = bisect.bisect_left(st, s - m)
        while i < len(v) and v[i][0] < e:
            if v[i][1] > s: out.add(v[i][2])
            i += 1
    return out

cat = collections.defaultdict(list)
for line in open(copies_tsv):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'family_id': h = {k: i for i, k in enumerate(f)}; continue
    cat[f[h['chrom']]].append((int(f[h['locus_start']]) - 50000, int(f[h['locus_end']]) + 50000, f[0]))
cat = {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in cat.items()}
def fam_at(chrom, s, e):
    if chrom not in cat: return set()
    v, st = cat[chrom]; i = bisect.bisect_left(st, s - 200000); out = set()
    while i < len(v) and v[i][0] < e:
        if v[i][1] > s: out.add(v[i][2])
        i += 1
    return out

def bin_of(r): return min(bisect.bisect_right(BINS, r) - 1, len(LAB) - 1)
use = collections.defaultdict(collections.Counter)   # bin -> class counts
asr = collections.defaultdict(list)                   # bin -> AS ratios (summary only)
edge = collections.defaultdict(collections.Counter)   # (gA, gB) -> Counter(max bin per molecule)
prim = collections.Counter()                          # gene -> primaries
n_mol = n_mult = 0
def flush(name, rows):
    global n_mol, n_mult
    n_mol += 1
    P = [r for r in rows if r[1] == '0']
    if len(P) != 1: return
    _, _, pc, pp, pe, pas, pms, _ = P[0]
    pms, pas = int(pms), int(pas); pp, pe = int(pp), int(pe)
    Pg = genes_at(pc, pp - 1, pe)
    for g in Pg: prim[g] += 1
    S = [r for r in rows if r[1] == '1']
    if not S or pms <= 0: return
    n_mult += 1
    Pf = fam_at(pc, pp - 1, pe)
    best = {}
    for _, _, sc, sp, se, sas, sms, _ in S:
        sp, se = int(sp), int(se); r = int(sms) / pms; b = bin_of(r)
        asr[b].append(int(sas) / pas if pas > 0 else 0.0)
        Sg = genes_at(sc, sp - 1, se); Sf = fam_at(sc, sp - 1, se)
        c = use[b]; c['records'] += 1
        par = any((min(p, s), max(p, s)) in truth for p in Pg for s in Sg if p != s)
        same = bool(Pf & Sf)
        if par: c['paralog'] += 1
        elif Sg - Pg: c['other_gene'] += 1
        elif Sg: c['same_gene'] += 1
        else: c['no_gene'] += 1
        if same: c['catalog_same_family'] += 1
        if par or same: c['useful'] += 1
        for p in Pg:
            for s in Sg:
                if p != s:
                    k = (min(p, s), max(p, s)); best[k] = max(best.get(k, -1), b)
    for k, b in best.items(): edge[k][b] += 1

with open(records) as fh:
    hdr = fh.readline()
    cur, rows = None, []
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[0] != cur:
            if cur is not None: flush(cur, rows)
            cur, rows = f[0], []
        rows.append(f)
    if cur is not None: flush(cur, rows)

with open(out + '.usefulness.tsv', 'w') as w:
    w.write('ms_ratio_bin\trecords\tparalog\tcatalog_same_family\tuseful\tuseful_frac\tother_gene\tsame_gene\tno_gene\tas_ratio_median\n')
    for b, lab in enumerate(LAB):
        c = use[b]; a = sorted(asr[b]); med = a[len(a) // 2] if a else float('nan')
        w.write(f"{lab}\t{c['records']}\t{c['paralog']}\t{c['catalog_same_family']}\t{c['useful']}\t{c['useful'] / c['records'] if c['records'] else 0:.4f}\t{c['other_gene']}\t{c['same_gene']}\t{c['no_gene']}\t{med:.3f}\n")

coding = lambda g: G[g]['biotype'] == 'protein_coding'
with open(out + '.graph.tsv', 'w') as w:
    w.write('min_bin\tedges\tgenes\tcomponents\tsilent_genes_joined\tsilent_coding_joined\tsilent_coding_truth_pair_in_component\tcopies_file\n')
    for b, lab in enumerate(LAB):
        pairs = [k for k, c in edge.items() if sum(n for bb, n in c.items() if bb >= b) >= K]
        comps = F.components(pairs)
        fn = f"{out}.R{lab.replace('>=', 'ge').replace('<', 'lt')}.copies.tsv"
        silent = silent_c = silent_ok = 0
        with open(fn, 'w') as cw:
            cw.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n')
            for i, comp in enumerate(sorted(comps, key=lambda s: sorted(s))):
                expressed = {g for g in comp if prim[g] > 0}
                for j, g in enumerate(sorted(comp)):
                    d = G[g]; ex = sorted(d['exons'])
                    cw.write(f"C{i}\t{j}\t{g}\t{d['chrom']}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{d['strand']}\t{prim[g]}\t{','.join(f'{s}-{e}' for s, e in ex)}\n")
                    if prim[g] == 0:
                        silent += 1
                        if coding(g):
                            silent_c += 1
                            if any((min(g, x), max(g, x)) in truth for x in expressed): silent_ok += 1
        genes = {g for k in pairs for g in k}
        w.write(f"{lab}\t{len(pairs)}\t{len(genes)}\t{len(comps)}\t{silent}\t{silent_c}\t{silent_ok}\t{fn}\n")
print(f"molecules {n_mol} with_secondaries {n_mult} genes_with_primaries {sum(1 for g in prim if prim[g] > 0)} silent_coding {sum(1 for g in G if coding(g) and prim[g] == 0)}")
