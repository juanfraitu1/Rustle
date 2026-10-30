#!/usr/bin/env python3
"""Sequence-defined family score (2026-10-30, user: "score the families based on what they really share", "Soto only works for up to
98% similarity of SDs"). No Compara. Truth at identity tier T: two EXPRESSED annotated genes (any biotype) whose spliced exon
sequences align (minimap2 -c -x asm20, all-vs-all) with the SHORTER covered >= --cov (0.5) at identity >= T. Predicted = gene
pairs inside one family (copies mapped to genes, famtruth.map_copy). Per tier T in --tiers: precision = predicted pairs that are
truth at T / predicted pairs that align at all or not (every predicted pair is judged: a pair below every tier is FALSE), recall =
truth pairs at T recovered. A predicted pair is 'ok' for precision at tier T when it meets tier T; so precision falls as T rises
and recall rises as T rises (fewer, closer truth pairs).
  seqscore.py --species human|gorilla --gtf BASE.families.gtf --fasta GENOME [--cov 0.5] [--tiers 0.85,0.90,0.95,0.98]
              [--cache FILE] COPIES_TSV...        (env as famtruth)"""
import sys, os, csv, collections, itertools, subprocess, tempfile, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
opt = lambda k, d: a[a.index(k) + 1] if k in a else d
species = opt('--species', 'gorilla'); gtf = opt('--gtf', None); FASTA = opt('--fasta', None); COV = float(opt('--cov', '0.5'))
TIERS = [float(x) for x in opt('--tiers', '0.85,0.90,0.95,0.98').split(',')]; CACHE = opt('--cache', None)
# Structural RNA arrays (rDNA RNA45S/RNA28S/RNA18S: 434 expressed copies at >= 98% on A119b, 94k pairs) are tandem repeat
# units, not gene families; they are dropped from the universe unless --keep-rna (2026-10-30).
RNA_BT = {'rRNA', 'rRNA_pseudogene', 'tRNA', 'Mt_tRNA', 'Mt_rRNA', 'snRNA', 'snoRNA', 'scaRNA', 'misc_RNA', 'miRNA', 'vault_RNA', 'sRNA'}
flags = {'--species', '--gtf', '--fasta', '--cov', '--tiers', '--cache'}
paths = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or a[i - 1] not in flags)]
if species == 'hg002':  # as pairscore.py: Ensembl rapid-release GFF3 + xref via guided_eval/guided3
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
    import guided3 as g3
    g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']; F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
    G = g3.load_genes('hg002')
else:
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
idx = F.exon_index(G); expressed = F.expressed_from_gtf(gtf, idx)
genes = sorted(g for g in expressed if G[g]['exons'] and ('--keep-rna' in a or G[g].get('biotype') not in RNA_BT))
best = None
if CACHE and os.path.exists(CACHE): best = pickle.load(open(CACHE, 'rb'))
if best is None:
    by_chrom = collections.defaultdict(list)
    for g in genes: by_chrom[G[g]['chrom']].append(g)
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'g.fa'); n = 0
        with open(fa, 'w') as w:
            for chrom, gs in by_chrom.items():
                raw = subprocess.run(['samtools', 'faidx', FASTA, chrom], capture_output=True, text=True).stdout
                seq = ''.join(raw.split('\n')[1:])
                for g in gs:
                    s = ''.join(seq[s0:e] for s0, e in sorted(G[g]['exons'])).upper()
                    if G[g]['strand'] == '-': s = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                    if len(s) >= 100: w.write(f">{g}\n{s}\n"); n += 1
        print(f"# {n} expressed genes aligned all-vs-all", file=sys.stderr)
        paf = subprocess.run(['minimap2', '-c', '-x', 'asm20', '-N', '500', '-p', '0.05', '--secondary=yes', '-t', '16', fa, fa], capture_output=True, text=True).stdout
    best = {}
    for line in paf.split('\n'):
        f = line.split('\t')
        if len(f) < 12 or f[0] == f[5]: continue
        q, t = f[0], f[5]; ql, tl = int(f[1]), int(f[6]); alen = int(f[3]) - int(f[2]); m, bl = int(f[9]), int(f[10])
        cov = alen / min(ql, tl); idn = m / bl if bl else 0
        k = (min(q, t), max(q, t))
        if cov >= COV and (k not in best or idn > best[k]): best[k] = idn
    if CACHE: pickle.dump(best, open(CACHE, 'wb'))
exp = set(genes)
def overlapping(x, y):  # genes sharing genomic sequence (nested / readthrough / antisense annotations) are not a family pair
    return G[x]['chrom'] == G[y]['chrom'] and G[x]['exons'][0][0] < G[y]['exons'][-1][1] and G[y]['exons'][0][0] < G[x]['exons'][-1][1]
best = {k: i for k, i in best.items() if k[0] in exp and k[1] in exp and not overlapping(*k)}
truth = {T: {k for k, i in best.items() if i >= T} for T in TIERS}
# GENE-AVERAGED precision / recall (2026-10-30): pair counts grow with the square of a family's size, so one 300-copy repeat family
# decides a pair score (A119b: 57,628 truth pairs at 98%, one component of hundreds of genes). Per gene g: recall_g = truth
# partners of g inside g's family / truth partners of g; precision_g = predicted partners of g that are truth / predicted
# partners; averaged over genes that have truth partners (recall) or predicted partners (precision).
def gene_avg(pred, tr):
    pp = collections.defaultdict(set); tp = collections.defaultdict(set)
    for x, y in pred: pp[x].add(y); pp[y].add(x)
    for x, y in tr: tp[x].add(y); tp[y].add(x)
    rec = [len(pp[g] & tp[g]) / len(tp[g]) for g in tp if tp[g]]
    prec = [len(pp[g] & tp[g]) / len(pp[g]) for g in pp if pp[g]]
    return (sum(prec) / len(prec) if prec else 0, sum(rec) / len(rec) if rec else 0, len(tp))
print("run\tpred_pairs\t" + "\t".join(f"P@{T:.2f}\tR@{T:.2f} (truth {len(truth[T])})\tgeneP@{T:.2f}\tgeneR@{T:.2f}" for T in TIERS))
for p in paths:
    fams, _, _ = F.load_catalog(p, G, idx)
    pred = {(x, y) for v in fams.values() for x, y in itertools.combinations(sorted(v), 2) if x in exp and y in exp and not overlapping(x, y)}
    cells = []
    for T in TIERS:
        ok = len(pred & truth[T]); gp, gr, ng = gene_avg(pred, truth[T])
        cells.append(f"{ok / len(pred) if pred else 0:.3f}\t{ok / len(truth[T]) if truth[T] else 0:.3f}\t{gp:.3f}\t{gr:.3f}")
    print(f"{os.path.basename(p)}\t{len(pred)}\t" + "\t".join(cells))
