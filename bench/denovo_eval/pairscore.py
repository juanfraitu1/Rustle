#!/usr/bin/env python3
"""Pair precision / recall only, no F, no family matching (2026-10-29, user: "we don't need bipartite matching anymore, the only thing
that is useful is precision and recall", "we dont need F anymore, just precision and recall", "pseudogenes don't need their own category"). Truth = Compara paralog pairs at --level among protein-coding genes with a de
novo transcript in the run (famtruth's expressed universe). Predicted = gene pairs inside one family after mapping copies to genes
(famtruth.map_copy). A predicted pair where either gene is NOT protein-coding (pseudogene, lncRNA) is IGNORED: neither hit nor false.
Also reports 'related' (Compara at any level / same HGNC group / same symbol) among the non-true predicted pairs, so a false pair means
an unrelated protein-coding gene pair.

  pairscore.py --species human|gorilla|hg002 --gtf RUN.families.gtf [--level simiiformes] [--fasta GENOME --noncoding] COPIES_TSV...

--noncoding (2026-10-29, user: "recover all protein coding, pseudogenes and lncRNAs but use Compara for the things that are validated as
protein coding, then the others can be judged based on discarding that they are just false positives due to domain sharing"): every
predicted pair with a pseudogene / lncRNA / other gene is judged by SEQUENCE: the two genes' spliced exon sequences (annotation exons,
not our loci) are aligned with minimap2 -c -x asm20; accepted when the shorter spliced sequence is covered >= 50% at >= 85% identity
(a paralog or a pseudogene of a family member), rejected as DOMAIN-ONLY otherwise (a single shared domain exon covers far less than
half of a gene). Reported: nc_pairs, nc_accepted, nc_domain_only, nc_precision. No recall for non-coding pairs (no truth).
  env as famtruth (+ DENOVO_HG002_GFF / DENOVO_HG002_XREF for hg002, DENOVO_HGNC)
"""
import sys, os, csv, collections, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species = a[a.index('--species') + 1]; gtf = a[a.index('--gtf') + 1]
level = a[a.index('--level') + 1] if '--level' in a else 'simiiformes'
paths = [x for i, x in enumerate(a) if not x.startswith('--') and a[i - 1] not in ('--species', '--gtf', '--level', '--fasta')]
if species == 'hg002':
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
    import guided3 as g3
    g3.HG002_XREF = os.environ['DENOVO_HG002_XREF']; F.GGO_GFF = os.environ['DENOVO_HG002_GFF']
    G = g3.load_genes('hg002')
else:
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
idx = F.exon_index(G)
expressed = F.expressed_from_gtf(gtf, idx)
truth, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS[level])
truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
anyp = set()
for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
    x, y = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
    if x and y: anyp.add((min(x, y), max(x, y)))
groups = collections.defaultdict(set)
hg = os.environ.get('DENOVO_HGNC')
if hg:
    for r in csv.DictReader(open(hg), delimiter='\t'):
        for g in r['gene_group_id'].split('|'):
            if g and r['ensembl_gene_id']: groups[r['ensembl_gene_id']].add(g)
def related(x, y):
    if G[x]['name'] == G[y]['name']: return True
    ex, ey = {e.split('.')[0] for e in G[x]['ens']}, {e.split('.')[0] for e in G[y]['ens']}
    return any((min(p, q), max(p, q)) in anyp for p in ex for q in ey) or any(groups[p] & groups[q] for p in ex for q in ey)
NONCODING = '--noncoding' in a; FASTA = a[a.index('--fasta') + 1] if '--fasta' in a else None
cov_cache = {}
def noncoding_judge(pairs_nc):
    """pairs_nc: set of (gene, gene). Returns {pair: (cov_shorter, identity)} from one minimap2 run of the spliced sequences."""
    import subprocess, tempfile
    need = {g for p in pairs_nc for g in p} - {g for p in cov_cache for g in p}
    genes_needed = sorted({g for p in pairs_nc if p not in cov_cache for g in p})
    if not genes_needed: return
    seqs = {}
    by_chrom = collections.defaultdict(list)
    for g in genes_needed: by_chrom[G[g]['chrom']].append(g)
    for chrom, gs in by_chrom.items():
        raw = subprocess.run(['samtools', 'faidx', FASTA, chrom], capture_output=True, text=True).stdout
        seq = ''.join(raw.split('\n')[1:])
        for g in gs:
            ex = sorted(G[g]['exons']); s_ = ''.join(seq[s0:e] for s0, e in ex).upper()
            if G[g]['strand'] == '-': s_ = s_.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            seqs[g] = s_
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'g.fa')
        with open(fa, 'w') as w:
            for g, s_ in seqs.items():
                if s_: w.write(f">{g}\n{s_}\n")
        paf = subprocess.run(['minimap2', '-c', '-x', 'asm20', '-N', '500', '-p', '0.05', '--secondary=yes', '-t', '8', fa, fa], capture_output=True, text=True).stdout
    best = {}
    for line in paf.split('\n'):
        f = line.split('\t')
        if len(f) < 12 or f[0] == f[5]: continue
        q, t = f[0], f[5]; ql, tl = int(f[1]), int(f[6]); alen = int(f[3]) - int(f[2]); m, bl = int(f[9]), int(f[10])
        cov = alen / min(ql, tl); idn = m / bl if bl else 0
        k = (min(q, t), max(q, t))
        if k not in best or cov > best[k][0]: best[k] = (cov, idn)
    for p in pairs_nc:
        if p not in cov_cache: cov_cache[p] = best.get(p, (0.0, 0.0))
print(f"# level {level}: truth pairs among expressed protein-coding genes {len(truth)} (universe {len(universe & expressed)})")
print("run\tpred_pairs\tignored_noncoding\thit\tfalse\tof_which_related\tprecision\trecall" + ("\tnc_pairs\tnc_accepted\tnc_domain_only\tnc_precision" if NONCODING else ""))
for p in paths:
    fams, _, _ = F.load_catalog(p, G, idx)
    pred = set()
    for v in fams.values():
        for x, y in itertools.combinations(sorted(v), 2): pred.add((x, y))
    coding = {(x, y) for x, y in pred if G[x]['biotype'] == 'protein_coding' and G[y]['biotype'] == 'protein_coding' and x in universe and y in universe}
    hit = len(coding & truth); false = [(x, y) for x, y in coding if (x, y) not in truth]
    rel = sum(1 for x, y in false if related(x, y))
    prec = hit / len(coding) if coding else 0.0; rec = hit / len(truth) if truth else 0.0
    extra = ''
    if NONCODING:
        nc = {(x, y) for x, y in pred if (x, y) not in coding and (G[x]['biotype'] != 'protein_coding' or G[y]['biotype'] != 'protein_coding')}
        noncoding_judge(nc)
        acc = sum(1 for q in nc if cov_cache[q][0] >= 0.5 and cov_cache[q][1] >= 0.85)
        extra = f"\t{len(nc)}\t{acc}\t{len(nc) - acc}\t{acc / len(nc) if nc else 0:.3f}"
        if os.environ.get('PAIRSCORE_DUMP'):  # per-pair verdicts of the non-coding judge: gene, gene, name, name, biotypes, cov, identity, verdict
            with open(p + '.ncpairs.tsv', 'w') as w:
                for x, y in sorted(nc):
                    c, i = cov_cache[(x, y)]
                    w.write(f"{x}\t{y}\t{G[x]['name']}\t{G[y]['name']}\t{G[x]['biotype']}\t{G[y]['biotype']}\t{c:.2f}\t{i:.3f}\t{'accepted' if c >= 0.5 and i >= 0.85 else 'domain_only'}\n")
    print(f"{p}\t{len(coding)}\t{len(pred) - len(coding)}\t{hit}\t{len(false)}\t{rel}\t{prec:.3f}\t{rec:.3f}{extra}")
