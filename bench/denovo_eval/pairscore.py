#!/usr/bin/env python3
"""Pair precision / recall only, no F, no family matching (2026-10-29, user: "we don't need bipartite matching anymore, the only thing
that is useful is precision and recall", "we dont need F anymore, just precision and recall", "pseudogenes don't need their own category"). Truth = Compara paralog pairs at --level among protein-coding genes with a de
novo transcript in the run (famtruth's expressed universe). Predicted = gene pairs inside one family after mapping copies to genes
(famtruth.map_copy). A predicted pair where either gene is NOT protein-coding (pseudogene, lncRNA) is IGNORED: neither hit nor false.
Also reports 'related' (Compara at any level / same HGNC group / same symbol) among the non-true predicted pairs, so a false pair means
an unrelated protein-coding gene pair.

  pairscore.py [--species human|gorilla|hg002|NAME --annotation GTF/GFF3 --xref TSV --ids symbol|own] [--paralogs TSV] --gtf RUN.families.gtf
               [--level simiiformes|TAXON,...] [--fasta GENOME --noncoding] [--universe-copies] COPIES_TSV...
--universe-copies (semi-guided, 2026-10-30): truth pairs among genes carrying any copy of THAT run (silent copies included) instead of
the GTF's expressed genes, so a silent copy of an unexpressed paralog can be a hit rather than a false pair.

--noncoding (2026-10-29, user: "recover all protein coding, pseudogenes and lncRNAs but use Compara for the things that are validated as
protein coding, then the others can be judged based on discarding that they are just false positives due to domain sharing"): every
predicted pair with a pseudogene / lncRNA / other gene is judged by SEQUENCE: the two genes' spliced exon sequences (annotation exons,
not our loci) are aligned with minimap2 -c -x asm20; accepted when the shorter spliced sequence is covered >= 50% at >= 85% identity
(a paralog or a pseudogene of a family member), rejected as DOMAIN-ONLY otherwise (a single shared domain exon covers far less than
half of a gene). Reported: nc_pairs, nc_accepted, nc_domain_only, nc_precision. No recall for non-coding pairs (no truth).
  env as famtruth (+ DENOVO_HG002_GFF / DENOVO_HG002_XREF for hg002, DENOVO_HGNC optional); genes: famtruth.load_genes
"""
import sys, os, csv, collections, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species = a[a.index('--species') + 1] if '--species' in a else 'human'; gtf = a[a.index('--gtf') + 1]
level = a[a.index('--level') + 1] if '--level' in a else 'simiiformes'
paths = [x for i, x in enumerate(a) if not x.startswith('--') and a[i - 1] not in F.SPECIES_FLAGS + ('--gtf', '--level', '--fasta', '--paralogs')]
# --universe-copies: see below
G = F.genes_from_argv(a, species)
idx = F.exon_index(G)
expressed = F.expressed_from_gtf(gtf, idx)
UNIVERSE_COPIES = '--universe-copies' in a  # semi-guided (2026-10-30): the expressed universe = genes carrying ANY copy of the run (silent copies included), per run
PARALOGS = F.need(a[a.index('--paralogs') + 1] if '--paralogs' in a else F.PARALOGS, 'DENOVO_PARALOGS or --paralogs')
truth, universe = F.truth_pairs(G, PARALOGS, F.levels(level))
truth_all = truth
truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
anyp = set()
for r in csv.reader(F._open(PARALOGS), delimiter='\t'):  # columns: gene, paralogue, type, last common ancestor (any species)
    if len(r) < 2 or r[0].startswith('Gene'): continue
    x, y = r[0].split('.')[0], r[1].split('.')[0]
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
print("run\tpred_pairs\tignored_noncoding\thit\tfalse\tof_which_related\tprecision\trecall" + ("\tignored_same_symbol" if os.environ.get('PAIRSCORE_SAME_SYMBOL') == 'ignore' else "") + ("\tnc_pairs\tnc_accepted\tnc_domain_only\tnc_precision" if NONCODING else ""))
for p in paths:
    fams, _, _ = F.load_catalog(p, G, idx)
    if UNIVERSE_COPIES:
        expressed = {g for v in fams.values() for g in v}
        truth = {tuple(sorted(q)) for q in truth_all if q[0] in expressed and q[1] in expressed}
    pred = set()
    for v in fams.values():
        for x, y in itertools.combinations(sorted(v), 2): pred.add((x, y))
    coding = {(x, y) for x, y in pred if G[x]['biotype'] == 'protein_coding' and G[y]['biotype'] == 'protein_coding' and x in universe and y in universe}
    same = {(x, y) for x, y in coding if G[x]['name'] == G[y]['name']} if os.environ.get('PAIRSCORE_SAME_SYMBOL') == 'ignore' else set()
    coding -= same  # r1530: same-symbol pairs are CAT duplicate ids of one Ensembl gene (213/213 share a source id); Compara cannot name them
    hit = len(coding & truth); false = [(x, y) for x, y in coding if (x, y) not in truth]
    if os.environ.get('PAIRSCORE_FAMDUMP'):  # per-family readout of the CODING truth (2026-10-30, user "how does it do in all the gene families"):
        import re as _re                      # family key = gene symbol with trailing digits/letters stripped; hit / missed truth pairs; false pairs
        key = lambda g: _re.sub(r'[0-9]+[A-Z]?$|[A-Z]$', '', G[g]['name'].split('-')[0]) or G[g]['name']
        fam = collections.defaultdict(lambda: [0, 0, 0])
        for x, y in truth:
            k = key(x) if key(x) == key(y) else f"{key(x)}/{key(y)}"; fam[k][0 if (x, y) in coding else 1] += 1
        for x, y in false:
            k = key(x) if key(x) == key(y) else f"{key(x)}/{key(y)}"; fam[k][2] += 1
        with open(os.environ['PAIRSCORE_FAMDUMP'] + '.' + os.path.basename(p) + '.famdump.tsv', 'w') as fh:
            fh.write('family\thit\tmissed\tfalse\n')
            for k, v in sorted(fam.items(), key=lambda kv: -(kv[1][0] + kv[1][1] + kv[1][2])): fh.write(f"{k}\t{v[0]}\t{v[1]}\t{v[2]}\n")
    if os.environ.get('PAIRSCORE_FALSEDUMP'):  # the false coding pairs of each run (gene ids + names), r1529 step 0
        with open(os.environ['PAIRSCORE_FALSEDUMP'] + '.' + os.path.basename(p) + '.false.tsv', 'w') as fh:
            for x, y in sorted(false): fh.write(f"{x}\t{y}\t{G[x]['name']}\t{G[y]['name']}\n")
    if os.environ.get('PAIRSCORE_HITDUMP'):  # the hit coding truth pairs of each run (gene ids + names), to diff two runs' recall
        with open(os.environ['PAIRSCORE_HITDUMP'] + '.' + os.path.basename(p) + '.hits.tsv', 'w') as fh:
            for x, y in sorted(coding & truth): fh.write(f"{x}\t{y}\t{G[x]['name']}\t{G[y]['name']}\n")
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
    print(f"{p}\t{len(coding)}\t{len(pred) - len(coding) - len(same)}\t{hit}\t{len(false)}\t{rel}\t{prec:.3f}\t{rec:.3f}{extra}" + (f"\t{len(same)}" if same else ""))
