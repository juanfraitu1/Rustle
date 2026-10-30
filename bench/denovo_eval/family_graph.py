#!/usr/bin/env python3
"""FAMILY VARIATION GRAPH from a sketch family (2026-10-30, user: the variation graph is the thesis model; the sketch only names the
members). For each requested family: whole member loci (introns included, strand-normalised to +) -> minimap2 all-vs-all (asm20)
-> seqwish induction -> GFA with one path per locus; odgi stats (length, nodes, edges, paths, components). vg msga is gone in
vg >= 1.76 (LEDGER_GENERAL §18 used it), seqwish is its replacement here.
  family_graph.py COPIES_TSV GENOME_FA OUT_DIR --family NAME=REGEX [--family ...] [--species human] [--pad 0]
  (run inside the pixi 'graph' environment: pixi run -e graph python3 bench/denovo_eval/family_graph.py ...; env as famtruth)
A family is chosen as the catalogue family holding most protein-coding genes whose name matches REGEX; every locus of that family
becomes a path. Writes OUT_DIR/NAME.loci.fa, NAME.members.tsv, NAME.paf, NAME.gfa, NAME.og and prints the stats line."""
import sys, os, re, csv, subprocess, collections, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
COPIES, GENOME, OUT = a[:3]
fams = [x.split('=', 1) for i, x in enumerate(a) if i > 0 and a[i - 1] == '--family']
species = a[a.index('--species') + 1] if '--species' in a else 'human'
PAD = int(a[a.index('--pad') + 1]) if '--pad' in a else 0
THREADS = os.environ.get('THREADS', '8')
os.makedirs(OUT, exist_ok=True)
for tool in ('minimap2', 'seqwish', 'odgi', 'samtools'):
    if not shutil.which(tool): sys.exit(f"{tool} not on PATH (use: pixi run -e graph ...)")
G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
rows = list(csv.DictReader(open(COPIES), delimiter='\t'))
idx = F.exon_index(G)
gene_of = {}  # locus -> annotated gene (guided: the tid is the gene; de novo: map the locus by exon overlap as the scorers do)
for r in rows:
    if r['tid'] in G: gene_of[r['tid']] = r['tid']
    else:
        g = F.map_copy(idx, r['chrom'], r['strand'] if r['strand'] in '+-' else '+', F.blocks(r['exons']))
        if g: gene_of[r['tid']] = g
fam_of = collections.defaultdict(set)  # gene -> families holding a locus of it
for r in rows:
    if r['tid'] in gene_of: fam_of[gene_of[r['tid']]].add(r['family_id'])
rc = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
print("family\tfamily_id\tloci\tnamed_coding\tseq_bp\tgraph_bp\tnodes\tedges\tpaths\tcomponents")
for name, pat in fams:
    rx = re.compile(pat)
    mem = [g for g, d in G.items() if rx.match(d['name']) and d['biotype'] == 'protein_coding' and g in fam_of]
    if not mem: print(f"{name}\tno member in a family"); continue
    fid = collections.Counter(f_ for g in mem for f_ in fam_of[g]).most_common(1)[0][0]
    loci = [r for r in rows if r['family_id'] == fid]
    total = 0
    with open(f'{OUT}/{name}.loci.fa', 'w') as w, open(f'{OUT}/{name}.members.tsv', 'w') as t:
        t.write('name\tbiotype\tregion\tstrand\tlength\n')
        for r in loci:
            g = gene_of.get(r['tid']); d = G.get(g) if g else None; nm = d['name'] if d else r['tid']; bt = d['biotype'] if d else 'unannotated'
            s, e = max(0, int(r['start']) - PAD), int(r['end']) + PAD; reg = f"{r['chrom']}:{s + 1}-{e}"
            seq = ''.join(subprocess.run(['samtools', 'faidx', GENOME, reg], capture_output=True, text=True).stdout.split('\n')[1:])
            if r['strand'] == '-': seq = seq.translate(rc)[::-1]
            total += len(seq)
            w.write(f">{nm}|{reg}|{r['strand']}\n{seq}\n"); t.write(f"{nm}\t{bt}\t{reg}\t{r['strand']}\t{len(seq)}\n")
    sh = lambda c: subprocess.run(c, shell=True, check=True, capture_output=True, text=True)
    sh(f"minimap2 -cx asm20 -X -t {THREADS} {OUT}/{name}.loci.fa {OUT}/{name}.loci.fa > {OUT}/{name}.paf")
    sh(f"seqwish -s {OUT}/{name}.loci.fa -p {OUT}/{name}.paf -g {OUT}/{name}.gfa -t {THREADS} -k 19")
    sh(f"odgi build -g {OUT}/{name}.gfa -o {OUT}/{name}.og")
    st = sh(f"odgi stats -i {OUT}/{name}.og -S").stdout.strip().split('\n')[-1].split('\t')
    par = {}
    def find(x):
        par.setdefault(x, x)
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for line in open(f"{OUT}/{name}.gfa"):
        if line[0] == 'S': find(line.split('\t')[1])
        elif line[0] == 'L': p = line.split('\t'); par[find(p[1])] = find(p[3])
    ncomp = len({find(x) for x in list(par)})
    print(f"{name}\t{fid}\t{len(loci)}\t{len(mem)}\t{total}\t{st[0]}\t{st[1]}\t{st[2]}\t{st[3]}\t{ncomp}")
