#!/usr/bin/env python3
"""r1533 prototype: loci from tied reads' own placements -> seqwish graph per connected locus set -> GraphAligner threading ->
one model per locus from its own assigned reads (additive GTF).  usage: thread_tied.py BAM AS_TABLE GENOME_FA OUT_PREFIX
Panel mode: PANEL_GENES=panel.genes.tsv restricts the tied reads to those with a record at the member spans of the families that
hold a base=0 multimapper member (r1532's 28 families); the spans choose reads only, never loci.
Needs on PATH: samtools, minimap2, seqwish, GraphAligner (pixi default + graph envs)."""
import sys, os, csv, re, subprocess, collections, pysam
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F
bam_path, as_table, genome, out = sys.argv[1:5]; TIE = 0.98; THREADS = os.environ.get('THREADS', '16')
os.makedirs(out + '.d', exist_ok=True); D = out + '.d'
def say(m): print(f'[thread_tied] {m}', flush=True)
best = {}
with open(as_table) as fh:
    for l in fh:
        f = l.rstrip('\n').split('\t')
        if len(f) > 3 and f[1].isdigit() and f[2].isdigit() and int(f[2]) >= TIE * int(f[1]): best[f[0]] = int(f[1])  # tied reads only
say(f'tied reads in the as-table: {len(best)}')
bam = pysam.AlignmentFile(bam_path)
if os.environ.get('PANEL_GENES'):
    G = F.load_cat(); fams = collections.defaultdict(list); want = set()
    for r in csv.DictReader(open(os.environ['PANEL_GENES']), delimiter='\t'):
        if int(r['tied_molecules']) < 2: continue
        fams[r['family']].append(r['gene'])
        if r['base'] == '0' and (int(r['primary_molecules']) < 2 or float(r['mapq0_share']) >= 0.5): want.add(r['family'])
    names = set()
    for f in want:
        for g in fams[f]:
            v = G[g]
            for rec in bam.fetch(v['chrom'], v['exons'][0][0], v['exons'][-1][1]):
                if rec.query_name in best and not rec.is_supplementary: names.add(rec.query_name)
    say(f'panel: {len(want)} families, {len(names)} tied reads with a record at a member span')
else:
    names = set(best)
with open(f'{D}/names.txt', 'w') as w: w.write('\n'.join(sorted(names)) + '\n')
subprocess.run(f'samtools view -@ {THREADS} -b -N {D}/names.txt -o {D}/tied.bam {bam_path} && samtools index {D}/tied.bam', shell=True, check=True)
# records at tied placements (AS >= TIE x best), blocks from the CIGAR (N = intron; D kept inside a block)
recs = collections.defaultdict(list); seqs = {}
for rec in pysam.AlignmentFile(f'{D}/tied.bam'):
    if rec.is_unmapped or rec.is_supplementary or rec.query_name not in best: continue
    if rec.query_sequence: seqs[rec.query_name] = rec.query_sequence
    if not rec.has_tag('AS') or rec.get_tag('AS') < TIE * best[rec.query_name]: continue
    blocks = []; pos = rec.reference_start
    for op, n in rec.cigartuples:
        if op in (0, 7, 8, 2):
            if blocks and blocks[-1][1] == pos: blocks[-1][1] = pos + n
            else: blocks.append([pos, pos + n])
            pos += n
        elif op == 3: pos += n
    strand = '-' if rec.is_reverse else '+'
    recs[(rec.reference_name, strand)].append((blocks[0][0], blocks[-1][1], rec.query_name, [tuple(b) for b in blocks]))
say(f'tied placements: {sum(len(v) for v in recs.values())} records of {len(seqs)} reads with sequence')
# loci = overlap clusters of placements per (chrom, strand)
loci = []  # (chrom, strand, start, end, exons[(s,e)], reads{name: blocks})
for (c, s), rs in recs.items():
    rs.sort()
    cur = None
    for st, en, name, blocks in rs:
        if cur and st < cur[3]: cur[3] = max(cur[3], en); cur[5][name] = blocks
        else:
            cur = [c, s, st, en, None, {name: blocks}]; loci.append(cur)
    for L in loci:
        if L[4] is None:
            bl = sorted(b for bs in L[5].values() for b in bs); ex = []
            for a, b in bl:
                if ex and a <= ex[-1][1]: ex[-1][1] = max(ex[-1][1], b)
                else: ex.append([a, b])
            L[4] = [tuple(x) for x in ex]
say(f'loci: {len(loci)}')
# connected locus sets by shared reads
parent = list(range(len(loci)))
def find(i):
    while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
    return i
seen = {}
for i, L in enumerate(loci):
    for n in L[5]:
        if n in seen: parent[find(i)] = find(seen[n])
        else: seen[n] = i
groups = collections.defaultdict(list)
for i in range(len(loci)): groups[find(i)].append(i)
fa = pysam.FastaFile(genome); RC = str.maketrans('ACGTacgt', 'TGCAtgca')
def path_seq(L):
    s = ''.join(fa.fetch(L[0], a, b) for a, b in L[4]).upper()
    return s.translate(RC)[::-1] if L[1] == '-' else s
models = []; stats = collections.Counter()
for gi, members in groups.items():
    if len(members) < 2: stats['single-locus group'] += 1; continue
    gd = f'{D}/g{gi}'; os.makedirs(gd, exist_ok=True)
    with open(f'{gd}/loci.fa', 'w') as w:
        for i in members: w.write(f'>L{i}\n{path_seq(loci[i])}\n')
    rn = {n for i in members for n in loci[i][5] if n in seqs}
    with open(f'{gd}/reads.fa', 'w') as w:
        for n in rn: w.write(f'>{n}\n{seqs[n]}\n')
    r = subprocess.run(f'minimap2 -cx asm20 -X -t 4 {gd}/loci.fa {gd}/loci.fa > {gd}/loci.paf 2>/dev/null && seqwish -s {gd}/loci.fa -p {gd}/loci.paf -g {gd}/g.gfa -t 4 -k 19 2>/dev/null && GraphAligner -g {gd}/g.gfa -f {gd}/reads.fa -a {gd}/aln.gaf -x vg -t 4 > {gd}/ga.log 2>&1', shell=True)
    if r.returncode or not os.path.exists(f'{gd}/aln.gaf'): stats['graph failed'] += 1; continue
    paths = {}
    for l in open(f'{gd}/g.gfa'):
        if l[0] == 'P': p = l.split('\t'); paths[int(p[1][1:])] = {int(x[:-1]) for x in p[2].split(',')}
    cnt = collections.Counter(n for s in paths.values() for n in s); priv = {i: {n for n in paths.get(i, ()) if cnt[n] == 1} for i in members}
    bestaln = {}
    for l in open(f'{gd}/aln.gaf'):
        c = l.rstrip('\n').split('\t'); m = int(c[9])
        if c[0] not in bestaln or m > bestaln[c[0]][0]: bestaln[c[0]] = (m, {int(x) for x in re.findall(r'[<>](\d+)', c[5])})
    assigned = collections.defaultdict(list)
    for n, (m, nodes) in bestaln.items():
        hits = [i for i in members if nodes & priv[i]]
        if len(hits) == 1 and n in loci[hits[0]][5]: assigned[hits[0]].append(n)   # one locus, and the read has a placement there
        else: stats['read unassigned'] += 1
    for i, ns in assigned.items():
        if len(ns) < 2: stats['locus < 2 reads'] += 1; continue
        L = loci[i]; bl = sorted(b for n in ns for b in L[5][n]); ex = []
        for a, b in bl:
            if ex and a <= ex[-1][1]: ex[-1][1] = max(ex[-1][1], b)
            else: ex.append([a, b])
        models.append((L[0], L[1], ex, len(ns), f'TT_{L[0]}_{L[2]}')); stats['models'] += 1
with open(out + '.gtf', 'w') as w:
    for c, s, ex, n, tid in models:
        attr = f'gene_id "{tid}"; transcript_id "{tid}"; reads "{n}";'
        w.write(f'{c}\tthread_tied\ttranscript\t{ex[0][0]+1}\t{ex[-1][1]}\t.\t{s}\t.\t{attr}\n')
        for a, b in ex: w.write(f'{c}\tthread_tied\texon\t{a+1}\t{b}\t.\t{s}\t.\t{attr}\n')
say(f'groups {len(groups)}; ' + ', '.join(f'{k} {v}' for k, v in sorted(stats.items())) + f'; wrote {out}.gtf')
