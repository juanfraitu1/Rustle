import sys, os, csv, collections, itertools, subprocess, statistics
sys.path.insert(0, os.path.expanduser('~/Rustle/bench/denovo_eval'))
import famtruth as F
p, d, arm = sys.argv[1:4]; K = 17
G = F.load_ggo(F.load_cat()); idx = F.exon_index(G)
expressed = F.expressed_from_gtf(os.path.expanduser(f'~/runs/{d}/{p}.families.gtf'), idx)
truth, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS['simiiformes'])
anyp = set()
for r in csv.DictReader(open(F.PARALOGS), delimiter='\t'):
    x, y = r['Gene stable ID'].split('.')[0], r['Human paralogue gene stable ID'].split('.')[0]
    if x and y: anyp.add((min(x, y), max(x, y)))
def related(x, y):
    if G[x]['name'] == G[y]['name']: return True
    ex, ey = {e.split('.')[0] for e in G[x]['ens']}, {e.split('.')[0] for e in G[y]['ens']}
    return any((min(a, b), max(a, b)) in anyp for a in ex for b in ey)
# families + loci from the arm copies
fam = collections.defaultdict(list)
for r in csv.DictReader(open(arm), delimiter='\t'):
    g = F.map_copy(idx, r['chrom'], r['strand'], F.blocks(r['exons']))
    if g and G[g]['biotype'] == 'protein_coding' and g in universe:
        fam[r['family_id']].append((g, r['chrom'], r['strand'], F.blocks(r['exons'])))
def kmers(seq):
    return {seq[i:i+K] for i in range(len(seq) - K + 1) if seq[i:i+K].isupper() and 'N' not in seq[i:i+K]}
def fetch(chrom, s, e):
    return subprocess.run(['samtools', 'faidx', os.environ['G'], f'{chrom}:{s+1}-{e}'], capture_output=True, text=True).stdout.split('\n', 1)[1].replace('\n', '')
cache = {}
def intron_kmers(chrom, strand, ex):
    key = (chrom, tuple(ex))
    if key in cache: return cache[key]
    ex = sorted(ex); ks = set(); bp = 0
    for (s0, e0), (s1, e1) in zip(ex, ex[1:]):
        if s1 - e0 < K: continue
        seq = fetch(chrom, e0, s1); bp += len(seq); ks |= kmers(seq)
        if strand == '-': ks |= kmers(seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1])
    cache[key] = (ks, bp); return cache[key]
rows = collections.defaultdict(list)
for fid, mem in fam.items():
    for (gx, cx, sx, ex), (gy, cy, sy, ey) in itertools.combinations(mem, 2):
        if gx == gy: continue
        pr = (min(gx, gy), max(gx, gy))
        cls = 'hit' if pr in truth else ('related' if related(gx, gy) else 'unrelated')
        kx, bx = intron_kmers(cx, sx, ex); ky, by = intron_kmers(cy, sy, ey)
        if min(len(kx), len(ky)) < 50: rows[cls].append(None); continue
        rows[cls].append(len(kx & ky) / min(len(kx), len(ky)))
for cls in ('hit', 'related', 'unrelated'):
    v = rows[cls]; na = sum(x is None for x in v); vv = sorted(x for x in v if x is not None)
    if not vv: print(cls, 'n', len(v), 'no-intron', na); continue
    q = lambda f: vv[min(len(vv)-1, int(f*len(vv)))]
    print(f"{p} {cls:9s} pairs {len(v):4d} no-intron/short {na:3d} | intron containment q10 {q(.1):.3f} q25 {q(.25):.3f} median {q(.5):.3f} q75 {q(.75):.3f} | >=0.05: {sum(x>=0.05 for x in vv)}/{len(vv)}  >=0.2: {sum(x>=0.2 for x in vv)}/{len(vv)}")
