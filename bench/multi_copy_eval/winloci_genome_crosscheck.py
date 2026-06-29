#!/usr/bin/env python3
"""Definitive portfolio filter: for each screen 'win', keep only gained RefSeq
transcripts that are BOTH (a) annotated at the CANDIDATE (secondary-dependent)
locus and (b) genuinely missed by WHOLE-GENOME StringTie (not a region artifact).
Reads:
  - the screen workflow result JSON (arg1)
  - candidates.json (arg2)
  - genome-wide StringTie gffcompare tmap (arg3)  [class '=' ref ids = ST genome-wide matches]
  - GGO_genomic.gff (arg4)  [for transcript coordinates]
"""
import json, re, sys
res_f, cand_f, gst_tmap, gff = sys.argv[1:5]

d = json.load(open(res_f)); res = d['result'] if 'result' in d else d
if isinstance(res, str): res = json.loads(res)
conf = res['confirmed']
cands = {c['gene_id']: c for c in json.load(open(cand_f))}

# genome-wide ST matched ref-tx set
st_genome = set()
for line in open(gst_tmap):
    f = line.rstrip('\n').split('\t')
    if len(f) > 2 and f[2] == '=':
        st_genome.add(f[1])

# coordinates of every gained transcript
want = set()
for w in conf:
    for g in w['gains']: want.add(g.replace('rna-', '').replace('gene-', ''))
txloc = {}
for line in open(gff):
    if line[0] == '#': continue
    f = line.split('\t')
    if len(f) < 9: continue
    m = re.search(r'ID=(?:rna-|gene-)?([^;]+)', f[8])
    if not m: continue
    tid = m.group(1)
    if tid in want and f[2] in ('mRNA', 'transcript', 'lnc_RNA', 'ncRNA'):
        txloc[tid] = (f[0], int(f[3]), int(f[4]))

rows = []
for w in conf:
    c = cands.get(w['gene']); v = w['verdict']
    cc = (c['chrom'], c['start'], c['end']) if c else None
    genuine = []   # at candidate locus AND genome-ST misses
    for g in w['gains']:
        tid = g.replace('rna-', '').replace('gene-', '')
        loc = txloc.get(tid)
        at_cand = bool(loc and cc and loc[0] == cc[0] and loc[1] <= cc[2] and loc[2] >= cc[1])
        st_has = (('rna-' + tid) in st_genome) or (tid in st_genome) or (g in st_genome)
        if at_cand and not st_has:
            genuine.append(g)
    rows.append(dict(gene=w['gene'], owner_frac=v.get('owner_frac', 0),
                     strict=v.get('n_strict_owner', 0), pure=v.get('n_pure_owner', 0),
                     reads=v.get('n_chain_reads', 0), ngains=len(w['gains']),
                     genuine=genuine, region=w['region']))

rows.sort(key=lambda r: (-len(r['genuine']), -r['owner_frac']))
ngenuine_tx = len(set(t for r in rows for t in r['genuine']))
nwin = sum(1 for r in rows if r['genuine'])
print(f"genome-wide ST matches {len(st_genome)} RefSeq transcripts")
print(f"\n=== GENUINE PORTFOLIO: {nwin} loci with >=1 gain that is at-candidate AND genome-ST-missed "
      f"({ngenuine_tx} unique transcripts) ===")
print(f"{'gene':22} {'own_frac':>8} {'pure':>5} {'reads':>6} {'genuine_tx':>10}  transcripts")
for r in rows:
    if not r['genuine']: continue
    print(f"{r['gene']:22} {r['owner_frac']:>8.3f} {r['pure']:>5} {r['reads']:>6} "
          f"{len(r['genuine']):>10}  {','.join(g.replace('rna-','') for g in r['genuine'])}")
print(f"\n=== dropped (no genuine gain after filters): "
      f"{sum(1 for r in rows if not r['genuine'])} of {len(rows)} screen-wins ===")
for r in rows:
    if not r['genuine']:
        print(f"  {r['gene']:22} own_frac={r['owner_frac']:.3f} ngains={r['ngains']} (all at-sibling or ST-has-genome-wide)")
