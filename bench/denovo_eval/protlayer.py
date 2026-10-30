#!/usr/bin/env python3
"""PROTEIN LAYER over sketch families (2026-10-30, user: "is there a way we can add an extra layer for the members that only work at
the protein level?"). Layer 1 = the sketch families (nucleotide; O2/O3 run there). Layer 2 = a graph whose nodes are layer-1
families and whose edges come from TRANSLATED alignment (mmseqs easy-search --search-type 2, TBLASTX-like) between the spliced
exon sequences of their protein-coding members: a hit counts when amino-acid identity >= --min-id and it covers >= --cov of the
shorter sequence. Two families are joined when the hits touch >= --min-frac of the protein-coding members of the SMALLER family
(the SD-merge rule of sd_merge.py, so one pseudogene cannot glue two 30-member families). Connected components of layer 2 are the
reconciled (super)families.
  protlayer.py COPIES_TSV GENOME_FA OUT [--min-id 0.4] [--cov 0.5] [--evalue 1e-5] [--min-frac 0.5[,0.3,...]] [--threads 16] [--species human]
Writes OUT.pc.fa (spliced exons of protein-coding members), OUT.m8 (reused if present), OUT.prot_edges.tsv (family pairs with hit
count, best identity, touched fraction of each side, touched gene names) and OUT_f{frac}.fam.copies.tsv per --min-frac value."""
import sys, os, csv, collections, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
opt = lambda k, d, f=float: f(a[a.index(k) + 1]) if k in a else d
MINID = opt('--min-id', 0.4); COV = opt('--cov', 0.5); EV = opt('--evalue', 1e-5); T = opt('--threads', 16, int)
FRACS = [float(x) for x in opt('--min-frac', '0.5', str).split(',')]
species = opt('--species', 'human', str)
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or not a[i - 1].startswith('--'))]
COPIES, GENOME, OUT = pos[:3]
G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
rows = list(csv.DictReader(open(COPIES), delimiter='\t'))
fam = {r['tid']: r['family_id'] for r in rows}
pc = [t for t in fam if t in G and G[t].get('biotype') == 'protein_coding' and G[t]['exons']]
members = collections.defaultdict(set)
for t in pc: members[fam[t]].add(t)
rc = str.maketrans('ACGT', 'TGCA')
fa = OUT + '.pc.fa'
if not os.path.exists(OUT + '.m8'):
    by = collections.defaultdict(list)
    for t in pc: by[G[t]['chrom']].append(t)
    with open(fa, 'w') as w:
        for chrom, ts in by.items():
            seq = ''.join(subprocess.run(['samtools', 'faidx', GENOME, chrom], capture_output=True, text=True).stdout.split('\n')[1:])
            for t in ts:
                s = ''.join(seq[s0:e] for s0, e in sorted(G[t]['exons'])).upper()
                if G[t]['strand'] == '-': s = s.translate(rc)[::-1]
                if len(s) >= 90: w.write(f">{t}\n{s}\n")
    print(f"{len(pc)} protein-coding members in {len(members)} families; translated all-vs-all", file=sys.stderr)
    subprocess.run(['mmseqs', 'easy-search', fa, fa, OUT + '.m8', OUT + '.tmp', '--search-type', '2', '-e', str(EV), '--threads', str(T),
                    '--max-seqs', '2000', '--format-output', 'query,target,evalue,pident,qcov,tcov,alnlen'], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
hits = {}  # (famA, famB) -> {(q, t): best pident}
for line in open(OUT + '.m8'):
    q, t, ev, pid, qc, tc, al = line.rstrip('\n').split('\t')
    if q == t or q not in fam or t not in fam or fam[q] == fam[t]: continue
    pid, qc, tc = float(pid), float(qc), float(tc)
    if pid > 1: pid /= 100
    if pid < MINID or max(qc, tc) < COV: continue
    k = (min(fam[q], fam[t]), max(fam[q], fam[t]))
    e = (min(q, t), max(q, t))
    hits.setdefault(k, {})
    if e not in hits[k] or pid > hits[k][e]: hits[k][e] = pid
nm = lambda t: G[t]['name']
edges = []
with open(OUT + '.prot_edges.tsv', 'w') as w:
    w.write('famA\tfamB\tpc_A\tpc_B\thits\tbest_id\tfrac_A\tfrac_B\ttouched_A\ttouched_B\n')
    for (fa_, fb_), d in sorted(hits.items()):
        ta = {x for e in d for x in e if fam[x] == fa_}; tb = {x for e in d for x in e if fam[x] == fb_}
        fra, frb = len(ta) / len(members[fa_]), len(tb) / len(members[fb_])
        edges.append((fa_, fb_, min(fra, frb)))
        w.write(f"{fa_}\t{fb_}\t{len(members[fa_])}\t{len(members[fb_])}\t{len(d)}\t{max(d.values()):.2f}\t{fra:.2f}\t{frb:.2f}\t"
                f"{','.join(sorted(nm(x) for x in ta))}\t{','.join(sorted(nm(x) for x in tb))}\n")
print(f"family pairs with a protein hit (id >= {MINID}, cov >= {COV}): {len(edges)}", file=sys.stderr)
for fr in FRACS:
    gp = {f_: f_ for f_ in set(fam.values())}
    def find(x):
        while gp[x] != x: gp[x] = gp[gp[x]]; x = gp[x]
        return x
    n = 0
    for fa_, fb_, f_ in edges:
        if f_ >= fr and find(fa_) != find(fb_): gp[find(fa_)] = find(fb_); n += 1
    with open(f"{OUT}_f{fr}.fam.copies.tsv", 'w') as w:
        w.write('\t'.join(rows[0].keys()) + '\n')
        for r in rows:
            r = dict(r); r['family_id'] = find(r['family_id']); w.write('\t'.join(r.values()) + '\n')
    print(f"min-frac {fr}: {n} merges, {len(set(find(f_) for f_ in fam.values()))} families", file=sys.stderr)
