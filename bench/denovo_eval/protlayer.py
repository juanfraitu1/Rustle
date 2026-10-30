#!/usr/bin/env python3
"""PROTEIN RECONCILER over sketch families (2026-10-30, user: "is there a way we can add an extra layer for the members that only
work at the protein level?" -> "can we try to be more specific" -> "lets implement the second rule"). Layer 1 = the sketch families
(nucleotide; O2/O3 run there). Layer 2 = a graph whose nodes are layer-1 families and whose edges come from TRANSLATED alignment
(mmseqs easy-search --search-type 2, TBLASTX-like) between the spliced exon sequences of their members: a hit counts when amino-acid
identity >= --min-id and it covers >= --cov of the shorter sequence. A family pair is an edge when (default gates of the second rule)
  --rbh         each family is the other's best-hit family (by best member identity), and
  --same-chrom  every touched member on both sides lies on one chromosome (SD families are local; superfamilies are not),
and two families are joined when the hits touch >= --min-frac of the members of the SMALLER family (the SD-merge rule of
sd_merge.py, so one pseudogene cannot glue two 30-member families). Connected components of layer 2 are the reconciled families.
  protlayer.py COPIES_TSV GENOME_FA OUT [--min-id 0.7] [--cov 0.5] [--evalue 1e-5] [--min-frac 0.5[,0.3,...]] [--threads 16]
               [--rbh] [--same-chrom] [--members all|coding --species human|gorilla] [--tag T]
--members all (default, ANNOTATION-FREE: the de novo form) takes every locus of the copies file, its spliced exons from the exons column;
--members coding (the guided form used for the A119b grid of 2026-10-30) keeps only annotated protein-coding members and needs the
catalogue (env as famtruth). --tag names the output OUT_<tag>_f{frac}.fam.copies.tsv.
Writes OUT.pc.fa (spliced exons), OUT.m8 (reused if present), OUT.prot_edges.tsv (family pairs with hit count, best identity, touched
fraction of each side, touched members, gate verdict) and OUT_f{frac}.fam.copies.tsv per --min-frac value."""
import sys, os, csv, collections, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
opt = lambda k, d, f=float: f(a[a.index(k) + 1]) if k in a else d
MINID = opt('--min-id', 0.7); COV = opt('--cov', 0.5); EV = opt('--evalue', 1e-5); T = opt('--threads', 16, int)
FRACS = [float(x) for x in opt('--min-frac', '0.5', str).split(',')]
species = opt('--species', 'human', str); RBH = '--rbh' in a; SAMECHR = '--same-chrom' in a; TAG = opt('--tag', '', str)
MEMBERS = opt('--members', 'all', str)
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or not a[i - 1].startswith('--'))]
COPIES, GENOME, OUT = pos[:3]
rows = list(csv.DictReader(open(COPIES), delimiter='\t'))
fam = {r['tid']: r['family_id'] for r in rows}
if MEMBERS == 'coding':
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
    pc = [t for t in fam if t in G and G[t].get('biotype') == 'protein_coding' and G[t]['exons']]
    seqs = {t: (G[t]['chrom'], G[t]['strand'], sorted(G[t]['exons'])) for t in pc}
    nm = lambda t: G[t]['name']
else:
    seqs = {r['tid']: (r['chrom'], r['strand'] if r['strand'] in '+-' else '+', sorted(F.blocks(r['exons']))) for r in rows}
    pc = list(seqs); nm = lambda t: t
chrom_of = {t: v[0] for t, v in seqs.items()}
members = collections.defaultdict(set)
for t in pc: members[fam[t]].add(t)
rc = str.maketrans('ACGT', 'TGCA')
fa = OUT + '.pc.fa'
if not os.path.exists(OUT + '.m8'):
    by = collections.defaultdict(list)
    for t in pc: by[seqs[t][0]].append(t)
    with open(fa, 'w') as w:
        for chrom, ts in by.items():
            seq = ''.join(subprocess.run(['samtools', 'faidx', GENOME, chrom], capture_output=True, text=True).stdout.split('\n')[1:])
            for t in ts:
                s = ''.join(seq[s0:e] for s0, e in seqs[t][2]).upper()
                if seqs[t][1] == '-': s = s.translate(rc)[::-1]
                if len(s) >= 90: w.write(f">{t}\n{s}\n")
    print(f"{len(pc)} members ({MEMBERS}) in {len(members)} families; translated all-vs-all", file=sys.stderr)
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
bestfam = {}  # family -> (best identity, partner family)
for (fa_, fb_), d in hits.items():
    m = max(d.values())
    for x, y in ((fa_, fb_), (fb_, fa_)):
        if x not in bestfam or m > bestfam[x][0]: bestfam[x] = (m, y)
edges = []
with open(OUT + '.prot_edges.tsv', 'w') as w:
    w.write('famA\tfamB\tn_A\tn_B\thits\tbest_id\tfrac_A\tfrac_B\ttouched_A\ttouched_B\tgate\n')
    for (fa_, fb_), d in sorted(hits.items()):
        ta = {x for e in d for x in e if fam[x] == fa_}; tb = {x for e in d for x in e if fam[x] == fb_}
        fra, frb = len(ta) / len(members[fa_]), len(tb) / len(members[fb_])
        ok = True
        if RBH and not (bestfam[fa_][1] == fb_ and bestfam[fb_][1] == fa_): ok = False
        if SAMECHR and len({chrom_of[x] for x in ta | tb}) > 1: ok = False
        if ok: edges.append((fa_, fb_, min(fra, frb)))
        w.write(f"{fa_}\t{fb_}\t{len(members[fa_])}\t{len(members[fb_])}\t{len(d)}\t{max(d.values()):.2f}\t{fra:.2f}\t{frb:.2f}\t"
                f"{','.join(sorted(nm(x) for x in ta))}\t{','.join(sorted(nm(x) for x in tb))}\t{'pass' if ok else 'fail'}\n")
print(f"family pairs with a protein hit (id >= {MINID}, cov >= {COV}): {len(hits)}, passing gates: {len(edges)}", file=sys.stderr)
for fr in FRACS:
    gp = {f_: f_ for f_ in set(fam.values())}
    def find(x):
        while gp[x] != x: gp[x] = gp[gp[x]]; x = gp[x]
        return x
    n = 0
    for fa_, fb_, f_ in edges:
        if f_ >= fr and find(fa_) != find(fb_): gp[find(fa_)] = find(fb_); n += 1
    with open(f"{OUT}_{TAG + '_' if TAG else ''}f{fr}.fam.copies.tsv", 'w') as w:
        w.write('\t'.join(rows[0].keys()) + '\n')
        for r in rows:
            r = dict(r); r['family_id'] = find(r['family_id']); w.write('\t'.join(r.values()) + '\n')
    print(f"min-frac {fr}: {n} merges, {len(set(find(f_) for f_ in fam.values()))} families", file=sys.stderr)
