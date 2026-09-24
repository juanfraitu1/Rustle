#!/usr/bin/env python3
"""O3 simulations with truth (docs/PREREG_o3_reference_bias_2026-09-23.md arm A, docs/PREREG_o3_rna_only_2026-09-23.md
control + addendum 2): reads from K annotated chr genes plus an EXTRA COPY that is absent from the reference, mapped to
the unmodified reference with the shipped minimap2 settings.

modes
  transcript  each transcript of the gene mutated independently at divergence d (arm A: where do the copy's reads go?)
  genomic     the gene SPAN mutated once, every transcript read off it (one real genomic copy: the o3_rna_flag positive
              control); also writes <out>.copies.fa/.mmi, the mutated spans, usable as an o3_rna_flag --confirm genome
  shuffled    the extra copy's transcript has exons 2 and 3 swapped (+ divergence d): how minimap2 represents an
              exon-order rearrangement (insertion vs clip vs supplementary), the structural detector's control

usage: o3_sim_copies.py MODE REF.gtf REF.fa CHROM K DIVERGENCE OUT_PREFIX SEED
Read names carry the truth: `template|gene|tx|i` / `extra|gene|tx|i` (`shuffled|gene|i` in shuffled mode)."""
import sys, re, random, collections, subprocess, statistics
import pysam
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from sim_reads import simulate_reads

MODE, gtf, fa_p, CHROM, K, div, out, seed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]), sys.argv[7], int(sys.argv[8])
assert MODE in ('transcript', 'genomic', 'shuffled')
fa = pysam.FastaFile(fa_p); rng = random.Random(seed)
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
MM2 = "minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4"

ex = collections.defaultdict(list); gene = {}; st = {}
for ln in open(gtf):
    if ln[0] == '#': continue
    r = ln.rstrip('\n').split('\t')
    if len(r) < 9 or r[0] != CHROM or r[2] != 'exon': continue
    tm = re.search(r'transcript_id "([^"]+)"', r[8]); gm = re.search(r'gene_id "([^"]+)"', r[8]) or re.search(r'gene_name "([^"]+)"', r[8])
    if not tm or not gm: continue
    ex[tm.group(1)].append((int(r[3]), int(r[4]))); gene[tm.group(1)] = gm.group(1); st[tm.group(1)] = r[6]
min_exons = 5 if MODE == 'shuffled' else 3
by_gene = collections.defaultdict(list)
for t, e in ex.items():
    if len(e) >= min_exons: by_gene[gene[t]].append(t)
genes = sorted(g for g, ts in by_gene.items() if (len(ts) == 1 if MODE == 'shuffled' else 1 <= len(ts) <= 6))
chosen = rng.sample(genes, K)

def mutate(seq, d, r):
    s = list(seq); n = int(round(d * len(s)))
    for i in r.sample(range(len(s)), n): s[i] = r.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)
def tx_exons(t):
    e = sorted(ex[t]); exs = [fa.fetch(CHROM, a - 1, b).upper() for a, b in e]
    if st[t] == '-': exs = [x.translate(COMP)[::-1] for x in exs][::-1]
    return e, exs

truth = {}; n = 0
copies_fa = open(out + '.copies.fa', 'w') if MODE == 'genomic' else None
with open(out + '.fq', 'w') as fq:
    for g in chosen:
        grng = random.Random(seed * 7919 + hash(g) % 100000); ts = by_gene[g]
        span = None
        if MODE == 'genomic':
            gs = min(a for t in ts for a, b in ex[t]); ge = max(b for t in ts for a, b in ex[t])
            span = mutate(fa.fetch(CHROM, gs - 1, ge).upper(), div, grng)
            copies_fa.write(f'>{g}_copy {CHROM}:{gs}-{ge} d={div}\n{span}\n')
        for t in (ts[:1] if MODE == 'shuffled' else ts):
            e, exs = tx_exons(t); body = ''.join(exs)
            if MODE == 'transcript': extra = mutate(body, div, grng)
            elif MODE == 'genomic':
                cp = [span[a - gs:b - gs + 1] for a, b in e]
                if st[t] == '-': cp = [x.translate(COMP)[::-1] for x in cp][::-1]
                extra = ''.join(cp)
            else:
                sh = exs[:]; sh[1], sh[2] = sh[2], sh[1]; extra = mutate(''.join(sh), div, random.Random(seed + hash(g) % 100000))
            kinds = (('template', body), ('shuffled' if MODE == 'shuffled' else 'extra', extra))
            for kind, seq in kinds:
                for i, (rd, q) in enumerate(simulate_reads(seq, 10, err=0.001, indel=0.0003, seed=seed * 31 + (1 if kind != 'template' else 0) + hash(t) % 100000)):
                    j5, j3 = rng.randint(0, 30), rng.randint(0, 30); rd = rd[j5:len(rd) - j3]; q = q[j5:len(q) - j3]
                    rid = f'{kind}|{g}|{t}|{i}' if MODE != 'shuffled' else f'{kind}|{g}|{i}'
                    fq.write(f'@{rid}\n{rd}\n+\n{q}\n'); truth[rid] = (kind, g, t, e[0][0], e[-1][1], len(e)); n += 1
if copies_fa: copies_fa.close()
with open(out + '.genes.txt', 'w') as fh:
    for g in chosen: fh.write(g + '\n')
subprocess.run(f"{MM2} {fa_p} {out}.fq 2>/dev/null | samtools sort -@2 -o {out}.bam - && samtools index {out}.bam", shell=True, check=True)
if MODE == 'genomic':
    subprocess.run(f"minimap2 -x splice:hq -d {out}.copies.mmi {out}.copies.fa 2>/dev/null", shell=True, check=True)
print(f'[o3_sim_copies {MODE} d={div}] {len(chosen)} genes, {n} reads -> {out}.bam', flush=True)

# --- classify the reads' placements
recs = collections.defaultdict(list)
for r in pysam.AlignmentFile(out + '.bam'):
    if r.is_supplementary and MODE != 'shuffled': continue
    recs[r.query_name].append(r)
S = collections.defaultdict(collections.Counter); de = collections.defaultdict(list); nint = collections.defaultdict(list); ins = collections.defaultdict(list)
for rid, (kind, g, t, gs, ge, ne) in truth.items():
    rs = recs.get(rid, []); prim = [r for r in rs if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
    if not prim: S[kind]['unmapped'] += 1; continue
    p = prim[0]; at = p.reference_start + 1 <= ge and p.reference_end >= gs
    S[kind]['absorbed at template locus' if at else 'primary elsewhere'] += 1
    secs = [r for r in rs if r.is_secondary]; pas = p.get_tag('AS') if p.has_tag('AS') else 0
    if any((r.get_tag('AS') if r.has_tag('AS') else 0) >= 0.98 * pas for r in secs) or p.mapping_quality == 0: S[kind]['AS-tied / MAPQ 0'] += 1
    if p.has_tag('de'): de[kind].append(p.get_tag('de'))
    if MODE == 'shuffled':
        S[kind]['with supplementary'] += any(r.is_supplementary for r in rs)
        nint[kind].append(sum(1 for op, l in p.cigartuples if op == 3)); ins[kind].append(max([l for op, l in p.cigartuples if op == 1] or [0]))
        S[kind]['insertion >= 50 bp'] += ins[kind][-1] >= 50
for kind in sorted(S):
    T = sum(1 for v in truth.values() if v[0] == kind); c = S[kind]
    line = f"[{MODE} d={div:.3f}] {kind:9s} reads {T:,}: " + ', '.join(f"{k} {v} ({100*v/T:.1f}%)" for k, v in c.items())
    if de[kind]: line += f" | median de {statistics.median(de[kind]):.4f}"
    if MODE == 'shuffled' and nint[kind]: line += f" | introns in primary median {statistics.median(nint[kind])} | largest insertion median {statistics.median(ins[kind])}"
    print(line)
