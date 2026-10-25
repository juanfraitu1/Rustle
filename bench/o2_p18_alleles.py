# P18: each catalog copy's allele on Jim's OTHER haplotype (the reference chromosome IS one haplotype's chromosome).
# Orthology: the copy's transcript hit on the other hap's matching chromosome lying between the primary asm5 hits of its 20-kb flanks.
import re, collections
cm = {l.split()[0]: l.rstrip().split('\t')[1:] for l in open('chrom_map.tsv')}
comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan'); rc = lambda s: s.translate(comp)[::-1]
seqs = {}; cur = None
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.fa'):
    if l[0] == '>': cur = l[1:].split()[0]; seqs[cur] = []
    else: seqs[cur].append(l.strip())
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
fl = collections.defaultdict(dict)
for h in ('mat', 'pat'):
    for l in open(f'fl.{h}.paf'):
        f = l.split('\t')
        if 'tp:A:P' not in l or int(f[11]) < 30: continue
        q = f[0].split('|'); key = (q[0] + '|' + q[1], h)
        if (int(f[3]) - int(f[2])) / int(f[1]) < 0.5: continue
        fl[key][q[2]] = (f[5], int(f[7]), int(f[8]))
def apply_cs(q, cs, qs, qe):
    out = [q[:qs]]; i = qs
    for op, val in re.findall(r'([:*+\-~])([0-9a-z]+)', cs):
        if op == ':': n = int(val); out.append(q[i:i + n]); i += n
        elif op == '*': out.append(val[0].upper()); i += 1
        elif op == '+': i += len(val)
        elif op == '-': out.append(val.upper())
    assert i == qe, (i, qe)
    out.append(q[qe:]); return ''.join(out)
hits = collections.defaultdict(list)
for h in ('mat', 'pat'):
    for l in open(f'tx.{h}.paf'):
        f = l.rstrip('\n').split('\t')
        cs = next(x[5:] for x in f[12:] if x.startswith('cs:Z:'))
        hits[(f[0], h)].append((f[5], int(f[7]), int(f[8]), f[4], int(f[2]), int(f[3]), int(f[1]), int(f[9]), int(f[10]), cs))
stat = collections.Counter(); out = open('alleles.fa', 'w'); tab = open('alleles.tsv', 'w')
print('copy\tref_chrom\tother_hap_chrom\tallele_edits\tlen_ref\tlen_allele\tstatus', file=tab)
for q, s in seqs.items():
    fam, idx, loc = q.split('|')[:3]; c = loc.split(':')[0]
    if c not in cm or cm[c][1] == '-': stat['sex chromosome / unplaced: no other allele'] += 1; continue
    oc = cm[c][1]; h = oc.split('_')[1]
    f2 = fl.get((fam + '|' + idx, h), {})
    anch = [v for v in f2.values() if v[0] == oc]
    if not anch: stat['no flank anchor on the other haplotype'] += 1; print(q, c, oc, '', len(s), '', 'no_anchor', sep='\t', file=tab); continue
    lo = min(v[1] for v in anch) - 50000; hi = max(v[2] for v in anch) + 50000
    cand = [x for x in hits[(q, h)] if x[0] == oc and x[1] < hi and x[2] > lo and (x[5] - x[4]) / x[6] >= 0.95]
    if len(cand) != 1:
        stat[f'{"no" if not cand else "several"} transcript hit inside the anchored window'] += 1
        print(q, c, oc, '', len(s), '', 'no_hit' if not cand else 'ambiguous', sep='\t', file=tab); continue
    t, ts, te, strand, qs, qe, ql, mt, bl, cs = cand[0]
    a = apply_cs(s, cs, qs, qe) if strand == '+' else rc(apply_cs(rc(s), cs, ql - qe, ql - qs))
    e = bl - mt
    stat['allele built: identical' if a == s else 'allele built: differs'] += 1
    print(q, c, oc, e, len(s), len(a), 'ok', sep='\t', file=tab)
    out.write(f'>{q}|{oc}:{ts}-{te}\n{a}\n')
for k, v in sorted(stat.items()): print(f'{k}: {v}')
