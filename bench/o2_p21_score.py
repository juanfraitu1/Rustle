# P21 scorer: copy_phasing O2 vs minimap2's primary on the P18 test reads
import csv, re, subprocess, collections
truth = {}
for l in open('/mnt/linuxdisk/tmp/p18/test.truth.tsv'):
    if l.startswith('read\t'): continue
    f = l.rstrip('\n').split('\t'); truth[f[0]] = f
prim = {}
for ln in subprocess.run(['samtools', 'view', '-F', '0x904', '/mnt/linuxdisk/tmp/p18/test.bam'], capture_output=True, text=True, check=True).stdout.splitlines():
    f = ln.split('\t'); t = truth[f[0]]; s = int(f[3]) - 1
    e = s + sum(int(n) for n, o in re.findall(r'(\d+)([MDN=X])', f[5]))
    prim[f[0]] = f[2] == t[3] and s < int(t[5]) and e > int(t[4])
p18 = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('/mnt/linuxdisk/tmp/p18/test.rows.tsv') if not l.startswith('read\t')}
conf_err = {r for r, x in p18.items() if x[5] == 'False' and x[6] == 'False' and float(x[8]) >= 0.99}
c = collections.Counter(); asg = {}
for r in csv.DictReader(open('O.o2.tsv'), delimiter='\t'):
    n = r['read']; t = truth[n]; g = t[6]
    for grp in ('all', g): c[(grp, 'tied rows')] += 1
    k = r['assigned']
    if k in ('ABSTAIN', 'HOPELESS'):
        for grp in ('all', g): c[(grp, k)] += 1
        continue
    right = k == f'{t[1]}:{t[2]}'; asg[n] = right
    for grp in ('all', g):
        c[(grp, 'assigned')] += 1; c[(grp, 'right')] += right
        c[(grp, 'fixed')] += right and not prim[n]; c[(grp, 'broken')] += (not right) and prim[n]
for k in sorted(c): print(k, c[k])
fixed40 = sum(1 for n in conf_err if asg.get(n))
tied40 = sum(1 for n in conf_err if n in {r['read'] for r in csv.DictReader(open('O.o2.tsv'), delimiter='\t')})
final = sum(asg[n] if n in asg else prim[n] for n in truth); M = sum(prim.values())
print(f'the 40 confident other-allele errors: in a tied row {tied40}, fixed {fixed40}')
a = c[('all', 'assigned')]; rt = c[('all', 'right')]
print(f"H1 precision {rt}/{a} = {rt / a if a else float('nan'):.4f} (>= .98): {'PASS' if a and rt / a >= .98 else 'FAIL'}")
print(f"H2 fixed {c[('all', 'fixed')]} > broken {c[('all', 'broken')]} and fixed40 {fixed40} >= 1: {'PASS' if c[('all', 'fixed')] > c[('all', 'broken')] and fixed40 >= 1 else 'FAIL'}")
print(f"H3 final {final} >= M {M}: {'PASS' if final >= M else 'FAIL'}")
fam = collections.Counter()
for r in csv.DictReader(open('O.families.tsv'), delimiter='\t'):
    fam['families'] += 1; fam['haplotypes'] += int(r['haplotypes']); fam['certified'] += int(r['certified'])
print(dict(fam), 'putative rows', sum(1 for _ in open('O.putative.tsv')) - 1)
