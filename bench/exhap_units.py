"""EXHAP units (LEDGER_O3): each family's exon units = distinct sets of reference loci of its autosomal exons (1416's -x sr mapping,
C 0.8, T 0.90), with GAIN truth (mat or pat holds more copies than ref) per unit. Writes units.tsv, copies.tsv, truth.tsv."""
import collections, os, sys
D = sys.argv[1]  # the sdcn directory
SEXREF = {'NC_073247.2', 'NC_073248.2', 'NC_011120.1'}
def sex(g, t): return t in SEXREF if g == 'ref' else t.split('_')[0] in ('chrX', 'chrY', 'chrM')
H = {g: collections.defaultdict(list) for g in ('ref', 'mat', 'pat')}
for g in H:
    for l in open(f'{D}/ex_sr_{g}.paf'):
        f = l.split('\t')
        if (int(f[3]) - int(f[2])) / int(f[1]) >= 0.8 and int(f[9]) / int(f[10]) >= 0.90 and not sex(g, f[5]):
            H[g][f[0]].append((f[5], int(f[7]), int(f[8]), f[4]))
def merged(xs):
    out = []
    for t, s, e, st in sorted(xs):
        if out and out[-1][0] == t and s < out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([t, s, e, st])
    return out
split = {l.split()[0]: l.split()[2] for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bpclust_truth.tsv')) if not l.startswith('family')}
units = collections.defaultdict(dict)
for u in H['ref']:
    fam = u.split('|')[1]
    loci = merged(H['ref'][u])
    key = tuple((t, s // 50, e // 50) for t, s, e, _ in loci)
    n = [len(loci), len(merged(H['mat'][u])), len(merged(H['pat'][u]))]
    gain = max(n[1], n[2]) > n[0]
    if key in units[fam]:
        units[fam][key][1] |= gain
    else:
        units[fam][key] = [loci, gain, n]
with open('units.tsv', 'w') as o, open('copies.tsv', 'w') as c, open('truth.tsv', 'w') as tr:
    o.write('unit\tfamily\tsplit\tn_ref\tn_mat\tn_pat\tgain\n')
    c.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\texons\n')
    tr.write('family\tsplit\tgain\tn_units\n')
    for fam in sorted(units, key=lambda x: int(x[5:])):
        for k, (key, (loci, gain, n)) in enumerate(sorted(units[fam].items())):
            uid = f'{fam}x{k}'
            o.write(f'{uid}\t{fam}\t{split[fam]}\t{n[0]}\t{n[1]}\t{n[2]}\t{int(gain)}\n')
            for i, (t, s, e, st) in enumerate(loci):
                c.write(f'{uid}\t{i}\t{uid}:{i}\t{t}\t{s}\t{e}\t1\t{st}\t{s}-{e}\n')
        tr.write(f"{fam}\t{split[fam]}\t{'GAIN' if any(v[1] for v in units[fam].values()) else 'NOGAIN'}\t{len(units[fam])}\n")
