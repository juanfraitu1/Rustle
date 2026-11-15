"""EXHAP units (bench/exhap_units.py): each family's exon units = distinct sets of reference loci of its autosomal exons, with
GAIN truth (mat or pat holds more copies than ref). No dev/held-out split (Amendment 1): split = all. The last column of
units.tsv names the unit's first exon (the one bench/exceil.py keeps). usage: exhap_units.py DIR SEX.tsv"""
import collections, sys
from common import sexset, load_hits, merged
D = sys.argv[1]; H = load_hits(D, sexset(sys.argv[2]))
units = collections.defaultdict(dict)
for u in H['ref']:
    fam = u.split('|')[1]
    loci = merged(H['ref'][u])
    key = tuple((t, s // 50, e // 50) for t, s, e, _ in loci)
    n = [len(loci), len(merged(H['mat'][u])), len(merged(H['pat'][u]))]
    gain = max(n[1], n[2]) > n[0]
    if key in units[fam]: units[fam][key][1] |= gain
    else: units[fam][key] = [loci, gain, n, u]
with open(f'{D}/units.tsv', 'w') as o, open(f'{D}/copies.tsv', 'w') as c, open(f'{D}/truth.tsv', 'w') as tr:
    o.write('unit\tfamily\tsplit\tn_ref\tn_mat\tn_pat\tgain\texon\n')
    c.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\texons\n')
    tr.write('family\tsplit\tgain\tn_units\n')
    for fam in sorted(units, key=lambda x: (len(x), x)):
        for k, (key, (loci, gain, n, ex)) in enumerate(sorted(units[fam].items(), key=lambda kv: kv[0])):
            uid = f'{fam}x{k}'
            o.write(f'{uid}\t{fam}\tall\t{n[0]}\t{n[1]}\t{n[2]}\t{int(gain)}\t{ex}\n')
            for i, (t, s, e, st) in enumerate(loci):
                c.write(f'{uid}\t{i}\t{uid}:{i}\t{t}\t{s}\t{e}\t1\t{st}\t{s}-{e}\n')
        tr.write(f"{fam}\tall\t{'GAIN' if any(v[1] for v in units[fam].values()) else 'NOGAIN'}\t{len(units[fam])}\n")
nu = sum(len(v) for v in units.values()); ng = sum(v[1] for f in units.values() for v in f.values())
print(f'units {nu} in {len(units)} families; gaining units {ng}; GAIN families {sum(any(v[1] for v in f.values()) for f in units.values())}')
