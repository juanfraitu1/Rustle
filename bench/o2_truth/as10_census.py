#!/usr/bin/env python3
"""The 10-AS filter against the ms near tie, over a whole sample (register 1954). Reads the BAM's score table.

    as10_census.py BAM.as_table.tsv [LABEL]

Per read: AS lead = AS of the primary - AS of its best rival (the best-AS non-primary alignment); the rule fires at
>= 10. Near tie on ms = a second alignment within 10 ms of the best, or the primary is below the best ms.
Columns of the table (src/as_table.rs): 2 best AS, 6 primary AS, 7 best ms, 8 second ms, 9 primary ms, 11 primary
spliced, 12 rival AS, 14 rival spliced."""
import sys, collections

def lead_bin(l): return '< 0' if l < 0 else '0' if l == 0 else '1-9' if l < 10 else '10-19' if l < 20 else '20-49' if l < 50 else '>= 50'
def near(best, second, prim): return second >= best - 10 or prim < best

def census(path):
    c = collections.Counter()
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if f[11] == 'NA': c['single'] += 1; continue
        lead = int(f[5]) - int(f[11]); t = near(int(f[6]), int(f[7]), int(f[8]))
        c[('fires' if lead >= 10 else 'abstains', 'near' if t else 'clear')] += 1; c[('bin', lead_bin(lead), 'near' if t else 'clear')] += 1
        if lead >= 10 and t: c[('intron', f[10], f[13])] += 1; c[('fires_near_bin', lead_bin(lead))] += 1
    return c

def selftest():
    assert lead_bin(-3) == '< 0' and lead_bin(9) == '1-9' and lead_bin(10) == '10-19' and lead_bin(50) == '>= 50'
    assert near(100, 95, 100) and not near(100, 80, 100) and near(100, 80, 99)
    print('as10_census selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest(); sys.exit()
    c = census(sys.argv[1]); lab = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
    multi = sum(v for k, v in c.items() if len(k) == 2 and k[0] in ('fires', 'abstains')); nt = c[('fires', 'near')] + c[('abstains', 'near')]
    print(f'## {lab}: {c["single"] + multi} reads, {c["single"]} with one alignment, {multi} with several')
    print('cell\treads\tshare of reads with several alignments')
    for a in ('fires', 'abstains'):
        for b in ('clear', 'near'): print(f'rule {a}, ms {b}\t{c[(a, b)]}\t{c[(a, b)] / max(1, multi):.4f}')
    print(f'near-tied on ms: {nt}; the rule fires on {c[("fires", "near")]} of them ({c[("fires", "near")] / max(1, nt):.3f})')
    print('AS lead\tms clear\tms near tie')
    for b in ('< 0', '0', '1-9', '10-19', '20-49', '>= 50'): print(f'{b}\t{c[("bin", b, "clear")]}\t{c[("bin", b, "near")]}')
    fn = c[('fires', 'near')]
    print('rule fires AND near tie on ms, by (primary spliced, rival spliced): ' + '; '.join(f'{k[1]}/{k[2]} {v} ({v / max(1, fn):.2f})' for k, v in sorted(c.items(), key=lambda kv: -kv[1]) if k[0] == 'intron'))
    print('  their AS lead: ' + '; '.join(f'{b} {c[("fires_near_bin", b)]}' for b in ('10-19', '20-49', '>= 50')))
