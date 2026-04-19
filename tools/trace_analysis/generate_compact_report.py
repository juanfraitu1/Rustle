import re
from pathlib import Path
from collections import Counter

loci = [
    ('44669742-44715973(+)', 'over_emission_44669742_L10__44406847-45079269.log'),
    ('17433800-17456446(+)', 'missing_STRG29_17433800__17190254-17521824.log'),
    ('110888244-110930143(+)', 'missing_110888244_L7__110886513-110932194.log'),
    ('57118036-57129805(+)', 'missing_L6_missing_57118036_L6__57112612-57238547.log'.replace('missing_L6_','')),
    ('20532687-20568095(-)', 'over_emission_20532687_L6__20126938-20733035.log'),
    ('40800807-40875796(-)', 'over_emission_40800807_L7__40800807-41473537.log'),
]

out = Path('/tmp/locus_traces/COMPACT_SUMMARY.md')
with open(out, 'w') as o:
    o.write('# Compact divergence summary per locus\n\n')
    o.write('Per locus: StringTie final output, kill reasons, and notable patterns.\n\n')
    for (label, fname) in loci:
        tpath = Path('/tmp/locus_traces') / fname
        o.write(f'## {label}\n\n')
        if not tpath.exists():
            o.write('TRACE NOT FOUND\n\n')
            continue
        # Kill reasons histogram
        reasons = Counter()
        final_keep = []
        final_drop = []
        with open(tpath) as f:
            for line in f:
                m = re.search(r'reason=(\w+)', line)
                if m:
                    reasons[m.group(1)] += 1
                if 'FINAL_FATE' in line and 'fate=KEEP' in line:
                    mm = re.search(r'n=(\d+) (\d+-\d+) cov=([\d.]+) strand=(\S+) exons=(\d+)', line)
                    if mm:
                        final_keep.append((int(mm.group(1)), mm.group(2), float(mm.group(3)), mm.group(4), int(mm.group(5))))
                if 'FINAL_FATE' in line and 'fate=DROP' in line:
                    mm = re.search(r'n=(\d+) (\d+-\d+) cov=([\d.]+) strand=(\S+) exons=(\d+)', line)
                    if mm:
                        final_drop.append((int(mm.group(1)), mm.group(2), float(mm.group(3)), mm.group(4), int(mm.group(5))))
        o.write('### StringTie kill-reason histogram\n\n')
        o.write('| reason | count |\n|---|---|\n')
        for r, c in sorted(reasons.items(), key=lambda x: -x[1]):
            o.write(f'| {r} | {c} |\n')
        o.write('\n')
        o.write(f'### FINAL decisions: {len(final_keep)} KEEP, {len(final_drop)} DROP\n\n')
        o.write('#### Kept\n\n')
        for (n, r, c, s, e) in final_keep:
            o.write(f'- pred[{n}] {r} cov={c:.2f} {s} exons={e}\n')
        o.write('\n#### Dropped (sample first 10)\n\n')
        for (n, r, c, s, e) in final_drop[:10]:
            o.write(f'- pred[{n}] {r} cov={c:.2f} {s} exons={e}\n')
        o.write('\n')
print(f'Wrote {out}')
