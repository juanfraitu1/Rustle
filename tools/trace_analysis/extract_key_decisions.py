import re
from pathlib import Path

loci = [
    ('44669742-44715973(+)', 'over_emission_44669742_L10__44406847-45079269.log'),
    ('17433800-17456446(+)', 'missing_STRG29_17433800__17190254-17521824.log'),
    ('110888244-110930143(+)', 'missing_110888244_L7__110886513-110932194.log'),
    ('57118036-57129805(+)', 'missing_57118036_L6__57112612-57238547.log'),
    ('20532687-20568095(-)', 'over_emission_20532687_L6__20126938-20733035.log'),
    ('40800807-40875796(-)', 'over_emission_40800807_L7__40800807-41473537.log'),
]

out_path = Path('/tmp/locus_traces/KEY_DECISIONS.md')
with open(out_path, 'w') as out:
    out.write('# StringTie key decisions per divergent locus\n\n')
    for (label, fname) in loci:
        tpath = Path('/tmp/locus_traces') / fname
        out.write(f'\n## {label}\n\n')
        if not tpath.exists():
            out.write('TRACE MISSING\n')
            continue
        # Extract FILTER decisions (who kills whom and why)
        out.write('### FILTER decisions (KILLED_BY lines)\n\n```\n')
        with open(tpath) as f:
            for line in f:
                if 'print_predcluster: FILTER' in line:
                    out.write(line)
        out.write('```\n\n')
        # Extract PRED_FATE summary
        out.write('### PRED_FATE summary\n\n```\n')
        with open(tpath) as f:
            for line in f:
                if 'print_predcluster: PRED_FATE' in line or 'print_predcluster: FINAL_FATE' in line:
                    out.write(line)
        out.write('```\n\n')
        # Extract FINAL_KEEP
        out.write('### FINAL_KEEP list\n\n```\n')
        with open(tpath) as f:
            for line in f:
                if 'print_predcluster: FINAL_KEEP' in line:
                    out.write(line)
        out.write('```\n\n')
        # Sample some BAD_JUNC / JUNC_DEMOTE
        out.write('### Junction kills (first 20)\n\n```\n')
        n = 0
        with open(tpath) as f:
            for line in f:
                if ('BAD_JUNC' in line or 'JUNC_DELETE' in line or 'JUNC_DEMOTE' in line) and '---' in line:
                    out.write(line)
                    n += 1
                    if n >= 20:
                        break
        out.write('```\n\n')

print(f'Wrote {out_path}')
import os
print(f'Size: {os.path.getsize(out_path)} bytes')
