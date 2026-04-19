import re
from pathlib import Path

# Target loci to slice
target_loci = [
    ('44406847-45079269', 'over_emission_44669742_L10'),
    ('17190254-17521824', 'missing_STRG29_17433800'),
    ('110886513-110932194', 'missing_110888244_L7'),
    ('57112612-57238547', 'missing_57118036_L6'),
    ('20126938-20733035', 'over_emission_20532687_L6'),
    ('40800807-41473537', 'over_emission_40800807_L7'),
]

log_path = 'stringtie_debug/trace_GGO_19_deep.log'
out_dir = Path('/tmp/locus_traces')
out_dir.mkdir(exist_ok=True)

# Map target bundle extents to line ranges
# Find PASS_ID lines and their ranges
pass_ranges = []  # (bundle_start_end, start_line, end_line)
last_pass = None
with open(log_path) as f:
    for ln, line in enumerate(f, 1):
        m = re.match(r'--- build_graphs: PASS_ID=\d+ bundle_start=(\d+) bundle_end=(\d+)', line)
        if m:
            if last_pass:
                pass_ranges.append((last_pass[0], last_pass[1], ln-1))
            last_pass = (f"{m.group(1)}-{m.group(2)}", ln)
if last_pass:
    pass_ranges.append((last_pass[0], last_pass[1], ln))

print(f'Indexed {len(pass_ranges)} bundles')

for (extent, label) in target_loci:
    for (ext, ls, le) in pass_ranges:
        if ext == extent:
            out_path = out_dir / f'{label}__{extent}.log'
            with open(log_path) as f:
                lines = []
                for ln, line in enumerate(f, 1):
                    if ls <= ln <= le:
                        lines.append(line)
                    elif ln > le:
                        break
            with open(out_path, 'w') as o:
                o.writelines(lines)
            print(f'{label} ({extent}): {len(lines)} lines -> {out_path}')
            break
    else:
        print(f'{label} ({extent}): NOT FOUND')
