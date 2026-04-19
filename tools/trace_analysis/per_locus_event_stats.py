import re
import sys
from collections import defaultdict

# Read StringTie bundles from trace log, get event counts per bundle
bundles = []  # (pass_id, start, end, lines)
cur_bundle = None
cur_events = defaultdict(int)

target_loci = [
    (44669742, 44715973, 2, 16, '+'),
    (110888244, 110930143, 13, 6, '+'),
    (40800807, 40875796, 8, 15, '-'),
    (57118036, 57129805, 6, 0, '+'),
    (20532687, 20568095, 1, 7, '-'),
    (17433800, 17456446, 5, 0, '+'),
]

def classify(line):
    line = line.strip()
    if 'JUNC_COLOR_BREAK' in line:
        m = re.search(r'strand_now=(-?\d+)', line)
        if m:
            return f'COLOR_BREAK_s{m.group(1)}'
    for kw in ['EJUNC_DELETE','BAD_JUNC','JUNC_DELETE','SRCSINK_COVBUILD_SOURCE','SRCSINK_DECISION','JUNC_DEMOTE']:
        if kw in line:
            return kw
    if 'good_junc: ACCEPTED' in line:
        return 'GJ_ACCEPT'
    if line.startswith('--- good_junc: KILL') or 'good_junc: KILL' in line:
        return 'GJ_KILL'
    if 'longtrim' in line.lower() and 'SPLIT' in line:
        return 'LONGTRIM_SPLIT'
    return None

with open('stringtie_debug/trace_GGO_19_deep.log') as f:
    for ln, line in enumerate(f, 1):
        m = re.match(r'--- build_graphs: PASS_ID=(\d+) bundle_start=(\d+) bundle_end=(\d+)', line)
        if m:
            if cur_bundle:
                bundles.append((cur_bundle, dict(cur_events)))
            cur_bundle = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
            cur_events = defaultdict(int)
            continue
        c = classify(line)
        if c:
            cur_events[c] += 1

if cur_bundle:
    bundles.append((cur_bundle, dict(cur_events)))

# Find target bundles
for (ls, le, st_n, ru_n, strand) in target_loci:
    print(f'\n=== Locus {ls}-{le} ({strand}) ST={st_n} RU={ru_n} ===')
    for ((pid, bs, be), ev) in bundles:
        if bs <= ls and be >= le:
            print(f'  Bundle PASS_ID={pid} {bs}-{be}')
            for k in sorted(ev.keys()):
                print(f'    {k}: {ev[k]}')
            break
