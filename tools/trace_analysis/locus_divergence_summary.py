import re
import subprocess
from pathlib import Path

loci = [
    ('44669742', '44715973', '+', 'over_emission', 'over_emission_44669742_L10__44406847-45079269.log'),
    ('17433800', '17456446', '+', 'STRG29_missing', 'missing_STRG29_17433800__17190254-17521824.log'),
    ('110888244', '110930143', '+', 'missing_L7', 'missing_110888244_L7__110886513-110932194.log'),
    ('57118036', '57129805', '+', 'missing_L6', 'missing_57118036_L6__57112612-57238547.log'),
    ('20532687', '20568095', '-', 'over_emission_L6', 'over_emission_20532687_L6__20126938-20733035.log'),
    ('40800807', '40875796', '-', 'over_emission_L7', 'over_emission_40800807_L7__40800807-41473537.log'),
]

for (ls, le, strand, label, tfile) in loci:
    print(f'\n========== {label} {ls}-{le} ({strand}) ==========')
    tpath = Path('/tmp/locus_traces') / tfile
    
    # Final StringTie transcripts at this locus
    st_txs = []
    with open('GGO_19_stringtie.gtf') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p)<9 or p[2]!='transcript' or p[0]!='NC_073243.2': continue
            if p[6] != strand: continue
            s, e = int(p[3]), int(p[4])
            if s >= int(ls) - 10 and e <= int(le) + 10:
                m = re.search(r'transcript_id "([^"]+)"', p[8])
                cov = re.search(r'cov "([^"]+)"', p[8])
                exons = 0
                with open('GGO_19_stringtie.gtf') as g:
                    for l2 in g:
                        if m.group(1) in l2 and '\texon\t' in l2:
                            exons += 1
                st_txs.append((m.group(1), float(cov.group(1)) if cov else 0, s, e))
    
    # Rustle transcripts
    ru_txs = []
    with open('/tmp/v.gtf') as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p)<9 or p[2]!='transcript' or p[0]!='NC_073243.2': continue
            if p[6] != strand: continue
            s, e = int(p[3]), int(p[4])
            if s >= int(ls) - 10 and e <= int(le) + 10:
                m = re.search(r'transcript_id "([^"]+)"', p[8])
                cov = re.search(r'cov "([^"]+)"', p[8])
                ru_txs.append((m.group(1), float(cov.group(1)) if cov else 0, s, e))
    
    print(f'StringTie: {len(st_txs)} tx')
    for (tid, cv, s, e) in st_txs[:5]:
        print(f'  {tid} cov={cv:.2f} {s}-{e}')
    print(f'Rustle:    {len(ru_txs)} tx')
    for (tid, cv, s, e) in ru_txs[:5]:
        print(f'  {tid} cov={cv:.2f} {s}-{e}')
    
    # Extract key decisions from the StringTie locus trace
    if tpath.exists():
        kills = {'EJUNC_DELETE':0, 'JUNC_DELETE':0, 'BAD_JUNC':0,
                 'COLOR_BREAK_s0':0, 'JUNC_DEMOTE':0, 'LONGTRIM_SPLIT':0,
                 'SRCSINK_COVBUILD_SOURCE':0, 'SRCSINK_DECISION':0,
                 'FINAL_KEEP':0, 'FINAL_DROP':0,
                 'RETAINED_INTRON':0, 'INCLUDED':0}
        with open(tpath) as f:
            for line in f:
                if 'JUNC_COLOR_BREAK' in line and 'strand_now=0' in line:
                    kills['COLOR_BREAK_s0'] += 1
                elif 'EJUNC_DELETE' in line:
                    kills['EJUNC_DELETE'] += 1
                elif 'JUNC_DELETE' in line and 'EJUNC' not in line:
                    kills['JUNC_DELETE'] += 1
                elif 'BAD_JUNC' in line:
                    kills['BAD_JUNC'] += 1
                elif 'JUNC_DEMOTE' in line:
                    kills['JUNC_DEMOTE'] += 1
                elif 'LONGTRIM_SPLIT' in line:
                    kills['LONGTRIM_SPLIT'] += 1
                elif 'SRCSINK_COVBUILD_SOURCE' in line:
                    kills['SRCSINK_COVBUILD_SOURCE'] += 1
                elif 'SRCSINK_DECISION' in line:
                    kills['SRCSINK_DECISION'] += 1
                elif 'FINAL_FATE' in line and 'fate=KEEP' in line:
                    kills['FINAL_KEEP'] += 1
                elif 'FINAL_FATE' in line and 'fate=DROP' in line:
                    kills['FINAL_DROP'] += 1
                elif 'reason=RETAINED_INTRON' in line:
                    kills['RETAINED_INTRON'] += 1
                elif 'reason=INCLUDED' in line:
                    kills['INCLUDED'] += 1
        print(f'StringTie events in bundle: {dict((k,v) for k,v in kills.items() if v>0)}')
