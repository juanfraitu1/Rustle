#!/usr/bin/env python3
"""L1/L2 scoring (PREREG docs/PREREG_loose_ends_L1_L2_2026-09-05.md).

usage: o2_l1l2_score.py <sweep_dir> [<sweep_dir_B>] [--families <list>]
Per sweep family: unit reads (a primary block inside a unit; truth-by-placement = unit with max block overlap),
status counts, origin-rejected, orphans, reads whose placement unit is a DROPPED member (their statuses and
whether an assigned read went to that unit), MAPQ-60 placement agreement. With a second dir the same families are
compared paired. For fam_MCL1_073242 the 62 audited junction anchors are scored (copies mapped by containment
from o2scale/fam_NPIPcore_073242).
"""
import csv, collections, glob, os, sys, pysam
M = '/mnt/linuxdisk/home/juanfraitu/mcl_ann'
BAM = pysam.AlignmentFile(os.environ.get('O2_BAM', '/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam'))
args = sys.argv[1:]; fams = None
if '--families' in args:
    i = args.index('--families'); fams = set(open(args[i + 1]).read().split()); del args[i:i + 2]
dirs = args

def score(d):
    cp = {r['copy_idx']: r for r in csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')}
    A = {r['read_name']: r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')}
    truth = {}; mapq = {}
    for i, r in cp.items():
        c, s, e = r['chrom'], int(r['start']), int(r['end'])
        for a in BAM.fetch(c, s, e):
            if a.flag & 2308 or a.query_name not in A: continue
            o = sum(min(b1, e) - max(b0, s) for b0, b1 in a.get_blocks() if b1 > s and b0 < e)
            if o > 0 and (a.query_name not in truth or o > truth[a.query_name][1]):
                truth[a.query_name] = (i, o); mapq[a.query_name] = a.mapping_quality
    st = collections.Counter(A[n]['status'] for n in truth)
    rej = sum(1 for n in truth if A[n].get('origin_rejected') == '1')
    orph = sum(1 for n in truth if A[n].get('n_candidates') == '0')
    q60 = [n for n in truth if mapq[n] >= 60]; asg60 = [n for n in q60 if A[n]['status'] == 'assigned']
    agree = sum(1 for n in asg60 if A[n]['catalog_copy_idx'] == truth[n][0])
    low = [n for n in truth if mapq[n] < 60]; low_asg = sum(1 for n in low if A[n]['status'] == 'assigned')
    dropped = {i for i, r in cp.items() if r.get('member_status') == 'dropped'}
    dr = [n for n in truth if truth[n][0] in dropped]
    dst = collections.Counter(A[n]['status'] for n in dr)
    d_to_it = sum(1 for n in dr if A[n]['status'] == 'assigned' and A[n]['catalog_copy_idx'] == truth[n][0])
    d_asg = dst['assigned']
    # reads assigned TO a dropped unit whose placement is elsewhere (the dropped unit stealing)
    stolen = sum(1 for n in truth if A[n]['status'] == 'assigned' and A[n]['catalog_copy_idx'] in dropped and truth[n][0] not in dropped)
    params = dict(l.rstrip('\n').split('\t', 1) for l in open(f'{d}/A.params.tsv') if '\t' in l)
    return dict(unit_reads=len(truth), assigned=st['assigned'], tied=st['tied'], amb=st['ambiguous'], rejected=rej, orphans=orph,
                q60=len(q60), q60_assigned=len(asg60), q60_agree=agree, low=len(low), low_assigned=low_asg,
                n_dropped_units=len(dropped), dropped_reads=len(dr), dropped_assigned=d_asg, dropped_to_it=d_to_it, dropped_tied=dst['tied'], dropped_amb=dst['ambiguous'], stolen_by_dropped=stolen,
                catalog_locus=params.get('read_star_catalog_locus', 'NA'), rejected_all=int(params.get('origin_rejected', 0)), orphans_all=int(params.get('orphans', 0)))

def anchors(d):
    core = {r['copy_idx']: (r['chrom'], int(r['start']), int(r['end'])) for r in csv.DictReader(open(f'{M}/o2scale/fam_NPIPcore_073242/copies.tsv'), delimiter='\t')}
    sw = {r['copy_idx']: (r['chrom'], int(r['start']), int(r['end'])) for r in csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')}
    ovl = lambda a, b: a[0] == b[0] and a[1] < b[2] and b[1] < a[2]
    m = {}
    for ci, cc in core.items():
        hits = [vi for vi, vv in sw.items() if ovl(cc, vv)]; m[ci] = hits[0] if len(hits) == 1 else None
    truth = [r for r in csv.DictReader(open(f'{M}/o2scale/truth_valid.tsv'), delimiter='\t')]
    A = {r['read_name']: r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')}
    st = collections.Counter(); wrong = []; right = 0
    for t in truth:
        r = A.get(t['read'])
        if r is None: st['absent'] += 1; continue
        st[r['status']] += 1
        if r['status'] == 'assigned':
            exp = m.get(t['anchor_copy'])
            if exp is None: st['unmappable'] += 1
            elif r['catalog_copy_idx'] == exp: right += 1
            else: wrong.append((t['read'][-8:], t['mapq'], exp, r['catalog_copy_idx'], r.get('min_p_value')))
    return dict(st), right, wrong

keys = ['unit_reads', 'assigned', 'tied', 'amb', 'rejected', 'orphans', 'q60', 'q60_assigned', 'q60_agree', 'low', 'low_assigned', 'n_dropped_units', 'dropped_reads', 'dropped_assigned', 'dropped_to_it', 'dropped_tied', 'dropped_amb', 'stolen_by_dropped']
tot = {d: collections.Counter() for d in dirs}
fam_list = sorted(os.path.basename(p) for p in glob.glob(f'{dirs[0]}/fam_*') if os.path.exists(f'{p}/A.done'))
if fams is not None: fam_list = [f for f in fam_list if f in fams]
for f in fam_list:
    if any(not os.path.exists(f'{d}/{f}/A.assignments.tsv') for d in dirs): continue
    rows = {d: score(f'{d}/{f}') for d in dirs}
    for d in dirs:
        for k in keys: tot[d][k] += rows[d][k]
    if rows[dirs[0]]['n_dropped_units'] or f in ('fam_MCL1_073242', 'fam_MCL7_073242'):
        for d in dirs:
            r = rows[d]; print(f"{f} [{os.path.basename(d)} catalog_locus={r['catalog_locus']}]: unit reads {r['unit_reads']} | assigned {r['assigned']} tied {r['tied']} amb {r['amb']} | rejected {r['rejected']} (all rows {r['rejected_all']}) orphans {r['orphans']} | q60 agree {r['q60_agree']}/{r['q60_assigned']} | MAPQ<60 {r['low']}: assigned {r['low_assigned']} | dropped units {r['n_dropped_units']}: reads {r['dropped_reads']} assigned {r['dropped_assigned']} (to it {r['dropped_to_it']}) tied {r['dropped_tied']} amb {r['dropped_amb']}; stolen by dropped {r['stolen_by_dropped']}")
    if f == 'fam_MCL1_073242':
        for d in dirs:
            st, right, wrong = anchors(f'{d}/{f}'); print(f"  NPIP 62 anchors [{os.path.basename(d)}]: {st} | right {right} wrong {len(wrong)} {wrong[:6]}")
print(f'=== totals over {len(fam_list)} families')
for d in dirs:
    t = tot[d]; print(f"[{os.path.basename(d)}] unit reads {t['unit_reads']} | assigned {t['assigned']} ({100*t['assigned']/max(t['unit_reads'],1):.1f}%) tied {t['tied']} amb {t['amb']} | rejected {t['rejected']} orphans {t['orphans']} | q60 agree {t['q60_agree']}/{t['q60_assigned']} = {t['q60_agree']/max(t['q60_assigned'],1):.4f} | MAPQ<60 {t['low']}: assigned {t['low_assigned']} | dropped units {t['n_dropped_units']}: reads {t['dropped_reads']} assigned {t['dropped_assigned']} to-it {t['dropped_to_it']} tied {t['dropped_tied']} amb {t['dropped_amb']}; stolen {t['stolen_by_dropped']}")
