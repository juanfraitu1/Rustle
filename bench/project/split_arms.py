#!/usr/bin/env python3
"""r1522: cut projected.gtf / base.gtf into the judged sets H (outside chr13:13-16 Mb) and G (target T >= 2).
usage: split_arms.py PANEL_DIR -> PANEL_DIR/{H,G}/{base,tier_own,tier_tied,tier_none}.gtf + arms.txt (reach, motifs)"""
import sys, re, csv, os, collections
d = sys.argv[1]
tg = {(r['family'], r['copy']): r for r in csv.DictReader(open(f'{d}/targets.tsv'), delimiter='\t')}
hot = lambda r: r['chrom'] == 'chr13' and 13_000_000 <= int(r['start']) < 16_000_000
keep = {'H': lambda r: not hot(r), 'G': lambda r: int(r['T']) >= 2}
tid2tg = {t: k for k, r in tg.items() if r['tids'] for t in r['tids'].split(',')}
out = open(f'{d}/arms.txt', 'w')
for arm, ok in keep.items():
    os.makedirs(f'{d}/{arm}', exist_ok=True)
    judged = {k for k, r in tg.items() if ok(r)}
    ws = {k: open(f'{d}/{arm}/{k}.gtf', 'w') for k in ('base', 'tier_own', 'tier_tied', 'tier_none')}
    gained = collections.defaultdict(set); mo = mn = 0
    for l in open(f'{d}/projected.gtf'):
        f = l.split('\t'); a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        fam, mi = a['gene_id'].rsplit('_', 1)
        if (fam, mi) not in judged: continue
        if f[2] == 'transcript':  # exon lines carry no tier: keep the transcript's
            tier = a['tier']; ok_, n = map(int, a['motif'].split('/')); mo += ok_; mn += n; gained[(fam, mi)].add(tier)
        ws['tier_' + tier].write(l)
    for l in open(f'{d}/base.gtf'):
        t = re.search(r'transcript_id "([^"]+)"', l.split('\t')[8]).group(1)
        if tid2tg.get(t) in judged: ws['base'].write(l)
    for h in ws.values(): h.close()
    low = [k for k in judged if int(tg[k]['reads']) <= 2]
    print(f"[{arm}] targets judged {len(judged)} of {len(tg)}; bar1 reach {sum(1 for k in low if gained[k])}/{len(low)} = {sum(1 for k in low if gained[k])/max(len(low),1):.3f}; "
          f"bar4 motifs {mo}/{mn} = {mo/max(mn,1):.3f}; tiers {collections.Counter(t for v in gained.values() for t in v)}", file=out)
out.close(); print(open(f'{d}/arms.txt').read())
