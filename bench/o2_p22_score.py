import sys, csv, collections
truth = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('/mnt/linuxdisk/tmp/p18/dev.truth.tsv') if not l.startswith('read\t')}
panel = {r['family']: r for r in csv.DictReader(open('panel.tsv'), delimiter='\t')}
for arm in sys.argv[1:]:
    reads = collections.defaultdict(list)
    for r in csv.DictReader(open(f'{arm}.reads.tsv'), delimiter='\t'): reads[r['cluster']].append(r['read'])
    tp = fp = 0; found = set(); detail = []
    for c in csv.DictReader(open(f'{arm}.clusters.tsv'), delimiter='\t'):
        if c['fate'] != 'new_copy': continue
        fam = c['family']; m = panel[fam]['masked_copy']
        src = collections.Counter((truth[r][1], truth[r][2]) for r in reads[c['cluster']])
        hit = sum(k for (f, cp), k in src.items() if f == fam and cp == m)
        ok = m != '-' and hit * 2 > len(reads[c['cluster']])
        tp += ok; fp += not ok
        if ok: found.add(fam)
        else: detail.append(f"{c['cluster']}({panel[fam]['arm']},{dict(src.most_common(2))})")
    print(f"{arm}: certified {tp+fp}  true {tp}  false {fp}  masked copies found {len(found)}/70  false: {' '.join(detail)}")
