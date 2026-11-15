#!/usr/bin/env python3
"""r1636: what StringTie did with the chains another assembly has and it lacks, from the modded build's decision log.

  st_fate.py RUN_JSONL CONTIG OURS_GTF STRINGTIE_GTF REFSEQ_GTF

Fate of each ours-only multi-exon chain on CONTIG, first that applies: killed:<reason> (pred_kill with that chain),
extracted_not_output (a flow path with that chain, no kill event), transfrag_only (a transfrag with that chain, never a path),
no_transfrag (StringTie never held that chain as one transfrag)."""
import sys, json, collections
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tool_match import read
log, contig = sys.argv[1:3]; ours, st, ref = (read(p) for p in sys.argv[3:6])
def key(s):
    try: return tuple((int(a) - 1, int(b)) for a, b in (x.split('-') for x in s.split(','))) if s else ()
    except ValueError: return ()  # single-exon or truncated field
killed = {}; path = set(); trf = collections.Counter(); bad = 0; good = set(); collapsed = {}
for l in open(log):
    if '"pred_kill"' in l or '"path_extracted"' in l or '"transfrag_define"' in l or '"good_junction"' in l or '"transfrag_collapse"' in l:
        try: e = json.loads(l, strict=False); p = e["payload"]
        except ValueError: bad += 1; continue  # a line cut by another thread's write
        if e['step'] == 'pred_kill': killed.setdefault(key(p.get('chain', '')), p['reason'])
        elif e['step'] == 'path_extracted': path.add(key(p.get('introns', '')))
        elif e['step'] == 'transfrag_define': trf[key(p.get('introns', ''))] += p.get('abund', 0)
        elif e['step'] == 'good_junction': good.add((e['start'], e['end'] - 1))
        elif e['step'] == 'transfrag_collapse':
            r = key(p.get('rep_introns', ''))
            for m in p.get('members', '').split(';' if ';' in p.get('members', '') else '|'):
                k = key(m.rsplit(':', 1)[0]) if ':' in m else ()
                if k and k != r: collapsed[k] = r
by_intron = collections.defaultdict(set)
for t in trf:
    for i in t: by_intron[i].add(t)
run = lambda a, b: any(b[k:k + len(a)] == a for k in range(len(b) - len(a) + 1))
def sub(ch):
    if ch in collapsed: return 'no_transfrag:collapsed_into_another'
    if any(i not in good for i in ch): return 'no_transfrag:junction_not_accepted'
    cand = set().union(*(by_intron.get(i, ()) for i in ch))
    if any(run(ch, t) for t in cand): return 'no_transfrag:inside_a_longer_transfrag'
    if any(run(t, ch) for t in cand): return 'no_transfrag:only_shorter_transfrags'
    return 'no_transfrag:other'
tab = collections.defaultdict(lambda: [0, 0]); ab = collections.defaultdict(list)
for c, s, ch in ours - st:
    if c != contig: continue
    f = 'killed:' + killed[ch] if ch in killed else 'extracted_not_output' if ch in path else 'transfrag_only' if ch in trf else sub(ch)
    tab[f][0] += 1; tab[f][1] += (c, s, ch) in ref
    if f == 'transfrag_only': ab[f].append(trf[ch])
n = sum(v[0] for v in tab.values()); print(f'unparsable log lines {bad}'); print(f'ours-only chains on {contig}: {n}\nfate\tchains\tshare\tcorrect\tcorrect_share')
for f, (k, ok) in sorted(tab.items(), key=lambda x: -x[1][0]): print(f'{f}\t{k}\t{k / n:.3f}\t{ok}\t{ok / k:.3f}')
a = sorted(ab['transfrag_only']); print('transfrag_only abundance quartiles', [a[len(a) * q // 4] for q in (1, 2, 3)] if a else [])
