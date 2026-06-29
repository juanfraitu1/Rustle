#!/usr/bin/env python3
"""Slice per-locus digests for the coupled-port understanding workflow.

Reads /tmp/divergent_loci.json + /tmp/{st,ru}_chr19.jsonl and writes one digest
per locus to /tmp/loci/<idx>.txt summarizing: divergent chains (RU-only/ST-only),
both tools' graph node counts, path_extracted (source/cov/introns), pred_kill, and
checktrf_result counts — everything needed to classify the divergence mechanism.
"""
import json, os, collections

loci = json.load(open('/tmp/divergent_loci.json'))

def load_events(path):
    evs = []
    for ln in open(path):
        ln = ln.strip()
        if not ln:
            continue
        try:
            e = json.loads(ln)
        except Exception:
            continue
        if e.get('step') == '_meta':
            continue
        evs.append(e)
    return evs

st_evs = load_events('/tmp/st_chr19.jsonl')
ru_evs = load_events('/tmp/ru_chr19.jsonl')

def overlap(e, s, en):
    es, ee = e.get('start', 0), e.get('end', 0)
    return es <= en and ee >= s

os.makedirs('/tmp/loci', exist_ok=True)

def fmt_path(e):
    p = e.get('payload', {})
    intr = p.get('introns', '')
    return f"  [{e['strand']}] {e.get('start')}-{e.get('end')} src={p.get('source')} cov={p.get('cov')} lc={p.get('longcov')} nintr={len(intr.split(',')) if intr else 0} introns={intr}"

def fmt_kill(e):
    p = e.get('payload', {})
    return f"  KILL [{e['strand']}] {e.get('start')}-{e.get('end')} reason={p.get('reason')} cov={p.get('cov')} killer_cov={p.get('killer_cov')} killer={p.get('killer_start')}-{p.get('killer_end')}"

manifest = []
for idx, L in enumerate(loci):
    s, en, strand = L['start'], L['end'], L['strand']
    # pad range slightly
    s2, en2 = s - 50, en + 50
    sub_st = [e for e in st_evs if e.get('strand') == strand and overlap(e, s2, en2)]
    sub_ru = [e for e in ru_evs if e.get('strand') == strand and overlap(e, s2, en2)]
    cat = ('MIXED' if L['n_ru'] and L['n_st'] else ('RU_PURE' if L['n_ru'] else 'ST_PURE'))
    lines = []
    lines.append(f"LOCUS {idx} [{strand}] {s}-{en}  category={cat}  RU_only={L['n_ru']} ST_only={L['n_st']}")
    lines.append("--- DIVERGENT CHAINS ---")
    for m in L['members']:
        lines.append(f"  {m['who']}-only [{m['who']}] {m['start']}-{m['end']} nintr={m['nintr']} introns={m['introns']}")
    # graph node counts
    def gnodes(evs):
        gl = [e for e in evs if e['step'] == 'graphnode_list']
        return [(e.get('start'), e.get('end'), e['payload'].get('n_nodes')) for e in gl]
    lines.append("--- GRAPHS (start,end,n_nodes) ---")
    lines.append(f"  ST: {gnodes(sub_st)}")
    lines.append(f"  RU: {gnodes(sub_ru)}")
    # path_extracted
    lines.append("--- ST path_extracted ---")
    for e in sub_st:
        if e['step'] == 'path_extracted':
            lines.append(fmt_path(e))
    lines.append("--- RU path_extracted ---")
    for e in sub_ru:
        if e['step'] == 'path_extracted':
            lines.append(fmt_path(e))
    # kills
    st_k = [e for e in sub_st if e['step'] == 'pred_kill']
    ru_k = [e for e in sub_ru if e['step'] == 'pred_kill']
    lines.append(f"--- ST pred_kill ({len(st_k)}) ---")
    for e in st_k:
        lines.append(fmt_kill(e))
    lines.append(f"--- RU pred_kill ({len(ru_k)}) ---")
    for e in ru_k:
        lines.append(fmt_kill(e))
    # checktrf summary
    def ck(evs):
        c = collections.Counter()
        for e in evs:
            if e['step'] == 'checktrf_result':
                c[str(e['payload'].get('stored', e['payload'].get('result', '?')))] += 1
        return dict(c)
    lines.append(f"--- checktrf_result (stored counts) ST={ck(sub_st)} RU={ck(sub_ru)} ---")
    open(f'/tmp/loci/{idx}.txt', 'w').write('\n'.join(lines) + '\n')
    manifest.append({'idx': idx, 'cat': cat, 'strand': strand, 'start': s, 'end': en,
                     'n_ru': L['n_ru'], 'n_st': L['n_st'],
                     'st_npaths': sum(1 for e in sub_st if e['step'] == 'path_extracted'),
                     'ru_npaths': sum(1 for e in sub_ru if e['step'] == 'path_extracted')})

json.dump(manifest, open('/tmp/loci/manifest.json', 'w'), indent=0)
cats = collections.Counter(m['cat'] for m in manifest)
print(f"wrote {len(manifest)} locus digests to /tmp/loci/  categories={dict(cats)}")
