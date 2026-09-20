#!/usr/bin/env python3
"""Classify rustle-only intron chains (vs StringTie) as TP/FP against the real RefSeq
annotation, split by extraction source (flow vs checktrf). Verifies the checktrf store-gate
fix impact. Usage: python3 bench/classify_checktrf_tpfp.py <ru.gtf> <st.gtf> <ru.jsonl> <annot.gff> <chrom>
"""
import re, sys, json, collections

def chains_gtf(path):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def chains_gff(path, chrom):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon': continue
        m = re.search(r'Parent=([^;]+)', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def safe(p):
    for l in open(p):
        try: yield json.loads(l)
        except ValueError: continue

ru_gtf, st_gtf, ru_jsonl, annot, chrom = sys.argv[1:6]
ru = chains_gtf(ru_gtf); st = chains_gtf(st_gtf); ann = chains_gff(annot, chrom)
ru_only = ru - st
parse = lambda s: tuple((int(a), int(b)) for a, b in re.findall(r'(\d+)-(\d+)', s))
src = collections.defaultdict(set)
for d in safe(ru_jsonl):
    if d.get('step') != 'path_extracted': continue
    pl = d['payload']; ch = tuple((a-1, b+1) for a, b in parse(pl.get('introns', '')))
    src[(d['strand'], ch)].add(pl.get('source', '?'))

res = collections.Counter()
for k in ru_only:
    s = src.get(k, set())
    has_ck = any('checktrf' in x for x in s); has_flow = any(x == 'flow' for x in s)
    bucket = 'checktrf_only' if (has_ck and not has_flow) else \
             'checktrf+flow' if (has_ck and has_flow) else 'flow_only'
    res[(bucket, 'TP' if k in ann else 'FP')] += 1
print(f"rustle-only chains: {len(ru_only)}  (annotation chains: {len(ann)})")
for b in ('flow_only', 'checktrf_only', 'checktrf+flow'):
    print(f"  {b:14s}: TP={res[(b,'TP')]:3d}  FP={res[(b,'FP')]:3d}")
