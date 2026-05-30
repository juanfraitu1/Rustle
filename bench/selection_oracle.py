#!/usr/bin/env python3
"""Selection-ceiling oracle: if Rustle kept exactly ST's winners on candidate-matching
clusters, what is the chain-level F1? Bounds sub-project-1's payoff before building.
Inputs: /tmp/ru.gtf /tmp/st.gtf ../GGO_19.gtf, /tmp/ru_pe.jsonl /tmp/st_pe.jsonl.

Coordinate convention (VERIFIED Step 1): path_extracted jsonl `introns` are the
actual intron boundaries (first..last intron base, inclusive), so jsonl (d,a) maps
to GTF (exon_end,next_start) = (d-1,a+1). Confirmed on two chains on NC_073243.2.
"""
import json, re, collections

def gtf_chains(p):
    tx=collections.defaultdict(list); strand={}
    for line in open(p):
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id "([^"]+)"',f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); strand[m.group(1)]=f[6]
    out={}
    for t,ex in tx.items():
        ex.sort()
        if len(ex)<2: continue
        out[(strand[t],tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1)))]=(ex[0][0],ex[-1][1])
    return out

def pe_chains(p):
    s=set()
    for line in open(p):
        try: e=json.loads(line)
        except Exception: continue
        if e.get("step")!="path_extracted": continue
        istr=e["payload"].get("introns","")
        if not istr: continue
        # VERIFIED: jsonl introns (d,a) == GTF (exon_end+1,next_start-1); map to GTF coords
        ch=tuple((int(d)-1,int(a)+1) for d,a in (tok.split("-") for tok in istr.split(",")))
        s.add((e["strand"],ch))
    return s

ref=gtf_chains('../GGO_19.gtf'); ru=gtf_chains('/tmp/ru.gtf'); st=gtf_chains('/tmp/st.gtf')
ru_pe=pe_chains('/tmp/ru_pe.jsonl'); st_pe=pe_chains('/tmp/st_pe.jsonl')

# NOTE: the path_extracted jsonls are chrom-filtered to NC_073243.2, while ru/st/ref GTFs
# are genome-wide. Baseline TP/FN/FP are reported genome-wide. A cluster can only become
# candidate-matching where BOTH jsonls have entries (i.e. on NC_073243.2), so the oracle's
# substitution is naturally confined to that chrom; off-chrom clusters keep Rustle's finals.

# cluster by overlapping span+strand over the union of all final+candidate chains
def cluster(chains_with_span):
    items=sorted((s,a,b,ch) for (s,ch),(a,b) in chains_with_span.items())
    cl={}; cur=None; lid=-1
    for s,a,b,ch in items:
        if cur is None or s!=cur[0] or a>cur[2]: lid+=1; cur=[s,a,b,lid]
        else: cur[2]=max(cur[2],b)
        cl[(s,ch)]=lid
    return cl
# build span map for all chains (finals)
allspan={**{(s,ch):v for (s,ch),v in ru.items()}, **{(s,ch):v for (s,ch),v in st.items()}, **{(s,ch):v for (s,ch),v in ref.items()}}
clid=cluster(allspan)

# candidate-matching clusters: ru_pe and st_pe identical within the cluster
pe_by_cl=collections.defaultdict(lambda:[set(),set()])
for k in ru_pe:
    if k in clid: pe_by_cl[clid[k]][0].add(k)
for k in st_pe:
    if k in clid: pe_by_cl[clid[k]][1].add(k)
matching_cl={c for c,(r,s) in pe_by_cl.items() if r==s and r}

# ORACLE output: on matching clusters keep ST's finals; elsewhere keep Rustle's finals
oracle=set()
for (s,ch),(a,b) in ru.items():
    c=clid.get((s,ch))
    if c in matching_cl:
        continue  # drop Rustle's winner; ST's will be added below
    oracle.add((s,ch))
for (s,ch),(a,b) in st.items():
    c=clid.get((s,ch))
    if c in matching_cl:
        oracle.add((s,ch))
refset=set(ref)
tp=len(oracle&refset); fn=len(refset-oracle); fp=len(oracle-refset)
f1=200*tp/(2*tp+fn+fp)

ru_tp=len(set(ru)&refset); ru_fn=len(refset-set(ru)); ru_fp=len(set(ru)-refset)
print(f"matching clusters: {len(matching_cl)}")
print(f"BASELINE Rustle: TP={ru_tp} FN={ru_fn} FP={ru_fp}")
print(f"ORACLE (ST winners on matching clusters): TP={tp} FN={fn} FP={fp} F1={f1:.2f}")
print(f"FP removed (baseline_FP - oracle_FP): {ru_fp-fp}")
print(f"FN recovered (baseline_FN - oracle_FN): {ru_fn-fn}")

# --- Sanity check: print a couple matching clusters and confirm candidates identical ---
print("\n--- SANITY: sample matching clusters (ru_pe == st_pe) ---")
shown=0
for c in sorted(matching_cl):
    r,s=pe_by_cl[c]
    assert r==s and r
    strand=next(iter(r))[0]
    print(f"cluster {c} strand={strand} #candidates={len(r)} (ru==st: {r==s})")
    for k in sorted(r):
        print(f"    cand introns: {k[1]}")
    shown+=1
    if shown>=2: break
