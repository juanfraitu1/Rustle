#!/usr/bin/env python3
"""Per-locus digests for the 24 StringTie single-exon transcripts, to diagnose what
seeds outcompete the SE candidate. Reads /tmp/ru_se_arch.jsonl (rustle parity log,
SE-flow-arch ON) + /tmp/st_all.gtf + /tmp/ru_se_arch.gtf. Writes /tmp/se_loci/<i>.txt."""
import re, json, os, collections

ST_SE = [
 (18162670,18166999,'+'),(18843129,18845947,'+'),(23303485,23306475,'-'),(25185146,25188394,'-'),
 (27405028,27412101,'+'),(30128450,30132405,'+'),(30184086,30193392,'+'),(30816172,30821379,'-'),
 (32312926,32315519,'+'),(33325221,33325548,'-'),(36671347,36673170,'+'),(40126487,40126701,'-'),
 (40716493,40719871,'-'),(69130713,69134243,'-'),(82898818,82904516,'+'),(86703361,86708995,'+'),
 (95384571,95387943,'+'),(101437484,101441874,'-'),(103053309,103057580,'+'),(103250946,103254957,'-'),
 (103639711,103642383,'+'),(110800078,110803739,'+'),(111761013,111767629,'+'),(111870669,111873292,'+'),
]

def load_gtf(path):
    tx=collections.defaultdict(list); strand={}
    for line in open(path):
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id "([^"]+)"', f[8]);
        if not m: continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); strand[m.group(1)]=f[6]
    return [(sorted(ex)[0][0],sorted(ex)[-1][1],strand[t],len(ex),t) for t,ex in tx.items()]

st_tx=load_gtf('/tmp/st_all.gtf'); ru_tx=load_gtf('/tmp/ru_se_arch.gtf')

evs=[]
for ln in open('/tmp/ru_se_arch.jsonl'):
    ln=ln.strip()
    if not ln: continue
    try: e=json.loads(ln)
    except: continue
    if e.get('step')=='_meta': continue
    evs.append(e)

def overlap(s,e,a,b): return min(e,b)-max(s,a)

os.makedirs('/tmp/se_loci', exist_ok=True)
for i,(ss,se,strd) in enumerate(ST_SE):
    pad=3000
    sub=[e for e in evs if e.get('strand')==strd and overlap(ss,se,e.get('start',0),e.get('end',0))> -pad and e.get('start',0)<se+pad and e.get('end',0)>ss-pad]
    L=[]
    L.append(f"ST SE #{i}: {ss}-{se} ({strd})  len={se-ss+1}")
    # rustle output at locus
    L.append("--- rustle FINAL output overlapping (>30%) ---")
    for (s,e,st2,nex,t) in ru_tx:
        if overlap(ss,se,s,e)>0.3*(se-ss): L.append(f"  RU {s}-{e} {st2} {nex}ex {t}")
    L.append("--- ST output overlapping ---")
    for (s,e,st2,nex,t) in st_tx:
        if overlap(ss,se,s,e)>0.3*(se-ss): L.append(f"  ST {s}-{e} {st2} {nex}ex {t}")
    # graph
    gl=[e for e in sub if e['step']=='graphnode_list']
    L.append("--- graph(s) covering locus (start,end,n_nodes) ---")
    for e in gl: L.append(f"  {e.get('start')}-{e.get('end')} n_nodes={e['payload'].get('n_nodes')}")
    # SEEDS overlapping the SE region (the competition)
    L.append("--- parse_trflong_seed overlapping SE region (the competing seeds) ---")
    for e in sub:
        if e['step']=='parse_trflong_seed' and overlap(ss,se,e.get('start',0),e.get('end',0))>0:
            p=e['payload']
            L.append(f"  seed {e['start']}-{e['end']} abund={p.get('entry_abund')} n_introns={p.get('n_introns')} introns={p.get('introns','')[:90]}")
    # path_extracted overlapping
    L.append("--- path_extracted overlapping SE region ---")
    for e in sub:
        if e['step']=='path_extracted' and overlap(ss,se,e.get('start',0),e.get('end',0))>0:
            p=e['payload']
            L.append(f"  path {e['start']}-{e['end']} src={p.get('source')} nexons={p.get('nexons')} cov={p.get('cov')} lc={p.get('longcov')}")
    # pred_kill overlapping
    L.append("--- pred_kill overlapping ---")
    for e in sub:
        if e['step']=='pred_kill' and overlap(ss,se,e.get('start',0),e.get('end',0))>0:
            p=e['payload']
            L.append(f"  KILL {e['start']}-{e['end']} reason={p.get('reason')} nexons={p.get('nexons')} cov={p.get('cov')}")
    open(f'/tmp/se_loci/{i}.txt','w').write('\n'.join(L)+'\n')
print(f"wrote {len(ST_SE)} SE-locus digests to /tmp/se_loci/")
