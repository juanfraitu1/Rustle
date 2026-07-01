import csv, json, pysam
from collections import defaultdict
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
bam=pysam.AlignmentFile(BAM,"rb"); refs=set(bam.references)
def mapq0frac(chrom,s,e,cap=3000):
    if chrom not in refs: return (0,0)
    n=m0=0
    for a in bam.fetch(chrom,max(0,s),min(e,s+400000)):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        n+=1
        if a.mapping_quality==0: m0+=1
        if n>=cap: break
    return (n,m0)

# ---- (A) MAPQ0 over KNOWN-ambiguous SUN co-located catalog copies ----
sun=list(csv.DictReader(open('bench/sun_identifiability.tsv'),delimiter='\t'))
tier_n=defaultdict(int); tier_m0=defaultdict(int)
tot_n=tot_m0=0
for r in sun:
    n,m0=mapq0frac(r['chrom'],int(r['start']),int(r['end']))
    tier_n[r['tier']]+=n; tier_m0[r['tier']]+=m0
    tot_n+=n; tot_m0+=m0
print("=== (A) MAPQ0 over SUN co-located catalog (known E_c-territory) ===")
print("ALL sun copies: reads=%d mapq0=%d = %.1f%% ambiguous"%(tot_n,tot_m0,100*tot_m0/max(1,tot_n)))
for t in sorted(tier_n):
    print("  tier%s: reads=%d mapq0=%d = %.1f%% ambiguous  (n_copies=%d)"%(
        t,tier_n[t],tier_m0[t],100*tier_m0[t]/max(1,tier_n[t]),sum(1 for r in sun if r['tier']==t)))

# ---- (B) E_c subset E_r : do co-located catalog families' copies co-group in denovo E_r? ----
# map each sun copy (chrom,start,end) to a denovo locus by overlap; check same denovo family
META="/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
dn=[]
for r in csv.DictReader(open(META),delimiter='\t'):
    dn.append((r['id'],r['chrom'],int(r['start']),int(r['end'])))
dn_by_chrom=defaultdict(list)
for i in dn: dn_by_chrom[i[1]].append(i)
loc2fam={}
for r in csv.DictReader(open('bench/denovo_families.tsv'),delimiter='\t'):
    for m in r['members'].split(','): loc2fam[m]=r['family_id']
def find_dn(chrom,s,e):
    best=None;bov=0
    for (i,c,ds,de) in dn_by_chrom.get(chrom,[]):
        ov=min(e,de)-max(s,ds)
        if ov>bov: bov=ov;best=i
    return best if bov>0 else None
sun_by_fam=defaultdict(list)
for r in sun: sun_by_fam[r['family']].append(r)
same=diff=unmapped=nofam=0; multi_copy_checked=0
no_er_family=0            # dnfams == 0 (no multi-copy E_r label -> UNEVALUABLE, not a drop)
split_details=[]          # dnfams >= 2 (genuine split across >=2 E_r families)
for fam,copies in sun_by_fam.items():
    if len(copies)<2: continue
    dnfams=defaultdict(list)   # dnfam_id -> [dn_locus,...] for the mapped copies
    genes=set()
    mapped=0
    for r in copies:
        d=find_dn(r['chrom'],int(r['start']),int(r['end']))
        if d is None: unmapped+=1; continue
        mapped+=1
        genes.add(r.get('gene',''))
        f=loc2fam.get(d)
        if f is None: nofam+=1; continue
        dnfams[f].append(d)
    if mapped>=2:
        multi_copy_checked+=1
        nf=len(dnfams)
        if nf==1: same+=1
        elif nf==0: no_er_family+=1
        else:
            diff+=1
            split_details.append((fam,sorted(genes),dict(dnfams)))
evaluable=same+diff        # families that DO carry a multi-copy E_r label
print("\n=== (B) E_c-territory (SUN co-located) families vs denovo E_r grouping ===")
print("multi-copy sun families with >=2 copies mapped to denovo loci: %d"%multi_copy_checked)
print("  all copies in ONE denovo E_r family (E_r contains E_c grouping): %d"%same)
print("  copies split across >=2 denovo E_r families (genuine split)   : %d"%diff)
print("  no multi-copy E_r label at all (UNEVALUABLE, not an E_c drop)  : %d"%no_er_family)
print("  => containment rate among evaluable: %d/%d = %.1f%%"%(same,evaluable,100*same/max(1,evaluable)))
print("  sun copies with no overlapping denovo locus (not de-novo expressed): %d"%unmapped)
print("  denovo locus not in any multi-copy E_r family: %d"%nofam)

# ---- (C) Correctly attribute the splits: EDGE_LINKED (a cross-family core_recip edge exists,
#          i.e. E_r component over-fragmentation) vs OPERATIONAL-SHARED-EXON leak (NO core_recip edge,
#          yet the loci de-tie -> the shipped symmetric core_recip>=0.13 oracle MISSES a genuine E_c edge;
#          absorbed only by the DEFINITIONAL permissive-local E_r^asym; R cannot fix it, R only splits). ----
DELTA=0.005; DE_MAX=0.05; MIN_READS=3
coord={i[0]:(i[1],i[2],i[3]) for i in dn}   # dn_locus -> (chrom,start,end)
# homology edges present in the operational E_r oracle (core_recip>=0.13 -> connected components = families)
er_edge=set()
try:
    for r in csv.DictReader(open('bench/denovo_family_edges.tsv'),delimiter='\t'):
        er_edge.add(frozenset((r['a'],r['b'])))
except FileNotFoundError:
    er_edge=None
def detie_reads(a_loc,b_loc,maxwin=400000,qcap=8000):
    """Faithful read_conflict.rs de_tied count between two DN loci (best/min de per read per locus)."""
    best=defaultdict(dict)
    for idx,loc in enumerate((a_loc,b_loc)):
        c,s,e=coord[loc]
        if c not in refs: continue
        n=0
        for al in bam.fetch(c,max(0,s),min(e,s+maxwin)):
            if al.is_unmapped: continue
            try: de=al.get_tag('de')
            except KeyError: continue
            cur=best[al.query_name].get(idx)
            if cur is None or de<cur: best[al.query_name][idx]=de
            n+=1
            if n>=qcap: break
    n=0
    for qn,dd in best.items():
        if 0 in dd and 1 in dd and abs(dd[0]-dd[1])<=DELTA and max(dd[0],dd[1])<=DE_MAX:
            n+=1
    return n
print("\n=== (C) attribution of the %d genuine splits (EDGE_LINKED vs operational-oracle shared-exon leak) ==="%len(split_details))
n_edge_linked=n_oracle_leak=0
for fam,genes,dnf in split_details:
    fam_ids=sorted(dnf)
    # any operational core_recip edge crossing the two E_r families?
    cross_edge=False
    if er_edge is not None:
        for i in range(len(fam_ids)):
            for j in range(i+1,len(fam_ids)):
                for a_loc in dnf[fam_ids[i]]:
                    for b_loc in dnf[fam_ids[j]]:
                        if frozenset((a_loc,b_loc)) in er_edge: cross_edge=True
    # de-tie reads between the representative loci of the split E_r families
    counts=[]
    for i in range(len(fam_ids)):
        for j in range(i+1,len(fam_ids)):
            a_loc=dnf[fam_ids[i]][0]; b_loc=dnf[fam_ids[j]][0]
            counts.append("%s~%s:%d"%(fam_ids[i],fam_ids[j],detie_reads(a_loc,b_loc)))
    label="EDGE_LINKED" if cross_edge else "OPERATIONAL-SHARED-EXON-LEAK"
    if cross_edge: n_edge_linked+=1
    else: n_oracle_leak+=1
    print("  sunfam %s genes=%s  E_r-families=%s  cross_core_recip_edge=%s  de_tie_reads=[%s]  -> %s"%(
        fam,",".join(g for g in genes if g),fam_ids,cross_edge,", ".join(counts),label))
print("  splits EDGE_LINKED (core_recip edge exists, R-fixable over-fragmentation): %d"%n_edge_linked)
print("  splits OPERATIONAL-SHARED-EXON-LEAK (no core_recip edge; miss of the 0.13 oracle,")
print("     absorbed only by definitional permissive-local E_r^asym; NOT R-fixable): %d"%n_oracle_leak)

# ---- (D) DECISIVE: do fully-Tier-1 SUN families actually form a de-tie edge (SURVIVE E_c),
#          or do they "vanish under E_c-as-O1"? A single distinguishing SUN moves de by ~1/read_len
#          (<= delta=0.005), so SUN-resolvable copies STILL de-tie. Test it directly: run the faithful
#          read_conflict.rs de_tied predicate over each fully-Tier-1 family's own co-located loci. ----
def family_ec_visible(loci,maxwin=400000,qcap=8000):
    """EC_VISIBLE iff some locus-pair accumulates >=MIN_READS de-tie reads (best/min de per read per locus)."""
    best=defaultdict(dict)
    for idx,loc in enumerate(loci):
        c,s,e=coord[loc]
        if c not in refs: continue
        n=0
        for al in bam.fetch(c,max(0,s),min(e,s+maxwin)):
            if al.is_unmapped: continue
            try: de=al.get_tag('de')
            except KeyError: continue
            cur=best[al.query_name].get(idx)
            if cur is None or de<cur: best[al.query_name][idx]=de
            n+=1
            if n>=qcap: break
    pair=Counter()
    for qn,dd in best.items():
        if len(dd)<2: continue
        it=list(dd.items())
        for x in range(len(it)):
            for y in range(x+1,len(it)):
                (li,di),(lj,dj)=it[x],it[y]
                if abs(di-dj)<=DELTA and max(di,dj)<=DE_MAX:
                    pair[(min(li,lj),max(li,lj))]+=1
    return any(v>=MIN_READS for v in pair.values())
from collections import Counter
SIZE_CAP=40
n_visible=n_dropped=n_skip=n_uncov=0
for fam,copies in sun_by_fam.items():
    if len(copies)<2: continue
    if not all(r['tier']=='1' for r in copies): continue   # fully-Tier-1 only
    loci=[]
    for r in copies:
        d=find_dn(r['chrom'],int(r['start']),int(r['end']))
        if d and d in coord and d not in loci: loci.append(d)
    if len(copies)>SIZE_CAP: n_skip+=1; continue
    if len(loci)<2: n_uncov+=1; continue
    if family_ec_visible(loci): n_visible+=1
    else: n_dropped+=1
ev=n_visible+n_dropped
print("\n=== (D) do fully-Tier-1 SUN families form a de-tie edge? (refutes 'fully-Tier-1 => edgeless-E_c => vanish') ===")
print("fully-Tier-1 SUN families: EC_VISIBLE(form a de-tie edge, SURVIVE E_c)=%d  EC_DROPPED(edgeless)=%d  skip(size>%d)=%d  uncov(<2 mapped loci)=%d"%(
    n_visible,n_dropped,SIZE_CAP,n_skip,n_uncov))
if ev: print("  => %d/%d = %.1f%% of evaluable fully-Tier-1 families are E_c-VISIBLE (a de-tie edge forms). EC_DROPPED=%d."%(
    n_visible,ev,100*n_visible/ev,n_dropped))
