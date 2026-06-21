import pysam, re, collections, json, random
random.seed(42)
bam=pysam.AlignmentFile("GGO.bam","rb")
DN=re.compile(r'DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$')

# Locus span: from the DN start, find read cluster extent. Reads at de-novo loci
# start near the given coordinate. We take reads overlapping [start-100, start+50kb]
# whose alignment START is within 5kb of the locus start (anchor), and take the
# span as the covered region. Cache spans.
span_cache={}
def locus_span(chrom, start):
    key=(chrom,start)
    if key in span_cache: return span_cache[key]
    ends=[]; starts=[]
    try:
        for r in bam.fetch(chrom, max(0,start-200), start+200000):
            if r.is_unmapped or r.is_supplementary: continue
            if abs(r.reference_start - start) <= 3000:
                starts.append(r.reference_start); ends.append(r.reference_end)
            if len(ends)>500: break
    except ValueError:
        span_cache[key]=None; return None
    if len(ends)<2:
        span_cache[key]=None; return None
    import statistics
    # robust end: use 90th percentile end to avoid runon
    ends.sort()
    e=ends[int(len(ends)*0.9)]
    s=min(starts)
    span_cache[key]=(s,e)
    return (s,e)

# For a pair, collect reads at locus A and locus B (within spans), keyed by query_name.
# A read "cross-maps" the pair if same query_name has alignments at BOTH loci
# (primary at one, secondary at other) AND at least one alignment is MAPQ0.
def fetch_locus(chrom, span):
    d=collections.defaultdict(list)
    for r in bam.fetch(chrom, span[0], span[1]):
        if r.is_unmapped or r.is_supplementary: continue
        # require substantial overlap with span
        ov=min(r.reference_end,span[1])-max(r.reference_start,span[0])
        if ov < 200: continue
        nm=r.get_tag("NM") if r.has_tag("NM") else None
        d[r.query_name].append((r.mapping_quality, nm, r.is_secondary, r.reference_start, r.reference_end))
    return d

def classify_pair(chrom, sA, sB):
    spanA=locus_span(chrom,sA); spanB=locus_span(chrom,sB)
    if spanA is None or spanB is None: return None
    # require non-overlapping loci
    if not (spanA[1] < spanB[0] or spanB[1] < spanA[0]): return None
    dA=fetch_locus(chrom,spanA); dB=fetch_locus(chrom,spanB)
    common=set(dA)&set(dB)
    # cross-mapping reads with a MAPQ0 alignment
    xmap=[]
    nm_same=0; nm_diff=0
    for q in common:
        a=dA[q]; b=dB[q]
        mqa=[x[0] for x in a]; mqb=[x[0] for x in b]
        if not (min(mqa)==0 or min(mqb)==0): continue  # must be ambiguous mapping
        # pick the alignment at each locus (best/primary): use the one with min NM available
        nma=[x[1] for x in a if x[1] is not None]
        nmb=[x[1] for x in b if x[1] is not None]
        if not nma or not nmb: continue
        # representative NM at each locus = min NM (best alignment there)
        na=min(nma); nb=min(nmb)
        xmap.append((q,na,nb))
        if na==nb: nm_same+=1
        else: nm_diff+=1
    n_x=len(xmap)
    if n_x==0: return {'chrom':chrom,'sA':sA,'sB':sB,'spanA':spanA,'spanB':spanB,'n_xmap':0,'nm_same':0,'nm_diff':0}
    return {'chrom':chrom,'sA':sA,'sB':sB,'spanA':spanA,'spanB':spanB,'n_xmap':n_x,'nm_same':nm_same,'nm_diff':nm_diff,
            'frac_same':nm_same/n_x}

# Sanity check on MAGEA pairs
for (s1,s2) in [(161251228,161458538),(164381222,164442447),(164397061,164426194)]:
    r=classify_pair("NC_073247.2",s1,s2)
    print("MAGEA",s1,s2,"->",{k:r[k] for k in('n_xmap','nm_same','nm_diff','frac_same')} if r and r['n_xmap']>0 else r)

# ---- Full run over co-located pairs that actually cross-map ----
pair_list=json.load(open('/tmp/coloc_pairs.json'))
print("total coloc pairs:",len(pair_list))
# Sample broadly but cap to keep runtime sane; iterate until we have >=300 cross-mapping pairs
random.shuffle(pair_list)
results=[]
xmap_pairs=[]   # pairs with >=1 cross-mapping MAPQ0 read (assignment-relevant)
checked=0
MIN_XMAP=3      # need at least 3 cross-mapping reads to classify
for (fid,chrom,sA,sB) in pair_list:
    if len(xmap_pairs)>=400: break
    checked+=1
    try:
        r=classify_pair(chrom,sA,sB)
    except Exception as e:
        continue
    if r is None: continue
    if r['n_xmap']>=MIN_XMAP:
        r['fid']=fid
        xmap_pairs.append(r)
    if checked%200==0:
        print(f"checked {checked}, xmap_pairs {len(xmap_pairs)}", flush=True)
print("CHECKED:",checked,"ASSIGNMENT-RELEVANT (>=%d xmap MAPQ0 reads):"%MIN_XMAP,len(xmap_pairs))
json.dump(xmap_pairs, open('/tmp/xmap_results.json','w'))
