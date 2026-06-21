import pysam, re, collections, json, random
random.seed(42)
bam=pysam.AlignmentFile("/home/juanfra/winloci_scratch/GGO.bam","rb")
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

import re as _re, os as _os, subprocess as _sub

FAM_TSV = "/mnt/c/Users/jfris/Desktop/Rustle/bench/denovo_families.tsv"
GFF     = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
OUT_TSV = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_resolution_census.tsv"

# Memoize the BAM-fetch hot path. classify_pair calls fetch_locus by global name (late binding),
# so rebinding the global here makes exhaustive pairing bounded by distinct loci, not pair count.
_fetch_cache = {}
_orig_fetch_locus = fetch_locus
def fetch_locus(chrom, span):
    key = (chrom, span)
    if key not in _fetch_cache:
        _fetch_cache[key] = _orig_fetch_locus(chrom, span)
    return _fetch_cache[key]

def _parse_member(m):
    mm = _re.match(r"DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$", m)
    return (mm.group(1), int(mm.group(2))) if mm else None

def _pair_class(p):
    if p["n_xmap"] < 3: return "low_support"
    if p["frac_same"] >= 1.0: return "k0_strict"
    if p["frac_same"] >= 0.95: return "k0"
    return "resolvable"

def per_family_census(win=2_000_000):
    """Per-family K=0 labels over ALL co-located copy pairs (same chrom, within `win`, threshold n_xmap>=1).
    low_support pairs (n_xmap<3) are counted in n_colocated_pairs but excluded from the family verdict."""
    rows, all_pairs = [], []
    for line in open(FAM_TSV):
        if line.startswith("family_id"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        fid, ncopies, members = parts[0], parts[1], parts[2]
        locs = [x for x in (_parse_member(m) for m in members.split(",")) if x]
        by_chrom = {}
        for c, s in locs: by_chrom.setdefault(c, []).append(s)
        classed = []
        for c, starts in by_chrom.items():
            starts = sorted(set(starts))
            for i in range(len(starts)):
                for j in range(i + 1, len(starts)):
                    if starts[j] - starts[i] > win: break  # sorted -> no closer pair beyond the window
                    try:
                        p = classify_pair(c, starts[i], starts[j])
                    except Exception:
                        continue
                    if p is None or p["n_xmap"] < 1: continue
                    p["fid"] = fid
                    classed.append(p); all_pairs.append(p)
        if not classed: continue
        cls   = [_pair_class(p) for p in classed]
        n_k0  = cls.count("k0") + cls.count("k0_strict")
        n_k0s = cls.count("k0_strict")
        n_res = cls.count("resolvable")
        n_ar  = n_k0 + n_res
        if   n_ar == 0:      verdict = "not_assignment_relevant"
        elif n_k0 == n_ar:   verdict = "k0"
        elif n_res == n_ar:  verdict = "resolvable"
        else:                verdict = "mixed"
        k0_pairs = ";".join(
            f"{p['chrom']}:{p['spanA'][0]}-{p['spanA'][1]}~{p['spanB'][0]}-{p['spanB'][1]}"
            for p, k in zip(classed, cls) if k in ("k0", "k0_strict"))
        rows.append(dict(family_id=fid, n_copies=int(ncopies), n_colocated_pairs=len(classed),
                         n_assignment_relevant=n_ar, n_resolvable=n_res, n_k0=n_k0, n_k0_strict=n_k0s,
                         family_verdict=verdict, k0_pairs=k0_pairs))
    cols = ["family_id","n_copies","n_colocated_pairs","n_assignment_relevant","n_resolvable",
            "n_k0","n_k0_strict","family_verdict","k0_pairs"]
    with open(OUT_TSV, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    return rows, all_pairs

def check_census():
    rows, pairs = per_family_census()
    ar   = [p for p in pairs if _pair_class(p) in ("k0", "k0_strict", "resolvable")]  # n_xmap>=3
    n    = len(ar)
    res  = sum(1 for p in ar if _pair_class(p) == "resolvable")
    k0   = sum(1 for p in ar if _pair_class(p) in ("k0", "k0_strict"))
    reads = sum(p["n_xmap"] for p in ar); diff = sum(p["nm_diff"] for p in ar)
    assert n >= 200, f"expected a broad assignment-relevant set, got {n}"
    assert 0.70 <= res / n <= 0.85, f"resolvable fraction off: {res/n:.2f}"
    assert 0.78 <= diff / reads <= 0.88, f"read-level resolvable off: {diff/reads:.2f}"
    dnfam1 = [r for r in rows if r["family_id"] == "DNFAM1"]
    assert dnfam1, "DNFAM1 missing from census"
    assert dnfam1[0]["family_verdict"] in ("k0", "mixed"), f"DNFAM1 verdict={dnfam1[0]['family_verdict']}"
    assert dnfam1[0]["n_k0"] >= 2, f"DNFAM1 must carry >=2 K0 pairs (pair2/pair3), got n_k0={dnfam1[0]['n_k0']}"
    for r in rows:
        assert r["n_resolvable"] + r["n_k0"] == r["n_assignment_relevant"], f"{r['family_id']} verdict-count mismatch"
        assert r["n_assignment_relevant"] <= r["n_colocated_pairs"]
        assert r["n_k0_strict"] <= r["n_k0"]
    print(f"OK  - census: {len(rows)} families; assignment-relevant {n} pairs -> {res} resolvable "
          f"({100*res/n:.0f}%), {k0} K0; reads {100*diff/reads:.0f}% resolvable; "
          f"DNFAM1 verdict={dnfam1[0]['family_verdict']} n_k0={dnfam1[0]['n_k0']}")
    return rows

if __name__ == "__main__":
    check_census()
