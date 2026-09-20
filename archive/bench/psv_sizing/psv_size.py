import pysam, re, collections
BAM="/home/juanfra/winloci_scratch/GGO.bam"
GFF="/home/juanfra/winloci_scratch/GGO_tx.gff"
REC="/tmp/gw/rec_NC_073235.2.gtf"
CHROM="NC_073235.2"
B=(31684678,31709055)        # XM_055380753.2 (Copy B, - strand)
SIS=(32710000,32745000)      # sister copies (~32.71-32.74M)
TX="XM_055380753.2"
TOL=4
def match(i,S): return any(abs(i[0]-j[0])<=TOL and abs(i[1]-j[1])<=TOL for j in S)
# XM_055380753.2 introns from GFF (exons -> introns), 0-based (donor=end_1b, acc=start_1b-1)
ex=[]
for l in open(GFF):
    f=l.split("\t")
    if len(f)<9 or f[2]!="exon": continue
    if f"rna-{TX}" in f[8] or f'Name={TX}' in f[8] or f"GenBank:{TX}" in f[8]:
        ex.append((int(f[3]),int(f[4])))
ex.sort()
xm=set((ex[i][1], ex[i+1][0]-1) for i in range(len(ex)-1))
# recovery RSTL.* introns at the B locus
rtx=collections.defaultdict(list)
for l in open(REC):
    f=l.rstrip().split("\t")
    if len(f)<9 or f[2]!="exon": continue
    s,e=int(f[3]),int(f[4])
    if not (s< B[1] and e> B[0]): continue
    m=re.search(r'transcript_id "([^"]+)"',f[8]); rtx[m.group(1)].append((s,e))
recov=set()
for t,exs in rtx.items():
    exs.sort()
    for i in range(len(exs)-1): recov.add((exs[i][1], exs[i+1][0]-1))
# decisive-B reads (own - strand) introns
bam=pysam.AlignmentFile(BAM,"rb")
def best_as(reg):
    d={}
    for r in bam.fetch(CHROM,reg[0],reg[1]):
        if r.is_unmapped: continue
        try: a=r.get_tag("AS")
        except KeyError: continue
        if r.query_name not in d or a>d[r.query_name][0]: d[r.query_name]=(a,r)
    return d
asB=best_as(B); asS={k:v[0] for k,v in best_as(SIS).items()}
bread_introns=set(); ndec=0
for name,(a,r) in asB.items():
    if not r.is_reverse: continue            # XM is - strand: own-strand only
    s=asS.get(name)
    if s is not None and a<=s: continue       # not decisive for B
    ndec+=1
    blocks=r.get_blocks()
    for i in range(len(blocks)-1):
        bread_introns.add((blocks[i][1], blocks[i+1][0]))
# tally
J=len(xm)
covered=[i for i in xm if match(i,bread_introns)]    # XM junctions present in B-reads
recovered=[i for i in xm if match(i,recov)]          # XM junctions the recovery got
dropped=[i for i in xm if match(i,bread_introns) and not match(i,recov)]  # in B-reads but recovery missed
datamiss=[i for i in xm if not match(i,bread_introns)]                    # not in any B-read
print(f"XM_055380753.2 has {J} annotated junctions; recovery RSTL has {len(recov)} junctions; decisive-B own-strand reads={ndec}")
print(f"  XM junctions COVERED by decisive-B reads:        {len(covered)}/{J}")
print(f"  XM junctions the RECOVERY captured:              {len(recovered)}/{J}")
print(f"  in B-reads but DROPPED by assembly (fixable):    {len(dropped)}  <-- PSV-aware assembly could recover these")
print(f"  NOT in any B-read (data/identifiability limit):  {len(datamiss)}  <-- PSV-awareness cannot help these")

# refinement: XM junctions present in ANY read at the locus (true data ceiling)
any_introns=set()
for r in bam.fetch(CHROM,B[0],B[1]):
    if r.is_unmapped: continue
    blocks=r.get_blocks()
    for i in range(len(blocks)-1):
        any_introns.add((blocks[i][1], blocks[i+1][0]))
in_any   =[i for i in xm if match(i,any_introns)]
in_anyB  =[i for i in xm if match(i,any_introns)]  # placeholder
true_datamiss=[i for i in xm if not match(i,any_introns)]
in_reads_not_recovered=[i for i in xm if match(i,any_introns) and not match(i,recov)]
print("\n--- refined (vs ANY read at locus) ---")
print(f"  XM junctions in ANY read at locus:               {len(in_any)}/{J}")
print(f"  TRUE data limit (in NO read at all):             {len(true_datamiss)}/{J}  <-- unrecoverable by any method")
print(f"  in reads but NOT recovered (assignment/assembly):{len(in_reads_not_recovered)}/{J}  <-- PSV-aware assembly's addressable ceiling")
