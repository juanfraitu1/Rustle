#!/usr/bin/env python3
"""Stage 1: extract EVERY BAM-derivable feature (except the PSV alleles) for DSFAM42 reads, labelled by
their PSV-assignment (assigned copy = ground-ish truth; tied = unassignable). Stage 2 (sklearn) asks:
can any of these features separate copies? Writes /tmp/dsfam42_features.tsv."""
import pysam, math

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
REGION = ("NC_073247.2", 59690000, 59790000)
# the 9 de-novo copy spans (from detect_and_assign stderr) for 5'/3'-offset features
COPIES = [(59690204,59723101),(59696701,59705136),(59706220,59715078),(59715729,59724620),
          (59725264,59733595),(59744368,59751694),(59753367,59761248),(59772963,59781277),(59782481,59789847)]

status, copy = {}, {}
for line in open("/tmp/dsfam42.assignments.tsv"):
    c = line.rstrip("\n").split("\t")
    if c[0] == "read_name": continue
    status[c[0]] = c[3]; copy[c[0]] = c[2]   # status, assigned_copy

def best_copy_idx(rs, re):
    best, bo = -1, 0
    for i, (s, e) in enumerate(COPIES):
        ov = min(re, e) - max(rs, s)
        if ov > bo: bo, best = ov, i
    return best

def seq_entropy(seq):
    if not seq: return 0.0
    from collections import Counter
    n = len(seq); c = Counter(seq)
    return -sum((v/n)*math.log2(v/n) for v in c.values() if v)

bam = pysam.AlignmentFile(BAM, "rb")
rows = {}
for a in bam.fetch(*REGION):
    if a.is_unmapped or a.query_name not in status: continue
    exonic = sum(l for op, l in (a.cigartuples or []) if op in (0,7,8))
    if a.query_name in rows and exonic <= rows[a.query_name][0]: continue   # keep longest aln per read
    ni = sum(1 for op,l in (a.cigartuples or []) if op==3)
    sc = sum(l for op,l in (a.cigartuples or []) if op in (4,5))
    sc5 = (a.cigartuples[0][1] if a.cigartuples and a.cigartuples[0][0] in (4,5) else 0)
    sc3 = (a.cigartuples[-1][1] if a.cigartuples and a.cigartuples[-1][0] in (4,5) else 0)
    tags = dict(a.get_tags())
    seq = a.query_sequence or ""
    q = a.query_qualities
    gc = (seq.count("G")+seq.count("C"))/len(seq) if seq else 0
    ci = best_copy_idx(a.reference_start, a.reference_end)
    s5, s3 = (COPIES[ci][0], COPIES[ci][1]) if ci >= 0 else (REGION[1], REGION[2])
    rows[a.query_name] = (exonic, ni+1, a.reference_end-a.reference_start, sc/max(1,a.infer_read_length() or len(seq) or 1),
        sc5, sc3, int(a.is_reverse), a.mapping_quality,
        float(tags.get("AS", 0)), int(tags.get("NM", 0)), float(tags.get("de", 0)),
        gc, (sum(q)/len(q) if q else 0), (min(q) if q else 0), seq_entropy(seq),
        a.reference_start - s5, s3 - a.reference_end, ci)

with open("/tmp/dsfam42_features.tsv", "w") as o:
    o.write("read\tstatus\tassigned_copy\texonic\tn_exons\tgspan\tsoftclip_frac\tsc5\tsc3\tis_rev\tmapq\t"
            "AS\tNM\tde\tgc\tqual_mean\tqual_min\tseq_entropy\toffset5\toffset3\tbest_copy\n")
    for r, f in rows.items():
        o.write(f"{r}\t{status[r]}\t{copy[r]}\t" + "\t".join(str(x) for x in f) + "\n")
from collections import Counter
print("wrote /tmp/dsfam42_features.tsv ; status:", dict(Counter(status[r] for r in rows)))
