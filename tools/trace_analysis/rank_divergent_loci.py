# Rank loci by Rustle/StringTie divergence
import re
from collections import defaultdict

# Per-chrom per-strand intervals: StringTie transcripts
st_tx = defaultdict(list)
with open('GGO_19_stringtie.gtf') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip().split('\t')
        if len(p)<9 or p[2]!='transcript': continue
        st_tx[(p[0], p[6])].append((int(p[3]), int(p[4])))

# Rustle transcripts
ru_tx = defaultdict(list)
with open('/tmp/v.gtf') as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip().split('\t')
        if len(p)<9 or p[2]!='transcript': continue
        ru_tx[(p[0], p[6])].append((int(p[3]), int(p[4])))

# For each StringTie locus (merged overlapping tx), count ST and RU tx fully or >80% inside
def merge(intervals):
    if not intervals: return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for s,e in intervals[1:]:
        if s <= out[-1][1]+500:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s,e])
    return [tuple(x) for x in out]

loci_info = []
for key, tx in st_tx.items():
    chrom, strand = key
    loci = merge(tx)
    ru_list = ru_tx.get(key, [])
    for (ls,le) in loci:
        st_count = sum(1 for s,e in tx if s>=ls and e<=le)
        ru_count = sum(1 for s,e in ru_list if min(e,le)-max(s,ls) > 0.8*(e-s))
        loci_info.append((chrom, strand, ls, le, st_count, ru_count, abs(st_count-ru_count), le-ls))

loci_info.sort(key=lambda x: (-x[6], -x[4]))
print("chrom\tstrand\tstart\tend\tst_tx\tru_tx\tabs_diff\tlen")
for x in loci_info[:20]:
    print('\t'.join(str(v) for v in x))
