import sys
from collections import defaultdict

# 1. Parse GGO_clusterd.gff to get transcript lengths and exon counts
ref_exons = defaultdict(int)
ref_chroms = {}
ref_strands = {}
ref_starts = {}
ref_ends = {}

with open('GGO_clusterd.gff') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        
        feature = parts[2]
        if feature == 'transcript':
            attrs = parts[8]
            t_id = [x for x in attrs.split(';') if 'transcript_id' in x][0].split('"')[1]
            ref_chroms[t_id] = parts[0]
            ref_strands[t_id] = parts[6]
            ref_starts[t_id] = int(parts[3])
            ref_ends[t_id] = int(parts[4])
        elif feature == 'exon':
            attrs = parts[8]
            t_id = [x for x in attrs.split(';') if 'transcript_id' in x][0].split('"')[1]
            ref_exons[t_id] += 1

# 2. Parse gffcmp.tracking to find missed transcripts
matched_refs = set()
with open('gffcmp.tracking') as f:
    for line in f:
        parts = line.strip().split('\t')
        class_code = parts[3]
        if class_code == '=':
            ref_id = parts[2].split('|')[1]
            matched_refs.add(ref_id)

# 3. Find missed multi-exon transcripts
missed_multi = []
for t_id, exons in ref_exons.items():
    if exons > 1 and t_id not in matched_refs:
        missed_multi.append({
            'id': t_id,
            'chrom': ref_chroms[t_id],
            'strand': ref_strands[t_id],
            'start': ref_starts[t_id],
            'end': ref_ends[t_id],
            'exons': exons
        })

# Sort by number of exons (descending) to find complex ones, and also print some simpler ones
missed_multi.sort(key=lambda x: x['exons'], reverse=True)

print(f"Total missed multi-exon transcripts: {len(missed_multi)}")
print("Top 5 complex missed transcripts:")
for m in missed_multi[:5]:
    print(f"ID: {m['id']}, Exons: {m['exons']}, Locus: {m['chrom']}:{m['start']}-{m['end']}({m['strand']})")

print("\nTop 5 simple (2-3 exons) missed transcripts:")
for m in [x for x in missed_multi if x['exons'] <= 3][:5]:
    print(f"ID: {m['id']}, Exons: {m['exons']}, Locus: {m['chrom']}:{m['start']}-{m['end']}({m['strand']})")
