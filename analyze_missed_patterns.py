import sys
from collections import defaultdict

# 1. Parse Reference GFF to get multi-exon status and metadata
ref_exons = defaultdict(int)
ref_info = {}

with open('GGO_clusterd.gff') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        
        feature = parts[2]
        if feature == 'transcript':
            attrs = parts[8]
            t_id = [x for x in attrs.split(';') if 'transcript_id' in x][0].split('"')[1]
            ref_info[t_id] = {'chrom': parts[0], 'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6]}
        elif feature == 'exon':
            attrs = parts[8]
            t_id = [x for x in attrs.split(';') if 'transcript_id' in x][0].split('"')[1]
            ref_exons[t_id] += 1

# 2. Parse gffcmp.tracking to find MATCHED multi-exon reference transcripts
matched_refs = set()
with open('gffcmp.tracking') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 4: continue
        class_code = parts[3]
        if class_code == '=':
            ref_id_full = parts[2]
            if '|' in ref_id_full:
                ref_id = ref_id_full.split('|')[1]
                matched_refs.add(ref_id)

# 3. Find MISSED multi-exon reference transcripts
missed_multi = []
for t_id, count in ref_exons.items():
    if count > 1 and t_id not in matched_refs:
        missed_multi.append(t_id)

# 4. Analyze class codes for these missed transcripts
# We need to see how they WERE matched if not with '='
missed_ref_to_best_code = {}
with open('gffcmp.tracking') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 4: continue
        ref_id_full = parts[2]
        if '|' in ref_id_full:
            ref_id = ref_id_full.split('|')[1]
            if ref_id in missed_multi:
                code = parts[3]
                # Prioritize 'near' matches over total misses
                current_best = missed_ref_to_best_code.get(ref_id, 'none')
                priority = {'=': 10, 'c': 8, 'j': 7, 'k': 6, 'm': 5, 'n': 4, 'none': 0}
                if priority.get(code, 1) > priority.get(current_best, 0):
                    missed_ref_to_best_code[ref_id] = code

code_stats = defaultdict(int)
for ref_id in missed_multi:
    code = missed_ref_to_best_code.get(ref_id, 'unmatched')
    code_stats[code] += 1

print(f"Total Multi-Exon Reference Transcripts: {len([t for t, c in ref_exons.items() if c > 1])}")
print(f"Total Missed Multi-Exon: {len(missed_multi)}")
print("\nMissed Multi-Exon Class Code Distribution (Best match found):")
for code, count in sorted(code_stats.items(), key=lambda x: x[1], reverse=True):
    print(f"  {code}: {count} ({(count/len(missed_multi))*100:.1f}%)")

print("\nSample of 'unmatched' (completely missing) transcripts:")
unmatched_ids = [r for r in missed_multi if missed_ref_to_best_code.get(r) is None or missed_ref_to_best_code.get(r) == 'unmatched']
for uid in unmatched_ids[:5]:
    info = ref_info[uid]
    print(f"  {uid} | {info['chrom']}:{info['start']}-{info['end']}({info['strand']}) | Exons: {ref_exons[uid]}")
