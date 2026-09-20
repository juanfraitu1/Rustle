import sys
from collections import defaultdict

# Parse gffcmp.tracking to find novel transcripts (class codes other than =)
novel_classes = defaultdict(int)
total_query = 0

with open('gffcmp.tracking') as f:
    for line in f:
        total_query += 1
        parts = line.strip().split('\t')
        class_code = parts[3]
        novel_classes[class_code] += 1

print(f"Total query transcripts: {total_query}")
print("Class codes distribution:")
for code, count in sorted(novel_classes.items(), key=lambda x: x[1], reverse=True):
    print(f"  {code}: {count} ({(count/total_query)*100:.1f}%)")
