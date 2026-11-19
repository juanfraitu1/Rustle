#!/bin/bash
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --job-name=r1829_sel
# PREREG r1829 layer 2, part A: HG002 residues 0-24, tier ms10 (as layer 1), first 30,000 by sha256 order; all their records
source /storage/group/kdm16/default/jxi21/fritz/env.sh
W=$FRITZ/ws/r1829_near_tie_census/layer2; mkdir -p $W; cd $W
B=$HOME/reads/HSA/HG002_kinnex.t2t.bam
pixi run -q --manifest-path $FRITZ/Rustle/pixi.toml python - "$B.as_table.tsv" names.txt <<'EOF'
import hashlib, sys
sel = []
for l in open(sys.argv[1]):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[4] in ('NA', '*') or f[11] == 'NA': continue
    if f[14] == f[4] and abs(int(f[16]) - int(f[15])) < 5000: continue
    if abs(int(f[8]) - int(f[12])) > 10: continue
    h = hashlib.sha256(f[0].encode()).hexdigest()
    if int(h[:8], 16) % 100 < 25: sel.append((h, f[0]))
sel.sort()
print(f"ms10 molecules in residues 0-24: {len(sel)}; taken {min(len(sel), 30000)}")
open(sys.argv[2], 'w').write(''.join(n + '\n' for _, n in sel[:30000]))
EOF
samtools view -@ 32 -b -N names.txt -o sel.bam $B && echo "extract rc $? $SECONDS s"
samtools view -c sel.bam
