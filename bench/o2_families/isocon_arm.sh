#!/bin/bash
# IsoCon on planted reads (register 1906): one pool per family, IsoCon pipeline, candidates aligned to the genome.
#   isocon_arm.sh PLANT_DIR SCENARIO GENOME.fa OUTDIR [CORES]     env ISOCON_DIR = directory with IsoCon's pixi.toml
# -> OUTDIR/isocon.gtf (primary alignments of final_candidates.fa, gene_id = family), OUTDIR/isocon.seq.tsv (per family:
#    candidates, equal to an isoform of an expressed copy, of silent copies only, of none)
set -uo pipefail
P=$1; SC=$2; FA=$3; O=$4; T=${5:-16}; H=$(cd "$(dirname "$0")" && pwd); mkdir -p "$O/pools" "$O/isocon" "$O/logs"
python3 - "$P" "$SC" "$O" <<'PY'
import sys, collections
p, sc, o = sys.argv[1:]
fam = {l.split('\t')[1]: l.split('\t')[0] for l in list(open(f'{p}/plan.tsv'))[1:]}
out = {}; name = None
for i, l in enumerate(open(f'{p}/{sc}.fq')):
    if i % 4 == 0: name = l[1:].strip()
    elif i % 4 == 1:
        f = fam[name.split('|')[1]]
        if f not in out: out[f] = open(f'{o}/pools/{f}.fa', 'w')
        out[f].write(f'>{name.replace("|", "_")}\n{l}')
PY
( cd "${ISOCON_DIR:?set ISOCON_DIR}" && ls "$O"/pools/*.fa | xargs -P 4 -I{} bash -c 't=$(basename {} .fa); pixi run IsoCon pipeline -fl_reads {} -outfolder '"$O"'/isocon/$t --nr_cores '"$T"' > '"$O"'/logs/$t.log 2>&1; echo "isocon $t exit $?"' )
cat /dev/null > "$O/cand.fa"; for f in "$O"/isocon/*/final_candidates.fa; do [ -s "$f" ] || continue; t=$(basename "$(dirname "$f")"); awk -v t="$t" '/^>/{print ">" t "__" substr($1,2); next}{print}' "$f" >> "$O/cand.fa"; done
minimap2 -ax splice:hq -uf --eqx -t "$T" "$FA" "$O/cand.fa" 2> "$O/mm.err" | samtools view -F 0x904 - > "$O/cand.sam"
python3 - "$P" "$SC" "$O" "$FA" "$H" <<'PY'
import sys, re, collections
p, sc, o, fa, h = sys.argv[1:]; sys.path.insert(0, h); import plant
with open(f'{o}/isocon.gtf', 'w') as w:
    for l in open(f'{o}/cand.sam'):
        f = l.split('\t'); pos = int(f[3]); ex = []; s = pos
        for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
            n = int(n)
            if op == 'N': ex.append((s, pos - 1)); pos += n; s = pos
            elif op in 'MD=X': pos += n
        ex.append((s, pos - 1)); st = '-' if int(f[1]) & 16 else '+'; fam = f[0].split('__')[0]
        at = f'gene_id "{fam}"; transcript_id "{f[0]}";'
        w.write(f'{f[2]}\tIsoCon\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{st}\t.\t{at}\n')
        for a, b in ex: w.write(f'{f[2]}\tIsoCon\texon\t{a}\t{b}\t.\t{st}\t.\t{at}\n')
FA = plant.Fai(fa); tx = collections.defaultdict(list); meta = {}
for l in open(f'{p}/loci.gtf'):
    f = l.rstrip('\n').split('\t')
    if f[2] == 'exon': t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tx[t].append((int(f[3]), int(f[4]))); meta[t] = (f[0], f[6])
hdr = None; expseq = collections.defaultdict(set); silseq = collections.defaultdict(set)
for l in open(f'{p}/plan.tsv'):
    f = l.rstrip('\n').split('\t')
    if hdr is None: hdr = f; continue
    v = f[hdr.index(sc)]
    if v == 'NA': continue
    s = ''.join(FA.fetch(meta[f[3]][0], a - 1, b) for a, b in sorted(tx[f[3]])).upper(); s = s[::-1].translate(plant.COMP) if meta[f[3]][1] == '-' else s
    (expseq if int(v) > 0 else silseq)[f[0]].add(s)
C = collections.defaultdict(collections.Counter); name = None
for l in open(f'{o}/cand.fa'):
    if l[0] == '>': name = l[1:].strip().split('__')[0]; continue
    s = l.strip().upper(); c = C[name]; c['candidates'] += 1
    c['expressed' if s in expseq[name] else 'silent_only' if s in silseq[name] else 'inside_expressed' if any(s in x for x in expseq[name]) else 'none'] += 1
with open(f'{o}/isocon.seq.tsv', 'w') as w:
    w.write('family\texpressed_sequences\tcandidates\tequal_expressed\tinside_expressed\tequal_silent_only\tnone\n')
    for f in sorted(set(expseq) | set(C)):
        c = C[f]; w.write(f"{f}\t{len(expseq[f])}\t{c['candidates']}\t{c['expressed']}\t{c['inside_expressed']}\t{c['silent_only']}\t{c['none']}\n")
print(open(f'{o}/isocon.seq.tsv').read(), end='')
PY
