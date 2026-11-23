#!/bin/bash
# IsoCon on the reads of Rustle's tied groups (register 1925): no annotation and no truth decide the pools.
#   isocon_groups.sh ASSEMBLY.gtf BAM AS_TABLE READS.fq GENOME.fa OUTDIR [CORES]     env ISOCON_DIR = directory with IsoCon's pixi.toml
# A group = loci joined by reads exactly tied between them on ms (sample_psv.tie_groups). Its tied reads make one pool
# (pools under 5 reads are skipped); IsoCon derives the distinct sequences of the pool; they are aligned to the genome.
# -> OUTDIR/groups.tsv (pool, loci, reads), OUTDIR/cand.fa, OUTDIR/cand.sam (primary alignments; name = POOL__candidate)
set -uo pipefail
GTF=$1; BAM=$2; TAB=$3; FQ=$4; FA=$5; O=$6; T=${7:-16}; H=$(cd "$(dirname "$0")" && pwd); mkdir -p "$O/pools" "$O/isocon" "$O/logs"
python3 - "$H" "$GTF" "$BAM" "$TAB" "$FQ" "$O" "$T" <<'PY' || exit 1
import sys
h, gtf, bam, tab, fq, o, t = sys.argv[1:]; sys.path.insert(0, h)
import sample_psv as S
meta, sets, groups, gr = S.tie_groups(gtf, bam, tab, int(t))
pool = {}
with open(f'{o}/groups.tsv', 'w') as w:
    w.write('pool\tloci\treads\n')
    for i, k in enumerate(sorted(groups, key=lambda k: (meta[k][0], meta[k][1])), 1):
        n = f'g{i:04d}'; w.write(f"{n}\t{','.join(sorted(groups[k]))}\t{len(gr[k])}\n")
        if len(gr[k]) >= 5:
            for q in gr[k]: pool[q] = n
out = {}; name = None
for i, l in enumerate(open(fq)):
    if i % 4 == 0: name = l[1:].split()[0]
    elif i % 4 == 1 and name in pool:
        f = pool[name]
        if f not in out: out[f] = open(f'{o}/pools/{f}.fa', 'w')
        out[f].write(f'>{name.replace("|", "_")}\n{l}')
for f in out.values(): f.close()
print(f'isocon_groups: {len(groups)} groups, {len(out)} pools of >= 5 tied reads, {len(pool)} reads')
PY
( cd "${ISOCON_DIR:?set ISOCON_DIR}" && ls "$O"/pools/*.fa | xargs -P 4 -I{} bash -c 't=$(basename {} .fa); pixi run IsoCon pipeline -fl_reads {} -outfolder '"$O"'/isocon/$t --nr_cores '"$T"' > '"$O"'/logs/$t.log 2>&1; echo "isocon $t exit $?"' )
cat /dev/null > "$O/cand.fa"; for f in "$O"/isocon/*/final_candidates.fa; do [ -s "$f" ] || continue; t=$(basename "$(dirname "$f")"); awk -v t="$t" '/^>/{print ">" t "__" substr($1,2); next}{print}' "$f" >> "$O/cand.fa"; done
minimap2 -ax splice:hq -uf --eqx -t "$T" "$FA" "$O/cand.fa" 2> "$O/mm.err" | samtools view -F 0x904 - > "$O/cand.sam"
echo "isocon_groups: $(grep -c '^>' "$O/cand.fa") candidates, $(wc -l < "$O/cand.sam") aligned"
