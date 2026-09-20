#!/bin/bash
# DNA GRAPH vs SPLICED-RNA GRAPH, on the SAME NODE SET.
#
# The comparison is only interpretable if the nodes are identical, so both graphs are built over the
# 27 loci of F(NPIPB11). What differs is the SEQUENCE each node contributes to the edge rule:
#
#   DNA node = the locus's genomic sequence (introns included)
#   RNA node = the locus's READ-SUPPORTED EXONS (>=3 IsoSeq reads), concatenated in genomic order
#
# Same tiers (asm20 UNION sensitive -k11 -w5), same rule (identity >=0.80, coverage >=0.50 on min).
# So any difference between the two graphs is attributable to the SUBSTRATE, not to the node set,
# the aligner settings, or the threshold.
#
# ⚠ Uses the FULL A119b BAM (96 GB, 68M reads, no -L in its @PG chain). Both Soto caches are -M -L
#   subsets and would silently deflate the RNA side.
# ⚠ -F 2308 (primary only) for the per-read CIGAR walk. N in an RNA CIGAR is an intron spliced OUT.
# ⚠ RNA node sequence is much SHORTER than the DNA node; coverage denominators therefore differ. That
#   is the substrate effect under measurement, and both lengths are reported per node.
# ⚠ A locus with no read support yields an EMPTY RNA node -- reported explicitly, never dropped silently.
set -uo pipefail
OUT=/home/juanfra/winloci_scratch/seedfam
HERE="$(cd "$(dirname "$0")" && pwd)"
HS=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam
T=4
MIN_READS=3

echo "=== building RNA nodes (read-supported exons per locus) ==="
: > "$OUT/gvg_rna.fa"
: > "$OUT/gvg_nodelen.tsv"
while IFS=$'\t' read -r c s e rest; do
  reg="$c:$((s+1))-$e"
  samtools view -F 2308 "$BAM" "$reg" 2>/dev/null \
    | python3 "$HERE/read_exons.py" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4 -o count \
    | awk -v C="$c" -v S="$s" -v E="$e" -v M="$MIN_READS" -F'\t' \
        '$4>=M && $3>S && $2<E {a=($2>S?$2:S); b=($3<E?$3:E); if(b-a>=20) print C"\t"a"\t"b}' \
    > "$OUT/.gvg_ex.bed"
  bp=$(awk '{s+=$3-$2}END{print s+0}' "$OUT/.gvg_ex.bed")
  printf "%s\t%d\t%d\n" "$reg" "$((e-s))" "$bp" >> "$OUT/gvg_nodelen.tsv"
  if [ "$bp" -gt 0 ]; then
    bedtools getfasta -fi "$HS" -bed "$OUT/.gvg_ex.bed" -fo "$OUT/.gvg_ex.fa"
    python3 -c "
s=''.join(l.strip() for l in open('$OUT/.gvg_ex.fa') if l[0]!='>')
open('$OUT/gvg_rna.fa','a').write('>$reg\n'+s+'\n')"
  fi
done < "$OUT/HS.loci.bed"
echo "  RNA nodes with sequence: $(grep -c '^>' "$OUT/gvg_rna.fa")  of $(wc -l < "$OUT/HS.loci.bed")"

echo "=== all-vs-all, both substrates, identical settings ==="
for tag in dna rna; do
  src="$OUT/HS.loci.fa"; [ "$tag" = rna ] && src="$OUT/gvg_rna.fa"
# ⚠ TIER CORRECTION 2026-08-10 (defect B1): the all-vs-all E_r pass below now runs the SHIPPED
# command -- `minimap2 -c -X --no-long-join -t N -k 11 -w 5`, ONE tier. It previously ran
# `-c --eqx -N 200 -p 0.02` with NO `-X` plus an `-x asm20` leg that the shipped default SKIPS
# (RUSTLE_ER_SENSITIVE_ONLY has defaulted true since 2026-08-07). `-N`/`-p` are inert here; `-X`
# is the operative difference (it implies --dual=no). See bench/er_tier.sh for the full rationale.
# GENOME passes in this script are deliberately UNCHANGED: there -X would be wrong.
  minimap2 -c -X --no-long-join -t $T -k 11 -w 5 "$src" "$src" > "$OUT/gvg_${tag}.paf" 2>/dev/null
  echo "  $tag records: $(wc -l < "$OUT/gvg_${tag}.paf")"
done

python3 "$HERE/graph_vs_graph_report.py" "$OUT"
echo "DONE"
