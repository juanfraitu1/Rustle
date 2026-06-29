#!/usr/bin/env bash
# Per-locus VG-vs-StringTie evaluation against RefSeq truth.
#
# For one candidate paralog gene, extract its region (+ its dominant sibling's
# region) from the full genome BAM, run StringTie / rustle-baseline / rustle-VG
# (win stack), score each against the RefSeq GFF, and emit a single JSON object
# describing which real RefSeq transcripts each tool matched and what VG gained
# over StringTie / baseline.
#
# Usage: winloci_eval.sh GENE CHROM START END STRAND [SIB_GENE]
# Env:   GI=<gene_introns.tsv>  PAD=60000  VGFLAGS="..."  (override VG env)
# Emits: one JSON object on stdout (everything else to stderr).
set -uo pipefail
GENE="${1:?gene}"; CHROM="${2:?chrom}"; START="${3:?start}"; END="${4:?end}"
STRAND="${5:?strand}"; SIB="${6:-}"
ROOT=/mnt/c/Users/jfris/Desktop/Rustle
SCRATCH="${SCRATCH:-/home/juanfra/winloci_scratch}"
GGO="${GGO:-$SCRATCH/GGO.bam}"
GFF="${GFF:-$SCRATCH/GGO_tx.gff}"
FASTA="${FASTA:-$SCRATCH/GGO.fasta}"
RUSTLE="$ROOT/target/release/rustle"
ST="$ROOT/tools/stringtie/stringtie"
GI="${GI:-$ROOT/bench/paralog_secondary_scan/scan_out/gene_introns.tsv}"
W="/tmp/winloci/$GENE"; rm -rf "$W"; mkdir -p "$W"

# Window policy: candidate gene +/- FLANK; merge with sibling. Same-chrom siblings
# within MERGE_DIST collapse into ONE CONTIGUOUS region (spans the copies between
# them -- fixes tandem arrays); far same-chrom or cross-chrom siblings -> two
# regions. (PAD kept as an alias for FLANK for back-compat.)
FLANK="${FLANK:-${PAD:-15000}}"; MERGE_DIST="${MERGE_DIST:-200000}"
cs=$(( START>FLANK ? START-FLANK : 1 )); ce=$(( END+FLANK ))
REGIONS="$CHROM:$cs-$ce"
if [ -n "$SIB" ] && [ -f "$GI" ]; then
  line=$(awk -F'\t' -v g="$SIB" '$1==g{print $2"\t"$3"\t"$4; exit}' "$GI")
  if [ -n "$line" ]; then
    sc=$(echo "$line"|cut -f1); ss=$(echo "$line"|cut -f2); se=$(echo "$line"|cut -f3)
    gap=$(( START>se ? START-se : (ss>END ? ss-END : 0) ))
    if [ "$sc" = "$CHROM" ] && [ "$gap" -le "$MERGE_DIST" ]; then
      lo=$(( START<ss ? START : ss )); hi=$(( END>se ? END : se ))
      lo=$(( lo>FLANK ? lo-FLANK : 1 )); hi=$(( hi+FLANK ))
      REGIONS="$CHROM:$lo-$hi"
    else
      ss2=$(( ss>FLANK ? ss-FLANK : 1 )); se2=$(( se+FLANK ))
      REGIONS="$CHROM:$cs-$ce $sc:$ss2-$se2"
    fi
  fi
fi
REGSTR="${REGIONS// /,}"
echo "[eval] $GENE regions='$REGIONS' sib=$SIB" >&2

samtools view -b "$GGO" $REGIONS 2>/dev/null | samtools sort -o "$W/locus.bam" - 2>/dev/null
samtools index "$W/locus.bam" 2>/dev/null || true
NREADS=$(samtools view -c -F 0x900 "$W/locus.bam" 2>/dev/null || echo 0)
NSEC=$(samtools view -c -f 0x100 "$W/locus.bam" 2>/dev/null || echo 0)
echo "[eval] $GENE reads(primary)=$NREADS secondary=$NSEC" >&2

"$ST" -L "$W/locus.bam" -o "$W/st.gtf" >/dev/null 2>&1 || true
"$RUSTLE" -L "$W/locus.bam" -o "$W/base.gtf" >/dev/null 2>&1 || true
env ${VGFLAGS:-RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1} \
  "$RUSTLE" --vg --vg-snp --genome-fasta "$FASTA" -L "$W/locus.bam" -o "$W/vg.gtf" >/dev/null 2>&1 || true

eqset() { local f; f=$(ls $1 2>/dev/null | head -1); [ -n "$f" ] && awk -F'\t' 'NR>1 && $3=="="{print $2}' "$f" | sort -u; }
for t in st base vg; do
  if [ -s "$W/$t.gtf" ]; then gffcompare -r "$GFF" "$W/$t.gtf" -o "$W/$t" >/dev/null 2>&1 || true; fi
done
eqset "$W/st.*.tmap"   > "$W/st.eq"   || true
eqset "$W/base.*.tmap" > "$W/base.eq" || true
eqset "$W/vg.*.tmap"   > "$W/vg.eq"   || true

python3 - "$GENE" "$REGSTR" "$STRAND" "$NREADS" "$NSEC" "$W" <<'PY'
import sys, json, re
gene, region, strand, nreads, nsec, W = sys.argv[1:7]
def load(p):
    try: return [x.strip() for x in open(p) if x.strip()]
    except FileNotFoundError: return []
st  = set(load(f"{W}/st.eq"))
base= set(load(f"{W}/base.eq"))
vg  = set(load(f"{W}/vg.eq"))
def txcount(p):
    try: return sum(1 for l in open(p) if "\ttranscript\t" in l)
    except FileNotFoundError: return 0
st_tx, base_tx, vg_tx = txcount(f"{W}/st.gtf"), txcount(f"{W}/base.gtf"), txcount(f"{W}/vg.gtf")
gain_st   = sorted(vg - st)
gain_base = sorted(vg - base)
loss_st   = sorted(st - vg)
if gain_st:
    cls = "win_vs_st"
elif gain_base:
    cls = "win_vs_base_only"
elif loss_st:
    cls = "regress"
else:
    cls = "tie"
print(json.dumps(dict(
    gene=gene, region=region, strand=strand,
    n_primary_reads=int(nreads), n_secondary=int(nsec),
    st_tx=st_tx, base_tx=base_tx, vg_tx=vg_tx,
    st_eq=sorted(st), base_eq=sorted(base), vg_eq=sorted(vg),
    vg_gain_vs_st=gain_st, vg_gain_vs_base=gain_base, vg_loss_vs_st=loss_st,
    classification=cls)))
PY
