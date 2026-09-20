#!/bin/bash
# gffcompare + SQANTI3 scoring driver for the human chr20 ordinary-chromosome assembler bakeoff
# (bench/CHR20_ASSEMBLER_COMPARISON.md). Runs AFTER bakeoff_chr20_ours.sh, bakeoff_chr20_stringtie.sh,
# bakeoff_chr20_flair.sh have each produced their GTF. Scores all three tools against the same chr20
# reference GTF (extracted+converted from chm13v2.0_RefSeq_full.gff.gz, see prep_chr20_ref.sh) and, for
# our own tool only, additionally scores the multicopy="false"-only subset of its GTF (the `multicopy`
# attribute is emitted by `--gtf`; see src/bin/copy_assign.rs:193-198).
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
REF_GTF="$W/chr20_ref.gtf"
REF_FA="$W/chr20.fa"
OURS_GTF="$W/ours/ours.gtf"
ST_GTF="$W/stringtie/st.gtf"
FLAIR_GTF="$W/flair/flair.isoforms.gtf"

# ---------------------------------------------------------------------------------------------------
# 1. gffcompare, all three tools, whole chr20.
# ---------------------------------------------------------------------------------------------------
mkdir -p "$W/gffcompare"
cd "$W/gffcompare"
gffcompare -r "$REF_GTF" -o ours "$OURS_GTF" > ours.gffcompare.log 2>&1
gffcompare -r "$REF_GTF" -o stringtie "$ST_GTF" > stringtie.gffcompare.log 2>&1
gffcompare -r "$REF_GTF" -o flair "$FLAIR_GTF" > flair.gffcompare.log 2>&1
echo "=== gffcompare .stats (whole chr20) ==="
for f in ours stringtie flair; do echo "--- $f ---"; cat "$f.stats"; done

# ---------------------------------------------------------------------------------------------------
# 2. Our tool's multicopy="false"-only subset (substitutes for "--non-hard-loci"; no new flag needed --
#    reuses the existing --gtf `multicopy` attribute). On this substrate copy_assign's own de novo family
#    detector found 0 co-located families on chr20 (see ours/ours.stderr.log: "0 co-located families"), so
#    EVERY one of the 976 transcripts is multicopy="false" -- the subset is byte-identical to the whole
#    "ours" GTF. We still materialize and score it explicitly rather than assert this, for the record.
# ---------------------------------------------------------------------------------------------------
python3 - "$OURS_GTF" "$W/ours/ours.multicopy_false.gtf" << 'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
keep_tids = set()
lines = open(src).readlines()
for l in lines:
    f = l.rstrip("\n").split("\t")
    if len(f) > 2 and f[2] == "transcript" and 'multicopy "false"' in f[8]:
        for attr in f[8].split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                keep_tids.add(attr.split('"')[1])
n_written = 0
with open(dst, "w") as out:
    for l in lines:
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] in ("transcript", "exon"):
            tid = None
            for attr in f[8].split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split('"')[1]
                    break
            if tid in keep_tids:
                out.write(l)
                n_written += 1
print(f"[chr20_score] multicopy=false subset: {len(keep_tids)} transcripts, {n_written} GTF rows -> {dst}")
PYEOF

gffcompare -r "$REF_GTF" -o ours_multicopy_false "$W/ours/ours.multicopy_false.gtf" > ours_multicopy_false.gffcompare.log 2>&1
echo "=== gffcompare .stats (ours, multicopy=false only) ==="
cat ours_multicopy_false.stats

# ---------------------------------------------------------------------------------------------------
# 3. SQANTI3 QC, all three tools, whole chr20. --report skip: only the *_classification.txt structural
#    category table is needed; --report html/pdf pulls in an R rendering path we don't need for this.
# ---------------------------------------------------------------------------------------------------
source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3
SQ=/mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3
cd "$SQ"
for pair in "ours:$OURS_GTF" "stringtie:$ST_GTF" "flair:$FLAIR_GTF"; do
  name="${pair%%:*}"; gtf="${pair#*:}"
  mkdir -p "$W/sqanti3/$name"
  python sqanti3_qc.py --isoforms "$gtf" --refGTF "$REF_GTF" --refFasta "$REF_FA" \
    -o "$name" -d "$W/sqanti3/$name" --report skip -t 4 \
    > "$W/sqanti3/$name.qc.log" 2>&1
  echo "SQANTI3 $name exit=$?"
done

echo "=== SQANTI3 structural category counts ==="
for name in ours stringtie flair; do
  f="$W/sqanti3/$name/${name}_classification.txt"
  echo "--- $name (total $(tail -n +2 "$f" | wc -l)) ---"
  tail -n +2 "$f" | cut -f6 | sort | uniq -c | sort -rn
done
