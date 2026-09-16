#!/bin/bash
# usage: chrom_score.sh CHROM LABEL GTF [--sqanti]
# gffcompare (+ optional SQANTI3) scoring driver for a human ordinary-chromosome assembler bakeoff
# (parameterized sibling of bench/chr20_score.sh). Writes gffcompare's own version and, when --sqanti is
# given, the SQANTI3 checkout's commit into $W/gffcompare/<label>.provenance.txt (hardening for
# docs/PREREG_gtf_refine_chr17_2026-09-16.md, "Run protocol (frozen)": tool-version drift between labels
# would confound E1-E5 the same way an unpinned env var would).
set -euo pipefail
CHROM=${1:?CHROM}; LABEL=${2:?LABEL}; GTF=${3:?GTF}; SQ_FLAG=${4:-}
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
SQANTI3_DIR=/mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3
mkdir -p "$W/gffcompare"; cd "$W/gffcompare"

{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "gffcompare_version: $(gffcompare --version 2>&1 | head -1)"
} > "$LABEL.provenance.txt"

# gffcompare, no flags beyond -r/-o (no -R/-Q -- those change precision and would confound E2). Output goes
# to files only (log + .stats) -- no scored metric (Query mRNAs / Transcript level / Intron chain level /
# etc.) is printed here, per the PREREG's no-peeking rule (docs/PREREG_gtf_refine_chr17_2026-09-16.md,
# "Failure policy"): nothing from .stats/.tmap/classification/junctions may be inspected before
# bench/gtf_refine_verdict.py computes the verdict.
gffcompare -r "$W/${CHROM}_ref.gtf" -o "$LABEL" "$GTF" > "$LABEL.gffcompare.log" 2>&1

if [ "$SQ_FLAG" = "--sqanti" ]; then
  echo "sqanti3_commit: $(git -C "$SQANTI3_DIR" rev-parse HEAD 2>/dev/null || echo unknown)" >> "$LABEL.provenance.txt"
  source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3
  mkdir -p "$W/sqanti3/$LABEL"; cd "$SQANTI3_DIR"
  python sqanti3_qc.py --isoforms "$GTF" --refGTF "$W/${CHROM}_ref.gtf" --refFasta "$W/${CHROM}.fa" \
    -o "$LABEL" -d "$W/sqanti3/$LABEL" --report skip -t 4 > "$W/sqanti3/$LABEL.qc.log" 2>&1
  echo "SQANTI3 $LABEL exit=$?"
fi
