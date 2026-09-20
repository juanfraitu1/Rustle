#!/usr/bin/env bash
# Run the full per-copy attribution pipeline for ONE identity level:
#   generate labeled synthetic tandem array -> align -> rustle --vg+tandem -> score.
#
# Usage: run_one.sh <identity> <out_dir> [seed] [copies] [reads_per_copy]
# Env:   RUSTLE (default target/release/rustle), MINIMAP2 (default minimap2)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

ID="${1:?usage: run_one.sh <identity> <out_dir> [seed] [copies] [reads_per_copy] [spacing]}"
OUT="${2:?need out_dir}"
SEED="${3:-1}"
COPIES="${4:-3}"
RPC="${5:-40}"
# Default to DISPERSED loci (200 kb apart): each copy forms its own bundle and
# normal VG discovery links them into ONE family via cross-copy multimappers
# (the DAZ-like path) — robust and free of tandem-bundling fragility. This tests
# the SAME fingerprint-EM that tandem mode feeds. Use a small spacing (~9500) to
# test the tandem-collapse path instead.
SPACING="${6:-200000}"

RUSTLE="${RUSTLE:-$REPO/target/release/rustle}"
MINIMAP2="${MINIMAP2:-minimap2}"
mkdir -p "$OUT"

# 1. generate labeled synthetic array
python3 "$HERE/gen_synthetic.py" --identity "$ID" --copies "$COPIES" \
    --reads-per-copy "$RPC" --out "$OUT" --seed "$SEED" --spacing "$SPACING" >&2

# 2. align spliced reads to the array; KEEP secondaries (near-identical copies
#    must multimap so the EM has cross-copy evidence to apportion).
#    -ax splice:hq  high-quality spliced; -uf  forward transcript strand;
#    -N 20 -p 0.5   retain up to 20 secondaries down to 50% of the primary score.
"$MINIMAP2" -ax splice:hq -uf -N 20 -p 0.5 --secondary=yes \
    "$OUT/genome.fa" "$OUT/reads.fa" 2>/dev/null \
    | samtools sort -o "$OUT/aligned.bam" - 2>/dev/null
samtools index "$OUT/aligned.bam"

# 3. rustle --vg + tandem, dumping per-read fingerprint-EM attribution
RUSTLE_VG_TANDEM=1 RUSTLE_VG_ANCHOR_PRIOR=1 \
RUSTLE_VG_FP_ATTR_TSV="$OUT/fp.tsv" \
    "$RUSTLE" -L "$OUT/aligned.bam" --vg --genome-fasta "$OUT/genome.fa" \
    -o "$OUT/out.gtf" 2>"$OUT/rustle.log"

# tandem fired?
grep -q "copy sub-bundles" "$OUT/rustle.log" && TANDEM_FIRED=1 || TANDEM_FIRED=0
echo "[run_one] identity=$ID tandem_fired=$TANDEM_FIRED fp_rows=$(wc -l <"$OUT/fp.tsv" 2>/dev/null || echo 0)" >&2

# 4. score attribution
if [[ -s "$OUT/fp.tsv" ]]; then
    python3 "$HERE/score_attribution.py" --tsv "$OUT/fp.tsv" \
        --reads "$OUT/reads.fa" --meta "$OUT/meta.json" --json "$OUT/score.json"
else
    echo "{\"error\":\"no fp.tsv — EM did not run (tandem_fired=$TANDEM_FIRED)\",\"identity_target\":$ID}" \
        | tee "$OUT/score.json"
fi
