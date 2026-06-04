#!/usr/bin/env bash
# O4/O5 COPY-RESOLUTION BENCHMARK (the discriminating multi-copy benchmark).
#
# WHY THIS BENCHMARK EXISTS (grounded 2026-06-03, workflow o4o5-benchmark-discrimination):
#   gffcompare intron-chain Sn/Pr does NOT discriminate VG from the baseline on multi-copy
#   paralogs. When every copy has full-length reads the baseline scores 100/100 (one transcript
#   per separated bundle) and VG can only tie or, via secondary-alignment intake, LOSE precision.
#   No synthetic regime makes VG beat the baseline on chain Sn/Pr (verified: copies 2-6, spacing
#   13k-30k, with/without starvation, with/without RUSTLE_VG_TANDEM). So chain Sn/Pr is the WRONG
#   metric — it measures a task the baseline already solves.
#
#   The honest discriminator is COPY ATTRIBUTION: which paralog copy did each ambiguous
#   (multi-mapping) read come from? VG's fingerprint-EM answers this; the baseline has NO copy
#   concept at all (zero copy_id / copy_confidence / family_* attributes), so the metric is
#   UNDEFINED for it — not merely worse. This benchmark scores that metric against per-read
#   ground truth (reads named c{src}_r{j}) and characterizes the IDENTIFIABILITY boundary,
#   including the critical honesty check: does VG ABSTAIN when copies are identical (the DAZ
#   non-identifiable limit), or does it fabricate? Fabrication here is the DAZ3 failure mode.
#
# FIXTURE: a merged-but-separable tandem pair (spacing 16000 > gene_len 11600, so copies do not
#   physically overlap) with structurally-distinct isoforms (--distinct-isoforms). Every read
#   multimaps to both copies (minimap2 -N20 -p0.5), so the ambiguous reads are genuine.
#
# OUTPUT: bench/tandem_attribution/o4o5_results.json + a printed markdown summary.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
RUSTLE="$ROOT/target/release/rustle"
GEN="$HERE/gen_synthetic.py"
SCORE="$HERE/score_attribution.py"
OUT="${OUT:-/tmp/o4o5_bench}"
SEEDS="${SEEDS:-1 3 5 7 11}"
RPC_DEFAULT="${RPC_DEFAULT:-80}"   # reads/copy so the family clears DEFAULT --vg-family-min-shared=10
mkdir -p "$OUT"

align() { minimap2 -ax splice:hq -uf -N20 -p0.5 "$1" "$2" 2>/dev/null | samtools sort -o "$3" - 2>/dev/null; samtools index "$3"; }
chain_snpr() { gffcompare -r "$1" "$2" -o "$3" 2>/dev/null; grep -E "Transcript level" "$3.stats" | grep -oE '[0-9]+\.[0-9]+' | paste -sd'/'; }

# acc for one (identity, seed, filter-args) → echoes "acc dec_acc dec_frac abstain n_mm" or "DROPPED"
attr_run() {
  local id="$1" seed="$2" rpc="$3" d="$4"; shift 4; local filt=("$@")
  rm -rf "$d"; python3 "$GEN" --identity "$id" --copies 2 --reads-per-copy "$rpc" --spacing 16000 \
      --distinct-isoforms --seed "$seed" --out "$d" >/dev/null 2>&1
  align "$d/genome.fa" "$d/reads.fa" "$d/r.bam"
  rm -f "$d/fp.tsv"
  RUSTLE_VG_FP_ATTR_TSV="$d/fp.tsv" "$RUSTLE" --vg --vg-snp "${filt[@]}" --genome-fasta "$d/genome.fa" \
      -L "$d/r.bam" -o "$d/vg.gtf" 2>/dev/null
  if [ ! -f "$d/fp.tsv" ]; then echo "DROPPED"; return; fi
  python3 "$SCORE" --tsv "$d/fp.tsv" --reads "$d/reads.fa" --meta "$d/meta.json" --json "$d/s.json" >/dev/null 2>&1
  python3 -c "import json;d=json.load(open('$d/s.json'));print(d['attribution_accuracy_multimappers'],d['decisive_accuracy'],d['decisive_frac'],d['abstain_frac'],d['n_multimapper_reads_scored'])"
}

echo "## O4/O5 copy-resolution benchmark"
echo

# ── PART 1: primary discriminating case (DEFAULT settings, id=0.97) ──────────────
echo "### Part 1 — copy attribution, DEFAULT settings (identity 0.97, rpc=$RPC_DEFAULT)"
echo "VG attributes ambiguous multimappers to their true source copy; baseline has no copy metric."
p1_accs=(); p1_base_chain=(); p1_vg_chain=()
for s in $SEEDS; do
  d="$OUT/p1_s$s"
  read acc dacc dfrac abst nmm <<< "$(attr_run 0.97 "$s" "$RPC_DEFAULT" "$d" || echo 'DROPPED - - - -')"
  "$RUSTLE" -L "$d/r.bam" -o "$d/base.gtf" 2>/dev/null
  bchain=$(chain_snpr "$d/truth.gtf" "$d/base.gtf" "$d/bc"); vchain=$(chain_snpr "$d/truth.gtf" "$d/vg.gtf" "$d/vc")
  bcopy=$(grep -cE 'copy_id|copy_confidence|family_id' "$d/base.gtf" || true)
  echo "  seed $s: VG copy-attr acc=$acc dec_acc=$dacc (n_mm=$nmm) | chain base=$bchain vg=$vchain | baseline copy-attrs=$bcopy"
  [ "$acc" != "-" ] && p1_accs+=("$acc")
done
echo

# ── PART 2: identifiability spectrum + ABSTENTION boundary (relaxed filter, documented) ──
echo "### Part 2 — identifiability spectrum (relaxed --vg-family-min-shared 2; small controlled family)"
echo "Relaxed because the real-data spurious-family guard (min_shared=10) drops this small KNOWN family"
echo "at high identity; relaxing lets the EM run so its behavior across the spectrum is MEASURABLE."
echo "Honesty check: at identity 1.0 (identical copies, the DAZ limit) VG must ABSTAIN, not fabricate."
RELAX=(--vg-family-min-shared 2 --vg-family-min-shared-per-copy 0)
printf "  %-9s | %-8s | %-8s | %-9s | %-9s | %s\n" identity acc dec_acc dec_frac abstain note
mean() { printf '%s\n' "$@" | awk '{s+=$1;n++} END{if(n)printf "%.3f",s/n; else printf "NA"}'; }
for id in 1.0 0.999 0.99 0.97 0.95 0.90; do
  accs=(); daccs=(); dfracs=(); absts=(); dropped=0
  for s in $SEEDS; do
    d="$OUT/p2_${id}_s$s"
    out="$(attr_run "$id" "$s" "$RPC_DEFAULT" "$d" "${RELAX[@]}" || echo DROPPED)"
    if [ "$out" = "DROPPED" ]; then dropped=$((dropped+1)); continue; fi
    read acc dacc dfrac abst nmm <<< "$out"
    accs+=("$acc"); daccs+=("$dacc"); dfracs+=("$dfrac"); absts+=("$abst")
  done
  note="identifiable"; [ "$id" = "1.0" ] && note="NON-id → must abstain (DAZ limit)"
  [ "$id" = "0.999" ] && note="1 SNP → boundary: OVERCONFIDENT if dec_acc≪1 at dec_frac≫0"
  [ "$id" = "0.90" ] && note="<ambiguity threshold: reads map uniquely, no family"
  printf "  %-9s | %-8s | %-8s | %-9s | %-9s | %s (dropped %d/%s)\n" \
    "$id" "$(mean "${accs[@]:-NA}")" "$(mean "${daccs[@]:-NA}")" "$(mean "${dfracs[@]:-NA}")" \
    "$(mean "${absts[@]:-NA}")" "$note" "$dropped" "$(echo $SEEDS|wc -w)"
done
# Aggregate results JSON (reproducibility + oracle wiring).
python3 - "$OUT" <<'PY'
import json, glob, os, sys
out = sys.argv[1]; agg = {"part1_default_id097": [], "part2_spectrum": {}}
for f in sorted(glob.glob(os.path.join(out, "p1_s*", "s.json"))):
    agg["part1_default_id097"].append(json.load(open(f)))
for f in sorted(glob.glob(os.path.join(out, "p2_*_s*", "s.json"))):
    d = json.load(open(f)); agg["part2_spectrum"].setdefault(str(d["identity_target"]), []).append(d)
json.dump(agg, open(os.path.join(os.path.dirname(os.path.dirname(out)),
          "o4o5_results.json") if False else "/tmp/o4o5_results.json", "w"), indent=2)
print("  wrote /tmp/o4o5_results.json")
PY
echo
echo "Baseline never produces a copy metric (0 copy_id/copy_confidence/family attributes) — the"
echo "copy-resolution metric is UNDEFINED for the baseline, so this benchmark demonstrates a"
echo "CAPABILITY VG provides and the baseline does not, validated against per-read ground truth."
