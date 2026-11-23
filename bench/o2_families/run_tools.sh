#!/bin/bash
# O2 family benchmark: align planted reads and run every assembler on the same BAM, then score against the planted truth.
#   run_tools.sh PLANT_DIR SCENARIO GENOME.fa OUTDIR [THREADS]
# Needs on PATH: minimap2, samtools, python3 (pysam), stringtie; Rustle binaries in RUSTLE_BIN (default ../../target/release).
# Optional arms, each run only when its variable is set (the paths are site-specific):
#   O2_FLAIR_PATH=DIR[:DIR]   directories holding flair and flair_bam2bed.py
#   O2_ISOQUANT_MANIFEST=pixi.toml   (pixi run --manifest-path ... isoquant)
#   O2_BAMBU_MANIFEST=pixi.toml      (pixi run --manifest-path ... Rscript bench/famcmp/bambu_run.R)
# MAP_P (default 0.1; 0.8 = the slim mapping mode of row 1902). Arms that fail are reported and scored as MISSING.
set -uo pipefail
P=$1; SC=$2; FA=$3; O=$4; T=${5:-32}
H=$(cd "$(dirname "$0")" && pwd); R=$(cd "$H/../.." && pwd); BIN=${RUSTLE_BIN:-$R/target/release}; DR=$R/tools/rustle_pipeline.sh
mkdir -p "$O" && cd "$O" || exit 1
say() { echo "[o2_families] $SC $*"; }
if [ ! -s x.bam ]; then
  minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p "${MAP_P:-0.1}" --secondary=yes -t "$T" "$FA" "$P/$SC.fq" 2> mm.err | samtools sort -@ 4 -o x.bam - && samtools index x.bam || { say "align FAILED"; exit 2; }
fi
"$BIN/as_table" --bam x.bam --out x.bam.as_table.tsv --threads 8 > as_table.log 2>&1 || say "as_table FAILED"
UNPOOL=${UNPOOL:-0} python3 "$H/tied_share.py" "$P/loci.gtf" x.bam x.bam.as_table.tsv "$P/$SC.truth.tsv" reads.txt 8 > reads.log 2>&1 || say "read provenance FAILED"
run() { local n=$1 s=$SECONDS; shift; "$@" > "$n.log" 2>&1; say "$n exit $? in $((SECONDS - s)) s"; }
# O2_RUSTLE_ARMS: which Rustle arms run (default all three; the mode needs the secondaries arm before it)
ARMS=" ${O2_RUSTLE_ARMS:-primary good mode} "
[[ $ARMS != *" primary "* ]] || run rustle_primary bash "$DR" assemble --seed-pool primary --bam "$O/x.bam" --fasta "$FA" --out "$O/rp" --bin "$BIN" --threads "$T"
[[ $ARMS != *" good "* ]] || run rustle_good bash "$DR" assemble --seed-pool good --bam "$O/x.bam" --fasta "$FA" --out "$O/rg" --bin "$BIN" --threads "$T"
[[ $ARMS != *" mode "* ]] || RUSTLE_MULTICOPY_EQ=0 run rustle_mode bash "$DR" multicopy --seed-pool good --bam "$O/x.bam" --fasta "$FA" --out "$O/rg" --bin "$BIN" --threads "$T"
# the O2 extension: the mode's assembly, a locus kept only when the EM leaves it >= 2 reads (register 1905 / 1906)
EMW=; [ "${O2_EM_WEIGHT:-1}" = 0 ] || EMW="--weighted ${O2_EM_WEIGHT:-1}"   # scores inside the EM, as stage o2 since 2026-11-23 (register 1921); O2_EM_WEIGHT=0 = counts-only
if [[ $ARMS == *" mode "* ]] && [ -s "$O/rg.multicopy.gtf" ]; then run rustle_mode_em python3 "$H/em_filter.py" "$O/rg.multicopy.gtf" "$O/x.bam" "$O/x.bam.as_table.tsv" "$O/rg.o2" --threads 8 $EMW; fi
if [[ $ARMS == *" primary "* ]] && [ -s "$O/rp.gtf" ]; then run rustle_primary_em python3 "$H/em_filter.py" "$O/rp.gtf" "$O/x.bam" "$O/x.bam.as_table.tsv" "$O/rp.o2" --threads 8 $EMW; fi
# opt-in arm `mono` (row 1919): the mode with single-exon models rescued on family evidence, in a copy of the secondaries run
if [[ $ARMS == *" mono "* ]] && [ -s "$O/rg.gtf" ]; then
  for x in gtf families.gtf molecules.tsv molecules.tsv.asbin; do [ ! -e "$O/rg.$x" ] || cp -p "$O/rg.$x" "$O/rgm.$x"; done
  RUSTLE_MULTICOPY_MONO=1 RUSTLE_MULTICOPY_EQ=0 run rustle_mode_mono bash "$DR" multicopy --seed-pool good --bam "$O/x.bam" --fasta "$FA" --out "$O/rgm" --bin "$BIN" --threads "$T"
  [ ! -s "$O/rgm.multicopy.gtf" ] || run rustle_mode_mono_em python3 "$H/em_filter.py" "$O/rgm.multicopy.gtf" "$O/x.bam" "$O/x.bam.as_table.tsv" "$O/rgm.o2" --threads 8 $EMW
fi
if command -v stringtie > /dev/null; then run stringtie stringtie -L -p 8 -o "$O/stringtie.gtf" "$O/x.bam"; else say "stringtie not on PATH: arm skipped"; fi
if [ -n "${O2_ISOQUANT_MANIFEST:-}" ]; then run isoquant pixi run --manifest-path "$O2_ISOQUANT_MANIFEST" isoquant -d pacbio_ccs --bam "$O/x.bam" --reference "$FA" -o "$O/iq" -p x -t 16; fi
if [ -n "${O2_FLAIR_PATH:-}" ]; then
  flair_arm() { export PATH="$O2_FLAIR_PATH:$PATH"; local tmp; tmp=$(mktemp -d "${TMPDIR:-/tmp}/fl.XXXXXX")
    samtools fastq -@ 4 -F 0x900 -0 "$tmp/reads.fastq" x.bam && flair_bam2bed.py -b x.bam -o "$O/flair" \
      && flair collapse -q "$O/flair.bed" -g "$FA" -r "$tmp/reads.fastq" -o "$O/fl" -t 16 --trust_ends --temp_dir "$tmp"; local rc=$?; rm -rf "$tmp"; return $rc; }
  run flair flair_arm
fi
if [ -n "${O2_BAMBU_MANIFEST:-}" ]; then run bambu pixi run --manifest-path "$O2_BAMBU_MANIFEST" Rscript "$R/bench/famcmp/bambu_run.R" "$O/x.bam" "$FA" "$O/bambu" 16 "${TMPDIR:-/tmp}/rc.$$"; fi
LABELS=(rustle_primary="$O/rp.gtf" primary+em="$O/rp.o2.gtf" rustle_secondaries="$O/rg.gtf" rustle_mode="$O/rg.multicopy.gtf" mode+em="$O/rg.o2.gtf" stringtie="$O/stringtie.gtf"
  flair="$O/fl.isoforms.gtf" isoquant="$O/iq/x/x.transcript_models.gtf" bambu="$O/bambu.gtf")
[[ $ARMS != *" mono "* ]] || LABELS+=(mode_mono="$O/rgm.multicopy.gtf" mode_mono+em="$O/rgm.o2.gtf")
python3 "$H/score.py" "$P" "$SC" --families --out "$O/score.tsv" "${LABELS[@]}" > /dev/null
# units of identical copies and divergence classes (register 1906); built once per planted set
[ -s "$P/units.tsv" ] || python3 "$H/units.py" "$P" "$FA" 8 > "$P/units.log" 2>&1 || say "units FAILED (see $P/units.log)"
[ ! -s "$P/units.tsv" ] || python3 "$H/units_score.py" "$P" "$SC" --out "$O/units.tsv" "${LABELS[@]}" > /dev/null || say "unit scoring FAILED"
grep -P "\tALL\t|^scenario" "$O/score.tsv"; [ ! -s "$O/units.tsv" ] || cut -f1-8 "$O/units.tsv"; cat reads.txt 2>/dev/null
