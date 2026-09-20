#!/bin/bash
# FLAIR arm of the human chr20 ordinary-chromosome assembler bakeoff (bench/CHR20_ASSEMBLER_COMPARISON.md).
# Sibling of bakeoff_flair.sh (gorilla NPIP substrate) -- DO NOT overwrite that one.
# Same pipeline: BAM -> FASTQ via `samtools fastq -F 2308` -> flair align -> flair correct (unguided,
# no --gtf) -> flair collapse. FLAIR at its documented defaults.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
BAM="$W/chr20.bam"
FA="$W/chr20.fa"
source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate flair
mkdir -p "$W/flair"; cd "$W/flair"

# PATH-ORDERING WORKAROUND (confirmed root cause, not a guess): flair's own `flair_cli.py:main()` calls
# `flair.set_unix_path()`, which prepends the flair PACKAGE directory itself
# (site-packages/flair/, where filter_transcriptome_align.py physically lives) to PATH on every
# invocation. That script's shebang is `#!/usr/bin/env python3`; the package dir has no python3 binary of
# its own, so `env` keeps searching the REST of PATH -- and on this machine /home/linuxbrew/.linuxbrew/bin
# (no `flair` package) sits ahead of the flair conda env's own bin in the inherited shell PATH, so
# `filter_transcriptome_align.py` was executing under the WRONG python3 and dying with
# "ModuleNotFoundError: No module named 'flair'" even though the flair package is fully installed.
# Forcing the flair env's own bin/ to the front of PATH here (still ahead of linuxbrew after flair
# prepends its package dir on top) fixes `env python3` resolution without touching the conda install.
export PATH="/home/juanfra/miniforge3/envs/flair/bin:$PATH"

# 1. reads. flair aligns from FASTQ, so the molecules must leave the BAM.
#    -F 2308 keeps primary/mapped/non-supplementary only, so all tools start from the SAME molecules.
if [ ! -s reads.fq ]; then
  samtools fastq -F 2308 -@ 4 "$BAM" > reads.fq 2> fastq.log
fi
echo "reads: $(( $(wc -l < reads.fq) / 4 ))"

# 2. flair align -- NOTE minimap2 runs with --secondary=no inside flair (flair_align.py:150).
[ -s flair.bed ] || flair align -g "$FA" -r reads.fq -o flair --threads 4 2>&1 | tail -5

# 3. flair correct -- SKIPPED. Confirmed empirically (`flair correct --help` + a live run) that the
#    installed FLAIR 3.0.0 hard-requires -f/--gtf, --junction_tab, or --junction_bed for this subcommand
#    ("FlairInputDataError: No junctions from GTF or junctionsBed to correct with. Exiting...") -- fully
#    unguided splice correction, which the original bakeoff_flair.sh relied on, no longer exists in this
#    version. Supplying our chr20 reference GTF here would make FLAIR's run *guided* while ours/StringTie
#    stay unguided, breaking the apples-to-apples comparison the task wants -- so instead we feed the RAW
#    `flair align` bed straight into `flair collapse`, staying unguided at the cost of skipping the
#    splice-site-correction step. Documented as a deviation in bench/CHR20_ASSEMBLER_COMPARISON.md.

# 4. flair collapse -- emits the isoform GTF that the common scorer reads.
[ -s flair.isoforms.gtf ] || flair collapse -g "$FA" -q flair.bed -r reads.fq \
     -o flair --threads 4 --generate_map 2>&1 | tail -5

ls -la "$W/flair"
