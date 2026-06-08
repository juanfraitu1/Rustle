#!/usr/bin/env bash
# Install SQANTI3 into a dedicated conda env. SQANTI3 ships its own env yml.
source "$(dirname "$0")/config.sh"
set -euo pipefail

if [ ! -d "$SQANTI3_DIR" ]; then
  mkdir -p "$(dirname "$SQANTI3_DIR")"
  git clone https://github.com/ConesaLab/SQANTI3.git "$SQANTI3_DIR"
fi
cd "$SQANTI3_DIR"
git fetch --tags
# Pinned to v5.5.4 (validated for this protocol 2026-06-08). Fall back only if absent.
git checkout v5.5.4 2>/dev/null || git checkout v5.5.3 2>/dev/null || true

# The solve + R/bioconductor downloads + pip builds (parasail/edlib compile from
# source) take many minutes; the env is ~8 GB. mamba env create aborts at the
# overwrite prompt when a partial env exists and stdin is not a tty, so skip if
# present (remove a broken env manually with `mamba env remove -n $SQANTI3_ENV`).
if mamba env list | grep -qE "^\s*${SQANTI3_ENV}\s"; then
  echo "[SQANTI3] env '$SQANTI3_ENV' already exists; skipping create"
else
  ENV_YML=$(ls SQANTI3.conda_env.yml environment.yml 2>/dev/null | head -1)
  mamba env create -n "$SQANTI3_ENV" -f "$ENV_YML"
fi

mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" --version 2>&1 | tail -3 \
  || mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" --help 2>&1 | head -5
echo "[SQANTI3] installed at $SQANTI3_DIR, env $SQANTI3_ENV"
