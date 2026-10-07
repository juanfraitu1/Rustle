#!/bin/bash
# Shim copied from the lab's benchmark_collapse/bin/collapse_isoforms_precise.py (flair invokes it by filename, the package ships no entry point); only the interpreter path is changed to the local flair environment.
exec /home/juanfra/miniforge3/envs/flair/bin/python -m flair.collapse_isoforms_precise "$@"
