#!/usr/bin/env python3
"""Non-coding-aware promotion of reference-absent copies: promote credible collapses on a
COLLAPSE-QUALITY bar (divergence band + genome coverage + own-locus + balanced co-segregation)
instead of a protein-BLASTx hit. Protein/ORF are LABELS, not gates. Purely additive: re-scores the
734 pre-ORF-gate consensuses already in gw_promoted/cons.fa and writes gw_noncoding_copies.json.
Every record is a FLAGGED reference-divergent candidate (copy_vs_allele='candidate-DNA-needed');
copy-vs-allele is not resolvable from RNA and needs DNA parCN. See
docs/superpowers/specs/2026-07-18-noncoding-promotion-design.md."""
import os
import sys
import json
import subprocess
from collections import Counter, defaultdict

sys.path.insert(0, "bench")
from copy_vs_allele_structural import load_sedef, sedef_partners  # stdlib-only module

# ---- paths (overridable in main via argparse) ----
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/gw_promoted"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)

# ---- the collapse-quality bar (single-edit recall/precision knob) ----
GID_LO, GID_HI = 0.60, 0.97   # genome_id band: <LO = repeat/chimera artifact, >=HI = ~reference allele
COV_MIN = 0.90                # consensus must map near full-length (not a chimeric fragment)
MIN_COLS = 8                  # co-segregation breadth (alt columns)
AF_LO, AF_HI = 0.15, 0.60     # overall alt-read fraction band (balanced collapse)
MIN_ALT = 5                   # minimum alt reads
MIN_DEPTH = 8                 # minimum primary reads (matches hidden_copy_scan.MIN_DEPTH)

# ---- honesty rail (set unconditionally on every promoted record) ----
CANDIDATE = "candidate-DNA-needed"
STATUS = "flagged-reference-divergent-candidate"


def classify_noncoding(rec):
    """Pure collapse-quality gate. rec keys: genome_id, genome_cov, own_locus (bool),
    alt_cols, alt_read_fraction, alt_reads, n_primary.
    Returns (promote: bool, call: str, reason: str). Only call=='noncoding-candidate' promotes."""
    gid = rec["genome_id"]
    cov = rec["genome_cov"]
    if gid >= GID_HI:
        return (False, "~REF", f"genome_id {gid:.3f} >= {GID_HI} (essentially the reference allele)")
    if gid < GID_LO:
        return (False, "artifact", f"genome_id {gid:.3f} < {GID_LO} (repeat/chimera, not a divergent copy)")
    if cov < COV_MIN:
        return (False, "chimera", f"genome_cov {cov:.2f} < {COV_MIN} (partial/chimeric consensus)")
    if not rec["own_locus"]:
        return (False, "not-own-locus", "best genome hit is not the candidate's own locus")
    if rec["alt_cols"] < MIN_COLS:
        return (False, "thin", f"alt_cols {rec['alt_cols']} < {MIN_COLS} (too few co-segregating columns)")
    if not (AF_LO <= rec["alt_read_fraction"] <= AF_HI):
        return (False, "thin", f"alt_read_fraction {rec['alt_read_fraction']:.3f} outside "
                               f"[{AF_LO},{AF_HI}] (unbalanced)")
    if rec["alt_reads"] < MIN_ALT or rec["n_primary"] < MIN_DEPTH:
        return (False, "thin", f"depth alt_reads={rec['alt_reads']} (min {MIN_ALT}) / "
                               f"n_primary={rec['n_primary']} (min {MIN_DEPTH})")
    return (True, "noncoding-candidate",
            f"divergent ({100*(1-gid):.1f}%) full-length (cov {cov:.2f}) balanced collapse "
            f"({rec['alt_cols']} cols, af {rec['alt_read_fraction']:.2f}) at own locus")
