#!/usr/bin/env python3
"""Validation E — rustle paralog counts vs Ensembl Compara gorilla truth.

Pulls Ensembl Compara homology TSVs for both gorilla (paralog counts) and
human (paralog counts + human→gorilla ortholog map), then for a panel of
canonical gene symbols produces a side-by-side count table.

See README.md in this directory for method, results, and caveats.

Usage:
    python3 run_validation.py
    -> reads gorilla.compara.tsv.gz, human.compara.tsv.gz
    -> writes compara_panel_validation.tsv

The two .tsv.gz files must be present in the same directory.  Download with:
    curl -O https://ftp.ensembl.org/pub/release-112/tsv/ensembl-compara/homologies/gorilla_gorilla/Compara.112.protein_default.homologies.tsv.gz
    mv Compara.112.protein_default.homologies.tsv.gz gorilla.compara.tsv.gz
    curl -O https://ftp.ensembl.org/pub/release-112/tsv/ensembl-compara/homologies/homo_sapiens/Compara.112.protein_default.homologies.tsv.gz
    mv Compara.112.protein_default.homologies.tsv.gz human.compara.tsv.gz
"""

from __future__ import annotations
import csv
import gzip
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
HUMAN_FILE   = HERE / "human.compara.tsv.gz"
GORILLA_FILE = HERE / "gorilla.compara.tsv.gz"
OUT_TSV      = HERE / "compara_panel_validation.tsv"

# Canonical human Ensembl gene IDs for the 9-family panel. Stable IDs that
# don't change across releases. Verify on the Ensembl gene page if updating.
HUMAN_ENSG = {
    "GOLGA6L1": "ENSG00000174450",
    "GOLGA8A":  "ENSG00000175265",
    "AMY1A":    "ENSG00000237763",
    "TBC1D3":   "ENSG00000274226",
    "NBPF1":    "ENSG00000219481",
    "OR1A1":    "ENSG00000172146",
    "MUC2":     "ENSG00000198788",
    "HBB":      "ENSG00000244734",
    "ZNF91":    "ENSG00000167232",
}

# Rustle counts per slide 28's multi_copy_eval numbers (gorilla data,
# expression-detected paralogs).  Re-derive from rustle's own
# family-discovery TSV when repeating this for paper artifacts.
RUSTLE_COUNTS = {
    "GOLGA6L1":  6, "GOLGA8A": 13, "AMY1A":   3,
    "TBC1D3":   11, "NBPF1":  23, "OR1A1":  26,
    "MUC2":     15, "HBB":     2, "ZNF91":  15,
}

# Compara TSV columns (release 112 schema):
COL_GENE      = 0   # gene_stable_id
COL_HOMOLOGY  = 4   # homology_type
COL_HOM_GENE  = 5   # homology_gene_stable_id
COL_HOM_SPEC  = 7   # homology_species


def index_human(targets: set[str]) -> tuple[dict, dict, dict]:
    """One pass through human.compara.tsv.gz.

    Returns three dicts keyed by ENSG:
      within_paralogs[ENSG]    -> set of paralog ENSG IDs (strict)
      other_paralogs[ENSG]     -> set of paralog ENSG IDs (other_paralog)
      gorilla_orthologs[ENSG]  -> set of gorilla ENSGGOG IDs
    """
    within  = {g: set() for g in targets}
    other   = {g: set() for g in targets}
    gorilla = {g: set() for g in targets}
    with gzip.open(HUMAN_FILE, "rt") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[COL_GENE] not in targets:
                continue
            spec = c[COL_HOM_SPEC]
            hom  = c[COL_HOMOLOGY]
            if spec == "homo_sapiens":
                if hom == "within_species_paralog":
                    within[c[COL_GENE]].add(c[COL_HOM_GENE])
                elif hom == "other_paralog":
                    other[c[COL_GENE]].add(c[COL_HOM_GENE])
            elif spec == "gorilla_gorilla" and "ortholog" in hom:
                gorilla[c[COL_GENE]].add(c[COL_HOM_GENE])
    return within, other, gorilla


def index_gorilla(targets: set[str]) -> tuple[dict, dict]:
    within = {g: set() for g in targets}
    other  = {g: set() for g in targets}
    with gzip.open(GORILLA_FILE, "rt") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[COL_GENE] not in targets:
                continue
            if c[COL_HOM_SPEC] != "gorilla_gorilla":
                continue
            hom = c[COL_HOMOLOGY]
            if hom == "within_species_paralog":
                within[c[COL_GENE]].add(c[COL_HOM_GENE])
            elif hom == "other_paralog":
                other[c[COL_GENE]].add(c[COL_HOM_GENE])
    return within, other


def main() -> None:
    if not HUMAN_FILE.exists() or not GORILLA_FILE.exists():
        raise SystemExit(
            "Missing Compara TSVs.  See module docstring for the curl commands."
        )
    targets = set(HUMAN_ENSG.values())
    h_within, h_other, gorilla_orthologs = index_human(targets)

    # Collect gorilla targets (orthologs of our human panel)
    g_targets: set[str] = set()
    for orth_set in gorilla_orthologs.values():
        g_targets.update(orth_set)
    g_within, g_other = index_gorilla(g_targets)

    print(
        f"{'symbol':<10}{'rustle':>8}"
        f"{'h_within':>10}{'h_broad':>9}"
        f"{'g_within':>10}{'g_broad':>9}"
        f"{'n_orth':>8}  notes"
    )
    print("-" * 92)

    rows = []
    for sym, hensg in HUMAN_ENSG.items():
        hw = len(h_within[hensg]) + 1
        hb = len(h_within[hensg] | h_other[hensg]) + 1
        orth = gorilla_orthologs[hensg]
        if orth:
            gw = max(len(g_within[o]) for o in orth) + 1
            gb = max(len(g_within[o] | g_other[o]) for o in orth) + 1
        else:
            gw = gb = 0
        note = (
            "no gorilla ortholog" if not orth
            else (f"{len(orth)} orthologs" if len(orth) > 1 else "")
        )
        rt = RUSTLE_COUNTS[sym]
        print(
            f"{sym:<10}{rt:>8}"
            f"{hw:>10}{hb:>9}"
            f"{gw:>10}{gb:>9}"
            f"{len(orth):>8}  {note}"
        )
        rows.append((sym, hensg, rt, hw, hb, gw, gb, len(orth), note))

    with open(OUT_TSV, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sym", "human_ensg", "rustle",
            "human_within", "human_broad",
            "gorilla_within", "gorilla_broad",
            "n_gorilla_orthologs", "note",
        ])
        for row in rows:
            w.writerow(row)
    print(f"\n✓ wrote {OUT_TSV}")


if __name__ == "__main__":
    main()
