#!/usr/bin/env python3
"""
Train a depth-3 decision tree on GGO_19 ML feature dump.

Usage:
    python3 bench/train_ml_filter.py \
        ml_training_features.jsonl \
        ml_training_all_transcripts.gtf \
        ml_training_cmp.ml_training_all_transcripts.gtf.tmap

Prints the trained tree as a Rust pub fn ml_predict body to stdout.
Paste the printed function into src/rustle/ml_filter.rs replacing the
placeholder ml_predict body.
"""

import sys
import json
import re
from pathlib import Path

import pandas as pd
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.model_selection import cross_val_score
import numpy as np


FEATURE_NAMES = [
    "cov_ratio", "longcov_ratio", "bpcov_ratio", "jct_mm_ratio",
    "nexons_k", "tlen_ratio", "longcov_abs", "is_intron_subset",
]


def parse_features(jsonl_path: str) -> dict:
    """Returns dict: intron_chain_str -> feature dict."""
    records = {}
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            chain = obj["intron_chain"]
            records[chain] = obj["features"]
    return records


def parse_gtf_intron_chains(gtf_path: str) -> dict:
    """Returns dict: tx_id -> intron_chain_str (1-indexed start, 0-indexed end format)."""
    tx_exons: dict = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" in line:
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            start, end = int(fields[3]), int(fields[4])
            attr = fields[8]
            m = re.search(r'transcript_id "([^"]+)"', attr)
            if not m:
                continue
            tx_id = m.group(1)
            tx_exons.setdefault(tx_id, []).append((start, end))

    chains = {}
    for tx_id, exons in tx_exons.items():
        exons.sort()
        if len(exons) < 2:
            chains[tx_id] = ""
            continue
        # intron: end_of_exon_i (GTF inclusive) + 1 for 1-based start,
        # start_of_next_exon (GTF 1-based = 0-based end of intron)
        parts = [f"{exons[i][1] + 1}-{exons[i+1][0]}" for i in range(len(exons) - 1)]
        chains[tx_id] = ",".join(parts)
    return chains


def parse_tmap(tmap_path: str) -> dict:
    """Returns dict: qry_id -> class_code."""
    cc = {}
    with open(tmap_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            # columns: ref_gene_id, ref_id, class_code, qry_gene_id, qry_id, ...
            cc[parts[4]] = parts[2]
    return cc


def tree_to_rust(clf: DecisionTreeClassifier, feature_names: list, cv_auc: float) -> str:
    """Convert a fitted DecisionTreeClassifier to a Rust if-else function body."""
    tree = clf.tree_

    def node_to_rust(node: int, indent: int) -> str:
        pad = "    " * indent
        if tree.children_left[node] == tree.children_right[node]:
            # Leaf: majority class
            values = tree.value[node][0]
            n_kill, n_keep = int(values[0]), int(values[1])
            keep = n_keep >= n_kill
            verdict = "true" if keep else "false"
            return f"{pad}return {verdict}; // keep={n_keep}, kill={n_kill}\n"

        feat = feature_names[tree.feature[node]]
        thresh = tree.threshold[node]
        left = node_to_rust(tree.children_left[node], indent + 1)
        right = node_to_rust(tree.children_right[node], indent + 1)
        return (
            f"{pad}if f.{feat} <= {thresh:.6f} {{\n"
            f"{left}"
            f"{pad}}}\n"
            f"{right}"
        )

    body = node_to_rust(0, 1)
    date = pd.Timestamp.now().strftime("%Y-%m-%d")
    return (
        f"/// Returns `true` → keep transcript, `false` → kill (longunder).\n"
        f"/// Trained on GGO_19 de novo, {date}, 5-fold CV AUC = {cv_auc:.3f}\n"
        f"pub fn ml_predict(f: &MlFeatures) -> bool {{\n"
        f"{body}"
        f"}}\n"
    )


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    jsonl_path, gtf_path, tmap_path = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"[train] Loading features from {jsonl_path}...", file=sys.stderr)
    feature_records = parse_features(jsonl_path)
    print(f"[train] {len(feature_records)} feature records loaded", file=sys.stderr)

    print(f"[train] Parsing intron chains from {gtf_path}...", file=sys.stderr)
    tx_chains = parse_gtf_intron_chains(gtf_path)

    print(f"[train] Parsing labels from {tmap_path}...", file=sys.stderr)
    tx_class = parse_tmap(tmap_path)

    # Build chain -> class_code mapping via tx_id join
    chain_to_class: dict = {}
    for tx_id, chain in tx_chains.items():
        if chain and tx_id in tx_class:
            chain_to_class[chain] = tx_class[tx_id]

    # Build training rows: join feature_records (keyed by chain) with labels
    rows = []
    for chain, feats in feature_records.items():
        if chain not in chain_to_class:
            continue  # no label — skip
        label = 1 if chain_to_class[chain] == "=" else 0
        row = {k: feats[k] for k in FEATURE_NAMES}
        row["label"] = label
        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"[train] {len(df)} labeled rows ({df['label'].sum():.0f} keep, {(1-df['label']).sum():.0f} kill)",
          file=sys.stderr)

    if len(df) < 20:
        print("[train] ERROR: fewer than 20 labeled examples — check paths or re-run dump.", file=sys.stderr)
        sys.exit(1)

    X = df[FEATURE_NAMES].values
    y = df["label"].values

    clf = DecisionTreeClassifier(max_depth=3, class_weight="balanced", random_state=42)
    cv_scores = cross_val_score(clf, X, y, cv=5, scoring="roc_auc")
    print(f"[train] 5-fold CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}", file=sys.stderr)

    clf.fit(X, y)

    text = export_text(clf, feature_names=FEATURE_NAMES)
    print("[train] Sklearn tree structure:", file=sys.stderr)
    print(text, file=sys.stderr)

    rust_fn = tree_to_rust(clf, FEATURE_NAMES, cv_scores.mean())
    print(rust_fn)


if __name__ == "__main__":
    main()
