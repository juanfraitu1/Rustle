#!/usr/bin/env python3
"""U4 CLI: read a SQANTI3 *_classification.txt + universe.tsv -> <arm>_recovery.tsv."""
import argparse, csv, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classification", required=True)
    ap.add_argument("--universe", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.universe) as f:
        uni_rows = list(csv.DictReader(f, delimiter="\t"))
    universe_tx = {r["transcript_id"] for r in uni_rows}
    fam_of = {r["transcript_id"]: r["family_id"] for r in uni_rows}

    with open(args.classification) as f:
        classif = list(csv.DictReader(f, delimiter="\t"))
    rec = lib_eval.recovery_set(classif, universe_tx)

    with open(args.out, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["ref_transcript", "family_id", "fsm", "ism"])
        for tx in sorted(rec):
            w.writerow([tx, fam_of.get(tx, "NA"),
                        "true" if rec[tx]["fsm"] else "false",
                        "true" if rec[tx]["ism"] else "false"])


if __name__ == "__main__":
    main()
