#!/usr/bin/env python3
"""U6 CLI: combine per-arm recovery.tsv + authenticity.tsv -> headtohead.tsv + REPORT.md."""
import argparse, csv, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


def _load_recovery(path):
    rec, fam = {}, {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rec[r["ref_transcript"]] = {"fsm": r["fsm"] == "true", "ism": r["ism"] == "true"}
            fam[r["ref_transcript"]] = r["family_id"]
    return rec, fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rustle-recovery", required=True)
    ap.add_argument("--stringtie-recovery", required=True)
    ap.add_argument("--authenticity", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    rustle_rec, fam = _load_recovery(args.rustle_recovery)
    stringtie_rec, fam2 = _load_recovery(args.stringtie_recovery)
    fam.update(fam2)
    authentic = {}
    with open(args.authenticity) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            authentic[r["ref_transcript"]] = r["authentic"] == "true"

    res = lib_eval.head_to_head(rustle_rec, stringtie_rec, authentic, fam)

    h2h = os.path.join(args.outdir, "headtohead.tsv")
    with open(h2h, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["category", "ref_transcript", "family_id"])
        for tx in res["rustle_only_fsm_authentic"]:
            w.writerow(["rustle_only_fsm_authentic", tx, fam.get(tx, "NA")])
        for tx in res["rustle_only_fsm_phantom"]:
            w.writerow(["rustle_only_fsm_phantom", tx, fam.get(tx, "NA")])

    with open(os.path.join(args.outdir, "REPORT.md"), "w") as o:
        o.write("# VG copy-recovery head-to-head\n\n")
        o.write(f"- Authentic rustle-VG-only FSM recoveries (StringTie2 misses): **{res['n_win']}**\n")
        o.write(f"- Phantom (secondary-only) rustle-VG-only FSM recoveries: **{res['n_phantom']}**\n")
        o.write(f"- Families won: {', '.join(res['families_won']) or 'none'}\n\n")
        o.write("See `headtohead.tsv` for per-transcript rows.\n")
    sys.stderr.write(f"[U6] win={res['n_win']} phantom={res['n_phantom']} -> {h2h}\n")


if __name__ == "__main__":
    main()
