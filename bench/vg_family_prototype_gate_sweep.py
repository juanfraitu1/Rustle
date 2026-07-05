#!/usr/bin/env python3
"""vg_family_prototype_gate_sweep.py — post-hoc sweep empirical FP gates.

Uses the feature table from vg_family_prototype_fp_characterize.py to filter
vg_family_prototype.tsv and runs vg_family_prototype_eval.py on each filtered
catalog.  Reports the P/R trade-off of each gate.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_gate_sweep.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import subprocess
from collections import defaultdict

BENCH = os.path.dirname(os.path.abspath(__file__))
CATALOG_TSV = os.path.join(BENCH, "vg_family_prototype.tsv")
FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_gate_sweep.json")
OUT_TXT = os.path.join(BENCH, "vg_family_prototype_gate_sweep.txt")
PYTHON = "/home/juanfra/miniforge3/bin/python"
EVAL_PY = os.path.join(BENCH, "vg_family_prototype_eval.py")


def load_features(path):
    rows = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows[int(row["family_id"])] = {k: (float(v) if "." in v else int(v)) for k, v in row.items()}
    return rows


def load_catalog(path):
    fams = defaultdict(list)
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            fams[int(row["fam_id"])].append(row["member"])
    return fams


def write_catalog(fams, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["fam_id", "member"])
        for fid in sorted(fams):
            for m in fams[fid]:
                w.writerow([fid, m])


def run_eval(catalog_path):
    base, _ = os.path.splitext(catalog_path)
    prefix = base + "_eval"
    json_path = prefix + ".json"
    cmd = [PYTHON, EVAL_PY, "--input", catalog_path, "--output-prefix", prefix]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"[!] eval failed for {catalog_path}", flush=True)
        print(proc.stderr[-2000:], flush=True)
        return None
    with open(json_path) as fh:
        return json.load(fh)


def main():
    features = load_features(FEATURES_TSV)
    catalog = load_catalog(CATALOG_TSV)

    # candidate gates (OR semantics across conditions)
    gates = [
        ("baseline", []),
        ("members>=80", [("n_members", 80)]),
        ("members>=70", [("n_members", 70)]),
        ("members>=60", [("n_members", 60)]),
        ("pair_hub_frac>=0.05", [("pair_hub_frac", 0.05)]),
        ("repeat_hub_frac>=0.05", [("repeat_hub_frac", 0.05)]),
        ("members>=80_OR_pair_hub_frac>=0.05", [("n_members", 80), ("pair_hub_frac", 0.05)]),
        ("members>=80_OR_repeat_hub_frac>=0.05", [("n_members", 80), ("repeat_hub_frac", 0.05)]),
        ("mean_pair_mult>=25", [("mean_pair_mult", 25)]),
        ("pair_hub_frac>=0.4", [("pair_hub_frac", 0.4)]),
        ("members>=60_AND_mean_pair_mult>=25", [("n_members", 60), ("mean_pair_mult", 25)]),
        ("members>=50_AND_mean_pair_mult>=20", [("n_members", 50), ("mean_pair_mult", 20)]),
        ("mean_pair_mult>=20_AND_pair_hub_frac>=0.3", [("mean_pair_mult", 20), ("pair_hub_frac", 0.3)]),
        ("members>=40_AND_mean_pair_mult>=18_AND_pair_hub_frac>=0.2", [("n_members", 40), ("mean_pair_mult", 18), ("pair_hub_frac", 0.2)]),
    ]

    results = []
    for name, conds in gates:
        print(f"[*] gate: {name}", flush=True)
        kept = {}
        dropped = 0
        for fid, members in catalog.items():
            feat = features[fid]
            drop = any(feat[k] >= v for k, v in conds)
            if drop:
                dropped += 1
            else:
                kept[fid] = members

        cat_path = os.path.join(BENCH, f"vg_family_prototype_gate_{name}.tsv")
        write_catalog(kept, cat_path)
        print(f"    kept {len(kept)} families, dropped {dropped}", flush=True)

        ev = run_eval(cat_path)
        if ev is None:
            continue

        ep = ev["truth1_protein_Ep"]
        dna = ev["truth2_dna_loose"]
        orc = ev["truth3_diploid_oracle"]
        res = dict(
            gate=name,
            conditions=str(conds),
            n_kept=len(kept),
            n_dropped=dropped,
            n_multicopy_families=ev["n_multicopy_families"],
            precision_Ep=ep["precision_Ep"],
            impure_blocks=ep["impure_blocks"],
            pure_blocks=ep["pure_blocks"],
            recall_oracle=orc["recall_oracle"],
            precision_oracle_task=orc["precision_oracle_taskformula"],
            precision_oracle_dedup=orc["precision_oracle_dedup"],
            n_allele=orc["n_allele"],
            n_oversize=orc["n_oversize"],
            n_multifam=orc["n_multifam"],
            distinct_fp_blocks=orc["distinct_fp_blocks"],
            component_recovery_recall=dna["component_recovery_recall"],
            pair_projection_recall_real_cdna=dna["pair_projection_recall_real_cdna"],
        )
        results.append(res)
        print(f"    P_Ep={res['precision_Ep']}  R_oracle={res['recall_oracle']}  "
              f"P_oracle(dedup)={res['precision_oracle_dedup']}  FP={res['n_allele']}/{res['n_oversize']}/{res['n_multifam']}", flush=True)

    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=False)
    print(f"[write] {OUT_JSON}", flush=True)

    lines = []
    lines.append("Gate sweep summary")
    def fmt(x):
        return f"{x:>7.4f}" if x is not None else "    —  "

    lines.append(f"{'gate':<45} {'kept':>6} {'drop':>5} {'multi':>6} {'P_Ep':>7} {'R_orac':>7} {'P_orac_dedup':>12} {'FP a/o/m':>12}")
    for r in results:
        lines.append(f"{r['gate']:<45} {r['n_kept']:>6} {r['n_dropped']:>5} {r['n_multicopy_families']:>6} "
                     f"{fmt(r['precision_Ep'])} {fmt(r['recall_oracle'])} {fmt(r['precision_oracle_dedup'])} "
                     f"{r['n_allele']}/{r['n_oversize']}/{r['n_multifam']}")
    txt = "\n".join(lines)
    print(txt)
    with open(OUT_TXT, "w") as fh:
        fh.write(txt + "\n")
    print(f"[write] {OUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
