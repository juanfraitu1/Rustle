#!/usr/bin/env python3
"""Amendment 43b: DNA labels of the flagged consensus sequences. Usage: varkmers_eval.py dev   (held-out values are not printed unless 'held' is given)."""
import collections
import json
import statistics as st
import sys

import numpy as np

sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/unmapped_rescue")
import varkmers as V  # noqa: E402

D2 = "/mnt/linuxdisk/tmp/dnaverify2"
O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"


def load():
    cnt = sum(np.fromfile(f"{D2}/chunk{i}.i32", dtype="<i4").astype(np.int64) for i in range(6))
    z = np.load(f"{D2}/var.npz")
    keys = list(z["keys"])
    feats = json.load(open(f"{O}/spec_features.json"))
    out = {}
    for i, k in enumerate(keys):
        a, p = cnt[z[f"all_{i}"]], cnt[z[f"pri_{i}"]]
        out[k] = dict(side=feats[k]["side"], verdict=feats[k]["verdict"], lab_all=V.label(a.tolist()), lab_pri=V.label(p.tolist()), n_all=len(a), n_pri=len(p),
                      f_all=float((a >= 2).mean()) if len(a) else None, f_pri=float((p >= 2).mean()) if len(p) else None,
                      med_all=float(np.median(a[a >= 2])) if (a >= 2).any() else None)
    return out


def main(side):
    lab = load()
    json.dump(lab, open(f"{O}/dna_labels.json", "w"))
    sel = {k: x for k, x in lab.items() if x["side"] == side}
    tr = [x for x in sel.values() if x["verdict"] == "TRUE" and x["n_pri"] >= 10]
    ok = sum(x["lab_pri"] == "DNA-SUPPORTED" for x in tr)
    print(f"[{side}] CALIBRATION: TRUE flags with >= 10 var_pri k-mers: {len(tr)}; DNA-SUPPORTED {ok} ({ok / max(1, len(tr)):.1%}) (bar >= 80%); "
          f"median share seen >= 2: {st.median(x['f_pri'] for x in tr):.3f}")
    for v in ("TRUE", "WRONG"):
        xs = [x for x in sel.values() if x["verdict"] == v]
        c = collections.Counter(x["lab_all"] for x in xs)
        meds = [x["med_all"] for x in xs if x["lab_all"] == "DNA-SUPPORTED" and x["med_all"] is not None]
        print(f"[{side}] {v:5s} flags {len(xs)}: var_all labels {dict(c)}; median share of var_all seen >= 2: "
              f"{st.median(x['f_all'] for x in xs if x['f_all'] is not None):.3f}; median WGS count of seen var_all k-mers (DNA-SUPPORTED): {st.median(meds) if meds else None}")


if __name__ == "__main__":
    main(sys.argv[1])
