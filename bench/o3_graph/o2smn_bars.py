#!/usr/bin/env python3
"""PREREG_o2_smn_showcase_2026-10-20: bars R1-R4 from o2loci.py score's score.json.   o2smn_bars.py SCORE.json OUT.json"""
import json
import sys

COPIES = ("S19SMN2.0", "S19SMN2.1")


def inc(n):
    k = n.get("FL", 0) + n.get("D7", 0)
    return n.get("FL", 0) / k if k else None


def bars(r):
    out = {}
    for c in COPIES:
        t = inc(r["truth"]["reads"][c])
        d = lambda a, k: abs(inc(r[a][k][c]) - t)
        n_o2, n_t = sum(r["o2"]["reads"][c].values()), sum(r["truth"]["reads"][c].values())
        b = dict(truth=round(t, 4), d_reads=dict(o2=round(d("o2", "reads"), 4), mm2=round(d("mm2", "reads"), 4)),
                 d_loci=dict(rustle_o2=round(d("o2", "rustle"), 4), rustle_mm2=round(d("mm2", "rustle"), 4),
                             stringtie_mm2=round(d("mm2", "stringtie"), 4), flair_mm2=round(d("mm2", "flair"), 4)),
                 rustle_o2_vs_truth=round(abs(inc(r["o2"]["rustle"][c]) - inc(r["truth"]["rustle"][c])), 4),
                 spanning=dict(o2=n_o2, truth=n_t))
        b["R1"] = b["d_reads"]["o2"] < b["d_reads"]["mm2"]
        b["R2"] = all(b["d_loci"]["rustle_o2"] < b["d_loci"][k] for k in ("rustle_mm2", "stringtie_mm2", "flair_mm2"))
        b["R3"] = b["rustle_o2_vs_truth"] <= 0.03
        b["R4"] = n_o2 >= 0.6 * n_t
        out[c] = b
    out["verdict"] = "PASS" if all(out[c][k] for c in COPIES for k in ("R1", "R2", "R3", "R4")) else "FAIL"
    return out


if __name__ == "__main__":
    res = bars(json.load(open(sys.argv[1])))
    json.dump(res, open(sys.argv[2], "w"), indent=1)
    print(json.dumps(res, indent=1))
