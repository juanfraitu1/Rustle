#!/usr/bin/env python3
"""The gates of docs/PREREG_seed_pool_real_reads_2026-10-07.md section 5, computed from the products of run.sh.

    gates.py --work W/SUB --sub human_chr16|human_chr17 --def /mnt/linuxdisk/tmp/rescore_2026-10-06/DEF --bin RELEASE_DIR

G0 (human_chr16): the G98 arm's assembled GTF, families GTF, clusters and loci GFF3 equal the default re-score products byte for byte.
G1 (human_chr16): own node 24, E-found within own nodes 8, locus-level 21, U2 F .645 on the G98 arm.
G1b: composition.py's own-node flags equal nodes.py's on the P and G98 arms (nodes.py stops on a shared span; then the gate is reported as not applicable).
G2: P_R1 clusters equal P's.   G3: copy_support's summary is the same under PYTHONHASHSEED 0 and 1.   G4: the P_SH (sharded) clusters equal P's.
Prints one line per gate and writes gates.json; exits 1 if a gate that applies fails.
"""
import argparse
import filecmp
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(p):
    with open(p) as fh:
        return json.load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--sub", required=True)
    ap.add_argument("--def", dest="defdir", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--copies", default="/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann/copies.hsa.tsv")
    ap.add_argument("--truth", default="/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann/truth.hsa.gtf")
    ap.add_argument("--family", default=None)
    ap.add_argument("--exons-json", default=None)
    a = ap.parse_args()
    S, out = a.work, {}

    def same(x, y):
        return os.path.exists(x) and os.path.exists(y) and filecmp.cmp(x, y, shallow=False)

    def report(gate, ok, note):
        out[gate] = dict(ok=ok, note=note)
        print(f"{gate} {'n/a ' if ok is None else ('PASS' if ok else 'FAIL')}  {note}")

    if a.sub == "human_chr16":
        g = f"{S}/G98/run"
        pairs = [(f"{g}.gtf", "gtf"), (f"{g}.families.gtf", "families.gtf"), (f"{g}.fam.clusters.tsv", "fam.clusters.tsv"), (f"{g}.fam.loci.gff3", "fam.loci.gff3")]
        res = {n: same(x, f"{a.defdir}/human_A119b.chr16.{n}") for x, n in pairs}
        report("G0", all(res.values()), "G98 vs the default re-score products: " + ", ".join(f"{n} {'identical' if v else 'DIFFERS'}" for n, v in res.items()))
        try:
            sup, fs = load(f"{S}/support.json")["arms"]["G98"], load(f"{S}/fs_G98.json")["truths"]["u2"]
            got = (sup["old_overlap_in_npip_nodes"], sup["tc_found_in_npip_nodes"], sup["locus_tc_found"], round(fs["f"], 3))
            report("G1", got == (24, 8, 21, 0.645), f"own node / E-found in own nodes / locus-level / U2 F = {got}, registered (24, 8, 21, 0.645)")
        except (OSError, KeyError) as e:
            report("G1", None, f"not scored yet ({e})")
    else:
        report("G0", None, "the re-score products exist for human chr16 only")
        report("G1", None, "the re-score numbers exist for human chr16 only")

    try:
        meta = load(f"{S}/comp.json")
        copies_tsv, truth = a.copies, a.truth
        exons = ["--exons-json", a.exons_json] if a.exons_json else []
        arms = [x for x in ("P", "G98") if x in meta["arms"]]
        cmd = [sys.executable, f"{HERE}/../default_rescore/nodes.py", "--copies", copies_tsv, "--truth", truth, "--family", meta["family"], *exons,
               *sum([["--arm", f"{x}={S}/{x}/run.fam.loci.gff3,{S}/{x}/run.fam.clusters.tsv"] for x in arms], []), "--out", f"{S}/nodes_check.json"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            report("G1b", None, "nodes.py stops (" + (r.stderr.strip().splitlines() or ["?"])[-1][:120] + ")")
        else:
            old = {x["cid"]: x["node"] for x in load(f"{S}/nodes_check.json")["rows"]}
            new = {x["cid"]: x["node"] for x in load(f"{S}/nodes.json")["rows"]}
            bad = [c for c in old for k in arms if old[c].get(k) != new[c].get(k)]
            report("G1b", not bad, f"composition.py vs nodes.py own-node flags on {'+'.join(arms)}: {len(bad)} differ")
    except (OSError, KeyError) as e:
        report("G1b", None, f"not scored yet ({e})")

    report("G2", same(f"{S}/P_R1/run.fam.clusters.tsv", f"{S}/P/run.fam.clusters.tsv") if os.path.exists(f"{S}/P_R1/run.fam.clusters.tsv") else None,
           "P_R1 clusters equal P's (Rule 1 over its own primaries is the identity)")
    try:
        report("G3", load(f"{S}/support_s0.json") == load(f"{S}/support_s1.json"), "copy_support summary under PYTHONHASHSEED 0 and 1")
    except OSError as e:
        report("G3", None, f"not scored yet ({e})")
    if os.path.exists(f"{S}/P_SH/run.fam.clusters.tsv"):
        report("G4", same(f"{S}/P_SH/run.fam.clusters.tsv", f"{S}/P/run.fam.clusters.tsv"), "sharded all-vs-all clusters equal the single-process run's")
    else:
        report("G4", None, "run.sh g4 not run")
    json.dump(out, open(f"{S}/gates.json", "w"), indent=1)
    if any(v["ok"] is False for v in out.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
