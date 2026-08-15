"""Third arm for the min_clusters false-positive survey: gate 1 fully OFF (min_clusters=1).

The 2-vs-3 comparison answers "does the knob cost anything at the value the ablation used". It
does not say how much margin is left. `min_clusters=1` is `n_clusters >= 1`, which every candidate
satisfies, so this arm is gate 1 DISABLED. Running it on the same alignments measures whether 2
sits on a cliff (many admissions at 1) or on a plateau (none), i.e. how much of gate 1's
protection is real versus already supplied by gates 2-5.

This is a MEASUREMENT of the gate, not a proposal: 1 is not a setting anything should ship or run
an experiment at. Reuses the reads.sorted.bam / ref_intact.fa already built by
adv_mincl_fp_survey.py -- only the copy_assign invocation is new.
"""
import argparse
import glob
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import v4b_3copy  # noqa: E402
from adv_mincl_fp_survey import PINNED, arm  # noqa: E402

if pathlib.Path(PINNED).exists():
    v4b_3copy.COPY_ASSIGN = PINNED


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/adv_mincl_fp")
    ap.add_argument("--limit", type=int, default=6)
    ap.add_argument("--only-candidate-loci", action="store_true", default=True)
    a = ap.parse_args()

    done = 0
    for p in sorted(glob.glob(f"{a.scratch}/*/result.json")):
        if done >= a.limit:
            break
        work = pathlib.Path(p).parent
        res = json.loads(pathlib.Path(p).read_text())
        if "default3" not in res or (work / "mincl1.json").exists():
            continue
        d3 = res["default3"]
        if a.only_candidate_loci and not (d3["n_dna_needs"] or sum(d3["collapsed_copies"])):
            continue  # gate never reached here; a third arm would be vacuous
        srt = work / "reads.sorted.bam"
        ref = work / "ref_intact.fa"
        if not srt.exists():
            continue
        chrom = res["locus"].split(":")[0]
        t0 = time.time()
        out = arm("mincl1", ref, f"{chrom}_local", res["local_ref_len"], srt, work, 1)
        out["seconds"] = round(time.time() - t0, 1)
        (work / "mincl1.json").write_text(json.dumps(out, indent=2))
        done += 1
        print(f"{res['family_id']}\tadmit3={d3['n_admitted_absent_copies']}"
              f"\tadmit2={res['mincl2']['n_admitted_absent_copies']}"
              f"\tadmit1={out['n_admitted_absent_copies']}\t{out['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
