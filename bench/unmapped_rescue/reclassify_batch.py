#!/usr/bin/env python3
"""Re-classify the finished runs of one animal with a consensus variant (Amendment 37): rebuild the consensus with the abPOA mode, align all runs' consensus sequences to each assembly in ONE
minimap2 call per assembly (one index load), split the PAFs back per run, then classify + rescore each run. Miniforge python. Resumable: exit 75 = run again.

    reclassify_batch.py <animal> <mode>"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import discover as D  # noqa: E402

SEP = "@@"


def main(animal, mode, suffixes=None):
    """suffixes: comma-separated run-name suffixes (default: the real run and the controls of seeds 5, 6, 7); a custom list gets its own pool and PAF names"""
    t0 = time.time()
    sufs = suffixes.split(",") if suffixes else ["", "_control", "_control_s6", "_control_s7"]
    tag = f"{mode}.{suffixes.replace(',', '+').replace('_', '')}" if suffixes else mode
    runs = [r for r in (animal + x for x in sufs) if os.path.exists(f"{D.W}/discover_{r}/clusters.tsv")]
    pool = f"{D.W}/batch_{animal}.{tag}.fa"
    if not os.path.exists(pool):
        with open(pool, "w") as o:
            for r in runs:
                d = f"{D.W}/discover_{r}"
                if not os.path.exists(D.variant_paths(d, mode)["cons"]):
                    D.reconsensus(r, mode)
                for ln in open(D.variant_paths(d, mode)["cons"]):
                    o.write(f">{r}{SEP}{ln[1:]}" if ln[0] == ">" else ln)
    assemblies = {"R": D.ANIMALS[animal]["R"], **D.ANIMALS[animal]["others"]}
    for which, (idx, _fa, _nm) in assemblies.items():
        flag = f"{D.W}/batch_{animal}.{tag}.{which}.paf.done"
        if os.path.exists(flag):
            continue
        if time.time() - t0 > D.TIME_BUDGET:
            print(f"{which}: paused after {time.time() - t0:.0f} s; run again")
            sys.exit(75)
        paf = f"{D.W}/batch_{animal}.{tag}.{which}.paf"
        subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {idx} {pool} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(flag, "w").write("ok")
        print(f"{which}: aligned in {time.time() - t0:.0f} s")
    for r in runs:
        d = f"{D.W}/discover_{r}"
        vp = D.variant_paths(d, mode)
        for which in assemblies:
            out = vp["paf"].format(which)
            with open(out, "w") as o:
                for ln in open(f"{D.W}/batch_{animal}.{tag}.{which}.paf"):
                    if ln.startswith(r + SEP):
                        o.write(ln[len(r) + len(SEP):])
            open(out + ".done", "w").write("ok")
    for r in runs:
        subprocess.run([sys.executable, f"{HERE}/discover.py", "classify", r, animal, mode], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([sys.executable, f"{HERE}/discover.py", "rescore", r, animal, mode], check=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
