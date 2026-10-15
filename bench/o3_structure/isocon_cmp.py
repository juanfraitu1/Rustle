#!/usr/bin/env python3
"""IsoCon 0.3.3 as the reconstruction comparator (docs/PREREG_o3_structure_projection_2026-10-14.md, "Reported beside").

    isocon_cmp.py WORLD VARIANT ROWS.tsv     # per family with an absent copy E: IsoCon on all its reads (A + B + E, FASTA), its best final
                                             # candidate against E's true transcript; one row appended per family

Run from the default pixi environment; IsoCon itself runs as `pixi run -e isocon IsoCon pipeline ... --nr_cores 4`.
"""
import csv
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proto as P  # noqa: E402

PIXI = os.path.expanduser("~/.pixi/bin/pixi")
MANIFEST = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "pixi.toml")


def main():
    world, v, rows_path = sys.argv[1:4]
    kinds = json.load(open(f"{world}/kinds.json"))
    scen = json.load(open(f"{world}/world.json"))
    T = P.read_fa(f"{world}/transcripts.fa")
    work = f"{world}/{v}/isocon"
    os.makedirs(work, exist_ok=True)
    reads = {}
    lines = open(f"{world}/reads.{v}.fq").read().splitlines()
    for i in range(0, len(lines), 4):
        name = lines[i][1:]
        reads.setdefault(name.split("|")[0][:-1], []).append((name, lines[i + 1]))   # family = copy id without its last letter
    out = []
    for fam in sorted(scen):
        copy = fam + "E"
        if copy not in kinds:
            continue
        d = f"{work}/{fam}"
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        with open(f"{d}/reads.fa", "w") as fh:
            for n, s in reads.get(fam, []):
                fh.write(f">{n}\n{s}\n")
        subprocess.run(f"{PIXI} run --manifest-path {MANIFEST} -e isocon IsoCon pipeline -fl_reads {d}/reads.fa -outfolder {d}/out "
                       f"--nr_cores 4 > {d}/log 2>&1", shell=True, check=False, executable="/bin/bash")
        fin = f"{d}/out/final_candidates.fa"
        row = dict(world=world.rstrip("/").split("/")[-1], variant=v, family=fam, scenario=scen[fam], n_final=0, best_idcov=0.0,
                   best_tcov=0.0, exact=0)
        if os.path.exists(fin) and os.path.getsize(fin) > 0:
            row["n_final"] = sum(1 for ln in open(fin) if ln.startswith(">"))
            with open(f"{d}/E.fa", "w") as fh:
                fh.write(f">{copy}\n{T[copy]}\n")
            for h in P.paf(P.mm2("-c --eqx -x map-hifi -N 5", f"{d}/E.fa", fin)):
                if h["tp"] != "P" or h["st"] != "+":
                    continue
                ic = (h["m"] / h["bl"]) * (h["qe"] - h["qs"]) / h["ql"]
                if ic > row["best_idcov"]:
                    row.update(best_idcov=round(ic, 5), best_tcov=round((h["te"] - h["ts"]) / len(T[copy]), 5))
            row["exact"] = int(row["best_idcov"] >= 0.999)
        out.append(row)
        shutil.rmtree(f"{d}/out/alignments", ignore_errors=True)
    new = not os.path.exists(rows_path)
    with open(rows_path, "a") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        if new:
            w.writeheader()
        w.writerows(out)
    print(f"{world.rstrip('/').split('/')[-1]}\t{v}\tIsoCon families {len(out)}\texact {sum(r['exact'] for r in out)}")


if __name__ == "__main__":
    main()
