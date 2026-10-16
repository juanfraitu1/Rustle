#!/usr/bin/env python3
"""A family set's worlds end to end in one command, fast and cheap (2026-10-16 goal: integrate O1, O2 and O3 so that whole
families are fast and cheap to test).

    fastset.py run OUTROOT FAMILIES SET SEED [--worlds INTACT,CROSS,OUTSIDE,REMOVED,BRANCH,SILENT] [--het 0.001]
                   [--local DIR] [--jobs 2] [--threads 4] [--sim-profile NAME]

1. the set's homology-closed local reference (localref.py build; once per family set, cached in --local or OUTROOT/local):
   the only step that loads the whole-genome index
2. REMOVED: the masked genome (realworld.py mask) and the masked local reference (localref.py mask, < 1 s)
3. every world simulated (realworld.py build)
4. every distinct read sequence of the full-genome worlds mapped ONCE to the local reference, its alignments handed to each
   world's read with that sequence (minimap2 maps a read on its sequence alone); REMOVED's reads mapped to the masked one
5. per world, `--jobs` at a time under the light lock (each < 2 GB): the current method (run.py cur), then the integrated
   O1 -> O2 -> O3 (o1seed.py, o2graph.py run W integ 2, detect.py)
6. summary.py e2e + bars; OUTROOT/timing.json holds every step's wall time and peak RSS
"""
import argparse
import collections
import concurrent.futures
import importlib.util
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
WORLDS = ("INTACT", "CROSS", "OUTSIDE", "REMOVED", "BRANCH", "SILENT")
CUR_BIN = "/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb"


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def union_reads(worlds):
    """{world: {read name: sequence}} -> (union {uname: sequence}, owners {uname: [(world, read name)]})."""
    by_seq, union, owners = {}, {}, collections.defaultdict(list)
    for w in sorted(worlds):
        for n, s in sorted(worlds[w].items()):
            u = by_seq.get(s)
            if u is None:
                u = by_seq[s] = f"u{len(by_seq)}"
                union[u] = s
            owners[u].append((w, n))
    return union, dict(owners)


def split_union(lines, owners):
    """SAM records of union reads -> {world: [the records, renamed to each world read with that sequence]}."""
    out = collections.defaultdict(list)
    for ln in lines:
        u, rest = ln.split("\t", 1)
        for w, n in owners.get(u, ()):
            out[w].append(f"{n}\t{rest}")
    return dict(out)


def read_fq(path):
    out, it = {}, open(path)
    for h in it:
        s = next(it).rstrip("\n")
        next(it)
        next(it)
        out[h[1:].split()[0].rstrip("\n")] = s
    return out


class Steps:
    """Runs and times commands under tools/rlock.sh; collects {world, step, wall s, peak RSS KB, rc}."""

    def __init__(self, root):
        self.root, self.rows = root, []

    def run(self, world, step, cmd, lock="light", timeout=590):
        tf = f"{self.root}/.time.{world}.{step}"
        log = f"{self.root}/{world}.log"
        t0 = time.time()
        rc = subprocess.run(f"RLOCK_TIMEOUT={timeout} bash tools/rlock.sh {lock} /usr/bin/time -o {tf} -f '%e %M' {cmd} "
                            f">> {log} 2>&1", shell=True, cwd=REPO, executable="/bin/bash").returncode
        wall, rss = (open(tf).read().split()[-2:] if os.path.exists(tf) else ("nan", "0"))
        self.rows.append(dict(world=world, step=step, wall=float(wall), rss_kb=int(rss), rc=rc, elapsed=round(time.time() - t0, 1)))
        if rc:
            raise RuntimeError(f"{world} {step} failed (exit {rc}); see {log}")


def map_union(root, worlds, local, threads):
    """Step 4: the full-genome worlds' distinct reads mapped once to the local reference, split into each world's E0.bam."""
    L = _load("localref")
    union, owners = union_reads({w: read_fq(f"{root}/{w}/reads.E0.fq") for w in worlds})
    os.makedirs(f"{root}/union", exist_ok=True)
    with open(f"{root}/union/reads.fq", "w") as fh:
        fh.writelines(f"@{u}\n{s}\n+\n{'~' * len(s)}\n" for u, s in union.items())
    fai = f"{root}/{worlds[0]}/genome.ref.fa.fai"
    header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in (ln.split("\t")[:2] for ln in open(fai)))
    sorts = {w: subprocess.Popen(f"samtools sort -o {root}/{w}/E0.bam -", shell=True, stdin=subprocess.PIPE, text=True,
                                 executable="/bin/bash") for w in worlds}
    for p in sorts.values():
        p.stdin.write(header)
    with open(f"{root}/union/mm2.log", "w") as log:
        mm = subprocess.Popen(f"{L.MM2} -t {threads} {local}/local.mmi {root}/union/reads.fq", shell=True,
                              stdout=subprocess.PIPE, stderr=log, text=True, executable="/bin/bash")
        for line in mm.stdout:
            if not line.startswith("@"):
                for w, recs in split_union([L.lift_line(line)], owners).items():
                    sorts[w].stdin.writelines(recs)
        if mm.wait():
            raise RuntimeError("union mapping failed")
    for w, p in sorts.items():
        p.stdin.close()
        if p.wait():
            raise RuntimeError(f"sorting {w}'s alignments failed")
        subprocess.run(["samtools", "index", f"{root}/{w}/E0.bam"], check=True)
        open(f"{root}/{w}/E0.mm2.log", "w").write(f"mapped as part of {root}/union (local reference {local})\n")
    return len(union), sum(len(v) for v in owners.values())


def cmd_run(a):
    root = os.path.abspath(a.outroot)
    os.makedirs(root, exist_ok=True)
    worlds = a.worlds.split(",")
    st = Steps(root)
    t0 = time.time()
    local = a.local or f"{root}/local"
    if not os.path.exists(f"{local}/local.mmi"):                                   # 1: once per family set
        st.run("set", "localref", f"python3 bench/o3_graph/localref.py build {a.families} {a.set} {local}", lock="heavy")
    het = f"--het {a.het}" if a.het > 0 else ""
    prof = f"--sim-profile {a.sim_profile}" if a.sim_profile else ""
    if "REMOVED" in worlds:                                                        # 2
        st.run("set", "mask_genome", f"python3 bench/o3_graph/realworld.py mask {root}/masked {a.families} --set {a.set} "
               f"--seed {a.seed}")
        st.run("set", "mask_local", f"python3 bench/o3_graph/localref.py mask {local} {root}/masked/mask.bed {root}/local_masked")
    for w in worlds:                                                               # 3
        g = f"--genome {root}/masked/masked.fa" if w == "REMOVED" else ""
        x = "--projection liftover" if w == "CROSS" else ""
        st.run(w, "simulate", f"python3 bench/o3_graph/realworld.py build {root}/{w} {a.families} --set {a.set} --scenario {w} "
               f"--seed {a.seed} {het} {prof} {g} {x}")
    full = [w for w in worlds if w != "REMOVED"]                                   # 4
    t = time.time()
    n_union, n_reads = map_union(root, full, local, a.threads) if full else (0, 0)
    st.rows.append(dict(world="set", step="map_union", wall=round(time.time() - t, 1), rss_kb=0, rc=0,
                        union_reads=n_union, world_reads=n_reads))
    if "REMOVED" in worlds:
        st.run("REMOVED", "map", f"python3 bench/o3_graph/localref.py map {root}/REMOVED {root}/local_masked --threads {a.threads}")

    def chain(w):                                                                  # 5
        W = f"{root}/{w}"
        st.run(w, "current", f"python3 bench/o3_graph/run.py cur {W} {CUR_BIN} --threads {a.threads}")
        st.run(w, "O1", f"python3 bench/o3_graph/o1seed.py {W}")
        st.run(w, "O2", f"python3 bench/o3_graph/o2graph.py run {W} integ 2")
        st.run(w, "O3", f"python3 bench/o3_graph/detect.py {W}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in [ex.submit(chain, w) for w in worlds]:
            f.result()
    st.run("set", "e2e", f"python3 bench/o3_graph/summary.py e2e {root}/e2e.json " + " ".join(f"{root}/{w}" for w in worlds))
    subprocess.run(f"python3 bench/o3_graph/summary.py bars {root}/e2e.json > {root}/bars.txt", shell=True, cwd=REPO, check=True)
    total = round(time.time() - t0, 1)
    json.dump(dict(total_wall=total, steps=st.rows), open(f"{root}/timing.json", "w"), indent=1)
    by = collections.defaultdict(float)
    for r in st.rows:
        by[r["step"]] += r["wall"]
    print(f"[fastset] {root}: {len(worlds)} worlds in {total} s; peak RSS {max(r['rss_kb'] for r in st.rows) / 1e6:.2f} GB; "
          + ", ".join(f"{k} {v:.0f} s" for k, v in by.items()) + f" | union {n_union} distinct of {n_reads} reads", flush=True)
    print(open(f"{root}/bars.txt").read(), end="")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("outroot")
    p.add_argument("families")
    p.add_argument("set")
    p.add_argument("seed", type=int)
    p.add_argument("--worlds", default=",".join(WORLDS))
    p.add_argument("--het", type=float, default=0.001)
    p.add_argument("--local")
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--sim-profile", default="")
    p.set_defaults(func=cmd_run)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
