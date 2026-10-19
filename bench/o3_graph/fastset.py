#!/usr/bin/env python3
"""A family set's worlds end to end in one command, fast and cheap (2026-10-16 goal: integrate O1, O2 and O3 so that whole
families are fast and cheap to test).

    fastset.py run OUTROOT FAMILIES SET SEED [--worlds INTACT,CROSS,OUTSIDE,REMOVED,BRANCH,SILENT] [--het 0.001]
                   [--local DIR] [--jobs 2] [--threads 4] [--sim-profile NAME] [--o2-distance clip|gapcomp]
                   [--member-genome-check] [--o2-genome off|fallback|combined] [--per-read-detectability]
                   [--trim-adapters ADAPTERS.json]

1. the set's homology-closed local reference (localref.py build; once per family set, cached in --local or OUTROOT/local):
   the only step that loads the whole-genome index
2. REMOVED: the masked genome (realworld.py mask) and the masked local reference (localref.py mask, < 1 s)
3. every world simulated (realworld.py build)
4. every distinct read sequence of the full-genome worlds mapped ONCE to the local reference, its alignments handed to each
   world's read with that sequence (minimap2 maps a read on its sequence alone); REMOVED's reads mapped to the masked one
5. per world, `--jobs` at a time under the light lock (each < 2 GB): the current method (run.py cur), then the integrated
   O1 -> O2 -> O3 (o1seed.py, o2graph.py run W integ 2 [+ integl|integc 2 with --o2-genome fallback|combined], detect.py)
6. summary.py e2e + bars (on that O2 table with --o2-genome; --per-read-detectability passed through); OUTROOT/timing.json holds every step's wall time and peak RSS
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
GENOME_STEP = {"fallback": "integl", "combined": "integc", "seeded": "integs", "resolve": "integr",   # --o2-genome
               "veto": "integv", "guard": "integt"}


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


def _time_fields(text):
    """The last line of /usr/bin/time -o output, split ('%e %M %U %S' since 2026-10-17; '%e %M' before)."""
    lines = [ln.split() for ln in text.splitlines() if ln.split()]
    return lines[-1] if lines else []


def parse_time(text):
    """/usr/bin/time -f '%e %M %U %S' (or '%e %M') output -> (wall s, peak RSS KB); (nan, 0) when the step left none (round
    12: rlock's timeout killed it before /usr/bin/time wrote)."""
    f = _time_fields(text)
    try:
        return float(f[0]), int(f[1])
    except (IndexError, ValueError):
        return float("nan"), 0


def parse_cpu(text):
    """The step's CPU seconds, user + system (its waited-for children's included), from '%e %M %U %S'; nan when absent."""
    f = _time_fields(text)
    try:
        return round(float(f[2]) + float(f[3]), 2)
    except (IndexError, ValueError):
        return float("nan")


class Steps:
    """Runs and times commands under tools/rlock.sh; collects {world, step, wall s, peak RSS KB, rc}."""

    def __init__(self, root, timeout=590):
        self.root, self.rows, self.timeout = root, [], timeout

    def run(self, world, step, cmd, lock="light", timeout=None):
        timeout = timeout or self.timeout
        tf = f"{self.root}/.time.{world}.{step}"
        log = f"{self.root}/{world}.log"
        t0 = time.time()
        rc = subprocess.run(f"RLOCK_TIMEOUT={timeout} bash tools/rlock.sh {lock} /usr/bin/time -o {tf} -f '%e %M %U %S' {cmd} "
                            f">> {log} 2>&1", shell=True, cwd=REPO, executable="/bin/bash").returncode
        text = open(tf).read() if os.path.exists(tf) else ""
        wall, rss = parse_time(text)
        self.rows.append(dict(world=world, step=step, wall=wall, rss_kb=rss, rc=rc, elapsed=round(time.time() - t0, 1),
                              cpu=parse_cpu(text)))
        if rc:
            raise RuntimeError(f"{world} {step} failed (exit {rc}); see {log}")


def _sig(families, set_name):
    import hashlib
    return hashlib.sha1(json.dumps(json.load(open(families))["sets"][set_name], sort_keys=True).encode()).hexdigest()


def write_local_source(local, families, set_name):
    """Record which family set a local reference was built from (round 9 review: a reused reference was never checked)."""
    os.makedirs(local, exist_ok=True)
    json.dump(dict(families=os.path.abspath(families), set=set_name, sha1=_sig(families, set_name)), open(f"{local}/source.json", "w"))


def local_matches(local, families, set_name):
    p = f"{local}/source.json"
    return os.path.exists(p) and json.load(open(p))["sha1"] == _sig(families, set_name)


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


def chain_steps(a, root, w):
    """One world's steps after mapping, [(step, command, lock)]: the current method, then the integrated O1 -> O2 (+ the read's own
    genome alignments, --o2-genome fallback|combined: o2graph.py run W integl|integc 2, which O3 does not read) -> O3."""
    W = f"{root}/{w}"
    steps = [("current", f"python3 bench/o3_graph/run.py cur {W} {CUR_BIN} --threads {a.set_threads if a.heavy_current else a.threads}",
              "heavy" if a.heavy_current else "light"),      # full families: copy_assign peaks ~5 GB, over the light limit
             ("O1", f"python3 bench/o3_graph/o1seed.py {W}" + (" --twins" if a.twins else "")
                    + ((" --ref-copies --local " + (f"{root}/local_masked" if w.startswith("REMOVED") else (a.local or f"{root}/local"))
                        + (" --ref-carried" if a.ref_carried else "")) if a.ref_copies else ""), "light"),
             ("O2", f"env O2_PROCS={a.o2_procs} O2_THREADS={a.threads}"
                    + (f" O2_ALN_CACHE={root}/.o2_aln_cache" if a.aln_cache else "")
                    + f" python3 bench/o3_graph/o2graph.py run {W} integ 2 --distance {a.o2_distance}"
                    + (" --union-exact" if a.union_exact else "") + (" --error-floor" if a.error_floor else "")
                    + (" --twin-abstain" if a.twin_abstain else "") + (" --o2-significance" if a.o2_significance else ""),
                    "light")]
    if a.o2_genome != "off":
        steps.append(("O2_genome", (f"env O2_ALN_CACHE={root}/.o2_aln_cache " if a.aln_cache and (a.struct_rescue or a.short_path_align)
                                    else "") + f"python3 bench/o3_graph/o2graph.py run {W} {GENOME_STEP[a.o2_genome]} 2"
                      + (" --error-floor" if a.error_floor else "") + (" --twin-abstain" if a.twin_abstain else "")
                      + (" --veto-abstain" if a.veto_abstain else "") + (" --o2-significance" if a.o2_significance else "")
                      + (" --locus-pairwise" if a.locus_pairwise else "") + (" --structure-check" if a.structure_check else "")
                      + (" --struct-rescue" if a.struct_rescue else "") + (" --short-path-align" if a.short_path_align else "")
                      + (" --twin-guard" if a.twin_guard else ""), "light"))
    steps.append(("O3", f"env O3_GRAPH_CACHE={root}/.graph_cache O3_PROCS={a.o3_procs} O3_THREADS={a.threads}"
                  + (f" O3_ALN_CACHE={root}/.o3_aln_cache" if a.aln_cache else "")
                  + f" python3 bench/o3_graph/detect.py {W}"
                  + (" --member-genome-check" if a.member_genome_check else "")
                  + (" --member-attribution" if a.member_attribution else "")
                  + (" --member-attribution-delta" if a.member_attribution_delta else "")
                  + (" --error-floor" if a.error_floor else "") + (" --o1-hosts" if a.o1_hosts else "")
                  + (" --mapped-sites" if a.mapped_sites else "") + (" --branch-genome-novel" if a.branch_genome_novel else "")
                  + (" --allele-discount" if a.allele_discount else "") + (" --mosaic-shared" if a.mosaic_shared else "")
                  + (" --haplotype-check" if a.haplotype_check else "")
                  + (" --haplotype-check-genome" if a.haplotype_check_genome else "")
                  + (" --branch-genome-segments" if a.branch_genome_segments else "")
                  + (" --branch-linear" if a.branch_linear else "")
                  + (f" --pool-table o2.{GENOME_STEP[a.o2_genome]}.m2.tsv" if a.o3_pool == "final" and a.o2_genome != "off" else "")
                  + ("" if a.all_arms else " --arms I"),
                  "light"))
    return steps


def e2e_command(a, root, worlds):
    """The set's end-to-end scoring: summary.py e2e on the integrated O2 table the run made (o2.integl.m2.tsv or o2.integc.m2.tsv
    with --o2-genome), per-read detectability with --per-read-detectability."""
    return (f"python3 bench/o3_graph/summary.py e2e {root}/e2e.json " + " ".join(f"{root}/{w}" for w in worlds)
            + (f" --o2 o2.{GENOME_STEP[a.o2_genome]}.m2.tsv" if a.o2_genome != "off" else "")
            + (" --per-read-detectability" if a.per_read_detectability else "")
            + (" --error-floor" if a.error_floor else ""))


def cmd_run(a):
    root = os.path.abspath(a.outroot)
    os.makedirs(root, exist_ok=True)
    worlds = a.worlds.split(",")
    st = Steps(root, a.step_timeout)
    t0 = time.time()
    local = a.local or f"{root}/local"
    if not os.path.exists(f"{local}/local.mmi"):                                   # 1: once per family set
        st.run("set", "localref", f"python3 bench/o3_graph/localref.py build {a.families} {a.set} {local}", lock="heavy")
        write_local_source(local, a.families, a.set)
    elif os.path.exists(f"{local}/source.json") and not local_matches(local, a.families, a.set):
        raise SystemExit(f"{local} was built for another family set; pass a fresh --local")
    het = f"--het {a.het}" if a.het > 0 else ""
    prof = (f"--sim-profile {a.sim_profile}" if a.sim_profile else "") + (" --het-draw" if a.sim_het_draw else "")
    if "REMOVED" in worlds:                                                        # 2
        st.run("set", "mask_genome", f"python3 bench/o3_graph/realworld.py mask {root}/masked {a.families} --set {a.set} "
               f"--seed {a.seed}")
        st.run("set", "mask_local", f"python3 bench/o3_graph/localref.py mask {local} {root}/masked/mask.bed {root}/local_masked")
    for w in worlds:                                                               # 3
        g = f"--genome {root}/masked/masked.fa" if w == "REMOVED" else ""
        x = "--projection liftover" if w == "CROSS" else ""
        st.run(w, "simulate", f"python3 bench/o3_graph/realworld.py build {root}/{w} {a.families} --set {a.set} --scenario {w} "
               f"--seed {a.seed} {het} {prof} {g} {x}")
    if a.trim_adapters:                                                            # 3b: an untrimmed library's adapters off
        for w in worlds:
            os.replace(f"{root}/{w}/reads.E0.fq", f"{root}/{w}/reads.E0.raw.fq")
            st.run(w, "trim", f"python3 bench/o3_graph/trim.py fastq {root}/{w}/reads.E0.raw.fq {root}/{w}/reads.E0.fq {a.trim_adapters}")
    full = [w for w in worlds if w != "REMOVED"]                                   # 4
    t = time.time()
    n_union, n_reads = map_union(root, full, local, a.set_threads) if full else (0, 0)   # alone: every core it may use
    st.rows.append(dict(world="set", step="map_union", wall=round(time.time() - t, 1), rss_kb=0, rc=0,
                        union_reads=n_union, world_reads=n_reads))
    if "REMOVED" in worlds:
        st.run("REMOVED", "map", f"python3 bench/o3_graph/localref.py map {root}/REMOVED {root}/local_masked --threads {a.set_threads}",
               lock="heavy")

    def chain(w):                                                                  # 5
        for step, cmd, lock in chain_steps(a, root, w):
            st.run(w, step, cmd, lock=lock)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in [ex.submit(chain, w) for w in worlds]:
            f.result()
    st.run("set", "e2e", e2e_command(a, root, worlds))
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
    p.add_argument("--sim-het-draw", action="store_true", help="the simulator draws each copy's het/hom site counts instead of "
                   "rounding (realworld.py --het-draw)")
    p.add_argument("--local")
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--set-threads", type=int, default=4, help="threads for steps that run alone: the shared mapping, REMOVED's "
                   "mapping, copy_assign under --heavy-current (output does not depend on it)")
    p.add_argument("--sim-profile", default="")
    p.add_argument("--o2-distance", choices=("clip", "gapcomp", "gapstruct"), default="clip")
    p.add_argument("--member-genome-check", action="store_true")
    p.add_argument("--union-exact", action="store_true", help="O2: paralog structures project despite a dropped terminal micro-exon "
                   "and stay apart from a copy's own structure unless their junctions match exactly (o2graph.py --union-exact)")
    p.add_argument("--error-floor", action="store_true", help="every read-level comparison with delta subtracts the world's "
                   "sequencing-error floor (median de of MAPQ-60 primaries): O2's unexplained gate, O3's genome check, and the "
                   "scorer's per-read detectability (delta is an error-free haplotype divergence)")
    p.add_argument("--step-timeout", type=int, default=590, help="seconds before rlock kills a step (output does not depend on it)")
    p.add_argument("--branch-genome-novel", action="store_true", help="O3's branch arm clusters only reads with >= MIN_NOVEL "
                   "genome-unexplained bases (detect.py --branch-genome-novel)")
    p.add_argument("--haplotype-check", action="store_true", help="arm I drops a PSV call whose haplotype -- its members' shared "
                   "substitutions, over every read carrying them -- stays within delta (detect.py --haplotype-check)")
    p.add_argument("--haplotype-check-genome", action="store_true", help="arm I also drops a PSV call whose members' haplotype, "
                   "in genome coordinates (their own primaries), stays within delta (detect.py --haplotype-check-genome)")
    p.add_argument("--twin-abstain", action="store_true",
                   help="O2: a read neither of two twins beats by the margin is abstained (o2graph.py --twin-abstain)")
    p.add_argument("--veto-abstain", action="store_true", help="O2's genome step: a path call the read's genome alignments "
                   "contradict on its aligned bases, with no confident winner to move it to, is abstained (o2graph.py "
                   "--veto-abstain)")
    p.add_argument("--o2-significance", action="store_true", help="O2: edits priced by their type's rate in the world "
                   "(integ: wdist.json), and the genome step abstains a path call that does not beat every rival by "
                   "ln(rivals / alpha) (o2graph.py --o2-significance)")
    p.add_argument("--locus-pairwise", action="store_true", help="O2's genome step: a locus-fallback call must beat every "
                   "rival by the margin counting, against each, the read bases the rival aligns and the call clips "
                   "(o2graph.py --locus-pairwise)")
    p.add_argument("--structure-check", action="store_true", help="O2's genome step: a path call is abstained when a rival "
                   "whose own locus explains the read better than its paths do is not beaten on the genome by the margin "
                   "(o2graph.py --structure-check)")
    p.add_argument("--struct-rescue", action="store_true", help="O2's genome step: a thin or tied path call decided again on "
                   "structure-aware path distances, a long indel at its length (o2graph.py --struct-rescue)")
    p.add_argument("--short-path-align", action="store_true", help="O2's genome step: reads the path step left without any "
                   "alignment (under ~200 bp: map-hifi's -s 200) aligned to the paths with -s 40 (o2graph.py --short-path-align)")
    p.add_argument("--twin-guard", action="store_true", help="O2's genome step: a call is abstained when a locus off the "
                   "catalogue (neither a seeded copy nor one O1 dropped as silent) explains the read as well "
                   "(o2graph.py --twin-guard)")
    p.add_argument("--ref-carried", action="store_true", help="with --ref-copies: the family rule's carried intron chain (extra "
                   "introns at the locus allowed) and every minimap2 chain kept (-P) (o1seed.py --ref-carried; register 1325)")
    p.add_argument("--ref-copies", action="store_true", help="O1: the catalogued copies' transcripts searched in the world's local "
                   "reference; an off-catalogue locus holding one by the family rule is seeded as a new copy (o1seed.py --ref-copies)")
    p.add_argument("--aln-cache", action="store_true",
                   help="exact alignment caches shared by the set's worlds: O2's path alignments per read sequence "
                   "(O2_ALN_CACHE) and O3's isoform-to-host alignments (O3_ALN_CACHE); outputs identical")
    p.add_argument("--branch-linear", action="store_true",
                   help="O3's branch step on the family's loci as linear sequences (minimap2), no family graph "
                   "(detect.py --branch-linear)")
    p.add_argument("--branch-genome-segments", action="store_true",
                   help="O3's branch step on the pieces each candidate's genome alignment leaves unexplained, no family graph "
                   "(detect.py --branch-genome-segments)")
    p.add_argument("--mosaic-shared", action="store_true", help="O3's mosaic distance only over positions >= 2 members cover "
                   "(detect.py --mosaic-shared)")
    p.add_argument("--o3-pool", choices=("paths", "final"), default="paths", help="arm I's pool: what O2's paths leave "
                   "unexplained (default) or what the --o2-genome step leaves (detect.py --pool-table)")
    p.add_argument("--allele-discount", action="store_true", help="O3's genome check discounts the donor's alleles, the edits "
                   "the reads O2 explains also carry beyond sequencing error (detect.py --allele-discount)")
    p.add_argument("--mapped-sites", action="store_true", help="O3 arm I clusters only on host positions some catalogued "
                   "copy's map covers (detect.py --mapped-sites)")
    p.add_argument("--o1-hosts", action="store_true", help="O3 arm I takes its hosts and their reads from O1 (seeded.json + "
                   "primary alignments), not from copy_assign's run: the integrated pipeline no longer needs the current method")
    p.add_argument("--twins", action="store_true", help="O1 marks a new copy within delta of a curated copy as its twin; O2 then "
                   "gives the twin a read only by the margin (o1seed.py --twins, o2graph.decide twins)")
    p.add_argument("--member-attribution-delta", action="store_true", help="O3 arm I: attribution whose margin is also more "
                   "than delta of each read's aligned length (detect.attribute delta)")
    p.add_argument("--member-attribution", action="store_true", help="O3 arm I: a call keeps only the reads its consensus explains "
                   "better than every catalogued copy by O2's margin (detect.attribute)")
    p.add_argument("--heavy-current", action="store_true", help="run copy_assign under the heavy lock (large families)")
    p.add_argument("--all-arms", action="store_true", help="O3: every legacy arm too (default: the integrated arm I only)")
    p.add_argument("--o3-procs", type=int, default=2, help="O3: families in parallel per world (output does not depend on it)")
    p.add_argument("--o2-procs", type=int, default=2, help="O2: families in parallel per world (output does not depend on it)")
    p.add_argument("--trim-adapters", default="", help="ADAPTERS.json: trim the library's adapters off every read before mapping")
    p.add_argument("--o2-genome", choices=("off",) + tuple(GENOME_STEP), default="off", help="O2 with the read's own genome "
                   "alignments: fallback = only the reads the paths leave unexplained (o2graph.py integl); combined = every read on "
                   "min(path, genome) per copy (integc); seeded = combined over the seeded copies, fallback for the unexplained (integs); resolve = "
                   "fallback plus O2's tied / thin reads decided on min(path, genome) (integr); the set is scored on that table")
    p.add_argument("--per-read-detectability", action="store_true", help="e2e: a REMOVED copy's read is detectable by its own `de`")
    p.set_defaults(func=cmd_run)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
