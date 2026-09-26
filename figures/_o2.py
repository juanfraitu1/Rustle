"""_o2 — shared data layer for Figures 4 and 5 (copy assignment, simulated reads with known source copies).

The experiment (docs/PREREG_o2_read_truth_2026-09-23.md; genome-wide on every sample:
docs/PREREG_genome_wide_copy_assignment_2026-09-25.md):

    copy table copies.tsv/.fa --bench/sim.py copies--> reads named `family|copy|i`, mapped genome-wide
        (shipped minimap2 settings) = SIM.bam, SIM.copies_used.tsv
    copy table minus the copies that received no reads  = CAT.copies.tsv/.fa (`--families` refuses a copy with
        no reads in its region)
    copy_assign --families CAT (default output: one row per read and family) = O2.assignments.tsv
    copy_assign --families CAT --union-certificate (union test)          = U2.assignments.tsv + U2.union_certificate.tsv
    bench/score.py reads --catalog CAT SIM O2|U2 --per-read      = the per-read scorer: its table (logged) and one
                                                                   row per MAPQ-0 read with its verdict under each
                                                                   reading (OWN / PRIMARY / ANY)

COPY TABLE (cfg `o2_copy_table`; user decision 2026-09-25: ONE default de novo family definition, and copy
assignment consumes the SAME families):
  families (default)  the families stage's copy table (run cache stage `families`: <sample>.fam.copies.tsv/.fa,
                      `mcl_families --from-gtf --emit-units`): one copy per member locus of every de novo family =
                      the locus representative transcript and its spliced exon sum. Its `max_family_identity` is the
                      best families-stage edge (genomic-span `-x asm20` identity) from the copy to another member.
  catalog             the LEGACY copy catalog (run cache stage `catalog`, gw_family_catalog), only when asked; its
                      `max_family_identity` is an exon-sum alignment identity. Work dirs get a `.catalog` suffix.

SCOPES (cfg `o2_scope`):
  genome (default)  every sample of figures/samples.tsv (cfg `o2_samples` narrows it), each with its own genome-wide
                    copy table (above). The simulation maps in parts sized so each call stays under ~8 min
                    (READS_PER_PART), and copy_assign runs in SHARDS of read-connected families (plan_shards), each
                    shard one call. Work dir ${work}/o2sim/<sample>/ (legacy catalog: <sample>.catalog/).
  dev               the development tables: human A119b chr16 catalog and gorilla OR6737 chr20 (NC_073244.2) catalog
                    (cfg `{human,gorilla}_o2_catalog_*`: LEGACY gw_family_catalog tables, named explicitly by the
                    inputs file), unsharded, work dirs ${work}/o2sim/{human,gorilla}/ (the runs behind the
                    2026-09-25 tables).
A call does at most `o2_heavy_steps_per_call` heavy steps (default 1: one simulation part, or one copy_assign shard)
and then stops with "re-run to continue" (Pending), so every command stays inside the machine's 10-minute rule. Wrap
each call in `flock /mnt/linuxdisk/tmp/rustle_heavy.lock` (the shared one-heavy-process lock).

`per_read()` joins the simulated BAM (truth, the aligner's primary placement) with the scorer's per-read verdicts
of both runs into one record per simulated read; nothing here re-implements the scorer's judge(). The tabulators
turn those records into the tidy tables of Figures 4 and 5, and `crosscheck()` requires the tallies to equal the
ALL rows the scorer printed. `per_read()` also measures, for every MAPQ-0 read, whether it has an NM-identical twin
(`twin_state`): a placement at another locus whose AS equals the best and whose NM equals the source-copy
placement's.

Identity of a source copy = the copy table's own `max_family_identity` (its identity to the most similar directly
aligned copy of its family; which alignment depends on the copy table, above). Samples and species are never pooled: every record and every table row carries both.
"""
from __future__ import annotations

import bisect
import collections
import csv
import hashlib
import json
import math
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import figlib

BENCH = figlib.REPO / "bench"
SCORE = BENCH / "score.py"
SIM = BENCH / "sim.py"
sys.path.insert(0, str(BENCH))
import score as _score  # noqa: E402  (bench/score.py: its locus rule `same_locus`, used for the aligner baseline)
import sim as _sim  # noqa: E402  (bench/sim.py: MM2, copies_selection, _fastq_key)

SCOPES = ("genome", "dev")
SPECIES_ORDER = ["human", "chimpanzee", "gorilla", "orangutan"]
# dev scope: the two species-keyed runs and the sample each belongs to
DEV_SPECIES = ["human", "gorilla"]
DEV_SAMPLE = {"human": "human_A119b", "gorilla": "gorilla_OR6737"}
DEV_CATALOG = {"human": "chr16", "gorilla": "chr20 (NC_073244.2)"}
# recorded runs (hash-seeded, 2026-09-23/24) in cfg['o2sim_dir']: human chr16 catalog, gorilla NC_073244.2 catalog
RECORDED_PREFIX = {"human": "h16", "gorilla": "g44"}
# display order of every per-sample panel: grouped by species, never pooled (figures/samples.tsv ids)
SAMPLE_ORDER = ["human_A119b", "human_testis", "chimp_PTR", "gorilla_OR6737", "gorilla_KB3781", "orangutan_PPY"]
SAMPLE_SPECIES = {"human_A119b": "human", "human_testis": "human", "chimp_PTR": "chimpanzee",
                  "gorilla_OR6737": "gorilla", "gorilla_KB3781": "gorilla", "orangutan_PPY": "orangutan"}
SPECIES_TITLE = {"human": "Human", "chimpanzee": "Chimpanzee", "gorilla": "Gorilla", "orangutan": "Orangutan"}
# small-multiple grids (rows x columns), species kept together: 2 x 3 for Fig. 5a-b, 3 x 2 for the Fig. 4 supplement
GRID_2x3 = [["human_A119b", "human_testis", "chimp_PTR"], ["gorilla_OR6737", "gorilla_KB3781", "orangutan_PPY"]]
GRID_3x2 = [["human_A119b", "human_testis"], ["gorilla_OR6737", "gorilla_KB3781"], ["chimp_PTR", "orangutan_PPY"]]
SAMPLE_LABEL = {"human_A119b": "Human A119b", "human_testis": "Human testis", "gorilla_OR6737": "Gorilla OR6737",
                "gorilla_KB3781": "Gorilla KB3781", "chimp_PTR": "Chimpanzee", "orangutan_PPY": "Orangutan"}
SAMPLE_DETAIL = {"human_A119b": "CHM13 v2.0", "human_testis": "testis, CHM13 v2.0",
                 "gorilla_OR6737": "testis, mGorGor1", "gorilla_KB3781": "fibroblast, mGorGor1",
                 "chimp_PTR": "mPanTro3", "orangutan_PPY": "mPonPyg2"}
# simulated reads per mapping part (PREREG_genome_wide_copy_assignment §2): ~8 min per part at the measured
# 0.021 s/read (human, CHM13 index) and 0.098 s/read (gorilla index); chimpanzee and orangutan take gorilla's rate
READS_PER_PART = {"human": 17000, "gorilla": 4400, "chimpanzee": 4400, "orangutan": 4400}
# copy_assign's read window around every supplied copy (src/bin/copy_assign.rs COPY_READ_PAD)
COPY_READ_PAD = 50_000
# shard cost model: (seconds per record inside the shard's read windows, fixed seconds per shard, seconds per record
# of a contig swept whole). Measured 2026-09-25 (check V3, human chr16 simulation, frozen binary sha1 22b3deb4):
# 287 families / 91,374 window records took 123 s (default) and 144 s (--union-certificate); 3 families / 312
# window records plus 69,817 records on the 24 contigs swept whole took 7-8 s. Real reads: the hard-locus run
# (--gtf, 415,854 records in the copy regions) 88 s.
SHARD_COST = {"sim": (1.6e-3, 10.0, 1e-4), "real": (2.2e-4, 30.0, 1e-4)}
SHARD_BUDGET_S = 480.0

try:  # the sample registry (figures/samples.py); absent only in a stripped checkout
    import samples  # noqa: E402
except ImportError:  # pragma: no cover
    samples = None


class Pending(Exception):
    """A bounded call did its share of heavy work; the same command must be run again to continue."""


class Budget:
    """At most `steps` heavy steps (one mapping part, one copy_assign shard) in this call."""

    def __init__(self, steps: int):
        self.steps = steps
        self.done: list[str] = []

    def take(self, what: str):
        if self.steps <= 0:
            raise Pending(f"heavy-step budget of this call used ({'; '.join(self.done)}); next: {what}")
        self.steps -= 1
        self.done.append(what)


def budget_from(cfg: dict) -> Budget:
    return Budget(int(cfg.get("o2_heavy_steps_per_call", "1")))


# ================================================================ samples and scopes
def o2_scope(cfg: dict) -> str:
    s = cfg.get("o2_scope", "genome")
    if s not in SCOPES:
        raise SystemExit(f"o2_scope={s!r}: expected one of {SCOPES}")
    return s


def sample_ids(cfg: dict) -> list[str]:
    """Samples of the build, registry order (cfg `o2_samples` = comma-separated ids or aliases narrows it)."""
    if o2_scope(cfg) == "dev":
        return [DEV_SAMPLE[sp] for sp in DEV_SPECIES]
    reg = list(samples.registry(cfg))
    want = [x.strip() for x in cfg.get("o2_samples", "").split(",") if x.strip()]
    return [samples.resolve(cfg, x) for x in want] if want else reg


def species_of(cfg: dict, sid: str) -> str:
    return samples.get(cfg, sid)["species"]


COPY_TABLES = {"families": "families", "catalog": "catalog"}   # cfg o2_copy_table -> run cache stage


def copy_table(cfg: dict) -> str:
    """Which copy table copy assignment consumes (cfg `o2_copy_table`): 'families' (default) or the legacy 'catalog'."""
    t = cfg.get("o2_copy_table", "families")
    if t not in COPY_TABLES:
        raise SystemExit(f"o2_copy_table={t!r}: expected one of {sorted(COPY_TABLES)} (families = the default de novo "
                         "family definition; catalog = the legacy gw_family_catalog)")
    return t


def genome_work_dir(cfg: dict, sid: str) -> Path:
    """${work}/o2sim/<sample>/ for the families copy table, <sample>.catalog/ for the legacy catalog (never mixed)."""
    return figlib.work_dir(cfg, "o2sim") / (sid if copy_table(cfg) == "families" else f"{sid}.catalog")


def catalog_paths(cfg: dict, sid: str) -> tuple[Path, Path]:
    """The sample's genome-wide copy table (copies.tsv, copies.fa): the families stage's by default, the legacy
    gw_family_catalog when cfg o2_copy_table=catalog; stops with the command that builds it."""
    stage = COPY_TABLES[copy_table(cfg)]
    tsv = samples.product(cfg, sid, stage, "copies")
    fa = samples.product(cfg, sid, stage, "copies_fa")
    if not (tsv.exists() and fa.exists()):
        raise SystemExit(f"[_o2] {sid}: no genome-wide copy table ({tsv}); run `python3 figures/make.py runs "
                         f"--sample {sid} --stage {stage}` first (under flock)")
    state, reason = samples.status(cfg, sid, stage)
    if state not in ("fresh", "adopt") and not cfg.get("o2_accept_stale_catalog"):
        raise SystemExit(f"[_o2] {sid}: the genome-wide copy table (stage {stage}) is {state} ({reason}); re-run "
                         f"`make.py runs --sample {sid} --stage {stage}`, or set o2_accept_stale_catalog=1 to use it "
                         "as it is")
    return tsv, fa


def splice_index(cfg: dict, sid: str) -> Path:
    row = samples.get(cfg, sid)
    mmi = Path(row["splice_mmi"]) if row["splice_mmi"] else samples.products(cfg, sid, "index")["mmi"]
    if not mmi.exists():
        raise SystemExit(f"[_o2] {sid}: splice index {mmi} absent; run `make.py runs --sample {sid} --stage index`")
    return mmi


def planned_parts(cfg: dict, species: str, cat_tsv, cat_fa) -> tuple[int, int, int]:
    """(reads the simulator will write, reads per part, parts) — from the catalog, before any read is simulated."""
    sel = _sim.copies_selection(str(cat_tsv), str(cat_fa))[3]
    n = sum(k for _, k, _ in sel)
    per = int(cfg.get(f"o2_reads_per_part_{species}", READS_PER_PART.get(species, 4400)))
    return n, per, max(1, math.ceil(n / per))


# ================================================================ runs: recorded or rebuilt
def recorded_runs(cfg: dict, species: str) -> dict:
    """The pre-2026-09-24 runs on disk (hash()-seeded simulation; the per-copy read streams are not reproducible)."""
    d = Path(cfg["o2sim_dir"])
    p = RECORDED_PREFIX[species]
    return {"species": species, "sample": DEV_SAMPLE[species], "catalog_scope": DEV_CATALOG[species],
            "sim": d / p, "catalog": d / f"{p}_cat.copies.tsv", "o2": d / f"{p}_o2", "u2": d / f"{p}_u2",
            "source": f"recorded run {d}/{p}* (hash-seeded, 2026-09-23/24)", "notes": []}


def _stamp_ok(path: Path, text: str) -> bool:
    try:
        return path.read_text() == text
    except OSError:
        return False


def _md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _write_if_changed(path: Path, text: str) -> bool:
    """Replace `path` only when its content differs (an unchanged product keeps its mtime, so the make-style checks
    of the steps downstream stay valid). Returns True when the file was (re)written."""
    try:
        if path.read_text() == text:
            return False
    except OSError:
        pass
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return True


def _clean_env_prefix() -> list[str]:
    """`env -u VAR ...` for every RUSTLE_* variable of the calling shell: copy_assign reads ~20 of them
    (RUSTLE_PSV_GENOMIC, RUSTLE_INTRON_PSV, RUSTLE_POSTERIOR_PRIOR, ...), the cache cannot see the environment,
    and the captions say "default settings"."""
    return ["env", *[x for k in sorted(os.environ) if k.startswith("RUSTLE_") for x in ("-u", k)]]


CLEAN_ENV = "RUSTLE_* unset (env -u)"


def _fa_records(path):
    """(family, copy_idx) -> list of FASTA lines (header included), in file order."""
    recs: dict = {}
    order = []
    key = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].split("|")
                key = (h[0], h[1].strip())
                recs[key] = [line]
                order.append(key)
            elif key is not None:
                recs[key].append(line)
    return recs, order


def _primary_count(bam, chrom: str, s: int, e: int) -> int:
    """`samtools view -c -F 2308 BAM chrom:s+1-e` (primary, mapped, non-supplementary records overlapping [s, e))."""
    return int(subprocess.run(["samtools", "view", "-c", "-F", "2308", str(bam), f"{chrom}:{s + 1}-{e}"],
                              capture_output=True, text=True, check=True).stdout.strip())


def derive_catalog(catalog_tsv, catalog_fa, sim_prefix, out_prefix, *, force=False) -> dict:
    """The catalog handed to `copy_assign --families`: the simulated catalog minus every copy that was NOT simulated
    (spliced sequence < 300 bp, or a single-copy family) AND carries no primary alignment (-F 2308) over its span in
    the simulated BAM — `--families` aborts on a copy with no reads in its region. Reproduces the recorded `h16_cat`
    (18 dropped, 2 unsimulated-but-covered copies kept) exactly."""
    out_prefix = Path(out_prefix)
    tsv, fa, dropped = (Path(f"{out_prefix}.copies.tsv"), Path(f"{out_prefix}.copies.fa"),
                        Path(f"{out_prefix}.dropped.tsv"))
    bam = f"{sim_prefix}.bam"
    used_p = f"{sim_prefix}.copies_used.tsv"
    if not force and all(figlib.fresh(p, catalog_tsv, catalog_fa, bam, used_p, __file__) for p in (tsv, fa, dropped)):
        return {"tsv": tsv, "fa": fa, "dropped": dropped}
    # (re)derived when this module changes too; each product is replaced only when its content differs, so an
    # unchanged catalog does not invalidate the copy_assign runs that read it
    used = {(r["family_id"], r["copy_idx"]) for r in csv.DictReader(open(used_p), delimiter="\t")}
    recs, order = _fa_records(catalog_fa)
    with open(catalog_tsv) as fh:
        header = fh.readline()
        rows = [l for l in fh if l.strip()]
    keep, drop = [], []
    for line in rows:
        f = line.rstrip("\n").split("\t")
        k = (f[0], f[1])
        if k in used:
            keep.append(line)
            continue
        chrom, s, e = f[3], int(f[4]), int(f[5])
        if _primary_count(bam, chrom, s, e) > 0:
            keep.append(line)
        else:
            seq_len = sum(len(x.strip()) for x in recs.get(k, [""])[1:])
            drop.append([f[0], f[1], chrom, s, e, seq_len, f[8], "not simulated, no primary read over its span"])
    kept = {tuple(l.split("\t")[:2]) for l in keep}
    _write_if_changed(tsv, header + "".join(keep))
    _write_if_changed(fa, "".join("".join(recs[k]) for k in order if k in kept))
    _write_if_changed(dropped, "family_id\tcopy_idx\tchrom\tstart\tend\tlen\tn_reads_real\treason\n"
                      + "".join("\t".join(map(str, r)) + "\n" for r in drop))
    return {"tsv": tsv, "fa": fa, "dropped": dropped}


def bam_contigs(bam) -> list[tuple[str, int]]:
    hdr = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True).stdout
    out = []
    for line in hdr.splitlines():
        if line.startswith("@SQ"):
            tags = dict(t.split(":", 1) for t in line.split("\t")[1:] if ":" in t)
            out.append((tags["SN"], int(tags["LN"])))
    return out


def _regions_text(contigs) -> str:
    return "".join(f"{c}:1-{n}\n" for c, n in contigs)


def ensure_runs(cfg: dict, sample: str, *, force=False, budget: Budget | None = None) -> dict:
    """Simulation plus both copy_assign runs of one sample (cached; HEAVY unless cached). Dev scope: `sample` is the
    species key ('human' / 'gorilla') or its sample id; genome scope: any sample id or alias."""
    if o2_scope(cfg) == "dev":
        sp = sample if sample in DEV_SPECIES else next(k for k, v in DEV_SAMPLE.items()
                                                        if v == samples.resolve(cfg, sample))
        return _ensure_runs_dev(cfg, sp, force=force)
    return _ensure_runs_genome(cfg, samples.resolve(cfg, sample), force=force, budget=budget or budget_from(cfg))


def _ensure_runs_dev(cfg: dict, species: str, *, force=False) -> dict:
    """Rebuild the development simulation and both assignment runs under `${work}/o2sim/<species>/` (cached).

    HEAVY: `sim.py copies` maps every read genome-wide with the shipped minimap2 settings, then two `copy_assign`
    runs; the measured wall time and peak RSS are in figures/README.md ("Build order and measured cost").

    Caching: the simulation reruns when its stamp (`sim.params`: catalog, index, seed, the md5 of bench/sim.py and
    the minimap2 command) changes or its inputs are newer than `sim.bam`; a rerun reuses a mapped part only when
    `sim.py`'s part key (FASTQ md5, index path/size/mtime, minimap2 command and version) still matches. `force`
    deletes the parts too. Each copy_assign run executes with every RUSTLE_* variable unset and records that in
    `<tag>.cmd`; the run reruns when its command changes."""
    d = figlib.work_dir(cfg, "o2sim") / species
    d.mkdir(parents=True, exist_ok=True)
    cat_tsv, cat_fa = cfg[f"{species}_o2_catalog_tsv"], cfg[f"{species}_o2_catalog_fa"]
    rec_dir = cfg.get("o2sim_dir")
    if rec_dir and Path(cat_tsv).resolve().parent == Path(rec_dir).resolve():
        # the recorded run directory holds DERIVED catalogs (<prefix>_cat.*); the simulation must start from the
        # catalog itself (gorilla: /mnt/linuxdisk/tmp/hom_c234.copies.*, byte-identical to g44_cat.*)
        print(f"[_o2] WARNING: {species}_o2_catalog_tsv = {cat_tsv} points into the recorded run directory "
              f"{rec_dir}; point it at the source catalog", file=sys.stderr)
    index, fasta = cfg[f"{species}_splice_mmi"], cfg[f"{species}_fasta"]
    seed, threads = cfg["o2_sim_seed"], cfg.get("threads", "4")
    binary = Path(cfg["bin"]) / "copy_assign"
    sim = d / "sim"
    bam = Path(f"{sim}.bam")
    stamp = d / "sim.params"
    notes = []
    base = f"catalog_tsv={cat_tsv}\ncatalog_fa={cat_fa}\nindex={index}\nseed={seed}\n"
    params = base + f"sim_py_md5={_md5(SIM)}\nminimap2={_sim.MM2}\n"
    complete = figlib.fresh(bam, cat_tsv, cat_fa, index) and Path(f"{sim}.copies_used.tsv").exists()
    # a stamp written before 2026-09-25's fingerprint lines holds `base` only: its simulation is accepted (the
    # recorded runs are not re-simulated to fill in a fingerprint) and the table notes say so
    legacy = not force and complete and _stamp_ok(stamp, base)
    if legacy:
        notes.append(f"simulation stamp {stamp} predates the simulator fingerprint (bench/sim.py md5, minimap2 "
                     f"command): accepted as built; bench/sim.py md5 now {_md5(SIM)}")
    elif force or not (complete and _stamp_ok(stamp, params)):
        # mapping in read-disjoint parts (identical records): o2_sim_parts_per_call > 0 bounds one call's wall time
        # (each part ~2-4 min incl. the ~1 min index load); re-run make.py data fig4 until the simulation completes
        stamp.unlink(missing_ok=True)
        Path(f"{sim}.sibling.tsv").unlink(missing_ok=True)  # written last: its presence marks a complete simulation
        if force:
            for q in list(d.glob("sim.part*.bam")) + [Path(f"{sim}.parts.key")]:
                q.unlink(missing_ok=True)
        figlib.run(["python3", str(SIM), "copies", cat_tsv, cat_fa, index, str(sim), str(seed), "--threads", threads,
                    "--parts", cfg.get("o2_sim_parts", "8"),
                    "--max-parts-per-call", cfg.get("o2_sim_parts_per_call", "0")],
                   log=d / "sim.log", cwd=d)
        if not Path(f"{sim}.sibling.tsv").exists():
            raise Pending(f"{species}: simulation mapping parts remain ({d})")
        stamp.write_text(params)
    _check_sim_bam(sim, notes)
    cat = derive_catalog(cat_tsv, cat_fa, sim, d / "cat", force=force)
    # every contig of the BAM header, one region each (what `tools/rustle_pipeline.sh assign` does; the recorded
    # runs used the same whole-contig lists hsa.regions.txt / ggo.regions.txt)
    regions = d / "regions.txt"
    if force or not figlib.fresh(regions, bam):
        _write_if_changed(regions, _regions_text(bam_contigs(bam)))
    for tag, extra in (("o2", []), ("u2", ["--union-certificate"])):
        out = d / tag
        cmd = [str(binary), "--bam", str(bam), "--fasta", fasta, "--regions", str(regions), "--families",
               str(cat["tsv"]), "--copies-fa", str(cat["fa"]), *extra, "--out", str(out)]
        cstamp = d / f"{tag}.cmd"
        want = " ".join(cmd) + f"\n{CLEAN_ENV}\n"
        built = figlib.fresh(Path(f"{out}.assignments.tsv"), bam, cat["tsv"], cat["fa"], regions, binary)
        # a run from before the RUSTLE_* guard (no .cmd stamp) is kept, and its environment is reported as unrecorded
        if force or not (built and (_stamp_ok(cstamp, want) or not cstamp.exists())):
            cstamp.unlink(missing_ok=True)
            figlib.run([*_clean_env_prefix(), *cmd], log=d / f"{tag}.log", cwd=d)
            cstamp.write_text(want)
        notes.append(f"copy_assign {tag}: " + (CLEAN_ENV + f" (stamp {cstamp.name})" if _stamp_ok(cstamp, want) else
                                               "run predates the RUSTLE_* guard; its environment is not recorded "
                                               "(the next forced rebuild runs it with RUSTLE_* unset)"))
    return {"species": species, "sample": DEV_SAMPLE[species], "catalog_scope": DEV_CATALOG[species], "sim": sim,
            "catalog": cat["tsv"], "source_catalog": cat_tsv, "o2": d / "o2", "u2": d / "u2",
            "source": f"rebuilt {d} (sim.py copies seed {seed}; copy_assign {binary})", "notes": notes}


def _check_sim_bam(sim: Path, notes: list):
    """The merged BAM holds exactly one primary-or-unmapped record per simulated read."""
    bam = Path(f"{sim}.bam")
    n_sim = sum(int(r["n_sim"]) for r in csv.DictReader(open(f"{sim}.copies_used.tsv"), delimiter="\t"))
    n_bam = int(subprocess.run(["samtools", "view", "-c", "-F", "2304", str(bam)], capture_output=True, text=True,
                               check=True).stdout.strip())
    if n_bam != n_sim:
        raise RuntimeError(f"{bam}: {n_bam} primary/unmapped records, but copies_used.tsv simulated {n_sim} reads")
    notes.append(f"{bam.name}: {n_bam} primary/unmapped records = the {n_sim} simulated reads of copies_used.tsv")


def _ensure_runs_genome(cfg: dict, sid: str, *, force=False, budget: Budget) -> dict:
    """Genome-wide simulation plus both copy_assign runs of one sample under `${work}/o2sim/<sample>/` (cached,
    bounded: at most `budget` heavy steps per call, then Pending).

    1. `sim.py copies ... --parts P --max-parts-per-call 1 --reuse-fastq --sibling none` (the first call adds
       `--simulate-only`: simulating ~1 M reads takes ~5 min by itself). P = planned_parts(). Complete when
       `sim.done` exists and `sim.params` (catalog, index, seed, parts, bench/sim.py md5, minimap2 command) matches.
    2. derive_catalog(), regions = every contig of the BAM header (as the dev runs).
    3. copy_assign default and --union-certificate in shards (plan_shards; cfg `o2_shard` 0 = one unsharded run,
       `o2_shard_union` 0 = the union run unsharded), merged into o2.assignments.tsv / u2.assignments.tsv /
       u2.union_certificate.tsv (rows of the shards in shard order; the scorer does not depend on row order)."""
    row = samples.get(cfg, sid)
    species = row["species"]
    d = genome_work_dir(cfg, sid)
    d.mkdir(parents=True, exist_ok=True)
    cat_tsv, cat_fa = catalog_paths(cfg, sid)
    index, fasta = splice_index(cfg, sid), row["fasta"]
    seed, threads = cfg["o2_sim_seed"], cfg.get("threads", "4")
    sim = d / "sim"
    bam = Path(f"{sim}.bam")
    n_expected, per_part, parts = planned_parts(cfg, species, cat_tsv, cat_fa)
    stamp = d / "sim.params"
    params = (f"catalog_tsv={cat_tsv}\ncatalog_fa={cat_fa}\nindex={index}\nseed={seed}\nparts={parts}\n"
              f"reads_per_part={per_part}\nsibling=none\nsim_py_md5={_md5(SIM)}\nminimap2={_sim.MM2}\n")
    notes = [f"simulation: {n_expected} reads planned before simulating, {parts} mapping parts of <= {per_part} reads"]
    done = Path(f"{sim}.done")
    complete = (done.exists() and figlib.fresh(bam, cat_tsv, cat_fa, index)
                and Path(f"{sim}.copies_used.tsv").exists() and _stamp_ok(stamp, params))
    if force or not complete:
        if force:
            for q in list(d.glob("sim.part*")) + [Path(f"{sim}.parts.key"), Path(f"{sim}.fq.key"), done]:
                q.unlink(missing_ok=True)
        stamp.unlink(missing_ok=True)
        done.unlink(missing_ok=True)
        fq_key = _sim._fastq_key(SimpleNamespace(copies_tsv=str(cat_tsv), copies_fa=str(cat_fa), seed=int(seed)))
        first = not _stamp_ok(Path(f"{sim}.fq.key"), fq_key)
        budget.take(f"{sid}: simulate {n_expected} reads" if first else f"{sid}: map one simulation part")
        cmd = ["python3", str(SIM), "copies", str(cat_tsv), str(cat_fa), str(index), str(sim), str(seed),
               "--threads", threads, "--parts", str(parts), "--max-parts-per-call", "1", "--reuse-fastq",
               "--sibling", "none"] + (["--simulate-only"] if first else [])
        figlib.run(cmd, log=d / "sim.log", cwd=d)
        if not done.exists():
            left = sum(not Path(f"{sim}.part{i}.bam").exists() for i in range(parts)) if parts > 1 else 1
            raise Pending(f"{sid}: {left} of {parts} simulation mapping parts remain ({d})")
        stamp.write_text(params)
    _check_sim_bam(sim, notes)
    cat = derive_catalog(cat_tsv, cat_fa, sim, d / "cat", force=force)
    regions = d / "regions.txt"
    contigs = bam_contigs(bam)
    if force or not figlib.fresh(regions, bam):
        _write_if_changed(regions, _regions_text(contigs))
    binary = Path(cfg["bin"]) / "copy_assign"
    notes.append(f"copy table: {copy_table(cfg)} ({cat_tsv})")
    runs = {"species": species, "sample": sid, "catalog_scope": "genome", "sim": sim, "catalog": cat["tsv"],
            "source_catalog": cat_tsv, "copy_table": copy_table(cfg), "o2": d / "o2", "u2": d / "u2", "notes": notes,
            "source": f"rebuilt {d} (sim.py copies seed {seed}, {parts} parts; copy_assign {binary})"}
    shard = cfg.get("o2_shard", "1") != "0"
    plan = None
    if shard:
        plan = plan_shards(bam, cat["tsv"], cat["fa"], contigs, d / "shards", kind="sim",
                           budget_s=float(cfg.get("o2_shard_budget_s", SHARD_BUDGET_S)))
        notes.append(plan["summary"])
    for tag, extra in (("o2", []), ("u2", ["--union-certificate"])):
        sharded = shard and not (tag == "u2" and cfg.get("o2_shard_union", "1") == "0")
        units = plan["shards"] if sharded else [{"dir": str(d), "tsv": str(cat["tsv"]), "fa": str(cat["fa"]),
                                                 "regions": str(regions), "name": "all"}]
        outs = []
        for u in units:
            out = Path(u["dir"]) / tag
            cmd = [str(binary), "--bam", str(bam), "--fasta", fasta, "--regions", u["regions"], "--families",
                   u["tsv"], "--copies-fa", u["fa"], *extra, "--out", str(out)]
            _run_copy_assign(cmd, out, [bam, u["tsv"], u["fa"], u["regions"], binary], budget,
                             f"{sid}: copy_assign {tag} shard {u['name']}", force=force)
            outs.append(out)
        if sharded:
            merge_outputs(outs, d / tag, union=(tag == "u2"))
        notes.append(f"copy_assign {tag}: {CLEAN_ENV}; " + (f"{len(outs)} shards merged (row order = shard order)"
                                                             if sharded else "one unsharded run"))
    return runs


def _run_copy_assign(cmd: list, out: Path, srcs: list, budget: Budget, what: str, *, force=False):
    """One copy_assign run with every RUSTLE_* unset, skipped when its products are fresh and its `.cmd` stamp
    matches; one heavy step."""
    cstamp = Path(f"{out}.cmd")
    want = " ".join(map(str, cmd)) + f"\n{CLEAN_ENV}\n"
    if not force and figlib.fresh(Path(f"{out}.assignments.tsv"), *srcs) and _stamp_ok(cstamp, want):
        return False
    budget.take(what)
    cstamp.unlink(missing_ok=True)
    timer = ["/usr/bin/time", "-f", "TIME elapsed=%e rss_kb=%M"] if os.path.exists("/usr/bin/time") else []
    figlib.run([*timer, *_clean_env_prefix(), *map(str, cmd)], log=Path(f"{out}.log"), cwd=out.parent)
    cstamp.write_text(want)
    return True


def merge_outputs(outs: list, merged: Path, *, union: bool):
    """Concatenate the shards' `assignments.tsv` (and `union_certificate.tsv`) under one header, in shard order."""
    for suffix in [".assignments.tsv"] + ([".union_certificate.tsv"] if union else []):
        header, body = None, []
        for o in outs:
            with open(f"{o}{suffix}") as fh:
                h = fh.readline()
                if header is None:
                    header = h
                elif h != header:
                    raise RuntimeError(f"{o}{suffix}: header differs from the first shard's")
                body.extend(l for l in fh if l.strip())
        Path(merged).parent.mkdir(parents=True, exist_ok=True)
        _write_if_changed(Path(f"{merged}{suffix}"), header + "".join(body))


# ================================================================ shards: read-connected components of families
class _UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, i):
        p = self.p
        while p[i] != i:
            p[i] = p[p[i]]
            i = p[i]
        return i

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[max(a, b)] = min(a, b)


def _name_hash(name: bytes) -> int:
    return int.from_bytes(hashlib.blake2b(name, digest_size=8).digest(), "little")


def _read_catalog(cat_tsv):
    with open(cat_tsv) as fh:
        header = fh.readline()
        rows = [l for l in fh if l.strip()]
    fams, copies = [], collections.defaultdict(list)
    for line in rows:
        f = line.rstrip("\n").split("\t")
        if f[0] not in copies:
            fams.append(f[0])
        copies[f[0]].append((f[3], int(f[4]), int(f[5])))
    return header, rows, fams, copies


def link_pieces(bam, clusters: dict, piece_dir: Path, *, deadline: float | None = None) -> tuple[list, list]:
    """Relation (b) of the shard plan: every mapped record (primary, secondary or supplementary) that overlaps a
    window cluster, as (64-bit hash of the read name, cluster id) pairs, plus the records per cluster. One cached
    piece per contig (`piece_dir/<contig>.npz`); with a `deadline` (time.time() value) no new contig starts after it
    and Pending is raised while pieces remain."""
    import numpy as np
    import pysam

    piece_dir.mkdir(parents=True, exist_ok=True)
    hashes, cids, counts = [], [], collections.Counter()
    todo = []
    for chrom, cl in clusters.items():
        p = piece_dir / f"{chrom}.npz"
        if p.exists():
            z = np.load(p)
            hashes.append(z["h"]); cids.append(z["c"])
            counts.update(dict(zip(z["cid"].tolist(), z["n"].tolist())))
        else:
            todo.append((chrom, cl, p))
    if todo:
        with pysam.AlignmentFile(str(bam)) as af:
            for chrom, cl, p in todo:
                if deadline is not None and time.time() > deadline:
                    raise Pending(f"record pass over {bam}: {len(todo)} contig pieces remain in {piece_dir}")
                his = [c[1] for c in cl]
                h_l, c_l, cnt = [], [], collections.Counter()
                for lo, hi, cid in cl:
                    for r in af.fetch(chrom, lo, hi):
                        if r.is_unmapped:
                            continue
                        rs, re_ = r.reference_start, r.reference_end or r.reference_start + 1
                        j = bisect.bisect_right(his, rs)
                        hv = None
                        while j < len(cl) and cl[j][0] < re_:
                            if hv is None:
                                hv = _name_hash(r.query_name.encode())
                            h_l.append(hv); c_l.append(cl[j][2]); cnt[cl[j][2]] += 1
                            j += 1
                h = np.array(h_l, dtype=np.uint64)
                c = np.array(c_l, dtype=np.uint32)
                ks = sorted(cnt)
                np.savez(p.with_suffix(".tmp.npz"), h=h, c=c, cid=np.array(ks, dtype=np.int64),
                         n=np.array([cnt[k] for k in ks], dtype=np.int64))
                p.with_suffix(".tmp.npz").replace(p)
                hashes.append(h); cids.append(c)
                counts.update(cnt)
    return hashes, cids, counts


def plan_shards(bam, cat_tsv, cat_fa, contigs: list, out_dir: Path, *, kind: str = "sim",
                budget_s: float = SHARD_BUDGET_S, n_shards: int | None = None, deadline: float | None = None,
                sample_frac: float | None = None, sample_seed: int = 20260925) -> dict:
    """Split the catalog into shards of whole READ-CONNECTED COMPONENTS of families, so that each shard's
    copy_assign run gives exactly the rows the one-run table gives for its families
    (docs/PREREG_genome_wide_copy_assignment_2026-09-25.md §4). Two families are joined when
      (a) their read windows (copy +/- COPY_READ_PAD, the records copy_assign loads) overlap or touch on a contig;
      (b) one read has records (any flag) inside a window of each (link_pieces over `bam`);
      (c) a family with copies on several contigs is joined, on each such contig that carries single-contig
          families, to the single-contig family nearest its copy (so no shard sweeps that contig whole).
    `contigs` = the one-run's regions [(name, length)]: contigs with no single-contig family are swept whole in
    every shard, as in the one run. Components are packed first-fit-decreasing into shards of estimated wall time
    <= budget_s (SHARD_COST[kind]: seconds per record inside windows, plus a fixed cost and the records of the
    whole-swept contigs); `n_shards` forces that many balanced shards (check V3). A component is never split; one
    above the budget is its own shard and is flagged. `sample_frac`: keep a seeded random fraction of the
    components, stratified by the contig of the component's first copy (stop rule of the pre-registration).
    Writes out_dir/shard_NNN/{cat.copies.tsv, cat.copies.fa, regions.txt}, components.tsv, plan.json; reused
    while plan.key (catalog md5s, BAM fingerprint, parameters, this module's md5) matches."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    st = os.stat(bam)
    key = json.dumps({"tsv": _md5(cat_tsv), "fa": _md5(cat_fa), "bam": [str(bam), st.st_size, int(st.st_mtime)],
                      "contigs": contigs, "kind": kind, "budget_s": budget_s, "n_shards": n_shards,
                      "pad": COPY_READ_PAD, "cost": SHARD_COST[kind], "sample_frac": sample_frac,
                      "sample_seed": sample_seed, "code": _md5(__file__)}, sort_keys=True)
    if _stamp_ok(out_dir / "plan.key", key) and (out_dir / "plan.json").exists():
        return json.loads((out_dir / "plan.json").read_text())
    header, rows, fams, copies = _read_catalog(cat_tsv)
    fidx = {f: i for i, f in enumerate(fams)}
    LN = dict(contigs)
    single = [len({c for c, _, _ in copies[f]}) == 1 for f in fams]
    single_contigs = {copies[f][0][0] for i, f in enumerate(fams) if single[i]}
    uf = _UF(len(fams))
    # (a) window clusters per contig
    wins = collections.defaultdict(list)
    for i, f in enumerate(fams):
        for c, s, e in copies[f]:
            if c not in LN:
                raise RuntimeError(f"{cat_tsv}: {f} has a copy on {c}, which is not among the regions")
            wins[c].append((max(0, s - COPY_READ_PAD), min(LN[c], e + COPY_READ_PAD), i))
    clusters, cl_members = {}, {}
    nid = 0
    for c in sorted(wins):
        cl = []
        for lo, hi, i in sorted(wins[c]):
            if cl and lo <= cl[-1][1]:
                cl[-1][1] = max(cl[-1][1], hi)
                uf.union(cl_members[cl[-1][2]][0], i)
                cl_members[cl[-1][2]].append(i)
            else:
                cl.append([lo, hi, nid])
                cl_members[nid] = [i]
                nid += 1
        clusters[c] = [tuple(x) for x in cl]
    # (c) anchors for families on several contigs
    for i, f in enumerate(fams):
        if single[i]:
            continue
        for c in sorted({c for c, _, _ in copies[f]}):
            if c not in single_contigs:
                continue
            s0 = min(s for cc, s, _ in copies[f] if cc == c)
            best = None
            for lo, hi, cid in clusters[c]:
                sm = [j for j in cl_members[cid] if single[j]]
                if not sm:
                    continue
                gap = 0 if lo <= s0 < hi else min(abs(s0 - lo), abs(s0 - hi))
                if best is None or gap < best[0]:
                    best = (gap, sm[0])
            uf.union(i, best[1])
    # (b) reads
    hashes, cids, counts = link_pieces(bam, clusters, out_dir / "links", deadline=deadline)
    rep = {cid: m[0] for cid, m in cl_members.items()}
    import numpy as np
    if hashes and sum(len(h) for h in hashes):
        h = np.concatenate(hashes)
        c = np.concatenate(cids)
        o = np.lexsort((c, h))
        h, c = h[o], c[o]
        brk = np.flatnonzero(np.diff(h) != 0) + 1
        starts = np.concatenate(([0], brk))
        ends = np.concatenate((brk, [len(h)]))
        multi = np.flatnonzero(c[starts] != c[ends - 1])   # sorted by cluster within a read: first != last
        for g in multi.tolist():
            cs = np.unique(c[starts[g]:ends[g]]).tolist()
            for x in cs[1:]:
                uf.union(rep[cs[0]], rep[x])
    # components
    comp = collections.defaultdict(list)
    for i in range(len(fams)):
        comp[uf.find(i)].append(i)
    comp_records = collections.Counter()
    for cid, n in counts.items():
        comp_records[uf.find(rep[cid])] += n
    whole = [c for c, _ in contigs if c not in single_contigs]   # swept whole in every shard, as in the one run
    whole_records = 0
    if whole:
        idx = {l.split("\t")[0]: int(l.split("\t")[2]) for l in subprocess.run(
            ["samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True).stdout.splitlines()}
        whole_records = sum(idx.get(c, 0) for c in whole)
    per_rec, fixed, per_whole = SHARD_COST[kind]
    fixed_s = fixed + per_whole * whole_records
    comps = sorted(comp.items(), key=lambda kv: (-comp_records[kv[0]], kv[1][0]))
    sampled_note = ""
    if sample_frac is not None:
        strata = collections.defaultdict(list)
        for root, members in comps:
            cs = {copies[fams[j]][0][0] for j in members}
            strata[copies[fams[members[0]]][0][0] if len(cs) == 1 else "multi-contig"].append((root, members))
        rng = random.Random(sample_seed)
        keep = []
        for k in sorted(strata):
            v = sorted(strata[k], key=lambda kv: kv[1][0])
            keep += rng.sample(v, max(1, math.ceil(sample_frac * len(v))))
        kept_roots = {r for r, _ in keep}
        sampled_note = (f"; SAMPLED {len(keep)} of {len(comps)} components (fraction {sample_frac}, seed "
                        f"{sample_seed}, stratified by contig)")
        comps = [kv for kv in comps if kv[0] in kept_roots]
    est = {root: per_rec * comp_records[root] for root, _ in comps}
    bins: list[list] = []  # [est_s, [roots]]
    if n_shards:
        bins = [[0.0, []] for _ in range(n_shards)]
        for root, _ in comps:
            b = min(bins, key=lambda x: x[0])
            b[0] += est[root]
            b[1].append(root)
    else:
        for root, _ in comps:
            for b in bins:
                if b[0] + est[root] + fixed_s <= budget_s:
                    b[0] += est[root]
                    b[1].append(root)
                    break
            else:
                bins.append([est[root], [root]])
    members_of = dict(comps)
    shards = []
    fa_recs, fa_order = _fa_records(cat_fa)
    for k, (e_s, roots) in enumerate(bins):
        if not roots:
            continue
        sd = out_dir / f"shard_{k:03d}"
        sd.mkdir(exist_ok=True)
        fset = {fams[j] for r in roots for j in members_of[r]}
        rows_k = [l for l in rows if l.split("\t", 1)[0] in fset]
        _write_if_changed(sd / "cat.copies.tsv", header + "".join(rows_k))
        _write_if_changed(sd / "cat.copies.fa", "".join("".join(fa_recs[q]) for q in fa_order if q[0] in fset))
        cset = {copies[f][0][0] for f in fset if single[fidx[f]]} | set(whole)
        _write_if_changed(sd / "regions.txt", _regions_text([(c, n) for c, n in contigs if c in cset]))
        shards.append({"name": f"{k:03d}", "dir": str(sd), "tsv": str(sd / "cat.copies.tsv"),
                       "fa": str(sd / "cat.copies.fa"), "regions": str(sd / "regions.txt"),
                       "n_families": len(fset), "n_copies": len(rows_k), "n_components": len(roots),
                       "records": sum(comp_records[r] for r in roots), "est_s": round(e_s + fixed_s, 1),
                       "over_budget": e_s + fixed_s > budget_s})
    with open(out_dir / "components.tsv", "w") as fh:
        fh.write("component\tn_families\tfamilies\trecords\n")
        for root, members in comps:
            fh.write(f"{fams[root]}\t{len(members)}\t{','.join(fams[j] for j in members)}\t{comp_records[root]}\n")
    n_over = sum(s["over_budget"] for s in shards)
    summary = (f"shard plan ({kind}): {len(fams)} families in {len(comp)} read-connected components -> "
               f"{len(shards)} shards; largest component {max((len(m) for _, m in comps), default=0)} families / "
               f"{max(comp_records.values(), default=0)} records; {len(whole)} contigs swept whole in every shard "
               f"({whole_records} records); estimated {sum(s['est_s'] for s in shards):.0f} s in total"
               + (f"; {n_over} shard(s) above the {budget_s:.0f}-s budget" if n_over else "") + sampled_note)
    plan = {"summary": summary, "shards": shards, "whole_contigs": whole, "n_components": len(comp),
            "sampled": sample_frac}
    (out_dir / "plan.json").write_text(json.dumps(plan, indent=1))
    (out_dir / "plan.key").write_text(key)
    print(f"[_o2] {summary}", file=sys.stderr)
    return plan


# ================================================================ scoring
def score_logs(runs: dict, log_dir: Path) -> dict:
    """`bench/score.py reads --per-read` on both assignment runs (light: one -F 2308 pass over the simulated BAM each).
    Returns {tag: log} and {tag + "_per_read": per-read TSV}: the log is the scorer's table, the TSV its per-read
    verdicts, which `per_read()` joins instead of re-implementing the scorer's judge()."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    key = runs.get("sample") or runs["species"]
    for tag in ("o2", "u2"):
        log = log_dir / f"score_reads_{runs['species']}_{tag}.txt"
        pr = log_dir / f"per_read_{runs['species']}_{tag}.tsv"
        if runs.get("catalog_scope") == "genome":
            log, pr = log_dir / f"score_reads_{key}_{tag}.txt", log_dir / f"per_read_{key}_{tag}.tsv"
        src = (f"{runs['sim']}.bam", f"{runs[tag]}.assignments.tsv", runs["catalog"], SCORE)
        if not (figlib.fresh(log, *src) and figlib.fresh(pr, *src)):
            figlib.run([sys.executable, str(SCORE), "reads", "--catalog", str(runs["catalog"]), str(runs["sim"]),
                        str(runs[tag]), "--per-read", str(pr)], log=log)
        out[tag] = log
        out[f"{tag}_per_read"] = pr
    return out


# ================================================================ per-read join
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")
# score.py reads' verdicts -> this module's outcome names ("lost" = no row at all, "no_own_row" = no row of the
# read's source family: both are score.py's "other" column and this module's "not_scored")
_VERDICT = {"lost": "not_scored", "no_own_row": "not_scored", "no_primary_row": "not_scored"}


def _ref_span(pos1: int, cigar: str):
    s = pos1 - 1
    return s, s + sum(int(n) for n, op in _CIG.findall(cigar) if op in "MDN=X")


def _load_per_read(path) -> dict:
    with open(path) as fh:
        return {r["read_name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def placements(bam, names) -> dict:
    """Every mapped, non-supplementary alignment (-F 2052: the primary and every secondary) of the named reads:
    name -> [(chrom, start, end, AS, NM)] (0-based half-open reference span)."""
    out = collections.defaultdict(list)
    p = subprocess.Popen(["samtools", "view", "-F", "2052", str(bam)], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        if line.split("\t", 1)[0] not in names:
            continue
        f = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in f[11:]}
        s, e = _ref_span(int(f[3]), f[5])
        out[f[0]].append((f[2], s, e, int(tags["AS"]), int(tags["NM"])))
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam}")
    return out


def twin_state(pls: list, source) -> str:
    """Why a MAPQ-0 read cannot be told apart genome-wide (TWIN_STATES), from all its placements `pls` and the source
    copy's catalog span `source` = (chrom, start, end). A source-copy placement overlaps the source copy; the rest
    are 'elsewhere'. best = the highest AS over all placements; the source-copy placement is the source-copy
    alignment with the highest AS (then the lowest NM).
      nm_twin               a placement elsewhere has AS = best and the source-copy placement's NM (an NM-identical
                            twin: the read is equally far from both loci)
      as_tie_nm_differs     placements elsewhere reach AS = best, none with that NM
      no_equal_as_elsewhere the source copy holds the best AS alone (MAPQ 0 from the aligner's chaining score)
      true_below_best       the best AS is elsewhere; the source-copy placement scores lower
      no_true_placement     no alignment overlaps the source copy"""
    true = [x for x in pls if source and x[0] == source[0] and min(x[2], source[2]) - max(x[1], source[1]) > 0]
    if not true:
        return "no_true_placement"
    best = max(x[3] for x in pls)
    t_as = max(x[3] for x in true)
    if t_as < best:
        return "true_below_best"
    t_nm = min(x[4] for x in true if x[3] == t_as)
    tied = [x for x in pls if x not in true and x[3] == best]
    if any(x[4] == t_nm for x in tied):
        return "nm_twin"
    return "as_tie_nm_differs" if tied else "no_equal_as_elsewhere"


def per_read(runs: dict, logs: dict) -> list[dict]:
    """One record per simulated read (primary or unmapped record of the simulated BAM).

    The aligner's placement and the identity band are computed here; everything about the assignment tables comes
    from `score.py reads --per-read` (`logs[o2_per_read]`, `logs[u2_per_read]`), which lists every MAPQ-0 read.
    The test fields are left empty for MAPQ > 0 reads: the AS-tied gate does not score them and no table uses them.
    `twin` = twin_state() for the MAPQ-0 reads ("na" otherwise)."""
    species, sample = runs["species"], runs.get("sample") or DEV_SAMPLE.get(runs["species"], runs["species"])
    cat, ident = {}, {}
    for r in csv.DictReader(open(runs["catalog"]), delimiter="\t"):
        k = (r["family_id"], r["copy_idx"])
        cat[k] = (r["chrom"], int(r["start"]), int(r["end"]))
        v = r.get("max_family_identity")
        ident[k] = float(v) if v not in (None, "", "NA") else float("nan")
    # the copy a primary falls on: largest raw overlap, ties -> first in catalog order (score.CopyIndex: the same rule
    # as the former scan of every copy of the contig, indexed so it stays linear genome-wide)
    idx = _score.CopyIndex(cat)
    pr_o2, pr_u2 = _load_per_read(logs["o2_per_read"]), _load_per_read(logs["u2_per_read"])
    pls = placements(f"{runs['sim']}.bam", set(pr_o2))  # pr_o2 lists every MAPQ-0 read
    out = []
    p = subprocess.Popen(["samtools", "view", "-F", "2304", f"{runs['sim']}.bam"], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        name, flag = f[0], int(f[1])
        t = tuple(name.split("|")[:2])
        rec = {"sample": sample, "species": species, "read": name, "family": t[0], "copy": t[1],
               "identity": ident.get(t, float("nan"))}
        rec["band"] = figlib.identity_band(rec["identity"]) or "unknown"
        # aligner primary: the catalog copy with the max raw overlap of the primary's reference span
        if flag & 4:
            rec["mapq"] = None
            rec["aligner"] = "unmapped"
        else:
            rec["mapq"] = int(f[4])
            s, e = _ref_span(int(f[3]), f[5])
            best = idx.copy_of(f[2], s, e)
            rec["aligner"] = ("outside_catalog" if best is None else
                              "true_copy" if (best == t or _score.same_locus(cat[best], cat.get(t))) else "other_copy")
        o2, u2 = pr_o2.get(name), pr_u2.get(name)
        if rec["mapq"] == 0 and (o2 is None or u2 is None):
            raise RuntimeError(f"{name}: MAPQ-0 read missing from score.py reads --per-read")
        if o2 is None:  # MAPQ > 0 or unmapped: not scored by the test
            o2 = {"n_rows": "0", "own_status": "no_row", "own_n_decisive": "", "OWN": "lost", "ANY": "lost",
                  "wrong_locus_rows": "0"}
            u2 = {"own_status": "no_row", "own_n_decisive": "", "ANY": "lost"}
        rec["in_table"] = int(o2["n_rows"]) > 0
        rec["own_row"] = o2["own_status"] != "no_row"
        rec["own_status"] = o2["own_status"] if rec["own_row"] else ""
        rec["family_psv"] = rec["own_row"] and int(o2["own_n_decisive"]) >= 1
        rec["own_outcome"] = _VERDICT.get(o2["OWN"], o2["OWN"])
        rec["assigned"] = rec["own_outcome"] in ("correct", "wrong", "conflict")
        rec["correct"] = rec["own_outcome"] == "correct"
        rec["any_outcome"] = _VERDICT.get(o2["ANY"], o2["ANY"])
        rec["foreign_claim"] = int(o2["wrong_locus_rows"]) > 0
        rec["identical"] = u2["own_status"] != "no_row" and int(u2["own_n_decisive"]) == 0
        rec["union_outcome"] = _VERDICT.get(u2["ANY"], u2["ANY"])
        rec["twin"] = twin_state(pls.get(name, []), cat.get(t)) if rec["mapq"] == 0 else "na"
        out.append(rec)
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(f"samtools view failed on {runs['sim']}.bam")
    return out


# ================================================================ cross-check against bench/score.py reads
_ROW = re.compile(r"^\s+ALL\s+n=\s*(\d+) correct\s+(\d+) wrong\s+(\d+) conflict\s+(\d+) abstain\s+(\d+) other\s+(\d+)")


def parse_score_reads(log) -> dict:
    """{'OWN'|'PRIMARY'|'ANY': (n, correct, wrong, conflict, abstain, other)} from the ALL rows of `score.py reads`."""
    out, view = {}, None
    for line in open(log):
        if line.startswith("== "):
            view = line[3:].strip()
            continue
        m = _ROW.match(line)
        if m and view:
            out[view] = tuple(int(x) for x in m.groups())
    return out


def crosscheck(reads: list[dict], logs: dict) -> list[str]:
    """Our MAPQ-0 tabulation must equal score.py reads' ALL rows (OWN on the default run, ANY on both)."""
    z = [r for r in reads if r["mapq"] == 0]

    def tally(key):
        c = collections.Counter(r[key] for r in z)
        return (len(z), c["correct"], c["wrong"], c["conflict"], c["abstain"], c["not_scored"])
    pairs = [("o2", "OWN", "own_outcome"), ("o2", "ANY", "any_outcome"), ("u2", "ANY", "union_outcome")]
    notes = []
    for tag, view, key in pairs:
        ref = parse_score_reads(logs[tag]).get(view)
        mine = tally(key)
        if ref is None and mine[0] == 0:  # score.py prints no ALL row when no read has MAPQ 0
            notes.append(f"cross-check: no MAPQ-0 read ({tag} {view}); score.py reads prints no ALL row")
            continue
        if ref != mine:
            raise RuntimeError(f"tabulation disagrees with score.py reads ({tag} {view}): score.py {ref} vs {mine}")
        notes.append(f"cross-check: score.py reads {tag} {view} ALL (n, correct, wrong, conflict, abstain, other) "
                     f"= {ref}, identical to this table's tally")
    return notes


# ================================================================ Figure 4 sets (UpSet) and per-read fates
SETS = ["unique", "mapq0", "not_scored", "identical", "family_psv", "assigned", "correct"]
SET_LABEL = {
    "unique": "MAPQ > 0",
    "mapq0": "MAPQ 0",
    "not_scored": "Not tested (no result for the read)",
    "identical": "No decisive site over all placements",
    "family_psv": "≥ 1 decisive site in the source family",
    "assigned": "Assigned in the source family*",
    "correct": "… to the source copy",
}
ALIGNER_STATES = ["true_copy", "other_copy", "outside_catalog", "unmapped"]
# twin_state() of a MAPQ-0 read (MAPQ > 0 and unmapped reads: "na")
TWIN_STATES = ["nm_twin", "as_tie_nm_differs", "no_equal_as_elsewhere", "true_below_best", "no_true_placement"]
# per-read fate of the MAPQ-0 reads, truth-selected reading (PREREG_genome_wide_copy_assignment §2), in bar order
FATES = ["correct", "wrong", "psv_unassigned", "no_site", "no_result"]
FATE_LABEL = {
    "correct": "Assigned to the source copy",
    "wrong": "Assigned to another copy",
    "psv_unassigned": "≥ 1 decisive site, left unassigned",
    "no_site": "No decisive site, left unassigned",
    "no_result": "No result in the source family",
}


def fate_of(r: dict) -> str | None:
    """The fate of one UpSet-table row's reads (None for MAPQ > 0 and unmapped reads). Tables built before
    2026-09-25's `own_row` column cannot tell 'no result in the source family' from 'no decisive site' when the read
    has rows in other families only; they count such reads under 'no_site' unless `not_scored` (no row at all)."""
    if r["mapq0"] != "1":
        return None
    if r["assigned"] == "1":
        return "correct" if r["correct"] == "1" else "wrong"
    own_row = r.get("own_row")
    if r["not_scored"] == "1" or own_row == "0":
        return "no_result"
    return "psv_unassigned" if r["family_psv"] == "1" else "no_site"


# ================================================================ tabulation
UPSET_HEADER = ["sample", "species", "catalog_scope", "unique", "mapq0", "not_scored", "identical", "family_psv",
                "assigned", "correct", "own_row", "foreign_claim", "any_outcome", "union_outcome", "twin",
                "aligner_primary", "identity_band", "n_reads"]


def upset_rows(reads: list[dict], catalog_scope: str) -> list[list]:
    c = collections.Counter()
    for r in reads:
        mq = r["mapq"]
        z = mq == 0
        key = (r["sample"], r["species"], catalog_scope, int(mq is not None and mq > 0), int(z),
               int(z and not r["in_table"]), int(z and r["identical"]), int(z and r["family_psv"]),
               int(z and r["assigned"]), int(z and r["correct"]), int(z and r["own_row"]),
               int(z and r["foreign_claim"]), r["any_outcome"] if z else "na", r["union_outcome"] if z else "na",
               r["twin"], r["aligner"], r["band"])
        c[key] += 1
    band_order = {b: i for i, b in enumerate(figlib.IDENTITY_BANDS + ["unknown"])}
    return [list(k) + [v] for k, v in sorted(c.items(), key=lambda kv: (kv[0][:-1], band_order[kv[0][-1]]))]


def table_sample(r: dict) -> str:
    """The sample of a table row; tables built before the `sample` column carry `species` only (dev scope)."""
    return r.get("sample") or DEV_SAMPLE.get(r["species"], r["species"])


def table_scope(r: dict) -> str:
    return r.get("catalog_scope") or DEV_CATALOG.get(r["species"], "")


def catalog_title(sample: str, scope: str, short: bool = False) -> str:
    """'Human A119b (CHM13 v2.0), genome-wide catalog' / '..., chr16 catalog (development)'; `short` drops the
    exposure note (the caption states it)."""
    base = f"{SAMPLE_LABEL.get(sample, sample)} ({SAMPLE_DETAIL.get(sample, '')})"
    if scope == "genome":
        return f"{base}, genome-wide catalog"
    if short:
        return f"{base}, {scope} catalog"
    note = ("development" if sample == "human_A119b" and scope == "chr16"
            else "not used to develop the copy-assignment rules")
    return f"{base}, {scope} catalog ({note})"


BAND_HEADER = ["sample", "species", "catalog_scope", "method", "identity_band", "n_reads", "correct", "wrong",
               "conflict", "abstain", "not_scored", "coverage", "accuracy", "n_copies", "n_copies_assigned", "acc_lo",
               "acc_hi", "acc_ci", "n_family_psv", "n_nm_twin"]
METHODS = ["family_certificate", "per_family_table", "union_certificate", "aligner_primary"]
METHOD_LABEL = {
    "family_certificate": "Copy-assignment test, scored in the source family*",
    "per_family_table": "Default output, any assigned result",
    "union_certificate": "Union test (one test over all placements)",
    "aligner_primary": "Aligner's primary alignment",
}
BOOT_B, BOOT_SEED = 2000, 20260925
MIN_CI_COPIES = 5  # fewer assigned source copies than this: no interval is drawn (acc_ci = too_few_copies)


def copy_level_ci(units: list[tuple[int, int]], B: int = BOOT_B, seed: int = BOOT_SEED):
    """95% interval for accuracy = sum(correct) / sum(assigned) when reads are clustered by source copy.

    units: one (correct, assigned) pair per source copy with >= 1 assigned read. Reads from one copy (10-100 per
    copy) are not independent, so a read-level Wilson interval is too narrow. With fewer than MIN_CI_COPIES copies
    no interval is given ((None, None, "too_few_copies")): a bootstrap over 2-3 copies is a point or a coin toss.
    Otherwise cluster bootstrap: resample the copies with replacement B times (fixed seed), percentile 2.5-97.5 of
    the pooled ratio. When every copy has the same ratio (e.g. every assigned read correct) the bootstrap is
    degenerate (a point), so the interval is the Wilson interval at the pooled ratio with the m COPIES as the sample
    size — the copy as the independent unit.
    Returns (lo, hi, method) or None when nothing is assigned."""
    m = len(units)
    if m == 0:
        return None
    if m < MIN_CI_COPIES:
        return None, None, "too_few_copies"
    k_all, n_all = sum(k for k, _ in units), sum(n for _, n in units)
    if len({k / n for k, n in units}) == 1:
        lo, hi = wilson(m * k_all / n_all, m)
        return lo, hi, "wilson_over_copies"
    rng = random.Random(seed)
    vals = []
    for _ in range(B):
        draw = rng.choices(units, k=m)
        vals.append(sum(k for k, _ in draw) / sum(n for _, n in draw))
    vals.sort()
    return vals[round(0.025 * (B - 1))], vals[round(0.975 * (B - 1))], "cluster_bootstrap"


def band_rows(reads: list[dict], catalog_scope: str) -> list[list]:
    """Per method x identity band, over the MAPQ-0 reads: coverage = (correct + wrong + conflict) / n,
    accuracy = correct / (correct + wrong + conflict) (score.py reads' convention; a conflict = two loci claimed).
    Aligner: correct = primary on the source copy's locus, wrong = on another catalog copy, not_scored = primary
    outside every catalog copy (the pre-registered baseline counts only placements on a catalog copy).
    n_copies = distinct source copies of the band's MAPQ-0 reads; n_copies_assigned = those with >= 1 assigned
    read; acc_lo / acc_hi / acc_ci = copy_level_ci() over the assigned reads grouped by source copy.
    n_family_psv (family_certificate rows only) = MAPQ-0 reads whose source family's result has >= 1 decisive site,
    the most that reading can assign by a decisive site; n_nm_twin (union_certificate rows only) = MAPQ-0 reads with
    an NM-identical twin elsewhere (twin_state), which no test can separate from their source."""
    z = [r for r in reads if r["mapq"] == 0]
    key_of = {"family_certificate": "own_outcome", "per_family_table": "any_outcome",
              "union_certificate": "union_outcome"}
    aligner_map = {"true_copy": "correct", "other_copy": "wrong", "outside_catalog": "not_scored",
                   "unmapped": "not_scored"}
    out = []
    sample = reads[0]["sample"] if reads else ""
    species = reads[0]["species"] if reads else ""
    for m in METHODS:
        for band in ["all"] + figlib.IDENTITY_BANDS:
            sel = [r for r in z if band == "all" or r["band"] == band]
            if m == "aligner_primary":
                outcome = [aligner_map[r["aligner"]] for r in sel]
            else:
                outcome = [r[key_of[m]] for r in sel]
            c = collections.Counter(outcome)
            per_copy = collections.defaultdict(lambda: [0, 0])
            for r, o in zip(sel, outcome):
                if o in ("correct", "wrong", "conflict"):
                    u = per_copy[(r["family"], r["copy"])]
                    u[0] += o == "correct"
                    u[1] += 1
            ci = copy_level_ci([tuple(u) for _, u in sorted(per_copy.items())])
            n = len(sel)
            a = c["correct"] + c["wrong"] + c["conflict"]
            out.append([sample, species, catalog_scope, m, band, n, c["correct"], c["wrong"], c["conflict"],
                        c["abstain"], c["not_scored"], a / n if n else float("nan"), c["correct"] / a if a else float("nan"),
                        len({(r["family"], r["copy"]) for r in sel}), len(per_copy),
                        ci[0] if ci else None, ci[1] if ci else None, ci[2] if ci else "",
                        sum(r["family_psv"] for r in sel) if m == "family_certificate" else None,
                        sum(r["twin"] == "nm_twin" for r in sel) if m == "union_certificate" else None])
    return out


def wilson(k: int, n: int, z: float = 1.96):
    """Wilson score interval for k successes of n (None when n == 0)."""
    if n <= 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def fnum(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# ================================================================ builds shared by Figures 4 and 5
def collect(cfg: dict, *, force=False, recorded=False, budget: Budget | None = None) -> list[tuple[dict, dict, list[dict]]]:
    """[(runs, logs, per-read records)] for every sample of the build, in registry order. Genome scope: stops with
    Pending after the call's heavy-step budget (re-run the same command); the samples done so far stay cached."""
    out = []
    budget = budget or budget_from(cfg)
    if recorded:
        for sp in DEV_SPECIES:
            runs = recorded_runs(cfg, sp)
            logs = score_logs(runs, figlib.work_dir(cfg, "o2sim") / "recorded" / "score")
            out.append((runs, logs, per_read(runs, logs)))
        return out
    scope = o2_scope(cfg)
    pending = []
    for sid in sample_ids(cfg):
        try:
            runs = ensure_runs(cfg, sid, force=force, budget=budget)
        except Pending as e:
            pending.append(str(e))
            break
        log_dir = (figlib.work_dir(cfg, "o2sim") / runs["species"] if scope == "dev" else genome_work_dir(cfg, sid)) / "score"
        logs = score_logs(runs, log_dir)
        out.append((runs, logs, per_read(runs, logs)))
    if pending:
        raise SystemExit(f"[_o2] not finished — {pending[0]}. Re-run the same command (under flock) to continue.")
    return out


def run_notes(runs: dict, reads: list[dict], logs: dict) -> list[str]:
    key = runs["sample"]
    return ([f"{key}: {n}" for n in crosscheck(reads, logs)] + [f"{key}: {runs['source']}"]
            + [f"{key}: {n}" for n in runs.get("notes", [])])


# ================================================================ the alignment-score margin rule (experiment C)
# docs/PREREG_genome_wide_copy_assignment_2026-09-25.md, Amendment 2: the Eichler lab's rule MR(T) (assign a read to its
# best alignment iff no other alignment of the read, anywhere in the genome, scores within T AS units; a read with no
# other alignment is assigned) against Rustle's per-read answer (MAPQ > 0: the aligner's primary; MAPQ 0: the
# copy-assignment result, union test = reading 'u', source family* = reading 's'). The rule and the join are computed
# by `bench/score.py eichler --sim` from the same sim.bam and copy_assign runs as Figures 4 and 5; this module checks the
# join read by read against those figures' own records and tabulates it.
MR_THRESHOLDS = list(_score.EICHLER_THRESHOLDS)
MR_HEADLINE_T = 10
MR_READINGS = ["u", "s"]
MR_READING_LABEL = {"u": "union test", "s": "test scored in the source family*"}
MR_STRATA_HEADER = ["sample", "species", "catalog_scope", "identity_band", "threshold", "reading", "stratum",
                    "n_reads", "n_copies", "n_margin0", "margin_rule_correct", "rustle_correct"]
MR_ACC_HEADER = ["sample", "species", "catalog_scope", "identity_band", "method", "threshold", "reading", "n_reads",
                 "assigned", "correct", "fraction_assigned", "fraction_correct", "n_copies_assigned", "acc_lo",
                 "acc_hi", "acc_ci"]
# fraction-correct rows: the margin rule at each T, Rustle per reading (T-independent), and the two "Rustle only"
# strata per T and reading (their fraction_assigned = the stratum's share of the row's reads)
MR_METHODS = ["margin_rule", "rustle", "rustle_only_aligner", "rustle_only_test"]
_ALIGNER_OF = {"correct": "true_copy", "wrong": "other_copy", "outside": "outside_catalog", "unmapped": "unmapped"}


def margin_rule_per_read(runs: dict, log_dir: Path) -> tuple[Path, Path]:
    """`score.py eichler --sim` on one sample's runs (light: one pass over sim.bam; cached on its inputs)."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    key = runs.get("sample") or runs["species"]
    out, log = log_dir / f"margin_rule_{key}.tsv", log_dir / f"margin_rule_{key}.txt"
    src = (f"{runs['sim']}.bam", f"{runs['o2']}.assignments.tsv", f"{runs['u2']}.assignments.tsv", runs["catalog"],
           SCORE)
    if not (figlib.fresh(out, *src) and figlib.fresh(log, *src)):
        figlib.run([sys.executable, str(SCORE), "eichler", "--sim", str(runs["sim"]), "--catalog", str(runs["catalog"]),
                    "--default", str(runs["o2"]), "--union", str(runs["u2"]), "--per-read", str(out)], log=log)
    return out, log


def margin_rule_join(path, reads: list[dict]) -> list[dict]:
    """The per-read rows of `score.py eichler --sim`, each checked against this module's per_read() record of the same
    read (the read set, the aligner's placement, and the MAPQ-0 verdicts of both readings must be identical) and given
    its identity band and source-copy key. Raises on any difference."""
    rec = {r["read"]: r for r in reads}
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    names = {r["read_name"] for r in rows}
    if names != set(rec) or len(names) != len(rows):
        raise RuntimeError(f"{path}: {len(rows)} rows / {len(names)} reads vs {len(rec)} simulated reads in per_read()")
    for r in rows:
        p = rec[r["read_name"]]
        mq = "" if p["mapq"] is None else str(p["mapq"])
        if r["mapq"] != mq or _ALIGNER_OF.get(r["primary_verdict"]) != p["aligner"]:
            raise RuntimeError(f"{path}: {r['read_name']}: MAPQ / aligner placement differ from per_read() "
                               f"({r['mapq']}, {r['primary_verdict']}) vs ({mq}, {p['aligner']})")
        if mq == "0" and (_VERDICT.get(r["rs_verdict"], r["rs_verdict"]) != p["own_outcome"]
                          or _VERDICT.get(r["ru_verdict"], r["ru_verdict"]) != p["union_outcome"]):
            raise RuntimeError(f"{path}: {r['read_name']}: test verdicts differ from score.py reads --per-read")
        r["band"] = p["band"]
        r["_copy"] = (p["family"], p["copy"])
    return rows


def margin_rule_tables(rows: list[dict], sample: str, species: str, scope: str) -> tuple[list[list], list[list]]:
    """(strata rows, fraction-correct rows) per identity band ('all' + figlib.IDENTITY_BANDS), threshold and reading.
    Every count is over the band's simulated reads; fraction_assigned = assigned / n_reads of the band; intervals =
    copy_level_ci over the assigned reads grouped by source copy."""
    bands = ["all"] + figlib.IDENTITY_BANDS
    n_band = collections.Counter()
    strata = collections.defaultdict(lambda: [0, set(), 0, 0, 0])   # n, copies, margin0, mr correct, rustle correct
    units = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))   # key -> copy -> [correct, n]
    for r in rows:
        bs = ("all", r["band"])
        rust = {x: (_score.rustle_assigns(r, x), _score.rustle_correct(r, x)) for x in MR_READINGS}
        for b in bs:
            n_band[b] += 1
            for x in MR_READINGS:
                a, c = rust[x]
                if a:
                    u = units[(b, "rustle", "", x)][r["_copy"]]
                    u[0] += c
                    u[1] += 1
        for T in MR_THRESHOLDS:
            ma = _score.mr_assigns(r, T)
            mc = ma and r["mr_verdict"] == "correct"
            for b in bs:
                if ma:
                    u = units[(b, "margin_rule", T, "")][r["_copy"]]
                    u[0] += mc
                    u[1] += 1
            for x in MR_READINGS:
                s = _score.mr_stratum(r, T, x)
                a, c = rust[x]
                for b in bs:
                    st = strata[(b, T, x, s)]
                    st[0] += 1
                    st[1].add(r["_copy"])
                    st[2] += r["margin"] == "0"
                    st[3] += mc
                    st[4] += a and c
                    if s in ("rustle_aligner", "rustle_test"):
                        u = units[(b, "rustle_only_" + s.split("_", 1)[1], T, x)][r["_copy"]]
                        u[0] += c
                        u[1] += 1
    srows, arows = [], []
    for b in bands:
        if not n_band[b]:
            continue
        for T in MR_THRESHOLDS:
            for x in MR_READINGS:
                for s in _score.MR_STRATA:
                    n, cps, m0, mc, rc = strata.get((b, T, x, s), [0, set(), 0, 0, 0])
                    srows.append([sample, species, scope, b, T, x, s, n, len(cps), m0, mc, rc])
        keys = ([("margin_rule", T, "") for T in MR_THRESHOLDS] + [("rustle", "", x) for x in MR_READINGS]
                + [(m, T, x) for m in ("rustle_only_aligner", "rustle_only_test") for T in MR_THRESHOLDS
                   for x in MR_READINGS])
        for m, T, x in keys:
            per = units.get((b, m, T, x), {})
            k = sum(v[0] for v in per.values())
            a = sum(v[1] for v in per.values())
            ci = copy_level_ci([tuple(v) for _, v in sorted(per.items())])
            n = n_band[b]
            arows.append([sample, species, scope, b, m, T, x, n, a, k, a / n if n else float("nan"),
                          k / a if a else float("nan"), len(per), ci[0] if ci else None, ci[1] if ci else None,
                          ci[2] if ci else ""])
    return srows, arows


def margin_rule_notes(sample: str, rows: list[dict]) -> list[str]:
    """Counts the caption quotes that are not sums of table rows: the exceptions to containment at T = 10 by kind."""
    out = []
    for x in MR_READINGS:
        ex = [r for r in rows if _score.mr_stratum(r, MR_HEADLINE_T, x) in ("both_differ", "margin_only")]
        kinds = collections.Counter(
            (_score.mr_stratum(r, MR_HEADLINE_T, x), _score.rustle_source(r), r["mr_verdict"],
             "rustle_correct" if _score.rustle_correct(r, x) else "rustle_not_correct") for r in ex)
        out.append(f"{sample}: containment exceptions at T = {MR_HEADLINE_T}, reading {MR_READING_LABEL[x]}: "
                   f"{len(ex)} reads; (stratum, Rustle source, margin-rule verdict, Rustle) = "
                   + ("; ".join(f"{k} {v}" for k, v in kinds.most_common()) or "none"))
    return out


def runs_as_built(cfg: dict, sample: str) -> dict:
    """The runs dict of one sample WITHOUT building anything (for tables derived from runs that already exist): dev
    scope = ${work}/o2sim/<species>/, genome scope = ${work}/o2sim/<sample>/. Raises when a product is missing."""
    if o2_scope(cfg) == "dev":
        sp = sample if sample in DEV_SPECIES else next(k for k, v in DEV_SAMPLE.items()
                                                        if v == samples.resolve(cfg, sample))
        d = figlib.work_dir(cfg, "o2sim") / sp
        runs = {"species": sp, "sample": DEV_SAMPLE[sp], "catalog_scope": DEV_CATALOG[sp]}
    else:
        sid = samples.resolve(cfg, sample)
        d = genome_work_dir(cfg, sid)
        runs = {"species": species_of(cfg, sid), "sample": sid, "catalog_scope": "genome",
                "copy_table": copy_table(cfg)}
        if not Path(f"{d}/sim.done").exists():
            raise SystemExit(f"[_o2] {sid}: no finished simulation in {d}")
    runs.update({"sim": d / "sim", "catalog": d / "cat.copies.tsv", "o2": d / "o2", "u2": d / "u2", "notes": [],
                 "source": f"runs as built in {d} (not rebuilt)"})
    for p in (f"{runs['sim']}.bam", runs["catalog"], f"{runs['o2']}.assignments.tsv", f"{runs['u2']}.assignments.tsv"):
        if not Path(p).exists():
            raise SystemExit(f"[_o2] {runs['sample']}: {p} missing; build the runs with `make.py data fig4`")
    t = max(Path(f"{runs[k]}.assignments.tsv").stat().st_mtime for k in ("o2", "u2"))
    runs["notes"].append(f"copy_assign runs used as built (last written {time.strftime('%Y-%m-%d %H:%M', time.localtime(t))}); "
                         f"the copy_assign binary now at {cfg['bin']}/copy_assign is "
                         f"{'newer' if (Path(cfg['bin']) / 'copy_assign').stat().st_mtime > t else 'older'}")
    return runs


# ================================================================ check V3 (shards == one run), CLI
def _sorted_rows(path) -> tuple[str, list[str]]:
    with open(path) as fh:
        h = fh.readline()
        return h, sorted(l for l in fh if l.strip())


def compare_tables(a, b) -> list[str]:
    """Differences between two TSVs compared as row multisets (same header required); [] when identical."""
    ha, ra = _sorted_rows(a)
    hb, rb = _sorted_rows(b)
    if ha != hb:
        return [f"headers differ: {a} vs {b}"]
    if ra == rb:
        return []
    ca, cb = collections.Counter(ra), collections.Counter(rb)
    only_a, only_b = list((ca - cb).elements()), list((cb - ca).elements())
    cols = ha.rstrip("\n").split("\t")
    diff_cols = collections.Counter()
    ka = {tuple(l.split("\t")[:2]): l for l in only_a}
    for l in only_b:
        k = tuple(l.split("\t")[:2])
        if k in ka:
            for col, x, y in zip(cols, ka[k].rstrip("\n").split("\t"), l.rstrip("\n").split("\t")):
                if x != y:
                    diff_cols[col] += 1
    return [f"{len(only_a)} rows only in {a}, {len(only_b)} only in {b}; differing columns (rows keyed by the first "
            f"two columns): {dict(diff_cols) or 'none paired'}"]


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="figures/_o2.py utilities: shard planning and merging (check V3)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan", help="plan shards of a catalog over a BAM (prints the plan; writes OUT)")
    p.add_argument("--bam", required=True)
    p.add_argument("--catalog", required=True, help="copies.tsv (copies.fa next to it, same prefix)")
    p.add_argument("--regions", required=True, help="the one-run's regions file (chrom:1-LEN per line)")
    p.add_argument("--out", required=True)
    p.add_argument("--kind", choices=sorted(SHARD_COST), default="sim")
    p.add_argument("--shards", type=int, help="force this many balanced shards (check V3)")
    p.add_argument("--budget-s", type=float, default=SHARD_BUDGET_S)
    p = sub.add_parser("merge", help="merge shard outputs PREFIX... into OUT (assignments + union_certificate)")
    p.add_argument("--out", required=True)
    p.add_argument("--union", action="store_true")
    p.add_argument("prefixes", nargs="+")
    p = sub.add_parser("compare", help="compare two assignment tables as row multisets")
    p.add_argument("a")
    p.add_argument("b")
    a = ap.parse_args(argv)
    if a.cmd == "plan":
        contigs = []
        for line in open(a.regions):
            c, _, rng = line.strip().rpartition(":")
            contigs.append((c, int(rng.split("-")[1])))
        fa = a.catalog[: -len(".tsv")] + ".fa" if a.catalog.endswith(".tsv") else a.catalog + ".fa"
        plan = plan_shards(a.bam, a.catalog, fa, contigs, Path(a.out), kind=a.kind, budget_s=a.budget_s,
                           n_shards=a.shards)
        for s in plan["shards"]:
            print("\t".join(str(s[k]) for k in ("name", "n_families", "n_copies", "n_components", "records", "est_s",
                                                "regions")))
    elif a.cmd == "merge":
        merge_outputs([Path(x) for x in a.prefixes], Path(a.out), union=a.union)
    elif a.cmd == "compare":
        d = compare_tables(a.a, a.b)
        print("\n".join(d) if d else f"IDENTICAL as row multisets: {a.a} == {a.b}")
        sys.exit(1 if d else 0)


if __name__ == "__main__":
    main()
