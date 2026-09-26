"""_liftoff — the Liftoff self-lift locus baseline and Liftoff's matching criteria (Figure 8).

Pre-registration (read it first): docs/PREREG_liftoff_loci_2026-09-25.md. User decision 2026-09-25: every locus
comparison is made in the Liftoff framework. Per species, `liftoff -copies` lifts the genome's own RefSeq annotation
(gene + pseudogene records) onto the same genome = the annotation-guided locus baseline; Rustle's guided-mode loci are
compared with it like for like, and Rustle's de novo loci, default de novo families (the driver's `families` stage
copy table; prereg amendment 5; the legacy copy catalog is scored by the same code as a labelled secondary set) and
missing-copy flags are scored against it as a reference. Its (record, extra copy) pairs are also the Liftoff family
reference of Figs 6s and 7 (copy_pairs / pair_families below).

UNITS OF WORK (every heavy unit is cached, <= 10 min, and budgeted: a call that runs out of budget raises
assembly.Pending = exit 75, "run the same command again"):

  per species (`${work}/liftoff/<species>/`)
    prepare   records.tsv (every gene/pseudogene record, its exon union) + one GFF per shard (a contig's records, or a
              block of them cut where no record spans; LIFTOFF_BLOCK_RECORDS)           two GFF passes, 1-4 min
    index     minimap2 -d genome.fa.mmi with Liftoff's own options (human: the §6iv index is reused, header checked)
    shard     liftoff -copies on one shard's records against the WHOLE genome (-sc 0.95, -f pseudogene)
    merge     liftoff_loci.tsv: classes in_place / moved / partial / extra_copy, the cross-shard overlap rule (M1-M3)
  per sample (`${work}/liftoff/support/<sample>.tsv`)
    support   >= 2 reads whose primary alignment (-F 2308) has an aligned block on the locus's exon union
  per species, optional (`fig8_guided_finder=1`; `${work}/liftoff/<species>/finder/`)
    finder    Rustle's guided candidate-locus search with every record as a seed (bench/guided_pipeline.py finders)

MATCHING (prereg §3). In one genome a coordinate match has no mismatch and no insertion, so Liftoff's sequence_ID of a
locus G placed on locus R's exons equals its coverage; both criteria reduce to
    cov(G | R) = |X(G) ∩ X(R)| / |X(G)|   (X = exon union, same contig, strand ignored)
G is found by R iff cov(G | R) >= 0.5 (Liftoff's -a) for some single R.

CLI (light unless it says so):
    python3 figures/_liftoff.py plan                      what is done / left per species and sample
    python3 figures/_liftoff.py run --species S [--budget-s 540] [--max-units N]     HEAVY (flock)
    python3 figures/_liftoff.py support --sample ID [--budget-s 540]                   moderate (flock)
    python3 figures/_liftoff.py finder --species S [--budget-s 540]                    HEAVY (flock)
    python3 figures/_liftoff.py validate [--budget-s 540]                              V-L1 (HEAVY, flock)
"""
from __future__ import annotations

import bisect
import collections
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import figlib

# ---------------------------------------------------------------- pre-registered parameters (prereg §1-§3)
LIFTOFF_ENV = "/home/juanfra/miniforge3/envs/liftoff"
LIFT_A = 0.5            # -a (default)
LIFT_S = 0.5            # -s (default)
LIFT_SC = 0.95          # -sc (pre-registered; default 1.0)
LIFT_OVERLAP = 0.1      # -overlap (default)
LIFT_THREADS = 4        # -p
LIFT_TYPES = ["pseudogene"]   # -f: parent types added to Liftoff's default "gene"
MM2_OPTIONS = "-a --end-bonus 5 --eqx -N 50 -p 0.5"   # Liftoff's default -mm2_options (its index is built with them)
SC_ROWS = [0.95, 0.98, 0.99, 1.0]                     # sensitivity rows: filters of the -sc 0.95 run
MATCH_COV = 0.5         # cov(G | R) >= 0.5 (Liftoff's -a)
SHORT_BP = 200          # exon union < 200 bp: reported, never scored for Rustle
SUPPORT_READS = 2       # >= 2 primary reads with an aligned block on the exon union (the rule of Fig. 6d)
BOOT_N, BOOT_SEED = 2000, 20260925
LIFTOFF_BLOCK_RECORDS = 2500   # a contig with more records is split into blocks (cfg fig8_liftoff_block_records)
HARD_LIMIT_S = 590             # a Liftoff call is killed after this (the 10-minute rule; cfg fig8_liftoff_hard_s)
DEV_CONTIGS = {"human": {"chr16"}, "gorilla": {"NC_073244.2"}}   # development contigs (prereg §5 exposure)
PARENT_TYPES = ("gene", "pseudogene")
LEAF_TYPES = {"exon", "CDS", "cDNA_match", "match", "start_codon", "stop_codon"}

GP_DIR = figlib.REPO / "bench"


# ---------------------------------------------------------------- small helpers
def log(msg: str):
    print(f"[liftoff] {msg}", file=sys.stderr, flush=True)


def fp(path) -> str:
    p = Path(path)
    try:
        st = p.stat()
        return f"{p.resolve()}|{st.st_size}|{int(st.st_mtime_ns)}"
    except OSError:
        return f"{p}|absent"


def attrs(col9: str) -> dict:
    out = {}
    for kv in col9.rstrip("\n").split(";"):
        k, sep, v = kv.partition("=")
        if sep:
            out[k.strip()] = v
    return out


def merge_iv(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def iv_str(iv) -> str:
    return ",".join(f"{s}-{e}" for s, e in iv)


def iv_parse(s: str):
    if not s or s in ("-", "NA"):
        return []
    out = []
    for x in s.split(","):
        a, _, b = x.partition("-")
        out.append((int(a), int(b)))
    return out


def iv_len(iv) -> int:
    return sum(e - s for s, e in iv)


def iv_inter(a, b) -> int:
    """|A ∩ B| for two merged, sorted interval lists."""
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def ov(a0, a1, b0, b1) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


# ---------------------------------------------------------------- species and targets
def species_list(cfg: dict) -> list[str]:
    import samples
    out = []
    for sid, row in samples.registry(cfg).items():
        if row["species"] not in out:
            out.append(row["species"])
    return out


def species_inputs(cfg: dict, species: str) -> dict:
    """The species' genome FASTA and RefSeq GFF from the sample registry; every sample of a species must agree."""
    import samples
    rows = [r for r in samples.registry(cfg).values() if r["species"] == species]
    if not rows:
        raise KeyError(f"no sample of species {species!r} in the registry")
    fa = {r["fasta"] for r in rows}
    gff = {r["annotation_gff"] for r in rows}
    if len(fa) != 1 or len(gff) != 1:
        raise ValueError(f"{species}: samples disagree on the genome or annotation: {fa} {gff}")
    return {"fasta": Path(fa.pop()), "gff": Path(gff.pop()), "samples": [r["id"] for r in rows]}


def liftoff_root(cfg: dict) -> Path:
    d = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "liftoff"
    d.mkdir(parents=True, exist_ok=True)
    return d


class Target:
    """One self-lift: a genome FASTA, its annotation, a work dir (a species, or the V-L1 mini-genome)."""

    def __init__(self, name: str, fasta: Path, gff: Path, wdir: Path, *, block_records: int, hard_s: int,
                 reuse_index: Path | None = None, contigs: list[str] | None = None, single: bool = False):
        self.name, self.fasta, self.gff, self.wdir = name, Path(fasta), Path(gff), Path(wdir)
        self.block_records, self.hard_s = block_records, hard_s
        self.reuse_index = reuse_index
        self.contigs = contigs          # restrict records to these contigs (V-L1)
        self.single = single            # one shard with every record (V-L1 arm A)
        self.wdir.mkdir(parents=True, exist_ok=True)

    # paths
    @property
    def genome(self) -> Path:
        return self.wdir / "genome.fa"

    @property
    def mmi(self) -> Path:
        return self.wdir / "genome.fa.mmi"

    @property
    def records(self) -> Path:
        return self.wdir / "records.tsv"

    @property
    def plan_path(self) -> Path:
        return self.wdir / "plan.tsv"

    @property
    def loci(self) -> Path:
        return self.wdir / "liftoff_loci.tsv"

    def shard_dir(self, label: str) -> Path:
        return self.wdir / "shards" / label


# Liftoff indexes already on this machine (built by Liftoff itself with its default options: k15 w10), reused when
# the header matches and they are newer than the FASTA (ensure_index checks); cfg `fig8_liftoff_index_<species>`
# overrides. human: §6iv (next to chm13v2.0.fa); gorilla / chimpanzee: the 2026-09-17 NPIP lift (ggo_npip/ref/).
KNOWN_LIFTOFF_INDEX = {
    "gorilla": "/mnt/linuxdisk/home/juanfraitu/ggo_npip/ref/GGO.fa.mmi",
    "chimpanzee": "/mnt/linuxdisk/home/juanfraitu/ggo_npip/ref/PTR.fa.mmi",
}


def species_target(cfg: dict, species: str) -> Target:
    inp = species_inputs(cfg, species)
    reuse = None
    for cand in (cfg.get(f"fig8_liftoff_index_{species}"), KNOWN_LIFTOFF_INDEX.get(species),
                 str(inp["fasta"]) + ".mmi"):
        if cand and Path(cand).exists():
            reuse = Path(cand)
            break
    return Target(species, inp["fasta"], inp["gff"], liftoff_root(cfg) / species,
                  block_records=int(cfg.get("fig8_liftoff_block_records", LIFTOFF_BLOCK_RECORDS)),
                  hard_s=int(cfg.get("fig8_liftoff_hard_s", HARD_LIMIT_S)), reuse_index=reuse)


def params_key() -> dict:
    return {"liftoff": "v1.6.3", "env": LIFTOFF_ENV, "a": LIFT_A, "s": LIFT_S, "sc": LIFT_SC, "overlap": LIFT_OVERLAP,
            "p": LIFT_THREADS, "types": LIFT_TYPES, "mm2": MM2_OPTIONS, "copies": True}


# ---------------------------------------------------------------- budget
class Budget:
    """Seconds left in this call (0/None = unlimited). `fits(pred)`: the first unit of a call always runs."""

    def __init__(self, seconds: float | None, max_units: int | None = None):
        self.limit = float(seconds or 0)
        self.t0 = time.time()
        self.units = 0
        self.max_units = max_units

    def remaining(self) -> float:
        return math.inf if self.limit <= 0 else self.limit - (time.time() - self.t0)

    def fits(self, pred_s: float) -> bool:
        if self.max_units is not None and self.units >= self.max_units:
            return False
        if self.units == 0:
            return self.remaining() > 0
        return pred_s <= self.remaining()

    def used(self):
        self.units += 1


def pending(msg: str):
    import assembly
    raise assembly.Pending(msg)


# ---------------------------------------------------------------- prepare: records + shard GFFs
REC_HEAD = ["idx", "id", "name", "rtype", "biotype", "contig", "start0", "end", "strand", "exons", "exonic_len", "shard"]


def _plan_blocks(recs: list[dict], cap: int) -> list[tuple[str, list[int]]]:
    """Shards per contig (natural order of first appearance): the contig's records, or blocks of <= cap records cut
    at positions no record spans (a block may exceed cap only when no record-free cut exists)."""
    by = collections.OrderedDict()
    for r in recs:
        by.setdefault(r["contig"], []).append(r)
    out = []
    for contig, rs in by.items():
        if len(rs) <= cap:
            out.append((contig, [r["idx"] for r in rs]))
            continue
        rs = sorted(rs, key=lambda r: (r["start0"], r["end"], r["idx"]))
        target = math.ceil(len(rs) / math.ceil(len(rs) / cap))   # balanced blocks of <= cap records
        blocks, cur, max_end = [], [], -1
        for r in rs:
            if cur and len(cur) >= target and r["start0"] >= max_end:
                blocks.append(cur)
                cur, max_end = [], -1
            cur.append(r)
            max_end = max(max_end, r["end"])
        if cur:
            blocks.append(cur)
        for i, b in enumerate(blocks, 1):
            out.append((f"{contig}.b{i}", sorted(x["idx"] for x in b)))
    return out


def _apply_splits(shards: list, recs: list[dict], splits: list[str]) -> list:
    """A shard that did not finish inside the hard limit is split in two at the record-free cut nearest its middle
    (labels L.s1 / L.s2; again if one of those does not finish). A shard of one record, or one with no record-free cut,
    cannot be split: prereg §6 stop rule (cluster)."""
    by_idx = {r["idx"]: r for r in recs}
    todo = list(shards)
    out = []
    while todo:
        lab, idxs = todo.pop(0)
        if lab not in splits:
            out.append((lab, idxs))
            continue
        rs = sorted((by_idx[i] for i in idxs), key=lambda r: (r["start0"], r["end"], r["idx"]))
        cuts, max_end = [], -1
        for k, r in enumerate(rs):
            if k and r["start0"] >= max_end:
                cuts.append(k)
            max_end = max(max_end, r["end"])
        if not cuts:
            raise RuntimeError(f"shard {lab} ({len(rs)} records) did not finish inside the hard limit and cannot be "
                               f"split (no record-free cut): run it on the cluster (prereg §6)")
        k = min(cuts, key=lambda c: abs(c - len(rs) / 2))
        todo[0:0] = [(f"{lab}.s1", sorted(r["idx"] for r in rs[:k])), (f"{lab}.s2", sorted(r["idx"] for r in rs[k:]))]
    return out


def prepare(t: Target, *, force=False) -> tuple[dict, bool]:
    """records.tsv (one row per gene/pseudogene record, with its exon union) and the GFF of every unfinished shard.
    Two streaming passes over the annotation (parents precede children in RefSeq GFF3, but the parent map is built in
    pass 1 so the order does not matter)."""
    splits = json.loads((t.wdir / "splits.json").read_text()) if (t.wdir / "splits.json").exists() else []
    key = {"gff": fp(t.gff), "types": list(PARENT_TYPES), "block_records": t.block_records, "contigs": t.contigs,
           "single": t.single, "v": 2, **({"splits": splits} if splits else {})}
    kpath = t.wdir / "prepare.key"
    recs = None
    if not force and kpath.exists() and json.loads(kpath.read_text()) == key and t.records.exists() \
            and t.plan_path.exists():
        plan = read_plan(t)
        need = [lab for lab in plan if not shard_done(t, lab) and not (t.shard_dir(lab) / "in.gff").exists()]
        if not need:
            return plan, False
    else:
        need = None
    t0 = time.time()
    # pass 1: records and the parent map of every non-leaf feature
    parent_of: dict[str, str] = {}
    recs, rec_of = [], {}
    contigs = set(t.contigs) if t.contigs else None
    with open(t.gff) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] in LEAF_TYPES:
                continue
            a = attrs(f[8])
            fid = a.get("ID")
            if not fid:
                continue
            par = a.get("Parent")
            if par:
                parent_of[fid] = par.split(",")[0]
            elif f[2] in PARENT_TYPES and (contigs is None or f[0] in contigs):
                if fid in rec_of:
                    raise ValueError(f"{t.gff}: duplicate record ID {fid}")
                rec_of[fid] = len(recs)
                recs.append({"idx": len(recs), "id": fid, "name": a.get("Name", fid), "rtype": f[2],
                             "biotype": a.get("gene_biotype", "-"), "contig": f[0], "start0": int(f[3]) - 1,
                             "end": int(f[4]), "strand": f[6], "ex": [], "cds": []})
    top_memo: dict[str, str | None] = {}

    def top(fid: str):
        if fid in top_memo:
            return top_memo[fid]
        chain, cur, res = [], fid, None
        for _ in range(64):
            if cur in rec_of:
                res = cur
                break
            chain.append(cur)
            nxt = parent_of.get(cur)
            if nxt is None:
                break
            cur = nxt
        for c in chain:
            top_memo[c] = res
        top_memo[fid] = res
        return res

    if t.single:
        shards = [("all", [r["idx"] for r in recs])]
    else:
        shards = _plan_blocks(recs, t.block_records)
    shards = _apply_splits(shards, recs, splits)
    shard_of = {}
    for lab, idxs in shards:
        for i in idxs:
            shard_of[i] = lab
    plan = collections.OrderedDict((lab, len(idxs)) for lab, idxs in shards)
    write_labels = set(plan) if need is None else set(need)
    write_labels = {lab for lab in write_labels if not shard_done(t, lab)}
    handles = {}
    try:
        for lab in write_labels:
            d = t.shard_dir(lab)
            d.mkdir(parents=True, exist_ok=True)
            handles[lab] = open(d / "in.gff.tmp", "w")
            handles[lab].write("##gff-version 3\n")
        # pass 2: route every line to its record's shard; collect exon/CDS unions
        with open(t.gff) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.split("\t", 8)
                if len(f) < 9 or f[2] in ("cDNA_match", "match"):
                    continue
                a = attrs(f[8])
                fid = a.get("ID")
                if fid in rec_of and "Parent" not in a:
                    r = fid
                else:
                    par = a.get("Parent")
                    if not par:
                        continue
                    r = top(par.split(",")[0])
                if r is None:
                    continue
                rec = recs[rec_of[r]]
                if f[0] != rec["contig"]:
                    continue
                if f[2] == "exon":
                    rec["ex"].append((int(f[3]) - 1, int(f[4])))
                elif f[2] == "CDS":
                    rec["cds"].append((int(f[3]) - 1, int(f[4])))
                lab = shard_of[rec["idx"]]
                if lab in handles:
                    handles[lab].write(line)
    finally:
        for h in handles.values():
            h.close()
    for lab in write_labels:
        d = t.shard_dir(lab)
        (d / "in.gff.tmp").replace(d / "in.gff")
    tmp = t.records.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fh:
        fh.write("\t".join(REC_HEAD) + "\n")
        for r in recs:
            iv = merge_iv(r["ex"] or r["cds"] or [(r["start0"], r["end"])])
            fh.write("\t".join(map(str, [r["idx"], r["id"], r["name"], r["rtype"], r["biotype"], r["contig"],
                                         r["start0"], r["end"], r["strand"], iv_str(iv), iv_len(iv),
                                         shard_of[r["idx"]]])) + "\n")
    tmp.replace(t.records)
    with open(t.plan_path, "w") as fh:
        fh.write("shard\tn_records\n")
        for lab, n in plan.items():
            fh.write(f"{lab}\t{n}\n")
    kpath.write_text(json.dumps(key))
    log(f"{t.name}: prepared {len(recs):,} records in {len(plan)} shards ({time.time() - t0:.0f} s; "
        f"{len(write_labels)} shard GFFs written)")
    return plan, True


def read_plan(t: Target) -> collections.OrderedDict:
    out = collections.OrderedDict()
    with open(t.plan_path) as fh:
        next(fh)
        for line in fh:
            lab, n = line.rstrip("\n").split("\t")
            out[lab] = int(n)
    return out


def read_records(t: Target) -> list[dict]:
    with open(t.records) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k in ("idx", "start0", "end", "exonic_len"):
            r[k] = int(r[k])
    return rows


# ---------------------------------------------------------------- index
def index_state(t: Target) -> str:
    return "done" if t.mmi.exists() and (t.wdir / "index.ok").exists() else "todo"


def _mmi_header(path) -> tuple:
    import struct
    with open(path, "rb") as fh:
        b = fh.read(24)
    if b[:4] != b"MMI\x02":
        raise ValueError(f"{path}: not a minimap2 index")
    return struct.unpack("<5I", b[4:24])   # w, k, b, n_seq, flag


def ensure_genome_link(t: Target):
    for src, dst in ((t.fasta, t.genome), (Path(str(t.fasta) + ".fai"), Path(str(t.genome) + ".fai"))):
        if src.exists() and not dst.exists():
            dst.symlink_to(src)
    if not t.genome.exists():
        raise FileNotFoundError(t.fasta)


def ensure_index(t: Target):
    """minimap2 -d genome.fa.mmi with Liftoff's options (Liftoff's own build_minimap2_index command), so every shard
    reuses one index. Human reuses the §6iv index when its header matches (k15 w10, flag 0) and names every contig."""
    ensure_genome_link(t)
    if index_state(t) == "done":
        return
    n_seq = sum(1 for _ in open(str(t.fasta) + ".fai"))
    if t.reuse_index is not None and not t.mmi.exists():
        w, k, b, n, flag = _mmi_header(t.reuse_index)
        if (w, k, flag) == (10, 15, 0) and n == n_seq and t.reuse_index.stat().st_mtime > t.fasta.resolve().stat().st_mtime:
            t.mmi.symlink_to(t.reuse_index)
            (t.wdir / "index.ok").write_text(f"reused {t.reuse_index} (w{w} k{k} n_seq {n})\n")
            log(f"{t.name}: reusing the index {t.reuse_index}")
            return
        log(f"{t.name}: {t.reuse_index} does not match (w{w} k{k} n{n} flag{flag}); building a new one")
    mm2 = f"{LIFTOFF_ENV}/bin/minimap2"
    cmd = [mm2, "-d", str(t.mmi), str(t.genome)] + MM2_OPTIONS.split(" ") + ["-t", str(LIFT_THREADS)]
    tmp_log = t.wdir / "index.log"
    log(f"{t.name}: building the Liftoff index ({' '.join(cmd)})")
    r = _timed(cmd, t.wdir / "index.time", tmp_log, cwd=t.wdir, hard_s=t.hard_s)
    if r != 0 or not t.mmi.exists():
        t.mmi.unlink(missing_ok=True)
        raise RuntimeError(f"{t.name}: index build failed ({r}); see {tmp_log}")
    w, k, b, n, flag = _mmi_header(t.mmi)
    if n != n_seq:
        t.mmi.unlink(missing_ok=True)
        raise RuntimeError(f"{t.name}: index has {n} sequences, genome {n_seq} (a multi-part index is not supported)")
    (t.wdir / "index.ok").write_text(f"built w{w} k{k} n_seq {n}\n")


TIME_FMT = "wall_s\t%e\npeak_rss_kb\t%M\nuser_s\t%U\nsys_s\t%S\nexit\t%x"

# Liftoff aligns the same gene sequences twice with the same options against the same index: once to place the
# records, once more in the -copies step (run_liftoff.map_extra_copies re-extracts every record's sequence into the
# same file). minimap2's output is a function of (index, query, options), so the second call can return the first
# call's SAM. This wrapper (passed as Liftoff's -m) does exactly that, keyed on the md5 of the query FASTA, every
# argument except -o, and the index's path/size/mtime; it keeps the cache inside the call's -dir, deleted with it.
MM2_CACHE_WRAPPER = r"""#!/bin/bash
# written by figures/_liftoff.py (MM2_CACHE_WRAPPER): minimap2 for Liftoff with a per-call result cache
REAL=__REAL__
if [ "$1" != "-o" ] || [ $# -lt 4 ]; then exec "$REAL" "$@"; fi
OUT=$2; IDX=$3; Q=$4
[ -f "$Q" ] && [ -f "$IDX" ] || exec "$REAL" "$@"
C="$(dirname "$OUT")/mm2cache"; mkdir -p "$C"
K=$( { md5sum < "$Q"; shift 2; printf '%s\n' "$@"; stat -L -c '%s %Y' "$IDX"; "$REAL" --version; } | md5sum | cut -c1-32)
if [ -s "$C/$K.sam" ] && [ -f "$C/$K.ok" ]; then
  cp "$C/$K.sam" "$OUT" && echo "[mm2cache] reused $C/$K.sam for $OUT" >&2 && exit 0
fi
"$REAL" -o "$OUT" "${@:3}"
rc=$?
if [ $rc -eq 0 ]; then cp "$OUT" "$C/$K.sam.tmp" && mv "$C/$K.sam.tmp" "$C/$K.sam" && : > "$C/$K.ok"; fi
exit $rc
"""


def mm2_wrapper(t_root: Path) -> Path:
    p = t_root / "mm2_cache.sh"
    body = MM2_CACHE_WRAPPER.replace("__REAL__", f"{LIFTOFF_ENV}/bin/minimap2")
    if not p.exists() or p.read_text() != body:
        p.write_text(body)
        p.chmod(0o755)
    return p


def _timed(cmd, time_file: Path, log_file: Path, *, cwd: Path, hard_s: int, env: dict | None = None) -> int:
    """Run cmd under /usr/bin/time in its own process group; kill the group after hard_s seconds (returns 124)."""
    e = dict(os.environ)
    e["PATH"] = f"{LIFTOFF_ENV}/bin:" + e.get("PATH", "")
    e.setdefault("TMPDIR", "/mnt/linuxdisk/tmp")
    e.update(env or {})
    full = ["/usr/bin/time", "-f", TIME_FMT, "-o", str(time_file)] + [str(c) for c in cmd]
    with open(log_file, "w") as lf:
        p = subprocess.Popen(full, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, env=e, start_new_session=True)
        try:
            return p.wait(timeout=hard_s)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGTERM)
            try:
                p.wait(timeout=15)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            return 124


def read_time(path: Path) -> dict:
    out = {}
    if path.exists():
        for line in open(path):
            k, _, v = line.strip().partition("\t")
            if v:
                try:
                    out[k] = float(v)
                except ValueError:
                    pass
    return out


# ---------------------------------------------------------------- shards
def shard_key(t: Target, label: str) -> dict:
    return {"params": params_key(), "genome": fp(t.fasta), "gff": fp(t.gff), "label": label,
            "block_records": t.block_records, "single": t.single, "contigs": t.contigs}


def shard_done(t: Target, label: str) -> bool:
    d = t.shard_dir(label)
    k = d / "done.key"
    return k.exists() and (d / "loci.tsv").exists() and json.loads(k.read_text()) == shard_key(t, label)


def predict_shard_s(t: Target, n_records: int, plan: dict | None = None) -> float:
    """Median measured seconds per record of this target's finished shards (x1.25); 0.2 s/record + 120 s before any
    shard has finished."""
    rates = []
    plan = plan or read_plan(t)
    for lab, n in plan.items():
        tm = read_time(t.shard_dir(lab) / "time.txt")
        if tm.get("wall_s") and n and shard_done(t, lab):
            rates.append(tm["wall_s"] / n)
    if not rates:
        return 0.2 * n_records + 120
    rates.sort()
    return rates[len(rates) // 2] * n_records * 1.25


def run_shard(t: Target, label: str):
    d = t.shard_dir(label)
    gff = d / "in.gff"
    if not gff.exists():
        raise RuntimeError(f"{t.name} {label}: {gff} missing (run prepare)")
    (d / "types.txt").write_text("\n".join(LIFT_TYPES) + "\n")
    inter = d / "inter"
    if inter.exists():
        shutil.rmtree(inter)
    for p in (d / "lifted.gff3", d / "unmapped.txt", Path(str(gff) + "_db")):
        p.unlink(missing_ok=True)
    cmd = [f"{LIFTOFF_ENV}/bin/liftoff", str(t.genome), str(t.genome), "-g", "in.gff", "-o", "lifted.gff3",
           "-u", "unmapped.txt", "-dir", "inter", "-copies", "-sc", str(LIFT_SC), "-a", str(LIFT_A), "-s",
           str(LIFT_S), "-overlap", str(LIFT_OVERLAP), "-p", str(LIFT_THREADS), "-f", "types.txt",
           "-m", str(mm2_wrapper(t.wdir))]
    t0 = time.time()
    log(f"{t.name} {label}: liftoff ({read_plan(t)[label]:,} records; hard limit {t.hard_s} s)")
    r = _timed(cmd, d / "time.txt", d / "liftoff.log", cwd=d, hard_s=t.hard_s)
    tm = read_time(d / "time.txt")
    if r == 124:
        (d / "toolong").write_text(f"killed after {t.hard_s} s at {read_plan(t)[label]} records\n")
        sp = t.wdir / "splits.json"
        splits = json.loads(sp.read_text()) if sp.exists() else []
        if label not in splits:
            sp.write_text(json.dumps(splits + [label]))
        shutil.rmtree(d / "inter", ignore_errors=True)
        pending(f"{t.name} {label}: Liftoff did not finish in {t.hard_s} s ({read_plan(t)[label]:,} records); the "
                f"shard is split in two at the next call")
    if r != 0 or not (d / "lifted.gff3").exists():
        raise RuntimeError(f"{t.name} {label}: liftoff failed ({r}); see {d / 'liftoff.log'}")
    n = parse_lifted(d / "lifted.gff3", d / "loci.tsv")
    with open(d / "lifted.gff3", "rb") as fi, gzip.open(d / "lifted.gff3.gz", "wb", compresslevel=4) as fo:
        shutil.copyfileobj(fi, fo)
    (d / "lifted.gff3").unlink()
    shutil.rmtree(inter, ignore_errors=True)
    Path(str(gff) + "_db").unlink(missing_ok=True)
    gff.unlink(missing_ok=True)
    (d / "done.key").write_text(json.dumps(shard_key(t, label)))
    log(f"{t.name} {label}: {n:,} placed features in {time.time() - t0:.0f} s "
        f"(peak {tm.get('peak_rss_kb', 0) / 1e6:.1f} GB)")


LOCI_HEAD = ["feature_id", "copy_num_id", "tag", "rtype", "contig", "start0", "end", "strand", "coverage",
             "sequence_id", "partial", "low_identity", "exons", "exonic_len"]


def parse_lifted(gff3: Path, out: Path) -> int:
    """Parent-level rows of one Liftoff output: its ID, copy tag, placed span, coverage / sequence_ID, flags, and the
    exon union of its children (CDS when it has no exon; else its span)."""
    parents, parent_of = {}, {}
    ex, cds = collections.defaultdict(list), collections.defaultdict(list)
    lines = []
    with open(gff3) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            fid = a.get("ID")
            if "copy_num_ID" in a and f[2] in PARENT_TYPES and "Parent" not in a:
                parents[fid] = (f, a)
                continue
            if fid and a.get("Parent"):
                parent_of[fid] = a["Parent"].split(",")[0]
            if f[2] in ("exon", "CDS"):
                lines.append((f[2], a.get("Parent", "").split(",")[0], int(f[3]) - 1, int(f[4])))

    def top(x):
        for _ in range(64):
            if x in parents:
                return x
            x = parent_of.get(x)
            if x is None:
                return None
        return None

    for typ, par, s, e in lines:
        p = top(par)
        if p is not None:
            (ex if typ == "exon" else cds)[p].append((s, e))
    tmp = out.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write("\t".join(LOCI_HEAD) + "\n")
        for fid, (f, a) in parents.items():
            iv = merge_iv(ex.get(fid) or cds.get(fid) or [(int(f[3]) - 1, int(f[4]))])
            fo.write("\t".join(map(str, [fid, a.get("copy_num_ID", "-"), a.get("extra_copy_number", "0"), f[2], f[0],
                                         int(f[3]) - 1, int(f[4]), f[6], a.get("coverage", "NA"),
                                         a.get("sequence_ID", "NA"), int(a.get("partial_mapping") == "True"),
                                         int(a.get("low_identity") == "True"), iv_str(iv), iv_len(iv)])) + "\n")
    tmp.replace(out)
    return len(parents)


# ---------------------------------------------------------------- merge (prereg §1 M1-M3, §2 classes)
MERGED_HEAD = ["species", "shard", "source_id", "name", "rtype", "biotype", "stratum", "short", "src_contig",
               "src_start0", "src_end", "src_strand", "cls", "tag", "contig", "start0", "end", "strand", "coverage",
               "sequence_id", "exons", "exonic_len", "note"]


def stratum_of(rtype: str, biotype: str) -> str:
    if rtype == "pseudogene" or "pseudogene" in (biotype or ""):
        return "pseudogene"
    if biotype == "protein_coding":
        return "protein_coding"
    if biotype == "lncRNA":
        return "lncRNA"
    return "other"


def merge(t: Target) -> Path:
    plan = read_plan(t)
    key = {lab: (t.shard_dir(lab) / "done.key").read_text() for lab in plan}
    kpath = t.wdir / "merge.key"
    if t.loci.exists() and kpath.exists() and json.loads(kpath.read_text()) == key:
        return t.loci
    recs = read_records(t)
    by_id = {r["id"]: r for r in recs}
    rows, unmapped, anomalies = [], set(), []
    for lab in plan:
        d = t.shard_dir(lab)
        if (d / "unmapped.txt").exists():
            unmapped.update(x.strip() for x in open(d / "unmapped.txt") if x.strip())
        with open(d / "loci.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                fid = row["feature_id"]
                tag = int(row["tag"])
                if fid in by_id and tag == 0:
                    src, is_copy = fid, False
                else:
                    src = fid.rsplit("_", 1)[0] if fid not in by_id else fid
                    if src not in by_id:
                        anomalies.append(f"{lab}: {fid} names no record")
                        continue
                    is_copy = True
                    if fid in by_id:
                        anomalies.append(f"{lab}: {fid} tag {tag} carries its record's own ID")
                r = by_id[src]
                row.update(shard=lab, source=src, is_copy=is_copy)
                rows.append(row)
    placed = {}   # source -> annotated placement row
    for row in rows:
        if not row["is_copy"]:
            if row["source"] in placed:
                anomalies.append(f"{row['source']}: placed twice")
            placed[row["source"]] = row
    # classes of annotated placements
    for row in rows:
        r = by_id[row["source"]]
        s0, e0 = int(row["start0"]), int(row["end"])
        if row["is_copy"]:
            row["cls"] = "extra_copy"
        elif row["partial"] == "1" or row["low_identity"] == "1":
            row["cls"] = "partial"
        elif row["contig"] == r["contig"] and ov(s0, e0, r["start0"], r["end"]) >= 0.5 * max(e0 - s0, r["end"] - r["start0"]):
            row["cls"] = "in_place"
        else:
            row["cls"] = "moved"
        row["note"] = ""
    # M1: an extra copy overlapping (same contig + strand, >= 1 bp) an annotated placement from ANOTHER shard whose
    # reference record does not overlap the copy's source record in the reference -> dropped
    ann_by = collections.defaultdict(list)
    for row in rows:
        if not row["is_copy"]:
            ann_by[(row["contig"], row["strand"])].append((int(row["start0"]), int(row["end"]), row))
    for v in ann_by.values():
        v.sort(key=lambda x: x[0])
    ann_starts = {k: [x[0] for x in v] for k, v in ann_by.items()}
    ann_maxlen = {k: max((x[1] - x[0] for x in v), default=0) for k, v in ann_by.items()}

    def ref_overlap(a, b) -> bool:
        ra, rb = by_id[a], by_id[b]
        return ra["contig"] == rb["contig"] and ra["strand"] == rb["strand"] and a != b and \
            ov(ra["start0"], ra["end"], rb["start0"], rb["end"]) > 0

    n_m1 = 0
    for row in rows:
        if row["cls"] != "extra_copy":
            continue
        k = (row["contig"], row["strand"])
        s0, e0 = int(row["start0"]), int(row["end"])
        v = ann_by.get(k, [])
        lo = bisect.bisect_left(ann_starts.get(k, []), s0 - ann_maxlen.get(k, 0))
        for i in range(lo, len(v)):
            a0, a1, arow = v[i]
            if a0 >= e0:
                break
            if a1 > s0 and arow["shard"] != row["shard"] and not ref_overlap(row["source"], arow["source"]):
                row["cls"], row["note"] = "dropped_M1", f"overlaps {arow['source']} ({arow['shard']})"
                n_m1 += 1
                break
    # M2: overlapping extra copies from different shards (same contig + strand): keep higher sequence_ID, then longer,
    # then the source record earlier in GFF order
    def fnum(x):
        try:
            return float(x)
        except ValueError:
            return -1.0
    cops = [row for row in rows if row["cls"] == "extra_copy"]
    cops.sort(key=lambda row: (-fnum(row["sequence_id"]), -(int(row["end"]) - int(row["start0"])),
                               by_id[row["source"]]["idx"]))
    kept = collections.defaultdict(list)
    n_m2 = 0
    for row in cops:
        k = (row["contig"], row["strand"])
        s0, e0 = int(row["start0"]), int(row["end"])
        hit = next((o for o in kept[k] if o["shard"] != row["shard"] and ov(s0, e0, int(o["start0"]), int(o["end"])) > 0), None)
        if hit is not None:
            row["cls"], row["note"] = "dropped_M2", f"overlaps copy of {hit['source']} ({hit['shard']})"
            n_m2 += 1
        else:
            kept[k].append(row)
    # M3: annotated placements of different shards overlapping on the same strand (counted, kept)
    n_m3 = 0
    for k, v in ann_by.items():
        act = []
        for a0, a1, arow in v:
            act = [x for x in act if x[1] > a0]
            for b0, b1, brow in act:
                if brow["shard"] != arow["shard"] and not ref_overlap(arow["source"], brow["source"]):
                    n_m3 += 1
                    arow["note"] = (arow["note"] + "; " if arow["note"] else "") + f"M3 overlaps {brow['source']}"
            act.append((a0, a1, arow))
    tmp = t.loci.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write("\t".join(MERGED_HEAD) + "\n")
        for row in sorted(rows, key=lambda x: (x["contig"], int(x["start0"]), x["source"])):
            r = by_id[row["source"]]
            st = stratum_of(r["rtype"], r["biotype"])
            fo.write("\t".join(map(str, [t.name, row["shard"], r["id"], r["name"], r["rtype"], r["biotype"], st,
                                         int(int(row["exonic_len"]) < SHORT_BP), r["contig"], r["start0"], r["end"],
                                         r["strand"], row["cls"], row["tag"], row["contig"], row["start0"],
                                         row["end"], row["strand"], row["coverage"], row["sequence_id"],
                                         row["exons"], row["exonic_len"], row["note"] or "-"])) + "\n")
        placed_ids = {row["source"] for row in rows if not row["is_copy"]}
        for r in recs:
            if r["id"] not in placed_ids:
                st = stratum_of(r["rtype"], r["biotype"])
                fo.write("\t".join(map(str, [t.name, r["shard"], r["id"], r["name"], r["rtype"], r["biotype"], st,
                                             int(r["exonic_len"] < SHORT_BP), r["contig"], r["start0"], r["end"],
                                             r["strand"], "unmapped", "-", "-", "-", "-", "-", "NA", "NA", "-", 0,
                                             "listed unmapped" if r["id"] in unmapped else "absent from the output"]))
                         + "\n")
    tmp.replace(t.loci)
    rep = {"records": len(recs), "placed": len(placed), "extra_copies_raw": sum(1 for r in rows if r["is_copy"]),
           "M1_dropped": n_m1, "M2_dropped": n_m2, "M3_overlaps": n_m3, "anomalies": anomalies[:50],
           "n_anomalies": len(anomalies), "shards": len(plan)}
    (t.wdir / "merge_report.json").write_text(json.dumps(rep, indent=1))
    kpath.write_text(json.dumps(key))
    log(f"{t.name}: merged {len(plan)} shards: {rep}")
    return t.loci


def read_loci(path: Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["iv"] = iv_parse(r["exons"])
        r["exonic_len"] = int(r["exonic_len"])
        r["short"] = r["short"] == "1"
        r["seqid_f"] = float(r["sequence_id"]) if r["sequence_id"] not in ("NA", "") else float("nan")
    return rows


# ---------------------------------------------------------------- driving one target
def target_state(t: Target) -> dict:
    st = {"prepare": "done" if t.plan_path.exists() and (t.wdir / "prepare.key").exists() else "todo",
          "index": index_state(t), "shards_done": 0, "shards": 0, "records_left": 0, "merge": "todo"}
    if st["prepare"] == "done":
        plan = read_plan(t)
        st["shards"] = len(plan)
        done = [lab for lab in plan if shard_done(t, lab)]
        st["shards_done"] = len(done)
        st["records_left"] = sum(n for lab, n in plan.items() if lab not in done)
        st["merge"] = "done" if t.loci.exists() and len(done) == len(plan) and (t.wdir / "merge.key").exists() else "todo"
        st["pred_left_s"] = sum(predict_shard_s(t, n, plan) for lab, n in plan.items() if lab not in done)
    return st


def ensure_target(t: Target, budget: Budget) -> Path:
    """Run prepare / index / every shard / merge of one target, within the budget. Returns liftoff_loci.tsv or raises
    Pending (exit 75) when work is left."""
    if not budget.fits(240):
        pending(f"{t.name}: prepare")
    plan, worked = prepare(t)
    if worked:
        budget.used()
    ensure_genome_link(t)
    if index_state(t) != "done":
        if not budget.fits(420):
            pending(f"{t.name}: index")
        ensure_index(t)
        budget.used()
    for lab, n in plan.items():
        if shard_done(t, lab):
            continue
        pred = predict_shard_s(t, n, plan)
        if not budget.fits(pred):
            pending(f"{t.name}: shard {lab} ({n:,} records, predicted {pred:.0f} s)")
        run_shard(t, lab)
        budget.used()
    return merge(t)


def ensure_species(cfg: dict, species: str, budget: Budget) -> Path:
    return ensure_target(species_target(cfg, species), budget)


def loci_path(cfg: dict, species: str) -> Path | None:
    t = species_target(cfg, species)
    st = target_state(t)
    return t.loci if st["merge"] == "done" else None


# ---------------------------------------------------------------- V-L1: one run vs per-contig shards (+ blocks)
VL1_CONTIGS = ["chr20", "chr21", "chr22"]


def validate(cfg: dict, budget: Budget) -> dict:
    """Mini-genome chr20+chr21+chr22 of CHM13: arm A one Liftoff run; arm B per-contig shards + merge; arm C
    per-contig shards with blocks (chr20 forced into >= 2 blocks) + merge. Bar: annotated placements identical,
    extra copies Jaccard >= 0.98 (prereg §1)."""
    inp = species_inputs(cfg, "human")
    root = liftoff_root(cfg) / "vl1"
    root.mkdir(parents=True, exist_ok=True)
    fa = root / "mini.fa"
    if not (fa.exists() and Path(str(fa) + ".fai").exists()):
        subprocess.run(f"samtools faidx {inp['fasta']} {' '.join(VL1_CONTIGS)} > {fa}.tmp && mv {fa}.tmp {fa} && "
                       f"samtools faidx {fa}", shell=True, check=True)
    arms = {"A_single": dict(single=True, block_records=10 ** 9),
            "B_contigs": dict(single=False, block_records=10 ** 9),
            "C_blocks": dict(single=False, block_records=600)}
    res = {}
    for arm, kw in arms.items():
        t = Target(f"vl1_{arm}", fa, inp["gff"], root / arm, contigs=VL1_CONTIGS,
                   hard_s=int(cfg.get("fig8_liftoff_hard_s", HARD_LIMIT_S)), **kw)
        if arm != "A_single":   # share one index (identical genome and options)
            a_mmi = root / "A_single" / "genome.fa.mmi"
            ensure_genome_link(t)
            if a_mmi.exists() and (root / "A_single" / "index.ok").exists() and not t.mmi.exists():
                t.mmi.symlink_to(a_mmi)
                (t.wdir / "index.ok").write_text("shared with arm A\n")
        res[arm] = read_loci(ensure_target(t, budget))
    def keyset(rows, copies):
        out = set()
        for r in rows:
            if copies and r["cls"] == "extra_copy":
                out.add((r["source_id"], r["contig"], r["start0"], r["end"], r["strand"]))
            if not copies and r["cls"] in ("in_place", "moved", "partial"):
                out.add((r["source_id"], r["contig"], r["start0"], r["end"], r["strand"]))
        return out
    report = {}
    A = res["A_single"]
    for arm in ("B_contigs", "C_blocks"):
        X = res[arm]
        ann_a, ann_x = keyset(A, False), keyset(X, False)
        cop_a, cop_x = keyset(A, True), keyset(X, True)
        jac = len(cop_a & cop_x) / max(1, len(cop_a | cop_x))
        report[arm] = {"annotated_A": len(ann_a), "annotated_arm": len(ann_x), "annotated_identical": ann_a == ann_x,
                       "annotated_only_A": sorted(ann_a - ann_x)[:20], "annotated_only_arm": sorted(ann_x - ann_a)[:20],
                       "copies_A": len(cop_a), "copies_arm": len(cop_x), "copies_both": len(cop_a & cop_x),
                       "copies_jaccard": round(jac, 4), "copies_only_A": sorted(cop_a - cop_x)[:20],
                       "copies_only_arm": sorted(cop_x - cop_a)[:20],
                       "pass": ann_a == ann_x and jac >= 0.98}
    (root / "vl1_report.json").write_text(json.dumps(report, indent=1))
    return report


# ---------------------------------------------------------------- read support per sample (prereg §4 C2)
def support_path(cfg: dict, sid: str) -> Path:
    d = liftoff_root(cfg) / "support"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{sid}.tsv"


def _support_key(lp: Path, bam) -> dict:
    return {"loci": fp(lp), "bam": fp(bam), "rule": f">= {SUPPORT_READS} primary reads (-F 2308) with a block on the "
                                                    f"exon union", "v": 1}


def ensure_support(cfg: dict, sid: str, species: str, budget: Budget) -> Path | None:
    """For every scorable Liftoff locus (classes in_place / moved / extra_copy, exon union >= 200 bp): the number of
    reads (capped at SUPPORT_READS) whose primary alignment (-F 2308) has an aligned block on the exon union. Cached
    per contig (`<sample>.parts/<contig>.tsv`), budgeted."""
    import pysam
    import samples
    lp = loci_path(cfg, species)
    if lp is None:
        return None
    out = support_path(cfg, sid)
    bam = samples.get(cfg, sid)["bam"]
    key = _support_key(lp, bam)
    kpath = Path(str(out) + ".key")
    if out.exists() and kpath.exists() and json.loads(kpath.read_text()) == key:
        return out
    parts = Path(str(out) + ".parts")
    parts.mkdir(exist_ok=True)
    kp = parts / "key.json"
    if not kp.exists() or json.loads(kp.read_text()) != key:
        shutil.rmtree(parts)
        parts.mkdir()
        kp.write_text(json.dumps(key))
    loci = [r for r in read_loci(lp) if r["cls"] in ("in_place", "moved", "extra_copy") and not r["short"]]
    by_contig = collections.defaultdict(list)
    for i, r in enumerate(loci):
        by_contig[r["contig"]].append(i)
    rate = []
    for tf in parts.glob("*.time"):
        a = json.loads(tf.read_text())
        if a["n"]:
            rate.append(a["s"] / a["n"])
    per = (sorted(rate)[len(rate) // 2] * 1.5) if rate else 0.01
    with pysam.AlignmentFile(bam) as bf:
        names = set(bf.references)
        for contig in sorted(by_contig, key=lambda c: -len(by_contig[c])):
            pf = parts / f"{contig}.tsv"
            if pf.exists():
                continue
            idxs = by_contig[contig]
            if not budget.fits(per * len(idxs)):
                pending(f"support {sid}: contig {contig} ({len(idxs):,} loci)")
            t0 = time.time()
            rows = []
            for i in idxs:
                r = loci[i]
                iv = r["iv"]
                n = 0
                if contig in names and iv:
                    for rd in bf.fetch(contig, iv[0][0], iv[-1][1]):
                        if rd.flag & 2308:
                            continue
                        if any(iv_inter([(bs, be)], iv) > 0 for bs, be in rd.get_blocks()):
                            n += 1
                            if n >= SUPPORT_READS:
                                break
                rows.append((r["source_id"], r["cls"], r["contig"], r["start0"], r["end"], n))
            tmp = pf.with_suffix(".tmp")
            with open(tmp, "w") as fo:
                for x in rows:
                    fo.write("\t".join(map(str, x)) + "\n")
            tmp.replace(pf)
            (parts / f"{contig}.time").write_text(json.dumps({"s": time.time() - t0, "n": len(idxs)}))
            budget.used()
    tmp = out.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write("source_id\tcls\tcontig\tstart0\tend\treads_capped\n")
        for contig in by_contig:
            with open(parts / f"{contig}.tsv") as fh:
                fo.write(fh.read())
    tmp.replace(out)
    kpath.write_text(json.dumps(key))
    return out


def read_support(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["source_id"], row["cls"], row["contig"], row["start0"], row["end"])] = int(row["reads_capped"])
    return out


def locus_key(r: dict) -> tuple:
    return (r["source_id"], r["cls"], r["contig"], r["start0"], r["end"])


# ---------------------------------------------------------------- Rustle loci (read-only; run-cache products)
def product_if_ready(cfg: dict, sid: str, stage: str, name: str) -> tuple[Path | None, str]:
    import samples
    state, reason = samples.status(cfg, sid, stage)
    p = samples.product(cfg, sid, stage, name)
    if state in ("fresh", "adopt") and p.exists():
        return p, state
    return None, f"{stage} {state}: run `python3 figures/make.py runs --sample {sid} --stage {stage}`"


def denovo_loci(gtf: Path, cache_dir: Path) -> list[dict]:
    """gene_id groups of an assembled GTF: contig, span, exon union of every transcript. Cached by the GTF's
    fingerprint (a genome-wide human GTF is ~0.5 GB)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cp = cache_dir / (hashlib.sha1(fp(gtf).encode()).hexdigest()[:16] + ".denovo_loci.tsv")
    if not cp.exists():
        ex = collections.defaultdict(list)
        meta = {}
        opener = gzip.open if str(gtf).endswith(".gz") else open
        with opener(gtf, "rt") as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.split("\t", 8)
                if len(f) < 9 or f[2] != "exon":
                    continue
                m = re.search(r'gene_id "([^"]+)"', f[8])
                if not m:
                    continue
                g = (f[0], m.group(1))
                ex[g].append((int(f[3]) - 1, int(f[4])))
                meta.setdefault(g, f[6])
        tmp = cp.with_suffix(".tmp")
        with open(tmp, "w") as fo:
            fo.write("contig\tlocus\tstrand\tstart0\tend\texons\n")
            for (contig, g), iv in ex.items():
                iv = merge_iv(iv)
                fo.write(f"{contig}\t{g}\t{meta[(contig, g)]}\t{iv[0][0]}\t{iv[-1][1]}\t{iv_str(iv)}\n")
        tmp.replace(cp)
    with open(cp) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["iv"] = iv_parse(r["exons"])
    return rows


def catalog_loci(copies_tsv: Path) -> list[dict]:
    rows = []
    with open(copies_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            iv = merge_iv(iv_parse(r["exons"])) if r.get("exons") else [(int(r["start"]), int(r["end"]))]
            rows.append({"contig": r["chrom"], "locus": f"{r['family_id']}:{r['copy_idx']}", "family": r["family_id"],
                         "strand": r.get("strand", "."), "start0": int(r["start"]), "end": int(r["end"]), "iv": iv})
    return rows


# ---------------------------------------------------------------- matching (prereg §3)
class Index:
    """Loci of one set by contig, binned for overlap queries."""
    BIN = 100_000

    def __init__(self, rows: list[dict]):
        self.rows = rows
        self.bins = collections.defaultdict(list)
        for i, r in enumerate(rows):
            if not r["iv"]:
                continue
            for b in range(r["iv"][0][0] // self.BIN, (r["iv"][-1][1] - 1) // self.BIN + 1):
                self.bins[(r["contig"], b)].append(i)

    def near(self, contig: str, iv) -> list[int]:
        if not iv:
            return []
        seen, out = set(), []
        for b in range(iv[0][0] // self.BIN, (iv[-1][1] - 1) // self.BIN + 1):
            for i in self.bins.get((contig, b), ()):
                if i not in seen:
                    seen.add(i)
                    out.append(i)
        return out

    def best_cov(self, contig: str, iv) -> tuple[float, int | None]:
        """max over R of cov(G | R) = |X(G) ∩ X(R)| / |X(G)| (G = iv), and that R's index."""
        L = iv_len(iv)
        best, bi = 0.0, None
        if L == 0:
            return best, bi
        for i in self.near(contig, iv):
            c = iv_inter(iv, self.rows[i]["iv"]) / L
            if c > best:
                best, bi = c, i
        return best, bi


def found_by(ref_rows: list[dict], query_rows: list[dict]) -> list[tuple[float, int | None]]:
    idx = Index(query_rows)
    return [idx.best_cov(r["contig"], r["iv"]) for r in ref_rows]


# ---------------------------------------------------------------- Liftoff copy pairs vs a family table (Figs 6s, 7, 8)
# The families Figures 6s, 7 and 8 score are the DEFAULT de novo families: the driver's `families` stage copy table
# (<id>.fam.copies.tsv; one copy per member locus, its representative's exons), read by catalog_loci (the same columns
# as the legacy gw_family_catalog copies.tsv, which the same code scores as a labelled secondary set).
def support_if_ready(cfg: dict, sid: str, species: str) -> tuple[list | None, dict | None, str]:
    """(Liftoff loci, the sample's read-support table, "") when the species' self-lift is merged and the sample's
    support table is complete for it; else (None, None, the command that builds what is missing). Runs nothing."""
    import samples
    lp = loci_path(cfg, species)
    if lp is None:
        return None, None, (f"Liftoff self-lift of {species} not merged: python3 figures/_liftoff.py run --species "
                            f"{species} (or make.py data fig8)")
    out = support_path(cfg, sid)
    kpath = Path(str(out) + ".key")
    if not (out.exists() and kpath.exists()
            and json.loads(kpath.read_text()) == _support_key(lp, samples.get(cfg, sid)["bam"])):
        return None, None, f"read support of {sid} on the {species} Liftoff loci absent: make.py data fig8"
    return read_loci(lp), read_support(out), ""


def copy_pairs(loci: list[dict], sc: float, sup: dict | None = None, keep=None) -> list[tuple[dict, dict]]:
    """Liftoff (source record's annotated placement, extra copy) pairs (prereg C3): the extra copy's sequence_ID >= sc,
    both exon unions >= SHORT_BP, the annotated placement in place or moved. With `sup`: both loci read-supported
    (>= SUPPORT_READS reads of the sample; the C2 rule). `keep(contig)`: both loci on kept contigs."""
    ann = {r["source_id"]: r for r in loci if r["cls"] in ("in_place", "moved") and not r["short"]}
    out = []
    for r in loci:
        if r["cls"] != "extra_copy" or r["short"] or not r["seqid_f"] >= sc - 1e-9:
            continue
        a = ann.get(r["source_id"])
        if a is None or (keep is not None and not (keep(a["contig"]) and keep(r["contig"]))):
            continue
        if sup is not None and (sup.get(locus_key(a), 0) < SUPPORT_READS or sup.get(locus_key(r), 0) < SUPPORT_READS):
            continue
        out.append((a, r))
    return out


def pair_families(pairs: list[tuple[dict, dict]], fam_loci: list[dict]) -> list[tuple[str, bool, set]]:
    """Per pair: (source_id, both loci covered, families shared). A Liftoff locus G is covered by a family member R
    when cov(G | R) >= MATCH_COV (Liftoff's -a on exon bases, prereg §3); the pair is in one family when a covering
    member of each locus share a family (`family` of catalog_loci rows)."""
    idx = Index(fam_loci)

    def fams(g):
        n = iv_len(g["iv"])
        return {fam_loci[i]["family"] for i in idx.near(g["contig"], g["iv"])
                if n and iv_inter(g["iv"], fam_loci[i]["iv"]) / n >= MATCH_COV}
    out = []
    for a, c in pairs:
        fa, fc = fams(a), fams(c)
        out.append((a["source_id"], bool(fa and fc), fa & fc))
    return out


# ---------------------------------------------------------------- bootstrap over source records
def boot_ci(clusters: list[str], hits: list[int], n=BOOT_N, seed=BOOT_SEED) -> tuple[float, float]:
    """95% interval of sum(hits) / len(hits), resampling source records (clusters) with replacement."""
    import numpy as np
    if not hits:
        return float("nan"), float("nan")
    ids = {}
    cn = collections.Counter()
    ch = collections.Counter()
    for c, h in zip(clusters, hits):
        k = ids.setdefault(c, len(ids))
        cn[k] += 1
        ch[k] += h
    K = len(ids)
    N = np.array([cn[k] for k in range(K)], dtype=float)
    H = np.array([ch[k] for k in range(K)], dtype=float)
    rng = np.random.default_rng(seed)
    vals = []
    step = max(1, min(n, 2_000_000 // max(K, 1)))
    done = 0
    while done < n:
        m = min(step, n - done)
        idx = rng.integers(0, K, size=(m, K))
        num, den = H[idx].sum(axis=1), N[idx].sum(axis=1)
        vals.append(num / np.maximum(den, 1))
        done += m
    v = np.concatenate(vals)
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


# ---------------------------------------------------------------- Rustle guided candidate search (prereg §4 C1)
# Seeds = EVERY gene/pseudogene record; otherwise bench/guided_pipeline.py unchanged: the record's unit (longest
# NM_/NR_ transcript, else any transcript, else its own exons, else its span) aligned with -x splice, its CDS envelope
# (else exon envelope, else span) with -x asm20, both `-c -N 100 -p 0.1` against the species' splice index;
# gp.transcript_hits (identity >= 0.80, >= 50% of the unit), gp.gene_body_chains (identity >= 0.80, >= 50% of
# min(query, extrapolated span)), gp.chain_first (hits overlapping any seed's span are blocked).
FINDER_TX_TYPES = ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript")   # guided_pipeline.Annotation
FINDER_FLAGS = ["-c", "-N", "100", "-p", "0.1"]
FINDER_UNITS_PER_SHARD = 10000
FINDER_ENV_MB_PER_SHARD = 25


def _gp():
    if str(GP_DIR) not in sys.path:
        sys.path.insert(0, str(GP_DIR))
    import guided_pipeline as gp
    return gp


def finder_dir(t: Target) -> Path:
    d = t.wdir / "finder"
    d.mkdir(parents=True, exist_ok=True)
    return d


def finder_key(cfg: dict, t: Target, mmi: Path) -> dict:
    mm2 = cfg.get("minimap2", "minimap2")
    ver = subprocess.run([mm2, "--version"], capture_output=True, text=True).stdout.strip()
    return {"gff": fp(t.gff), "genome": fp(t.fasta), "mmi": fp(mmi), "minimap2": f"{mm2} {ver}", "flags": FINDER_FLAGS,
            "units_per_shard": int(cfg.get("fig8_finder_units", FINDER_UNITS_PER_SHARD)),
            "env_mb": int(cfg.get("fig8_finder_env_mb", FINDER_ENV_MB_PER_SHARD)), "v": 1}


def finder_prepare(cfg: dict, t: Target, fd: Path, key: dict) -> list[tuple[str, str, Path]]:
    """Unit and envelope FASTAs of every record (guided_pipeline.Annotation semantics, keyed by record ID), split into
    query shards; returns [(label, preset, fasta)]."""
    import pysam
    gp = _gp()
    kpath = fd / "prepare.key"
    lpath = fd / "shards.tsv"
    if kpath.exists() and json.loads(kpath.read_text()) == key and lpath.exists():
        return [tuple(l.rstrip("\n").split("\t")[:2]) + (fd / l.rstrip("\n").split("\t")[2],)
                for l in open(lpath)]
    recs = read_records(t)
    rid = {r["id"]: r for r in recs}
    tx_parent, tx_name = {}, {}
    exons, cds = collections.defaultdict(list), collections.defaultdict(list)
    with open(t.gff) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            if f[2] in FINDER_TX_TYPES:
                a = attrs(f[8])
                if a.get("Parent", "").split(",")[0] in rid:
                    tx_parent[a["ID"]] = a["Parent"].split(",")[0]
                    tx_name[a["ID"]] = a.get("Name", a["ID"])
            elif f[2] in ("exon", "CDS"):
                a = attrs(f[8])
                par = a.get("Parent", "").split(",")[0]
                (exons if f[2] == "exon" else cds)[par].append((int(f[3]) - 1, int(f[4])))
    txs_of = collections.defaultdict(list)
    for tx, g in tx_parent.items():
        if exons.get(tx):
            txs_of[g].append(tx)
    genome = pysam.FastaFile(str(t.fasta))
    units, envs = [], []
    for r in recs:
        txs = txs_of.get(r["id"], [])
        span = lambda tx: sum(e - s for s, e in gp.merge(exons[tx]))
        pool = [tx for tx in txs if tx_name[tx].startswith(("NM_", "NR_"))] or txs
        if pool:
            tx = max(pool, key=lambda x: (span(x), x))
            blocks = gp.merge(exons[tx])
        else:
            blocks = gp.merge(exons.get(r["id"], []))
        if blocks:
            seq = "".join(genome.fetch(r["contig"], s, e) for s, e in blocks).upper()
        else:
            seq = genome.fetch(r["contig"], r["start0"], r["end"]).upper()
        units.append((r["id"], gp.rc(seq) if r["strand"] == "-" else seq))
        cd = [c for tx in txs for c in cds.get(tx, [])]
        if cd:
            env = (min(s for s, _ in cd), max(e for _, e in cd))
        elif blocks:
            env = (blocks[0][0], blocks[-1][1])
        else:
            env = (r["start0"], r["end"])
        eseq = genome.fetch(r["contig"], env[0], env[1]).upper()
        envs.append((r["id"], gp.rc(eseq) if r["strand"] == "-" else eseq))
    shards = []
    per = key["units_per_shard"]
    for i in range(0, len(units), per):
        shards.append((f"u{len(shards):03d}", "splice", units[i:i + per]))
    cap, cur, bp = key["env_mb"] * 1_000_000, [], 0
    ne = 0
    for q in envs:
        if cur and bp + len(q[1]) > cap:
            shards.append((f"e{ne:03d}", "asm20", cur))
            ne, cur, bp = ne + 1, [], 0
        cur.append(q)
        bp += len(q[1])
    if cur:
        shards.append((f"e{ne:03d}", "asm20", cur))
    out = []
    with open(lpath, "w") as lf:
        for lab, preset, seqs in shards:
            fa = fd / f"{lab}.fa"
            with open(fa, "w") as fo:
                for n, sq in seqs:
                    fo.write(f">{n}\n{sq}\n")
            lf.write(f"{lab}\t{preset}\t{fa.name}\t{len(seqs)}\t{sum(len(x[1]) for x in seqs)}\n")
            out.append((lab, preset, fa))
    kpath.write_text(json.dumps(key))
    log(f"{t.name}: finder queries: {len(units):,} units, {len(envs):,} envelopes in {len(shards)} shards")
    return out


def _finder_shard_done(fd: Path, lab: str, key: dict) -> bool:
    k = fd / f"{lab}.done"
    return k.exists() and json.loads(k.read_text()) == key


def finder_shard(cfg: dict, t: Target, fd: Path, lab: str, preset: str, fa: Path, mmi: Path, key: dict):
    gp = _gp()
    mm2 = cfg.get("minimap2", "minimap2")
    paf = fd / f"{lab}.paf"
    cmd = [mm2] + FINDER_FLAGS[:1] + ["-x", preset] + FINDER_FLAGS[1:] + ["-t", str(LIFT_THREADS), "-o", str(paf) + ".tmp",
                                                                       str(mmi), str(fa)]
    r = _timed(cmd, fd / f"{lab}.time", fd / f"{lab}.log", cwd=fd, hard_s=t.hard_s)
    if r != 0:
        Path(str(paf) + ".tmp").unlink(missing_ok=True)
        raise RuntimeError(f"{t.name} finder {lab}: minimap2 exit {r} (124 = over {t.hard_s} s: lower "
                           f"fig8_finder_units / fig8_finder_env_mb); see {fd / (lab + '.log')}")
    Path(str(paf) + ".tmp").replace(paf)
    out = fd / f"{lab}.hits.tsv"
    with open(out, "w") as fo:
        if preset == "splice":
            for h in gp.transcript_hits(str(paf)):
                bl = gp.tx_exon_blocks(h)
                fo.write("\t".join(map(str, ["tx", h["q"], h["chrom"], h["strand"], h["s"], h["e"], h["nm"], h["bl"],
                                             "-", "-", iv_str(merge_iv(bl))])) + "\n")
        else:
            for c in gp.gene_body_chains(str(paf)):
                bl = sum(x["bl"] for x in c["recs"])
                fo.write("\t".join(map(str, ["chain", c["q"], c["chrom"], c["strand"], c["s"], c["e"], c["nm"], bl,
                                             c["xs"], c["xe"], "-"])) + "\n")
    paf.unlink()
    (fd / f"{lab}.done").write_text(json.dumps(key))
    fa.unlink(missing_ok=True)   # regenerated by finder_prepare if the key ever changes


def ensure_finder(cfg: dict, species: str, budget: Budget) -> Path:
    """Rustle guided candidates of one species: finder/candidates.tsv (contig, start0, end, source record, kind,
    identity of the leading hit, exons). HEAVY (splice index ~13 GB per mapping call); budgeted per query shard."""
    import samples
    t = species_target(cfg, species)
    if not t.records.exists():
        pending(f"{species}: the Liftoff prepare step (records.tsv) comes first")
    sid = next(s for s, r in samples.registry(cfg).items() if r["species"] == species)
    mmi = samples.product(cfg, sid, "index", "mmi")
    if not Path(mmi).exists():
        pending(f"{species}: splice index {mmi} missing (make.py runs --sample {sid} --stage index)")
    fd = finder_dir(t)
    key = finder_key(cfg, t, Path(mmi))
    cand = fd / "candidates.tsv"
    if cand.exists() and (fd / "candidates.key").exists() and json.loads((fd / "candidates.key").read_text()) == key:
        return cand
    if not budget.fits(300):
        pending(f"{species}: finder prepare")
    shards = finder_prepare(cfg, t, fd, key)
    rates = collections.defaultdict(list)
    sizes = {l.split("\t")[0]: int(l.rstrip("\n").split("\t")[4]) for l in open(fd / "shards.tsv")}
    for lab, preset, fa in shards:
        tm = read_time(fd / f"{lab}.time")
        if _finder_shard_done(fd, lab, key) and tm.get("wall_s"):
            rates[preset].append(tm["wall_s"] / max(1, sizes[lab]))
    for lab, preset, fa in shards:
        if _finder_shard_done(fd, lab, key):
            continue
        r = sorted(rates[preset])
        pred = (r[len(r) // 2] * sizes[lab] * 1.25) if r else 600.0
        if not budget.fits(pred):
            pending(f"{species}: finder shard {lab} ({preset}, {sizes[lab] / 1e6:.0f} Mb, predicted {pred:.0f} s)")
        log(f"{species}: finder shard {lab} ({preset}, {sizes[lab] / 1e6:.0f} Mb)")
        finder_shard(cfg, t, fd, lab, preset, fa, Path(mmi), key)
        budget.used()
    # candidates: gp.chain_first per contig, seeds = every record, blocked_by on an interval index
    gp = _gp()
    recs = read_records(t)
    rec = {r["id"]: {"chrom": r["contig"], "start0": r["start0"], "end": r["end"], "family": r["id"]} for r in recs}
    spans = Index([{"contig": r["contig"], "iv": [(r["start0"], r["end"])]} for r in recs])

    def blocked(seeds, rec_, h):
        return any(ov(spans.rows[i]["iv"][0][0], spans.rows[i]["iv"][0][1], h["s"], h["e"]) > 0
                   for i in spans.near(h["chrom"], [(h["s"], h["e"])]))
    gp.blocked_by = blocked
    tx_by, ch_by = collections.defaultdict(list), collections.defaultdict(list)
    for lab, preset, fa in shards:
        for line in open(fd / f"{lab}.hits.tsv"):
            k, q, chrom, strand, s0, e0, nm, bl, xs, xe, ex = line.rstrip("\n").split("\t")
            h = {"q": q, "chrom": chrom, "strand": strand, "s": int(s0), "e": int(e0), "nm": int(nm), "bl": int(bl),
                 "blocks": ex}
            if k == "tx":
                tx_by[chrom].append(h)
            else:
                h.update(xs=int(xs), xe=int(xe))
                ch_by[chrom].append(h)
    seeds = set(rec)
    rows = []
    for chrom in sorted(set(tx_by) | set(ch_by)):
        for c in gp.chain_first(tx_by.get(chrom, []), ch_by.get(chrom, []), seeds, rec):
            lead = c["tx"] or c["chain"]
            ex = c["tx"]["blocks"] if c["tx"] else iv_str([(c["s"], c["e"])])
            rows.append((c["chrom"], c["s"], c["e"], c["family"], "tx" if c["tx"] else "chain",
                         f"{lead['nm'] / max(1, lead['bl']):.4f}", ex))
    tmp = cand.with_suffix(".tmp")
    with open(tmp, "w") as fo:
        fo.write("contig\tstart0\tend\tsource_id\tkind\tidentity\texons\n")
        for r in sorted(rows, key=lambda x: (x[0], x[1])):
            fo.write("\t".join(map(str, r)) + "\n")
    tmp.replace(cand)
    (fd / "candidates.key").write_text(json.dumps(key))
    log(f"{species}: {len(rows):,} guided candidate loci")
    return cand


def read_candidates(path: Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["iv"] = merge_iv(iv_parse(r["exons"]))
        r["identity"] = float(r["identity"])
    return rows


# ---------------------------------------------------------------- plan / CLI
def plan_report(cfg: dict) -> list[str]:
    import samples
    out = []
    vl1 = liftoff_root(cfg) / "vl1" / "vl1_report.json"
    if vl1.exists():
        rep = json.loads(vl1.read_text())
        out.append("V-L1 " + "; ".join(f"{a}: {'PASS' if r['pass'] else 'FAIL'} (annotated identical "
                                        f"{r['annotated_identical']}, copies Jaccard {r['copies_jaccard']})"
                                        for a, r in rep.items()))
    else:
        out.append("V-L1 not done: python3 figures/_liftoff.py validate --budget-s 540 (repeat while exit 75)")
    for sp in species_list(cfg):
        t = species_target(cfg, sp)
        st = target_state(t)
        pred = st.get("pred_left_s")
        sp_json = t.wdir / "splits.json"
        splits = json.loads(sp_json.read_text()) if sp_json.exists() else []
        fd = t.wdir / "finder" / "candidates.tsv"
        out.append(f"{sp:11s} prepare {st['prepare']:4s} index {st['index']:4s} shards {st['shards_done']}/"
                   f"{st['shards'] or '?'} records left {st['records_left']:,} merge {st['merge']}"
                   + (f" (~{pred / 60:.0f} min left)" if pred else "")
                   + (f"; split after the hard limit: {', '.join(splits)}" if splits else "")
                   + f"; guided finder {'done' if fd.exists() else 'not run'}")
    for sid, row in samples.registry(cfg).items():
        sp = row["species"]
        sp_done = loci_path(cfg, sp) is not None
        sup = support_path(cfg, sid)
        out.append(f"  {sid:15s} support {'done' if sup.exists() else ('todo' if sp_done else 'waits for ' + sp)}")
    return out


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "run", "support", "validate", "finder"):
        p = sub.add_parser(name)
        p.add_argument("--inputs")
        p.add_argument("--budget-s", type=float, default=0)
        p.add_argument("--max-units", type=int)
        p.add_argument("--species")
        p.add_argument("--sample")
    a = ap.parse_args(argv)
    cfg = figlib.load_inputs(a.inputs)
    budget = Budget(a.budget_s, a.max_units)
    if a.cmd == "plan":
        print("\n".join(plan_report(cfg)))
        return
    if a.cmd == "run":
        for sp in ([a.species] if a.species else species_list(cfg)):
            ensure_species(cfg, sp, budget)
        return
    if a.cmd == "support":
        import samples
        sid = samples.resolve(cfg, a.sample)
        p = ensure_support(cfg, sid, samples.get(cfg, sid)["species"], budget)
        print(p or "the species' Liftoff table is not finished")
        return
    if a.cmd == "validate":
        print(json.dumps(validate(cfg, budget), indent=1))
        return
    if a.cmd == "finder":
        print(ensure_finder(cfg, a.species, budget))
        return


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
