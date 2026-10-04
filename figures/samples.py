"""samples — the sample registry and the genome-wide run cache (`make.py samples`, `make.py runs`).

REGISTRY. `figures/samples.tsv` (inputs key `samples`) holds one row per sample. Its cells are `${key}` references into the
inputs file, which holds every machine path, or literal text; `-` or an empty cell means "none". Columns:

    id           sample id (<species>_<individual or library>), e.g. gorilla_OR6737
    alias        the legacy species key the figure modules pass ('gorilla', 'human'), or '-'
    species      human | gorilla | chimpanzee | orangutan  (numbers from different species are never pooled)
    individual   the animal / library, as far as it is known
    tissue       tissue or cell type; 'unknown' when no record says so
    tissue_source  where the tissue statement comes from (BAM header, docs/, file names)
    genome       reference assembly the BAM is aligned to
    bam fasta annotation_gff annotation_gtf splice_mmi   inputs of the pipeline stages (gff: uncompressed; flag needs it)
    flag_confirm comma-separated NAME=PATH.mmi genomes the flag stage confirms against (driver --confirm), or '-'
    stringtie_gtf flair_gtf isoseq_gff   the lab's ANNOTATION-FREE (de novo) transcript sets (only where the lab ran
                 them), or '-': StringTie -L without -G; FLAIR bam2bed + collapse without -f/--gtf (flair correct
                 skipped); IsoSeq collapse (benchmark_collapse/run_stringtie.sbatch, run_flair.sbatch)
    stringtie_guided_gtf flair_guided_gtf   (optional columns) ANNOTATION-GUIDED runs of the same tools (StringTie -G,
                 FLAIR with the annotation), supplied by the user; '-' for every sample until then. They feed only the
                 separate guided path of figures 1-3 (assembly.guided_samples), never a panel with the de novo methods.
                 There is no Rustle counterpart column: Rustle has no annotation-guided transcript assembly
                 (docs/PREREG_guided_transcript_comparison_2026-09-25.md).

    registry(cfg)               {id: row} with paths expanded, in file order
    resolve(cfg, key)           sample id for an id or an alias ('gorilla' -> 'gorilla_OR6737', 'human' -> 'human_A119b')
    get(cfg, key)               the row of one sample
    baseline(cfg, key, tool)    the lab's stringtie | flair | isoseq transcript set (annotation-free), or the user's
                                stringtie_guided | flair_guided one (annotation-guided), or None
    verify(cfg, key)            [(level, message)]: files, indexes, BAM @SQ vs FASTA .fai vs annotation seqids vs .mmi

RUN CACHE. Each stage of `tools/rustle_pipeline.sh` runs GENOME-WIDE on the whole BAM of a sample, into
`${work}/runs/<id>/` with PREFIX `<id>` (assemble_primary: `<id>.primary`; the driver's PREFIX.cache/ stays on):

    stage             driver call                                   products (name: path)
    assemble          assemble                                      gtf <id>.gtf, molecules <id>.molecules.tsv
    assemble_primary  assemble --no-seed-secondaries                gtf <id>.primary.gtf
    families          families (de novo mode, needs assemble)       clusters <id>.fam.clusters.tsv, loci, params,
                      = THE default family definition               copies <id>.fam.copies.tsv, copies_fa (the copy table
                                                                    copy assignment consumes)
    families_primary  families on <id>.primary.gtf (needs assemble_primary; PREFIX <id>.primary)
                                                                    clusters <id>.primary.fam.clusters.tsv, loci, params
    catalog           catalog (LEGACY copy catalog)                 copies <id>.cat.copies.tsv, copies_fa, families, pairs
    assign            assign (needs families: the driver's assign   assignments <id>.assign.assignments.tsv
                      reads the families' copy table since
                      2026-10-02; the opt-in candidates stage is not run)
    index             minimap2 -x splice -d (only when the registry's splice_mmi file is absent)   mmi
    flag              flag --index splice_mmi --gff annotation_gff [--confirm ...]   scan, calls <id>.flag.missing_copy.tsv

    product(cfg, key, stage, name=None)   a product's path (runs nothing; None-safe for a figure that only reads)
    status(cfg, key, stage)               (state, reason): fresh | adopt | run | stale | blocked
    ensure(cfg, key, stage, force=False)  run the stage if it is not fresh (HEAVY; foreground), return {name: path};
                                          raises StagePending when a bounded call made progress and must be repeated
    restamp(cfg, key, stage, proof)       re-stamp a STALE stage without running it, from a proof file of
                                          cmp-identical products (the binaries changed, the outputs provably did not)

BOUNDED CALLS. `catalog` runs `--piecewise` with a per-call budget (exit 75 while pieces or the merge remain), and the
stages with a genome-wide minimap2 all-vs-all (`families`, `families_primary`, `catalog`) run it through the shard wrapper
`tools/mm2_shard.sh` (RUSTLE_MINIMAP2 in the child's environment, inputs key `runs_mm2_wrapper`; `-` = plain minimap2):
MM2_SHARD_BUDGET_S (`runs_mm2_budget_s`, default 480) and MM2_SHARD_DEADLINE (call start + `runs_call_budget_s`,
default 560) bound one call unless the caller's environment sets them. The wrapper's shards concatenate to the single
run byte for byte (cmp-checked 2026-09-25 on chr16 for mcl_families and gw_family_catalog), so the wrapper is NOT part of
the stage key: a stamp written with or without it stays valid (older stamps that recorded it in `env` are compared
with it removed). A call whose all-vs-all ran out of budget after making progress raises StagePending (`make.py runs`
exits 75: repeat the same command); a call that made no progress fails, so a resume loop cannot spin.

A stage is FRESH when its products exist and its stamp `<stage>.key` matches the current key: the BAM / FASTA / index /
annotation fingerprints (path, size, mtime), the sha1 of every binary the stage runs, the driver's CODE hash
of that stage (driver_stage_code: the stage's own function plus the shared code, comment edits and other stages'
edits do not re-run; stamps from before it are mapped through LEGACY_DRIVER_CODE), the minimap2 version where the stage aligns, the RUSTLE_* variables
of the environment, and the key of the upstream stage. Products without a stamp are ADOPTED once when they are newer than
every input and binary (the rule assembly.py always used). Every run records wall time and peak RSS
(`/usr/bin/time`) in `<stage>.time` and the driver output in `<stage>.driver.log`.

The two legacy samples keep their existing assemblies: `${work}/assembly/<alias>/rustle.genome.*` and
`rustle_primary.genome.*` are MOVED (hard link + atomic symlink swap: same inode, zero bytes copied) to
`${work}/runs/<id>/<id>.*` / `<id>.primary.*` the first time a run-cache function looks at that sample; the old names stay
as symlinks, so any reader of the old paths (fig_secondary's molecules table) keeps working.
"""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

import figlib

HERE = Path(__file__).resolve().parent
COLUMNS = ["id", "alias", "species", "individual", "tissue", "tissue_source", "genome", "bam", "fasta", "annotation_gff",
           "annotation_gtf", "splice_mmi", "flag_confirm", "stringtie_gtf", "flair_gtf", "isoseq_gff",
           "stringtie_guided_gtf", "flair_guided_gtf"]
# columns a registry file may leave out (added 2026-09-25; absent = '-' for every sample)
OPTIONAL_COLUMNS = {"stringtie_guided_gtf", "flair_guided_gtf"}
PATH_COLUMNS = ["bam", "fasta", "annotation_gff", "annotation_gtf", "splice_mmi", "stringtie_gtf", "flair_gtf",
                "isoseq_gff", "stringtie_guided_gtf", "flair_guided_gtf"]
# annotation-free (de novo) lab baselines: the tool comparisons of figures 1-3
BASELINES = {"stringtie": "stringtie_gtf", "flair": "flair_gtf", "isoseq": "isoseq_gff"}
# annotation-guided runs of the same tools (user-supplied): the separate guided path of figures 1-3 only
GUIDED_BASELINES = {"stringtie_guided": "stringtie_guided_gtf", "flair_guided": "flair_guided_gtf"}


# ---------------------------------------------------------------- registry
def _expand(v: str, cfg: dict) -> str:
    for k, val in cfg.items():
        if isinstance(val, str):
            v = v.replace("${" + k + "}", val)
    return v


def registry(cfg: dict) -> dict:
    path = Path(cfg.get("samples") or HERE / "samples.tsv")
    with open(path) as fh:
        rows = list(csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t"))
    out = {}
    for r in rows:
        missing = [c for c in COLUMNS if c not in r and c not in OPTIONAL_COLUMNS]
        if missing:
            raise ValueError(f"{path}: row {r.get('id')!r} lacks columns {missing}")
        row = {"_unresolved": []}
        for c in COLUMNS:
            v = _expand((r.get(c) or "").strip(), cfg)
            if "${" in v:  # a key this inputs file does not define: the cell is unusable, verify() reports it
                row["_unresolved"].append(f"{c}={v}")
                v = ""
            row[c] = None if v in ("", "-") else v
        out[row["id"]] = row
    return out


def resolve(cfg: dict, key: str) -> str:
    reg = registry(cfg)
    if key in reg:
        return key
    for sid, row in reg.items():
        if row["alias"] == key:
            return sid
    raise KeyError(f"unknown sample {key!r}; known: {', '.join(reg)} (aliases: "
                   f"{', '.join(r['alias'] for r in reg.values() if r['alias'])})")


def get(cfg: dict, key: str) -> dict:
    return registry(cfg)[resolve(cfg, key)]


def baseline(cfg: dict, key: str, tool: str):
    """A tool's transcript set for one sample, or None: `stringtie` / `flair` / `isoseq` = the lab's annotation-free
    runs; `stringtie_guided` / `flair_guided` = the user's annotation-guided runs (never compared with the others)."""
    return get(cfg, key)[BASELINES.get(tool) or GUIDED_BASELINES[tool]]


def alias_of(cfg: dict, key: str) -> str | None:
    return get(cfg, key)["alias"]


# ---------------------------------------------------------------- verification
def _bam_sq(bam) -> list:
    h = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True).stdout
    out = []
    for line in h.splitlines():
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            out.append((f["SN"], int(f["LN"])))
    return out


def _fai(fasta) -> list:
    p = Path(str(fasta) + ".fai")
    if not p.exists():
        return None
    return [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(p) if l.strip()]


def _cache_dir(cfg) -> Path:
    d = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "runs" / "_meta"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fp(path) -> str:
    """Fingerprint of an input file: resolved path, size, mtime (ns)."""
    p = Path(path)
    try:
        st = p.stat()
        return f"{p.resolve()}|{st.st_size}|{st.st_mtime_ns}"
    except OSError:
        return f"{p}|absent"


def _cached(cfg, tag: str, path, compute):
    """Memoise compute(path) on disk, keyed by the file's fingerprint (seqid scans of GB-sized annotations)."""
    k = hashlib.sha1(f"{tag}\t{_fp(path)}".encode()).hexdigest()[:16]
    f = _cache_dir(cfg) / f"{tag}.{k}.json"
    if f.exists():
        return json.loads(f.read_text())
    v = compute(path)
    f.write_text(json.dumps(v))
    return v


def annotation_seqids(cfg, path) -> dict:
    """{'seqids': [...] in file order, 'lengths': {seqid: len} from ##sequence-region pragmas}; gz-aware."""
    def compute(p):
        cat = f"zcat {shq(p)}" if str(p).endswith(".gz") else f"cat {shq(p)}"
        cmd = cat + " | awk -F'\\t' '/^##sequence-region/{split($0,a,\" \"); print \"L\\t\"a[2]\"\\t\"a[4]; next} " \
                    "/^#/{next} $1!=last{print \"S\\t\"$1; last=$1}'"
        seen, order, lengths = set(), [], {}
        out = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            t = line.split("\t")
            if t[0] == "L" and len(t) == 3 and t[2].isdigit():
                lengths[t[1]] = int(t[2])
            elif t[0] == "S" and t[1] not in seen:
                seen.add(t[1])
                order.append(t[1])
        return {"seqids": order, "lengths": lengths}
    return _cached(cfg, "seqids", path, compute)


def shq(p) -> str:
    return "'" + str(p).replace("'", "'\\''") + "'"


def mmi_header(path) -> dict:
    """minimap2 index header: k, w, flag and the (name, length) of every sequence (mm_idx_dump layout)."""
    with open(path, "rb") as fh:
        if fh.read(4) != b"MMI\x02":
            raise ValueError(f"{path}: not a minimap2 index")
        w, k, b, n_seq, flag = struct.unpack("<5I", fh.read(20))
        seqs = []
        for _ in range(n_seq):
            (ln,) = struct.unpack("<B", fh.read(1))
            name = fh.read(ln).decode()
            (length,) = struct.unpack("<I", fh.read(4))
            seqs.append((name, length))
    return {"k": k, "w": w, "flag": flag, "seqs": seqs}


def bam_records(cfg, bam) -> int:
    """Mapped alignment records (samtools idxstats; index-only, seconds)."""
    def compute(p):
        out = subprocess.run(["samtools", "idxstats", str(p)], capture_output=True, text=True, check=True).stdout
        return sum(int(l.split("\t")[2]) for l in out.splitlines() if l.strip())
    return _cached(cfg, "records", bam, compute)


def _cmp_names(a: list, b: list, la: str, lb: str) -> list:
    """Messages for names/lengths of `a` vs `b` (lists of (name, len) or names)."""
    msgs = []
    da = dict(a) if a and isinstance(a[0], tuple) else {n: None for n in a}
    db = dict(b) if b and isinstance(b[0], tuple) else {n: None for n in b}
    only_a = [n for n in da if n not in db]
    only_b = [n for n in db if n not in da]
    diff_len = [n for n in da if n in db and da[n] is not None and db[n] is not None and da[n] != db[n]]
    if only_a:
        msgs.append(f"{len(only_a)} {la} contig(s) absent from {lb}: {', '.join(only_a[:6])}{' ...' if len(only_a) > 6 else ''}")
    if only_b:
        msgs.append(f"{len(only_b)} {lb} seqid(s) absent from {la}: {', '.join(only_b[:6])}{' ...' if len(only_b) > 6 else ''}")
    if diff_len:
        msgs.append(f"{len(diff_len)} contig(s) with different lengths in {la} vs {lb}: {', '.join(diff_len[:6])}")
    return msgs


def verify(cfg: dict, key: str) -> list:
    """[(OK|WARN|FAIL, message)] for one sample. Light: headers, indexes, one seqid scan per annotation (cached)."""
    row = get(cfg, key)
    res = []
    ok = lambda m: res.append(("OK", m))  # noqa: E731
    warn = lambda m: res.append(("WARN", m))  # noqa: E731
    fail = lambda m: res.append(("FAIL", m))  # noqa: E731
    for u in row["_unresolved"]:
        fail(f"unresolved cell {u} (key missing from {cfg.get('_inputs_file', 'the inputs file')})")
    bam, fasta = row["bam"], row["fasta"]
    if not bam or not Path(bam).exists():
        fail(f"bam missing: {bam}")
        return res
    idx = [p for p in (bam + ".bai", bam + ".csi", bam[:-4] + ".bai") if Path(p).exists()]
    (ok if idx else fail)(f"bam {bam} ({Path(bam).stat().st_size / 1e9:.1f} GB, {bam_records(cfg, bam):,} mapped records); "
                          f"index {idx[0] if idx else 'MISSING'}")
    sq = _bam_sq(bam)
    if not fasta or not Path(fasta).exists():
        fail(f"fasta missing: {fasta}")
    else:
        fai = _fai(fasta)
        if fai is None:
            fail(f"fasta {fasta}: no .fai")
        else:
            m = _cmp_names(sq, fai, "BAM @SQ", "FASTA .fai")
            (warn if m else ok)(f"BAM @SQ ({len(sq)}) vs {Path(fasta).name}.fai ({len(fai)}): "
                                + ("; ".join(m) if m else "identical names and lengths"))
    for col in ("annotation_gff", "annotation_gtf"):
        p = row[col]
        if not p:
            (warn if col == "annotation_gff" else ok)(f"{col}: none")
            continue
        if not Path(p).exists():
            fail(f"{col} missing: {p}")
            continue
        a = annotation_seqids(cfg, p)
        names = a["seqids"]
        m = _cmp_names([n for n, _ in sq], names, "BAM @SQ", Path(p).name)
        lm = [n for n, L in sq if n in a["lengths"] and a["lengths"][n] != L]
        if lm:
            m.append(f"{len(lm)} ##sequence-region length(s) differ from the BAM: {', '.join(lm[:6])}")
        if col == "annotation_gff" and p.endswith(".gz"):
            m.append("gzip-compressed: the flag stage (missing_copy_flag --gff) reads plain text only")
        (warn if m else ok)(f"{col} {p}: {len(names)} seqids; " + ("; ".join(m) if m else "all match the BAM @SQ"))
    mmi = row["splice_mmi"]
    if not mmi:
        warn("splice_mmi: none (the flag stage needs one; `make.py runs --stage index` builds it)")
    elif not Path(mmi).exists():
        warn(f"splice_mmi {mmi}: not built yet (`make.py runs --sample {row['id']} --stage index`)")
    else:
        h = mmi_header(mmi)
        m = _cmp_names(sq, h["seqs"], "BAM @SQ", Path(mmi).name)
        preset = "splice (k15 w5)" if (h["k"], h["w"]) == (15, 5) else f"NOT a splice preset (k{h['k']} w{h['w']})"
        (warn if m or "NOT" in preset else ok)(f"splice_mmi {mmi} ({Path(mmi).stat().st_size / 1e9:.1f} GB): {preset}; "
                                               + ("; ".join(m) if m else "sequences match the BAM @SQ"))
    for spec in (row["flag_confirm"] or "").split(","):
        if spec.strip():
            name, _, p = spec.strip().partition("=")
            (ok if Path(p).exists() else fail)(f"flag_confirm {name}: {p}{'' if Path(p).exists() else ' MISSING'}")
    for tool, col in {**BASELINES, **GUIDED_BASELINES}.items():
        p = row[col]
        if not p:
            continue
        mode = "annotation-guided" if tool in GUIDED_BASELINES else "annotation-free"
        if not Path(p).exists():
            fail(f"{tool} baseline ({mode}) missing: {p}")
            continue
        names = annotation_seqids(cfg, p)["seqids"]
        extra = [n for n in names if n not in dict(sq)]
        (warn if extra else ok)(f"{tool} baseline ({mode}) {p}: {len(names)} seqids"
                                + (f"; {len(extra)} absent from the BAM: {', '.join(extra[:6])}" if extra
                                   else ", all BAM contigs"))
    return res


# ---------------------------------------------------------------- run cache
STAGE_ORDER = ["assemble", "assemble_primary", "families", "families_primary", "catalog", "assign", "index", "flag"]
STAGES = {
    "assemble": {"driver": "assemble", "suffix": "", "extra": [], "bins": ["copy_assign", "as_table"], "mm2": False,
                 "products": {"gtf": ".gtf", "molecules": ".molecules.tsv"}, "needs": []},
    "assemble_primary": {"driver": "assemble", "suffix": ".primary", "extra": ["--no-seed-secondaries"],
                         "bins": ["copy_assign"], "mm2": False, "products": {"gtf": ".gtf"}, "needs": []},
    # THE default de novo family definition (user decision 2026-09-25): the families AND their copy table (one copy
    # per member locus = its representative transcript; the gw_family_catalog copies contract), which copy
    # assignment consumes (figures/_o2.py); `catalog` below is the legacy copy catalog
    "families": {"driver": "families", "suffix": "", "extra": [], "bins": ["mcl_families"], "mm2": True,
                 "products": {"clusters": ".fam.clusters.tsv", "loci": ".fam.loci.tsv", "params": ".fam.params.tsv",
                              "copies": ".fam.copies.tsv", "copies_fa": ".fam.copies.fa"},
                 "needs": ["assemble"]},
    # the same de novo family run on the primaries-only assembly (fig. 6d compares the two seeding configurations)
    "families_primary": {"driver": "families", "suffix": ".primary", "extra": [], "bins": ["mcl_families"], "mm2": True,
                         "products": {"clusters": ".fam.clusters.tsv", "loci": ".fam.loci.tsv",
                                      "params": ".fam.params.tsv"},
                         "needs": ["assemble_primary"]},
    # LEGACY copy catalog (gw_family_catalog; kept, not the default definition; _o2 reads it only when asked:
    # cfg o2_copy_table=catalog). --piecewise: the representatives are built one contig per call (cached in PREFIX.cache/reps/), then merged
    # and the catalog continues as one run — the same products (cmp-checked, rust agent 2026-09-25). The per-call
    # budget is added by `ensure` (CATALOG_BUDGET_S), not here, so changing it does not invalidate the stamp.
    "catalog": {"driver": "catalog", "suffix": "", "extra": ["--piecewise"], "bins": ["gw_family_catalog"], "mm2": True,
                "products": {"copies": ".cat.copies.tsv", "copies_fa": ".cat.copies.fa", "families": ".cat.families.tsv",
                             "pairs": ".cat.pairs.tsv"}, "needs": []},
    # the driver's assign reads the families' copy table (<id>.fam.copies.*) since 2026-10-02 (the legacy catalog only with
    # --legacy-catalog, not passed here); the opt-in candidates stage (ruling R14) is not a run-cache stage
    "assign": {"driver": "assign", "suffix": "", "extra": [], "bins": ["copy_assign"], "mm2": True,
               "products": {"assignments": ".assign.assignments.tsv"}, "needs": ["families"]},
    "index": {"driver": None, "bins": [], "mm2": True, "needs": []},
    "flag": {"driver": "flag", "suffix": "", "extra": [], "bins": ["missing_copy_flag"], "mm2": True,
             "products": {"scan": ".flag_scan.scan.tsv", "consensus": ".flag_scan.consensus.fa",
                          "calls": ".flag.missing_copy.tsv"}, "needs": ["index"]},
}
# processes that mean "a heavy run is already going" (the machine rule: one heavy run at a time)
HEAVY = {"copy_assign", "as_table", "gw_family_catalog", "mcl_families", "missing_copy_flag", "minimap2", "gffcompare",
         "sqanti3_qc.py", "sqanti3_filter.py"}


def run_dir(cfg: dict, key: str) -> Path:
    d = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "runs" / resolve(cfg, key)
    d.mkdir(parents=True, exist_ok=True)
    return d


def prefix(cfg: dict, key: str, stage: str) -> Path:
    sid = resolve(cfg, key)
    return run_dir(cfg, sid) / (sid + STAGES[stage].get("suffix", ""))


def products(cfg: dict, key: str, stage: str) -> dict:
    if stage == "index":
        mmi = get(cfg, key)["splice_mmi"]
        return {"mmi": Path(mmi) if mmi else run_dir(cfg, key) / f"{resolve(cfg, key)}.splice.mmi"}
    p = str(prefix(cfg, key, stage))
    return {n: Path(p + ext) for n, ext in STAGES[stage]["products"].items()}


def product(cfg: dict, key: str, stage: str, name: str | None = None) -> Path:
    """Path of a genome-wide product (first product of the stage when `name` is None). Runs nothing."""
    _migrate_legacy(cfg, key)
    prods = products(cfg, key, stage)
    return prods[name] if name else next(iter(prods.values()))


_SHA_MEMO: dict = {}


def _sha1(path) -> str:
    p = Path(path)
    fp = _fp(p)
    if fp not in _SHA_MEMO:
        h = hashlib.sha1()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _SHA_MEMO[fp] = h.hexdigest()[:16]
    return _SHA_MEMO[fp]


_MM2: list = []
ALL_VS_ALL_STAGES = ("families", "families_primary", "catalog")   # one genome-wide minimap2 all-vs-all each
DEFAULT_WRAPPER = figlib.REPO / "tools" / "mm2_shard.sh"


def shard_wrapper(cfg: dict) -> str | None:
    """The all-vs-all shard wrapper the families/catalog stages run through (inputs key `runs_mm2_wrapper`; default
    tools/mm2_shard.sh; '-' / 'none' = plain minimap2), or None."""
    v = str(cfg.get("runs_mm2_wrapper", DEFAULT_WRAPPER)).strip()
    return None if v.lower() in ("", "-", "none", "off", "0") else v


def _is_wrapper(path: str | None, cfg: dict | None = None) -> bool:
    """True when `path` is the shard wrapper (proven output-neutral), by resolved path or file name."""
    if not path:
        return False
    cands = {str(DEFAULT_WRAPPER.resolve())}
    w = shard_wrapper(cfg or {}) if cfg is not None else None
    if w:
        cands.add(str(Path(w).resolve()))
    try:
        rp = str(Path(shutil.which(path) or path).resolve())
    except OSError:
        rp = path
    return rp in cands or Path(path).name == "mm2_shard.sh"


def _real_minimap2() -> str:
    """The minimap2 that actually aligns: RUSTLE_MINIMAP2 unless it is the shard wrapper, whose real minimap2 is
    MM2_SHARD_MINIMAP2 (default `minimap2` on PATH)."""
    exe = os.environ.get("RUSTLE_MINIMAP2", "minimap2")
    if _is_wrapper(exe):
        exe = os.environ.get("MM2_SHARD_MINIMAP2", "minimap2")
    return exe


def _minimap2_version(cfg) -> str:
    """'<path of the real minimap2> <version>' (the shard wrapper is transparent: see BOUNDED CALLS)."""
    if not _MM2:
        exe = _real_minimap2()
        try:
            v = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
        except OSError:
            v = "absent"
        _MM2.append(f"{shutil.which(exe) or exe} {v}")
    return _MM2[0]


def _normalize_key(k: dict) -> dict:
    """A stamp key with the shard wrapper removed (stamps written before 2026-09-25 14:30 recorded a caller's
    RUSTLE_MINIMAP2=tools/mm2_shard.sh in `env` and as the minimap2 path): the same products either way."""
    k = json.loads(json.dumps(k))
    if "driver_code" in k:   # a whole-script hash of a stamp written before driver_stage_code
        k["driver_code"] = _legacy_driver_code(k["driver_code"], k.get("driver_stage"))
    env = k.get("env") or {}
    if _is_wrapper(env.get("RUSTLE_MINIMAP2")):
        env.pop("RUSTLE_MINIMAP2")
    mm = k.get("minimap2")
    if isinstance(mm, str):
        path, sep, ver = mm.partition(" ")
        if _is_wrapper(path):
            real = _real_minimap2()
            k["minimap2"] = f"{shutil.which(real) or real}{sep}{ver}"
    return k


def _stage_inputs(cfg, row, stage) -> tuple[dict, list, str | None]:
    """(input files {label: path}, extra driver args, blocking reason or None)."""
    inputs = {"bam": row["bam"], "fasta": row["fasta"]}
    extra = list(STAGES[stage].get("extra", []))
    if stage == "index":
        return {"fasta": row["fasta"]}, [], None
    if stage == "flag":
        mmi = products(cfg, row["id"], "index")["mmi"]
        inputs["splice_mmi"] = str(mmi)
        extra += ["--index", str(mmi)]
        if not row["annotation_gff"]:
            return inputs, extra, "no annotation_gff (the driver would scan the de novo GTF instead; not wired here)"
        if row["annotation_gff"]:
            if row["annotation_gff"].endswith(".gz"):
                return inputs, extra, "annotation_gff is gzip-compressed; missing_copy_flag reads plain text"
            inputs["annotation_gff"] = row["annotation_gff"]
            extra += ["--gff", row["annotation_gff"]]
        for spec in (row["flag_confirm"] or "").split(","):
            if spec.strip():
                name, _, p = spec.strip().partition("=")
                inputs[f"confirm_{name}"] = p
                extra += ["--confirm", f"{name}={p}"]
    for label, p in inputs.items():
        if label == "splice_mmi":
            continue
        if not p or not Path(p).exists():
            return inputs, extra, f"{label} missing ({p})"
    return inputs, extra, None


_STAGE_FN = re.compile(r"^stage_(\w+)\(\)\s*\{")
_VAR_REF = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)")
_ASSIGN = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\+?=")


def driver_stage_code(path, stage: str) -> str:
    """sha1 (16 hex) of the driver CODE one stage runs (2026-09-25; replaces the whole-script hash as the stamp's
    `driver_code`, so an edit of one stage no longer re-runs every stage of every sample). Non-comment, non-blank
    lines, as assembly._code_hash, WITHOUT: the bodies of the other `stage_*()` functions, the final stage dispatch
    (`case "$STAGE" in` .. `esac` at column 0), and, of the argument parser (`while [ $# -gt 0 ]; do` .. `done`) and
    of the lines made only of assignments (the defaults), every arm / assignment of a variable that the kept code
    never references (a new flag of another stage leaves this stage's hash alone; a changed default it reads, e.g.
    SEED_SEC for assemble, does not). Shared code (helpers, the checks, the cache setup) stays in every stage's hash."""
    lines = [l.rstrip("\n") for l in open(path)]
    kept, defaults, parser = [], [], []   # (index, text); defaults: (index, [(var, token)]); parser: (index, [(vars, arm)])
    fn, in_dispatch, in_parser = None, False, False
    for i, raw in enumerate(lines):
        t = raw.strip()
        if not t or t.startswith("#"):
            continue
        m = _STAGE_FN.match(raw)
        if m:
            fn = m.group(1)
        if fn is not None:
            if fn == stage:
                kept.append((i, t))
            if raw == "}":
                fn = None
            continue
        if raw.startswith('case "$STAGE" in'):
            in_dispatch = True
        if in_dispatch:
            in_dispatch = raw != "esac"
            continue
        if raw.startswith("while [ $# -gt 0 ]"):
            in_parser = True
            continue
        if in_parser:
            if raw == "done":
                in_parser = False
                continue
            arms = []
            for arm in t.split(";;"):
                arm = arm.strip()
                head, sep, body = arm.partition(")")
                if not sep:
                    continue
                vs = {a.group(1) for tok in body.split(";") if (a := _ASSIGN.match(tok.strip()))}
                if vs:
                    arms.append((vs, arm))
            parser.append((i, arms))
            continue
        toks = [x.strip() for x in t.split(";") if x.strip()]
        if toks and all(_ASSIGN.match(x) for x in toks):
            defaults.append((i, [(_ASSIGN.match(x).group(1), x) for x in toks]))
        else:
            kept.append((i, t))
    refs = {v for _, t in kept for v in _VAR_REF.findall(t)}
    while True:   # defaults a kept default reads are kept too
        more = {v for _, toks in defaults for var, x in toks if var in refs for v in _VAR_REF.findall(x)} - refs
        if not more:
            break
        refs |= more
    out = kept + [(i, "; ".join(x for var, x in toks if var in refs)) for i, toks in defaults] \
        + [(i, ";; ".join(arm for vs, arm in arms if vs & refs)) for i, arms in parser]
    h = hashlib.sha1()
    for _, t in sorted(out):
        if t:
            h.update(t.encode() + b"\n")
    return h.hexdigest()[:16]


# Stamps written before driver_stage_code (2026-09-25) carry the WHOLE-script hash (assembly._code_hash). For each
# such hash, the per-stage hash of that same script (computed from a saved copy of it), so their stages stay fresh
# while their own code is unchanged. 53ae2f7d79c5dcaa = tools/rustle_pipeline.sh as of 2026-09-25 16:15 (every
# run-cache stamp at that time).
LEGACY_DRIVER_CODE = {
    "53ae2f7d79c5dcaa": {"assemble": "e52f85367deaaa00", "families": "5a9f8f561448466c", "catalog": "2cabf97e6c987bbe",
                         "assign": "52e1874fa465d165", "flag": "ad9d6ed6c3756c14"},
}


def _legacy_driver_code(code: str | None, driver_stage: str | None) -> str | None:
    """A stamp's `driver_code` in the per-stage scheme (a whole-script hash of LEGACY_DRIVER_CODE -> that script's
    hash for `driver_stage`); anything else unchanged."""
    return LEGACY_DRIVER_CODE.get(code or "", {}).get(driver_stage or "", code)


def _key(cfg, row, stage, inputs, extra, upstream: dict) -> dict:
    spec = STAGES[stage]
    k = {"stage": stage, "driver_stage": spec["driver"], "args": extra,
         "inputs": {lab: _fp(p) for lab, p in sorted(inputs.items())},
         "bins": {b: _sha1(Path(cfg["bin"]) / b) for b in spec["bins"]},
         "env": {e: v for e, v in sorted(os.environ.items()) if e.startswith("RUSTLE_") and e != "RUSTLE_CACHE_DIR"
                 and not (e == "RUSTLE_MINIMAP2" and _is_wrapper(v, cfg))},
         "upstream": upstream}
    if spec["driver"]:
        k["driver_code"] = driver_stage_code(cfg["driver"], spec["driver"])
    if spec["mm2"]:
        k["minimap2"] = _minimap2_version(cfg)
    return k


def _digest(k: dict) -> str:
    return hashlib.sha1(json.dumps(k, sort_keys=True).encode()).hexdigest()[:16]


def _stamp_path(cfg, key, stage) -> Path:
    return run_dir(cfg, key) / f"{stage}.key"


def _read_stamp(cfg, key, stage) -> dict | None:
    p = _stamp_path(cfg, key, stage)
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def _write_stamp(cfg, key, stage, k: dict, adopted: bool, extra: dict | None = None):
    """Write `<stage>.key`. `extra` (not part of the digest) records how the products were made or proven: the
    all-vs-all wrapper of the run (`run_env`) or a re-stamp's proof (`restamped`)."""
    p = _stamp_path(cfg, key, stage)
    tmp = p.with_suffix(".key.tmp")
    tmp.write_text(json.dumps({"digest": _digest(k), "adopted": adopted,
                               "written": _dt.datetime.now().isoformat(timespec="seconds"), "key": k,
                               **(extra or {})},
                              indent=1, sort_keys=True) + "\n")
    tmp.replace(p)
    if STAGES[stage]["driver"]:  # the per-product driver-code stamp assembly.py and _o1.py have always read
        Path(str(prefix(cfg, key, stage)) + ".driver_code").write_text(k["driver_code"] + "\n")


def _plan(cfg, key, stage, _memo=None) -> dict:
    """{'state', 'reason', 'key', 'digest', 'products'} for one stage, propagating upstream re-runs."""
    memo = {} if _memo is None else _memo
    sid = resolve(cfg, key)
    if (sid, stage) in memo:
        return memo[(sid, stage)]
    _migrate_legacy(cfg, sid)
    row = get(cfg, sid)
    prods = products(cfg, sid, stage)
    inputs, extra, blocked = _stage_inputs(cfg, row, stage)
    ups = {u: _plan(cfg, sid, u, memo) for u in STAGES[stage]["needs"]}
    out = {"products": prods, "stage": stage, "sample": sid}
    if blocked:
        out.update(state="blocked", reason=blocked, key=None, digest=None)
    elif stage == "index":
        mmi = prods["mmi"]
        k = _key(cfg, row, stage, inputs, extra, {})
        out.update(key=k, digest=_fp(mmi) if mmi.exists() else _digest(k))
        out.update(state="fresh", reason="index present") if mmi.exists() else \
            out.update(state="run", reason="splice index absent")
    else:
        bad_up = [u for u, p in ups.items() if p["state"] == "blocked"]
        if bad_up:
            out.update(state="blocked", reason=f"upstream {', '.join(bad_up)} blocked", key=None, digest=None)
        else:
            k = _key(cfg, row, stage, inputs, extra, {u: p["digest"] for u, p in ups.items()})
            d = _digest(k)
            out.update(key=k, digest=d)
            stamp = _read_stamp(cfg, sid, stage)
            have = all(p.exists() for p in prods.values())
            pending_up = [u for u, p in ups.items() if p["state"] in ("run", "stale")]
            old = _normalize_key(stamp.get("key", {})) if stamp else {}
            if not have:
                out.update(state="run", reason="products absent" + (f" (after {', '.join(pending_up)})" if pending_up else ""))
            elif stamp and (stamp.get("digest") == d or _digest(old) == d):
                out.update(state="fresh", reason=f"stamp {d}" + (" (adopted)" if stamp.get("adopted") else "")
                           + (" (re-stamped from a proof)" if stamp.get("restamped") else ""))
            elif stamp:
                changed = sorted(f for f in set(k) | set(old) if k.get(f) != old.get(f))
                out.update(state="stale", reason=f"key changed: {', '.join(changed)}", changed=changed,
                           old_digest=stamp.get("digest"))
            elif pending_up:
                out.update(state="stale", reason=f"no stamp and upstream {', '.join(pending_up)} will re-run")
            else:
                srcs = [p for p in inputs.values() if p] + [Path(cfg["bin"]) / b for b in STAGES[stage]["bins"]] \
                    + [q for u in ups.values() for q in u["products"].values()]
                old_code = Path(str(prefix(cfg, sid, stage)) + ".driver_code")
                code_ok = (not old_code.exists()
                           or _legacy_driver_code(old_code.read_text().strip(), STAGES[stage]["driver"]) == k.get("driver_code"))
                if all(figlib.fresh(p, *srcs) for p in prods.values()) and code_ok:
                    out.update(state="adopt", reason="products newer than every input and binary; no stamp yet")
                else:
                    out.update(state="stale", reason="no stamp; products older than an input/binary"
                                                     + ("" if code_ok else " or made by other driver code"))
    memo[(sid, stage)] = out
    return out


def status(cfg: dict, key: str, stage: str) -> tuple[str, str]:
    p = _plan(cfg, key, stage)
    return p["state"], p["reason"]


def _busy() -> list:
    me = os.getpid()
    out = []
    for d in Path("/proc").iterdir():
        if not d.name.isdigit() or int(d.name) == me:
            continue
        try:
            comm = (d / "comm").read_text().strip()
            if comm in HEAVY or comm[:15] in {h[:15] for h in HEAVY}:
                out.append(f"{d.name}:{comm}")
        except OSError:
            pass
    return out


def _lock(cfg, key, stage):
    lock = run_dir(cfg, key) / f"{stage}.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            pid = int(lock.read_text().strip() or 0)
            os.kill(pid, 0)
            raise RuntimeError(f"{lock}: stage already running (pid {pid})")
        except (ProcessLookupError, ValueError):
            lock.unlink(missing_ok=True)
            return _lock(cfg, key, stage)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    return lock


TIME_FMT = "wall_s\t%e\npeak_rss_kb\t%M\nuser_s\t%U\nsys_s\t%S\nexit\t%x"
# catalog --piecewise: one call does at most this much (seconds; the first contig piece of a call always runs) and
# exits 75 while pieces or the merge remain. Override with `catalog_budget_s` in the inputs file.
CATALOG_BUDGET_S = "420"
# ... and contigs of more than this many BAM records are cut into read-free pieces of about this size (still exact;
# human A119b chr1, 5.5 M records, exceeded 10 min as one piece). Override with `catalog_piece_records`.
CATALOG_PIECE_RECORDS = "1000000"
PIECES_PENDING_EXIT = 75


class StagePending(RuntimeError):
    """A bounded stage call (catalog --piecewise, or an all-vs-all through the shard wrapper) made progress but has
    more to do: run the same stage again."""


# ---- the all-vs-all shard wrapper (families, catalog): child environment and the outcome of a bounded call
MM2_SHARD_ROOT_DEFAULT = "/mnt/linuxdisk/tmp/mm2_shard_cache"   # = tools/mm2_shard.sh's MM2_SHARD_DIR default
# wrapper.log messages (tools/mm2_shard.sh `log` calls), classified; the wrapper exits 75/76 on the budget ones
_WRAP_BUDGET = ("budget: shard", "hit the deadline", "budget exhausted before the index build", "is gone: stopping",
                "stopping the parent")
_WRAP_PROGRESS = ("index: built", "layout:")
_WRAP_FAIL = ("failed", "REFUSED", "verification failed", "no completion line", "differs from its .ok")


def _aligner_env(cfg: dict, stage: str, t0: float) -> tuple[dict, dict]:
    """(environment overlay for the stage's child, provenance for the stamp). The families and catalog stages run
    their all-vs-all through the shard wrapper unless the caller's RUSTLE_MINIMAP2 names another aligner; the
    wrapper's per-call budget is set here unless the caller's environment sets it."""
    if stage not in ALL_VS_ALL_STAGES:
        return {}, {}
    env: dict = {}
    use = os.environ.get("RUSTLE_MINIMAP2")
    if not use:
        use = shard_wrapper(cfg)
        if not use:
            return {}, {}
        env["RUSTLE_MINIMAP2"] = use   # the same string every call: it is part of the Rust caches' keys
    if not _is_wrapper(use, cfg):
        return env, {"aligner": use}
    if "MM2_SHARD_BUDGET_S" not in os.environ:
        env["MM2_SHARD_BUDGET_S"] = str(int(float(cfg.get("runs_mm2_budget_s", "480"))))
    if "MM2_SHARD_DEADLINE" not in os.environ:
        env["MM2_SHARD_DEADLINE"] = str(int(t0 + float(cfg.get("runs_call_budget_s", "560"))))
    wpath = Path(shutil.which(use) or use)
    shard_env = {k: v for k, v in sorted({**os.environ, **env}.items()) if k.startswith("MM2_SHARD_")
                 and k not in ("MM2_SHARD_DEADLINE",)}
    return env, {"aligner": use, "aligner_sha1": _sha1(wpath) if wpath.exists() else "absent",
                 "real_minimap2": _minimap2_version(cfg), "shard_env": shard_env,
                 "note": "tools/mm2_shard.sh: shards concatenate to the single run byte for byte (cmp, chr16); "
                         "not part of the stage key"}


def _wrapper_activity(t0: float) -> list[tuple[Path, list[str]]]:
    """[(shard dir, wrapper.log messages written since t0)] for every wrapper key dir touched since t0, oldest
    first (the wrapper's own log, `%F %T message`, in $MM2_SHARD_DIR/<md5>/a<sha>/q<md5>/t<threads>.<spec>/)."""
    root = Path(os.environ.get("MM2_SHARD_DIR") or MM2_SHARD_ROOT_DEFAULT)
    out = []
    for wl in root.glob("*/a*/q*/t*/wrapper.log"):
        try:
            if wl.stat().st_mtime < t0 - 2:
                continue
            msgs = []
            for line in open(wl, errors="replace"):
                try:
                    t = _dt.datetime.strptime(line[:19], "%Y-%m-%d %H:%M:%S").timestamp()
                except ValueError:
                    continue
                if t >= int(t0):   # the log stamps whole seconds: the call's own second onward
                    msgs.append(line[20:].rstrip("\n"))
            if msgs:
                out.append((wl.stat().st_mtime, wl.parent, msgs))
        except OSError:
            continue
    return [(d, m) for _, d, m in sorted(out, key=lambda x: x[0])]


def all_vs_all_outcome(t0: float) -> tuple[str, str]:
    """How the shard wrapper's work in a call that started at t0 ended: ('pending' | 'stuck' | 'failed' | 'complete'
    | 'none' | 'unknown', last message). 'pending' = it stopped on its budget after making progress (a shard mapped,
    the index built or the query split) -> call again; 'stuck' = it stopped on its budget without progress (one step
    does not fit a call: exit 76 territory) -> do not loop."""
    acts = _wrapper_activity(t0)
    if not acts:
        return "none", "no shard-wrapper activity in this call"
    msgs = [m for _, ms in acts for m in ms]
    last = acts[-1][1][-1]
    progress = any((m.startswith("shard ") and " done: " in m) or m.startswith(_WRAP_PROGRESS) for m in msgs)
    if any(x in last for x in _WRAP_FAIL):
        return "failed", f"{last} ({acts[-1][0]})"
    if last.startswith("complete:"):
        return "complete", f"{last} ({acts[-1][0]})"
    if any(x in last for x in _WRAP_BUDGET):
        return ("pending" if progress else "stuck"), f"{last} ({acts[-1][0]})"
    return "unknown", f"{last} ({acts[-1][0]})"


def _cleanup_refine_tmp(t0: float) -> list[str]:
    """Remove the catalog's `rustle_refine_<pid>_*` temporary FASTAs written since t0 by a process that no longer
    exists (the shard wrapper stops the catalog with SIGTERM on its budget, and Drop never runs)."""
    d = Path(os.environ.get("TMPDIR") or "/mnt/linuxdisk/tmp")   # = figlib.run's default TMPDIR
    gone = []
    for f in d.glob("rustle_refine_*"):
        m = re.match(r"rustle_refine_(\d+)_", f.name)
        try:
            if not m or not f.is_file() or Path(f"/proc/{m.group(1)}").exists() or f.stat().st_mtime < t0:
                continue
            f.unlink()
            gone.append(f.name)
        except OSError:
            continue
    return gone


def ensure(cfg: dict, key: str, stage: str, *, force: bool = False, ignore_busy: bool = False) -> dict:
    """Run one genome-wide stage if it is not fresh (upstream stages first). HEAVY: foreground, one process."""
    sid = resolve(cfg, key)
    for u in STAGES[stage]["needs"]:
        ensure(cfg, sid, u, ignore_busy=ignore_busy)
    plan = _plan(cfg, sid, stage)
    if plan["state"] == "blocked":
        raise RuntimeError(f"{sid} {stage}: blocked — {plan['reason']}")
    if plan["state"] == "fresh" and not force:
        return plan["products"]
    if plan["state"] == "adopt" and not force:
        _write_stamp(cfg, sid, stage, plan["key"], adopted=True)
        print(f"[runs] {sid} {stage}: adopted existing products ({plan['reason']})", file=sys.stderr)
        return plan["products"]
    busy = _busy()
    if busy and not ignore_busy:
        raise RuntimeError(f"{sid} {stage}: another heavy process is running ({', '.join(busy)}); one heavy run at a "
                           "time (pass --ignore-busy to override)")
    row = get(cfg, sid)
    d = run_dir(cfg, sid)
    inputs, extra, _ = _stage_inputs(cfg, row, stage)
    if stage == "index":
        mmi = plan["products"]["mmi"]
        threads = str(min(3, int(cfg.get("threads", "4"))))  # minimap2 indexes with at most 3 threads
        cmd = [os.environ.get("RUSTLE_MINIMAP2", "minimap2"), "-x", "splice", "-t", threads, "-d", str(mmi) + ".tmp",
               row["fasta"]]
    else:
        cmd = ["bash", cfg["driver"], STAGES[stage]["driver"], "--bam", row["bam"], "--fasta", row["fasta"],
               "--out", str(prefix(cfg, sid, stage)), "--bin", cfg["bin"], "--threads", cfg.get("threads", "4")] + extra
        if stage == "catalog":  # bounded call; not part of the stage key (see STAGES["catalog"])
            cmd += ["--budget-s", str(cfg.get("catalog_budget_s", CATALOG_BUDGET_S)),
                    "--piece-records", str(cfg.get("catalog_piece_records", CATALOG_PIECE_RECORDS))]
    timef = d / f"{stage}.time"
    lock = _lock(cfg, sid, stage)
    try:
        t0 = time.time()
        env, run_env = _aligner_env(cfg, stage, t0)
        with open(timef, "w") as fh:
            fh.write(f"sample\t{sid}\nstage\t{stage}\nstarted\t{_dt.datetime.now().isoformat(timespec='seconds')}\n"
                     f"command\t{' '.join(cmd)}\n"
                     + (f"env\t{' '.join(f'{k}={v}' for k, v in env.items())}\n" if env else ""))
        print(f"[runs] {sid} {stage}: {' '.join(f'{k}={v}' for k, v in env.items())}{' ' if env else ''}"
              f"{' '.join(cmd)}", file=sys.stderr)
        rc = figlib.run(["/usr/bin/time", "-a", "-o", str(timef), "-f", TIME_FMT] + cmd,
                        log=d / f"{stage}.driver.log", env=env, check=False)
        with open(timef, "a") as fh:
            fh.write(f"finished\t{_dt.datetime.now().isoformat(timespec='seconds')}\n")
        if rc == PIECES_PENDING_EXIT and stage == "catalog":
            raise StagePending(f"{sid} {stage}: {time.time() - t0:.0f} s of bounded work done, more pieces or the merge "
                               f"remain — run the same stage again (progress: {d / (stage + '.driver.log')}, "
                               f"{prefix(cfg, sid, stage)}.catalog.log)")
        if rc != 0 and run_env.get("aligner_sha1"):   # the all-vs-all ran through the shard wrapper
            outcome, msg = all_vs_all_outcome(t0)
            gone = _cleanup_refine_tmp(t0)
            with open(timef, "a") as fh:
                fh.write(f"all_vs_all\t{outcome}: {msg}\n" + (f"removed_tmp\t{','.join(gone)}\n" if gone else ""))
            if outcome == "pending":
                raise StagePending(f"{sid} {stage}: {time.time() - t0:.0f} s; the all-vs-all shard wrapper stopped on its "
                                   f"budget after making progress ({msg}) — run the same stage again")
            if outcome == "stuck":
                raise RuntimeError(f"{sid} {stage}: the all-vs-all made NO progress in a whole call ({msg}): one step (the "
                                   "index or a single shard) does not fit MM2_SHARD_BUDGET_S. Raise runs_mm2_budget_s / "
                                   "runs_call_budget_s for one call, or set MM2_SHARD_BP smaller (a new shard layout; "
                                   "finished shards of the old layout are not reused)")
        if rc != 0:
            raise RuntimeError(f"{sid} {stage} failed (exit {rc}) after {time.time() - t0:.0f} s — see "
                               f"{d / (stage + '.driver.log')} and the driver's PREFIX.*.log")
        if stage == "index":
            Path(str(mmi) + ".tmp").replace(mmi)
        missing = [f"{n} ({p})" for n, p in plan["products"].items() if not Path(p).exists()]
        if missing:   # never stamp a run that did not make its products (e.g. a binary older than the stage's products)
            raise RuntimeError(f"{sid} {stage}: the run exited 0 but did not write {', '.join(missing)} — see "
                               f"{d / (stage + '.driver.log')} (families: a {cfg['bin']}/mcl_families older than the "
                               "copy table writes no PREFIX.fam.copies.tsv; rebuild it)")
        # the key the run was started under; how the all-vs-all ran is recorded beside it, not in it
        _write_stamp(cfg, sid, stage, plan["key"], adopted=False, extra={"run_env": run_env} if run_env else None)
    finally:
        lock.unlink(missing_ok=True)
    return plan["products"]


# ---------------------------------------------------------------- re-stamp from a proof (no run)
# key fields a proof may bridge: the code that made the products changed, not what they were made from
RESTAMP_FIELDS = {"bins", "driver_code", "minimap2", "upstream"}


def _file_sha1_full(path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_proof(path) -> dict:
    """A re-stamp proof: TSV, '#' lines are comments. Rows:
         product<TAB>NAME<TAB>CACHED_PATH<TAB>FRESH_PATH   one per product of the stage (every product, exactly once):
                                                           CACHED = the run cache's product, FRESH = the same product
                                                           made by the CURRENT binaries/driver (e.g. in a scratch
                                                           PREFIX), cmp-identical to CACHED
         command<TAB>TEXT                                  how FRESH was made (optional, recorded)
       Any other first field is refused."""
    out = {"products": {}, "command": []}
    for n, line in enumerate(open(path), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "product" and len(f) == 4:
            if f[1] in out["products"]:
                raise ValueError(f"{path}:{n}: product {f[1]!r} listed twice")
            out["products"][f[1]] = (f[2], f[3])
        elif f[0] == "command" and len(f) >= 2:
            out["command"].append("\t".join(f[1:]))
        else:
            raise ValueError(f"{path}:{n}: expected 'product<TAB>NAME<TAB>CACHED<TAB>FRESH' or 'command<TAB>TEXT'")
    return out


def restamp(cfg: dict, key: str, stage: str, proof_path) -> dict:
    """Re-stamp a STALE stage under the current key WITHOUT running it, when `proof_path` lists, for every product of
    the stage, a FRESH copy made by the current code that is byte-identical to the cached product. Checked here, not
    trusted: (1) the stage is stale only through code fields (RESTAMP_FIELDS: binaries, driver code, minimap2 build,
    and the upstream key when every upstream stage is fresh) — never inputs, arguments or RUSTLE_* settings; (2) every
    product is listed once, CACHED is the stage's product path; (3) FRESH exists, is another file, and is newer than
    every binary of the stage and the driver (made after the current code existed); (4) CACHED and FRESH have equal
    size and sha1 (computed now). The proof (its path, sha1 and text, the per-product sha1, the old digest and the
    changed fields) is recorded in the new stamp under `restamped`. Returns that record."""
    sid = resolve(cfg, key)
    if stage == "index":
        raise ValueError("the index stage is keyed by the index file itself; nothing to re-stamp")
    plan = _plan(cfg, sid, stage)
    if plan["state"] == "fresh":
        raise ValueError(f"{sid} {stage} is already fresh ({plan['reason']}); nothing to re-stamp")
    if plan["state"] != "stale" or "changed" not in plan:
        raise ValueError(f"{sid} {stage} is {plan['state']} ({plan['reason']}): only a stale stage with a stamp can be "
                         "re-stamped")
    bad = [f for f in plan["changed"] if f not in RESTAMP_FIELDS]
    if bad:
        raise ValueError(f"{sid} {stage}: key fields {bad} changed; a proof can bridge only {sorted(RESTAMP_FIELDS)} "
                         "(changed inputs, arguments or settings mean the products may legitimately differ: re-run)")
    if "upstream" in plan["changed"]:
        ups = {u: _plan(cfg, sid, u)["state"] for u in STAGES[stage]["needs"]}
        if any(s != "fresh" for s in ups.values()):
            raise ValueError(f"{sid} {stage}: upstream key changed and upstream stages are not all fresh ({ups}); "
                             "re-stamp or re-run them first")
    proof = read_proof(proof_path)
    prods = plan["products"]
    if set(proof["products"]) != set(prods):
        raise ValueError(f"proof lists products {sorted(proof['products'])}; {stage} has {sorted(prods)} (every one, "
                         "exactly once)")
    code_files = [Path(cfg["bin"]) / b for b in STAGES[stage]["bins"]] + ([Path(cfg["driver"])]
                                                                          if STAGES[stage]["driver"] else [])
    code_mtime = max((p.stat().st_mtime for p in code_files if p.exists()), default=0.0)
    record = {"proof": str(Path(proof_path).resolve()), "proof_sha1": _file_sha1_full(proof_path),
              "proof_text": Path(proof_path).read_text(), "from_digest": plan.get("old_digest"),
              "changed": plan["changed"], "command": proof["command"], "products": {},
              "checked": _dt.datetime.now().isoformat(timespec="seconds")}
    for name, (cached, fresh_p) in sorted(proof["products"].items()):
        want = prods[name]
        if Path(cached).resolve() != want.resolve():
            raise ValueError(f"product {name}: proof names {cached}, the stage's product is {want}")
        fp = Path(fresh_p)
        if not fp.is_file():
            raise ValueError(f"product {name}: fresh copy {fp} does not exist")
        if fp.resolve() == want.resolve() or os.path.samefile(fp, want):
            raise ValueError(f"product {name}: the fresh copy IS the cached product (a proof needs a second file)")
        if fp.stat().st_mtime < code_mtime:
            raise ValueError(f"product {name}: fresh copy {fp} is older than the current binaries/driver (made by "
                             "older code?)")
        sa, sb = want.stat().st_size, fp.stat().st_size
        ha = _file_sha1_full(want)
        hb = _file_sha1_full(fp) if sa == sb else None
        if sa != sb or ha != hb:
            raise ValueError(f"product {name}: {want} and {fp} differ (size {sa} vs {sb}); not cmp-identical")
        record["products"][name] = {"cached": str(want), "fresh": str(fp.resolve()), "bytes": sa, "sha1": ha,
                                    "fresh_mtime": _dt.datetime.fromtimestamp(fp.stat().st_mtime).isoformat(
                                        timespec="seconds")}
    _write_stamp(cfg, sid, stage, plan["key"], adopted=False, extra={"restamped": record})
    with open(run_dir(cfg, sid) / f"{stage}.restamp.log", "a") as fh:
        fh.write(json.dumps({"stage": stage, "digest": plan["digest"], **record}, sort_keys=True) + "\n")
    return record


# ---------------------------------------------------------------- legacy assemblies (figures 1-3's cache)
LEGACY = [("rustle.genome", ""), ("rustle_primary.genome", ".primary")]
LEGACY_LOGS = {"rustle.assemble.driver.log": "assemble.driver.log",
               "rustle_primary.assemble.driver.log": "assemble_primary.driver.log"}


def legacy_moves(cfg, key) -> list:
    """[(old, new)] for the legacy assembly files of an alias sample that are still regular files."""
    row = get(cfg, key)
    if not row["alias"]:
        return []
    old_dir = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "assembly" / row["alias"]
    if not old_dir.is_dir():
        return []
    new_dir = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "runs" / row["id"]
    moves = []
    for f in sorted(old_dir.iterdir()):
        if f.is_symlink() or not f.is_file():
            continue
        for old_stem, new_suffix in LEGACY:
            if f.name.startswith(old_stem + "."):
                moves.append((f, new_dir / (row["id"] + new_suffix + f.name[len(old_stem):])))
        if f.name in LEGACY_LOGS:
            moves.append((f, new_dir / LEGACY_LOGS[f.name]))
    return moves


def _migrate_legacy(cfg, key, verbose=False) -> list:
    """Move the legacy assembly files (hard link, then an atomic swap of the old name for a symlink). Idempotent."""
    done = []
    for old, new in legacy_moves(cfg, key):
        if new.exists():
            if not new.samefile(old):  # both exist and differ: never clobber, leave for the analyst
                print(f"[runs] legacy {old} and {new} both exist and differ; left as is", file=sys.stderr)
            continue
        new.parent.mkdir(parents=True, exist_ok=True)
        st = old.stat()
        os.link(old, new)
        tmp = old.with_name(old.name + ".symlink.tmp")
        tmp.unlink(missing_ok=True)
        os.symlink(os.path.relpath(new, old.parent), tmp)
        os.replace(tmp, old)
        nst = new.stat()
        assert (nst.st_ino, nst.st_size, nst.st_mtime_ns) == (st.st_ino, st.st_size, st.st_mtime_ns)
        done.append((old, new, st.st_ino))
        if verbose:
            print(f"[runs] moved {old} -> {new} (inode {st.st_ino})", file=sys.stderr)
    return done


# ---------------------------------------------------------------- cost model (dry-run estimates)
# Recorded measurements the estimates rest on (each is quoted in the dry-run's `basis` column):
#   assemble          gorilla_OR6737 10.71 M records 243 s; human_A119b 68.03 M 1,159 s (driver logs, 2026-09-25, current
#                     binaries, as_table included) -> 72 s + 16.0 s per M records. Peak RSS never recorded; r1086 measured
#                     the streaming assembler at 1.93 GB on the whole human genome (primaries only) + the best-AS table.
#   assemble_primary  91 s / 464 s on the same two BAMs -> 21 s + 6.5 s per M records; ~1.9 GB (r1086).
#   families          never run genome-wide. Per contig (fig7, 2026-09-25): locus FASTA 25-163 MB -> 34-416 s, time
#                     roughly QUADRATIC in the locus FASTA (55 MB 46 s, 163 MB 416 s); < 4 GB per contig. The copy
#                     table (--emit-units) adds one load of the family contigs (~3 GB for a whole human genome) and a
#                     spliced exon sum per member: human chr16 +0.1 s / +30 MB (2026-09-25).
#   catalog           current binary: human chr16 slice 1.79 M records 294 s / 7 GB (fig6, 2026-09-25). Genome-wide only
#                     with the pre-wave-7 binary: GGO 26 min, PTR 38, PPY 20, HSA 41 (winloci_data/*_gwcat.log, 07-17),
#                     peak 25.5 GB (o1_replicate/fibro_gwcat.time) and 25.7 GB (o1_gw/ggo_gw.time).
#   assign            gorilla NC_073244.2 slice (tes44, 473 K records, 2026-09-24, pre-wave-7) 83 s -> 175 s per M records;
#                     peak never recorded per stage (whole tes44 driver run 15.7 GB, the flag stage's index).
#   index             CHM13 splice index 403 s / 19.2 GB (memory note 09-17).
#   flag              KB3781 genome-wide scan ~8 min, testis ~25 min (r1090/r1094 notes); align + verdict with three
#                     13 GB indexes 3:46 / 17.5 GB (gw22/o3/fib_align.log); tes44 whole run 15.7 GB with one index.
# genome-wide de novo locus FASTA (sum of gene spans of the assembled GTF), measured 2026-09-25 on human_A119b (1,980 MB)
# and gorilla_OR6737 (1,010 MB); a sample whose GTF is not assembled yet borrows its species' (apes: the gorilla's)
LOCUS_FASTA_MB = {"human": 1980, "gorilla": 1010, "chimpanzee": 1010, "orangutan": 1010}


def _measured(cfg, sid, stage):
    f = run_dir(cfg, sid) / f"{stage}.time"
    if not f.exists():
        return None
    kv = dict(l.rstrip("\n").split("\t", 1) for l in open(f) if "\t" in l)
    if kv.get("exit") != "0" or "wall_s" not in kv:
        return None
    return float(kv["wall_s"]), float(kv["peak_rss_kb"]) / 1e6


def _locus_fasta_mb(cfg, sid, stage: str = "families") -> tuple[float, str]:
    """Genome-wide de novo locus FASTA size: from this sample's assembled GTF when present (sum of gene spans), else
    the same species' estimate."""
    gtf = products(cfg, sid, STAGES[stage]["needs"][0])["gtf"]
    if gtf.exists():
        def compute(p):
            spans = {}
            with open(p) as fh:
                for line in fh:
                    f = line.split("\t", 9)
                    if len(f) < 9 or f[2] != "transcript":
                        continue
                    m = re.search(r'gene_id "([^"]+)"', f[8])
                    g = (f[0], m.group(1) if m else f[8])
                    s, e = int(f[3]), int(f[4])
                    a = spans.get(g)
                    spans[g] = (min(a[0], s), max(a[1], e)) if a else (s, e)
            return sum(e - s + 1 for s, e in spans.values()) / 1e6
        return _cached(cfg, "locusmb", gtf, compute), "gene spans of this sample's GTF"
    sp = get(cfg, sid)["species"]
    return LOCUS_FASTA_MB.get(sp, 1500), f"{sp} guess (GTF not assembled yet)"


def estimate(cfg, key, stage) -> dict:
    """{'wall_s': (lo, hi), 'rss_gb': (lo, hi) or None, 'basis': str, 'measured': bool}."""
    sid = resolve(cfg, key)
    m = _measured(cfg, sid, stage)
    if m:
        return {"wall_s": (m[0], m[0]), "rss_gb": (m[1], m[1]), "basis": "measured (" + stage + ".time)",
                "measured": True}
    row = get(cfg, sid)
    n = bam_records(cfg, row["bam"]) / 1e6 if row["bam"] and Path(row["bam"]).exists() else 0.0
    if stage == "assemble":
        t = 72 + 16.0 * n
        return {"wall_s": (0.8 * t, 1.3 * t), "rss_gb": (2.0, 6.0), "measured": False,
                "basis": f"72 s + 16.0 s/M x {n:.1f} M records (OR6737 243 s, A119b 1159 s); RSS unmeasured "
                         "(1.9 GB streaming + best-AS table)"}
    if stage == "assemble_primary":
        t = 21 + 6.5 * n
        return {"wall_s": (0.8 * t, 1.3 * t), "rss_gb": (1.5, 3.0), "measured": False,
                "basis": f"21 s + 6.5 s/M x {n:.1f} M (OR6737 91 s, A119b 464 s); RSS ~1.9 GB (r1086)"}
    if stage in ("families", "families_primary"):
        mb, src = _locus_fasta_mb(cfg, sid, stage)
        r = mb / 163.0
        return {"wall_s": (416 * r, 416 * r * r), "rss_gb": None, "measured": False,
                "basis": f"NEVER RUN GENOME-WIDE; locus FASTA ~{mb:.0f} MB ({src}) vs chr2 163 MB / 416 s, "
                         "linear..quadratic; RSS unknown (per contig < 4 GB; genome-wide PAF may be GBs)"}
    if stage == "catalog":
        return {"wall_s": (164 * n, 2 * 164 * n + 600), "rss_gb": (7.0, 26.0), "measured": False,
                "basis": f"164 s/M x {n:.1f} M (chr16 294 s / 7 GB, current binary); old genome-wide binary peaked "
                         "25.5-25.7 GB"}
    if stage == "assign":
        return {"wall_s": (175 * n * 0.5, 175 * n * 1.5), "rss_gb": None, "measured": False,
                "basis": f"175 s/M x {n:.1f} M (tes44 slice 83 s, pre-wave-7); RSS never recorded"}
    if stage == "index":
        if products(cfg, sid, "index")["mmi"].exists():
            return {"wall_s": (0, 0), "rss_gb": (0, 0), "basis": "index present", "measured": True}
        return {"wall_s": (350, 500), "rss_gb": (18.0, 21.0), "measured": False,
                "basis": "CHM13 splice index 403 s / 19.2 GB"}
    if stage == "flag":
        confirm = len([s for s in (row["flag_confirm"] or "").split(",") if s.strip()])
        mmi = products(cfg, sid, "index")["mmi"]
        gb = (mmi.stat().st_size / 1e9 if mmi.exists() else 13.0) + 2.5 + (1.5 if confirm else 0)
        return {"wall_s": (8 * 60 + 120, 25 * 60 + 60 * (1 + confirm)), "rss_gb": (gb - 1.5, gb + 1.0),
                "measured": False,
                "basis": f"scan 8-25 min (KB3781 / testis), align {1 + confirm} index(es) ~1 min each; RSS = index "
                         "size + ~2.5 GB (fib_align 17.5 GB with 3 indexes, tes44 15.7 GB with 1)"}
    raise KeyError(stage)


def _fmt_s(s: float) -> str:
    return f"{s / 3600:.1f} h" if s >= 5400 else f"{s / 60:.0f} min" if s >= 90 else f"{s:.0f} s"


# ---------------------------------------------------------------- CLI (make.py samples / make.py runs)
def cli_samples(cfg, verify_: bool, only: str | None = None):
    reg = registry(cfg)
    sids = [resolve(cfg, only)] if only else list(reg)
    for sid in sids:
        row = reg[sid]
        print(f"{sid}  ({row['species']}, {row['tissue']}; alias {row['alias'] or '-'})")
        for c in COLUMNS[3:]:
            if c in ("tissue",):
                continue
            print(f"    {c:15s} {row[c] or '-'}")
        if verify_:
            for level, msg in verify(cfg, sid):
                print(f"    {level:4s} {msg}")
    return 0


def queue(cfg, sids, stages) -> list:
    memo: dict = {}
    out = []
    for sid in sids:
        for st in stages:
            out.append(_plan(cfg, sid, st, memo))
    return out


def cli_runs(cfg, sample: str, stage: str, dry_run: bool, max_stages: int | None, force: bool, ignore_busy: bool):
    reg = registry(cfg)
    sids = list(reg) if sample == "all" else [resolve(cfg, sample)]
    asked = STAGE_ORDER if stage == "all" else [stage]
    for s in asked:
        if s not in STAGES:
            sys.exit(f"unknown stage {s!r}; stages: {', '.join(STAGE_ORDER)}, all")
    need = set(asked)  # the queue shows (and runs) the upstream stages a requested stage needs
    while True:
        more = {u for s in need for u in STAGES[s]["needs"]} - need
        if not more:
            break
        need |= more
    stages = [s for s in STAGE_ORDER if s in need]
    for sid in sids:
        for old, new, ino in _migrate_legacy(cfg, sid, verbose=True):
            pass
    plans = queue(cfg, sids, stages)
    todo = [p for p in plans if p["state"] in ("run", "stale") or (force and p["state"] != "blocked")]
    if dry_run:
        hdr = f"{'sample':16s} {'stage':17s} {'state':7s} {'est. wall':>17s} {'est. peak RSS':>14s}  flags / reason / basis"
        print(hdr)
        tot_lo = tot_hi = 0.0
        for p in plans:
            e = estimate(cfg, p["sample"], p["stage"])
            q = p["state"] in ("run", "stale")
            wall = "-" if not q else (f"{_fmt_s(e['wall_s'][0])}" if e["wall_s"][0] == e["wall_s"][1]
                                      else f"{_fmt_s(e['wall_s'][0])}-{_fmt_s(e['wall_s'][1])}")
            rss = "-" if not q else ("unknown" if e["rss_gb"] is None else
                                     f"{e['rss_gb'][0]:.0f}-{e['rss_gb'][1]:.0f} GB")
            flags = []
            if q:
                tot_lo += e["wall_s"][0]
                tot_hi += e["wall_s"][1]
                if e["wall_s"][1] > 600:
                    flags.append(">10min")
                if e["rss_gb"] is None or e["rss_gb"][1] > 20:
                    flags.append(">20GB?" if e["rss_gb"] is None else ">20GB")
            print(f"{p['sample']:16s} {p['stage']:17s} {p['state']:7s} {wall:>17s} {rss:>14s}  "
                  f"{'[' + ' '.join(flags) + '] ' if flags else ''}{p['reason']}"
                  + (f" | {e['basis']}" if q else ""))
        print(f"# queued {len(todo)} stage(s); estimated {_fmt_s(tot_lo)} .. {_fmt_s(tot_hi)} of serial wall time; "
              f"'adopt' = stamp written on the next real call, no work")
        return 0
    n = 0
    for p in plans:
        if p["state"] == "blocked":
            print(f"[runs] {p['sample']} {p['stage']}: blocked — {p['reason']}", file=sys.stderr)
            continue
        if p["state"] in ("fresh",) and not force:
            continue
        if max_stages is not None and n >= max_stages and p["state"] != "adopt":
            print(f"[runs] --max-stages {max_stages} reached; re-run to continue (exit {PIECES_PENDING_EXIT})",
                  file=sys.stderr)
            return PIECES_PENDING_EXIT
        # re-plan: an upstream stage run earlier in this call changes this one's key
        cur = _plan(cfg, p["sample"], p["stage"])
        if cur["state"] == "fresh" and not force:
            continue
        try:
            ensure(cfg, p["sample"], p["stage"], force=force, ignore_busy=ignore_busy)
        except StagePending as e:  # one bounded call per stage run; its downstream stages must wait
            print(f"[runs] {e} (exit {PIECES_PENDING_EXIT})", file=sys.stderr)
            return PIECES_PENDING_EXIT
        if cur["state"] != "adopt":
            n += 1
            e = _measured(cfg, p["sample"], p["stage"])
            if e:
                print(f"[runs] {p['sample']} {p['stage']}: {_fmt_s(e[0])}, peak {e[1]:.1f} GB", file=sys.stderr)
    return 0


def cli_restamp(cfg, sample: str, stage: str, proof: str) -> int:
    """make.py runs --sample S --stage ST --restamp PROOF: re-stamp one stale stage from a proof (no run)."""
    if sample == "all" or stage == "all":
        sys.exit("--restamp needs one --sample and one --stage")
    if stage not in STAGES:
        sys.exit(f"unknown stage {stage!r}; stages: {', '.join(STAGE_ORDER)}")
    try:
        rec = restamp(cfg, sample, stage, proof)
    except (ValueError, OSError) as e:
        print(f"[runs] restamp refused: {e}", file=sys.stderr)
        return 1
    sid = resolve(cfg, sample)
    print(f"[runs] {sid} {stage}: re-stamped without running (key fields bridged: {', '.join(rec['changed'])}; "
          f"{len(rec['products'])} product(s) cmp-identical: "
          + ", ".join(f"{k} {v['bytes']:,} B sha1 {v['sha1'][:12]}" for k, v in rec["products"].items())
          + f"); proof recorded in {_stamp_path(cfg, sid, stage)}", file=sys.stderr)
    return 0
