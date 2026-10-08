"""assembly — shared provisioning for the transcript-assembly figures (1-3).

`species` below is a sample key of figures/samples.tsv: a sample id (`gorilla_OR6737`, `chimp_PTR`, ...) or a legacy alias
(`gorilla` = gorilla_OR6737, `human` = human_A119b; the figure modules pass these). For one sample it produces (cached):

    ${work}/runs/<id>/<id>.gtf, <id>.primary.gtf
                                 our genome-wide assemblies, from the run cache (samples.ensure, stages `assemble` and
                                 `assemble_primary`; rustle: driver default, loci seeded with secondaries within 2% of the
                                 molecule's best alignment score; rustle_primary: --no-seed-secondaries) — the expensive
                                 step, recomputed only with `force_arms=True` or when its run-cache key changes
                                 (the old `${work}/assembly/<alias>/rustle*.genome.*` names are symlinks to these)
    ${work}/assembly/<alias or id>/eval_<contigs>/ref.gtf
                                 the annotation restricted to the evaluation contigs (the directory is keyed by
                                 the evaluation scope, so a changed scope never reuses a stale restriction)
    eval_<contigs>/<tool>.gtf    each arm's transcripts restricted to the same contigs
                                 tools: rustle, rustle_primary, stringtie, flair, isoseq (all ANNOTATION-FREE: MODE_DENOVO);
                                 and, only where the user registered them, stringtie_guided, flair_guided
                                 (ANNOTATION-GUIDED: the separate guided path; never in a panel with the others)
    eval_<contigs>/gc_<tool>.stats / .tmap / .refmap   gffcompare of <tool>.gtf against ref.gtf
                                 (`force=True` recomputes the restrictions and gffcompare, not the assemblies;
                                 gffcompare's .annotated.gtf/.loci/.tracking are deleted — nothing reads them)

Evaluation contigs (every sample, genome-wide by default): every contig the sample's annotation covers. When the
annotation covers every contig of the genome (gorilla, chimpanzee, orangutan) there is no restriction at all
(`eval_all/`); human CHM13 RefSeq leaves out chrM, so the human samples are scored on the 24 annotated contigs
(`eval_annotated/`) and chrM transcripts are left out for every method alike (they could only count as unmatched).
`eval_contigs_<sample id>` (comma list; `all` = the default) narrows one sample for development. The old key
`human_regions` (chr20-22, chosen when the human FLAIR and IsoSeq outputs were scored on the cluster) is no longer
read: genome-wide gffcompare of the 3.2 M-transcript human IsoSeq GTF takes 107 s / 5.1 GB here. Every method is
restricted to the SAME contigs, so a comparison is always like for like. The lab baselines (stringtie, flair, isoseq)
exist only for human_A119b and gorilla_OR6737 (samples.baseline; benchmark_samples()).

Annotation GTF of a sample (annotation_gtf): the registry's `annotation_gtf`, or, where only a RefSeq GFF3 exists
(chimpanzee, orangutan), a GTF derived from it with the `gff_to_gtf` binary, one contig per call, concatenated in the
GFF's contig order -- the same converter (a stand-in for `gffread -T`, every transcript type kept) that built the
human and gorilla reference GTFs (benchmark_collapse/score_vs_annotation.sh). Cached under
`${work}/assembly/<id>/annotation.gff_to_gtf.gtf` with a key sidecar (GFF fingerprint + converter sha1).

Budgeted calls (figures 1-3): `figs_budget_s` (or `<fig>_budget_s`, e.g. `fig2_budget_s`), seconds; 0 = no limit.
A build checks it before each heavy unit (one gffcompare, one contig's BAM pass, one SQANTI3 call) and, when it is
used up, stops with exit code 75 (Pending) before writing any table; every finished unit is cached, so the same
command continues where it stopped. `figs_plan=1` makes the figure 1-3 builds print the heavy units they would run,
with an estimated cost, and run nothing.

Our assembler always runs genome-wide through the pipeline driver (`tools/rustle_pipeline.sh assemble`), then
the GTF is restricted, exactly as the external tools' genome-wide outputs are.

Modes (apples to apples). The tool comparison is annotation-free on every side: Rustle `assemble` (reads + genome),
StringTie -L without -G, FLAIR bam2bed + collapse without -f/--gtf (flair correct skipped), IsoSeq collapse
(MODE_DENOVO_METHODS; verified in benchmark_collapse/run_stringtie.sbatch and run_flair.sbatch). Annotation-guided
tool runs (registry columns stringtie_guided_gtf / flair_guided_gtf, GUIDED_TOOLS) are scored by the same code into
separate tables and figures (fig1_guided, fig2_guided_*, fig3_guided_bins) that build only when such a GTF exists;
Rustle has no annotation-guided transcript assembly, so those have no Rustle row (RUSTLE_NO_GUIDED;
docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md).
"""
from __future__ import annotations

import gzip
import hashlib
import math
import re
import sys
import time
from pathlib import Path

import figlib
import samples

TOOLS = ["rustle", "rustle_primary", "stringtie", "flair", "isoseq"]
SPECIES = ["gorilla", "human"]  # legacy aliases of gorilla_OR6737 / human_A119b (figures/samples.tsv)

# ---------------------------------------------------------------- comparison modes (apples to apples)
# Every tool comparison of figures 1-3 is ANNOTATION-FREE on every side (verified in benchmark_collapse/
# run_stringtie.sbatch and run_flair.sbatch). The wording below is printed on every tool-comparison panel, table
# note and caption (GLOSSARY "Transcript modes").
MODE_DENOVO = "annotation-free (de novo)"
MODE_DENOVO_METHODS = ("annotation-free (de novo): Rustle assemble (reads + genome), StringTie -L without -G, FLAIR "
                       "collapse without annotation (flair correct skipped), IsoSeq collapse")
MODE_GUIDED = "annotation-guided"
# the guided runs of the same tools (user-supplied registry columns stringtie_guided_gtf / flair_guided_gtf)
GUIDED_TOOLS = ["stringtie_guided", "flair_guided"]
GUIDED_NA = "guided comparison: not available (guided StringTie/FLAIR GTFs not supplied)"
RUSTLE_NO_GUIDED = ("Rustle has no annotation-guided transcript assembly (its guided mode defines families and loci, "
                    "Figs 7-8), so the guided comparison has no Rustle row")
GUIDED_CAVEAT = ("the guided tools were given the annotation they are scored against: their numbers describe those "
                 "runs and are never compared with an annotation-free number")


def guided_tools(cfg: dict, key: str) -> list[str]:
    """The annotation-guided tool runs registered for one sample (registry stringtie_guided_gtf / flair_guided_gtf)."""
    return [t for t in GUIDED_TOOLS if samples.baseline(cfg, key, t)]


def guided_samples(cfg: dict) -> list[str]:
    """Sample keys (sample_key) with at least one annotation-guided tool run; [] until the user supplies them."""
    return [sample_key(cfg, sid) for sid in samples.registry(cfg) if guided_tools(cfg, sid)]


def guided_status(cfg: dict) -> str:
    """One line for build logs and table notes: GUIDED_NA, or which samples and tools the guided path covers."""
    keys = guided_samples(cfg)
    if not keys:
        return GUIDED_NA
    return ("guided comparison (separate table and figure, never drawn with the annotation-free methods): "
            + "; ".join(f"{k}: {', '.join(guided_tools(cfg, k))}" for k in keys))


def species_dir(cfg: dict, species: str) -> Path:
    """`${work}/assembly/<alias>` for the two legacy samples (their eval caches live there), else `<sample id>`."""
    row = samples.get(cfg, species)
    d = figlib.work_dir(cfg, "assembly") / (row["alias"] or row["id"])
    d.mkdir(parents=True, exist_ok=True)
    return d


def eval_dir(cfg: dict, species: str) -> Path:
    """Restricted GTFs and gffcompare outputs, keyed by the evaluation contig set."""
    d = species_dir(cfg, species) / f"eval_{scope_key(cfg, species)}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _natural(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def genome_contigs(cfg: dict, species: str) -> list[str]:
    """The sample's genome contigs, in .fai order (the BAM @SQ equals the .fai: `make.py samples --verify`)."""
    fasta = samples.get(cfg, species)["fasta"]
    return [l.split("\t", 1)[0] for l in open(str(fasta) + ".fai") if l.strip()]


def annotated_contigs(cfg: dict, species: str) -> list[str]:
    """Contigs the sample's annotation covers, in genome order (cached seqid scan, samples.annotation_seqids)."""
    row = samples.get(cfg, species)
    src = row["annotation_gtf"] or row["annotation_gff"]
    if not src:
        raise KeyError(f"no annotation for sample {row['id']} (figures/samples.tsv)")
    have = set(samples.annotation_seqids(cfg, src)["seqids"])
    return [c for c in genome_contigs(cfg, species) if c in have]


_WARNED: set = set()


def evaluation_contigs(cfg: dict, species: str) -> set[str] | None:
    """None = every contig (the annotation covers the whole genome); else the annotated contigs (human: all but
    chrM), or `eval_contigs_<sample id>` when that development key narrows the sample."""
    sid = samples.resolve(cfg, species)
    if "human_regions" in cfg and "human_regions" not in _WARNED and \
            cfg["human_regions"].strip().lower() not in ("all", "*", ""):
        _WARNED.add("human_regions")
        print(f"[assembly] WARNING: inputs key human_regions={cfg['human_regions']} is no longer read (every sample "
              f"is scored genome-wide); set eval_contigs_<sample id>=... to narrow a sample", file=sys.stderr)
    raw = (cfg.get(f"eval_contigs_{sid}") or "all").strip()
    if raw.lower() not in ("all", "*", ""):
        return {c.strip() for c in raw.split(",") if c.strip()}
    annotated = annotated_contigs(cfg, species)
    if set(annotated) >= set(genome_contigs(cfg, species)):
        return None
    return set(annotated)


def scope_key(cfg: dict, species: str) -> str:
    """Directory key of the evaluation scope: 'all', 'annotated' (every annotated contig), or the contig list."""
    contigs = evaluation_contigs(cfg, species)
    if contigs is None:
        return "all"
    if (cfg.get(f"eval_contigs_{samples.resolve(cfg, species)}") or "all").strip().lower() in ("all", "*", ""):
        return "annotated"
    return "-".join(sorted(contigs))


def is_genome_wide(cfg: dict, species: str) -> bool:
    return scope_key(cfg, species) in ("all", "annotated")


def unannotated_contigs(cfg: dict, species: str) -> list[str]:
    """Genome contigs outside the evaluation scope because the annotation does not cover them (human: chrM)."""
    if not is_genome_wide(cfg, species):
        return []
    have = set(annotated_contigs(cfg, species))
    return [c for c in genome_contigs(cfg, species) if c not in have]


def scope_label(cfg: dict, species: str) -> str:
    """'genome-wide' (every annotated contig) or the contig list, as printed in the tables and on the figures."""
    if is_genome_wide(cfg, species):
        return "genome-wide"
    return ",".join(sorted(evaluation_contigs(cfg, species), key=_natural))


def sample_key(cfg: dict, key: str) -> str:
    """The key the figure tables use for a sample: its legacy alias ('gorilla', 'human') if it has one, else its id."""
    row = samples.get(cfg, key)
    return row["alias"] or row["id"]


def benchmark_samples(cfg: dict) -> list[str]:
    """Sample keys (sample_key) of every registry sample with all three lab baselines (StringTie, FLAIR, IsoSeq
    collapse): the full comparison of figures 1-3. Gorilla first, then human (the panel order), then any other."""
    keys = [sample_key(cfg, sid) for sid in samples.registry(cfg)
            if all(samples.baseline(cfg, sid, t) for t in samples.BASELINES)]
    order = {"gorilla": 0, "human": 1}
    return sorted(keys, key=lambda k: order.get(k, 2))


def all_samples(cfg: dict) -> list[str]:
    """Sample keys of every registry sample, in registry order (the Rustle-only supplementary results)."""
    return [sample_key(cfg, sid) for sid in samples.registry(cfg)]


def sample_label(cfg: dict, key: str) -> str:
    """Plain label of a sample for figures and tables, e.g. 'Gorilla OR6737 (testis)'."""
    row = samples.get(cfg, key)
    ind = row["id"].split("_", 1)[1] if "_" in row["id"] else row["id"]
    tissue = row["tissue"] if row["tissue"] and row["tissue"] != "unknown" else "tissue not recorded"
    if ind == tissue:   # human_testis: the library is named by its tissue
        return f"{row['species'].capitalize()}, {tissue}"
    return f"{row['species'].capitalize()} {ind} ({tissue})"


# ---------------------------------------------------------------- annotation GTF (derived from GFF3 when needed)
def _src_fp(path) -> str:
    p = Path(path)
    try:
        st = p.stat()
        return f"{p.resolve()}|{st.st_size}|{st.st_mtime_ns}"
    except OSError:
        return f"{p}|absent"


def _file_sha1(path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def derived_annotation_path(cfg: dict, species: str) -> Path:
    return figlib.work_dir(cfg, "assembly") / samples.resolve(cfg, species) / "annotation.gff_to_gtf.gtf"


def annotation_gtf(cfg: dict, species: str, *, plan_only=False) -> Path | None:
    """The sample's annotation as GTF: the registry's `annotation_gtf`, else a GTF derived from `annotation_gff` with
    `gff_to_gtf` (one call per contig, concatenated in the GFF's contig order; cached). `plan_only` returns the
    path without building it (None when it would have to be built)."""
    row = samples.get(cfg, species)
    if row["annotation_gtf"]:
        return Path(row["annotation_gtf"])
    gff = row["annotation_gff"]
    if not gff:
        raise KeyError(f"no annotation for sample {row['id']} (figures/samples.tsv)")
    dst = derived_annotation_path(cfg, species)
    conv = Path(cfg["bin"]) / "gff_to_gtf"
    contigs = samples.annotation_seqids(cfg, gff)["seqids"]
    key = f"{_src_fp(gff)}\tgff_to_gtf={_file_sha1(conv)}\t{','.join(contigs)}"
    side = Path(str(dst) + ".key")
    if dst.exists() and side.exists() and side.read_text().rstrip("\n") == key:
        return dst
    if plan_only:
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(dst) + ".tmp")
    part = Path(str(dst) + ".part.gtf")
    with open(tmp, "w") as fo:
        for c in contigs:
            figlib.run([str(conv), str(gff), c, str(part)], log=dst.parent / "gff_to_gtf.log")
            with open(part) as fi:
                for line in fi:
                    fo.write(line)
    part.unlink(missing_ok=True)
    tmp.replace(dst)
    side.write_text(key + "\n")
    return dst


# ---------------------------------------------------------------- per-call budget and plan mode (figures 1-3)
PENDING_EXIT = 75


class Pending(SystemExit):
    """A budgeted build call stopped before its next heavy unit (every finished unit is cached): run it again."""

    def __init__(self, msg: str):
        print(f"[pending] {msg}; run the same command again to continue", file=sys.stderr)
        super().__init__(PENDING_EXIT)


class Budget:
    """Wall-time budget of one build call: cfg['<fig>_budget_s'], else cfg['figs_budget_s']; 0 or unset = none."""

    def __init__(self, cfg: dict, fig: str):
        self.fig = fig
        self.limit = float(cfg.get(f"{fig}_budget_s") or cfg.get("figs_budget_s") or 0)
        self.t0 = time.time()

    def remaining(self) -> float:
        return math.inf if self.limit <= 0 else self.limit - (time.time() - self.t0)

    def check(self, what: str):
        if self.remaining() <= 0:
            raise Pending(f"{self.fig}: the {self.limit:.0f} s budget of this call is used up; next unit: {what}")


def plan_only(cfg: dict) -> bool:
    return str(cfg.get("figs_plan", "")).strip().lower() in ("1", "true", "yes")


def restrict_gtf(src, dst: Path, contigs: set[str] | None, *, force=False) -> Path:
    """Copy GTF/GFF lines whose seqname is in `contigs` (all when None); comment lines dropped."""
    dst = Path(dst)
    if not force and figlib.fresh(dst, src):
        return dst
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with _open(src) as fi, open(tmp, "w") as fo:
        for line in fi:
            if not line or line[0] == "#":
                continue
            if contigs is None or line.split("\t", 1)[0] in contigs:
                fo.write(line)
    tmp.replace(dst)
    return dst


def ensure_rustle(cfg: dict, species: str, tool: str, *, force=False) -> Path:
    """Genome-wide assembly of a sample through the run cache (samples.ensure): `rustle` = stage `assemble`,
    `rustle_primary` = stage `assemble_primary` (--no-seed-secondaries). Fresh products are reused (run-cache key:
    BAM/FASTA fingerprints, binary hashes, driver code hash); `force` re-assembles."""
    assert tool in ("rustle", "rustle_primary")
    stage = "assemble" if tool == "rustle" else "assemble_primary"
    return samples.ensure(cfg, species, stage, force=force)["gtf"]


def _code_hash(path) -> str:
    """sha1 of a shell script's non-comment, non-blank lines."""
    import hashlib
    h = hashlib.sha1()
    for line in open(path):
        t = line.strip()
        if t and not t.startswith("#"):
            h.update(t.encode() + b"\n")
    return h.hexdigest()[:16]


RUSTLE_STAGE = {"rustle": "assemble", "rustle_primary": "assemble_primary"}


def rustle_product(cfg: dict, species: str, tool: str) -> Path:
    """The genome-wide Rustle GTF of a sample from the run cache WITHOUT assembling: fresh -> its path; adoptable
    (products newer than every input, no stamp yet) -> stamped and returned; stale or missing -> RuntimeError naming
    the `make.py runs` command (figures 1-3 never start a genome-wide assembly implicitly)."""
    stage = RUSTLE_STAGE[tool]
    sid = samples.resolve(cfg, species)
    state, reason = samples.status(cfg, sid, stage)
    if state == "fresh":
        return samples.product(cfg, sid, stage, "gtf")
    if state == "adopt":
        return samples.ensure(cfg, sid, stage)["gtf"]   # writes the stamp only, runs nothing
    raise RuntimeError(f"{sid} {stage} is {state} ({reason}): run `python3 figures/make.py runs --sample {sid} --stage "
                       f"{stage}` first (figure builds never re-assemble implicitly)")


def arm_source(cfg: dict, species: str, tool: str, *, check=True) -> Path:
    """Genome-wide GTF of one method for a sample: Rustle from the run cache (rustle_product; `check=False` only
    names the path), the others the lab's annotation-free baselines (samples.baseline), or, for GUIDED_TOOLS, the
    user's annotation-guided runs (the separate guided path only)."""
    if tool in RUSTLE_STAGE:
        return rustle_product(cfg, species, tool) if check else samples.product(cfg, species, RUSTLE_STAGE[tool], "gtf")
    src = samples.baseline(cfg, species, tool)
    if not src:
        if tool in GUIDED_TOOLS:
            raise KeyError(f"no {tool} GTF for sample {samples.resolve(cfg, species)} ({GUIDED_NA})")
        raise KeyError(f"no {tool} baseline for sample {samples.resolve(cfg, species)} (the lab ran the tools on "
                       "human_A119b and gorilla_OR6737 only)")
    return Path(src)


def ensure_arm(cfg: dict, species: str, tool: str, *, force=False, force_arms=False) -> Path:
    """`<tool>.gtf` restricted to the evaluation contigs (`force` re-restricts; `force_arms` re-assembles through
    the run cache; otherwise a stale or missing assembly raises, see rustle_product)."""
    d = eval_dir(cfg, species)
    contigs = evaluation_contigs(cfg, species)
    if tool in RUSTLE_STAGE and force_arms:
        src = ensure_rustle(cfg, species, tool, force=True)
    else:
        src = arm_source(cfg, species, tool)
    return restrict_gtf(src, d / f"{tool}.gtf", contigs, force=force)


def outside_scope(cfg: dict, species: str, tool: str) -> dict:
    """{contig: transcripts} of one method's genome-wide GTF on the contigs the evaluation leaves out (human: chrM,
    which the annotation does not cover); {} when nothing is left out. Cached in the evaluation dir (keyed by the
    source's path, size and mtime), so the source is read once."""
    contigs = evaluation_contigs(cfg, species)
    if contigs is None:
        return {}
    src = arm_source(cfg, species, tool, check=False)
    f = eval_dir(cfg, species) / f"{tool}.outside_scope.tsv"
    key = f"# {_src_fp(src)}\t{','.join(sorted(contigs))}"
    if f.exists():
        lines = f.read_text().splitlines()
        if lines and lines[0] == key:
            return {c: int(n) for c, n in (l.split("\t") for l in lines[1:] if l)}
    ids: dict = {}
    tid = re.compile(r'transcript_id[ =]"?([^";]+)')
    with _open(src) as fh:
        for line in fh:
            c = line.split("\t", 1)[0]
            if c in contigs or not line or line[0] == "#":
                continue
            m = tid.search(line)
            if m:
                ids.setdefault(c, set()).add(m.group(1))
    out = {c: len(v) for c, v in sorted(ids.items(), key=lambda kv: _natural(kv[0]))}
    f.write_text(key + "\n" + "".join(f"{c}\t{n}\n" for c, n in out.items()))
    return out


def gffcompare_state(cfg: dict, species: str, tool: str) -> str:
    """'' when the cached gffcompare of one method is current, else the work gffcompare() would do (plan / budget).
    Reads file times only."""
    d = eval_dir(cfg, species)
    ann = annotation_gtf(cfg, species, plan_only=True)
    if ann is None:
        return "derive the annotation GTF from the GFF3, restrict, gffcompare"
    src = arm_source(cfg, species, tool, check=False)
    if not Path(src).exists():
        return f"no {tool} GTF yet ({src})"
    ref, q = d / "ref.gtf", d / f"{tool}.gtf"
    if not (figlib.fresh(ref, ann) and figlib.fresh(q, src)):
        return "restrict, gffcompare"
    stats = d / f"gc_{tool}.stats"
    if not (figlib.fresh(stats, q, ref) and (d / f"gc_{tool}.{q.name}.tmap").exists()
            and (d / f"gc_{tool}.{q.name}.refmap").exists()):
        return "gffcompare"
    return ""


def ensure_ref(cfg: dict, species: str, *, force=False) -> Path:
    d = eval_dir(cfg, species)
    return restrict_gtf(annotation_gtf(cfg, species), d / "ref.gtf", evaluation_contigs(cfg, species), force=force)


def gffcompare(cfg: dict, species: str, tool: str, *, force=False, force_arms=False) -> dict:
    """Run gffcompare of the arm against the restricted annotation; return {'stats','tmap','refmap'} paths.

    `force` recomputes the restrictions and gffcompare; the genome-wide assemblies are recomputed only with
    `force_arms` (they take minutes to tens of minutes each). gffcompare writes `<prefix>.<query basename>.tmap/
    .refmap` NEXT TO THE QUERY, so the query lives in the evaluation dir and the prefix is `gc_<tool>`."""
    d = eval_dir(cfg, species)
    ref = ensure_ref(cfg, species, force=force)
    q = ensure_arm(cfg, species, tool, force=force, force_arms=force_arms)
    prefix = d / f"gc_{tool}"
    stats = Path(str(prefix) + ".stats")
    tmap = d / f"gc_{tool}.{q.name}.tmap"
    refmap = d / f"gc_{tool}.{q.name}.refmap"
    if force or not (figlib.fresh(stats, q, ref) and tmap.exists() and refmap.exists()):
        figlib.run([cfg.get("gffcompare", "gffcompare"), "-r", str(ref), "-o", str(prefix), str(q)],
                   log=d / f"gc_{tool}.log", cwd=d)
        # older gffcompare versions write the summary to the bare prefix
        if not stats.exists() and prefix.exists():
            prefix.replace(stats)
        for ext in (".annotated.gtf", ".loci", ".tracking", ".combined.gtf"):
            Path(str(prefix) + ext).unlink(missing_ok=True)
    return {"stats": stats, "tmap": tmap, "refmap": refmap, "query": q, "ref": ref}


_LEVEL = re.compile(r"^\s*(Base|Exon|Intron|Intron chain|Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)")
_MATCH = re.compile(r"^\s*Matching (intron chains|transcripts|loci):\s+(\d+)")
_QUERY = re.compile(r"Query mRNAs\s*:\s*(\d+) in\s+(\d+) loci\s+\((\d+) multi-exon")
_REF = re.compile(r"Reference mRNAs\s*:\s*(\d+) in\s+(\d+) loci\s+\((\d+) multi-exon")


def parse_stats(path) -> dict:
    """gffcompare .stats -> {'<level>_sn','<level>_pr', 'matching_<what>', 'query_mrnas', 'ref_mrnas', ...}."""
    out: dict = {}
    for line in open(path):
        m = _LEVEL.match(line)
        if m:
            k = m.group(1).lower().replace(" ", "_")
            out[f"{k}_sn"], out[f"{k}_pr"] = float(m.group(2)), float(m.group(3))
            continue
        m = _MATCH.match(line)
        if m:
            out["matching_" + m.group(1).replace(" ", "_")] = int(m.group(2))
            continue
        m = _QUERY.search(line)
        if m:
            out["query_mrnas"], out["query_loci"], out["query_multiexon"] = map(int, m.groups())
            continue
        m = _REF.search(line)
        if m:
            out["ref_mrnas"], out["ref_loci"], out["ref_multiexon"] = map(int, m.groups())
    return out


def read_tmap(path):
    """gffcompare .tmap rows as dicts (query-centred: class_code, ref_id, qry_id, num_exons, ...)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def exact_matched_refs(tmap_path) -> set[str]:
    """Reference transcript ids matched with class code '=' (identical intron chain) by at least one query."""
    return {r["ref_id"] for r in read_tmap(tmap_path) if r.get("class_code") == "=" and r.get("ref_id") not in ("", "-")}


def ref_transcripts(ref_gtf) -> dict:
    """{transcript_id: {'gene', 'chrom', 'strand', 'exons': [(s0, e)], 'introns': ((d, a), ...)}} (0-based)."""
    tx: dict = {}
    for line in _open(ref_gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8])
        gid = re.search(r'gene_id "([^"]+)"', f[8])
        if not tid:
            continue
        t = tx.setdefault(tid.group(1), {"gene": gid.group(1) if gid else "", "chrom": f[0], "strand": f[6],
                                         "exons": []})
        t["exons"].append((int(f[3]) - 1, int(f[4])))
    for t in tx.values():
        t["exons"].sort()
        t["introns"] = tuple((a[1], b[0]) for a, b in zip(t["exons"], t["exons"][1:]))
    return tx
