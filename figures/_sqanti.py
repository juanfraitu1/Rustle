"""_sqanti — private helpers for fig2: SQANTI3 QC + rules filter on one method, and the tallies the figure plots.

Protocol (bench/SQANTI3_POLISH.md, 2026-09-19), nothing customised; SQANTI3 5.5.4:
    sqanti3_qc.py --isoforms ARM.gtf --refGTF REF.gtf --refFasta GENOME.fa --report skip -t N -d DIR -o PREFIX
    sqanti3_filter.py rules --sqanti_class DIR/PREFIX_classification.txt --filter_gtf DIR/PREFIX_corrected.gtf
                            --skip_report -d FILTDIR -o PREFIX            (the default rules JSON)

SQANTI3 runs ONE CONTIG AT A TIME (as the recorded per-chromosome runs did): the genome, the annotation and the
method's GTF are cut to that contig. A structural category is decided per isoform against the annotation of its own
contig and every rule of the default filter JSON (perc_A_downstream_TTS, all_canonical, RTS_stage, min_cov) is per
isoform, so the per-contig tallies add up to the tally of one run over all the contigs; each run stays short and is
cached on its own, so an interrupted build resumes at the next contig. A method with more than `fig2_chunk_tx`
transcripts on a contig runs in parts of whole transcripts (same argument). Contigs: GENOME-WIDE by default, i.e.
every contig of the sample's evaluation scope (assembly.evaluation_contigs: every contig the annotation covers; the
human annotation leaves out chrM); `fig2_contigs_<sample id>` (comma list) narrows one sample for development. In the
reference every record with an empty gene_id (gene-level RefSeq records without an mRNA child, e.g. 6.6% of the
gorilla GTF) gets its own transcript id as gene_id, so SQANTI3 does not merge all of them into one gene called "".

Per-contig products (the genome and annotation cut, the method's cut) are REPLACED ONLY WHEN THEIR CONTENT CHANGES:
a re-split that writes the same bytes keeps the old file and its time stamp, so an unchanged contig never re-runs
SQANTI3 after, e.g., a re-assembly that reproduces the same transcripts or a wider contig set.
"""
from __future__ import annotations

import filecmp
import os
import re
import sys
from pathlib import Path

import figlib

DEFAULT_ARMS = ["rustle", "rustle_primary", "stringtie", "flair", "isoseq"]  # = fig1's methods; cfg['fig2_arms'] narrows it
RUSTLE_ARMS = ["rustle", "rustle_primary"]   # the supplementary samples (no lab baselines)
LEGACY_KEYS = ("gorilla_sqanti_contigs", "human_regions")
_WARNED: set = set()


def arms(cfg: dict, key: str | None = None) -> list[str]:
    """Methods SQANTI3 classifies for a sample: the five of fig. 1 for a sample with lab baselines (narrowed by
    cfg['fig2_arms']), Rustle's two configurations for every other sample."""
    import assembly
    raw = cfg.get("fig2_arms", "")
    sel = [a.strip() for a in raw.split(",") if a.strip()] if raw else DEFAULT_ARMS
    if key is not None and key not in assembly.benchmark_samples(cfg):
        sel = [a for a in sel if a in RUSTLE_ARMS] or RUSTLE_ARMS
    return [a for a in figlib.TOOL_ORDER if a in sel]


def sqanti_contigs(cfg: dict, key: str) -> list[str]:
    """Contigs SQANTI3 runs on for a sample: its evaluation scope (genome-wide: every annotated contig), in genome
    order, unless `fig2_contigs_<sample id>` narrows it (development)."""
    import assembly
    import samples
    for k in LEGACY_KEYS:
        if k in cfg and k not in _WARNED and cfg[k].strip().lower() not in ("all", "*", ""):
            _WARNED.add(k)
            print(f"[fig2] WARNING: inputs key {k}={cfg[k]} is no longer read (SQANTI3 runs genome-wide); set "
                  "fig2_contigs_<sample id>=... to narrow a sample", file=sys.stderr)
    sid = samples.resolve(cfg, key)
    raw = (cfg.get(f"fig2_contigs_{sid}") or "").strip()
    if raw and raw.lower() not in ("all", "*"):
        return [c.strip() for c in raw.split(",") if c.strip()]
    ev = assembly.evaluation_contigs(cfg, key)
    return [c for c in assembly.genome_contigs(cfg, key) if ev is None or c in ev]


def is_genome_wide(cfg: dict, key: str) -> bool:
    import assembly
    import samples
    raw = (cfg.get(f"fig2_contigs_{samples.resolve(cfg, key)}") or "").strip().lower()
    return raw in ("", "all", "*") and assembly.is_genome_wide(cfg, key)


def reference_home(cfg: dict, key: str) -> str:
    """The sample key whose per-contig genome and annotation cuts this sample shares (the first registry sample with
    the same FASTA and annotation: human_testis -> human, gorilla_KB3781 -> gorilla)."""
    import assembly
    import samples
    row = samples.get(cfg, key)
    for sid, r in samples.registry(cfg).items():
        if r["fasta"] == row["fasta"] and (r["annotation_gtf"], r["annotation_gff"]) == \
                (row["annotation_gtf"], row["annotation_gff"]):
            return assembly.sample_key(cfg, sid)
    return assembly.sample_key(cfg, key)


def sqanti_env(cfg: dict) -> dict:
    env_bin = str(Path(cfg["sqanti3_env"]) / "bin")
    return {"PATH": env_bin + os.pathsep + os.environ.get("PATH", "")}


# ---- cache keys: a cached product is reused only while the mtime check passes AND its `.key` sidecar (what it was
# made from: source path + size + mtime, the SQANTI3 install + version) still matches. A product cached before the
# sidecars existed is adopted once with the current key (it passed the mtime check, which was the old rule).
def _src_key(src) -> str:
    p = Path(src)
    try:
        st = p.stat()
        return f"{p.resolve()}\t{st.st_size}\t{int(st.st_mtime)}"
    except OSError:
        return f"{p}\tabsent"


def sqanti_version(cfg: dict) -> str:
    """'<sqanti3_dir> <version>' from the install's src/config.py (a changed install re-runs every cached call)."""
    d = Path(cfg["sqanti3_dir"])
    try:
        m = re.search(r"__version__\s*=\s*['\"]([^'\"]+)", (d / "src" / "config.py").read_text())
        return f"{d.resolve()} {m.group(1) if m else '?'}"
    except OSError:
        return f"{d} ?"


def _key_ok(dst: Path, key: str) -> bool:
    side = Path(str(dst) + ".key")
    if not side.exists():
        side.write_text(key + "\n")   # adopt a product cached before the sidecars existed
        return True
    return side.read_text().rstrip("\n") == key


def _key_write(dst: Path, key: str):
    Path(str(dst) + ".key").write_text(key + "\n")


def _cached_ok(dst: Path, src, key: str) -> bool:
    """A per-contig cut is current when its key sidecar (source path, size, mtime) matches; a cut made before the
    sidecars existed is adopted by the old time-stamp rule. (Time stamps alone are not enough here: a re-cut that
    wrote the same bytes keeps the older file on purpose, see _replace_if_changed.)"""
    dst = Path(dst)
    if not dst.exists():
        return False
    side = Path(str(dst) + ".key")
    if side.exists():
        return side.read_text().rstrip("\n") == key
    return figlib.fresh(dst, src) and _key_ok(dst, key)


def _replace_if_changed(tmp: Path, dst: Path) -> bool:
    """Move `tmp` over `dst` unless `dst` already holds the same bytes (then `tmp` is dropped and `dst` keeps its
    time stamp, so the SQANTI3 runs cached on it stay valid). Returns True when `dst` changed."""
    if dst.exists() and filecmp.cmp(tmp, dst, shallow=False):
        tmp.unlink()
        return False
    tmp.replace(dst)
    return True


def subset_genome(fasta, contigs, dst: Path, *, force=False) -> Path:
    """The genome cut to `contigs` (+ .fai), via pysam (light: reads only those contigs)."""
    import pysam

    dst = Path(dst)
    key = _src_key(fasta) + "\t" + ",".join(contigs)
    if not force and figlib.fresh(dst, fasta) and Path(str(dst) + ".fai").exists() and _key_ok(dst, key):
        return dst
    tmp = Path(str(dst) + ".tmp")
    with pysam.FastaFile(str(fasta)) as fa, open(tmp, "w") as fo:
        for c in contigs:
            seq = fa.fetch(c)
            fo.write(f">{c}\n")
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i + 80] + "\n")
    tmp.replace(dst)
    pysam.faidx(str(dst))
    _key_write(dst, key)
    return dst


_TID = re.compile(r'transcript_id "([^"]*)"')


def sqanti_references(src_gtf, dst_of: dict, *, force=False) -> dict:
    """{contig: reference GTF of that contig} in one pass over the annotation; empty gene_id -> the record's
    transcript_id (see module doc). Cached: rewritten only when a file is missing or older than the source."""
    import gzip

    dst_of = {c: Path(d) for c, d in dst_of.items()}
    key = _src_key(src_gtf) + "\treference"
    if not force and all(_cached_ok(d, src_gtf, key) for d in dst_of.values()):
        return dst_of
    opener = gzip.open if str(src_gtf).endswith(".gz") else open
    tmp = {c: Path(str(d) + ".tmp") for c, d in dst_of.items()}
    fos = {c: open(t, "w") for c, t in tmp.items()}
    try:
        with opener(src_gtf, "rt") as fi:
            for line in fi:
                if not line or line[0] == "#":
                    continue
                fo = fos.get(line.split("\t", 1)[0])
                if fo is None:
                    continue
                if 'gene_id ""' in line:
                    m = _TID.search(line)
                    if m and m.group(1):
                        line = line.replace('gene_id ""', f'gene_id "{m.group(1)}"', 1)
                fo.write(line)
    finally:
        for fo in fos.values():
            fo.close()
    for c in dst_of:
        _replace_if_changed(tmp[c], dst_of[c])
        _key_write(dst_of[c], key)
    return dst_of


def split_gtf(src_gtf, dst_of: dict, *, force=False) -> dict:
    """{contig: the arm's records on that contig} in one pass (comment lines dropped); cached on the source (its
    path, size and mtime: pointing an arm at another GTF re-splits even when that GTF is older than the cache)."""
    import gzip

    dst_of = {c: Path(d) for c, d in dst_of.items()}
    key = _src_key(src_gtf) + "\tsplit"
    if not force and all(_cached_ok(d, src_gtf, key) for d in dst_of.values()):
        return dst_of
    opener = gzip.open if str(src_gtf).endswith(".gz") else open
    tmp = {c: Path(str(d) + ".tmp") for c, d in dst_of.items()}
    fos = {c: open(t, "w") for c, t in tmp.items()}
    try:
        with opener(src_gtf, "rt") as fi:
            for line in fi:
                if not line or line[0] == "#":
                    continue
                fo = fos.get(line.split("\t", 1)[0])
                if fo is not None:
                    fo.write(line)
    finally:
        for fo in fos.values():
            fo.close()
    for c in dst_of:
        _replace_if_changed(tmp[c], dst_of[c])
        _key_write(dst_of[c], key)
    return dst_of


def n_transcripts(gtf) -> int:
    """Distinct transcript_id values of a GTF (SQANTI3 is not run on a contig where an arm has none)."""
    ids = set()
    with open(gtf) as fh:
        for line in fh:
            m = _TID.search(line)
            if m:
                ids.add(m.group(1))
    return len(ids)


def chunk_gtf(gtf, max_tx: int, *, force=False) -> list[Path]:
    """[gtf] when it has <= max_tx transcripts; otherwise equal parts `<gtf>.partK.gtf` of whole transcripts (lines
    grouped by transcript_id, first-appearance order; records without a transcript_id are dropped). SQANTI3's
    categories and the default rules filter are per isoform against the isoform's own contig, so per-part tallies
    sum to the whole-contig tally; parts only keep each SQANTI3 call inside the foreground time budget."""
    gtf = Path(gtf)
    n = n_transcripts(gtf)
    if n <= max_tx:
        return [gtf]
    k = -(-n // max_tx)
    parts = [gtf.with_name(f"{gtf.stem}.part{i}.gtf") for i in range(k)]
    # reuse only a partition that is newer than the GTF AND covers it exactly: after fig2_chunk_tx changes k, the
    # first parts of the old partition are still newer than the GTF but hold only part of its transcripts
    if not force and all(figlib.fresh(p, gtf) for p in parts) and sum(n_transcripts(p) for p in parts) == n:
        return parts
    order, lines = [], {}
    with open(gtf) as fh:
        for line in fh:
            m = _TID.search(line)
            if not m:
                continue
            t = m.group(1)
            if t not in lines:
                order.append(t)
                lines[t] = []
            lines[t].append(line)
    per = -(-len(order) // k)
    for i, p in enumerate(parts):
        tmp = Path(str(p) + ".tmp")
        with open(tmp, "w") as fo:
            for t in order[i * per:(i + 1) * per]:
                fo.writelines(lines[t])
        tmp.replace(p)
    return parts


def _clear(d: Path):
    import shutil

    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)


def qc_cached(cfg: dict, isoforms, ref_gtf, genome, outdir: Path, prefix: str) -> bool:
    """True when run_qc would reuse its cached classification (reads time stamps and the key sidecar only)."""
    outdir = Path(outdir)
    cls = outdir / f"{prefix}_classification.txt"
    side = Path(str(cls) + ".key")
    key = "\t".join([sqanti_version(cfg), _src_key(isoforms), _src_key(ref_gtf), _src_key(genome)])
    return (figlib.fresh(cls, isoforms, ref_gtf, genome) and (outdir / f"{prefix}_corrected.gtf").exists()
            and (not side.exists() or side.read_text().rstrip("\n") == key))


def filter_cached(cfg: dict, classification, outdir: Path, prefix: str) -> bool:
    res = Path(outdir) / f"{prefix}_RulesFilter_result_classification.txt"
    side = Path(str(res) + ".key")
    key = "\t".join([sqanti_version(cfg), _src_key(classification)])
    return figlib.fresh(res, classification) and (not side.exists() or side.read_text().rstrip("\n") == key)


def run_qc(cfg: dict, isoforms, ref_gtf, genome, outdir: Path, prefix: str, *, force=False) -> dict:
    """SQANTI3 QC (HEAVY; foreground). Returns {'classification', 'corrected_gtf'}; cached on the inputs."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cls = outdir / f"{prefix}_classification.txt"
    gtf = outdir / f"{prefix}_corrected.gtf"
    log = outdir.parent / f"{outdir.name}.log"
    key = "\t".join([sqanti_version(cfg), _src_key(isoforms), _src_key(ref_gtf), _src_key(genome)])
    if force or not (figlib.fresh(cls, isoforms, ref_gtf, genome) and gtf.exists() and _key_ok(cls, key)):
        # SQANTI3 silently REUSES corrected FASTA / ORFs / genePreds found in its output dir, so a re-run after the
        # arm GTF changed would classify the OLD isoforms: start from an empty directory (ours, under ${work}).
        _clear(outdir)
        cmd = [str(Path(cfg["sqanti3_env"]) / "bin" / "python"), str(Path(cfg["sqanti3_dir"]) / "sqanti3_qc.py"),
               "--isoforms", str(isoforms), "--refGTF", str(ref_gtf), "--refFasta", str(genome),
               "--report", "skip", "-t", str(cfg.get("threads", "4")), "-d", str(outdir), "-o", prefix]
        figlib.run(cmd, log=log, env=sqanti_env(cfg), cwd=outdir.parent)
        if cls.exists():
            _key_write(cls, key)
    if not cls.exists():
        raise RuntimeError(f"SQANTI3 QC wrote no classification: {cls} (see {log})")
    return {"classification": cls, "corrected_gtf": gtf}


def run_filter(cfg: dict, qc: dict, outdir: Path, prefix: str, *, force=False) -> Path:
    """SQANTI3 rules filter with the default JSON; returns the RulesFilter classification (has `filter_result`)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    res = outdir / f"{prefix}_RulesFilter_result_classification.txt"
    log = outdir.parent / f"{outdir.name}.log"
    key = "\t".join([sqanti_version(cfg), _src_key(qc["classification"])])
    if force or not (figlib.fresh(res, qc["classification"]) and _key_ok(res, key)):
        _clear(outdir)
        cmd = [str(Path(cfg["sqanti3_env"]) / "bin" / "python"), str(Path(cfg["sqanti3_dir"]) / "sqanti3_filter.py"),
               "rules", "--sqanti_class", str(qc["classification"]), "--filter_gtf", str(qc["corrected_gtf"]),
               "--skip_report", "-d", str(outdir), "-o", prefix, "-c", str(cfg.get("threads", "4"))]
        figlib.run(cmd, log=log, env=sqanti_env(cfg), cwd=outdir.parent)
        if res.exists():
            _key_write(res, key)
    if not res.exists():
        raise RuntimeError(f"SQANTI3 rules filter wrote no result: {res} (see {log})")
    return res


def _rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def fold(category: str) -> str:
    return figlib.SQANTI_CATEGORY_MAP.get(category, "Other")


def tally(classification, filtered=None) -> dict:
    """{'n': total isoforms, 'fine': {structural_category: n}, 'sub': {(structural_category, subcategory): [n,
    n_multiexon]}, 'n_multi' (isoforms with > 1 exon), 'n_fsm', 'n_fsm_multi', 'n_pass', 'n_fsm_pass'} of one arm.

    PASS = filter_result 'Isoform' in the RulesFilter classification (the other value is 'Artifact')."""
    fine: dict = {}
    sub: dict = {}
    n = n_multi = n_fsm_multi = 0
    for r in _rows(classification):
        n += 1
        c = r["structural_category"]
        multi = int(r["exons"]) > 1
        fine[c] = fine.get(c, 0) + 1
        cell = sub.setdefault((c, r.get("subcategory", "")), [0, 0])
        cell[0] += 1
        cell[1] += multi
        n_multi += multi
        n_fsm_multi += multi and c == "full-splice_match"
    out = {"n": n, "fine": fine, "sub": sub, "n_multi": n_multi, "n_fsm": fine.get("full-splice_match", 0),
           "n_fsm_multi": n_fsm_multi, "n_pass": None, "n_fsm_pass": None}
    if filtered is not None:
        n_pass = n_fsm_pass = n_f = 0
        for r in _rows(filtered):
            n_f += 1
            if r.get("filter_result") == "Isoform":
                n_pass += 1
                if r["structural_category"] == "full-splice_match":
                    n_fsm_pass += 1
        if n_f != n:
            raise RuntimeError(f"filter result has {n_f} isoforms, QC classification {n}: {filtered}")
        out.update(n_pass=n_pass, n_fsm_pass=n_fsm_pass)
    return out


def merge_tallies(ts: list[dict]) -> dict:
    """Sum of per-contig tallies (same keys as tally())."""
    out = {"n": 0, "fine": {}, "sub": {}, "n_multi": 0, "n_fsm": 0, "n_fsm_multi": 0, "n_pass": 0, "n_fsm_pass": 0}
    for t in ts:
        for k in ("n", "n_multi", "n_fsm", "n_fsm_multi"):
            out[k] += t[k]
        for k, v in t["fine"].items():
            out["fine"][k] = out["fine"].get(k, 0) + v
        for k, (a, b) in t["sub"].items():
            cell = out["sub"].setdefault(k, [0, 0])
            cell[0] += a
            cell[1] += b
        for k in ("n_pass", "n_fsm_pass"):
            out[k] = None if out[k] is None or t[k] is None else out[k] + t[k]
    return out


EMPTY_TALLY = {"n": 0, "fine": {}, "sub": {}, "n_multi": 0, "n_fsm": 0, "n_fsm_multi": 0, "n_pass": 0,
               "n_fsm_pass": 0}


CATEGORY_HEADER = ["species", "scope", "tool", "structural_category", "category", "n", "n_total", "frac"]
FILTER_HEADER = ["species", "scope", "tool", "n_total", "n_pass", "pass_frac", "n_fsm", "fsm_frac", "n_fsm_pass",
                 "n_multiexon", "n_fsm_multiexon", "fsm_multiexon_frac"]
SUBCATEGORY_HEADER = ["species", "scope", "tool", "structural_category", "subcategory", "n", "n_multiexon",
                      "n_total", "n_total_multiexon", "frac"]


def category_rows(species: str, scope: str, tool: str, t: dict) -> list[list]:
    rows = []
    for cat, k in sorted(t["fine"].items(), key=lambda kv: (figlib.SQANTI_ORDER.index(fold(kv[0])), -kv[1])):
        rows.append([species, scope, tool, cat, fold(cat), k, t["n"], k / t["n"] if t["n"] else None])
    return rows


def subcategory_rows(species: str, scope: str, tool: str, t: dict) -> list[list]:
    """SQANTI3 subcategory counts (all isoforms and multi-exon isoforms) of every structural category."""
    rows = []
    order = lambda kv: (figlib.SQANTI_ORDER.index(fold(kv[0][0])), kv[0][0], -kv[1][0], kv[0][1])
    for (cat, subcat), (k, k_multi) in sorted(t["sub"].items(), key=order):
        rows.append([species, scope, tool, cat, subcat, k, k_multi, t["n"], t["n_multi"],
                     k / t["n"] if t["n"] else None])
    return rows


def filter_row(species: str, scope: str, tool: str, t: dict) -> list:
    n, nm = t["n"], t["n_multi"]
    return [species, scope, tool, n, t["n_pass"], (t["n_pass"] / n) if n and t["n_pass"] is not None else None,
            t["n_fsm"], (t["n_fsm"] / n) if n else None, t["n_fsm_pass"], nm, t["n_fsm_multi"],
            (t["n_fsm_multi"] / nm) if nm else None]
