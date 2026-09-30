"""_lrc — %LRC (LRGASP long-read coverage) of transcript models: a REPORTING metric only (no filter, no default, no
decision rule, no tool ranking; docs/PREREG_lrc_metric_2026-09-29.md).

Definition (LRGASP, Pardo-Palacios et al. 2024, Nat. Methods 21:1349, Box 1): "%LRC — Fraction of the transcript model
sequence length mapped by one or more long reads." The paper fixes no procedure, so the prereg does:
    reads     the sample's registry BAM (the one every arm was assembled from), PRIMARY alignments only
              (samtools -F 2308: no unmapped, secondary or supplementary record); no MAPQ filter, so a multi-mapped
              read counts once, at its primary placement
    covered   a reference base under an aligned read base: CIGAR M, = or X (pysam get_blocks; `samtools depth`
              without -J). Introns (N), deletions (D) and clips are not covered; strand is ignored
    model     one GTF transcript (transcript_id on one contig); its exon records merged, L = exonic bases
    %LRC      covered exonic bases / L. Exon level, not span level: a read spliced across a retained intron does
              not cover it
Classes (LRGASP Extended Data Fig. 2, the printed labels): > 0.98, 0.75-0.98 (both ends included), < 0.75; decided
in integers (covered * 50 > 49 * L; covered * 4 < 3 * L), so a model at exactly 0.98 or 0.75 is in the middle class.

Products (cached under ${work}/lrc/<sample id>/):
    union/<contig>.npz, union/manifest.tsv
                     the union of the primary aligned blocks of one contig (sorted, disjoint, 0-based half-open),
                     one BAM pass per sample, contig by contig; the manifest's first line is the key (BAM
                     fingerprint + read filter): a changed BAM starts the union over. With a budget the pass stops
                     before the next contig once it is used up (assembly.Pending, exit 75 = run the same command
                     again), so every call stays a light job
    <tool>.lrc.tsv   per model on the sample's evaluation contigs (assembly.evaluation_contigs; human: all but chrM):
                     contig, transcript_id, n_exons (exon records), length, covered, lrc; <tool>.lrc.npz holds the
                     numeric columns for the tables; both keyed (.key) on the union key and the GTF's fingerprint
Tables (`tables`; figlib.write_table with provenance, to --data, default ${work}/lrc/tables; no figure reads them):
    lrc_summary      species, sample, scope, tool, subset (all | multi-exon), n_models, mean_lrc, median_lrc,
                     frac_gt98, frac_75to98, frac_lt75, frac_full (= 1), frac_zero (= 0)
    lrc_sqanti       species, sample, scope, tool, category (fig. 2 fold FSM/ISM/NIC/NNC/Other), n_models, mean_lrc,
                     median_lrc, frac_gt98, frac_75to98, frac_lt75: the per-contig GTFs SQANTI3 classified for
                     fig. 2 (${work}/fig2/<sample>/<contig>/<tool>.gtf, only where cached), joined by isoform id; a
                     classification whose input changed since SQANTI3 ran is skipped and counted in the notes

Arms: rustle and rustle_primary from the run cache on every sample; the lab's annotation-free StringTie / FLAIR /
IsoSeq collapse GTFs where the registry has them (human_A119b, gorilla_OR6737). Samples and species are never pooled.

    python3 figures/_lrc.py union  --sample ID [--budget-s 150] [--threads 2]   the BAM pass (exit 75 = run again)
    python3 figures/_lrc.py score  --sample ID [--tools rustle,...]             per-model %LRC of every arm (light)
    python3 figures/_lrc.py tables [--samples ID,...] [--data DIR]              summary tables of the scored samples
"""
from __future__ import annotations

import gzip
import itertools
import re
import sys
import time
from array import array
from pathlib import Path

import numpy as np

import figlib

EXCLUDE_FLAGS = 0x4 | 0x100 | 0x800   # 2308 (repo invariant): unmapped, secondary, supplementary
READ_FILTER = "primary only (-F 2308), no MAPQ filter; covered = CIGAR M/=/X (pysam get_blocks); strand ignored"
COMPACT_AT = 1 << 23                  # flat block coordinates (int64) held before folding them into the union
ARMS = ["rustle", "rustle_primary", "stringtie", "flair", "isoseq"]
CLASS_COLUMNS = ["frac_gt98", "frac_75to98", "frac_lt75"]
SUMMARY_HEADER = ["species", "sample", "scope", "tool", "subset", "n_models", "mean_lrc", "median_lrc"] \
    + CLASS_COLUMNS + ["frac_full", "frac_zero"]
SQANTI_HEADER = ["species", "sample", "scope", "tool", "category", "n_models", "mean_lrc", "median_lrc"] \
    + CLASS_COLUMNS


# ---------------------------------------------------------------- the metric (pure; figures/test_lrc.py)
def merge_intervals(starts, ends):
    """Union of half-open intervals -> (starts, ends) sorted and disjoint (touching intervals are joined)."""
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    if s.size == 0:
        return s.copy(), e.copy()
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    run = np.maximum.accumulate(e)
    new = np.flatnonzero(s[1:] > run[:-1]) + 1      # an interval starting past every earlier end opens a new one
    return s[np.concatenate(([0], new))], run[np.concatenate((new - 1, [s.size - 1]))]


def covered_bases(starts, ends, ustarts, uends):
    """Bases of each interval [starts[i], ends[i]) that lie inside the union (ustarts, uends) (sorted, disjoint)."""
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    if len(ustarts) == 0:
        return np.zeros(s.size, dtype=np.int64)
    lens = uends - ustarts
    before = np.concatenate(([0], np.cumsum(lens)))  # union bases left of interval j

    def upto(x):   # union bases in (-inf, x)
        j = np.searchsorted(ustarts, x, side="right") - 1
        jj = np.maximum(j, 0)
        return np.where(j >= 0, before[jj] + np.clip(x - ustarts[jj], 0, lens[jj]), 0)

    return upto(e) - upto(s)


def model_lrc(tidx, starts, ends, n_models: int, ustarts, uends):
    """(n_exons, length, covered) per model from its exon records (`tidx` = model index of each record). A model's
    records are merged before measuring, so an overlapping or duplicated exon counts once; n_exons counts records."""
    t = np.asarray(tidx, dtype=np.int64)
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    n_exons = np.bincount(t, minlength=n_models)
    if t.size == 0:
        z = np.zeros(n_models, dtype=np.int64)
        return n_exons, z, z.copy()
    span = int(e.max()) + 1          # each model on its own stretch of the number line: merging never crosses models
    ms, me = merge_intervals(s + t * span, e + t * span)
    mt = ms // span
    ms, me = ms - mt * span, me - mt * span
    cov = covered_bases(ms, me, ustarts, uends)
    length = np.bincount(mt, weights=me - ms, minlength=n_models).astype(np.int64)
    covered = np.bincount(mt, weights=cov, minlength=n_models).astype(np.int64)
    return n_exons, length, covered


def classes(length, covered) -> dict:
    """Boolean masks of the LRGASP classes, in integers (exactly 0.98 and 0.75 fall in the middle class)."""
    L = np.asarray(length, dtype=np.int64)
    c = np.asarray(covered, dtype=np.int64)
    gt98 = c * 50 > 49 * L
    lt75 = c * 4 < 3 * L
    return {"frac_gt98": gt98, "frac_75to98": ~gt98 & ~lt75, "frac_lt75": lt75}


def summarize(length, covered) -> list:
    """[n_models, mean, median, frac_gt98, frac_75to98, frac_lt75, frac_full, frac_zero] of a set of models."""
    L = np.asarray(length, dtype=np.int64)
    c = np.asarray(covered, dtype=np.int64)
    n = int(L.size)
    if n == 0:
        return [0] + [None] * 7
    lrc = c / L
    cls = classes(L, c)
    return [n, float(lrc.mean()), float(np.median(lrc))] + [float(cls[k].mean()) for k in CLASS_COLUMNS] \
        + [float((c == L).mean()), float((c == 0).mean())]


# ---------------------------------------------------------------- the BAM pass (union of primary aligned blocks)
def _fold(us, ue, flat: array):
    b = np.frombuffer(flat, dtype=np.int64).reshape(-1, 2) if len(flat) else np.empty((0, 2), dtype=np.int64)
    return merge_intervals(np.concatenate((us, b[:, 0])), np.concatenate((ue, b[:, 1])))


def contig_union(bam, contig: str, *, threads: int = 2):
    """(starts, ends, n_primary): the union of the aligned blocks of the primary alignments on one contig."""
    import pysam

    us = ue = np.empty(0, dtype=np.int64)
    flat, n = array("q"), 0
    chain = itertools.chain.from_iterable
    with pysam.AlignmentFile(str(bam), threads=threads) as fh:
        for r in fh.fetch(contig):
            if r.flag & EXCLUDE_FLAGS:
                continue
            n += 1
            flat.extend(chain(r.get_blocks()))
            if len(flat) >= COMPACT_AT:
                us, ue = _fold(us, ue, flat)
                flat = array("q")
    us, ue = _fold(us, ue, flat)
    return us, ue, n


def _fp(path) -> str:
    p = Path(path)
    try:
        st = p.stat()
        return f"{p.resolve()}|{st.st_size}|{st.st_mtime_ns}"
    except OSError:
        return f"{p}|absent"


def sample_dir(cfg: dict, key: str) -> Path:
    import samples
    d = figlib.work_dir(cfg, "lrc") / samples.resolve(cfg, key)
    d.mkdir(parents=True, exist_ok=True)
    return d


def union_key(cfg: dict, key: str) -> str:
    import samples
    return f"{_fp(samples.get(cfg, key)['bam'])}\t{READ_FILTER}"


def union_contigs(cfg: dict, key: str) -> list[str]:
    """The sample's evaluation contigs (assembly.evaluation_contigs) in BAM header order."""
    import pysam
    import assembly
    import samples
    ev = assembly.evaluation_contigs(cfg, key)
    with pysam.AlignmentFile(str(samples.get(cfg, key)["bam"])) as fh:
        return [c for c in fh.references if ev is None or c in ev]


def read_manifest(cfg: dict, key: str) -> dict:
    """{contig: (n_primary, n_intervals, bases)} of the finished contigs; {} when the union's key changed."""
    f = sample_dir(cfg, key) / "union" / "manifest.tsv"
    if not f.exists():
        return {}
    lines = f.read_text().splitlines()
    if not lines or lines[0] != "# " + union_key(cfg, key):
        return {}
    return {c: (int(a), int(b), int(x)) for c, a, b, x in (l.split("\t") for l in lines[1:] if l)}


def build_union(cfg: dict, key: str, *, budget_s: float = 0.0, threads: int = 2) -> dict:
    """The union on every evaluation contig, one contig at a time (cached, resumable; read_manifest's dict). With
    budget_s > 0 it stops before the next contig once that many seconds are used (assembly.Pending, exit 75)."""
    import pysam
    import assembly
    import samples

    d = sample_dir(cfg, key) / "union"
    d.mkdir(parents=True, exist_ok=True)
    man = d / "manifest.tsv"
    done = read_manifest(cfg, key)
    if not done:
        man.write_text("# " + union_key(cfg, key) + "\n")
    bam = samples.get(cfg, key)["bam"]
    with pysam.AlignmentFile(str(bam)) as fh:
        mapped = {s.contig: s.mapped for s in fh.get_index_statistics()}
    t0 = time.time()
    for c in union_contigs(cfg, key):
        if c in done:
            continue
        if not mapped.get(c):            # no record at all: an empty union, nothing to store
            us = ue = np.empty(0, dtype=np.int64)
            n = 0
        else:
            if budget_s > 0 and time.time() - t0 > budget_s:
                raise assembly.Pending(f"lrc union {samples.resolve(cfg, key)}: {budget_s:.0f} s used; next contig {c}")
            us, ue, n = contig_union(bam, c, threads=threads)
            tmp = d / f"{c}.npz.tmp"
            with open(tmp, "wb") as fo:
                np.savez(fo, starts=us, ends=ue)
            tmp.replace(d / f"{c}.npz")
        done[c] = (n, int(us.size), int((ue - us).sum()))
        with open(man, "a") as fo:
            fo.write(f"{c}\t{n}\t{us.size}\t{done[c][2]}\n")
    return done


def load_union(cfg: dict, key: str, contig: str, done: dict):
    """(starts, ends) of one contig's union; a contig outside the finished union raises."""
    if contig not in done:
        raise RuntimeError(f"lrc union of {key} lacks {contig}: run `python3 figures/_lrc.py union --sample {key}`")
    if done[contig][1] == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    z = np.load(sample_dir(cfg, key) / "union" / f"{contig}.npz")
    return z["starts"], z["ends"]


# ---------------------------------------------------------------- models (GTF exon records)
_TID = 'transcript_id "'


def read_models(gtf, contigs: set | None = None) -> dict:
    """{contig: (ids, tidx, starts, ends)} of a GTF's exon records (0-based half-open), only on `contigs` when given;
    records without a transcript_id are skipped. Gzip-aware."""
    per: dict = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon" or (contigs is not None and f[0] not in contigs):
                continue
            i = f[8].find(_TID)
            if i < 0:
                continue
            j = f[8].find('"', i + len(_TID))
            tid = f[8][i + len(_TID):j]
            slot = per.get(f[0])
            if slot is None:
                slot = per[f[0]] = ({}, [], array("q"), array("q"), array("q"))
            index, ids, tix, st, en = slot
            k = index.get(tid)
            if k is None:
                k = index[tid] = len(ids)
                ids.append(tid)
            tix.append(k)
            st.append(int(f[3]) - 1)
            en.append(int(f[4]))
    return {c: (ids, np.frombuffer(tix, dtype=np.int64), np.frombuffer(st, dtype=np.int64),
                np.frombuffer(en, dtype=np.int64)) for c, (_, ids, tix, st, en) in per.items()}


def score_models(models: dict, union_of) -> dict:
    """{contig: (ids, n_exons, length, covered)}; union_of(contig) -> (ustarts, uends)."""
    out = {}
    for c, (ids, tix, st, en) in models.items():
        us, ue = union_of(c)
        out[c] = (ids,) + model_lrc(tix, st, en, len(ids), us, ue)
    return out


# ---------------------------------------------------------------- per sample x arm (cached)
def arms(cfg: dict, key: str, only: list[str] | None = None) -> list[str]:
    """Rustle's two arms on every sample, plus the lab baselines the registry has for it."""
    import samples
    have = ["rustle", "rustle_primary"] + [t for t in samples.BASELINES if samples.baseline(cfg, key, t)]
    return [t for t in ARMS if t in have and (not only or t in only)]


def score_arm(cfg: dict, key: str, tool: str, *, force: bool = False) -> Path:
    """<tool>.lrc.tsv / .npz of one arm on the sample's evaluation contigs (the union must be complete)."""
    import assembly

    d = sample_dir(cfg, key)
    tsv, npz = d / f"{tool}.lrc.tsv", d / f"{tool}.lrc.npz"
    src = assembly.arm_source(cfg, key, tool, check=False)
    contigs = union_contigs(cfg, key)
    done = read_manifest(cfg, key)
    missing = [c for c in contigs if c not in done]
    if missing:
        raise RuntimeError(f"lrc union of {key} is not finished ({len(missing)} contigs left): run "
                           f"`python3 figures/_lrc.py union --sample {key}` until it exits 0")
    k = f"{union_key(cfg, key)}\t{_fp(src)}\t{len(contigs)} contigs"
    side = Path(str(npz) + ".key")
    if not force and npz.exists() and tsv.exists() and side.exists() and side.read_text().rstrip("\n") == k:
        return npz
    scored = score_models(read_models(src, set(contigs)), lambda c: load_union(cfg, key, c, done))
    tmp = Path(str(tsv) + ".tmp")
    cols = {"n_exons": [], "length": [], "covered": []}
    with open(tmp, "w") as fo:
        fo.write("contig\ttranscript_id\tn_exons\tlength\tcovered\tlrc\n")
        for c in (c for c in contigs if c in scored):
            ids, nx, L, cov = scored[c]
            for i, tid in enumerate(ids):
                fo.write(f"{c}\t{tid}\t{nx[i]}\t{L[i]}\t{cov[i]}\t{cov[i] / L[i]:.6f}\n")
            cols["n_exons"].append(nx)
            cols["length"].append(L)
            cols["covered"].append(cov)
    tmp.replace(tsv)
    cat = {k2: (np.concatenate(v) if v else np.empty(0, dtype=np.int64)) for k2, v in cols.items()}
    tmpz = Path(str(npz) + ".tmp")
    with open(tmpz, "wb") as fo:
        np.savez(fo, **cat)
    tmpz.replace(npz)
    side.write_text(k + "\n")
    return npz


def summary_rows(cfg: dict, key: str, tool: str) -> list[list]:
    import assembly
    import samples
    row = samples.get(cfg, key)
    z = np.load(sample_dir(cfg, key) / f"{tool}.lrc.npz")
    lead = [row["species"], row["id"], assembly.scope_label(cfg, key), tool]
    multi = z["n_exons"] > 1
    return [lead + ["all"] + summarize(z["length"], z["covered"]),
            lead + ["multi-exon"] + summarize(z["length"][multi], z["covered"][multi])]


# ---------------------------------------------------------------- per SQANTI3 category (fig. 2 cache, where present)
def sqanti_rows(cfg: dict, key: str, tool: str) -> tuple[list[list], dict]:
    """Rows of lrc_sqanti for one sample x arm, and {'classified', 'joined', 'stale', 'contigs'} for the notes."""
    import assembly
    import samples
    import _sqanti as S

    skey = assembly.sample_key(cfg, key)
    root = figlib.work_dir(cfg, "fig2") / skey
    href = figlib.work_dir(cfg, "fig2") / S.reference_home(cfg, key)
    done = read_manifest(cfg, key)
    part = re.compile(rf"^qc_{re.escape(tool)}(\.part\d+)?$")
    by_cat: dict = {}
    info = {"classified": 0, "joined": 0, "stale": 0, "contigs": []}
    for cdir in sorted(p for p in root.iterdir() if p.is_dir()) if root.exists() else []:
        c = cdir.name
        gtf = cdir / f"{tool}.gtf"
        qdirs = sorted(q for q in cdir.iterdir() if q.is_dir() and part.match(q.name))
        if not gtf.exists() or not qdirs or c not in done:
            continue
        cats: dict = {}
        for q in qdirs:
            sfx = q.name[len(f"qc_{tool}"):]
            cls = q / f"{tool}_classification.txt"
            if not cls.exists():
                continue
            n_cls = sum(1 for _ in open(cls)) - 1
            if not S.qc_cached(cfg, cdir / f"{tool}{sfx}.gtf", href / c / "ref.gtf", href / c / "genome.fa", q, tool):
                info["stale"] += n_cls
                continue
            for r in S._rows(cls):
                cats[r["isoform"]] = S.fold(r["structural_category"])
        if not cats:
            continue
        info["contigs"].append(c)
        info["classified"] += len(cats)
        ids, nx, L, cov = score_models(read_models(gtf, {c}), lambda cc: load_union(cfg, key, cc, done)).get(
            c, ([], [], [], []))
        for i, tid in enumerate(ids):
            cat = cats.get(tid)
            if cat is None:
                continue
            info["joined"] += 1
            cell = by_cat.setdefault(cat, ([], []))
            cell[0].append(int(L[i]))
            cell[1].append(int(cov[i]))
    row = samples.get(cfg, key)
    scope = ",".join(info["contigs"])
    rows = [[row["species"], row["id"], scope, tool, cat] + summarize(*by_cat[cat])[:6]
            for cat in figlib.SQANTI_ORDER if cat in by_cat]
    return rows, info


# ---------------------------------------------------------------- CLI
def cmd_tables(cfg: dict, keys: list[str], data_dir: Path):
    import assembly
    import samples
    summ, sq, notes_sq, inputs = [], [], [], {}
    for key in keys:
        sid = samples.resolve(cfg, key)
        for tool in arms(cfg, sid):
            npz = sample_dir(cfg, sid) / f"{tool}.lrc.npz"
            if not npz.exists():
                print(f"[lrc] {sid} {tool}: not scored (python3 figures/_lrc.py score --sample {sid})", file=sys.stderr)
                continue
            inputs[f"{sid} {tool} gtf"] = assembly.arm_source(cfg, sid, tool, check=False)
            summ += summary_rows(cfg, sid, tool)
            rows, info = sqanti_rows(cfg, sid, tool)
            sq += rows
            if info["classified"] or info["stale"]:
                notes_sq.append(f"{sid} {tool}: {info['joined']} of {info['classified']} classified isoforms joined "
                                f"on {','.join(info['contigs'])}; {info['stale']} skipped (classification older than "
                                "its input GTF)")
        inputs[f"{sid} bam"] = samples.get(cfg, sid)["bam"]
    base = ["%LRC (LRGASP, Pardo-Palacios et al. 2024, Box 1): fraction of the transcript model sequence length mapped "
            "by one or more long reads; docs/PREREG_lrc_metric_2026-09-29.md", "reads: " + READ_FILTER,
            "model = GTF transcript, exons merged; classes > 0.98, 0.75-0.98 (inclusive), < 0.75 (LRGASP Extended "
            "Data Fig. 2); reporting only (no filter, no ranking); samples and species never pooled"]
    gen = "figures/_lrc.py tables"
    figlib.write_table("lrc_summary", SUMMARY_HEADER, summ, generator=gen, inputs=inputs, data_dir=data_dir,
                       notes=base + ["scope: each sample's evaluation contigs (assembly.evaluation_contigs); "
                                     "multi-exon = more than one exon record"])
    figlib.write_table("lrc_sqanti", SQANTI_HEADER, sq, generator=gen, inputs=inputs, data_dir=data_dir,
                       notes=base + ["the per-contig GTFs SQANTI3 classified for fig. 2 (development scope), joined "
                                     "by isoform id"] + notes_sq)
    return data_dir


def main(argv=None):
    import argparse
    import samples
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("union", "score", "tables"):
        p = sub.add_parser(name)
        p.add_argument("--inputs")
        p.add_argument("--sample" if name != "tables" else "--samples", required=name != "tables")
        if name == "union":
            p.add_argument("--budget-s", type=float, default=150.0)
            p.add_argument("--threads", type=int, default=2)
        if name == "score":
            p.add_argument("--tools", default="")
            p.add_argument("--force", action="store_true")
        if name == "tables":
            p.add_argument("--data")
    a = ap.parse_args(argv)
    cfg = figlib.load_inputs(a.inputs)
    if a.cmd == "union":
        t0 = time.time()
        done = build_union(cfg, a.sample, budget_s=a.budget_s, threads=a.threads)
        print(f"[lrc] union {samples.resolve(cfg, a.sample)}: {len(done)} contigs, "
              f"{sum(v[0] for v in done.values())} primary alignments, {sum(v[2] for v in done.values())} bases "
              f"covered ({time.time() - t0:.0f} s this call)")
        return
    if a.cmd == "score":
        only = [t.strip() for t in a.tools.split(",") if t.strip()]
        for tool in arms(cfg, a.sample, only):
            t0 = time.time()
            print(f"[lrc] {tool}: {score_arm(cfg, a.sample, tool, force=a.force)} ({time.time() - t0:.0f} s)")
        return
    keys = [k.strip() for k in (a.samples or "").split(",") if k.strip()] or list(samples.registry(cfg))
    out = Path(a.data) if a.data else figlib.work_dir(cfg, "lrc") / "tables"
    print(f"[lrc] tables in {cmd_tables(cfg, keys, out)}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
