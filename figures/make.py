#!/usr/bin/env python3
"""make.py — build and render the publication figures.

    python3 figures/make.py list                       what each figure shows and which tables it reads
    python3 figures/make.py data  FIG|all [--force]    regenerate figures/data/*.tsv from the raw inputs (HEAVY:
                                                        runs the pipeline, gffcompare, SQANTI3, simulations;
                                                        foreground, one heavy process at a time)
                          [--recorded]                  figures that support it tabulate the recorded development
                                                        runs instead (light; tables marked provisional; written to
                                                        figures/data_recorded/, never over figures/data/)
                          [--set KEY=VALUE ...]         override one input key for this call (e.g. fig2_arms=...)
    python3 figures/make.py plot  FIG|all [--data DIR] render figures/out/<fig>.{pdf,png,svg} from the tables
    python3 figures/make.py check                      every table present, with provenance; every figure renders
    python3 figures/make.py samples [--verify] [--sample ID]   the sample registry (figures/samples.tsv), resolved;
                                                        --verify checks files, indexes and contig names (light)
    python3 figures/make.py runs [--sample ID|all] [--stage STAGE|all] [--dry-run] [--max-stages N] [--force]
                                                        GENOME-WIDE pipeline runs per sample, cached under
                                                        ${work}/runs/<sample>/ (HEAVY; foreground; one stage at a time;
                                                        --dry-run prints the queue with estimated cost). Stages:
                                                        assemble assemble_primary families families_primary catalog
                                                        assign index flag.
                                                        Exit 75 = a bounded call (catalog pieces, a sharded all-vs-all)
                                                        or --max-stages left work: run the same command again
    python3 figures/make.py runs --sample ID --stage STAGE --restamp PROOF.tsv
                                                        re-stamp a STALE stage without running it, only from a proof of
                                                        cmp-identical products made by the current code (samples.restamp;
                                                        the proof is recorded in the stamp)

Machine-local paths come from `figures/inputs.local.tsv` (copy `inputs.example.tsv`), or `--inputs FILE`.
Each figure is a module `figures/fig_*.py` exposing:

    META  = {"id": "fig1", "title": ..., "claim": ..., "tables": [...]}
    build(cfg, data_dir, force)   -> writes its tables (figlib.write_table), caching heavy steps in work_dir
    plot(data_dir, out_dir)       -> renders the figure (figlib.use_style + figlib.save)
"""
from __future__ import annotations

import argparse
import importlib
import inspect
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figlib  # noqa: E402

RECORDED_DIR = HERE / "data_recorded"

# registration order = figure numbering in the paper
MODULES = [
    "fig_intron_chain",      # 1  gffcompare intron-chain SN/PR vs StringTie / FLAIR / IsoSeq
    "fig_sqanti",            # 2  SQANTI3 structural categories + rules filter
    "fig_secondary",         # 3  transcripts built with secondary alignments in multi-mapping loci
    "fig_assignability",     # 4  UpSet: which reads can be assigned to a copy, by copy identity
    "fig_assign_accuracy",   # 5  copy assignment accuracy / coverage by divergence, vs aligner and tools
    "fig_family_spectrum",   # 6  family definition across the identity spectrum (O1)
    "fig_family_recovery",   # 7  family recovery, de novo vs guided: sensitivity, precision, bipartite matching (O1)
    "fig_loci",              # 8  loci in the Liftoff framework: self-lift baseline, guided search, de novo loci
]


def _commented_keys(path) -> set:
    """Optional keys documented as `# key<TAB>...` comment lines of the inputs file."""
    out = set()
    for line in open(path):
        if line.startswith("# ") and "\t" in line:
            out.add(line[2:].split("\t", 1)[0].strip())
    return out


def _pdf_width_mm(path):
    import re
    m = re.search(rb"/MediaBox \[\s*[\d.]+ [\d.]+ ([\d.]+) ([\d.]+)", Path(path).read_bytes())
    return float(m.group(1)) / 72 * 25.4 if m else None


def modules():
    out = []
    for name in MODULES:
        try:
            out.append(importlib.import_module(name))
        except ModuleNotFoundError as e:
            if e.name == name:
                print(f"[make] {name}: not implemented yet", file=sys.stderr)
                continue
            raise
    return out


def select(which: str):
    mods = modules()
    if which == "all":
        return mods
    sel = [m for m in mods if m.META["id"] == which or m.__name__ == which]
    if not sel:
        sys.exit(f"unknown figure {which!r}; try: {', '.join(m.META['id'] for m in mods)}")
    return sel


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    d = sub.add_parser("data")
    d.add_argument("fig")
    d.add_argument("--force", action="store_true", help="recompute cached intermediates too")
    d.add_argument("--inputs", help="inputs TSV (default figures/inputs.local.tsv)")
    d.add_argument("--recorded", action="store_true",
                   help="tabulate the recorded development runs (only figures whose build() takes recorded=)")
    d.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override an input key")
    p = sub.add_parser("plot")
    p.add_argument("fig")
    p.add_argument("--out", default=str(figlib.OUT_DIR))
    p.add_argument("--data", default=str(figlib.DATA_DIR), help="tables to render (e.g. figures/data_recorded)")
    sub.add_parser("check")
    sm = sub.add_parser("samples")
    sm.add_argument("--verify", action="store_true", help="check files, indexes, contig names (light; seconds-minutes)")
    sm.add_argument("--sample", help="one sample id or alias (default: all)")
    sm.add_argument("--inputs", help="inputs TSV (default figures/inputs.local.tsv)")
    r = sub.add_parser("runs")
    r.add_argument("--sample", default="all", help="sample id, alias, or all")
    r.add_argument("--stage", default="all",
                   help="assemble|assemble_primary|families|families_primary|catalog|assign|index|flag|all")
    r.add_argument("--dry-run", action="store_true", help="print the queue with estimated wall time and peak RSS")
    r.add_argument("--max-stages", type=int, help="run at most N stages in this call (re-run to continue)")
    r.add_argument("--force", action="store_true", help="re-run the selected stages even when fresh")
    r.add_argument("--ignore-busy", action="store_true", help="run even if another heavy process is running")
    r.add_argument("--restamp", metavar="PROOF.tsv",
                   help="re-stamp one stale --sample/--stage WITHOUT running it, from a proof file listing every "
                        "product as cmp-identical to a copy made by the current code (see samples.read_proof)")
    r.add_argument("--inputs", help="inputs TSV (default figures/inputs.local.tsv)")
    a = ap.parse_args(argv)

    if a.cmd in ("samples", "runs"):
        import samples
        cfg = figlib.load_inputs(a.inputs)
        if a.cmd == "samples":
            sys.exit(samples.cli_samples(cfg, a.verify, a.sample))
        if a.restamp:
            sys.exit(samples.cli_restamp(cfg, a.sample, a.stage, a.restamp))
        sys.exit(samples.cli_runs(cfg, a.sample, a.stage, a.dry_run, a.max_stages, a.force, a.ignore_busy))

    if a.cmd == "list":
        for m in modules():
            M = m.META
            print(f"{M['id']}  {M['title']}\n      claim:  {M['claim']}\n      tables: {', '.join(M['tables'])}")
        return
    if a.cmd == "data":
        cfg = figlib.load_inputs(a.inputs)
        example = figlib.load_inputs(figlib.HERE / "inputs.example.tsv")
        documented = set(example) | _commented_keys(figlib.HERE / "inputs.example.tsv")
        for kv in a.set:
            k, sep, v = kv.partition("=")
            if not sep:
                sys.exit(f"--set expects KEY=VALUE, got {kv!r}")
            # a documented key `name_<placeholder>` stands for every key starting with `name_`
            prefixes = tuple(d.split("<", 1)[0] for d in documented if "<" in d)
            if k not in cfg and k not in documented and not (prefixes and k.startswith(prefixes)):
                print(f"[make] WARNING: --set {k}: not a key of inputs.example.tsv (typo?)", file=sys.stderr)
            for kk, vv in list(cfg.items()):
                v = v.replace("${" + kk + "}", vv)
            cfg[k] = v
        for m in select(a.fig):
            kw = {}
            if a.recorded:
                if "recorded" not in inspect.signature(m.build).parameters:
                    print(f"[make] {m.META['id']}: no recorded mode; skipped", file=sys.stderr)
                    continue
                kw["recorded"] = True
            # recorded (development) tables never overwrite the committed ones: they go to data_recorded/
            out = RECORDED_DIR if kw else figlib.DATA_DIR
            out.mkdir(parents=True, exist_ok=True)
            print(f"[make] data {m.META['id']} ({m.__name__}){' --recorded -> ' + str(out) if kw else ''}",
                  file=sys.stderr)
            m.build(cfg, out, a.force, **kw)
        return
    if a.cmd == "plot":
        figlib.use_style()
        for m in select(a.fig):
            print(f"[make] plot {m.META['id']} ({m.__name__})", file=sys.stderr)
            m.plot(Path(a.data), Path(a.out))
        return
    if a.cmd == "check":
        ok = True
        figlib.use_style()
        tmp_out = Path("/tmp/rustle_figures_check")
        for m in modules():
            supp = list(m.META.get("supplementary_tables", []))
            for t in list(m.META["tables"]) + supp:
                path = figlib.DATA_DIR / f"{t}.tsv"
                if not path.exists():
                    if t in supp:   # drawn only once its build writes it (e.g. the guided path, the all-sample tables)
                        print(f"absent   {m.META['id']}: {t} (supplementary; appears with make.py data {m.META['id']})")
                        continue
                    print(f"MISSING  {m.META['id']}: {path}")
                    ok = False
                    continue
                meta = figlib.table_meta(t)
                commit = meta.get("commit", ["?"])[0]
                n_rows = len(figlib.read_table(t))
                flags = []
                if n_rows == 0:
                    flags.append("EMPTY")
                    ok = False
                if commit.endswith("-dirty") or commit in ("unknown", "?"):
                    flags.append("uncommitted generator")
                if figlib.provisional_notes([t]):
                    flags.append("PROVISIONAL")
                print(f"ok       {m.META['id']}: {t} ({n_rows} rows; commit {commit}, {meta.get('date', ['?'])[0]})"
                      + (f"  [{'; '.join(flags)}]" if flags else ""))
            try:
                paths = m.plot(figlib.DATA_DIR, tmp_out) or []
            except Exception as e:  # noqa: BLE001
                print(f"FAILS    {m.META['id']}: {e!r}")
                ok = False
                continue
            print(f"renders  {m.META['id']}")
            for p in paths:
                if p.suffix != ".pdf":
                    continue
                w = _pdf_width_mm(p)
                if w is not None and w > 183.0:
                    print(f"WIDE     {m.META['id']}: {p.name} is {w:.1f} mm (> 183 mm double column)")
                    ok = False
                shipped = figlib.OUT_DIR / p.name
                if not shipped.exists() or shipped.read_bytes() != p.read_bytes():
                    print(f"STALE    {m.META['id']}: {shipped} differs from a fresh render (run make.py plot {m.META['id']})")
                    ok = False
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
