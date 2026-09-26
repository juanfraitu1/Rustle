"""figlib — shared infrastructure for the publication figures (style, palette, tidy tables, UpSet, saving).

Every figure module in this directory uses these pieces so the set reads as one system:

* **palette** — one colour per tool, fixed order, validated colour-vision-deficiency safe on ALL pairs
  (dataviz validator, 2026-09-25: worst all-pairs CVD ΔE 9.2, normal-vision 16.3; FLAIR aqua is below 3:1 on
  white, so tools are always direct-labelled or legended, never colour-alone). Rustle primaries-only (no
  secondary-alignment seeding) is the Rustle hue with a hatch / open marker, not a new hue.
* **tidy tables** — every plotted number lives in `figures/data/<table>.tsv` with a `#`-comment provenance
  header (generator, git commit, date, inputs with size+mtime). `make.py plot` needs only these tables, so the
  figures re-render anywhere; `make.py data` regenerates the tables from the raw inputs.
* **species are never pooled** — a panel shows one species; the column `species` is mandatory in every table
  whose numbers come from reads.
"""
from __future__ import annotations

import csv
import datetime as _dt
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA_DIR = HERE / "data"
OUT_DIR = HERE / "out"

# ---------------------------------------------------------------- palette (validated; fixed order, never cycled)
TOOL_ORDER = ["rustle", "rustle_primary", "stringtie", "flair", "isoseq"]
TOOL_LABEL = {
    "rustle": "Rustle",
    "rustle_primary": "Rustle (primaries only)",
    "stringtie": "StringTie",
    "flair": "FLAIR",
    "isoseq": "IsoSeq collapse",
}
TOOL_COLOR = {
    "rustle": "#2a78d6",
    "rustle_primary": "#2a78d6",  # same hue; distinguished by HATCH / open marker (secondary encoding)
    "stringtie": "#eb6834",
    "flair": "#1baf7a",
    "isoseq": "#4a3aa7",
}
TOOL_HATCH = {"rustle_primary": "////"}
TOOL_MARKER = {"rustle": "o", "rustle_primary": "o", "stringtie": "s", "flair": "^", "isoseq": "D"}
TOOL_MARKER_FILLED = {"rustle_primary": False}
# ANNOTATION-GUIDED runs of the same tools (the separate guided path of figures 1-3; never drawn in a panel with the
# annotation-free methods above): the tool's hue, open marker / cross-hatch, and "(guided)" in every label
GUIDED_TOOL_ORDER = ["stringtie_guided", "flair_guided"]
TOOL_LABEL.update({"stringtie_guided": "StringTie (guided, -G)", "flair_guided": "FLAIR (guided, with annotation)"})
TOOL_COLOR.update({"stringtie_guided": TOOL_COLOR["stringtie"], "flair_guided": TOOL_COLOR["flair"]})
TOOL_HATCH.update({"stringtie_guided": "xxxx", "flair_guided": "xxxx"})
TOOL_MARKER.update({"stringtie_guided": "s", "flair_guided": "^"})
TOOL_MARKER_FILLED.update({"stringtie_guided": False, "flair_guided": False})
# on-figure mode lines (the wording of assembly.MODE_DENOVO_METHODS / GUIDED_NA; plot() reads no inputs file)
MODE_DENOVO_LINE = ("Annotation-free (de novo) comparison: Rustle assemble (reads + genome), StringTie -L without -G, "
                    "FLAIR collapse without annotation (flair correct skipped), IsoSeq collapse.")
GUIDED_NA_LINE = "Guided comparison: not available (guided StringTie/FLAIR GTFs not supplied)."
GUIDED_SEPARATE_LINE = "Guided comparison: separate figure ({name}), never drawn with the annotation-free methods."


def mode_lines(guided_figure: str, guided_present: bool) -> str:
    """The two-line mode statement every tool-comparison figure prints (figures 1-3)."""
    return MODE_DENOVO_LINE + "\n" + (GUIDED_SEPARATE_LINE.format(name=guided_figure) if guided_present
                                      else GUIDED_NA_LINE)

SPECIES_LABEL = {"gorilla": "Gorilla (testis, OR6737)", "human": "Human (A119b)"}

# ink (text wears these, never a series colour)
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8983"
GRID = "#e4e3de"
SURFACE = "#ffffff"

# sequential blue ramp (reference instance; ordinal use starts no lighter than step 250)
BLUE = {100: "#cde2fb", 250: "#86b6ef", 300: "#6da7ec", 350: "#5598e7", 400: "#3987e5", 450: "#2a78d6",
        500: "#256abf", 550: "#1c5cab", 600: "#184f95", 650: "#104281", 700: "#0d366b"}

# SQANTI3 structural categories: an ORDINAL ramp from "matches the annotation" to "novel", plus neutral
# "artifact-prone" categories. Tools are identified by x position + label in that figure, not colour.
SQANTI_ORDER = ["FSM", "ISM", "NIC", "NNC", "Other"]
SQANTI_LABEL = {
    "FSM": "Full splice match",
    "ISM": "Incomplete splice match",
    "NIC": "Novel in catalog",
    "NNC": "Novel not in catalog",
    "Other": "Antisense / intergenic / genic / genic intron / fusion",
}
SQANTI_COLOR = {"FSM": BLUE[650], "ISM": BLUE[450], "NIC": BLUE[300], "NNC": BLUE[250], "Other": "#bdbcb6"}
SQANTI_CATEGORY_MAP = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
}  # everything else -> "Other"

# copy-identity bands for read assignability (closest-sibling identity of the read's source copy)
IDENTITY_BANDS = ["identical", "99.5-100%", "99-99.5%", "98-99%", "<98%"]
IDENTITY_COLOR = {"identical": INK_2, "99.5-100%": BLUE[650], "99-99.5%": BLUE[500], "98-99%": BLUE[350],
                  "<98%": BLUE[250]}


def identity_band(identity: float | None) -> str | None:
    """Band of a copy's identity to its closest sibling (a fraction in [0, 1]); None if unknown."""
    if identity is None or identity != identity:  # NaN
        return None
    if identity >= 1.0:
        return "identical"
    if identity >= 0.995:
        return "99.5-100%"
    if identity >= 0.99:
        return "99-99.5%"
    if identity >= 0.98:
        return "98-99%"
    return "<98%"


# ---------------------------------------------------------------- matplotlib style
MM = 1 / 25.4
WIDTH_SINGLE = 89 * MM   # single column (inches)
WIDTH_DOUBLE = 183 * MM  # double column


# Arial is not shipped with matplotlib; register it from the usual places so every machine that has it lays the text
# out with the same metrics (the SVGs name Arial, so a DejaVu layout would shift when viewed on a machine with Arial).
FONT_FILES = ["arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"]
FONT_DIRS = [HERE / "fonts", Path("/mnt/c/Windows/Fonts"), Path("C:/Windows/Fonts"), Path("/usr/share/fonts/truetype/msttcorefonts"),
             Path("/Library/Fonts"), Path("/System/Library/Fonts/Supplemental")]


def _register_fonts():
    from matplotlib import font_manager
    for d in FONT_DIRS:
        found = [d / f for f in FONT_FILES if (d / f).exists()]
        if found:
            for f in found:
                font_manager.fontManager.addfont(str(f))
            return str(d)
    return None


def use_style():
    import matplotlib as mpl

    mpl.use("Agg")
    _register_fonts()
    mpl.rcParams.update({
        "svg.hashsalt": "rustle-figures",  # deterministic SVG ids: unchanged data re-renders byte-identical files
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.edgecolor": INK_2,
        "axes.labelcolor": INK,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "text.color": INK,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "lines.linewidth": 1.5,
        "lines.markersize": 4.5,
        "legend.frameon": False,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "pdf.fonttype": 42,   # editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "hatch.linewidth": 0.6,
        "hatch.color": "#ffffff",
    })


def panel_label(ax, letter: str, x: float = -0.18, y: float = 1.06):
    """Bold lower-case panel letter at the top-left of an axes (Nature style)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def tool_bar_kwargs(tool: str) -> dict:
    kw = {"color": TOOL_COLOR[tool], "edgecolor": SURFACE, "linewidth": 0.8}
    if tool in TOOL_HATCH:
        kw.update(hatch=TOOL_HATCH[tool], color=TOOL_COLOR[tool], alpha=0.55)
    return kw


def tool_marker_kwargs(tool: str) -> dict:
    filled = TOOL_MARKER_FILLED.get(tool, True)
    c = TOOL_COLOR[tool]
    return {
        "marker": TOOL_MARKER[tool],
        "markerfacecolor": c if filled else SURFACE,
        "markeredgecolor": c,
        "markeredgewidth": 1.2,
        "color": c,
        "linestyle": "none",
    }


def save(fig, name: str, out_dir: Path | None = None, formats=("pdf", "png", "svg")) -> list[Path]:
    out_dir = Path(out_dir or OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        p = out_dir / f"{name}.{fmt}"
        # write a sibling temp file, then rename: a file held open by a Windows viewer (/mnt/c) cannot be
        # truncated in place (EINVAL) but can be replaced
        tmp = out_dir / f".{name}.{fmt}.tmp"
        fig.savefig(tmp, format=fmt, dpi=300 if fmt == "png" else None, bbox_inches="tight", pad_inches=0.02,
                    metadata=_save_metadata(fmt))
        os.replace(tmp, p)
        paths.append(p)
    return paths


def _save_metadata(fmt):
    # deterministic outputs: no creation date in PDF/SVG, so re-rendering unchanged data gives identical files
    if fmt == "pdf":
        return {"CreationDate": None, "Creator": "Rustle figures", "Producer": "matplotlib"}
    if fmt == "svg":
        return {"Date": None, "Creator": "Rustle figures"}
    return {}


# ---------------------------------------------------------------- provenance + tidy tables
def git_commit() -> str:
    try:
        sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "src", "bench", "figures", "tools", "Cargo.toml", "Cargo.lock"],
                               capture_output=True, text=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def file_fingerprint(path) -> str:
    p = Path(path)
    try:
        st = p.stat()
        return f"{p.resolve()}\tsize={st.st_size}\tmtime={int(st.st_mtime)}"
    except OSError:
        return f"{p}\tabsent"


def write_table(name: str, header: list[str], rows, *, generator: str, inputs: dict | None = None,
                notes: list[str] | None = None, data_dir: Path | None = None) -> Path:
    """Write `figures/data/<name>.tsv` with a provenance header. Rows are sequences aligned to `header`."""
    data_dir = Path(data_dir or DATA_DIR)
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / f"{name}.tsv"
    tmp = path.with_suffix(".tsv.tmp")  # written whole, then renamed: an interrupted build never leaves half a table
    with open(tmp, "w", newline="") as fh:
        fh.write(f"# table: {name}\n# generator: {generator}\n# commit: {git_commit()}\n")
        fh.write(f"# date: {_dt.date.today().isoformat()}\n")
        for k, v in (inputs or {}).items():
            fh.write(f"# input {k}: {file_fingerprint(v) if v else 'none'}\n")
        for n in notes or []:
            fh.write(f"# note: {n}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if x is None else (f"{x:.6g}" if isinstance(x, float) else x) for x in r])
    tmp.replace(path)
    return path


def read_table(name: str, data_dir: Path | None = None) -> list[dict]:
    """Rows of `figures/data/<name>.tsv` as dicts (strings); `#` lines are provenance and skipped."""
    path = Path(data_dir or DATA_DIR) / f"{name}.tsv"
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def table_meta(name: str, data_dir: Path | None = None) -> dict:
    path = Path(data_dir or DATA_DIR) / f"{name}.tsv"
    meta = {}
    with open(path) as fh:
        for l in fh:
            if not l.startswith("#"):
                break
            k, _, v = l[1:].strip().partition(": ")
            meta.setdefault(k, []).append(v)
    return meta


def provisional_notes(tables, data_dir: Path | None = None) -> list[str]:
    """The `# note: provisional...` lines of the given tables (absent tables are skipped)."""
    out = []
    for t in tables:
        try:
            out += [n for n in table_meta(t, data_dir).get("note", []) if n.lower().startswith("provisional")]
        except FileNotFoundError:
            pass
    return out


def stamp_provisional(fig, tables, data_dir: Path | None = None) -> bool:
    """Print one 'PROVISIONAL' line at the bottom-right of `fig` when any of `tables` carries a provisional note
    (a figure built from recorded development runs must say so on its face, not only in the caption)."""
    if not provisional_notes(tables, data_dir):
        return False
    fig.text(0.995, 0.002, "PROVISIONAL tables (development runs or subsets; see caption)", fontsize=5.5,
             color=INK_3, ha="right", va="bottom")
    return True


# ---------------------------------------------------------------- machine-local inputs
def load_inputs(path: Path | None = None) -> dict:
    """`figures/inputs.local.tsv` (key<TAB>value; `${key}` expands earlier keys and environment variables).

    Falls back to `figures/inputs.example.tsv`, which documents every key the data builders read."""
    path = Path(path) if path else (HERE / "inputs.local.tsv" if (HERE / "inputs.local.tsv").exists()
                                    else HERE / "inputs.example.tsv")
    cfg: dict[str, str] = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        k, _, v = line.partition("\t")
        v = v.strip()
        for kk, vv in list(cfg.items()) + list(os.environ.items()):
            v = v.replace("${" + kk + "}", vv)
        cfg[k.strip()] = v
    cfg["_inputs_file"] = str(path)
    return cfg


def work_dir(cfg: dict, fig: str) -> Path:
    """Scratch for heavy intermediates (never inside the repo): `${work}/<fig>/`."""
    d = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / fig
    d.mkdir(parents=True, exist_ok=True)
    return d


def run(cmd, *, log: Path | None = None, env: dict | None = None, cwd=None, check=True):
    """Run a command in the foreground (one heavy process at a time — the machine rule), logging to `log`."""
    e = dict(os.environ)
    e.setdefault("TMPDIR", "/mnt/linuxdisk/tmp")
    e.update(env or {})
    if log:
        with open(log, "w") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=e, cwd=cwd, shell=isinstance(cmd, str))
    else:
        r = subprocess.run(cmd, env=e, cwd=cwd, shell=isinstance(cmd, str))
    if check and r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}" + (f" — see {log}" if log else ""))
    return r.returncode


def fresh(target: Path, *sources) -> bool:
    """True if `target` exists and is newer than every existing source (cheap make-style caching)."""
    t = Path(target)
    if not t.exists():
        return False
    tm = t.stat().st_mtime
    return all(not Path(s).exists() or Path(s).stat().st_mtime <= tm for s in sources)


# ---------------------------------------------------------------- UpSet (no external dependency)
def upset(fig, membership_counts: dict, set_names: list[str], *, stack: dict | None = None,
          stack_order: list[str] | None = None, stack_colors: dict | None = None, set_labels: dict | None = None,
          min_size: int = 1, max_columns: int = 20, gridspec=None, bar_ylabel: str = "Reads",
          set_xlabel: str = "Set size", break_at: tuple[float, float] | None = None):
    """Draw an UpSet plot into `fig` (or into `gridspec`, a SubplotSpec) and return its axes dict.

    membership_counts: {frozenset(set names): count} — EXCLUSIVE intersections (each item counted once, in the
      exact combination of sets it belongs to; the empty set may be included and is drawn as "none").
    stack: optional {frozenset: {category: count}} — splits each intersection bar by a category (e.g. identity
      band); categories drawn in `stack_order` with `stack_colors`.
    Columns are sorted by size (descending); columns smaller than `min_size` or beyond `max_columns` are folded
    into one "other" column so nothing is silently dropped.
    break_at: optional (low_top, high_bottom) — a true broken y axis for one dominant column: the lower axes shows
      0..low_top, the upper axes high_bottom..max, with break marks between them (every bar keeps its real
      height; nothing is compressed). Ignored when no column exceeds low_top. The returned "bar" is the lower
      axes; "bar_top" is the upper one (or None).
    """
    import matplotlib.gridspec as gs
    import numpy as np

    items = [(k, v) for k, v in membership_counts.items() if v > 0]
    items.sort(key=lambda kv: (-kv[1], sorted(kv[0])))
    shown = [kv for kv in items if kv[1] >= min_size][:max_columns]
    folded = [kv for kv in items if kv not in shown]
    cols = [k for k, _ in shown]
    sizes = [v for _, v in shown]
    if folded:
        cols.append("__other__")
        sizes.append(sum(v for _, v in folded))
    n_sets = len(set_names)
    outer = gridspec if gridspec is not None else gs.GridSpec(1, 1, figure=fig)[0]
    g = gs.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer, height_ratios=[2.2, 0.28 * n_sets + 0.2],
                                   width_ratios=[0.9, max(3.0, 0.32 * len(cols))], hspace=0.05, wspace=0.03)
    broken = break_at is not None and max(sizes, default=0) > break_at[0]
    if broken:
        gb = gs.GridSpecFromSubplotSpec(2, 1, subplot_spec=g[0, 1], height_ratios=[1, 2.6], hspace=0.08)
        ax_top = fig.add_subplot(gb[0])
        ax_bar = fig.add_subplot(gb[1], sharex=ax_top)
    else:
        ax_top = None
        ax_bar = fig.add_subplot(g[0, 1])
    ax_mat = fig.add_subplot(g[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(g[1, 0], sharey=ax_mat)
    x = np.arange(len(cols))
    bar_axes = [ax_bar] + ([ax_top] if broken else [])
    # intersection bars (optionally stacked); drawn on both halves of a broken axis
    for ax in bar_axes:
        if stack and stack_order:
            bottom = np.zeros(len(cols))
            for cat in stack_order:
                vals = []
                for c in cols:
                    if c == "__other__":
                        vals.append(sum((stack.get(k, {}) or {}).get(cat, 0) for k, _ in folded))
                    else:
                        vals.append((stack.get(c, {}) or {}).get(cat, 0))
                vals = np.array(vals, dtype=float)
                ax.bar(x, vals, bottom=bottom, width=0.72, color=(stack_colors or {}).get(cat, INK_3),
                       edgecolor=SURFACE, linewidth=0.5, label=cat if ax is ax_bar else None)
                bottom += vals
        else:
            ax.bar(x, sizes, width=0.72, color=INK_2, edgecolor=SURFACE, linewidth=0.5)
    for xi, s in zip(x, sizes):
        ax = ax_top if broken and s > break_at[0] else ax_bar
        ax.annotate(f"{s:,}", (xi, s), xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
                    fontsize=5.5, color=INK_2, rotation=90)
    if broken:
        top = max(sizes)
        ax_bar.set_ylim(0, break_at[0])
        ax_top.set_ylim(break_at[1], top + 0.35 * (top - break_at[1]))
        ax_top.spines["bottom"].set_visible(False)
        ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_top.yaxis.set_major_locator(__import__("matplotlib.ticker", fromlist=["MaxNLocator"]).MaxNLocator(2))
        # break marks: short diagonals on the y spine of both halves
        kw = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=5, linestyle="none", color=INK_2, mec=INK_2, mew=0.8,
                  clip_on=False)
        ax_top.plot([0], [0], transform=ax_top.transAxes, **kw)
        ax_bar.plot([0], [1], transform=ax_bar.transAxes, **kw)
    else:
        ax_bar.set_ylim(0, max(sizes, default=1) * 1.18)
    ax_bar.set_ylabel(bar_ylabel)
    ax_bar.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_bar.spines["bottom"].set_visible(False)
    # matrix
    ax_mat.set_ylim(-0.5, n_sets - 0.5)
    for i in range(n_sets):
        if i % 2 == 0:
            ax_mat.axhspan(i - 0.5, i + 0.5, color="#f4f3ef", zorder=0, linewidth=0)
    for xi, c in enumerate(cols):
        if c == "__other__":
            ax_mat.text(xi, (n_sets - 1) / 2, "other", rotation=90, ha="center", va="center", fontsize=5.5,
                        color=INK_3)
            continue
        on = sorted(n_sets - 1 - set_names.index(s) for s in c if s in set_names)  # frozenset order is hash order
        ax_mat.scatter([xi] * n_sets, range(n_sets), s=14, color="#d6d5cf", zorder=2, linewidths=0)
        if on:
            ax_mat.scatter([xi] * len(on), on, s=16, color=INK, zorder=3, linewidths=0)
            if len(on) > 1:
                ax_mat.plot([xi, xi], [min(on), max(on)], color=INK, linewidth=1.0, zorder=2)
    ax_mat.set_yticks(range(n_sets))
    ax_mat.set_yticklabels([(set_labels or {}).get(s, s) for s in reversed(set_names)])
    ax_mat.tick_params(axis="both", length=0)
    ax_mat.tick_params(axis="x", labelbottom=False)
    ax_mat.grid(False)
    for sp in ax_mat.spines.values():
        sp.set_visible(False)
    # set sizes (left, growing leftwards)
    set_size = {s: 0 for s in set_names}
    for k, v in items:
        for s in k:
            if s in set_size:
                set_size[s] += v
    ys = [n_sets - 1 - set_names.index(s) for s in set_names]
    ax_set.barh(ys, [set_size[s] for s in set_names], height=0.55, color=INK_3, edgecolor=SURFACE)
    ax_set.invert_xaxis()
    ax_set.set_xlabel(set_xlabel)
    ax_set.tick_params(axis="y", left=False, labelleft=False)
    ax_set.spines["left"].set_visible(False)
    ax_set.grid(axis="x", color=GRID, linewidth=0.5)
    ax_set.grid(axis="y", visible=False)
    return {"bar": ax_bar, "bar_top": ax_top, "matrix": ax_mat, "sets": ax_set, "columns": cols}
