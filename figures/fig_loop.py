"""Figure 9 — the closed loop: tied reads assigned to one copy, given to that copy, and re-assembled.

Pre-registration: docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md (user decision 2026-09-25 16:00, item 3). Pipeline:
tools/rustle_reassemble.sh (union -> home -> pass2 -> g0) after the driver's assemble and families stages; the home
table is bench/loop_home.py, the pass-2 read filter is RUSTLE_READ_HOME_TABLE (src/rustle/vg_family/denovo_assemble.rs).

  (a) G0, per sample: records the loop moves (removed away from home, R_drop; admitted at home beyond the default
      read pool, R_add) against the pre-registered floor of 20 records (below it the sample is bounded, not scored).
  (b) Simulation (truth known): simulated copies whose exact intron chain is recovered, per arm (P, default, LOOP-S,
      LOOP-P, ORACLE = every read given to its true copy).
  (c) Annotation-free intron-chain sensitivity and precision against the annotation (as Fig. 1), per sample and arm.
  (d) Reference chains gained minus lost by LOOP-S against the default, by tie fraction (Fig. 3's bins), with the
      number of read-sharing groups that carry the gains.
  (e) Loci in the Liftoff framework (Fig. 8 C2): read-supported Liftoff loci found by the arm's de novo loci, on the
      primary-read universe U1 and the candidate-placement universe U2 (which counts tied secondaries).

Tables (each read only when present; the panel says "not built" otherwise):
  fig9_loop_g0       per sample: A, H, home-table molecules, kept at home, R_add, R_drop, G0 status
  fig9_loop_sim      per simulation and arm: M5 counts
  fig9_loop_chains   per sample and arm: M1 (and M2 paired against the arm's base)
  fig9_loop_ties     per sample, tie bin: M3
  fig9_loop_liftoff  per sample, arm, universe: M4

Build: `python3 figures/make.py data fig9` tabulates what the loop runs have produced (runs nothing heavy);
`--recorded` tabulates the development runs (key fig9_dev_dir) into data_recorded/, marked provisional.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import figlib

TABLES = ["fig9_loop_g0", "fig9_loop_sim", "fig9_loop_chains", "fig9_loop_ties", "fig9_loop_liftoff"]
ARMS = ["P", "GOOD", "LOOP-S", "LOOP-P", "ORACLE"]
ARM_LABEL = {"P": "Rustle, primary alignments only", "GOOD": "Rustle (default)",
             "LOOP-S": "Loop: assigned reads at their copy (default otherwise)",
             "LOOP-P": "Loop: assigned reads at their copy (primaries otherwise)",
             "ORACLE": "Every read at its true copy (simulation only)"}
ARM_COLOR = {"P": figlib.BLUE[250], "GOOD": figlib.BLUE[450], "LOOP-S": figlib.BLUE[650], "LOOP-P": figlib.BLUE[350],
             "ORACLE": figlib.INK_2}
ARM_HATCH = {"P": "////", "LOOP-P": "////"}
SAMPLE_ORDER = ["human_testis", "gorilla_KB3781", "chimp_PTR", "orangutan_PPY", "human_A119b", "gorilla_OR6737"]
VERDICT = {"human_testis", "gorilla_KB3781", "chimp_PTR", "orangutan_PPY"}
G0_FLOOR = 20
GEN = "figures/fig_loop.py"


def _existing_tables(data_dir: Path = figlib.DATA_DIR) -> list[str]:
    return [t for t in TABLES if (Path(data_dir) / f"{t}.tsv").exists()]


META = {
    "id": "fig9",
    "title": "The closed loop: tied reads assigned to one copy, given to that copy and re-assembled",
    "claim": ("Pre-registered (docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md): after pass 1 (the default assembly and "
              "families), reads with equally good alignments at several copies are assigned by one test over every "
              "candidate copy; each assigned read is then taken only at its copy and the sample is re-assembled. "
              "Judged like the annotation-free tools (Fig. 1, Fig. 3) and in the Liftoff framework (Fig. 8), with "
              "the verdict on four held-out samples (claims LOOP1-LOOP3)."),
    "tables": _existing_tables(),
    "supplementary_tables": [t for t in TABLES if t not in _existing_tables()],
}


# ================================================================ build (tabulates; runs nothing heavy)
def _read_kv(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        for line in fh:
            k, _, v = line.rstrip("\n").partition("\t")
            if k and k != "item":
                out[k] = v
    return out


def _g0_row(sid: str, species: str, prefix: Path) -> list | None:
    summ, g0 = Path(f"{prefix}.loop.home.summary.tsv"), Path(f"{prefix}.pass2.g0.tsv")
    if not g0.exists():
        return None
    s = _read_kv(summ) if summ.exists() else {}
    g = _read_kv(g0)
    moved = int(g.get("R_drop", 0)) + int(g.get("R_add", 0))
    return [sid, species, "verdict" if sid in VERDICT else "development", s.get("A_assigned", ""),
            s.get("H_no_record_at_home", ""), g.get("home_table_molecules", ""), g.get("kept_at_home", ""),
            g.get("R_add", ""), g.get("R_drop", ""), "bounded" if moved < G0_FLOOR else "open"]


G0_HEADER = ["sample", "species", "exposure", "A_assigned", "H_no_record_at_home", "home_table_molecules",
             "kept_at_home", "R_add", "R_drop", "g0"]
SIM_HEADER = ["simulation", "arm", "n_transcripts", "n_multi_exon_copies", "multi_exon_chain_found",
              "n_single_exon_copies", "single_exon_covered", "n_transcripts_at_copies", "n_transcripts_matching_a_copy"]


def _sim_rows(name: str, m5: Path) -> list[list]:
    with open(m5) as fh:
        return [[name] + [r[h] for h in SIM_HEADER[1:]] for r in csv.DictReader(fh, delimiter="\t")]


def build(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    import samples
    notes = ["pre-registration docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md; union-certificate verdicts only (r1092); "
             "families frozen from pass 1 (r395); tie widths: assignment exact tie (1.0), seeding 0.98"]
    if recorded:
        dev = Path(cfg.get("fig9_dev_dir") or "/mnt/linuxdisk/tmp/rustle_figures_dev/loop")
        rows = []
        for m5 in sorted(dev.glob("*/m5.tsv")):
            rows += _sim_rows(m5.parent.name, m5)
        if rows:
            figlib.write_table("fig9_loop_sim", SIM_HEADER, rows, generator=GEN, data_dir=data_dir,
                               inputs={p.parent.name: p for p in sorted(dev.glob("*/m5.tsv"))},
                               notes=["provisional: development simulation (legacy catalog; not a verdict)"] + notes)
        META["tables"] = _existing_tables(data_dir)
        return
    reg = samples.registry(cfg)
    g0, inputs = [], {}
    for sid in SAMPLE_ORDER:
        if sid not in reg:
            continue
        prefix = samples.run_dir(cfg, sid) / sid
        r = _g0_row(sid, reg[sid]["species"], prefix)
        if r:
            g0.append(r)
            inputs[f"{sid}_g0"] = Path(f"{prefix}.pass2.g0.tsv")
    if g0:
        figlib.write_table("fig9_loop_g0", G0_HEADER, g0, generator=GEN, inputs=inputs, notes=notes, data_dir=data_dir)
    sims = [s.partition("=") for s in str(cfg.get("fig9_sims", "")).split(",") if "=" in s]
    rows = []
    for name, _, d in sims:
        m5 = Path(d) / "m5.tsv"
        if m5.exists():
            rows += _sim_rows(name, m5)
    if rows:
        figlib.write_table("fig9_loop_sim", SIM_HEADER, rows, generator=GEN, notes=notes, data_dir=data_dir,
                           inputs={n: Path(d) / "m5.tsv" for n, _, d in sims})
    print("[fig9] M1-M4 (fig9_loop_chains/_ties/_liftoff) are tabulated from <sample>.pass2.gtf / .pass2p.gtf with "
          "figures/assembly.py (gffcompare, Fig. 1 and 3 rules) and figures/_liftoff.py (C2, universes U1 and U2) once "
          "the pass-2 assemblies exist; until then those panels say 'not built'.", file=sys.stderr)
    META["tables"] = _existing_tables(data_dir)


# ================================================================ plot
def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _na(ax, text):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=6, color=figlib.INK_3,
            wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.grid(False)


def _panel_g0(ax, rows):
    rows = sorted(rows, key=lambda r: SAMPLE_ORDER.index(r["sample"]) if r["sample"] in SAMPLE_ORDER else 99)
    y = list(range(len(rows)))[::-1]
    drop = [int(r["R_drop"] or 0) for r in rows]
    add = [int(r["R_add"] or 0) for r in rows]
    ax.barh(y, drop, color=figlib.BLUE[650], height=0.6, label="Removed away from the assigned copy")
    ax.barh(y, add, left=drop, color=figlib.BLUE[250], height=0.6, label="Admitted at the assigned copy")
    ax.axvline(G0_FLOOR, color=figlib.INK_2, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['sample']}{'' if r['exposure'] == 'verdict' else ' (dev.)'}" for r in rows])
    ax.set_xlabel("Alignment records moved (floor: 20)")
    ax.grid(axis="x")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=5.4, ncol=1)


def _panel_sim(ax, rows):
    sims = sorted({r["simulation"] for r in rows})
    width = 0.8 / len(ARMS)
    for i, arm in enumerate(ARMS):
        xs, vals = [], []
        for j, s in enumerate(sims):
            r = next((r for r in rows if r["simulation"] == s and r["arm"] == arm), None)
            if r and int(r["n_multi_exon_copies"]):
                xs.append(j + (i - (len(ARMS) - 1) / 2) * width)
                vals.append(int(r["multi_exon_chain_found"]) / int(r["n_multi_exon_copies"]))
        if xs:
            ax.bar(xs, vals, width=width, color=ARM_COLOR[arm], hatch=ARM_HATCH.get(arm), edgecolor=figlib.SURFACE,
                   linewidth=0.5, label=ARM_LABEL[arm])
    ax.set_xticks(range(len(sims)))
    ax.set_xticklabels(sims)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Simulated multi-exon copies with\ntheir exact intron chain")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.2, ncol=1)


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    figlib.use_style()
    data_dir = Path(data_dir)
    have = set(_existing_tables(data_dir))
    tab = {t: (figlib.read_table(t, data_dir) if t in have else []) for t in TABLES}
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE, 150 * figlib.MM))
    gs = fig.add_gridspec(3, 2, left=0.14, right=0.97, top=0.94, bottom=0.08, hspace=1.0, wspace=0.55)
    axes = {"a": fig.add_subplot(gs[0, 0]), "b": fig.add_subplot(gs[0, 1]), "c": fig.add_subplot(gs[1, 0]),
            "d": fig.add_subplot(gs[1, 1]), "e": fig.add_subplot(gs[2, :])}
    pending = "not built (tools/rustle_reassemble.sh, then make.py data fig9)"
    _panel_g0(axes["a"], tab["fig9_loop_g0"]) if tab["fig9_loop_g0"] else _na(axes["a"], "G0 counts: " + pending)
    _panel_sim(axes["b"], tab["fig9_loop_sim"]) if tab["fig9_loop_sim"] else _na(axes["b"], "Simulation: " + pending)
    _na(axes["c"], "Intron-chain sensitivity and precision per arm: " + pending)
    _na(axes["d"], "Chains gained / lost by tie fraction: " + pending)
    _na(axes["e"], "Liftoff loci found, universes U1 (primary reads) and U2 (candidate placements): " + pending)
    titles = {"a": "Records the loop moves, per sample (G0)",
              "b": "Simulation: copies recovered, per arm",
              "c": "Annotation-free intron chains vs the annotation",
              "d": "Chains gained minus lost vs the default, by tie fraction",
              "e": "Loci in the Liftoff framework (reference, not competitor)"}
    for k, ax in axes.items():
        figlib.panel_label(ax, k, x=-0.02, y=1.08)
        ax.text(0.03, 1.08, titles[k], transform=ax.transAxes, fontsize=6.8, fontweight="bold", va="bottom",
                ha="left")
    fig.text(0.01, 0.995, "Annotation-free (de novo) re-assembly; Rustle-internal arms, not a tool comparison. "
             "Species never pooled.", fontsize=5.8, color=figlib.INK_2, ha="left", va="top")
    figlib.stamp_provisional(fig, TABLES, data_dir)
    paths = figlib.save(fig, "fig9_loop", out_dir)
    plt.close(fig)
    return paths


# ================================================================ caption numbers
def summary(data_dir: Path = figlib.DATA_DIR):
    """Print G0 per sample and the M5 rows; LOOP1 needs the M1-M4 tables (not built yet)."""
    have = set(_existing_tables(data_dir))
    if "fig9_loop_g0" in have:
        for r in figlib.read_table("fig9_loop_g0", data_dir):
            print(f"G0 {r['sample']} ({r['exposure']}): A {r['A_assigned']}, H {r['H_no_record_at_home']}, R_drop "
                  f"{r['R_drop']}, R_add {r['R_add']} -> {r['g0']}")
    if "fig9_loop_sim" in have:
        for r in figlib.read_table("fig9_loop_sim", data_dir):
            print(f"M5 {r['simulation']} {r['arm']}: multi-exon {r['multi_exon_chain_found']} of "
                  f"{r['n_multi_exon_copies']}, single-exon {r['single_exon_covered']} of {r['n_single_exon_copies']}, "
                  f"transcripts matching a copy {r['n_transcripts_matching_a_copy']} of {r['n_transcripts_at_copies']}")
    for t in TABLES:
        if t not in have:
            print(f"{t}: not built")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    if sys.argv[1:2] == ["summary"]:
        summary(Path(sys.argv[2]) if len(sys.argv) > 2 else figlib.DATA_DIR)
    else:
        sys.exit("usage: python3 figures/fig_loop.py summary [DATA_DIR]")
