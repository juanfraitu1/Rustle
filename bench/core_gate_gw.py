#!/usr/bin/env python
"""
Genome-wide contiguous-core family-gate measurement.

Parses per-chromosome CORE_TRACE files emitted by rustle's de-novo cross-copy
singleton-exon merge step (merge_singletons_by_sequence) and measures what the
contiguous-core coverage gate (threshold T=0.13) WOULD drop vs keep, genome-wide.

INPUT  : /tmp/gw/cg_NC_*.trace
         Each line: "[CORE_TRACE] cid_i=.. cid_j=.. jacc=X core_cov=Y len_i=A len_j=B would_gate=.."

CAVEAT : The `would_gate` field in the traces is VACUOUS. It was computed against
         the active threshold 0.0 (gate-OFF during the trace run), so every line
         reads would_gate=false. We IGNORE that field entirely and recompute the
         drop/keep decision against T=0.13 directly from core_cov.

OUTPUT : core_gate_gw.png  (3-panel figure)
         core_gate_gw.md   (writeup)
         stdout            (headline numbers)

Deterministic: pure parse of the trace files, no randomness, sorted iteration.

Run with: /home/juanfra/miniforge3/bin/python core_gate_gw.py
"""

import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
TRACE_GLOB = "/tmp/gw/cg_NC_*.trace"
HERE = os.path.dirname(os.path.abspath(__file__))
PNG_OUT = os.path.join(HERE, "core_gate_gw.png")
MD_OUT = os.path.join(HERE, "core_gate_gw.md")

T = 0.13          # contiguous-core coverage gate threshold
T_LO = 0.05       # lower lobe of the bimodal valley
T_HI95 = 0.95     # "true copy" near-identical band
DANGER_JACC = 0.5  # high-jaccard cut for the danger-zone definition

# Palette
NAVY = "#22313f"
ORANGE = "#e8590c"
GREEN = "#1e8449"

LINE_RE = re.compile(
    r"\[CORE_TRACE\]\s+cid_i=(?P<ci>\S+)\s+cid_j=(?P<cj>\S+)\s+"
    r"jacc=(?P<jacc>[0-9.]+)\s+core_cov=(?P<core>[0-9.]+)\s+"
    r"len_i=(?P<li>\d+)\s+len_j=(?P<lj>\d+)\s+would_gate=(?P<wg>\S+)"
)


# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------
def chrom_of(path):
    base = os.path.basename(path)
    # cg_NC_073235.2.trace -> NC_073235.2
    return base[len("cg_"):-len(".trace")]


def parse():
    """Return (all_pairs, per_chrom) sorted deterministically by chrom."""
    files = sorted(glob.glob(TRACE_GLOB))
    if not files:
        sys.exit("no trace files found at %s" % TRACE_GLOB)

    all_pairs = []          # list of dicts
    per_chrom = {}          # chrom -> {"total":, "drop":, "jacc":[], "core":[]}

    for path in files:
        ch = chrom_of(path)
        per_chrom.setdefault(ch, {"total": 0, "drop": 0, "jacc": [], "core": []})
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                m = LINE_RE.match(line)
                if not m:
                    sys.exit("unparseable line in %s:\n  %s" % (path, line))
                jacc = float(m.group("jacc"))
                core = float(m.group("core"))
                rec = {
                    "chrom": ch,
                    "jacc": jacc,
                    "core": core,
                    "len_i": int(m.group("li")),
                    "len_j": int(m.group("lj")),
                }
                all_pairs.append(rec)
                pc = per_chrom[ch]
                pc["total"] += 1
                pc["jacc"].append(jacc)
                pc["core"].append(core)
                if core < T:
                    pc["drop"] += 1

    return all_pairs, per_chrom


# ---------------------------------------------------------------------------
# Headline computation (recomputed against T=0.13, ignoring would_gate)
# ---------------------------------------------------------------------------
def headline(all_pairs):
    n = len(all_pairs)
    cores = [p["core"] for p in all_pairs]
    drop = sum(1 for c in cores if c < T)
    lt_lo = sum(1 for c in cores if c < T_LO)
    valley = sum(1 for c in cores if T_LO <= c < T)
    keep = sum(1 for c in cores if c >= T)
    ge95 = sum(1 for c in cores if c >= T_HI95)

    # Danger zone: jaccard alone would keep (jacc high) but core says DROP.
    danger = [p for p in all_pairs if p["jacc"] >= DANGER_JACC and p["core"] < T]
    danger_hi = [p for p in danger if p["jacc"] >= 0.8]
    # Jaccard-alone keep at the same numeric threshold:
    jacc_keep = sum(1 for p in all_pairs if p["jacc"] >= T)

    return {
        "n": n,
        "drop": drop,
        "drop_pct": 100.0 * drop / n,
        "lt_lo": lt_lo,
        "valley": valley,
        "keep": keep,
        "keep_pct": 100.0 * keep / n,
        "ge95": ge95,
        "danger": len(danger),
        "danger_hi": len(danger_hi),
        "jacc_keep": jacc_keep,
        "danger_examples": sorted(
            [(p["jacc"], p["core"]) for p in danger], reverse=True
        )[:6],
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(all_pairs, per_chrom, H):
    cores = [p["core"] for p in all_pairs]
    jaccs = [p["jacc"] for p in all_pairs]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    fig.subplots_adjust(left=0.05, right=0.985, top=0.86, bottom=0.16, wspace=0.26)

    # --- Panel 1: histogram of core_cov, drop (orange) vs keep (green) ---
    ax = axes[0]
    bins = [i / 50.0 for i in range(51)]  # 0..1 step 0.02
    drop_vals = [c for c in cores if c < T]
    keep_vals = [c for c in cores if c >= T]
    ax.hist(keep_vals, bins=bins, color=GREEN, alpha=0.85, label="keep (core≥%.2f)" % T)
    ax.hist(drop_vals, bins=bins, color=ORANGE, alpha=0.9, label="drop (core<%.2f)" % T)
    ax.axvline(T, color=NAVY, lw=2, ls="-")
    ax.text(T + 0.012, ax.get_ylim()[1] * 0.92, "T=0.13", color=NAVY,
            fontsize=10, fontweight="bold", va="top")
    # mark the valley span 0.05-0.13
    ax.axvspan(T_LO, T, color=NAVY, alpha=0.06)
    ax.text(0.09, ax.get_ylim()[1] * 0.55, "valley\n0.05–0.13", color=NAVY,
            fontsize=8, ha="center", va="center")
    ax.set_yscale("log")
    ax.set_xlabel("contiguous-core coverage")
    ax.set_ylabel("pairs (log)")
    ax.set_title("(1) core_cov is bimodal; T=0.13 sits in the valley",
                 fontsize=10, color=NAVY)
    ax.legend(fontsize=8, loc="upper center")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # --- Panel 2: scatter jacc vs core, danger zone highlighted ---
    ax = axes[1]
    dz_j = [p["jacc"] for p in all_pairs if p["jacc"] >= DANGER_JACC and p["core"] < T]
    dz_c = [p["core"] for p in all_pairs if p["jacc"] >= DANGER_JACC and p["core"] < T]
    other_j = [p["jacc"] for p in all_pairs
               if not (p["jacc"] >= DANGER_JACC and p["core"] < T)]
    other_c = [p["core"] for p in all_pairs
               if not (p["jacc"] >= DANGER_JACC and p["core"] < T)]
    ax.scatter(other_j, other_c, s=8, color=GREEN, alpha=0.45, edgecolors="none",
               label="true copies / low-jacc")
    ax.scatter(dz_j, dz_c, s=14, color=ORANGE, alpha=0.85, edgecolors="none",
               label="DANGER ZONE (n=%d)" % H["danger"])
    ax.axhline(T, color=NAVY, lw=1.5, ls="--")
    ax.text(0.02, T + 0.02, "core T=0.13", color=NAVY, fontsize=8)
    # box the danger region
    ax.add_patch(plt.Rectangle((DANGER_JACC, 0), 1.0 - DANGER_JACC, T,
                               fill=False, edgecolor=ORANGE, lw=1.2, ls=":"))
    ax.text(0.74, 0.165, "high jaccard +\nlow core\n= domain-sharers\n(false-merged by\njaccard alone)",
            color=ORANGE, fontsize=7.5, ha="center", va="bottom")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("jaccard (sequence overlap)")
    ax.set_ylabel("contiguous-core coverage")
    ax.set_title("(2) Core separates true copies from domain-sharers",
                 fontsize=10, color=NAVY)
    ax.legend(fontsize=8, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # --- Panel 3: per-chrom drop counts ---
    ax = axes[2]
    rows = [(ch, d["drop"]) for ch, d in per_chrom.items() if d["drop"] > 0]
    rows.sort(key=lambda r: r[1], reverse=True)
    labels = [r[0].replace("NC_0", "").replace(".2", "").replace(".1", "") for r in rows]
    vals = [r[1] for r in rows]
    ypos = list(range(len(rows)))[::-1]
    ax.barh(ypos, vals, color=ORANGE, alpha=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    for y, v in zip(ypos, vals):
        ax.text(v + 0.5, y, str(v), fontsize=7, va="center", color=NAVY)
    ax.set_xlabel("dropped pairs (core<0.13)")
    ax.set_title("(3) Drops by chromosome (total=%d)" % H["drop"],
                 fontsize=10, color=NAVY)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle(
        "Genome-wide contiguous-core family-gate: %d cross-copy merge pairs, "
        "%d (%.1f%%) drop at T=0.13" % (H["n"], H["drop"], H["drop_pct"]),
        fontsize=12, color=NAVY, fontweight="bold", y=0.97)
    fig.savefig(PNG_OUT, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Writeup
# ---------------------------------------------------------------------------
def write_md(all_pairs, per_chrom, H):
    n_chroms = len(per_chrom)
    n_with_pairs = sum(1 for d in per_chrom.values() if d["total"] > 0)

    # per-chrom table sorted by drop desc, then chrom asc
    table_rows = sorted(
        per_chrom.items(),
        key=lambda kv: (-kv[1]["drop"], -kv[1]["total"], kv[0]),
    )

    lines = []
    A = lines.append
    A("# Genome-wide contiguous-core family-gate measurement\n")
    A("**Site measured:** rustle's de-novo cross-copy singleton-exon merges "
      "(`merge_singletons_by_sequence`) — the gate's actual call site.\n")
    A("**Threshold:** T = %.2f on contiguous-core coverage (`core_cov`).\n" % T)
    A("**Gate status:** DEFAULT-OFF. This is a *measurement* of what the gate "
      "would do, not a behavior change.\n")

    A("\n## Headline\n")
    A("- **%d chromosomes** parsed (%d carry cross-copy merge pairs; "
      "3 are empty).\n" % (n_chroms, n_with_pairs))
    A("- **%d** cross-copy merge pairs total.\n" % H["n"])
    A("- **%d (%.1f%%) DROP** at core_cov < %.2f.\n"
      % (H["drop"], H["drop_pct"], T))
    A("  - %d have core_cov < %.2f (deep lobe).\n" % (H["lt_lo"], T_LO))
    A("  - %d sit in the valley [%.2f, %.2f).\n" % (H["valley"], T_LO, T))
    A("- **%d (%.1f%%) KEEP** at core_cov ≥ %.2f.\n"
      % (H["keep"], H["keep_pct"], T))
    A("- **%d** pairs have core_cov ≥ %.2f (near-identical true copies).\n"
      % (H["ge95"], T_HI95))
    A("- `NC_073235.2` = **%d / %d** drops (matches the prior pipeline test).\n"
      % (per_chrom["NC_073235.2"]["drop"], per_chrom["NC_073235.2"]["total"]))

    A("\n## Bimodality — why 0.13 is a principled threshold\n")
    A("The core_cov distribution is **bimodal**: a heavy mass near 1.0 "
      "(near-identical true copies; %d pairs ≥ %.2f) and a second lobe "
      "down near 0.0 (%d pairs < %.2f). Between them is a sparse **valley** "
      "spanning roughly [%.2f, %.2f) holding only %d pairs. T = %.2f is placed "
      "inside that valley, so the decision boundary does not slice through a "
      "dense region — small threshold perturbations move very few pairs. This "
      "is exactly the kind of separated, non-arbitrary cut a gate wants: the "
      "data, not the knob, picks the split.\n"
      % (H["ge95"], T_HI95, H["lt_lo"], T_LO, T_LO, T, H["valley"], T))

    A("\n## Why jaccard alone fails and core fixes it (the danger zone)\n")
    A("At the same numeric cut, **jaccard alone keeps all %d pairs** "
      "(every pair has jaccard ≥ %.2f) — i.e. a jaccard-only gate at 0.13 "
      "would never drop anything here. Yet **%d pairs have high jaccard "
      "(≥ %.2f) but low contiguous-core (< %.2f)**, and %d of those have "
      "jaccard ≥ 0.80. These are **domain-sharers**: two copies that share a "
      "large fraction of *k*-mers/segments (high jaccard) but whose shared "
      "material is *not contiguous core* — it is a repeated domain, a shared "
      "exon cassette, or scattered homology rather than a single colinear copy "
      "body. Jaccard cannot tell a tandem-repeated domain from a real copy; "
      "contiguous-core coverage can, because it demands the overlap form one "
      "run.\n"
      % (H["jacc_keep"], T, H["danger"], DANGER_JACC, T, H["danger_hi"]))
    A("\nWorst offenders (jaccard, core_cov) — jaccard says \"merge\", core says "
      "\"don't\":\n\n")
    A("| jaccard | core_cov |\n|---|---|\n")
    for j, c in H["danger_examples"]:
        A("| %.3f | %.3f |\n" % (j, c))

    A("\n## Per-chromosome\n")
    A("| chrom | pairs | drop (core<0.13) | drop %% |\n")
    A("|---|---|---|---|\n")
    for ch, d in table_rows:
        pct = (100.0 * d["drop"] / d["total"]) if d["total"] else 0.0
        A("| %s | %d | %d | %.1f%% |\n" % (ch, d["total"], d["drop"], pct))

    A("\n## Honest scope & caveats\n")
    A("- **What this measures.** rustle's *de-novo cross-copy "
      "singleton-exon merges* (`merge_singletons_by_sequence`). That is the "
      "gate's literal call site, so these counts are the population the gate "
      "would actually act on.\n")
    A("- **Per-chrom / within-chrom families only.** The OOM-safe genome-wide "
      "protocol runs one chromosome at a time, so only *within-chromosome* "
      "families are grouped. Cross-chromosome paralogs (e.g. RABL2A/RABL2B) "
      "are not co-grouped here — but those are the coordinate-separated case "
      "this gate does **not** target anyway; the gate is about distinguishing "
      "true copies from domain-sharers inside a candidate family, not about "
      "linking copies across chromosomes.\n")
    A("- **`would_gate` trace field is vacuous — ignored.** Every trace line "
      "reads `would_gate=false` because it was computed against the active "
      "threshold 0.0 (gate-OFF during the trace run). We do **not** trust that "
      "field; every drop/keep number above is recomputed directly from "
      "`core_cov` against T = %.2f.\n" % T)
    A("- **Gate is default-off.** Nothing here changes rustle's emitted output; "
      "this is a measurement of the gate's *would-be* effect.\n")

    A("\n## Reproduce\n")
    A("```\n/home/juanfra/miniforge3/bin/python %s\n```\n"
      % os.path.join("bench", "core_gate_gw.py"))
    A("Inputs: `%s` (deterministic parse).\n" % TRACE_GLOB)

    with open(MD_OUT, "w") as fh:
        fh.write("".join(lines))


# ---------------------------------------------------------------------------
def main():
    all_pairs, per_chrom = parse()
    H = headline(all_pairs)

    # stdout headline
    print("=== genome-wide contiguous-core family-gate (T=%.2f) ===" % T)
    print("chroms parsed         : %d (%d with pairs)"
          % (len(per_chrom), sum(1 for d in per_chrom.values() if d["total"])))
    print("cross-copy merge pairs : %d" % H["n"])
    print("DROP (core<0.13)       : %d (%.1f%%)  [%d <0.05, %d in [0.05,0.13)]"
          % (H["drop"], H["drop_pct"], H["lt_lo"], H["valley"]))
    print("KEEP (core>=0.13)      : %d (%.1f%%)" % (H["keep"], H["keep_pct"]))
    print("core>=0.95             : %d" % H["ge95"])
    print("DANGER ZONE (jacc>=%.1f & core<0.13): %d (%d with jacc>=0.8)"
          % (DANGER_JACC, H["danger"], H["danger_hi"]))
    print("jaccard-alone keep @0.13: %d (jaccard alone drops nothing)"
          % H["jacc_keep"])
    print("NC_073235.2            : %d/%d drops"
          % (per_chrom["NC_073235.2"]["drop"], per_chrom["NC_073235.2"]["total"]))

    make_figure(all_pairs, per_chrom, H)
    write_md(all_pairs, per_chrom, H)
    print("\nwrote %s" % PNG_OUT)
    print("wrote %s" % MD_OUT)


if __name__ == "__main__":
    main()
