#!/usr/bin/env python3
"""
Generate figures for the LOC115930085 ZNF383 paralog S2 falsifier.

Outputs (in same directory):
  fig_loc115930085_reads.png     — alignment pattern at the locus
  fig_loc115930085_recovery.png  — per-pipeline transcript recovery + exon match
  fig_family1_refinement.png     — bundle mega-region vs per-gene sub-loci

Inputs (paths configurable at top):
  GGO.bam (full genome BAM, indexed)
  GGO_genomic.gff (NCBI annotation)
  /tmp/ks_*.gtf (rustle/stringtie outputs from the S2 run)
  /tmp/ks_vg_refined.tsv (refined family report)
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pysam

ROOT = Path("/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer")
BAM = ROOT / "GGO_KRABchr.bam"
GFF = ROOT / "GGO_genomic.gff"
OUT = ROOT / "docs/figures"

CHROM = "NC_073244.2"
LOC_S, LOC_E = 48_772_798, 48_792_019   # LOC115930085 = ZNF383
PARALOG_S, PARALOG_E = 48_384_759, 48_415_507  # LOC101153727 = ZNF383 paralog
CLUSTER_S, CLUSTER_E = 48_100_000, 49_000_000  # KRAB-ZNF cluster window


# ── Fig 1: read alignment at LOC115930085 ───────────────────────────────────
def read_track(ax, bam_path, chrom, start, end):
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    primary, secondary = [], []
    for r in bam.fetch(chrom, start, end):
        if r.is_unmapped:
            continue
        if r.is_secondary or r.is_supplementary:
            secondary.append(r)
        else:
            primary.append(r)
    primary.sort(key=lambda r: r.reference_start)
    secondary.sort(key=lambda r: r.reference_start)
    y = 0
    for r in primary:
        ax.add_patch(Rectangle((r.reference_start, y - 0.35), r.reference_length, 0.7,
                                facecolor="#2b6fa8", edgecolor="none"))
        y += 1
    gap = 2
    for r in secondary:
        ax.add_patch(Rectangle((r.reference_start, y - 0.35 + gap), r.reference_length, 0.7,
                                facecolor="#d45e3e", edgecolor="none"))
        y += 1
    ax.set_xlim(start, end)
    ax.set_ylim(-1, y + gap)
    ax.set_yticks([len(primary) / 2, len(primary) + gap + len(secondary) / 2])
    ax.set_yticklabels([f"primary\nn={len(primary)}", f"secondary\nn={len(secondary)}"])
    ax.set_xlabel(f"{chrom} (bp)")
    ax.set_title(
        f"LOC115930085 (ZNF383) — {len(primary)+len(secondary)} reads total, "
        f"{len(secondary)}/{len(primary)+len(secondary)} secondary")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return len(primary), len(secondary)


def fig1_reads():
    fig, ax = plt.subplots(figsize=(10, 4))
    p, s = read_track(ax, BAM, CHROM, LOC_S, LOC_E)
    ax.text(0.01, 0.98,
            "StringTie -L discards secondary\n→ sees only 4 reads, cannot assemble",
            transform=ax.transAxes, fontsize=9, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="#f5f5f5", edgecolor="gray"))
    fig.tight_layout()
    fig.savefig(OUT / "fig_loc115930085_reads.png", dpi=150)
    fig.savefig(OUT / "fig_loc115930085_reads.pdf")
    plt.close(fig)
    return p, s


# ── Fig 2: per-pipeline recovery ────────────────────────────────────────────
def count_tx(gtf, chrom, start, end):
    if not Path(gtf).exists():
        return 0, []
    tx, exons_by_tx = [], {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            s, e = int(parts[3]), int(parts[4])
            if e < start or s > end:
                continue
            attrs = parts[8]
            tid = ""
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id "):
                    tid = kv.split('"')[1]
            if parts[2] == "transcript":
                tx.append(tid)
            elif parts[2] == "exon" and tid:
                exons_by_tx.setdefault(tid, []).append((s, e))
    return len(tx), exons_by_tx


def ref_exons(gff, chrom, gene_name):
    exons = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            if parts[2] != "exon":
                continue
            if f"gene={gene_name};" not in parts[8] and f"gene={gene_name}" not in parts[8]:
                continue
            exons.append((int(parts[3]), int(parts[4])))
    # Dedupe
    exons = sorted(set(exons))
    return exons


def fig2_recovery():
    pipelines = [
        ("StringTie -L",           "/tmp/ks_stringtie.gtf", "#808080"),
        ("Rustle default",         "/tmp/ks_rustle.gtf",    "#9a9a9a"),
        ("Rustle --vg none",       "/tmp/ks_vg_none.gtf",   "#e7a54e"),
        ("Rustle --vg (EM)",       "/tmp/ks_vg_refined.gtf","#2b8a4e"),
    ]
    rows = []
    for name, gtf, color in pipelines:
        n, exmap = count_tx(gtf, CHROM, LOC_S, LOC_E)
        # Only keep tx whose exons actually fall within LOC_S..LOC_E
        filtered = {tid: exs for tid, exs in exmap.items()
                    if any(s <= LOC_E and e >= LOC_S for s, e in exs)}
        rows.append((name, len(filtered), filtered, color))

    ref = ref_exons(str(GFF), CHROM, "LOC115930085")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 6),
        gridspec_kw=dict(height_ratios=[1, 2.2], hspace=0.4))

    # Top: bar chart of tx count
    names = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    colors = [r[3] for r in rows]
    bars = ax1.bar(names, counts, color=colors, edgecolor="black", linewidth=0.5)
    for b, c in zip(bars, counts):
        ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.15, str(c),
                 ha="center", fontsize=11)
    ax1.set_ylabel("Transcripts at\nLOC115930085")
    ax1.set_title("Transcript recovery per pipeline")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_ylim(0, max(counts) + 2)

    # Bottom: exon diagrams (ref + Rustle-VG tx)
    y = 0
    track_h = 0.6
    if ref:
        for (s, e) in ref:
            ax2.add_patch(Rectangle((s, y - track_h/2), e - s, track_h,
                                    facecolor="#2b6fa8", edgecolor="black"))
        ax2.plot([min(s for s, _ in ref), max(e for _, e in ref)], [y, y],
                 "-", color="#2b6fa8", lw=0.6)
        ax2.text(LOC_S - 1000, y, "NCBI ref\nXM_031006335.3",
                 ha="right", va="center", fontsize=9)
        y += 1.2
    vg_exons = rows[3][2]  # Rustle --vg (EM) transcripts
    # show up to 2
    for tid in list(vg_exons.keys())[:2]:
        exs = sorted(vg_exons[tid])
        for (s, e) in exs:
            ax2.add_patch(Rectangle((s, y - track_h/2), e - s, track_h,
                                    facecolor="#2b8a4e", edgecolor="black"))
        ax2.plot([min(s for s, _ in exs), max(e for _, e in exs)], [y, y],
                 "-", color="#2b8a4e", lw=0.6)
        ax2.text(LOC_S - 1000, y, tid, ha="right", va="center", fontsize=9)
        y += 1.0
    ax2.set_xlim(LOC_S - 1200, LOC_E + 1000)
    ax2.set_ylim(-0.7, y)
    ax2.set_xlabel(f"{CHROM} (bp)")
    ax2.set_yticks([])
    ax2.set_title("Exon structure: reference vs Rustle-VG recovery")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig_loc115930085_recovery.png", dpi=150)
    fig.savefig(OUT / "fig_loc115930085_recovery.pdf")
    plt.close(fig)


# ── Fig 3: family decomposition before/after ────────────────────────────────
def fig3_family_refinement():
    before_regions = [
        ("NC_073244.2", 47407487, 47528439, "-"),
        ("NC_073244.2", 48184934, 48614722, "-"),
        ("NC_073244.2", 48184934, 48614722, "+"),
        ("NC_073244.2", 48620414, 48999333, "-"),
        ("NC_073244.2", 48620414, 48999333, "+"),
    ]
    after = []
    tsv = Path("/tmp/ks_vg_refined.tsv")
    if tsv.exists():
        with open(tsv) as f:
            f.readline()  # header
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) < 4:
                    continue
                regs = parts[3].split(";")
                # Find the family that contains LOC115930085 (48,772,798).
                if not any(f"48772" in r or f"4877" in r for r in regs):
                    continue
                for r in regs:
                    chrom_s, strand = r.rsplit(":", 1)
                    chrom, span = chrom_s.rsplit(":", 1)
                    s, e = span.split("-")
                    after.append((chrom, int(s), int(e), strand))
                break

    # Reference gene names in cluster
    genes = []
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9 or parts[0] != CHROM or parts[2] != "gene":
                continue
            s, e = int(parts[3]), int(parts[4])
            if e < CLUSTER_S or s > CLUSTER_E:
                continue
            if "zinc finger" not in parts[8]:
                continue
            name = ""
            for kv in parts[8].split(";"):
                if kv.startswith("gene="):
                    name = kv[5:]
            genes.append((s, e, parts[6], name))

    fig, (ax_before, ax_genes, ax_after) = plt.subplots(
        3, 1, figsize=(12, 6),
        gridspec_kw=dict(height_ratios=[1, 1, 1.4], hspace=0.5))

    track_h = 0.7
    for (_c, s, e, st) in before_regions:
        y = 0.5 if st == "+" else -0.5
        ax_before.add_patch(Rectangle((s, y - track_h/2), e - s, track_h,
                                      facecolor="#d45e3e", edgecolor="black", alpha=0.7))
    ax_before.set_xlim(CLUSTER_S, CLUSTER_E)
    ax_before.set_ylim(-1.5, 1.5)
    ax_before.set_yticks([-0.5, 0.5])
    ax_before.set_yticklabels(["− strand", "+ strand"])
    ax_before.set_title(f"BEFORE: {len(before_regions)} bundle-level regions (wide mega-bundles)")
    ax_before.spines["top"].set_visible(False)
    ax_before.spines["right"].set_visible(False)

    for (s, e, st, name) in genes:
        y = 0.5 if st == "+" else -0.5
        ax_genes.add_patch(Rectangle((s, y - track_h/2), e - s, track_h,
                                     facecolor="#2b6fa8", edgecolor="black", alpha=0.8))
        ax_genes.text((s + e) / 2, y + (0.6 if st == "+" else -0.8),
                      name, fontsize=7, ha="center", rotation=30)
    ax_genes.set_xlim(CLUSTER_S, CLUSTER_E)
    ax_genes.set_ylim(-1.8, 1.8)
    ax_genes.set_yticks([-0.5, 0.5])
    ax_genes.set_yticklabels(["− strand", "+ strand"])
    ax_genes.set_title("Reference: annotated ZNF genes in cluster")
    ax_genes.spines["top"].set_visible(False)
    ax_genes.spines["right"].set_visible(False)

    for (_c, s, e, st) in after:
        y = 0.5 if st == "+" else -0.5
        ax_after.add_patch(Rectangle((s, y - track_h/2), e - s, track_h,
                                     facecolor="#2b8a4e", edgecolor="black", alpha=0.8))
    # Highlight LOC115930085
    ax_after.add_patch(Rectangle((LOC_S, -0.5 - track_h/2), LOC_E - LOC_S, track_h,
                                 facecolor="none", edgecolor="red", linewidth=2))
    ax_after.annotate("LOC115930085\n(ZNF383)",
                      xy=((LOC_S + LOC_E)/2, -0.5 - track_h/2),
                      xytext=((LOC_S + LOC_E)/2, -1.6),
                      fontsize=8, ha="center",
                      arrowprops=dict(arrowstyle="->", color="red"))
    # Highlight LOC101153727
    ax_after.add_patch(Rectangle((PARALOG_S, 0.5 - track_h/2), PARALOG_E - PARALOG_S, track_h,
                                 facecolor="none", edgecolor="red", linewidth=2))
    ax_after.annotate("LOC101153727\n(ZNF383 paralog)",
                      xy=((PARALOG_S + PARALOG_E)/2, 0.5 + track_h/2),
                      xytext=((PARALOG_S + PARALOG_E)/2, 1.6),
                      fontsize=8, ha="center",
                      arrowprops=dict(arrowstyle="->", color="red"))
    ax_after.set_xlim(CLUSTER_S, CLUSTER_E)
    ax_after.set_ylim(-2.2, 2.2)
    ax_after.set_yticks([-0.5, 0.5])
    ax_after.set_yticklabels(["− strand", "+ strand"])
    ax_after.set_xlabel(f"{CHROM} (bp)")
    ax_after.set_title(f"AFTER: {len(after)} per-gene sub-loci (RUSTLE_VG_REPORT_GAP refinement)")
    ax_after.spines["top"].set_visible(False)
    ax_after.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT / "fig_family1_refinement.png", dpi=150)
    fig.savefig(OUT / "fig_family1_refinement.pdf")
    plt.close(fig)


if __name__ == "__main__":
    p, s = fig1_reads()
    print(f"fig1 reads: primary={p} secondary={s}")
    fig2_recovery()
    print("fig2 recovery: written")
    fig3_family_refinement()
    print("fig3 refinement: written")
    print(f"outputs in {OUT}")
