#!/usr/bin/env python3
"""Controlled in-silico validation of the read-coherence + PSV copy-split in the
co-located / collapsed-tandem regime, and an empirical map of the IDENTIFIABILITY
BOUNDARY (the thesis identifiability theorem).

MODEL
-----
A family of N paralog copies that SHARE ONE intron chain (co-located / collapsed
tandem). This is the adversarial regime: read-coherence (the structural axis) alone
collapses the whole family to a SINGLE transcript, because every copy traverses the
same ordered intron chain. The only separating signal is the COPY axis: P paralog-
sequence-variant (PSV) columns, where each copy carries a distinct fixed allele
vector that differs from every other copy at >= K columns.

Reads are simulated per copy. Each read covers a random contiguous sub-interval of
the transcript (models partial spanning / finite read length L): it observes only
the PSV columns inside its span, every other column is None. Each OBSERVED allele is
independently corrupted with per-base error e (HiFi ~ 0.003) to one of the other
three bases.

The split logic is a FAITHFUL port of `split_readchain_by_psv` (and its helpers
`select_identifiable_copies`, `read_consistent_with`, `distinct_columns`,
`consensus_haplotype`) from src/rustle/vg_family/copy_split.rs. Group by chain;
candidate copies = distinct allele vectors with >= min_reads support; split iff
>= 2 pairwise >= K-distinct copies survive the greedy identifiable-set selection.

IDENTIFIABILITY THEOREM (empirical form)
----------------------------------------
A copy is recovered iff enough error-FREE reads that span all of its >= K
distinguishing columns accumulate on its EXACT true allele vector to clear the
min_reads threshold, AND it is >= K-distinct from every already-admitted copy.
Per-base error e fragments support off the true vector at rate ~ (1 - (1-e)^s)
where s = #columns the read spans. The boundary therefore sits where the expected
count of error-free spanning reads per copy crosses min_reads: recovery -> 1.0 once
coverage * P(read spans >= K cols) * (1-e)^(#spanned) clears threshold, -> merged
below it.

DETERMINISM
-----------
random.seed and numpy.random.seed are set from FIXED INTEGER CONSTANTS only (never
from the wall clock or any nondeterministic source). The whole pipeline is re-seeded
per configuration so reruns are byte-identical.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# Deterministic seeds: FIXED INTEGER CONSTANTS ONLY (never wall-clock).
# ---------------------------------------------------------------------------
MASTER_SEED = 20260616          # global master constant
DEMO_SEED = 11                  # single-config demo
SWEEP_SEED_BASE = 700           # base for per-cell sweep seeds (cell offset added)

BASES = b"ACGT"                 # allele alphabet


# ===========================================================================
# FAITHFUL PORT of src/rustle/vg_family/copy_split.rs
# ReadObs / CopyIsoform and split_readchain_by_psv + helpers.
# In the Rust, an allele is a u8 base or None; here Some(base:int) / None.
# ===========================================================================

@dataclass(frozen=True)
class ReadObs:
    # read-coherence STRUCTURAL key: ordered intron chain (donor, acceptor).
    intron_chain: tuple[tuple[int, int], ...]
    # COPY axis: allele at each family PSV column (None if read does not span it).
    psv_alleles: tuple[int | None, ...]


@dataclass
class CopyIsoform:
    intron_chain: tuple[tuple[int, int], ...]
    allele_vector: tuple[int | None, ...]
    read_count: int
    identifiable: bool


def read_consistent_with(read, copy) -> bool:
    """True if read's observed (non-None) columns all agree with copy's alleles.
    Mirrors copy_split.rs::read_consistent_with."""
    for r, c in zip(read, copy):
        if r is not None and c is not None:
            if r != c:
                return False
        elif r is not None and c is None:
            return False  # read observes a base where copy haplotype is unobserved
        # (r is None): no constraint
    return True


def distinct_columns(a, b) -> int:
    """Count columns where both observe a base and the bases differ.
    Mirrors copy_split.rs::distinct_columns."""
    n = 0
    for x, y in zip(a, b):
        if x is not None and y is not None and x != y:
            n += 1
    return n


def select_identifiable_copies(candidates, min_psv_k):
    """Greedy maximal set pairwise >= K-distinct, candidates in sorted order.
    Mirrors copy_split.rs::select_identifiable_copies."""
    chosen = []
    for cand in candidates:
        if all(distinct_columns(cand, c) >= min_psv_k for c in chosen):
            chosen.append(cand)
    return chosen


def consensus_haplotype(group, n_psv_columns):
    """Per-column majority (lowest base wins ties), None where unobserved.
    Mirrors copy_split.rs::consensus_haplotype (BTreeMap tally, lowest base on tie)."""
    out = []
    for col in range(n_psv_columns):
        tally: dict[int, int] = {}
        for r in group:
            b = r.psv_alleles[col]
            if b is not None:
                tally[b] = tally.get(b, 0) + 1
        if not tally:
            out.append(None)
            continue
        # highest count; ties broken by lowest base.
        best = None
        best_count = -1
        for base in sorted(tally):  # ascending base => lowest wins ties
            c = tally[base]
            if c > best_count:
                best_count = c
                best = base
        out.append(best)
    return tuple(out)


def _vec_sort_key(vec):
    """BTreeMap<Vec<Option<u8>>> ordering: Rust orders None < Some, then by base.
    Encode None as (0,0) and Some(b) as (1,b) per column so the tuple sorts like
    Rust's Option<u8> Ord."""
    return tuple((0, 0) if a is None else (1, a) for a in vec)


def split_readchain_by_psv(reads, n_psv_columns, min_psv_k, min_reads_per_copy):
    """Faithful port of copy_split.rs::split_readchain_by_psv."""
    # 1. group by EXACT intron chain (structural axis); BTreeMap order = sorted chains.
    by_chain: dict[tuple, list[ReadObs]] = {}
    for r in reads:
        by_chain.setdefault(r.intron_chain, []).append(r)

    out: list[CopyIsoform] = []

    for chain in sorted(by_chain.keys()):
        group = by_chain[chain]

        # 2. tally distinct non-all-None allele vectors.
        vector_counts: dict[tuple, int] = {}
        for r in group:
            if all(a is None for a in r.psv_alleles):
                continue  # fully ambiguous: can never define a copy
            vector_counts[r.psv_alleles] = vector_counts.get(r.psv_alleles, 0) + 1

        # candidate copies = vectors with >= min_reads support, in BTreeMap order.
        cand_keys = sorted(vector_counts.keys(), key=_vec_sort_key)
        candidates = [v for v in cand_keys if vector_counts[v] >= min_reads_per_copy]

        if n_psv_columns >= min_psv_k:
            split_copies = select_identifiable_copies(candidates, min_psv_k)
        else:
            split_copies = []

        if len(split_copies) >= 2:
            # 2a. SPLIT: apportion each read to its UNIQUE consistent copy.
            counts = [0] * len(split_copies)
            for r in group:
                consistent = None
                unique = True
                for i, copy in enumerate(split_copies):
                    if read_consistent_with(r.psv_alleles, copy):
                        if consistent is not None:
                            unique = False
                            break
                        consistent = i
                if unique and consistent is not None:
                    counts[consistent] += 1
            for copy, count in zip(split_copies, counts):
                out.append(CopyIsoform(chain, copy, count, True))
        else:
            # 3. MERGED, not identifiable.
            allele_vector = consensus_haplotype(group, n_psv_columns)
            out.append(CopyIsoform(chain, allele_vector, len(group), False))

    out.sort(key=lambda c: (c.intron_chain, _vec_sort_key(c.allele_vector)))
    return out


# ===========================================================================
# SIMULATION: co-located / collapsed tandem family, partial-spanning reads.
# ===========================================================================

@dataclass
class FamilyTruth:
    n_copies: int
    n_psv: int
    shared_chain: tuple[tuple[int, int], ...]
    allele_vectors: list[tuple[int, ...]]   # one fixed FULL allele vector per copy


def make_family(n_copies, n_psv, min_k, rng) -> FamilyTruth:
    """N copies sharing ONE intron chain, each a distinct full PSV allele vector
    differing from every other by >= min_k columns. Rejection-sampled deterministically
    from rng so any two copies are >= min_k-distinct (the identifiability premise)."""
    shared_chain = ((1000, 2000), (3000, 4000))  # ONE shared chain for the whole family
    vectors: list[tuple[int, ...]] = []
    attempts = 0
    while len(vectors) < n_copies:
        attempts += 1
        if attempts > 100000:
            raise RuntimeError("cannot place copies >= K-distinct; raise n_psv")
        cand = tuple(BASES[rng.randrange(4)] for _ in range(n_psv))
        if all(sum(a != b for a, b in zip(cand, v)) >= min_k for v in vectors):
            vectors.append(cand)
    return FamilyTruth(n_copies, n_psv, shared_chain, vectors)


def simulate_reads(fam: FamilyTruth, coverage, span_frac, err, rng) -> list[ReadObs]:
    """Per copy, `coverage` reads. Each read covers a contiguous block of PSV columns
    of length ~ span_frac * P (>=1), at a random start (partial spanning / read length).
    Observed alleles inside the span are error-corrupted at per-base rate `err`."""
    P = fam.n_psv
    span_len = max(1, int(round(span_frac * P)))
    reads: list[ReadObs] = []
    for true_vec in fam.allele_vectors:
        for _ in range(coverage):
            # contiguous span [start, start+span_len)
            if span_len >= P:
                start = 0
                slen = P
            else:
                start = rng.randrange(0, P - span_len + 1)
                slen = span_len
            alleles: list[int | None] = [None] * P
            for col in range(start, start + slen):
                base = true_vec[col]
                if rng.random() < err:
                    # corrupt to one of the OTHER three bases (deterministic via rng)
                    others = [b for b in BASES if b != base]
                    base = others[rng.randrange(3)]
                alleles[col] = base
            reads.append(ReadObs(fam.shared_chain, tuple(alleles)))
    return reads


def _clean_support_per_copy(fam: FamilyTruth, reads: list[ReadObs], K: int) -> float:
    """Order parameter of the identifiability theorem: the largest distinct error-free
    allele-vector support, spanning >= K columns, available PER TRUE COPY. A copy is
    recoverable when this support reaches min_reads (the split's candidate threshold).
    For each true copy we count reads whose observed (non-None) columns (a) all agree
    with that copy, (b) span >= K columns, then take the count of the single most-
    supported exact observed vector among those (the candidate the split would form)."""
    per_copy = []
    for tv in fam.allele_vectors:
        vec_counts: dict[tuple, int] = {}
        for r in reads:
            obs_cols = [i for i, a in enumerate(r.psv_alleles) if a is not None]
            if len(obs_cols) < K:
                continue
            if all(r.psv_alleles[i] == tv[i] for i in obs_cols):  # error-free & consistent
                vec_counts[r.psv_alleles] = vec_counts.get(r.psv_alleles, 0) + 1
        per_copy.append(max(vec_counts.values()) if vec_counts else 0)
    return sum(per_copy) / len(per_copy)


def recovered_copies(out: list[CopyIsoform]):
    """Identifiable copies emitted by the split (the recovered set)."""
    return [c for c in out if c.identifiable]


def fraction_recovered(fam: FamilyTruth, out: list[CopyIsoform]) -> float:
    """Fraction of TRUE copies correctly recovered: an emitted identifiable copy
    whose (non-None) allele vector is consistent with exactly one true copy AND
    matches it. We map each emitted copy to the unique true vector it agrees with."""
    emitted = recovered_copies(out)
    matched_truth = set()
    for c in emitted:
        hits = [
            ti for ti, tv in enumerate(fam.allele_vectors)
            if all(a is None or a == tv[i] for i, a in enumerate(c.allele_vector))
        ]
        if len(hits) == 1:
            matched_truth.add(hits[0])
    return len(matched_truth) / fam.n_copies


# ===========================================================================
# (1) SINGLE-CONFIG DEMO
# ===========================================================================

def run_demo():
    N, P, C, e, K = 3, 8, 10, 0.003, 2
    MIN_READS = 2
    rng = random.Random(DEMO_SEED)
    np.random.seed(DEMO_SEED)  # seeded for byte-identical reruns even though unused below
    fam = make_family(N, P, K, rng)
    reads = simulate_reads(fam, coverage=C, span_frac=1.0, err=e, rng=rng)

    # Read-coherence ALONE: group by chain only -> one transcript.
    chains = {r.intron_chain for r in reads}
    rc_only = len(chains)

    out = split_readchain_by_psv(reads, P, K, MIN_READS)
    rec = recovered_copies(out)
    frac = fraction_recovered(fam, out)

    print("=" * 64)
    print("DEMO  (N=%d copies, P=%d PSVs, C=%d reads/copy, e=%.3f, K=%d)" % (N, P, C, e, K))
    print("  shared intron chain (collapsed tandem):", fam.shared_chain)
    print("-" * 64)
    print("  read-coherence ALONE -> %d transcript(s) (the whole family collapses)" % rc_only)
    print("  read-coherence + PSV -> %d identifiable copies" % len(rec))
    print("  truth = %d copies; fraction recovered = %.3f" % (N, frac))
    print("-" * 64)
    for ci, c in enumerate(rec):
        vec = "".join(chr(b) if b is not None else "." for b in c.allele_vector)
        print("    copy %d: alleles=%s  reads=%d  identifiable=%s"
              % (ci, vec, c.read_count, c.identifiable))
    print("  TRUTH allele vectors:")
    for ti, tv in enumerate(fam.allele_vectors):
        print("    copy %d: %s" % (ti, "".join(chr(b) for b in tv)))
    print("=" * 64)
    return dict(N=N, P=P, C=C, e=e, K=K, rc_only=rc_only,
                recovered=len(rec), frac=frac)


# ===========================================================================
# (2) SWEEP: identifiability curve over (span fraction, coverage).
# ===========================================================================

def run_sweep():
    N, P, K = 3, 8, 2
    MIN_READS = 2
    err = 0.003

    span_fracs = [0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.875, 1.0]
    coverages = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
    REPLICATES = 40  # independent seeded families+reads per cell, averaged

    heat = np.zeros((len(span_fracs), len(coverages)))
    # The continuous governing quantity (the theorem's order parameter): the expected
    # number of error-FREE reads per copy that span >= K of that copy's columns AND
    # whose exact observed vector reaches >= min_reads support (i.e. the support that
    # can actually pin a true copy). When this crosses min_reads, recovery -> 1.0.
    support_axis = np.zeros((len(span_fracs), len(coverages)))

    for si, sf in enumerate(span_fracs):
        for ci, cov in enumerate(coverages):
            fr_sum = 0.0
            sup_sum = 0.0
            for rep in range(REPLICATES):
                # deterministic per-cell seed from fixed constants only.
                seed = SWEEP_SEED_BASE + si * 100000 + ci * 1000 + rep
                rng = random.Random(seed)
                fam = make_family(N, P, K, rng)
                reads = simulate_reads(fam, coverage=cov, span_frac=sf, err=err, rng=rng)
                out = split_readchain_by_psv(reads, P, K, MIN_READS)
                fr_sum += fraction_recovered(fam, out)
                # max over true copies of (count of reads exactly matching a true-copy
                # sub-span and spanning >= K of its columns), averaged per copy -> the
                # effective per-copy "clean spanning support" that the split can use.
                sup_sum += _clean_support_per_copy(fam, reads, K)
            heat[si, ci] = fr_sum / REPLICATES
            support_axis[si, ci] = sup_sum / REPLICATES

    # error-floor sweep: fix a generous geometry, vary e, to show the floor.
    errs = [0.0, 0.001, 0.003, 0.01, 0.03, 0.05, 0.08, 0.12, 0.2]
    cov_fixed, sf_fixed = 12, 1.0
    err_curve = []
    for ei, e in enumerate(errs):
        fr_sum = 0.0
        for rep in range(REPLICATES):
            seed = SWEEP_SEED_BASE + 9_000_000 + ei * 1000 + rep
            rng = random.Random(seed)
            fam = make_family(N, P, K, rng)
            reads = simulate_reads(fam, coverage=cov_fixed, span_frac=sf_fixed, err=e, rng=rng)
            out = split_readchain_by_psv(reads, P, K, MIN_READS)
            fr_sum += fraction_recovered(fam, out)
        err_curve.append(fr_sum / REPLICATES)

    return dict(
        N=N, P=P, K=K, err=err, min_reads=MIN_READS,
        span_fracs=span_fracs, coverages=coverages, replicates=REPLICATES,
        heat=heat, support_axis=support_axis,
        errs=errs, err_curve=err_curve, cov_fixed=cov_fixed, sf_fixed=sf_fixed,
    )


# ===========================================================================
# PLOT
# ===========================================================================

NAVY = "#22313f"
ORANGE = "#e8590c"
GREEN = "#1e8449"


def make_figure(demo, sweep, path):
    heat = sweep["heat"]
    span_fracs = sweep["span_fracs"]
    coverages = sweep["coverages"]

    # navy -> orange -> green colormap (merged -> boundary -> recovered).
    cmap = LinearSegmentedColormap.from_list("split", [NAVY, ORANGE, GREEN])

    fig = plt.figure(figsize=(13, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.34)

    # --- Panel A: heatmap fraction recovered vs (coverage, span fraction) ---
    axA = fig.add_subplot(gs[0, 0])
    im = axA.imshow(heat, origin="lower", aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
    axA.set_xticks(range(len(coverages)))
    axA.set_xticklabels(coverages, fontsize=8)
    axA.set_yticks(range(len(span_fracs)))
    axA.set_yticklabels(["%.0f%%" % (s * 100) for s in span_fracs], fontsize=8)
    axA.set_xlabel("coverage  (reads / copy)", fontsize=9)
    axA.set_ylabel("read span  (fraction of PSVs observed)", fontsize=9)
    axA.set_title("A. copies recovered  (N=%d, K=%d, e=%.3f)" % (sweep["N"], sweep["K"], sweep["err"]),
                  fontsize=10, color=NAVY)
    # boundary contour at 0.99 (recovery -> 1.0).
    axA.contour(heat, levels=[0.99], colors="white", linewidths=1.6)
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.03)
    cb.set_label("fraction recovered", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # --- Panel B: identifiability curve vs the theorem's order parameter ---
    #     x = per-copy clean (error-free, >= K-spanning) support for the candidate vector.
    #     The split forms a candidate copy at >= min_reads; recovery snaps to 1.0 there.
    axB = fig.add_subplot(gs[0, 1])
    x = sweep["support_axis"].ravel()
    y = heat.ravel()
    order = np.argsort(x)
    # bin and average to draw the empirical boundary curve.
    xb = x[order]
    yb = y[order]
    nb = 14
    edges = np.linspace(xb.min(), xb.max() + 1e-9, nb + 1)
    cx, cy = [], []
    for k in range(nb):
        m = (xb >= edges[k]) & (xb < edges[k + 1])
        if m.any():
            cx.append(xb[m].mean())
            cy.append(yb[m].mean())
    axB.scatter(xb, yb, s=12, color=ORANGE, alpha=0.35, edgecolors="none")
    axB.plot(cx, cy, "-", color=NAVY, lw=2.0, label="empirical boundary")
    axB.axvline(sweep["min_reads"], color=GREEN, lw=1.4, ls="--", alpha=0.85)
    axB.text(sweep["min_reads"] + 0.15, 0.30, " min_reads = %d\n (candidate threshold)" % sweep["min_reads"],
             color=GREEN, fontsize=7.5, va="center")
    axB.set_xlabel("clean $\\geq$K-spanning support / copy", fontsize=9)
    axB.set_ylabel("fraction of copies recovered", fontsize=9)
    axB.set_title("B. identifiability boundary", fontsize=10, color=NAVY)
    axB.set_ylim(-0.05, 1.08)
    axB.text(0.03, 0.10, "merged", color=NAVY, fontsize=8, transform=axB.transAxes)
    axB.text(0.60, 0.92, "recovered", color=GREEN, fontsize=8, transform=axB.transAxes)
    for s in axB.spines.values():
        s.set_color("#cccccc")

    # --- Panel C: error floor ---
    axC = fig.add_subplot(gs[0, 2])
    axC.plot(sweep["errs"], sweep["err_curve"], "-o", color=ORANGE, ms=5, lw=1.8)
    axC.axvline(0.003, color=GREEN, lw=1.2, ls="--", alpha=0.8)
    axC.text(0.003, 0.12, " HiFi e=0.003", color=GREEN, fontsize=8, rotation=90, va="bottom")
    axC.set_xlabel("per-base error e", fontsize=9)
    axC.set_ylabel("fraction recovered", fontsize=9)
    axC.set_title("C. error floor  (C=%d, span=full)" % sweep["cov_fixed"], fontsize=10, color=NAVY)
    axC.set_ylim(-0.05, 1.08)
    for s in axC.spines.values():
        s.set_color("#cccccc")

    fig.suptitle(
        "Read-coherence + PSV copy-split in the collapsed-tandem regime  "
        "(read-coherence alone -> %d transcript; +PSV -> %d copies)"
        % (demo["rc_only"], demo["recovered"]),
        fontsize=11, color=NAVY, y=1.0,
    )
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ===========================================================================
# Boundary summary
# ===========================================================================

def summarize_boundary(sweep):
    heat = sweep["heat"]
    cov = sweep["coverages"]
    sf = sweep["span_fracs"]
    # smallest coverage achieving full recovery at full span.
    full_span_row = heat[-1]  # span_frac == 1.0
    cov_full = None
    for j, v in enumerate(full_span_row):
        if v >= 0.99:
            cov_full = cov[j]
            break
    # smallest span achieving full recovery at high coverage.
    high_cov_col = heat[:, -1]  # coverage == max
    sf_full = None
    for i, v in enumerate(high_cov_col):
        if v >= 0.99:
            sf_full = sf[i]
            break
    return cov_full, sf_full


if __name__ == "__main__":
    # Master seed (constant). Per-config seeds are derived from fixed constants too.
    random.seed(MASTER_SEED)
    np.random.seed(MASTER_SEED)

    demo = run_demo()
    sweep = run_sweep()
    make_figure(demo, sweep, "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_split_sim.png")

    cov_full, sf_full = summarize_boundary(sweep)
    print()
    print("SWEEP boundary (N=%d, P=%d, K=%d, e=%.3f, min_reads=%d, %d reps/cell):"
          % (sweep["N"], sweep["P"], sweep["K"], sweep["err"], sweep["min_reads"], sweep["replicates"]))
    print("  full recovery at full span needs coverage >= %s reads/copy" % cov_full)
    print("  full recovery at max coverage needs span >= %s of PSVs"
          % ("%.0f%%" % (sf_full * 100) if sf_full is not None else "n/a"))
    print("  error floor: recovered at e=0.003 -> %.3f ; collapses by e=%.2f -> %.3f"
          % (sweep["err_curve"][sweep["errs"].index(0.003)],
             sweep["errs"][-1], sweep["err_curve"][-1]))
    print("  figure -> bench/copy_split_sim.png")
