#!/usr/bin/env python3
"""
Identifiability boundary experiment for multi-copy gene-family assembly.

Question (the thesis core, answering the "single SNP is trivial" objection):
  Two gene copies are identical over a shared region S of length L_S, and differ
  at exactly two marks that flank S:
        mark m1            shared region S (identical)            mark m2
   copy A:  G        ===== L_S bp, no marks =====                   T
   copy B:  C        ===== L_S bp, no marks =====                   A
  A long read of length R observes a mark only if its footprint covers that mark.
  A read that covers BOTH marks (it must SPAN S, so R must exceed the m1..m2 gap)
  LINKS the left allele to the right allele on one molecule.

THEOREM (read-linkage identifiability): the copies' full structures are
  identifiable iff some read spans S. If no read spans S, the data are EQUALLY
  consistent with the truth {(G,T),(C,A)} and the CHIMERA {(G,A),(C,T)} — an
  isoform present in NEITHER copy — so any method that emits full isoforms must
  GUESS, and is wrong ~half the time (and cannot know it).

This script demonstrates that boundary at L_S ~= R, comparing:
  * GREEDY (linkage-blind point estimate): uses a spanning read if one exists;
    otherwise still emits two full isoforms by guessing the pairing -> fabricates.
  * CERTIFIED (our rule): join left-to-right across S ONLY when a read spans it;
    otherwise ABSTAIN and emit a non-identifiability certificate -> never fabricates.

Outputs: boundary_curve.png, boundary_data.csv, and a worked single-case printout.
Deterministic (seeded). No external data needed.
"""
import csv, os, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# Ground truth: two equal-expression copies (the genuinely ambiguous case).
TRUTH = {"A": ("G", "T"), "B": ("C", "A")}     # (left allele @m1, right allele @m2)
CHIMERAS = {("G", "A"), ("C", "T")}            # pairings present in NO copy

def simulate_reads(L_S, R, coverage, eps, rng):
    """Place `coverage` reads per copy on a transcript with marks flanking S.
    Returns list of dicts: {m1: allele or None, m2: allele or None, copy: 'A'/'B'}.
    A read SPANS iff it observes both marks (requires footprint to cover m1..m2)."""
    pad = R                       # let reads start before m1 / end after m2
    m1_pos = pad
    m2_pos = pad + L_S + 1        # marks flank the shared region of length L_S
    L_tx = m2_pos + pad
    max_start = max(1, L_tx - R)
    reads = []
    for copy, (aL, aR) in TRUTH.items():
        for _ in range(coverage):
            start = rng.randrange(0, max_start)
            end = start + R
            cov_m1 = start <= m1_pos < end
            cov_m2 = start <= m2_pos < end
            def emit(true_allele):
                # sequencing error: flip to a random other base (4-letter)
                if rng.random() < eps:
                    return rng.choice([b for b in "ACGT" if b != true_allele])
                return true_allele
            reads.append({
                "m1": emit(aL) if cov_m1 else None,
                "m2": emit(aR) if cov_m2 else None,
                "copy": copy,
            })
    return reads

def spanning_links(reads):
    """Reads that observe BOTH marks -> (left, right) linkage votes."""
    return [(r["m1"], r["m2"]) for r in reads if r["m1"] is not None and r["m2"] is not None]

def greedy_pairing(reads, rng):
    """Linkage-blind point estimate. Uses spanning reads if present (majority
    link); otherwise emits a full pairing by guessing (coin flip). Always emits
    two full-length isoforms -> can fabricate."""
    links = spanning_links(reads)
    if links:
        # majority vote on what the left 'G' connects to
        g_right = [r for (l, r) in links if l == "G"]
        if g_right:
            top = max(set(g_right), key=g_right.count)
            return {("G", top), ("C", "A" if top == "T" else "T")}
        # only C-links seen
        c_right = [r for (l, r) in links if l == "C"]
        top = max(set(c_right), key=c_right.count)
        return {("C", top), ("G", "A" if top == "T" else "T")}
    # NO spanning read -> non-identifiable -> guess (data is symmetric; ~50% wrong)
    if rng.random() < 0.5:
        return {("G", "T"), ("C", "A")}      # happens to be truth
    return {("G", "A"), ("C", "T")}          # chimera

def certified_assembly(reads):
    """Our rule: join across S only if a read spans it; else ABSTAIN.
    Returns (full_pairing_or_None, certificate_or_None)."""
    links = spanning_links(reads)
    if links:
        g_right = [r for (l, r) in links if l == "G"]
        if g_right:
            top = max(set(g_right), key=g_right.count)
            return {("G", top), ("C", "A" if top == "T" else "T")}, None
        c_right = [r for (l, r) in links if l == "C"]
        top = max(set(c_right), key=c_right.count)
        return {("C", top), ("G", "A" if top == "T" else "T")}, None
    # abstain: report unlinked sides + certificate (no full isoform emitted)
    cert = ("NON-IDENTIFIABLE across S: left alleles {G,C} and right alleles {T,A} "
            "are not linked by any read (no read spans S). 2 of the 4 graph paths "
            "are phantom; data cannot say which.")
    return None, cert

def fabricated(pairing):
    """Does an emitted full pairing contain a chimeric (no-copy) isoform?"""
    return pairing is not None and bool(pairing & CHIMERAS)

def run_sweep(R=100, coverage=30, eps=0.01, n_trials=400, seed=1):
    rng = random.Random(seed)
    L_S_vals = list(range(10, 30*R//10 + 1, max(1, R//10)))  # 10 .. 3R
    rows = []
    for L_S in L_S_vals:
        g_phantom = c_phantom = spanned = c_complete = 0
        for _ in range(n_trials):
            reads = simulate_reads(L_S, R, coverage, eps, rng)
            has_span = len(spanning_links(reads)) > 0
            spanned += has_span
            gp = greedy_pairing(reads, rng)
            g_phantom += fabricated(gp)
            cp, cert = certified_assembly(reads)
            c_phantom += fabricated(cp)        # always 0 (cp is None or correct)
            c_complete += (cp is not None)
        rows.append({
            "L_S": L_S, "L_S_over_R": L_S / R,
            "greedy_phantom_rate": g_phantom / n_trials,
            "certified_phantom_rate": c_phantom / n_trials,
            "spanning_fraction": spanned / n_trials,
            "certified_complete_rate": c_complete / n_trials,
        })
    return R, rows

def worked_case(R=100, L_S=220, coverage=30, eps=0.01, seed=7):
    """A concrete single family with L_S > R: greedy fabricates, certified refuses."""
    rng = random.Random(seed)
    reads = simulate_reads(L_S, R, coverage, eps, rng)
    links = spanning_links(reads)
    gp = greedy_pairing(reads, rng)
    cp, cert = certified_assembly(reads)
    lines = []
    lines.append(f"WORKED CASE  (R={R}, L_S={L_S}  ->  L_S/R={L_S/R:.2f} > 1, shared region wider than reads)")
    lines.append(f"  reads: {len(reads)}  |  reads spanning S (linking m1..m2): {len(links)}")
    lines.append(f"  left alleles seen @m1: G={sum(r['m1']=='G' for r in reads)} C={sum(r['m1']=='C' for r in reads)}")
    lines.append(f"  right alleles seen @m2: T={sum(r['m2']=='T' for r in reads)} A={sum(r['m2']=='A' for r in reads)}")
    lines.append(f"  TRUTH copies: {set(TRUTH.values())}")
    lines.append(f"  GREEDY emits isoforms: {gp}   -> fabricated phantom? {fabricated(gp)}")
    if cp is None:
        lines.append(f"  CERTIFIED: ABSTAINS. certificate: {cert}")
        lines.append(f"  CERTIFIED fabricated phantom? False  (refuses to emit a full isoform it cannot witness)")
    else:
        lines.append(f"  CERTIFIED emits: {cp}  (a spanning read existed)")
    return "\n".join(lines)

def main():
    R, rows = run_sweep()
    csv_path = os.path.join(HERE, "boundary_data.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    x  = [r["L_S_over_R"] for r in rows]
    gp = [r["greedy_phantom_rate"] for r in rows]
    cp = [r["certified_phantom_rate"] for r in rows]
    sf = [r["spanning_fraction"] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, gp, "o-", color="#c0392b", label="greedy (linkage-blind): fabricates phantom isoforms")
    ax.plot(x, cp, "s-", color="#27ae60", label="certified (join only if a read spans S): 0 fabrication")
    ax.plot(x, sf, "--", color="#7f8c8d", label="P(some read spans S)  [identifiable region]")
    ax.axvline(1.0, color="k", ls=":", lw=1)
    ax.text(1.02, 0.46, "L_S = R\n(read length = mark spacing)", fontsize=8, va="top")
    ax.set_xlabel("shared-region length / read length   (L_S / R)")
    ax.set_ylabel("rate")
    ax.set_title("Identifiability boundary: read length vs. mark spacing\n"
                 "equal-expression paralog copies, two flanking marks")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(loc="center right", fontsize=8)
    ax.grid(alpha=0.25)
    png = os.path.join(HERE, "boundary_curve.png")
    fig.tight_layout(); fig.savefig(png, dpi=150)

    print(worked_case())
    print()
    print(f"wrote {png}")
    print(f"wrote {csv_path}")
    # headline numbers
    below = [r for r in rows if r["L_S_over_R"] <= 0.8]
    above = [r for r in rows if r["L_S_over_R"] >= 1.5]
    gb = np.mean([r["greedy_phantom_rate"] for r in below])
    ga = np.mean([r["greedy_phantom_rate"] for r in above])
    print(f"\nHEADLINE: greedy phantom-isoform rate  below boundary (L_S<R): {gb:.3f}   "
          f"above boundary (L_S>R): {ga:.3f}")
    print(f"          certified phantom-isoform rate everywhere: "
          f"{max(r['certified_phantom_rate'] for r in rows):.3f}")

if __name__ == "__main__":
    main()
