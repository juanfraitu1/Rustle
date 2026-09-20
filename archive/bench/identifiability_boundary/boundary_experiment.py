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

# ----------------------------------------------------------------------------
# EXTENSION 1: structure as a phasing channel (skip-junction bridging).
# Two copies differ by SPLICING, not (only) SNPs:
#     copy A:  e1 - S - e2   (INCLUDES the shared exon S)
#     copy B:  e1 ---- e2    (SKIPS S; uses an e1->e2 skip junction)
# e1 and e2 can be nucleotide-identical between copies — the ONLY signal is the
# junction. A copy-B read that covers the e1->e2 skip junction simultaneously
# (a) IDENTIFIES copy B and (b) LINKS the e1 and e2 flanks on one molecule WITHOUT
# spanning S's length (it skips S). So junction-linkage is independent of L_S,
# whereas SNP-only linkage needs a read to span all L_S bp. This recovers
# identifiability exactly where the SNP-only channel collapses.
# (No existing tool uses junctions as an INPUT linking channel — longcallR tests
#  allele-specific junctions only DOWNSTREAM, conditioned on a SNP haplotype.)
def run_structural_sweep(R=100, coverage=30, n_trials=400, seed=2, anchor=8, flank=120):
    """Compare P(family identifiable) for SNP-only vs structural(skip) channels,
    sweeping L_S/R. flank = length of e1 (and e2) around the junction."""
    rng = random.Random(seed)
    L_S_vals = list(range(10, 3*R + 1, max(1, R//10)))
    rows = []
    for L_S in L_S_vals:
        # SNP-only: identifiable iff some read spans the L_S-wide shared region.
        snp_id = 0
        # structural: identifiable iff some copy-B read covers the e1->e2 skip
        # junction (footprint has >= anchor bp in e1 AND in e2). On the SPLICED
        # B transcript (length 2*flank), this is independent of L_S.
        struct_id = 0
        m_gap = L_S + 1                      # m1..m2 distance for the SNP channel
        max_start_snp = max(1, (2*R + m_gap) - R)
        spliced_len = 2 * flank              # B's transcript with S excised
        max_start_str = max(1, spliced_len - R)
        jpos = flank                         # e1->e2 junction at the e1/e2 boundary
        for _ in range(n_trials):
            # SNP channel: 2*coverage reads; any spanning S?
            span = any(rng.randrange(0, max_start_snp) <= R - m_gap
                       for _ in range(2*coverage)) if R >= m_gap else False
            snp_id += span
            # structural channel: copy-B reads; any covering the skip junction
            # with >= anchor bp on each side?
            jcov = False
            for _ in range(coverage):
                s = rng.randrange(0, max_start_str)
                if s + anchor <= jpos and s + R - anchor >= jpos:
                    jcov = True; break
            struct_id += jcov
        rows.append({"L_S": L_S, "L_S_over_R": L_S / R,
                     "snp_only_identifiable": snp_id / n_trials,
                     "structural_identifiable": struct_id / n_trials})
    return rows

# ----------------------------------------------------------------------------
# EXTENSION 2: unknown copy number K. In the CUT regime (L_S > R, no spanning
# read), with K copies you observe K left alleles and K right alleles as
# marginals but cannot pair them. A linkage-blind method commits to one of K!
# pairings; the expected number of CORRECT pairs is 1 (fixed points of a random
# permutation), so it emits ~K-1 chimeric isoforms -> phantom FRACTION (K-1)/K,
# rising to 1 as K grows. The certified rule abstains -> 0, and reports only that
# ">= K copies share S" (the count itself is non-identifiable without spanning reads).
def run_k_sweep(K_vals=range(2, 9), n_trials=2000, seed=3):
    rng = random.Random(seed)
    rows = []
    for K in K_vals:
        greedy_phantom_frac = 0.0
        for _ in range(n_trials):
            perm = list(range(K)); rng.shuffle(perm)            # greedy's guessed pairing
            correct = sum(perm[i] == i for i in range(K))       # fixed points
            greedy_phantom_frac += (K - correct) / K
        rows.append({"K": K,
                     "greedy_phantom_fraction": greedy_phantom_frac / n_trials,
                     "analytic_(K-1)/K": (K - 1) / K,
                     "certified_phantom_fraction": 0.0})
    return rows

def make_extended_figure(R, rows_boundary, rows_struct, rows_k, path):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    # Panel A: the boundary (SNP-only)
    ax = axes[0]
    x = [r["L_S_over_R"] for r in rows_boundary]
    ax.plot(x, [r["greedy_phantom_rate"] for r in rows_boundary], "o-", color="#c0392b", ms=4, label="greedy: fabricates")
    ax.plot(x, [r["certified_phantom_rate"] for r in rows_boundary], "s-", color="#27ae60", ms=4, label="certified: 0 fabrication")
    ax.plot(x, [r["spanning_fraction"] for r in rows_boundary], "--", color="#7f8c8d", label="P(read spans S)")
    ax.axvline(1.0, color="k", ls=":", lw=1)
    ax.set_title("A. The boundary (SNP marks)\nidentifiability collapses at L_S = R")
    ax.set_xlabel("L_S / R"); ax.set_ylabel("rate"); ax.set_ylim(-0.03, 1.03); ax.grid(alpha=0.25); ax.legend(fontsize=7, loc="center right")
    # Panel B: structure rescues identifiability
    ax = axes[1]
    xs = [r["L_S_over_R"] for r in rows_struct]
    snp = [r["snp_only_identifiable"] for r in rows_struct]
    st  = [r["structural_identifiable"] for r in rows_struct]
    ax.plot(xs, snp, "--", color="#7f8c8d", label="SNP-only channel")
    ax.plot(xs, st,  "^-", color="#2980b9", ms=4, label="+ structural channel (skip junction)")
    ax.fill_between(xs, snp, st, where=[b>=a for a,b in zip(snp,st)], color="#2980b9", alpha=0.15)
    ax.axvline(1.0, color="k", ls=":", lw=1)
    ax.set_title("B. Structure as a phasing channel\na skip-read bridges S without spanning it")
    ax.set_xlabel("L_S / R"); ax.set_ylabel("P(family identifiable)"); ax.set_ylim(-0.03, 1.03); ax.grid(alpha=0.25); ax.legend(fontsize=7, loc="center right")
    # Panel C: unknown K makes the cut catastrophic
    ax = axes[2]
    ks = [r["K"] for r in rows_k]
    ax.plot(ks, [r["greedy_phantom_fraction"] for r in rows_k], "o-", color="#c0392b", ms=5, label="greedy (cut regime)")
    ax.plot(ks, [r["analytic_(K-1)/K"] for r in rows_k], ":", color="#c0392b", alpha=0.6, label="analytic (K-1)/K")
    ax.plot(ks, [r["certified_phantom_fraction"] for r in rows_k], "s-", color="#27ae60", ms=5, label="certified")
    ax.set_title("C. Unknown copy number K\nthe more copies share S, the worse the cut")
    ax.set_xlabel("number of copies K"); ax.set_ylabel("fraction of emitted isoforms that are phantom"); ax.set_ylim(-0.03, 1.03); ax.grid(alpha=0.25); ax.legend(fontsize=7, loc="center right")
    fig.suptitle("Identifiability of multi-copy gene-family assembly: the boundary, the structural channel that extends it, and how unknown K sharpens it", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=150)

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

    # --- extensions ---
    rows_struct = run_structural_sweep(R=R)
    rows_k = run_k_sweep()
    ext_png = os.path.join(HERE, "boundary_extended.png")
    make_extended_figure(R, rows, rows_struct, rows_k, ext_png)
    with open(os.path.join(HERE, "structural_data.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_struct[0].keys())); w.writeheader(); w.writerows(rows_struct)
    with open(os.path.join(HERE, "k_data.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_k[0].keys())); w.writeheader(); w.writerows(rows_k)
    print(f"wrote {ext_png}")

    # structural headline: identifiability above the SNP boundary (L_S > R)
    sa = [r for r in rows_struct if r["L_S_over_R"] >= 1.5]
    snp_above = np.mean([r["snp_only_identifiable"] for r in sa])
    str_above = np.mean([r["structural_identifiable"] for r in sa])
    print(f"\nEXT-1 (structure as a channel), at L_S > R: "
          f"P(identifiable) SNP-only={snp_above:.3f}  vs  +structural(skip)={str_above:.3f}  "
          f"-> structure recovers identifiability the SNP channel loses.")
    # K headline
    print("EXT-2 (unknown K), cut regime: greedy phantom fraction by K: " +
          ", ".join(f"K={r['K']}:{r['greedy_phantom_fraction']:.2f}" for r in rows_k) +
          "  | certified: 0 for all K.")

if __name__ == "__main__":
    main()
