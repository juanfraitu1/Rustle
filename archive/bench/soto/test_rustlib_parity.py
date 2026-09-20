#!/usr/bin/env python3
"""Regression test for rustlib.py — run this before trusting any bench number.

TWO JOBS.

(1) CROSS-IMPLEMENTATION AGREEMENT. The Python layer exists so core principles can be checked quickly,
    which is only useful if it computes the SAME rule as the shipped engine.
    TRUE RUST PARITY IS NOW TESTED (2026-08-07). `RUSTLE_ER_EDGE_DUMP=<prefix> gw_family_catalog ...`
    writes the returned E_r edge set, the node list, the effective parameters and the Rust's OWN
    minimap2 PAF; the last section below re-derives that edge set with rustlib and demands an EXACT
    match. Running rustlib over the Rust's own PAF is what makes the comparison about the EDGE RULE
    rather than about the alignment invocation.
    ⚠ Exact agreement on a fixture is only as strong as the fixture's separating power. Each dump is
    therefore also run through a WRONG-DENOMINATOR counterfactual, and the test prints whether that
    counterfactual actually moved the edge set. A fixture that agrees under both denominators has
    proven nothing about the denominator.

(2) REGRESSION ON THE FIVE MISTAKES that produced quotable wrong numbers in the 2026-08-06/07 session.
    Each has a test that FAILS if the mistake is reintroduced:
      * single-record vs aggregated coverage      (PAF figure contradicted the graph figure)
      * uniform vs two-tier identity floor        (retracted RNA densities 0.547 / 0.698)
      * argmax best-hit vs qualifying-edge count  ("6/19 orthologs" when it was 19/19)
      * node-name convention                      (a false "no alignment records at all")
      * reads overlapping vs contained in a window(neighbour's transcripts became the representative)
      * the SINGLE-PATH CEILING denominator       (five denominators gave 0.33-1.00 for one model)

Stdlib only. Exits non-zero on failure. Skips cleanly when a fixture is absent.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rustlib as R  # noqa: E402

FAILED = []
SKIPPED = []


def check(name, cond, detail=""):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}   {detail}")
        FAILED.append(name)


def skip(name, why):
    print(f"  skip  {name}   ({why})")
    SKIPPED.append(name)


# ---------------------------------------------------------------- node naming
print("node naming (mistake 4: a false 'no alignment records at all')")
check("canon_node is 0-based", R.canon_node("chr16:15117337-15126130") == "L~chr16_15117336_15126130",
      R.canon_node("chr16:15117337-15126130"))
check("parse_node round-trips both forms",
      R.parse_node("L~chr16_15117336_15126130") == ("chr16", 15117336, 15126130)
      and R.parse_node("chr16:15117337-15126130") == ("chr16", 15117336, 15126130))

# ---------------------------------------------------------------- edge rule
print("\nedge rule (mistakes 1 and 2)")
PAF = """A\t1000\t0\t900\t+\tB\t1000\t0\t900\t880\t900\t60\tde:f:0.02
A\t1000\t0\t300\t+\tC\t1000\t0\t300\t295\t300\t60\tde:f:0.02
A\t1000\t400\t700\t+\tC\t1000\t400\t700\t295\t300\t60\tde:f:0.02
A\t1000\t0\t900\t+\tD\t1000\t0\t900\t600\t900\t60\tde:f:0.35
"""
with tempfile.NamedTemporaryFile("w", suffix=".paf", delete=False) as fh:
    fh.write(PAF)
    p = fh.name
try:
    single = R.er_edges([(p, 0.60)], aggregate=False)
    agg = R.er_edges([(p, 0.60)], aggregate=True)
    # A-B: one record at coverage 0.90 -> edge under both
    check("single-record admits a clean pair", ("A", "B") in single)
    # A-C: two records of 0.30 each. Single-record REJECTS (0.30 < 0.50); aggregate ACCEPTS (0.60).
    check("single-record REJECTS scattered blocks (shipped)", ("A", "C") not in single,
          "coverage must not be summed under the shipped rule")
    check("aggregate ACCEPTS the same scattered blocks", ("A", "C") in agg,
          "aggregate=True must union the intervals")
    check("aggregate is not the default", single != agg)
    # A-D: identity 0.65 -> passes a 0.60 floor, fails a 0.80 floor. The floors must not collapse.
    lo = R.er_edges([(p, 0.60)])
    hi = R.er_edges([(p, 0.80)])
    check("identity floor is honoured per tier", ("A", "D") in lo and ("A", "D") not in hi,
          f"lo={sorted(lo)} hi={sorted(hi)}")
    two_tier = R.er_edges([(p, 0.80), (p, 0.60)])
    check("multi-tier UNIONs rather than collapsing", ("A", "D") in two_tier)
    check("SHIPPED_TIERS is sensitive-only at 0.60 (asm20 retired 2026-08-07)",
          R.SHIPPED_TIERS == (("sensitive", 0.60),), str(R.SHIPPED_TIERS))
finally:
    os.unlink(p)

# ---------------------------------------------------------------- components / scoring
print("\ncomponents and scoring (mistake 3: argmax best-hit)")
nodes = ["chr1:1-100", "chr1:201-300", "chr1:401-500", "chr1:601-700"]
edges = {("chr1:1-100", "chr1:201-300"), ("chr1:201-300", "chr1:401-500")}
C = R.components(nodes, edges)
check("components largest-first with singletons kept", [len(c) for c in C] == [3, 1], str(C))
check("density of a 3-node path", abs(R.density(C[0], edges) - 2 / 3) < 1e-9)
check("density of a singleton is NaN not 0", R.density(C[1], edges) != R.density(C[1], edges))
members = [("chr1", 0, 100, "G1"), ("chr1", 200, 300, "G2"), ("chr1", 400, 500, "G3")]
sc = R.score_family(C[0], members, edges)
check("recall counts every qualifying member, not the best hit", sc["recall"] == 1.0, str(sc))
check("purity is loci-with-a-member / loci", sc["purity"] == 1.0, str(sc))

# ---------------------------------------------------------------- BAM-backed (skip if absent)
print("\nBAM primitives (mistake 5: overlapping vs contained)")
BAM = os.environ.get("RUSTLE_TEST_BAM", "/home/juanfra/winloci_scratch/chr16_sub.bam")
if not os.path.exists(BAM):
    skip("sec_frac / read_exons", f"no BAM at {BAM}")
else:
    sf = R.sec_frac(BAM, "chr16", 29663350, 29688504)          # NPIPB11
    check("sec_frac keeps secondaries (NPIPB11 is duplicated)", sf["sec_frac"] > 0.8, str(sf))
    check("sec_frac counts primaries too", sf["primary"] > 0, str(sf))
    inside = R.read_exons(BAM, "chr16", 29016080, 29053452, contained=True)
    outside = R.read_exons(BAM, "chr16", 29016080, 29053452, contained=False)
    check("contained=True stays inside the window",
          all(a >= 29016080 and b <= 29053452 for a, b in inside), str(inside[:3]))
    check("contained=False admits more than contained=True (the neighbour leak)",
          len(outside) >= len(inside), f"{len(outside)} vs {len(inside)}")

# ---------------------------------------------------------------- cross-implementation check
print("\ncross-implementation check")
# ⚠ THIS IS NOT RUST PARITY. edges.tsv is emitted by bench's OWN python harness
# (dna_seed_family_bench.py), so comparing against it is python-vs-python. True parity needs a
# gw_family_catalog run whose edge set is dumped; the Rust currently emits copies/families, not edges.
# Marked explicitly so a passing run is never mistaken for verification against the shipped engine.
HARNESS_PAF = os.environ.get("RUSTLE_PARITY_PAF",
                             "/home/juanfra/winloci_scratch/seedfam/dnapr/meeting/loci_ava_sens.paf")
HARNESS_EDGES = os.environ.get("RUSTLE_PARITY_EDGES",
                               "/home/juanfra/winloci_scratch/seedfam/dnapr/meeting/edges.tsv")
if not (os.path.exists(HARNESS_PAF) and os.path.exists(HARNESS_EDGES)):
    skip("harness cross-check", "no harness edges.tsv / PAF fixture on disk")
else:
    other = set()
    with open(HARNESS_EDGES) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 2 or f[0] == "u":
                continue
            other.add((f[0], f[1]) if f[0] < f[1] else (f[1], f[0]))
    nodes = {n for e in other for n in e}
    py = R.er_edges([(HARNESS_PAF, 0.60)], nodes=nodes)
    inter = py & other
    j = len(inter) / max(len(py | other), 1)
    # The harness AGGREGATES identity and coverage across records; rustlib uses the SHIPPED
    # single-record rule. They are expected to differ slightly, and the direction is informative:
    # aggregation should admit MORE, never fewer. A large gap either way means one of them drifted.
    check("the two implementations agree to within 1%", j >= 0.99,
          f"jaccard={j:.4f} python-only={len(py - other)} harness-only={len(other - py)}")
    print(f"        harness={len(other)} rustlib={len(py)} shared={len(inter)} jaccard={j:.4f}")
    print(f"        (differences are the aggregated-vs-single-record rule, mistake 1 — "
          f"rustlib follows the SHIPPED single-record rule)")

# ---------------------------------------------------------------- per-record floors (2026-08-07 bug)
print("\nper-record floors (paf_pairs used to argmax by coverage BEFORE testing the floors)")
# rec1 has the highest coverage but FAILS identity; rec2 clears BOTH floors. The Rust unions over
# records (denovo_pipeline.rs:3687) so this is an edge. Selecting the argmax-coverage record first
# and only then testing the floors loses it. Row 3 carries NO de:f: at all: the Rust falls back to
# nmatch/blocklen (denovo_pipeline.rs:3635), so dropping the record loses that edge too.
PR = """0\t1000\t0\t900\t+\t1\t1000\t0\t900\t550\t900\t60\tde:f:0.4500
0\t1000\t0\t700\t+\t1\t1000\t0\t700\t665\t700\t60\tde:f:0.0500
2\t1000\t0\t800\t+\t3\t1000\t0\t800\t784\t800\t60
"""
with tempfile.NamedTemporaryFile("w", suffix=".paf", delete=False) as fh:
    fh.write(PR)
    p = fh.name
try:
    e = R.edges_from_paf(p, 0.60, 0.50)
    check("a lower-coverage record that clears both floors still makes an edge", ("0", "1") in e,
          f"got {sorted(e)}; the argmax-by-coverage record fails identity 0.55 < 0.60")
    check("a record with no de:f: falls back to nmatch/blocklen, not dropped", ("2", "3") in e,
          f"got {sorted(e)}")
    pp = R.paf_pairs(p, 0.60, 0.50)
    check("the reported exemplar is the highest-coverage PASSING record",
          abs(pp[("0", "1")]["coverage"] - 0.70) < 1e-9
          and abs(pp[("0", "1")]["identity"] - 0.95) < 1e-9, str(pp[("0", "1")]))
finally:
    os.unlink(p)

# ---------------------------------------------------------------- TRUE Rust parity
print("\nTRUE Rust parity (RUSTLE_ER_EDGE_DUMP vs rustlib, on the Rust's OWN PAF)")
# Each entry is the <prefix> passed to RUSTLE_ER_EDGE_DUMP. Override with a colon-separated list in
# RUSTLE_PARITY_DUMPS. Missing stems skip; a stem that exists but is EMPTY is a FAILURE, never a
# "0 edges" finding (RUSTLE_SHARED_EXON=1 returns before the dump, and a non-zero minimap2 exit
# yields an empty edge set — the Rust prints an explicit stderr warning for the latter).
DUMPS = os.environ.get("RUSTLE_PARITY_DUMPS", ":".join([
    "/home/juanfra/winloci_scratch/erdump/f78/f78",        # 78 genomic loci, 336,061 PAF records
    "/home/juanfra/winloci_scratch/erdump/mix",            # 24 genomic loci
    "/home/juanfra/winloci_scratch/erdump/rna",            # --homology-primary, spliced substrate
    "/home/juanfra/winloci_scratch/erdump/sumcov/sc",      # RUSTLE_ER_SUM_COVERAGE=1
    "/home/juanfra/winloci_scratch/erdump/longer/lg",      # RUSTLE_ER_COVERAGE_LONGER=1
])).split(":")


def _pairkey(a, b):
    """⚠ Rep ids in the dump are INTEGER strings. Sorting them LEXICOGRAPHICALLY ('10' < '2')
    disagrees with the Rust's numeric (rep_i, rep_j) ordering and manufactures a fake divergence —
    that alone produced 318 rust-only and 318 python-only edges on the 78-locus fixture before it
    was spotted. Sort numerically whenever both ids parse as integers."""
    try:
        x, y = int(a), int(b)
    except ValueError:
        return (a, b) if a < b else (b, a)
    return (x, y) if x < y else (y, x)


ran_any = False
for stem in [s for s in DUMPS if s]:
    tag = os.path.basename(stem)
    if not os.path.exists(stem + ".edges.tsv"):
        skip(f"Rust parity [{tag}]", f"no dump at {stem}.edges.tsv")
        continue
    prm = dict(ln.split("\t", 1)
               for ln in open(stem + ".params.tsv").read().splitlines()[1:] if "\t" in ln)
    rows = open(stem + ".edges.tsv").read().splitlines()[1:]
    nrows = open(stem + ".nodes.tsv").read().splitlines()[1:]
    pafs = [f for f in os.listdir(os.path.dirname(stem) or ".")
            if f.startswith(os.path.basename(stem) + ".er.") and f.endswith(".paf")]
    pafs = [os.path.join(os.path.dirname(stem) or ".", f) for f in pafs]
    if not rows or not nrows or not pafs:
        check(f"Rust parity [{tag}] dump is non-empty", False,
              f"{len(rows)} edges / {len(nrows)} nodes / {len(pafs)} PAFs — an empty dump is a BUG, "
              f"not a zero-edge result")
        continue
    ran_any = True
    rust = {_pairkey(r.split("\t")[0], r.split("\t")[1]) for r in rows}
    floor = float(prm["sensitive_identity"])
    mincov = float(prm["min_coverage"])
    denom = "max" if prm["coverage_denominator"] == "span(max)" else "min"
    agg = prm.get("env.RUSTLE_ER_SUM_COVERAGE", "<unset>") not in ("<unset>", "0", "")
    tiers = [(p, floor) for p in sorted(pafs)]
    py = {_pairkey(*k) for k in R.er_edges(tiers, min_coverage=mincov, denom=denom, aggregate=agg)}
    inter = rust & py
    j = len(inter) / max(len(rust | py), 1)
    check(f"Rust parity [{tag}] edge sets are IDENTICAL", rust == py,
          f"rust={len(rust)} python={len(py)} shared={len(inter)} rust_only={len(rust - py)} "
          f"python_only={len(py - rust)} jaccard={j:.4f} "
          f"first_rust_only={sorted(rust - py)[:3]} first_python_only={sorted(py - rust)[:3]}")
    # Separating power: does the WRONG denominator actually move this fixture?
    wrong = {_pairkey(*k) for k in R.er_edges(tiers, min_coverage=mincov,
                                              denom=("min" if denom == "max" else "max"),
                                              aggregate=agg)}
    disc = "DISCRIMINATES" if wrong != rust else "⚠ does NOT discriminate"
    print(f"        {tag}: nodes={len(nrows)} rust={len(rust)} python={len(py)} jaccard={j:.4f} "
          f"| denom={denom} agg={agg} | wrong-denominator {disc} ({len(wrong)} edges)")

check("at least one Rust dump was actually diffed", ran_any,
      "every dump stem was missing — this section proved NOTHING; regenerate with "
      "RUSTLE_ER_EDGE_DUMP=<prefix> gw_family_catalog ...")

# ---------------------------------------------------------------- single-path ceiling (mistake 6)
print("\nsingle-path ceiling (mistake 6: five denominators gave 0.33-1.00 for the SAME model)")

# (a) HAND TOY. Six junctions, worked out by hand and confirmed against full enumeration below.
#     j1=(10,20) j2=(15,25) j3=(20,30) j4=(30,40) j5=(35,50) j6=(60,70)
#     touch_ok=True  conflicts j1-j2, j2-j3, j4-j5; j6 isolated.
#       ceiling 4 = {j1,j3,j4,j6}; cover 2 = {j1,j3,j5}+{j2,j4,j6}; max antichain 2.
#     touch_ok=False j1 and j3 now TOUCH at 20 and therefore conflict, so {j1,j2,j3} is a pairwise-
#       conflicting TRIPLE: ceiling 3 = {j1,j4,j6}, cover 3, antichain 3.
#     ⚠ The strict-mode cover was first hand-predicted as 2. That hand value was WRONG and brute force
#       caught it — which is exactly why the toy is here rather than a bare assertion.
TOY = [(10, 20), (15, 25), (20, 30), (30, 40), (35, 50), (60, 70)]
check("toy: ceiling is 4 under the half-open convention", R.single_path_ceiling(TOY, True) == 4,
      str(R.max_chain(TOY, True)))
check("toy: ceiling drops to 3 when a zero-length exon is disallowed",
      R.single_path_ceiling(TOY, False) == 3, str(R.max_chain(TOY, False)))
check("toy: min_path_cover is 2 / 3 (NOT the ceiling — it is the DUAL number)",
      (R.min_path_cover(TOY, True), R.min_path_cover(TOY, False)) == (2, 3),
      f"{R.min_path_cover(TOY, True)} / {R.min_path_cover(TOY, False)}")
check("toy: nesting conflicts", R.junction_conflict((10, 100), (40, 50)))
check("toy: shared donor conflicts", R.junction_conflict((10, 50), (10, 90)))
check("toy: shared acceptor conflicts", R.junction_conflict((10, 90), (50, 90)))
check("toy: touching is compatible iff touch_ok",
      not R.junction_conflict((10, 20), (20, 30), True)
      and R.junction_conflict((10, 20), (20, 30), False))

# (b) DILWORTH, live, on randomised inputs. min chain cover (matching) == max antichain (sweep) are
#     computed by completely unrelated code, so agreement checks both. Also brute-forced at small n.
import itertools  # noqa: E402
import random  # noqa: E402


def _brute_chain(J, tk):
    js = sorted(set(J))
    for k in range(len(js), 0, -1):
        for c in itertools.combinations(js, k):
            if all(not R.junction_conflict(a, b, tk) for a, b in itertools.combinations(c, 2)):
                return k
    return 0


def _brute_antichain(J, tk):
    js = sorted(set(J))
    for k in range(len(js), 0, -1):
        for c in itertools.combinations(js, k):
            if all(R.junction_conflict(a, b, tk) for a, b in itertools.combinations(c, 2)):
                return k
    return 0


random.seed(20260807)
bad_chain = bad_ac = bad_dil = 0
N_CASES = int(os.environ.get("RUSTLE_CEILING_CASES", "600"))
for _ in range(N_CASES):
    J = []
    for _i in range(random.randint(1, 8)):
        d = random.randint(0, 40)
        J.append((d, d + random.randint(1, 15)))
    for tk in (True, False):
        if R.single_path_ceiling(J, tk) != _brute_chain(J, tk):
            bad_chain += 1
        ac = R.max_antichain(J, tk)
        if ac != _brute_antichain(J, tk):
            bad_ac += 1
        if R.min_path_cover(J, tk) != ac:
            bad_dil += 1
check(f"max_chain == brute force ({N_CASES} random cases x 2 modes)", bad_chain == 0, f"{bad_chain} bad")
check("max_antichain == brute force", bad_ac == 0, f"{bad_ac} bad")
check("DILWORTH: n - max_matching == max_antichain (two unrelated implementations)",
      bad_dil == 0, f"{bad_dil} bad")

# (c) THE RECORDED NPIP NUMBERS, denominators included. These are the published pairs; if the relation or
#     the min_len filter drifts, they move. NPIPB2 is 9/14 only WITH introns_of_transcripts' min_len=20
#     default — at min_len=0 a 1-bp transcript-to-genome alignment gap at chr16:11,969,658 makes it 10/15.
CEIL_REC = {"NPIPB6": (12, 21), "NPIPA2": (12, 24), "NPIPB9": (12, 21),
            "NPIPB15": (8, 14), "NPIPA5": (10, 18), "NPIPB2": (9, 14)}
if not os.path.exists(R.GFF):
    skip("recorded NPIP ceilings", f"no GFF at {R.GFF}")
else:
    for g, (ceil, union) in CEIL_REC.items():
        tx = R.introns_of_transcripts("chr16", g)
        U = sorted({iv for t in tx.values() for iv in t})
        rec = R.locus_ceiling(U)
        check(f"{g} ceiling is {ceil}/{union}",
              (rec["ceiling"], rec["n_junctions"]) == (ceil, union),
              f"got {rec['ceiling']}/{rec['n_junctions']} over {len(tx)} transcripts")
    # NPIPB6 is the worked example: 16 transcripts, 21 junctions, ceiling 12, and FIVE reps needed.
    tx = R.introns_of_transcripts("chr16", "NPIPB6")
    U = sorted({iv for t in tx.values() for iv in t})
    rec = R.locus_ceiling(U)
    check("NPIPB6 needs 5 single-path reps (min_paths is NOT the ceiling)", rec["min_paths"] == 5,
          str(rec["min_paths"]))
    check("NPIPB6 has 16 annotated transcripts", len(tx) == 16, str(len(tx)))
    check("the min_len=20 filter is what makes NPIPB2 9/14 and not 10/15",
          R.locus_ceiling(sorted({iv for t in R.introns_of_transcripts("chr16", "NPIPB2", min_len=0).values()
                                  for iv in t}))["n_junctions"] == 15,
          "min_len=0 should keep the 1-bp alignment gap")

# (d) ⚠ THE BANNED UNIVERSE. `read_exons` returns MERGED exon blocks, so the only junctions derivable from
#     it are the gaps BETWEEN consecutive blocks — which are pairwise disjoint BY CONSTRUCTION. Ceiling
#     then equals |J| at every locus and the normalised score is 1.000 everywhere: the killed-metric
#     pattern (prediction ⊆ its own truth by construction). This test asserts the degeneracy so that
#     anyone who reintroduces the route fails loudly instead of reporting a perfect score.
if not os.path.exists(BAM):
    skip("banned read_exons-gap universe stays banned", f"no BAM at {BAM}")
else:
    blocks = R.read_exons(BAM, "chr16", 28623032, 28645151)          # NPIPB6
    gaps = [(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1)]
    rec = R.locus_ceiling(gaps) if gaps else {"ceiling": 0, "n_junctions": 0}
    check("read_exons gaps are degenerate (ceiling == |J|) — DO NOT use them as the universe",
          rec["ceiling"] == rec["n_junctions"],
          f"{rec['ceiling']}/{rec['n_junctions']} — if this ever differs, re-read why the route was banned")
    print(f"        NPIPB6 read_exons route: {len(blocks)} merged blocks -> "
          f"{rec['n_junctions']} gap-junctions, ceiling {rec['ceiling']} (ratio would be 1.000)")

# ---------------------------------------------------------------- the persisted rep-quality record
print("\nrep_quality record (the ceiling wired to a real catalog)")
import rep_quality as Q  # noqa: E402

# The exons -> junctions conversion. copies.tsv `exons` is 0-based half-open genomic BLOCKS; the junctions
# are the GAPS between consecutive blocks, so n_junctions == n_exon - 1 always.
row = {"exons": "100-200,300-400,500-600"}
check("rep_junctions turns exon blocks into the gaps between them",
      Q.rep_junctions(row) == [(200, 300), (400, 500)], str(Q.rep_junctions(row)))
check("a single-exon (DNA-substrate) rep has ZERO junctions, not NA",
      Q.rep_junctions({"exons": "100-200"}) == [])
check("a rep with no exons column is NA, not zero", Q.rep_junctions({}) is None)

# END-TO-END on the worked example. NPIPB6's representative carries 7 junctions that are in the annotated
# union of 21, whose single-path ceiling is 12 -> raw 0.333, normalised 0.583. These are the RECORDED
# numbers; if either the ceiling or the rep-join drifts, this fails.
QCOP = os.environ.get("RUSTLE_REPQ_COPIES", "/home/juanfra/winloci_scratch/hp/NPIPREG.copies.tsv")
if not (os.path.exists(QCOP) and os.path.exists(R.GFF)):
    skip("NPIPB6 end-to-end 7/21 and 7/12", f"no catalog at {QCOP}")
else:
    ns = argparse.Namespace(copies=QCOP, bam=None, chrom="chr16", genes="NPIPB6", gene_prefix=None,
                            universe="annotation", min_reads=3, min_intron_len=20,
                            touch_ok=True, contained=True)
    r = Q.build_rows(ns)[0]
    check("NPIPB6 record is (on 7, union 21, ceiling 12)",
          (r["rep_on_universe"], r["n_junctions_union"], r["ceiling"]) == (7, 21, 12), str(r))
    check("NPIPB6 raw 0.3333 / normalised 0.5833 (the recorded pair)",
          (r["raw"], r["normalised"]) == ("0.3333", "0.5833"), f"{r['raw']} / {r['normalised']}")
    check("NPIPB6 needs 5 single-path reps", r["min_paths"] == 5, str(r["min_paths"]))
    check("the record joins to copies.tsv by node_key",
          r["node_key"] == R.canon_node(f"chr16:{r['rep_start'] + 1}-{r['rep_end']}"), str(r["node_key"]))

# ⚠ 0/0 MUST NOT BE WRITTEN AS 0.0. A locus with no junction clearing the floor has ceiling 0; a silent
# zero there reads as "the model explained nothing" when the truth is "there was nothing to explain".
empty = R.locus_ceiling([])
check("an empty junction set gives ceiling 0 and min_paths 0, not a crash",
      (empty["ceiling"], empty["min_paths"], empty["n_junctions"]) == (0, 0, 0), str(empty))

print(f"\n{len(FAILED)} failed, {len(SKIPPED)} skipped")
sys.exit(1 if FAILED else 0)
