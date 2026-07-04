#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the shared exon-graph
LIBRARY of `bench/recombination_bridge_detector.py` (RBD) -> `vg_family/bridge_detector.rs`
(the O1 over-merge-gate migration, RETIREMENT_AND_MIGRATION.md Tier-2 #6, full-parity path).

It imports the SHIPPED RBD and runs its REAL library functions (aln_id / revcomp /
load_skeletons / family_exons / exon_match_tensor / build_graph / colinear_cov /
pair_coverage / cross_component_maxid / subfams_colocated / subfam_dom_gene /
subfam_label) to emit a SELF-CONTAINED fixture to
`src/rustle/vg_family/testdata/bridge_detector_fixture.json` with four parts:

  (a) "aln_id": >= 80 (q,t) pairs -> the REAL edlib id.  PROVES the pure-Rust HW-DP
      == edlib.  Classes: exact-infix, q-longer-than-t, indels(1-3), substitutions,
      revcomp-only, both/one-empty, N-containing, periodic/repeats (multiple
      equal-distance infixes), q==t, case differences.
  (b) "load_skeletons": real sampled skeleton rows + synthetic edge cases (unsorted,
      leading/trailing, whole-span, zero-length, adjacent) -> the REAL RBD exon lists
      (RBD's 1-based-inclusive complement, NO intron sort, b>=a filter).
      "family_exons": >= 20 families -> per member (skel exon coords + RAW case-carrying
      genome slices) and the REAL family_exons recs (spliced, strand-ORDER-reversed for
      '-', uppercased).  Ships the genome slices so no 3.4G genome is needed at test time.
      Includes '-' strand multi-exon members (order reversal) + a dropped (no-skeleton)
      member.
  (c) "tensor_families": >= 15 families -> the REAL exon_match_tensor / build_graph
      adjacency+edge_exons / colinear_cov (all ordered pairs) / pair_coverage / connected
      components / cross_component_maxid / subfams_colocated / subfam_dom_gene /
      subfam_label.  Cost-bounded (small exons) so the Rust HW-DP test runs fast.
  (d) "subfam_cases": hand-built coord+partition cases that FORCE subfams_colocated
      True/False + subfam_dom_gene/label tie-breaking (first-seen), via REAL RBD.

DEFERRED (NOT ported -- the gates call the LIBRARY above, not these): analyze_family,
main, _pure_controls / _control_precompute / control_recall_cost, load_families /
target_families / load_strand (loaders wired by the gate driver, not the library math),
CARD_ID / N_CONTROL (used only by the deferred control-recall).

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_bridge_detector_fixture.py
Deterministic (sorted selection, PYTHONHASHSEED=0).
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")   # repo root (RBD uses relative paths)

import json
import struct
import tempfile

sys.path.insert(0, "bench")
import recombination_bridge_detector as RBD   # SHIPPED library (ground truth)
import pysam
import networkx as nx

OUT = os.path.join("src", "rustle", "vg_family", "testdata", "bridge_detector_fixture.json")


# ------------------------------------------------------------------ helpers
def rc(s):
    return RBD.revcomp(s)


def fb(x):
    """Raw IEEE-754 u64 bit pattern of a float -- shipped alongside every bit-compared
    value so the Rust test reads EXACT bits (serde_json's decimal float parser is 1-ULP
    imprecise; ints round-trip exactly).  Tolerance NONE."""
    return struct.unpack("<Q", struct.pack("<d", float(x)))[0]


def comps_of(G):
    """connected components as sorted lists, sorted by first element (deterministic)."""
    return sorted((sorted(c) for c in nx.connected_components(G)), key=lambda c: c[0])


# ================================================================== part (a) aln_id
def build_aln_id_cases():
    cases = []

    def add(cls, q, t):
        idv = RBD.aln_id(q, t)
        cases.append(dict(cls=cls, q=q, t=t, id=idv, bits=fb(idv)))

    # exact infix (q is a substring of t)
    add("exact_infix", "ACGT", "TTACGTGG")
    add("exact_infix", "GATTACA", "CCCGATTACAAAA")
    add("exact_infix", "ACGTACGT", "ACGTACGT")            # q == t
    add("exact_infix_prefix", "ACGT", "ACGTGGGG")
    add("exact_infix_suffix", "ACGT", "GGGGACGT")
    add("q_equals_t", "ACGTACGTAC", "ACGTACGTAC")
    add("q_equals_t", "TTTTGGGGCCCCAAAA", "TTTTGGGGCCCCAAAA")

    # q longer than t (must still equal edlib HW)
    add("q_longer", "ACGTACGTACGT", "ACGT")
    add("q_longer", "ACGTACGTACGT", "ACG")
    add("q_longer", "GATTACAGATTACA", "GATTACA")
    add("q_longer", "AAAAAAAAAAAA", "AAA")
    add("q_longer_mismatch", "ACGTACGTTTTT", "ACGT")

    # single substitutions
    add("sub1", "ACGTACGT", "ACGAACGT")
    add("sub1", "GATTACA", "CCCGATCACAAA")
    add("sub2", "ACGTACGTAC", "TTACCTACGACGG")
    add("sub3", "ACGTACGTACGT", "ACGAACTTACGA")

    # indels 1-3
    add("ins1", "ACGTACGT", "TTACGTAACGTGG")
    add("del1", "ACGTACGT", "TTACGACGTGG")
    add("indel2", "ACGTACGTAC", "TTACGTGTACGG")
    add("indel3", "ACGTACGTACGT", "ACGTCGTACT")

    # revcomp-only matches: forward is poor, revcomp(q) matches t well
    for base in ["ACGTACGTAC", "GATTACAGGA", "TTGGCCAATT", "ACACACGTGT", "CCGGTTAACC"]:
        add("revcomp_only", rc(base), "TTT" + base + "GGG")
    add("revcomp_only_exact", rc("GATTACA"), "GATTACA")
    add("revcomp_noise", rc("ACGTACGTAC"), "AA" + "ACGTACGTAC" + "T")

    # empties
    add("both_empty", "", "")
    add("q_empty", "", "ACGT")
    add("t_empty", "ACGT", "")
    add("t_empty2", "GATTACA", "")

    # N-containing (edlib treats N as a normal mismatchable symbol)
    add("n_in_q", "ACNTACGT", "TTACGTACGTGG")
    add("n_in_t", "ACGTACGT", "TTACNTACGTGG")
    add("n_both", "ACNTACNT", "GGACNTACNTAA")
    add("n_run", "NNNN", "ACNNNNGT")
    add("all_n", "NNNNN", "NNNNNNN")

    # periodic / tandem repeats (multiple equal-distance infixes -> tie handling)
    add("periodic", "ACAC", "ACACACACACAC")
    add("periodic", "ATAT", "TATATATATATA")
    add("tandem", "GATTACA", "GATTACAGATTACAGATTACA")
    add("homopolymer", "AAAA", "AAAAAAAAAAAA")
    add("homopolymer_sub", "AAAA", "AAGAAAAAAA")
    add("periodic_offset", "CAGCAG", "GCAGCAGCAGCAGC")

    # case differences (revcomp table maps lowercase too; aln_id is case-sensitive on match)
    add("lowercase_q", "acgt", "TTACGTGG")
    add("lowercase_t", "ACGT", "ttacgtgg")
    add("mixed_case", "AcGtAcGt", "acgtACGTacgt")
    add("lower_exact", "gattaca", "cccgattacaaaa")
    add("lower_revcomp", rc("acgtacgt"), "acgtacgt")

    # longer realistic-ish random-but-fixed pairs (varied identity)
    seqs = [
        ("AAGGCCTTAGCTAGCTAGGCTA", "TTAAGGCCTTAGCTAGCTAGGCTAGG"),
        ("ATGCATGCATGCATGCATGC", "ATGCATGGATGCATGCTTGC"),
        ("GGGGATCCGGGGATCCGGGG", "AAGGGGATCCGGGGATCCGGGGTT"),
        ("TACGTTACGGATCCAATTGG", "TTACGTTACGCATCCAATTGGA"),
        ("CATGGCATGGCATGGCATGG", "GGCATGGCATGACATGGCATGG"),
        ("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT"),   # no similarity
        ("AAACCCGGGTTTAAACCCGG", "AAACCCGGGTTTAAACCCGG"),    # q==t medium
    ]
    for q, t in seqs:
        add("realish", q, t)
        add("realish_rc", rc(q), t)      # revcomp orientation

    # a few more mismatch-heavy / worst-case
    add("orthogonal", "AAAAAAAAAA", "CCCCCCCCCC")
    add("orthogonal", "ACGTACGTAC", "TGCATGCATG")
    add("single_base_q", "A", "TTTTAGGGG")
    add("single_base_q_miss", "A", "TTTTTTTTT")
    add("two_base", "AC", "GGACGG")

    # deterministic pseudo-random pairs across the identity spectrum (fixed seed).
    import random
    rng = random.Random(0)
    bases = "ACGT"
    for _ in range(30):
        lt = rng.randint(20, 60)
        t = "".join(rng.choice(bases) for _ in range(lt))
        # carve q from a random infix of t, then apply a few random edits
        li = rng.randint(8, max(8, lt))
        st = rng.randint(0, max(0, lt - min(li, lt)))
        q = list(t[st:st + li]) if lt else []
        nedit = rng.randint(0, 4)
        for _ in range(nedit):
            if not q:
                break
            op = rng.randint(0, 2)
            p = rng.randint(0, len(q) - 1)
            if op == 0:                       # substitution
                q[p] = rng.choice(bases)
            elif op == 1:                     # deletion
                del q[p]
            else:                             # insertion
                q.insert(p, rng.choice(bases))
        qs = "".join(q) if q else "A"
        # half the time test the revcomp orientation
        if rng.random() < 0.5:
            qs = rc(qs)
        add("rand", qs, t)

    return cases


# ================================================================== part (b) load_skeletons
def build_load_skeletons_cases(skel):
    cases = []
    seen = 0
    # (b.1) REAL rows: read the raw skeleton file, look up RBD's parsed exon list.
    with open(RBD.SKEL_TSV) as fh:
        next(fh)
        picked_multi = 0
        picked_single = 0
        for ln in fh:
            t = ln.rstrip("\n").split("\t")
            if len(t) < 6:
                t = (t + [""] * 6)[:6]
            chrom, start, end, introns = t[0], int(t[1]), int(t[2]), t[5]
            key = (chrom, start, end)
            exons = skel.get(key)
            if exons is None:
                continue
            n_iv = introns.count("-") if introns else 0
            if n_iv >= 2 and picked_multi < 18:
                picked_multi += 1
            elif n_iv == 0 and picked_single < 6:
                picked_single += 1
            elif n_iv == 1 and picked_multi < 24:
                picked_multi += 1
            else:
                continue
            cases.append(dict(chrom=chrom, start=start, end=end, introns=introns,
                              exons=[list(e) for e in exons], kind="real"))
            seen += 1
            if picked_multi >= 24 and picked_single >= 6:
                break

    # (b.2) SYNTHETIC rows via a temp skeleton file + REAL RBD.load_skeletons monkeypatch.
    synth_lines = [
        ("chrA", 100, 300, "150-200"),                 # one interior intron
        ("chrB", 100, 300, ""),                         # blank -> single exon [start,end]
        ("chrC", 100, 300, "100-150"),                  # leading intron a==start
        ("chrD", 100, 300, "200-300"),                  # trailing intron b==end
        ("chrE", 100, 300, "100-300"),                  # whole span intron -> []
        ("chrF", 100, 300, "250-280;120-150"),          # UNSORTED (RBD does NOT sort!)
        ("chrG", 100, 300, "120-150;150-200"),          # adjacent introns
        ("chrH", 100, 300, "150-150"),                  # zero-length intron (a==b)
        ("chrI", 100, 400, "120-140;200-230;330-360"),  # three sorted introns
        ("chrJ", 100, 400, "330-360;120-140;200-230"),  # three UNSORTED introns
    ]
    tmpf = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    tmpf.write("chrom\tstart\tend\tn_exon\tn_reads\tintrons\n")
    for chrom, start, end, introns in synth_lines:
        tmpf.write(f"{chrom}\t{start}\t{end}\t0\t0\t{introns}\n")
    tmpf.close()
    old = RBD.SKEL_TSV
    RBD.SKEL_TSV = tmpf.name
    syn_skel = RBD.load_skeletons()
    RBD.SKEL_TSV = old
    os.unlink(tmpf.name)
    for chrom, start, end, introns in synth_lines:
        exons = syn_skel[(chrom, start, end)]
        cases.append(dict(chrom=chrom, start=start, end=end, introns=introns,
                          exons=[list(e) for e in exons], kind="synthetic"))
    return cases


# ================================================================== part (b) family_exons
def capture_family_exons(fid, members, skel, strand, fa):
    """Run the REAL family_exons and capture the RAW case-carrying genome slices it fetched
    (keyed (chrom, s-1, e)) so the Rust test can reproduce recs with a mock fetcher."""
    recs = RBD.family_exons(members, skel, strand, fa)
    skel_out = []
    slices = {}
    for m in members:
        key = (m["chrom"], m["start"], m["end"])
        exons = skel.get(key)
        if exons is None:
            continue
        skel_out.append(dict(chrom=m["chrom"], start=m["start"], end=m["end"],
                             exons=[list(e) for e in exons]))
        for (s, e) in exons:
            raw = fa.fetch(m["chrom"], s - 1, e)      # case-carrying, BEFORE .upper()
            slices[f"{m['chrom']}:{s-1}:{e}"] = raw
    recs_out = [dict(dn=r["dn"], gene=r["gene"], chrom=r["chrom"], start=r["start"],
                     end=r["end"], strand=r["strand"], exons=r["exons"], exlen=r["exlen"])
                for r in recs]
    return dict(fid=fid,
                members=[dict(dn=m["dn"], gene=m["gene"], chrom=m["chrom"],
                              start=m["start"], end=m["end"]) for m in members],
                strand={m["dn"]: strand.get(m["dn"], "+") for m in members},
                skel=skel_out,
                slices=[dict(chrom=k.split(":")[0], start0=int(k.split(":")[1]),
                             end0=int(k.split(":")[2]), bytes=v) for k, v in slices.items()],
                recs=recs_out)


# ================================================================== part (c) tensor
def capture_tensor_family(fid, recs):
    best = RBD.exon_match_tensor(recs)
    G, edge_exons = RBD.build_graph(recs, best)
    n = len(recs)
    # best in INSERTION order (a outer, b middle, i inner); trailing element = id bits.
    best_list = [[a, i, b, idv, j, fb(idv)] for (a, i, b), (idv, j) in best.items()]
    edges = sorted({tuple(sorted((u, v))) for u, v in G.edges()})
    ee = [[a, b, list(v)] for (a, b), v in edge_exons.items()]
    colinear = [[b, nn, (cv := RBD.colinear_cov(recs, best, b, nn)), fb(cv)]
                for b in range(n) for nn in range(n) if b != nn]
    paircov = [[a, b, (pv := RBD.pair_coverage(recs, best, a, b)), fb(pv)]
               for a in range(n) for b in range(n) if a != b]
    comps = comps_of(G)
    xcid = RBD.cross_component_maxid(recs, best, comps)
    coloc = RBD.subfams_colocated(recs, comps)
    dom = [RBD.subfam_dom_gene(recs, c) for c in comps]
    lab = [RBD.subfam_label(recs, c) for c in comps]
    return dict(
        fid=fid,
        recs=[dict(dn=r["dn"], gene=r["gene"], chrom=r["chrom"], start=r["start"],
                   end=r["end"], strand=r["strand"], exons=r["exons"]) for r in recs],
        best=best_list, edges=[list(e) for e in edges], edge_exons=ee,
        colinear_cov=colinear, pair_coverage=paircov,
        components=[list(c) for c in comps], cross_component_maxid=xcid,
        cross_component_maxid_bits=fb(xcid),
        subfams_colocated=coloc, subfam_dom_gene=dom, subfam_label=lab,
    )


def dp_cost(recs):
    """estimate HW-DP cell-ops of exon_match_tensor for a family (for cost-bounded pick)."""
    n = len(recs)
    cost = 0
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            for qa in recs[a]["exons"]:
                if len(qa) < RBD.MIN_EXON:
                    continue
                for tb in recs[b]["exons"]:
                    if len(tb) < RBD.MIN_EXON:
                        continue
                    lo, hi = (len(qa), len(tb)) if len(qa) <= len(tb) else (len(tb), len(qa))
                    if lo / hi < 0.5:
                        continue
                    cost += 2 * len(qa) * len(tb)     # x2 for fwd+revcomp
    return cost


# ================================================================== part (d) subfam cases
def build_subfam_cases():
    """Hand-built recs (coords+gene only) + partitions -> REAL RBD outputs; forces
    subfams_colocated True/False + dom-gene FIRST-SEEN tie-breaking."""
    cases = []

    def mk(recs_coords, partition):
        recs = [dict(chrom=c, start=s, end=e, gene=g) for (c, s, e, g) in recs_coords]
        return dict(
            recs_coords=[list(x) for x in recs_coords],
            partition=[list(p) for p in partition],
            subfams_colocated=RBD.subfams_colocated(recs, partition),
            subfam_dom_gene=[RBD.subfam_dom_gene(recs, p) for p in partition],
            subfam_label=[RBD.subfam_label(recs, p) for p in partition],
        )

    # same chrom, OVERLAPPING intervals across parts -> colocated True
    cases.append(mk([("c1", 100, 200, "GA"), ("c1", 150, 250, "GB")], [[0], [1]]))
    # same chrom, NON-overlapping -> False
    cases.append(mk([("c1", 100, 200, "GA"), ("c1", 300, 400, "GB")], [[0], [1]]))
    # different chrom -> False
    cases.append(mk([("c1", 100, 200, "GA"), ("c2", 150, 250, "GB")], [[0], [1]]))
    # touching boundary (start==end) counts as overlap (<= / >=) -> True
    cases.append(mk([("c1", 100, 200, "GA"), ("c1", 200, 300, "GB")], [[0], [1]]))
    # dom-gene tie: comp iterated in given order; first-seen among tied-max wins
    cases.append(mk([("c1", 1, 2, "ZED"), ("c1", 3, 4, "ALPHA"),
                     ("c1", 5, 6, "ZED"), ("c1", 7, 8, "ALPHA")], [[0, 1, 2, 3]]))
    cases.append(mk([("c1", 1, 2, "ALPHA"), ("c1", 3, 4, "ZED"),
                     ("c1", 5, 6, "ALPHA"), ("c1", 7, 8, "ZED")], [[0, 1, 2, 3]]))
    # clear majority
    cases.append(mk([("c1", 1, 2, "MAJ"), ("c1", 3, 4, "MAJ"), ("c1", 5, 6, "MIN")], [[0, 1, 2]]))
    # three-part partition, one overlapping pair among parts 0 & 2
    cases.append(mk([("cX", 10, 20, "A"), ("cY", 10, 20, "B"), ("cX", 15, 25, "C")],
                    [[0], [1], [2]]))
    return cases


# ================================================================== main
def main():
    print("[load] skeletons / strand / families ...", flush=True)
    skel = RBD.load_skeletons()
    strand = RBD.load_strand()
    fams = RBD.load_families()
    targets = RBD.target_families()
    fa = pysam.FastaFile(RBD.GENOME)

    # ---- part (a)
    aln_cases = build_aln_id_cases()
    print(f"       aln_id cases: {len(aln_cases)}")

    # ---- part (b) load_skeletons
    ls_cases = build_load_skeletons_cases(skel)
    print(f"       load_skeletons cases: {len(ls_cases)} "
          f"({sum(1 for c in ls_cases if c['kind']=='real')} real + "
          f"{sum(1 for c in ls_cases if c['kind']=='synthetic')} synthetic)")

    # ---- part (b) family_exons: pick >=20 families, ensure '-' multi-exon + a dropped member.
    fe_out = []
    have_minus_multiexon = False
    have_dropped = False
    picked_fe = []
    for fid in sorted(fams):
        members = [dict(dn=m["dn"], gene=m["gene"], chrom=m["chrom"],
                        start=m["start"], end=m["end"]) for m in fams[fid]]
        recs = RBD.family_exons(members, skel, strand, fa)
        if len(recs) < 2:
            continue
        # bound the per-member exon length so the fixture stays small
        totlen = sum(sum(len(s) for s in r["exons"]) for r in recs)
        if totlen > 9000 or len(members) > 6:
            continue
        # does it have a '-' strand multi-exon member?  a dropped member?
        minus_multi = any(r["strand"] == "-" and len(r["exons"]) >= 2 for r in recs)
        dropped = any((m["chrom"], m["start"], m["end"]) not in skel for m in members)
        picked_fe.append((fid, members, minus_multi, dropped))

    # greedily assemble ~24 covering the two required properties
    chosen = []
    for fid, members, mm, dr in picked_fe:
        want = False
        if len(chosen) < 24:
            want = True
        if mm and not have_minus_multiexon:
            want = True
        if dr and not have_dropped:
            want = True
        if want:
            chosen.append(fid)
            have_minus_multiexon = have_minus_multiexon or mm
            have_dropped = have_dropped or dr
        if len(chosen) >= 24 and have_minus_multiexon and have_dropped:
            break
    # add a guaranteed synthetic dropped-member family if none found
    for fid in chosen:
        members = [dict(dn=m["dn"], gene=m["gene"], chrom=m["chrom"],
                        start=m["start"], end=m["end"]) for m in fams[fid]]
        fe_out.append(capture_family_exons(fid, members, skel, strand, fa))
    if not have_dropped:
        # append a synthetic family whose 2nd member has no skeleton (guaranteed drop)
        base = fams[chosen[0]][0]
        members = [dict(dn=base["dn"], gene=base["gene"], chrom=base["chrom"],
                        start=base["start"], end=base["end"]),
                   dict(dn="DN_SYN_NOSKEL", gene="SYN", chrom="chrNOPE",
                        start=1, end=999)]
        fe_out.append(capture_family_exons(-1, members, skel, strand, fa))
        have_dropped = True
    print(f"       family_exons families: {len(fe_out)} "
          f"(minus_multiexon={have_minus_multiexon}, dropped={have_dropped})")

    # ---- part (c) tensor: cost-bounded pick of >=15 families with edges; prefer coverage.
    BUDGET = 120_000_000     # total HW-DP cell-ops across the fixture
    cand = []
    for fid in sorted(fams):
        members = [dict(dn=m["dn"], gene=m["gene"], chrom=m["chrom"],
                        start=m["start"], end=m["end"]) for m in fams[fid]]
        if len(members) > 8:
            continue
        recs = RBD.family_exons(members, skel, strand, fa)
        if len(recs) < 2:
            continue
        best = RBD.exon_match_tensor(recs)
        G, _ = RBD.build_graph(recs, best)
        ncomp = nx.number_connected_components(G)
        arts = len(set(nx.articulation_points(G))) if G.number_of_edges() else 0
        c = dp_cost(recs)
        cand.append((c, fid, recs, G.number_of_edges(), ncomp, arts))
    cand.sort(key=lambda x: x[0])

    MAX_FAMS = 20            # keep the fixture small + the Rust HW-DP test snappy
    tensor_out = []
    spent = 0
    have_multicomp = 0
    have_art = 0
    chosen_t = set()
    # first pass: cheapest families with edges (hard cap on count + cost)
    for c, fid, recs, edges, ncomp, arts in cand:
        if fid in chosen_t:
            continue
        if edges == 0:
            continue
        if len(tensor_out) >= MAX_FAMS:
            break
        if spent + c > BUDGET and len(tensor_out) >= 15:
            break
        tensor_out.append(capture_tensor_family(fid, recs))
        chosen_t.add(fid)
        spent += c
        have_multicomp += (ncomp >= 2)
        have_art += (arts >= 1)
    # ensure multi-component coverage (over-merged: >=2 comps) even if a bit costlier
    if have_multicomp < 2:
        for c, fid, recs, edges, ncomp, arts in cand:
            if fid in chosen_t or ncomp < 2:
                continue
            tensor_out.append(capture_tensor_family(fid, recs))
            chosen_t.add(fid)
            spent += c
            have_multicomp += 1
            if have_multicomp >= 2:
                break
    if have_art < 1:
        for c, fid, recs, edges, ncomp, arts in cand:
            if fid in chosen_t or arts < 1:
                continue
            tensor_out.append(capture_tensor_family(fid, recs))
            chosen_t.add(fid)
            spent += c
            have_art += 1
            break
    print(f"       tensor families: {len(tensor_out)} (est HW-DP cells={spent:,}; "
          f"multicomp={have_multicomp}, articulation={have_art})")

    # ---- part (d)
    subfam_cases = build_subfam_cases()
    print(f"       subfam_cases: {len(subfam_cases)}")

    fa.close()

    fixture = dict(
        ID_THRESH=RBD.ID_THRESH, MIN_EXON=RBD.MIN_EXON, COLINEAR_COV=RBD.COLINEAR_COV,
        aln_id=aln_cases,
        load_skeletons=ls_cases,
        family_exons=fe_out,
        tensor_families=tensor_out,
        subfam_cases=subfam_cases,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)
    print(f"[fixture] wrote {os.path.normpath(OUT)} ({os.path.getsize(OUT):,} bytes)")


if __name__ == "__main__":
    main()
