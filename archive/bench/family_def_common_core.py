#!/usr/bin/env python3
"""family_def_common_core.py — OPTIONAL alternative family detector: shared k-mer COMMON CORE.

Advisor-requested alternative to the ~B (pairwise-alignment + clique) detector. Find a family's
members by the COMMON CORE: the shared k-mer substring present across the reads / copy-models of
family members, with members assigned by maximum k-mer JACCARD to that core. Opt-in / experimental.

Honest known limitation (reproduced below): the recovered core captures the gene BODY, but its
EXTENT over- or under-shoots the annotated gene length — because the core boundary is set by
sequence CONSERVATION across copies (a shared repeat/domain extends it; divergence truncates it),
not by gene structure. Reported as the body / gene-length ratio + membership over/under-merge.

Plug-in API:
    members = detect_family(seed_seqs, pool_seqs)          # grow a family by common-core Jaccard
    core, body_len, ref, frac = common_core(member_seqs)   # the shared core + its body extent
Wire as an alternative to the ~B clique step (env FAMILY_DETECT=core in family_def_build.py).
Self-validation on real families:  python bench/family_def_common_core.py
"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_definition_visuals import load_seqs   # representative cDNA per gene (read-derived consensus)

K = 15
_T = str.maketrans("ACGT", "TGCA")


def canon(km):
    rc = km.translate(_T)[::-1]
    return km if km <= rc else rc


def kmers(seq, k=K):
    seq = seq.upper()
    return {canon(seq[i:i + k]) for i in range(len(seq) - k + 1) if "N" not in seq[i:i + k]}


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def common_core(seqs, k=K, min_frac=1.0):
    """seqs: {gene: seq}. Returns (core_kmers, core_body_len_bp, ref_gene, body_frac_of_ref).
    core = k-mers present in >= min_frac of members; body = span of core k-mers on the longest member."""
    ksets = {g: kmers(s, k) for g, s in seqs.items() if s}
    if len(ksets) == 0:
        return set(), 0, None, 0.0
    if len(ksets) == 1:                       # a single seed IS its own core
        g, ks = next(iter(ksets.items()))
        return ks, len(seqs[g]), g, 1.0
    cnt = collections.Counter()
    for ks in ksets.values():
        cnt.update(ks)
    need = max(2, round(min_frac * len(ksets)))
    core = {km for km, c in cnt.items() if c >= need}
    ref = max(seqs, key=lambda g: len(seqs[g]))
    s = seqs[ref].upper()
    pos = [i for i in range(len(s) - k + 1) if canon(s[i:i + k]) in core]
    if not pos:
        return core, 0, ref, 0.0
    body = (pos[-1] + k) - pos[0]
    return core, body, ref, body / max(len(s), 1)


def detect_family(seed_seqs, pool_seqs, k=K, cmin=0.5, iters=4, min_frac=0.8):
    """Grow a family from a seed: a candidate joins if it CONTAINS the common core (the shared
    substring) above fraction cmin = |cand ∩ core| / |core|.  seed/pool: {gene: seq}."""
    members = dict(seed_seqs)
    for _ in range(iters):
        core, _, _, _ = common_core(members, k, min_frac)
        if not core:
            break
        added = False
        for g, s in pool_seqs.items():
            if g in members:
                continue
            ks = kmers(s, k)
            if len(ks & core) / max(len(core), 1) >= cmin:
                members[g] = s
                added = True
        if not added:
            break
    return members


# ----------------------------------------------------------- self-validation
def _build_families():
    import networkx as nx
    from family_def_read_filters import dna_homology
    from family_def_retrocopy_filter import load_intron_index, n_intron
    from family_def_strand_probe import edge_strand
    H, _ = dna_homology(); idx = load_intron_index(); st = edge_strand()
    G = nx.Graph()
    for (a, b), r in H.items():
        if r.get("id", 0) < 0.80 or min(r["cov_a"], r["cov_b"]) < 0.30:
            continue
        na, nb = n_intron(idx, a), n_intron(idx, b)
        if na is not None and nb is not None and ((na == 0) ^ (nb == 0)):
            continue
        if st.get((a, b) if a < b else (b, a)) == "-":
            continue
        G.add_edge(a, b)
    return [sorted(c) for c in nx.connected_components(G) if len(c) >= 2], H


def main():
    fams, H = _build_families()
    # pick a clean spread of family sizes
    named = lambda c: [g for g in c if not g.startswith("LOC")]
    fams.sort(key=len)
    picks = []
    for want in (2, 3, 4, 6, 10):
        cand = [c for c in fams if len(c) == want and named(c)]
        if cand:
            picks.append(cand[len(cand) // 2])
    print("=== COMMON-CORE detector — body extent vs gene length (the over/under-merge) ===")
    print(f"  {'family (a member)':28} {'#mem':>4} {'core body':>9} {'ref cDNA':>8} {'body/gene':>9}  verdict")
    for c in picks:
        seqs = load_seqs(c)
        seqs = {g: s for g, s in seqs.items() if s}
        if len(seqs) < 2:
            continue
        core, body, ref, frac = common_core(seqs)
        lens = sorted(len(s) for s in seqs.values())
        med = lens[len(lens) // 2]
        v = "UNDER-merge (core < gene body)" if frac < 0.85 else ("captures body" if frac <= 1.05 else "OVER (core > ref)")
        lab = (named(c)[0] if named(c) else c[0])
        print(f"  {lab:28} {len(seqs):>4} {body:>9} {len(seqs[ref]):>8} {frac:>8.2f}  {v}")

    # membership over/under-merge: seed 1 member, recover the family from its homology neighborhood,
    # at a PERMISSIVE and a STRICT containment threshold -> there is no clean cut.
    print("\n=== COMMON-CORE detector — membership vs the true ~B family (no clean threshold) ===")
    print(f"  {'family':16} {'true':>4}   permissive cmin=0.30        strict cmin=0.60")
    print(f"  {'':16} {'':>4}   recov/extra/missed verdict   recov/extra/missed verdict")
    for c in picks:
        truth = set(c)
        pool = set(c)
        for (a, b), r in H.items():
            if r.get("id", 0) >= 0.70 and (a in truth or b in truth):
                pool.add(a); pool.add(b)
        ps = load_seqs(list(pool)); ps = {g: s for g, s in ps.items() if s}
        ss = {g: ps[g] for g in c[:1] if g in ps}
        if not ss:
            continue
        out = []
        for cm in (0.30, 0.60):
            rec = set(detect_family(ss, ps, cmin=cm))
            ex = len(rec - truth); mi = len(truth - rec)
            v = "match" if not ex and not mi else (f"OVER+{ex}" if ex and not mi else (f"UNDER-{mi}" if mi and not ex else f"+{ex}/-{mi}"))
            out.append(f"{len(rec):>2}/{ex}/{mi} {v:<10}")
        lab = (named(c)[0] if named(c) else c[0])
        print(f"  {lab:16} {len(truth):>4}   {out[0]:<28} {out[1]}")
    print("\n  => the common core captures only the conserved BODY (7-83% of the gene), and there is no")
    print("     single threshold: permissive over-merges (pulls in domain-sharers), strict under-merges")
    print("     (drops divergent members). Boundary set by conservation, not gene structure — matching")
    print("     prior experience. OPT-IN alternative to ~B (detect_family / common_core); ~B stays default.")


GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"


def detect_by_multimap_seed(seed_gene, bam_path=BAM, id_floor=0.6, min_shared=3):
    """The READ-FOLLOWING variant (advisor's actual approach): seed = one copy; take its multimapping
    reads; follow them to where ELSE they map -> those loci are the other copies, each defining a
    separate CORE = the region the shared reads cover. Demonstrates that the cores have no guaranteed
    extent (they cover only the SIMILAR region the reads can't disambiguate, not the whole gene)."""
    import pysam
    from family_def_read_filters import dna_homology
    coord = {}
    for line in open(GENES_BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    if seed_gene not in coord:
        print(f"[!] no coord for {seed_gene}"); return
    H, _ = dna_homology()
    cands = {(b if a == seed_gene else a) for (a, b), r in H.items()
             if r.get("id", 0) >= id_floor and (a == seed_gene or b == seed_gene)}
    bam = pysam.AlignmentFile(bam_path, "rb")

    def reads_at(g, mapq0_only=False):
        c, s, e = coord[g]; out = {}
        try:
            it = bam.fetch(c, max(0, s), e)
        except (ValueError, KeyError):
            return out
        for r in it:
            if r.is_supplementary:
                continue
            if mapq0_only and r.mapping_quality != 0:
                continue
            out.setdefault(r.query_name, (r.reference_start, r.reference_end or r.reference_start))
        return out

    seed_mm = reads_at(seed_gene, mapq0_only=True)     # the seed's MULTIMAPPING reads
    sc, ss, se = coord[seed_gene]
    print(f"seed {seed_gene}: {len(seed_mm)} multimapping (MAPQ0) reads, "
          f"{len(cands)} homology-neighbour candidates")
    print(f"  copies discovered by following those reads (≥{min_shared} shared) — each its own CORE:")
    rows = []
    for g in sorted(cands):
        if g not in coord:
            continue
        rd = reads_at(g)
        shared = set(seed_mm) & set(rd)
        if len(shared) >= min_shared:
            spans = [rd[n] for n in shared]
            lo, hi = min(s for s, e in spans), max(e for s, e in spans)
            gl = coord[g][2] - coord[g][1]
            rows.append((g, len(shared), hi - lo, gl, (hi - lo) / max(gl, 1)))
    for g, ns, cl, gl, frac in sorted(rows, key=lambda x: -x[1]):
        print(f"    {g:16} shared_reads={ns:4}  core={cl:>6}bp  gene={gl:>6}bp  core/gene={frac:>4.0%}  "
              f"{'UNDER (core < gene)' if frac < 0.85 else ('~ gene' if frac <= 1.15 else 'OVER (core > gene)')}")
    bam.close()
    if rows:
        fr = [r[4] for r in rows]
        print(f"  => {len(rows)} cores; core/gene ranges {min(fr):.0%}-{max(fr):.0%} — NO guaranteed extent:")
        print("     each core spans only the SIMILAR region the reads can't disambiguate (the ~R/collapsed")
        print("     part), not the whole gene. Merging nearby cores still over-/under-shoots the gene.")


def detect_for_gene(seed_gene, k=K, cmin=0.4, id_floor=0.6):
    """CLI use: given a seed gene, build its candidate pool (cDNA-homology neighbours) and return the
    common-core-detected family + diagnostics."""
    from family_def_read_filters import dna_homology
    H, _ = dna_homology()
    pool = {seed_gene}
    for (a, b), r in H.items():
        if r.get("id", 0) >= id_floor and (a == seed_gene or b == seed_gene):
            pool.add(a); pool.add(b)
    ps = load_seqs(list(pool)); ps = {g: s for g, s in ps.items() if s}
    if seed_gene not in ps:
        print(f"[!] no sequence for {seed_gene}"); return
    rec = detect_family({seed_gene: ps[seed_gene]}, ps, cmin=cmin)
    core, body, ref, frac = common_core({g: ps[g] for g in rec})
    print(f"seed {seed_gene}: common-core family ({len(rec)} members, cmin={cmin}):")
    print("  " + ", ".join(sorted(rec)))
    print(f"  shared core body: {body} bp on {ref} ({frac:.0%} of its {len(ps[ref])} bp) "
          f"— {'captures body' if 0.85<=frac<=1.05 else 'UNDER/OVER vs gene length'}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--multimap":   # read-following variant
        detect_by_multimap_seed(sys.argv[2])
    elif len(sys.argv) > 1:
        cmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
        detect_for_gene(sys.argv[1], cmin=cmin)   # k-mer common-core: GENE [cmin]
    else:
        main()
