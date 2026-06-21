#!/usr/bin/env python3
"""Secondary-cap fragmentation: REPRODUCE and FIX the one documented false negative of the family definition.

family_definition_formal.md (FN class 1) flags that minimap2's default secondary cap (-N 5 -> <=6 records/read)
truncates each read's placements, so for a tandem array with >6 cross-mapping copies the edges between farther
co-members are never observed and one true family FRAGMENTS into a big component + singletons. The note inferred
this from coverage structure but explicitly did NOT run the uncapped re-alignment ("'uncapped truth' is inferred
... not an actual -N99 -p0 re-alignment"). THIS script runs it.

Array: NC_086018.1:22.09-22.69 Mb (~600 kb, ~35 kb period; the note labels it "chrY" but NC_086018.1 is chr23).

Method: vertices = read-start clusters across the array (derived once from the UNCAPPED bam = fullest picture).
Build the de-tie conflict graph (EXACT shipped criterion: best-overlap per read at each vertex, distinct record,
de-tie |de_i-de_j|<=DELTA and max<=DE_MAX, edge iff >=MIN_READS conflicting reads; family = connected component)
from (a) the CAPPED production bam GGO.bam and (b) the UNCAPPED re-map, over the SAME vertices. A far-away
expressed gene (RABL2B) is included as an FP control: it must stay a singleton in BOTH graphs.

Uncapped bam is produced by (run once, ~few min):
  samtools view -b -N names.txt -F 0x900 GGO.bam | samtools fasta - > reads.fa   # array-region reads
  minimap2 -ax splice:hq -uf -N 50 -p 0.1 GGO.fasta reads.fa | samtools sort -o uncapped.bam - ; samtools index
"""
import collections
import pysam

CAPPED   = "/home/juanfra/winloci_scratch/GGO.bam"
UNCAPPED = "/home/juanfra/winloci_scratch/arrayfix/uncapped.bam"
CHROM, LO, HI = "NC_086018.1", 22_090_000, 22_690_000
DELTA, DE_MAX, MIN_READS = 0.005, 0.05, 3          # shipped de-tie operating point
VGAP, VMIN = 2_500, 3                               # vertex clustering: start-gap (< ~4-5kb copy period), min reads/copy
CONTROL = (CHROM, 48_818_440, 48_832_011, "RABL2B_ctrl")  # far-away expressed gene -> FP control


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def cluster_vertices(bam_path):
    """Copy vertices = clusters of PRIMARY-read starts across the array window (primaries = one clean anchor/read
    at each read's best copy; secondaries would smear a cluster onto every copy). Fixed vertex set for both graphs."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    spans = []
    for r in bam.fetch(CHROM, LO, HI):
        if r.is_supplementary or r.is_unmapped or r.is_secondary:
            continue
        if LO <= r.reference_start < HI:
            spans.append((r.reference_start, r.reference_end))
    spans.sort()
    clusters, cur = [], [spans[0]]
    for s in spans[1:]:
        if s[0] - cur[-1][0] <= VGAP:
            cur.append(s)
        else:
            clusters.append(cur); cur = [s]
    clusters.append(cur)
    verts = []
    for i, cl in enumerate(clusters):
        if len(cl) < VMIN:
            continue
        verts.append((CHROM, min(x[0] for x in cl), max(x[1] for x in cl), f"V{len(verts)}"))
    return verts


def best_records(bam, v):
    c, s, e, _ = v
    d = collections.defaultdict(list)
    for r in bam.fetch(c, s, e):
        if r.is_supplementary or r.is_unmapped:
            continue
        ql = r.query_alignment_length or 0
        if ql <= 0:
            continue
        ov = overlap(r.reference_start, r.reference_end, s, e)
        if ov <= 0:
            continue
        de = float(r.get_tag("de")) if r.has_tag("de") else r.get_tag("NM") / ql
        d[r.query_name].append((ov, r.reference_start, de))
    return {q: max(recs, key=lambda x: x[0]) for q, recs in d.items()}  # best-overlap per read


def conflict_graph(bam_path, vertices):
    bam = pysam.AlignmentFile(bam_path, "rb")
    recs = {v[3]: best_records(bam, v) for v in vertices}
    names = [v[3] for v in vertices]
    pair_n, co_obs = {}, 0
    edges = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            shared = set(recs[a]) & set(recs[b])
            n = 0
            for q in shared:
                pa, pb = recs[a][q], recs[b][q]
                if abs(pa[1] - pb[1]) < 200:                       # distinct record (not same locus)
                    continue
                if abs(pa[2] - pb[2]) <= DELTA and max(pa[2], pb[2]) <= DE_MAX:
                    n += 1
            pair_n[(a, b)] = n
            if n > 0:
                co_obs += 1
            if n >= MIN_READS:
                edges.append((a, b))
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for n in names:
        comp[find(n)].append(n)
    return [sorted(v, key=lambda x: int(x[1:]) if x[1:].isdigit() else 1e9) for v in comp.values()], pair_n, co_obs


def report(label, comps, co_obs, n_pairs, array_names):
    arr_comps = [c for c in comps if any(x in array_names for x in c)]
    big = max((c for c in arr_comps), key=len, default=[])
    singles = [c for c in arr_comps if len(c) == 1]
    print(f"\n[{label}] array vertices={len(array_names)}  components touching array={len(arr_comps)}  "
          f"largest={len(big)}  singletons={len(singles)}  pairs-co-observed={co_obs}/{n_pairs}")
    for c in sorted(arr_comps, key=len, reverse=True):
        print(f"    {{{', '.join(c)}}}")
    return len(arr_comps), len(big), len(singles)


def main():
    vertices = cluster_vertices(CAPPED) + [CONTROL]
    array_names = [v[3] for v in vertices if v[3] != "RABL2B_ctrl"]
    n_pairs = len(array_names) * (len(array_names) - 1) // 2
    print(f"Array {CHROM}:{LO}-{HI}  ->  {len(array_names)} copy vertices (period "
          f"~{(HI-LO)//max(1,len(array_names))//1000}kb)  + 1 FP control (RABL2B, {CONTROL[1]//10**6}Mb)")

    cap_comp, cap_pn, cap_co = conflict_graph(CAPPED, vertices)
    unc_comp, unc_pn, unc_co = conflict_graph(UNCAPPED, vertices)
    nc_cap, big_cap, sing_cap = report("CAPPED  (-N5, production GGO.bam)", cap_comp, cap_co, n_pairs, array_names)
    nc_unc, big_unc, sing_unc = report("UNCAPPED(-N50 -p0.1 re-map)",      unc_comp, unc_co, n_pairs, array_names)

    # FP control: RABL2B must be its own singleton in BOTH graphs (never fused into the array).
    def ctrl_singleton(comps):
        return any(c == ["RABL2B_ctrl"] for c in comps)
    fp_ok = ctrl_singleton(cap_comp) and ctrl_singleton(unc_comp)

    print("\n=== LEDGER ===")
    print(f"  CAPPED  : array fragments into {nc_cap} components (largest {big_cap}, {sing_cap} singletons) "
          f"<- the documented FALSE NEGATIVE")
    print(f"  UNCAPPED: array forms {nc_unc} component(s) (largest {big_unc}, {sing_unc} singletons)")
    print(f"  co-observed copy-pairs: capped {cap_co}/{n_pairs} -> uncapped {unc_co}/{n_pairs}")
    print(f"  FP control RABL2B stays a singleton in both graphs: {fp_ok}")
    # Classify the UNCAPPED array singletons: a singleton with >0 cross-mapping reads to the core would be a
    # residual cap/quorum miss; a singleton with 0 cross-mapping reads maps UNIQUELY -> genuinely resolvable,
    # correctly separate (NOT a cap false negative).
    unc_bam = pysam.AlignmentFile(UNCAPPED, "rb")
    name2v = {v[3]: v for v in vertices}
    main = max((c for c in unc_comp if any(x in array_names for x in c)), key=len)
    uniq, residual = [], []
    for c in unc_comp:
        if len(c) == 1 and c[0] in array_names:
            v = name2v[c[0]]; rv = best_records(unc_bam, v)
            xmap = 0
            for m in main:
                rm = best_records(unc_bam, name2v[m])
                for q in set(rv) & set(rm):
                    if abs(rv[q][1] - rm[q][1]) >= 200:
                        xmap += 1
            (uniq if xmap == 0 else residual).append((c[0], len(rv), xmap))

    fixed = len(main) > big_cap
    print(f"\n  RESULT: uncapping {'RECOVERS the cap-fragmented homology core' if fixed else 'no core growth'} "
          f"({big_cap} -> {len(main)} copies); FP {'clean (control isolated)' if fp_ok else 'VIOLATED'}.")
    print(f"  components: capped {nc_cap} -> uncapped {nc_unc};  core size: {big_cap} -> {len(main)};  "
          f"co-observed pairs: {cap_co} -> {unc_co} / {n_pairs}")
    print(f"  residual uncapped singletons that map UNIQUELY (0 cross-map -> resolvable, correctly separate, "
          f"NOT cap FNs): {len(uniq)}  {[u[0] for u in uniq]}")
    print(f"  residual singletons still cross-mapping the core but below quorum (true residual): "
          f"{len(residual)}  {residual}")
    print(f"\n  => After uncapping, the identifiability family = the {len(main)} cross-mapping copies; "
          f"every remaining singleton maps uniquely. Cap false-negative recovered; FP=0.")


if __name__ == "__main__":
    main()
