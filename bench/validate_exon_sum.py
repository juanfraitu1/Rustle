#!/usr/bin/env python3
"""validate_exon_sum.py — ANNOTATION-FREE, de-circularising validation of the genome-wide conflict
family catalog by FLNC exon-sum alignment.

The family DEFINITION is the de-tie read-conflict graph (loci among which reads are genuinely confused).
That fires on ANY shared region — including a single shared exon/domain — so it can transitively merge a
real paralog pair with a copy that only shares a domain (the GSTM2+SEC22B contamination seen when copies
were mapped to RefSeq gene names). This validator does NOT use annotation and does NOT use the conflict
graph. It takes each copy's exon-sum (the spliced, FLNC-derived consensus sequence `<prefix>.copies.fa`)
and all-vs-all aligns a family's copies with minimap2. A copy is CONFIRMED iff its spliced sequence aligns
FULL-LENGTH (identity >= ID, coverage-of-shorter >= COV) to a sibling — a genuinely independent test,
because read-conflict can fire on a fragment while full-length exon-sum alignment cannot.

  purity(family) = (largest mutually-aligning component) / (#copies)   [1.0 = pure paralog family]

Run: python bench/validate_exon_sum.py [prefix] [minimap2] [ID] [COV]
"""
import sys, subprocess, collections, os, tempfile, shutil

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/gw_xchrom_catalog"
MM2 = sys.argv[2] if len(sys.argv) > 2 else "minimap2"
ID = float(sys.argv[3]) if len(sys.argv) > 3 else 0.80   # min gap-compressed identity (1 - de)
COV = float(sys.argv[4]) if len(sys.argv) > 4 else 0.50  # min aligned fraction of the SHORTER sequence


def abort(msg):
    """Hard, loud, nonzero. A validator that cannot fail certifies nothing (o1_ledger.md §6am)."""
    print(f"VALIDATE_EXON_SUM ABORT: {msg}", file=sys.stderr)
    sys.exit(2)


# GUARD (§6am): a missing or present-but-EMPTY input used to yield an empty evidence set that still
# printed purity rates and the OVER-MERGES verdict. Abort BEFORE any alignment and BEFORE
# {PREFIX}.refined.tsv is truncated further down.
for _p in (f"{PREFIX}.copies.fa", f"{PREFIX}.families.tsv"):
    if not os.path.isfile(_p):
        abort(f"required input missing: {_p}")
    if os.path.getsize(_p) == 0:
        abort(f"required input is EMPTY (zero-length evidence scores as a pass): {_p}")
# GUARD (§6am): the aligner IS the evidence. If it is not executable every family returns 0 edges,
# which reads as purity 1/n and prints OVER-MERGES from zero alignments. Resolved once, up front.
if shutil.which(MM2) is None and not (os.path.isfile(MM2) and os.access(MM2, os.X_OK)):
    abort(f"minimap2 not found / not executable: {MM2!r} (pass an explicit path as argv[2])")

# --- load exon-sum sequences, grouped by family ---  header: >fid|ci|chrom:start-end|strand|nexon=N
fam_seqs = collections.defaultdict(dict)   # fid -> {ci: seq}
fam_locus = collections.defaultdict(dict)  # fid -> {ci: (chrom, start, end, strand)}
fid = ci = None
seq = []
def flush():
    if fid is not None and seq:
        fam_seqs[fid][ci] = "".join(seq)
for line in open(f"{PREFIX}.copies.fa"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        flush()
        h = line[1:].split("|")
        fid, ci = h[0], h[1]
        chrom, span = h[2].rsplit(":", 1)
        s, e = span.split("-")
        fam_locus[fid][ci] = (chrom, int(s), int(e), h[3])
        seq = []
    else:
        seq.append(line)
flush()

fams = {}
for line in open(f"{PREFIX}.families.tsv"):
    if line.startswith("family_id"):
        continue
    f = line.rstrip("\n").split("\t")
    fams[f[0]] = {"n_copies": int(f[1]), "cross_chrom": f[4] == "true"}

# GUARD (§6am): the section-3 flagship used to be skipped in SILENCE when absent, and align_family()
# returns an empty set below 2 copies — so a 0- or 1-copy family printed "CLIQUE" (0 edges trivially
# equals the 0 edges a clique of that size requires). The flagship claim is load-bearing (quoted in
# bench/COPY_ASSIGNMENT_AND_GATE.md), so it is checked HERE, before any alignment and before
# {PREFIX}.refined.tsv is truncated below.
RABL2_FID = os.environ.get("RABL2_FID", "GWFAM50")
if RABL2_FID not in fam_seqs:
    abort(f"flagship family {RABL2_FID} absent from {PREFIX}.copies.fa — the clique test cannot be run "
          f"(set RABL2_FID=<family_id> if this catalog numbers it differently)")
if len(fam_seqs[RABL2_FID]) < 2:
    abort(f"flagship family {RABL2_FID} has {len(fam_seqs[RABL2_FID])} copy — the clique test is VACUOUS "
          f"below 2 copies (0 edges trivially equals the 0 edges a clique requires)")

def align_family(seqs):
    """seqs: {ci: sequence}. Return set of edges {(ci_a, ci_b)} where the two spliced sequences align
    full-length (identity>=ID, coverage-of-shorter>=COV). minimap2 tries both strands."""
    cis = list(seqs)
    if len(cis) < 2:
        return set()
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        path = fh.name
        for c in cis:
            fh.write(f">{c}\n{seqs[c]}\n")
    try:
        # asm20: tuned for similar assemblies (<=~20% divergence) — these are de-tied multimappers, so
        # the copies are by construction similar enough that reads confused between them. -X = all-vs-all
        # skipping self & dual; -c emits the de:f gap-compressed-divergence tag.
        res = subprocess.run([MM2, "-cx", "asm20", "-X", "--no-long-join", path, path],
                             capture_output=True, text=True)
    finally:
        os.unlink(path)
    # GUARD (§6am): the aligner ran unchecked — a per-family minimap2 failure returned empty stdout,
    # hence 0 edges, hence purity 1/n and a printed OVER-MERGES verdict produced by NO alignment.
    if res.returncode != 0:
        abort(f"minimap2 exited {res.returncode} on a {len(cis)}-copy family "
              f"(a failed aligner silently scores as 0 edges): {res.stderr.strip()[-500:]}")
    out = res.stdout
    edges = set()
    lens = {c: len(seqs[c]) for c in cis}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
        if q == t:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else (int(f[9]) / max(1, int(f[10])))
        cov = (qe - qs) / max(1, min(lens.get(q, ql), lens.get(t, tl)))
        if ident >= ID and cov >= COV:
            edges.add(tuple(sorted((q, t))))
    return edges

def distinct_loci(comp, locus):
    """Collapse copies at the SAME locus, return the number of spatially-DISTINCT loci in `comp`.
    Principled & threshold-free: distinct paralog copies occupy DISJOINT genomic spans, so two copies are
    the SAME locus iff their spans OVERLAP on the same (chrom, strand) (a gene + its own nested fragment /
    a second isoform). This is the distinct-LOCUS guarantee that exon-sum homology alone cannot provide
    (two isoforms of one gene trivially align). Returns (n_distinct_loci, [representative_ci per locus])."""
    parent = {c: c for c in comp}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for i in range(len(comp)):
        for j in range(i + 1, len(comp)):
            a, b = comp[i], comp[j]
            ca, sa, ea, ta = locus[a]
            cb, sb, eb, tb = locus[b]
            if ca == cb and ta == tb and min(ea, eb) > max(sa, sb):  # same chrom+strand, spans overlap
                union(a, b)
    groups = collections.defaultdict(list)
    for c in comp:
        groups[find(c)].append(c)
    # representative per locus = the widest-span copy (most complete)
    reps = [max(g, key=lambda c: locus[c][2] - locus[c][1]) for g in groups.values()]
    return len(groups), reps

def components(nodes, edges):
    """Connected components of the exon-sum-homology graph (each = a refined paralog family)."""
    adj = {n: set() for n in nodes}
    for a, b in edges:
        if a in adj and b in adj:
            adj[a].add(b); adj[b].add(a)
    seen, comps = set(), []
    for s in nodes:
        if s in seen:
            continue
        comp, stack = [], [s]; seen.add(s)
        while stack:
            u = stack.pop(); comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); stack.append(v)
        comps.append(comp)
    return comps

def largest_component(nodes, edges):
    cs = components(nodes, edges)
    return max((len(c) for c in cs), default=0)

def pct(a, b):
    return f"{a}/{b} = {100*a/b:.1f}%" if b else "n/a"
def med(xs):
    return sorted(xs)[len(xs)//2] if xs else 0.0

n = pure = pure_xc = n_xc = total_edges = 0
purities, pur_same, pur_xc = [], [], []
contaminated = []
refined = []  # (orig_fid, comp_index, [ci...]) for each >=2-copy homology component
for f, seqs in fam_seqs.items():
    nodes = list(seqs)
    if len(nodes) < 2:
        continue
    edges = align_family(seqs)
    total_edges += len(edges)
    comps = components(nodes, edges)
    lc = max((len(c) for c in comps), default=0)
    pur = lc / len(nodes)
    cross = fams.get(f, {}).get("cross_chrom", False)
    n += 1
    n_xc += cross
    purities.append(pur)
    (pur_xc if cross else pur_same).append(pur)
    if pur >= 0.999:
        pure += 1
        pure_xc += cross
    elif pur < 0.6:
        contaminated.append((f, cross, len(nodes), lc, round(pur, 2)))
    for k, comp in enumerate(comps):
        if len(comp) < 2:
            continue
        nloci, reps = distinct_loci(comp, fam_locus[f])   # require >=2 spatially-DISTINCT loci
        if nloci >= 2:                                     # homology-validated AND multi-LOCUS
            refined.append((f, k, reps))

# GUARD (§6am): guard the EVIDENCE, not just the files. With no scorable family there is nothing to
# be pure or impure about, and with zero alignments genome-wide every family scores purity 1/n by
# construction — both print rates and the OVER-MERGES verdict from an empty computation. Abort here,
# BEFORE {PREFIX}.refined.tsv (a committed artefact) is truncated below.
if n == 0:
    abort(f"no family in {PREFIX}.copies.fa has >=2 copies — nothing to validate")
if total_edges == 0:
    abort(f"0 exon-sum alignment edges across all {n} scored families — the alignment evidence set is "
          f"EMPTY, so every purity is 1/n by construction (check minimap2 and {PREFIX}.copies.fa)")

print(f"=== EXON-SUM (FLNC) validation of {PREFIX} ===")
print(f"   annotation-free, de-circularised; minimap2 asm20, identity>={ID}, coverage-of-shorter>={COV}\n")
print(f"--- 1. PURITY of the raw conflict-graph catalog (how often its copies are really all homologous) ---")
print(f"families with >=2 copies: {n}  (same-chrom: {n - n_xc}, cross-chrom: {n_xc})")
print(f"PURE (all copies one homology component): all {pct(pure, n)} | "
      f"same-chrom {pct(pure - pure_xc, n - n_xc)} | cross-chrom {pct(pure_xc, n_xc)}")
print(f"median purity: all={med(purities):.2f}  same-chrom={med(pur_same):.2f}  cross-chrom={med(pur_xc):.2f}")
# GUARD (§6am): this verdict used to print unconditionally, i.e. it could not come out the other way.
# It now states only what the purities above actually show.
if pure < n:
    print(f"\n  => the raw cross-chrom conflict graph OVER-MERGES: read-conflict fires on shared repeats/domains,")
    print(f"     so many cross-chrom 'families' have copies whose full spliced sequences do NOT mutually align.")
else:
    print(f"\n  => NO over-merge detected: all {n} families are a single exon-sum-homology component,")
    print(f"     i.e. this run does NOT support the OVER-MERGES claim.")
print(f"\nmost-contaminated (purity<0.6) — read-conflict edges the exon-sum alignment rejects:")
for f, cross, nc, lc, pur in sorted(contaminated, key=lambda x: x[4])[:12]:
    print(f"  {f} xc={cross} n_copies={nc} largest_homologous={lc} purity={pur}")

# --- 2. the homology-REFINED catalog (split each family into exon-sum-homology components) ---
ref_xc = 0
sizes = collections.Counter()
ref_rows = []
for orig, k, comp in refined:
    chroms = set(fam_locus[orig][c][0] for c in comp)
    cross = len(chroms) > 1
    ref_xc += cross
    sizes[len(comp)] += 1
    ref_rows.append((orig, k, comp, sorted(chroms), cross))
print(f"\n--- 2. REFINED catalog (exon-sum-homology component AND >=2 spatially-DISTINCT loci) ---")
print(f"refined families: {len(refined)}  (cross-chrom: {ref_xc})   [raw conflict families: {n}]")
print(f"size distribution: " + ", ".join(f"{s}cp×{sizes[s]}" for s in sorted(sizes)))

with open(f"{PREFIX}.refined.tsv", "w") as fh:
    fh.write("refined_id\torig_family\tn_copies\tn_chroms\tchroms\tcross_chrom\tcopies\n")
    for i, (orig, k, comp, chroms, cross) in enumerate(sorted(ref_rows, key=lambda r: -len(r[2]))):
        loc = ";".join(f"{fam_locus[orig][c][0]}:{fam_locus[orig][c][1]}-{fam_locus[orig][c][2]}" for c in comp)
        fh.write(f"REF{i}\t{orig}\t{len(comp)}\t{len(chroms)}\t{','.join(chroms)}\t{cross}\t{loc}\n")
print(f"wrote {PREFIX}.refined.tsv")

# --- 3. RABL2 flagship ---
print("\n=== RABL2 flagship (cross-chrom) ===")
rabl2_fid = RABL2_FID   # presence and >=2 copies already hard-checked at load time (GUARD §6am above)
seqs = fam_seqs[rabl2_fid]
edges = align_family(seqs)
comps = components(list(seqs), edges)
lc = max((len(c) for c in comps), default=0)
print(f"  {rabl2_fid}: {len(seqs)} copies, largest homology component = {lc}, purity = {lc/len(seqs):.2f}")
print(f"  loci: {[f'{fam_locus[rabl2_fid][c][0]}:{fam_locus[rabl2_fid][c][1]}' for c in seqs]}")
print(f"  all-pairs aligned (clique): {len(edges)} edges among {len(seqs)} copies "
      f"({'CLIQUE' if len(edges) == len(seqs)*(len(seqs)-1)//2 else 'not full clique'})")
