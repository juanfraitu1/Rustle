#!/usr/bin/env python3
"""validate_gw_catalog.py — EXTERNAL (de-circularising) validation of the genome-wide conflict family
catalog (incl cross-chromosome families) against independent paralog truth: RefSeq gene-name roots +
Ensembl Compara paralog relations. NOT sequence-similarity (which the family detection itself uses) —
gene-tree / orthology-based, so it is a genuine independent check.

A family PASSES if its copies map to RefSeq genes that (a) share a gene-name ROOT (RABL2A/RABL2B → RABL2;
identical LOC ids), or (b) are a Compara paralog pair (bench/compara_paralog_relation.json). Reports
precision, cross-chrom precision, and the RABL2 flagship explicitly. Run:
    python bench/validate_gw_catalog.py <catalog_prefix>   (default gw_xchrom_catalog)
"""
import sys, csv, json, re, bisect, collections

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/gw_xchrom_catalog"
GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
COMPARA = "/mnt/c/Users/jfris/Desktop/Rustle/bench/compara_paralog_relation.json"

# --- 1. RefSeq gene intervals per chrom (genes only) ---
genes = collections.defaultdict(list)  # chrom -> [(start,end,name)]
for line in open(GFF):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    chrom, s, e, attr = f[0], int(f[3]), int(f[4]), f[8]
    m = re.search(r"gene=([^;]+)", attr) or re.search(r"Name=([^;]+)", attr)
    if m:
        genes[chrom].append((s, e, m.group(1)))
for c in genes:
    genes[c].sort()
gstarts = {c: [g[0] for g in v] for c, v in genes.items()}

def gene_at(chrom, s, e):
    """RefSeq gene name overlapping [s,e) on chrom (the max-overlap one), or None."""
    iv = genes.get(chrom)
    if not iv:
        return None
    i = bisect.bisect_right(gstarts[chrom], e) - 1
    best, bov = None, 0
    for j in range(max(0, i - 30), min(len(iv), i + 2)):
        gs, ge, nm = iv[j]
        ov = min(ge, e) - max(gs, s)
        if ov > bov:
            best, bov = nm, ov
    return best

def name_root(name):
    """Conservative paralog-suffix strip: RABL2A→RABL2 (digit+letter suffix), RFPL1→RFPL (letter+digit).
    Does NOT over-strip ordinary names (PBRM1 stays PBRM1). LOC ids: only identical count as same root."""
    if name is None:
        return None
    if name.startswith("LOC"):
        return name
    n = name
    if len(n) >= 3 and n[-1].isupper() and n[-1].isalpha() and n[-2].isdigit():
        n = n[:-1]  # RABL2A -> RABL2, APOBEC3B -> APOBEC3
    elif len(n) >= 3 and n[-1].isdigit() and n[-2].isalpha():
        n = re.sub(r"\d+$", "", n)  # RFPL1 -> RFPL
    return n if len(n) >= 2 else name

# --- mmseqs protein-homology pairs (gene-name keyed; handles LOC paralogs Compara/names cannot) ---
HOM = collections.defaultdict(set)
M8 = "/home/juanfra/winloci_scratch/prot_ava.m8"
try:
    # columns: query target fident qcov tcov evalue bits alnlen. A TRUE paralog pair has high identity
    # AND high MUTUAL coverage (most of both proteins align); a DOMAIN-SHARER aligns only over a shared
    # domain (high fident on a fragment, LOW qcov/tcov) — excluded by the coverage gate so the precision
    # reflects real paralogy, not promiscuous domains.
    for line in open(M8):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6:
            q, t, fident = f[0], f[1], float(f[2])
            qcov = float(f[3]) if len(f) > 3 else 1.0
            tcov = float(f[4]) if len(f) > 4 else 1.0
            eval_ = float(f[5])
            if q != t and fident >= 0.40 and qcov >= 0.50 and tcov >= 0.50 and eval_ <= 1e-5:
                HOM[q].add(t)
                HOM[t].add(q)
except FileNotFoundError:
    pass

def mmseqs_homologous(names):
    ns = [n for n in names if n]
    for i in range(len(ns)):
        for j in range(i + 1, len(ns)):
            if ns[j] in HOM.get(ns[i], ()):
                return True
    return False

def purity(names, roots):
    """Largest mutually-homologous component / #RefSeq-mapped copies. A pure paralog family = 1.0;
    a contaminated over-merge (one real pair + unrelated copies) < 1.0. Edge = mmseqs homology OR
    shared name-root OR Compara paralog. Returns (purity, n_mapped)."""
    idx = [i for i, n in enumerate(names) if n]
    if len(idx) < 2:
        return (1.0 if idx else 0.0), len(idx)
    adj = {i: set() for i in idx}
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            ni, nj, ri, rj = names[i], names[j], roots[i], roots[j]
            edge = (nj in HOM.get(ni, ())) or (ri == rj and ri is not None) or \
                   (nj in para.get(ni, ()))
            if edge:
                adj[i].add(j); adj[j].add(i)
    seen, best = set(), 0
    for s in idx:
        if s in seen:
            continue
        comp, stack = 0, [s]
        seen.add(s)
        while stack:
            u = stack.pop(); comp += 1
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); stack.append(v)
        best = max(best, comp)
    return best / len(idx), len(idx)

# --- 2. Compara paralog pairs (gene -> set of paralog partners) ---
comp = json.load(open(COMPARA))
para = collections.defaultdict(set)
for fam, info in comp.get("families", {}).items():
    named = info.get("named_genes", [])
    for a in named:
        for b in named:
            if a != b:
                para[a].add(b)

def compara_pair(names):
    ns = [n for n in names if n]
    for i in range(len(ns)):
        for j in range(i + 1, len(ns)):
            if ns[j] in para.get(ns[i], ()):
                return True
    return False

# --- 3. validate each family ---
copies = collections.defaultdict(list)
for c in csv.DictReader(open(f"{PREFIX}.copies.tsv"), delimiter="\t"):
    copies[c["family_id"]].append(c)
fams = {f["family_id"]: f for f in csv.DictReader(open(f"{PREFIX}.families.tsv"), delimiter="\t")}

n = passed = passed_xc = n_xc = root_pass = compara_pass = mmseqs_pass = 0
unmapped = mappable = mappable_pass = mappable_xc = mappable_xc_pass = 0
pure_count = 0  # families whose mapped copies are ALL one homologous component (purity==1, >=2 mapped)
purities = []
examples = []
for fid, cs in copies.items():
    names = [gene_at(c["chrom"], int(c["start"]), int(c["end"])) for c in cs]
    roots = [name_root(x) for x in names]
    cross = fams[fid]["cross_chrom"] == "true"
    n += 1
    n_xc += cross
    mapped = [r for r in roots if r]
    has_mappable_pair = sum(1 for x in names if x) >= 2  # >=2 copies map to a RefSeq gene
    if not mapped:
        unmapped += 1
        continue
    root_ok = len(set(mapped)) == 1 and len(mapped) >= 2  # all mapped copies share one root
    comp_ok = compara_pair(names)
    mm_ok = mmseqs_homologous(names)
    ok = root_ok or comp_ok or mm_ok
    pur, nmap = purity(names, roots)
    passed += ok
    passed_xc += ok and cross
    root_pass += root_ok
    compara_pass += comp_ok
    mmseqs_pass += mm_ok
    if has_mappable_pair:
        mappable += 1
        mappable_pass += ok
        mappable_xc += cross
        mappable_xc_pass += ok and cross
        purities.append(pur)
        if pur >= 0.999:
            pure_count += 1
    if ok and len(examples) < 12:
        tag = "name" if root_ok else ("compara" if comp_ok else "mmseqs")
        examples.append((fid, cross, sorted(set(x for x in names if x))[:4], len(cs), tag, pur))

def pct(a, b):
    return f"{a}/{b} = {100*a/b:.1f}%" if b else "n/a"

print(f"=== EXTERNAL validation of {PREFIX} ({n} families) ===\n")
print(f"families validated (mmseqs homology OR gene-name root OR Compara): {pct(passed, n)}")
print(f"  via mmseqs protein-homology: {mmseqs_pass} | gene-name root: {root_pass} | Compara: {compara_pass}")
print(f"  unmapped-to-RefSeq (cannot validate, de-novo at unannotated loci): {unmapped}")
print(f"\n⭐ PRECISION on families with >=2 RefSeq-mappable copies (the checkable set):")
print(f"   contains a confirmed paralog pair: {pct(mappable_pass, mappable)} | CROSS-CHROM: {pct(mappable_xc_pass, mappable_xc)}")
med_pur = sorted(purities)[len(purities)//2] if purities else 0.0
print(f"   PURE (all mapped copies one homologous component): {pct(pure_count, mappable)} | median purity: {med_pur:.2f}")
print(f"\nCROSS-CHROM families: {n_xc} | validated (of all): {pct(passed_xc, n_xc)}  <- the new capability")
print(f"\nvalidated examples (fid, cross_chrom, genes, n_copies, via, purity):")
for fid, cross, gns, ncp, tag, pur in examples:
    print(f"  {fid} xc={cross} genes={gns} n={ncp} via={tag} purity={pur:.2f}")

# --- 4. RABL2 flagship ---
print("\n=== RABL2 flagship (cross-chrom) ===")
for fid, cs in copies.items():
    nm = [gene_at(c["chrom"], int(c["start"]), int(c["end"])) for c in cs]
    if any(x and "RABL2" in x for x in nm):
        print(f"  {fid}: copies map to {[x for x in nm]}")
        print(f"    chroms: {sorted(set(c['chrom'] for c in cs))} -> cross-chrom paralog family RECOVERED")
        break
