#!/usr/bin/env python3
"""KEY=dupfamboundary, Phase R2 (T1-d): the gorilla SEDEF arm. Descriptive, no verdict class.

Spec: docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md section 5 and Amendment 1 (git 17ad07cd).

    t1_gorilla_pairs.py prepare   build the retained SEDEF tables, the atoms and the atom classes at tau 0.90, 0.95, 0.98
    t1_gorilla_pairs.py counts    the input-only counts of Gate 6 / section 5 and the positive control (no S-rate, no report table)
    t1_gorilla_pairs.py gate0     size, row count and sha256 of the R2 inputs against the registered prefixes
    t1_gorilla_pairs.py report    the S-rate / A_viol / B_viol tables (REFUSES to run without a release file)

Reading decisions (the prereg text is followed literally; where it allows two readings the conservative one is used, see the build report):
  * mitochondrial rows are dropped first, then identical (contig1, start1, end1, contig2, start2, end2) tuples are collapsed (largest column 21, then
    larger column 12, then first in file), then rows are retained at fracMatch >= tau. (Filtering first and deduplicating second gives the same set.)
  * a class edge is an edge of edges.tsv, whose coverage is the minimum over both atoms (the pre-registered two-sided form), with coverage >= 0.5.
  * a gene touches a class when the SUM of its merged-exon bp over all atoms of the class is >= 100 (not >= 100 in one atom).
  * eligibility of a pair: both genes have a non-empty U, they differ, and their merged exon unions share 0 bp; genes on different contigs share 0 bp.
  * distance is the gap between the half-open gene spans (start1 - 1, end), 0 if they overlap; < 100 kb, 100 kb to < 1 Mb, >= 1 Mb, or cross-contig.
  * depth-matched: a pairs.tsv row (any family status present in the file) with identity (not identity_cs) >= tau, both genes in the 1to1 layer and
    in the row's family, the pair eligible, and the family in the verdict set; rows are resolved before the identity filter, so a row naming a gene
    outside the layer stops the run whatever its identity.

Everything below the CLI is pure functions on in-memory tables so that the unit tests can drive them with made-up data.
The atoms themselves come from bench/dna_sd_atoms.py (restored unchanged, sha256 9d96d976...) run as a subprocess in cigar mode.
"""
import argparse
import bisect
import collections
import csv
import hashlib
import json
import math
import os
import platform
import random
import re
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
ATOMS_SCRIPT = os.path.join(REPO, "bench", "dna_sd_atoms.py")
ATOMS_SHA256 = "9d96d976eb52dca6e056b1bc8256f006835d72b265a8697c925d15cfcb562bce"

MITO_CONTIG = "NC_011120.1"
DEV_CONTIGS = ("NC_073241.2", "NC_073242.2", "NC_073244.2", "NC_073234.2")  # the SD x family development slice and the most-junctions contig
TAUS = (0.90, 0.95, 0.98)
TAU_PRIMARY = 0.90
CLASS_COVERAGE = 0.5      # a class edge needs two-sided (min over both atoms) CIGAR coverage >= 0.5
TOUCH_BP = 100            # a gene touches a class iff >= 100 of its exonic bp lie in atoms of the class
MIN_PAIRS = 20            # a table is UNDERPOWERED with fewer than 20 depth-matched pairs ...
MIN_FAMILIES = 8          # ... or fewer than 8 families
BOOTSTRAP_B = 2000
SEED = 20260930
N_COLUMNS = 34
COL_ALN_LEN = 11          # 0-based index of column 12
COL_FRAC = 20             # 0-based index of column 21


# ---------------------------------------------------------------------------------------------------------------------------
# the SEDEF table
# ---------------------------------------------------------------------------------------------------------------------------
class SedefRow(collections.namedtuple("SedefRow", "order line key frac aln_len")):
    __slots__ = ()


def read_sedef(lines):
    """Parse SEDEF text lines. Returns (rows, stats). Header lines (leading '#') and blank lines are skipped; rows that involve the
    mitochondrial contig NC_011120.1 on either side are dropped and counted; a data row without exactly 34 columns stops the run."""
    rows = []
    stats = {"data_rows": 0, "header_lines": 0, "mito_rows_dropped": 0}
    for ln in lines:
        ln = ln.rstrip("\n")
        if not ln.strip():
            continue
        if ln.startswith("#"):
            stats["header_lines"] += 1
            continue
        f = ln.split("\t")
        if len(f) != N_COLUMNS:
            raise ValueError(f"SEDEF row with {len(f)} columns (expected {N_COLUMNS}): {ln[:100]}")
        stats["data_rows"] += 1
        if f[0] == MITO_CONTIG or f[3] == MITO_CONTIG:
            stats["mito_rows_dropped"] += 1
            continue
        key = (f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]))
        rows.append(SedefRow(len(rows), ln, key, float(f[COL_FRAC]), int(f[COL_ALN_LEN])))
    return rows, stats


def dedupe(rows):
    """One row per identical (contig1, start1, end1, contig2, start2, end2): the largest fracMatch (column 21), ties by the larger aln_len
    (column 12), then the first in file. Survivors keep file order. Returns (kept, number removed)."""
    best = {}
    for r in rows:
        b = best.get(r.key)
        if b is None or (r.frac, r.aln_len) > (b.frac, b.aln_len):
            best[r.key] = r
    kept = sorted(best.values(), key=lambda r: r.order)
    return kept, len(rows) - len(kept)


def retain(rows, tau):
    """Rows with fracMatch >= tau."""
    return [r for r in rows if r.frac >= tau]


def write_sedef(rows, path):
    with open(path, "w") as fh:
        for r in rows:
            fh.write(r.line + "\n")


# ---------------------------------------------------------------------------------------------------------------------------
# atoms and atom classes
# ---------------------------------------------------------------------------------------------------------------------------
def check_atoms_script(expected_sha256=ATOMS_SHA256, path=ATOMS_SCRIPT):
    """The restored script must be the registered one, byte for byte."""
    with open(path, "rb") as fh:
        h = hashlib.sha256(fh.read()).hexdigest()
    if h != expected_sha256:
        raise RuntimeError(f"{path} has sha256 {h}, registered {expected_sha256}")
    return h


def run_atoms(sd_path, contigs_path, out_prefix, python=None):
    """bench/dna_sd_atoms.py in cigar mode, format 'gorilla' (34 columns). Returns its stderr text."""
    check_atoms_script()
    p = subprocess.run([python or sys.executable, "-B", ATOMS_SCRIPT, sd_path, "gorilla", "@" + contigs_path, out_prefix, "cigar"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"dna_sd_atoms.py failed ({p.returncode}): {p.stderr[-500:]}")
    return p.stderr


def read_nodes(prefix):
    """[(contig, start, end)] in atom-index order."""
    out = []
    with open(prefix + ".nodes.tsv") as fh:
        next(fh)
        for k, ln in enumerate(fh):
            idx, c, s, e = ln.rstrip("\n").split("\t")
            if int(idx) != k:
                raise ValueError(f"atom index {idx} at position {k}")
            out.append((c, int(s), int(e)))
    return out


def read_edges(prefix):
    """[(i, j, identity, coverage)]; coverage is the minimum over both atoms (the pre-registered two-sided form)."""
    out = []
    with open(prefix + ".edges.tsv") as fh:
        next(fh)
        for ln in fh:
            i, j, ident, cov = ln.rstrip("\n").split("\t")
            out.append((int(i), int(j), float(ident), float(cov)))
    return out


def build_classes(n_atoms, edges, min_cov=CLASS_COVERAGE):
    """Union-find over the atom edges with two-sided coverage >= min_cov. The class id of an atom is the smallest atom index of its class.
    Returns (class_of, number of edges used)."""
    parent = list(range(n_atoms))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    used = 0
    for i, j, _ident, cov in edges:
        if cov >= min_cov:
            used += 1
            a, b = find(i), find(j)
            if a != b:
                if a < b:
                    parent[b] = a
                else:
                    parent[a] = b
    return [find(i) for i in range(n_atoms)], used


def class_stats(class_of):
    sizes = collections.Counter(class_of)
    return {"classes": len(sizes), "largest": max(sizes.values()) if sizes else 0, "singletons": sum(1 for v in sizes.values() if v == 1)}


class AtomIndex:
    """Overlap queries on atoms (half-open intervals) per contig."""

    def __init__(self, atoms):
        by = collections.defaultdict(list)
        for idx, (c, s, e) in enumerate(atoms):
            by[c].append((s, e, idx))
        self.by = {c: sorted(v) for c, v in by.items()}
        self.starts = {c: [x[0] for x in v] for c, v in self.by.items()}
        self.maxlen = {c: max(e - s for s, e, _ in v) for c, v in self.by.items()}

    def overlaps(self, contig, s, e):
        """[(atom index, overlap bp)] for atoms overlapping [s, e) by at least 1 bp, in atom-index order."""
        if contig not in self.by or e <= s:
            return []
        lo = bisect.bisect_left(self.starts[contig], s - self.maxlen[contig])
        hi = bisect.bisect_left(self.starts[contig], e)
        out = []
        for a_s, a_e, idx in self.by[contig][lo:hi]:
            ov = min(e, a_e) - max(s, a_s)
            if ov > 0:
                out.append((idx, ov))
        out.sort()
        return out


def touch_classes(blocks, contig, atom_index, class_of):
    """{class id: exonic bp of the gene in atoms of that class} for exon blocks [(start, end)) on `contig`."""
    t = collections.Counter()
    for s, e in blocks:
        for idx, ov in atom_index.overlaps(contig, s, e):
            t[class_of[idx]] += ov
    return dict(t)


def gene_U(touch, min_bp=TOUCH_BP):
    """U(g): the classes with at least `min_bp` exonic bp of the gene."""
    return frozenset(k for k, v in touch.items() if v >= min_bp)


def tau_tag(tau):
    return f"{int(round(tau * 100)):03d}"


def prepare_tau(rows, tau, outdir, reuse=False):
    """Retain rows at tau, write sd_geNNN.bed, run the atoms script, build the classes. `rows` should be deduplicated.
    With reuse=True the atoms are not recomputed when sd_geNNN.bed is byte-identical to what would be written and the outputs exist.
    Returns a summary dict (and leaves out<NNN>.{nodes,edges,edges_shorter}.tsv in outdir)."""
    tag = tau_tag(tau)
    kept = retain(rows, tau)
    os.makedirs(outdir, exist_ok=True)
    bed = os.path.join(outdir, f"sd_ge{tag}.bed")
    prefix = os.path.join(outdir, "out" + tag)
    text = "".join(r.line + "\n" for r in kept)
    done = os.path.exists(prefix + ".nodes.tsv") and os.path.exists(prefix + ".edges.tsv") and os.path.exists(bed)
    same = False
    if reuse and done:
        with open(bed) as fh:
            same = fh.read() == text
    contigs = sorted({k for r in kept for k in (r.key[0], r.key[3])})
    cf = os.path.join(outdir, f"contigs{tag}.txt")
    log = ""
    if not same:
        with open(bed, "w") as fh:
            fh.write(text)
        with open(cf, "w") as fh:
            fh.write("\n".join(contigs) + ("\n" if contigs else ""))
        log = run_atoms(bed, cf, prefix)
    atoms = read_nodes(prefix)
    edges = read_edges(prefix)
    class_of, used = build_classes(len(atoms), edges)
    st = class_stats(class_of)
    return {"tau": tau, "rows_retained": len(kept), "contigs": len(contigs), "atoms": len(atoms), "edges_total": len(edges), "edges_used": used,
            "classes": st["classes"], "largest_class_atoms": st["largest"], "singleton_classes": st["singletons"], "prefix": prefix,
            "reused": same, "log": log.strip().splitlines()}


# ---------------------------------------------------------------------------------------------------------------------------
# genes, the Compara layer, the pair universe
# ---------------------------------------------------------------------------------------------------------------------------
Gene = collections.namedtuple("Gene", "contig name type biotype start1 end blocks")
Layer = collections.namedtuple("Layer", "families gene_family gff_to_gene status n_families_all")
DIST_CLASSES = ("same<100kb", "same100kb-1Mb", "same>1Mb", "cross-contig")


def merge_blocks(blocks):
    """Sorted blocks with overlapping blocks merged (adjacent blocks stay apart; bp counts do not depend on it)."""
    out = []
    for s, e in sorted(blocks):
        if out and s < out[-1][1]:
            if e > out[-1][1]:
                out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return tuple(out)


def read_genes(lines):
    """families_gw/species/gorilla/genes.tsv. Returns (genes, index, duplicate keys); the join key is (contig, start1, end).
    Exon blocks are 0-based half-open, merged. A duplicated key keeps its first record in the index and is reported."""
    it = iter(lines)
    header = next(it).rstrip("\n").split("\t")
    col = {h: i for i, h in enumerate(header)}
    for need in ("contig", "name", "type", "biotype", "start1", "end", "exons"):
        if need not in col:
            raise ValueError(f"genes table lacks column {need}")
    genes, index, dups = [], {}, []
    for ln in it:
        ln = ln.rstrip("\n")
        if not ln.strip():
            continue
        f = ln.split("\t")
        f += [""] * (len(header) - len(f))
        blocks = []
        for tok in f[col["exons"]].split(","):
            if tok:
                a, b = tok.split("-")
                blocks.append((int(a), int(b)))
        g = Gene(f[col["contig"]], f[col["name"]], f[col["type"]], f[col["biotype"]], int(f[col["start1"]]), int(f[col["end"]]), merge_blocks(blocks))
        key = (g.contig, g.start1, g.end)
        if key in index:
            dups.append(key)
        else:
            index[key] = len(genes)
        genes.append(g)
    return genes, index, dups


def read_layer(lines, index):
    """gorilla_truth.tsv: the 1to1 rows, keyed by (gorilla_contig, start1, end1) into the gene table. A 1to1 row whose key is absent from
    the gene table, or a gene projected into two families, stops the run."""
    rd = csv.DictReader(lines, delimiter="\t")
    families = collections.defaultdict(list)
    gene_family, gff_to_gene, status, missing = {}, {}, {}, []
    for r in rd:
        status.setdefault(r["family_id"], r["family_status"])
        if r["member_status"] != "1to1":
            continue
        key = (r["gorilla_contig"], int(r["start1"]), int(r["end1"]))
        g = index.get(key)
        if g is None:
            missing.append(key)
            continue
        if g in gene_family and gene_family[g] != r["family_id"]:
            raise ValueError(f"gene {key} is projected into two families")
        gene_family[g] = r["family_id"]
        gff_to_gene[r["gorilla_gff_id"]] = g
        if g not in families[r["family_id"]]:
            families[r["family_id"]].append(g)
    if missing:
        raise ValueError(f"{len(missing)} 1to1 rows have no gene record, e.g. {missing[:3]}")
    return Layer({f: sorted(v) for f, v in families.items()}, gene_family, gff_to_gene, status, len(status))


def read_pairs(lines):
    """pairs.tsv rows as dicts; identity and cov_short become floats, NA or empty identity becomes None."""
    out = []
    for r in csv.DictReader(lines, delimiter="\t"):
        r = dict(r)
        r["identity"] = None if r["identity"] in ("", "NA") else float(r["identity"])
        r["cov_short"] = None if r.get("cov_short", "") in ("", "NA") else float(r["cov_short"])
        out.append(r)
    return out


def verdict_families(layer, genes, dev_contigs=DEV_CONTIGS):
    """Families none of whose 1to1 projected genes, eligible or not, lies on a development contig."""
    return frozenset(f for f, members in layer.families.items() if not any(genes[g].contig in dev_contigs for g in members))


def exon_overlap_bp(a, b):
    """Shared bp of two sorted, merged block lists."""
    i = j = tot = 0
    while i < len(a) and j < len(b):
        ov = min(a[i][1], b[j][1]) - max(a[i][0], b[j][0])
        if ov > 0:
            tot += ov
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def shared_exonic_bp(g, h):
    """Shared exonic bp of two genes: 0 when they lie on different contigs (coordinates are only comparable on one contig)."""
    if g.contig != h.contig:
        return 0
    return exon_overlap_bp(g.blocks, h.blocks)


def distance_class(g, h):
    """Gap between the nearest ends of the two gene spans (half-open, 0 if they overlap); cross-contig is its own class."""
    if g.contig != h.contig:
        return "cross-contig"
    gap = max(0, max(g.start1 - 1, h.start1 - 1) - min(g.end, h.end))
    if gap < 100_000:
        return "same<100kb"
    if gap < 1_000_000:
        return "same100kb-1Mb"
    return "same>1Mb"


def pair_eligible(g, h, genes, U):
    """Both genes eligible (U non-empty), g differs from h, and the exon unions share 0 bp."""
    return g != h and bool(U.get(g)) and bool(U.get(h)) and shared_exonic_bp(genes[g], genes[h]) == 0


def same_family_pairs(layer, genes, U, verdict=None):
    """[(family, g, h, distance class)] for every eligible pair of 1to1 genes of one family (whole or partial), g < h."""
    out = []
    for fam in sorted(layer.families):
        if verdict is not None and fam not in verdict:
            continue
        m = layer.families[fam]
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                if pair_eligible(m[i], m[j], genes, U):
                    out.append((fam, m[i], m[j], distance_class(genes[m[i]], genes[m[j]])))
    return out


def control_pairs(layer, genes, U, verdict=None):
    """[(family g, family h, g, h, distance class)] for every eligible pair of labelled genes of DIFFERENT families, g < h."""
    pool = sorted((g, fam) for fam in layer.families if verdict is None or fam in verdict for g in layer.families[fam] if U.get(g))
    out = []
    for i in range(len(pool)):
        for j in range(i + 1, len(pool)):
            (g, fg), (h, fh) = pool[i], pool[j]
            if fg != fh and pair_eligible(g, h, genes, U):
                out.append((fg, fh, g, h, distance_class(genes[g], genes[h])))
    return out


def _resolve_rows(rows, layer):
    out = []
    for r in rows:
        try:
            g, h = layer.gff_to_gene[r["a_gene"]], layer.gff_to_gene[r["b_gene"]]
        except KeyError as e:
            raise ValueError(f"pairs.tsv names a gene outside the Compara layer: {e}")
        if layer.gene_family[g] != r["family_id"] or layer.gene_family[h] != r["family_id"]:
            raise ValueError(f"pairs.tsv row of family {r['family_id']} names genes of {layer.gene_family[g]} and {layer.gene_family[h]}")
        out.append((r, g, h))
    return out


def eligible_rows(rows, layer, genes, U, verdict=None):
    """Rows of pairs.tsv (whole families) whose pair is eligible, with any identity, inside the verdict set when one is given."""
    out = []
    for r, g, h in _resolve_rows(rows, layer):
        if verdict is not None and r["family_id"] not in verdict:
            continue
        if pair_eligible(g, h, genes, U):
            out.append(dict(r, g=g, h=h, dist_class=distance_class(genes[g], genes[h])))
    return sorted(out, key=lambda d: (d["family_id"], d["a_gene"], d["b_gene"]))


def depth_matched(rows, layer, genes, U, tau, verdict=None):
    """Eligible rows of pairs.tsv with identity >= tau (NA excluded)."""
    return [d for d in eligible_rows(rows, layer, genes, U, verdict) if d["identity"] is not None and d["identity"] >= tau]


# ---------------------------------------------------------------------------------------------------------------------------
# the report statistics (called on made-up tables by the tests; on the real verdict set only after the independent code review)
# ---------------------------------------------------------------------------------------------------------------------------
Z95 = 1.959963984540054


def wilson(k, n, z=Z95):
    """Wilson score interval of k out of n; (None, None) for an empty table."""
    if n == 0:
        return (None, None)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def s_rate(pairs, U):
    """pairs: [(family, g, h)]. S(g, h) = 1 iff U(g) and U(h) intersect. Pair-weighted rate, family-weighted rate (mean over families of
    the family's pair share), A_viol = 1 - pair-weighted rate and the Wilson interval of the pair-weighted rate."""
    by = collections.defaultdict(list)
    for fam, g, h in pairs:
        by[fam].append(1 if (U[g] & U[h]) else 0)
    n = sum(len(v) for v in by.values())
    shared = sum(sum(v) for v in by.values())
    fam_rates = [sum(by[f]) / len(by[f]) for f in sorted(by)]
    pw = shared / n if n else None
    return {"n_pairs": n, "n_families": len(by), "shared": shared, "pair_weighted": pw,
            "family_weighted": (sum(fam_rates) / len(fam_rates)) if fam_rates else None,
            "a_viol": (1 - pw) if pw is not None else None, "wilson": wilson(shared, n), "by_family": {f: by[f] for f in sorted(by)}}


def b_viol(control, U):
    """(pairs of different families whose genes share a class, pairs): a ceiling-type quantity."""
    k = sum(1 for _fg, _fh, g, h, _c in control if U[g] & U[h])
    return (k, len(control))


def b_viol_by_class(control, U):
    """The same count as b_viol, split by the registered distance classes: {class: (pairs sharing a class, pairs)}, descriptive."""
    by = {}
    for _fg, _fh, g, h, c in control:
        k, n = by.get(c, (0, 0))
        by[c] = (k + (1 if U[g] & U[h] else 0), n + 1)
    return {c: by[c] for c in sorted(by)}


def underpowered(n_pairs, n_families):
    return n_pairs < MIN_PAIRS or n_families < MIN_FAMILIES


def percentile_interval(values, alpha_permille=50):
    """Percentile interval at 1 - alpha: sorted[floor(B * alpha/2)] and sorted[ceil(B * (1 - alpha/2)) - 1] (integer arithmetic)."""
    v = sorted(values)
    B = len(v)
    lo = (B * (alpha_permille // 2)) // 1000
    hi = -((-B * (1000 - alpha_permille // 2)) // 1000) - 1
    return (v[lo], v[hi])


def bootstrap_values(by_family, B=BOOTSTRAP_B, seed=SEED):
    """Family-cluster bootstrap: B resamples of the families with replacement (families sorted by id, one random.Random(seed) stream).
    Returns (pair-weighted rates, family-weighted rates), one value per resample."""
    fams = sorted(by_family)
    rng = random.Random(seed)
    pw, fw = [], []
    for _ in range(B):
        draw = [fams[rng.randrange(len(fams))] for _ in fams]
        n = sum(len(by_family[f]) for f in draw)
        pw.append(sum(sum(by_family[f]) for f in draw) / n)
        fw.append(sum(sum(by_family[f]) / len(by_family[f]) for f in draw) / len(draw))
    return pw, fw


def bootstrap_rates(by_family, B=BOOTSTRAP_B, seed=SEED):
    pw, fw = bootstrap_values(by_family, B, seed)
    return {"B": B, "seed": seed, "n_families": len(by_family), "pair_weighted": percentile_interval(pw), "family_weighted": percentile_interval(fw)}


def pair_rows(depth_pairs, U, genes):
    """The per-pair table behind a report row: names, distance class, identity, aligned fraction, the number of classes of each gene,
    S (1 iff U(g) and U(h) intersect) and the 'path' flag (either gene touches two or more classes)."""
    rows = []
    for d in depth_pairs:
        g, h = d["g"], d["h"]
        rows.append({"family_id": d["family_id"], "a_gene": d["a_gene"], "b_gene": d["b_gene"], "contig_a": genes[g].contig, "contig_b": genes[h].contig,
                     "dist_class": d["dist_class"], "identity": d["identity"], "cov_short": d["cov_short"], "n_classes_a": len(U[g]), "n_classes_b": len(U[h]),
                     "S": 1 if (U[g] & U[h]) else 0, "path": 1 if (len(U[g]) >= 2 or len(U[h]) >= 2) else 0})
    return rows


def report_table(depth_pairs, control, U, tau):
    """One row of the report at tau from depth-matched pairs (dicts with family_id, g, h, cov_short) and the control pool."""
    r = s_rate([(d["family_id"], d["g"], d["h"]) for d in depth_pairs], U)
    boot = bootstrap_rates(r["by_family"]) if r["n_families"] else None
    return {"tau": tau, "n_pairs": r["n_pairs"], "n_families": r["n_families"], "shared": r["shared"], "pair_weighted": r["pair_weighted"],
            "family_weighted": r["family_weighted"], "a_viol": r["a_viol"], "wilson": r["wilson"], "bootstrap": boot,
            "b_viol": list(b_viol(control, U)), "b_viol_by_class": {c: list(v) for c, v in b_viol_by_class(control, U).items()},
            "underpowered": underpowered(r["n_pairs"], r["n_families"]),
            "short_aligned_pairs": sum(1 for d in depth_pairs if d.get("cov_short") is not None and d["cov_short"] < 0.5),
            "path_pairs": sum(1 for d in depth_pairs if len(U[d["g"]]) >= 2 or len(U[d["h"]]) >= 2)}


# ---------------------------------------------------------------------------------------------------------------------------
# the positive control, Gate 0 fingerprints
# ---------------------------------------------------------------------------------------------------------------------------
def lrpap_copies(lines, genes):
    """The 8 full-length LRPAP1 copies: the rows of the gorilla self-lift (liftoff_loci.tsv) whose source_id is gene-LRPAP1, each joined to
    the annotated record of its landing site (the in_place row is LRPAP1 itself; the others name the record in `note`, 'overlaps gene-X (contig)')."""
    by = collections.defaultdict(list)
    for i, g in enumerate(genes):
        by[(g.contig, g.name)].append(i)
    out = []
    for r in csv.DictReader(lines, delimiter="\t"):
        if r["source_id"] != "gene-LRPAP1":
            continue
        if r["cls"] == "in_place":
            name = r["name"]
        else:
            m = re.match(r"overlaps gene-(\S+) \(", r["note"])
            if not m:
                raise ValueError(f"cannot read the overlapped record from note {r['note']!r}")
            name = m.group(1)
        cand = by.get((r["contig"], name), [])
        if len(cand) != 1:
            raise ValueError(f"LRPAP1 copy {name} on {r['contig']}: {len(cand)} gene records")
        out.append({"name": name, "gene": cand[0], "cls": r["cls"], "contig": r["contig"], "sequence_id": float(r["sequence_id"])})
    if len(out) != 8:
        raise ValueError(f"{len(out)} gene-LRPAP1 rows in the self-lift table, expected 8")
    return out


def positive_control(copies, U):
    """Each copy touches at least one class, and the graph whose vertices are the copies and whose edges join copies with intersecting U is
    connected. `passed` is False otherwise (the arm is then INVALID)."""
    parent = list(range(len(copies)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(copies)):
        for j in range(i + 1, len(copies)):
            if U.get(copies[i]["gene"], frozenset()) & U.get(copies[j]["gene"], frozenset()):
                a, b = find(i), find(j)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    comps = collections.defaultdict(list)
    for i, c in enumerate(copies):
        comps[find(i)].append(c["name"])
    without = [c["name"] for c in copies if not U.get(c["gene"])]
    ids = [c["sequence_id"] for c in copies if c["cls"] != "in_place"]
    connected = len(comps) == 1
    return {"n_copies": len(copies), "without_class": without, "connected": connected, "passed": (not without) and connected,
            "components": sorted(sorted(v) for v in comps.values()), "classes": {c["name"]: sorted(U.get(c["gene"], frozenset())) for c in copies},
            "identity_range": (min(ids), max(ids)) if ids else None}


def fingerprint(path):
    """size, newline count and sha256 of a file (read in 1 MiB chunks)."""
    h, size, lines = hashlib.sha256(), 0, 0
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
            size += len(b)
            lines += b.count(b"\n")
    return {"path": path, "size": size, "lines": lines, "sha256": h.hexdigest()}


def gate0_rows(paths, expected_prefix):
    """Gate 0 rows: the fingerprint of every input and whether its sha256 starts with the registered prefix (None: nothing registered)."""
    rows = []
    for name in sorted(paths):
        f = fingerprint(paths[name])
        pre = expected_prefix.get(name)
        f["name"], f["expected_prefix"] = name, pre
        f["ok"] = True if pre is None else f["sha256"].startswith(pre)
        rows.append(f)
    return rows


# ---------------------------------------------------------------------------------------------------------------------------
# the counts (Gate 6 / section 5, input-only) and the report (released only after an independent code review)
# ---------------------------------------------------------------------------------------------------------------------------
WINDATA = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
DEFAULT_PATHS = {
    "sedef": f"{WINDATA}/GGO_sedef_final.bed",
    "genes": "/mnt/linuxdisk/tmp/rustle_figures/families_gw/species/gorilla/genes.tsv",
    "truth": "/mnt/linuxdisk/tmp/fewcopy_gorilla_2026-10-07/gorilla_truth.tsv",
    "pairs": "/mnt/linuxdisk/tmp/fewcopy_gorilla_2026-10-07/pairs.tsv",
    "liftoff": "/mnt/linuxdisk/tmp/rustle_figures/liftoff/gorilla/liftoff_loci.tsv",
}
# registered in section 13 of the prereg (sha256, first 16 hex)
EXPECTED_PREFIX = {"sedef": "41aa1c9c53a18d68", "genes": "6f9bfbd776e8ecc8", "truth": "9e6cb80ad518f528", "pairs": "fb12608ae431a105", "liftoff": "af503bf2713caf54"}
DEFAULT_OUT = "/mnt/linuxdisk/tmp/hier_run/r2"
RELEASE_TOKEN = "RELEASE R2 REPORT AFTER INDEPENDENT CODE REVIEW"


class ReleaseRefused(PermissionError):
    """The report was asked for without the release file that is written only after the independent code review."""


def require_release(release_file):
    """Raise ReleaseRefused unless the first line of the release file is the registered token. Called by cmd_report AND by compute_report, so that
    importing the module does not bypass it."""
    ok = False
    if release_file and os.path.exists(release_file):
        with open(release_file) as fh:
            ok = fh.readline().strip() == RELEASE_TOKEN
    if not ok:
        raise ReleaseRefused("report: refused. The S-rate, A_viol and B_viol tables are computed only after the independent code review; "
                             "a release file with the registered token is required (--release-file).")


def default_registry():
    """The registered sha256 prefixes of section 13 (the five inputs) and the restored atoms script."""
    return dict(EXPECTED_PREFIX, **{"dna_sd_atoms.py": ATOMS_SHA256[:16]})


def read_registry(path):
    """A registry file: one `name<TAB>sha256 prefix` per line. The tests hand it the hashes of their made-up inputs; the real run uses default_registry()."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                out[f[0]] = f[1]
    return out


REQUIRED_INPUTS = ("sedef", "genes", "truth", "pairs", "liftoff", "dna_sd_atoms.py")


def validity(inputs, control, control_error=None):
    """(valid, reasons): the report is INVALID when a required input has no fingerprint row, no registered sha256 prefix or a mismatching one, or when the
    positive control failed, was not run, or could not be computed (control_error)."""
    by = {r["name"]: r for r in inputs}
    reasons = []
    for n in REQUIRED_INPUTS:
        r = by.get(n)
        if r is None:
            reasons.append(f"input {n}: no fingerprint row")
        elif r.get("expected_prefix") is None:
            reasons.append(f"input {n}: no registered sha256 prefix")
        elif not r["ok"]:
            reasons.append(f"input {n}: sha256 {r['sha256'][:16]} does not start with the registered {r['expected_prefix']}")
    if control_error is not None:
        reasons.append(f"positive control could not be computed ({control_error})")
    elif control is None:
        reasons.append("positive control not run")
    elif not control["passed"]:
        reasons.append(f"positive control failed (copies without a class: {control['without_class']}; copy graph connected: {control['connected']})")
    return (not reasons, reasons)


def environment_record():
    """The interpreter and the libraries a Gate 0 record names (python3, numpy, scipy, pysam; None when not installed)."""
    import importlib.metadata as md
    mods = {}
    for name in ("numpy", "scipy", "pysam"):
        try:
            mods[name] = md.version(name)
        except md.PackageNotFoundError:
            mods[name] = None
    return {"python_version": sys.version.split()[0], "executable": sys.executable, "platform": platform.platform(), "modules": mods}


REGISTERED_DEPTH_MATCHED = {"0.90": (38, 33), "0.95": (15, 15), "0.98": (5, 5)}   # Gate 6 tripwires of section 5 (pairs, families in the verdict set)
B_VIOL_CAVEAT = ("B_viol is a ceiling-type quantity: the control pool is not identity-matched and mostly cross-contig, so it bounds, and does not "
                 "estimate, how often genes of different families share a class.")


def _f(x):
    return "NA" if x is None else f"{x:.3f}"


def format_report(rep):
    """The printed report: one block per tau, every number of report.json that the review asked to see, and INVALID with its reasons when a gate failed."""
    if not rep["valid"]:
        return ["INVALID: the R2 report is withheld and no rate is printed."] + [f"  {r}" for r in rep["reasons"]]
    L = []
    for k in sorted(rep["tau"]):
        t = rep["tau"][k]
        reg = REGISTERED_DEPTH_MATCHED.get(k)
        tag = "" if reg is None else f" (registered {reg[0]}/{reg[1]}: {'match' if (t['n_pairs'], t['n_families']) == reg else 'DIFFERS'})"
        L.append(f"tau {k}: {t['n_pairs']} depth-matched pairs in {t['n_families']} families{tag}; UNDERPOWERED {t['underpowered']}")
        w = t["wilson"]
        L.append(f"  S-rate pair-weighted {_f(t['pair_weighted'])} (Wilson {_f(w[0])}-{_f(w[1])}), family-weighted {_f(t['family_weighted'])}; A_viol {_f(t['a_viol'])}")
        b = t["bootstrap"]
        if b:
            L.append(f"  bootstrap ({b['B']} family resamples, seed {b['seed']}): pair-weighted {_f(b['pair_weighted'][0])}-{_f(b['pair_weighted'][1])}, "
                     f"family-weighted {_f(b['family_weighted'][0])}-{_f(b['family_weighted'][1])}")
        else:
            L.append("  bootstrap: NA (no families)")
        L.append(f"  path pairs (a gene touches two or more classes): {t['path_pairs']} of {t['n_pairs']}; "
                 f"pairs aligned < 0.5 of the shorter transcript: {t['short_aligned_pairs']}")
        bc = ", ".join(f"{c} {v[0]}/{v[1]}" for c, v in t["b_viol_by_class"].items())
        L.append(f"  B_viol: {t['b_viol'][0]} of {t['b_viol'][1]} eligible different-family pairs share a class (by distance class: {bc or 'none'})")
        L.append(f"  {B_VIOL_CAVEAT}")
    return L


def tau_key(tau):
    return f"{tau:.2f}"


def class_counts(items):
    c = collections.Counter(x[-1] for x in items)
    return {k: c[k] for k in sorted(c)}


def universe_for_genes(info, genes, gene_indices):
    """U(g) for the given genes at one tau, from the atoms that prepare_tau left in info['prefix']."""
    atoms = read_nodes(info["prefix"])
    edges = read_edges(info["prefix"])
    class_of, _used = build_classes(len(atoms), edges)
    index = AtomIndex(atoms)
    return {g: gene_U(touch_classes(genes[g].blocks, genes[g].contig, index, class_of)) for g in sorted(gene_indices)}


def universe_for_tau(info, genes, layer):
    """U(g) for the labelled genes at one tau."""
    return universe_for_genes(info, genes, layer.gene_family)


def counts_for_tau(info, raw_rows, genes, layer, pairs, verdict, U):
    sfp_v = same_family_pairs(layer, genes, U, verdict)
    sfp_a = same_family_pairs(layer, genes, U, None)
    dm_v = depth_matched(pairs, layer, genes, U, info["tau"], verdict)
    dm_a = depth_matched(pairs, layer, genes, U, info["tau"], None)
    fam_v = len({d["family_id"] for d in dm_v})
    return {
        "rows_retained_raw": len(retain(raw_rows, info["tau"])), "rows_retained": info["rows_retained"], "contigs": info["contigs"], "atoms": info["atoms"],
        "edges_total": info["edges_total"], "edges_used": info["edges_used"], "classes": info["classes"], "largest_class_atoms": info["largest_class_atoms"],
        "singleton_classes": info["singleton_classes"],
        "labelled_genes": len(U), "eligible_labelled_genes": sum(1 for g in U if U[g]), "path_genes": sum(1 for g in U if len(U[g]) >= 2),
        "same_family_pairs": {"verdict": class_counts(sfp_v), "all": class_counts(sfp_a),
                              "partial_verdict": sum(1 for f, *_x in sfp_v if layer.status[f] == "partial"),
                              "whole_verdict": sum(1 for f, *_x in sfp_v if layer.status[f] == "whole")},
        "control_pairs": {"verdict": class_counts(control_pairs(layer, genes, U, verdict)), "all": class_counts(control_pairs(layer, genes, U, None))},
        "pairs_tsv_eligible_verdict": len(eligible_rows(pairs, layer, genes, U, verdict)),
        "pairs_tsv_eligible_all": len(eligible_rows(pairs, layer, genes, U, None)),
        "depth_matched": {"verdict": {"pairs": len(dm_v), "families": fam_v}, "all": {"pairs": len(dm_a), "families": len({d["family_id"] for d in dm_a})}},
        "underpowered_verdict": underpowered(len(dm_v), fam_v),
    }


def load_inputs(a):
    with open(a.sedef) as fh:
        raw, sedef_stats = read_sedef(fh)
    deduped, removed = dedupe(raw)
    with open(a.genes) as fh:
        genes, index, dups = read_genes(fh)
    with open(a.truth) as fh:
        layer = read_layer(fh, index)
    with open(a.pairs) as fh:
        pairs = read_pairs(fh)
    return raw, sedef_stats, deduped, removed, genes, dups, layer, pairs


def cmd_prepare(a):
    raw, stats, deduped, removed, *_ = load_inputs(a)
    out = os.path.join(a.outdir, "atoms")
    summary = {}
    for tau in TAUS:
        info = prepare_tau(deduped, tau, out, reuse=not a.rebuild)
        summary[tau_key(tau)] = {k: v for k, v in info.items() if k not in ("prefix", "log")}
        print(f"tau {tau_key(tau)}: rows {info['rows_retained']}, atoms {info['atoms']}, edges {info['edges_total']} (used {info['edges_used']}), classes {info['classes']}"
              + (" (reused)" if info["reused"] else ""))
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "prepare.json"), "w") as fh:
        json.dump({"sedef": dict(stats, duplicates_removed=removed, unique_tuples=len(deduped)), "tau": summary}, fh, indent=1, sort_keys=True)
    return 0


# the recon values of section 5 (tau 0.90), taken on un-deduplicated rows with two-sided coverage > 0.5; Gate 6 prints the counts under the registered rules
RECON_TAU090 = {"atoms": 40631, "classes": 12483, "edges": 372793}


def format_counts(c):
    L = ["== R2 counts (input-only; denominators first) =="]
    s = c["sedef"]
    L.append(f"SEDEF: data rows {s['data_rows']}, header lines {s['header_lines']}, rows on {MITO_CONTIG} dropped {s['mito_rows_dropped']}, rows kept {s['rows_after_mito']}, "
             f"identical tuples removed {s['duplicates_removed']}, unique tuples {s['unique_tuples']}")
    y = c["layer"]
    L.append(f"layer: 1to1 genes {y['genes_1to1']} in {y['families_with_1to1_genes']} families of {y['families_all']} (status {y['status_counts']}); "
             f"families with >= 2 projected genes {y['families_with_two_or_more_genes']}; verdict-set families {y['families_in_verdict_set']} "
             f"({y['genes_in_verdict_set']} genes); gene table rows {y['gene_table_rows']}, duplicated join keys {y['gene_table_duplicate_keys']}")
    p = c["pairs_tsv"]
    L.append(f"pairs.tsv: rows {p['rows']} from {p['families']} families, {p['with_identity']} with an identity; identity >= 0.90: {p['identity_ge_0.90']}, "
             f"of which aligned fraction (cov_short) < 0.5: {p['short_aligned_lt_0.5_among_ge_0.90']}")
    for k in sorted(c["tau"]):
        t = c["tau"][k]
        L.append(f"tau {k}: rows retained {t['rows_retained_raw']} (before dedupe), {t['rows_retained']} (after); atoms {t['atoms']}, atom edges {t['edges_total']} "
                 f"({t['edges_used']} with two-sided coverage >= 0.5), classes {t['classes']} (largest {t['largest_class_atoms']} atoms, {t['singleton_classes']} singletons)")
        if k == tau_key(TAU_PRIMARY):
            L.append(f"   registered recon values at tau 0.90 (un-deduplicated rows, two-sided coverage > 0.5): {RECON_TAU090['atoms']} atoms, {RECON_TAU090['classes']} classes, "
                     f"{RECON_TAU090['edges']} atom edges; this run (identical tuples collapsed, coverage >= 0.5, as registered): {t['atoms']} atoms "
                     f"({t['atoms'] - RECON_TAU090['atoms']:+d}), {t['classes']} classes ({t['classes'] - RECON_TAU090['classes']:+d}), {t['edges_total']} atom edges "
                     f"({t['edges_total'] - RECON_TAU090['edges']:+d}); the differences are the two registered rule changes")
        L.append(f"   labelled genes {t['labelled_genes']}, eligible (U non-empty) {t['eligible_labelled_genes']}, flagged 'path' (|U| >= 2) {t['path_genes']}")
        sf = t["same_family_pairs"]
        L.append(f"   eligible same-family pairs, verdict set: {sf['verdict']} (whole families {sf['whole_verdict']}, partial {sf['partial_verdict']}); all families: {sf['all']}")
        L.append(f"   eligible different-family pairs (control pool), verdict set: {t['control_pairs']['verdict']}; all: {t['control_pairs']['all']}")
        d = t["depth_matched"]
        L.append(f"   pairs.tsv rows that are eligible: verdict set {t['pairs_tsv_eligible_verdict']}, all {t['pairs_tsv_eligible_all']}; depth-matched at identity >= {k}: "
                 f"verdict set {d['verdict']['pairs']} pairs in {d['verdict']['families']} families (UNDERPOWERED: {t['underpowered_verdict']}); all {d['all']['pairs']} in {d['all']['families']}")
    pc = c.get("positive_control")
    if pc is None:
        L.append("positive control: not run (no self-lift table given)")
    else:
        L.append(f"positive control (tau 0.90): {pc['n_copies']} copies, without a class {pc['without_class']}, connected {pc['connected']}, passed {pc['passed']}, "
                 f"components {pc['components']}, identity range of the non-reference copies {pc['identity_range']}")
    return "\n".join(L)


def compute_counts(a):
    raw, stats, deduped, removed, genes, dups, layer, pairs = load_inputs(a)
    verdict = verdict_families(layer, genes)
    out = os.path.join(a.outdir, "atoms")
    c = {"sedef": dict(stats, rows_after_mito=len(raw), duplicates_removed=removed, unique_tuples=len(deduped))}
    fam_sizes = [len(m) for m in layer.families.values()]
    c["layer"] = {"genes_1to1": len(layer.gene_family), "families_all": layer.n_families_all, "families_with_1to1_genes": len(layer.families),
                  "families_with_two_or_more_genes": sum(1 for n in fam_sizes if n >= 2), "families_in_verdict_set": len(verdict),
                  "genes_in_verdict_set": sum(len(layer.families[f]) for f in verdict), "status_counts": dict(sorted(collections.Counter(layer.status.values()).items())),
                  "gene_table_rows": len(genes), "gene_table_duplicate_keys": len(dups)}
    num = [r for r in pairs if r["identity"] is not None]
    ge90 = [r for r in num if r["identity"] >= 0.90]
    c["pairs_tsv"] = {"rows": len(pairs), "families": len({r["family_id"] for r in pairs}), "with_identity": len(num), "identity_ge_0.90": len(ge90),
                      "short_aligned_lt_0.5_among_ge_0.90": sum(1 for r in ge90 if r["cov_short"] is not None and r["cov_short"] < 0.5)}
    c["tau"] = {}
    primary_info = None
    for tau in TAUS:
        info = prepare_tau(deduped, tau, out, reuse=not a.rebuild)
        U = universe_for_tau(info, genes, layer)
        c["tau"][tau_key(tau)] = counts_for_tau(info, raw, genes, layer, pairs, verdict, U)
        if tau == TAU_PRIMARY:
            primary_info = info
    c["positive_control"] = None
    if a.liftoff:
        with open(a.liftoff) as fh:
            copies = lrpap_copies(fh, genes)
        Uc = universe_for_genes(primary_info, genes, [x["gene"] for x in copies])
        pc = positive_control(copies, Uc)
        pc["identity_range"] = list(pc["identity_range"]) if pc["identity_range"] else None
        c["positive_control"] = pc
    return c


def cmd_counts(a):
    c = compute_counts(a)
    c["inputs"] = [fingerprint(p) for p in (a.sedef, a.genes, a.truth, a.pairs)]
    os.makedirs(a.outdir, exist_ok=True)
    print(format_counts(c))
    with open(os.path.join(a.outdir, "counts.json"), "w") as fh:
        json.dump(c, fh, indent=1, sort_keys=True)
    print("wrote counts.json")
    return 0


def _registry_for(a):
    reg = read_registry(a.registry) if getattr(a, "registry", None) else default_registry()
    need = list(REQUIRED_INPUTS)
    missing = [n for n in need if n not in reg]
    if missing:
        raise ValueError(f"registry lacks {missing}: every input of the run must have a registered sha256 prefix")
    short = sorted(n for n in need if len(reg[n]) < 16)
    if short:
        raise ValueError(f"registry prefix too short for {short}: a registered sha256 prefix has at least 16 hex characters")
    return reg


def cmd_gate0(a):
    paths = {"sedef": a.sedef, "genes": a.genes, "truth": a.truth, "pairs": a.pairs, "liftoff": a.liftoff or DEFAULT_PATHS["liftoff"],
             "dna_sd_atoms.py": ATOMS_SCRIPT}
    try:
        reg = _registry_for(a)
    except ValueError as e:
        print(f"gate0: {e}", file=sys.stderr)
        return 2
    rows = gate0_rows(paths, reg)
    env = environment_record()
    os.makedirs(a.outdir, exist_ok=True)
    bad = 0
    with open(os.path.join(a.outdir, "gate0_r2.tsv"), "w") as fh:
        fh.write("name\tpath\tsize\tlines\tsha256\texpected_prefix\tok\n")
        for r in rows:
            fh.write(f"{r['name']}\t{r['path']}\t{r['size']}\t{r['lines']}\t{r['sha256']}\t{r['expected_prefix']}\t{r['ok']}\n")
            print(f"{r['name']:16s} {r['size']:>12d} bytes {r['lines']:>9d} lines sha256 {r['sha256'][:16]} registered {r['expected_prefix']} {'OK' if r['ok'] else 'MISMATCH'}")
            bad += 0 if r["ok"] else 1
    with open(os.path.join(a.outdir, "gate0_r2_env.tsv"), "w") as fh:
        fh.write(f"python_version\t{env['python_version']}\nexecutable\t{env['executable']}\nplatform\t{env['platform']}\n")
        for m, v in env["modules"].items():
            fh.write(f"{m}\t{v if v is not None else 'not installed'}\n")
    print(f"interpreter: python {env['python_version']} at {env['executable']}; " + ", ".join(f"{m} {v or 'not installed'}" for m, v in env["modules"].items()))
    print("Gate 0: " + ("PASS" if not bad else f"FAIL ({bad} mismatches)"))
    return 0 if not bad else 1


def compute_report(a):
    """The release-gated report. Refuses without the token (also when called directly), validates the registry, then checks the Gate 0 fingerprints of its
    own inputs and the positive control BEFORE any rate is computed; a failed gate returns {'valid': False, 'reasons': [...]} with no table."""
    require_release(getattr(a, "release_file", None))
    registry = _registry_for(a)                       # before anything is built: a bad registry leaves no trace
    raw, stats, deduped, removed, genes, dups, layer, pairs = load_inputs(a)
    verdict = verdict_families(layer, genes)
    out = os.path.join(a.outdir, "atoms")
    infos = {tau: prepare_tau(deduped, tau, out, reuse=not a.rebuild) for tau in TAUS}
    paths = {"sedef": a.sedef, "genes": a.genes, "truth": a.truth, "pairs": a.pairs, "liftoff": a.liftoff, "dna_sd_atoms.py": ATOMS_SCRIPT}
    inputs = gate0_rows(paths, registry)
    pc, pc_error = None, None
    try:
        with open(a.liftoff) as fh:
            copies = lrpap_copies(fh, genes)
        pc = positive_control(copies, universe_for_genes(infos[TAU_PRIMARY], genes, [x["gene"] for x in copies]))
        pc["identity_range"] = list(pc["identity_range"]) if pc["identity_range"] else None
    except (ValueError, KeyError, OSError, csv.Error) as e:
        pc_error = f"{type(e).__name__}: {e}"
    valid, reasons = validity(inputs, pc, pc_error)
    with open(os.path.abspath(__file__), "rb") as fh:
        scorer_sha256 = hashlib.sha256(fh.read()).hexdigest()
    rep = {"release": RELEASE_TOKEN, "valid": valid,
           "gates": {"inputs": inputs, "positive_control": pc, "positive_control_error": pc_error, "environment": environment_record(),
                     "registry": getattr(a, "registry", None) or "default (the registered constants of section 13)", "scorer_sha256": scorer_sha256}}
    if not valid:
        rep["reasons"] = reasons
        return rep
    rep["tau"], rep["pair_rows"] = {}, {}
    for tau in TAUS:
        U = universe_for_tau(infos[tau], genes, layer)
        dm = depth_matched(pairs, layer, genes, U, tau, verdict)
        ctl = control_pairs(layer, genes, U, verdict)
        t = report_table(dm, ctl, U, tau)
        t["wilson"] = list(t["wilson"])
        if t["bootstrap"]:
            t["bootstrap"]["pair_weighted"] = list(t["bootstrap"]["pair_weighted"])
            t["bootstrap"]["family_weighted"] = list(t["bootstrap"]["family_weighted"])
        rep["tau"][tau_key(tau)] = t
        rep["pair_rows"][tau_key(tau)] = pair_rows(dm, U, genes)
    return rep


def cmd_report(a):
    try:
        require_release(a.release_file)
    except ReleaseRefused as e:
        print(e, file=sys.stderr)
        return 2
    if not a.liftoff or not os.path.exists(a.liftoff):
        print(f"report: refused. The positive control needs the gorilla self-lift table (--liftoff); {a.liftoff!r} does not exist.", file=sys.stderr)
        return 2
    try:
        rep = compute_report(a)
    except ValueError as e:
        if str(e).startswith("registry "):
            print(f"report: refused. {e}", file=sys.stderr)
            return 2
        raise
    os.makedirs(a.outdir, exist_ok=True)
    pair_tables = rep.pop("pair_rows", {})
    if not rep["valid"]:                              # an INVALID run leaves no table of an earlier run next to it (exact names, our own outputs)
        for tau in TAUS:
            stale = os.path.join(a.outdir, f"report_pairs_{tau_key(tau)}.tsv")
            if os.path.exists(stale):
                os.remove(stale)
    with open(os.path.join(a.outdir, "report.json"), "w") as fh:
        json.dump(rep, fh, indent=1, sort_keys=True)
    cols = ["family_id", "a_gene", "b_gene", "contig_a", "contig_b", "dist_class", "identity", "cov_short", "n_classes_a", "n_classes_b", "S", "path"]
    for k in sorted(pair_tables):
        with open(os.path.join(a.outdir, f"report_pairs_{k}.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in pair_tables[k]:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    env = rep["gates"]["environment"]
    pcr = rep["gates"]["positive_control"]
    print(f"interpreter: python {env['python_version']}; positive control passed {pcr['passed'] if pcr else None}; "
          f"inputs matching their registered sha256 prefix: {sum(1 for r in rep['gates']['inputs'] if r['ok'])} of {len(rep['gates']['inputs'])}")
    print("\n".join(format_report(rep)))
    return 0 if rep["valid"] else 3


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("prepare", cmd_prepare), ("counts", cmd_counts), ("gate0", cmd_gate0), ("report", cmd_report)):
        sp = sub.add_parser(name)
        sp.add_argument("--sedef", default=DEFAULT_PATHS["sedef"])
        sp.add_argument("--genes", default=DEFAULT_PATHS["genes"])
        sp.add_argument("--truth", default=DEFAULT_PATHS["truth"])
        sp.add_argument("--pairs", default=DEFAULT_PATHS["pairs"])
        sp.add_argument("--liftoff", default=DEFAULT_PATHS["liftoff"] if name == "report" else None,
                        help="the gorilla self-lift table (liftoff_loci.tsv); enables the positive control in `counts`; `report` always runs the control and defaults to the registered table")
        sp.add_argument("--outdir", default=DEFAULT_OUT)
        sp.add_argument("--rebuild", action="store_true", help="recompute the atoms even when the retained table is unchanged")
        if name in ("gate0", "report"):
            sp.add_argument("--registry", default=None, help="TSV `name<TAB>sha256 prefix` replacing the registered prefixes (the tests use it for made-up inputs)")
        if name == "report":
            sp.add_argument("--release-file", default=None, help="file whose first line is the release token (written only after the independent code review)")
        sp.set_defaults(func=fn)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
