#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 RNA-only
EDGE-ORACLE loaders + the allele-demote decision (`vg_family/edge_oracle.rs`),
the migration increment (RETIREMENT_AND_MIGRATION.md Tier-1 #4) that ports what
the family-definition driver consumes from `bench/rna_only_edge_oracle.py`.

It reuses the SHIPPED Python loaders verbatim -- `rna_only_edge_oracle.{load_pair_core,
load_universal_aln,load_allele}` + `demote_gene`'s exact expression -- run over
REDUCED, self-contained inputs, and emits to
`src/rustle/vg_family/testdata/edge_oracle_fixture.json`:

  (a) pair_core : a reduced EDGES_TSV + reduced gene_of_dn covering
      - multi-DN-edge gene-pairs that FORCE the MAX aggregation (incl. pairs where
        the max is NOT the first edge -> proves order-invariance),
      - single-edge pairs,
      - an edge whose endpoint is dropped from gene_of_dn (missing -> SKIP),
      - an edge whose two DN endpoints project to the SAME gene (ga==gb -> SKIP),
      with the Python golden `load_pair_core(reduced_gene_of_dn)` on that reduced TSV.

  (b) universal_aln : a reduced SHAREDLEN_UNIV_TSV subset (numeric + zero aln_frac)
      PLUS synthetic empty-aln_frac rows (-> None), with the Python golden.

  (c) allele : a reduced A1_TSV subset (incl. the known allele-as-copy genes), with
      the Python golden `load_allele`.

  (d) demote : a diverse set of demote_gene inputs (real demoted/not + synthetic
      boundary cases at each threshold) -> the Python demote decision.

Every golden is produced by the REAL shipped function pointed at the reduced input
(via a temporary module-constant swap), so the fixture is byte-faithful to ship.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_edge_oracle_fixture.py
Writes: src/rustle/vg_family/testdata/edge_oracle_fixture.json
Deterministic (sorted selection; PYTHONHASHSEED=0).
"""
import os
import sys

# Match the shipped module's own PYTHONHASHSEED=0 re-exec (importing RO does it
# anyway; we set it here so our sorting/selection is deterministic regardless).
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
import tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rna_only_edge_oracle as RO   # noqa: E402  (loaders + demote_gene spec)

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "edge_oracle_fixture.json")


def pair_key(a, b):
    """Sorted [gene_a, gene_b] = Rust `gene_pair_key` / Python `frozenset`."""
    return [a, b] if a <= b else [b, a]


def run_with_const(const_name, tmp_text, fn, *args):
    """Write tmp_text to a temp file, point RO.<const_name> at it, run the shipped
    loader fn(*args), restore, and return the result. Faithful: it runs the REAL
    shipped function on the reduced input."""
    old = getattr(RO, const_name)
    fd, path = tempfile.mkstemp(suffix=".tsv")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(tmp_text)
        setattr(RO, const_name, path)
        return fn(*args)
    finally:
        setattr(RO, const_name, old)
        os.unlink(path)


# ---------------------------------------------------------------- pair_core
def build_pair_core(gene_of_dn):
    # Read every edge row (a, b, core_recip str) in file order.
    rows = []
    with open(RO.EDGES_TSV) as fh:
        next(fh)
        for ln in fh:
            a, b, c = ln.rstrip("\n").split("\t")
            rows.append((a, b, c))

    # Group by gene-pair (only the non-skip edges), preserving file order.
    pair_rows = defaultdict(list)   # frozenset(gene-pair) -> [(a,b,c_str), ...]
    same_gene_rows = []             # rows where gene_of_dn[a] == gene_of_dn[b]
    for a, b, c in rows:
        ga, gb = gene_of_dn.get(a), gene_of_dn.get(b)
        if not (ga and gb):
            continue
        if ga == gb:
            same_gene_rows.append((a, b, c))
            continue
        pair_rows[frozenset((ga, gb))].append((a, b, c))

    # Classify pairs (sorted iteration for determinism).
    def keyf(k):
        return tuple(sorted(k))
    multi_not_first, multi_is_first, singles = [], [], []
    for k in sorted(pair_rows, key=keyf):
        v = pair_rows[k]
        cvs = [float(c) for _, _, c in v]
        if len(v) == 1:
            singles.append(k)
        elif cvs[0] != max(cvs):
            multi_not_first.append(k)      # MAX is a NON-first edge (order matters)
        else:
            multi_is_first.append(k)

    # Select a representative, diverse subset.
    sel_pairs = (multi_not_first[:3] + multi_is_first[:2] + singles[:3])
    sel_pairs = list(dict.fromkeys(sel_pairs))     # de-dup, keep order

    # A dedicated single-edge pair for the MISSING-endpoint skip (drop one endpoint
    # so the whole pair vanishes) -- choose one NOT already selected and whose DNs
    # are not reused by the selection.
    used_dn = {a for k in sel_pairs for (a, b, c) in pair_rows[k]} | \
              {b for k in sel_pairs for (a, b, c) in pair_rows[k]}
    missing_pair = None
    for k in singles:
        if k in sel_pairs:
            continue
        (a, b, c) = pair_rows[k][0]
        if a not in used_dn and b not in used_dn:
            missing_pair = k
            break
    assert missing_pair is not None, "no clean single-edge pair for the missing case"
    (ma, mb, mc) = pair_rows[missing_pair][0]
    drop_dn = ma                                   # drop endpoint a -> pair skipped

    # Pick one same-gene skip row whose DNs are not reused.
    same_row = None
    for (a, b, c) in same_gene_rows:
        if a not in used_dn and b not in used_dn and a != drop_dn and b != drop_dn:
            same_row = (a, b, c)
            break
    assert same_row is not None, "no clean same-gene (ga==gb) skip row"

    # Assemble reduced edge rows: selected pairs (file order) + missing row + same row.
    reduced_edge_rows = []
    for k in sel_pairs:
        reduced_edge_rows.extend(pair_rows[k])
    reduced_edge_rows.append((ma, mb, mc))
    reduced_edge_rows.append(same_row)

    # Reduced gene_of_dn = mapping for every endpoint in the reduced rows, MINUS the
    # deliberately-dropped endpoint (drop_dn).
    reduced_dns = set()
    for a, b, c in reduced_edge_rows:
        reduced_dns.add(a)
        reduced_dns.add(b)
    reduced_god = {dn: gene_of_dn[dn] for dn in sorted(reduced_dns)
                   if dn in gene_of_dn and dn != drop_dn}

    edges_tsv = "a\tb\tcore_recip\n" + \
        "".join(f"{a}\t{b}\t{c}\n" for a, b, c in reduced_edge_rows)

    # GOLDEN via the shipped load_pair_core on the reduced TSV + reduced gene_of_dn.
    golden = run_with_const("EDGES_TSV", edges_tsv, RO.load_pair_core, reduced_god)

    # Assertions that the skips actually fired (fixture demonstrates the paths).
    gk_missing = frozenset((gene_of_dn[ma], gene_of_dn[mb]))
    assert gk_missing not in golden, "missing-endpoint pair unexpectedly present"
    gk_same = frozenset((gene_of_dn[same_row[0]],))
    assert gk_same not in golden and \
        frozenset((gene_of_dn[same_row[0]], gene_of_dn[same_row[1]])) not in golden, \
        "same-gene pair unexpectedly present"

    golden_list = sorted(
        [pair_key(*sorted(k)) + [v] for k, v in golden.items()],
        key=lambda e: (e[0], e[1]))

    diversity = dict(
        n_pairs=len(golden),
        multi_edge_pairs=sum(1 for k in sel_pairs if len(pair_rows[k]) > 1),
        max_not_first=sum(1 for k in sel_pairs if k in set(multi_not_first)),
        skipped_missing=1,
        skipped_same_gene=1,
    )
    return dict(edges_tsv=edges_tsv, gene_of_dn=reduced_god,
                golden=golden_list, diversity=diversity)


# ---------------------------------------------------------------- universal_aln
UNIV_HDR = "gene_a\tgene_b\tdn_a\tdn_b\tcls\tin_ep\tcore\taln_len\taln_frac"


def build_universal_aln():
    with open(RO.SHAREDLEN_UNIV_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(header)}
        data = [ln.rstrip("\n").split("\t") for ln in fh]
    # deterministic selection: sort by (gene_a, gene_b); take a spread of zero +
    # nonzero aln_frac rows.
    data.sort(key=lambda f: (f[I["gene_a"]], f[I["gene_b"]]))
    nonzero = [f for f in data if f[I["aln_frac"]] not in ("", "0.0")]
    zero = [f for f in data if f[I["aln_frac"]] == "0.0"]
    sel = nonzero[:15] + zero[:15]

    def line_of(f):
        return "\t".join(f)

    # synthetic empty-aln_frac rows -> None (the real cache has none; the loader
    # handles "" -> None, so we exercise that path explicitly).
    syn_empty = [
        ["SYN_EMPTYB", "SYN_EMPTYA", "dn1", "dn2", "arb", "0", "0.0", "0", ""],
        ["SYN_EMPTYC", "SYN_EMPTYD", "dn3", "dn4", "arb", "0", "0.0", "0", ""],
    ]
    reduced_rows = sel + syn_empty
    tsv = UNIV_HDR + "\n" + "".join(line_of(f) + "\n" for f in reduced_rows)

    golden = run_with_const("SHAREDLEN_UNIV_TSV", tsv, RO.load_universal_aln)
    golden_list = sorted(
        [pair_key(*sorted(k)) + [(None if v is None else v)] for k, v in golden.items()],
        key=lambda e: (e[0], e[1]))
    return dict(tsv=tsv, golden=golden_list,
                n_empty=sum(1 for _, _, v in golden_list if v is None))


# ---------------------------------------------------------------- allele
def build_allele():
    with open(RO.A1_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(header)}
        data = [ln.rstrip("\n").split("\t") for ln in fh]
    genes = sorted({f[I["gene"]] for f in data})
    # ensure the named allele-as-copy genes are present, plus a deterministic spread.
    must = [g for g in ("DHRSX", "LOC115932259", "LOC129530050", "FOXO1") if g in genes]
    rest = [g for g in genes if g not in must][:12]
    sel_genes = list(dict.fromkeys(must + rest))
    sel_rows = [f for f in data if f[I["gene"]] in set(sel_genes)]
    tsv = "\t".join(header) + "\n" + "".join("\t".join(f) + "\n" for f in sel_rows)

    golden = run_with_const("A1_TSV", tsv, RO.load_allele)
    golden_json = {g: dict(balanced_frac=v["balanced_frac"], het_risk=v["het_risk"],
                           copy_like=v["copy_like"], med_maf=v["med_maf"])
                   for g, v in golden.items()}
    return dict(tsv=tsv, golden=golden_json), golden


# ---------------------------------------------------------------- demote
def demote_decision(allele, g):
    """EXACT replica of the nested `rna_only_edge_oracle.demote_gene` (lines
    1235-1237): `a is not None and a["balanced_frac"] >= 0.90 and
    a["copy_like"] <= 0.10`. Thresholds cited: DEMOTE_BAL_MIN=0.90 /
    DEMOTE_COPY_MAX=0.10 (family_rna_refine.py:181-182)."""
    a = allele.get(g)
    return a is not None and a["balanced_frac"] >= 0.90 and a["copy_like"] <= 0.10


def build_demote(real_allele):
    # allele dict = real genes + synthetic BOUNDARY genes at each threshold.
    def sig(bal, het, cop, maf):
        return dict(balanced_frac=bal, het_risk=het, copy_like=cop, med_maf=maf)

    allele = {g: dict(balanced_frac=v["balanced_frac"], het_risk=v["het_risk"],
                      copy_like=v["copy_like"], med_maf=v["med_maf"])
              for g, v in real_allele.items()}
    syn = {
        "SYN_BAL_EXACT": sig(0.90, 1.0, 0.10, 0.5),      # both boundaries inclusive -> TRUE
        "SYN_BAL_JUST_BELOW": sig(0.8999, 1.0, 0.05, 0.5),   # bal < 0.90 -> FALSE
        "SYN_COPY_JUST_ABOVE": sig(0.95, 1.0, 0.1001, 0.5),  # copy > 0.10 -> FALSE
        "SYN_COPY_EXACT": sig(0.95, 1.0, 0.10, 0.5),     # copy == 0.10 (inclusive) -> TRUE
        "SYN_BAL_HIGH_COPY_HIGH": sig(0.99, 1.0, 0.50, 0.5),  # copy high -> FALSE
        "SYN_BOTH_LOW": sig(0.50, 0.0, 0.50, 0.4),       # both low -> FALSE
    }
    allele.update(syn)

    # diverse case genes: real demoted, real not-demoted, all synthetic, and ABSENT.
    case_genes = []
    for g in ("DHRSX", "LOC115932259", "LOC129530050", "FOXO1"):
        if g in allele:
            case_genes.append(g)
    case_genes += list(syn.keys())
    case_genes.append("NOT_A_REAL_GENE")            # absent -> a is None -> FALSE
    case_genes = list(dict.fromkeys(case_genes))

    cases = [dict(gene=g, expect=bool(demote_decision(allele, g))) for g in case_genes]
    return dict(allele=allele, cases=cases)


# ---------------------------------------------------------------- main
def main():
    print("[load] gene_of_dn projection (shipped family_er_pr) ...", flush=True)
    gene_of_dn = RO.load_gene_of_dn()

    pair_core = build_pair_core(gene_of_dn)
    universal_aln = build_universal_aln()
    allele_fx, real_allele = build_allele()
    demote = build_demote(real_allele)

    fixture = dict(
        paths=dict(edges_tsv=os.path.basename(RO.EDGES_TSV),
                   sharedlen_univ_tsv=os.path.basename(RO.SHAREDLEN_UNIV_TSV),
                   a1_tsv=os.path.basename(RO.A1_TSV)),
        thresholds=dict(demote_bal_min=0.90, demote_copy_max=0.10),
        pair_core=pair_core,
        universal_aln=universal_aln,
        allele=allele_fx,
        demote=demote,
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)

    d = pair_core["diversity"]
    n_true = sum(1 for c in demote["cases"] if c["expect"])
    print(f"[fixture] wrote {os.path.normpath(OUT)}")
    print(f"          pair_core: {d['n_pairs']} pairs "
          f"(multi-edge={d['multi_edge_pairs']} max-not-first={d['max_not_first']} "
          f"skip-missing={d['skipped_missing']} skip-same-gene={d['skipped_same_gene']}); "
          f"reduced gene_of_dn={len(pair_core['gene_of_dn'])}")
    print(f"          universal_aln: {len(universal_aln['golden'])} rows "
          f"({universal_aln['n_empty']} None)")
    print(f"          allele: {len(allele_fx['golden'])} genes")
    print(f"          demote: {len(demote['cases'])} cases "
          f"({n_true} demoted / {len(demote['cases']) - n_true} not); "
          f"allele dict={len(demote['allele'])} genes")


if __name__ == "__main__":
    main()
