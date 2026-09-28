#!/usr/bin/env python3
"""Soto 2025 gene-family replication: one module, one subcommand per step.

WHAT THIS MEASURES (register T15 / row 1085 / row 858). The chain uses Soto's own famCN (Table S1C), Soto's
CAT v4 genes and Soto's gene universe. It measures CONCORDANCE with Soto et al. 2025 (Cell 188:5363), not an
independent replication. The famCN leg is circular: Soto built their families with famCN, so any score that
splits on famCN and is then evaluated against Soto is partly circular and must say so (register 858).

THE HEADLINE CHAIN (ledger §6ie-§6ip; inputs under winloci_data/soto_replication/, see docs/DATA.md):

    S=bench/soto/soto_replication.py
    python3 $S genesets --out-eligible g1793.tsv --out-full g2334.tsv
    python3 $S edges --sedef final_human_clean.bed --geneset g2334.tsv --cat-bed cat_v4.bed \\
        --extra-anchors bench/soto/acro_extra_anchors.tsv --out-shared shared.tsv        # 4,192 edges
    #   (or skip: bench/soto/shared_exons_2334_finalhuman.tsv is this step's frozen output)
    python3 $S cluster --shared shared.tsv --geneset g1793.tsv --full-geneset g2334.tsv \\
        --famcn bench/soto/soto_famCN_S1C.tsv --mad-statistic median --out rep_median.tsv
    python3 $S score --predicted rep_median.tsv                     # --truth defaults to S1C

Expected (§6ip, median / mean MAD): ARI 0.6959 / 0.6862, exact 241/491 (49.1%) / 264/491 (53.8%), pair
P/R/F1 0.841/0.595/0.697 / 0.906/0.554/0.687; bipartite MICRO 0.784/0.709 / 0.812/0.721, MACRO 0.731/0.718 /
0.770/0.739, undetected 99/491 / 88/491. `cluster` output is byte-identical to the frozen
`replicated_families_2334_{median,mean}_finalhuman.tsv`.

NATIVE v1.0 ARM (2026-09-28): `edges --native-v1 --sedef final_v1.bed` (the user's CHM13 v1.0 SEDEF with CIGARs,
151,759 rows, trailing header row stripped; identity lift, all 3,628 rows >= 0.98 used, 0 dropped) -> 4,384 edges
(frozen: shared_exons_2334_finalv1_native.tsv). Same cluster/score chain, median / mean MAD: ARI 0.6985 / 0.6894,
exact 241 / 264 (unchanged set), pair P/R/F1 0.837/0.601/0.700 / 0.900/0.560/0.691, bipartite MICRO 0.780/0.715 /
0.807/0.728, undetected 98 / 87. The gain is the acrocentrics the v2.0->v1.0 guard dropped (27 genes newly placed,
all with a truth-correct partner; ID_328 detected but fused with the CN-degenerate ID_321). Union with the v2.0
edges scores identically. Earlier v1.0 arms without SEDEF's own CIGAR: 0.6316 (linear), 0.6618 (mappy CIGAR).

CURATION ARM (2026-09-28, KEY=rule; figs/soto_cur_rule.md): `curate` recomputes Soto's SD98 gene set from SD98
regions + CAT v4 (their `-f 1` rule: UCSC v1.0 track -> 5,154 genes / 1,864 eligible = their 1,793 + the 71 they
removed by hand; `--sedef final_v1.bed` -> 5,364 / 1,999) and applies an explicit version of their manual
curation. Rule search (scratchpad/soto_cur/rule_search.py, 44 vocabulary atoms, every <= 3-atom and/or form,
selected on a seeded random half, scored on the other): the ONE-condition rule "drop a gene when >= 0.3 of its
exonic bp lies inside a same-strand SD98 partner that outranks it (protein_coding > transcribed > unprocessed
pseudogene; HGNC-named > readthrough > clone-named; more transcripts; longer)" is what every split picks
(36/50 seeds); held-out P/R 0.78/0.80 (50-seed mean; fixed a priori 0.80/0.82), and 2-3 atoms, depth-2/3 trees
and a 300-tree forest do NOT beat it out of sample (F1 .77/.76/.72/.79 vs .78). Plugged in:
    python3 $S curate --cat-bed cat_v4.bed --sd98-bed sd98_v1.bed --out-eligible ge.tsv --out-full gf.tsv
    python3 $S edges --native-v1 --sedef final_v1_clean.bed --geneset gf.tsv --cat-bed cat_v4.bed --out-shared s.tsv
    python3 $S cluster --shared s.tsv --geneset ge.tsv --full-geneset gf.tsv --famcn ... --out f.tsv; score
recomputed+uncurated 1,864: ARI 0.6946 / 0.6822 (median / mean), exact 234 / 257; recomputed+rule (72 dropped,
1,792): 0.6983 / 0.6886, exact 240 / 263 -- vs their hand-curated 1,793: 0.6985 / 0.6894, 241 / 264. Fully native
(`--sedef final_v1.bed` regions too, 1,921 eligible): 0.6981 / 0.6884. Off by default: no other subcommand calls
any of it (headline chain re-verified byte-identical).

SUBCOMMANDS
    genesets   the 2,334-gene universe (every S1C Gene ID) and the 1,793 family-eligible genes
               (S1C `In Table S1 (SD98 gene set)` = Yes), as 2-column `gene_id biotype` TSVs sorted by id.
               New: the hand-made files it replaces had no generator. They drive `cluster` byte-identically.
    curate     OPT-IN: SD98 regions (+ CAT v4) -> SD98 gene set -> explicit curation rule -> eligible + full
               genesets for `edges`/`cluster` (see CURATION ARM above; `--rule none` = uncurated).
    edges      steps 1-4: SEDEF rows >= 0.98 identity -> lift both sides to CHM13 v1.0 -> walk the SEDEF CIGAR
               -> project CAT v4 exons >= 0.99 covered across the pair -> `gene_a gene_b` edge TSV (sorted).
    cluster    steps 5-6: connected components -> famCN MAD split -> family call; --full-geneset adds the
               single-eligible-seed islands (§6im) and attaches non-eligible members (§6ii).
    dennislab  the Dennis-lab notebook's algorithm (§6if Finding 2: ARI 0.541 mean / 0.564 median, worse;
               PARKED, kept as evidence that "their real algorithm" was run).
    famcn      WSSD read-depth famCN at arbitrary CHM13 v2.0 intervals (needs pyBigWig or bigBedToBed). Kept
               as infrastructure; the §6io CN weak-edge lever that used it was NOT shipped.
    score      fixed-universe ARI / exact-family / pair P-R-F1 (`--only pairs`), then Hungarian 1:1 family
               matching (`--only bipartite`); default `--only all` prints both, in that order.

IN-REPO INPUTS (resolved relative to this file, not the CWD)
    soto_famCN_S1C.tsv                 Soto Table S1C (truth, famCN, biotypes; pinned by REPRODUCE.md)
    soto_parCN_S1E.tsv                 Soto Table S1E (dual v1.0/v2.0 coordinates -> the liftover); restored
                                       from cd37ccb0^ (wave 3 archived it and never restored it)
    acro_extra_anchors.tsv             the §6il/§6in acrocentric liftover anchors; restored from cd37ccb0^
    shared_exons_2334_finalhuman.tsv   frozen `edges` output (§6ip, 4,192 edges), so steps 5-6 and the
                                       scorers reproduce from a clone without the 54 MB SEDEF BED and the
                                       84 MB CAT BED, which are not committed

ENVIRONMENT. `famcn --tool pybigwig` needs pyBigWig (/home/juanfra/miniforge3/bin/python3 has it; the
linuxbrew python3 first on PATH does not). `score` needs scikit-learn, numpy and scipy. genesets, edges,
cluster and dennislab are stdlib-only. Heavy imports happen inside the functions that need them.

OLD -> NEW (wave 7, 2026-09-24; the old scripts are at tag notebook-2026-09-24, flags unchanged unless noted)
    bench/soto/soto_replicate_from_sedef.py ...        -> soto_replication.py edges ...
                                                          (--s1e now optional, default bench/soto/soto_parCN_S1E.tsv)
    bench/soto/soto_cluster_from_shared.py ...         -> soto_replication.py cluster ...
                                                          (default --mad-statistic mean kept: the paper's prose)
    bench/soto/soto_cluster_dennislab_algorithm.py ... -> soto_replication.py dennislab ...
                                                          (default --mad-statistic median kept)
    bench/soto/famcn_from_wssd.py ...                  -> soto_replication.py famcn ...
                                                          (--s1e default now module-relative, same file)
    bench/soto/soto_score_against_truth.py ...         -> soto_replication.py score --only pairs ...
    bench/soto/soto_bipartite_match_score.py ...       -> soto_replication.py score --only bipartite ...
    (both scorers back to back)                        -> soto_replication.py score ...
    (--truth was required in both scorers)             -> optional, default bench/soto/soto_famCN_S1C.tsv
    bench/soto/soto_attach_noncoding_members.py ...    -> DROPPED; use `cluster --full-geneset` (§6ii verified
                                                          the same partition, ARI 0.6820). attach() is kept.
    bench/soto/soto_replicate_clustering.py ...        -> DROPPED (minimap2 map-back path, superseded by the
                                                          SEDEF path, §6ie; register 946). load_exons() is kept.
    bench/soto/rustlib.py                              -> REMOVED from the tree (0 importers, not Soto-specific);
                                                          `git show notebook-2026-09-24:bench/soto/rustlib.py`
    bench/soto_vs_us_referee.py                        -> bench/score.py referee (wave 7, step P1)

    functions:
    famcn_from_wssd.{build_liftover, lift, cn_for_interval, GAP_GUARD, BASE_URL}       -> same names
    soto_replicate_clustering.load_exons                                                -> load_exons
    soto_replicate_clustering.{_cigar, _project, parse_region, open_maybe_gz, main}     -> dropped
                                                                                           (_cigar == cigar_ops)
    soto_replicate_from_sedef.{cigar_ops, build_blocks, bwalk_to_genomic, genomic_to_bwalk,
                               project_a_to_bwalk, project_bwalk_to_a, find_shared_exons} -> same names
    soto_cluster_from_shared.{ELIGIBLE, mad_mean, mad_median, step5_step6}              -> same names
    inline islands block of soto_cluster_from_shared.main                              -> single_seed_islands
    inline attach/rewrite block of soto_cluster_from_shared.main                        -> attach_extra_members
    soto_attach_noncoding_members.attach                                                -> attach (also takes
                                                                                           an edge list)
    soto_cluster_dennislab_algorithm.{raw_components, dennislab_families}               -> same names
                                                                                           (over components())
    soto_score_against_truth.{load_truth, pairs_of, score}                              -> same names
    soto_bipartite_match_score.{build_families, main body}                              -> build_families,
                                                                                           bipartite_score
    the 3 famCN loaders, 6 geneset loaders, 4 edge-file loaders, 2 anchor loaders and 2 predicted-assignment
    loaders                                         -> load_famcn, load_geneset, read_edges + adjacency /
                                                       load_shared, load_extra_anchors, load_predicted +
                                                       restrict_universe
    the 4 DFS component loops                       -> components() (discovery order preserved: SEDEFFAM ids
                                                       depend on it)
"""
import argparse, csv, os, re, statistics as st, subprocess, sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
S1C = os.path.join(HERE, "soto_famCN_S1C.tsv")
S1E = os.path.join(HERE, "soto_parCN_S1E.tsv")
ANCHORS = os.path.join(HERE, "acro_extra_anchors.tsv")
FROZEN_EDGES = os.path.join(HERE, "shared_exons_2334_finalhuman.tsv")

BASE_URL = "http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/wssd"
GAP_GUARD = 50_000  # v2 bp around a regime switch treated as unmappable

ELIGIBLE = {"protein_coding", "unprocessed_pseudogene",
            "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene"}


# ---------------------------------------------------------------------------------------------------------
# Loaders (each replaces 2-6 byte-for-byte copies in the old scripts)
# ---------------------------------------------------------------------------------------------------------

def load_geneset(path):
    """A geneset TSV (`gene_id`, optional `biotype`) -> (set of gene ids, {gene: biotype or ''})."""
    genes, biotype = set(), {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
            biotype[r["gene_id"]] = r.get("biotype", "")
    return genes, biotype


def load_famcn(path):
    """gene -> famCN from any of the famCN tables used here (S1C `Gene ID`/`Median famCN`, or a
    `gene_id` + `famCN`/`famCN_median` table). Rows without a numeric value are skipped."""
    famcn = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r.get("gene_id") or r.get("Gene ID")
            v = r.get("famCN") or r.get("famCN_median") or r.get("Median famCN")
            try:
                famcn[gid] = float(v)
            except (TypeError, ValueError):
                pass
    return famcn


def read_edges(path):
    """The `gene_a gene_b` shared-exon edge TSV as a list of (a, b) pairs, in file order."""
    with open(path) as fh:
        return [(r["gene_a"], r["gene_b"]) for r in csv.DictReader(fh, delimiter="\t")]


def adjacency(edges, keep=None):
    """gene -> set(partner genes), keys in first-seen order. With `keep`, only edges whose BOTH ends are in
    `keep` (the cluster backbone); without, every edge (the dennislab path reads the file unfiltered)."""
    shared = defaultdict(set)
    for ga, gb in edges:
        if keep is None or (ga in keep and gb in keep):
            shared[ga].add(gb)
            shared[gb].add(ga)
    return shared


def load_shared(path, keep=None):
    return adjacency(read_edges(path), keep)


def load_extra_anchors(path):
    """--extra-anchors TSV (chrom, v2_pos, offset) -> list of triples, or None when no path is given."""
    if not path:
        return None
    extra_anchors = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            extra_anchors.append((r["chrom"], int(r["v2_pos"]), int(r["offset"])))
    return extra_anchors


def load_predicted(path):
    """A predicted assignment TSV (gene_id, family_id; empty family_id = unplaced) -> {gene: family_id}."""
    predicted = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("family_id"):
                predicted[r["gene_id"]] = r["family_id"]
    return predicted


def restrict_universe(clean_truth, eligible_only_universe=None):
    universe = set(clean_truth)
    if eligible_only_universe:
        with open(eligible_only_universe) as fh:
            restrict = {r["gene_id"] for r in csv.DictReader(fh, delimiter="\t")}
        universe &= restrict
    return universe


def load_truth(path):
    """Returns (gene_to_family: clean single-family genes only, ambiguous: set of excluded gene ids)."""
    families = defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fid = r["Family ID"]
            if fid:
                families[r["Gene ID"]].add(fid)
    ambiguous = {g for g, f in families.items() if len(f) > 1}
    clean = {g: next(iter(f)) for g, f in families.items() if len(f) == 1}
    return clean, ambiguous


# ---------------------------------------------------------------------------------------------------------
# genesets
# ---------------------------------------------------------------------------------------------------------

def soto_genesets(truth_path):
    """From S1C: (full, eligible) as sorted lists of (gene_id, biotype). full = every distinct `Gene ID`
    (2,334); eligible = the genes with `In Table S1 (SD98 gene set)` = Yes (1,793). The biotype is the
    gene's `Biotype`; S1C has no gene whose rows disagree on either column (checked 2026-09-24), and a gene
    that did would take its first row's value."""
    biotype, eligible = {}, set()
    with open(truth_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["Gene ID"]
            if not g:
                continue
            biotype.setdefault(g, r["Biotype"])
            if r["In Table S1 (SD98 gene set)"] == "Yes":
                eligible.add(g)
    full = sorted(biotype.items())
    return full, [(g, b) for g, b in full if g in eligible]


def write_geneset(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype"])
        w.writerows(rows)


# ---------------------------------------------------------------------------------------------------------
# SD98 gene set from SD98 regions + the (reverse-engineered) manual curation step (2026-09-28, KEY=rule)
# ---------------------------------------------------------------------------------------------------------
# Soto's step 2 (A_SD98_regions.md, their released code): SD98 gene = CAT v4 gene with >= 1 exon FULLY
# contained (`bedtools intersect -f 1`) in a merged SD98 region (chrY/chrM excluded; chrX dropped at the gene-
# count step); family-eligible = biotype grep "protein_coding\|unprocessed_pseudogene". Then redundancy
# candidates = self-intersect of those genes' transcripts, `bedtools intersect -wao -s -f 0.9`, different
# gene id; then MANUAL curation ("71 genes ... removed"; "kept only fusion genes that encoded an alternative
# protein"). `genes_to_remove.txt` was never released. The functions below recompute the gene set and the
# per-gene vocabulary the curation could have used, and `curation_drop` applies an explicit rule fitted on
# one random half of the genes and validated on the other (scratchpad/soto_cur/rule_search.py; report
# figs/soto_cur_rule.md). Everything here is NEW code reached only through `curate`; no existing subcommand
# calls it (byte-identical when unused).

ELIGIBLE_GREP = ("protein_coding", "unprocessed_pseudogene")   # their literal grep (substring match)
BIOTYPE_RANK = {"unprocessed_pseudogene": 0, "transcribed_unprocessed_pseudogene": 1,
                "translated_unprocessed_pseudogene": 1, "protein_coding": 2}
NAME_RANK = {"clone": 0, "MSTRG": 0, "LOC": 0, "ENSG": 0, "readthrough": 1, "named": 2}
CURATION_RULES = ("none", "redundant", "contained", "contained_or_flag", "contained_outranked",
                  "contained_or_flag_outranked")


def name_class(name):
    """Ensembl clone-accession names (AC012345.1, AL627309.3 ...) -> 'clone'; GENE1-GENE2 readthrough names ->
    'readthrough' (antisense/divergent -AS1/-DT/-IT/-OT suffixes and X-1 style are 'named'); MSTRG/LOC/ENSG
    placeholders keep their own class; everything else is 'named' (HGNC-style)."""
    if name.startswith("MSTRG"):
        return "MSTRG"
    if name.startswith("LOC"):
        return "LOC"
    if name.startswith("ENSG"):
        return "ENSG"
    if re.match(r"^(A[CLPF]|BX|CR|CU|FP|Z)\d{5,6}\.\d+$", name) or re.match(r"^[A-Z]{1,2}\d{6}\.\d+$", name):
        return "clone"
    if "-" in name and not re.search(r"-(AS\d*|DT|IT\d*|OT\d*)$", name):
        a, b = name.split("-", 1)
        if len(a) >= 2 and len(b) >= 2 and not re.fullmatch(r"[0-9]+|[A-Z]", b):
            return "readthrough"
    return "named"


def merge_intervals(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def overlap_bp(u1, u2):
    """Total overlap between two merged interval lists."""
    i = j = t = 0
    while i < len(u1) and j < len(u2):
        a, b = max(u1[i][0], u2[j][0]), min(u1[i][1], u2[j][1])
        if b > a:
            t += b - a
        if u1[i][1] < u2[j][1]:
            i += 1
        else:
            j += 1
    return t


def load_cat_genes(cat_bed):
    """CAT v4 BED (37 columns) -> gene_id -> dict(name, biotype, chrom, strand, start, end, source, tx (list of
    (tx_start, tx_end, exons, cds_exons)), exon_set (distinct exons), exons (merged union), cds (merged union))."""
    genes = {}
    with open(cat_bed) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            gid = f[18]
            s, e, ts, te = int(f[1]), int(f[2]), int(f[6]), int(f[7])
            sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
            exs = [(s + so, s + so + sz) for sz, so in zip(sizes, starts)]
            cds = [(max(a, ts), min(b, te)) for a, b in exs if te > ts and min(b, te) > max(a, ts)]
            g = genes.get(gid)
            if g is None:
                g = genes[gid] = {"name": f[12], "biotype": f[19], "chrom": f[0], "strand": f[5], "start": s,
                                  "end": e, "source": f[32] if len(f) > 32 else "", "tx": [],
                                  "_ex": [], "_cds": []}
            g["start"], g["end"] = min(g["start"], s), max(g["end"], e)
            g["tx"].append((s, e, exs, cds))
            g["_ex"].extend(exs)
            g["_cds"].extend(cds)
    for g in genes.values():
        g["exon_set"] = sorted(set(g.pop("_ex")))       # distinct exons (their exons BED; the `-f 1` test is per exon)
        g["exons"] = merge_intervals(g["exon_set"])     # exonic union (bp-level containment features)
        g["cds"] = merge_intervals(g.pop("_cds"))
    return genes


def load_regions_bed(path):
    """BED of regions -> chrom -> merged, sorted interval list."""
    per = defaultdict(list)
    with open(path) as fh:
        for ln in fh:
            f = ln.split("\t")
            if len(f) < 3 or ln.startswith("#"):
                continue
            per[f[0]].append((int(f[1]), int(f[2])))
    return {c: merge_intervals(v) for c, v in per.items()}


def sd98_regions_from_sedef(sedef_path, min_identity=0.98):
    """Their step 1 on a native SEDEF table: every row with identity (field 21) >= min_identity contributes
    BOTH sides; merged per chromosome (`bedtools merge` equivalent)."""
    per = defaultdict(list)
    with open(sedef_path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            try:
                ident = float(f[20])
            except ValueError:
                continue
            if ident < min_identity:
                continue
            per[f[0]].append((int(f[1]), int(f[2])))
            per[f[3]].append((int(f[4]), int(f[5])))
    return {c: merge_intervals(v) for c, v in per.items()}


def _fully_contained(regions, s, e):
    """True if [s, e) lies inside one merged region (bisect on the sorted merged list)."""
    import bisect
    i = bisect.bisect_right(regions, (s, float("inf"))) - 1
    return i >= 0 and regions[i][0] <= s and e <= regions[i][1]


def sd98_gene_set(genes, regions, autosomal=True):
    """Their step 2: gene ids with >= 1 exon fully contained in an SD98 region (`-f 1`); chrY/chrM never count;
    chrX excluded when autosomal (their gene-count step). Returns {gene_id: n_contained_exons}."""
    out = {}
    for gid, g in genes.items():
        if g["chrom"] in ("chrY", "chrM") or (autosomal and g["chrom"] == "chrX"):
            continue
        regs = regions.get(g["chrom"])
        if not regs:
            continue
        n = sum(1 for s, e in g["exon_set"] if _fully_contained(regs, s, e))
        if n:
            out[gid] = n
    return out


def is_eligible_biotype(bt):
    return any(k in bt for k in ELIGIBLE_GREP)


def curation_features(genes, elig):
    """Per gene in `elig` (the eligible SD98 gene ids), the vocabulary of their curation step, computed against
    PARTNERS = other eligible SD98 genes on the same chromosome and strand whose gene bodies overlap:
      flag                their literal candidate flag: some transcript of g lies >= 90% (of its own length) inside a
                          transcript of a different gene (`bedtools intersect -wao -s -f 0.9`, self-intersect)
      flag_by             the set of container genes behind `flag`
      partners            [(h, frac_of_g_exonic_bp_inside_h, g_covers_h_90)] for partners sharing >= 1 exonic bp
      exon_contain_frac   max over partners of exonic bp of g inside h / exonic bp of g
      exon_novel_frac     exonic bp of g outside every partner / exonic bp of g (1.0 = no partner)
      cds_novel_frac      coding bp of g outside every partner's CDS / coding bp of g (None = no CDS) --
                          "encodes an alternative protein" proxy
      n_covered           partners whose exon union is >= 90% covered by g's (fusion/readthrough container)
      body_contained      gene body >= 90% inside a partner's body
      name_class, biotype, n_tx, n_exons, span, has_cds, source (annotation facts, for the rank)
    """
    by_key = defaultdict(list)
    for gid in elig:
        g = genes[gid]
        by_key[(g["chrom"], g["strand"])].append(gid)
    feats = {}
    for gid in elig:
        g = genes[gid]
        s, e = g["start"], g["end"]
        ps = [h for h in by_key[(g["chrom"], g["strand"])]
              if h != gid and genes[h]["start"] < e and genes[h]["end"] > s]
        ex_len = sum(b - a for a, b in g["exons"])
        cds_len = sum(b - a for a, b in g["cds"])
        partners, flag_by, body_contained = [], set(), False
        for h in ps:
            gh = genes[h]
            ov = overlap_bp(g["exons"], gh["exons"])
            if ov:
                covers = ov >= 0.9 * sum(b - a for a, b in gh["exons"])
                partners.append((h, ov / ex_len, covers))
            if min(e, gh["end"]) - max(s, gh["start"]) >= 0.9 * (e - s):
                body_contained = True
            for ts, te, _, _ in g["tx"]:
                if any(min(te, ue) - max(ts, us) >= 0.9 * (te - ts) for us, ue, _, _ in gh["tx"]):
                    flag_by.add(h)
                    break
        shared = [h for h, _, _ in partners]
        ex_union = merge_intervals([iv for h in shared for iv in genes[h]["exons"]])
        cds_union = merge_intervals([iv for h in shared for iv in genes[h]["cds"]])
        feats[gid] = {
            "flag": bool(flag_by), "flag_by": flag_by, "partners": partners,
            "exon_contain_frac": max([fr for _, fr, _ in partners] + [0.0]),
            "exon_novel_frac": 1 - overlap_bp(g["exons"], ex_union) / ex_len if ex_len else 1.0,
            "cds_novel_frac": (1 - overlap_bp(g["cds"], cds_union) / cds_len) if cds_len else None,
            "n_covered": sum(1 for _, _, c in partners if c), "n_partners": len(partners),
            "body_contained": body_contained, "name_class": name_class(g["name"]), "biotype": g["biotype"],
            "n_tx": len(g["tx"]), "n_exons": len(g["exon_set"]), "span": e - s, "has_cds": cds_len > 0,
            "source": g["source"],
        }
    return feats


def curation_rank(f):
    """The annotation-quality order the curation appears to follow: protein_coding > transcribed_unprocessed
    > unprocessed pseudogene; HGNC-named > readthrough-named > clone/placeholder-named; more transcripts;
    longer span."""
    return (BIOTYPE_RANK.get(f["biotype"], 0), NAME_RANK.get(f["name_class"], 0), f["n_tx"], f["span"])


def curation_drop(feats, rule, contain=0.5):
    """Gene ids the rule removes. Rules (all restricted to genes with >= 1 same-strand exon-sharing partner):
      none                          nothing (the uncurated recomputed set)
      redundant                     their literal candidate flag alone
      contained                     >= `contain` of exonic bp inside one partner
      contained_or_flag             the above OR their flag
      contained_outranked           >= `contain` of exonic bp inside a partner that OUTRANKS it (curation_rank)
      contained_or_flag_outranked   (>= `contain` inside h OR flagged by h) AND h outranks it
    """
    if rule == "none":
        return set()
    drop = set()
    for gid, f in feats.items():
        if rule == "redundant":
            if f["flag"]:
                drop.add(gid)
            continue
        rk = curation_rank(f)
        for h, fr, _ in f["partners"]:
            hit = fr >= contain
            if rule in ("contained_or_flag", "contained_or_flag_outranked"):
                hit = hit or (h in f["flag_by"])
            if hit and (rule in ("contained", "contained_or_flag") or curation_rank(feats[h]) > rk):
                drop.add(gid)
                break
    return drop


# ---------------------------------------------------------------------------------------------------------
# Liftover CHM13 v2.0 -> v1.0 from S1E's dual coordinates, and WSSD famCN (was famcn_from_wssd.py)
# ---------------------------------------------------------------------------------------------------------

def build_liftover(s1e_path, extra_anchors=None):
    """Fit a piecewise-constant v2.0 -> v1.0 map from Soto's dual-coordinate anchors.

    `extra_anchors`: OPTIONAL iterable of (chrom, v2_pos, offset) triples to merge in alongside S1E's own
    anchors, exactly as if they were additional dual-coordinate rows -- e.g. from a direct minimap2
    realignment of a specific gene's own sequence against both assemblies, used where Soto's own SD98-
    paralog anchors happen to be sparse enough that a real, single-breakpoint acrocentric region gets
    treated as one wide "uncertain span" spanning both sides of the true (and already correctly known)
    breakpoint (docs/o1_ledger.md §6il). Omit for the original behaviour, byte-identical to before this
    parameter existed -- this is validated per-locus evidence, not a blanket relaxation of the guard.

    Returns {chrom: [(v2_from, offset), ...]} sorted by v2_from, plus the anchor count per chromosome.
    """
    anchors = defaultdict(list)
    with open(s1e_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m1 = re.match(r"(chr[\w]+):(\d+)-(\d+)", r.get("SD98_v1.0", "") or "")
            m2 = re.match(r"(chr[\w]+):(\d+)-(\d+)", r.get("SD98_v2.0", "") or "")
            if not (m1 and m2) or m1.group(1) != m2.group(1):
                continue
            anchors[m1.group(1)].append((int(m2.group(2)), int(m1.group(2)) - int(m2.group(2))))
    for chrom, v2_pos, offset in (extra_anchors or ()):
        anchors[chrom].append((v2_pos, offset))
    table, spans = {}, {}
    for c, pts in anchors.items():
        pts.sort()
        regimes, prev = [], None
        for v2, off in pts:
            if off != prev:
                regimes.append((v2, off))
                prev = off
        table[c] = regimes
        # the v2 window where a regime switch happens is uncertain: between the last anchor of one
        # regime and the first of the next.
        bounds = []
        for i in range(1, len(regimes)):
            lo = max(v2 for v2, o in pts if o == regimes[i - 1][1])
            bounds.append((lo, regimes[i][0]))
        spans[c] = bounds
    return table, spans


def lift(table, spans, chrom, start, end):
    """v2.0 interval -> v1.0 interval, or None when the mapping is not trustworthy here."""
    regimes = table.get(chrom)
    if not regimes:
        return None
    for lo, hi in spans.get(chrom, ()):
        if start < hi + GAP_GUARD and end > lo - GAP_GUARD:
            return None  # straddles / sits inside a regime switch
    off = regimes[0][1]
    for v2_from, o in regimes:
        if start >= v2_from:
            off = o
    return (start + off, end + off)


def cn_for_interval(bb, chrom, start, end, tool):
    """Length-weighted mean CN over the WSSD windows covering [start, end).

    `tool="pybigwig"`: OPT-IN alternate reader added when `bigBedToBed` was not installed on this machine.
    Reads the SAME bigBed files in-process via the `pyBigWig` package instead of shelling out --
    `.entries()` returns the identical (start, end, tab-separated-rest) rows `bigBedToBed` would print,
    just with chrom/start/end already split out; the CN value is the LAST tab field either way (column 10
    of bigBedToBed's 1-based output == `rest.split("\t")[-1]` here). Verified to reproduce
    `bigBedToBed`-derived famCN (via the pre-existing famcn_ours_all.tsv, same 10 samples) within
    floating-point rounding before being trusted for anything new (docs/o1_ledger.md, follow-up to §6in).
    """
    if tool == "pybigwig":
        import pyBigWig
        try:
            bw = pyBigWig.open(bb)
        except RuntimeError:
            return None
        try:
            entries = bw.entries(chrom, start, end) or []
        except RuntimeError:
            entries = []
        finally:
            bw.close()
        tot = n = 0.0
        for s, e, rest in entries:
            try:
                cn = float(rest.rsplit("\t", 1)[-1])
            except ValueError:
                continue
            w = min(e, end) - max(s, start)
            if w > 0:
                tot += cn * w
                n += w
        return tot / n if n > 0 else None

    try:
        out = subprocess.run(
            [tool, f"-chrom={chrom}", f"-start={start}", f"-end={end}", bb, "/dev/stdout"],
            capture_output=True, text=True, timeout=300,
        ).stdout
    except subprocess.TimeoutExpired:
        return None
    tot = n = 0.0
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 10:
            continue
        try:
            s, e, cn = int(f[1]), int(f[2]), float(f[9])
        except ValueError:
            continue
        w = min(e, end) - max(s, start)
        if w > 0:
            tot += cn * w
            n += w
    return tot / n if n > 0 else None


# ---------------------------------------------------------------------------------------------------------
# Steps 1-4: SEDEF CIGAR -> shared-exon edges (was soto_replicate_from_sedef.py + load_exons)
# ---------------------------------------------------------------------------------------------------------

def load_exons(cat_bed, keep_genes):
    """gene -> list of (chrom, start, end) exons, from the CAT bigBed-derived BED."""
    ex = defaultdict(set)
    meta = {}
    with open(cat_bed) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            gid = f[18]
            if gid not in keep_genes:
                continue
            meta[gid] = (f[12], f[19])
            c, st = f[0], int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
            for sz, so in zip(sizes, starts):
                ex[gid].add((c, st + so, st + so + sz))
    return ex, meta


def cigar_ops(cg):
    n = 0
    for ch in cg:
        if ch.isdigit():
            n = n * 10 + (ord(ch) - 48)
        else:
            yield n, ch
            n = 0


def build_blocks(start1, ops):
    """Walk a SEDEF CIGAR (D consumes side1 only, I consumes side2 only, M/=/X consumes both).

    Returns (blocks, p1_end, bwalk_end). `blocks` is a list of (a_from, a_to, bwalk_from, bwalk_to)
    for every matched run: `a` is side1's plain genomic position (side1 is always '+'). `bwalk` is a
    strand-agnostic walk counter for side2 -- 0 at the alignment's own start, increasing monotonically
    with the CIGAR regardless of strand2 -- converted to a genomic position only via
    bwalk_to_genomic/genomic_to_bwalk, so the block list itself never needs to know strand2.
    p1_end/bwalk_end let the caller verify the CIGAR reconciles with the row's own interval lengths
    (p1_end must equal end1, bwalk_end must equal end2-start2) before trusting the row.
    """
    p1 = start1
    bw = 0
    blocks = []
    for n, op in ops:
        if op in "M=X":
            blocks.append((p1, p1 + n, bw, bw + n))
            p1 += n
            bw += n
        elif op == "D":
            p1 += n
        elif op == "I":
            bw += n
        # S/H/P not expected in an internal SEDEF self-alignment CIGAR; ignored if present.
    return blocks, p1, bw


def bwalk_to_genomic(bw, start2, end2, strand2):
    return start2 + bw if strand2 == "+" else end2 - bw


def genomic_to_bwalk(pos, start2, end2, strand2):
    return pos - start2 if strand2 == "+" else end2 - pos


def project_a_to_bwalk(blocks, a_lo, a_hi):
    """Project a-axis interval [a_lo,a_hi) through the matched blocks to a bwalk interval, or None."""
    lo = hi = None
    for a_from, a_to, bw_from, _bw_to in blocks:
        if a_to <= a_lo:
            continue
        if a_from >= a_hi:
            break
        if lo is None:
            lo = bw_from + max(a_lo, a_from) - a_from
        hi = bw_from + min(a_hi, a_to) - a_from
    return (lo, hi) if lo is not None and hi is not None and hi > lo else None


def project_bwalk_to_a(blocks, bw_lo, bw_hi):
    """Project a bwalk interval [bw_lo,bw_hi) through the matched blocks to an a-axis interval, or None."""
    lo = hi = None
    for a_from, _a_to, bw_from, bw_to in blocks:
        if bw_to <= bw_lo:
            continue
        if bw_from >= bw_hi:
            break
        if lo is None:
            lo = a_from + max(bw_lo, bw_from) - bw_from
        hi = a_from + min(bw_hi, bw_to) - bw_from
    return (lo, hi) if lo is not None and hi is not None and hi > lo else None


def find_shared_exons(blocks, side2_start_v2, side2_end_v2, side2_strand,
                       own_chrom, own_lo_v1, own_hi_v1, own_offset, own_axis_is_bwalk,
                       other_chrom, other_lo_v1, other_hi_v1, other_offset, other_axis_is_bwalk,
                       per_chrom, min_cov, shared):
    """For genes on the OTHER side (v1.0 space, restricted to [other_lo_v1,other_hi_v1)) whose exon
    is covered >=min_cov by the whole aligned span, project the exon through `blocks` onto the OWN
    side and link to whichever OWN-side gene's exon it lands in. Exactly one of
    own_axis_is_bwalk/other_axis_is_bwalk must be True: `blocks`' "a" axis is always side1 (plain),
    its "bwalk" axis is always side2 (strand-aware) -- call once with own=side1/other=side2 and once
    with the roles swapped. Returns the number of successful exon projections (for reporting only).
    """
    hits = per_chrom.get(other_chrom)
    if not hits:
        return 0
    n_proj = 0
    for s, e, g_other in hits:
        if e <= other_lo_v1:
            continue
        if s >= other_hi_v1:
            break
        if (min(e, other_hi_v1) - max(s, other_lo_v1)) / (e - s) < min_cov:
            continue
        v2_s, v2_e = max(s, other_lo_v1) - other_offset, min(e, other_hi_v1) - other_offset
        if other_axis_is_bwalk:
            bw_lo, bw_hi = sorted((
                genomic_to_bwalk(v2_s, side2_start_v2, side2_end_v2, side2_strand),
                genomic_to_bwalk(v2_e, side2_start_v2, side2_end_v2, side2_strand),
            ))
            proj = project_bwalk_to_a(blocks, bw_lo, bw_hi)
            if proj is None:
                continue
            p_v2_lo, p_v2_hi = proj
        else:
            proj = project_a_to_bwalk(blocks, v2_s, v2_e)
            if proj is None:
                continue
            g_lo = bwalk_to_genomic(proj[0], side2_start_v2, side2_end_v2, side2_strand)
            g_hi = bwalk_to_genomic(proj[1], side2_start_v2, side2_end_v2, side2_strand)
            p_v2_lo, p_v2_hi = min(g_lo, g_hi), max(g_lo, g_hi)
        p_lo, p_hi = p_v2_lo + own_offset, p_v2_hi + own_offset
        # NOTE: unlike the retired minimap2 map-back path (soto_replicate_clustering.py, tag
        # notebook-2026-09-24), where a region's own trivial self-alignment to its OWN originating location
        # must be excluded, there is no analogous
        # "self-mapping" artifact here: side1 and side2 are two DISTINCT loci SEDEF itself reported as
        # a real duplication pair, so a projection from side2 landing inside side1's own span is the
        # EXPECTED, correct outcome (that is literally where side1's gene lives), not a self-hit. The
        # only real self-link risk -- a gene ending up linked to itself -- is guarded below by
        # `g_own != g_other`. (An earlier version of this function incorrectly rejected every
        # same-chromosome pair here, since a projection landing inside its own OWN span is always
        # true; caught via a same-chromosome case, ID_2/CHM13_G0020704 x CHM13_G0020810 on chr16,
        # that should have linked and did not until this was removed.)
        n_proj += 1
        for s2, e2, g_own in per_chrom.get(own_chrom, ()):
            if e2 <= p_lo:
                continue
            if s2 >= p_hi:
                break
            if g_own != g_other:
                shared[g_own].add(g_other)
                shared[g_other].add(g_own)
    return n_proj


# ---------------------------------------------------------------------------------------------------------
# Steps 5-6: components -> famCN MAD split -> family call (was soto_cluster_from_shared.py + attach)
# ---------------------------------------------------------------------------------------------------------

def components(adj):
    """Connected components of `adj` (gene -> set of partners, symmetric) by iterative DFS, in DISCOVERY
    order: keys are scanned in insertion order and each component lists genes in visit order. The order
    matters -- step5_step6 numbers families by it (a stable sort by size keeps ties in this order), so do
    not substitute a union-find or a sorted variant, or SEDEFFAM ids renumber."""
    seen, comps = set(), []
    for g in adj:
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


def mad_mean(vals):
    """Mean absolute deviation about the mean -- the paper's own METHODS-text wording ("mean absolute
    deviation")."""
    if len(vals) < 2:
        return 0.0
    m = sum(vals) / len(vals)
    return sum(abs(v - m) for v in vals) / len(vals)


def mad_median(vals):
    """Median absolute deviation about the median, UNSCALED (scale=1.0) -- what their actual released
    code computes (B_SD98_families.ipynb: `stats.median_abs_deviation(wssd_clust_median)`, no `scale=`
    override, so scipy's own default of 1.0 applies -- NOT the "mean absolute deviation" the paper's own
    prose describes). Confirmed via two independent reads of the notebook's raw source. Robust to the
    extreme-CN outliers this project has already documented in famCN distributions (e.g. BET1L=702),
    which a mean-based statistic is not.
    """
    if len(vals) < 2:
        return 0.0
    s = sorted(vals)
    n = len(s)
    med = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    dev = sorted(abs(v - med) for v in vals)
    return dev[n // 2] if n % 2 else (dev[n // 2 - 1] + dev[n // 2]) / 2


def step5_step6(shared, genes, biotype, famcn, mad_threshold, out_path, mad_fn=mad_mean):
    """shared: gene -> set(partner genes). genes: set of all SD98 gene ids (for singleton bookkeeping).
    biotype: gene -> biotype string. famcn: gene -> float famCN (genes absent are treated as un-splittable).
    Writes `out_path` (gene_id, biotype, family_id, n_members, famCN, status). The split is the retired
    soto_replicate_clustering.py's inline steps 5-6 (same MAD formula, same biotype eligibility set, same
    singleton bookkeeping); that script also wrote a gene_name column and REPFAM ids.
    """
    comps = components(shared)
    seen = {x for comp in comps for x in comp}

    final = []
    for comp in comps:
        vals = [(famcn[g], g) for g in comp if g in famcn]
        if len(vals) < 2 or mad_fn([v for v, _ in vals]) < mad_threshold:
            final.append(comp)
            continue
        vals.sort()
        cur, groups = [], []
        for v, g in vals:
            if cur and mad_fn([x for x, _ in cur] + [v]) >= mad_threshold:
                groups.append([g2 for _, g2 in cur])
                cur = []
            cur.append((v, g))
        if cur:
            groups.append([g2 for _, g2 in cur])
        nocn = [g for g in comp if g not in famcn]
        if groups and nocn:
            groups[0].extend(nocn)
        final.extend(groups)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        fam_i = n_fam = n_single = 0
        for comp in sorted(final, key=lambda c: -len(c)):
            eligible = any(biotype.get(g) in ELIGIBLE for g in comp)
            if len(comp) >= 2 and eligible:
                fid = f"SEDEFFAM{fam_i}"
                fam_i += 1
                n_fam += 1
                status = "family"
            else:
                fid = ""
                status = "singleton" if len(comp) < 2 else "no_coding_member"
                n_single += len(comp)
            for g in sorted(comp):
                w.writerow([g, biotype.get(g, ""), fid, len(comp),
                            f"{famcn[g]:.2f}" if g in famcn else "", status])
        for g in sorted(genes - seen):
            w.writerow([g, biotype.get(g, ""), "", 1, f"{famcn[g]:.2f}" if g in famcn else "", "singleton"])
            n_single += 1
    return n_fam, n_single, len(comps)


def single_seed_islands(shared, full_adj, genes):
    """SINGLE-ELIGIBLE-SEED ISLANDS (docs/o1_ledger.md §6im): every raw (full, unfiltered) connected
    component of `full_adj` that touches EXACTLY ONE gene of `genes` has all its edges admitted into
    `shared` (mutated in place) before clustering. Returns (n_islands, n_island_genes)."""
    n_islands, n_island_genes = 0, 0
    for comp in components(full_adj):
        if sum(1 for x in comp if x in genes) == 1:
            n_islands += 1
            n_island_genes += len(comp)
            for x in comp:
                for y in full_adj[x]:
                    shared[x].add(y)
                    shared[y].add(x)
    return n_islands, n_island_genes


def attach(gene_family, extra_genes, shared_edges):
    """For each extra (non-eligible) gene with >=1 shared-exon edge to a gene already in gene_family,
    attach it to whichever family it has the MOST such edges with (a disclosed tie-break for the rare
    case of edges to more than one family -- Soto's own text does not specify one; an exact tie goes to
    the family seen first in edge-file order). `shared_edges` is the edge TSV path or its read_edges()
    list. Returns {gene: family_id} for attached genes only.
    """
    if isinstance(shared_edges, str):
        shared_edges = read_edges(shared_edges)
    edge_count = defaultdict(lambda: defaultdict(int))
    for a, b in shared_edges:
        if a in extra_genes and b in gene_family:
            edge_count[a][gene_family[b]] += 1
        if b in extra_genes and a in gene_family:
            edge_count[b][gene_family[a]] += 1
    return {g: max(fams.items(), key=lambda kv: kv[1])[0] for g, fams in edge_count.items() if fams}


def attach_extra_members(out_path, genes, full_genes_all, full_biotype, famcn, edges):
    """--full-geneset, second half: attach non-eligible members onto the backbone families step5_step6
    just wrote to `out_path`, then rewrite it with the combined result (same 6-column shape, extra genes
    marked in `status`). Returns (n_extra_considered, n_attached)."""
    gene_family = {}
    rows = []
    with open(out_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append(r)
            if r["family_id"]:
                gene_family[r["gene_id"]] = r["family_id"]

    # exclude genes already placed by the single-eligible-seed-island step (gene_family already has their
    # real family_id from the rows just read back) -- otherwise attach() would reconsider them and the
    # write loop below would emit a duplicate row for the same gene_id.
    extra_genes = full_genes_all - genes - set(gene_family)
    attached = attach(gene_family, extra_genes, edges)

    fam_size = defaultdict(int)
    for r in rows:
        if r["family_id"]:
            fam_size[r["family_id"]] += 1
    for f in attached.values():
        fam_size[f] += 1

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        for r in rows:
            n = fam_size[r["family_id"]] if r["family_id"] else int(r["n_members"])
            w.writerow([r["gene_id"], r["biotype"], r["family_id"], n, r["famCN"], r["status"]])
        for g in sorted(attached):
            f = attached[g]
            w.writerow([g, full_biotype.get(g, ""), f, fam_size[f],
                        f"{famcn[g]:.2f}" if g in famcn else "", "attached_noncoding_member"])
        for g in sorted(extra_genes - set(attached)):
            w.writerow([g, full_biotype.get(g, ""), "", 1,
                        f"{famcn[g]:.2f}" if g in famcn else "", "extra_gene_no_attachment"])
    return len(extra_genes), len(attached)


# ---------------------------------------------------------------------------------------------------------
# The Dennis-lab notebook algorithm (was soto_cluster_dennislab_algorithm.py) -- PARKED arm
# ---------------------------------------------------------------------------------------------------------

def raw_components(shared):
    """Connected components of the shared-exon graph -- the SAME initial object their pipeline calls
    `clusters` (loaded from their data/SD98_exon_clusters.txt), before any MAD filtering."""
    return [frozenset(comp) for comp in components(shared)]


def dennislab_families(comps, biotype, famcn, mad_threshold, mad_fn):
    """Faithful port of the algorithm quoted in cmd_dennislab's docstring. Returns a list of frozensets
    (deduplicated families) plus the set of genes that appear in some LOW-dispersion cluster but never
    became part of any final family (i.e. non-coding-only clusters with no coding bridge -- these keep
    membership information but never get a Family ID, matching "SD98 genes associated with other gene
    features... were also assigned a gene family ID" ONLY when they actually join one via a coding bridge)
    and the set of genes whose ONLY cluster was high-dispersion (discarded outright).
    """
    low, high_genes = [], set()
    for comp in comps:
        vals = [famcn[g] for g in comp if g in famcn]
        is_low = len(vals) < 2 or mad_fn(vals) < mad_threshold
        if is_low:
            low.append(comp)
        else:
            high_genes.update(comp)

    seeds = [g for c in low for g in c if biotype.get(g) in ELIGIBLE]

    families_raw = []
    for seed in seeds:
        cluster_list = [seed]
        seen_local = {seed}
        i = 0
        while i < len(cluster_list):
            if biotype.get(cluster_list[i]) in ELIGIBLE:
                for c in low:
                    if cluster_list[i] in c:
                        for g in c:
                            if g not in seen_local:
                                seen_local.add(g)
                                cluster_list.append(g)
            i += 1
        families_raw.append(frozenset(cluster_list))

    families = sorted(set(families_raw), key=lambda f: sorted(f))
    in_family = set().union(*families) if families else set()
    low_genes = set().union(*low) if low else set()
    orphaned_low = low_genes - in_family  # non-coding-only low-dispersion clusters: never bridged in
    return families, orphaned_low, high_genes - in_family


# ---------------------------------------------------------------------------------------------------------
# Scoring against S1C (was soto_score_against_truth.py + soto_bipartite_match_score.py)
# ---------------------------------------------------------------------------------------------------------

def pairs_of(label_map, genes):
    by_label = defaultdict(list)
    for g in genes:
        by_label[label_map[g]].append(g)
    pairs = set()
    for lbl, members in by_label.items():
        if lbl.startswith("__"):
            continue
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                pairs.add(frozenset((members[i], members[j])))
    return pairs


def score(truth, predicted, universe_genes):
    from sklearn.metrics import adjusted_rand_score

    scored = sorted(universe_genes)
    truth_labels, pred_labels = [], []
    for i, g in enumerate(scored):
        t = truth.get(g, "")
        truth_labels.append(t if (t and not t.startswith("Unassigned")) else f"__s{i}")
        p = predicted.get(g, "")
        pred_labels.append(p if p else f"__o{i}")
    ari = adjusted_rand_score(truth_labels, pred_labels)

    truth_map, pred_map = dict(zip(scored, truth_labels)), dict(zip(scored, pred_labels))
    truth_fam = defaultdict(set)
    for g, t in truth_map.items():
        if not t.startswith("__s"):
            truth_fam[t].add(g)
    pred_fam = defaultdict(set)
    for g, p in pred_map.items():
        if not p.startswith("__o"):
            pred_fam[p].add(g)
    pred_sets = {frozenset(v) for v in pred_fam.values()}
    n_exact = sum(1 for v in truth_fam.values() if frozenset(v) in pred_sets)

    truth_pairs, pred_pairs = pairs_of(truth_map, scored), pairs_of(pred_map, scored)
    tp = len(truth_pairs & pred_pairs)
    prec = tp / len(pred_pairs) if pred_pairs else 0.0
    rec = tp / len(truth_pairs) if truth_pairs else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return dict(n_genes=len(scored), ari=ari, n_exact=n_exact, n_truth_fam=len(truth_fam),
                n_pred_fam=len(pred_fam), precision=prec, recall=rec, f1=f1)


def build_families(label_map, universe):
    fam = defaultdict(set)
    for g in universe:
        lbl = label_map.get(g, "")
        if lbl:
            fam[lbl].add(g)
    return fam


def bipartite_score(clean_truth, ambiguous, predicted, universe, show_worst=5):
    """Hungarian 1:1 matching of predicted to true families (max total gene overlap); prints the MICRO /
    MACRO precision-recall and the undetected true families (the old soto_bipartite_match_score.py
    stdout). Returns the summary numbers as a dict."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    # "Unassigned_*" truth labels are Soto's own true singletons -- exclude them as TRUE families (a
    # singleton has no meaningful bipartite partner), but leave the genes in the universe for the
    # predicted side's precision accounting (a predicted family claiming one of these genes still pays a
    # precision cost if that gene isn't real overlap with any other true family it's matched against).
    true_fam = defaultdict(set)
    for g in universe:
        t = clean_truth.get(g, "")
        if t and not t.startswith("Unassigned"):
            true_fam[t].add(g)
    pred_fam = build_families(predicted, universe)

    true_ids = sorted(true_fam)
    pred_ids = sorted(pred_fam)
    n_true, n_pred = len(true_ids), len(pred_ids)
    print(f"scored universe: {len(universe)} genes ({len(ambiguous)} ambiguous multi-family genes excluded)")
    print(f"true families: {n_true}   predicted families: {n_pred}")

    # overlap[i][j] = |pred_fam[pred_ids[i]] & true_fam[true_ids[j]]|
    overlap = np.zeros((n_pred, n_true), dtype=int)
    for i, p in enumerate(pred_ids):
        pf = pred_fam[p]
        for j, t in enumerate(true_ids):
            overlap[i, j] = len(pf & true_fam[t])

    # linear_sum_assignment requires a square-ish cost matrix; pad with zero-overlap dummies so every
    # true family gets a (possibly dummy, zero-overlap) predicted partner and vice versa, then match on
    # -overlap (maximize overlap == minimize its negative).
    n = max(n_pred, n_true)
    padded = np.zeros((n, n), dtype=int)
    padded[:n_pred, :n_true] = overlap
    row_ind, col_ind = linear_sum_assignment(-padded)

    matches = []  # (pred_id_or_None, true_id_or_None, overlap)
    for i, j in zip(row_ind, col_ind):
        p = pred_ids[i] if i < n_pred else None
        t = true_ids[j] if j < n_true else None
        matches.append((p, t, int(padded[i, j])))

    micro_num_p = micro_den_p = micro_num_r = micro_den_r = 0
    macro_p_list, macro_r_list = [], []
    zero_overlap_true = []
    for p, t, ov in matches:
        if t is None:
            continue  # a predicted family matched to a dummy (more predicted than true families)
        tsize = len(true_fam[t])
        psize = len(pred_fam[p]) if p is not None else 0
        micro_num_r += ov
        micro_den_r += tsize
        if p is not None:
            micro_num_p += ov
            micro_den_p += psize
            macro_p_list.append(ov / psize if psize else 0.0)
        else:
            macro_p_list.append(0.0)
        macro_r_list.append(ov / tsize if tsize else 0.0)
        if ov == 0:
            zero_overlap_true.append((t, tsize, p))

    micro_prec = micro_num_p / micro_den_p if micro_den_p else 0.0
    micro_rec = micro_num_r / micro_den_r if micro_den_r else 0.0
    macro_prec = sum(macro_p_list) / len(macro_p_list) if macro_p_list else 0.0
    macro_rec = sum(macro_r_list) / len(macro_r_list) if macro_r_list else 0.0

    print()
    print("=== bipartite-matched (Hungarian algorithm, max total overlap) ===")
    print(f"MICRO (gene-weighted)  precision={micro_prec:.3f}  recall/sensitivity={micro_rec:.3f}")
    print(f"MACRO (family-weighted) precision={macro_prec:.3f}  recall/sensitivity={macro_rec:.3f}")
    print(f"true families with ZERO overlap in their matched predicted partner "
          f"(structurally undetected): {len(zero_overlap_true)}/{n_true} "
          f"({100*len(zero_overlap_true)/n_true:.1f}%)")

    if show_worst and zero_overlap_true:
        by_size = sorted(zero_overlap_true, key=lambda x: -x[1])[: show_worst]
        print(f"\nlargest undetected true families (up to {show_worst}):")
        for t, size, p in by_size:
            print(f"  {t}: {size} members, matched to predicted family {p!r} (0 shared genes)")
    return dict(n_true=n_true, n_pred=n_pred, micro_precision=micro_prec, micro_recall=micro_rec,
                macro_precision=macro_prec, macro_recall=macro_rec, n_undetected=len(zero_overlap_true))


# ---------------------------------------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------------------------------------

def cmd_genesets(a):
    """Derive the two genesets the chain reads from S1C. `--out-full`: every S1C Gene ID (2,334).
    `--out-eligible`: S1C `In Table S1 (SD98 gene set)` = Yes (1,793). Both 2-column (`gene_id`, `biotype`),
    sorted by gene_id. Before wave 7 these were hand-made files under winloci_data/soto_replication/
    (soto_2334_geneset.tsv, soto_1793_geneset.tsv) with no generator; `--out-full` reproduces
    soto_2334_geneset.tsv byte for byte, and the eligible file has the same (gene_id, biotype) set as
    soto_1793_geneset.tsv (which carried two more columns, in an unrecorded order) and drives `cluster`
    byte-identically."""
    full, eligible = soto_genesets(a.truth)
    if a.out_full:
        write_geneset(a.out_full, full)
    if a.out_eligible:
        write_geneset(a.out_eligible, eligible)
    print(f"[genesets] {len(full)} genes in the S1C universe, {len(eligible)} family-eligible (In Table S1 = Yes)"
          + (f" -> {a.out_full}" if a.out_full else "")
          + (f", {a.out_eligible}" if a.out_eligible else ""), file=sys.stderr)


def cmd_curate(a):
    """Recompute Soto's SD98 gene set from SD98 regions + CAT v4 and apply an explicit, automatic version of
    their manual curation step (2026-09-28, KEY=rule; report figs/soto_cur_rule.md). OPT-IN: nothing else in
    this module calls it, and the headline chain (`genesets` from S1C) is byte-identical without it.

    Regions: `--sd98-bed` (already-merged SD98 regions, e.g. the UCSC v1.0 sedefSegDups track at $24 >= 0.98
    -- their input, 97.8 Mbp) or `--sedef` (a native SEDEF table: rows with field-21 identity >= --min-identity,
    both sides, merged -- the same file `edges --native-v1` reads, so the chain is fully native to it).
    Genes: >= 1 exon fully contained in a region (`-f 1`), chrY/chrM never, chrX dropped from the output (their
    gene-count step); eligible = biotype grep "protein_coding|unprocessed_pseudogene". With the UCSC track this
    gives 5,154 genes / 1,864 eligible = Soto's 1,793 + the 71 they removed by hand.

    Curation (`--rule`, see curation_drop): the candidate partners are other eligible SD98 genes on the same
    chromosome and strand that share exonic bp (their transcript self-intersect universe, chrX included at this
    step as in their code). `--rule none` writes the uncurated set. The default rule is the one that survived
    the held-out test in the report; `--contain` is its exonic-containment threshold.

    Outputs: `--out-eligible` (curated family-eligible genes, `gene_id biotype`, sorted); `--out-full` the
    universe for `edges --geneset` / `cluster --full-geneset`: `--full-source s1c` (default) = S1C's 2,334 genes
    UNION the curated eligible set, so the only thing that changes vs the headline is the eligible backbone;
    `--full-source sd98` = every SD98 gene of any biotype (their 5,154-style set) minus the dropped ones.
    `--out-dropped` lists the removed genes with the feature that fired; `--out-features` dumps the per-gene
    vocabulary for every eligible gene (what scratchpad/soto_cur/rule_search.py consumed).
    """
    genes = load_cat_genes(a.cat_bed)
    if a.sedef:
        regions = sd98_regions_from_sedef(a.sedef, a.min_identity)
    else:
        regions = load_regions_bed(a.sd98_bed)
    bp = sum(e - s for c, v in regions.items() if c not in ("chrX", "chrY", "chrM") for s, e in v)
    print(f"[regions] {sum(len(v) for v in regions.values())} merged SD98 regions, {bp:,} autosomal bp",
          file=sys.stderr)
    sd_auto = sd98_gene_set(genes, regions, autosomal=True)
    sd_x = sd98_gene_set(genes, regions, autosomal=False)
    elig_x = {g for g in sd_x if is_eligible_biotype(genes[g]["biotype"])}
    elig_auto = {g for g in elig_x if g in sd_auto}
    print(f"[genes] {len(sd_auto)} autosomal SD98 genes (all biotypes), {len(elig_auto)} family-eligible "
          f"({len(elig_x) - len(elig_auto)} more on chrX, used only as curation partners)", file=sys.stderr)

    feats = curation_features(genes, elig_x)
    drop = curation_drop(feats, a.rule, a.contain) & elig_auto
    n_univ = sum(1 for g in elig_auto if feats[g]["n_partners"] >= 1)
    print(f"[curate] rule={a.rule} contain={a.contain}: {n_univ} eligible genes have a same-strand exon-sharing "
          f"partner (curation universe); {len(drop)} dropped -> {len(elig_auto) - len(drop)} eligible",
          file=sys.stderr)

    eligible_rows = sorted((g, genes[g]["biotype"]) for g in elig_auto - drop)
    if a.out_eligible:
        write_geneset(a.out_eligible, eligible_rows)
    if a.out_full:
        if a.full_source == "s1c":
            full = dict(soto_genesets(a.truth)[0])
            for g, b in eligible_rows:
                full.setdefault(g, b)
        else:
            full = {g: genes[g]["biotype"] for g in sd_auto if g not in drop}
        write_geneset(a.out_full, sorted(full.items()))
        print(f"[full] {len(full)} genes ({a.full_source}) -> {a.out_full}", file=sys.stderr)
    if a.out_dropped:
        with open(a.out_dropped, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene_id", "gene_name", "biotype", "name_class", "n_tx", "flag", "exon_contain_frac",
                        "best_partner", "best_partner_name"])
            for g in sorted(drop):
                f = feats[g]
                best = max(f["partners"], key=lambda p: p[1]) if f["partners"] else None
                w.writerow([g, genes[g]["name"], f["biotype"], f["name_class"], f["n_tx"], int(f["flag"]),
                            f"{f['exon_contain_frac']:.4f}", best[0] if best else "",
                            genes[best[0]]["name"] if best else ""])
    if a.out_features:
        cols = ["gene_id", "gene_name", "chrom", "strand", "biotype", "name_class", "source", "n_tx", "n_exons",
                "span", "has_cds", "flag", "n_flag_by", "body_contained", "n_partners", "n_covered",
                "exon_contain_frac", "exon_novel_frac", "cds_novel_frac", "partners"]
        with open(a.out_features, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for g in sorted(elig_x):
                f = feats[g]
                w.writerow([g, genes[g]["name"], genes[g]["chrom"], genes[g]["strand"], f["biotype"],
                            f["name_class"], f["source"], f["n_tx"], f["n_exons"], f["span"], int(f["has_cds"]),
                            int(f["flag"]), len(f["flag_by"]), int(f["body_contained"]), f["n_partners"],
                            f["n_covered"], f"{f['exon_contain_frac']:.4f}", f"{f['exon_novel_frac']:.4f}",
                            "" if f["cds_novel_frac"] is None else f"{f['cds_novel_frac']:.4f}",
                            ";".join(f"{h}:{fr:.4f}:{int(c)}:{int(h in f['flag_by'])}"
                                     for h, fr, c in sorted(f["partners"]))])
    print(f"[done] {len(eligible_rows)} eligible genes" + (f" -> {a.out_eligible}" if a.out_eligible else ""),
          file=sys.stderr)


def cmd_edges(a):
    """Redo Soto's SD98/shared-exon/famCN replication using a fresh, unmerged, native CHM13 v2.0 SEDEF
    output, in place of the merged UCSC SD98 track that caused the merge-artifact bug documented in the
    archived tile_sd98_regions.py: bedtools-merging ~11k SD units into 817 blocks (mean 119.7 kb, max 4.25 Mb)
    made each block's own perfect self-alignment the minimap2 primary, so `-p 0.5` discarded every true
    paralog hit and made 69 entire Soto families structurally invisible. Tiling (20 kb / 10 kb windows) was
    a proxy fix; this step uses the real thing: "Soto's unmerged SD98 unit BED", per that file's own
    closing note.

    DEVIATION FROM SOTO'S LITERAL RECIPE, DISCLOSED. Their step 3 is "extract SD98 region FASTA, map back
    to the genome with minimap2". This step skips that second alignment pass entirely and instead uses
    SEDEF's OWN already-computed pairwise CIGAR (one row = one real duplication call between two specific
    regions) to project exons from one side to the other -- the same "project through the CIGAR" principle
    the retired map-back path (soto_replicate_clustering.py, tag notebook-2026-09-24) used, just fed by
    SEDEF's alignment instead of a second minimap2 run. This avoids re-introducing ANY risk of the exact
    `-p 0.5` bug above, since minimap2 is not invoked here at all. Report this as a deliberate, disclosed
    choice.

    INPUT FORMAT (verified empirically against this file, not assumed from generic PAF/SAM convention):
    34-column native SEDEF output, identical schema to the gorilla side's GGO_sedef_final.bed.
      field(1-idx)  1    2      3    4      5      6    7  8       9        ...  21        23           33
      content       chr1 start1 end1 chrom2 start2 end2 sc strand1 strand2  ...  identity  divergence   CIGAR
    strand1 is always '+' (checked: 36143/36143). CIGAR semantics confirmed by direct arithmetic against
    this file's own field 11/12 (M+D total == end1-start1; M+I total == end2-start2): **D consumes side1,
    I consumes side2** -- the reverse of typical minimap2/PAF convention (D=target-only, I=query-only) --
    so this file's CIGAR must NOT be walked with the usual SAM assumption.

    COORDINATES. Input is CHM13 v2.0 (chrY present; Soto's own v1.0 SD track has none). Lifted to v1.0
    using the SAME per-chromosome constant-offset table already built and validated for famCN
    (build_liftover, fitted from soto_parCN_S1E.tsv's dual v1.0/v2.0 coordinates) -- reused verbatim, not
    re-derived. A pair is dropped (and counted) if EITHER side fails to lift (straddles a regime switch,
    or its chromosome has no anchors) -- never guessed, matching the famcn step's own stated policy.

    --native-v1 (2026-09-28): the SEDEF input is ALREADY CHM13 v1.0 (the user's `final_v1.bed`, same 34-column
    schema and CIGAR convention, verified: 3,628/3,628 rows >= 0.98 reconcile), so the v2.0->v1.0 lift is the
    identity and S1E/--extra-anchors are not read. Every row keeps both sides; offsets are 0. Omit for the
    original behaviour (byte-identical to before this flag existed, checked against the frozen
    shared_exons_2334_finalhuman.tsv).
    """
    if a.native_v1:
        lift_fn = lambda chrom, start, end: (start, end)
        print("[liftover] --native-v1: input is CHM13 v1.0 already, identity map (S1E not read)", file=sys.stderr)
    else:
        extra_anchors = load_extra_anchors(a.extra_anchors)
        table, spans = build_liftover(a.s1e, extra_anchors=extra_anchors)
        lift_fn = lambda chrom, start, end: lift(table, spans, chrom, start, end)
        print(f"[liftover] {len(table)} chromosomes with anchors", file=sys.stderr)

    genes, _biotype = load_geneset(a.geneset)
    exons, _meta = load_exons(a.cat_bed, genes)
    per_chrom = defaultdict(list)
    for g, evs in exons.items():
        for c, s, e in evs:
            per_chrom[c].append((s, e, g))
    for c in per_chrom:
        per_chrom[c].sort()
    print(f"[genes] {len(genes)} SD98 genes, {sum(len(v) for v in exons.values())} exons", file=sys.stderr)

    n_rows = n_ident = n_lift_ok = n_cigar_bad = 0
    n_proj_total = 0
    shared = defaultdict(set)

    with open(a.sedef) as fh:
        for line in fh:
            n_rows += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 33:
                continue
            chrom1, start1, end1 = f[0], int(f[1]), int(f[2])
            chrom2, start2, end2 = f[3], int(f[4]), int(f[5])
            strand2 = f[9]
            try:
                identity = float(f[20])
            except ValueError:
                continue
            if identity < a.min_identity:
                continue
            n_ident += 1

            lift1 = lift_fn(chrom1, start1, end1)
            lift2 = lift_fn(chrom2, start2, end2)
            if lift1 is None or lift2 is None:
                continue
            v1_start1, v1_end1 = lift1
            v1_start2, v1_end2 = lift2
            offset1 = v1_start1 - start1
            offset2 = v1_start2 - start2
            n_lift_ok += 1

            blocks, p1_end, bw_end = build_blocks(start1, cigar_ops(f[32]))
            if p1_end != end1 or bw_end != (end2 - start2):
                n_cigar_bad += 1
                continue

            n_proj_total += find_shared_exons(
                blocks, start2, end2, strand2,
                own_chrom=chrom1, own_lo_v1=v1_start1, own_hi_v1=v1_end1, own_offset=offset1,
                own_axis_is_bwalk=False,
                other_chrom=chrom2, other_lo_v1=v1_start2, other_hi_v1=v1_end2, other_offset=offset2,
                other_axis_is_bwalk=True,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )
            n_proj_total += find_shared_exons(
                blocks, start2, end2, strand2,
                own_chrom=chrom2, own_lo_v1=v1_start2, own_hi_v1=v1_end2, own_offset=offset2,
                own_axis_is_bwalk=True,
                other_chrom=chrom1, other_lo_v1=v1_start1, other_hi_v1=v1_end1, other_offset=offset1,
                other_axis_is_bwalk=False,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )

            if a.limit and n_ident >= a.limit:
                break

    print(f"[rows] {n_rows} total, {n_ident} >= identity {a.min_identity}, "
          f"{n_lift_ok} both-sides lifted (dropped {n_ident - n_lift_ok}), "
          f"{n_cigar_bad} CIGAR-length mismatches rejected", file=sys.stderr)
    print(f"[link] {n_proj_total} exon projections -> {len(shared)} genes with >=1 shared exon",
          file=sys.stderr)

    # sorted, not dict/set iteration order: hash-order-dependent output is a reproducibility defect, not a
    # cosmetic one -- the underlying edge set was already confirmed identical run-to-run, but row ORDER
    # previously varied with Python's per-run string-hash randomization, which a byte-identity check would
    # wrongly flag.
    edge_pairs = sorted({tuple(sorted((g, p))) for g, partners in shared.items() for p in partners})
    with open(a.out_shared, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene_a", "gene_b"])
        w.writerows(edge_pairs)
    print(f"[done] {len(edge_pairs)} unique shared-exon edges -> {a.out_shared}", file=sys.stderr)


def cmd_cluster(a):
    """Steps 5-6 of Soto's clustering (connected components -> famCN MAD split -> family/singleton call)
    on an `edges` output. ONE implementation of the famCN-split rule, shared by the SEDEF path and (before
    wave 7) the minimap2 map-back path, instead of two copies that could drift.

    Add --full-geneset (+ build --shared over that same full universe, see `edges --geneset`) to ALSO admit
    non-family-eligible genes (lncRNA, processed_pseudogene, etc.) as MEMBERS of an already-formed family --
    Soto's own definition includes them ("SD98 genes associated with other gene features... were also
    assigned a gene family ID [when they join one]"), and the advisor's own framing is specifically
    "replicate Soto, whose definition includes pseudogenes/lncRNAs". Opt-in, off by default: the
    --full-geneset omitted case is BYTE-IDENTICAL to before this flag existed (docs/o1_ledger.md §6ih).

    This does NOT add those genes as full graph nodes able to bridge two families together -- that was
    tried and rejected (§6ih: precision 0.925->0.730, the same promiscuous-bridge-gene failure as §6ie's own
    mega-component bug). It attaches each extra gene to whichever already-formed family it shares an exon
    with (ties broken by edge count), via attach() -- the retired soto_attach_noncoding_members.py's
    function, which §6ii verified gives the EXACT same partition as that script (390 fams / 1,817 genes,
    same ARI). Its rationale: Soto's S1C places ~541 non-eligible genes (mostly lncRNA/processed_pseudogene)
    in a family, and scoring only over the 1,793 eligible genes is the shrinking-denominator metric trap
    retracted on 2026-08-02; letting those genes act as graph nodes instead collapses precision (§6ih).

    SINGLE-ELIGIBLE-SEED ISLANDS (docs/o1_ledger.md §6im), also gated by --full-geneset: a raw (full,
    unfiltered) connected component that touches EXACTLY ONE eligible-biotype gene cannot be a promiscuous
    bridge between two eligible families -- bridging requires the component to reach >=2 eligible genes, and
    if it did, it would (by definition of "connected component") already be ONE component containing both,
    not two separate ones. So such a component's full edges are admitted wholesale before clustering, with
    none of the precision risk §6ih measured (that risk is specifically a promiscuous gene reaching INTO a
    second eligible-anchored component). This recovers the §6ik `no_eligible_seed_isolated` case (3/5 of the
    largest fully-missed true families: a lone eligible gene with only non-eligible siblings, correctly
    edge-connected but stranded because the filtered backbone graph has no node for its partners).
    """
    mad_fn = mad_median if a.mad_statistic == "median" else mad_mean

    genes, biotype = load_geneset(a.geneset)

    # the BACKBONE clustering step must only ever see eligible-eligible edges -- an edge touching a
    # --full-geneset-only gene must NOT let step5_step6 treat that gene as a graph node (that is exactly
    # the naive, rejected approach: docs/o1_ledger.md §6ih measured it collapsing precision 0.925->0.730
    # by letting non-eligible genes bridge two components together). The unfiltered edge list is still
    # used as-is for the attach() call below, which needs the extra genes' own edges.
    edges = read_edges(a.shared)
    shared = adjacency(edges, keep=genes)

    full_biotype = {}
    if a.full_geneset:
        full_genes_all, full_biotype = load_geneset(a.full_geneset)
        biotype.update(full_biotype)

        full_adj = adjacency(edges, keep=full_genes_all)
        n_islands, n_island_genes = single_seed_islands(shared, full_adj, genes)
        print(f"[single-eligible-seed] {n_islands} raw component(s) ({n_island_genes} genes total) with "
              f"exactly one eligible gene -- admitted wholesale", file=sys.stderr)

    famcn = load_famcn(a.famcn)

    n_fam, n_single, n_comps = step5_step6(shared, genes, biotype, famcn, a.mad, a.out, mad_fn=mad_fn)
    print(f"[step4-input] {len(genes)} SD98 genes, {len(shared)} genes with >=1 shared exon, "
          f"{n_comps} raw components", file=sys.stderr)
    print(f"[done] {n_fam} families, {n_single} singleton/ineligible genes -> {a.out}", file=sys.stderr)

    if not a.full_geneset:
        return

    n_extra, n_attached = attach_extra_members(a.out, genes, full_genes_all, full_biotype, famcn, edges)
    print(f"[attach] {n_extra} extra genes from --full-geneset considered, "
          f"{n_attached} attached to an existing family -> {a.out} (rewritten)", file=sys.stderr)


def cmd_dennislab(a):
    """The ACTUAL Dennis-lab family-clustering algorithm, reverse-engineered from their own released code
    (github.com/mydennislab/HSD_brain_evolution, section_I&II/B_SD98_families.ipynb -- fetched and read
    directly, 2026-09-11), not from the paper's own natural-language METHODS description. Every prior
    reimplementation in this project (the retired soto_replicate_clustering.py, and step5_step6 behind
    `cluster`) modeled their step 5 as "build one shared-exon graph, take connected components, split any
    component whose famCN MAD is too high into smaller coherent sub-groups" -- a literal reading of
    "groupings where the MAD of CN was less than one were selected". Their real code does something
    structurally different, confirmed by reading it directly (verbatim quotes below), not by re-guessing
    from prose:

        get_mad(elements) = stats.median_abs_deviation(median(wssd rows for elements))   # MEDIAN, not MEAN
        low_dispersion_clusters  = [c for c in raw_clusters if get_mad(c) <  1]
        high_dispersion_clusters = [c for c in raw_clusters if get_mad(c) >= 1]           # DISCARDED outright,
                                                                                           # never split further
        genes = [protein_coding/unprocessed_pseudogene elements of any low_dispersion_cluster]
        for gene in genes:
            gene_cluster = [gene]; i = 0
            while True:
                if gene_cluster[i] is protein_coding/unprocessed_pseudogene:
                    for cluster in low_dispersion_clusters:
                        if gene_cluster[i] in cluster:
                            gene_cluster = list(set(gene_cluster) | cluster)   # merge the WHOLE cluster in
                i += 1
                if i == len(gene_cluster):
                    families.append(gene_cluster); break
        families = dedup(sorted(families))   # the same underlying group is found once per coding seed it holds

    So: (1) MAD filtering happens FIRST, as a pass/fail gate on each RAW shared-exon component, not as a
    post-hoc split of an over-large one -- a component that fails is dropped from family formation entirely,
    not partitioned into smaller MAD-coherent pieces (the "greedy sorted agglomeration" every prior
    reimplementation invented for "how do you split it" was answering a question their code never asks).
    (2) Two low-dispersion clusters that share ONLY a non-coding gene (lncRNA / processed pseudogene) are
    NEVER merged -- expansion propagates only through protein-coding/unprocessed-pseudogene bridge genes.
    This is a real, structural difference from plain connected-components over the whole shared-exon graph
    (which would merge on ANY shared gene, coding or not), independent of the MAD statistic question
    (mean vs median, see mad_mean/mad_median) and worth testing separately.

    Result (§6if Finding 2): ARI 0.541 (mean) / 0.564 (median) on the 1,793 universe -- worse than
    `cluster`; PARKED, not adopted. Reads --shared UNFILTERED (unlike `cluster`, which keeps only edges
    inside --geneset).
    """
    mad_fn = mad_median if a.mad_statistic == "median" else mad_mean

    genes, biotype = load_geneset(a.geneset)
    shared = load_shared(a.shared)
    famcn = load_famcn(a.famcn)

    comps = raw_components(shared)
    families, orphaned_low, orphaned_high = dennislab_families(comps, biotype, famcn, a.mad, mad_fn)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        assigned = set()
        for i, fam in enumerate(families):
            fid = f"DENNISFAM{i}"
            for g in sorted(fam):
                w.writerow([g, biotype.get(g, ""), fid, len(fam),
                            f"{famcn[g]:.2f}" if g in famcn else "", "family"])
                assigned.add(g)
        for g in sorted(genes - assigned):
            status = "low_dispersion_no_coding_bridge" if g in orphaned_low else \
                     "high_dispersion_discarded" if g in orphaned_high else "singleton"
            w.writerow([g, biotype.get(g, ""), "", 1, f"{famcn[g]:.2f}" if g in famcn else "", status])

    n_low = sum(1 for c in comps
                if len([famcn[g] for g in c if g in famcn]) < 2
                or mad_fn([famcn[g] for g in c if g in famcn]) < a.mad)
    print(f"[step4-input] {len(genes)} SD98 genes, {len(shared)} genes with >=1 shared exon, "
          f"{len(comps)} raw components ({n_low} low-dispersion, {len(comps) - n_low} high-dispersion)",
          file=sys.stderr)
    print(f"[done] {len(families)} families (mad-statistic={a.mad_statistic}), "
          f"{len(orphaned_low)} low-dispersion-but-no-coding-bridge, "
          f"{len(orphaned_high)} high-dispersion-discarded -> {a.out}", file=sys.stderr)


def cmd_famcn(a):
    """Compute famCN (WSSD read-depth copy number) at ARBITRARY coordinates, replicating Soto's CN leg.

    WHY THIS EXISTS. Soto define a gene family as shared-exon genes whose family copy number agrees
    (mean-absolute-deviation < 1), so famCN is their SPLITTER. Their supplementary table S1C publishes famCN,
    but only at THEIR annotated genes -- and our over-merged families contain predicted copies with no Soto
    gene, which is exactly where we need the number. This computes famCN ourselves from the per-sample WSSD
    tracks, so it can be evaluated at any interval we predict.

    VALIDATION (3 samples vs Soto's published famCN, n=269): SYCE1 2.19 vs 2.01, PMCHL1 3.93 vs 3.99,
    OR4K7P 5.64 vs 5.98, USP32P3 8.37 vs 8.60, CBWD6 13.67 vs 13.70, AC133919.3 24.80 vs 22.76,
    NPIPB15 46.67 vs 48.25 -- 7/8 within 9% over a 24x range. The 8th (BET1L, famCN 702, subtelomeric) needs
    far more samples; treat CN > ~100 as unresolved at low sample counts.

    COORDINATES. The WSSD tracks are T2T-CHM13 **v1.0**; we work in **v2.0**. Rather than a chain file, the
    liftover is derived from Soto's own table S1E, which carries BOTH `SD98_v1.0` and `SD98_v2.0` for 1833
    paralogs: per chromosome the difference is a CONSTANT (a few hundred bp), except the five acrocentrics,
    which have exactly one breakpoint each in the short arm that v2.0 rebuilt. That gives a piecewise-constant
    map fitted to real anchor pairs. Intervals landing in the uncertain zone between two regimes are reported
    as unmapped rather than guessed -- a wrong offset would silently return another locus's copy number, which
    is worse than an absent value.

    CIRCULARITY WARNING. famCN is CONSTITUTIVE of the Soto ground truth: they built those families with it. Any
    score that uses famCN to split and is then evaluated against Soto is partly circular and must say so.
    The honest uses are (a) characterising what a CN signal would buy ("CN would resolve N of our M residual
    over-merges"), and (b) splitting when evaluated on an INDEPENDENT truth set.
    """
    extra_anchors = load_extra_anchors(a.extra_anchors)
    table, spans = build_liftover(a.s1e, extra_anchors=extra_anchors)
    print(f"[liftover] fitted over {len(table)} chromosomes; "
          f"{sum(len(v) for v in spans.values())} regime switches guarded", file=sys.stderr)

    if a.wssd_dir:
        bbs = sorted(os.path.join(a.wssd_dir, f) for f in os.listdir(a.wssd_dir) if f.endswith("_wssd.bb"))
        if not bbs:
            sys.exit(f"no *_wssd.bb in {a.wssd_dir}")
    else:
        listing = subprocess.run(["curl", "-s", "--max-time", "60", BASE_URL + "/"],
                                 capture_output=True, text=True).stdout
        bbs = [f"{BASE_URL}/{m}" for m in re.findall(r'href="([^"]+_wssd\.bb)"', listing)]
    bbs = bbs[: a.samples]
    print(f"[wssd] {len(bbs)} sample track(s)", file=sys.stderr)

    rows = list(csv.DictReader(open(a.intervals), delimiter="\t"))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([a.id_col, "chrom", "start", "end", "v1_start", "v1_end",
                    "famCN_median", "famCN_mad", "n_samples", "status"])
        done = unmapped = 0
        for r in rows:
            if not r.get("chrom"):
                continue
            chrom, s, e = r["chrom"], int(r["start"]), int(r["end"])
            lifted = lift(table, spans, chrom, s, e)
            if lifted is None:
                w.writerow([r.get(a.id_col, ""), chrom, s, e, "", "", "", "", 0, "UNMAPPED_v1"])
                unmapped += 1
                continue
            v1s, v1e = lifted
            # One subprocess per (interval, sample); with the full 271-sample panel that is ~10.8k calls
            # for 40 intervals, so they run on a thread pool (each worker just waits on bigBedToBed).
            with ThreadPoolExecutor(max_workers=a.jobs) as ex:
                vals = [v for v in ex.map(lambda bb: cn_for_interval(bb, chrom, v1s, v1e, a.tool), bbs)
                        if v is not None]
            if not vals:
                w.writerow([r.get(a.id_col, ""), chrom, s, e, v1s, v1e, "", "", 0, "NO_CN"])
                continue
            med = st.median(vals)
            mad = st.median([abs(v - med) for v in vals])
            w.writerow([r.get(a.id_col, ""), chrom, s, e, v1s, v1e,
                        f"{med:.3f}", f"{mad:.3f}", len(vals), "OK"])
            done += 1
            if done % 50 == 0:
                print(f"  {done} intervals done", file=sys.stderr)
    print(f"[done] {done} with CN, {unmapped} unmappable -> {a.out}", file=sys.stderr)


def cmd_score(a):
    """Score a predicted gene->family_id assignment against Soto's own published truth (S1C).

    `--only pairs` (was soto_score_against_truth.py): the FIXED-UNIVERSE methodology this project settled on
    after retracting an earlier shrinking-set version of this exact metric (docs/o1_ledger.md,
    project_soto_full_replication.md 2026-08-02 section): every gene in the truth universe is scored,
    including ones the prediction never placed (given a unique singleton label rather than being silently
    dropped from the comparison). Prints ARI, exact family matches and pair precision/recall/F1.

    Genes with MORE THAN ONE distinct Family ID across their S1C rows (their own table has 149 such genes,
    mostly among the ~541 non-eligible-biotype "extra" genes outside their 1,793-gene family-eligible set --
    docs/o1_ledger.md §6ih) are EXCLUDED from scoring: a partition-comparison metric needs one ground-truth
    label per gene, and picking one of several real, simultaneously-true Family IDs for such a gene would be
    arbitrary. This is a disclosed exclusion (reported in the output), not a silent one.

    By default scores over ALL 2,185 genes with unambiguous ground truth (Soto's real family-membership
    universe, including non-eligible-biotype members) -- the fair, complete comparison. Pass
    --eligible-only-universe to instead reproduce the narrower 1,793-gene comparison this project's earlier
    sections used (valid as a DIFFERENT, more limited question -- "how well do we cluster the family-eligible
    genes" -- not a substitute for the complete one).

    `--only bipartite` (was soto_bipartite_match_score.py): OPTIMAL BIPARTITE MATCHING between predicted and
    true families (Hungarian algorithm, scipy.optimize.linear_sum_assignment), maximizing total gene overlap
    over a 1:1 assignment of predicted families to true families. SCOPE NOTE (standing project rule): this is
    a bipartite match used ONLY as an evaluation/scoring technique to compare two ALREADY-COMPUTED partitions
    after the fact. It plays no role in how families are defined or how genes are assigned to them anywhere
    in this pipeline -- the standing rule ("no bipartite matching or facility-location step" in the
    family-DEFINITION method) is about modeling choices, not about how a finished result gets graded against
    an external benchmark. Reports, per matched (predicted, true) family pair, PRECISION (overlap /
    |predicted family|) and RECALL (overlap / |true family|), aggregated MICRO (gene-weighted: sum over all
    matched pairs, dominated by large families) and MACRO (family-weighted mean), plus the true families
    that matched NOTHING (structurally undetected). Same load_truth() universe as `--only pairs`, so the
    numbers are directly comparable.

    `--only all` (default): the pairs output, then the bipartite output.
    """
    clean_truth, ambiguous = load_truth(a.truth)
    predicted = load_predicted(a.predicted)
    universe = restrict_universe(clean_truth, a.eligible_only_universe)

    if a.only in ("all", "pairs"):
        r = score(clean_truth, predicted, universe)
        print(f"scored universe: {r['n_genes']} genes ({len(ambiguous)} ambiguous multi-family genes "
              f"excluded from all scoring, {'restricted to --eligible-only-universe' if a.eligible_only_universe else 'full truth universe'})")
        print(f"ARI: {r['ari']:.4f}")
        print(f"exact family match: {r['n_exact']}/{r['n_truth_fam']} = {100*r['n_exact']/r['n_truth_fam']:.1f}%")
        print(f"pair precision/recall/F1: {r['precision']:.3f}/{r['recall']:.3f}/{r['f1']:.3f}")
        print(f"predicted families: {r['n_pred_fam']}")
    if a.only in ("all", "bipartite"):
        bipartite_score(clean_truth, ambiguous, predicted, universe, a.show_worst)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="soto_replication.py",
        description="Soto 2025 gene-family replication (CONCORDANCE with Soto, not independent: register "
                    "T15/858). Headline chain: genesets -> edges -> cluster -> score (ledger §6ie-§6ip).",
        epilog="Run `<subcommand> --help` for each step; the module docstring has the recipe, the §6ip numbers "
               "and the old-script -> subcommand map.")
    sub = ap.add_subparsers(dest="cmd", required=True, metavar="SUBCOMMAND")

    p = sub.add_parser("genesets", help="derive the 2,334-gene and 1,793-gene genesets from S1C (new)",
                       description=cmd_genesets.__doc__)
    p.add_argument("--truth", default=S1C, help="Soto Table S1C (default: %(default)s)")
    p.add_argument("--out-eligible", help="write the 1,793 family-eligible genes here (gene_id, biotype)")
    p.add_argument("--out-full", help="write the 2,334-gene S1C universe here (gene_id, biotype)")
    p.set_defaults(func=cmd_genesets)

    p = sub.add_parser("curate", help="OPT-IN (2026-09-28): recompute the SD98 gene set from SD98 regions + CAT v4 "
                                      "and apply an explicit version of Soto's manual curation (figs/soto_cur_rule.md)",
                       description=cmd_curate.__doc__)
    p.add_argument("--cat-bed", required=True, help="CAT v4 BED (CHM13 v1.0, 37 columns)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sd98-bed", help="merged SD98 region BED (their input: UCSC v1.0 sedefSegDups $24>=0.98)")
    src.add_argument("--sedef", help="native SEDEF table (field 21 identity); both sides of rows >= --min-identity, merged")
    p.add_argument("--min-identity", type=float, default=0.98, help="with --sedef (default %(default)s)")
    p.add_argument("--rule", choices=CURATION_RULES, default="contained_outranked",
                   help="curation rule (see curation_drop; `none` = uncurated). Default: %(default)s = the "
                        "one-condition rule that survived the held-out test (figs/soto_cur_rule.md)")
    p.add_argument("--contain", type=float, default=0.3,
                   help="exonic-bp containment threshold of the contained* rules (default %(default)s; the "
                        "held-out plateau is 0.30-0.35, a depth-1 tree picks 0.32)")
    p.add_argument("--truth", default=S1C, help="S1C, for --full-source s1c (default: %(default)s)")
    p.add_argument("--full-source", choices=["s1c", "sd98"], default="s1c",
                   help="--out-full = S1C's 2,334 genes UNION the curated eligible set (s1c, default) or every "
                        "SD98 gene of any biotype minus the dropped ones (sd98)")
    p.add_argument("--out-eligible", help="curated family-eligible geneset TSV (gene_id, biotype)")
    p.add_argument("--out-full", help="full-universe geneset TSV for `edges --geneset` / `cluster --full-geneset`")
    p.add_argument("--out-dropped", help="TSV of the removed genes")
    p.add_argument("--out-features", help="TSV of the per-gene curation vocabulary (every eligible gene)")
    p.set_defaults(func=cmd_curate)

    p = sub.add_parser("edges", help="steps 1-4: SEDEF CIGAR -> shared-exon edge TSV "
                                     "(was soto_replicate_from_sedef.py)")
    p.add_argument("--sedef", required=True, help="native CHM13 v2.0 SEDEF output (34 columns)")
    p.add_argument("--min-identity", type=float, default=0.98, help="SD98 floor (field 21, 1-indexed)")
    p.add_argument("--s1e", default=S1E, help="soto_parCN_S1E.tsv (builds the v2.0->v1.0 liftover; "
                                              "default: %(default)s)")
    p.add_argument("--geneset", required=True, help="geneset TSV (gene_id column; the headline uses the "
                                                    "2,334-gene `genesets --out-full` file)")
    p.add_argument("--cat-bed", required=True, help="CAT v4 BED (v1.0)")
    p.add_argument("--min-cov", type=float, default=0.99, help="bedtools -f equivalent")
    p.add_argument("--out-shared", required=True, help="TSV: gene_a<TAB>gene_b (one row per edge, deduped)")
    p.add_argument("--limit", type=int, default=0, help="stop after N qualifying rows (0 = no limit; for smoke tests)")
    p.add_argument("--extra-anchors",
                   help="OPT-IN: TSV (chrom, v2_pos, offset columns) of extra, individually-validated "
                        "liftover anchors -- e.g. from directly aligning a specific gene's own sequence "
                        "against both genome versions (docs/o1_ledger.md §6il) -- merged into the S1E "
                        "anchor table before fitting regimes. The headline passes "
                        "bench/soto/acro_extra_anchors.tsv. Omit for the original behaviour "
                        "(byte-identical to before this flag existed).")
    p.add_argument("--native-v1", action="store_true",
                   help="OPT-IN: --sedef is already CHM13 v1.0 (e.g. final_v1.bed): the v2.0->v1.0 lift is "
                        "the identity, --s1e/--extra-anchors are not read, no row is dropped for lifting. "
                        "Omit for the original v2.0 behaviour (byte-identical).")
    p.set_defaults(func=cmd_edges)

    p = sub.add_parser("cluster", help="steps 5-6: components -> famCN MAD split -> families "
                                       "(was soto_cluster_from_shared.py)")
    p.add_argument("--shared", required=True, help="gene_a<TAB>gene_b TSV (`edges` output)")
    p.add_argument("--geneset", required=True)
    p.add_argument("--famcn", required=True)
    p.add_argument("--mad", type=float, default=1.0)
    p.add_argument("--mad-statistic", choices=["mean", "median"], default="mean",
                   help="mean = the paper's own METHODS-text wording; median = what their released "
                        "code (B_SD98_families.ipynb) actually computes (scipy median_abs_deviation, "
                        "unscaled) -- default stays 'mean' so this flag is opt-in, not a silent change")
    p.add_argument("--full-geneset",
                   help="OPT-IN: also admit genes in this (larger) geneset that are NOT in --geneset "
                        "as MEMBERS of an already-formed family, via a shared exon -- never as a way to "
                        "found a family or merge two together. --shared must have been built over this "
                        "SAME full geneset (`edges --geneset <this file>`), or the "
                        "extra genes will have no edges to attach through. Omit for the original, "
                        "eligible-only behaviour (byte-identical to before this flag existed).")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_cluster)

    p = sub.add_parser("dennislab", help="the Dennis-lab notebook algorithm, PARKED arm "
                                         "(was soto_cluster_dennislab_algorithm.py)")
    p.add_argument("--shared", required=True)
    p.add_argument("--geneset", required=True)
    p.add_argument("--famcn", required=True)
    p.add_argument("--mad", type=float, default=1.0)
    p.add_argument("--mad-statistic", choices=["mean", "median"], default="median",
                   help="default 'median' here (unlike `cluster`'s default 'mean'), "
                        "since this step specifically exists to test their REAL algorithm faithfully")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_dennislab)

    p = sub.add_parser("famcn", help="WSSD read-depth famCN at arbitrary v2.0 intervals (was famcn_from_wssd.py)")
    p.add_argument("--intervals", required=True,
                   help="TSV with chrom/start/end (+ any id columns); v2.0 coordinates")
    p.add_argument("--wssd-dir", help="directory of *_wssd.bb; omit to stream from UCSC")
    p.add_argument("--samples", type=int, default=8, help="how many samples to median over")
    p.add_argument("--jobs", type=int, default=8, help="parallel bigBedToBed calls per interval")
    p.add_argument("--s1e", default=S1E, help="soto_parCN_S1E.tsv (default: %(default)s)")
    p.add_argument("--tool", default="bigBedToBed",
                   help="bigBedToBed (default, needs the UCSC binary) or 'pybigwig' (reads the same "
                        "bigBed files in-process via the pyBigWig package, no subprocess/binary needed)")
    p.add_argument("--id-col", default="family_id")
    p.add_argument("--extra-anchors",
                   help="OPT-IN: TSV (chrom, v2_pos, offset columns) merged into the S1E anchor table, "
                        "same mechanism/file format as `edges --extra-anchors` "
                        "(docs/o1_ledger.md §6il/§6in). Omit for the original behaviour.")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_famcn)

    p = sub.add_parser("score", help="ARI / exact / pair P-R-F1 and bipartite family matching vs S1C "
                                     "(was soto_score_against_truth.py + soto_bipartite_match_score.py)")
    p.add_argument("--predicted", required=True,
                   help="TSV with gene_id, family_id columns (empty family_id = unplaced)")
    p.add_argument("--truth", default=S1C, help="soto_famCN_S1C.tsv (default: %(default)s)")
    p.add_argument("--eligible-only-universe", help="restrict scoring to this geneset's gene_id column "
                   "(reproduces the narrower 1,793-gene comparison; omit for the full, honest universe)")
    p.add_argument("--only", choices=["all", "pairs", "bipartite"], default="all",
                   help="pairs = the old soto_score_against_truth.py stdout; bipartite = the old "
                        "soto_bipartite_match_score.py stdout; all (default) = both, in that order")
    p.add_argument("--show-worst", type=int, default=5,
                   help="bipartite: print this many worst-matched true families (by size) for inspection")
    p.set_defaults(func=cmd_score)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
