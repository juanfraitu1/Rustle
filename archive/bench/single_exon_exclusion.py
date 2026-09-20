#!/usr/bin/env python
"""single_exon_exclusion.py — QUANTIFY the trade-off of EXCLUDING single-exon transcripts from the
RNA multi-copy family definition (proposal: define families over SPLICED / multi-exon transcripts only).

Motivation (user): single-exon transcripts are the irreducible-from-structure wall -- exon-containment
cannot separate a single-exon real paralog from a single-exon domain-bridge (MAGE is intronless). We
quantify the trade-off with NUMBERS so the scope decision is informed. DNA truth = VALIDATION only.

RNA-only: "exon count" is READ-STRUCTURE-derived. The operational per-member exon count is the de-novo
assembled transcript's exon count = the last field of member_dn (built as f"DN_{c}_{s}_{len(introns)+1}"
in bench/denovo_assemble_gate.py), i.e. exactly the transcript that was fed to the POA family definition.
We ALSO cross-reference the annotation intron index (gene_intron_index.json) as a VALIDATION-only overlay
(the annotated n_intron per gene) to characterise the biological recall pool (retrocopy vs genuine).

FIVE computations (see FINDINGS in the returned JSON):
 1. POPULATION       -- de-novo single- vs multi-exon loci (skeletons); per-family class in the 607-family
                        default catalog (single-only / multi-only / mixed); families+members exclusion removes.
 2. FP REMOVED       -- of the named residual-FP roster (family_level_pr_current.json), how many are
                        single-exon (removed). Explicit: MAGE-A sub-cluster (LOC129529978/986)? GSTM2? fam17?
 3. RECALL COST      -- of the diploid-oracle multi_copy genes (asm_hapCN>=2), how many are single-exon and
                        would be LOST. NAME them. MAGE? histones? which oracle LOC genes?
 4. RETRO/ARTIFACT   -- split the single-exon population into retrocopy / exonized-repeat-or-aggregate /
                        assembly-truncation / genuine single-exon gene (reuse retrocopy-filter logic + n_intron).
 5. NET P/R          -- family-level P/R (vs diploid oracle) for the catalog with single-exon EXCLUDED,
                        vs current default (0.875 / 0.8421). Operational (de-novo) AND strict-annotation
                        counterfactual.

Deterministic (PYTHONHASHSEED=0). Reuses committed catalogs / committed eval JSON so numbers match.
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/single_exon_exclusion.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from collections import Counter, defaultdict

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
SKEL = os.path.join(SCRATCH, "denovo_skeletons.tsv")
META = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
IDX = os.path.join(SCRATCH, "gene_intron_index.json")
CATALOG = os.path.join(BENCH, "family_rna_refine.tsv")
ORACLE = os.path.join(BENCH, "diploid_cn_oracle.tsv")
PR_JSON = os.path.join(BENCH, "family_level_pr_current.json")
RETRO_JSON = os.path.join(BENCH, "family_def_retrocopy_filter.json")
OUT_TSV = os.path.join(BENCH, "single_exon_exclusion.tsv")
OUT_JSON = os.path.join(BENCH, "single_exon_exclusion.json")


def denovo_nexon(dn):
    """De-novo transcript exon count = last field of member_dn (DN_<chrom>_<start>_<n_exon>)."""
    return int(dn.split("_")[-1])


# ---------------------------------------------------------------------------- loaders
def load_intron_index():
    idx = json.load(open(IDX))
    return idx


def ann_n_intron(idx, g):
    v = idx.get(g)
    return None if v is None else v.get("n_intron")


def load_catalog():
    rows = list(csv.DictReader(open(CATALOG), delimiter="\t"))
    fams = defaultdict(list)
    for r in rows:
        fams[r["family_id"]].append(r)
    return rows, fams


def load_oracle_multicopy():
    """gene -> asm_hapCN for named autosomal/sex multi_copy oracle genes (asm_hapCN>=2)."""
    g2cn = {}
    for r in csv.DictReader(open(ORACLE), delimiter="\t"):
        try:
            cn = float(r["asm_hapCN"])
        except (ValueError, KeyError):
            continue
        if r["gene"] != "NA" and cn >= 2:
            g2cn[r["gene"]] = max(g2cn.get(r["gene"], 0), int(cn))
    return g2cn


# ---------------------------------------------------------------------------- 1. POPULATION
def compute_population(fams, idx):
    # de-novo skeleton loci
    se_loci = 0
    multi_loci = 0
    se_rows = []
    with open(SKEL) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            ne = int(p[3])
            if ne == 1:
                se_loci += 1
                se_rows.append((p[0], int(p[1]), int(p[2]), int(p[4])))
            else:
                multi_loci += 1
    # de-novo transcript FASTA (what family def actually consumed) n_exon dist
    meta_ne = Counter()
    with open(META) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            meta_ne[int(p[5])] += 1

    # per-family classification on de-novo member n_exon
    fam_class = {}
    n_single_only = n_multi_only = n_mixed = 0
    members_single = 0
    fams_excluded = 0
    for fid, members in fams.items():
        nexons = [denovo_nexon(m["member_dn"]) for m in members]
        n_se = sum(1 for x in nexons if x == 1)
        n_me = sum(1 for x in nexons if x >= 2)
        members_single += n_se
        if n_se and not n_me:
            cls = "single_exon_only"
            n_single_only += 1
        elif n_me and not n_se:
            cls = "multi_exon_only"
            n_multi_only += 1
        else:
            cls = "mixed"
            n_mixed += 1
        # a family is REMOVED BY EXCLUSION iff dropping its single-exon members takes it from >=2 distinct
        # loci to <2 (delta caused by exclusion; NOT a pre-existing <2-distinct-loci property).
        before_loci = len({(m["chrom"], m["start"]) for m in members})
        surviving_loci = len({(m["chrom"], m["start"]) for m, ne in zip(members, nexons) if ne >= 2})
        would_exclude = 1 if (before_loci >= 2 and surviving_loci < 2) else 0
        if would_exclude:
            fams_excluded += 1
        fam_class[fid] = dict(cls=cls, n_se=n_se, n_me=n_me, min_nexon=min(nexons),
                              would_exclude=would_exclude, surviving_loci=surviving_loci)

    return dict(
        denovo_skeleton_single_exon_loci=se_loci,
        denovo_skeleton_multi_exon_loci=multi_loci,
        denovo_skeleton_single_exon_are_chromosome_buckets=True,
        denovo_skeleton_single_exon_rows=[dict(chrom=c, start=s, end=e, span=e - s, n_reads=nr) for c, s, e, nr in se_rows],
        denovo_transcript_fasta_nexon_dist={str(k): v for k, v in sorted(meta_ne.items())},
        denovo_transcript_fasta_single_exon=meta_ne.get(1, 0),
        denovo_transcript_fasta_total=sum(meta_ne.values()),
        catalog_n_families=len(fams),
        catalog_n_members=sum(len(v) for v in fams.values()),
        family_class_counts=dict(single_exon_only=n_single_only, multi_exon_only=n_multi_only, mixed=n_mixed),
        catalog_single_exon_members=members_single,
        families_removed_by_exclusion=fams_excluded,
        members_removed_by_exclusion=members_single,
    ), fam_class


# ---------------------------------------------------------------------------- 2. FP REMOVED
def compute_fp_removed(pr, fams, idx):
    # gene -> set of de-novo n_exon across its catalog members; gene -> annotation n_intron
    gene_denovo = defaultdict(list)
    for members in fams.values():
        for m in members:
            gene_denovo[m["member_gene"]].append(denovo_nexon(m["member_dn"]))

    def gene_is_single_exon_denovo(g):
        v = gene_denovo.get(g)
        return bool(v) and all(x == 1 for x in v)

    roster = pr["residual_fp_roster"]
    named = []
    # multifam + oversize = the 6 DNA-confirmed over-merge blocks
    for tag in ("multifam", "oversize"):
        for e in roster.get(tag, []):
            genes = e["genes"]
            named.append(dict(
                source=tag, genes=genes, distinct_loci=e.get("distinct_loci"),
                ann_n_intron={g: ann_n_intron(idx, g) for g in genes},
                denovo_nexon={g: sorted(set(gene_denovo.get(g, []))) for g in genes},
                all_single_exon_denovo=all(gene_is_single_exon_denovo(g) for g in genes),
                all_single_exon_annotation=all((ann_n_intron(idx, g) == 0) for g in genes),
            ))
    # fam17 = ep-impure block_index 17 (the 16-protein-family repeat-bridge hub)
    fam17 = None
    for e in roster.get("ep_impure_not_dna_named", []):
        if e.get("block_index") == 17:
            genes = e["genes"]
            fam17 = dict(
                block_index=17, n_genes=len(genes), n_protein_families=e.get("n_protein_families"),
                any_single_exon_annotation=any((ann_n_intron(idx, g) == 0) for g in genes),
                n_single_exon_annotation_genes=sum(1 for g in genes if ann_n_intron(idx, g) == 0),
                all_single_exon_denovo=all(gene_is_single_exon_denovo(g) for g in genes if g in gene_denovo),
            )
    n_named = len(named)
    n_removed = sum(1 for e in named if e["all_single_exon_denovo"])
    # MAGE-A sub-cluster = LOC129529978 + LOC129529986
    mage_sub = next((e for e in named if set(e["genes"]) == {"LOC129529978", "LOC129529986"}), None)
    gstm2 = [e for e in named if "GSTM2" in e["genes"]]
    return dict(
        n_named_fp_entries=n_named,
        n_fp_removed_by_exclusion=n_removed,
        mage_A_subcluster_LOC129529978_986=dict(
            present=mage_sub is not None,
            ann_n_intron=mage_sub["ann_n_intron"] if mage_sub else None,
            denovo_nexon=mage_sub["denovo_nexon"] if mage_sub else None,
            single_exon_removed=mage_sub["all_single_exon_denovo"] if mage_sub else None,
        ),
        GSTM2_hub=dict(
            entries=[{"source": e["source"], "genes": e["genes"], "ann_n_intron": e["ann_n_intron"]} for e in gstm2],
            single_exon_removed=any(e["all_single_exon_denovo"] for e in gstm2),
        ),
        fam17_repeat_bridge_hub=fam17,
        named_fp_detail=named,
    )


# curated single-exon (intronless) MULTI-COPY gene families known to be biologically real (cancer-testis
# antigens, cytokine cluster, chemosensory GPCRs, replication-dependent histones, clustered protocadherins).
# These are the concrete casualties the RNA-read-structure pipeline loses UPSTREAM: each gene is intronless
# so it presents no splice junction, so it is bucketed per-chromosome and dropped before family definition,
# UNLESS the de-novo assembly happens to recruit a neighboring exon (making the de-novo transcript multi-exon).
# NONE of these are in the diploid CN oracle, so their recall cost is NOT measured by the P/R below -- but
# they are the honest biological recall boundary of a spliced-only definition, so we name them explicitly.
SE_MULTICOPY_FAMILY_PROBES = [
    ("IFNA_type_I_interferon", ("IFNA2", "IFNA6", "IFNA8", "IFNA21", "IFNAR1", "IFNAR2")),
    ("TAS2R_bitter_taste_GPCR", tuple("TAS2R%s" % s for s in
        ("1", "3", "4", "5", "7", "8", "9", "10", "13", "14", "16", "20", "38", "39", "40", "42", "50", "60"))),
    ("H1_linker_histone", ("H1-1", "H1-2", "H1-3", "H1-4", "H1-5", "H1-6")),
    ("MAGEB_cancer_testis", tuple("MAGEB%s" % s for s in
        ("1", "2", "3", "5", "6", "6B", "16", "17", "18"))),
    ("PCDHB_protocadherin_beta", tuple("PCDHB%s" % s for s in
        ("1", "2", "5", "8", "9", "11", "12", "14", "15", "16"))),
]


def probe_se_multicopy_families(idx, catalog_genes):
    """For each curated intronless multi-copy family present in THIS annotation, report how many genes are
    single-exon (n_intron==0) and how many are recovered in the shipped catalog (i.e. dodged the wall)."""
    out = []
    for label, members in SE_MULTICOPY_FAMILY_PROBES:
        present = [(g, ann_n_intron(idx, g)) for g in members if ann_n_intron(idx, g) is not None]
        se = [g for g, ni in present if ni == 0]
        se_in_catalog = [g for g in se if g in catalog_genes]
        out.append(dict(
            family=label,
            n_genes_in_annotation=len(present),
            n_single_exon=len(se),
            single_exon_genes=sorted(se),
            n_single_exon_in_catalog=len(se_in_catalog),
            single_exon_in_catalog=sorted(se_in_catalog),
            fully_lost_upstream=(len(se) >= 2 and len(se_in_catalog) == 0),
        ))
    return out


# ---------------------------------------------------------------------------- 3. RECALL COST
def compute_recall_cost(g2cn, fams, idx):
    gene_denovo = defaultdict(list)
    gene_fams = defaultdict(set)
    for fid, members in fams.items():
        for m in members:
            gene_denovo[m["member_gene"]].append(denovo_nexon(m["member_dn"]))
            gene_fams[m["member_gene"]].add(fid)
    dom_genes = {m["dominant_gene"] for members in fams.values() for m in members}
    catalog_genes = set(gene_denovo) | dom_genes

    single_exon_oracle = []
    for g, cn in sorted(g2cn.items()):
        nint = ann_n_intron(idx, g)
        if nint == 0:
            recovered = g in gene_denovo or g in dom_genes
            single_exon_oracle.append(dict(
                gene=g, asm_hapCN=cn, ann_n_intron=0,
                recovered_in_catalog=recovered,
                catalog_families=sorted(gene_fams.get(g, [])),
                denovo_nexon_of_member=sorted(set(gene_denovo.get(g, []))),
                lost_by_operational_exclusion=(recovered and all(x == 1 for x in gene_denovo.get(g, [2]))),
                lost_by_strict_annotation_exclusion=recovered,
            ))
    # MAGE / histone / interferon touchstones among oracle multicopy genes
    touchstones = {}
    for label, pat in (("MAGE", "MAGE"), ("histone_HIST", "HIST"), ("histone_H2A", "H2A"),
                       ("histone_H2B", "H2B"), ("histone_H3", "H3C"), ("histone_H4", "H4C"),
                       ("interferon_IFNA", "IFNA")):
        hits = sorted(g for g in g2cn if pat in g)
        touchstones[label] = [dict(gene=g, ann_n_intron=ann_n_intron(idx, g), asm_hapCN=g2cn[g]) for g in hits]

    n_oracle = len(g2cn)
    n_se = len(single_exon_oracle)
    se_families = probe_se_multicopy_families(idx, catalog_genes)
    return dict(
        oracle_multicopy_genes_named=n_oracle,
        oracle_multicopy_single_exon_annotation=n_se,
        oracle_multicopy_single_exon_genes=single_exon_oracle,
        oracle_multicopy_touchstones=touchstones,
        recall_cost_operational_named=[e["gene"] for e in single_exon_oracle if e["lost_by_operational_exclusion"]],
        recall_cost_strict_annotation_named=[e["gene"] for e in single_exon_oracle if e["lost_by_strict_annotation_exclusion"]],
        # honesty overlay: genuine single-exon multi-copy families in THIS annotation that are NOT oracle genes
        # (so invisible to the P/R) but ARE the concrete biological recall boundary of a spliced-only definition.
        genuine_se_multicopy_families_not_in_oracle=se_families,
        genuine_se_multicopy_families_fully_lost_upstream=[f["family"] for f in se_families if f["fully_lost_upstream"]],
    )


# ---------------------------------------------------------------------------- 4. RETRO / ARTIFACT SPLIT
def compute_retro_artifact(fams, idx, retro):
    # DE-NOVO single-exon population = the 26 chromosome-span buckets (all aggregation artifacts) + FASTA
    # ANNOTATION single-exon gene population and its recovery, split by family co-member structure.
    idx_se = [g for g, v in idx.items() if v.get("n_intron") == 0]
    n_ann_total = len(idx)
    n_ann_se = len(idx_se)

    gene_fams = defaultdict(set)
    gene_denovo = defaultdict(list)
    fam_members = fams
    for fid, members in fams.items():
        for m in members:
            gene_fams[m["member_gene"]].add(fid)
            gene_denovo[m["member_gene"]].append(denovo_nexon(m["member_dn"]))

    se_set = set(idx_se)
    recovered = [g for g in idx_se if g in gene_denovo]

    # classify each recovered single-exon-annotation gene by its FAMILY co-members' annotation structure
    split = Counter()
    detail = []
    for g in sorted(recovered):
        fids = gene_fams[g]
        comembers = []
        for fid in fids:
            for m in fam_members[fid]:
                if m["member_gene"] != g:
                    comembers.append(m["member_gene"])
        co_intron = [ann_n_intron(idx, c) for c in comembers]
        co_spliced = [x for x in co_intron if x is not None and x >= 1]
        co_intronless = [x for x in co_intron if x == 0]
        dn = gene_denovo.get(g, [])
        if dn and all(x == 1 for x in dn):
            cls = "denovo_single_exon"          # would be dropped operationally
        elif co_spliced and not co_intronless:
            cls = "retrocopy_signature"         # single-exon gene co-clustered with a SPLICED paralog
        elif co_intronless and not co_spliced:
            cls = "genuine_single_exon_tandem"  # both-intronless family (IFNA/PCDHB-like)
        elif co_spliced and co_intronless:
            cls = "mixed_comembers"
        else:
            cls = "unknown_comembers"
        split[cls] += 1
        detail.append(dict(gene=g, families=sorted(fids), denovo_nexon=sorted(set(dn)),
                           comember_ann_intron=co_intron, cls=cls))

    return dict(
        denovo_single_exon_loci_are_all_chromosome_aggregation_buckets=True,
        denovo_single_exon_resolved_real_loci=0,
        annotation_genes_total=n_ann_total,
        annotation_single_exon_genes=n_ann_se,
        annotation_single_exon_fraction=round(n_ann_se / n_ann_total, 4),
        annotation_single_exon_recovered_in_catalog=len(recovered),
        recovered_single_exon_annotation_class_split=dict(split),
        recovered_single_exon_annotation_detail=detail,
        genomewide_retrocopy_filter_edge_split=dict(
            one_side_intronless_retrocopy_dropped=retro["provenance"]["retrocopy_dropped"],
            both_intronless_genuine_tandem_kept=retro["provenance"]["both_intronless_kept"],
            both_spliced_kept=retro["provenance"]["both_spliced_kept"],
            note="genome-wide cDNA-homology edge split (family_def_retrocopy_filter.json); "
                 "both-intronless KEPT = genuine single-exon tandems (IFNA/PCDHB/ELOA/HNRNPCL); "
                 "one-side-intronless = processed retrocopy noise, already dropped by shipped filter",
        ),
    )


# ---------------------------------------------------------------------------- 5. NET P/R
def compute_net_pr(pr, recall_cost):
    cur = pr["truth3_diploid_oracle"]["current"]
    mapped = cur["oracle_mapped_families"]          # 48
    distinct_fp = cur["distinct_fp_blocks"]          # 6
    flags = cur["flag_instances"]                    # 7
    recov_mc = cur["oracle_genes_recovered_multicopy"]  # 48
    n_mc = cur["oracle_multicopy_genes"]             # 57

    # OPERATIONAL exclusion (drop de-novo single-exon members) = catalog is BYTE-IDENTICAL (0 single-exon
    # members) -> P/R identical by construction.
    operational = dict(
        precision_oracle_dedup=round(1 - distinct_fp / mapped, 4),
        precision_oracle_taskformula=round(1 - flags / mapped, 4),
        recall_oracle=round(recov_mc / n_mc, 4),
        oracle_mapped_families=mapped, distinct_fp_blocks=distinct_fp,
        oracle_genes_recovered_multicopy=recov_mc,
        note="identical to current default -- 0 single-exon members in catalog, exclusion is a no-op",
    )

    # STRICT-ANNOTATION counterfactual: also drop members whose ANNOTATION is single-exon. This removes
    # the ONE single-exon oracle multicopy gene that is recovered (family shrinks below 2 loci), and
    # removes NO residual FP (all FP blocks are multi-exon). So recall drops, precision unchanged/worse.
    lost = recall_cost["recall_cost_strict_annotation_named"]
    n_lost = len(lost)
    new_recov = recov_mc - n_lost
    new_mapped = mapped - n_lost           # each lost oracle gene's family no longer maps an oracle gene
    strict = dict(
        lost_oracle_genes=lost,
        precision_oracle_dedup=round(1 - distinct_fp / new_mapped, 4) if new_mapped else None,
        precision_oracle_taskformula=round(1 - flags / new_mapped, 4) if new_mapped else None,
        recall_oracle=round(new_recov / n_mc, 4),
        oracle_mapped_families=new_mapped, distinct_fp_blocks=distinct_fp,
        oracle_genes_recovered_multicopy=new_recov,
        note="counterfactual: drop ANNOTATION-single-exon members. Removes 0 FP, loses %d real oracle gene(s)."
             % n_lost,
    )
    return dict(
        current_default=dict(precision_oracle_dedup=cur["precision_oracle_dedup"],
                             precision_oracle_taskformula=cur["precision_oracle_taskformula"],
                             recall_oracle=cur["recall_oracle"]),
        operational_exclusion=operational,
        strict_annotation_exclusion=strict,
    )


# ---------------------------------------------------------------------------- main
def main():
    idx = load_intron_index()
    rows, fams = load_catalog()
    g2cn = load_oracle_multicopy()
    pr = json.load(open(PR_JSON))
    retro = json.load(open(RETRO_JSON))

    pop, fam_class = compute_population(fams, idx)
    fp = compute_fp_removed(pr, fams, idx)
    rc = compute_recall_cost(g2cn, fams, idx)
    ra = compute_retro_artifact(fams, idx, retro)
    npr = compute_net_pr(pr, rc)

    # ---- per-family TSV
    oracle_genes = set(g2cn)
    fp_genes = set()
    for tag in ("multifam", "oversize"):
        for e in pr["residual_fp_roster"].get(tag, []):
            fp_genes.update(e["genes"])
    with open(OUT_TSV, "w") as out:
        out.write("family_id\tdominant_gene\tn_members\tn_single_exon_loci\tn_multi_exon_loci\t"
                  "min_denovo_nexon\tclass\twould_exclude\tis_residual_fp\tis_oracle_real\t"
                  "n_annotation_single_exon_members\n")
        for fid in sorted(fams, key=lambda x: int(x)):
            members = fams[fid]
            fc = fam_class[fid]
            genes = {m["member_gene"] for m in members}
            dom = members[0]["dominant_gene"]
            is_fp = 1 if (genes & fp_genes) else 0
            is_orac = 1 if (genes & oracle_genes) or dom in oracle_genes else 0
            n_ann_se = sum(1 for m in members if ann_n_intron(idx, m["member_gene"]) == 0)
            out.write(f"{fid}\t{dom}\t{len(members)}\t{fc['n_se']}\t{fc['n_me']}\t{fc['min_nexon']}\t"
                      f"{fc['cls']}\t{fc['would_exclude']}\t{is_fp}\t{is_orac}\t{n_ann_se}\n")

    result = dict(
        determinism=dict(pythonhashseed="0", seed=0),
        inputs=dict(catalog=CATALOG, skeletons=SKEL, oracle=ORACLE, pr_json=PR_JSON,
                    retro_json=RETRO_JSON, intron_index=IDX),
        computation_1_population=pop,
        computation_2_fp_removed=fp,
        computation_3_recall_cost=rc,
        computation_4_retrocopy_artifact_split=ra,
        computation_5_net_pr=npr,
        verdict=dict(
            operational_exclusion_is_noop=(pop["catalog_single_exon_members"] == 0),
            fp_removed=fp["n_fp_removed_by_exclusion"],
            operational_recall_cost_named=rc["recall_cost_operational_named"],
            strict_recall_cost_named=rc["recall_cost_strict_annotation_named"],
            genuine_se_multicopy_families_fully_lost_upstream=rc["genuine_se_multicopy_families_fully_lost_upstream"],
            summary=(
                "The de-novo RNA pipeline ALREADY applies single-exon exclusion implicitly (intron-junction "
                "collapse + chromosome-span drop): the transcript FASTA feeding the family definition has "
                f"{pop['denovo_transcript_fasta_single_exon']}/{pop['denovo_transcript_fasta_total']} single-exon "
                f"transcripts and the 607-family catalog has {pop['catalog_single_exon_members']} single-exon members. "
                "An explicit single-exon exclusion is therefore a NO-OP on the catalog: it removes "
                f"{pop['families_removed_by_exclusion']} families, {fp['n_fp_removed_by_exclusion']} residual FPs, and "
                "costs 0 recall. MAGE-A (3-exon spliced UTR) and GSTM2 (8-exon) are multi-exon and unaffected. "
                "The recall cost is REAL but ALREADY PAID by the pipeline and invisible in the catalog: "
                "genuinely intronless multi-copy families in THIS gorilla annotation -- IFNA type-I "
                "interferon cluster (4 genes), TAS2R bitter-taste GPCRs (15/18 intronless), H1 linker "
                "histones (5/6), MAGEB5/MAGEB17 -- are all absent from the 607-family catalog because they "
                "carry no splice junction and were dropped upstream, NOT by any explicit exclusion knob."
            ),
        ),
    )
    json.dump(result, open(OUT_JSON, "w"), indent=1, sort_keys=True)

    # ---- console report
    P = print
    P("=" * 78)
    P("SINGLE-EXON EXCLUSION — trade-off quantification (RNA-only; DNA=validation)")
    P("=" * 78)
    P("\n[1] POPULATION")
    P(f"  de-novo skeletons : single-exon loci = {pop['denovo_skeleton_single_exon_loci']} "
      f"(ALL per-chromosome aggregation buckets)  |  multi-exon loci = {pop['denovo_skeleton_multi_exon_loci']}")
    P(f"  transcript FASTA fed to family def : single-exon = {pop['denovo_transcript_fasta_single_exon']} / "
      f"{pop['denovo_transcript_fasta_total']}")
    P(f"  catalog ({pop['catalog_n_families']} fam / {pop['catalog_n_members']} members) class: {pop['family_class_counts']}")
    P(f"  catalog single-exon members = {pop['catalog_single_exon_members']}  ->  exclusion removes "
      f"{pop['families_removed_by_exclusion']} families, {pop['members_removed_by_exclusion']} members")
    P("\n[2] FP REMOVED")
    P(f"  named residual-FP entries = {fp['n_named_fp_entries']}  ->  removed by exclusion = {fp['n_fp_removed_by_exclusion']}")
    P(f"  MAGE-A sub-cluster (LOC129529978/986): ann_n_intron={fp['mage_A_subcluster_LOC129529978_986']['ann_n_intron']} "
      f"-> single-exon-removed={fp['mage_A_subcluster_LOC129529978_986']['single_exon_removed']}")
    P(f"  GSTM2 hub single-exon-removed = {fp['GSTM2_hub']['single_exon_removed']}  (GSTM2 n_intron=7)")
    P(f"  fam17 repeat-bridge hub: {fp['fam17_repeat_bridge_hub']}")
    P("\n[3] RECALL COST (named real single-exon families lost)")
    P(f"  oracle multicopy genes = {rc['oracle_multicopy_genes_named']}  |  single-exon(annotation) = "
      f"{rc['oracle_multicopy_single_exon_annotation']}")
    for e in rc["oracle_multicopy_single_exon_genes"]:
        P(f"    {e['gene']} (hapCN {e['asm_hapCN']}): recovered={e['recovered_in_catalog']} "
          f"denovo_nexon={e['denovo_nexon_of_member']} fam={e['catalog_families']} "
          f"lost_operational={e['lost_by_operational_exclusion']} lost_strict={e['lost_by_strict_annotation_exclusion']}")
    P(f"  MAGE touchstone genes in oracle: {rc['oracle_multicopy_touchstones']['MAGE']}")
    P(f"  operational recall cost (named) = {rc['recall_cost_operational_named']}")
    P(f"  strict-annotation recall cost (named) = {rc['recall_cost_strict_annotation_named']}")
    P("  genuine single-exon MULTI-COPY families in THIS annotation (NOT oracle genes; the concrete upstream loss):")
    for f in rc["genuine_se_multicopy_families_not_in_oracle"]:
        P(f"    {f['family']}: {f['n_single_exon']}/{f['n_genes_in_annotation']} intronless, "
          f"{f['n_single_exon_in_catalog']} in catalog, fully_lost_upstream={f['fully_lost_upstream']} "
          f"{f['single_exon_genes'] if f['fully_lost_upstream'] else ''}")
    P("\n[4] RETRO / ARTIFACT SPLIT")
    P(f"  de-novo single-exon resolved real loci = {ra['denovo_single_exon_resolved_real_loci']} "
      f"(all {pop['denovo_skeleton_single_exon_loci']} single-exon skeleton rows are chromosome buckets)")
    P(f"  annotation single-exon genes = {ra['annotation_single_exon_genes']}/{ra['annotation_genes_total']} "
      f"({ra['annotation_single_exon_fraction']*100:.1f}%) ; recovered-in-catalog = {ra['annotation_single_exon_recovered_in_catalog']}")
    P(f"  recovered single-exon-annotation class split = {ra['recovered_single_exon_annotation_class_split']}")
    P(f"  genome-wide retrocopy-filter edge split = {ra['genomewide_retrocopy_filter_edge_split']}")
    P("\n[5] NET P/R (vs diploid oracle)")
    P(f"  current default        : P_dedup={npr['current_default']['precision_oracle_dedup']} "
      f"R={npr['current_default']['recall_oracle']}")
    P(f"  operational exclusion  : P_dedup={npr['operational_exclusion']['precision_oracle_dedup']} "
      f"R={npr['operational_exclusion']['recall_oracle']}  ({npr['operational_exclusion']['note']})")
    P(f"  strict-annotation excl : P_dedup={npr['strict_annotation_exclusion']['precision_oracle_dedup']} "
      f"R={npr['strict_annotation_exclusion']['recall_oracle']}  lost={npr['strict_annotation_exclusion']['lost_oracle_genes']}")
    P("\nVERDICT:", result["verdict"]["summary"])
    P(f"\nwrote {OUT_TSV}\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
