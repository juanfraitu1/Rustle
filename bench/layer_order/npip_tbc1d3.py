#!/usr/bin/env python3
"""NPIP/TBC1D3 on human CHM13 (2026-09-16): the layer-order integration and the nested edge-test lattice. There is one
subcommand per stage. The library is lattice_common.py in this directory. Reports: bench/LAYER_ORDER_NPIP_TBC1D3.md
and bench/NESTED_LATTICE_NPIP_TBC1D3.md.

Wave 7 (2026-09-24) folded the 12 files of bench/layer_order/ into this CLI and lattice_common.py. The old files are
at git tag notebook-2026-09-24 (`git show notebook-2026-09-24:bench/layer_order/<old>.py`). Each stage's body is the
old script's code, moved verbatim; only paths, imports and the helpers shared through lattice_common changed. Every
output is byte-identical to the old scripts' output, with the repo-root path fixed and PYTHONHASHSEED=0 (see
Determinism below).

Old command -> new command. Run from any directory, in the foreground, in this order. IS = ROOT/integrate_slim,
LAT = ROOT/lattice.
  lo_expr_recount.py OUT [--ignore IDS] EXTRA...  -> npip_tbc1d3.py expr-recount OUT [--ignore IDS] EXTRA...    (26 s)
  lo_corrected_tables.py                          -> npip_tbc1d3.py corrected-tables                               (3 s)
  (IS/P_N2_clusters.tsv is a FROZEN input of layer-order. The archived lo_p_variants.py wrote it; that script needs the
   off-repo light/scripts/layer_protein_bounded.py. See tag notebook-2026-09-20:archive/bench/layer_order/.)
  lo_analysis.py                                  -> npip_tbc1d3.py layer-order                                    (4 s)
  lattice_edges.py                                -> npip_tbc1d3.py lattice-edges                       (64 s, 2.5 GB)
  lattice_expr.py                                 -> npip_tbc1d3.py lattice-expr                                  (43 s)
  lattice_levels.py                               -> npip_tbc1d3.py lattice-levels                                (44 s)
  lattice_truth.py                                -> npip_tbc1d3.py lattice-truth                                  (8 s)
  lattice_filtration.py [--l1 c2_loose]           -> npip_tbc1d3.py lattice-filtration [--l1 c2_loose]        (5 s each)
  lattice_check_c2.py                             -> npip_tbc1d3.py lattice-check-c2                              (30 s)
  lattice_report_tables.py                        -> npip_tbc1d3.py lattice-report                                (26 s)
  (both reproduce blocks, in dependency order)    -> npip_tbc1d3.py all [--with-check-c2]                      (~5 min)
Old library names map to lattice_common.* (the table is in its docstring). Only layer-order uses the layer-order
operators below (containment, verdict, tournament, join_P_over_D*, refine, expr_layer), so they stay in this file.

--root DIR (or $LO_ROOT) is the directory that holds light/, heavy/, integrate_slim/ and lattice/. It may be given
before or after the subcommand. The default is /mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3, the FROZEN
2026-09-16 results. Every stage OVERWRITES its outputs under ROOT, so to rerun, copy the four directories and point
--root at the copy. Read-only inputs outside ROOT: o1_falsemerge/ (E1 catalog PAFs, regions, nodes, E0/E1),
winloci_data/{hgnc, soto_replication, Reference/chm13v2.0_RefSeq_full.gff.gz} and human_testis.t2t.bam. Tools:
samtools (PATH), numpy and scipy (the scorers), and the mcl_port Rust bin through bench/lib.py (layer-order REFINE;
$RUSTLE_MCL_PORT_BIN). docs/DATA.md lists the substrate.

Determinism (found in wave 7): in the old scripts, ties in 7 outputs were broken by set iteration order, which follows
Python's per-process string-hash seed. The 7 outputs are IS/expr_groups.tsv, IS/expr_sweep.tsv, IS/analysis.out,
LAT/member_groups.tsv, LAT/groups.tsv, LAT/triangle_drops.tsv and LAT/report_tables.md; one stdout line of
corrected-tables is affected too. The frozen 2026-09-16 files therefore hold one random seed's tie order: the same rows
and the same numbers, some in another order. This CLI re-executes itself with PYTHONHASHSEED=0 (HASH_SEED below), so
two runs give identical bytes. Those bytes equal the old scripts' output under PYTHONHASHSEED=0.
"""
import argparse
import collections
import csv
import itertools
import os
import random
import re
import sys
import time
from collections import defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import lattice_common as LC  # noqa: E402  (also puts bench/ on sys.path; paths are read as LC.X at call time)
import lib  # noqa: E402  (bench/lib.py: mcl -> the mcl_port Rust bin)
from lattice_common import (  # noqa: E402,F401
    CAT_CHROMS, CIG, ID_COL, ID_L3, LEVELS, PAF, SEF_MIN, SOTO_FIELDS, UF, IntervalIndex, Log, bip_jaccard, c2,
    c_tree, catalog_context, catalog_keys, chain_groups, chain_stats, components, count_reads, exonic_bases_in, f1,
    fmt, fmtv, fnum, gff_exon_index, groups, hgnc_all, hgnc_lookup, hgnc_tables, join_labels, key_of, load_exons,
    load_genes, membership, merge, refines, score_counts, score_lo, score_member_anchored, shortest_path, soto_flag,
    soto_load, soto_map_gene, soto_ok, split_counts, tests, truss3, tsv, uf_components, write, write_str,
    write_windows_bed)

HASH_SEED = "0"
FAMS = ("NPIP", "TBC1D3", "pooled")
MODES = ("any", "unique", "unique_mr")
# the RefSeq readthrough records that overlap a member on the same strand (expr-recount --ignore in the reproduce block)
READTHROUGH_OVER_MEMBERS = ("PKD1P3-NPIPA1", "LOC131696449", "PKD1P4-NPIPA8", "PKD1P5-LOC105376752", "PDXDC2P-NPIPB14P",
                            "TBC1D3P1-DHX40P1")
# universe genes heavy/EXPR.counts.tsv lacks (expr-recount EXTRA ids in the reproduce block)
EXPR_RECOUNT_EXTRA = ("USP6NL", "LOC100420408", "LOC100420311", "TBC1D29P", "LOC124905656", "LOC100420289")
# corrected member rule (was lo_corrected_tables.py)
DESC_NPIP = re.compile(r"nuclear pore complex[- ]interacting protein", re.I)
DESC_TBC = re.compile(r"^TBC1 domain family member 3( |[A-Z]|$)")
FAMNAME = re.compile(r"NPIP|TBC1D3")
SYM = re.compile(r"^(NPIP.*|TBC1D3(?:$|[A-Z]|P\d|-).*)$")  # the original symbol rule (members.py), for the audit line

# layer-order tables (lo_analysis.py module level); set by load_layer_order(), used only by the layer-order stage
U = NAME = MEMB = SIDE = LAYERS = C_MID_AS_BUILT = P_EDGES = P_W = D_EDGE_ROWS = D_EDGES = READS = None
CT_LAB = CT_ROWS = CT_LIT = CT_CLUSTERS = None


# ==================================================================================================================== stages
def cmd_expr_recount(args):
    """was bench/layer_order/lo_expr_recount.py OUT_TSV [--ignore ID,ID,...] EXTRA_ID ...

    NPIP/TBC1D3 layer order (integration) — testis read counts for EXPR on the corrected universe.

    Rule = heavy/scripts/expr_counts.py as documented in heavy/README.md (EXPR) and re-derived independently by
    verify_slim_expr/recount_all.py (0 of 321 genes differ): human_testis.t2t.bam, primary reads only (samtools -F 2308),
    read blocks split at N and D (M/=/X consume both), a read counts for gene G if >= 1 block overlaps >= 1 bp of any
    exon of G (strand ignored); unique = the read's blocks hit exons of exactly ONE RefSeq gene/pseudogene genome-wide
    and that gene is G; exon-less gene/pseudogene records count on a gene-body exon. (One function with lattice-expr:
    lattice_common.count_reads.)

    This stage is the recount_all.py logic with the gene set = heavy/EXPR.counts.tsv (321) + extra ids given on the
    command line (universe genes EXPR lacks). It prints any difference against EXPR.counts.tsv (expected: none).
    Writes OUT_TSV and OUT_TSV.windows.bed.

    --ignore (audit 2026-09-16): also count n_reads_unique_mr = reads whose exon hits, after removing the listed records,
    are exactly {G} (for a listed record G itself, the other listed records are removed). Used with the RefSeq
    readthrough records that overlap a member on the same strand, so 'unique' agrees with the member rule's
    one-record-per-copy intent.
    """
    out_tsv = args.out
    ignore = set(args.ignore.split(",")) if args.ignore is not None else set()
    extra = set(args.extra)
    expr = {r["gene_id"]: r for r in csv.DictReader(open(f"{LC.HEAVY}/EXPR.counts.tsv"), delimiter="\t")}
    want = set(expr) | extra
    genes, exons, index = gff_exon_index(strip_prefix=True)
    print(f"[gff] gene/pseudogene ids {len(genes)}; distinct exons {sum(len(v) for v in exons.values())}")
    spans = defaultdict(list)
    for g in want:
        c, s, e = genes["gene-" + g]
        spans[c].append((s, e))
    bed = out_tsv + ".windows.bed"
    nwin, _ = write_windows_bed(bed, spans)
    n_any, n_uni, n_uni_mr, n_rec, n_hit = count_reads(bed, index, want, ignore)
    print(f"[bam] windows {nwin}; primary records {n_rec}; on >= 1 exon {n_hit}")
    diff = 0
    with open(out_tsv, "w") as out:
        out.write("gene_id\tn_reads_any\tn_reads_unique\tin_EXPR_table" + ("\tn_reads_unique_mr" if ignore else "") + "\n")
        for g in sorted(want):
            a, u = n_any[g], n_uni[g]
            if g in expr:
                if (a, u) != (int(expr[g]["n_reads_any"]), int(expr[g]["n_reads_unique"])):
                    diff += 1
                    print("DIFF", g, (a, u), (expr[g]["n_reads_any"], expr[g]["n_reads_unique"]))
            else:
                print(f"EXTRA {g} any {a} unique {u}")
            out.write(f"{g}\t{a}\t{u}\t{'yes' if g in expr else 'no'}" + (f"\t{n_uni_mr[g]}" if ignore else "") + "\n")
    print(f"EXPR.counts.tsv genes recounted {len(expr)}; differing {diff}; extra genes {len(extra - set(expr))}")
    if ignore:
        print(f"unique_mr: ignored records {sorted(ignore)}; genes whose unique count rises: "
              + ", ".join(f"{g} {n_uni[g]}->{n_uni_mr[g]}" for g in sorted(want) if n_uni_mr[g] != n_uni[g]))


def cmd_corrected_tables(args):
    """was bench/layer_order/lo_corrected_tables.py

    NPIP/TBC1D3 layer order (integration) — apply the verification fixes to the light/heavy builds.

    Writes corrected tables NEXT TO the originals (suffix .corrected.tsv; originals untouched):
      light/members.corrected.tsv            member rule applied the same way to both families (fix: members asymmetry)
      light/P.groups.corrected.tsv           in_P genes + every other U gene with a §6ko protein as 'P|other:<gene>' (audit
                                             2026-09-16: PKD1, DHX40); non-coding members are NOT P singletons (fix: P encoding)
      light/P.members_status.corrected.tsv   every corrected member, group_id 'P|NA' when not in P
      light/D.groups.corrected.tsv           layer_dna.py rerun on the corrected members, duplicate-coordinate node fix,
                                             catalog-absent members -> 'D|NA' (fix: per-chromosome catalogs, cheapest route);
                                             + every other U gene that is a catalog node, with its catalog cluster (audit
                                             2026-09-16: TBC1D26 -> MCL24, status clustered_nonmember_group)
      light/D.edges.corrected.tsv            same edge rule as layer_dna.py (pre-MCL --dump-graph edges touching member groups)
      light/C.groups.corrected.tsv           TBC1D3 positional clusters moved out of clade_L1 (fix: positional != clade);
                                             members absent from the reference trees -> 'NA' (not assessed)
      light/universe.corrected.tsv           U = members + P + D + C (S1 dropped: deferred layer), per-layer universe flags
      light/truth_hgnc.corrected.tsv, truth_soto.corrected.tsv, truth_guided.corrected.tsv (E0/E1 filled for every row),
      light/truth_literature_subfamilies.corrected.tsv
      heavy/EXPR.counts.corrected.tsv        + rows for universe genes EXPR lacked (integrate_slim/expr_recount.tsv)

    Corrected member rule (both families):
      description  RefSeq description matches /nuclear pore complex[- ]interacting protein/i  (NPIP)
                   or /^TBC1 domain family member 3( |[A-Z]|$)/                                (TBC1D3; excludes TBC1D30-32)
      readthrough  a record whose description says 'readthrough' and whose name/description names NPIP or TBC1D3 is a
                   member only if its span overlaps NO same-strand description-rule record (one record per copy: the
                   readthrough is kept only when its family part has no gene record of its own).
    """
    from truth import excluded  # §6ko r2 filter (was protein_families.excluded)
    LIGHT, HEAVY, INT = LC.LIGHT, LC.HEAVY, LC.INT
    genes, by_coord = load_genes()

    # ------------------------------------------------------------------ A. members
    old_members = {r["gene_id"]: r for r in tsv(f"{LIGHT}/members.tsv")}
    desc = {}
    for gid, g in genes.items():
        if DESC_NPIP.search(g["description"]):
            desc[gid] = "NPIP"
        elif DESC_TBC.search(g["description"]):
            desc[gid] = "TBC1D3"
    rt_cand = {gid: ("NPIP" if "NPIP" in (g["name"] + " " + g["description"]) else "TBC1D3")
               for gid, g in genes.items()
               if "readthrough" in g["description"] and FAMNAME.search(g["name"] + " " + g["description"])}
    members, rt_log = {}, []
    for gid, fam in desc.items():
        members[gid] = (fam, "description" + ("; symbol" if SYM.match(genes[gid]["name"]) else ""))
    for gid, fam in sorted(rt_cand.items()):
        g = genes[gid]
        s, e = int(g["start0"]), int(g["end"])
        ov = sorted(genes[d]["name"] for d in desc if genes[d]["chrom"] == g["chrom"] and genes[d]["strand"] == g["strand"]
                    and int(genes[d]["start0"]) < e and s < int(genes[d]["end"]))
        rt_log.append(f"{g['name']} ({g['description']}): overlaps same-strand description-rule records {ov or '-'} -> "
                      f"{'NOT member' if ov else 'member'}")
        if not ov:
            members[gid] = (fam, "readthrough (family part has no own gene record)")
    sym_only = [genes[g]["name"] for g in genes if SYM.match(genes[g]["name"]) and g not in members and g not in rt_cand]
    lit = {r["name"]: r for r in tsv(LC.LIT)}
    db = soto_load()
    exons = load_exons()
    mrows = []
    for gid, (fam, basis) in members.items():
        g = genes[gid]
        m = soto_map_gene(db, g["name"], g["chrom"], g["strand"], exons.get(gid) or [(int(g["start0"]), int(g["end"]))])
        mrows.append({"gene_id": gid, "name": g["name"], "biotype": g["biotype"], "chrom": g["chrom"], "start": g["start0"],
                      "end": g["end"], "strand": g["strand"], "family": fam, "member_basis": basis,
                      "in_lit_truth_31": "yes" if g["name"] in lit else "no", "description": g["description"],
                      "in_original_members": "yes" if gid in old_members else "no", **m})
    mrows.sort(key=lambda x: (x["family"], x["chrom"], int(x["start"])))
    # audit (2026-09-16, minor): the readthrough restriction is not a general 'one record per copy' rule. Flag every member
    # whose span lies inside another same-strand member record (non-readthrough containment), and report the member count
    # under the alternative 'fold the enclosed record into its encloser' option. No layer is changed by this column.
    for x in mrows:
        s, e = int(x["start"]), int(x["end"])
        enc = sorted(y["name"] for y in mrows if y is not x and y["chrom"] == x["chrom"] and y["strand"] == x["strand"]
                     and int(y["start"]) <= s and e <= int(y["end"]) and (int(y["end"]) - int(y["start"])) > (e - s))
        x["span_inside_member_record"] = ",".join(enc)
    n_enc = collections.Counter(x["family"] for x in mrows if x["span_inside_member_record"])
    print(f"[members] records whose span lies inside another same-strand member record: "
          f"{[(x['name'], x['span_inside_member_record']) for x in mrows if x['span_inside_member_record']]}; member count if "
          f"they are folded into the encloser: NPIP {sum(x['family'] == 'NPIP' for x in mrows) - n_enc['NPIP']}, TBC1D3 "
          f"{sum(x['family'] == 'TBC1D3' for x in mrows) - n_enc['TBC1D3']}")
    mcols = ["gene_id", "name", "biotype", "chrom", "start", "end", "strand", "family", "member_basis", "in_lit_truth_31",
             "description", "in_original_members", "span_inside_member_record"] + SOTO_FIELDS
    write_str(f"{LIGHT}/members.corrected.tsv", mrows, mcols)
    M = {r["gene_id"]: r for r in mrows}
    fam_of = {g: r["family"] for g, r in M.items()}
    print(f"[members] corrected {len(M)}: NPIP {sum(v == 'NPIP' for v in fam_of.values())}, TBC1D3 "
          f"{sum(v == 'TBC1D3' for v in fam_of.values())}; added {sorted(genes[g]['name'] for g in set(M) - set(old_members))}; "
          f"removed {sorted(genes[g]['name'] for g in set(old_members) - set(M))}; symbol-rule names outside the rule {sym_only or '-'}")
    for x in rt_log:
        print("   readthrough:", x)

    # ------------------------------------------------------------------ B. P
    pidx = tsv(f"{LIGHT}/work/P/proteins.index.tsv")
    P_universe = {r["gene_id"] for r in pidx}
    Prow = [r for r in tsv(f"{LIGHT}/P.groups.tsv") if r["p_status"] == "in_P"]
    for r in Prow:
        r["is_member"] = "yes" if r["gene_id"] in M else "no"
        r["member_family"] = fam_of.get(r["gene_id"], "")
    pcols = list(tsv(f"{LIGHT}/P.groups.tsv")[0].keys())
    write_str(f"{LIGHT}/P.groups.corrected.tsv", Prow, pcols)
    Pg = {r["gene_id"]: r["group_id"] for r in Prow}
    old_status = {r["gene_id"]: r for r in tsv(f"{LIGHT}/P.members_status.tsv")}
    cds = {r["gene_id"] for r in tsv(f"{LIGHT}/work/refseq/cds.tsv")}
    srows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["name"])):
        if gid in Pg:
            st, grp = "in_P", Pg[gid]
        elif gid in old_status:
            st, grp = old_status[gid]["p_status"], "P|NA"
        else:
            st = ";".join(x for x in ("excluded_r2_biotype" if excluded(r["biotype"], 2) else "",
                                      "no_CDS_in_RefSeq" if gid not in cds else "") if x) or "CHECK"
            grp = "P|NA"
        srows.append({"gene_id": gid, "name": r["name"], "biotype": r["biotype"], "family": r["family"], "p_status": st,
                      "group_id": grp, "in_P_universe": "yes" if gid in P_universe else "no"})
    write_str(f"{LIGHT}/P.members_status.corrected.tsv",
          srows, ["gene_id", "name", "biotype", "family", "p_status", "group_id", "in_P_universe"])
    assert all((r["group_id"] != "P|NA") == (r["in_P_universe"] == "yes") for r in srows), "member in P universe but not in_P"
    print(f"[P] in_P rows {len(Prow)} (members {sum(r['is_member'] == 'yes' for r in Prow)}); members not in P "
          f"{sum(r['group_id'] == 'P|NA' for r in srows)}: {dict(collections.Counter(r['p_status'] for r in srows if r['group_id'] == 'P|NA'))}")

    # ------------------------------------------------------------------ C. D (layer_dna.py logic, corrected members)
    D_all, E0_all, E1_all, D_catalog_of = {}, {}, {}, {}
    drows, erows = [], []
    for tag, c in LC.CATALOGS.items():
        k2g = catalog_keys(by_coord, c["kind"], c["path"])
        mem = membership(k2g, c["D"])
        e0 = membership(k2g, c["E0"])
        e1 = membership(k2g, c["E1"])
        assert all(mem[g][0] == e1[g][0] for g in mem), f"D != E1 on {tag}"
        for g in mem:
            D_catalog_of[g] = tag
            D_all[g] = mem[g]
            E0_all[g] = e0[g]
            E1_all[g] = e1[g]
        mclus = {mem[m][0] for m in M if m in mem and mem[m][0]}
        in_groups = sorted((g for g in mem if mem[g][0] in mclus),
                           key=lambda x: (mem[x][0], genes[x]["chrom"], int(genes[x]["start0"]), x))
        for gid in in_groups:
            g = genes[gid]
            cl, k, f = mem[gid]
            drows.append({"gene_id": gid, "group_id": f"D|{tag}|{cl}", "name": g["name"], "biotype": g["biotype"],
                          "chrom": g["chrom"], "start": g["start0"], "end": g["end"], "strand": g["strand"],
                          "is_member": "yes" if gid in M else "no", "member_family": fam_of.get(gid, ""),
                          "catalog": tag, "cluster_id": cl, "record_key": k, "folded_into": f, "status": "clustered"})
        for m in M:
            if m in mem and not mem[m][0]:
                g = genes[m]
                drows.append({"gene_id": m, "group_id": f"D|{tag}|singleton:{m}", "name": g["name"], "biotype": g["biotype"],
                              "chrom": g["chrom"], "start": g["start0"], "end": g["end"], "strand": g["strand"],
                              "is_member": "yes", "member_family": fam_of[m], "catalog": tag, "cluster_id": "",
                              "record_key": mem[m][1], "folded_into": mem[m][2], "status": "in_catalog_unclustered"})
        keyset = {key_of(mem[g][1]) for g in in_groups}
        grp = {g: f"D|{tag}|{mem[g][0]}" for g in mem if mem[g][0]}
        for line in open(c["D"] + ".graph.tsv"):
            f = line.rstrip("\n").split("\t")
            if len(f) != 3:
                continue
            a, b = key_of(f[0]), key_of(f[1])
            if a not in keyset and b not in keyset:
                continue
            for ga in k2g.get(a, []) or [{"gene_id": "?", "name": "?"}]:
                for gb in k2g.get(b, []) or [{"gene_id": "?", "name": "?"}]:
                    erows.append({"catalog": tag, "u_gene_id": ga["gene_id"], "u_name": ga["name"],
                                  "v_gene_id": gb["gene_id"], "v_name": gb["name"], "weight": f[2],
                                  "u_group": grp.get(ga["gene_id"], ""), "v_group": grp.get(gb["gene_id"], ""),
                                  "u_key": f[0], "v_key": f[1]})
        print(f"[D] {tag}: keys {len(k2g)} (genes mapped {len(mem)}); member clusters {sorted(mclus)}; genes in them "
              f"{len(in_groups)} (members {sum(g in M for g in in_groups)})")
    for m in M:
        if m not in D_all:
            g = genes[m]
            drows.append({"gene_id": m, "group_id": "D|NA", "name": g["name"], "biotype": g["biotype"], "chrom": g["chrom"],
                          "start": g["start0"], "end": g["end"], "strand": g["strand"], "is_member": "yes",
                          "member_family": fam_of[m], "catalog": "", "cluster_id": "", "record_key": "", "folded_into": "",
                          "status": "not_in_any_catalog (outside D's universe)"})
    dcols = ["gene_id", "group_id", "name", "biotype", "chrom", "start", "end", "strand", "is_member", "member_family",
             "catalog", "cluster_id", "record_key", "folded_into", "status"]
    write_str(f"{LIGHT}/D.groups.corrected.tsv", drows, dcols)
    write_str(f"{LIGHT}/D.edges.corrected.tsv",
          erows, ["catalog", "u_gene_id", "u_name", "v_gene_id", "v_name", "weight", "u_group", "v_group", "u_key", "v_key"])
    old_d = {(r["gene_id"], r["group_id"]) for r in tsv(f"{LIGHT}/D.groups.tsv") if r["status"] == "clustered"}
    new_d = {(r["gene_id"], r["group_id"]) for r in drows if r["status"] == "clustered"}
    old_e = {tuple(r[c] for c in ("u_gene_id", "v_gene_id", "weight")) for r in tsv(f"{LIGHT}/D.edges.tsv")}
    new_e = {tuple(r[c] for c in ("u_gene_id", "v_gene_id", "weight")) for r in erows}
    print(f"[D] clustered (gene, group) rows original {len(old_d)} corrected {len(new_d)} identical {old_d == new_d}; "
          f"edge rows original {len(old_e)} corrected {len(new_e)} identical {old_e == new_e}; not in any catalog: "
          f"{sorted(genes[r['gene_id']]['name'] for r in drows if r['group_id'] == 'D|NA')}")
    Dg = {r["gene_id"]: r["group_id"] for r in drows if r["group_id"] != "D|NA"}

    # ------------------------------------------------------------------ D. C
    oc = {r["gene_id"]: r for r in tsv(f"{LIGHT}/C.groups.tsv")}
    crows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["chrom"], int(kv[1]["start"]))):
        o = oc.get(gid)
        if o is None or o["in_reference_trees"] != "yes":
            crows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_reference_trees": "no",
                          "C_fine": "C|NA", "C_mid": "Cmid|NA", "C_L1": "CL1|NA", "positional_cluster": "",
                          "evidence_fine": "", "evidence_mid": "", "evidence_L1": "",
                          "note": "not in the 31-record reference trees: not assessed (outside C's universe)"})
            continue
        l1, ev1, pos, note = o["clade_L1"], o["evidence_L1"], "", ""
        if r["family"] == "TBC1D3" and not l1.startswith("CL1|singleton"):
            pos = f"{l1[4:]} (intron tree {ev1.split('SH-aLRT/UFBoot ')[-1]}; POSITIONAL, not a clade: o1_ledger §6jp/§6jr/§6js, L14917-14923)"
            l1, ev1, note = f"CL1|singleton:{gid}", "", "clade_L1 cleared: TBC1D3 cluster1/2 are positional"
        crows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_reference_trees": "yes",
                      "C_fine": o["clade_id"], "C_mid": o["clade_mid"], "C_L1": l1, "positional_cluster": pos,
                      "evidence_fine": o["evidence_fine"], "evidence_mid": o["evidence_mid"], "evidence_L1": ev1,
                      "note": note})
    ccols = ["gene_id", "name", "family", "in_reference_trees", "C_fine", "C_mid", "C_L1", "positional_cluster",
             "evidence_fine", "evidence_mid", "evidence_L1", "note"]
    write_str(f"{LIGHT}/C.groups.corrected.tsv", crows, ccols)
    Cr = {r["gene_id"]: r for r in crows if r["in_reference_trees"] == "yes"}
    print(f"[C] rows {len(crows)}; in trees {len(Cr)}; TBC1D3 clade_L1 cleared "
          f"{sum(1 for r in crows if r['positional_cluster'])}")

    # ------------------------------------------------------------------ E. universe
    memP = {Pg[g] for g in Pg if g in M}
    memD = {Dg[g] for g in Dg if g in M and "singleton" not in Dg[g]}
    U = set(M) | {g for g in Pg if Pg[g] in memP} | {g for g in Dg if Dg[g] in memD} | set(Cr)
    side = {}
    for g in U:
        fams = set()
        if g in M:
            fams.add(fam_of[g])
        for Lg in (Pg, Dg):
            if g in Lg:
                fams |= {fam_of[m] for m in M if Lg.get(m) == Lg[g]}
        assert len(fams) == 1, (g, fams)
        side[g] = fams.pop()
    expr = {r["gene_id"]: r for r in tsv(f"{HEAVY}/EXPR.counts.tsv")}
    rec = {r["gene_id"]: r for r in tsv(f"{INT}/expr_recount.tsv")}
    hgt = hgnc_tables()
    old_soto = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_soto.tsv")}
    urows, hrows, srows2, grows = [], [], [], []
    order = sorted(U, key=lambda g: (genes[g]["chrom"], int(genes[g]["start0"]), g))
    for g in order:
        r = genes[g]
        name = r["name"]
        ex = rec.get(name)
        assert ex is not None, name
        h, how = hgnc_lookup(hgt, g, name)
        hrows.append({"gene_id": g, "name": name, "is_member": "yes" if g in M else "no", "match": how,
                      "hgnc_id": h["hgnc_id"] if h else "", "hgnc_symbol": h["symbol"] if h else "",
                      "locus_group": h["locus_group"] if h else "", "gene_group": h["gene_group"] if h else "",
                      "gene_group_id": h["gene_group_id"] if h else ""})
        if g in old_soto:
            s = dict(old_soto[g])
            s["is_member"] = "yes" if g in M else "no"
        else:
            m = soto_map_gene(db, name, r["chrom"], r["strand"], exons.get(g) or [(int(r["start0"]), int(r["end"]))])
            flag = soto_flag(m)
            s = {"gene_id": g, "name": name, "is_member": "yes" if g in M else "no", "flag": flag, **m}
        srows2.append(s)
        c = Cr.get(g)
        layers = [L for L, ok in (("P", Pg.get(g) in memP), ("D", Dg.get(g) in memD), ("C", g in Cr)) if ok]
        urows.append({"gene_id": g, "name": name, "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                      "end": r["end"], "strand": r["strand"], "is_member": "yes" if g in M else "no",
                      "member_family": fam_of.get(g, ""), "family_side": side[g],
                      "in_P_universe": "yes" if g in P_universe else "no",
                      "P_group": Pg.get(g, f"P|other:{g}" if g in P_universe else "P|NA"),
                      "in_D_universe": "yes" if g in D_all else "no",
                      "D_group": (f"D|{D_catalog_of[g]}|{D_all[g][0]}" if D_all[g][0] else f"D|{D_catalog_of[g]}|singleton:{g}")
                      if g in D_all else "D|NA",
                      "D_folded_into": D_all[g][2] if g in D_all else "",
                      "in_C_universe": "yes" if c else "no",
                      "C_L1": c["C_L1"] if c else "CL1|NA", "C_mid": c["C_mid"] if c else "Cmid|NA",
                      "C_fine": c["C_fine"] if c else "C|NA", "positional_cluster": c["positional_cluster"] if c else "",
                      "layers_placing_gene_with_member": ",".join(layers),
                      "E0_group": (f"E0|{D_catalog_of[g]}|{E0_all[g][0]}" if E0_all[g][0] else f"E0|{D_catalog_of[g]}|singleton:{g}")
                      if g in D_all else "",
                      "hgnc_gene_group_id": h["gene_group_id"] if h else "", "hgnc_gene_group": h["gene_group"] if h else "",
                      "soto_families": s["soto_families"], "soto_flag": s["flag"],
                      "n_reads_any": ex["n_reads_any"], "n_reads_unique": ex["n_reads_unique"],
                      "description": r["description"]})
    ucols = list(urows[0].keys())
    write_str(f"{LIGHT}/universe.corrected.tsv", urows, ucols)
    # audit (2026-09-16, important): the group tables must carry every U gene of the layer's universe, not only the genes
    # of member groups, so that containment computed from P.groups / D.groups equals containment computed from the
    # universe table. P universe rule (spec, 'P defined on coding genes'): every U gene with a §6ko protein
    # (work/P/proteins.index.tsv), including coding genes P does not place with a member ('P|other:<gene>').
    # D universe: every U gene that is a node (record) of one of the two catalogs, with its catalog cluster.
    plen = {r["gene_id"]: r["length_aa"] for r in pidx}
    have_p = {r["gene_id"] for r in Prow}
    extra_p = []
    for g in order:
        if g in P_universe and g not in have_p:
            r = genes[g]
            extra_p.append({"gene_id": g, "group_id": f"P|other:{g}", "name": r["name"], "biotype": r["biotype"],
                            "chrom": r["chrom"], "start": r["start0"], "end": r["end"], "strand": r["strand"],
                            "protein_len_aa": plen[g], "is_member": "yes" if g in M else "no",
                            "member_family": fam_of.get(g, ""), "p_status": "in_P_universe_not_with_member"})
    write_str(f"{LIGHT}/P.groups.corrected.tsv", Prow + extra_p, pcols)
    have_d = {r["gene_id"] for r in drows}
    extra_d = []
    for g in order:
        if g in D_all and g not in have_d:
            r = genes[g]
            cl, k, f = D_all[g]
            extra_d.append({"gene_id": g, "group_id": f"D|{D_catalog_of[g]}|{cl}" if cl else f"D|{D_catalog_of[g]}|singleton:{g}",
                            "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                            "end": r["end"], "strand": r["strand"], "is_member": "yes" if g in M else "no",
                            "member_family": fam_of.get(g, ""), "catalog": D_catalog_of[g], "cluster_id": cl,
                            "record_key": k, "folded_into": f,
                            "status": "clustered_nonmember_group" if cl else "in_catalog_unclustered"})
    write_str(f"{LIGHT}/D.groups.corrected.tsv", drows + extra_d, dcols)
    print(f"[P/D group tables] rows added for U genes of the layer universe outside member groups: P "
          f"{[(genes[x['gene_id']]['name'], x['group_id']) for x in extra_p]}; D "
          f"{[(genes[x['gene_id']]['name'], x['group_id'], x['folded_into']) for x in extra_d]}")
    pu_tab = {r["gene_id"]: r["group_id"] for r in Prow + extra_p}
    du_tab = {r["gene_id"]: r["group_id"] for r in drows + extra_d if r["group_id"] != "D|NA"}
    assert pu_tab == {u["gene_id"]: u["P_group"] for u in urows if u["in_P_universe"] == "yes"}, "P table != universe"
    assert du_tab == {u["gene_id"]: u["D_group"] for u in urows if u["in_D_universe"] == "yes"}, "D table != universe"
    write_str(f"{LIGHT}/truth_hgnc.corrected.tsv", hrows, list(hrows[0].keys()))
    write_str(f"{LIGHT}/truth_soto.corrected.tsv", srows2, ["gene_id", "name", "is_member", "flag"] + SOTO_FIELDS)
    # guided truth: every U gene in a catalog + every gene sharing an E0/E1 cluster with a member, cluster ALWAYS filled
    mcl0 = {(D_catalog_of[m], E0_all[m][0]) for m in M if m in E0_all and E0_all[m][0]}
    mcl1 = {(D_catalog_of[m], E1_all[m][0]) for m in M if m in E1_all and E1_all[m][0]}
    gl = [g for g in D_all if g in U or (D_catalog_of[g], E0_all[g][0]) in mcl0 or (D_catalog_of[g], E1_all[g][0]) in mcl1]
    for g in sorted(gl, key=lambda x: (genes[x]["chrom"], int(genes[x]["start0"]), x)):
        r = genes[g]
        grows.append({"gene_id": g, "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"], "start": r["start0"],
                      "end": r["end"], "is_member": "yes" if g in M else "no", "in_universe": "yes" if g in U else "no",
                      "catalog": D_catalog_of[g],
                      "E0_cluster": f"{D_catalog_of[g]}|{E0_all[g][0]}" if E0_all[g][0] else "unclustered",
                      "E0_folded_into": E0_all[g][2],
                      "E1_cluster": f"{D_catalog_of[g]}|{E1_all[g][0]}" if E1_all[g][0] else "unclustered",
                      "E1_folded_into": E1_all[g][2], "record_key": E0_all[g][1]})
    write_str(f"{LIGHT}/truth_guided.corrected.tsv", grows, list(grows[0].keys()))
    # literature
    lrows = []
    for gid, r in sorted(M.items(), key=lambda kv: (kv[1]["family"], kv[1]["chrom"], int(kv[1]["start"]))):
        t = lit.get(r["name"])
        lrows.append({"gene_id": gid, "name": r["name"], "family": r["family"], "in_literature_truth": "yes" if t else "no",
                      "level1": t["level1"] if t else "", "level2": t["level2"] if t else "",
                      "level1_kind": ("" if not t else "phylogenetic subfamily (Dishuck 2025)" if r["family"] == "NPIP"
                                      else "POSITIONAL genomic cluster (Guitart 2024; o1_ledger L14917-14923)"),
                      "npipb_named_subfamily": ("yes" if r["name"] in {"NPIPB3", "NPIPB4", "NPIPB5", "NPIPB11", "NPIPB12",
                                                                        "NPIPB13"} else "no") if r["family"] == "NPIP" and t else ""})
    write_str(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv", lrows, list(lrows[0].keys()))
    # EXPR corrected
    xrows = []
    names_U = {genes[g]["name"]: g for g in U}
    for n, r in sorted(rec.items()):
        o = expr.get(n)
        gid = "gene-" + n
        row = dict(o) if o else {"gene_id": n, "name": n, "chrom": genes[gid]["chrom"], "n_reads_any": r["n_reads_any"],
                                 "n_reads_unique": r["n_reads_unique"],
                                 "expressed_ge3_any": "yes" if int(r["n_reads_any"]) >= 3 else "no",
                                 "expressed_ge3_unique": "yes" if int(r["n_reads_unique"]) >= 3 else "no",
                                 "seed_family": "", "biotype": genes[gid]["biotype"], "in_S2": "no", "S2_group_id": ""}
        row["member_family_corrected"] = fam_of.get(gid, "")
        row["in_universe_corrected"] = "yes" if n in names_U else "no"
        row["row_source"] = "EXPR.counts.tsv" if o else "integrate_slim/expr_recount.tsv (same rule)"
        xrows.append(row)
    write_str(f"{HEAVY}/EXPR.counts.corrected.tsv", xrows, list(xrows[0].keys()))
    print(f"[U] {len(U)} genes (members {len(M)}); per layer universe: P {sum(u['in_P_universe'] == 'yes' for u in urows)}, "
          f"D {sum(u['in_D_universe'] == 'yes' for u in urows)}, C {sum(u['in_C_universe'] == 'yes' for u in urows)}; "
          f"placed with a member by P {sum('P' in u['layers_placing_gene_with_member'].split(',') for u in urows)}, "
          f"D {sum('D' in u['layers_placing_gene_with_member'].split(',') for u in urows)}; family side "
          f"{dict(collections.Counter(side.values()))}; chromosomes {dict(collections.Counter(genes[g]['chrom'] for g in U))}")
    print(f"[truths] HGNC match {dict(collections.Counter(h['match'] for h in hrows))}; Soto flags "
          f"{dict(collections.Counter(s['flag'] for s in srows2))}; guided rows {len(grows)}; EXPR corrected rows {len(xrows)}")


# ============================================================================================================== layer order (was lo_analysis.py)
def load_layer_order():
    """The module-level tables of lo_analysis.py (U, layers, P/D edges, READS, C_tree), which only the
    layer-order stage needs (lattice-edges computes c_tree itself)."""
    global U, NAME, MEMB, SIDE, LAYERS, C_MID_AS_BUILT, P_EDGES, P_W, D_EDGE_ROWS, D_EDGES, READS
    global CT_LAB, CT_ROWS, CT_LIT, CT_CLUSTERS
    U = {r["gene_id"]: r for r in tsv(f"{LC.LIGHT}/universe.corrected.tsv")}
    NAME = {g: r["name"] for g, r in U.items()}
    MEMB = {g for g, r in U.items() if r["is_member"] == "yes"}
    SIDE = {g: r["family_side"] for g, r in U.items()}

    LAYERS = {"P": layer_from("P_group", "in_P_universe"), "D": layer_from("D_group", "in_D_universe"),
              "C_L1": layer_from("C_L1", "in_C_universe"), "C_mid": layer_from("C_mid", "in_C_universe"),
              "C_fine": layer_from("C_fine", "in_C_universe")}
    C_MID_AS_BUILT = dict(LAYERS["C_mid"])
    # C_mid as built holds ONE recovered group (the named NPIPB subfamily) and singletons elsewhere, so it is not a level
    # between C_L1 and C_fine (A6-9 and B6-9 are singletons in it). Hierarchical mid level = C_mid ∨ C_fine (by construction
    # c(C_mid ⊇ C_fine) = 1). The as-built layer is kept as C_mid_ab in every table.
    LAYERS["C_mid"] = join_labels(C_MID_AS_BUILT, LAYERS["C_fine"])
    LAYERS["C_mid_ab"] = C_MID_AS_BUILT

    P_EDGES = {}
    for r in tsv(f"{LC.LIGHT}/P.edges.tsv"):
        P_EDGES[frozenset((r["u_gene_id"], r["v_gene_id"]))] = (float(r["weight"]), float(r["identity"]),
                                                                 float(r["coverage_longer"]))
    P_W = {k: v[0] for k, v in P_EDGES.items()}
    D_EDGE_ROWS = tsv(f"{LC.LIGHT}/D.edges.corrected.tsv")
    D_EDGES = {}
    for r in D_EDGE_ROWS:
        if "?" in (r["u_gene_id"], r["v_gene_id"]) or r["u_gene_id"] == r["v_gene_id"]:
            continue
        k = frozenset((r["u_gene_id"], r["v_gene_id"]))
        D_EDGES[k] = max(D_EDGES.get(k, 0.0), float(r["weight"]))
    _rec = {r["gene_id"]: r for r in tsv(f"{LC.INT}/expr_recount.tsv")}
    READS = {g: {"unique": int(r["n_reads_unique"]), "any": int(r["n_reads_any"]),
                 "unique_mr": int(_rec[NAME[g]]["n_reads_unique_mr"])} for g, r in U.items()}
    CT_LAB, CT_ROWS, CT_LIT, CT_CLUSTERS = c_tree(U, NAME)
    LAYERS.update({"Ctree_top": CT_LAB["Ctree_top"], "Ctree_min": CT_LAB["Ctree_min"]})


def fam_genes(fam, side=None):
    side = side or SIDE
    return set(side) if fam == "pooled" else {g for g in side if side[g] == fam}


def layer_from(col, uni_col):
    return {g: r[col] for g, r in U.items() if r[uni_col] == "yes"}


def groups_of(lab, genes):
    out = collections.defaultdict(set)
    for g in genes:
        out[lab[g]].add(g)
    return out


def pairs_of(lab, genes):
    s = set()
    for G in groups_of(lab, genes).values():
        for a, b in itertools.combinations(sorted(G), 2):
            s.add((a, b))
    return s


# ---------------------------------------------------------------------------------------------------------- catalogs
def catalog_all():
    """gene_id -> (catalog, cluster_id or '', record key, folded representative or '') for every gene of both E1 catalogs
    (lattice_common.catalog_keys/membership), plus gene rows."""
    genes, by_coord = load_genes()
    out = {}
    for tag, c in LC.CATALOGS.items():
        k2g = catalog_keys(by_coord, c["kind"], c["path"])
        for g, (cl, k, f) in membership(k2g, c["D"]).items():
            out[g] = (tag, cl, k, f)
    return out, genes


# ---------------------------------------------------------------------------------------------------------- containment
def containment(X, Y, fam, layers, side=None):
    LX, LY = layers[X], layers[Y]
    S = set(LX) & set(LY) & fam_genes(fam, side)
    pY, pX = pairs_of(LY, S), pairs_of(LX, S)
    c = Fraction(len(pY & pX), len(pY)) if pY else None
    gY = [G for G in groups_of(LY, S).values() if len(G) >= 2]
    inside = sum(1 for G in gY if len({LX[g] for g in G}) == 1)
    return {"family": fam, "X": X, "Y": Y, "n_genes_in_both": len(S), "pairs_Y": len(pY), "pairs_Y_in_X": len(pY & pX),
            "c_X_contains_Y": "vacuous (0 pairs)" if c is None else fmt(c), "_c": c,
            "groups_Y_ge2": len(gY), "groups_Y_inside_one_X": inside,
            "nesting": "NA" if not gY else fmt(Fraction(inside, len(gY)))}


def verdict(a, b):
    """X above Y iff c(X ⊇ Y) > c(Y ⊇ X); NA when either containment is over 0 pairs."""
    if a["_c"] is None or b["_c"] is None:
        return "NA (0 pairs on one side)"
    if a["_c"] > b["_c"]:
        return f"{a['X']} above {a['Y']}"
    if b["_c"] > a["_c"]:
        return f"{a['Y']} above {a['X']}"
    return "tie"


def tournament(names, fam, layers, side=None):
    rows = []
    for X, Y in itertools.combinations(names, 2):
        a, b = containment(X, Y, fam, layers, side), containment(Y, X, fam, layers, side)
        rows.append({"family": fam, "X": X, "Y": Y, "c(X>=Y)": a["c_X_contains_Y"], "n_pairs_Y": a["pairs_Y"],
                     "groups_Y_inside_X": f"{a['groups_Y_inside_one_X']}/{a['groups_Y_ge2']}",
                     "c(Y>=X)": b["c_X_contains_Y"], "n_pairs_X": b["pairs_Y"],
                     "groups_X_inside_Y": f"{b['groups_Y_inside_one_X']}/{b['groups_Y_ge2']}",
                     "n_genes": a["n_genes_in_both"], "verdict": verdict(a, b)})
    return rows


def join_P_over_D(LP, LD):
    genes = set(LP) | set(LD)
    vert = {g: ("D", LD[g]) if g in LD else ("Ponly", g) for g in genes}
    links = []
    for G in groups_of(LP, set(LP)).values():
        G = sorted(G)
        links += [(vert[G[0]], vert[h]) for h in G[1:]]
    comps = uf_components(set(vert.values()), links)
    lab = {}
    for comp in comps:
        members = sorted(g for g in genes if vert[g] in comp)
        name = "Pjoin|" + ("+".join(sorted({v[1] for v in comp if v[0] == "D"})) or members[0])
        for g in members:
            lab[g] = name
    return lab


def join_P_over_D_whole(cat, n2, pidx_genes):
    """JOIN on whole catalog groups. cat: gene -> (catalog, cluster, key, folded) for all catalog genes; n2: gene -> k=2 MCL
    cluster for N_2 genes (member clusters = PC1/PC2); pidx_genes: genes with a protein. Returns labels for every gene of
    a member-containing component, and the P label used per gene."""
    Pl = {}
    for g in pidx_genes:
        Pl[g] = n2.get(g, f"P|own:{g}")
    vert = {}
    for g in set(cat) | set(Pl):
        if g in cat and cat[g][1]:
            vert[g] = ("D", f"D|{cat[g][0]}|{cat[g][1]}")
        elif g in cat:
            vert[g] = ("Dsingle", g)
        else:
            vert[g] = ("Ponly", g)
    links = []
    for G in groups_of(Pl, set(Pl)).values():
        if len(G) < 2:
            continue
        G = sorted(G)
        links += [(vert[G[0]], vert[h]) for h in G[1:]]
    comps = uf_components(set(vert.values()), links)
    by_vert = collections.defaultdict(list)
    for g, v in vert.items():
        by_vert[v].append(g)
    lab = {}
    for comp in comps:
        genes_c = sorted(g for v in comp for g in by_vert[v])
        if not set(genes_c) & MEMB:
            continue
        name = "PjoinW|" + "+".join(sorted(v[1] for v in comp if v[0] == "D"))
        for g in genes_c:
            lab[g] = name
    return lab, Pl


def refine(LX, LY, edge_w, use_mcl, outside="attach"):
    """X's grouping recomputed inside each Y group (see module docstring). outside: 'attach' | 'single' | 'together'
    (the X group's genes outside Y's universe form one part of their own)."""
    lab = {}
    for xg, G in groups_of(LX, set(LX)).items():
        inside = collections.defaultdict(set)
        out_genes = []
        for g in G:
            if g in LY:
                inside[LY[g]].add(g)
            else:
                out_genes.append(g)
        parts = []
        for yg, part in inside.items():
            if len(part) == 1 or not use_mcl:
                parts.append(set(part))
                continue
            E = {tuple(sorted(k)): w for k, w in edge_w.items() if k <= part}
            cl = lib.mcl(E) if E else []
            got = set().union(*map(set, cl)) if cl else set()
            parts += [set(c) for c in cl] + [{g} for g in part - got]
        if outside == "together" and out_genes:
            parts.append(set(out_genes))
            out_genes = []
        for g in sorted(out_genes):
            if outside == "single":
                parts.append({g})
            elif use_mcl:
                best = max(parts, key=lambda p: (sum(edge_w.get(frozenset((g, h)), 0.0) for h in p), len(p)), default=None)
                if best is not None and sum(edge_w.get(frozenset((g, h)), 0.0) for h in best) > 0:
                    best.add(g)
                else:
                    parts.append({g})
            else:
                best = max(parts, key=len, default=None)
                if best is not None:
                    best.add(g)
                else:
                    parts.append({g})
        for i, p in enumerate(sorted(parts, key=lambda p: sorted(p))):
            for g in p:
                lab[g] = f"{xg}#r{i}" if len(parts) > 1 else xg
    return lab


def expr_layer(lab, expressed, edges):
    """EXPR(L): {L group: [components >= 2]}, plus dropped expressed singletons per group."""
    out, dropped = {}, {}
    for grp, G in groups_of(lab, set(lab)).items():
        Ex = G & expressed
        links = [tuple(k) for k in edges(Ex)]
        comps = [c for c in uf_components(Ex, links)]
        out[grp] = [c for c in comps if len(c) >= 2]
        dropped[grp] = sorted(g for c in comps if len(c) == 1 for g in c)
    return out, dropped


def edge_set(kind):
    return {"P": [P_EDGES], "P_ref": [P_EDGES], "D": [D_EDGES], "P_join": [P_EDGES, D_EDGES]}.get(kind)


def edge_fn(kind):
    pools = edge_set(kind)
    if kind == "comembership" or pools is None:  # C layers: clade co-membership
        return lambda S: [frozenset(p) for p in itertools.combinations(sorted(S), 2)]
    return lambda S: [k for pool in pools for k in pool if k <= S]


def expr_labels(comps_by_group):
    lab = {}
    for grp, comps in comps_by_group.items():
        for i, c in enumerate(sorted(comps, key=lambda c: sorted(c))):
            for g in c:
                lab[g] = f"{grp}#e{i}"
    return lab


# ---------------------------------------------------------------------------------------------------------- DNA pair stats
def dna_pair_stats(pairs):
    """identity / coverage for D edges (catalog, key_u, key_v) from the catalog PAF, mcl_families graph rule
    (records >= 300 bp, identity >= 0.70, pooled identity, union of aligned intervals on the longer gene / its exon-union
    length). Exon-union lengths: c16_19_20 nodes.tsv; c15_17_22 light/work/refseq/exons.tsv (checked via weight)."""
    exlen = {}
    for r in tsv(LC.NODES_C16):
        exlen[("c16_19_20", f"{r['chrom']}:{int(r['start']) + 1}-{r['end']}")] = sum(
            int(b) - int(a) for a, b in (x.split("-") for x in r["exons"].split(",")))
    for r in tsv(f"{LC.LIGHT}/work/refseq/exons.tsv"):
        exlen.setdefault(("c15_17_22", f"{r['chrom']}:{int(r['start0']) + 1}-{r['end']}"),
                         sum(int(b) - int(a) for a, b in (x.split("-") for x in r["exons"].split(","))))
    want = collections.defaultdict(set)
    for cat, a, b in pairs:
        want[cat].add(frozenset((a, b)))
    acc = {}
    for cat, S in want.items():
        keys = set().union(*S)
        with open(LC.PAF[cat]) as fh:
            for line in fh:
                f = line.split("\t", 12)
                q, t = f[0], f[5]
                if q == t or q not in keys or t not in keys or frozenset((q, t)) not in S:
                    continue
                nm, bl = int(f[9]), int(f[10])
                if bl < 300 or nm / max(bl, 1) < 0.70:
                    continue
                k = (cat, frozenset((q, t)))
                e = acc.setdefault(k, {"nm": 0, "bl": 0, "iv": collections.defaultdict(list)})
                e["nm"] += nm
                e["bl"] += bl
                e["iv"][q].append((int(f[2]), int(f[3])))
                e["iv"][t].append((int(f[7]), int(f[8])))
    out = {}
    for (cat, pr), e in acc.items():
        a, b = sorted(pr)
        la, lb = exlen.get((cat, a)), exlen.get((cat, b))
        if la is None or lb is None:
            continue
        longer = a if la >= lb else b
        merged = merge(sorted(e["iv"][longer]))
        cov = min(1.0, sum(y - x for x, y in merged) / max(la if longer == a else lb, 1))
        out[(cat, pr)] = (e["nm"] / e["bl"], cov)
    return out


def cmd_layer_order(args):
    """was bench/layer_order/lo_analysis.py

    NPIP/TBC1D3 layer order (integration, revised after the 2026-09-16 audit): containment, tournament, enforcement cost,
    truth agreement, EXPR operator, disagreement lists — layers P (protein), D (= the RefSeq E1 guided catalog, §6kl; NOT
    §0★★ clause 2/4), C (subfamily clades) on the CORRECTED tables (bench/layer_order/lo_corrected_tables.py).
    Spec: docs/superpowers/specs/2026-09-16-family-layer-order-design.md (SCOPE AMENDMENT + Definitions, incl. the
    2026-09-16 audit amendments). S1/S2/S3 are deferred (user scope 2026-09-16 15:03) and not read.

    Conventions (all stated in bench/LAYER_ORDER_NPIP_TBC1D3.md):
      layer universe  P: U gene with a §6ko protein (work/P/proteins.index.tsv), including coding genes P does not place with a
                      member ('P|other:<gene>' = its own group). D: U gene that is a node of one of the two E1 catalogs.
                      C layers: U gene that is a leaf of the §6js reference trees (literature C: 31 leaves; C_tree: the 30 leaves
                      common to the exon and intron trees). Outside its universe a gene is 'not in layer', never a singleton.
      C layers        literature-anchored (CIRCULAR reference): C_L1, C_mid (:= as-built C_mid ∨ C_fine), C_mid_as_built,
                      C_fine. Clause-5 (§0★★) split-system layers: all SH-aLRT > 75 splits of either tree restricted to the
                      common leaves, smaller side = cluster, kept iff compatible with every other kept split (Buneman).
                      Ctree_top = maximal clusters (the partition whose pairs are 'share >= 1 subfamily'); Ctree_min = minimal.
      c(X ⊇ Y)        |pairs(Y) ∩ pairs(X)| / |pairs(Y)| over genes in U_X ∩ U_Y (∩ family); 0 pairs -> vacuous (verdict NA).
      nesting         fraction of Y groups (>= 2 genes after restriction to U_X ∩ U_Y) inside one X group.
      JOIN(P over D)  components of the graph whose vertices are D groups (+ P-universe genes outside D as own vertices), two
                      vertices joined when two of their genes share a P group. Two variants: 'U' (D groups cut to U, as first
                      reported) and 'whole' (whole catalog groups; P groups for genes outside U from the k = 2 MCL on N_2,
                      integrate_slim/P_N2_clusters.tsv; genes outside N_2 are their own P group), U re-closed afterwards.
      REFINE(X in Y)  X's grouping recomputed inside each Y group: P -> mcl_port.mcl (inflation 2.8) on P.edges induced on
                      X group ∩ Y group; C -> clade co-membership ∩ Y group. A gene outside Y's universe is a FREE CHOICE,
                      variants reported: 'attach' (P: to the refined part with the largest summed X-edge weight; C: kept
                      with the largest part of its clade), 'single' (P: its own group; C: isolated), and for C 'together'
                      (the clade's outside genes form one part of their own).
      EXPR(L)         within each L group, connected components (>= 2 genes) of the subgraph induced on expressed genes.
                      Expressed = reads >= t, t = 3 in the main tables (sweep t = 1..5). Read modes: any (PRIMARY, spec rule),
                      unique (the read hits exons of one RefSeq record genome-wide), unique_mr (same, ignoring the 6 readthrough
                      records that overlap a member on the same strand). Edges: P -> P.edges; D -> D.edges.corrected; C ->
                      clade co-membership; P_join -> P.edges ∪ D.edges; P_ref -> P.edges. EXPR(L) ⊆ L holds by construction;
                      with co-membership edges no L group can split, so T2 holds trivially there (reported as guaranteed).
      truth scores    (A) on U (recall CONDITIONED ON THE PREDICTION: U is the closure of the scored layers);
                      (B) member-anchored, layer-independent: pairs with >= 1 member endpoint over every gene of a truth group
                      that holds a member (HGNC gene group; Soto family), restricted to the layer's universe genome-wide
                      (P: genes with a protein; D: catalog nodes). Bipartite F only when both sides have >= 1 pair.
    """
    LIGHT, INT = LC.LIGHT, LC.INT
    load_layer_order()
    say = Log()

    LP, LD = LAYERS["P"], LAYERS["D"]
    say(f"[universe] U {len(U)}; members {len(MEMB)}; universes P {len(LP)} D {len(LD)} C_lit {len(LAYERS['C_L1'])} "
        f"C_tree {len(LAYERS['Ctree_top'])}")
    tab_p = {r["gene_id"]: r["group_id"] for r in tsv(f"{LIGHT}/P.groups.corrected.tsv")}
    tab_d = {r["gene_id"]: r["group_id"] for r in tsv(f"{LIGHT}/D.groups.corrected.tsv") if r["group_id"] != "D|NA"}
    say(f"[check] P.groups.corrected == universe P labels: {tab_p == LP}; D.groups.corrected == universe D labels: "
        f"{tab_d == LD}")

    # ------------------------------------------------ containment + tournament
    main_layers = ["P", "D", "Ctree_top", "Ctree_min"]
    ref_layers = ["C_L1", "C_mid", "C_mid_ab", "C_fine"]
    allL = main_layers + ref_layers
    crow, trow = [], []
    for fam in FAMS:
        for X in allL:
            for Y in allL:
                if X != Y:
                    r = containment(X, Y, fam, LAYERS)
                    r.pop("_c")
                    crow.append(r)
        trow += tournament(allL, fam, LAYERS)
    write_str(f"{INT}/containment.tsv", crow)
    write_str(f"{INT}/tournament.tsv", trow)
    for r in trow:
        say(f"[tournament] {r['family']:6s} {r['X']:>9s} vs {r['Y']:<9s}: c(X⊇Y) {r['c(X>=Y)']} [{r['n_pairs_Y']}] "
            f"{r['groups_Y_inside_X']} | c(Y⊇X) {r['c(Y>=X)']} [{r['n_pairs_X']}] {r['groups_X_inside_Y']} | genes "
            f"{r['n_genes']} -> {r['verdict']}")
    # one common gene set for every main layer
    common = set.intersection(*(set(LAYERS[x]) for x in main_layers + ["C_L1"]))
    lay_c = {x: {g: LAYERS[x][g] for g in common} for x in allL}
    crow2 = []
    for fam in FAMS:
        crow2 += tournament(allL, fam, lay_c)
    write_str(f"{INT}/tournament_common_genes.tsv", crow2)
    say(f"[common] genes in every layer universe (P ∩ D ∩ C_tree ∩ C_lit): {len(common)} "
        f"(NPIP {len(common & fam_genes('NPIP'))}, TBC1D3 {len(common & fam_genes('TBC1D3'))})")
    for r in crow2:
        if r["family"] == "pooled":
            say(f"[common] pooled {r['X']:>9s} vs {r['Y']:<9s}: {r['c(X>=Y)']} [{r['n_pairs_Y']}] {r['groups_Y_inside_X']} | "
                f"{r['c(Y>=X)']} [{r['n_pairs_X']}] {r['groups_X_inside_Y']} -> {r['verdict']}")
    # all compatible clusters (not only the partition levels) inside one P / D group
    for fam, (L, K) in CT_CLUSTERS.items():
        for X in ("P", "D"):
            gid = {NAME[g]: g for g in LAYERS[X]}
            ok = [s for s in K if all(n in gid for n in s)]
            inside = sum(1 for s in ok if len({LAYERS[X][gid[n]] for n in s}) == 1)
            say(f"[C_tree hierarchy] {fam}: compatible supported clusters {len(K)}; inside one {X} group {inside}/{len(ok)} "
                f"(clusters with every leaf in {X}'s universe)")
    write_str(f"{INT}/ctree_clusters.tsv", CT_ROWS)
    write_str(f"{INT}/ctree_literature_groups.tsv", CT_LIT)

    # ------------------------------------------------ enforcement
    cat, genes = catalog_all()
    pidx_genes = {r["gene_id"] for r in tsv(f"{LIGHT}/work/P/proteins.index.tsv")}
    n2 = {}
    for r in tsv(f"{INT}/P_N2_clusters.tsv"):
        n2[r["gene_id"]] = r["k2_cluster"]
    # the k = 2 MCL member clusters must equal P's groups
    for grp, G in groups_of(LP, {g for g in LP if not LP[g].startswith("P|other")}).items():
        assert len({n2[g] for g in G}) == 1 and sum(1 for x in n2.values() if x == n2[next(iter(G))]) == len(G), grp
    Pj = join_P_over_D(LP, LD)
    Pjw, Pl_w = join_P_over_D_whole(cat, n2, pidx_genes)
    Pr = refine(LP, LD, P_W, use_mcl=True, outside="attach")
    Prs = refine(LP, LD, P_W, use_mcl=True, outside="single")
    enf = dict(LAYERS)
    enf.update({"P_join": Pj, "P_join_whole": Pjw, "P_ref": Pr, "P_ref_single": Prs})
    cl_names = ("C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min")
    for cl in cl_names:
        for Yn, LY in (("D", LD), ("P", LP)):
            enf[f"{cl}_ref{Yn}"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="attach")
            enf[f"{cl}_ref{Yn}_single"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="single")
            enf[f"{cl}_ref{Yn}_together"] = refine(LAYERS[cl], LY, {}, use_mcl=False, outside="together")
    extra_w = sorted(genes[g]["name"] for g in set(Pjw) - set(U))
    say(f"[JOIN whole] member-containing components: " + " | ".join(
        f"{k}: {len(v)} genes" for k, v in sorted(groups_of(Pjw, set(Pjw)).items())) + f"; genes outside U {len(extra_w)}: "
        f"{extra_w}")
    for g in sorted(set(Pjw) - set(U), key=lambda x: genes[x]["name"]):
        say(f"   {genes[g]['name']} biotype {genes[g]['biotype']} catalog group {cat[g][1] if g in cat else '-'} P label "
            f"{Pl_w.get(g, 'no protein')}")
    changes = []
    plan = [("P_join", Pj, LP), ("P_join_whole", Pjw, LP), ("P_ref", Pr, LP), ("P_ref_single", Prs, LP)] + [
        (f"{c}_{s}", enf[f"{c}_{s}"], LAYERS[c]) for c in cl_names for s in ("refD", "refD_single", "refD_together", "refP", "refP_single",
                                                                "refP_together")]
    for name, lab, ref in plan:
        S = set(ref)
        before, after = pairs_of(ref, S), pairs_of(lab, S & set(lab))
        added = sorted((NAME[a], NAME[b]) for a, b in after - before)
        removed = sorted((NAME[a], NAME[b]) for a, b in before - after)
        extra = sorted(genes[g]["name"] for g in set(lab) - S)
        changes.append({"enforced": name, "pairs_on_original_universe_before": len(before), "after": len(after),
                        "pairs_added": len(added), "pairs_removed": len(removed),
                        "genes_added_to_universe": len(extra),
                        "added_examples": "; ".join(f"{a}-{b}" for a, b in added[:12]),
                        "removed_examples": "; ".join(f"{a}-{b}" for a, b in removed[:12])})
        say(f"[enforce] {name}: pairs on original universe {len(before)} -> {len(after)} (+{len(added)} / -{len(removed)}); "
            f"universe +{len(extra)} genes")
    write_str(f"{INT}/enforcement_changes.tsv", changes)
    for X, Y in (("P_join", "D"), ("P_join_whole", "D"), ("D", "P_ref"), ("D", "P_ref_single"), ("D", "C_L1_refD"),
                 ("P", "C_L1_refP"), ("P", "C_L1_refP_single"), ("D", "Ctree_top_refD"), ("P", "Ctree_top_refP")):
        r = containment(X, Y, "pooled", enf, side={**SIDE, **{g: "x" for g in set(enf[X]) | set(enf[Y]) if g not in SIDE}})
        say(f"[nesting] c({X} ⊇ {Y}) pooled = {r['c_X_contains_Y']} (pairs {r['pairs_Y']}); groups inside one "
            f"{r['groups_Y_inside_one_X']}/{r['groups_Y_ge2']}")
    erow = []
    for g in sorted(set(U) | set(Pjw), key=lambda g: (SIDE.get(g, "~"), genes[g]["name"])):
        erow.append({"gene_id": g, "name": genes[g]["name"], "family_side": SIDE.get(g, "outside U (whole-group JOIN)"),
                     "is_member": "yes" if g in MEMB else "no",
                     **{k: enf[k].get(g, "NA") for k in ("P", "D", "P_join", "P_join_whole", "P_ref", "P_ref_single",
                                                         "C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min")}})
    write_str(f"{INT}/enforced_partitions.tsv", erow)

    # ------------------------------------------------ truths (A): on U
    hg = {g: r["hgnc_gene_group_id"] for g, r in U.items() if r["hgnc_gene_group_id"]}
    so = {g: r["soto_families"] for g, r in U.items() if r["soto_flag"] == "ok" and r["soto_families"]}
    e0 = {g: r["E0_group"] for g, r in U.items() if r["E0_group"]}
    lit = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
    litL1 = {g: r["level1"] for g, r in lit.items() if r["level1"] and r["family"] == "NPIP"}
    litL2 = {g: r["level2"] for g, r in lit.items() if r["level2"]}
    named = {g: ("named" if r["npipb_named_subfamily"] == "yes" else f"other:{g}") for g, r in lit.items()
             if r["family"] == "NPIP" and r["in_literature_truth"] == "yes"}
    all_hg = hgnc_all(genes)
    db = soto_load()
    exons = load_exons()
    for g in set(Pjw) - set(U):
        if g in all_hg:
            hg[g] = all_hg[g]
        fs = soto_ok(db, exons, genes, g)
        if fs:
            so[g] = fs
    say(f"[truth A] whole-group JOIN genes outside U: HGNC {sorted((genes[g]['name'], hg[g]) for g in set(Pjw) - set(U) if g in hg)}"
        f"; Soto ok {sorted((genes[g]['name'], so[g]) for g in set(Pjw) - set(U) if g in so)}")
    PU = set(LP)
    PU_join_w = {g for g in Pjw if g in pidx_genes} | PU
    trows = []
    plan = [("P", "as built", LP, PU), ("P", "JOIN over D, U-restricted (P universe)", Pj, PU),
            ("P", "JOIN over D, U-restricted (own universe)", Pj, set(Pj)),
            ("P", "JOIN over D, whole groups (P universe, re-closed)", Pjw, PU_join_w),
            ("P", "JOIN over D, whole groups (own universe, re-closed)", Pjw, set(Pjw)),
            ("P", "REFINE in D, USP6NL attached", Pr, PU), ("P", "REFINE in D, USP6NL single", Prs, PU),
            ("D", "as built = after (D unchanged)", LD, set(LD))]
    for layer, variant, lab, uni in plan:
        for tname, truth in (("HGNC gene_group_id (superfamily-level for TBC1D3)", hg), ("Soto family (flag ok)", so)):
            for fam in FAMS:
                side = {g: SIDE.get(g, "TBC1D3") for g in uni}  # whole-group JOIN genes are all on the TBC1D3 side
                r = score_lo(lab, truth, fam_genes(fam, side) & uni)
                trows.append({"layer": layer, "variant": variant, "truth": tname, "family": fam, **r})
    for layer, variant, lab, tname, truth in (
            ("D", "as built", LD, "E0 guided catalog (construction sensitivity, not truth)", e0),
            ("C_L1", "as built", LAYERS["C_L1"], "literature L1 NPIPA|NPIPB (CIRCULAR)", litL1),
            ("C_mid_ab", "as built", C_MID_AS_BUILT, "literature named NPIPB subfamily (CIRCULAR)", named),
            ("C_fine", "as built", LAYERS["C_fine"], "literature L2 paralog groups (CIRCULAR)", litL2),
            ("Ctree_top (clause 5)", "as built", LAYERS["Ctree_top"], "literature L1 (NPIP only)", litL1),
            ("Ctree_top (clause 5)", "as built", LAYERS["Ctree_top"], "literature L2 paralog groups", litL2),
            ("Ctree_min (clause 5)", "as built", LAYERS["Ctree_min"], "literature L2 paralog groups", litL2),
            ("Ctree_root (clause 5, rooted variant)", "as built", CT_LAB["Ctree_root"], "literature L1 (NPIP only)", litL1)):
        for fam in FAMS:
            r = score_lo(lab, truth, fam_genes(fam) & set(lab))
            trows.append({"layer": layer, "variant": variant, "truth": tname, "family": fam, **r})
    trows.append({"layer": "C_mid (= C_mid_ab ∨ C_fine)", "variant": "not scored",
                  "truth": "named NPIPB ∨ L2: identical by construction (both sides are the same join)", "family": "-"})
    cols = ["layer", "variant", "truth", "family", "n_genes", "truth_pairs", "pred_pairs", "tp_pairs", "pair_precision",
            "pair_recall", "bip_F_count(§6ks)", "bip_R", "bip_P", "bip_F_jaccard"]
    write_str(f"{INT}/truth_agreement.tsv", trows, cols)
    for r in trows:
        if r["family"] in ("TBC1D3", "pooled", "NPIP") and r.get("n_genes"):
            say(f"[truth A] {r['layer']} | {r['variant']} | {r['truth'][:30]} | {r['family']} n={r['n_genes']} "
                f"P {r['pair_precision']} R {r['pair_recall']} bipF {r['bip_F_count(§6ks)']} (pairs pred {r['pred_pairs']} "
                f"truth {r['truth_pairs']})")

    # ------------------------------------------------ truths (B): member-anchored, layer-independent gene sets
    hg_parts = {g: frozenset(v.split("|")) for g, v in all_hg.items()}
    mem_hg = set().union(*(hg_parts.get(m, frozenset()) for m in MEMB))
    hg_group_genes = {g for g, p in hg_parts.items() if p & mem_hg}
    n_2227_sym = sum(1 for r in tsv(LC.HGNC) if "2227" in r["gene_group_id"].split("|"))
    say(f"[truth B] HGNC groups holding a member: {sorted(mem_hg)}; HGNC symbols in those groups {n_2227_sym}; RefSeq "
        f"genes in them {len(hg_group_genes)}, with a §6ko protein {len(hg_group_genes & pidx_genes)}, catalog nodes "
        f"{len(hg_group_genes & set(cat))}")
    so_parts = {g: frozenset(v.split(";")) for g, v in so.items() if g in U}
    mem_so = set().union(*(so_parts.get(m, frozenset()) for m in MEMB))
    cand = {r["best_refseq_gene_id"] for r in tsv(f"{LIGHT}/truth_soto_families.tsv")
            if r["family_id"] in mem_so and r["best_refseq_gene_id"]}
    n_soto_genes = sum(1 for r in tsv(f"{LIGHT}/truth_soto_families.tsv") if r["family_id"] in mem_so)
    for g in cand - set(U):
        fs = soto_ok(db, exons, genes, g)
        if fs and set(fs.split(";")) & mem_so:
            so_parts[g] = frozenset(fs.split(";"))
    so_group_genes = {g for g, p in so_parts.items() if p & mem_so}
    say(f"[truth B] Soto families holding a member (flag ok): {sorted(mem_so)}; Soto genes in them {n_soto_genes}; RefSeq "
        f"genes forward-mapped into them (flag ok) {len(so_group_genes)} (outside U "
        f"{sorted(genes[g]['name'] for g in so_group_genes - set(U))})")
    brows = []
    for tname, parts, grp_genes in (("HGNC gene group", hg_parts, hg_group_genes), ("Soto family", so_parts, so_group_genes)):
        for layer, variant, lab, uni in (
                ("P", "as built", LP, pidx_genes), ("P", "JOIN U-restricted", Pj, pidx_genes),
                ("P", "JOIN whole groups", Pjw, pidx_genes), ("P", "REFINE attach", Pr, pidx_genes),
                ("P", "REFINE single", Prs, pidx_genes), ("D", "as built", LD, set(cat))):
            for fam in ("NPIP", "TBC1D3", "pooled"):
                mem_f = {m for m in MEMB if fam == "pooled" or SIDE[m] == fam}
                gset = ((set(U) | grp_genes) & uni)
                # a gene outside U is outside every member group of every layer here, except whole-group JOIN genes
                labx = dict(lab)
                if layer == "D":
                    for g in gset - set(labx):
                        if g in cat:
                            labx[g] = f"D|{cat[g][0]}|{cat[g][1]}" if cat[g][1] else f"D|single:{g}"
                r = score_member_anchored(labx, parts, gset, mem_f & gset)
                brows.append({"truth": tname, "layer": layer, "variant": variant, "family": fam,
                              "gene_set": "U ∪ all genes of member-holding truth groups, ∩ layer universe genome-wide", **r})
                say(f"[truth B] {tname} | {layer} {variant} | {fam}: genes {r['n_genes']} (members {r['n_members']}) "
                    f"P {r['pair_precision']} R {r['pair_recall']} (pred {r['pred_pairs']} truth {r['truth_pairs']} tp "
                    f"{r['tp_pairs']})")
    write_str(f"{INT}/truth_member_anchored.tsv", brows)

    # ------------------------------------------------ EXPR
    xrows, t2rows, grows, sweep = [], [], [], []
    layer_set = {"P": LP, "D": LD, "P_join": Pj, "P_ref": Pr, "C_L1": LAYERS["C_L1"], "C_mid": LAYERS["C_mid"],
                 "C_fine": LAYERS["C_fine"], "Ctree_top": LAYERS["Ctree_top"], "Ctree_min": LAYERS["Ctree_min"]}
    # T2 precondition E_M ⊆ E_L inside each L group (layer edges)
    pre = {}
    for L, M in itertools.permutations(layer_set, 2):
        S = set(layer_set[L]) & set(layer_set[M])
        gM = [G for G in groups_of(layer_set[M], S).values() if len(G) >= 2]
        if not gM or any(len({layer_set[L][g] for g in G}) != 1 for G in gM):
            continue
        fM, fL = edge_fn(M), edge_fn(L)
        miss = 0
        tot = 0
        for G in groups_of(layer_set[M], S).values():
            if len(G) < 2:
                continue
            eM, eL = set(fM(G)), set(fL(G))
            tot += len(eM)
            miss += len(eM - eL)
        pre[(L, M)] = (tot, miss)
    for (L, M), (tot, miss) in sorted(pre.items()):
        say(f"[T2 precondition] {L} ⊇ {M} on data: M-edges inside M groups {tot}; not in E_L {miss} -> "
            f"{'holds' if miss == 0 else 'FAILS'}")
    for mode in MODES:
        for t in (1, 2, 3, 4, 5):
            expressed = {g for g in U if READS[g][mode] >= t}
            mem_ex = {f: sum(1 for m in MEMB if SIDE[m] == f and m in expressed) for f in ("NPIP", "TBC1D3")}
            for emode in ("layer_edges", "comembership"):
                EX, splits, distinct = {}, [], set()
                for L, lab in layer_set.items():
                    kind = "comembership" if emode == "comembership" else L
                    comps, dropped = expr_layer(lab, expressed, edge_fn(kind))
                    EX[L] = comps
                    for grp, cs in comps.items():
                        for c in cs:
                            distinct.add(frozenset(c))
                        if len(cs) > 1:
                            splits.append(f"{L}:{grp}: " + " | ".join(",".join(sorted(NAME[g] for g in c)) for c in cs))
                    if t == 3:
                        for grp, G in groups_of(lab, set(lab)).items():
                            ex = sorted(G & expressed)
                            if len(G) < 2 or not (G & MEMB):
                                continue
                            xrows.append({"expr_mode": mode, "edges": emode, "layer": L, "group": grp, "n_genes": len(G),
                                          "n_expressed": len(ex), "n_EXPR_groups": len(comps[grp]),
                                          "EXPR_groups": " | ".join(",".join(sorted(NAME[g] for g in c))
                                                                    for c in sorted(comps[grp], key=lambda c: -len(c))),
                                          "expressed_dropped_singletons": ",".join(NAME[g] for g in dropped[grp]),
                                          "split": "yes" if len(comps[grp]) > 1 else "no"})
                        bad = sum(1 for grp, cs in comps.items() for c in cs if len({lab[g] for g in c}) != 1)
                        grows.append({"expr_mode": mode, "edges": emode, "layer": L,
                                      "EXPR_groups": sum(len(cs) for cs in comps.values()),
                                      "not_inside_one_L_group (0 by construction)": bad})
                nchk = nviol = 0
                for L, M in itertools.permutations(layer_set, 2):
                    if (L, M) not in pre:
                        continue
                    S = set(layer_set[L]) & set(layer_set[M])
                    eL = expr_labels(EX[L])
                    viol, checked = [], 0
                    for c in (c for cs in EX[M].values() for c in cs):
                        cS = c & S
                        if len(cS) < 2:
                            continue
                        checked += 1
                        if len({eL.get(g, f"none:{g}") for g in cS}) != 1:
                            viol.append(",".join(sorted(NAME[g] for g in cS)))
                    nchk += checked
                    nviol += len(viol)
                    if t == 3:
                        t2rows.append({"expr_mode": mode, "edges": emode, "coarse_L": L, "fine_M": M,
                                       "precondition_E_M_in_E_L": ("guaranteed (co-membership)" if emode == "comembership"
                                                                   else "holds" if pre[(L, M)][1] == 0 else
                                                                   f"fails ({pre[(L, M)][1]} of {pre[(L, M)][0]} M edges)"),
                                       "EXPR_M_groups_checked": checked, "violations": len(viol),
                                       "detail": " || ".join(viol)})
                sweep.append({"expr_mode": mode, "t": t, "edges": emode, "NPIP_members_expressed": mem_ex["NPIP"],
                              "TBC1D3_members_expressed": mem_ex["TBC1D3"], "U_genes_expressed": len(expressed),
                              "distinct_EXPR_groups": len(distinct), "L_groups_split": len(splits),
                              "T2_checks": nchk, "T2_violations": nviol,
                              "splits": " || ".join(splits)})
                say(f"[EXPR sweep] {mode:9s} t>={t} {emode:12s}: members NPIP {mem_ex['NPIP']}/27 TBC1D3 "
                    f"{mem_ex['TBC1D3']}/19; U expressed {len(expressed)}; distinct EXPR groups {len(distinct)}; L groups "
                    f"split {len(splits)}; T2 checks {nchk} violations {nviol}" + (f"; splits: {splits}" if splits else ""))
    write_str(f"{INT}/expr_groups.tsv", xrows)
    write_str(f"{INT}/expr_nesting.tsv", grows)
    write_str(f"{INT}/expr_T2.tsv", t2rows)
    write_str(f"{INT}/expr_sweep.tsv", sweep)
    um = []
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        rd = READS[g]
        um.append({"gene_id": g, "name": NAME[g], "family": U[g]["member_family"], "biotype": U[g]["biotype"],
                   "n_reads_any": rd["any"], "n_reads_unique": rd["unique"], "n_reads_unique_mr": rd["unique_mr"],
                   "expressed_any_ge3": "yes" if rd["any"] >= 3 else "no",
                   "expressed_unique_ge3": "yes" if rd["unique"] >= 3 else "no",
                   "expressed_unique_mr_ge3": "yes" if rd["unique_mr"] >= 3 else "no",
                   "layers": ",".join(L for L in ("P", "D", "C_L1", "Ctree_top") if g in LAYERS[L])})
    write_str(f"{INT}/member_expression.tsv", um)
    for fam in ("NPIP", "TBC1D3"):
        sub = [r for r in um if r["family"] == fam]
        say(f"[EXPR] {fam} members {len(sub)}: >=3 any {sum(r['expressed_any_ge3'] == 'yes' for r in sub)}, unique "
            f"{sum(r['expressed_unique_ge3'] == 'yes' for r in sub)}, unique_mr "
            f"{sum(r['expressed_unique_mr_ge3'] == 'yes' for r in sub)}; any<3: "
            f"{[r['name'] for r in sub if r['expressed_any_ge3'] == 'no']}")

    # ------------------------------------------------ member reconciliation vs §6jg truth (22 NPIP + 9 TBC1D3 records)
    mrec = []
    mrows = {r["gene_id"]: r for r in tsv(f"{LIGHT}/members.corrected.tsv")}
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        r = mrows[g]
        inlit = r["in_lit_truth_31"] == "yes"
        if inlit:
            why = "in §6jg truth"
        elif "readthrough" in r["member_basis"]:
            why = "readthrough record whose family part has no gene record of its own (member rule)"
        elif r["biotype"] == "protein_coding":
            why = f"coding RefSeq LOC record ('{r['description']}'); §6jg's 22 NPIP records predate it"
        else:
            why = f"{r['biotype']} ('{r['description']}'); §6jg's TBC1D3 truth is the 9 protein-coding copies" \
                if r["family"] == "TBC1D3" else f"{r['biotype']} ('{r['description']}'); not among §6jg's 22 records"
        if r["span_inside_member_record"]:
            why += f"; span lies inside member record {r['span_inside_member_record']} (same strand)"
        mrec.append({"member": r["name"], "family": r["family"], "biotype": r["biotype"], "chrom": r["chrom"],
                     "in_6jg_truth": "yes" if inlit else "no", "reason": why})
    write_str(f"{INT}/member_reconciliation.tsv", mrec)
    for fam in ("NPIP", "TBC1D3"):
        say(f"[members] {fam}: {sum(1 for x in mrec if x['family'] == fam)} members; in §6jg truth "
            f"{sum(1 for x in mrec if x['family'] == fam and x['in_6jg_truth'] == 'yes')}; extra: "
            f"{[x['member'] for x in mrec if x['family'] == fam and x['in_6jg_truth'] == 'no']}")

    # ------------------------------------------------ disagreement lists
    dP, dD, dPD, dC = [], [], [], []
    for g in sorted(U, key=lambda g: (SIDE[g], NAME[g])):
        if g in LP:
            mates = [m for m in MEMB if m != g and LP.get(m) == LP[g]]
            if mates:
                not_d = [m for m in mates if not (g in LD and m in LD and LD[g] == LD[m])]
                if not_d and len(not_d) == len(mates):
                    best = max(mates, key=lambda m: P_EDGES.get(frozenset((g, m)), (0, 0, 0))[0])
                    w, i, c = P_EDGES.get(frozenset((g, best)), (float("nan"),) * 3)
                    ws = [P_EDGES[frozenset((g, m))] for m in mates if frozenset((g, m)) in P_EDGES]
                    dP.append({"gene": NAME[g], "biotype": U[g]["biotype"], "chrom": U[g]["chrom"],
                               "is_member": U[g]["is_member"], "P_group": LP[g],
                               "D_status": LD.get(g, "not in D universe (outside both catalogs)"),
                               "D_folded_into": U[g]["D_folded_into"],
                               "member_D_group": ";".join(sorted({LD[m] for m in mates if m in LD})),
                               "n_member_mates": len(mates), "n_member_P_edges": len(ws),
                               "best_member": NAME[best], "blastp_identity": fmt(i), "coverage_longer": fmt(c),
                               "weight": fmt(w),
                               "identity_range": f"{min(x[1] for x in ws):.3f}-{max(x[1] for x in ws):.3f}" if ws else "",
                               "coverage_range": f"{min(x[2] for x in ws):.3f}-{max(x[2] for x in ws):.3f}" if ws else ""})
        if g in LD:
            mates = [m for m in MEMB if m != g and LD.get(m) == LD[g]]
            if mates and not all(g in LP and m in LP and LP[g] == LP[m] for m in mates):
                pm = [m for m in mates if g in LP and m in LP and LP[g] == LP[m]]
                if pm:
                    continue  # P also groups g with some member: not a D-only co-membership
                es = [(D_EDGES[frozenset((g, m))], m) for m in mates if frozenset((g, m)) in D_EDGES]
                best = max(es)[1] if es else None
                es_any = [(D_EDGES[frozenset((g, m))], m) for m in MEMB if m != g and frozenset((g, m)) in D_EDGES]
                best_any = max(es_any) if es_any else None
                grp_genes = [h for h in LD if h != g and LD[h] == LD[g]]
                ea = [(D_EDGES[frozenset((g, h))], h) for h in grp_genes if frozenset((g, h)) in D_EDGES]
                bany = max(ea)[1] if ea else None
                rt = lambda h: bool(U[h]["D_folded_into"]) or "readthrough" in U[h]["description"]  # noqa: E731
                flags = []
                if U[g]["D_folded_into"]:
                    flags.append("gene folded into another locus")
                if "readthrough" in U[g]["description"]:
                    flags.append("gene is a readthrough record")
                if all(rt(m) for m in mates):
                    flags.append("its only member mates are readthrough/folded records")
                if not es and bany is not None and rt(bany):
                    flags.append(f"no direct member edge; strongest group edge is to readthrough/folded {NAME[bany]}")
                dD.append({"gene": NAME[g], "gene_id": g, "biotype": U[g]["biotype"], "chrom": U[g]["chrom"],
                           "is_member": U[g]["is_member"], "D_group": LD[g], "folded_into": U[g]["D_folded_into"],
                           "P_status": LP.get(g, "not in P universe (non-coding / r2-excluded)"),
                           "n_member_mates": len(mates), "n_member_D_edges_same_group": len(es),
                           "best_member_same_D_group": NAME[best] if best else "", "best_member_id": best or "",
                           "D_weight_same_group": fmt(max(es)[0]) if es else "",
                           "best_member_any_D_group": f"{NAME[best_any[1]]} ({best_any[0]:.3f}, {LD.get(best_any[1])})"
                           if best_any else "",
                           "best_group_partner_if_no_member_edge": "" if es or bany is None else
                           f"{NAME[bany]} ({max(ea)[0]:.3f})",
                           "readthrough_fold_flags": "; ".join(flags)})
    keyrow = {}
    for r in D_EDGE_ROWS:
        keyrow[frozenset((r["u_gene_id"], r["v_gene_id"]))] = (r["catalog"], r["u_key"], r["v_key"], r["u_gene_id"])
    want = []
    for r in dD:
        k = frozenset((r["gene_id"], r["best_member_id"]))
        if r["best_member_id"] and k in keyrow:
            c_, uk, vk, _ = keyrow[k]
            want.append((c_, uk, vk))
    st = dna_pair_stats(want)
    n_ok = n_chk = 0
    for r in dD:
        k = frozenset((r["gene_id"], r["best_member_id"]))
        r["DNA_identity"], r["DNA_cov_longer_exonic"], r["id_x_cov_equals_weight"] = "", "", ""
        if r["best_member_id"] and k in keyrow:
            c_, uk, vk, _ = keyrow[k]
            s = st.get((c_, frozenset((uk, vk))))
            if s:
                r["DNA_identity"], r["DNA_cov_longer_exonic"] = fmt(s[0]), fmt(s[1])
                eq = abs(s[0] * s[1] - float(r["D_weight_same_group"])) < 5e-4
                r["id_x_cov_equals_weight"] = "yes" if eq else f"no ({s[0] * s[1]:.3f})"
                n_chk += 1
                n_ok += eq
    say(f"[disagree] D-only genes {len(dD)}; DNA identity x coverage (best member edge in the same D group) reproduces the "
        f"D weight for {n_ok}/{n_chk}")
    for g in sorted(MEMB, key=lambda g: (SIDE[g], NAME[g])):
        if g in LP and g in LD:
            mp = {NAME[m] for m in U if m != g and m in LP and m in LD and LP[m] == LP[g]}
            md = {NAME[m] for m in U if m != g and m in LP and m in LD and LD[m] == LD[g]}
            if mp != md:
                dPD.append({"member": NAME[g], "P_group": LP[g], "D_group": LD[g],
                            "P_only_comates": ",".join(sorted(mp - md)), "D_only_comates": ",".join(sorted(md - mp)),
                            "member_comates_P_only": ",".join(sorted(x for x in mp - md if "gene-" + x in MEMB)),
                            "member_comates_D_only": ",".join(sorted(x for x in md - mp if "gene-" + x in MEMB))})
    for fam in ("NPIP", "TBC1D3"):
        cnt = collections.Counter(LD[m] for m in MEMB if SIDE[m] == fam and m in LD)
        main_g = cnt.most_common(1)[0][0]
        for m in sorted(MEMB, key=lambda x: NAME[x]):
            if SIDE[m] == fam and m in LD and LD[m] != main_g:
                dPD.append({"member": NAME[m], "P_group": LP.get(m, "not in P universe"), "D_group": LD[m],
                            "P_only_comates": "", "D_only_comates": f"outside the family's main D group {main_g} "
                                                                    f"(folded into {U[m]['D_folded_into'] or '-'})",
                            "member_comates_P_only": "", "member_comates_D_only": ""})
    for cl in ("C_L1", "C_mid", "C_fine", "Ctree_top", "Ctree_min"):
        LCL = LAYERS[cl]
        for other_name, LO in (("P", LP), ("D", LD)):
            S = set(LCL) & set(LO)
            bad = sorted(pairs_of(LCL, S) - pairs_of(LO, S))
            dC.append({"clade_level": cl, "vs": other_name, "clade_pairs_on_shared_genes": len(pairs_of(LCL, S)),
                       "clade_pairs_split_by_other": len(bad),
                       "examples": "; ".join(f"{NAME[a]}-{NAME[b]}" for a, b in bad[:10]),
                       "clade_genes_outside_other_universe": ",".join(sorted(NAME[g] for g in set(LCL) - set(LO)))})
    write_str(f"{INT}/disagree_P_not_D.tsv", dP)
    write_str(f"{INT}/disagree_D_not_P.tsv", [{k: v for k, v in r.items() if k not in ("gene_id", "best_member_id")}
                                          for r in dD])
    write_str(f"{INT}/disagree_members_P_vs_D.tsv", dPD)
    write_str(f"{INT}/disagree_clades.tsv", dC)
    say(f"[disagree] P-with-member-not-D {len(dP)}: {[r['gene'] for r in dP]}")
    say(f"[disagree] members with different P vs D co-mates / off-main D group {len(dPD)}")
    for r in dC:
        say(f"[disagree] {r['clade_level']} vs {r['vs']}: clade pairs {r['clade_pairs_on_shared_genes']}, split "
            f"{r['clade_pairs_split_by_other']}; clade genes outside {r['vs']}: {r['clade_genes_outside_other_universe']}")
    say.dump(f"{INT}/analysis.out")


# ============================================================================================================== nested edge-test lattice (was lattice_*.py)
def cmd_lattice_edges(args):
    """was bench/layer_order/lattice_edges.py

    Nested edge-test lattice, step 1: the unified edge table and the L0 closure.

    Outputs (lattice/):
      nodes.tsv          one row per gene of V = U (68, integrate_slim universe) closed under L0 edges
      edges.tsv          one row per gene pair inside V with ANY evidence (blastp HSP pair, catalog PAF record, S2 edge)
      edges_all_c2.tsv   every clause-2-approximation DNA edge (primary or loose) in both catalogs (L1 before the closure cut)
      closure.tsv        BFS rounds of the L0 closure; protein-only and DNA-only closures of U
      edges_build.out    validation lines (E1 dump reproduction, S1 dump reproduction, §6ko edge-set reproduction)

    Evidence and provenance per edge (column prefix):
      p_*   protein, light/work/P/blastp.tsv (outfmt 'qseqid sseqid nident length qstart qend sstart send bitscore', searched
            with -evalue 1e-5 against the 20,088-protein database; 4,430 proteins searched). Per ordered (query, subject) the
            shipped bench/protein_families.edges_from rule (greedy by bitscore, non-overlapping HSPs on the LONGER protein;
            identity = sum nident / sum length; coverage = merged HSP span on the longer protein / its length). The pair keeps
            the best-weight qualifying direction (coverage >= 0.30), else the best-coverage direction (reported, not
            qualifying). e-value: not saved per HSP (only the search cutoff 1e-5 is known). p_cov_union_longer /
            p_qualifies_union: the §0★★★.1 union-cover form (all HSP intervals on the longer protein merged; >= 0.30 in either
            direction); used only for a closure row, not by any level.
      d_*   DNA, the two E1 catalogs' all-vs-all gene-body PAFs (minimap2 -x asm20 -c -X -N 50 -p 0.1; -X = one direction per
            pair). Keys = gene-body spans (1-based inclusive). Exon unions: c15_17_22 light/work/refseq/exons.tsv, c16_19_20
            o1_falsemerge/lit/aj_ho/refseq/nodes.tsv; a record without exons is one span exon (mcl_families --exonless-span).
            d_e1_*: src/rustle/vg_family/annotation_families.rs graph_from_paf_loci re-implemented (records >= 300 bp and
            identity >= 0.70; pooled identity; cov_longer = merged aligned span on the gene with the longer exon union / that
            exon-union length, capped at 1; gates: >= 1 exonic bp covered on the longer gene and >= 1 exon-to-exon bp on one
            record; edge iff cov_longer >= 0.30). d_shared_exon_frac = max over those records of min(exonic bp of A in the
            record, exonic bp of B in the record) / min(exon-union length A, B)  (= --min-shared-exon-frac numerator and
            denominator). d_e1_identity_gapexcl = sum matches / sum CIGAR M over the same records (indels excluded, the SEDEF
            fracMatch analogue).
            d_c2_*: clause 2 (§0★★) APPROXIMATED on the same PAF: gene-body = chains built exactly as
            bench/guided_pipeline.gene_body_chains on each direction (records of the pair grouped by strand; gap <= query body,
            span <= 2 x query body, query order), a chain passes iff identity >= 0.80 and aligned query bp >= 0.50 x
            min(body u, body v) (clause 2's literal 'shorter body'; nodes have fixed extents); exon = merged query intervals of
            records with identity >= 0.80 cover >= 0.50 of the query gene's exon union (proxy for 'spliced transcript aligns
            at >= 0.80 over >= 0.50'; record identity includes introns).
            PRIMARY (correction pass 2026-09-16): both disjuncts also carry the two target-side requirements of the shipped
            clause 2 (bench/denovo_shared_def.py cmd_families; §0★★★.1 'its target overlaps v's exons'): (1) the passing
            chain's target span, or each exon-proxy record's target interval, overlaps >= 1 exonic base of v; (2) the shipped
            strand check: when both u and v are spliced (exon union >= 2 blocks), v's strand must equal u's strand for a '+'
            record/chain and the opposite strand for a '-' one (gene bodies are extracted on the genomic + strand).
            Variants: d_c2nostrand_* (requirement 1 only), d_c2exontgt_* (requirement 1 on the exon proxy only, gene-body chain
            unrestricted), d_c2loose_* (neither requirement: the 17:03 build's primary).
            d_c2x_*: the chains with the guided finder's own denominator min(query body, extrapolated target span), plus
            requirements 1-2 — on a gene-body PAF the extrapolation is clipped at the target BODY end, so a few-hundred-bp
            overlap at a body edge passes (e.g. NPIPA8-PKD1P1, 600 bp at identity 1.0); variant only. d_c2x_gb_chain_raw = a
            chain exists under that denominator without requirements 1-2, checked against the shipped
            gene_body_chains by lattice_check_c2.py. d_c2xloose_approx (raw chain OR loose exon proxy) is used only to define V.
            d_w98_gapexcl / d_w98_gapincl: §0★★★.1 single-record w_98 = max identity over ALL PAF records r of the pair with
            sx(r) = min(exonic bp of A in r, exonic bp of B in r) >= 0.30 x min(exon-union A, B); NA when no record witnesses.
            d_shared_exon_frac_allrec: f_ex as a maximum over ALL records (the definition's form; primary f_ex keeps the
            shipped E1-record restriction).
      s2_*  heavy/S2.edges.tsv (SD98 map-back; max_identity = max PAF col10/col11 over mappings; n_shared_exons).
      ctree_* clause-5 split system (lo_analysis.c_tree on light/C.supported_clades.tsv): annotation only.
    """
    LIGHT, HEAVY, OUT = LC.LIGHT, LC.HEAVY, LC.OUT
    DUMP, DUMP_S1 = LC.DUMP, LC.DUMP_S1
    T0 = time.time()
    say = Log()

    # ---------------------------------------------------------------------------------------- genes, catalogs
    ctx = catalog_context(say)
    genes, by_coord, exons_tsv, U = ctx.genes, ctx.by_coord, ctx.exons_tsv, ctx.U
    key2genes, key_blocks, key_cat, key_exlen, gene_key = (ctx.key2genes, ctx.key_blocks, ctx.key_cat, ctx.key_exlen,
                                                           ctx.gene_key)

    FLIP = {"+": "-", "-": "+"}

    pair_recs = collections.defaultdict(list)
    n_lines = 0
    for cat, path in PAF.items():
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 11 or f[0] == f[5]:
                    continue
                n_lines += 1
                qk, tk = key_of(f[0]), key_of(f[5])
                cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
                mlen = sum(int(n) for n, op in CIG.findall(cg) if op in "M=X") if cg else None
                rec = {"ql": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "tl": int(f[6]), "ts": int(f[7]),
                       "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10]), "mlen": mlen}
                if qk <= tk:
                    pair_recs[(qk, tk)].append((rec, True))
                else:
                    pair_recs[(tk, qk)].append((rec, False))
    say(f"[paf] records (non-self) {n_lines}; key pairs {len(pair_recs)}; {time.time() - T0:.0f}s")

    def dna_attrs(a, b, recs, sa, sb):
        """a <= b keys; recs: list of (rec, a_is_query); sa, sb: strands of the genes on keys a and b (the strand check of
        clause 2 depends on them, everything else does not)."""
        blocks_a, blocks_b = key_blocks.get(a), key_blocks.get(b)
        out = {"d_paf_records": len(recs)}
        # ---- E1 (graph_from_paf_loci, deferred-pair branch)
        aiv, biv = [], []
        nm = bl = ml = 0
        ml_ok = True
        exon_exon = 0
        n_e1 = 0
        for r, aq in recs:
            if r["bl"] < 300 or r["nm"] / max(r["bl"], 1) < 0.70:
                continue
            n_e1 += 1
            (as_, ae, bs, be) = (r["qs"], r["qe"], r["ts"], r["te"]) if aq else (r["ts"], r["te"], r["qs"], r["qe"])
            aiv.append((as_, ae))
            biv.append((bs, be))
            nm += r["nm"]
            bl += r["bl"]
            if r["mlen"] is None:
                ml_ok = False
            else:
                ml += r["mlen"]
            ax = exonic_bases_in(blocks_a, a[1], as_, ae) if blocks_a else 0
            bx = exonic_bases_in(blocks_b, b[1], bs, be) if blocks_b else 0
            exon_exon = max(exon_exon, min(ax, bx))
        la = recs[0][0]["ql"] if recs[0][1] else recs[0][0]["tl"]
        lb = recs[0][0]["tl"] if recs[0][1] else recs[0][0]["ql"]
        da = key_exlen.get(a, la)
        db = key_exlen.get(b, lb)
        # ---- §0★★★.1 single-record w_98 and all-record f_ex (every PAF record of the pair, no E1 filter)
        den_x = min(da, db)
        w98_gi = w98_ge = None
        w98_gi_e1 = w98_ge_e1 = None
        sx_all = 0
        for r, aq in recs:
            (as_, ae, bs, be) = (r["qs"], r["qe"], r["ts"], r["te"]) if aq else (r["ts"], r["te"], r["qs"], r["qe"])
            ax = exonic_bases_in(blocks_a, a[1], as_, ae) if blocks_a else 0
            bx = exonic_bases_in(blocks_b, b[1], bs, be) if blocks_b else 0
            sx = min(ax, bx)
            sx_all = max(sx_all, sx)
            if sx >= 0.30 * den_x:
                gi = r["nm"] / max(r["bl"], 1)
                ge = r["nm"] / r["mlen"] if r["mlen"] else None
                w98_gi = gi if w98_gi is None else max(w98_gi, gi)
                if ge is not None:
                    w98_ge = ge if w98_ge is None else max(w98_ge, ge)
                if r["bl"] >= 300 and gi >= 0.70:
                    w98_gi_e1 = gi if w98_gi_e1 is None else max(w98_gi_e1, gi)
                    if ge is not None:
                        w98_ge_e1 = ge if w98_ge_e1 is None else max(w98_ge_e1, ge)
        out["d_w98_gapexcl"] = w98_ge
        out["d_w98_gapincl"] = w98_gi
        out["_w98_e1"] = (w98_ge_e1, w98_gi_e1)  # E1-record-restricted w_98, for the build log only (not written)
        out["d_shared_exon_frac_allrec"] = sx_all / max(1, min(da, db))
        out["d_e1_records"] = n_e1
        out["d_exon_union_bp_a"] = da
        out["d_exon_union_bp_b"] = db
        out["d_body_bp_a"] = la
        out["d_body_bp_b"] = lb
        out["d_shared_exon_bp"] = exon_exon
        out["d_shared_exon_frac"] = exon_exon / max(1, min(da, db))
        if n_e1:
            longer_is_a = da >= db
            gk, iv, den, blocks = (a, aiv, da, blocks_a) if longer_is_a else (b, biv, db, blocks_b)
            m = merge(iv)
            covered = sum(exonic_bases_in(blocks, gk[1], s, e) for s, e in m) if blocks else 0
            numer = sum(e - s for s, e in m)
            cov = min(1.0, numer / max(den, 1))
            ident = nm / max(bl, 1)
            gate = blocks is not None and covered >= 1 and exon_exon >= 1
            out.update({"d_e1_identity": ident, "d_e1_identity_gapexcl": (nm / ml if ml_ok and ml else None),
                        "d_e1_cov_longer": cov, "d_e1_cov_gene": "A" if longer_is_a else "B",
                        "d_e1_exonic_bp_covered_longer": covered, "d_e1_gate_exonic": gate,
                        "d_e1_edge": gate and cov >= 0.30, "d_e1_weight": ident * cov if gate and cov >= 0.30 else None})
        else:
            out.update({"d_e1_identity": None, "d_e1_identity_gapexcl": None, "d_e1_cov_longer": None, "d_e1_cov_gene": "",
                        "d_e1_exonic_bp_covered_longer": None, "d_e1_gate_exonic": False, "d_e1_edge": False,
                        "d_e1_weight": None})
        # ---- clause 2 approximation
        # chain candidates: (passes, chain identity >= 0.80, aligned frac, identity, direction); '' direction = no eligible chain
        NEG = (False, False, -1.0, None, "")
        best = {n: NEG for n in ("gb", "gx", "gb_ns", "gb_loose", "gx_loose")}
        best_ex = {n: (False, 0.0, "") for n in ("ex", "ex_ns", "ex_loose")}
        both_spliced = len(blocks_a or ()) >= 2 and len(blocks_b or ()) >= 2
        for direction in ("A->B", "B->A"):
            u_strand, v_strand = (sa, sb) if direction == "A->B" else (sb, sa)
            tkey, tblocks = (b, blocks_b) if direction == "A->B" else (a, blocks_a)

            def strand_ok(rec_strand, u_strand=u_strand, v_strand=v_strand):
                # bench/denovo_shared_def.py cmd_families: orient = u strand for a '+' hit, flipped for '-'; skip when both
                # nodes are spliced and v's strand differs from orient
                if not both_spliced:
                    return True
                return v_strand == (u_strand if rec_strand == "+" else FLIP.get(u_strand, u_strand))

            view = []
            for r, aq in recs:
                q_is_a = aq
                if direction == "A->B":
                    if q_is_a:
                        v = dict(qs=r["qs"], qe=r["qe"], ts=r["ts"], te=r["te"], Lq=r["ql"], clen=r["tl"])
                    else:
                        v = dict(qs=r["ts"], qe=r["te"], ts=r["qs"], te=r["qe"], Lq=r["tl"], clen=r["ql"])
                else:
                    if q_is_a:
                        v = dict(qs=r["ts"], qe=r["te"], ts=r["qs"], te=r["qe"], Lq=r["tl"], clen=r["ql"])
                    else:
                        v = dict(qs=r["qs"], qe=r["qe"], ts=r["ts"], te=r["te"], Lq=r["ql"], clen=r["tl"])
                v.update(strand=r["strand"], nm=r["nm"], bl=r["bl"])
                v["touch"] = bool(tblocks) and exonic_bases_in(tblocks, tkey[1], v["ts"], v["te"]) > 0
                view.append(v)
            Lq, clen = view[0]["Lq"], view[0]["clen"]
            for strand in ("+", "-"):
                rs = [v for v in view if v["strand"] == strand]
                if not rs:
                    continue
                sok = strand_ok(strand)
                for grp in chain_groups(rs, Lq):
                    okx, idn, fracx, ok, frac, cts, cte = chain_stats(grp, Lq, clen)
                    touch = bool(tblocks) and exonic_bases_in(tblocks, tkey[1], cts, cte) > 0
                    cands = [("gb_loose", (ok, idn >= 0.80, frac, idn, direction)),
                             ("gx_loose", (okx, idn >= 0.80, fracx, idn, direction))]
                    if touch:
                        cands.append(("gb_ns", (ok, idn >= 0.80, frac, idn, direction)))
                        if sok:
                            cands += [("gb", (ok, idn >= 0.80, frac, idn, direction)),
                                      ("gx", (okx, idn >= 0.80, fracx, idn, direction))]
                    for name, cand in cands:
                        if cand[:3] > best[name][:3]:
                            best[name] = cand
            qkey = a if direction == "A->B" else b
            qblocks = blocks_a if direction == "A->B" else blocks_b
            qlen_ex = da if direction == "A->B" else db
            for name, keep in (("ex_loose", lambda v: True), ("ex_ns", lambda v: v["touch"]),
                               ("ex", lambda v: v["touch"] and strand_ok(v["strand"]))):
                m = merge([(v["qs"], v["qe"]) for v in view if v["nm"] / v["bl"] >= 0.80 and keep(v)])
                exb = sum(exonic_bases_in(qblocks, qkey[1], s, e) for s, e in m) if qblocks else 0
                fr = exb / max(1, qlen_ex)
                if (fr >= 0.50, fr) > (best_ex[name][0], best_ex[name][1]):
                    best_ex[name] = (fr >= 0.50, fr, direction)

        def frac_or_na(x):
            return None if x < 0 else x
        out.update({"d_both_spliced": both_spliced,
                    "d_c2_genebody": best["gb"][0], "d_c2_gb_best_frac": frac_or_na(best["gb"][2]),
                    "d_c2_gb_chain_identity": best["gb"][3], "d_c2_gb_direction": best["gb"][4],
                    "d_c2_exon": best_ex["ex"][0], "d_c2_exon_best_frac": frac_or_na(best_ex["ex"][1]),
                    "d_c2_exon_direction": best_ex["ex"][2], "d_c2_approx": best["gb"][0] or best_ex["ex"][0],
                    "d_c2x_genebody": best["gx"][0], "d_c2x_gb_best_frac": frac_or_na(best["gx"][2]),
                    "d_c2x_approx": best["gx"][0] or best_ex["ex"][0], "d_c2x_gb_chain_raw": best["gx_loose"][0],
                    "d_c2nostrand_approx": best["gb_ns"][0] or best_ex["ex_ns"][0],
                    "d_c2exontgt_approx": best["gb_loose"][0] or best_ex["ex_ns"][0],
                    "d_c2loose_genebody": best["gb_loose"][0], "d_c2loose_gb_best_frac": frac_or_na(best["gb_loose"][2]),
                    "d_c2loose_gb_chain_identity": best["gb_loose"][3], "d_c2loose_exon": best_ex["ex_loose"][0],
                    "d_c2loose_exon_best_frac": frac_or_na(best_ex["ex_loose"][1]),
                    "d_c2loose_approx": best["gb_loose"][0] or best_ex["ex_loose"][0],
                    "d_c2xloose_approx": best["gx_loose"][0] or best_ex["ex_loose"][0]})
        return out

    dna = {}
    for (a, b), recs in pair_recs.items():
        dna[(a, b)] = dna_attrs(a, b, recs, genes[key2genes[a][0]]["strand"], genes[key2genes[b][0]]["strand"])

    def _n(f):
        return sum(1 for v in dna.values() if v[f])

    say(f"[dna] key pairs evaluated {len(dna)} (strands of each key's first gene; gene pairs on the {sum(1 for v in key2genes.values() if len(v) > 1)} "
        f"multi-gene keys are re-evaluated with their own strands below); E1 edges {_n('d_e1_edge')}")
    say(f"[dna] clause-2 approx, PRIMARY (v-exon overlap + strand check, shorter-body denominator): {_n('d_c2_approx')} "
        f"(gene-body {_n('d_c2_genebody')}; exon proxy {_n('d_c2_exon')}); v-exon overlap without strand check "
        f"{_n('d_c2nostrand_approx')}; v-exon overlap on the exon proxy only {_n('d_c2exontgt_approx')}; neither requirement "
        f"(17:03 primary) {_n('d_c2loose_approx')} (gene-body {_n('d_c2loose_genebody')}; exon proxy {_n('d_c2loose_exon')})")
    say(f"[dna] finder denominator: with both requirements {_n('d_c2x_approx')} (gene-body {_n('d_c2x_genebody')}); raw chains "
        f"{_n('d_c2x_gb_chain_raw')}; raw chains OR loose exon proxy (17:03 c2x) {_n('d_c2xloose_approx')}")
    _w = collections.Counter()
    for v in dna.values():
        for i, f in enumerate(("d_w98_gapexcl", "d_w98_gapincl")):
            allr = v[f] is not None and v[f] >= 0.98
            e1r = v["_w98_e1"][i] is not None and v["_w98_e1"][i] >= 0.98
            _w[(f, allr, e1r)] += 1
    say(f"[dna] w_98 >= 0.98 key pairs (all records / E1 records only / decisions differing): gap-excl "
        f"{_w[('d_w98_gapexcl', True, True)] + _w[('d_w98_gapexcl', True, False)]} / "
        f"{_w[('d_w98_gapexcl', True, True)] + _w[('d_w98_gapexcl', False, True)]} / "
        f"{_w[('d_w98_gapexcl', True, False)] + _w[('d_w98_gapexcl', False, True)]}; gap-incl "
        f"{_w[('d_w98_gapincl', True, True)] + _w[('d_w98_gapincl', True, False)]} / "
        f"{_w[('d_w98_gapincl', True, True)] + _w[('d_w98_gapincl', False, True)]} / "
        f"{_w[('d_w98_gapincl', True, False)] + _w[('d_w98_gapincl', False, True)]}; {time.time() - T0:.0f}s")

    # ---- validation against the dumps (E1 = D layer graph; S1 = shipped shared-exon 0.30)
    for tag, dumps, rule in (("E1", DUMP, lambda v: v["d_e1_edge"]),
                             ("S1 (E1 + shared-exon >= 0.30)", DUMP_S1, lambda v: v["d_e1_edge"] and v["d_shared_exon_frac"] >= 0.30)):
        for cat, path in dumps.items():
            dump = {}
            for line in open(path):
                u, v, w = line.rstrip("\n").split("\t")
                ku, kv = key_of(u), key_of(v)
                dump[(min(ku, kv), max(ku, kv))] = float(w)
            mine = {k: v for k, v in dna.items() if key_cat.get(k[0]) == cat and rule(v)}
            both = set(dump) & set(mine)
            wdiff = sum(1 for k in both if abs(round(mine[k]["d_e1_weight"], 6) - dump[k]) > 1.5e-6)
            say(f"[validate {tag} {cat}] dump edges {len(dump)}; recomputed {len(mine)}; common {len(both)}; dump-only "
                f"{len(set(dump) - set(mine))}; recomputed-only {len(set(mine) - set(dump))}; weight mismatches (>1.5e-6) {wdiff}")
            for k in sorted(set(dump) - set(mine))[:3]:
                say(f"   dump-only example {k} w {dump[k]} recomputed {dna.get(k)}")
            for k in sorted(set(mine) - set(dump))[:3]:
                say(f"   recomputed-only example {k} {mine[k]}")

    # ------------------------------------------------------------------------------------------------ gene-level DNA pairs
    def gpair(x, y):
        return (x, y) if x < y else (y, x)

    def same_locus(x, y):
        gx, gy = genes[x], genes[y]
        return gx["chrom"] == gy["chrom"] and int(gx["start0"]) < int(gy["end"]) and int(gy["start0"]) < int(gx["end"])

    gdna = {}
    n_restrand = 0
    for (a, b), v0 in dna.items():
        for x in key2genes.get(a, []):
            for y in key2genes.get(b, []):
                if x != y:
                    sx_, sy_ = genes[x]["strand"], genes[y]["strand"]
                    if (sx_, sy_) != (genes[key2genes[a][0]]["strand"], genes[key2genes[b][0]]["strand"]):
                        v = dna_attrs(a, b, pair_recs[(a, b)], sx_, sy_)  # strand check with this gene pair's strands
                        n_restrand += 1
                    else:
                        v = v0
                    k = gpair(x, y)
                    vv = {f: val for f, val in v.items() if not f.startswith("_")}
                    if k[0] != x:  # A/B oriented by key; relabel to gene order
                        for f1, f2 in (("d_exon_union_bp_a", "d_exon_union_bp_b"), ("d_body_bp_a", "d_body_bp_b")):
                            vv[f1], vv[f2] = v[f2], v[f1]
                        vv["d_e1_cov_gene"] = {"A": "B", "B": "A"}.get(v["d_e1_cov_gene"], "")
                        for f in ("d_c2_gb_direction", "d_c2_exon_direction"):  # (d_c2x_* carry no direction)
                            vv[f] = {"A->B": "B->A", "B->A": "A->B"}.get(v[f], "")
                    vv["d_catalog"] = key_cat[a]
                    gdna[k] = vv
    say(f"[dna] gene pairs with PAF records {len(gdna)}; gene pairs re-evaluated with their own strands {n_restrand}")

    # ------------------------------------------------------------------------------------------------ protein
    pidx = {r["pid"]: r for r in tsv(f"{LIGHT}/work/P/proteins.index.tsv")}
    plen = {p: int(r["length_aa"]) for p, r in pidx.items()}
    p2g = {p: r["gene_id"] for p, r in pidx.items()}
    searched = {x.strip() for x in open(f"{LIGHT}/work/P/searched.txt") if x.strip()}
    say(f"[protein] proteins {len(pidx)} (genes {len(set(p2g.values()))}); searched {len(searched)}")
    hs = collections.defaultdict(list)
    for line in open(f"{LIGHT}/work/P/blastp.tsv"):
        q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip("\n").split("\t")
        if q != s:
            hs[(q, s)].append((float(bits), int(nid), int(ln), int(q0) - 1, int(q1), int(s0) - 1, int(s1)))
    say(f"[protein] ordered HSP pairs {len(hs)}; {time.time() - T0:.0f}s")
    prot = {}
    for (q, s), rows in hs.items():
        longer_is_q = plen[q] >= plen[s]
        taken, L, N = [], 0, 0
        for bits, nid, ln, q0, q1, s0, s1 in sorted(rows, reverse=True):  # shipped greedy order
            iv = (q0, q1) if longer_is_q else (s0, s1)
            if any(iv[0] < y and x < iv[1] for x, y in taken):
                continue
            taken.append(iv)
            L += ln
            N += nid
        cov = sum(y - x for x, y in merge(taken)) / max(plen[q], plen[s])
        idn = N / L if L else 0.0
        qual = cov >= 0.30
        w = idn * min(cov, 1.0)
        # §0★★★.1 union-cover form: every HSP interval on the longer protein, merged (monotone in the HSP set)
        cov_u = sum(y - x for x, y in merge([(q0, q1) if longer_is_q else (s0, s1) for _b, _n, _l, q0, q1, s0, s1 in rows])) / max(plen[q], plen[s])
        k = gpair(p2g[q], p2g[s])
        cur = prot.get(k)
        mx = max(r[0] for r in rows)
        cu = max(cov_u, cur["p_cov_union_longer"]) if cur else cov_u
        cand = {"p_qualifies_6ko": qual, "p_aa_identity": idn, "p_cov_longer": cov, "p_weight": w if qual else None,
                "p_max_bitscore": max(mx, cur["p_max_bitscore"]) if cur else mx,
                "p_directions": ",".join(sorted(set((cur["p_directions"].split(",") if cur else []) + [f"{pidx[q]['name']}>{pidx[s]['name']}"]))),
                "p_cov_union_longer": cu, "p_qualifies_union": cu >= 0.30}
        if cur is None:
            prot[k] = cand
        else:
            better = (qual, w if qual else cov) > (cur["p_qualifies_6ko"], (cur["p_weight"] or 0) if cur["p_qualifies_6ko"] else cur["p_cov_longer"])
            if better:
                prot[k] = cand
            else:
                cur["p_max_bitscore"], cur["p_directions"] = cand["p_max_bitscore"], cand["p_directions"]
                cur["p_cov_union_longer"], cur["p_qualifies_union"] = cu, cu >= 0.30
    del hs
    say(f"[protein] gene pairs with HSPs {len(prot)}; §6ko-qualifying (shipped greedy cover) "
        f"{sum(1 for v in prot.values() if v['p_qualifies_6ko'])}; union-cover qualifying {sum(1 for v in prot.values() if v['p_qualifies_union'])} "
        f"(union-only {sum(1 for v in prot.values() if v['p_qualifies_union'] and not v['p_qualifies_6ko'])}, greedy-only "
        f"{sum(1 for v in prot.values() if v['p_qualifies_6ko'] and not v['p_qualifies_union'])}); {time.time() - T0:.0f}s")
    from truth import edges_from, pair_hsps  # (was protein_families; wave 7)

    E_ship = edges_from(pair_hsps(f"{LIGHT}/work/P/blastp.tsv", plen), plen, 0.0)
    ship_g = {gpair(p2g[u], p2g[v]) for u, v in E_ship}
    mine_g = {k for k, v in prot.items() if v["p_qualifies_6ko"]}
    say(f"[validate protein] shipped edges_from {len(E_ship)} edges; recomputed qualifying gene pairs {len(mine_g)}; "
        f"symmetric difference {len(ship_g ^ mine_g)}")
    pe = {gpair(r["u_gene_id"], r["v_gene_id"]): r for r in tsv(f"{LIGHT}/P.edges.tsv")}
    pdiff = sum(1 for k, r in pe.items() if k not in prot or abs(prot[k]["p_aa_identity"] - float(r["identity"])) > 6e-5
                or abs(prot[k]["p_cov_longer"] - float(r["coverage_longer"])) > 6e-5)
    say(f"[validate protein] light/P.edges.tsv rows {len(pe)}; identity/coverage mismatches (> 6e-5) {pdiff}")

    # ------------------------------------------------------------------------------------------------ S2
    s2g = {}
    for r in tsv(f"{HEAVY}/S2.genes.tsv"):
        cand = [g["gene_id"] for g in by_coord.get((r["chrom"], int(r["start"]) + 1, int(r["end"])), []) if g["name"] == r["name"]]
        if len(cand) == 1:
            s2g[r["gene_id"]] = cand[0]
    s2 = {}
    for r in tsv(f"{HEAVY}/S2.edges.tsv"):
        if r["gene_a"] in s2g and r["gene_b"] in s2g:
            s2[gpair(s2g[r["gene_a"]], s2g[r["gene_b"]])] = r
    say(f"[s2] genes mapped {len(s2g)} of {len(tsv(f'{HEAVY}/S2.genes.tsv'))}; edges mapped {len(s2)}")

    # ------------------------------------------------------------------------------------------------ closure over L0
    ADJ = {n: collections.defaultdict(set) for n in ("U", "0", "P", "D", "0loose", "Dloose", "0union")}

    def link(n, k):
        ADJ[n][k[0]].add(k[1])
        ADJ[n][k[1]].add(k[0])

    for k, v in prot.items():
        sl = same_locus(*k)
        if v["p_qualifies_6ko"]:
            link("U", k)  # V: same-locus links included
            if not sl:
                for n in ("0", "P", "0loose"):
                    link(n, k)
        if v["p_qualifies_union"] and not sl:
            link("0union", k)
    for k, v in gdna.items():
        sl = same_locus(*k)
        if v["d_c2loose_approx"] or v["d_c2xloose_approx"] or v["d_e1_edge"]:  # every strict variant is a subset of these
            link("U", k)
        if sl:
            continue
        if v["d_c2_approx"]:
            for n in ("0", "D", "0union"):
                link(n, k)
        if v["d_c2loose_approx"]:
            for n in ("0loose", "Dloose"):
                link(n, k)

    def bfs(adj, seeds):
        dist = {s: 0 for s in seeds}
        frontier = list(seeds)
        rounds = []
        while frontier:
            nxt = []
            for x in frontier:
                for y in adj.get(x, ()):
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        nxt.append(y)
            if nxt:
                rounds.append(len(nxt))
            frontier = nxt
        return dist, rounds

    CL = {n: bfs(ADJ[n], set(U)) for n in ADJ}
    V, roundsU = CL["U"]
    V0, rounds = CL["0"]
    g2p = {g: p for p, g in p2g.items()}
    prev_V = {r["gene_id"] for r in tsv(f"{OUT}/pre_correction_1712/nodes.tsv")} if os.path.exists(
        f"{OUT}/pre_correction_1712/nodes.tsv") else None

    def crow_of(name, n, extra=True):
        dist, rr = CL[n]
        row = {"closure": name, "genes": len(dist), "new_genes_per_hop": ",".join(map(str, rr)),
               "never_searched_proteins": sum(1 for g in dist if g in g2p and g2p[g] not in searched)}
        if extra:
            row.update({"no_protein": sum(1 for g in dist if g not in g2p),
                        "outside_E1_catalogs": sum(1 for g in dist if g not in gene_key)})
        row["outside_V"] = sum(1 for g in dist if g not in V)
        return row

    crow = [crow_of("V: union of all reported L0 operationalisations (protein §6ko; clause-2 approx with and without the "
                    "v-exon/strand requirements, both denominators; E1; same-locus links included)", "U"),
            crow_of("primary L0 (protein §6ko greedy cover OR clause-2 approx with v-exon overlap + strand check), same-locus "
                    "links excluded", "0"),
            crow_of("protein only (§6ko greedy cover)", "P"),
            crow_of("DNA only (primary clause-2 approx)", "D"),
            crow_of("L0 with the union-cover protein test (§0★★★.1 form) OR primary clause-2 approx", "0union"),
            crow_of("pre-correction primary L0 (protein OR clause-2 approx without v-exon/strand requirements)", "0loose"),
            crow_of("DNA only, clause-2 approx without v-exon/strand requirements (pre-correction)", "Dloose")]
    for r in crow:
        say(f"[closure] {r['closure']}: {r['genes']} genes; hops {r['new_genes_per_hop']}; never-searched proteins "
            f"{r['never_searched_proteins']}; no §6ko protein {r.get('no_protein', '-')}; outside both E1 catalogs "
            f"{r.get('outside_E1_catalogs', '-')}; outside V {r['outside_V']}")
    say(f"[closure] V minus primary closure {len(set(V) - set(V0))}; V identical to the 17:03 build's V: "
        f"{None if prev_V is None else set(V) == prev_V}")
    write(f"{OUT}/closure.tsv", crow)

    # ------------------------------------------------------------------------------------------------ C_tree annotation
    CT_LAB, CT_ROWS, CT_LIT, CT_CLUSTERS = c_tree(U, {g: r["name"] for g, r in U.items()})  # (was lo_analysis)

    ct_sup = {}
    for r in CT_ROWS:
        ct_sup[(r["family"], frozenset(r["cluster_smaller_side"].split(",")))] = r["support"]
    ct_top = CT_LAB["Ctree_top"]

    def ctree_pair(x, y):
        nx_, ny_ = genes[x]["name"], genes[y]["name"]
        for fam, (L, K) in CT_CLUSTERS.items():
            if nx_ in L and ny_ in L:
                cs = sorted((s for s in K if nx_ in s and ny_ in s), key=len)
                if cs:
                    return fam, ",".join(sorted(cs[0])), ct_sup.get((fam, cs[0]), ""), ct_top.get(x) == ct_top.get(y)
                return fam, "none (only the whole family)", "", ct_top.get(x) == ct_top.get(y) and "singleton" not in ct_top.get(x, "")
        return "", "", "", ""

    # ------------------------------------------------------------------------------------------------ write
    cols = ["gene_a", "gene_b", "name_a", "name_b", "chrom_a", "chrom_b", "same_locus",
            "p_evidence", "p_searched_a", "p_searched_b", "p_aa_identity", "p_cov_longer", "p_weight", "p_max_bitscore",
            "p_evalue", "p_qualifies_6ko", "p_aa50", "p_directions", "p_cov_union_longer", "p_qualifies_union",
            "d_evidence", "d_catalog", "d_paf_records", "d_e1_records", "d_body_bp_a", "d_body_bp_b", "d_exon_union_bp_a",
            "d_exon_union_bp_b", "d_e1_identity", "d_e1_identity_gapexcl", "d_e1_cov_longer", "d_e1_cov_gene",
            "d_e1_cov_denominator", "d_e1_exonic_bp_covered_longer", "d_e1_gate_exonic", "d_e1_edge", "d_e1_weight",
            "d_shared_exon_bp", "d_shared_exon_frac", "d_shared_exon_denominator", "d_shared_exon_frac_allrec",
            "d_w98_gapexcl", "d_w98_gapincl", "d_both_spliced",
            "d_c2_genebody", "d_c2_gb_best_frac", "d_c2_gb_chain_identity", "d_c2_gb_direction", "d_c2_exon",
            "d_c2_exon_best_frac", "d_c2_exon_direction", "d_c2_approx", "d_c2nostrand_approx", "d_c2exontgt_approx",
            "d_c2loose_genebody", "d_c2loose_gb_best_frac", "d_c2loose_gb_chain_identity", "d_c2loose_exon",
            "d_c2loose_exon_best_frac", "d_c2loose_approx",
            "d_c2x_genebody", "d_c2x_gb_best_frac", "d_c2x_approx", "d_c2x_gb_chain_raw", "d_c2xloose_approx",
            "d_clause2_evaluation",
            "s2_edge", "s2_max_identity", "s2_n_shared_exons", "s2_n_projected_exon_pairs", "s2_same_locus",
            "ctree_family", "ctree_smallest_common_cluster", "ctree_support", "ctree_same_top_cluster"]

    def row_for(k):
        x, y = k
        r = {"gene_a": x, "gene_b": y, "name_a": genes[x]["name"], "name_b": genes[y]["name"], "chrom_a": genes[x]["chrom"],
             "chrom_b": genes[y]["chrom"], "same_locus": same_locus(x, y)}
        p = prot.get(k)
        r["p_evidence"] = p is not None
        r["p_searched_a"] = (g2p[x] in searched) if x in g2p else "no_protein"
        r["p_searched_b"] = (g2p[y] in searched) if y in g2p else "no_protein"
        if p:
            r.update(p)
            r["p_evalue"] = "<=1e-5 (search cutoff; per-HSP value not saved)"
            r["p_aa50"] = p["p_qualifies_6ko"] and p["p_aa_identity"] >= 0.50
        else:
            r.update({"p_qualifies_6ko": False, "p_aa50": False, "p_evalue": "NA", "p_qualifies_union": False})
        d = gdna.get(k)
        r["d_evidence"] = d is not None
        if d:
            r.update(d)
            r["d_e1_cov_denominator"] = (f"exon-union bp of gene {d['d_e1_cov_gene']} (longer exon union): "
                                         f"{d['d_exon_union_bp_a'] if d['d_e1_cov_gene'] == 'A' else d['d_exon_union_bp_b']}"
                                         if d["d_e1_cov_gene"] else "NA")
            r["d_shared_exon_denominator"] = min(d["d_exon_union_bp_a"], d["d_exon_union_bp_b"])
            r["d_clause2_evaluation"] = ("approximate (gene-body chain on gene-body PAF; exon proxy; both require >= 1 exonic "
                                         "base of v on the target side and the shipped strand check for spliced pairs)")
        else:
            cats = {gene_key.get(x) and key_cat[gene_key[x]], gene_key.get(y) and key_cat[gene_key[y]]}
            why = ("no PAF record for the pair" if len(cats) == 1 and None not in cats
                   else "not in one E1 catalog (cross-catalog or outside both)")
            r.update({"d_catalog": "NA", "d_e1_edge": False, "d_c2_approx": False, "d_c2_genebody": False, "d_c2_exon": False,
                      "d_c2nostrand_approx": False, "d_c2exontgt_approx": False, "d_c2loose_genebody": False,
                      "d_c2loose_exon": False, "d_c2loose_approx": False,
                      "d_c2x_genebody": False, "d_c2x_approx": False, "d_c2x_gb_chain_raw": False, "d_c2xloose_approx": False,
                      "d_shared_exon_frac": "NA", "d_shared_exon_frac_allrec": "NA", "d_w98_gapexcl": "NA",
                      "d_w98_gapincl": "NA", "d_clause2_evaluation": f"not satisfiable: {why}"})
        s = s2.get(k)
        r["s2_edge"] = s is not None
        if s:
            r.update({"s2_max_identity": s["max_identity"], "s2_n_shared_exons": s["n_shared_exons"],
                      "s2_n_projected_exon_pairs": s["n_projected_exon_pairs"], "s2_same_locus": s["same_locus"]})
        else:
            r["s2_max_identity"] = "NA"
        fam, cl, sup, same_top = ctree_pair(x, y)
        r.update({"ctree_family": fam, "ctree_smallest_common_cluster": cl, "ctree_support": sup,
                  "ctree_same_top_cluster": same_top})
        return r

    keys = set()
    for k in prot:
        if k[0] in V and k[1] in V:
            keys.add(k)
    for k in gdna:
        if k[0] in V and k[1] in V:
            keys.add(k)
    for k in s2:
        if k[0] in V and k[1] in V:
            keys.add(k)
    rows = [row_for(k) for k in sorted(keys)]
    write(f"{OUT}/edges.tsv", rows, cols)
    say(f"[write] edges.tsv rows {len(rows)} (protein HSP pairs {sum(1 for r in rows if r['p_evidence'])}, PAF pairs "
        f"{sum(1 for r in rows if r['d_evidence'])}, S2 pairs {sum(1 for r in rows if r['s2_edge'])}); {time.time() - T0:.0f}s")
    # every clause-2 DNA edge of both catalogs (for DNA-level components outside V, e.g. chaining checks)
    c2rows = [{"gene_a": k[0], "gene_b": k[1], "name_a": genes[k[0]]["name"], "name_b": genes[k[1]]["name"],
               "same_locus": same_locus(*k), "d_c2_approx": v["d_c2_approx"], "d_c2loose_approx": v["d_c2loose_approx"],
               "d_w98_gapexcl": v["d_w98_gapexcl"], "d_e1_identity_gapexcl": v["d_e1_identity_gapexcl"],
               "d_e1_identity": v["d_e1_identity"], "d_shared_exon_frac": v["d_shared_exon_frac"], "d_e1_edge": v["d_e1_edge"]}
              for k, v in sorted(gdna.items()) if v["d_c2_approx"] or v["d_c2loose_approx"]]
    write(f"{OUT}/edges_all_c2.tsv", c2rows)

    # ------------------------------------------------------------------------------------------------ nodes
    all_hg = hgnc_all(genes)
    db = soto_load()
    lit = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
    nrows = []
    for g in sorted(V):
        r = genes[g]
        ex = exons_tsv.get(g) or [(int(r["start0"]), int(r["end"]))]
        m = soto_map_gene(db, r["name"], r["chrom"], r["strand"], ex)
        flag = soto_flag(m)
        u = U.get(g, {})
        lr = lit.get(g, {})
        nrows.append({"gene_id": g, "name": r["name"], "biotype": r["biotype"], "chrom": r["chrom"], "start0": r["start0"],
                      "end": r["end"], "strand": r["strand"], "readthrough": "readthrough" in r["description"],
                      "in_U": g in U, "is_member": u.get("is_member", "no"), "member_family": u.get("member_family", ""),
                      "family_side": u.get("family_side", ""), "l0_hops_from_U": V[g], "in_primary_L0_closure": g in V0, "catalog": key_cat[gene_key[g]] if g in gene_key else "none",
                      "has_protein": g in g2p, "protein_searched": (g2p[g] in searched) if g in g2p else "no_protein",
                      "hgnc_gene_group_id": all_hg.get(g, ""), "soto_families": m["soto_families"], "soto_flag": flag,
                      "lit_level1": lr.get("level1", ""), "lit_level2": lr.get("level2", ""),
                      "lit_named_npipb": lr.get("npipb_named_subfamily", ""), "lit_in_truth": lr.get("in_literature_truth", ""),
                      "P_group": u.get("P_group", ""), "in_P_universe": u.get("in_P_universe", ""),
                      "D_group": u.get("D_group", ""), "in_D_universe": u.get("in_D_universe", ""),
                      "C_L1": u.get("C_L1", ""), "C_fine": u.get("C_fine", ""), "in_C_universe": u.get("in_C_universe", ""),
                      "Ctree_top": CT_LAB["Ctree_top"].get(g, ""), "Ctree_min": CT_LAB["Ctree_min"].get(g, ""),
                      "description": r["description"]})
    write(f"{OUT}/nodes.tsv", nrows)
    say(f"[write] nodes.tsv rows {len(nrows)}; {time.time() - T0:.0f}s")
    say.dump(f"{OUT}/edges_build.out")


def cmd_lattice_expr(args):
    """was bench/layer_order/lattice_expr.py

    Nested edge-test lattice, step 2: testis read counts for every node of V (lattice/nodes.tsv).

    Rule = expr-recount (was lo_expr_recount.py; itself = heavy/scripts/expr_counts.py; one function,
    lattice_common.count_reads): human_testis.t2t.bam, primary reads only (samtools -F 2308), read blocks split at N and
    D, a read counts for gene G if >= 1 block overlaps >= 1 bp of an exon of G (strand ignored); 'unique' = the read's
    blocks hit exons of exactly one RefSeq record genome-wide. Exon-less records count on their gene body. Only the gene
    set differs (all of V instead of U). Check: the 68 U genes must reproduce integrate_slim/expr_recount.tsv (any and
    unique).

    -> lattice/expr_counts.tsv, lattice/expr_counts.windows.bed, lattice/expr_counts.out
    """
    INT, OUT = LC.INT, LC.OUT
    T0 = time.time()
    say = Log()
    want = {r["gene_id"] for r in tsv(f"{OUT}/nodes.tsv")}
    genes, exons, index = gff_exon_index(strip_prefix=False)
    say(f"[gff] gene/pseudogene ids {len(genes)}; V genes {len(want)} (missing from GFF {len(want - set(genes))}); "
        f"{time.time() - T0:.0f}s")
    spans = defaultdict(list)
    for g in want:
        c, s, e = genes[g]
        spans[c].append((s, e))
    bed = f"{OUT}/expr_counts.windows.bed"
    nwin, nbp = write_windows_bed(bed, spans)
    n_any, n_uni, _, n_rec, n_hit = count_reads(bed, index, want)
    say(f"[bam] windows {nwin} ({nbp} bp); primary records {n_rec}; on >= 1 exon {n_hit}; {time.time() - T0:.0f}s")
    prev = {"gene-" + r["gene_id"]: r for r in tsv(f"{INT}/expr_recount.tsv")}
    diff = [g for g in prev if g in want and (n_any[g], n_uni[g]) != (int(prev[g]["n_reads_any"]), int(prev[g]["n_reads_unique"]))]
    say(f"[check] genes also in integrate_slim/expr_recount.tsv: {sum(1 for g in prev if g in want)}; differing any/unique: "
        f"{len(diff)} {diff[:10]}")
    with open(f"{OUT}/expr_counts.tsv", "w") as out:
        out.write("gene_id\tn_reads_any\tn_reads_unique\n")
        for g in sorted(want):
            out.write(f"{g}\t{n_any[g]}\t{n_uni[g]}\n")
    for t in (1, 3):
        say(f"[summary] V genes with any-overlap reads >= {t}: {sum(1 for g in want if n_any[g] >= t)}")
    say.dump(f"{OUT}/expr_counts.out")


def cmd_lattice_levels(args):
    """was bench/layer_order/lattice_levels.py

    Nested edge-test lattice, step 3: levels G_0..G_3 (primary tests in lattice_common.tests), their 3-truss (triangle)
    variant, T1/T2 sanity checks, chaining evidence, operationalisation variants.

    Inputs: lattice/nodes.tsv, lattice/edges.tsv (lattice_edges.py), lattice/expr_counts.tsv (lattice_expr.py).
    Outputs (lattice/):
      groups.tsv          per gene: component label (gene_id of the component representative) per level x variant
      levels.tsv          per level x variant: edges, components, largest, member-holding groups
      member_groups.tsv   every member-holding group per level x variant: size, members, pulled-in non-members
      sanity.tsv          T1 (nesting across levels) and T2 (expression views) violation counts (must be 0), with the
                          non-vacuity counts: coarse blocks with >= 2 genes, how many the finer partition splits, member-holding
      expr_views.tsv      member-holding expressed components per expression set x variant x level, including both 3-truss
                          views comp(truss(G_k[X])) and comp(truss(G_k)[X])
      chaining.tsv        genes of interest (PKD1 / readthrough / DHX40 / RNFT1 / TBC-domain neighbours): group per level and
                          the shortest connecting path to the family anchor with its edge attributes
      triangle_drops.tsv  member-holding groups of G_k that lose genes in the 3-truss variant (2-copy groups included)
      levels.out          log
    """
    OUT = LC.OUT
    T0 = time.time()
    say = Log()

    nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
    NAME = {g: r["name"] for g, r in nodes.items()}
    MEM = {g for g, r in nodes.items() if r["is_member"] == "yes"}
    FAM = {g: r["member_family"] for g, r in nodes.items() if r["is_member"] == "yes"}
    _ec = tsv(f"{OUT}/expr_counts.tsv")
    reads = {r["gene_id"]: int(r["n_reads_any"]) for r in _ec}
    reads_u = {r["gene_id"]: int(r["n_reads_unique"]) for r in _ec}
    V = set(nodes)

    # ---------------------------------------------------------------------------------------------- edges (compact)
    NEED = ["gene_a", "gene_b", "same_locus", "p_qualifies_6ko", "p_aa_identity", "d_c2_approx", "d_e1_edge", "d_e1_identity",
            "d_e1_cov_longer", "d_shared_exon_frac", "d_e1_identity_gapexcl", "s2_max_identity", "d_c2_genebody", "d_c2_exon",
            "p_cov_longer", "d_c2_gb_chain_identity", "d_c2_gb_best_frac", "d_c2_exon_best_frac", "d_c2x_approx",
            "d_c2nostrand_approx", "d_c2exontgt_approx", "d_c2loose_approx", "d_w98_gapexcl", "d_w98_gapincl"]
    ROWS = []
    with open(f"{OUT}/edges.tsv") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ix = {c: hdr.index(c) for c in NEED}
        for f in rd:
            ROWS.append({c: f[i] for c, i in ix.items()})
    say(f"[load] nodes {len(V)}; edge rows {len(ROWS)}; {time.time() - T0:.0f}s")
    EDGE = {(r["gene_a"], r["gene_b"]): r for r in ROWS}

    VARIANTS = {  # primary: L1 = clause-2 approx with v-exon overlap + strand check; L3 = single-record w_98, gap-excluded
        "primary": dict(),
        "with_same_locus": dict(with_same_locus=True),
        "L0=(P_and_aa>=0.50)_or_D": dict(p_aa_min=0.50),
        "L1=c2_no_strand_check": dict(l1="c2_nostrand"),
        "L1=c2_vexon_on_exon_proxy_only": dict(l1="c2_exontgt"),
        "L1=c2_no_vexon_no_strand": dict(l1="c2_loose"),
        "L1=c2x_extrapolated": dict(l1="c2x"),
        "L1=E1_as_built(D graph)": dict(l1="e1"),
        "L1=E1_at_0.80/0.50": dict(l1="e1c2"),
        "L3=w98_gap-inclusive": dict(id_which="w98_gapincl"),
        "L3=pooled_gap-excluded": dict(id_which="pooled_gapexcl"),
        "L3=pooled_gap-inclusive": dict(id_which="pooled_gapincl"),
        "L3=S2_SD98_mapback": dict(id_which="s2"),
        "17:03_tests_exact": dict(l1="c2_loose", id_which="pooled_gapexcl"),  # the 17:12 report's tests, unrounded thresholds
    }

    def level_edges(**kw):
        E = {k: [] for k in LEVELS}
        for r in ROWS:
            t = tests(r, **kw)
            for k, ok in zip(LEVELS, t):
                if ok:
                    E[k].append((r["gene_a"], r["gene_b"]))
        return E

    def summarize(tag, lab, E, nodeset, rows_levels, rows_groups):
        G = groups(lab)
        nonsingle = [S for S in G.values() if len(S) >= 2]
        mem_groups = [S for S in G.values() if S & MEM]
        rows_levels.append({"variant": tag[0], "level": tag[1], "nodes": len(nodeset), "edges": len(E),
                            "components_ge2": len(nonsingle), "singletons": sum(1 for S in G.values() if len(S) == 1),
                            "largest": max(len(S) for S in G.values()),
                            "member_groups": len(mem_groups),
                            "member_groups_sizes": ";".join(f"{len(S)}(NPIP {sum(1 for g in S & MEM if FAM[g] == 'NPIP')},TBC1D3 "
                                                            f"{sum(1 for g in S & MEM if FAM[g] == 'TBC1D3')})"
                                                            for S in sorted(mem_groups, key=lambda s: (-len(s), sorted(s))))})
        for S in sorted(mem_groups, key=lambda s: (-len(s), sorted(s))):
            mems = sorted(NAME[g] for g in S & MEM)
            non = sorted(S - MEM, key=lambda g: (int(nodes[g]["l0_hops_from_U"]), NAME[g]))
            bt = collections.Counter(nodes[g]["biotype"] for g in S)
            ch = collections.Counter(nodes[g]["chrom"] for g in S)
            rows_groups.append({"variant": tag[0], "level": tag[1], "group_rep": NAME[min(S)], "size": len(S),
                                "n_NPIP_members": sum(1 for g in S & MEM if FAM[g] == "NPIP"),
                                "n_TBC1D3_members": sum(1 for g in S & MEM if FAM[g] == "TBC1D3"),
                                "members": ",".join(mems), "n_nonmembers": len(non),
                                "nonmembers_first80": ",".join(NAME[g] for g in non[:80]),
                                "n_outside_E1_catalogs": sum(1 for g in S if nodes[g]["catalog"] == "none"),
                                "n_protein_never_searched": sum(1 for g in S if nodes[g]["protein_searched"] == "no"),
                                "biotypes": ";".join(f"{k}:{v}" for k, v in bt.most_common()),
                                "chroms": ";".join(f"{k}:{v}" for k, v in ch.most_common())})

    rows_levels, rows_groups, gl_rows = [], [], {g: {"gene_id": g, "name": NAME[g], "is_member": nodes[g]["is_member"],
                                                       "member_family": nodes[g]["member_family"],
                                                       "n_reads_any": reads.get(g, 0)} for g in V}
    LAB = {}
    EDGES = {}
    for vname, kw in VARIANTS.items():
        E = level_edges(**kw)
        EDGES[vname] = E
        for k in LEVELS:
            lab = components(V, E[k])
            LAB[(vname, k)] = lab
            summarize((vname, k), lab, E[k], V, rows_levels, rows_groups)
            for g in V:
                gl_rows[g][f"{vname}|{k}"] = lab[g]
        say(f"[levels] {vname}: edges " + ", ".join(f"{k} {len(E[k])}" for k in LEVELS) + f"; {time.time() - T0:.0f}s")

    # ---- triangle (3-truss) variant of the primary levels (and of the 17:12 report's tests, unrounded)
    TRI = {}
    for base, tag in (("primary", "triangle"), ("17:03_tests_exact", "17:03_tests_exact_triangle")):
        for k in LEVELS:
            kept, dropped, depth = truss3(EDGES[base][k])
            TRI[(tag, k)] = kept
            lab = components(V, kept)
            LAB[(tag, k)] = lab
            summarize((tag if tag != "triangle" else "triangle(3-truss)", k), lab, kept, V, rows_levels, rows_groups)
            for g in V:
                gl_rows[g][f"{tag}|{k}"] = lab[g]
            say(f"[{tag}] {k}: edges {len(EDGES[base][k])} -> {len(kept)} (dropped {dropped}, peel depth {depth}); "
                f"{time.time() - T0:.0f}s")

    # ---- node variant: readthrough records removed (clause 1: readthrough spans must not be nodes)
    RT = {g for g in V if nodes[g]["readthrough"] == "yes"}
    V_nrt = V - RT
    for k in LEVELS:
        E = [e for e in EDGES["primary"][k] if e[0] in V_nrt and e[1] in V_nrt]
        lab = components(V_nrt, E)
        LAB[("no_readthrough_nodes", k)] = lab
        summarize(("no_readthrough_nodes", k), lab, E, V_nrt, rows_levels, rows_groups)
        for g in V:
            gl_rows[g][f"no_readthrough_nodes|{k}"] = lab[g] if g in lab else "removed"
    say(f"[no-readthrough] readthrough records in V removed: {len(RT)} ({sorted(NAME[g] for g in RT & MEM)} are members)")

    write(f"{OUT}/levels.tsv", rows_levels)
    write(f"{OUT}/member_groups.tsv", rows_groups)
    gcols = ["gene_id", "name", "is_member", "member_family", "n_reads_any"] + [c for c in next(iter(gl_rows.values())) if "|" in c]
    write(f"{OUT}/groups.tsv", [gl_rows[g] for g in sorted(V, key=lambda g: NAME[g])], gcols)

    # ---------------------------------------------------------------------------------------------- sanity T1 / T2
    srows = []

    def add(check, variant, detail, fine, coarse, nodes_, viol):
        """one sanity row; groups_checked = all fine groups (singletons included, the 17:12 count); non-vacuity columns count
        coarse blocks with >= 2 genes that the fine partition actually splits."""
        nb, ns, nm = split_counts(fine, coarse, nodes_, MEM)
        srows.append({"check": check, "variant": variant, "detail": detail,
                      "groups_checked": len(groups({n: fine[n] for n in nodes_})), "violations": len(viol),
                      "coarse_blocks_ge2": nb, "coarse_blocks_split": ns, "member_holding_blocks_split": nm,
                      "example": ";".join(",".join(sorted(NAME[g] for g in S))[:200] for S in viol[:2])})

    T1_VARIANTS = list(VARIANTS) + ["triangle", "17:03_tests_exact_triangle", "no_readthrough_nodes"]
    for vname in T1_VARIANTS:
        for i in range(len(LEVELS)):
            for j in range(i + 1, len(LEVELS)):
                fine, coarse = LAB[(vname, LEVELS[j])], LAB[(vname, LEVELS[i])]
                viol = refines(fine, coarse)
                add("T1 nesting: components of G_j refine G_i", vname, f"{LEVELS[j]} in {LEVELS[i]}", fine, coarse, set(fine), viol)
    # triangle inside plain components (truss(E_k) subset of E_k)
    for tag, base in (("triangle", "primary"), ("17:03_tests_exact_triangle", "17:03_tests_exact")):
        for k in LEVELS:
            viol = refines(LAB[(tag, k)], LAB[(base, k)])
            add("triangle components refine plain components", f"{tag} vs {base}", k, LAB[(tag, k)], LAB[(base, k)],
                set(LAB[(tag, k)]), viol)

    XSETS = (("any-overlap reads>=3", {g for g in V if reads.get(g, 0) >= 3}),
             ("any-overlap reads>=1", {g for g in V if reads.get(g, 0) >= 1}),
             ("unique reads>=3", {g for g in V if reads_u.get(g, 0) >= 3}))
    erows = []

    def expr_row(xname, vname, view, k, lab):
        G = groups(lab)
        mg = sorted([S for S in G.values() if S & MEM and len(S) >= 2], key=lambda S: (-len(S), sorted(S)))
        erows.append({"expression_set": xname, "variant": vname, "view": view, "level": k,
                      "member_components": ";".join(f"{len(S)}({len(S & MEM)})" for S in mg),
                      "members": " | ".join(",".join(sorted(NAME[g] for g in S & MEM)) for S in mg),
                      "nonmembers": " | ".join(",".join(sorted(NAME[g] for g in S - MEM)) for S in mg),
                      "member_singletons": ",".join(sorted(NAME[g] for S in G.values() if len(S) == 1 for g in S & MEM))})
        say(f"[EXPR {xname}] {vname} {view} {k}: member-holding expressed components (>= 2 genes): "
            + " | ".join(f"{len(S)}: members {sorted(NAME[g] for g in S & MEM)} + {len(S - MEM)} non-members" for S in mg))

    for xname, X in XSETS:
        say(f"[T2] expression set X: {xname}: {len(X)} of {len(V)} nodes (members {len(X & MEM)})")
        for vname in ("primary", "with_same_locus", "L1=E1_as_built(D graph)", "17:03_tests_exact"):
            labX = {}
            for k in LEVELS:
                EX = [e for e in EDGES[vname][k] if e[0] in X and e[1] in X]
                labX[k] = components(X, EX)
                viol = refines(labX[k], LAB[(vname, k)], X)
                add(f"T2b: G_k[X] components refine G_k restricted to X ({xname})", vname, k, labX[k], LAB[(vname, k)], X, viol)
                if vname in ("primary", "17:03_tests_exact"):
                    expr_row(xname, vname, "comp(G_k[X])", k, labX[k])
            for i in range(len(LEVELS) - 1):
                viol = refines(labX[LEVELS[i + 1]], labX[LEVELS[i]])
                add(f"T2a: G_(k+1)[X] components refine G_k[X] ({xname})", vname, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                    labX[LEVELS[i + 1]], labX[LEVELS[i]], X, viol)
        # triangle: two expression views, both nest (truss(G_k[X]) subset of truss(G_k)[X] subset of G_k[X])
        for base, tag in (("primary", "triangle"), ("17:03_tests_exact", "17:03_tests_exact_triangle")):
            labTX, labTX2 = {}, {}
            for k in LEVELS:
                EX = [e for e in EDGES[base][k] if e[0] in X and e[1] in X]
                keptX, _, _ = truss3(EX)
                labTX[k] = components(X, keptX)                                               # comp(truss(G_k[X]))
                labTX2[k] = components(X, [e for e in TRI[(tag, k)] if e[0] in X and e[1] in X])  # comp(truss(G_k)[X])
                viol = refines(labTX[k], LAB[(tag, k)], X)
                add(f"T2b: truss(G_k[X]) components refine truss(G_k) restricted to X ({xname})", tag, k, labTX[k],
                    LAB[(tag, k)], X, viol)
                viol = refines(labTX[k], labTX2[k], X)
                add(f"truss views: comp(truss(G_k[X])) refines comp(truss(G_k)[X]) ({xname})", tag, k, labTX[k], labTX2[k], X,
                    viol)
                viol = refines(labTX2[k], LAB[(tag, k)], X)
                add(f"T2b: comp(truss(G_k)[X]) refines truss(G_k) restricted to X ({xname})", tag, k, labTX2[k],
                    LAB[(tag, k)], X, viol)
                expr_row(xname, base, "comp(truss(G_k[X]))", k, labTX[k])
                expr_row(xname, base, "comp(truss(G_k)[X])", k, labTX2[k])
            for i in range(len(LEVELS) - 1):
                viol = refines(labTX[LEVELS[i + 1]], labTX[LEVELS[i]])
                add(f"T2a: truss(G_(k+1)[X]) refine truss(G_k[X]) ({xname})", tag, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                    labTX[LEVELS[i + 1]], labTX[LEVELS[i]], X, viol)
                viol = refines(labTX2[LEVELS[i + 1]], labTX2[LEVELS[i]])
                add(f"T2a: truss(G_(k+1))[X] refine truss(G_k)[X] ({xname})", tag, f"{LEVELS[i + 1]} in {LEVELS[i]}",
                    labTX2[LEVELS[i + 1]], labTX2[LEVELS[i]], X, viol)
        say(f"[T2] {xname} done; {time.time() - T0:.0f}s")
    write(f"{OUT}/sanity.tsv", srows)
    write(f"{OUT}/expr_views.tsv", erows)
    say(f"[sanity] checks {len(srows)}; total violations {sum(r['violations'] for r in srows)}")

    # ---------------------------------------------------------------------------------------------- chaining
    INTEREST = ["PKD1", "PKD1P1", "PKD1P2", "PKD1P3", "PKD1P6", "PKD1P3-NPIPA1", "LOC131696449", "PKD1P4-NPIPA8",
                "PKD1P5-LOC105376752", "PKD1P6-NPIPP1", "PDXDC2P-NPIPB14P", "NPIPB1P", "NPIPB14P", "LOC100505915",
                "DHX40", "DHX40P1", "RNFT1", "RNFT1-DT", "RNFT1P3", "TBC1D3P1-DHX40P1", "TBC1D3P1", "USP6", "USP6NL", "TBC1D26",
                "TBC1D29P", "TBC1D28", "LOC100420408", "TBC1D3P5", "TBC1D3P7", "LOC124905656", "TBC1D3P6", "LOC100420289"]
    by_name = collections.defaultdict(list)
    for g in V:
        by_name[NAME[g]].append(g)
    ANCHOR = {"NPIP": by_name["NPIPB2"][0], "TBC1D3": by_name["TBC1D3"][0]}

    def side_of(name):
        return "NPIP" if any(x in name for x in ("PKD1", "NPIP", "PDXDC2P", "LOC131696449", "LOC100505915")) else "TBC1D3"

    ADJ = {}
    for _k in LEVELS:
        ADJ[_k] = collections.defaultdict(set)
        for _a, _b in EDGES["primary"][_k]:
            ADJ[_k][_a].add(_b)
            ADJ[_k][_b].add(_a)

    def edge_desc(a, b):
        r = EDGE.get((a, b)) or EDGE.get((b, a))
        return (f"{NAME[a]}-{NAME[b]}[p {r['p_qualifies_6ko']} aa {fmt3(r['p_aa_identity'])}; c2 {r['d_c2_approx']} "
                f"(gb {r['d_c2_genebody']}, ex {r['d_c2_exon']}); sef {fmt3(r['d_shared_exon_frac'])}; w98 {fmt3(r['d_w98_gapexcl'])}; "
                f"pooled id {fmt3(r['d_e1_identity_gapexcl'])}; same_locus {r['same_locus']}]")

    def fmt3(s):
        v = fnum(s)
        return "NA" if v is None else f"{v:.3f}"

    crow = []
    for nm in INTEREST:
        for g in by_name.get(nm, []):
            fam = side_of(nm)
            anc = ANCHOR[fam]
            row = {"gene": nm, "is_member": nodes[g]["is_member"], "catalog": nodes[g]["catalog"], "anchor": NAME[anc],
                   "readthrough": nodes[g]["readthrough"], "n_reads_any": reads.get(g, 0)}
            for vname in ("primary", "triangle", "no_readthrough_nodes"):
                for k in LEVELS:
                    lab = LAB[(vname, k)]
                    if g not in lab:
                        row[f"{vname}|{k}"] = "removed"
                        continue
                    row[f"{vname}|{k}"] = "with anchor" if lab[g] == lab[anc] else f"apart ({len(groups(lab)[lab[g]])})"
            for k in LEVELS:
                p = shortest_path(ADJ[k], anc, g) if LAB[("primary", k)][g] == LAB[("primary", k)][anc] else None
                row[f"path|{k}"] = " > ".join(edge_desc(p[i], p[i + 1]) for i in range(len(p) - 1)) if p and len(p) <= 8 else (
                    f"path of {len(p) - 1} edges" if p else "")
            crow.append(row)
    write(f"{OUT}/chaining.tsv", crow)

    # ---------------------------------------------------------------------------------------------- triangle drops
    trows = []
    for k in LEVELS:
        G = groups(LAB[("primary", k)])
        for S in G.values():
            if not (S & MEM) or len(S) < 2:
                continue
            labT = LAB[("triangle", k)]
            parts = groups({g: labT[g] for g in S})
            big = max(parts.values(), key=lambda P: (len(P & MEM), len(P), sorted(P)))  # the part holding most members
            lost = S - big
            trows.append({"level": k, "group_rep": NAME[min(S)], "size": len(S), "members": ",".join(sorted(NAME[g] for g in S & MEM)),
                          "n_members": len(S & MEM), "triangle_parts": len(parts),
                          "triangle_singletons": sum(1 for P in parts.values() if len(P) == 1),
                          "main_part_size": len(big), "main_part_n_members": len(big & MEM),
                          "largest_part_size": max(len(P) for P in parts.values()),
                          "genes_outside_main_part": len(lost),
                          "members_outside_main_part": ",".join(sorted(NAME[g] for g in lost & MEM)),
                          "nonmembers_outside_main_part_first40": ",".join(sorted(NAME[g] for g in lost - MEM)[:40]),
                          "is_2copy_group": len(S) == 2})
    write(f"{OUT}/triangle_drops.tsv", trows)
    say.dump(f"{OUT}/levels.out")
    say(f"[done] {time.time() - T0:.0f}s")


def cmd_lattice_truth(args):
    """was bench/layer_order/lattice_truth.py

    Nested edge-test lattice, step 4: agreement of every level with the truths, side by side with the prior study's
    report-only groupings (P MCL, D = E1 MCL, literature C, clause-5 C_tree) on the SAME genes.

    Metrics (same definitions as bench/layer_order/lo_analysis.score, re-implemented with group counts for speed and checked
    against integrate_slim/truth_agreement.tsv): genes = those with a truth label (and in the prediction's universe); pairwise
    precision / recall over same-group pairs; bipartite F with ONE-TO-ONE JACCARD matching (lo_analysis.bip_jaccard: Hungarian
    assignment maximising summed Jaccard; R = matched genes / genes, P = matched genes / genes of matched predicted groups);
    F is NA when the truth or the prediction has 0 same-group pairs.

    Truths (per gene, lattice/nodes.tsv): HGNC gene_group_id (superfamily-level for TBC1D3); Soto family (flag ok: matched,
    not weak, not ambiguous — exon-overlap mapping of light/scripts/soto_map.py); literature L1 NPIPA|NPIPB (NPIP only);
    literature L2 paralog groups (NPIP A/B groups; TBC1D3 Guitart Fig 6B/6C groups M, CDKL — phylogenetic, read post hoc;
    not positional. CORRECTED 2026-09-16, was: "TBC1D3 positional AE/CDKL as in the prior study" — the name-mapped
    AE/CDKL truth was wrong for 7/9 copies; see bench/TBC1D3_GUITART_TRUTH_CORRECTION.md. Comment only, no behaviour
    change: this module reads level2 from whatever truth TSV it is pointed at.).

    Gene sets:
      U     the 68-gene universe of the prior study, per family side (NPIP / TBC1D3 / pooled).
            'U all'   : every U gene with a truth label (a lattice node without edges is its own group)
            'U ∩ X'   : the report-only layer X's universe (P: coding; D: E1 catalog nodes; C: tree leaves)
      V     the closure V (8,070 genes), pooled; lattice levels only (recall conditioned on V).
      group the member-holding group of each family anchor (NPIPB2, TBC1D3) at each level: pair precision of ALL its genes
            with a truth label (members and pulled-in non-members), with the truth-family composition. U-scope scores only see
            members; this is the in-group precision that must be read next to them.
    Circularity: Soto families are SD98 (>= 98% identity) duplications with a shared-exon map-back, the conventions L2
    (shared-exon >= 0.30) and L3 (identity >= 0.98) test, so Soto agreement at L2/L3 is partly by construction. Clause 5
    (C_tree) was developed on NPIP (§6jp-§6js) and is not an independent comparator on these families.
    Outputs: lattice/truth.tsv, lattice/truth_ingroup.tsv, lattice/truth.out
    """
    INT, OUT = LC.INT, LC.OUT
    say = Log()

    nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
    grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
    V = set(nodes)
    U = {g for g, r in nodes.items() if r["in_U"] == "yes"}
    SIDE = {g: nodes[g]["family_side"] for g in U}

    LATTICE = {  # label -> groups.tsv column
        "L0": "primary|L0", "L1": "primary|L1", "L2": "primary|L2", "L3": "primary|L3",
        "L0Δ": "triangle|L0", "L1Δ": "triangle|L1", "L2Δ": "triangle|L2", "L3Δ": "triangle|L3",
        "L1[E1@.80/.50]": "L1=E1_at_0.80/0.50|L1", "L2[E1@.80/.50]": "L1=E1_at_0.80/0.50|L2", "L3[E1@.80/.50]": "L1=E1_at_0.80/0.50|L3",
        "L1[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L1", "L2[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L2",
        "L3[no v-exon/strand]": "L1=c2_no_vexon_no_strand|L3",
        "L3[w98 gap-incl]": "L3=w98_gap-inclusive|L3", "L3[pooled gap-excl]": "L3=pooled_gap-excluded|L3",
        "L3[pooled gap-incl]": "L3=pooled_gap-inclusive|L3", "L3[id S2]": "L3=S2_SD98_mapback|L3",
        "17:03 L0": "17:03_tests_exact|L0", "17:03 L1": "17:03_tests_exact|L1", "17:03 L2": "17:03_tests_exact|L2",
        "17:03 L3": "17:03_tests_exact|L3", "17:03 L3Δ": "17:03_tests_exact_triangle|L3",
    }
    LAT = {k: {g: grp[g][c] for g in V} for k, c in LATTICE.items()}
    REPORT = {
        "P (§6ko MCL)": ({g: nodes[g]["P_group"] for g in U if nodes[g]["in_P_universe"] == "yes"}),
        "D (E1 MCL)": ({g: nodes[g]["D_group"] for g in U if nodes[g]["in_D_universe"] == "yes"}),
        "C_L1 (lit, circular)": ({g: nodes[g]["C_L1"] for g in U if nodes[g]["in_C_universe"] == "yes"}),
        "C_fine (lit, circular)": ({g: nodes[g]["C_fine"] for g in U if nodes[g]["in_C_universe"] == "yes"}),
        "Ctree_top (clause 5)": ({g: nodes[g]["Ctree_top"] for g in U if nodes[g]["Ctree_top"]}),
        "Ctree_min (clause 5)": ({g: nodes[g]["Ctree_min"] for g in U if nodes[g]["Ctree_min"]}),
    }
    TRUTH = {
        "HGNC gene group": {g: r["hgnc_gene_group_id"] for g, r in nodes.items() if r["hgnc_gene_group_id"]},
        "Soto family (flag ok)": {g: r["soto_families"] for g, r in nodes.items() if r["soto_flag"] == "ok" and r["soto_families"]},
        "literature L1 (NPIPA|NPIPB)": {g: r["lit_level1"] for g, r in nodes.items() if r["lit_level1"] and r["member_family"] == "NPIP"},
        "literature L2": {g: r["lit_level2"] for g, r in nodes.items() if r["lit_level2"]},
    }
    for t, d in TRUTH.items():
        say(f"[truth] {t}: genes labelled in V {len(d)}, in U {len(set(d) & U)}")

    rows = []

    def fam_genes(fam):
        return set(U) if fam == "pooled" else {g for g in U if SIDE[g] == fam}

    # ---- check the scorer against the prior study (integrate_slim/truth_agreement.tsv, P as built / D as built, Soto)
    prior = {(r["layer"], r["variant"], r["truth"], r["family"]): r for r in tsv(f"{INT}/truth_agreement.tsv")}
    for lay, var, name in (("P", "as built", "P (§6ko MCL)"), ("D", "as built = after (D unchanged)", "D (E1 MCL)")):
        for tname, ptname in (("Soto family (flag ok)", "Soto family (flag ok)"),
                              ("HGNC gene group", "HGNC gene_group_id (superfamily-level for TBC1D3)")):
            for fam in ("NPIP", "TBC1D3", "pooled"):
                pr = prior[(lay, var, ptname, fam)]
                lab = REPORT[name]
                s = score_counts(lab, TRUTH[tname], fam_genes(fam) & set(lab))
                fj = s.get("bip_F_jaccard")
                fj = fj if isinstance(fj, str) else ("NA" if fj is None else f"{fj:.3f}")
                mine = (s["n_genes"], s.get("pred_pairs"), s.get("tp_pairs"), fj)
                ok = (str(s["n_genes"]) == pr["n_genes"] and str(s.get("pred_pairs", "")) == pr.get("pred_pairs", "")
                      and str(s.get("tp_pairs", "")) == pr.get("tp_pairs", "")
                      and (fj == pr["bip_F_jaccard"] or (fj.startswith("NA") and pr["bip_F_jaccard"] == "NA")))
                say(f"[check vs prior] {name} {tname} {fam}: mine n {mine[0]} pred {mine[1]} tp {mine[2]} F_jac {mine[3]} | prior n "
                    f"{pr['n_genes']} pred {pr.get('pred_pairs')} tp {pr.get('tp_pairs')} F_jac {pr['bip_F_jaccard']} -> {'same' if ok else 'DIFFERENT'}")

    # ---- (A) on U, same genes
    for tname, truth in TRUTH.items():
        fams = ("NPIP",) if tname.startswith("literature L1") else ("NPIP", "TBC1D3", "pooled")
        for fam in fams:
            F = fam_genes(fam)
            sets = {"U all": F}
            for rname, rlab in REPORT.items():
                sets[f"U ∩ {rname}"] = F & set(rlab)
            for sname, S in sets.items():
                for lname, lab in list(LAT.items()) + ([(sname[4:], REPORT[sname[4:]])] if sname != "U all" else []):
                    s = score_counts(lab, truth, S)
                    rows.append({"scope": "U", "truth": tname, "family": fam, "gene_set": sname, "layer": lname,
                                 "kind": "report-only" if lname in REPORT else "lattice", **s})
    # ---- (B) on V (closure), pooled, lattice only
    for tname in ("HGNC gene group", "Soto family (flag ok)"):
        for lname, lab in LAT.items():
            s = score_counts(lab, TRUTH[tname], V)
            rows.append({"scope": "V (L0 closure)", "truth": tname, "family": "pooled", "gene_set": "V", "layer": lname,
                         "kind": "lattice", **s})
    cols = ["scope", "truth", "family", "gene_set", "layer", "kind", "n_genes", "truth_pairs", "pred_pairs", "tp_pairs",
            "pair_precision", "pair_recall", "bip_R_jaccard", "bip_P_jaccard", "bip_F_jaccard"]
    write(f"{OUT}/truth.tsv", rows, cols)
    say(f"[write] truth.tsv rows {len(rows)}")

    # ---- (C) in-group precision: every labelled gene of the anchor's group (members and pulled-in non-members)
    NAME = {g: r["name"] for g, r in nodes.items()}
    ANCH = {"NPIP": next(g for g in V if NAME[g] == "NPIPB2"), "TBC1D3": next(g for g in V if NAME[g] == "TBC1D3")}
    MEMS = {g for g, r in nodes.items() if r["is_member"] == "yes"}
    grows = []
    for lname, lab in LAT.items():
        for fam, anc in ANCH.items():
            S = {g for g in V if lab[g] == lab[anc]}
            for tname in ("Soto family (flag ok)", "HGNC gene group"):
                truth = TRUTH[tname]
                lg = sorted(g for g in S if g in truth)
                fams = collections.Counter(truth[g] for g in lg)
                pairs = c2(len(lg))
                same = sum(c2(n) for n in fams.values())
                comp = []
                for f, n in fams.most_common():
                    gs = sorted(NAME[g] for g in lg if truth[g] == f)
                    comp.append(f"{f}:{n}[{','.join(gs[:12])}{',...' if len(gs) > 12 else ''}]")
                grows.append({"layer": lname, "anchor_family": fam, "group_size": len(S), "group_members": len(S & MEMS),
                              "truth": tname, "labelled_genes": len(lg), "labelled_members": len(set(lg) & MEMS),
                              "labelled_pairs": pairs, "same_family_pairs": same,
                              "in_group_pair_precision": same / pairs if pairs else "NA", "n_truth_families": len(fams),
                              "truth_family_composition": " ".join(comp)})
    write(f"{OUT}/truth_ingroup.tsv", grows)
    say(f"[write] truth_ingroup.tsv rows {len(grows)}")
    say.dump(f"{OUT}/truth.out")


def cmd_lattice_filtration(args):
    """was bench/layer_order/lattice_filtration.py (--l1 c2_loose -> lattice-filtration --l1 c2_loose)

    Nested edge-test lattice, step 5: identity filtration inside L1 (a threshold filtration = single linkage on FIXED
    evidence; it is not a 'more information' view: T3(b) is not tested here).

    For each member-holding L1 group (NPIP, TBC1D3) and each threshold t in (L1 itself, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99,
    1.00): components of the L1 edges inside the group whose identity field w >= t, with the shared-exon clause OFF
    (L1 AND w >= t) and ON (L2 AND w >= t). An edge whose field is missing passes no threshold (unsatisfiable, not imputed).
    Fields (all compared unrounded):
      w98_gapexcl     §0★★★.1 single-record w_98, gap-excluded (the primary L3 field: ON at 0.98 = L3 exactly). Monotone.
      w98_gapincl     single-record w_98, gap-inclusive.
      pooled_gapexcl  identity pooled over the E1 records, gap-excluded (the 17:12 report's field). NOT monotone.
      pooled_gapincl  pooled, gap-inclusive. NOT monotone.
    w_98 needs a record witnessing shared-exon >= 0.30 of the smaller exon union, so its OFF sweep is not shared-exon-free;
    the nested listing therefore shows OFF with pooled_gapexcl and ON with w98_gapexcl (and ON with pooled_gapexcl).

    Single-linkage view: for genes x, y of the group, b(x, y) = max over paths of the minimum edge identity (bottleneck; from a
    maximum spanning forest). Components at threshold t are exactly the classes of b >= t, so the thresholds form a dendrogram
    (nested by construction). A reference group G (literature subfamily, or a clause-5 split) inside its reference set R
    APPEARS at t iff it is whole (min over pairs in G of b >= t) and separated from R \\ G (max over x in G, y in R \\ G of
    b < t); the appearance interval is (sep, whole], empty when sep >= whole (G fragments before it separates).
    Reference sets: NPIP = literature records (Dishuck 2025 truth) that are E1 catalog nodes (NPIPB1P excluded: outside both
    catalogs); TBC1D3 = the 9 clause-5 C_tree leaves.

    Outputs: lattice/filtration.txt (nested listing + appearance table), lattice/filtration_groups.tsv,
    lattice/filtration_appearance.tsv.
    usage: lattice-filtration                 primary L1 (clause-2 approx with v-exon overlap + strand check)
           lattice-filtration --l1 c2_loose   the 17:12 report's L1 (no v-exon/strand requirements), unrounded fields;
                                                 outputs get the suffix .17_03_tests_exact
    """
    OUT = LC.OUT
    L1_OPT = args.l1
    out = Log(flush=False)

    FIELDS = ("w98_gapexcl", "w98_gapincl", "pooled_gapexcl", "pooled_gapincl")

    SUF = "" if L1_OPT == "c2" else ".17_03_tests_exact"
    GRID = [None, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99, 1.00]
    nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
    grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
    GCOL = {"c2": "primary|L1", "c2_loose": "L1=c2_no_vexon_no_strand|L1"}[L1_OPT]
    NAME = {g: r["name"] for g, r in nodes.items()}
    BYNAME = collections.defaultdict(list)
    for g, n in NAME.items():
        BYNAME[n].append(g)
    MEM = {g for g, r in nodes.items() if r["is_member"] == "yes"}
    FAM = {g: nodes[g]["member_family"] for g in MEM}
    # L1 edges with identity fields / shared-exon attribute: (a, b, {field: value}, sef)
    E1 = []
    with open(f"{OUT}/edges.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = tests(r, l1=L1_OPT)
            if t[1]:
                E1.append((r["gene_a"], r["gene_b"], {f: fnum(r[ID_COL[f]]) for f in FIELDS}, fnum(r["d_shared_exon_frac"])))
                assert (t[3] == (t[2] and E1[-1][2]["w98_gapexcl"] is not None and E1[-1][2]["w98_gapexcl"] >= 0.98))
    L1lab = {g: grp[g][GCOL] for g in nodes}
    anchor = {"NPIP": BYNAME["NPIPB2"][0], "TBC1D3": BYNAME["TBC1D3"][0]}
    REF = {
        "NPIP": {g for g in nodes if nodes[g]["lit_in_truth"] == "yes" and FAM.get(g) == "NPIP" and nodes[g]["catalog"] != "none"},
        "TBC1D3": {g for g in nodes if nodes[g]["Ctree_top"] and FAM.get(g) == "TBC1D3"},
    }

    def names(S):
        return sorted(NAME[g] for g in S)

    def short(n):
        return n.replace("NPIP", "").replace("TBC1D3", "T3") if n.startswith(("NPIP", "TBC1D3")) else n

    def edges_at(S, t, sef_on, which):
        for a, b, w, sef in E1:
            if a in S and b in S:
                if sef_on and (sef is None or sef < SEF_MIN):
                    continue
                if t is not None:
                    v = w[which]
                    if v is None or v < t:
                        continue
                yield a, b

    def bottleneck(S, sef_on, which):
        """b(x, y) for all x, y in S via Kruskal on the field (descending); returns dict of dicts restricted to ref ∪ members."""
        es = sorted(((e[2][which], e[0], e[1]) for e in E1 if e[0] in S and e[1] in S and e[2][which] is not None
                     and (not sef_on or (e[3] is not None and e[3] >= SEF_MIN))), reverse=True)
        uf = UF(S)
        members = {g: {g} for g in S}
        b = collections.defaultdict(dict)
        keep = set(REF["NPIP"]) | set(REF["TBC1D3"]) | MEM
        for w, x, y in es:
            rx, ry = uf.find(x), uf.find(y)
            if rx == ry:
                continue
            A = [g for g in members[rx] if g in keep]
            B = [g for g in members[ry] if g in keep]
            for p in A:
                for q in B:
                    b[p][q] = w
                    b[q][p] = w
            uf.union(rx, ry)
            r = uf.find(rx)
            members[r] = members[rx] | members[ry]
        return b

    def appear(G, R, b):
        G = set(G) & R
        rest = R - G
        whole = min((b[x].get(y, float("-inf")) for x in G for y in G if x < y), default=float("inf"))
        sep = max((b[x].get(y, float("-inf")) for x in G for y in rest), default=float("-inf"))
        return sep, whole

    def fmt(v):
        if v == float("inf"):
            return "+inf"
        if v == float("-inf"):
            return "-inf"
        return f"{v:.6f}"

    rows = []
    for fam in ("NPIP", "TBC1D3"):
        S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
        fam_members = S & MEM
        out(f"=== {fam}: member-holding L1 group ({GCOL}), {len(S)} genes ({len(fam_members)} members); reference set for appearance: "
            f"{len(REF[fam])} genes {names(REF[fam])}")
        for sef_on, field in ((False, "pooled_gapexcl"), (True, "w98_gapexcl"), (True, "pooled_gapexcl")):
            out(f"--- {fam}, shared-exon clause {'ON (L2 AND w >= t)' if sef_on else 'OFF (L1 AND w >= t)'}; field = {field}")
            for t in GRID:
                lab = components(S, list(edges_at(S, t, sef_on, field)))
                G = groups(lab)
                memg = sorted((C for C in G.values() if C & MEM), key=lambda C: (-len(C & MEM), -len(C), names(C)))
                nonmem_only = [C for C in G.values() if not (C & MEM) and len(C) >= 2]
                tag = "L1" if t is None and not sef_on else ("L2" if t is None else f"w>={t:.2f}")
                if sef_on and t == 0.98 and field == "w98_gapexcl":
                    tag += " (= L3)"
                venn = " ".join("{" + " ".join(short(n) for n in names(C & MEM)) + (f" +{len(C - MEM)}" if C - MEM else "") + "}"
                                for C in memg)
                out(f"  [{tag}] groups holding members: {len(memg)}; non-member groups (>=2): {len(nonmem_only)}")
                out(f"     {venn}")
                for C in memg:
                    if C - MEM and len(C - MEM) <= 25:
                        out(f"       non-members with {{{' '.join(short(n) for n in names(C & MEM))[:60]}...}}: {', '.join(names(C - MEM))}")
                    elif C - MEM:
                        out(f"       non-members with {{{' '.join(short(n) for n in names(C & MEM))[:60]}...}}: {len(C - MEM)} genes, "
                            f"e.g. {', '.join(names(C - MEM)[:12])}")
                for C in memg:
                    rows.append({"family": fam, "shared_exon": "on" if sef_on else "off", "field": field, "threshold": tag,
                                 "group_members": ",".join(names(C & MEM)), "n_members": len(C & MEM), "size": len(C),
                                 "nonmembers": ",".join(names(C - MEM)) if len(C - MEM) <= 200 else f"{len(C - MEM)} genes"})
            out()

    write(f"{OUT}/filtration_groups{SUF}.tsv", rows)

    # ---------------------------------------------------------------------------------------------- appearance thresholds
    lit = {g: nodes[g] for g in REF["NPIP"]}
    REFGROUPS = [
        ("NPIP", "NPIPA (lit L1)", {g for g in REF["NPIP"] if lit[g]["lit_level1"] == "NPIPA"}),
        ("NPIP", "NPIPB (lit L1)", {g for g in REF["NPIP"] if lit[g]["lit_level1"] == "NPIPB"}),
        ("NPIP", "A6-9 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "A6-9"}),
        ("NPIP", "B3-5 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B3-5"}),
        ("NPIP", "B6-9 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B6-9"}),
        ("NPIP", "B12/13 (lit L2)", {g for g in REF["NPIP"] if lit[g]["lit_level2"] == "B12/13"}),
        ("NPIP", "named NPIPB {B3,B4,B5,B11,B12,B13}", {g for g in REF["NPIP"] if lit[g]["lit_named_npipb"] == "yes"}),
        ("TBC1D3", "clause-5 {B,F,G,H}", {BYNAME[n][0] for n in ("TBC1D3B", "TBC1D3F", "TBC1D3G", "TBC1D3H")}),
        ("TBC1D3", "clause-5 {TBC1D3,D,E,K}", {BYNAME[n][0] for n in ("TBC1D3", "TBC1D3D", "TBC1D3E", "TBC1D3K")}),
    ]
    arows = []
    out("=== Appearance thresholds (G appears at t iff whole(G) >= t > sep(G); grid first hit and exact interval)")
    for which, sef_on in (("pooled_gapexcl", False), ("pooled_gapexcl", True), ("pooled_gapincl", False),
                         ("w98_gapexcl", True), ("w98_gapincl", True)):
        if True:
            B = {}
            for fam in ("NPIP", "TBC1D3"):
                S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
                B[fam] = bottleneck(S, sef_on, which)
            for fam, gname, G in REFGROUPS:
                sep, whole = appear(G, REF[fam], B[fam])
                first = next((t for t in GRID[1:] if whole >= t > sep), None)
                arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": fam, "group": gname,
                              "genes": ",".join(names(G)), "sep_max_bottleneck_to_rest": fmt(sep),
                              "whole_min_bottleneck_inside": fmt(whole),
                              "appears_interval": f"({fmt(sep)}, {fmt(whole)}]" if whole > sep else "never (fragments before it separates)",
                              "first_grid_threshold": "none" if first is None else f"{first:.2f}"})
            # A vs B as a split (both sides separated from each other)
            A = next(G for f, n, G in REFGROUPS if n.startswith("NPIPA"))
            Bb = next(G for f, n, G in REFGROUPS if n.startswith("NPIPB"))
            ab = max((B["NPIP"][x].get(y, float("-inf")) for x in A for y in Bb), default=float("-inf"))
            arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": "NPIP",
                          "group": "NPIPA | NPIPB split (no A-B pair together)", "genes": "",
                          "sep_max_bottleneck_to_rest": fmt(ab), "whole_min_bottleneck_inside": "",
                          "appears_interval": f"t > {fmt(ab)}",
                          "first_grid_threshold": next((f"{t:.2f}" for t in GRID[1:] if t > ab), "none")})
            tb = [G for f, n, G in REFGROUPS if f == "TBC1D3"]
            x = max((B["TBC1D3"][p].get(q, float("-inf")) for p in tb[0] for q in tb[1]), default=float("-inf"))
            arows.append({"identity": which, "shared_exon": "on" if sef_on else "off", "family": "TBC1D3",
                          "group": "{B,F,G,H} | {TBC1D3,D,E,K} split (no cross pair together)", "genes": "",
                          "sep_max_bottleneck_to_rest": fmt(x), "whole_min_bottleneck_inside": "", "appears_interval": f"t > {fmt(x)}",
                          "first_grid_threshold": next((f"{t:.2f}" for t in GRID[1:] if t > x), "none")})
    for r in arows:
        out(f"  [{r['identity']}, shared-exon {r['shared_exon']}] {r['family']} {r['group']}: sep {r['sep_max_bottleneck_to_rest']}, "
            f"whole {r['whole_min_bottleneck_inside']} -> appears {r['appears_interval']}; first grid threshold {r['first_grid_threshold']}")
    write(f"{OUT}/filtration_appearance{SUF}.tsv", arows)

    # ---------------------------------------------------------------------------------------------- member dendrogram (merge heights)
    out()
    out("=== Single-linkage merge heights among members (pooled gap-excluded identity with shared-exon OFF; w98 gap-excluded with "
        "shared-exon ON): each line = a merge of two member sets at bottleneck h (paths may pass through non-members of the L1 group)")
    for fam in ("NPIP", "TBC1D3"):
        S = {g for g in nodes if L1lab[g] == L1lab[anchor[fam]]}
        M = sorted(S & MEM, key=lambda g: NAME[g])
        for sef_on, field in ((False, "pooled_gapexcl"), (True, "w98_gapexcl")):
            b = bottleneck(S, sef_on, field)
            pairs = sorted(((b[x].get(y, float("-inf")), x, y) for i, x in enumerate(M) for y in M[i + 1:]), reverse=True)
            uf = UF(M)
            sets = {g: {g} for g in M}
            out(f"--- {fam} members ({len(M)}), shared-exon {'ON' if sef_on else 'OFF'}, field {field}")
            for h, x, y in pairs:
                rx, ry = uf.find(x), uf.find(y)
                if rx == ry:
                    continue
                A, Bs = sets[rx], sets[ry]
                if h == float("-inf"):
                    roots = {uf.find(g) for g in M}
                    out("   never joined by any identity-bearing path: " + " | ".join(
                        "{" + " ".join(short(NAME[g]) for g in sorted(sets[r0], key=lambda g: NAME[g])) + "}" for r0 in sorted(roots)))
                    break
                out(f"   h={h:.6f}: {{{' '.join(short(NAME[g]) for g in sorted(A, key=lambda g: NAME[g]))}}} + "
                    f"{{{' '.join(short(NAME[g]) for g in sorted(Bs, key=lambda g: NAME[g]))}}}")
                uf.union(rx, ry)
                r = uf.find(rx)
                sets[r] = A | Bs
    out.dump(f"{OUT}/filtration{SUF}.txt")


def cmd_lattice_check_c2(args):
    """was bench/layer_order/lattice_check_c2.py

    Independent check of lattice-edges' clause-2 GENE-BODY re-implementation: for every PAF pair touching the NPIP and
    TBC1D3 member-holding L1 groups (plus a random sample of other catalog pairs), write the pair's records in BOTH
    orientations to a scratch PAF and run the SHIPPED bench/guided_pipeline.gene_body_chains on it; a pair passes iff
    >= 1 chain is returned (chains are only emitted when identity >= 0.80 and aligned >= 0.50 x min(query, extrapolated
    span)). Compare with edges.tsv d_c2x_gb_chain_raw (the shipped denominator, before the v-exon-overlap and strand
    requirements that bench/denovo_shared_def.py applied to the returned chains; the primary d_c2_genebody additionally
    uses min(body u, body v)). Pairs checked: every PAF pair touching a member-holding L1 group of the primary, the c2x
    variant or the no-v-exon/no-strand variant (the largest L1 groups), plus 3,000 random other catalog pairs.
    Output: lattice/check_c2.out

    Second check (correction pass 2026-09-16): the target-side requirements. On the SAME shipped chains, apply the
    shipped bench/denovo_shared_def.py cmd_families edge loop verbatim in body coordinates (ExonIndex on the two nodes'
    exon unions, orient = u strand for a '+' chain and flipped for '-', skip v when both nodes are spliced and v's strand
    differs from orient) and compare with edges.tsv d_c2x_genebody (shipped denominator + v-exon overlap + strand check).
    The primary d_c2_genebody runs the same target/strand code with the shorter-body denominator.

    Wave 7: the old script had not run since 2026-09-19. It imported denovo_shared_def, archived in cleanup wave 1
    (567e092e), for ExonIndex alone; and it exec'd the head of lattice_edges.py, which has needed __file__ since §6r9
    (a94ff6d1). ExonIndex is now lattice_common.IntervalIndex (the same code) and the head is
    lattice_common.catalog_context (the same code), so nothing is exec'd. denovo_shared_def.py itself is at tag
    notebook-2026-09-19:archive/bench/denovo_shared_def.py.
    """
    import guided_pipeline as gp  # the SHIPPED gene_body_chains is what this stage checks (bench/guided_pipeline.py)
    OUT = LC.OUT
    # node exon unions per body key, exactly as lattice-edges builds them
    ctx = catalog_context(lambda s: print(s, flush=True))
    KEY_BLOCKS = {f"{k[0]}:{k[1]}-{k[2]}": v for k, v in ctx.key_blocks.items()}

    rng = random.Random(20260916)  # = the old random.seed(20260916) + random.sample
    grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
    anch = {r["name"]: r["gene_id"] for r in grp.values()}
    L1 = {g: r["primary|L1"] for g, r in grp.items()}
    L1x = {g: r["L1=c2x_extrapolated|L1"] for g, r in grp.items()}
    L1l = {g: r["L1=c2_no_vexon_no_strand|L1"] for g, r in grp.items()}
    fam_groups = {L1[anch["NPIPB2"]], L1[anch["TBC1D3"]]}
    fam_groups_x = {L1x[anch["NPIPB2"]], L1x[anch["TBC1D3"]]}
    fam_groups_l = {L1l[anch["NPIPB2"]], L1l[anch["TBC1D3"]]}
    genes = ctx.genes  # light/work/refseq/genes.tsv (the old script hardcoded the default root's copy; audit F9)
    key_of_gene = {g: f"{r['chrom']}:{int(r['start0']) + 1}-{r['end']}" for g, r in genes.items()}

    want = {}
    want_t = {}
    other = []
    with open(f"{OUT}/edges.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["d_evidence"] != "yes":
                continue
            k = frozenset((key_of_gene[r["gene_a"]], key_of_gene[r["gene_b"]]))
            if (L1[r["gene_a"]] in fam_groups or L1[r["gene_b"]] in fam_groups or L1x[r["gene_a"]] in fam_groups_x
                    or L1x[r["gene_b"]] in fam_groups_x or L1l[r["gene_a"]] in fam_groups_l or L1l[r["gene_b"]] in fam_groups_l):
                want[k] = r["d_c2x_gb_chain_raw"]
                want_t[k] = (r["d_c2x_genebody"], genes[r["gene_a"]]["strand"], genes[r["gene_b"]]["strand"], key_of_gene[r["gene_a"]])
            else:
                other.append((k, r["d_c2x_gb_chain_raw"], (r["d_c2x_genebody"], genes[r["gene_a"]]["strand"],
                                                            genes[r["gene_b"]]["strand"], key_of_gene[r["gene_a"]])))
    for k, v, vt in rng.sample(other, min(3000, len(other))):
        want[k] = v
        want_t[k] = vt
    lines = collections.defaultdict(list)
    for cat, path in PAF.items():
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == f[5]:
                    continue
                k = frozenset((f[0], f[5]))
                if k in want:
                    lines[k].append(f)
    tmp = f"{OUT}/check_c2.tmp.paf"
    agree = disagree = 0
    agree_t = disagree_t = 0
    examples = []
    FLIPS = {"+": "-", "-": "+"}
    for k, exp in want.items():
        with open(tmp, "w") as out:
            for f in lines[k]:
                out.write("\t".join(f) + "\n")
                sw = [f[5], f[6], f[7], f[8], f[4], f[0], f[1], f[2], f[3]] + f[9:]
                out.write("\t".join(sw) + "\n")
        chains = gp.gene_body_chains(tmp)
        got = "yes" if chains else "no"
        if got == exp:
            agree += 1
        else:
            disagree += 1
            if len(examples) < 10:
                examples.append((sorted(k), exp, got))
        # shipped cmd_families edge loop on these chains (body coordinates; chrom field = the target body key)
        exp_t, strand_a, strand_b, key_a = want_t[k]
        ka, kb = sorted(k)
        strand = {key_a: strand_a, (kb if key_a == ka else ka): strand_b}
        nd = [{"idx": i, "chrom": kk, "strand": strand[kk], "exons": [(s0 - (int(kk.rsplit(":", 1)[1].split("-")[0]) - 1),
                                                                        e0 - (int(kk.rsplit(":", 1)[1].split("-")[0]) - 1))
                                                                       for s0, e0 in KEY_BLOCKS[kk]]} for i, kk in enumerate((ka, kb))]
        idx = IntervalIndex.from_nodes(nd)  # (was denovo_shared_def.ExonIndex(nd))
        spliced = {n["idx"]: len(n["exons"]) >= 2 for n in nd}
        by_key = {n["chrom"]: n for n in nd}
        edge = False
        for c in chains:
            u = by_key[c["q"]]
            orient = u["strand"] if c["strand"] == "+" else FLIPS.get(u["strand"], u["strand"])
            for v in idx.hits(c["chrom"], c["s"], c["e"]):
                if v == u["idx"]:
                    continue
                if spliced[u["idx"]] and spliced[v] and nd[v]["strand"] != orient:
                    continue
                edge = True
        got_t = "yes" if edge else "no"
        if got_t == exp_t:
            agree_t += 1
        else:
            disagree_t += 1
            if len(examples) < 20:
                examples.append(("target/strand", sorted(k), exp_t, got_t))
    os.remove(tmp)
    msg = [f"pairs checked {len(want)} (member-holding L1 groups: {len(want) - min(3000, len(other))}; random other catalog pairs "
           f"{min(3000, len(other))}); raw chain exists (d_c2x_gb_chain_raw): agree {agree}; disagree {disagree}",
           f"shipped cmd_families v-exon overlap + strand check on those chains vs d_c2x_genebody: agree {agree_t}; disagree "
           f"{disagree_t}"] + [f"  example {e}" for e in examples]
    print("\n".join(msg))
    open(f"{OUT}/check_c2.out", "w").write("\n".join(msg) + "\n")


def cmd_lattice_report(args):
    """was bench/layer_order/lattice_report_tables.py

    Nested edge-test lattice, step 6: every table quoted in bench/NESTED_LATTICE_NPIP_TBC1D3.md, regenerated from the
    lattice/ result files (no new computation except shortest paths, simple counts and component labels over edges.tsv).
    All numbers are formatted here from UNROUNDED values (edges.tsv and truth.tsv store full precision).
    Output: lattice/report_tables.md
    """
    OUT = LC.OUT

    L = []

    def p(s=""):
        L.append(s)

    def f3(x, nd=3):
        if isinstance(x, str) and (x.startswith("NA") or x == ""):
            return "NA"
        v = fnum(x) if isinstance(x, str) else x
        return "NA" if v is None else f"{v:.{nd}f}"

    nodes = {r["gene_id"]: r for r in tsv(f"{OUT}/nodes.tsv")}
    NAME = {g: r["name"] for g, r in nodes.items()}
    BYN = {r["name"]: g for g, r in nodes.items()}
    MEM = {g for g, r in nodes.items() if r["is_member"] == "yes"}
    V = set(nodes)

    NEED = ["gene_a", "gene_b", "name_a", "name_b", "chrom_a", "chrom_b", "same_locus", "p_evidence", "p_aa_identity",
            "p_cov_longer", "p_qualifies_6ko", "p_aa50", "p_qualifies_union", "d_evidence", "d_catalog", "d_e1_records",
            "d_body_bp_a", "d_body_bp_b", "d_exon_union_bp_a", "d_exon_union_bp_b", "d_e1_identity", "d_e1_identity_gapexcl",
            "d_e1_cov_longer", "d_e1_edge", "d_shared_exon_bp", "d_shared_exon_frac", "d_shared_exon_frac_allrec",
            "d_w98_gapexcl", "d_w98_gapincl", "d_both_spliced", "d_c2_genebody", "d_c2_gb_best_frac", "d_c2_gb_chain_identity",
            "d_c2_exon", "d_c2_exon_best_frac", "d_c2_approx", "d_c2nostrand_approx", "d_c2exontgt_approx",
            "d_c2loose_genebody", "d_c2loose_gb_best_frac", "d_c2loose_gb_chain_identity", "d_c2loose_exon",
            "d_c2loose_exon_best_frac", "d_c2loose_approx", "d_c2x_approx", "s2_edge", "s2_max_identity", "ctree_family"]
    rows = []
    with open(f"{OUT}/edges.tsv") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ix = [(c, hdr.index(c)) for c in NEED]
        for f in rd:
            rows.append({c: f[i] for c, i in ix})
    EDG = {frozenset((r["name_a"], r["name_b"])): r for r in rows}
    T = [tests(r) for r in rows]  # primary tests per row
    T_OLD = [tests(r, l1="c2_loose", id_which="pooled_gapexcl") for r in rows]  # the 17:12 report's tests, unrounded

    def edges_where(pred):
        return [(r["gene_a"], r["gene_b"]) for r in rows if pred(r)]

    def lab_of(edges):
        return components(V, edges)

    def grp_of(lab, name):
        return {g for g in V if lab[g] == lab[BYN[name]]}

    # ------------------------------------------------------------------ edge table summary
    cnt = collections.Counter()
    for r, t, to in zip(rows, T, T_OLD):
        cnt["rows"] += 1
        cnt["same_locus"] += r["same_locus"] == "yes"
        cnt["protein HSP pair"] += r["p_evidence"] == "yes"
        cnt["protein §6ko-qualifying (shipped greedy cover)"] += r["p_qualifies_6ko"] == "yes"
        cnt["protein union-cover qualifying (§0★★★.1 form; closure row only)"] += r["p_qualifies_union"] == "yes"
        cnt["protein §6ko at aa >= 0.50"] += r["p_aa50"] == "yes"
        cnt["PAF pair"] += r["d_evidence"] == "yes"
        cnt["PAF pair with both genes spliced (strand check applies)"] += r["d_both_spliced"] == "yes"
        cnt["E1 edge (as built)"] += r["d_e1_edge"] == "yes"
        cnt["clause-2 approx, PRIMARY (v-exon overlap + strand check)"] += r["d_c2_approx"] == "yes"
        cnt["  gene-body clause (primary)"] += r["d_c2_genebody"] == "yes"
        cnt["  exon proxy clause (primary)"] += r["d_c2_exon"] == "yes"
        cnt["clause-2 approx, v-exon overlap without strand check"] += r["d_c2nostrand_approx"] == "yes"
        cnt["clause-2 approx, v-exon overlap on the exon proxy only"] += r["d_c2exontgt_approx"] == "yes"
        cnt["clause-2 approx, no v-exon overlap, no strand check (17:03 primary)"] += r["d_c2loose_approx"] == "yes"
        cnt["  gene-body clause (loose)"] += r["d_c2loose_genebody"] == "yes"
        cnt["  exon proxy clause (loose)"] += r["d_c2loose_exon"] == "yes"
        cnt["clause-2 approx, finder denominator + v-exon overlap + strand check (variant)"] += r["d_c2x_approx"] == "yes"
        cnt["PAF pair without E1-qualifying record"] += r["d_evidence"] == "yes" and r["d_e1_records"] == "0"
        cnt["w_98 gap-excluded defined (a record witnesses sx >= 0.30)"] += fnum(r["d_w98_gapexcl"]) is not None
        cnt["S2 edge"] += r["s2_edge"] == "yes"
        cnt["clause-5 C_tree pair annotated"] += r["ctree_family"] != ""
        for k, ok in zip(LEVELS, t):
            cnt[f"passes t_{k[1]} (primary)"] += ok
        cnt["passes t_1 AND f_ex over ALL records >= 0.30 (L2 with the definition's f_ex)"] += t[1] and (fnum(r["d_shared_exon_frac_allrec"]) or 0) >= 0.30
        cnt["passes t_1 with d_catalog NA (cross-trio DNA evidence)"] += t[1] and r["d_catalog"] == "NA"
        for k, ok in zip(LEVELS, to):
            cnt[f"passes t_{k[1]} under the 17:12 report's tests, unrounded (L1 loose, pooled gap-excluded L3)"] += ok
        cnt["17:12 tests: L3 at 4-dp-rounded pooled identity (the 17:12 count)"] += to[2] and (fnum(r["d_e1_identity_gapexcl"]) is not None and round(fnum(r["d_e1_identity_gapexcl"]), 4) >= 0.98)
    p("## T-edges: edge table lattice/edges.tsv")
    p("| quantity | pairs |")
    p("|---|---|")
    for k, v in cnt.items():
        p(f"| {k} | {v:,} |")
    p()
    p("Edges whose pooled gap-excluded identity is in [0.97995, 0.98) and that pass t_2 under the 17:12 tests "
      "(admitted at L3 by 4-dp rounding):")
    for r, to in zip(rows, T_OLD):
        v = fnum(r["d_e1_identity_gapexcl"])
        if to[2] and v is not None and 0.97995 <= v < 0.98:
            p(f"- {r['name_a']}–{r['name_b']}: pooled gap-excl {v:.7f}; w_98 gap-excl {f3(r['d_w98_gapexcl'], 7)}")
    p()

    # ------------------------------------------------------------------ closure
    p("## T-closure: lattice/closure.tsv")
    p("| closure | genes | new genes per hop | never-searched proteins | no §6ko protein | outside both E1 catalogs | outside V |")
    p("|---|---|---|---|---|---|---|")
    for r in tsv(f"{OUT}/closure.tsv"):
        p(f"| {r['closure']} | {r['genes']} | {r['new_genes_per_hop']} | {r.get('never_searched_proteins', '')} | {r.get('no_protein', '')} | "
          f"{r.get('outside_E1_catalogs', '')} | {r.get('outside_V', '')} |")
    p()

    # ------------------------------------------------------------------ levels: member-holding groups
    mg = tsv(f"{OUT}/member_groups.tsv")
    p("## T-levels: member-holding groups (lattice/member_groups.tsv); cells = size (members / non-members); member singletons "
      "listed separately; every multi-gene member-holding group is listed")
    p("| variant | level | edges | NPIP group(s) | TBC1D3 group(s) | member singletons |")
    p("|---|---|---|---|---|---|")
    lv = {(r["variant"], r["level"]): r for r in tsv(f"{OUT}/levels.tsv")}
    by = collections.defaultdict(list)
    for r in mg:
        by[(r["variant"], r["level"])].append(r)
    for (v, k), rs in by.items():
        npip = [r for r in rs if int(r["n_NPIP_members"]) > 0 and int(r["size"]) > 1]
        tbc = [r for r in rs if int(r["n_TBC1D3_members"]) > 0 and int(r["size"]) > 1]
        single = [r["members"] for r in rs if r["size"] == "1"]
        both = [r for r in rs if int(r["n_NPIP_members"]) > 0 and int(r["n_TBC1D3_members"]) > 0]

        def cell(lst):
            return "; ".join(f"{r['size']} ({int(r['n_NPIP_members']) + int(r['n_TBC1D3_members'])} / {r['n_nonmembers']})"
                             + (f" [{r['members']}]" if int(r['size']) <= 3 or (int(r['n_NPIP_members']) + int(r['n_TBC1D3_members'])) <= 2 else "")
                             for r in lst)
        p(f"| {v} | {k} | {lv[(v, k)]['edges']} | {cell(npip)} | {cell([r for r in tbc if r not in both]) or ('same group' if both else '')} | "
          f"{len(single)}: {', '.join(sorted(single))} |")
    p()

    # ------------------------------------------------------------------ non-members of the primary member-holding groups
    p("## T-nonmembers: primary member-holding groups (L1-L3): chromosomes, biotypes, name classes, non-members")
    for k in ("L0", "L1", "L2", "L3"):
        lab = lab_of([(r["gene_a"], r["gene_b"]) for r, t in zip(rows, T) if t[LEVELS.index(k)]])
        for anc in ("NPIPB2", "TBC1D3"):
            S = grp_of(lab, anc)
            non = S - MEM
            ch = collections.Counter(nodes[g]["chrom"] for g in sorted(S))  # sorted: deterministic tie order in most_common
            bt = collections.Counter(nodes[g]["biotype"] for g in sorted(S))
            cls = collections.Counter()
            for g in sorted(non):
                n = NAME[g]
                for pre in ("ZNF", "BNIP3P", "SMG1", "PLA2G10", "BOLA2", "PKD1", "USP", "DHX40", "GOLGA", "LOC", "TBC1D", "VN1R", "SLC7A5", "PDXDC"):
                    if n.startswith(pre):
                        cls[pre] += 1
                        break
            p(f"- {k} {anc}: {len(S)} genes ({len(S & MEM)} members); chrom {dict(ch.most_common(6))}; biotypes {dict(bt.most_common(4))}; "
              f"never-searched proteins {sum(1 for g in S if nodes[g]['protein_searched'] == 'no')}; outside both E1 catalogs "
              f"{sum(1 for g in S if nodes[g]['catalog'] == 'none')}; name classes of non-members {dict(cls.most_common())}")
            if len(non) <= 160:
                p(f"  - non-members: {', '.join(sorted(NAME[g] for g in non))}")
    p()

    # ------------------------------------------------------------------ L0 variants and path NPIP -> TBC1D3
    prot = lambda r: r["p_qualifies_6ko"] == "yes"
    aa50 = lambda r: prot(r) and (fnum(r["p_aa_identity"]) or 0) >= 0.50
    uni = lambda r: r["p_qualifies_union"] == "yes"
    d = lambda r: r["d_c2_approx"] == "yes"
    dl = lambda r: r["d_c2loose_approx"] == "yes"
    gb = lambda r: r["d_c2_genebody"] == "yes"
    exo = lambda r: r["d_c2_exon"] == "yes"
    e1 = lambda r: r["d_e1_edge"] == "yes"
    e1c2 = lambda r: e1(r) and (fnum(r["d_e1_identity"]) or 0) >= 0.80 and (fnum(r["d_e1_cov_longer"]) or 0) >= 0.50
    nsl = lambda r: r["same_locus"] != "yes"
    p("## T-L0variants: which L0 forms join NPIP (NPIPB2) and TBC1D3 (components on V, same-locus pairs excluded)")
    p("| L0 edge form | edges | joined | NPIPB2 group | TBC1D3 group |")
    p("|---|---|---|---|---|")
    for tag, pred in (("primary: P ∨ D", lambda r: prot(r) or d(r)),
                      ("(P ∧ aa ≥ 0.50) ∨ D", lambda r: aa50(r) or d(r)),
                      ("P_union-cover ∨ D", lambda r: uni(r) or d(r)),
                      ("P ∨ D gene-body disjunct only", lambda r: prot(r) or gb(r)),
                      ("P ∨ D exon-proxy disjunct only", lambda r: prot(r) or exo(r)),
                      ("P ∨ E1 as built (0.70/0.30 with exon-to-exon gate)", lambda r: prot(r) or e1(r)),
                      ("P ∨ E1 at 0.80/0.50", lambda r: prot(r) or e1c2(r)),
                      ("(P ∧ aa ≥ 0.50) ∨ E1 at 0.80/0.50", lambda r: aa50(r) or e1c2(r)),
                      ("P only", prot), ("D only (= L1)", d),
                      ("17:12: P ∨ D_loose", lambda r: prot(r) or dl(r)),
                      ("17:12: (P ∧ aa ≥ 0.50) ∨ D_loose", lambda r: aa50(r) or dl(r))):
        E = edges_where(lambda r: nsl(r) and pred(r))
        lab = lab_of(E)
        A, B = grp_of(lab, "NPIPB2"), grp_of(lab, "TBC1D3")
        p(f"| {tag} | {len(E):,} | {'yes' if lab[BYN['NPIPB2']] == lab[BYN['TBC1D3']] else 'no'} | {len(A):,} | {len(B):,} |")
    p()
    def edge_line(r):
        return (f"protein §6ko {r['p_qualifies_6ko']} (aa {f3(r['p_aa_identity'])}, cov {f3(r['p_cov_longer'])}); clause-2 {r['d_c2_approx']} "
                f"(gene-body {r['d_c2_genebody']} frac {f3(r['d_c2_gb_best_frac'])} chain id {f3(r['d_c2_gb_chain_identity'])}, exon proxy "
                f"{r['d_c2_exon']} frac {f3(r['d_c2_exon_best_frac'])}); loose clause-2 {r['d_c2loose_approx']}; bodies {r['d_body_bp_a']}/"
                f"{r['d_body_bp_b']} bp; exon unions {r['d_exon_union_bp_a']}/{r['d_exon_union_bp_b']} bp; E1 {r['d_e1_edge']}; shared-exon "
                f"{f3(r['d_shared_exon_frac'])} ({r['d_shared_exon_bp'] or 'NA'} bp); w98 gap-excl {f3(r['d_w98_gapexcl'], 4)}; pooled "
                f"gap-excl {f3(r['d_e1_identity_gapexcl'], 4)}; catalog {r['d_catalog']}; chroms {r['chrom_a']}/{r['chrom_b']}")

    for tag, pred in (("primary L0", lambda i, r: T[i][0]), ("(P ∧ aa ≥ 0.50) ∨ D", lambda i, r: nsl(r) and (aa50(r) or d(r)))):
        adj = collections.defaultdict(dict)
        for i, r in enumerate(rows):
            if pred(i, r):
                adj[r["gene_a"]][r["gene_b"]] = r
                adj[r["gene_b"]][r["gene_a"]] = r
        pp = shortest_path(adj, BYN["NPIPB2"], BYN["TBC1D3"])
        p(f"## T-L0path ({tag}): one shortest path NPIPB2 -> TBC1D3 ({len(pp) - 1 if pp else 'no'} edges; BFS with sorted neighbours)")
        for i in range(len(pp) - 1 if pp else 0):
            r = adj[pp[i]][pp[i + 1]]
            p(f"- {NAME[pp[i]]} – {NAME[pp[i + 1]]}: {edge_line(r)}")
        p()

    # ------------------------------------------------------------------ sanity
    san = tsv(f"{OUT}/sanity.tsv")
    agg = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    for r in san:
        key = r["check"]
        a = agg[key]
        a[0] += 1
        a[1] += int(r["groups_checked"])
        a[2] += int(r["violations"])
        a[3] += int(r["coarse_blocks_ge2"])
        a[4] += int(r["coarse_blocks_split"])
        a[5] += int(r["member_holding_blocks_split"])
    p("## T-sanity: lattice/sanity.tsv (non-vacuity: coarse blocks with >= 2 genes that the finer partition actually splits)")
    p("| check | checks | groups checked (all, singletons included) | violations | coarse blocks >= 2 genes | of which split | member-holding split |")
    p("|---|---|---|---|---|---|---|")
    for k, (n, g, v, nb, ns, nm) in agg.items():
        p(f"| {k} | {n} | {g:,} | {v} | {nb:,} | {ns:,} | {nm} |")
    p(f"| **total** | {len(san)} | {sum(int(r['groups_checked']) for r in san):,} | {sum(int(r['violations']) for r in san)} | | | |")
    p()
    p("Primary only, per level pair (T1) and per level (T2b, any-overlap reads >= 3):")
    for r in san:
        if r["variant"] == "primary" and (r["check"].startswith("T1") or r["check"].startswith("T2b: G_k[X] components refine G_k restricted to X (any-overlap reads>=3")):
            p(f"- {r['check'][:40]} {r['detail']}: coarse blocks >= 2: {r['coarse_blocks_ge2']}, split {r['coarse_blocks_split']}, "
              f"member-holding split {r['member_holding_blocks_split']}")
    p()

    # ------------------------------------------------------------------ chaining genes
    ch = tsv(f"{OUT}/chaining.tsv")
    p("## T-chaining: lattice/chaining.tsv (with / apart from the family anchor NPIPB2 or TBC1D3; apart (n) = size of its own group)")
    p("| gene | member | readthrough | L0 | L1 | L2 | L3 | L1Δ | L2Δ | L3Δ | L2 no-readthrough nodes | L3 no-readthrough nodes |")
    p("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in ch:
        def c(x):
            return "with" if x == "with anchor" else x.replace("apart ", "apart")
        p(f"| {r['gene']} | {r['is_member']} | {r['readthrough']} | " + " | ".join(c(r[f"primary|{k}"]) for k in LEVELS) + " | "
          + " | ".join(c(r[f"triangle|{k}"]) for k in LEVELS[1:]) + f" | {c(r['no_readthrough_nodes|L2'])} | {c(r['no_readthrough_nodes|L3'])} |")
    p()
    p("Connecting paths (primary):")
    for r in ch:
        if r["gene"] in ("PKD1", "PKD1P1", "DHX40", "USP6", "TBC1D26", "TBC1D29P", "TBC1D3P7", "LOC100505915"):
            for k in ("L1", "L2", "L3"):
                if r[f"path|{k}"]:
                    p(f"- {r['gene']} {k}: {r['path|' + k]}")
    p()

    # ------------------------------------------------------------------ Soto families of neighbours of interest
    p("## T-soto-neighbours: Soto family (nodes.tsv soto_families, flag) of genes discussed in the chaining sections")
    for n in ("NPIPA1", "NPIPA5", "NPIPA6", "NPIPA9", "NPIPB2", "NPIPB9", "PKD1P6-NPIPP1", "PKD1", "PKD1P1", "PKD1P2", "PKD1P3", "PKD1P6",
              "LOC131696449", "PKD1P3-NPIPA1", "DHX40", "DHX40P1", "TBC1D3P1-DHX40P1", "RNFT1-DT", "TBC1D3", "TBC1D3P1", "TBC1D3P2",
              "TBC1D3P3", "TBC1D3P4", "USP6", "USP32", "USP32P1", "USP32P2", "USP32P3", "USP32P4", "TBC1D26", "TBC1D28", "CA4", "TBC1D29P"):
        if n in BYN:
            g = BYN[n]
            p(f"- {n}: {nodes[g]['soto_families'] or '-'} ({nodes[g]['soto_flag']}); {nodes[g]['chrom']}; member {nodes[g]['is_member']}")
    p()

    # ------------------------------------------------------------------ truth pivot
    tr = tsv(f"{OUT}/truth.tsv")
    lat = ["L0", "L1", "L2", "L3", "L1Δ", "L2Δ", "L3Δ", "L3[w98 gap-incl]", "L3[pooled gap-excl]", "L3[pooled gap-incl]", "L3[id S2]",
           "L1[no v-exon/strand]", "17:03 L3"]
    p("## T-truth-U: bipartite F (one-to-one Jaccard) on the same genes; last column = the prior study's report-only grouping "
      "(F; pair P / R) on those genes (lattice/truth.tsv). U-scope scores only members; read T-truth-ingroup next to them.")
    for truth in ("HGNC gene group", "Soto family (flag ok)", "literature L1 (NPIPA|NPIPB)", "literature L2"):
        p(f"### {truth}")
        p("| family | gene set | n | " + " | ".join(lat) + " | report-only grouping: F (P / R) |")
        p("|---|---|---|" + "---|" * len(lat) + "---|")
        for fam in ("NPIP", "TBC1D3", "pooled"):
            for gs in ("U all", "U ∩ P (§6ko MCL)", "U ∩ D (E1 MCL)", "U ∩ C_L1 (lit, circular)", "U ∩ C_fine (lit, circular)",
                       "U ∩ Ctree_top (clause 5)", "U ∩ Ctree_min (clause 5)"):
                sel = {r["layer"]: r for r in tr if r["scope"] == "U" and r["truth"] == truth and r["family"] == fam and r["gene_set"] == gs}
                if not sel or sel[lat[0]]["n_genes"] in ("0", "1"):
                    continue
                rep = [r for r in sel.values() if r["kind"] == "report-only"]
                repc = (f"{rep[0]['layer']}: {f3(rep[0]['bip_F_jaccard'])} ({f3(rep[0]['pair_precision'])} / {f3(rep[0]['pair_recall'])})"
                        if rep else "")
                p(f"| {fam} | {gs} | {sel[lat[0]]['n_genes']} | " + " | ".join(f3(sel[x]["bip_F_jaccard"]) for x in lat) + f" | {repc} |")
        p()
    p("## T-truth-U-pairs: pair precision / recall of the lattice levels on 'U all' (lattice/truth.tsv)")
    p("| truth | family | n | " + " | ".join(lat[:7]) + " |")
    p("|---|---|---|" + "---|" * 7)
    for truth in ("HGNC gene group", "Soto family (flag ok)", "literature L1 (NPIPA|NPIPB)", "literature L2"):
        for fam in ("NPIP", "TBC1D3", "pooled"):
            sel = {r["layer"]: r for r in tr if r["scope"] == "U" and r["truth"] == truth and r["family"] == fam and r["gene_set"] == "U all"}
            if not sel or sel["L0"]["n_genes"] in ("0", "1"):
                continue
            p(f"| {truth} | {fam} | {sel['L0']['n_genes']} | " + " | ".join(f"{f3(sel[x]['pair_precision'])} / {f3(sel[x]['pair_recall'])}" for x in lat[:7]) + " |")
    p()
    latV = ["L0", "L1", "L2", "L3", "L1Δ", "L2Δ", "L3Δ", "L1[E1@.80/.50]", "L3[w98 gap-incl]", "L3[pooled gap-excl]", "L3[id S2]",
            "L1[no v-exon/strand]", "17:03 L3"]
    p("## T-truth-V: lattice levels on the closure V, pooled (recall conditioned on V; lattice/truth.tsv)")
    p("| truth | n | " + " | ".join(latV) + " |")
    p("|---|---|" + "---|" * len(latV))
    for truth in ("HGNC gene group", "Soto family (flag ok)"):
        sel = {r["layer"]: r for r in tr if r["scope"].startswith("V") and r["truth"] == truth}
        p(f"| {truth}: F | {sel['L0']['n_genes']} | " + " | ".join(f3(sel[x]["bip_F_jaccard"]) for x in latV) + " |")
        p(f"| {truth}: pair P / R | | " + " | ".join(f"{f3(sel[x]['pair_precision'])} / {f3(sel[x]['pair_recall'])}" for x in latV) + " |")
    p()
    ig = tsv(f"{OUT}/truth_ingroup.tsv")
    p("## T-truth-ingroup: in-group pair precision of the anchor's group (all labelled genes of the group, members and "
      "non-members; lattice/truth_ingroup.tsv)")
    p("| layer | family | group size (members) | Soto: labelled genes (members) | Soto pairs same / all = precision | Soto families | HGNC: labelled genes | HGNC precision |")
    p("|---|---|---|---|---|---|---|---|")
    igd = {(r["layer"], r["anchor_family"], r["truth"]): r for r in ig}
    for layer in ("L0", "L1", "L2", "L3", "L1Δ", "L2Δ", "L3Δ", "L3[w98 gap-incl]", "L3[pooled gap-excl]", "L1[no v-exon/strand]",
                  "L2[no v-exon/strand]", "L3[no v-exon/strand]", "17:03 L2", "17:03 L3"):
        for fam in ("NPIP", "TBC1D3"):
            s_, h_ = igd[(layer, fam, "Soto family (flag ok)")], igd[(layer, fam, "HGNC gene group")]
            p(f"| {layer} | {fam} | {s_['group_size']} ({s_['group_members']}) | {s_['labelled_genes']} ({s_['labelled_members']}) | "
              f"{s_['same_family_pairs']} / {s_['labelled_pairs']} = {f3(s_['in_group_pair_precision'])} | {s_['n_truth_families']} | "
              f"{h_['labelled_genes']} | {f3(h_['in_group_pair_precision'])} |")
    p()
    p("Soto family composition of the primary L2 / L3 groups:")
    for layer in ("L2", "L3", "L3Δ", "17:03 L2", "17:03 L3"):
        for fam in ("NPIP", "TBC1D3"):
            p(f"- {layer} {fam}: {igd[(layer, fam, 'Soto family (flag ok)')]['truth_family_composition']}")
    p()

    # ------------------------------------------------------------------ filtration appearance
    for suf, title in (("", "primary L1 groups"), (".17_03_tests_exact", "the 17:12 report's L1 groups (no v-exon/strand requirements), unrounded")):
        p(f"## T-appearance{suf}: lattice/filtration_appearance{suf}.tsv ({title})")
        p("| field | shared-exon | family | reference group | separated above (max bottleneck to the rest) | whole up to (min bottleneck inside) | appears | first grid threshold |")
        p("|---|---|---|---|---|---|---|---|")
        for r in tsv(f"{OUT}/filtration_appearance{suf}.tsv"):
            p(f"| {r['identity']} | {r['shared_exon']} | {r['family']} | {r['group']} | {r['sep_max_bottleneck_to_rest']} | "
              f"{r['whole_min_bottleneck_inside'] or '-'} | {r['appears_interval']} | {r['first_grid_threshold']} |")
        p()

    # ------------------------------------------------------------------ triangle drops
    p("## T-triangle: lattice/triangle_drops.tsv (member-holding groups of G_k under the 3-truss)")
    p("| level | group size | members | triangle parts | of which singletons | part holding the members: size (members) | members outside it | non-members outside it (first 40) | Soto families of a 2-gene group |")
    p("|---|---|---|---|---|---|---|---|---|")
    for r in tsv(f"{OUT}/triangle_drops.tsv"):
        soto2 = ""
        if r["is_2copy_group"] == "yes":
            soto2 = "; ".join(f"{n}: {nodes[BYN[n]]['soto_families'] or '-'} ({nodes[BYN[n]]['soto_flag']})" for n in r["members"].split(","))
        p(f"| {r['level']} | {r['size']} | {r['n_members']} | {r['triangle_parts']} | {r['triangle_singletons']} | {r['main_part_size']} "
          f"({r['main_part_n_members']}) | {r['members_outside_main_part'] or '-'} | {r['nonmembers_outside_main_part_first40'][:300]} | {soto2} |")
    p()

    # ------------------------------------------------------------------ expression
    p("## T-expr: member-holding expressed components (>= 2 genes) by expression set, variant, view and level (lattice/expr_views.tsv)")
    p("Any-overlap = primary reads (-F 2308) of any MAPQ with a block on an exon of the gene; unique = the read's blocks hit exons of "
      "exactly one RefSeq record genome-wide (lattice_expr.py); neither is the §0★★★.1 u >= 3 (MAPQ >= 1) convention.")
    p("| expression set | variant | view | level | components: size (members) | member singletons |")
    p("|---|---|---|---|---|---|")
    for r in tsv(f"{OUT}/expr_views.tsv"):
        p(f"| {r['expression_set']} | {r['variant']} | {r['view']} | {r['level']} | {r['member_components'] or '-'} | {r['member_singletons'] or '-'} |")
    p()
    p("Members and non-members of the reads >= 3 and unique >= 3 components (primary and 17:12 tests):")
    for r in tsv(f"{OUT}/expr_views.tsv"):
        if r["expression_set"] in ("any-overlap reads>=3", "unique reads>=3") and r["level"] in ("L0", "L3"):
            p(f"- {r['expression_set']} {r['variant']} {r['view']} {r['level']}: members [{r['members']}]; non-members [{r['nonmembers']}]")
    p()
    reads = {r["gene_id"]: (int(r["n_reads_any"]), int(r["n_reads_unique"])) for r in tsv(f"{OUT}/expr_counts.tsv")}
    # class (c) check at reads >= 1 on L0: do NPIP and TBC1D3 separate, and do all connecting paths pass through unexpressed copies?
    E0 = [(r["gene_a"], r["gene_b"]) for r, t in zip(rows, T) if t[0]]
    for tt in (1, 3):
        X = {g for g in V if reads[g][0] >= tt}
        labX = components(X, [e for e in E0 if e[0] in X and e[1] in X])
        a, b = BYN["NPIPB2"], BYN["TBC1D3"]
        sep = a not in X or b not in X or labX[a] != labX[b]
        p(f"- L0, any-overlap reads >= {tt}: NPIPB2 expressed {a in X}, TBC1D3 expressed {b in X}; separated in G_0[X]: {sep}; "
          f"joined in G_0: {lab_of(E0)[a] == lab_of(E0)[b]} (so every G_0 path between them uses a copy with < {tt} reads)")
    p()

    # ------------------------------------------------------------------ hubs inside member-holding groups
    body, span_exon = {}, {}
    for r in rows:
        if r["d_evidence"] == "yes":
            for side, g in (("a", r["gene_a"]), ("b", r["gene_b"])):
                body[g] = int(r[f"d_body_bp_{side}"])
                span_exon[g] = r[f"d_exon_union_bp_{side}"] == r[f"d_body_bp_{side}"]
    p("## T-hubs: L1 member-holding groups, internal edges and highest-degree nodes (primary and the 17:12 loose L1)")
    p("| L1 form | group (anchor) | genes | internal L1 edges | E1 edge | shared-exon bp > 0 | gene-body only | exon proxy only | exon proxy only with 0 shared exonic bp | top-degree nodes: name (degree, body bp, exon union = body?) |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    L1FORMS = (("primary", "d_c2_approx", "d_c2_genebody", "d_c2_exon"), ("17:12 loose", "d_c2loose_approx", "d_c2loose_genebody", "d_c2loose_exon"))
    for tag, col, gcol, ecol in L1FORMS:
        E = [r for r in rows if nsl(r) and r[col] == "yes"]
        lab = lab_of([(r["gene_a"], r["gene_b"]) for r in E])
        for anc in ("NPIPB2", "TBC1D3"):
            S = grp_of(lab, anc)
            es = [r for r in E if r["gene_a"] in S and r["gene_b"] in S]
            deg = collections.Counter()
            for r in es:
                deg[r["gene_a"]] += 1
                deg[r["gene_b"]] += 1
            top = sorted(S, key=lambda g: (-deg[g], NAME[g]))[:8]
            gbo = sum(1 for r in es if r[gcol] == "yes" and r[ecol] != "yes")
            exo_ = [r for r in es if r[ecol] == "yes" and r[gcol] != "yes"]
            p(f"| {tag} | {anc} | {len(S)} | {len(es)} | {sum(r['d_e1_edge'] == 'yes' for r in es)} | "
              f"{sum((fnum(r['d_shared_exon_bp']) or 0) > 0 for r in es)} | {gbo} | {len(exo_)} | {sum(1 for r in exo_ if (fnum(r['d_shared_exon_bp']) or 0) == 0)} | "
              + ", ".join(f"{NAME[g]} ({deg[g]}, {body.get(g, 'NA')}, {'yes' if span_exon.get(g) else 'no'})" for g in top) + " |")
    p()
    p("Hub attribution: degree of the 17:12 loose-L1 top-degree nodes inside their loose L1 group, under each clause-2 form "
      "(edges restricted to that group's genes):")
    p("| hub | group | loose (17:12) | v-exon overlap on exon proxy only | v-exon overlap, no strand check | primary | loose edges passing only by the exon proxy | of those with 0 shared exonic bp | loose edges passing by the gene-body chain |")
    p("|---|---|---|---|---|---|---|---|---|")
    El = [r for r in rows if nsl(r) and r["d_c2loose_approx"] == "yes"]
    labl = lab_of([(r["gene_a"], r["gene_b"]) for r in El])
    for anc in ("NPIPB2", "TBC1D3"):
        S = grp_of(labl, anc)
        es = [r for r in El if r["gene_a"] in S and r["gene_b"] in S]
        deg = collections.Counter()
        for r in es:
            deg[r["gene_a"]] += 1
            deg[r["gene_b"]] += 1
        hubs = sorted(S, key=lambda g: (-deg[g], NAME[g]))[:6] + [BYN[n] for n in ("VN1R91P", "BNIP3P16", "PDXDC1", "LGALS9B")
                                                                   if n in BYN and BYN[n] in S and BYN[n] not in sorted(S, key=lambda g: (-deg[g], NAME[g]))[:6]]
        for h in hubs:
            mine = [r for r in rows if nsl(r) and h in (r["gene_a"], r["gene_b"]) and r["gene_a"] in S and r["gene_b"] in S]
            cnts = [sum(1 for r in mine if r[c] == "yes") for c in ("d_c2loose_approx", "d_c2exontgt_approx", "d_c2nostrand_approx", "d_c2_approx")]
            exonly = [r for r in mine if r["d_c2loose_exon"] == "yes" and r["d_c2loose_genebody"] != "yes"]
            p(f"| {NAME[h]} ({body.get(h, 'NA')} bp) | {anc} | " + " | ".join(map(str, cnts))
              + f" | {len(exonly)} | {sum(1 for r in exonly if (fnum(r['d_shared_exon_bp']) or 0) == 0)} | {sum(1 for r in mine if r['d_c2loose_genebody'] == 'yes')} |")
    p()

    # ------------------------------------------------------------------ component size profile, primary vs triangle
    grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
    p("## T-sizes: component-size profile per level (all of V), primary vs 3-truss, and the member-holding groups (from groups.tsv)")
    p("| level | primary: components >= 2 | primary: 2-gene components | primary: singletons | 3-truss: components >= 2 | 3-truss: 2-gene | 3-truss: singletons | NPIP group primary -> 3-truss (% removed) | TBC1D3 group primary -> 3-truss (% removed) | 2-gene primary components holding a member |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for base, tri, tag in (("primary", "triangle", ""), ("17:03_tests_exact", "17:03_tests_exact_triangle", " (17:12 tests)")):
        for k in LEVELS:
            out = []
            for v in (base, tri):
                c = collections.Counter(r[f"{v}|{k}"] for r in grp.values())
                sizes = collections.Counter(c.values())
                out.append((sum(n for sz, n in sizes.items() if sz >= 2), sizes.get(2, 0), sizes.get(1, 0)))
            c = collections.Counter(r[f"{base}|{k}"] for r in grp.values())
            two_mem = sorted({grp[g][f"{base}|{k}"] for g in MEM if c[grp[g][f"{base}|{k}"]] == 2})
            two_mem_names = ["+".join(sorted(NAME[g] for g in grp if grp[g][f"{base}|{k}"] == lab)) for lab in two_mem]
            cells = []
            for anc in ("NPIPB2", "TBC1D3"):
                a = sum(1 for g in grp if grp[g][f"{base}|{k}"] == grp[BYN[anc]][f"{base}|{k}"])
                b = sum(1 for g in grp if grp[g][f"{tri}|{k}"] == grp[BYN[anc]][f"{tri}|{k}"])
                cells.append(f"{a} -> {b} ({100 * (a - b) / a:.0f}%)")
            p(f"| {k}{tag} | {out[0][0]} | {out[0][1]} | {out[0][2]} | {out[1][0]} | {out[1][1]} | {out[1][2]} | {cells[0]} | {cells[1]} | {', '.join(two_mem_names) or '-'} |")
    p()

    # ------------------------------------------------------------------ L2 paths to co-duplicated neighbours + span-exon diagnostic
    adj2 = collections.defaultdict(dict)
    for r, t in zip(rows, T):
        if t[2]:
            adj2[r["gene_a"]][r["gene_b"]] = r
            adj2[r["gene_b"]][r["gene_a"]] = r
    p("## T-L2paths: shortest primary L2 paths from NPIPB2 / TBC1D3 (edge: shared-exon fraction; w98 gap-excl; exon union = body flags)")
    for src, tgts in (("NPIPB2", ("ZNF429", "SMG1", "BOLA2", "PLA2G10CP", "PDXDC1", "PKD1", "PKD1P1")), ("TBC1D3", ("USP32", "USP6", "CA4", "DHX40"))):
        for tgt in tgts:
            pp = shortest_path(adj2, BYN[src], BYN[tgt]) if tgt in BYN else None
            if not pp:
                p(f"- {src} -> {tgt}: not in the {src} L2 group")
                continue
            parts = []
            for i in range(len(pp) - 1):
                r = adj2[pp[i]][pp[i + 1]]
                a_, b_ = (pp[i], pp[i + 1]) if r["gene_a"] == pp[i] else (pp[i + 1], pp[i])
                parts.append(f"{NAME[pp[i]]}-{NAME[pp[i + 1]]} (sef {f3(r['d_shared_exon_frac'])}, w98 {f3(r['d_w98_gapexcl'])}, "
                             f"span-exon {NAME[a_]} {'yes' if span_exon.get(a_) else 'no'} / {NAME[b_]} {'yes' if span_exon.get(b_) else 'no'})")
            p(f"- {src} -> {tgt}: " + " > ".join(parts))
    p()
    p("Diagnostic (not a level of the lattice): L2 with the shared-exon clause made unsatisfiable on any edge touching a pseudogene "
      "record whose exon union equals its whole gene body (the exon-less-record convention), primary L1 otherwise.")
    E2d = []
    for r, t in zip(rows, T):
        if t[2]:
            bad = any(span_exon.get(g) and "pseudogene" in nodes[g]["biotype"] for g in (r["gene_a"], r["gene_b"]))
            if not bad:
                E2d.append((r["gene_a"], r["gene_b"]))
    lab2d = components(V, E2d)
    for anc in ("NPIPB2", "TBC1D3"):
        S = {g for g in nodes if lab2d[g] == lab2d[BYN[anc]]}
        znf = sum(1 for g in S if NAME[g].startswith("ZNF"))
        p(f"- {anc}: group {len(S)} genes, members {len(S & MEM)}, ZNF* {znf}, PKD1* {sum(1 for g in S if NAME[g].startswith('PKD1'))}, "
          f"non-members: {', '.join(sorted(NAME[g] for g in S - MEM)[:80])}")
    p()

    # ------------------------------------------------------------------ L3 non-member routes and named edges
    E3 = [r for r, t in zip(rows, T) if t[3]]
    lab3 = lab_of([(r["gene_a"], r["gene_b"]) for r in E3])
    p("## T-L3routes: non-members of the primary L3 member-holding groups and their L3 edges")
    for anc in ("NPIPB2", "TBC1D3"):
        S = grp_of(lab3, anc)
        for g in sorted(S - MEM, key=lambda g: NAME[g]):
            es = [r for r in E3 if g in (r["gene_a"], r["gene_b"])]
            p(f"- {anc} group: {NAME[g]} ({nodes[g]['biotype']}): " + "; ".join(
                f"{r['name_b'] if r['gene_a'] == g else r['name_a']} (sef {f3(r['d_shared_exon_frac'])}, w98 {f3(r['d_w98_gapexcl'], 4)})"
                for r in es[:8]))
    p()
    p("## T-named-edges: attributes of edges quoted in the report")
    for a_, b_ in (("USP31", "USP6"), ("NPIPB2", "NPIPB9"), ("NPIPB9", "BNIP3P16"), ("NPIPB2", "LOC131696449"), ("LOC131696449", "PKD1"),
                   ("TBC1D3", "TBC1D3P1-DHX40P1"), ("TBC1D3", "LGALS9B"), ("TBC1D3D", "LOC105371848"), ("TBC1D3D", "LOC105371853"),
                   ("LOC105371848", "LOC105371853"), ("NPIPA2", "NPIPB2"), ("LOC112268174", "NPIPB9"), ("NPIPA8", "PKD1P1"),
                   ("LOC100190986", "NPIPB5"), ("TBC1D3P3", "TBC1D3P4")):
        r = EDG.get(frozenset((a_, b_)))
        if r is None:
            p(f"- {a_}–{b_}: no edge row")
            continue
        i = rows.index(r)
        p(f"- {a_}–{b_}: primary tests {T[i]}; 17:12 tests {T_OLD[i]}; same-locus {r['same_locus']}; {edge_line(r)}; loose gene-body "
          f"frac {f3(r['d_c2loose_gb_best_frac'])} chain id {f3(r['d_c2loose_gb_chain_identity'], 4)}; loose exon frac "
          f"{f3(r['d_c2loose_exon_best_frac'])}; E1 pooled gap-incl {f3(r['d_e1_identity'], 4)}; w98 gap-incl {f3(r['d_w98_gapincl'], 4)}; "
          f"both spliced {r['d_both_spliced']}")
    p()
    p("## T-L3variants: member partition of each L3 identity form (member_groups.tsv), and w_98 vs pooled decisions")
    for v in ("primary", "L3=w98_gap-inclusive", "L3=pooled_gap-excluded", "L3=pooled_gap-inclusive", "L3=S2_SD98_mapback",
              "L1=c2_no_strand_check", "17:03_tests_exact"):
        rs = [r for r in mg if r["variant"] == v and r["level"] == "L3"]
        multi = [f"{{{r['members']}}}+{r['n_nonmembers']}" for r in rs if int(r["size"]) > 1]
        single = sorted(r["members"] for r in rs if r["size"] == "1")
        p(f"- {v} (L3 edges {lv[(v, 'L3')]['edges']}): " + " · ".join(multi) + f"; member singletons: {', '.join(single)}")
    for tag, l1 in (("primary L1", "c2"), ("17:12 loose L1", "c2_loose")):
        c = collections.Counter()
        for r in rows:
            t2 = tests(r, l1=l1)[2]
            if not t2:
                continue
            w = fnum(r["d_w98_gapexcl"])
            q = fnum(r["d_e1_identity_gapexcl"])
            c[(w is not None and w >= 0.98, q is not None and q >= 0.98)] += 1
        p(f"- {tag}: t_2 edges {sum(c.values())}; w_98 gap-excl >= 0.98 and pooled gap-excl >= 0.98: {c[(True, True)]}; w_98 only: "
          f"{c[(True, False)]}; pooled only: {c[(False, True)]}; neither: {c[(False, False)]}")
    p()
    p("## T-soto-members: Soto family (flag) and literature labels of every member")
    for g in sorted(MEM, key=lambda g: (nodes[g]["member_family"], NAME[g])):
        p(f"- {NAME[g]} ({nodes[g]['member_family']}): Soto {nodes[g]['soto_families'] or '-'} ({nodes[g]['soto_flag']}); lit L1 "
          f"{nodes[g]['lit_level1'] or '-'}; lit L2 {nodes[g]['lit_level2'] or '-'}; C_tree_top {nodes[g]['Ctree_top'] or '-'}; "
          f"reads any/unique {reads[g][0]}/{reads[g][1]}")
    p()
    p("## T-crosstrio: catalogs of the members; DNA evidence across chromosome trios")
    p(f"- members by (family, catalog, chromosome): {dict(collections.Counter((nodes[g]['member_family'], nodes[g]['catalog'], nodes[g]['chrom']) for g in sorted(MEM)))}")
    p(f"- primary L1 edges with no PAF catalog (cross-trio): {sum(1 for r, t in zip(rows, T) if t[1] and r['d_catalog'] == 'NA')}; "
      f"PAF pairs whose genes lie in different chromosome trios: {sum(1 for r in rows if r['d_evidence'] == 'yes' and nodes[r['gene_a']]['catalog'] != nodes[r['gene_b']]['catalog'])}")
    p()
    with open(f"{OUT}/report_tables.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {OUT}/report_tables.md ({len(L)} lines)")


# ==================================================================================================================== CLI
class _Tee:
    """stdout copy for `all`: expr-recount and corrected-tables write no log of their own (the reproduce block
    redirected their stdout to IS/expr_recount.out and IS/corrected_tables.out)."""

    def __init__(self, path):
        self.fh = open(path, "w")
        self.term = sys.stdout

    def write(self, s):
        self.term.write(s)
        self.fh.write(s)

    def flush(self):
        self.term.flush()
        self.fh.flush()

    def close(self):
        self.fh.close()


def cmd_all(args):
    """Both reproduce blocks in dependency order (LAYER_ORDER §11 without the archived variant/probe scripts, then
    NESTED_LATTICE §11). Stage names go to stderr; stdout is the stages' own output, in order."""
    ns = argparse.Namespace
    steps = [("expr-recount", cmd_expr_recount,
              ns(out=f"{LC.INT}/expr_recount.tsv", ignore=",".join(READTHROUGH_OVER_MEMBERS), extra=list(EXPR_RECOUNT_EXTRA)),
              f"{LC.INT}/expr_recount.out"),
             ("corrected-tables", cmd_corrected_tables, ns(), f"{LC.INT}/corrected_tables.out"),
             ("layer-order", cmd_layer_order, ns(), None),
             ("lattice-edges", cmd_lattice_edges, ns(), None),
             ("lattice-expr", cmd_lattice_expr, ns(), None),
             ("lattice-levels", cmd_lattice_levels, ns(), None),
             ("lattice-truth", cmd_lattice_truth, ns(), None),
             ("lattice-filtration", cmd_lattice_filtration, ns(l1="c2"), None),
             ("lattice-filtration --l1 c2_loose", cmd_lattice_filtration, ns(l1="c2_loose"), None)]
    if args.with_check_c2:
        steps.append(("lattice-check-c2", cmd_lattice_check_c2, ns(), None))
    steps.append(("lattice-report", cmd_lattice_report, ns(), None))
    for name, fn, a, log in steps:
        t0 = time.time()
        print(f"== {name} (root {LC.ROOT})", file=sys.stderr, flush=True)
        if log:
            tee = _Tee(log)
            sys.stdout = tee
            try:
                fn(a)
            finally:
                sys.stdout = tee.term
                tee.close()
        else:
            fn(a)
        print(f"== {name} done, {time.time() - t0:.0f} s", file=sys.stderr, flush=True)


def _pin_hash_seed():
    """Re-execute with PYTHONHASHSEED=HASH_SEED unless already set so (see Determinism in the module docstring)."""
    if os.environ.get("PYTHONHASHSEED") != HASH_SEED:
        os.environ["PYTHONHASHSEED"] = HASH_SEED
        os.execv(sys.executable, [sys.executable] + sys.argv)


def main(argv=None):
    p = argparse.ArgumentParser(prog="npip_tbc1d3.py", description=__doc__.split("\n\n")[0],
                                epilog="Old -> new commands, paths and determinism: see the module docstring "
                                       "(python3 -c 'import npip_tbc1d3; help(npip_tbc1d3)').",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    root_help = ("directory holding light/ heavy/ integrate_slim/ lattice/ (default $LO_ROOT or "
                 f"{LC.DEFAULT_ROOT}, the frozen results; stages OVERWRITE their outputs there)")
    p.add_argument("--root", help=root_help)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--root", default=argparse.SUPPRESS, help=root_help)
    sub = p.add_subparsers(dest="cmd", required=True, metavar="STAGE")

    def add(name, fn, help_):
        s = sub.add_parser(name, parents=[common], help=help_, description=help_)
        s.set_defaults(func=fn)
        return s

    s = add("expr-recount", cmd_expr_recount,
            "testis read counts for heavy/EXPR.counts.tsv genes + EXTRA ids (was lo_expr_recount.py)")
    s.add_argument("out", metavar="OUT_TSV", help="output table (the reproduce block uses IS/expr_recount.tsv); "
                                                  "OUT_TSV.windows.bed is written next to it")
    s.add_argument("--ignore", metavar="ID,ID,...",
                   help="also count n_reads_unique_mr, ignoring these records (the reproduce block: the 6 readthroughs "
                        f"{','.join(READTHROUGH_OVER_MEMBERS)})")
    s.add_argument("extra", nargs="*", metavar="EXTRA_ID",
                   help=f"gene names EXPR.counts.tsv lacks (the reproduce block: {' '.join(EXPR_RECOUNT_EXTRA)})")
    add("corrected-tables", cmd_corrected_tables,
        "light/*.corrected.tsv + heavy/EXPR.counts.corrected.tsv (was lo_corrected_tables.py)")
    add("layer-order", cmd_layer_order,
        "containment, tournaments, JOIN/REFINE, truths A/B, EXPR, disagreements -> IS/ (was lo_analysis.py)")
    add("lattice-edges", cmd_lattice_edges,
        "LAT/{edges,nodes,closure,edges_all_c2}.tsv, edges_build.out (was lattice_edges.py; 64 s, 2.5 GB)")
    add("lattice-expr", cmd_lattice_expr, "LAT/expr_counts.tsv, expr_counts.out (was lattice_expr.py)")
    add("lattice-levels", cmd_lattice_levels,
        "LAT/{levels,member_groups,groups,sanity,expr_views,chaining,triangle_drops}.tsv, levels.out "
        "(was lattice_levels.py)")
    add("lattice-truth", cmd_lattice_truth, "LAT/truth.tsv, truth_ingroup.tsv, truth.out (was lattice_truth.py)")
    s = add("lattice-filtration", cmd_lattice_filtration,
            "LAT/filtration.txt, filtration_{groups,appearance}.tsv (was lattice_filtration.py)")
    s.add_argument("--l1", choices=("c2", "c2_loose"), default="c2",
                   help="c2 = primary L1; c2_loose = the 17:12 report's L1, outputs suffixed .17_03_tests_exact")
    add("lattice-check-c2", cmd_lattice_check_c2,
        "LAT/check_c2.out: shipped gene_body_chains + v-exon/strand loop vs edges.tsv (was lattice_check_c2.py)")
    add("lattice-report", cmd_lattice_report, "LAT/report_tables.md (was lattice_report_tables.py)")
    s = add("all", cmd_all, "every stage above in dependency order (both reproduce blocks)")
    s.add_argument("--with-check-c2", action="store_true", help="also run lattice-check-c2 (30 s)")
    # `expr-recount OUT --ignore IDS EXTRA...` (the documented order) mixes a positional, an option and a
    # nargs='*' positional, which argparse only accepts intermixed on Python >= 3.13; move `--ignore X` to just
    # after the stage name so every Python 3 parses it the same way (the old lo_expr_recount.py parsed by hand).
    argv = list(sys.argv[1:] if argv is None else argv)
    if "expr-recount" in argv:
        k = argv.index("expr-recount")
        for j in range(k + 1, len(argv)):
            if argv[j] == "--ignore" and j + 1 < len(argv):
                pair = argv[j:j + 2]
                del argv[j:j + 2]
                argv[k + 1:k + 1] = pair
                break
            if argv[j].startswith("--ignore="):
                argv.insert(k + 1, argv.pop(j))
                break
    a = p.parse_args(argv)
    if getattr(a, "root", None):
        LC.set_root(os.path.abspath(a.root))
    a.func(a)


if __name__ == "__main__":
    _pin_hash_seed()
    main()
