#!/usr/bin/env python3
"""Nested edge-test lattice, step 6: every table quoted in bench/NESTED_LATTICE_NPIP_TBC1D3.md, regenerated from the
lattice/ result files (no new computation except shortest paths, simple counts and component labels over edges.tsv).
All numbers are formatted here from UNROUNDED values (edges.tsv and truth.tsv store full precision).
Output: lattice/report_tables.md
"""
import collections
import csv
import sys

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
from lattice_common import LEVELS, OUT, components, fnum, tests, tsv  # noqa: E402

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


def bfs_path(src, dst, adjd):
    prev = {src: None}
    q = collections.deque([src])
    while q:
        x = q.popleft()
        if x == dst:
            break
        for y in sorted(adjd[x]):
            if y not in prev:
                prev[y] = x
                q.append(y)
    if dst not in prev:
        return None
    out = [dst]
    while prev[out[-1]] is not None:
        out.append(prev[out[-1]])
    return out[::-1]


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
    pp = bfs_path(BYN["NPIPB2"], BYN["TBC1D3"], adj)
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
        pp = bfs_path(BYN[src], BYN[tgt], adj2) if tgt in BYN else None
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
