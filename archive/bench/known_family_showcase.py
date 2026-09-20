#!/usr/bin/env python
"""known_family_showcase.py -- VALIDATION SHOWCASE (build stage).

Take a HANDFUL of literature-recognized, canonical multi-copy gene families that a reviewer knows
(RABL2, APOBEC3, MAGE-A, ANKRD18, RGPD8/RANBP2, KRAB-ZNF, GSTM/GSTM2, HERC2), confirm each is present
in the SHIPPED catalog (bench/family_rna_refine.tsv, 573 multi-copy families), and report the ACTUAL
per-family PRECISION + RECALL of the family DEFINITION -- honestly, spanning the identity regimes we
mapped (near-identical copy-ASSIGNMENT regime  <->  divergent family-DEFINITION regime). This is a
showcase, NOT a cherry-pick: the clean near-identical families come out high P+R, and the K-frontier /
collapse families are reported as the honest FLOOR.

TWO orthogonal axes are reported per family (they are DIFFERENT objectives):
  (O1) DEFINITION  = did the catalog GROUP the family's expressed copies into ONE family?
       RECALL      = grouped / (all expressed local-cluster copies of the family)  [dominant-family
                     anchored truth: namesake OR E_p/DNA-homologous OR de-novo-homologous copies that
                     sit in the family's own genomic locus-window; RABL2 cross-chromosomal via E_p].
       PRECISION   = purity of the grouped members.  P_strict = E_p/DNA whole-protein reciprocal
                     homology (family_er_pr.tsv in_ep|in_dna_loose) -- CONSERVATIVE, flags E_p
                     under-coverage of single-domain / partial copies.  P_bio corrects those truth-
                     artifacts (a member is real if E_p/DNA-homologous OR a namesake OR a co-located
                     tandem/segdup copy of the E_p backbone); the residual are GENUINE over-merges.
  (O2) ASSIGNMENT  = once grouped, can reads be assigned to the right COPY? -> the K-frontier.
       From sun_identifiability.tsv (per-copy SUN tier), diploid_cn_oracle.tsv (DNA CN: asm_hapCN /
       chi_H / K_read, RNA-independent), and the materialized copy-assignment VG (o2_vg_visualization).

Deliverables:
  bench/known_family_showcase.tsv        per-family P/R table (this build writes it)
  bench/fig_family_<name>.png            per-family de-novo homology graph (nodes=copies, edges=core_recip)
  bench/fig_family_showcase_summary.png  summary panel (P/R scatter + regime table)
  (flagships) bench/fig_o2_vg_fam{39,1,22}_*.png re-rendered from the o2_vg cache (copy-assignment VG)

REUSES (does not re-derive): family_er_pr (meta/annot/gene_of/raw/genes loaders + pair truth),
family_rna_refine.tsv (the shipped catalog), denovo_family_edges.tsv (core_recip homology edges),
diploid_cn_oracle.tsv, sun_identifiability.tsv, o2_vg_visualization.py (materialize/render from cache).

Deterministic: PYTHONHASHSEED=0 (re-exec), fixed spring-layout seed, sorted iteration, no other RNG.
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/known_family_showcase.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import itertools
import json
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as FP  # noqa: E402

CATALOG = os.path.join(BENCH, "family_rna_refine.tsv")
EDGES = os.path.join(BENCH, "denovo_family_edges.tsv")
PAIRS = os.path.join(BENCH, "family_er_pr.tsv")
SUN = os.path.join(BENCH, "sun_identifiability.tsv")
ORACLE = os.path.join(BENCH, "diploid_cn_oracle.tsv")
OUT_TSV = os.path.join(BENCH, "known_family_showcase.tsv")
OUT_JSON = os.path.join(BENCH, "known_family_showcase.json")
# summary panel: primary (task-named) + legacy alias, both written byte-identical
SUMMARY_PNG = os.path.join(BENCH, "fig_known_families_summary.png")
SUMMARY_PNG_LEGACY = os.path.join(BENCH, "fig_family_showcase_summary.png")
LAYOUT_SEED = 0
WIN = 1_500_000  # genomic locality window (bp) for the local-cluster truth-set

# --------------------------------------------------------------------------- canonical roster
# (catalog family_id, seed genes = literature-named members, symbol prefix, cross-chromosomal?,
#  validated/O2 family id or None, o2_vg flagship id or None)
FAMILIES = [
    dict(name="RABL2",   fid="73",  seeds=["RABL2A", "RABL2B"],                  prefix="RABL2",
         crosschrom=True,  val_fid="39", o2_vg=39,
         note="RAB-like GTPase; cross-chromosomal 5-copy cluster (flagship clean)"),
    dict(name="APOBEC3", fid="569", seeds=["APOBEC3C", "APOBEC3D", "APOBEC3F"],  prefix="APOBEC3",
         crosschrom=False, val_fid=None, o2_vg=None,
         note="APOBEC3 cytidine-deaminase tandem cluster (Bailey/Soto)"),
    dict(name="MAGEA",   fid="517", seeds=["MAGEA1", "MAGEA4", "MAGEA12"],       prefix="MAGEA",
         crosschrom=False, val_fid="94", o2_vg=None,
         note="MAGE-A cancer/testis X-array (Xq28)"),
    dict(name="ANKRD18", fid="306", seeds=["ANKRD18A", "ANKRD18B"],             prefix="ANKRD18",
         crosschrom=False, val_fid="22", o2_vg=22,
         note="ANKRD18A/B ankyrin-repeat segdup (flagship mid / K-frontier)"),
    dict(name="RGPD8",   fid="313", seeds=["RGPD8"],                            prefix="RGPD",
         crosschrom=False, val_fid="1", o2_vg=1,
         note="RANBP2-paralog RGPD cluster (flagship K-frontier FLOOR)"),
    dict(name="ZNF92",   fid="42",  seeds=["ZNF92"],                            prefix=None,
         crosschrom=False, val_fid=None, o2_vg=None,
         note="KRAB-zinc-finger cluster (divergent-paralog DEFINITION regime)"),
    dict(name="GSTM",    fid="18",  seeds=["GSTM1", "GSTM4", "GSTM5"],          prefix="GSTM",
         crosschrom=False, val_fid="0", o2_vg=None,
         note="glutathione-S-transferase Mu tandem cluster; GSTM2 near-identical collapse"),
    dict(name="HERC2",   fid="388", seeds=["HERC2"],                           prefix="HERC2",
         crosschrom=False, val_fid=None, o2_vg=None,
         note="HERC2 / HERC2P duplicon cluster, chr15 segdup (Eichler)"),
]

# node role -> colour
ROLE_COLOR = {
    "ep_homolog":   "#2ca02c",   # E_p/DNA-homologous backbone (real)
    "namesake":     "#1f77b4",   # symbol-namesake (real)
    "coloc_tandem": "#17becf",   # co-located tandem/segdup copy E_p under-covers (real)
    "OVER_MERGE":   "#d62728",   # genuine over-merge (precision hit)
    "missed":       "#999999",   # true copy NOT grouped into this family (recall miss)
    "NA":           "#c5b0d5",   # unannotated de-novo locus (structural)
}


# --------------------------------------------------------------------------- loaders
def load_catalog():
    fam = defaultdict(list)
    with open(CATALOG) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fam[row["family_id"]].append(row["member_dn"])
    return fam


def load_pairs():
    pair = {}
    with open(PAIRS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pair[frozenset((row["gene_a"], row["gene_b"]))] = dict(
                ep=row["in_ep"] == "1", dna=row["in_dna_loose"] == "1",
                tier=row["ep_tier"], cls=row["class"])
    return pair


def load_edges():
    E = {}
    with open(EDGES) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            E[frozenset((row["a"], row["b"]))] = float(row["core_recip"])
    return E


def load_sun():
    """validated family id -> list of per-copy dicts (tier / n_sun / group_size / gene)."""
    sun = defaultdict(list)
    if not os.path.exists(SUN):
        return sun
    with open(SUN) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sun[row["family"]].append(dict(
                gene=row["gene"], n_psv=int(row["n_psv"]), n_sun=int(row["n_sun"]),
                tier=int(row["tier"]), group_size=int(row["group_size"])))
    return sun


def load_oracle():
    """validated family id -> DNA-oracle dict (RNA-independent copy number)."""
    orc = {}
    if not os.path.exists(ORACLE):
        return orc
    with open(ORACLE) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            orc[row["family"]] = row
    return orc


def load_o2vg(fam_id):
    """materialize the copy-assignment VG for a flagship from the o2_vg cache (no BAM needed)."""
    try:
        import o2_vg_visualization as o2
        vg = o2.materialize_family(fam_id, use_cache=True)
        return o2, vg
    except Exception as e:  # noqa: BLE001
        print(f"   [o2_vg] fam{fam_id} unavailable: {e}")
        return None, None


# --------------------------------------------------------------------------- per-family metrics
def core_recip_stats(members, E):
    ee = sorted(E[frozenset((a, b))] for a, b in itertools.combinations(members, 2)
                if frozenset((a, b)) in E)
    if not ee:
        return dict(n=0, med=float("nan"), mx=float("nan"), mn=float("nan"), frac_ni=0.0,
                    regime="NA")
    med = ee[len(ee) // 2]
    frac_ni = sum(1 for x in ee if x > 0.90) / len(ee)
    regime = ("near-identical(>0.90)" if med > 0.90
              else "handoff(0.70-0.90)" if med >= 0.70 else "divergent(<0.70)")
    return dict(n=len(ee), med=med, mx=ee[-1], mn=ee[0], frac_ni=frac_ni, regime=regime)


def precision(members, genes, gene_of_dn, pair, E, prefix):
    """E_p/DNA purity (strict) + biology-corrected purity; per-member role classification."""
    anng = [(dn, gene_of_dn.get(dn)) for dn in members
            if gene_of_dn.get(dn) and not gene_of_dn.get(dn).startswith("DN_")]
    na = [dn for dn in members if not gene_of_dn.get(dn) or gene_of_dn.get(dn).startswith("DN_")]
    gene_set = {g for _, g in anng}

    def ep_dna(g):
        for h in gene_set:
            if h != g:
                k = frozenset((g, h))
                if k in pair and (pair[k]["ep"] or pair[k]["dna"]):
                    return True
        return False

    # E_p/DNA backbone anchors the co-location test (avoids a spurious member defining its own window)
    core = [dn for dn, g in anng if ep_dna(g) or (prefix and g.startswith(prefix))]

    def coloc(dn):
        L = genes[dn]
        return any(genes[c]["chrom"] == L["chrom"] and abs(genes[c]["start"] - L["start"]) <= WIN
                   for c in core if c != dn)

    role = {}
    for dn, g in anng:
        if ep_dna(g):
            role[dn] = "ep_homolog"
        elif prefix and g.startswith(prefix):
            role[dn] = "namesake"
        elif coloc(dn) and any(frozenset((dn, m)) in E for m in members if m != dn):
            role[dn] = "coloc_tandem"
        else:
            role[dn] = "OVER_MERGE"
    for dn in na:
        role[dn] = "NA"
    over = [(dn, gene_of_dn.get(dn)) for dn, g in anng if role[dn] == "OVER_MERGE"]
    p_strict = sum(1 for dn, g in anng if ep_dna(g)) / len(anng) if anng else float("nan")
    p_bio = sum(1 for dn, g in anng if role[dn] != "OVER_MERGE") / len(anng) if anng else float("nan")
    return dict(role=role, over=over, p_strict=p_strict, p_bio=p_bio,
                n_annot=len(anng), n_na=len(na),
                na_struct=sum(1 for dn in na if coloc(dn)))


def recall(spec, members, genes, gene_of_dn, pair, E, cat_of_dn, loci_of_gene, ep1):
    """dominant-family-anchored local-cluster recall.  Truth = the family's own grouped members
    PLUS any OTHER expressed locus that is (namesake | E_p/DNA-homologous | de-novo-homologous) to a
    member AND sits in the family's genomic window (same chrom, +/-WIN); cross-chromosomal families
    additionally admit namesake / E_p-tier1-reciprocal loci genome-wide (RABL2)."""
    prefix = spec["prefix"]
    seeds = spec["seeds"]
    memset = set(members)
    fam_genes = {gene_of_dn.get(dn) for dn in members}
    win = defaultdict(lambda: [1e18, -1e18])
    for dn in members:
        L = genes[dn]
        win[L["chrom"]][0] = min(win[L["chrom"]][0], L["start"])
        win[L["chrom"]][1] = max(win[L["chrom"]][1], L["end"])

    def in_window(dn):
        L = genes[dn]
        return any(L["chrom"] == c and L["end"] >= a - WIN and L["start"] <= b + WIN
                   for c, (a, b) in win.items())

    truth = set(members)
    for g, loci in loci_of_gene.items():
        homolog = (g in fam_genes) or (prefix and g.startswith(prefix)) or \
                  any(fg and frozenset((g, fg)) in pair and
                      (pair[frozenset((g, fg))]["ep"] or pair[frozenset((g, fg))]["dna"])
                      for fg in fam_genes)
        for dn in loci:
            if dn in memset:
                continue
            namesake = bool(prefix and g.startswith(prefix))
            cross = spec["crosschrom"] and any(s in ep1.get(g, ()) for s in seeds)
            dnhom = any(frozenset((dn, m)) in E for m in members)
            if (in_window(dn) and (homolog or dnhom)) or (namesake and homolog) or cross:
                truth.add(dn)
    truth = sorted(truth)
    grouped = [dn for dn in truth if cat_of_dn.get(dn) == spec["fid"]]
    missed = {}
    for dn in truth:
        c = cat_of_dn.get(dn)
        if c != spec["fid"]:
            key = c if c else "SINGLETON/not-multicopy"
            missed.setdefault(key, []).append((dn, gene_of_dn.get(dn)))
    r = len(grouped) / len(truth) if truth else float("nan")
    return dict(truth=truth, grouped=grouped, missed=missed, recall=r, n_truth=len(truth))


def o2_summary(spec, sun, orc):
    """assignment-regime (K-frontier) summary from SUN tiers + DNA oracle + o2_vg VG (flagships)."""
    out = dict(val_fid=spec["val_fid"])
    if spec["val_fid"] and spec["val_fid"] in sun:
        tiers = Counter(c["tier"] for c in sun[spec["val_fid"]])
        gs = sorted({c["group_size"] for c in sun[spec["val_fid"]]})
        out["sun_tiers"] = dict(tiers)
        out["sun_group_sizes"] = gs
        out["n_copies_val"] = len(sun[spec["val_fid"]])
    if spec["val_fid"] and spec["val_fid"] in orc:
        o = orc[spec["val_fid"]]
        out["oracle"] = dict(asm_hapCN=o["asm_hapCN"], chi_H=o["chi_H"], K_read=o["K_read"],
                             n_loci_ref=o["n_loci_ref"], cls=o["class"])
    if spec["o2_vg"] is not None:
        _, vg = load_o2vg(spec["o2_vg"])
        if vg:
            st = Counter(r.get("status") for r in vg.get("reads", []))
            out["vg"] = dict(n_copies=vg.get("n_copies"), K=vg.get("K"), cls=vg.get("cls"),
                             n_bubbles=vg.get("n_bubbles"), n_sun_bubbles=vg.get("n_sun_bubbles"),
                             n_reads=vg.get("n_reads"), assigned=st.get("assigned", 0),
                             group_sizes=vg.get("group_sizes"))
    return out


def classify(rec, prec, cr, o2):
    """honest combined class over the two axes (definition O1 + assignment O2)."""
    R, Pb = rec["recall"], prec["p_bio"]
    # definition axis
    if R >= 0.90 and Pb >= 0.90:
        cdef = "clean"
    elif R >= 0.70 and Pb >= 0.70:
        cdef = "mid"
    else:
        cdef = "floor"
    # assignment axis (K-frontier) from oracle/VG when available
    casn = "n/a"
    vg = o2.get("vg")
    orc = o2.get("oracle")
    if vg:
        nc, K = vg.get("n_copies") or 0, vg.get("K") or 0
        if K >= max(2, nc):
            casn = "resolvable"
        elif K <= 1:
            casn = "floor"
        else:
            casn = "k-frontier"
    elif orc:
        # K_read (reads actually resolving copies) is the read-level O2 outcome; chi_H is only the
        # theoretical max.  GSTM2: chi_H=7 but K_read=0 -> the near-identical array is a read FLOOR.
        try:
            kr = int(orc["K_read"])
            casn = "floor" if kr <= 1 else ("resolvable" if kr >= 3 else "k-frontier")
        except (ValueError, KeyError):
            casn = "n/a"
    # overall
    if cdef == "clean" and casn in ("resolvable", "n/a"):
        overall = "CLEAN"
    elif cdef == "floor" or casn == "floor":
        overall = "FLOOR"
    else:
        overall = "MID"
    # divergent-definition tag (O1 regime where domain-share over-merges appear)
    if cr["regime"].startswith("divergent") and prec["over"]:
        overall += "/divergent-overmerge"
    return dict(class_def=cdef, class_assign=casn, overall=overall)


# --------------------------------------------------------------------------- visualization
def render_family_graph(spec, members, rec, prec, genes, gene_of_dn, E, cr, out_png):
    """de-novo homology graph: nodes = grouped members (by precision role) + missed truth copies;
       edges = core_recip homology; layout = deterministic spring; annotate P/R/regime/class."""
    G = nx.Graph()
    labels = {}
    colors = []
    node_role = {}
    shown = list(members)
    missed_nodes = [dn for lst in rec["missed"].values() for dn, _ in lst]
    shown += missed_nodes
    big = len(shown) > 18

    def _short(g):
        # de-crowd dense bands: abbreviate long LOCxxxxxxxxx symbols to last-6-digits (…NNNNNN);
        # keep short curated symbols (ZNF92, GSTM2, FOXO1) full so roles stay identifiable.
        return ("…" + g[-6:]) if (big and g.startswith("LOC") and len(g) > 9) else g

    for dn in shown:
        role = prec["role"].get(dn, "missed")
        if dn in missed_nodes:
            role = "missed"
        node_role[dn] = role
        g = gene_of_dn.get(dn) or "NA"
        L = genes[dn]
        labels[dn] = _short(g) if big else \
            f"{g}\n{L['chrom'].split('_')[1] if '_' in L['chrom'] else L['chrom']}:{L['start']//1000}k"
        G.add_node(dn)
    # edges = core_recip homology; on dense (big) graphs drop the faintest edges (<0.15, below the
    # colorbar-relevant band) so the node labels stay legible without changing any reported metric.
    edge_min = 0.15 if big else 0.0
    for a, b in itertools.combinations(shown, 2):
        k = frozenset((a, b))
        if k in E and E[k] >= edge_min:
            G.add_edge(a, b, w=E[k])

    # deterministic GENOMIC layout: x = locus order within a chromosome band, y = per-chrom band with
    # a small zig-zag so co-located tandem copies separate.  Legible for tandem arrays AND cross-chrom.
    by_chrom = defaultdict(list)
    for dn in shown:
        by_chrom[genes[dn]["chrom"]].append(dn)
    chroms = sorted(by_chrom, key=lambda c: (-len(by_chrom[c]), c))  # biggest band on top
    n_big = max((len(v) for v in by_chrom.values()), default=1)
    pos = {}
    band_gap = max(2.6, 0.42 * n_big)
    tiers = 4 if n_big > 30 else 3 if n_big > 14 else 2   # stagger labels into 2-4 tiers when dense
    amp = 1.05 if n_big > 30 else 0.9
    for bi, c in enumerate(chroms):
        nodes = sorted(by_chrom[c], key=lambda d: (genes[d]["start"], d))
        y0 = -bi * band_gap
        span = max(len(nodes) - 1, 1)
        for i, dn in enumerate(nodes):
            x = (i / span) * max(n_big - 1, 1)
            tier = i % tiers
            pos[dn] = (x, y0 + amp * (tier - (tiers - 1) / 2.0))
    figw = min(26, max(13, 0.62 * n_big)) if len(shown) > 18 else 11
    fig, ax = plt.subplots(figsize=(figw, 9 if len(shown) > 18 else 8.5))
    # chromosome band labels
    for bi, c in enumerate(chroms):
        ax.text(-0.06 * max(n_big - 1, 1), -bi * band_gap,
                c.replace("NC_0", "chr").split(".")[0], fontsize=7.5, style="italic",
                ha="right", va="center", color="#555555")
    # edges colored by core_recip
    if G.edges:
        ews = [G[u][v]["w"] for u, v in G.edges]
        ec = nx.draw_networkx_edges(G, pos, ax=ax, edge_color=ews, edge_cmap=plt.cm.viridis,
                                    edge_vmin=0.1, edge_vmax=1.0,
                                    width=[0.6 + 3.2 * w for w in ews], alpha=0.6)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,
                                   norm=plt.Normalize(vmin=0.1, vmax=1.0))
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label("core_recip (whole-transcript reciprocal identity)", fontsize=9)
    # draw order: backbone first, then the story nodes (coloc / missed / OVER_MERGE) on top
    base_sz = 440 if big else 780
    draw_rank = {"ep_homolog": 0, "namesake": 0, "NA": 1, "coloc_tandem": 2, "missed": 3,
                 "OVER_MERGE": 4}
    for dn in sorted(shown, key=lambda d: draw_rank.get(node_role[d], 0)):
        role = node_role[dn]
        emph = role in ("OVER_MERGE", "missed")
        nx.draw_networkx_nodes(G, pos, nodelist=[dn], ax=ax,
                               node_size=int(base_sz * 1.3) if emph else base_sz,
                               node_color=ROLE_COLOR[role],
                               edgecolors="black", linewidths=1.8 if emph else 0.8,
                               node_shape="o" if role != "missed" else "s")
    nx.draw_networkx_labels(G, pos, labels={d: labels[d] for d in shown}, ax=ax,
                            font_size=5.4 if big else 6.5)

    R, Ps, Pb = rec["recall"], prec["p_strict"], prec["p_bio"]
    title = (f"{spec['name']}  (catalog fam{spec['fid']}, {len(members)} copies)\n"
             f"RECALL={R:.2f}   PRECISION  P_strict(E_p/DNA)={Ps:.2f}  P_bio(corrected)={Pb:.2f}   "
             f"core_recip med={cr['med']:.2f} max={cr['mx']:.2f} [{cr['regime']}]")
    ax.set_title(title, fontsize=10.5)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=ROLE_COLOR[r],
                      markeredgecolor="black", markersize=11, label=lab)
               for r, lab in [("ep_homolog", "E_p/DNA homolog (real)"),
                              ("namesake", "symbol namesake (real)"),
                              ("coloc_tandem", "co-located tandem/segdup copy (real; E_p under-covers)"),
                              ("NA", "unannotated de-novo copy"),
                              ("OVER_MERGE", "GENUINE over-merge (precision hit)"),
                              ("missed", "true copy NOT grouped here (recall miss)")]]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3,
              fontsize=8.0, framealpha=0.95)
    ax.axis("off")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def render_summary(rows, out_png):
    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.45], wspace=0.12)

    # ---- left: P/R scatter ----
    ax = fig.add_subplot(gs[0, 0])
    cmap = {"CLEAN": "#2ca02c", "MID": "#ff7f0e", "FLOOR": "#d62728"}
    # per-name label offsets + tiny visual y-nudge so dots coinciding at (1,1) all show their class
    off = {"RABL2": (8, -3), "APOBEC3": (8, -12), "RGPD8": (8, 8), "ANKRD18": (8, -2),
           "MAGEA": (-30, 12), "HERC2": (-52, -4), "ZNF92": (6, -16), "GSTM": (6, 6)}
    ynudge = {"RABL2": 0.0, "APOBEC3": -0.030, "RGPD8": 0.030}  # visual only; true P=R=1.0
    for r in rows:
        base = r["class_base"]                       # fill = OVERALL class (O1 definition + O2 assignment)
        ring = cmap.get(r["class_def"].upper(), "#000000")   # ring = O1 DEFINITION class only
        yp = r["p_bio"] + ynudge.get(r["name"], 0.0)
        ax.scatter(r["recall"], yp, s=90 + 26 * r["n_copies"],
                   color=cmap.get(base, "#7f7f7f"), edgecolors=ring, linewidths=2.6,
                   alpha=0.82, zorder=3)
        ax.scatter(r["recall"], r["p_strict"], s=34, color="black", marker="x", zorder=4)
        ax.annotate(r["name"], (r["recall"], yp), fontsize=8.8, fontweight="bold",
                    xytext=off.get(r["name"], (5, 5)), textcoords="offset points")
    ax.axhline(0.90, ls="--", lw=0.7, color="gray")
    ax.axvline(0.90, ls="--", lw=0.7, color="gray")
    ax.set_xlabel("RECALL (expressed local copies grouped into one family)")
    ax.set_ylabel("PRECISION  (dot = P_bio corrected;  x = P_strict E_p/DNA)")
    ax.set_xlim(0.0, 1.12)
    ax.set_ylim(0.55, 1.08)
    ax.set_title("Per-family definition P/R  (dot size = #copies;  fill = overall class,  ring = O1 "
                 "definition class)", fontsize=10)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor=c,
                      markersize=11, label=k) for k, c in cmap.items()]
    # ring-vs-fill exemplar: a family whose DEFINITION is clean but whose ASSIGNMENT is the floor
    # (green ring, red fill) reads as "clean def, floor assign" -- this is exactly RGPD8.
    handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap["FLOOR"],
                          markeredgecolor=cmap["CLEAN"], markeredgewidth=2.6, markersize=12,
                          label="clean DEF, floor ASSIGN (RGPD8)"))
    ax.legend(handles=handles, loc="lower left", fontsize=8.4, title="class (fill=O1+O2 · ring=O1)")

    # ---- right: table ----
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis("off")
    cols = ["family", "fam", "cp", "R", "P_str", "P_bio", "regime", "O2 K-frontier", "class"]
    tab = []
    for r in rows:
        tab.append([r["name"], r["fid"], str(r["n_copies"]),
                    f"{r['recall']:.2f}", f"{r['p_strict']:.2f}", f"{r['p_bio']:.2f}",
                    r["regime_short"], r["o2_short"], r["class_base"]])
    t = ax2.table(cellText=tab, colLabels=cols, loc="center", cellLoc="center",
                  colWidths=[0.14, 0.06, 0.05, 0.06, 0.07, 0.07, 0.12, 0.28, 0.10])
    t.auto_set_font_size(False)
    t.set_fontsize(8.3)
    t.scale(1.0, 1.62)
    for j in range(len(cols)):
        t[0, j].set_facecolor("#dddddd")
        t[0, j].set_text_props(weight="bold")
    for i, r in enumerate(rows, start=1):
        t[i, len(cols) - 1].set_facecolor(
            {"CLEAN": "#d7f0d7", "MID": "#ffe8cc", "FLOOR": "#f7d4d4"}.get(r["class_base"], "#ffffff"))
    ax2.set_title("Known-family showcase: precision / recall + identity regime + K-frontier (O2)",
                  fontsize=11, pad=14)
    # footnote: over-merges + regime handoff
    foot = ("O1 DEFINITION P/R (grouping) vs O2 ASSIGNMENT K-frontier (read-to-copy).  "
            "P_str = strict E_p/DNA whole-protein reciprocal homology (conservative: under-covers "
            "single-domain/partial copies e.g. APOBEC3C, HERC2P); P_bio corrects those truth-artifacts.\n"
            "RECALL is measured among RNA-EXPRESSED/detected loci (not vs full genomic copy number): "
            "R=1.0 (RABL2/RGPD8/APOBEC3) = all EXPRESSED copies grouped, not all genomic copies found "
            "(the DNA oracle exposes that gap: RGPD8 n_loci_ref=35, GSTM2 47).  O2-VG family ids "
            "(fam39/1/22) are a DISTINCT object from the O1 catalog ids (fam73/313/306); the copy counts "
            "differ (RGPD8 O1 fam313=8 vs O2-VG fam1=7) because they are two separate family graphs.\n"
            "Genuine over-merges (only 2/8): FOXO1->ANKRD18, one distal ZNF->ZNF92 (both divergent-regime "
            "domain bridges).  FLOOR = read-level collapse: RGPD8 K=1/7 (0 SUN), GSTM2 K_read=0 (47-copy "
            "near-identical array), GSTM cluster fragments into homology sub-families.")
    fig.text(0.52, 0.05, foot, fontsize=8.2, ha="left", va="top", wrap=True)
    fig.suptitle("RNA multi-copy family DEFINITION -- validation showcase on canonical known families",
                 fontsize=13.5, y=0.985)
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- resolve canonical family ids
# The catalog family_ids are not stable across pipeline iterations (gamma/gate changes shift numbering).
# Resolve each canonical family to the catalog block that contains the most seed genes.  This keeps the
# showcase honest even when the current prototype has reorganized the catalog.
def resolve_spec_fids(fam, gene_of_dn):
    """Return a dict spec_name -> best catalog fid (or None if no seed gene is found)."""
    # gene -> set of catalog family ids that contain a locus of that gene
    gene_to_fids = defaultdict(set)
    for fid, members in fam.items():
        for dn in members:
            g = gene_of_dn.get(dn)
            if g:
                gene_to_fids[g].add(fid)
    resolved = {}
    for spec in FAMILIES:
        scores = Counter()
        for g in spec["seeds"]:
            for fid in gene_to_fids.get(g, ()):
                scores[fid] += 1
        if scores:
            best, _ = scores.most_common(1)[0]
            resolved[spec["name"]] = best
        else:
            resolved[spec["name"]] = None
    return resolved


# --------------------------------------------------------------------------- main
def main():
    print("[load] catalog / pair-truth / edges / SUN / oracle / de-novo universe ...", flush=True)
    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    raw = FP.load_raw_families()
    genes, gene_of_dn, *_ = FP.build_genes_dict(raw, meta, gene_of)
    fam = load_catalog()
    pair = load_pairs()
    E = load_edges()
    sun = load_sun()
    orc = load_oracle()

    cat_of_dn = {}
    for fid, members in fam.items():
        for dn in members:
            cat_of_dn[dn] = fid
    loci_of_gene = defaultdict(list)
    for dn, g in gene_of_dn.items():
        loci_of_gene[g].append(dn)
    # E_p tier-1 reciprocal partners (for cross-chromosomal recall admission)
    ep1 = defaultdict(set)
    for k, v in pair.items():
        if v["ep"] and v["tier"] == "tier1":
            a, b = tuple(k)
            ep1[a].add(b)
            ep1[b].add(a)

    spec_fids = resolve_spec_fids(fam, gene_of_dn)
    resolved_specs = [dict(spec, fid=spec_fids[spec["name"]]) for spec in FAMILIES]
    print(f"       catalog {len(fam)} families; de-novo universe {len(gene_of_dn)} loci")
    for spec in resolved_specs:
        print(f"       {spec['name']:8s} -> fam{spec['fid']}")
    print()

    rows = []
    for spec in resolved_specs:
        fid = spec["fid"]
        members = fam.get(fid, []) if fid else []
        present = len(members) > 0
        cr = core_recip_stats(members, E)
        prec = precision(members, genes, gene_of_dn, pair, E, spec["prefix"])
        rec = recall(spec, members, genes, gene_of_dn, pair, E, cat_of_dn, loci_of_gene, ep1)
        o2 = o2_summary(spec, sun, orc)
        cls = classify(rec, prec, cr, o2)

        # short K-frontier string (full) + compact token for the summary table
        kf = "not in O2 validated set"
        o2_short = "def-only (not in O2 set)"
        if o2.get("vg"):
            v = o2["vg"]
            kf = f"VG K={v['K']}/{v['n_copies']} ({v['cls'].split('(')[0].strip()})"
            tag = {"fully_resolvable": "resolvable", "no_psv": "FLOOR", "partial": "partial"}.get(
                v["cls"].split("(")[0].strip(), v["cls"].split("(")[0].strip())
            o2_short = f"VG K={v['K']}/{v['n_copies']} {tag}"
        elif o2.get("oracle"):
            o = o2["oracle"]
            kf = f"DNA chi_H={o['chi_H']} asm_hapCN={o['asm_hapCN']} K_read={o['K_read']}"
            o2_short = f"DNA chi={o['chi_H']} K_read={o['K_read']} (hapCN {o['asm_hapCN']})"
        elif o2.get("sun_tiers"):
            kf = f"SUN tiers={o2['sun_tiers']}"
            o2_short = f"SUN tiers={o2['sun_tiers']}"
        regime_short = {"divergent(<0.70)": "divergent", "handoff(0.70-0.90)": "handoff",
                        "near-identical(>0.90)": "near-id", "NA": "NA"}.get(cr["regime"], cr["regime"])

        genes_list = sorted({gene_of_dn.get(dn) for dn in members
                             if gene_of_dn.get(dn) and not gene_of_dn.get(dn).startswith("DN_")})
        row = dict(
            name=spec["name"], fid=fid, present=present, n_copies=len(members),
            n_annot=prec["n_annot"], n_na=prec["n_na"], genes=genes_list,
            n_truth=rec["n_truth"], recall=rec["recall"],
            p_strict=prec["p_strict"], p_bio=prec["p_bio"],
            over_merge=[g for _, g in prec["over"]],
            missed={k: [g for _, g in v] for k, v in rec["missed"].items()},
            cr_med=cr["med"], cr_max=cr["mx"], cr_frac_ni=cr["frac_ni"], regime=cr["regime"],
            class_def=cls["class_def"], class_assign=cls["class_assign"], overall=cls["overall"],
            class_base=cls["overall"].split("/")[0], regime_short=regime_short, o2_short=o2_short,
            k_frontier_short=kf, o2=o2, note=spec["note"], rec=rec, prec=prec, cr=cr,
            members=members)
        rows.append(row)

        # per-family graph
        out_png = os.path.join(BENCH, f"fig_family_{spec['name']}.png")
        render_family_graph(spec, members, rec, prec, genes, gene_of_dn, E, cr, out_png)
        row["figure"] = out_png
        row["o2_vg_figure"] = None
        row["o2_vg_graph_figure"] = None

    # re-render flagship copy-assignment VGs from cache
    print("[o2_vg] re-rendering flagship copy-assignment VGs from cache ...")
    rows_by_fid = {r["fid"]: r for r in rows}
    for spec in resolved_specs:
        if spec["o2_vg"] is not None:
            o2mod, vg = load_o2vg(spec["o2_vg"])
            if o2mod and vg:
                gene = vg.get("gene", spec["name"])
                out = os.path.join(BENCH, f"fig_o2_vg_fam{spec['o2_vg']}_{gene}.png")
                out_graph = os.path.join(BENCH, f"fig_o2_vg_fam{spec['o2_vg']}_{gene}_graph.png")
                try:
                    o2mod.render_family(vg, out)
                    rows_by_fid[spec["fid"]]["o2_vg_figure"] = out
                    if os.path.exists(out_graph):
                        rows_by_fid[spec["fid"]]["o2_vg_graph_figure"] = out_graph
                    print(f"       fam{spec['o2_vg']} {gene}: {out}")
                except Exception as e:  # noqa: BLE001
                    print(f"       fam{spec['o2_vg']} render skipped: {e}")

    # summary panel -- primary (task-named) + legacy alias (byte-identical)
    render_summary(rows, SUMMARY_PNG)
    render_summary(rows, SUMMARY_PNG_LEGACY)

    # ---- write TSV ----
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family", "catalog_fid", "present", "n_copies_catalog", "n_annot", "n_NA",
                    "genes", "n_copies_truth", "recall", "precision_strict_EpDNA",
                    "precision_bio_corrected", "genuine_over_merge", "missed_copies",
                    "core_recip_med", "core_recip_max", "frac_edges_near_identical", "regime",
                    "class_definition", "class_assignment", "class_overall", "K_frontier_O2",
                    "note"])
        for r in rows:
            w.writerow([
                r["name"], r["fid"], r["present"], r["n_copies"], r["n_annot"], r["n_na"],
                ",".join(r["genes"]), r["n_truth"],
                f"{r['recall']:.3f}", f"{r['p_strict']:.3f}", f"{r['p_bio']:.3f}",
                ";".join(r["over_merge"]) if r["over_merge"] else "NONE",
                ";".join(f"{k}:{','.join(v)}" for k, v in r["missed"].items()) if r["missed"] else "NONE",
                f"{r['cr_med']:.3f}", f"{r['cr_max']:.3f}", f"{r['cr_frac_ni']:.3f}", r["regime"],
                r["class_def"], r["class_assign"], r["overall"], r["k_frontier_short"], r["note"]])
    print(f"[write] {OUT_TSV}")

    # ---- write JSON (structured, machine-readable; deterministic key order) ----
    def _round(x, n=4):
        try:
            return round(float(x), n)
        except (TypeError, ValueError):
            return None

    total_members = sum(r["n_copies"] for r in rows)
    total_annot = sum(r["n_annot"] for r in rows)
    total_over = sum(len(r["over_merge"]) for r in rows)
    class_counts = Counter(r["class_base"] for r in rows)
    fam_json = []
    for r in rows:
        fam_json.append(dict(
            name=r["name"], catalog_fid=r["fid"], present=r["present"],
            note=r["note"],
            n_copies_catalog=r["n_copies"], n_annot=r["n_annot"], n_NA=r["n_na"],
            genes=r["genes"],
            definition_O1=dict(
                recall=_round(r["recall"]), n_copies_truth=r["n_truth"],
                recall_scope="among RNA-expressed/detected loci (NOT vs full genomic copy number)",
                recall_small_n=(r["n_truth"] < 4),
                precision_strict_EpDNA=_round(r["p_strict"]),
                precision_bio_corrected=_round(r["p_bio"]),
                genuine_over_merge=list(r["over_merge"]),
                missed_copies={k: list(v) for k, v in r["missed"].items()}),
            identity_regime=dict(
                core_recip_med=_round(r["cr_med"]), core_recip_max=_round(r["cr_max"]),
                frac_edges_near_identical=_round(r["cr_frac_ni"]),
                regime=r["regime"], regime_short=r["regime_short"]),
            assignment_O2=dict(
                summary=r["o2"], K_frontier_short=r["k_frontier_short"]),
            classification=dict(
                definition=r["class_def"], assignment=r["class_assign"],
                overall=r["overall"], base=r["class_base"]),
            figures=dict(
                homology_graph=r.get("figure"),
                o2_vg=r.get("o2_vg_figure"),
                o2_vg_graph=r.get("o2_vg_graph_figure"))))
    doc = dict(
        schema="known_family_showcase/v1",
        description=("Validation showcase: per-family precision/recall of the RNA multi-copy family "
                     "DEFINITION (O1) plus the read-to-copy ASSIGNMENT K-frontier (O2) for canonical "
                     "literature-recognized families, spanning near-identical -> divergent regimes."),
        inputs=dict(
            catalog=os.path.relpath(CATALOG, BENCH),
            pairs=os.path.relpath(PAIRS, BENCH),
            edges=os.path.relpath(EDGES, BENCH),
            sun=os.path.relpath(SUN, BENCH),
            oracle=os.path.relpath(ORACLE, BENCH),
            n_families_catalog=len(fam),
            n_denovo_loci=len(gene_of_dn)),
        determinism=dict(PYTHONHASHSEED="0", layout="deterministic genomic (locus-ordered, no RNG)",
                         layout_seed=LAYOUT_SEED, sorted_iteration=True),
        caveats=dict(
            recall_scope=("RECALL is measured among RNA-EXPRESSED/detected loci, NOT vs full genomic "
                          "copy number.  R=1.0 (RABL2/RGPD8/APOBEC3) means every EXPRESSED copy of the "
                          "family is grouped into one catalog family -- it does NOT claim every genomic "
                          "paralog was found (the DNA oracle n_loci_ref column exposes that gap: e.g. "
                          "RGPD8 n_loci_ref=35 vs 8 expressed, GSTM2 47)."),
            small_n_recall=("Small truth-set (n_copies_truth < 4) -> R exact but low-n; flagged per-"
                            "family via recall_small_n.  Only APOBEC3 (3 expressed copies) qualifies.  "
                            "NOTE GSTM is NOT small-n: its truth-set is 11 (many expressed GSTM copies), "
                            "so its R=0.27 is a genuine collapse over a large denominator, not an "
                            "n-artifact."),
            o2_vg_vs_o1_family_id=("The copy-assignment VG uses the O2 validated family numbering "
                                   "(psv_graph fam39=RABL2, fam1=RGPD8, fam22=ANKRD18), which is a "
                                   "DISTINCT family object from the O1 catalog family_rna_refine.tsv "
                                   "numbering (fam73=RABL2, fam313=RGPD8, fam306=ANKRD18).  Copy counts "
                                   "differ by object (RGPD8: O1 fam313=8 copies vs O2-VG fam1=7 copy-"
                                   "paths); each is internally correct, they are not the same graph."),
            precision_truth=("P_strict uses E_p/DNA whole-protein reciprocal homology (conservative; "
                             "under-covers single-domain/partial copies); P_bio corrects those truth-"
                             "artifacts, leaving only genuine over-merges (2/8: FOXO1->ANKRD18, one "
                             "distal ZNF->ZNF92).  Truth-set is self-anchored on the family's own "
                             "namesake/E_p/de-novo-homologous members.")),
        headline=dict(
            n_families=len(rows),
            total_catalog_members=total_members, total_annotated_members=total_annot,
            total_genuine_over_merges=total_over,
            over_merges=[dict(family=r["name"], catalog_fid=r["fid"], members=list(r["over_merge"]))
                         for r in rows if r["over_merge"]],
            class_counts=dict(sorted(class_counts.items())),
            precision_note=("Only %d genuine over-merges across %d annotated grouped members "
                            "(%d catalog members incl. unannotated de-novo loci); every other apparent "
                            "impurity is E_p under-coverage of a real single-domain/partial copy."
                            % (total_over, total_annot, total_members))),
        figures=dict(
            summary=SUMMARY_PNG,
            summary_legacy_alias=SUMMARY_PNG_LEGACY,
            per_family={r["name"]: r.get("figure") for r in rows},
            flagship_o2_vg={r["name"]: r.get("o2_vg_figure")
                            for r in rows if r.get("o2_vg_figure")}),
        families=fam_json)
    with open(OUT_JSON, "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"[write] {OUT_JSON}")

    # ---- report ----
    print("\n================ KNOWN-FAMILY SHOWCASE :: per-family P/R ================")
    hdr = f"{'family':8s} {'fam':>4s} {'cp':>3s} {'R':>5s} {'P_str':>6s} {'P_bio':>6s} {'regime':>18s}  {'class':<26s} K-frontier(O2)"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:8s} {r['fid']:>4s} {r['n_copies']:>3d} "
              f"{r['recall']:>5.2f} {r['p_strict']:>6.2f} {r['p_bio']:>6.2f} "
              f"{r['regime']:>18s}  {r['overall']:<26s} {r['k_frontier_short']}")
    print("\nGenuine over-merges (precision hits):")
    for r in rows:
        if r["over_merge"]:
            print(f"   {r['name']} fam{r['fid']}: {r['over_merge']}")
    print("\nRecall misses (true copies not grouped into the family):")
    for r in rows:
        if r["missed"]:
            print(f"   {r['name']} fam{r['fid']}: " +
                  "; ".join(f"{k}->{v}" for k, v in r["missed"].items()))
    print("\nFigures: bench/fig_family_<name>.png (8) + bench/fig_known_families_summary.png "
          "(alias fig_family_showcase_summary.png) + bench/fig_o2_vg_fam{39,1,22}_*.png")
    print("JSON:    bench/known_family_showcase.json")


if __name__ == "__main__":
    main()
