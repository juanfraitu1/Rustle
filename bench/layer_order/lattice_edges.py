#!/usr/bin/env python3
"""Nested edge-test lattice, step 1: the unified edge table and the L0 closure.

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
import collections
import re
import sys
import time

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/scripts")
from lattice_common import (CAT_CHROMS, DUMP, DUMP_S1, HEAVY, INT, LIGHT, OUT, PAF, tsv, write)  # noqa: E402

T0 = time.time()
LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def key_of(s):
    c, r = s.rsplit(":", 1)
    a, b = r.split("-")
    return (c, int(a), int(b))


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def exonic_bases_in(blocks, gene_start1, s, e):
    off = gene_start1 - 1
    acc = 0
    for bs, be in blocks:
        lo, hi = max(bs - off, s), min(be - off, e)
        if hi > lo:
            acc += hi - lo
    return acc


# ------------------------------------------------------------------------------------------------ genes
genes = {r["gene_id"]: r for r in tsv(f"{LIGHT}/work/refseq/genes.tsv")}
by_coord = collections.defaultdict(list)
for g in genes.values():
    by_coord[(g["chrom"], int(g["start0"]) + 1, int(g["end"]))].append(g)
exons_tsv = {r["gene_id"]: [tuple(map(int, x.split("-"))) for x in r["exons"].split(",")]
             for r in tsv(f"{LIGHT}/work/refseq/exons.tsv")}
U = {r["gene_id"]: r for r in tsv(f"{LIGHT}/universe.corrected.tsv")}
say(f"[genes] RefSeq gene/pseudogene records {len(genes)}; with exon rows {len(exons_tsv)}; U {len(U)}")

# ------------------------------------------------------------------------------------------------ catalogs: keys, exons
import lo_corrected_tables as CT  # noqa: E402

key2genes, key_blocks, key_cat = {}, {}, {}
for cat, c in CT.CATALOGS.items():
    k2g = CT.catalog_keys(by_coord, c["kind"], c["path"])
    node_ex = {}
    if c["kind"] == "nodes":
        names = {r["idx"]: r["name"] for r in tsv(c["path"] + ".names.tsv")}
        for r in tsv(c["path"]):
            k = (r["chrom"], int(r["start"]) + 1, int(r["end"]))
            bl = [tuple(map(int, x.split("-"))) for x in r["exons"].split(",")] if r["exons"] else []
            node_ex.setdefault(k, []).extend(bl)
    for k, gs in k2g.items():
        key2genes[k] = [g["gene_id"] for g in gs]
        key_cat[k] = cat
        if c["kind"] == "nodes":
            bl = node_ex.get(k, [])
        else:
            bl = [b for g in gs for b in exons_tsv.get(g["gene_id"], [])]
        if not bl:
            bl = [(k[1] - 1, k[2])]
        key_blocks[k] = merge(bl)
key_exlen = {k: sum(e - s for s, e in b) for k, b in key_blocks.items()}
gene_key = {}
for k, gs in key2genes.items():
    for g in gs:
        gene_key[g] = k
say(f"[catalogs] keys {len(key2genes)} (c15_17_22 {sum(1 for k in key_cat if key_cat[k] == 'c15_17_22')}, c16_19_20 "
    f"{sum(1 for k in key_cat if key_cat[k] == 'c16_19_20')}); keys with 0 genes {sum(1 for v in key2genes.values() if not v)}; "
    f"keys with >1 gene {sum(1 for v in key2genes.values() if len(v) > 1)}; genes with a key {len(gene_key)}")


# ------------------------------------------------------------------------------------------------ DNA pair attributes
def chain_groups(recs, Lq):
    """The record groups (chains) of bench/guided_pipeline.gene_body_chains on one (query, target, strand) record list,
    in target order. recs: dicts qs qe ts te nm bl strand."""
    rs = sorted(recs, key=lambda r: (r["ts"], r["te"]))
    out, cur = [], []
    strand = rs[0]["strand"]
    for r in rs + [None]:
        ok = False
        if r is not None and cur:
            gap = r["ts"] - max(x["te"] for x in cur)
            span = r["te"] - min(x["ts"] for x in cur)
            order = r["qs"] >= cur[-1]["qs"] if strand == "+" else r["qs"] <= cur[-1]["qs"]
            ok = gap <= Lq and span <= 2 * Lq and order
        if r is None or (cur and not ok):
            out.append(cur)
            cur = []
        if r is not None:
            cur.append(r)
    return out


def chain_stats(cur, Lq, clen):
    """(passes with the finder's denominator, identity, aligned frac (finder), passes with the shorter-body denominator,
    aligned frac (shorter body), target span start, target span end) of one chain."""
    strand = cur[0]["strand"]
    nm, bl = sum(x["nm"] for x in cur), sum(x["bl"] for x in cur)
    qiv = merge([(x["qs"], x["qe"]) for x in cur])
    aligned = sum(e - s for s, e in qiv)
    ts, te = min(x["ts"] for x in cur), max(x["te"] for x in cur)
    if strand == "+":
        xs, xe = ts - qiv[0][0], te + (Lq - qiv[-1][1])
    else:
        xs, xe = ts - (Lq - qiv[-1][1]), te + qiv[0][0]
    xs, xe = max(0, xs), min(clen, xe)
    den = min(Lq, xe - xs)
    lit = min(Lq, clen)
    return (nm / bl >= 0.80 and aligned >= 0.50 * den, nm / bl, aligned / den if den > 0 else float("inf"),
            nm / bl >= 0.80 and aligned >= 0.50 * lit, aligned / lit, ts, te)


def chain_eval(recs, Lq, clen):
    """bench/guided_pipeline.gene_body_chains on one (query, target, strand) record list. recs: dicts qs qe ts te nm bl.
    Returns list of (passes (finder denominator), identity, aligned_frac (finder), passes (shorter body), aligned_frac
    (shorter body)) — the 17:03 interface, kept for the audit scripts."""
    return [chain_stats(c, Lq, clen)[:5] for c in chain_groups(recs, Lq)]


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
from protein_families import edges_from, pair_hsps  # noqa: E402

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
prev_V = {r["gene_id"] for r in tsv(f"{OUT}/pre_correction_1712/nodes.tsv")} if __import__("os").path.exists(
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
import lo_analysis as LA  # noqa: E402

ct_sup = {}
for r in LA.CT_ROWS:
    ct_sup[(r["family"], frozenset(r["cluster_smaller_side"].split(",")))] = r["support"]
ct_top = LA.CT_LAB["Ctree_top"]


def ctree_pair(x, y):
    nx_, ny_ = genes[x]["name"], genes[y]["name"]
    for fam, (L, K) in LA.CT_CLUSTERS.items():
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
all_hg = LA.hgnc_all(genes)
import soto_map  # noqa: E402

db = soto_map.load()
lit = {r["gene_id"]: r for r in tsv(f"{LIGHT}/truth_literature_subfamilies.corrected.tsv")}
nrows = []
for g in sorted(V):
    r = genes[g]
    ex = exons_tsv.get(g) or [(int(r["start0"]), int(r["end"]))]
    m = soto_map.map_gene(db, r["name"], r["chrom"], r["strand"], ex)
    flag = ("unmatched" if not m["soto_gene_id"] else "weak_match" if m["soto_match_quality"] == "weak"
            else "ambiguous_multi_family" if m["soto_ambiguous"] == "yes" else "ok")
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
                  "Ctree_top": LA.CT_LAB["Ctree_top"].get(g, ""), "Ctree_min": LA.CT_LAB["Ctree_min"].get(g, ""),
                  "description": r["description"]})
write(f"{OUT}/nodes.tsv", nrows)
say(f"[write] nodes.tsv rows {len(nrows)}; {time.time() - T0:.0f}s")
with open(f"{OUT}/edges_build.out", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
