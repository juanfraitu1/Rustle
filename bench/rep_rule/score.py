#!/usr/bin/env python3
"""bench/rep_rule/score.py — the readouts and the registered decision of docs/PREREG_locus_representative_rule_2026-10-04.md:
the de novo locus representative R_M (most-reads, the shipped rule) vs R_J (most-junctions). Products of bench/rep_rule/run.sh.

    score.py h3-inputs --species S --contig C --out PREFIX   H3 copies: the contig's annotated protein-coding genes (human CAT/Liftoff
                                                             v2.0 slim GFF, gorilla RefSeq) -> PREFIX.copies.tsv (bench/copy_support.py's
                                                             columns, family ALL) + PREFIX.truth.gtf (their transcripts, gene_id = cid)
    score.py score --species S --contig C                    H1 + H2 + H3 + the two-gene representative count -> W/<S>_<C>/score.json
    score.py decide [--md PATH]                              per-contig tables + the decision of the prereg -> W/decision.json (+ markdown)

H1 reuses Figure 7's own scorers and references (no re-implementation): figures/_o1_recovery.py `family_score` (the family_score
binary of the run's bin dir), `score_arm` (its Python re-derivation, per-family rows) and `check_against_family_score`, against
Compara families at Primates restricted to the contig (the genome-wide table cached by `compara_families`, restricted by
`_o1.filter_rows` as `compara_contig_truth` does: the raw Compara download that function checks first is no longer on disk), Soto 2025,
the NPIP reference set (`npip_u2_truth`, chr16) and the contig's protein-homology families (the Figure 7 cache); figures/_liftoff.py
`support_if_ready` / `copy_pairs(loci, 0.95, support, both loci on the contig)` / `pair_families` / `catalog_loci` on each arm's copy
table (Figure 7's Liftoff rows, per contig). The family_score GFF is the Figure 7 GFF slice of the contig.
Environment (as bench/rep_rule/run.sh): REP_WORK (products, default /mnt/linuxdisk/tmp/rep_rule), REP_FIG7 (the read-only Figure 7
cache), REP_BIN (the binary dir of the runs; family_score runs from it), REP_CAT_GFF (the human H3 annotation).
H2: each arm's PREFIX.fam.copies.tsv: junctions per copy = gaps >= 50 bp between consecutive exons of its `exons` column (the rule's
own junction, the Task 1 review's F5), `n_exon - 1` (the prereg's wording) beside; loci whose representative differs between the arms
(loci.gff3 joined by gene Name); and a check that a Python re-derivation of both rules from the GTF reproduces each loci.gff3.
H3: bench/copy_support.py's JSON / per-copy table (strict FOUND = spliced-expressed and a same-strand locus whose representative
carries >= k of the gene's read-supported junctions, k = min(2, annotated introns)).
Two-gene representatives (the fused-locus caveat; reported beside, no rule): a representative whose exons overlap the exon unions of
>= 2 same-strand protein-coding genes of the GFF slice that are exon-disjoint from each other (RefSeq readthrough genes excluded).
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "figures"))

# Paths: the same environment variables as bench/rep_rule/run.sh (which passes its resolved values down), defaults = the run of
# 2026-10-04. REP_BIN is the ONE binary dir both arms ran with (its mcl_families / family_score sha1 are checked against the run
# logs, and family_score runs from it); the inputs file's `bin` is overridden by it.
W = Path(os.environ.get("REP_WORK", "/mnt/linuxdisk/tmp/rep_rule"))
FIG7 = Path(os.environ.get("REP_FIG7", "/mnt/linuxdisk/tmp/rustle_figures/fig7/current"))          # read-only
BIN = Path(os.environ.get("REP_BIN", "/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release"))
CAT_GFF = Path(os.environ.get(
    "REP_CAT_GFF", "/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz"))
SAMPLE = {"human": "human_A119b", "gorilla": "gorilla_OR6737"}          # the BAMs of figures/inputs.local.tsv human_bam / gorilla_bam
ARMS = ("R_M", "R_J")
ARM_RULE = {"R_M": "most-reads (shipped)", "R_J": "most-junctions"}
DEV = [("human", "chr16"), ("gorilla", "NC_073244.2")]
HELD_OUT = [("human", "chr2"), ("human", "chr6"), ("human", "chr8"), ("human", "chr10"), ("gorilla", "NC_073234.2")]
STATUS = {("human", "chr16"): "development", ("gorilla", "NC_073244.2"): "development",
          ("human", "chr2"): "held out (reused verdict set)", ("human", "chr6"): "held out (untouched)",
          ("human", "chr8"): "held out (reused verdict set)", ("human", "chr10"): "held out (reused verdict set)",
          ("gorilla", "NC_073234.2"): "held out (untouched)"}
PRIMARY = {"human": "compara", "gorilla": "liftoff"}
MIN_JUNCTION_GAP = 50          # the rule's junction (mcl_families MIN_JUNCTION_GAP = bench/copy_support.py MIN_INTRON)
F_TOL, P_TOL = 0.005, 0.01     # the prereg's clause (a)
LIFTOFF_SC = 0.95              # Figure 7's copy_pairs threshold


def _attrs(col9: str) -> dict:
    return dict(kv.split("=", 1) for kv in col9.strip().split(";") if "=" in kv)


def _merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _inter(a, b) -> int:
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def _sha1(p: Path) -> str:
    h = hashlib.sha1()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy_ro(src: Path, dst: Path) -> Path:
    """A read-only Figure 7 product copied into the scratch dir (re-copied when the source differs)."""
    if not dst.exists() or _sha1(dst) != _sha1(src):
        shutil.copyfile(src, dst)
    return dst


# ================================================================ H3 inputs
def h3_inputs(species: str, contig: str, out: str):
    """Protein-coding genes on `contig` -> copy_support.py's copies table + a truth GTF of their transcripts (gene_id = cid)."""
    import figlib
    cfg = figlib.load_inputs()
    if species == "human":
        src, opener = CAT_GFF, (lambda p: gzip.open(p, "rt"))
    else:
        src, opener = Path(cfg["gorilla_ref_gff"]), (lambda p: open(p))
    genes, parent_of, exons = {}, {}, collections.defaultdict(list)
    pref = contig + "\t"
    with opener(src) as fh:
        for ln in fh:
            if not ln.startswith(pref):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = _attrs(f[8])
            if f[2] == "gene":
                if a.get("gene_biotype") == "protein_coding":
                    name = a.get("gene_name") or a.get("Name") or a["ID"]
                    genes[a["ID"]] = dict(name=name, start=int(f[3]), end=int(f[4]), strand=f[6])
            elif f[2] == "exon":
                for p in a.get("Parent", "").split(","):
                    if p:
                        exons[p].append((int(f[3]), int(f[4])))
            elif f[2] not in ("CDS", "start_codon", "stop_codon", "five_prime_UTR", "three_prime_UTR"):
                if "ID" in a and "Parent" in a:
                    parent_of[a["ID"]] = a["Parent"].split(",")
    tx = collections.defaultdict(dict)          # gene ID -> {transcript ID: exons}
    for t, ps in parent_of.items():
        for p in ps:
            if p in genes and exons.get(t):
                tx[p][t] = sorted(set(exons[t]))
    for g in genes:                              # exons parented by the gene itself (no transcript record)
        if exons.get(g) and not tx.get(g):
            tx[g][g] = sorted(set(exons[g]))
    order = sorted(genes, key=lambda g: (genes[g]["start"], genes[g]["end"], g))
    cp, tg = Path(out + ".copies.tsv"), Path(out + ".truth.gtf")
    with open(cp, "w") as fo:
        fo.write("cid\tfamily\tname\tchrom\tterr_lo0\tterr_hi\tstrand\tisoform_gene\tbiotype\tn_tx\tsource\n")
        for g in order:
            v = genes[g]
            fo.write(f"{g}\tALL\t{v['name']}\t{contig}\t{v['start'] - 1}\t{v['end']}\t{v['strand']}\t{g}\tprotein_coding\t"
                     f"{len(tx.get(g, {}))}\t{src.name}\n")
    with open(tg, "w") as fo:
        for g in order:
            v = genes[g]
            for t in sorted(tx.get(g, {})):
                ex = tx[g][t]
                at = f'gene_id "{g}"; transcript_id "{t}"; gene_name "{v["name"]}";'
                fo.write(f"{contig}\ttruth\ttranscript\t{ex[0][0]}\t{max(e for _, e in ex)}\t.\t{v['strand']}\t.\t{at}\n")
                for s, e in ex:
                    fo.write(f"{contig}\ttruth\texon\t{s}\t{e}\t.\t{v['strand']}\t.\t{at}\n")
    n_no_tx = sum(1 for g in genes if not tx.get(g))
    print(f"[h3-inputs] {species} {contig}: {len(genes)} protein-coding genes ({n_no_tx} without exon records), "
          f"{sum(len(v) for v in tx.values())} transcripts from {src} -> {cp}, {tg}")


# ================================================================ H1
def compara_truth(cfg: dict, contig: str, d: Path) -> tuple[Path, str]:
    """Compara families at Primates restricted to the contig, as `_o1_recovery.compara_contig_truth` builds them: the cached
    genome-wide table (`compara_families`' product) filtered by `_o1.filter_rows`. `compara_contig_truth` itself refuses today
    (NotBuilt: the raw Compara download `compara_gw` is gone from disk); the cached table is used read-only, never rebuilt."""
    import _o1
    import _o1_recovery as R
    fams = _o1.species_dir(cfg, "human") / f"compara.{R.COMPARA_HEADLINE}.families.tsv"
    if not fams.exists():
        raise FileNotFoundError(f"Compara families table {fams} absent (not rebuilt here)")
    dst = _o1.filter_rows(fams, d / f"human_{contig}_compara.tsv", "Contig", keep={contig})
    fig = FIG7 / f"human_{contig}_compara.tsv"
    note = (f"restricted from {fams} ({_sha1(fams)[:12]}); identical to Figure 7's cached {fig.name}: "
            f"{fig.exists() and _sha1(fig) == _sha1(dst)}")
    return dst, note


def h1(cfg: dict, species: str, contig: str, d: Path) -> dict:
    import _liftoff as L
    import _o1_recovery as R
    genes = _copy_ro(FIG7 / "gff" / f"{species}_{contig}.gff", d / f"{species}_{contig}.gff")
    truths, notes = {}, {}
    if species == "human":
        truths["compara"], notes["compara"] = compara_truth(cfg, contig, d)
        truths["soto"] = R.bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
        if contig == "chr16":
            truths["npip_u2"] = R.npip_u2_truth(cfg)
    truths["referee"] = _copy_ro(FIG7 / f"{species}_{contig}_ref.families.tsv", d / f"{species}_{contig}_ref.families.tsv")
    spans = R.gene_spans(genes, contig)
    out = {"truth_files": {t: str(p) for t, p in truths.items()}, "truth_notes": notes, "genes_gff": str(genes)}
    loci, sup, why = L.support_if_ready(cfg, SAMPLE[species], species)
    if loci is None:
        out["liftoff_missing"] = why
    else:
        pairs = L.copy_pairs(loci, LIFTOFF_SC, sup, lambda c: c == contig)
    for arm in ARMS:
        clusters = d / f"{arm}.fam.clusters.tsv"
        res = {}
        for t, tp in truths.items():
            label = f"{species}_{contig}_{arm}_{t}"
            fs = R.family_score(cfg, clusters, genes, tp, contig, label, d / f"{label}.family_score.txt")
            mine = R.score_arm(clusters, spans, tp, contig)
            R.check_against_family_score(mine, fs, label)
            keep = ("truth_families", "truth_genes", "clusters_scored", "clusters_total", "loci_total", "sens", "prec", "f",
                    "matched", "pred_members", "pair_tp", "truth_pairs", "pred_pairs", "pair_sens", "pair_prec", "exact",
                    "touched")
            res[t] = {k: mine[k] for k in keep} | {"collapsed": fs.get("collapsed"), "no_locus": fs.get("no_locus"),
                                                   "per_family": mine["per_family"]}
        if loci is not None:
            pf = L.pair_families(pairs, L.catalog_loci(d / f"{arm}.fam.copies.tsv"))
            k = sum(1 for _, _, shared in pf if shared)
            res["liftoff"] = {"pairs": len(pf), "recovered": k, "both_covered": sum(1 for _, b, _ in pf if b),
                              "sens": k / len(pf) if pf else None,
                              "recovered_ids": sorted(sid for sid, _, shared in pf if shared),
                              "pair_ids": sorted(sid for sid, _, _ in pf)}
        out[arm] = res
    return out


# ================================================================ H2
def _exons_col(s: str):
    return sorted((int(a), int(b)) for a, b in (x.split("-") for x in s.split(",")))


def j50(exons) -> int:
    """Junctions = gaps >= 50 bp between consecutive exons (0-based half-open, coordinate-sorted)."""
    ex = sorted(exons)
    return sum(1 for a, b in zip(ex, ex[1:]) if b[0] - a[1] >= MIN_JUNCTION_GAP)


def copies_profile(path: Path) -> dict:
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    jj, jn, bp = [], [], 0
    for r in rows:
        ex = _exons_col(r["exons"])
        jj.append(j50(ex))
        jn.append(int(r["n_exon"]) - 1)
        bp += sum(e - s for s, e in ex)
    n = len(rows)
    return {"copies": n, "families": len({r["family_id"] for r in rows}),
            "median_j50": float(statistics.median(jj)) if jj else None,
            "median_nexon_m1": float(statistics.median(jn)) if jn else None,
            "mean_j50": sum(jj) / n if n else None, "mean_nexon_m1": sum(jn) / n if n else None,
            "frac_ge2_j50": sum(1 for x in jj if x >= 2) / n if n else None,
            "frac_ge2_nexon_m1": sum(1 for x in jn if x >= 2) / n if n else None,
            "copies_j50_ne_nexon_m1": sum(1 for a, b in zip(jj, jn) if a != b),
            "total_exon_bp": bp, "gene_ids": sorted(r["gene_id"] for r in rows)}


def loci_gff3(path: Path) -> dict:
    """gene Name -> (chrom, strand, exons tuple (GFF 1-based closed, file order))."""
    loci = {}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = _attrs(f[8])
        if f[2] == "gene":
            loci[a["Name"]] = [f[0], f[6], []]
        elif f[2] == "exon":
            loci[a["gene"]][2].append((int(f[3]), int(f[4])))
    return {k: (v[0], v[1], tuple(v[2])) for k, v in loci.items()}


def gtf_reps(gtf: Path) -> dict:
    """Python re-derivation of both rules from the GTF (a CHECK of the binary, mcl_families `gtf_loci`): gene_id ->
    {rule: (rep transcript_id, its exons sorted by start (stable), its junction count)}."""
    def attr(s, key):
        m = re.search(key + r' "([^"]*)"', s)
        return m.group(1) if m else None
    gene_of, reads, strand, exons, order, seen = {}, {}, {}, collections.defaultdict(list), [], set()
    for ln in open(gtf):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        t = attr(f[8], "transcript_id")
        if t is None:
            continue
        if f[2] == "transcript":
            g = attr(f[8], "gene_id") or t
            if g not in seen:
                seen.add(g)
                order.append(g)
            gene_of[t] = g
            strand[t] = f[6]
            r = attr(f[8], "reads")
            try:
                reads[t] = int(r) if r is not None else 0
            except ValueError:
                reads[t] = 0
        elif f[2] == "exon":
            exons[t].append((int(f[3]), int(f[4])))
    txs = collections.defaultdict(list)
    for t, g in gene_of.items():
        txs[g].append(t)

    def junctions(ex):
        iv = sorted(ex)
        return sum(1 for a, b in zip(iv, iv[1:]) if max(0, b[0] - (a[1] + 1)) >= MIN_JUNCTION_GAP)

    def span(t):
        ex = exons.get(t)
        return (max(e for _, e in ex) - min(s for s, _ in ex)) if ex else 0
    out = {}
    for g in order:
        ts = sorted(txs.get(g, []))
        if not any(exons.get(t) for t in ts):
            continue
        res = {}
        for rule in ARMS:
            best, bk = None, None
            for t in ts:                      # last maximum wins (Rust max_by_key)
                key = (reads.get(t, 0), span(t)) if rule == "R_M" else (junctions(exons.get(t, [])), reads.get(t, 0), span(t))
                if bk is None or key >= bk:
                    best, bk = t, key
            ex = sorted(exons.get(best, []), key=lambda x: x[0])
            res[rule] = (best, tuple(ex), junctions(ex))
        out[g] = res
    return out


def h2(d: Path, gtf: Path) -> dict:
    out = {arm: copies_profile(d / f"{arm}.fam.copies.tsv") for arm in ARMS}
    lg = {arm: loci_gff3(d / f"{arm}.fam.loci.gff3") for arm in ARMS}
    names = set(lg["R_M"]) & set(lg["R_J"])
    changed = sorted(n for n in names if lg["R_M"][n][2] != lg["R_J"][n][2])
    reps = gtf_reps(gtf)
    mism = {arm: sorted(g for g in lg[arm] if g in reps and reps[g][arm][1] != lg[arm][g][2]) for arm in ARMS}
    missing = {arm: sorted(set(lg[arm]) ^ set(reps)) for arm in ARMS}
    rep_tx_changed = sorted(g for g in reps if reps[g]["R_M"][0] != reps[g]["R_J"][0])
    changed_set = set(changed)
    copies_changed = {arm: sum(1 for g in out[arm]["gene_ids"] if g in changed_set) for arm in ARMS}
    for arm in ARMS:
        out[arm].pop("gene_ids")
    jm = [reps[g]["R_M"][2] for g in reps]
    jj = [reps[g]["R_J"][2] for g in reps]
    out["loci"] = {"loci_R_M": len(lg["R_M"]), "loci_R_J": len(lg["R_J"]), "loci_rep_exons_changed": len(changed),
                   "loci_rep_transcript_changed": len(rep_tx_changed), "copies_with_changed_rep": copies_changed,
                   "rederivation_mismatches": {a: len(v) for a, v in mism.items()},
                   "rederivation_mismatch_examples": {a: v[:10] for a, v in mism.items()},
                   "loci_not_in_both": {a: len(v) for a, v in missing.items()},
                   "all_loci_median_j50": {"R_M": statistics.median(jm), "R_J": statistics.median(jj)},
                   "all_loci_frac_ge2_j50": {"R_M": sum(1 for x in jm if x >= 2) / len(jm), "R_J": sum(1 for x in jj if x >= 2) / len(jj)}}
    return out, lg


# ================================================================ two-gene representatives (the fused-locus caveat)
def pc_genes_slice(gff: Path, contig: str) -> list:
    """Protein-coding genes of a RefSeq GFF slice (readthrough genes excluded): (start0, end, strand, name, merged exon union)."""
    genes, par, ex = {}, {}, collections.defaultdict(list)
    for ln in open(gff):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != contig:
            continue
        a = _attrs(f[8])
        if f[2] == "gene":
            if a.get("gene_biotype") == "protein_coding" and "readthrough" not in a.get("description", "").lower():
                genes[a["ID"]] = (f[6], a.get("Name", a["ID"]))
        elif f[2] == "exon":
            for p in a.get("Parent", "").split(","):
                ex[p].append((int(f[3]) - 1, int(f[4])))
        elif "ID" in a and "Parent" in a:
            par[a["ID"]] = a["Parent"].split(",")
    union = collections.defaultdict(list)
    for p, iv in ex.items():
        for g in ([p] if p in genes else [q for q in par.get(p, []) if q in genes]):
            union[g].extend(iv)
    out = []
    for g, iv in union.items():
        m = _merge(iv)
        out.append((m[0][0], m[-1][1], genes[g][0], genes[g][1], m))
    out.sort()
    return out


def two_gene_reps(lg: dict, genes: list) -> dict:
    """{locus Name: [gene names]} of the representatives whose exons overlap >= 2 mutually exon-disjoint same-strand genes."""
    import bisect
    starts = [g[0] for g in genes]
    maxlen = max((g[1] - g[0] for g in genes), default=0)
    out = {}
    for name, (chrom, strand, exons) in lg.items():
        ex = _merge([(s - 1, e) for s, e in exons])
        lo, hi = ex[0][0], ex[-1][1]
        i0 = bisect.bisect_left(starts, lo - maxlen)
        hit = []
        for g in genes[i0:bisect.bisect_left(starts, hi)]:
            if g[1] <= lo or g[2] != strand:
                continue
            if _inter(ex, g[4]) > 0:
                hit.append(g)
        if len(hit) < 2:
            continue
        disjoint = any(_inter(a[4], b[4]) == 0 for x, a in enumerate(hit) for b in hit[x + 1:])
        if disjoint:
            out[name] = sorted({g[3] for g in hit})
    return out


def changed_two_gene(d: Path, h3: dict, lg: dict, tg: dict) -> dict:
    """For each H3 gene that changes: does a two-gene representative of the arm that FINDS it overlap the gene's exon union on its
    strand? (An upper bound on the gains / losses a fused representative can explain: overlap, not the junction-carrying locus.)"""
    ids = {x["cid"] for x in h3["changed"]}
    union, strand = collections.defaultdict(list), {}
    for ln in open(d / "h3.truth.gtf"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        if g in ids:
            union[g].append((int(f[3]) - 1, int(f[4])))
            strand[g] = f[6]
    out = {"gains": 0, "gains_with_two_gene_R_J_locus": 0, "losses": 0, "losses_with_two_gene_R_M_locus": 0}
    for x in h3["changed"]:
        u = _merge(union.get(x["cid"], []))
        gain = bool(x["R_J_found"]) and not x["R_M_found"]
        arm = "R_J" if gain else "R_M"
        hit = [n for n in tg[arm] if lg[arm][n][1] == strand.get(x["cid"])
               and _inter(_merge([(s - 1, e) for s, e in lg[arm][n][2]]), u) > 0]
        x["two_gene_locus_in_finding_arm"] = ";".join(f"{n}:{'+'.join(tg[arm][n])}" for n in hit)
        if gain:
            out["gains"] += 1
            out["gains_with_two_gene_R_J_locus"] += bool(hit)
        else:
            out["losses"] += 1
            out["losses_with_two_gene_R_M_locus"] += bool(hit)
    return out


# ================================================================ H3
def h3_read(prefix: Path) -> dict | None:
    js, tsv = Path(f"{prefix}.json"), Path(f"{prefix}.copies.tsv")
    if not js.exists():
        return None
    summ = json.load(open(js))
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    out = {"copies": summ["copies"], "spliced_expressed": summ["spliced_expressed"], "arms": summ["arms"], "by_k": {},
           "changed": []}
    for k in ("0", "1", "2"):
        sub = [r for r in rows if r["k"] == k]
        out["by_k"][k] = {"genes": len(sub), "spliced_expressed": sum(int(r["spliced_expressed"]) for r in sub),
                          **{f"{a}_found": sum(int(r[f"{a}_strict_found"]) for r in sub) for a in ARMS}}
    for r in rows:
        if r["R_M_strict_found"] != r["R_J_strict_found"]:
            out["changed"].append({"cid": r["cid"], "name": r["name"], "k": int(r["k"]), "support_reads": int(r["support_reads"]),
                                   "supported_junctions": int(r["supported_junctions"]),
                                   **{f"{a}_found": int(r[f"{a}_strict_found"]) for a in ARMS},
                                   **{f"{a}_rep_supported_junctions_max": int(r[f"{a}_rep_supported_junctions_max"]) for a in ARMS},
                                   **{f"{a}_locus_supported_junctions_max": int(r.get(f"{a}_locus_supported_junctions_max", -1) or -1)
                                      for a in ARMS}})
    out["rows"] = [{k: r[k] for k in ("cid", "name", "k", "spliced_expressed", "support_reads", "supported_junctions",
                                       "R_M_strict_found", "R_J_strict_found", "R_M_rep_supported_junctions_max",
                                       "R_J_rep_supported_junctions_max", "R_M_locus_level_found", "R_J_locus_level_found")}
                   for r in rows] if len(rows) <= 50 else None
    return out


# ================================================================ provenance
def run_log(p: Path) -> dict:
    out = {}
    for ln in open(p):
        ln = ln.rstrip("\n")
        if "\t" in ln and not ln.startswith(("\t", " ")):
            k, _, v = ln.partition("\t")
            out.setdefault(k, v)
        m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", ln)
        if m:
            parts = [float(x) for x in m.group(1).split(":")]
            out["wall_s"] = sum(x * 60 ** i for i, x in enumerate(reversed(parts)))
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", ln)
        if m:
            out["max_rss_kb"] = int(m.group(1))
    return out


def score(species: str, contig: str):
    import figlib
    cfg = figlib.load_inputs()
    cfg["bin"] = str(BIN)                  # family_score (R.family_score) runs from REP_BIN
    d = W / f"{species}_{contig}"
    logs = {arm: run_log(d / f"{arm}.run.log") for arm in ARMS}
    bin_sha = {}
    for b in ("mcl_families", "family_score"):
        seen = {logs[a].get(p + f"{b}_sha1") for a in ARMS for p in ("", "after_")}
        if len(seen) != 1 or None in seen:
            raise RuntimeError(f"{d}: the two arms did not run one {b} (sha1 before/after, both arms: {seen})")
        bin_sha[b] = seen.pop()
    for a in ARMS:
        if logs[a].get("exit") != "0":
            raise RuntimeError(f"{d}/{a}.run.log: exit {logs[a].get('exit')}")
    if logs["R_J"].get("params_representative_row") != "present" or logs["R_M"].get("params_representative_row") != "absent":
        raise RuntimeError(f"{d}: params.tsv representative row not as registered")
    cur = {k: _sha1(BIN / k) for k in bin_sha}
    if cur != bin_sha:
        raise RuntimeError(f"the binary changed since the runs: {bin_sha} -> {cur}")
    gtf = d / f"{species}_{contig}.denovo.gtf"
    res = {"species": species, "contig": contig, "status": STATUS[(species, contig)], "binary_sha1": bin_sha,
           "gtf": str(gtf), "gtf_sha1": _sha1(gtf),
           "paf_identical": _sha1(d / "R_M.fam.loci.paf") == _sha1(d / "R_J.fam.loci.paf"),
           "loci_fa_identical": _sha1(d / "R_M.fam.loci.fa") == _sha1(d / "R_J.fam.loci.fa"),
           "clusters_identical": _sha1(d / "R_M.fam.clusters.tsv") == _sha1(d / "R_J.fam.clusters.tsv"),
           "runs": {a: {k: logs[a].get(k) for k in ("wall_s", "max_rss_kb", "date", "rule", "bridge_regroup")} for a in ARMS}}
    res["H1"] = h1(cfg, species, contig, d)
    res["H2"], lg = h2(d, gtf)
    genes = pc_genes_slice(d / f"{species}_{contig}.gff", contig)
    copies = {a: {r["gene_id"] for r in csv.DictReader(open(d / f"{a}.fam.copies.tsv"), delimiter="\t")} for a in ARMS}
    tg = {a: two_gene_reps(lg[a], genes) for a in ARMS}
    res["two_gene_reps"] = {
        "genes_used": len(genes),
        **{a: {"all_loci": len(tg[a]), "copies": sum(1 for n in tg[a] if n in copies[a])} for a in ARMS},
        "only_R_J": sorted(set(tg["R_J"]) - set(tg["R_M"])), "only_R_M": sorted(set(tg["R_M"]) - set(tg["R_J"])),
        "only_R_J_genes": {n: tg["R_J"][n] for n in sorted(set(tg["R_J"]) - set(tg["R_M"]))[:40]}}
    h3 = h3_read(d / "h3.support")
    if h3 is not None:
        h3["run"] = {k: v for k, v in run_log(d / "h3.run.log").items() if k in ("wall_s", "max_rss_kb", "exit")}
        h3["two_gene"] = changed_two_gene(d, h3, lg, tg)
    res["H3"] = h3
    if species == "human" and contig == "chr16":
        res["H3_npip"] = h3_read(d / "npip.support")
        res["H3_npip"]["run"] = {k: v for k, v in run_log(d / "npip.run.log").items() if k in ("wall_s", "max_rss_kb", "exit")}
    js = d / "score.json"
    js.write_text(json.dumps(res, indent=1, default=str))
    print(f"[score] {species} {contig} -> {js}")
    _print_contig(res)


def _f(x, n=3):
    return "NA" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))


def _print_contig(res: dict):
    h1r = res["H1"]
    for t in [t for t in ("compara", "soto", "npip_u2", "referee", "liftoff") if t in h1r["R_M"]]:
        a, b = h1r["R_M"][t], h1r["R_J"][t]
        if t == "liftoff":
            print(f"  H1 {t:8s} R_M {a['recovered']}/{a['pairs']}  R_J {b['recovered']}/{b['pairs']}")
        else:
            print(f"  H1 {t:8s} R_M F {_f(a['f'])} s {_f(a['sens'])} p {_f(a['prec'])} | R_J F {_f(b['f'])} s {_f(b['sens'])} "
                  f"p {_f(b['prec'])}")
    h2r = res["H2"]
    for a in ARMS:
        x = h2r[a]
        print(f"  H2 {a} copies {x['copies']} fam {x['families']} median j50 {x['median_j50']} (n_exon-1 {x['median_nexon_m1']}) "
              f">=2 {_f(x['frac_ge2_j50'])} bp {x['total_exon_bp']}")
    print(f"  H2 loci {h2r['loci']}")
    if res.get("H3"):
        h = res["H3"]
        print(f"  H3 genes {h['copies']} spliced-expressed {h['spliced_expressed']} found R_M {h['arms']['R_M']['strict_found']} "
              f"R_J {h['arms']['R_J']['strict_found']} changed {len(h['changed'])}")
    print(f"  two-gene reps {res['two_gene_reps']['R_M']} vs {res['two_gene_reps']['R_J']}")


# ================================================================ decision
def write_changed(allres: dict, path: str):
    """One row per H3 gene whose strict FOUND differs between the arms, every contig."""
    cols = ["species", "contig", "status", "cid", "name", "k", "support_reads", "supported_junctions", "R_M_found", "R_J_found",
            "R_M_rep_supported_junctions_max", "R_J_rep_supported_junctions_max", "R_M_locus_supported_junctions_max",
            "R_J_locus_supported_junctions_max", "two_gene_locus_in_finding_arm"]
    with open(path, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for (sp, c), r in allres.items():
            for x in (r.get("H3") or {}).get("changed", []):
                row = {"species": sp, "contig": c, "status": r["status"], **x}
                fo.write("\t".join(str(row.get(k, "")) for k in cols) + "\n")


def write_npip(allres: dict, path: str):
    n = (allres.get(("human", "chr16")) or {}).get("H3_npip")
    if not n or not n.get("rows"):
        return
    cols = list(n["rows"][0].keys())
    with open(path, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for x in n["rows"]:
            fo.write("\t".join(str(x[k]) for k in cols) + "\n")


def decide(md: str | None, changed_tsv: str | None = None, npip_tsv: str | None = None):
    allres = {}
    for sp, c in DEV + HELD_OUT:
        p = W / f"{sp}_{c}" / "score.json"
        if p.exists():
            allres[(sp, c)] = json.load(open(p))
    if changed_tsv:
        write_changed(allres, changed_tsv)
    if npip_tsv:
        write_npip(allres, npip_tsv)
    clauses, verdict = [], True
    for sp, c in HELD_OUT:
        r = allres.get((sp, c))
        if r is None or r.get("H3") is None:
            clauses.append({"species": sp, "contig": c, "missing": True})
            verdict = False
            continue
        prim = PRIMARY[sp]
        a, b = r["H1"]["R_M"].get(prim), r["H1"]["R_J"].get(prim)
        row = {"species": sp, "contig": c, "status": r["status"], "primary": prim}
        if prim == "liftoff":
            if a is None or not a["pairs"]:
                row["a"] = {"pass": None, "why": f"Liftoff pairs on the contig: {a['pairs'] if a else 'absent'} (sensitivity undefined)"}
            else:
                ok = b["sens"] >= a["sens"] - F_TOL
                row["a"] = {"pass": ok, "R_M_sens": a["sens"], "R_J_sens": b["sens"], "pairs": a["pairs"]}
        else:
            ok_f = b["f"] >= a["f"] - F_TOL
            ok_p = b["prec"] >= a["prec"] - P_TOL
            row["a"] = {"pass": ok_f and ok_p, "R_M_f": a["f"], "R_J_f": b["f"], "R_M_prec": a["prec"], "R_J_prec": b["prec"],
                        "f_pass": ok_f, "prec_pass": ok_p}
        fm, fj = r["H3"]["arms"]["R_M"]["strict_found"], r["H3"]["arms"]["R_J"]["strict_found"]
        row["b"] = {"pass": fj >= fm, "R_M_found": fm, "R_J_found": fj}
        mm, mj = r["H2"]["R_M"]["median_j50"], r["H2"]["R_J"]["median_j50"]
        nm, nj = r["H2"]["R_M"]["median_nexon_m1"], r["H2"]["R_J"]["median_nexon_m1"]
        row["c"] = {"pass": mj >= mm, "R_M_median_j50": mm, "R_J_median_j50": mj, "R_M_median_nexon_m1": nm,
                    "R_J_median_nexon_m1": nj, "pass_nexon_m1": nj >= nm}
        for k in ("a", "b", "c"):
            if row[k]["pass"] is False:
                verdict = False
        clauses.append(row)
    out = {"adopt_R_J_as_default": verdict, "clauses": clauses,
           "rule": ("R_J becomes the default iff on every held-out contig (a) F vs the primary reference >= R_M - 0.005 and "
                    "precision >= R_M - 0.01 where defined; (b) H3 found >= R_M; (c) H2 median junctions per copy >= R_M")}
    (W / "decision.json").write_text(json.dumps(out, indent=1, default=str))
    text = tables_md(allres, out)
    if md:
        Path(md).write_text(text)
    print(text)


def tables_md(allres: dict, dec: dict) -> str:
    L = []
    order = [k for k in DEV + HELD_OUT if k in allres]
    L.append("### H1 — families (bipartite F / sensitivity / precision; Liftoff = pairs recovered / pairs)\n")
    L.append("| contig | status | arm | Compara F (s / p) | Soto F (s / p) | NPIP set F | protein homology F (s / p) | Liftoff pairs |")
    L.append("|---|---|---|---|---|---|---|---|")
    for k in order:
        r = allres[k]
        for arm in ARMS:
            h = r["H1"][arm]

            def fsp(t):
                if t not in h:
                    return "—"
                x = h[t]
                return f"{_f(x['f'])} ({_f(x['sens'])} / {_f(x['prec'])})"
            lf = h.get("liftoff")
            L.append(f"| {k[1]} | {r['status']} | {arm} | {fsp('compara')} | {fsp('soto')} | "
                     f"{_f(h['npip_u2']['f']) if 'npip_u2' in h else '—'} | {fsp('referee')} | "
                     f"{(str(lf['recovered']) + '/' + str(lf['pairs'])) if lf else 'absent'} |")
    L.append("\n### H2 — copy table structure (junction = gap >= 50 bp between exons; `n_exon - 1` in brackets)\n")
    L.append("| contig | arm | copies | families | median junctions / copy | copies with >= 2 junctions | total exon bp | loci whose rep changed |")
    L.append("|---|---|---|---|---|---|---|---|")
    for k in order:
        r = allres[k]
        for arm in ARMS:
            x = r["H2"][arm]
            ch = r["H2"]["loci"]
            L.append(f"| {k[1]} | {arm} | {x['copies']} | {x['families']} | {_f(x['median_j50'], 1)} [{_f(x['median_nexon_m1'], 1)}] | "
                     f"{_f(x['frac_ge2_j50'])} [{_f(x['frac_ge2_nexon_m1'])}] | {x['total_exon_bp']:,} | "
                     f"{(str(ch['loci_rep_exons_changed']) + ' / ' + str(ch['loci_R_M'])) if arm == 'R_J' else ''} |")
    L.append("\n### H3 — annotated protein-coding genes FOUND (strict rule, bench/copy_support.py)\n")
    L.append("| contig | genes | spliced-expressed | found R_M | found R_J | R_J gains (a two-gene R_J rep overlaps) | R_J losses (a two-gene R_M rep overlaps) | k=2 genes: found R_M / R_J | locus level R_M / R_J | two-gene reps R_M / R_J, all loci (in the copy table) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for k in order:
        r = allres[k]
        h = r.get("H3")
        t = r["two_gene_reps"]
        tg = f"{t['R_M']['all_loci']} / {t['R_J']['all_loci']} ({t['R_M']['copies']} / {t['R_J']['copies']})"
        if not h:
            L.append(f"| {k[1]} | — | — | — | — | — | — | — | — | {tg} |")
            continue
        g = sum(1 for x in h["changed"] if x["R_J_found"] and not x["R_M_found"])
        lo = sum(1 for x in h["changed"] if x["R_M_found"] and not x["R_J_found"])
        k2 = h["by_k"]["2"]
        tw = h.get("two_gene", {})
        L.append(f"| {k[1]} | {h['copies']} | {h['spliced_expressed']} | {h['arms']['R_M']['strict_found']} | "
                 f"{h['arms']['R_J']['strict_found']} | {g} ({tw.get('gains_with_two_gene_R_J_locus', '?')}) | "
                 f"{lo} ({tw.get('losses_with_two_gene_R_M_locus', '?')}) | {k2['R_M_found']} / {k2['R_J_found']} | "
                 f"{h['arms']['R_M']['locus_level_found']} / {h['arms']['R_J']['locus_level_found']} | {tg} |")
    L.append("\n### Run times (`/usr/bin/time -v`; families = one driver call per arm, H3 = one copy_support.py call per contig)\n")
    L.append("| contig | families R_M wall s (max RSS GB) | families R_J wall s (max RSS GB) | H3 wall s (max RSS GB) |")
    L.append("|---|---|---|---|")
    for k in order:
        r = allres[k]

        def wt(x):
            return f"{float(x['wall_s']):.0f} ({int(x['max_rss_kb']) / 1048576:.2f})" if x and x.get("wall_s") is not None else "—"
        L.append(f"| {k[1]} | {wt(r['runs']['R_M'])} | {wt(r['runs']['R_J'])} | {wt((r.get('H3') or {}).get('run'))} |")
    if allres.get(("human", "chr16"), {}).get("H3_npip"):
        n = allres[("human", "chr16")]["H3_npip"]
        L.append("\n### H3 on chr16 — the 25 NPIP copies (strict found; locus level beside)\n")
        L.append("| arm | strict found | locus-level found | any same-strand overlapping locus |")
        L.append("|---|---|---|---|")
        for arm in ARMS:
            x = n["arms"][arm]
            L.append(f"| {arm} | {x['strict_found']} | {x['locus_level_found']} | {x['old_overlap']} |")
    L.append("\n### Decision (held-out contigs only)\n")
    L.append("| contig | (a) primary reference | (b) H3 found | (c) H2 median junctions / copy |")
    L.append("|---|---|---|---|")
    for row in dec["clauses"]:
        if row.get("missing"):
            L.append(f"| {row['contig']} | missing | missing | missing |")
            continue
        a = row["a"]
        if row["primary"] == "liftoff":
            at = (f"n/a — {a['why']}" if a["pass"] is None else
                  f"{'PASS' if a['pass'] else 'FAIL'}: Liftoff sens {_f(a['R_J_sens'])} vs {_f(a['R_M_sens'])} ({a['pairs']} pairs)")
        else:
            at = (f"{'PASS' if a['pass'] else 'FAIL'}: Compara F {_f(a['R_J_f'])} vs {_f(a['R_M_f'])} (bar {_f(a['R_M_f'] - F_TOL)}), "
                  f"prec {_f(a['R_J_prec'])} vs {_f(a['R_M_prec'])} (bar {_f(a['R_M_prec'] - P_TOL)})")
        b, c = row["b"], row["c"]
        L.append(f"| {row['contig']} | {at} | {'PASS' if b['pass'] else 'FAIL'}: {b['R_J_found']} vs {b['R_M_found']} | "
                 f"{'PASS' if c['pass'] else 'FAIL'}: {_f(c['R_J_median_j50'], 1)} vs {_f(c['R_M_median_j50'], 1)} "
                 f"[n_exon-1: {_f(c['R_J_median_nexon_m1'], 1)} vs {_f(c['R_M_median_nexon_m1'], 1)}] |")
    L.append(f"\n**Decision: {'R_J becomes the default' if dec['adopt_R_J_as_default'] else 'R_M stays the default; R_J stays opt-in'}.**\n")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("h3-inputs")
    a.add_argument("--species", required=True, choices=["human", "gorilla"])
    a.add_argument("--contig", required=True)
    a.add_argument("--out", required=True)
    s = sub.add_parser("score")
    s.add_argument("--species", required=True, choices=["human", "gorilla"])
    s.add_argument("--contig", required=True)
    dd = sub.add_parser("decide")
    dd.add_argument("--md", help="also write the markdown tables here")
    dd.add_argument("--changed-tsv", help="write the H3 genes whose FOUND differs between the arms (every contig) here")
    dd.add_argument("--npip-tsv", help="write the 25 NPIP copies' H3 rows (chr16) here")
    x = ap.parse_args(argv)
    if x.cmd == "h3-inputs":
        h3_inputs(x.species, x.contig, x.out)
    elif x.cmd == "score":
        score(x.species, x.contig)
    else:
        decide(x.md, x.changed_tsv, x.npip_tsv)


if __name__ == "__main__":
    main()
