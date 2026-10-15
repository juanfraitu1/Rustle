#!/usr/bin/env python3
"""De novo analogue of Liftoff's `-copies` (prereg docs/PREREG_denovo_copies_2026-10-14.md).

    copies.py --species S --sample ID --base BASE.gtf --run DIR --out J [--liftoff ROOT]

Liftoff `-copies` lifts every annotated gene and keeps extra placements at sequence identity >= -sc. Here the templates
of a template-copies run (DIR/templates.tsv, one representative transcript per base locus) are the genes, and their
projections (DIR/batches/*.paf, minimap2 -x splice -c -N 50 -p 0.1) are the alignments. A projected copy is an alignment
passing Liftoff's thresholds of figure 8's run, computed as Liftoff computes them (find_best_mapping.
convert_all_children_coords) over the template's bases: coverage = aligned template bases / template length >= -a (0.5);
identity = matching bases / (template length + target-only bases) >= -sc (0.95).

recall  Liftoff extra copies (not short). Reachable: the source gene's placement (in_place / moved, not short) is
        covered >= 50% by a base locus; the own templates are those loci's templates (the analogue of lifting that
        gene). Found (own): a projected copy of an own template covers >= 50% of the copy's exons; found (any): a
        projected copy of any template does. Misses: an own alignment covers the copy but fails a threshold
        (below_threshold); an own template reached minimap2's cap of 51 alignments (capped); neither (no_alignment).
sites   projected copies off their own locus (the template's base locus covers < 50% of the hit), merged into sites
        (one contig, either strand, overlap >= 50% of the smaller exon union), classed by the best hit's exons:
        Liftoff extra copy, Liftoff gene (in_place / moved), assembled base locus only, none. Retro-like: the best hit
        is unspliced and its template has >= 2 exons. Unexplained: class none and not retro-like.
C1      recall (own) >= 0.80 with >= 30 reachable copies (fewer: not powered, no verdict).
C2      unexplained sites <= 0.10 of sites.
"""
import argparse
import collections
import glob
import json
import os
import re
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import guided_pipeline as GP  # noqa: E402
import _liftoff as LF  # noqa: E402
import templates as T  # noqa: E402

LIFTOFF = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
CAP = 51                 # minimap2 -N 50: one primary and at most 50 secondary alignments
SELF_GATE = 0.90         # share of templates with an alignment on their own exons (else index and base disagree)
C1_BAR, C1_MIN, C2_BAR = 0.80, 30, 0.10
PLACED = ("in_place", "moved", "extra_copy")
CLASSES = ("liftoff_copy", "liftoff_gene", "assembled", "none")


def metrics(qlen, matches, cg):
    """Liftoff's (coverage, sequence identity) of one alignment over the template's bases: coverage = aligned template
    bases (M/=/X) / template length; identity = matches / (template length + target-only bases, CIGAR D)."""
    m = d = 0
    for n, op in CIGAR.findall(cg):
        if op in "M=X":
            m += int(n)
        elif op == "D":
            d += int(n)
    return m / qlen, matches / (qlen + d)


def cov(a, b):
    """|a ∩ b| / |a| for exon unions a, b."""
    n = LF.iv_len(a)
    return LF.iv_inter(a, b) / n if n else 0.0


def batch_pafs(run):
    """Every batch's PAF, in batch order; refuses a run with an unprojected batch."""
    fas = sorted(glob.glob(os.path.join(run, "batches", "b*.fa")))
    missing = [os.path.basename(f) for f in fas if not os.path.exists(f[:-3] + ".done")]
    if not fas or missing:
        raise SystemExit(f"{run}: every batch must be projected first; missing: {missing}")
    return [f[:-3] + ".paf" for f in fas]


def scan(pafs, templates, keep):
    """One pass over every alignment. -> projected: the alignments passing Liftoff's thresholds; kept: every alignment
    of the templates in `keep` (for the miss decomposition); alignments: count per template; self: the templates with
    an alignment covering >= 50% of their own exons. Rows are {template, contig, strand, iv, cov, seq}."""
    out = dict(projected=[], kept=[], alignments=collections.Counter(), self=set())
    for path in pafs:
        for h in GP.parse_paf(path):
            if not h["cg"] or not h["qlen"]:
                continue
            q = h["q"]
            out["alignments"][q] += 1
            c, s = metrics(h["qlen"], h["nm"], h["cg"])
            row = dict(template=q, contig=h["chrom"], strand=h["strand"], iv=GP.tx_exon_blocks(h), cov=c, seq=s)
            t = templates.get(q)
            if t is not None and row["contig"] == t.contig and cov(t.exons, row["iv"]) >= LF.MATCH_COV:
                out["self"].add(q)
            if c >= LF.LIFT_A and s >= LF.LIFT_SC:
                out["projected"].append(row)
            if q in keep:
                out["kept"].append(row)
    out["alignments"] = dict(out["alignments"])
    return out


def self_rate(sc, templates):
    return len(sc["self"] & set(templates)) / len(templates) if templates else 0.0


def check_self(sc, templates):
    r = self_rate(sc, templates)
    if r < SELF_GATE:
        raise SystemExit(f"only {r:.3f} of templates align to their own exons (< {SELF_GATE}): the index and the base "
                         "GTF are not the same assembly")
    return r


def own_templates(loci, base_rows, templates):
    """-> ({source_id: the templates of the base loci covering its placement >= 50%}, the union of those sets)."""
    tid_of = {t.gene: tid for tid, t in templates.items()}
    bidx = LF.Index(base_rows)
    own = collections.defaultdict(set)
    for r in loci:
        if r["cls"] in ("in_place", "moved") and not r["short"]:
            for i in bidx.near(r["contig"], r["iv"]):
                b = base_rows[i]
                if cov(r["iv"], b["iv"]) >= LF.MATCH_COV and b["locus"] in tid_of:
                    own[r["source_id"]].add(tid_of[b["locus"]])
    placed = {r["source_id"] for r in loci if r["cls"] in ("in_place", "moved") and not r["short"]}
    return {s: own.get(s, set()) for s in placed}, set().union(*own.values()) if own else set()


def _covering(idx, r, pred=lambda h: True):
    """The rows of `idx` that cover >= 50% of locus r's exons."""
    return [idx.rows[i] for i in idx.near(r["contig"], r["iv"])
            if pred(idx.rows[i]) and cov(r["iv"], idx.rows[i]["iv"]) >= LF.MATCH_COV]


def _bin(n):
    return "1" if n == 1 else ("2-5" if n <= 5 else ">=6")


def recall(loci, base_rows, templates, sc):
    own, _ = own_templates(loci, base_rows, templates)
    copies = [r for r in loci if r["cls"] == "extra_copy" and not r["short"]]
    per_source = collections.Counter(r["source_id"] for r in copies)
    pidx, kidx = LF.Index(sc["projected"]), LF.Index(sc["kept"])
    rec = dict(copies=len(copies), reachable=0, unreachable_unplaced=0, unreachable_unassembled=0, found_own=0,
               found_any=0, misses=dict(below_threshold=0, capped=0, no_alignment=0),
               by_copy_count={k: [0, 0] for k in ("1", "2-5", ">=6")}, below_threshold_best=[])
    for r in copies:
        s = r["source_id"]
        if s not in own:
            rec["unreachable_unplaced"] += 1
            continue
        if not own[s]:
            rec["unreachable_unassembled"] += 1
            continue
        rec["reachable"] += 1
        hit = _covering(pidx, r, lambda h: h["template"] in own[s])
        rec["found_own"] += bool(hit)
        rec["found_any"] += bool(hit or _covering(pidx, r))
        b = rec["by_copy_count"][_bin(per_source[s])]
        b[0], b[1] = b[0] + bool(hit), b[1] + 1
        if hit:
            continue
        near = _covering(kidx, r, lambda h: h["template"] in own[s])
        if near:
            rec["misses"]["below_threshold"] += 1
            best = max(near, key=lambda h: (h["seq"], h["cov"]))
            rec["below_threshold_best"].append([round(best["seq"], 4), round(best["cov"], 4)])
        elif any(sc["alignments"].get(t, 0) >= CAP for t in own[s]):
            rec["misses"]["capped"] += 1
        else:
            rec["misses"]["no_alignment"] += 1
    rec["below_threshold_best"].sort()
    return rec


def site_groups(off):
    """Copies whose exon unions overlap by >= 50% of the smaller one, transitively -> [[copy, ...], ...].

    Swept by exon block. The open blocks are kept per component, so a block is compared only with other components'
    open blocks; the exonic bases a pair shares so far are summed, and the pair joins as soon as the sum reaches the
    bar (the sum only grows, so this is the rule itself). A component found to share enough with the block is not
    scanned further. A span sweep is quadratic here: one assembly put ~5,000 loci in a tandem array whose copies all
    overlap (A119b chr13 13.9-15.5 Mb)."""
    n = len(off)
    size = [LF.iv_len(h["iv"]) for h in off]
    last = [h["iv"][-1][1] for h in off]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        if len(opened.get(a, ())) < len(opened.get(b, ())):
            a, b = b, a
        parent[b] = a
        if b in opened:
            opened.setdefault(a, []).extend(opened.pop(b))
        return a

    blocks = sorted((h["contig"], s, e, i) for i, h in enumerate(off) for s, e in h["iv"])
    opened, shared, contig = {}, {}, None
    for tick, (c, s, e, i) in enumerate(blocks):
        if c != contig:
            opened, shared, contig = {}, {}, c
        if tick % 4096 == 0:
            opened = {r: v for r, v in ((r, [x for x in v if x[0] > s]) for r, v in opened.items()) if v}
            shared = {k: v for k, v in shared.items() if last[k[0]] > s and last[k[1]] > s}
        ri = find(i)
        for r in list(opened):
            if r == ri or r not in opened:
                continue
            for e2, j in opened[r]:
                if e2 <= s:
                    continue
                k = (j, i) if j < i else (i, j)
                v = shared.get(k, 0) + min(e, e2) - s
                if v >= 0.5 * min(size[i], size[j]):
                    shared.pop(k, None)
                    ri = union(ri, r)
                    break
                shared[k] = v
        opened.setdefault(ri, []).append((e, i))
    groups = collections.defaultdict(list)
    for i, h in enumerate(off):
        groups[find(i)].append(h)
    return list(groups.values())


def sites(loci, base_rows, templates, projected):
    """Projected copies off their own locus, merged into sites and classed (module docstring)."""
    locus_iv = {(b["contig"], b["locus"]): b["iv"] for b in base_rows}
    off = []
    for h in projected:
        t = templates.get(h["template"])
        own = locus_iv.get((t.contig, t.gene), t.exons) if t is not None else []
        if not (t is not None and h["contig"] == t.contig and cov(h["iv"], own) >= LF.MATCH_COV):
            off.append(h)
    groups = site_groups(off)
    extra = LF.Index([r for r in loci if r["cls"] == "extra_copy"])
    genes = LF.Index([r for r in loci if r["cls"] in ("in_place", "moved")])
    bidx = LF.Index(base_rows)
    classes = {k: [0, 0] for k in CLASSES}
    for g in groups:
        best = min(g, key=lambda h: (-h["seq"], -h["cov"], h["template"], h["strand"]))
        k = next((name for name, idx in (("liftoff_copy", extra), ("liftoff_gene", genes), ("assembled", bidx))
                  if idx.best_cov(best["contig"], best["iv"])[0] >= LF.MATCH_COV), "none")
        retro = len(best["iv"]) == 1 and len(templates[best["template"]].exons) >= 2
        classes[k][0] += 1
        classes[k][1] += retro
    n = len(groups)
    unexplained = classes["none"][0] - classes["none"][1]
    return dict(n=n, classes=classes, unexplained=unexplained, unexplained_share=unexplained / n if n else 0.0,
                templates_with_sites=len({h["template"] for h in off}))


def verdicts(rec, st):
    r = rec["found_own"] / rec["reachable"] if rec["reachable"] else float("nan")
    powered = rec["reachable"] >= C1_MIN
    return (dict(bar=C1_BAR, recall_own=r, reachable=rec["reachable"], powered=powered,
                 passes=bool(r >= C1_BAR) if powered else None),
            dict(bar=C2_BAR, unexplained_share=st["unexplained_share"], passes=bool(st["unexplained_share"] <= C2_BAR)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--species", "--sample", "--base", "--run", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--liftoff", default=LIFTOFF)
    a = ap.parse_args(argv)
    pafs = batch_pafs(a.run)
    templates = T.read_templates(os.path.join(a.run, "templates.tsv"))
    loci = LF.read_loci(Path(a.liftoff) / a.species / "liftoff_loci.tsv")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    base_rows = LF.denovo_loci(Path(a.base), Path(os.path.dirname(os.path.abspath(a.out))) / "denovo_cache")
    sc = scan(pafs, templates, own_templates(loci, base_rows, templates)[1])
    rate = check_self(sc, templates)
    rec = recall(loci, base_rows, templates, sc)
    rec["recall_own"] = rec["found_own"] / rec["reachable"] if rec["reachable"] else None
    rec["recall_any"] = rec["found_any"] / rec["reachable"] if rec["reachable"] else None
    st = sites(loci, base_rows, templates, sc["projected"])
    c1, c2 = verdicts(rec, st)
    res = dict(sample=a.sample, species=a.species, run=a.run, base=a.base, templates=len(templates),
               alignments=sum(sc["alignments"].values()), projected=len(sc["projected"]),
               capped_templates=sum(v >= CAP for v in sc["alignments"].values()), self_rate=rate,
               recall=rec, sites=st, c1=c1, c2=c2)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps(dict(sample=a.sample, c1=c1, c2=c2, recall_any=rec["recall_any"]), sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
