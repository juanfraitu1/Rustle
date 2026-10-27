#!/usr/bin/env python3
"""Score the nodes of one arm (an ideal-expression run, or a real-library product) against the Soto genes (docs/PREREG_ideal_arm_pilot_2026-10-07.md and its amendment). Descriptive: it prints
tables and one template sentence, no verdict, and no cut that moves with the arm.

    score_soto_nodes.py --truth SIMPREFIX --out OUT  (--unit node --loci-gff3 PREFIX.fam.loci.gff3 | --unit transcript|locus --gtf ANY.gtf)
                        (--freeze-single-cut CUT.json | --single-cut-file CUT.json | --no-single-cut)
                        [--compete all|soto] [--fused-rule span|exon] [--bam BAM] [--assemble-log LOG] [--real-library [--reads-per-gene TSV --reads-column COL]]
                        [--subset IDS [--subset-label L]] [--compare LABEL=SUMMARY.json ...]

POPULATION: the listed (Soto) genes that can be simulated, by S1C Gene ID (never by name): the denominator is the same in every arm. The genes no rule lets us simulate are listed apart as
`not_simulable`. The gene UNIVERSE (competition for nodes, fusion) is every gene with an eligible transcript in the truth tables: it comes from the ideal arm and does not shrink with an
arm's reads. Truth = sim_genes.py tables (targets.tsv, transcripts.tsv); a gene's exons are the union of the canonical blocks of its eligible transcripts and its chains are theirs.

UNIT: `node` = the families stage's representative of a locus (PREFIX.fam.loci.gff3), the primary unit; `transcript` = a transcript of any GTF; `locus` = a gene_id of the GTF (union of its
transcripts). The unit moves the headline by ten points on one run, so it is printed on every table.

MATCH: candidates are (gene, unit) pairs whose EXONS overlap; they are taken one-to-one, greedily, by exonic overlap bp (descending), then exonic Jaccard (descending), then gene id and unit id
(ascending): that is the tie policy, and `tie sensitive genes` counts the scored genes whose outcome would change under the reverse order. A match is a recovery only if the shared exon bp are
>= 0.5 of the SMALLER of the two exon unions (two-sided); a one-bp overlap is `weak_overlap`. Strand is ignored (recorded as strand_agree).

CLASSES (ordered, mutually exclusive; `not_simulable` is listed apart):
  no_reads            the gene has no reads in the arm (simulated arms: the truth table; real arms: the reads table). A node on it is counted apart (`shared`).
  reads_not_on_gene   it has reads but none has a primary alignment on its exons (needs --bam; reads are matched to their source by name `transcript|gene|fl|k`)
  no_locus            no unit overlaps its exons (reported separately for single-exon and multi-exon genes)
  outcompeted         units overlap it but every one was taken by a gene that overlaps it more
  weak_overlap        matched, but the shared exon bp are < 0.5 of the smaller exon union
  fused_excess        the matched unit covers a gene that the gene's OWN annotation (its exon union taken as a unit) does not cover. Perfect nodes show 21-27% fused by the span rule and
                      9-11% by the exon rule, so only the excess is a loss; both rules are always computed, `--fused-rule` picks the one that decides the class
  extent_wrong        multi-exon gene: the unit's intron chain is not the chain of any annotated isoform of the gene (the end exons are the chain's; the outer ends are jittered and not compared).
                      single-exon gene: the unit is spliced, or its exonic Jaccard is below the single-exon cut FROZEN from the ideal arm's dev half (--freeze-single-cut writes it, the other
                      arms read it with --single-cut-file; --no-single-cut leaves single-exon extent unjudged)
  node_ok             none of the above
The genes with a node at all (recovered, whatever the class) are split into own-read nodes (the gene has primary reads on it) and shared nodes (no reads of its own: the node exists
because other genes' reads, or secondary alignments, reach it).

Writes OUT.per_gene.tsv (every listed gene), OUT.by_family.tsv, OUT.summary.json, OUT.report.txt (the text printed on stdout). Coordinates 0-based half-open. Standard library (+ pysam for --bam).
"""
import argparse
import bisect
import collections
import csv
import json
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BASELINE_DOC = os.path.join(HERE, "..", "SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md")
BIOTYPES = ["protein_coding", "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "lncRNA"]
ORDER = ["no_reads", "reads_not_on_gene", "no_locus", "outcompeted", "weak_overlap", "fused_excess", "extent_wrong", "node_ok"]
TRUTH_COLUMNS = ("simulable", "single_exon", "half", "source", "sd98", "s1c_families", "n_reads")


def die(msg):
    print(msg, file=sys.stderr)
    sys.exit(2)


# ------------------------------------------------------------------------------------------------- small helpers
def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def blocks_of(s):
    return [tuple(map(int, b.split("-"))) for b in s.split(",")] if s else []


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def bp(iv):
    return sum(e - s for s, e in iv)


def ov_bp(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def chain_of(exons):
    return frozenset((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))


def percentile(vals, p):
    """The p-th percentile by linear interpolation between order statistics (numpy's default), in exact integer arithmetic for the position."""
    v = sorted(vals)
    if not v:
        return None
    lo, rem = divmod(p * (len(v) - 1), 100)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (rem / 100.0)


def median(vals):
    return float(statistics.median(vals)) if vals else None


def family_key(f):
    m = re.search(r"(\d+)$", f)
    return (f[:m.start()] if m else f, int(m.group(1)) if m else 0, f)


def pct(a, b):
    return 100.0 * a / b if b else None


def fmt(x, nd=3):
    return "n/a" if x is None else f"{x:.{nd}f}"


def fmt_pct(x):
    return "n/a" if x is None else f"{x:.1f}%"


# ------------------------------------------------------------------------------------------------- inputs
def load_truth(prefix):
    """(universe, scored, not_simulable). universe[gene] = dict(chrom, strand, exons, bp, span, chains, single, n_reads) for every gene with an ELIGIBLE transcript (the ideal arm's universe);
    scored[gene] = the listed gene's targets row if it can be simulated; not_simulable = the other listed genes' rows."""
    tx_blocks, chains, n_reads, meta = collections.defaultdict(list), collections.defaultdict(set), collections.Counter(), {}
    rows = read_tsv(prefix + ".transcripts.tsv")
    if rows and "eligible" not in rows[0]:
        die(f"{prefix}.transcripts.tsv has no `eligible` column: it was not written by this version of sim_genes.py")
    for r in rows:
        n_reads[r["gene"]] += int(r["n_reads"])
        if not int(r["eligible"]):
            continue
        tx_blocks[r["gene"]].extend(blocks_of(r["canon_blocks"]))
        chains[r["gene"]].add(frozenset(tuple(map(int, x.split("-"))) for x in r["chain"].split(",")) if r["chain"] else frozenset())
        meta[r["gene"]] = (r["chrom"], r["strand"])
    universe = {}
    for g, bl in tx_blocks.items():
        ex = merge(bl)
        universe[g] = dict(chrom=meta[g][0], strand=meta[g][1], exons=ex, bp=bp(ex), span=(ex[0][0], ex[-1][1]), chains=chains[g], single=all(not c for c in chains[g]), n_reads=n_reads[g])
    targets = read_tsv(prefix + ".targets.tsv")
    if targets and any(c not in targets[0] for c in TRUTH_COLUMNS):
        die(f"{prefix}.targets.tsv lacks columns {[c for c in TRUTH_COLUMNS if c not in targets[0]]}: it was not written by this version of sim_genes.py")
    scored, not_simulable = collections.OrderedDict(), []
    for r in targets:
        if int(r["simulable"]) and r["gene"] in universe:
            scored[r["gene"]] = r
        else:
            not_simulable.append(r)
    return universe, scored, not_simulable


def parse_attr(text, key):
    m = re.search(key + r' "([^"]+)"', text)
    return m.group(1) if m else None


def load_units(path, unit, chroms):
    """uid -> dict(chrom, strand, exons, bp, span). unit `node`: PREFIX.fam.loci.gff3 (gene line + exon lines of the representative); `transcript` / `locus`: any GTF (transcript_id, gene_id).
    `chroms`: keep only these contigs (None = all)."""
    units = {}
    if unit == "node":
        span, strand, chrom, exons = {}, {}, {}, collections.defaultdict(list)
        with open(path) as fh:
            for ln in fh:
                if not ln or ln[0] == "#":
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) < 9 or (chroms is not None and f[0] not in chroms):
                    continue
                at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                if f[2] == "gene":
                    nid = at.get("Name") or at["ID"].replace("gene-", "", 1)
                    span[nid], strand[nid], chrom[nid] = (int(f[3]) - 1, int(f[4])), f[6], f[0]
                elif f[2] == "exon":
                    nid = at.get("gene") or at["Parent"].replace("gene-", "", 1)
                    exons[nid].append((int(f[3]) - 1, int(f[4])))
                    chrom.setdefault(nid, f[0])
                    strand.setdefault(nid, f[6])
        for nid, ex in exons.items():
            m = merge(ex)
            units[nid] = dict(chrom=chrom[nid], strand=strand[nid], exons=m, bp=bp(m), span=span.get(nid, (m[0][0], m[-1][1])))
        return units
    tx_span, tx_exons, tx_meta, tx_gene = {}, collections.defaultdict(list), {}, {}
    with open(path) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon") or (chroms is not None and f[0] not in chroms):
                continue
            t = parse_attr(f[8], "transcript_id")
            if t is None:
                continue
            tx_gene.setdefault(t, parse_attr(f[8], "gene_id") or t)
            tx_meta.setdefault(t, (f[0], f[6]))
            if f[2] == "transcript":
                tx_span[t] = (int(f[3]) - 1, int(f[4]))
            else:
                tx_exons[t].append((int(f[3]) - 1, int(f[4])))
    if unit == "transcript":
        for t, ex in tx_exons.items():
            m = merge(ex)
            units[t] = dict(chrom=tx_meta[t][0], strand=tx_meta[t][1], exons=m, bp=bp(m), span=tx_span.get(t, (m[0][0], m[-1][1])))
    else:
        grouped = collections.defaultdict(list)
        for t in sorted(tx_exons):
            grouped[tx_gene[t]].append(t)
        for g, ts in grouped.items():
            m = merge([b for t in ts for b in tx_exons[t]])
            spans = [tx_span.get(t, (min(s for s, _ in tx_exons[t]), max(e for _, e in tx_exons[t]))) for t in ts]
            units[g] = dict(chrom=tx_meta[ts[0]][0], strand=tx_meta[ts[0]][1], exons=m, bp=bp(m), span=(min(s for s, _ in spans), max(e for _, e in spans)))
    return units


def read_baseline(path):
    """The 71.0% reference line of the 2026-09-22 fidelity document (None if it is not there). It is a RefSeq-NAME join of 1,078 genes and an older pipeline: printed once, as such."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        text = fh.read()
    n = re.search(r"Soto genes carrying a Family ID:\s*\**([\d,]+)\**", text)
    m = re.search(r"\*\*(\d+)\s*=\s*([\d.]+)%\*\*", text)
    f = re.search(r"\*\*(\d+)\s*\(([\d.]+)%\)\s*sit on a locus", text)
    if not (n and m and f):
        return None
    return dict(n=int(n.group(1).replace(",", "")), matched=int(m.group(1)), matched_share=float(m.group(2)), fused=int(f.group(1)), fused_share=float(f.group(2)),
                source=f"{os.path.basename(path)} (section 6w5, 2026-09-22)")


def read_ids(path):
    """Gene ids from a table with a `Gene ID` / `gene_id` column, or a plain list (one id per line). A table with only name columns is refused: names are not keys."""
    with open(path, newline="") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    if lines and ("\t" in lines[0] or lines[0].strip() in ("Gene ID", "gene_id")):
        head = lines[0].split("\t")
        for col in ("Gene ID", "gene_id"):
            if col in head:
                i = head.index(col)
                return [ln.split("\t")[i].strip() for ln in lines[1:] if len(ln.split("\t")) > i]
        die(f"{path}: a subset is joined on gene ids ('Gene ID' or 'gene_id' column); this table has {head} (names are not keys: 1,774 names for 2,334 ids)")
    return [ln.split()[0] for ln in lines]


def read_reads_table(path, column):
    out = {}
    for r in read_tsv(path):
        g = r.get("gene_id") or r.get("Gene ID") or r.get("gene")
        if g is not None and column in r:
            out[g] = int(float(r[column]))
    return out


# ------------------------------------------------------------------------------------------------- matching
class SpanIndex:
    """Intervals sorted by start with a running maximum of ends: `overlapping(s, e)` yields the items whose [start, end) overlaps [s, e) by >= 1 bp."""

    def __init__(self, items):   # items: (start, end, key)
        self.items = sorted(items)
        self.starts = [it[0] for it in self.items]
        self.maxend = []
        m = -1
        for it in self.items:
            m = max(m, it[1])
            self.maxend.append(m)

    def overlapping(self, s, e):
        i = bisect.bisect_left(self.starts, e) - 1
        while i >= 0 and self.maxend[i] > s:
            if self.items[i][1] > s:
                yield self.items[i]
            i -= 1


def candidates(genes, units, compete):
    """[(overlap bp, jaccard, gene, unit)] for the competing genes whose exons overlap a unit's, and the number of overlapping units per gene."""
    by_chrom = collections.defaultdict(list)
    for u, d in units.items():
        by_chrom[d["chrom"]].append((d["span"][0], d["span"][1], u))
    idx = {c: SpanIndex(v) for c, v in by_chrom.items()}
    out, n_cand = [], collections.Counter()
    for g in sorted(compete):
        gd = genes[g]
        if gd["chrom"] not in idx:
            continue
        for _s, _e, u in idx[gd["chrom"]].overlapping(*gd["span"]):
            o = ov_bp(gd["exons"], units[u]["exons"])
            if o > 0:
                out.append((o, o / max(gd["bp"] + units[u]["bp"] - o, 1), g, u))
                n_cand[g] += 1
    return out, n_cand


def greedy(cands, reverse_ties=False):
    """One-to-one by exonic overlap descending, then Jaccard descending, then gene id and unit id (ascending, or descending with reverse_ties)."""
    order = sorted(cands, key=lambda c: c[3], reverse=reverse_ties)
    order.sort(key=lambda c: c[2], reverse=reverse_ties)
    order.sort(key=lambda c: (-c[0], -c[1]))
    taken_g, taken_u, match = set(), set(), {}
    for o, _j, g, u in order:
        if g in taken_g or u in taken_u:
            continue
        taken_g.add(g)
        taken_u.add(u)
        match[g] = u
    return match


def covers(unit, gd, rule):
    """Does the unit cover the gene? Span rule: the unit's span overlaps >= half of the gene's span. Exon rule: the unit's exons overlap >= half of the gene's exon bp."""
    if rule == "span":
        (us, ue), (s, e) = unit["span"], gd["span"]
        ov = min(ue, e) - max(us, s)
        return ov > 0 and 2 * ov >= (e - s)
    return gd["bp"] > 0 and 2 * ov_bp(unit["exons"], gd["exons"]) >= gd["bp"]


def covered(unit, universe, gidx, rule, skip):
    """The genes of the universe other than `skip` that the unit covers."""
    if unit["chrom"] not in gidx:
        return set()
    return {h for _s, _e, h in gidx[unit["chrom"]].overlapping(*unit["span"]) if h != skip and covers(unit, universe[h], rule)}


def bam_read_stats(bam, scored, universe):
    """per scored gene: the reads simulated from it that the BAM holds (primary or unmapped records) and the number with a PRIMARY alignment overlapping its exons. Reads are tied to their
    source by name `transcript|gene|fl|k`."""
    import pysam
    sim, on = collections.Counter(), collections.Counter()
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        for rec in fh.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            parts = rec.query_name.split("|")
            if len(parts) != 4 or parts[1] not in scored:
                continue
            g = parts[1]
            sim[g] += 1
            if rec.is_unmapped or rec.reference_name != universe[g]["chrom"]:
                continue
            ex = universe[g]["exons"]
            for s, e in rec.get_blocks():
                if any(x < e and y > s for x, y in ex):
                    on[g] += 1
                    break
    return sim, on


# ------------------------------------------------------------------------------------------------- scoring
def build_rows(universe, scored, units, compete_mode):
    compete = set(universe) if compete_mode == "all" else set(scored)
    cands, n_cand = candidates(universe, units, compete)
    match, rmatch = greedy(cands), greedy(cands, reverse_ties=True)

    def signature(g, m):   # a gene's outcome apart from the id of its unit
        u = m.get(g)
        return None if u is None else (ov_bp(universe[g]["exons"], units[u]["exons"]), units[u]["span"], units[u]["bp"], units[u]["strand"])

    tie_sensitive = sum(1 for g in scored if signature(g, match) != signature(g, rmatch))
    by_chrom = collections.defaultdict(list)
    for g, d in universe.items():
        by_chrom[d["chrom"]].append((d["span"][0], d["span"][1], g))
    gidx = {c: SpanIndex(v) for c, v in by_chrom.items()}
    rows = collections.OrderedDict()
    for g in scored:
        t, gd = scored[g], universe[g]
        row = dict(gene=g, gene_name=t.get("s1c_name") or t.get("gene_name") or "", biotype=t.get("s1c_biotype") or t.get("biotype") or "", source=t.get("source") or "other",
                   sd98=t.get("sd98", ""), families=t.get("s1c_families", ""), half=t.get("half", ""), single_exon=int(gd["single"]), n_reads=gd["n_reads"], n_candidates=n_cand[g], matched=0,
                   unit="", overlap_bp="", overlap_frac="", jaccard=None, d5=None, d3=None, width_ratio=None, strand_agree="", n_genes_span="", n_genes_exon="", cov_span="", cov_exon="",
                   own_cov_span="", own_cov_exon="", excess_span=0, excess_exon=0, chain_equal="", unit_spliced="", locus_exists=0)
        u = match.get(g)
        if u is not None:
            ud = units[u]
            o = ov_bp(gd["exons"], ud["exons"])
            gs, ge = gd["span"]
            us, ue = ud["span"]
            own = dict(chrom=gd["chrom"], span=gd["span"], exons=gd["exons"])
            cov = {r: covered(ud, universe, gidx, r, g) for r in ("span", "exon")}
            own_cov = {r: covered(own, universe, gidx, r, g) for r in ("span", "exon")}
            row.update(matched=1, unit=u, overlap_bp=o, overlap_frac=o / max(min(gd["bp"], ud["bp"]), 1), jaccard=o / max(gd["bp"] + ud["bp"] - o, 1), width_ratio=ud["bp"] / max(gd["bp"], 1),
                       d5=(us - gs) if gd["strand"] == "+" else (ge - ue), d3=(ge - ue) if gd["strand"] == "+" else (us - gs),
                       strand_agree=(int(ud["strand"] == gd["strand"]) if ud["strand"] in "+-" and gd["strand"] in "+-" else ""),
                       n_genes_span=len(cov["span"]) + int(covers(ud, gd, "span")), n_genes_exon=len(cov["exon"]) + int(covers(ud, gd, "exon")),
                       cov_span=len(cov["span"]), cov_exon=len(cov["exon"]), own_cov_span=len(own_cov["span"]), own_cov_exon=len(own_cov["exon"]),
                       excess_span=len(cov["span"] - own_cov["span"]), excess_exon=len(cov["exon"] - own_cov["exon"]),
                       chain_equal=int(chain_of(ud["exons"]) in gd["chains"]), unit_spliced=int(len(ud["exons"]) > 1), locus_exists=int(2 * o >= min(gd["bp"], ud["bp"])))
        rows[g] = row
    return rows, tie_sensitive


def judge_extent(row, cut):
    """True when the matched unit's extent is wrong. Multi-exon gene: the unit's chain is not an annotated isoform's. Single-exon gene: the unit is spliced, or its Jaccard is below the frozen cut."""
    if row["single_exon"]:
        return bool(row["unit_spliced"]) or (cut is not None and row["jaccard"] < cut)
    return not row["chain_equal"]


def classify(row, rule, cut, real, have_bam):
    if real:
        if row["reads_on_gene"] is not None and row["reads_on_gene"] == 0:
            return "no_reads"
    else:
        if row["n_reads"] == 0:
            return "no_reads"
        if have_bam and row["reads_on_gene"] == 0:
            return "reads_not_on_gene"
    if not row["matched"]:
        return "no_locus" if row["n_candidates"] == 0 else "outcompeted"
    if not row["locus_exists"]:
        return "weak_overlap"
    if row["excess_" + rule] > 0:
        return "fused_excess"
    return "extent_wrong" if judge_extent(row, cut) else "node_ok"


def freeze_cut(rows, rule, percentile_p):
    """The single-exon cut: the lowest p% of the exonic Jaccard of the dev half's single-exon genes that reach the extent test (matched, non-weak, no excess fusion, unspliced node, with reads)."""
    vals = [r["jaccard"] for r in rows.values() if r["single_exon"] and r["half"] == "dev" and r["locus_exists"] and r["excess_" + rule] == 0 and not r["unit_spliced"]
            and r["n_reads"] > 0 and (r["reads_on_gene"] is None or r["reads_on_gene"] > 0)]
    return percentile(vals, percentile_p), len(vals)


def counts(rows, key="class"):
    c = collections.Counter(r[key] for r in rows)
    return {k: c[k] for k in ORDER}


def biotype_group(bt):
    return bt if bt in BIOTYPES else "other"


def group_classes(rows, keyfn, keys=None):
    groups = collections.OrderedDict((k, []) for k in keys) if keys else collections.OrderedDict()
    for r in rows:
        groups.setdefault(keyfn(r), []).append(r)
    return collections.OrderedDict((k, counts(v)) for k, v in groups.items())


def descriptors(rows):
    return dict(n=len(rows), median_jaccard=median([r["jaccard"] for r in rows]), median_abs_d5=median([abs(r["d5"]) for r in rows]), median_abs_d3=median([abs(r["d3"]) for r in rows]),
                median_width_ratio=median([r["width_ratio"] for r in rows]))


def family_rows(rows):
    fam = collections.defaultdict(list)
    for r in rows:
        for f in (r["families"].split(";") if r["families"] else ["(none)"]):
            fam[f].append(r)
    out = []
    for f in sorted(fam, key=family_key):
        rs = fam[f]
        c = counts(rs)
        out.append(dict(family=f, n_scored=len(rs), **c, all_node_ok=int(c["node_ok"] == len(rs))))
    return out


# ------------------------------------------------------------------------------------------------- output
PER_GENE_COLUMNS = ["gene", "gene_name", "biotype", "source", "sd98", "families", "half", "single_exon", "n_reads", "reads_on_gene", "n_candidates", "matched", "unit", "overlap_bp",
                    "overlap_frac", "jaccard", "d5", "d3", "width_ratio", "strand_agree", "n_genes_span", "n_genes_exon", "cov_span", "cov_exon", "own_cov_span", "own_cov_exon",
                    "excess_span", "excess_exon", "chain_equal", "unit_spliced", "locus_exists", "support", "class_span", "class_exon", "class", "reason"]


def write_per_gene(path, rows, not_sim):
    by_gene = {t["gene"]: t for t in not_sim}
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PER_GENE_COLUMNS)
        for g in sorted(set(rows) | set(by_gene)):
            if g in rows:
                r = rows[g]
                w.writerow([("%.4f" % r[c]) if c in ("jaccard", "width_ratio", "overlap_frac") and r.get(c) not in (None, "") else ("" if r.get(c) is None else r.get(c, "")) for c in PER_GENE_COLUMNS])
            else:
                t = by_gene[g]
                row = dict(gene=g, gene_name=t.get("s1c_name", ""), biotype=t.get("s1c_biotype") or t.get("biotype", ""), source=t.get("source", ""), sd98=t.get("sd98", ""),
                           families=t.get("s1c_families", ""), half=t.get("half", ""), class_span="not_simulable", class_exon="not_simulable", reason=t.get("reason", ""))
                row["class"] = "not_simulable"
                w.writerow([row.get(c, "") for c in PER_GENE_COLUMNS])


def write_by_family(path, fam):
    cols = ["family", "n_scored"] + ORDER + ["all_node_ok"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in fam:
            w.writerow([r[c] for c in cols])


def class_table(title, ctx, groups):
    lines = [f"TABLE {title}  {ctx}", " | ".join(["class"] + list(groups))]
    for k in ORDER:
        lines.append(" | ".join([k] + [str(groups[g][k]) for g in groups]))
    lines.append(" | ".join(["n"] + [str(sum(groups[g].values())) for g in groups]))
    lines.append(" | ".join(["node_ok share"] + [fmt_pct(pct(groups[g]["node_ok"], sum(groups[g].values()))) for g in groups]))
    return lines


def template(summ, ctx):
    c, pop = summ["classes"], summ["population"]
    n = pop["simulable"]
    return (f"{ctx} Of {n} simulable Soto genes (by gene id), {c['node_ok']} ({fmt_pct(pct(c['node_ok'], n))}) are node_ok; {c['no_locus']} have no node on their exons "
            f"({summ['classes_by_stratum']['single_exon']['no_locus']} of them single-exon), {c['outcompeted']} lost their node to a gene that overlaps it more, {c['weak_overlap']} share too few exon bp "
            f"with their node, {c['fused_excess']} sit in a node that covers a gene their own annotation does not, {c['extent_wrong']} have another chain or extent; "
            f"{c['no_reads']} have no reads and {c['reads_not_on_gene']} have reads that never align on them. {pop['not_simulable']['genes']} listed genes cannot be simulated "
            f"({pop['not_simulable']['sd98_eligible']} in the SD98 set) and are outside the denominator.")


def build_report(summ, ctx, args):
    L = [f"score_soto_nodes.py  truth={args.truth}  units={args.loci_gff3 if args.unit == 'node' else args.gtf}", f"context {ctx}"]
    pop = summ["population"]
    ns = pop["not_simulable"]
    L.append(f"tie policy: units are taken by exonic overlap bp (desc), then exonic Jaccard (desc), then gene id and unit id (asc); tie sensitive genes (outcome changes under the reverse order): {summ['tie']['sensitive_genes']}")
    L.append("")
    L.append(f"TABLE 1  population  {ctx}")
    L.append(f"listed genes {pop['listed']}; simulable (the denominator, by gene id) {pop['simulable']}; not simulable {ns['genes']} ({ns['sd98_eligible']} in the SD98 set): "
             + ", ".join(f"{k}: {v}" for k, v in sorted(ns["by_reason"].items())))
    L.append("simulable by annotation source: " + ", ".join(f"{k} {v}" for k, v in pop["by_source"].items()) + f"; single-exon genes {summ['classes_by_stratum']['single_exon_n']}")
    L.append("")
    L += class_table("2  ordered classes, all genes and by exon count", ctx, collections.OrderedDict(
        [("ALL", summ["classes"]), ("single-exon", summ["classes_by_stratum"]["single_exon"]), ("multi-exon", summ["classes_by_stratum"]["multi_exon"])]))
    L.append("")
    L += class_table("3  ordered classes by biotype", ctx, summ["classes_by_biotype"])
    L.append("")
    other_rule = "exon" if args.fused_rule == "span" else "span"
    L += class_table(f"4  the same classes under the other fusion rule ({other_rule}); the classes above use {args.fused_rule}", ctx,
                     collections.OrderedDict([(f"rule={args.fused_rule}", summ["classes"]), (f"rule={other_rule}", summ["classes_by_rule"][other_rule])]))
    L.append("")
    fu = summ["fusion"]
    L.append(f"TABLE 5  fusion among the {fu['matched']} matched genes: the annotation's own overlap versus the excess  {ctx}")
    L.append(f"span rule: raw fused {fu['raw_span']} ({fmt_pct(pct(fu['raw_span'], fu['matched']))}), excess {fu['excess_span']}; "
             f"exon rule: raw fused {fu['raw_exon']} ({fmt_pct(pct(fu['raw_exon'], fu['matched']))}), excess {fu['excess_exon']}")
    L.append("")
    rc = summ["recovered"]
    L.append(f"TABLE 6  nodes with a real overlap, own-read versus shared  {ctx}")
    L.append(f"genes with a node (matched, overlap >= half of the smaller exon union): {rc['nodes_with_overlap']}; own-read nodes (the gene has primary reads on it) {rc['own_read']}; "
             f"shared nodes (no reads of its own on it: the node exists through other genes' reads or secondary alignments) {rc['shared_only']}; read provenance unknown {rc['unknown']}")
    nr = summ["no_reads"]
    L.append(f"genes without reads in the arm: {nr['genes']}, of which {nr['with_node']} have a node anyway (shared)")
    L.append("")
    L.append(f"TABLE 7  descriptors of nodes (conditioned on a node existing: descriptors, not estimands)  {ctx}")
    for k in ("node_ok", "recovered"):
        d = summ["descriptors"][k]
        L.append(f"{k}: n={d['n']}, median Jaccard {fmt(d['median_jaccard'])}, median |d5| {fmt(d['median_abs_d5'], 0)}, median |d3| {fmt(d['median_abs_d3'], 0)}, median width ratio {fmt(d['median_width_ratio'], 2)}")
    L.append("")
    L += class_table("8  by annotation source", ctx, summ["classes_by_source"])
    L.append("")
    L += class_table("9  by half of the Soto split (dev / heldout)", ctx, summ["classes_by_half"])
    sc = summ["single_cut"]
    if sc is None:
        L.append("single-exon extent: not judged (--no-single-cut): a single-exon gene is extent_wrong only if its node is spliced")
    elif sc.get("frozen_now"):
        L.append(f"single-exon extent: this run FROZE the cut {sc['value']:.4f} ({sc['rule']}; n={sc['n']}); by construction about {sc['percentile']}% of this run's dev-half single-exon genes lie below it, read the held-out half")
    else:
        L.append(f"single-exon extent: frozen cut {sc['value']:.4f} from {os.path.basename(sc['file'])} ({sc['rule']}; n={sc['n']}); the same cut in every arm")
    fam = summ["families"]
    L.append("")
    L.append(f"TABLE 10  by Soto family  {ctx}")
    L.append(f"{fam['n']} families with a scored member; {fam['all_node_ok']} have every scored member node_ok; {fam['none_node_ok']} have none (per-family rows: by_family.tsv)")
    for label, sub in summ["subsets"].items():
        L.append("")
        L += class_table(f"11  subset '{label}' ({sub['n']} scored of {sub['n_ids']} ids)", ctx, collections.OrderedDict([("subset", sub["classes"])]))
    if summ["comparators"]:
        L.append("")
        L.append("COMPARATORS  (each row is its own scorer run; the unit, the competition scope and the fusion rule are those printed)")
        L.append("label | context | n | " + " | ".join(ORDER))
        rows = [dict(label="this arm", context=ctx, n=pop["simulable"], classes=summ["classes"])] + summ["comparators"]
        for c in rows:
            L.append(f"{c['label']} | {c['context']} | n={c['n']} | " + " | ".join(f"{k} {fmt_pct(pct(c['classes'][k], c['n']))}" for k in ORDER))
    if summ["assembly_polish"]:
        L.append("")
        L.append("ASSEMBLY POLISH (from the assemble log):")
        L += ["  " + x for x in summ["assembly_polish"]]
    ref = summ["reference_line_not_like_for_like"]
    L.append("")
    if ref:
        L.append(f"reference line, NOT like-for-like (a RefSeq-name join of {ref['n']} genes, an older pipeline, another library): {ref['matched']} of {ref['n']} matched = {ref['matched_share']:.1f}% ({ref['source']}). "
                 f"The like-for-like comparators are the real-library products scored by this same scorer (COMPARATORS).")
    L.append("")
    L.append("INTERPRETATION TEMPLATE: " + summ["interpretation_template"])
    return "\n".join(L) + "\n"


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, help="sim_genes.py prefix (targets.tsv, transcripts.tsv): the ideal arm's truth, or the arm's own")
    ap.add_argument("--out", required=True)
    ap.add_argument("--unit", choices=["node", "transcript", "locus"], default="node")
    ap.add_argument("--loci-gff3", help="the families stage's PREFIX.fam.loci.gff3 (--unit node)")
    ap.add_argument("--gtf", help="any assembled GTF (--unit transcript or locus)")
    ap.add_argument("--compete", choices=["all", "soto"], default="all", help="genes that compete for units: every gene of the universe, or the listed genes only")
    ap.add_argument("--fused-rule", choices=["span", "exon"], default="span", help="the rule that decides the fused_excess class (both are always computed)")
    ap.add_argument("--freeze-single-cut", metavar="JSON", help="calibration arm: compute the single-exon Jaccard cut from this run's dev half, write it here and apply it")
    ap.add_argument("--cut-percentile", type=int, default=5, help="with --freeze-single-cut: the percentile of the dev-half Jaccards that the cut is")
    ap.add_argument("--single-cut-file", metavar="JSON", help="read the frozen cut")
    ap.add_argument("--no-single-cut", action="store_true", help="leave single-exon extent unjudged")
    ap.add_argument("--bam", help="the arm's BAM: reads_not_on_gene and own-read versus shared nodes (simulated arms only)")
    ap.add_argument("--real-library", action="store_true", help="the arm is a real library: the truth tables' reads are not its reads")
    ap.add_argument("--reads-per-gene", help="real library: TSV gene_id, <column>: reads with a primary alignment on the gene")
    ap.add_argument("--reads-column", default="prim_on_exons")
    ap.add_argument("--assemble-log", help="the driver's PREFIX.assemble.log: its ASSEMBLY POLISH lines are printed")
    ap.add_argument("--subset", action="append", default=[], help="a file of gene ids (a `Gene ID` column or one id per line) to report apart")
    ap.add_argument("--subset-label", action="append", default=[])
    ap.add_argument("--compare", action="append", default=[], metavar="LABEL=SUMMARY.json", help="another run's summary to print beside this one")
    ap.add_argument("--baseline-doc", default=DEFAULT_BASELINE_DOC)
    a = ap.parse_args(argv)
    if a.unit == "node" and not a.loci_gff3:
        ap.error("--unit node needs --loci-gff3 (the families stage's PREFIX.fam.loci.gff3)")
    if a.unit != "node" and not a.gtf:
        ap.error(f"--unit {a.unit} needs --gtf")
    if sum(map(bool, (a.freeze_single_cut, a.single_cut_file, a.no_single_cut))) != 1:
        ap.error("give exactly one of --single-cut-file CUT.json (a frozen cut), --freeze-single-cut CUT.json (the calibration arm) or --no-single-cut")
    if a.freeze_single_cut and os.path.exists(a.freeze_single_cut):
        ap.error(f"{a.freeze_single_cut} already exists: a frozen cut is not recomputed (delete it deliberately, or read it with --single-cut-file)")
    if a.real_library and a.freeze_single_cut:
        ap.error("the single-exon cut is frozen on the simulated ideal arm, not on a real library")
    if a.reads_per_gene and not a.real_library:
        ap.error("--reads-per-gene is for --real-library")
    if len(a.subset_label) > len(a.subset):
        ap.error("more --subset-label than --subset")
    return a


def main(argv=None):
    a = parse_args(argv)
    universe, scored, not_sim = load_truth(a.truth)
    units = load_units(a.loci_gff3 if a.unit == "node" else a.gtf, a.unit, {d["chrom"] for d in universe.values()})
    rows, tie_sensitive = build_rows(universe, scored, units, a.compete)
    # reads: a simulated arm knows its own (the truth tables) and the BAM says where they align; a real library has a table or nothing
    have_bam = bool(a.bam) and not a.real_library
    if a.real_library:
        table = read_reads_table(a.reads_per_gene, a.reads_column) if a.reads_per_gene else {}
        for g, r in rows.items():
            r["n_reads"] = r["reads_on_gene"] = table.get(g)
    else:
        _sim, on_bam = bam_read_stats(a.bam, scored, universe) if have_bam else ({}, {})
        for g, r in rows.items():
            r["reads_on_gene"] = (0 if r["n_reads"] == 0 else on_bam[g]) if (have_bam or r["n_reads"] == 0) else None
    cut, single_cut = None, None
    if a.single_cut_file:
        with open(a.single_cut_file) as fh:
            single_cut = dict(json.load(fh), frozen_now=False, file=a.single_cut_file)
        cut = single_cut["value"]
    elif a.freeze_single_cut:
        cut, n = freeze_cut(rows, a.fused_rule, a.cut_percentile)
        if cut is None or n < 2:
            die("--freeze-single-cut: fewer than two dev-half single-exon genes reach the extent test in this run, nothing to calibrate on")
        payload = dict(value=cut, percentile=a.cut_percentile, n=n, half="dev",
                       rule=f"the lowest {a.cut_percentile}% of the exonic Jaccard of the dev half's single-exon genes matched one-to-one (non-weak, no excess fusion, unspliced node, with reads) in the "
                            f"calibration arm; unit={a.unit}, compete={a.compete}, fused-rule={a.fused_rule}, truth={os.path.basename(a.truth)}")
        with open(a.freeze_single_cut, "w") as fh:
            json.dump(payload, fh, indent=1)
        single_cut = dict(payload, frozen_now=True, file=a.freeze_single_cut)
    for r in rows.values():
        for rule in ("span", "exon"):
            r["class_" + rule] = classify(r, rule, cut, a.real_library, have_bam)
        r["class"] = r["class_" + a.fused_rule]
        r["support"] = "" if not r["locus_exists"] else ("unknown" if r["reads_on_gene"] is None else ("own_read" if r["reads_on_gene"] > 0 else "shared_only"))
        r["reason"] = ""
    lst = list(rows.values())
    ctx = f"[unit={a.unit} compete={a.compete} fused-rule={a.fused_rule} match=two-sided>=0.5 arm={'real' if a.real_library else 'sim'}]"
    ns_reason = collections.Counter(re.sub(r"\s*\(.*\)\s*$", "", t.get("reason", "")) or "unspecified" for t in not_sim)
    pop = dict(listed=len(scored) + len(not_sim), simulable=len(scored),
               not_simulable=dict(genes=len(not_sim), sd98_eligible=sum(1 for t in not_sim if t.get("sd98") == "Yes"), by_reason=dict(ns_reason)),
               by_source=dict(collections.Counter(r["source"] for r in lst)))
    matched = [r for r in lst if r["matched"]]
    summ = dict(tool="score_soto_nodes.py",
                inputs=dict(truth=a.truth, unit=a.unit, loci_gff3=a.loci_gff3, gtf=a.gtf, compete=a.compete, fused_rule=a.fused_rule, real_library=a.real_library, bam=a.bam,
                            baseline_doc=os.path.basename(a.baseline_doc)),
                context=ctx, population=pop, n_universe_genes=len(universe), n_units=len(units), classes=counts(lst),
                classes_by_rule={rule: counts(lst, "class_" + rule) for rule in ("span", "exon")},
                classes_by_stratum=dict(single_exon=counts([r for r in lst if r["single_exon"]]), multi_exon=counts([r for r in lst if not r["single_exon"]]),
                                        single_exon_n=sum(1 for r in lst if r["single_exon"])),
                classes_by_biotype=group_classes(lst, lambda r: biotype_group(r["biotype"]), BIOTYPES + ["other"]),
                classes_by_source=group_classes(lst, lambda r: r["source"], ["CAT", "Liftoff"]),
                classes_by_half=group_classes([r for r in lst if r["half"]], lambda r: r["half"], ["dev", "heldout"]),
                fusion=dict(matched=len(matched), raw_span=sum(1 for r in matched if r["n_genes_span"] >= 2), raw_exon=sum(1 for r in matched if r["n_genes_exon"] >= 2),
                            excess_span=sum(1 for r in matched if r["excess_span"] > 0), excess_exon=sum(1 for r in matched if r["excess_exon"] > 0)),
                recovered=dict(nodes_with_overlap=sum(r["locus_exists"] for r in lst), **{k: sum(1 for r in lst if r["support"] == k) for k in ("own_read", "shared_only", "unknown")}),
                no_reads=dict(genes=sum(1 for r in lst if r["class"] == "no_reads"), with_node=sum(1 for r in lst if r["class"] == "no_reads" and r["locus_exists"])),
                descriptors=dict(node_ok=descriptors([r for r in lst if r["class"] == "node_ok"]), recovered=descriptors([r for r in lst if r["locus_exists"]])),
                single_cut=single_cut, tie=dict(sensitive_genes=tie_sensitive, policy="exonic overlap bp desc, exonic Jaccard desc, gene id asc, unit id asc"), comparators=[], subsets={},
                assembly_polish=[], reference_line_not_like_for_like=read_baseline(a.baseline_doc))
    fam = family_rows(lst)
    summ["families"] = dict(n=len(fam), all_node_ok=sum(r["all_node_ok"] for r in fam), none_node_ok=sum(1 for r in fam if r["node_ok"] == 0))
    for i, path in enumerate(a.subset):
        label = a.subset_label[i] if i < len(a.subset_label) else os.path.basename(path)
        ids = set(read_ids(path))
        sub = [r for r in lst if r["gene"] in ids]
        summ["subsets"][label] = dict(n_ids=len(ids), n=len(sub), classes=counts(sub))
    for spec in a.compare:
        label, _, path = spec.partition("=")
        with open(path) as fh:
            other = json.load(fh)
        summ["comparators"].append(dict(label=label, context=other.get("context", ""), n=other["population"]["simulable"], classes=other["classes"], file=path))
    if a.assemble_log:
        with open(a.assemble_log, encoding="utf-8", errors="replace") as fh:
            summ["assembly_polish"] = [ln.rstrip("\n") for ln in fh if "ASSEMBLY POLISH" in ln]
    summ["interpretation_template"] = template(summ, ctx)
    write_per_gene(a.out + ".per_gene.tsv", rows, not_sim)
    write_by_family(a.out + ".by_family.tsv", fam)
    with open(a.out + ".summary.json", "w") as fh:
        json.dump(summ, fh, indent=1)
    report = build_report(summ, ctx, a)
    with open(a.out + ".report.txt", "w") as fh:
        fh.write(report)
    sys.stdout.write(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
