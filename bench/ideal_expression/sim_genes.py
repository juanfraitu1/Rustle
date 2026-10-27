#!/usr/bin/env python3
"""Ideal-arm read simulation for a GENE LIST: sim_windows.py's recipe applied to the 2,334 Soto et al. 2025 genes (docs/PREREG_ideal_arm_pilot_2026-10-07.md).

    sim_genes.py --out PREFIX [--genes FILE] [--pad 25000] [--reads-per-tx 10] [--seed 20261008] [--err 0.001] [--jitter 30] [--procs 1] [--no-reads]
                 [--weights FILE --total-reads N [--neighbour-weight median|median-all|X] [--missing-weight 0]]
    sim_genes.py gates   --out PREFIX --genome FA [--genes FILE]    G1 manifest (+ canonical introns re-read from the genome) and G2 read ends; exit 1 on FAIL
    sim_genes.py compare PREFIX_A PREFIX_B                          G3 (simulation half): the two simulations are byte-identical file by file
    sim_genes.py bamqc   --out PREFIX --bam BAM                     G4 mapping sanity from the BAM alone (one primary-or-unmapped record per read; shares printed, no bar)
    sim_genes.py depth-subset SHALLOW DEEP                          G-depth: the shallow arm is the first reads of the deep one, read for read (map the deep pool once, cut the rest)
    sim_genes.py readnames --out PREFIX --to FILE                   the name of every read, for `samtools view -N`
    sim_genes.py control --out PREFIX --to CTRL [--shift N]         the annotation as an assembled GTF and as one node per gene (positive control; --shift 50000: negative control)

Genes: the unique `Gene ID` values of bench/soto/soto_famCN_S1C.tsv (2,334 CAT/Liftoff v2.0 gene ids) or a plain list (one id per line). Windows: every listed gene's span +- pad, merged
per contig and clipped to the contig. Transcripts: every transcript of every gene (any biotype) whose span overlaps a window, so expressed neighbours are present as in the earlier ideal
runs. Canonicalisation, the skip rules (spliced length < 120 bp, molecule > 30 kb, a kept intron > 200 kb) and the read model are those of bench/ideal_expression/sim_windows.py and
bench/sim.py (`sim_windows.canon`, `sim_windows.load_cat`, `sim.simulate_reads`: error .001, indel error/3, no truncation, 0 to --jitter bp trimmed from each end, read k of
transcript t drawn from `sim.stable_seed(seed, t, k, ...)`), imported, not copied; the two decisions sim_windows.main inlines (the skip rule and the read loop) are restated here with
sim_windows' own constants and are pinned to it by the parity test of test_soto_pilot.py. Reads are named `transcript|gene|fl|k`: the truth label is the SOURCE copy, whatever the aligner does.

IDEAL ARM (default): --reads-per-tx reads for every simulable transcript. WEIGHTED ARM (--weights FILE --total-reads N): FILE is a TSV `gene_id<TAB>weight`; N reads are allocated over
the genes that have a simulable transcript and a positive weight, proportionally to weight (exact largest-remainder rounding, ties to the smaller gene id), and a gene's reads are split
equally over its simulable transcripts (the remainder to the first transcripts in transcript-id order). A listed gene absent from FILE gets --missing-weight (default 0); every non-listed
neighbour gets the one fixed --neighbour-weight (default: the median of the positive weights of the listed genes). A gene with weight 0 gets no reads and is reported as skipped.

Writes PREFIX.fq, PREFIX.transcripts.tsv (one row per transcript of the selected genes), PREFIX.targets.tsv (one row per listed gene: status simulated|skipped and the reason),
PREFIX.windows.tsv, PREFIX.skipped.tsv (the transcripts the skip rules reject), PREFIX.sim.json. Coordinates 0-based half-open; an intron is (donor_end, acceptor_start).
Output depends only on the arguments: no hash order, clock, path or process count enters any file.
"""
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import multiprocessing
import os
import random
import re
import statistics
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sim_windows as sw  # noqa: E402  (puts bench/ on sys.path: lib, sim)
import lib  # noqa: E402
import sim  # noqa: E402
import pysam  # noqa: E402

DEFAULT_GENES = os.path.join(HERE, "..", "soto", "soto_famCN_S1C.tsv")
DEFAULT_SPLIT = os.path.join(HERE, "..", "soto", "soto_split_2026-09-29.tsv")
SEED = 20261008   # the pilot's replicate: 20261006 and 20261007 are the two replicates of the NPIP / TBC1D3 ideal runs, whose reads it must not reuse
MIN_SPLICED = 120   # sim_windows.main inlines this literal
SUFFIXES = (".fq", ".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv", ".identical_seq.tsv", ".sim.json")
CHUNK_BP = 2_000_000   # read bases per work unit (the unit is also the unit of the process pool; the output does not depend on it)
ANCHOR = 20   # G2: bases of a read end matched against the transcript
TX_COLUMNS = ["transcript", "gene", "gene_name", "biotype", "is_soto", "chrom", "strand", "raw_blocks", "canon_blocks", "chain", "n_raw_introns", "n_noncanonical", "n_short",
              "spliced_len", "seq_sha1", "eligible", "simulated", "skip_reason", "n_reads"]
TARGET_COLUMNS = ["gene", "gene_name", "s1c_name", "biotype", "s1c_biotype", "s1c_families", "n_families", "half", "source", "sd98", "chrom", "strand", "start0", "end",
                  "n_tx_annotated", "n_tx_eligible", "n_tx_simulated", "n_reads", "exon_bp", "n_exons", "single_exon", "simulable", "identical_seq_genes", "status", "reason"]


# ------------------------------------------------------------------------------------------------- inputs
def read_gene_list(path):
    """Unique gene ids in first-seen order: the `Gene ID` column of an S1C-format table, else one id per line ('#' lines and blanks ignored)."""
    with open(path, newline="") as fh:
        lines = fh.read().splitlines()
    body = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    ids = []
    if body and "Gene ID" in body[0].split("\t"):
        col = body[0].split("\t").index("Gene ID")
        ids = [ln.split("\t")[col].strip() for ln in body[1:] if len(ln.split("\t")) > col]
    else:
        ids = [ln.split()[0] for ln in body]
    seen, out = set(), []
    for g in ids:
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def source_of(gene_id):
    """The annotation a gene id comes from: CAT (CHM13_G...) or Liftoff (LOFF_G...); '' for any other id."""
    return "CAT" if gene_id.startswith("CHM13_") else "Liftoff" if gene_id.startswith("LOFF_") else ""


def family_key(f):
    m = re.search(r"(\d+)$", f)
    return (f[:m.start()] if m else f, int(m.group(1)) if m else 0, f)


def load_s1c(path):
    """gene id -> dict(name, biotype, sd98, families): the S1C columns `Gene Name`, `Biotype`, `In Table S1 (SD98 gene set)`, `Family ID` (a gene can sit in several families; ids sorted naturally)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = (r.get("Gene ID") or "").strip()
            if not g:
                continue
            d = out.setdefault(g, dict(name=(r.get("Gene Name") or "").strip(), biotype=(r.get("Biotype") or "").strip(), sd98=(r.get("In Table S1 (SD98 gene set)") or "").strip(), families=set()))
            f = (r.get("Family ID") or "").strip()
            if f and f != "N/A":
                d["families"].add(f)
    for d in out.values():
        d["families"] = sorted(d["families"], key=family_key)
    return out


def load_split(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["gene_id"]] = r["half"]
    return out


def read_weights(path):
    """gene id -> weight from a TSV `gene_id<TAB>weight` (an optional header line, '#' comments); a repeated gene id, a negative or a non-finite weight is an error."""
    out, first = {}, True
    with open(path, newline="") as fh:
        for n, ln in enumerate(fh, 1):
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 2:
                f = ln.split()
            is_first, first = first, False
            try:
                w = float(f[1])
            except (ValueError, IndexError):
                if is_first and len(f) >= 2:
                    continue   # a header: the first content line, two or more fields, the second not a number
                raise ValueError(f"{path} line {n}: no numeric weight in {ln.strip()!r}")
            if not math.isfinite(w) or w < 0:
                raise ValueError(f"{path} line {n}: weight {f[1]!r} is negative or not finite")
            if f[0] in out:
                raise ValueError(f"{path} line {n}: gene {f[0]} appears twice")
            out[f[0]] = w
    return out


# ------------------------------------------------------------------------------------------------- the plan
def merge_windows(spans, pad, limit=None):
    """Spans +- pad, merged (touching windows merge, as in sim_windows), clipped to [0, limit]."""
    out = []
    for s, e in sorted(spans):
        lo, hi = max(0, s - pad), e + pad
        if limit is not None:
            hi = min(hi, limit)
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [tuple(w) for w in out]


def skip_reason(spliced, chain):
    """The simulate/skip rule of sim_windows.main, with its constants and its reason strings ('' = simulable)."""
    if spliced < MIN_SPLICED:
        return f"under {MIN_SPLICED} bp"
    if spliced > sw.MAX_MOLECULE:
        return f"molecule {spliced} bp > {sw.MAX_MOLECULE}"
    if any(x - d > sw.MAX_INTRON for d, x in chain):
        return "kept intron > 200 kb"
    return ""


def reason_class(reason):
    return f"molecule > {sw.MAX_MOLECULE} bp" if reason.startswith("molecule ") else reason


def allocate_weighted(weights, total):
    """Largest-remainder allocation of `total` reads over genes proportionally to `weights` (exact rational arithmetic; remainder ties go to the smaller gene id)."""
    items = sorted(weights.items())
    for g, w in items:
        if not math.isfinite(w) or w < 0:
            raise ValueError(f"weight of {g} is {w!r}: weights must be finite and >= 0")
    fr = {g: Fraction(w) for g, w in items}
    s = sum(fr.values())
    if s == 0:
        raise ValueError("every weight is zero: nothing to allocate")
    quota = {g: f * total / s for g, f in fr.items()}
    base = {g: q.numerator // q.denominator for g, q in quota.items()}
    left = total - sum(base.values())
    for g in sorted(quota, key=lambda g: (-(quota[g] - base[g]), g))[:left]:
        base[g] += 1
    return base


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def split_equal(n, k):
    return [n // k + (1 if i < n % k else 0) for i in range(k)]


def blocks_of(s):
    return [tuple(map(int, b.split("-"))) for b in s.split(",")] if s else []


def merge_blocks(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def plan(a, listed, fa):
    """Windows, selected genes and one row per transcript of the selected genes, contig by contig (FASTA order). Returns (rows, gene_info, windows, n_genes_in_windows, found)."""
    sw.CAT_GFF = a.gff   # sim_windows.load_cat reads this module global
    listed_set = set(listed)
    tf = pysam.TabixFile(a.gff)
    gff_contigs = list(tf.contigs)
    tf.close()
    order = [c for c in fa.references if c in gff_contigs]
    rows, gene_info, windows, found, n_sel = [], {}, [], set(), 0
    for chrom in order:
        genes, tx = sw.load_cat(chrom)
        mine = sorted(g for g in genes if g in listed_set)
        if not mine:
            continue
        found.update(mine)
        wins = merge_windows([genes[g]["span"] for g in mine], a.pad, fa.get_reference_length(chrom))
        windows += [(chrom, lo, hi) for lo, hi in wins]
        sel = {g for g, v in genes.items() if any(v["span"][0] < w[1] and v["span"][1] > w[0] for w in wins)}
        n_sel += len(sel)
        n_tx = collections.Counter(g for g, _st, _bl in tx.values() if g in sel)
        for g in sorted(sel):
            v = genes[g]
            gene_info[g] = dict(chrom=chrom, strand=v["strand"], name=v["name"], biotype=v["biotype"], start0=v["span"][0], end=v["span"][1], n_tx=n_tx[g])
        for t, (g, st, blocks) in sorted(tx.items(), key=lambda kv: (kv[1][2][0][0], kv[0])):
            if g not in sel:
                continue
            cb, chain, n_nc, n_short = sw.canon(fa, chrom, st, blocks)
            spliced = sum(e - s for s, e in cb)
            reason = skip_reason(spliced, chain)
            seq_sha1 = "" if reason else hashlib.sha1(lib.spliced1(fa, chrom, [(s + 1, e) for s, e in cb], st).encode("ascii")).hexdigest()
            rows.append(dict(transcript=t, gene=g, gene_name=genes[g]["name"], biotype=genes[g]["biotype"], is_soto=int(g in listed_set), chrom=chrom, strand=st,
                             raw_blocks=",".join(f"{s}-{e}" for s, e in blocks), canon_blocks=",".join(f"{s}-{e}" for s, e in cb), chain=",".join(f"{d}-{x}" for d, x in chain),
                             n_raw_introns=len(blocks) - 1, n_noncanonical=n_nc, n_short=n_short, spliced_len=spliced, seq_sha1=seq_sha1, eligible=int(not reason), simulated=0,
                             skip_reason=reason, n_reads=0))
    for chrom in gff_contigs:
        if chrom not in fa.references:
            stray = sorted(g for g in sw.load_cat(chrom)[0] if g in listed_set)
            if stray:
                raise SystemExit(f"{len(stray)} listed gene(s) ({stray[0]}, ...) lie on GFF contig {chrom}, which is not in the genome FASTA")
    return rows, gene_info, windows, n_sel, found


def assign_reads(rows, a, listed_set, weights):
    """Fill n_reads / simulated per row. Ideal arm: reads-per-tx for every eligible row. Weighted arm: largest-remainder over genes, equal split over a gene's eligible rows."""
    info = None
    elig = collections.defaultdict(list)
    for r in rows:
        if r["eligible"]:
            elig[r["gene"]].append(r)
    if weights is None:
        for r in rows:
            if r["eligible"]:
                r["n_reads"], r["simulated"] = a.reads_per_tx, 1
        return info
    soto_w = {g: weights[g] for g in sorted(listed_set) if g in weights}
    pos = [w for w in soto_w.values() if w > 0]
    mode = a.neighbour_weight
    if mode == "median":
        if not pos:
            raise SystemExit("--neighbour-weight median: no listed gene has a positive weight")
        nb, rule = float(statistics.median(pos)), "median of the positive Soto weights"
    elif mode == "median-all":
        if not soto_w:
            raise SystemExit("--neighbour-weight median-all: no listed gene appears in the weights file")
        nb, rule = float(statistics.median(soto_w.values())), "median of the Soto weights in the file (zeros included)"
    else:
        try:
            nb = float(mode)
        except ValueError:
            nb = float("nan")
        if not math.isfinite(nb) or nb < 0:
            raise SystemExit(f"--neighbour-weight {mode}: expected median, median-all or a number >= 0")
        rule = "fixed value"
    gw = {}
    for g in sorted(elig):
        w = weights.get(g, a.missing_weight) if g in listed_set else nb
        if w > 0:
            gw[g] = w
    if not gw:
        raise SystemExit("no gene with a simulable transcript has a positive weight")
    alloc = allocate_weighted(gw, a.total_reads)
    dropped = sum(1 for g in listed_set if weights.get(g, a.missing_weight) > 0 and g not in elig)
    for g, n in alloc.items():
        txs = sorted(elig[g], key=lambda r: r["transcript"])
        for r, k in zip(txs, split_equal(n, len(txs))):
            r["n_reads"], r["simulated"] = k, int(k > 0)
    info = dict(file=os.path.basename(a.weights), file_sha256=file_sha256(a.weights), total_reads=a.total_reads,
                lines=len(weights), soto_with_weight=len(soto_w), soto_missing=len(listed_set) - len(soto_w), weights_of_non_listed_genes_ignored=len([g for g in weights if g not in listed_set]),
                missing_weight=a.missing_weight, neighbour_weight=nb, neighbour_weight_rule=rule, genes_with_reads=sum(1 for n in alloc.values() if n > 0),
                listed_with_positive_weight_and_no_simulable_transcript=dropped)
    return info


def identical_sequences(rows):
    """{sequence hash: set of genes} for the sequences that two or more genes share (byte-identical spliced transcript sequences of eligible transcripts)."""
    by_seq = collections.defaultdict(set)
    for r in rows:
        if r["seq_sha1"]:
            by_seq[r["seq_sha1"]].add(r["gene"])
    return {h: gs for h, gs in by_seq.items() if len(gs) >= 2}


def target_rows(listed, gene_info, rows, s1c, split, zero_weight=None):
    """One row per listed gene. The truth columns (simulable, exon_bp, n_exons, single_exon, n_tx_eligible, identical_seq_genes) come from the ELIGIBLE transcripts, so they do not depend on the
    arm's reads; status / reason / n_reads / n_tx_simulated are the arm's. `zero_weight`: None in the ideal arm, else the set of listed genes whose effective weight is 0."""
    by_gene = collections.defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    other = collections.defaultdict(set)
    for gs in identical_sequences(rows).values():
        for g in gs:
            other[g] |= gs - {g}
    out = []
    for g in listed:
        info = gene_info.get(g)
        s = s1c.get(g, {})
        row = dict(gene=g, gene_name="", s1c_name=s.get("name", ""), biotype="", s1c_biotype=s.get("biotype", ""), s1c_families=";".join(s.get("families", [])),
                   n_families=len(s.get("families", [])), half=split.get(g, ""), source=source_of(g), sd98=s.get("sd98", ""), chrom="", strand="", start0="", end="", n_tx_annotated=0,
                   n_tx_eligible=0, n_tx_simulated=0, n_reads=0, exon_bp=0, n_exons=0, single_exon="", simulable=0, identical_seq_genes=0, status="skipped", reason="")
        if info is None:
            row["reason"] = "not in annotation"
            out.append(row)
            continue
        txs = by_gene[g]
        elig = [r for r in txs if r["eligible"]]
        sim_tx = [r for r in txs if r["simulated"]]
        union = merge_blocks([b for r in elig for b in blocks_of(r["canon_blocks"])])
        row.update(gene_name=info["name"], biotype=info["biotype"], chrom=info["chrom"], strand=info["strand"], start0=info["start0"], end=info["end"], n_tx_annotated=info["n_tx"],
                   n_tx_eligible=len(elig), n_tx_simulated=len(sim_tx), n_reads=sum(r["n_reads"] for r in txs), exon_bp=sum(e - s for s, e in union), n_exons=len(union),
                   single_exon=(int(all(not r["chain"] for r in elig)) if elig else ""), simulable=int(bool(elig)), identical_seq_genes=len(other[g]))
        if not txs:
            row["reason"] = "no transcript with exons"
        elif not elig:
            c = collections.Counter(reason_class(r["skip_reason"]) for r in txs)
            row["reason"] = "no simulable transcript (" + ", ".join(f"{k}: {c[k]}" for k in sorted(c)) + ")"
        elif not sim_tx:
            row["reason"] = "weight 0" if zero_weight is not None and g in zero_weight else "0 reads allocated"
        else:
            row["status"] = "simulated"
        out.append(row)
    return out


def not_simulable_manifest(listed, gene_info, rows, targets):
    """The listed genes that no rule lets us simulate: how many, which of them are in the SD98 gene set, and why."""
    by_gene = collections.defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    reasons, n, sd98 = collections.Counter(), 0, 0
    for t in targets:
        if t["simulable"]:
            continue
        n += 1
        sd98 += int(t["sd98"] == "Yes")
        if t["gene"] not in gene_info:
            reasons["not in annotation"] += 1
        elif not by_gene[t["gene"]]:
            reasons["no transcript with exons"] += 1
        else:
            reasons[" + ".join(sorted({reason_class(r["skip_reason"]) for r in by_gene[t["gene"]]}))] += 1
    return dict(genes=n, sd98_eligible=sd98, by_reason=dict(reasons))


# ------------------------------------------------------------------------------------------------- reads
_FA = None   # per-process genome handle of the pool workers


def _open_worker(genome):
    global _FA
    _FA = pysam.FastaFile(genome)


def _reads_text(fa, work, seed, err, jitter):
    """FASTQ text of the reads of a list of (transcript, gene, chrom, strand, canon_blocks, n) in the order given: sim_windows' read loop."""
    buf = io.StringIO()
    for t, g, chrom, st, canon, n in work:
        cb = blocks_of(canon)
        body0 = lib.spliced1(fa, chrom, [(s + 1, e) for s, e in cb], st)
        for k in range(n):
            rng = random.Random(sim.stable_seed(str(seed), t, str(k), "jitter"))
            lo, hi = rng.randint(0, jitter), rng.randint(0, jitter)
            body = body0[lo:len(body0) - hi] if len(body0) - hi > lo + 100 else body0
            (rd, q), = sim.simulate_reads(body, 1, err=err, indel=err / 3, seed=sim.stable_seed(str(seed), t, str(k), "errors"), trunc_frac=0.0)
            sim.write_fastq(buf, f"{t}|{g}|fl|{k}", (rd, q))
    return buf.getvalue()


def _chunk_job(job):
    work, seed, err, jitter = job
    return _reads_text(_FA, work, seed, err, jitter)


def chunks(rows):
    """Contiguous work units of about CHUNK_BP read bases, in row order (so concatenating their outputs in order is the single-process output whatever the process count)."""
    cur, bp = [], 0
    for r in rows:
        if not r["n_reads"]:
            continue
        cur.append((r["transcript"], r["gene"], r["chrom"], r["strand"], r["canon_blocks"], r["n_reads"]))
        bp += r["n_reads"] * r["spliced_len"]
        if bp >= CHUNK_BP:
            yield cur
            cur, bp = [], 0
    if cur:
        yield cur


def write_fastq(a, rows):
    h, n, bases = hashlib.sha256(), 0, 0
    jobs = ((w, a.seed, a.err, a.jitter) for w in chunks(rows))
    with open(a.out + ".fq", "wb") as fh:
        def put(text):
            nonlocal n, bases
            data = text.encode("ascii")
            fh.write(data)
            h.update(data)
            n += data.count(b"\n") // 4   # four lines per record, no embedded newlines
            bases += sum(len(ln) for ln in text.split("\n")[1::4])
        if a.procs == 1:
            fa = pysam.FastaFile(a.genome)
            for job in jobs:
                put(_reads_text(fa, *job))
            fa.close()
        else:
            ctx = multiprocessing.get_context("fork" if "fork" in multiprocessing.get_all_start_methods() else None)
            with ctx.Pool(a.procs, initializer=_open_worker, initargs=(a.genome,)) as pool:
                for text in pool.imap(_chunk_job, jobs):
                    put(text)
    return n, h.hexdigest(), bases


# ------------------------------------------------------------------------------------------------- writing
def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--genes", default=DEFAULT_GENES, help="S1C-format table with a `Gene ID` column, or one gene id per line (default: bench/soto/soto_famCN_S1C.tsv)")
    ap.add_argument("--gff", default=sw.CAT_GFF, help="tabix-indexed CAT/Liftoff GFF3 (default: the v2.0 slim GFF)")
    ap.add_argument("--genome", default=sw.GENOME, help="genome FASTA with .fai")
    ap.add_argument("--s1c", default=DEFAULT_GENES, help="table with `Gene Name`, `Biotype`, `Family ID` columns, for the targets table")
    ap.add_argument("--split", default=DEFAULT_SPLIT, help="soto_split table (gene_id, half), for the targets table")
    ap.add_argument("--pad", type=int, default=25_000)
    ap.add_argument("--reads-per-tx", type=int, default=10)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--err", type=float, default=0.001)
    ap.add_argument("--jitter", type=int, default=30)
    ap.add_argument("--procs", type=int, default=1, choices=[1, 2, 3, 4], help="worker processes for the read simulation (the output does not depend on it)")
    ap.add_argument("--first", type=int, default=None, help="use only the first N genes of the list (smoke runs, the determinism gate)")
    ap.add_argument("--no-reads", action="store_true", help="write the tables only (no FASTQ)")
    ap.add_argument("--weights", default=None, help="WEIGHTED ARM: TSV gene_id<TAB>weight")
    ap.add_argument("--total-reads", type=int, default=None, help="WEIGHTED ARM: reads to allocate over the genes")
    ap.add_argument("--neighbour-weight", default="median", help="weight of every non-listed neighbour gene: median (of the positive listed weights), median-all, or a number")
    ap.add_argument("--missing-weight", type=float, default=0.0, help="weight of a listed gene that is absent from --weights")
    a = ap.parse_args(argv)
    if (a.weights is None) != (a.total_reads is None):
        ap.error("--weights and --total-reads go together")
    if a.reads_per_tx < 1 or a.pad < 0 or a.jitter < 0 or (a.total_reads is not None and a.total_reads < 1):
        ap.error("--reads-per-tx and --total-reads must be >= 1, --pad and --jitter >= 0")
    return a


def simulate(a):
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    listed = read_gene_list(a.genes)
    if a.first is not None:
        listed = listed[:a.first]
    listed_set = set(listed)
    weights = read_weights(a.weights) if a.weights else None
    s1c, split = load_s1c(a.s1c), load_split(a.split)
    fa = pysam.FastaFile(a.genome)
    rows, gene_info, windows, n_sel, found = plan(a, listed, fa)
    fa.close()
    winfo = assign_reads(rows, a, listed_set, weights)
    zero = None if weights is None else {g for g in listed if weights.get(g, a.missing_weight) == 0}
    targets = target_rows(listed, gene_info, rows, s1c, split, zero)
    write_tsv(a.out + ".transcripts.tsv", TX_COLUMNS, rows)
    write_tsv(a.out + ".targets.tsv", TARGET_COLUMNS, targets)
    with open(a.out + ".windows.tsv", "w") as fh:
        fh.write("chrom\tstart0\tend\n")
        for c, lo, hi in windows:
            fh.write(f"{c}\t{lo}\t{hi}\n")
    skipped = [dict(transcript=r["transcript"], gene=r["gene"], gene_name=r["gene_name"], is_soto=r["is_soto"], reason=r["skip_reason"], spliced_len=r["spliced_len"]) for r in rows if r["skip_reason"]]
    write_tsv(a.out + ".skipped.tsv", ["transcript", "gene", "gene_name", "is_soto", "reason", "spliced_len"], skipped)
    shared = identical_sequences(rows)
    tx_of = collections.defaultdict(list)
    for r in rows:
        if r["seq_sha1"] in shared:
            tx_of[r["seq_sha1"]].append(r["transcript"])
    write_tsv(a.out + ".identical_seq.tsv", ["seq_sha1", "n_transcripts", "n_genes", "genes", "transcripts"],
              [dict(seq_sha1=h, n_transcripts=len(tx_of[h]), n_genes=len(shared[h]), genes=";".join(sorted(shared[h])), transcripts=";".join(sorted(tx_of[h]))) for h in sorted(shared)])
    n_reads_planned = sum(r["n_reads"] for r in rows)
    n_written, digest, read_bases = (None, None, None) if a.no_reads else write_fastq(a, rows)
    if n_written is not None and n_written != n_reads_planned:
        raise RuntimeError(f"wrote {n_written} reads, planned {n_reads_planned}")
    st = collections.Counter(t["status"] for t in targets)
    meta = dict(tool="sim_genes.py", prereg="docs/PREREG_ideal_arm_pilot_2026-10-07.md", genes_listed=len(listed), genes_in_annotation=len(found), genes_missing=len(listed) - len(found),
                pad=a.pad, windows=len(windows), window_bp=sum(hi - lo for _, lo, hi in windows), genes_in_windows=n_sel, neighbour_genes=n_sel - len(found),
                transcripts=len(rows), eligible_transcripts=sum(r["eligible"] for r in rows), simulated_transcripts=sum(r["simulated"] for r in rows),
                skipped=dict(collections.Counter(reason_class(r["skip_reason"]) for r in rows if r["skip_reason"])),
                longest_molecule=max((r["spliced_len"] for r in rows if r["simulated"]), default=0), reads=n_reads_planned, seed=a.seed, err=a.err, jitter=a.jitter,
                reads_per_tx=(a.reads_per_tx if weights is None else None), read_bases=read_bases, weights=winfo,
                not_simulable=not_simulable_manifest(listed, gene_info, rows, targets),
                identical_sequences=dict(sequences_in_2plus_genes=len(shared), genes=len({g for gs in shared.values() for g in gs}),
                                         listed_genes=len({g for gs in shared.values() for g in gs if g in listed_set}),
                                         sequences_in_2plus_listed_genes=sum(1 for gs in shared.values() if len(gs & listed_set) >= 2),
                                         listed_genes_in_those=len({g for gs in shared.values() if len(gs & listed_set) >= 2 for g in gs & listed_set})),
                targets=dict(simulated=st["simulated"], skipped=st["skipped"],
                reasons=dict(collections.Counter(re.sub(r"\(.*", "", t["reason"]).strip() for t in targets if t["reason"]))),
                fastq_sha256=digest, genome=os.path.basename(a.genome), gff=os.path.basename(a.gff))
    with open(a.out + ".sim.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta))
    return meta


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in MODES:
        return MODES[argv[0]](argv[1:])
    simulate(parse_args(argv))
    return 0


# ------------------------------------------------------------------------------------------------- gates (written apart from the generator: they re-read the files)
def tsv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def json_file(path):
    with open(path) as fh:
        return json.load(fh)


READ_NAME = re.compile(r"^([^|]+)\|([^|]+)\|fl\|(\d+)$")
GATE_MOTIFS = {"GTAG", "GCAG", "ATAC"}   # re-stated here, not imported from sim_windows: the check must not share the generator's table


def fastq_iter(path):
    """(header without '@', sequence, quality) per record; a malformed record raises ValueError."""
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq, plus, qual = fh.readline().rstrip("\n"), fh.readline(), fh.readline().rstrip("\n")
            if not h.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"{path}: malformed record at {h.strip()!r}")
            yield h.rstrip("\n")[1:], seq, qual


def gate_manifest(prefix, genes=None, genome=None, need_reads=True):
    """G1: the FASTQ matches the per-transcript counts and names, every listed gene is simulated or listed with a reason, and (with a genome) every kept intron is canonical.
    Returns (ok, lines); the first line starts with `G1 PASS` or `G1 FAIL`."""
    fails = []
    rows = tsv_rows(prefix + ".transcripts.tsv")
    targets = tsv_rows(prefix + ".targets.tsv")
    meta = json_file(prefix + ".sim.json")
    by_t = {r["transcript"]: r for r in rows}
    total = sum(int(r["n_reads"]) for r in rows)
    n = 0
    if need_reads:
        if not os.path.exists(prefix + ".fq"):
            fails.append("G1 FAIL reads: there is no FASTQ")
        else:
            count, ks, bad_name, wrong_gene = collections.Counter(), collections.defaultdict(set), [], []
            try:
                for name, _seq, _q in fastq_iter(prefix + ".fq"):
                    n += 1
                    m = READ_NAME.match(name)
                    if not m or m.group(1) not in by_t or not int(by_t[m.group(1)]["n_reads"]):
                        bad_name.append(name)
                        continue
                    t, g, k = m.group(1), m.group(2), int(m.group(3))
                    if g != by_t[t]["gene"]:
                        wrong_gene.append(name)
                    count[t] += 1
                    ks[t].add(k)
            except ValueError as e:
                fails.append(f"G1 FAIL reads: {e}")
            if n != total:
                fails.append(f"G1 FAIL reads: the FASTQ holds {n} records, the per-transcript counts sum to {total}")
            if n != meta.get("reads"):
                fails.append(f"G1 FAIL reads: the FASTQ holds {n} records, sim.json says {meta.get('reads')}")
            if bad_name:
                fails.append(f"G1 FAIL names: {len(bad_name)} read name(s) are not a listed transcript with reads (first: {bad_name[0]})")
            if wrong_gene:
                fails.append(f"G1 FAIL names: {len(wrong_gene)} read(s) carry a gene that is not their transcript's (first: {wrong_gene[0]})")
            off = [t for t, r in by_t.items() if int(r["n_reads"]) != count[t] or ks[t] != set(range(int(r["n_reads"])))]
            if off:
                fails.append(f"G1 FAIL names: {len(off)} transcript(s) whose reads are not exactly k = 0..n_reads-1 (first: {off[0]})")
    listed = read_gene_list(genes) if genes else [r["gene"] for r in targets]
    trow = collections.defaultdict(list)
    for r in targets:
        trow[r["gene"]].append(r)
    for g in listed:
        if len(trow[g]) != 1:
            fails.append(f"G1 FAIL targets: listed gene {g} has {len(trow[g])} rows in the targets table")
            continue
        r = trow[g][0]
        if r["status"] == "skipped" and not r["reason"].strip():
            fails.append(f"G1 FAIL targets: {g} is skipped without a reason")
        elif r["status"] == "simulated" and not (int(r["n_reads"]) > 0 and int(r["n_tx_simulated"]) > 0):
            fails.append(f"G1 FAIL targets: {g} is marked simulated but has no reads")
        elif r["status"] not in ("simulated", "skipped"):
            fails.append(f"G1 FAIL targets: {g} has status {r['status']!r}")
    extra = sorted(set(trow) - set(listed))
    if extra:
        fails.append(f"G1 FAIL targets: {len(extra)} gene(s) in the targets table are not in the gene list (first: {extra[0]})")
    n_int = 0
    if genome:
        fa = pysam.FastaFile(genome)
        comp = str.maketrans("ACGT", "TGCA")
        for r in rows:
            if not int(r["n_reads"]):
                continue
            if not 120 <= int(r["spliced_len"]) <= 30_000:
                fails.append(f"G1 FAIL introns: {r['transcript']} has a spliced length of {r['spliced_len']} bp")
            for x in (r["chain"].split(",") if r["chain"] else []):
                d, acc = map(int, x.split("-"))
                n_int += 1
                donor, accep = fa.fetch(r["chrom"], d, d + 2).upper(), fa.fetch(r["chrom"], acc - 2, acc).upper()
                motif = donor + accep if r["strand"] == "+" else accep[::-1].translate(comp) + donor[::-1].translate(comp)
                if motif not in GATE_MOTIFS or acc - d < 50 or acc - d > 200_000:
                    fails.append(f"G1 FAIL introns: {r['transcript']} intron {d}-{acc} on {r['strand']} is not a canonical intron of 50 bp to 200 kb (motif {motif})")
        fa.close()
    if fails:
        return False, fails
    return True, [f"G1 PASS manifest: {n} reads = the per-transcript counts, every read name is a listed transcript with its gene, {len(listed)} listed genes simulated or skipped with a reason"
                  + (f", {n_int} kept introns canonical" if genome else "")]


def _unique_offset(window, anchor, first, last):
    """Offsets p in [first, last] with window[p:p+len(anchor)] == anchor; the offset if it is unique, else None."""
    hits, p = [], window.find(anchor, first)
    while p != -1 and p <= last:
        hits.append(p)
        p = window.find(anchor, p + 1)
    return hits[0] if len(hits) == 1 else None


def gate_read_ends(prefix, genome, min_reads=8):
    """G2: the 5' and 3' trims of the reads (read off the transcript sequence with a 20-base anchor) are not constant. Pooled over the reads each end takes >= 2 values, none exceeds the
    declared jitter, and no transcript that is long enough for the jitter to apply (spliced length >= 2 jitter + 101) with >= `min_reads` readable reads has one trim on all of them
    (with 31 possible trims chance gives that probability 31^-(n-1) < 1e-10 for n = 8: any such transcript is a defect, not luck)."""
    meta = json_file(prefix + ".sim.json")
    jitter = int(meta["jitter"])
    rows = {r["transcript"]: r for r in tsv_rows(prefix + ".transcripts.tsv") if int(r["n_reads"])}
    fa = pysam.FastaFile(genome)
    cache = collections.OrderedDict()

    def body(t):
        if t not in cache:
            r = rows[t]
            cache[t] = lib.spliced1(fa, r["chrom"], [(s + 1, e) for s, e in blocks_of(r["canon_blocks"])], r["strand"])
            if len(cache) > 32:
                cache.popitem(last=False)
        return cache[t]

    o5, o3 = collections.defaultdict(list), collections.defaultdict(list)
    n_read = 0
    for name, seq, _q in fastq_iter(prefix + ".fq"):
        m = READ_NAME.match(name)
        if not m or m.group(1) not in rows:
            continue
        t, b = m.group(1), body(m.group(1))
        n_read += 1
        lo = _unique_offset(b[:jitter + ANCHOR], seq[:ANCHOR], 0, jitter)
        start = max(0, len(b) - jitter - ANCHOR)
        q = _unique_offset(b[start:], seq[-ANCHOR:], 0, len(b) - start - ANCHOR)
        if lo is not None:
            o5[t].append(lo)
        if q is not None:
            o3[t].append(len(b) - (start + q + ANCHOR))
    fa.close()
    fails, report = [], {}
    for end, o in (("5'", o5), ("3'", o3)):
        vals = collections.Counter(v for lst in o.values() for v in lst)
        if len(vals) < 2:
            fails.append(f"G2 FAIL read ends: the {end} trim takes {len(vals)} distinct value(s) over {sum(vals.values())} readable reads")
        if vals and max(vals) > jitter:
            fails.append(f"G2 FAIL read ends: a {end} trim of {max(vals)} bp exceeds the declared jitter of {jitter} bp")
        tested = 0
        for t, lst in sorted(o.items()):
            if len(lst) >= min_reads and int(rows[t]["spliced_len"]) >= 2 * jitter + 101:
                tested += 1
                if len(set(lst)) == 1:
                    fails.append(f"G2 FAIL read ends: {t} has a constant {end} trim ({lst[0]} bp) on all {len(lst)} readable reads")
        report[end] = (len(vals), max(vals, default=0), tested)
    if fails:
        return False, fails
    (k5, m5, t5), (k3, m3, t3) = report["5'"], report["3'"]
    return True, [f"G2 PASS read ends vary: the 5' trim takes {k5} values (max {m5} bp), the 3' trim {k3} (max {m3} bp); {t5} / {t3} transcripts tested for a constant 5' / 3' trim, none constant"]


def gate_compare(prefix_a, prefix_b, suffixes=SUFFIXES):
    """G3 (simulation half): the two simulations agree byte for byte on every output file that either has."""
    differ, same = [], []
    for suf in suffixes:
        ea, eb = os.path.exists(prefix_a + suf), os.path.exists(prefix_b + suf)
        if not ea and not eb:
            continue
        if ea != eb or file_sha256(prefix_a + suf) != file_sha256(prefix_b + suf):
            differ.append(suf)
        else:
            same.append(suf)
    if differ or not same:
        return False, ["G3 FAIL: " + ("differs: " + ", ".join(differ) if differ else "no output files to compare")]
    return True, [f"G3 PASS: {len(same)} output files byte-identical ({', '.join(same)})"]


def gate_bamqc(prefix, bam):
    """G4: mapping sanity from the BAM alone. PASS = the BAM holds exactly one primary-or-unmapped record per FASTQ read, under the FASTQ's own names. Printed beside (no bar):
    the share of reads with a primary alignment and the share of those whose primary alignment overlaps an exon of the SIMULATED source gene (paralog ambiguity), all / listed / neighbour."""
    rows = tsv_rows(prefix + ".transcripts.tsv")
    tx = {r["transcript"]: r for r in rows if int(r["n_reads"])}
    gene_blocks = collections.defaultdict(list)
    for r in tx.values():
        gene_blocks[r["gene"]].extend(blocks_of(r["canon_blocks"]))
    gene_blocks = {g: merge_blocks(b) for g, b in gene_blocks.items()}
    gene_chrom = {r["gene"]: r["chrom"] for r in tx.values()}
    soto = {r["transcript"]: r["is_soto"] == "1" for r in tx.values()}
    fq_names = collections.Counter(name for name, _s, _q in fastq_iter(prefix + ".fq"))
    st = dict(reads=sum(fq_names.values()), primary_or_unmapped_records=0, mapped=0, unmapped=0, secondary_records=0, primary_on_source_gene=0,
              soto=dict(reads=0, mapped=0, primary_on_source_gene=0), neighbour=dict(reads=0, mapped=0, primary_on_source_gene=0))
    seen, fails, unknown = collections.Counter(), [], 0
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        for rec in fh.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                st["secondary_records"] += 1
                continue
            q = rec.query_name
            seen[q] += 1
            st["primary_or_unmapped_records"] += 1
            m = READ_NAME.match(q)
            if not m or m.group(1) not in tx:
                unknown += 1
                continue
            grp = st["soto" if soto[m.group(1)] else "neighbour"]
            grp["reads"] += 1
            if rec.is_unmapped:
                st["unmapped"] += 1
                continue
            st["mapped"] += 1
            grp["mapped"] += 1
            g = m.group(2)
            if rec.reference_name == gene_chrom.get(g) and any(s < rec.reference_end and e > rec.reference_start for s, e in gene_blocks.get(g, [])):
                st["primary_on_source_gene"] += 1
                grp["primary_on_source_gene"] += 1
    dup = sorted(q for q, c in seen.items() if c > 1)
    missing = sorted(set(fq_names) - set(seen))
    extra = sorted(set(seen) - set(fq_names))
    if st["primary_or_unmapped_records"] != st["reads"]:
        fails.append(f"G4 FAIL bam: {st['primary_or_unmapped_records']} primary-or-unmapped records for {st['reads']} FASTQ reads")
    if dup:
        fails.append(f"G4 FAIL bam: {len(dup)} read(s) with more than one primary-or-unmapped record (first: {dup[0]})")
    if missing:
        fails.append(f"G4 FAIL bam: {len(missing)} FASTQ read(s) with no record (first: {missing[0]})")
    if extra or unknown:
        fails.append(f"G4 FAIL bam: {len(extra)} record name(s) not in the FASTQ, {unknown} not a listed transcript")
    share = lambda a, b: (f"{100.0 * a / b:.2f}%" if b else "n/a")
    info = [f"    reads with a primary alignment: {st['mapped']}/{st['reads']} = {share(st['mapped'], st['reads'])}; primary alignment on an exon of its source gene: "
            f"{st['primary_on_source_gene']}/{st['mapped']} = {share(st['primary_on_source_gene'], st['mapped'])} (no bar; secondary records: {st['secondary_records']})",
            f"    listed-gene reads: on source {st['soto']['primary_on_source_gene']}/{st['soto']['mapped']} = {share(st['soto']['primary_on_source_gene'], st['soto']['mapped'])}; "
            f"neighbour reads: {st['neighbour']['primary_on_source_gene']}/{st['neighbour']['mapped']} = {share(st['neighbour']['primary_on_source_gene'], st['neighbour']['mapped'])}"]
    if fails:
        return False, fails + info, st
    return True, [f"G4 PASS mapping sanity: one primary-or-unmapped record per read ({st['reads']})"] + info, st


def gate_depth_subset(small, large):
    """G-depth: the shallow arm is exactly the first reads of the deep arm. Read k of transcript t depends only on (seed, t, k), so a pool simulated at depth N contains the depth-n arm
    (n < N) read for read; this checks it on the files: same transcript tables apart from the read counts, same simulation parameters, and every read of the shallow FASTQ byte-identical to the
    read of the same name in the deep one, with nothing else of the deep one below the shallow count. Mapping the deep pool once and cutting the shallow arm by read name is then exact."""
    fails = []
    ms, ml = json_file(small + ".sim.json"), json_file(large + ".sim.json")
    for k in ("seed", "err", "jitter", "pad", "genes_listed", "transcripts", "eligible_transcripts"):
        if ms.get(k) != ml.get(k):
            fails.append(f"G-depth FAIL parameters: {k} is {ms.get(k)} in the shallow arm and {ml.get(k)} in the deep one")
    rs, rl = tsv_rows(small + ".transcripts.tsv"), tsv_rows(large + ".transcripts.tsv")
    cols = [c for c in TX_COLUMNS if c not in ("n_reads", "simulated")]
    if len(rs) != len(rl) or any(a[c] != b[c] for a, b in zip(rs, rl) for c in cols):
        fails.append("G-depth FAIL tables: the transcript tables differ in more than the read counts")
    n_small = {r["transcript"]: int(r["n_reads"]) for r in rs}
    n_large = {r["transcript"]: int(r["n_reads"]) for r in rl}
    if any(n > n_large.get(t, 0) for t, n in n_small.items()):
        fails.append("G-depth FAIL tables: a transcript has more reads in the shallow arm than in the deep one")
    digest = {}
    for name, seq, q in fastq_iter(small + ".fq"):
        digest[name] = hashlib.sha1((seq + "\n" + q).encode("ascii")).digest()
    if len(digest) != sum(n_small.values()):
        fails.append(f"G-depth FAIL reads: the shallow FASTQ holds {len(digest)} distinct names for {sum(n_small.values())} planned reads")
    seen = differ = 0
    for name, seq, q in fastq_iter(large + ".fq"):
        m = READ_NAME.match(name)
        if m and int(m.group(3)) < n_small.get(m.group(1), 0):
            d = digest.get(name)
            if d is None:
                fails.append(f"G-depth FAIL reads: {name} is among the first reads of the deep arm but not in the shallow one")
            elif d != hashlib.sha1((seq + "\n" + q).encode("ascii")).digest():
                differ += 1
            else:
                seen += 1
    if differ:
        fails.append(f"G-depth FAIL reads: {differ} read(s) of the shallow arm differ from the read of the same name in the deep one")
    if seen != len(digest) and not differ:
        fails.append(f"G-depth FAIL reads: {len(digest) - seen} read(s) of the shallow arm are not in the deep one")
    if fails:
        return False, fails
    return True, [f"G-depth PASS: all {seen} reads of the shallow arm (depth {ms.get('reads_per_tx')}) are byte-identical to the first reads of each transcript in the deep arm (depth {ml.get('reads_per_tx')})"]


def mode_depth_subset(argv):
    ap = argparse.ArgumentParser(description="G-depth")
    ap.add_argument("--json", default=None, help="also write {ok, lines} here")
    ap.add_argument("small")
    ap.add_argument("large")
    a = ap.parse_args(argv)
    ok, lines = gate_depth_subset(a.small, a.large)
    print("\n".join(lines))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(dict(ok=ok, lines=lines), fh, indent=1)
    return 0 if ok else 1


def mode_readnames(argv):
    ap = _arg_parser("write the name of every read of a simulation, in file order (the list `samtools view -N` cuts a shallower arm out of a deeper pool with)")
    ap.add_argument("--to", required=True)
    a = ap.parse_args(argv)
    with open(a.to, "w") as fh:
        for name, _s, _q in fastq_iter(a.out + ".fq"):
            fh.write(name + "\n")
    return 0


def mode_control(argv):
    """The annotation written as an assembled product, for the scorer's controls: CTRL.gtf (every ELIGIBLE transcript, gene_id = the gene) and CTRL.loci.gff3 (one node per gene, the families-stage
    format: the gene's longest eligible isoform, ties to the smaller transcript id). --shift moves every coordinate (the negative control: a perfect annotation in the wrong place)."""
    ap = _arg_parser("annotation-as-loci control")
    ap.add_argument("--to", required=True, help="control prefix")
    ap.add_argument("--shift", type=int, default=0, help="add this many bp to every coordinate")
    a = ap.parse_args(argv)
    rows = [r for r in tsv_rows(a.out + ".transcripts.tsv") if int(r["eligible"])]
    order = {}
    for r in rows:
        order.setdefault(r["chrom"], len(order))
    rows.sort(key=lambda r: (order[r["chrom"]], blocks_of(r["canon_blocks"])[0][0], r["transcript"]))
    lines = []
    for r in rows:
        bl = [(x + a.shift, y + a.shift) for x, y in blocks_of(r["canon_blocks"])]
        lines.append(f'{r["chrom"]}\trustle\ttranscript\t{bl[0][0] + 1}\t{bl[-1][1]}\t.\t{r["strand"]}\t.\tgene_id "{r["gene"]}"; transcript_id "{r["transcript"]}"; reads "{r["n_reads"]}"; cov "10.0"; TPM "1.0";')
        for i, (x, y) in enumerate(bl, 1):
            lines.append(f'{r["chrom"]}\trustle\texon\t{x + 1}\t{y}\t.\t{r["strand"]}\t.\tgene_id "{r["gene"]}"; transcript_id "{r["transcript"]}"; exon_number "{i}";')
    with open(a.to + ".gtf", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    rep = {}
    for r in rows:
        if r["gene"] not in rep or (-int(r["spliced_len"]), r["transcript"]) < (-int(rep[r["gene"]]["spliced_len"]), rep[r["gene"]]["transcript"]):
            rep[r["gene"]] = r
    out = ["##gff-version 3"]
    for g in sorted(rep, key=lambda g: (order[rep[g]["chrom"]], blocks_of(rep[g]["canon_blocks"])[0][0], g)):
        r = rep[g]
        bl = [(x + a.shift, y + a.shift) for x, y in blocks_of(r["canon_blocks"])]
        out.append(f'{r["chrom"]}\t.\tgene\t{bl[0][0] + 1}\t{bl[-1][1]}\t.\t{r["strand"]}\t.\tID=gene-{g};Name={g}')
        out += [f'{r["chrom"]}\t.\texon\t{x + 1}\t{y}\t.\t{r["strand"]}\t.\tParent=gene-{g};gene={g}' for x, y in bl]
    with open(a.to + ".loci.gff3", "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"control: {len(rows)} transcripts in {a.to}.gtf, {len(rep)} nodes in {a.to}.loci.gff3 (shift {a.shift})")
    return 0


def _arg_parser(desc, with_genes=False):
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--out", required=True, help="the simulation prefix")
    return ap


def mode_gates(argv):
    ap = _arg_parser("G1 and G2")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--genes", default=None, help="the gene list the simulation was run with (default: the genes of the targets table)")
    ap.add_argument("--json", default=None, help="where to write the verdicts (default: PREFIX.gates.json)")
    a = ap.parse_args(argv)
    ok1, l1 = gate_manifest(a.out, genes=a.genes, genome=a.genome)
    ok2, l2 = gate_read_ends(a.out, a.genome) if os.path.exists(a.out + ".fq") else (False, ["G2 FAIL read ends: there is no FASTQ"])
    for ln in l1 + l2:
        print(ln)
    with open(a.json or a.out + ".gates.json", "w") as fh:
        json.dump(dict(G1=dict(ok=ok1, lines=l1), G2=dict(ok=ok2, lines=l2)), fh, indent=1)
    return 0 if ok1 and ok2 else 1


def mode_compare(argv):
    ap = argparse.ArgumentParser(description="G3, simulation half")
    ap.add_argument("--only", default=None, help="comma-separated output suffixes to compare (default: every output file either simulation has)")
    ap.add_argument("--json", default=None, help="also write {ok, lines} here")
    ap.add_argument("a")
    ap.add_argument("b")
    a = ap.parse_args(argv)
    ok, lines = gate_compare(a.a, a.b, suffixes=tuple(a.only.split(",")) if a.only else SUFFIXES)
    print("\n".join(lines))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(dict(ok=ok, lines=lines), fh, indent=1)
    return 0 if ok else 1


def mode_bamqc(argv):
    ap = _arg_parser("G4")
    ap.add_argument("--bam", required=True)
    a = ap.parse_args(argv)
    ok, lines, st = gate_bamqc(a.out, a.bam)
    print("\n".join(lines))
    with open(a.out + ".g4.json", "w") as fh:
        json.dump(dict(ok=ok, lines=lines, stats=st), fh, indent=1)
    return 0 if ok else 1


MODES = {"gates": mode_gates, "compare": mode_compare, "bamqc": mode_bamqc, "depth-subset": mode_depth_subset, "readnames": mode_readnames, "control": mode_control}


if __name__ == "__main__":
    sys.exit(main())
