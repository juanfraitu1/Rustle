#!/usr/bin/env python3
"""readthrough_eval — the pre-registered scorer of the read-end readthrough junction filter.

Binding: docs/PREREG_readthrough_ends_representatives_2026-09-25.md (the "prereg"): §2 arms (BASE, R, RQ1, NULL),
§3 metrics (a)-(e), §4 substrates and power floors, §5 bar and verdicts. Nothing here chooses a threshold: every
tolerance, floor and bar below is copied from the prereg (constants block). Rule source: docs/READTHROUGH_G50K_AND_
LAST_EXON_2026-09-25.md and bench/mechanism/readthrough_rules.py (the instrument whose read-strand and read-end
conventions are reused for the read-derived clusters and the NULL junction census).

An ARM is an assembled GTF (one `gene_id` = one de novo locus) plus, optionally, its families (the driver's
`families` stage products `<prefix>.clusters.tsv` / `<prefix>.copies.tsv`; auto-detected as `<gtf stem>.fam.*` next
to the GTF, the `tools/rustle_pipeline.sh` layout). A SUBSTRATE is a sample of figures/samples.tsv (its BAM, genome,
RefSeq GFF, annotation GTF) restricted to a contig set. Every universe (read-supported genes, read-derived clusters,
reference loci) depends on the BAM and the annotation only, so it is IDENTICAL across arms; it is cached per
(sample, contig) under --work.

SUBCOMMANDS

  score    every metric of every arm on one substrate + the bar's per-clause verdict, as ONE tidy table.
      python3 bench/mechanism/readthrough_eval.py score --sample human_A119b --contigs chr20 \\
          --arm BASE=<base.gtf> --arm R=<r.gtf> --arm RQ1=<rq1.gtf> --arm NULL=<null.gtf> \\
          [--families ARM=<prefix>] [--base BASE] [--null NULL] [--out DIR] [--work DIR] \\
          [--fs-bin DIR] [--heldout] [--budget-s S]
      Substrate: --contigs c1,c2 (exact set) or --drop-contigs c1,c2 (every contig the sample's annotation covers,
      minus these; e.g. V1 = --sample human_A119b --drop-contigs chr16,chr20). Neither = every annotated contig.
      Output: DIR/readthrough_eval.<sample>.<substrate key>.tsv (+ .json provenance). A substrate that is not a
      development substrate (human_A119b chr16/chr20, gorilla_OR6737 NC_073244.2) is REFUSED unless --heldout is
      given: held-out substrates are scored only after the binary's sha1 is recorded (prereg §8).
      NULL gets metrics (a)-(d) only (prereg §2). Arms named on the command line in any order; BASE is the
      comparator of every clause, NULL the comparator of A1's null part (one NULL, matched to arm R, serves every arm).
      --budget-s: the BAM pass is cached per contig; a call that runs out of budget exits 75 before the next contig
      (run the same command again).

  null     the NULL arm's junction list (prereg §2 as realised by its Amendment 2): per contig, junctions drawn at
           random (seed 20260925, one RNG per sample:contig) from the canonical S >= 2 junctions that arm R does not
           flag (RQ1 ⊆ R), matched to arm R's flagged junctions by floor(log2 S), until the DISTINCT primary
           alignments carrying them reach those carrying an R-flagged junction.
      python3 bench/mechanism/readthrough_eval.py null --sample human_A119b --contigs chr20 \\
          --flags <arm R's <out>.readthrough_junctions.tsv> [--target <contig TAB alignments>] --out <null.tsv>
      --flags: every row is taken as flagged by R. Formats read: the assembler's dump (prereg Amendment 1 item 6:
      contig, donor, acceptor = 1-based inclusive intron, strand, S, U, V1, rule r|rq1), the design's rt_all table
      (start, end, strand, ..., RQ1; needs --contigs with ONE contig), or a NULL list. --target overrides the
      per-contig target (default: distinct primaries, -F 2308, carrying an R-flagged junction, from the BAM). The
      assembler's NULL arm removes the alignments whose intron chain holds a listed junction, as it removes flagged
      ones; <out>.summary.tsv gives target vs achieved per contig.

  synth    a SYNTHETIC arm for testing the scorer (never an arm of the prereg: transcripts are dropped AFTER assembly,
           the arms remove reads BEFORE pass 1): the GTF restricted to --contigs, every transcript carrying a listed
           junction dropped, then every gene_id that lost a transcript re-split into components of transcripts sharing
           a junction (single-exon transcripts join the first component whose exons they overlap).
      python3 bench/mechanism/readthrough_eval.py synth --gtf <base.gtf> --contigs chr20 --junctions <tsv> \\
          [--rq1-only] --out <synthetic.gtf>

  verdict  §5's verdict per arm from the per-substrate tables of `score` (development tables are ignored). Readings
           fixed here: a verdict substrate is a held-out one whose BASE FUSED >= 50; any "not_measured" clause on a
           verdict substrate caps the arm at "keep opt-in"; a clause below its power floor is reported, not judged;
           "G1 or G3 is lower on more than half" is read per clause (G1 on > half, or G3 on > half); the NPIP cap
           (§6) is an input (default unknown = capped), this scorer does not measure NPIP.
      python3 bench/mechanism/readthrough_eval.py verdict <table.tsv> ... [--npip-cap none|capped|unknown]

  selftest unit fixtures (tiny GFF / GTF / genome / BAM written to a temp dir): every metric's rule, the verdict
           logic and the NULL draw's determinism.

METRICS (prereg §3; the table's `metric` column)
  a.*  FUSED (primary) = loci holding >= 1 spliced transcript whose exon union overlaps (>= 1 bp) the exon unions of
       >= 2 annotated genes on the representative's strand whose spans do not overlap each other; rep_fused (the
       representative's exons), span_cover (the span holds >= 50% of the exon union of >= 2 such genes),
       fused_junction (a transcript's junction joins the exons of two such genes), strata, absorbed genes.
  b.*  TES / TSS recovered genes on the fixed read-supported gene universe (±25 / ±250 bp; span ends, and the
       representative's own ends); read-derived start / end clusters (fixed per sample) hit by a locus span end.
  c.*  gffcompare 0.12.10 (-r annotation restricted to the substrate contigs): intron-chain precision = query
       multi-exon transcripts with class '=' / query multi-exon transcripts (.tmap); matching reference intron chains
       (.stats); lost / gained matched reference chains vs BASE; multi-gene transcripts.
  d.*  loci in the Liftoff framework (cov(G | R) >= 0.5, strand ignored): the merged Liftoff table when it exists,
       else the annotation's own gene / pseudogene records ("annotation-only"; extra copies then not measured).
  e.*  families: human = Ensembl Compara families at Primates (family_score --chrom ALL --per-family --pairwise,
       the dev build); Soto 2025 descriptive; apes = Liftoff copy-pair recall. "not measured" when an input is absent.

Conventions: GTF/GFF coordinates 1-based inclusive; internal intervals 0-based half-open; read strand = (`ts`
absent -> '+') XOR FLAG 0x10 (readthrough_rules.py); primary = FLAG & 0x904 == 0 (-F 2308).
"""
from __future__ import annotations

import argparse
import array
import bisect
import collections
import hashlib
import json
import math
import os
import pickle
import random
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
FIGURES = REPO / "figures"
if str(FIGURES) not in sys.path:
    sys.path.insert(0, str(FIGURES))

import figlib  # noqa: E402
import samples  # noqa: E402
import assembly  # noqa: E402
import _liftoff as L  # noqa: E402

# ---------------------------------------------------------------- frozen by the prereg (never tuned here)
PREREG = "docs/PREREG_readthrough_ends_representatives_2026-09-25.md"
SEED = 20260925
TES_TOL, TSS_TOL = 25, 250               # §3b
TES_GAP, TSS_GAP, CLUSTER_MIN = 25, 100, 3   # §3b read-derived clusters (the doc's §5 / §7 rules)
PRIMING_WIN, PRIMING_MAX_A = 20, 12     # §3b internal priming: < 12 A in the 20 genomic bases downstream
SUPPORT_READS = 2                        # §3: >= 2 primary reads with an aligned block on the exon union
SPAN_COVER_FRAC = 0.5                    # §3a SPAN-COVER
MATCH_COV = L.MATCH_COV                  # §3d cov(G | R) >= 0.5
SHORT_BP = L.SHORT_BP                    # §3d exon union >= 200 bp
X_A1 = 0.10                              # §5 A1
Y = 0.01                                 # §5 G2 / G4
G5_TOL, G5_REFUTE = 0.005, 0.010         # §5 G5 / refute
FLOORS = {"a": 50, "b": 1000, "c": 500, "d": 1000, "e_human": 30, "e_ape": 30}   # §4 power floors
DEV = {"human_A119b": {"chr16", "chr20"}, "gorilla_OR6737": {"NC_073244.2"}}      # §4 development substrates
HUMAN = "human"
DEFAULT_WORK = Path("/mnt/linuxdisk/tmp/rustle_figures_dev/readthrough_ends/eval")
DEFAULT_FS_BIN = Path("/mnt/linuxdisk/home/juanfraitu/rustle_target_dev/release")
GFFCOMPARE_VERSION = "0.12.10"
CACHE_VERSION = "1"
PENDING = 75
CANON = {"+": {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}, "-": {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}}
COMP = str.maketrans("ACGTN", "TGCAN")
BIN = 10_000


def log(msg: str):
    print(f"[readthrough_eval] {msg}", file=sys.stderr, flush=True)


def fp(path) -> str:
    return L.fp(path)


def sha1_file(path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def key16(*parts) -> str:
    return hashlib.sha1("\n".join(map(str, parts)).encode()).hexdigest()[:16]


class Budget:
    def __init__(self, seconds: float):
        self.limit = seconds or 0
        self.t0 = time.time()

    def check(self, what: str):
        if self.limit > 0 and time.time() - self.t0 > self.limit:
            log(f"budget of {self.limit:.0f} s used up; next unit: {what}; run the same command again")
            sys.exit(PENDING)


# ================================================================ substrate
class Substrate:
    def __init__(self, cfg: dict, sample: str, contigs: str | None, drop: str | None):
        row = samples.get(cfg, sample)
        self.sid, self.species = row["id"], row["species"]
        self.bam, self.fasta, self.gff = row["bam"], row["fasta"], row["annotation_gff"]
        genome = assembly.genome_contigs(cfg, self.sid)
        if contigs:
            want = [c.strip() for c in contigs.split(",") if c.strip()]
            bad = [c for c in want if c not in genome]
            if bad:
                raise SystemExit(f"contigs not in the {self.sid} genome: {bad}")
            self.contigs = [c for c in genome if c in set(want)]
            self.key = "-".join(self.contigs) if len(self.contigs) <= 4 else f"{len(self.contigs)}contigs_{key16(*self.contigs)[:8]}"
        else:
            dropped = {c.strip() for c in (drop or "").split(",") if c.strip()}
            ann = assembly.annotated_contigs(cfg, self.sid)
            self.contigs = [c for c in ann if c not in dropped]
            self.key = "annotated" + ("_minus_" + "_".join(sorted(dropped)) if dropped else "")
        if not self.contigs:
            raise SystemExit("empty substrate")
        self.dev = self.sid in DEV and set(self.contigs) <= DEV[self.sid]
        self.cset = set(self.contigs)
        self.cfg = cfg

    def label(self) -> str:
        return f"{self.sid}:{self.key}"


# ================================================================ annotation (one pass over the RefSeq GFF)
_RE = {k: re.compile(rf"(?:^|;){k}=([^;\n]*)") for k in ("ID", "Parent", "Name", "gene_biotype", "description")}


class Rec:
    """A gene / pseudogene record: exon union by the Liftoff framework's rule (merged exon children; CDS if none;
    else the span), 0-based half-open; `ends` = (5' end, 3' end) of every annotated transcript (1-based)."""
    __slots__ = ("idx", "id", "name", "type", "biotype", "contig", "strand", "start1", "end1", "rt", "iv", "src",
                 "ends", "support")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)
        self.support = 0

    @property
    def pc(self) -> bool:
        return self.biotype == "protein_coding" and not re.match(r"^LOC\d+$", self.name)

    @property
    def lnc_or_loc(self) -> bool:
        return self.biotype == "lncRNA" or bool(re.match(r"^LOC\d+$", self.name))


def parse_gff(gff: str, contigs: set) -> dict:
    """{contig: [Rec]} for the gene / pseudogene records on `contigs` (record order = file order)."""
    ftype, parent, recs = {}, {}, {}
    exons, cds = collections.defaultdict(list), collections.defaultdict(list)
    with open(gff) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            t = line.find("\t")
            if line[:t] not in contigs:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = f[8]
            m = _RE["ID"].search(a)
            mid = m.group(1) if m else None
            m = _RE["Parent"].search(a)
            mp = m.group(1).split(",")[0] if m else None
            if mid:
                ftype[mid] = f[2]
                if mp:
                    parent[mid] = mp
            if f[2] in ("gene", "pseudogene") and mid:
                d = _RE["description"].search(a)
                desc = urllib.parse.unquote(d.group(1)).lower() if d else ""
                n = _RE["Name"].search(a)
                b = _RE["gene_biotype"].search(a)
                recs[mid] = Rec(idx=None, id=mid, name=n.group(1) if n else mid, type=f[2],
                                biotype=b.group(1) if b else "", contig=f[0], strand=f[6], start1=int(f[3]),
                                end1=int(f[4]), rt="readthrough" in desc, iv=None, src=None, ends=None)
            elif f[2] == "exon" and mp:
                exons[mp].append((int(f[3]) - 1, int(f[4])))
            elif f[2] == "CDS" and mp:
                cds[mp].append((int(f[3]) - 1, int(f[4])))
    memo: dict = {}

    def top(x):
        if x in memo:
            return memo[x]
        chain, cur, res = [], x, None
        for _ in range(64):
            if cur in recs:
                res = cur
                break
            chain.append(cur)
            cur = parent.get(cur)
            if cur is None:
                break
        for c in chain:
            memo[c] = res
        return res
    rex, rcds, rends = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(set)
    for p, ex in exons.items():
        r = top(p)
        if r is None:
            continue
        rex[r].extend(ex)
        s1, e1 = min(a for a, _ in ex) + 1, max(b for _, b in ex)
        rends[r].add((s1, e1) if recs[r].strand != "-" else (e1, s1))
    for p, ex in cds.items():
        r = top(p)
        if r is not None:
            rcds[r].extend(ex)
    out = collections.defaultdict(list)
    for rid, r in recs.items():
        if rex.get(rid):
            r.iv, r.src = L.merge_iv(rex[rid]), "exon"
        elif rcds.get(rid):
            r.iv, r.src = L.merge_iv(rcds[rid]), "cds"
        else:
            r.iv, r.src = [(r.start1 - 1, r.end1)], "span"
        r.ends = sorted(rends.get(rid, ()))
        r.idx = len(out[r.contig])
        out[r.contig].append(r)
    return dict(out)


def load_annotation(sub: Substrate, work: Path) -> dict:
    d = work / "annotation"
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{key16(CACHE_VERSION, fp(sub.gff), *sub.contigs)}.pkl"
    fields = [f for f in Rec.__slots__ if f != "support"]
    if p.exists():
        with open(p, "rb") as fh:
            raw = pickle.load(fh)      # plain tuples: loadable whatever module name the script runs under
        return {c: [Rec(**dict(zip(fields, t))) for t in v] for c, v in raw.items()}
    log(f"annotation: one pass over {sub.gff} for {len(sub.contigs)} contig(s)")
    t0 = time.time()
    ann = parse_gff(sub.gff, sub.cset)
    for c in sub.contigs:
        ann.setdefault(c, [])
    tmp = p.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump({c: [tuple(getattr(r, f) for f in fields) for r in v] for c, v in ann.items()}, fh,
                    protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(p)
    log(f"annotation: {sum(len(v) for v in ann.values()):,} records ({time.time() - t0:.0f} s)")
    return ann


# ================================================================ BAM pass: record support + read-derived clusters
def read_strand_ends(rd):
    """(strand, 5' end, 3' end, introns 1-based inclusive) of a spliced alignment, or None (readthrough_rules.py)."""
    pos, introns = rd.reference_start, []
    for op, n in rd.cigartuples:
        if op == 3:
            introns.append((pos + 1, pos + n))
            pos += n
        elif op in (0, 2, 7, 8):
            pos += n
    if not introns:
        return None
    ts = rd.get_tag("ts") if rd.has_tag("ts") else "+"
    st = "+" if (ts == "+") != rd.is_reverse else "-"
    e5, e3 = (rd.reference_start + 1, rd.reference_end) if st == "+" else (rd.reference_end, rd.reference_start + 1)
    return st, e5, e3, introns


def clusters(pos: list, gap: int) -> list:
    """[(mode, n, lo, hi)] of sorted positions; a new cluster when the gap to the previous position exceeds `gap`;
    mode = most frequent position (ties: the smallest)."""
    out, grp = [], []
    for x in pos:
        if grp and x - grp[-1] > gap:
            out.append(grp)
            grp = []
        grp.append(x)
    if grp:
        out.append(grp)
    res = []
    for g in out:
        c = collections.Counter(g)
        mode = max(c.items(), key=lambda kv: (kv[1], -kv[0]))[0]
        res.append((mode, len(g), g[0], g[-1]))
    return res


def internal_priming(fa, contig: str, mode: int, strand: str) -> bool:
    """True when >= 12 A in the 20 genomic bases downstream of the 3' end `mode` (1-based), strand-aware
    (readthrough_rules.py's fetch windows)."""
    if strand == "+":
        down = fa.fetch(contig, mode, mode + PRIMING_WIN).upper()
    else:
        down = fa.fetch(contig, max(0, mode - PRIMING_WIN - 1), mode - 1).upper().translate(COMP)[::-1]
    return down.count("A") >= PRIMING_MAX_A


def read_pass(sub: Substrate, ann: dict, work: Path, budget: Budget) -> dict:
    """Per contig (cached): capped read support of every record, read-derived TES / TSS clusters."""
    import pysam
    d = work / "reads" / sub.sid / key16(CACHE_VERSION, fp(sub.bam), fp(sub.gff), fp(sub.fasta))
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    bam = fa = None
    for c in sub.contigs:
        p = d / f"{c}.json"
        if p.exists():
            out[c] = json.loads(p.read_text())
            continue
        budget.check(f"BAM pass {sub.sid} {c}")
        if bam is None:
            bam, fa = pysam.AlignmentFile(sub.bam), pysam.FastaFile(sub.fasta)
        t0 = time.time()
        recs = ann.get(c, [])
        bins = collections.defaultdict(list)
        for r in recs:
            for s, e in r.iv:
                for b in range(s // BIN, (e - 1) // BIN + 1):
                    bins[b].append((s, e, r.idx))
        cnt = [0] * len(recs)
        e5s, e3s = {"+": [], "-": []}, {"+": [], "-": []}
        n_prim = 0
        if c in bam.references:
            for rd in bam.fetch(c):
                if rd.flag & 0x904:
                    continue
                n_prim += 1
                hit = set()
                for bs, be in rd.get_blocks():
                    for b in range(bs // BIN, (be - 1) // BIN + 1):
                        for s, e, i in bins.get(b, ()):
                            if s < be and e > bs:
                                hit.add(i)
                for i in hit:
                    if cnt[i] < SUPPORT_READS:
                        cnt[i] += 1
                x = read_strand_ends(rd)
                if x:
                    e5s[x[0]].append(x[1])
                    e3s[x[0]].append(x[2])
        tes, tss = [], []
        for st in ("+", "-"):
            for mode, n, lo, hi in clusters(sorted(e3s[st]), TES_GAP):
                primed = internal_priming(fa, c, mode, st) if n >= CLUSTER_MIN else None
                tes.append([st, mode, n, lo, hi, primed])
            for mode, n, lo, hi in clusters(sorted(e5s[st]), TSS_GAP):
                tss.append([st, mode, n, lo, hi])
        res = {"support": cnt, "tes": tes, "tss": tss, "n_primary": n_prim,
               "n_spliced": len(e3s["+"]) + len(e3s["-"])}
        tmp = p.with_suffix(".tmp")
        tmp.write_text(json.dumps(res))
        tmp.replace(p)
        out[c] = res
        log(f"BAM pass {c}: {n_prim:,} primaries, {res['n_spliced']:,} spliced ({time.time() - t0:.0f} s)")
    return out


# ================================================================ arm GTF (gtf_loci's reading, mcl_families.rs)
def _attr(s: str, key: str):
    """mcl_families.rs gtf_loci's `attr`: the FIRST occurrence of `key "`."""
    pat = f'{key} "'
    i = s.find(pat)
    if i < 0:
        return None
    i += len(pat)
    j = s.find('"', i)
    return s[i:j] if j >= 0 else None


class Tx:
    __slots__ = ("id", "gene", "contig", "strand", "reads", "exons", "iv")

    def __init__(self, tid, gene, strand, reads):
        self.id, self.gene, self.strand, self.reads = tid, gene, strand, reads
        self.contig, self.exons, self.iv = None, [], []

    @property
    def spliced(self) -> bool:
        return len(self.exons) > 1

    def introns(self):
        """(left exon's last base, right exon's first base), 1-based, per junction."""
        return [(a[1], b[0] + 1) for a, b in zip(self.exons, self.exons[1:])]


class Locus:
    __slots__ = ("id", "contig", "txs", "rep", "strand", "iv", "rep_iv", "start1", "end1", "spliced")

    def five(self, rep=False):
        """5' end (1-based) of the span, or of the representative (None when it has no exon)."""
        if rep:
            iv = self.rep.iv
            return None if not iv else (iv[0][0] + 1 if self.strand == "+" else iv[-1][1])
        return self.start1 if self.strand == "+" else self.end1

    def three(self, rep=False):
        if rep:
            iv = self.rep.iv
            return None if not iv else (iv[-1][1] if self.strand == "+" else iv[0][0] + 1)
        return self.end1 if self.strand == "+" else self.start1


def load_arm(gtf: str, contigs: set) -> tuple[dict, dict]:
    """({tid: Tx}, {gene_id: Locus}) on `contigs`; representative = most reads, then the longer span, then the last
    transcript id (gtf_loci: sorted ids, Rust max_by_key keeps the last maximum)."""
    txs: dict = {}
    ex = collections.defaultdict(list)
    opener = assembly._open
    with opener(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in contigs:
                continue
            t = _attr(f[8], "transcript_id")
            if t is None:
                continue
            if f[2] == "transcript":
                g = _attr(f[8], "gene_id") or t
                rv = _attr(f[8], "reads")
                txs[t] = Tx(t, g, f[6], int(rv) if rv is not None and rv.isdigit() else 0)
            elif f[2] == "exon":
                ex[t].append((f[0], int(f[3]) - 1, int(f[4])))
    loci: dict = {}
    by_gene = collections.defaultdict(list)
    for t in txs.values():
        e = sorted(ex.get(t.id, []), key=lambda x: (x[1], x[2]))
        if e:
            t.contig = e[0][0]
            t.exons = [(a, b) for _, a, b in e]
            t.iv = L.merge_iv(t.exons)
        by_gene[t.gene].append(t)
    for g, ts in by_gene.items():
        withx = [t for t in ts if t.exons]
        if not withx:
            continue
        rep, best = None, None
        for t in sorted(ts, key=lambda t: t.id.encode()):
            span = (max(b for _, b in t.exons) - min(a for a, _ in t.exons)) if t.exons else 0
            k = (t.reads, span)
            if best is None or k >= best:
                rep, best = t, k
        lc = Locus()
        lc.id, lc.contig, lc.txs, lc.rep = g, withx[0].contig, withx, rep
        lc.strand = rep.strand if rep.strand in ("+", "-") else "."
        lc.iv = L.merge_iv([x for t in withx for x in t.exons])
        lc.rep_iv = rep.iv
        lc.start1, lc.end1 = lc.iv[0][0] + 1, lc.iv[-1][1]
        lc.spliced = any(t.spliced for t in withx)
        loci[g] = lc
    return txs, loci


# ================================================================ helpers on intervals and genes
class GeneIndex:
    def __init__(self, recs: list):
        self.recs = recs
        self.idx = L.Index([{"contig": r.contig, "iv": r.iv} for r in recs])

    def overlapping(self, contig: str, iv, strand: str | None):
        out = []
        for i in self.idx.near(contig, iv):
            r = self.recs[i]
            if strand is not None and r.strand != strand:
                continue
            o = L.iv_inter(iv, r.iv)
            if o > 0:
                out.append((r, o))
        return out

    def at_base(self, contig: str, p1: int, strand: str):
        """Records on `strand` whose exon union contains the 1-based base p1."""
        out = set()
        for i in self.idx.near(contig, [(p1 - 1, p1)]):
            r = self.recs[i]
            if r.strand == strand and any(s < p1 <= e for s, e in r.iv):
                out.add(r.id)
        return out


def disjoint_pairs(genes: list) -> list:
    """Pairs of records whose spans do not overlap (1-based closed)."""
    out = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            a, b = genes[i], genes[j]
            if a.end1 < b.start1 or b.end1 < a.start1:
                out.append((a, b))
    return out


def within(sorted_pos: list, x: int, tol: int) -> bool:
    i = bisect.bisect_left(sorted_pos, x - tol)
    return i < len(sorted_pos) and sorted_pos[i] <= x + tol


# ================================================================ metrics
class Row:
    COLS = ["kind", "substrate", "sample", "species", "dev", "arm", "metric", "value", "k", "n", "base_value",
            "null_value", "bar", "qualifies", "result", "refute_trigger", "note"]


def fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if math.isnan(v):
            return "NA"
        return f"{v:.6f}".rstrip("0").rstrip(".") if v != int(v) else f"{v:.1f}"
    return str(v)


class Scorer:
    def __init__(self, sub: Substrate, ann: dict, reads: dict, work: Path):
        self.sub, self.ann, self.reads, self.work = sub, ann, reads, work
        self.all_recs = [r for c in sub.contigs for r in ann.get(c, [])]
        for c in sub.contigs:
            sup = reads[c]["support"]
            for r in ann.get(c, []):
                r.support = sup[r.idx] if r.idx < len(sup) else 0
        self.genes = [r for r in self.all_recs if not r.rt]            # §3 gene set (readthrough-described out)
        self.gidx = GeneIndex(self.genes)
        self.universe = [r for r in self.genes if r.support >= SUPPORT_READS]   # fixed read-supported genes
        self.n_span_src = sum(1 for r in self.all_recs if r.src != "exon")

    # ------------------------------------------------------------ (a) + owners
    def owners(self, loci: dict) -> dict:
        own = {}
        for lc in loci.values():
            if lc.strand == "." or not lc.rep_iv:
                own[lc.id] = None
                continue
            best, bk = None, None
            for r, o in self.gidx.overlapping(lc.contig, lc.rep_iv, lc.strand):
                k = (o, L.iv_inter(lc.iv, r.iv))
                if bk is None or k > bk or (k == bk and r.id < best.id):
                    best, bk = r, k
            own[lc.id] = best
        return own

    def fused_class(self, lc: Locus, iv) -> list:
        if not lc.spliced or lc.strand == ".":
            return []
        gs = [r for r, _ in self.gidx.overlapping(lc.contig, iv, lc.strand)]
        return disjoint_pairs(gs)

    def junction_fused(self, lc: Locus) -> bool:
        for t in lc.txs:
            if not t.spliced or t.strand not in ("+", "-"):
                continue
            for p_left, p_right in t.introns():
                A = self.gidx.at_base(lc.contig, p_left, t.strand)
                B = self.gidx.at_base(lc.contig, p_right, t.strand)
                if not A or not B or (A & B):
                    continue
                for a in A:
                    for b in B:
                        ra, rb = self._rec(a), self._rec(b)
                        if ra.end1 < rb.start1 or rb.end1 < ra.start1:
                            return True
        return False

    def _rec(self, rid: str):
        if not hasattr(self, "_byid"):
            self._byid = {r.id: r for r in self.genes}
        return self._byid[rid]

    def metrics_a(self, loci: dict, own: dict, txs: dict, base_loci: dict | None) -> dict:
        m = collections.OrderedDict()
        fused = rep_fused = span_cover = jf = pc = lnc = other = unstranded = 0
        for lc in loci.values():
            if lc.spliced and lc.strand == ".":
                unstranded += 1
            pairs = self.fused_class(lc, lc.iv)
            if pairs:
                fused += 1
                if self.junction_fused(lc):
                    jf += 1
                if any(a.pc and b.pc for a, b in pairs):
                    pc += 1
                elif any(a.lnc_or_loc or b.lnc_or_loc for a, b in pairs):
                    lnc += 1
                else:
                    other += 1
            if self.fused_class(lc, lc.rep_iv):
                rep_fused += 1
            if lc.spliced and lc.strand != ".":
                span = [(lc.start1 - 1, lc.end1)]
                cov = [r for r, o in self.gidx.overlapping(lc.contig, span, lc.strand)
                       if o >= SPAN_COVER_FRAC * L.iv_len(r.iv)]
                if disjoint_pairs(cov):
                    span_cover += 1
        m["a.fused"] = (fused, None, None, "loci with a spliced transcript whose exon union overlaps >= 2 same-strand "
                        "(rep strand) genes with non-overlapping spans")
        m["a.fused_junction"] = (jf, None, None, "FUSED loci with a junction joining the exons of two such genes")
        m["a.rep_fused"] = (rep_fused, None, None, "representative's exons only")
        m["a.span_cover"] = (span_cover, None, None, "span holds >= 50% of the exon union of >= 2 such genes")
        m["a.fused_both_protein_coding"] = (pc, None, None, "a qualifying pair both protein_coding, non-LOC")
        m["a.fused_any_lncRNA_or_LOC"] = (lnc, None, None, "else a qualifying pair with a lncRNA or LOC gene")
        m["a.fused_other"] = (other, None, None, "every other FUSED locus")
        m["a.spliced_loci_unstranded_rep"] = (unstranded, None, None,
                                              "spliced loci whose representative has no strand (never FUSED)")
        owned = collections.defaultdict(list)
        for lid, r in own.items():
            if r is not None:
                owned[r.id].append(loci[lid])
        uni = set(r.id for r in self.universe)
        absorbed = absorbed_any = 0
        lidx = L.Index([{"contig": lc.contig, "iv": lc.iv, "lc": lc} for lc in loci.values()])
        for g in self.universe:
            if owned.get(g.id):
                continue
            same = anyst = False
            for i in lidx.near(g.contig, g.iv):
                lc = lidx.rows[i]["lc"]
                o = own.get(lc.id)
                if o is None or o.id == g.id or L.iv_inter(g.iv, lc.iv) <= 0:
                    continue
                anyst = True
                if lc.strand == g.strand:
                    same = True
                    break
            absorbed += same
            absorbed_any += anyst
        m["a.absorbed_genes"] = (absorbed, absorbed, len(uni), "read-supported genes owning no locus, exons overlapped "
                                 "by a same-strand locus owned by another gene (fixed universe)")
        m["a.absorbed_genes_any_strand"] = (absorbed_any, absorbed_any, len(uni), "same without the strand condition")
        m["a.loci"] = (len(loci), None, None, "")
        m["a.transcripts"] = (sum(len(lc.txs) for lc in loci.values()), None, None, "")
        if base_loci is not None:
            m["a.loci_lost_vs_base"] = (overlap_none(base_loci, loci), None, None,
                                        "BASE loci whose exons no arm locus overlaps")
            m["a.loci_gained_vs_base"] = (overlap_none(loci, base_loci), None, None,
                                          "arm loci whose exons no BASE locus overlaps")
        return m, owned

    # ------------------------------------------------------------ (b)
    def metrics_b(self, loci: dict, owned: dict) -> dict:
        m = collections.OrderedDict()
        n = len(self.universe)
        res = {}
        for lab, rep in (("span", False), ("rep", True)):
            tes = tss = 0
            for g in self.universe:
                ls = owned.get(g.id, [])
                if not ls or not g.ends:
                    continue
                t3 = [e for _, e in g.ends]
                t5 = [s for s, _ in g.ends]
                e3 = [lc.three(rep) for lc in ls if lc.three(rep) is not None]
                e5 = [lc.five(rep) for lc in ls if lc.five(rep) is not None]
                if any(abs(y - x) <= TES_TOL for y in e3 for x in t3):
                    tes += 1
                if any(abs(y - x) <= TSS_TOL for y in e5 for x in t5):
                    tss += 1
            res[lab] = (tes, tss)
        m["b.universe"] = (n, None, None, "read-supported annotated genes (>= 2 primary reads with a block on the exon "
                           "union; readthrough-described records out)")
        m["b.tes_recovered"] = (res["span"][0] / n if n else None, res["span"][0], n,
                                "owning a locus whose SPAN 3' end is within 25 bp of an annotated 3' end of the gene")
        m["b.tss_recovered"] = (res["span"][1] / n if n else None, res["span"][1], n, "span 5' end within 250 bp")
        m["b.tes_recovered_rep"] = (res["rep"][0] / n if n else None, res["rep"][0], n, "representative's own 3' end")
        m["b.tss_recovered_rep"] = (res["rep"][1] / n if n else None, res["rep"][1], n, "representative's own 5' end")
        # read-derived clusters (fixed per sample) vs locus span ends
        three = collections.defaultdict(list)
        five = collections.defaultdict(list)
        for lc in loci.values():
            if lc.strand in ("+", "-"):
                three[(lc.contig, lc.strand)].append(lc.three())
                five[(lc.contig, lc.strand)].append(lc.five())
        for v in list(three.values()) + list(five.values()):
            v.sort()
        tes_c, tss_c = collections.defaultdict(list), collections.defaultdict(list)
        n_tes = n_tes_hit = n_tss = n_tss_hit = 0
        for c in self.sub.contigs:
            for st, mode, k, lo, hi, primed in self.reads[c]["tes"]:
                if k >= CLUSTER_MIN and primed is False:
                    n_tes += 1
                    tes_c[(c, st)].append(mode)
                    n_tes_hit += within(three.get((c, st), []), mode, TES_TOL)
            for st, mode, k, lo, hi in self.reads[c]["tss"]:
                if k >= CLUSTER_MIN:
                    n_tss += 1
                    tss_c[(c, st)].append(mode)
                    n_tss_hit += within(five.get((c, st), []), mode, TSS_TOL)
        for v in list(tes_c.values()) + list(tss_c.values()):
            v.sort()
        m["b.read_tes_clusters_hit"] = (n_tes_hit / n_tes if n_tes else None, n_tes_hit, n_tes,
                                        "read-derived 3' end clusters (gap > 25 bp, >= 3 reads, not internal priming) "
                                        "whose mode is within 25 bp of a same-strand locus span 3' end")
        m["b.read_tss_clusters_hit"] = (n_tss_hit / n_tss if n_tss else None, n_tss_hit, n_tss,
                                        "5' start clusters (gap > 100 bp, >= 3 reads), mode within 250 bp of a span 5' end")
        st_loci = [lc for lc in loci.values() if lc.strand in ("+", "-")]
        h3 = sum(1 for lc in st_loci if within(tes_c.get((lc.contig, lc.strand), []), lc.three(), TES_TOL))
        h5 = sum(1 for lc in st_loci if within(tss_c.get((lc.contig, lc.strand), []), lc.five(), TSS_TOL))
        m["b.loci_span3_at_read_tes"] = (h3 / len(st_loci) if st_loci else None, h3, len(st_loci),
                                         "descriptive, prediction-conditioned")
        m["b.loci_span5_at_read_tss"] = (h5 / len(st_loci) if st_loci else None, h5, len(st_loci),
                                         "descriptive, prediction-conditioned")
        return m

    # ------------------------------------------------------------ (c)
    def multi_gene_transcripts(self, txs: dict) -> int:
        n = 0
        for t in txs.values():
            if not t.spliced or t.strand not in ("+", "-") or t.contig is None:
                continue
            gs = [r for r, _ in self.gidx.overlapping(t.contig, t.iv, t.strand)]
            if disjoint_pairs(gs):
                n += 1
        return n

    # ------------------------------------------------------------ (d)
    def reference_loci(self, cfg: dict) -> tuple[str, list, list, list]:
        """(mode, supported annotated, supported extra copies or None, every placed locus for `located`)."""
        lp = None
        try:
            lp = L.loci_path(cfg, self.sub.species)
        except Exception as e:   # noqa: BLE001 — an absent target is "not built"
            log(f"Liftoff target {self.sub.species}: {e}")
        if lp is not None:
            loci, sup, why = L.support_if_ready(cfg, self.sub.sid, self.sub.species)
            if loci is not None:
                keep = [r for r in loci if r["contig"] in self.sub.cset]
                placed = [r for r in keep if r["cls"] in ("in_place", "moved", "extra_copy")]
                ok = [r for r in placed if not r["short"] and sup.get(L.locus_key(r), 0) >= SUPPORT_READS]
                return ("liftoff", [r for r in ok if r["cls"] in ("in_place", "moved")],
                        [r for r in ok if r["cls"] == "extra_copy"], placed)
            log(f"Liftoff table merged but {why}; falling back to annotation-only")
        rows = [{"contig": r.contig, "iv": r.iv, "id": r.id} for r in self.all_recs]
        ann = [x for x, r in zip(rows, self.all_recs) if L.iv_len(r.iv) >= SHORT_BP and r.support >= SUPPORT_READS]
        return "annotation-only", ann, None, rows

    def metrics_d(self, ref: tuple, loci: dict) -> dict:
        mode, ann, extra, placed = ref
        m = collections.OrderedDict()
        rows = [{"contig": lc.contig, "iv": lc.iv, "id": lc.id} for lc in loci.values()]
        ridx = L.Index(rows)
        note = ("Liftoff table" if mode == "liftoff" else
                "annotation-only: the annotation's own gene/pseudogene records (Liftoff table not merged)")
        for lab, ref_rows in (("annotated", ann), ("extra_copy", extra)):
            if ref_rows is None:
                m[f"d.found_{lab}"] = (None, None, None, "not measured: Liftoff table not merged (extra copies)")
                continue
            k = sum(1 for g in ref_rows if ridx.best_cov(g["contig"], g["iv"])[0] >= MATCH_COV)
            m[f"d.found_{lab}"] = (k / len(ref_rows) if ref_rows else None, k, len(ref_rows),
                                   f"read-supported reference loci (exon union >= 200 bp) covered >= 50% by one "
                                   f"locus; {note}")
        gidx = L.Index(placed)
        k = sum(1 for r in rows if gidx.best_cov(r["contig"], r["iv"])[0] >= MATCH_COV)
        m["d.located"] = (k / len(rows) if rows else None, k, len(rows),
                          "loci covered >= 50% by one reference locus (descriptive)")
        sup = list(ann) + list(extra or [])
        pairs_g, pairs_r = collections.defaultdict(set), collections.defaultdict(set)
        for gi, g in enumerate(sup):
            lg = L.iv_len(g["iv"])
            for i in ridx.near(g["contig"], g["iv"]):
                r = rows[i]
                o = L.iv_inter(g["iv"], r["iv"])
                if lg and o / lg >= MATCH_COV and o / L.iv_len(r["iv"]) >= MATCH_COV:
                    pairs_g[gi].add(i)
                    pairs_r[i].add(gi)
        one = sum(1 for gi, rs in pairs_g.items() if len(rs) == 1 and len(pairs_r[next(iter(rs))]) == 1)
        m["d.reciprocal_one_to_one"] = (one, one, len(sup), "reference locus and de novo locus each >= 50% of the "
                                        "other, unique both ways (descriptive)")
        return m


def overlap_none(a: dict, b: dict) -> int:
    idx = L.Index([{"contig": lc.contig, "iv": lc.iv} for lc in b.values()])
    n = 0
    for lc in a.values():
        if not any(L.iv_inter(lc.iv, idx.rows[i]["iv"]) > 0 for i in idx.near(lc.contig, lc.iv)):
            n += 1
    return n


# ================================================================ (c) gffcompare
def gffcompare_version() -> str:
    try:
        r = subprocess.run(["gffcompare", "--version"], capture_output=True, text=True)
        return (r.stdout + r.stderr).strip().split()[-1].lstrip("v")
    except OSError:
        return "absent"


def run_gffcompare(arm: str, gtf_src: str, sub: Substrate, ref: Path, wdir: Path) -> dict:
    wdir.mkdir(parents=True, exist_ok=True)
    q = assembly.restrict_gtf(gtf_src, wdir / f"{arm}.gtf", sub.cset, force=True)
    prefix = wdir / f"gc_{arm}"
    stats = Path(str(prefix) + ".stats")
    tmap = wdir / f"gc_{arm}.{q.name}.tmap"
    figlib.run(["gffcompare", "-r", str(ref), "-o", str(prefix), str(q)], log=wdir / f"gc_{arm}.log", cwd=wdir)
    if not stats.exists() and prefix.exists():
        prefix.replace(stats)
    for ext in (".annotated.gtf", ".loci", ".tracking", ".combined.gtf"):
        Path(str(prefix) + ext).unlink(missing_ok=True)
    st = assembly.parse_stats(stats)
    n = k = 0
    matched = set()
    for r in assembly.read_tmap(tmap):
        if int(r.get("num_exons", "0") or 0) > 1:
            n += 1
            if r.get("class_code") == "=":
                k += 1
                if r.get("ref_id") not in ("", "-"):
                    matched.add(r["ref_id"])
    q.unlink(missing_ok=True)       # the restricted copy of the arm (hundreds of MB genome-wide); .tmap/.stats kept
    return {"stats": st, "prec_k": k, "prec_n": n, "matched_refs": matched}


def metrics_c(gc: dict, chains: dict, base_gc: dict | None, n_multi: int) -> dict:
    m = collections.OrderedDict()
    st = gc["stats"]
    k, n = gc["prec_k"], gc["prec_n"]
    m["c.chain_precision"] = (k / n if n else None, k, n, "query multi-exon transcripts with class '=' / query "
                              "multi-exon transcripts (.tmap)")
    mic = st.get("matching_intron_chains")
    m["c.matching_intron_chains"] = (mic, mic, st.get("ref_multiexon"), "matching reference intron chains (.stats)")
    for lev in ("intron_chain", "transcript", "locus"):
        for sp in ("sn", "pr"):
            v = st.get(f"{lev}_{sp}")
            m[f"c.gffcompare_{lev}_{sp}"] = (v / 100 if v is not None else None, None, None, "gffcompare .stats (%/100)")
    m["c.query_transcripts"] = (st.get("query_mrnas"), None, None, "")
    m["c.multi_gene_transcripts"] = (n_multi, None, None, "spliced transcripts overlapping the exons of >= 2 same-"
                                     "strand genes with non-overlapping spans")
    if base_gc is not None:
        mine = {chains.get(r) for r in gc["matched_refs"]} - {None}
        base = {chains.get(r) for r in base_gc["matched_refs"]} - {None}
        m["c.ref_chains_lost_vs_base"] = (len(base - mine), None, len(base), "distinct reference intron chains matched "
                                          "('=') in BASE and not in the arm")
        m["c.ref_chains_gained_vs_base"] = (len(mine - base), None, len(base), "")
    return m


# ================================================================ (e) families
def families_paths(prefix: str | None, gtf: str) -> dict:
    if prefix is None:
        stem = gtf[:-4] if gtf.endswith(".gtf") else gtf
        prefix = stem + ".fam"
    cl, co = Path(prefix + ".clusters.tsv"), Path(prefix + ".copies.tsv")
    return {"prefix": prefix, "clusters": cl if cl.exists() else None, "copies": co if co.exists() else None}


def hub_profile(clusters: Path, cset: set) -> tuple[int, int, int]:
    sizes = collections.Counter()
    with open(clusters) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci, ki = hdr.index("cluster_id"), hdr.index("chrom")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ki] in cset:
                sizes[f[ci]] += 1
    return len(sizes), max(sizes.values()) if sizes else 0, sum(sizes.values())


def metrics_e(cfg: dict, fs_bin: Path, sub: Substrate, arm: str, fam: dict, n_loci: int, wdir: Path) -> dict:
    m = collections.OrderedDict()
    if fam["clusters"] is None:
        m["e.status"] = (None, None, None, f"not measured: families absent ({fam['prefix']}.clusters.tsv)")
        return m
    nf, big, inf = hub_profile(fam["clusters"], sub.cset)
    m["e.families"] = (nf, None, None, "families with >= 1 locus on the substrate")
    m["e.largest_family"] = (big, None, None, "hub check (loci on the substrate)")
    m["e.loci_in_families"] = (inf, None, None, "")
    m["e.loci_entering_graph"] = (n_loci, None, None, "loci of the arm GTF on the substrate (gtf_loci)")
    import _o1_recovery as R
    cfg_fs = dict(cfg, bin=str(fs_bin))
    if sub.species == HUMAN:
        refs = []
        try:
            refs.append(("compara", R.compara_families(cfg)[0]))
        except Exception as e:  # noqa: BLE001
            m["e.compara_status"] = (None, None, None, f"not measured: Compara families: {e}")
        try:
            refs.append(("soto", R.soto_gw(cfg)))
        except Exception as e:  # noqa: BLE001
            m["e.soto_status"] = (None, None, None, f"not measured: Soto: {e}")
        genes_gff = R._o1.annotation_cache(cfg, HUMAN)["genes_gff"]
        for ref, tpath in refs:
            r = R.fs_score(cfg_fs, fam["clusters"], genes_gff, tpath, f"{arm}_{ref}", wdir / "family_score",
                           keep=sub.cset)
            note = "Compara families at Primates" if ref == "compara" else "Soto 2025 (not independent; descriptive)"
            m[f"e.{ref}_families_scored"] = (r["truth_families"], None, None, note)
            if not r["truth_families"]:
                m[f"e.{ref}_f"] = (None, None, None, f"not measured: no {ref} family on the substrate")
                continue
            if r.get("empty"):   # families scored, no cluster overlaps any: a measured zero, precision undefined
                m[f"e.{ref}_f"] = (0.0, None, None, "one-to-one bipartite F (pooled); no truth/prediction overlap")
                m[f"e.{ref}_sens"] = (0.0, 0, r["truth_genes"], "bipartite sensitivity")
                m[f"e.{ref}_prec"] = (None, 0, 0, "undefined: no predicted member matched")
                pl = [ln for ln in open(r["log"]) if "| pairwise |" in ln]
                mm = R._FS_PAIR.search(pl[-1]) if pl else None
                if mm:
                    pt, pp, tp = (int(mm.group(i)) for i in (1, 2, 3))
                    m[f"e.{ref}_pair_sens"] = (tp / pt if pt else None, tp, pt, "pairwise")
                    m[f"e.{ref}_pair_prec"] = (tp / pp if pp else None, tp, pp, "pairwise")
                continue
            m[f"e.{ref}_f"] = (r["f"], None, None, "one-to-one bipartite F (pooled)")
            m[f"e.{ref}_sens"] = (r["sens"], r["matched"], r["truth_genes"], "bipartite sensitivity")
            m[f"e.{ref}_prec"] = (r["prec"], r["matched"], r["pred_members"],
                                  "bipartite precision (upper bound; not tie-invariant, r1045)")
            m[f"e.{ref}_pair_sens"] = (r["pair_sens"], r["pair_tp"], r["pair_truth"], "pairwise")
            m[f"e.{ref}_pair_prec"] = (r["pair_prec"], r["pair_tp"], r["pair_pred"], "pairwise")
    loci, sup, why = L.support_if_ready(cfg, sub.sid, sub.species)
    if loci is None or fam["copies"] is None:
        m["e.liftoff_pair_recall"] = (None, None, None, "not measured: " + (why or f"families copy table absent "
                                                                            f"({fam['prefix']}.copies.tsv)"))
    else:
        res = L.pair_families(L.copy_pairs(loci, 0.95, sup, lambda c: c in sub.cset), L.catalog_loci(fam["copies"]))
        k = sum(1 for _, _, shared in res if shared)
        m["e.liftoff_pair_recall"] = (k / len(res) if res else None, k, len(res),
                                      "Liftoff (record, extra copy) pairs, sequence_ID >= 0.95, both read-supported, "
                                      "covering copies share a family")
    return m


# ================================================================ bar (prereg §5)
def _val(M: dict, arm: str, key: str):
    x = M.get(arm, {}).get(key)
    return None if x is None else x[0]


def _kn(M: dict, arm: str, key: str):
    x = M.get(arm, {}).get(key)
    return (None, None) if x is None else (x[1], x[2])


def clauses(M: dict, species: str, arm: str, base: str, null: str | None) -> list:
    """[(clause, value, base_value, null_value, bar, qualifies, result, refute_trigger, note, k, n)]; A1's value is
    the arm's relative FUSED reduction (k = arm FUSED, n = BASE FUSED) and its null_value NULL's reduction; every
    other clause compares the arm's measure (value) with BASE's (base_value)."""
    out = []
    # A1
    b, a = _val(M, base, "a.fused"), _val(M, arm, "a.fused")
    nv = _val(M, null, "a.fused") if null else None
    q = b is not None and b >= FLOORS["a"]
    if b is None or a is None:
        out.append(("A1", None, None, None, "", None, "not_measured", "", "FUSED absent", a, b))
    else:
        red = (b - a) / b if b else None
        nred = (b - nv) / b if (b and nv is not None) else None
        if red is None:
            res, note, q = "not_measured", "BASE FUSED = 0", None
        elif red < X_A1:
            res, note = "fail", f"reduction {red:.4f} < {X_A1}"
        elif nred is None:
            res, note = "not_measured", f"reduction {red:.4f} >= {X_A1}; NULL arm absent"
            q = q if q else None
        elif nred < 0.5 * red:
            res, note = "pass", f"reduction {red:.4f}; NULL {nred:.4f} < half"
        else:
            res, note = "fail", f"reduction {red:.4f}; NULL {nred:.4f} >= half"
        note += f"; FUSED arm {a} / BASE {b} / NULL {nv if nv is not None else 'absent'}"
        out.append(("A1", red, None, nred, f"reduction >= {X_A1} and NULL < reduction/2", q,
                    res if (q or res == "not_measured") else f"below_floor({res})", "yes" if res == "fail" else "no",
                    note, a, b))
    # G1 precision (raw counts)
    ka, na = _kn(M, arm, "c.chain_precision")
    kb, nb = _kn(M, base, "c.chain_precision")
    mb = _val(M, base, "c.matching_intron_chains")
    qc = mb is not None and mb >= FLOORS["c"]
    if None in (ka, na, kb, nb) or not na or not nb:
        out.append(("G1", None, None, None, "", None, "not_measured", "", "gffcompare absent"))
    else:
        ok = ka * nb >= kb * na
        out.append(("G1", ka / na, kb / nb, None, "arm >= BASE", qc, ("pass" if ok else "fail") if qc else
                    f"below_floor({'pass' if ok else 'fail'})", "no" if ok else "yes", "lower = refute count"))
    # G2 matching reference chains
    ma = _val(M, arm, "c.matching_intron_chains")
    if ma is None or mb is None:
        out.append(("G2", ma, mb, None, "", None, "not_measured", "", "gffcompare absent"))
    else:
        ok = 100 * ma >= 99 * mb
        trig = 100 * ma < 98 * mb
        out.append(("G2", ma, mb, None, f"arm >= BASE x (1 - {Y})", qc,
                    ("pass" if ok else "fail") if qc else f"below_floor({'pass' if ok else 'fail'})",
                    "yes" if trig else "no", "refute trigger: loss > 2Y"))
    # G3 TES-recovered genes
    ta, tb = _kn(M, arm, "b.tes_recovered")[0], _kn(M, base, "b.tes_recovered")[0]
    nu = _val(M, base, "b.universe")
    qb = nu is not None and nu >= FLOORS["b"]
    if ta is None or tb is None:
        out.append(("G3", ta, tb, None, "", None, "not_measured", "", ""))
    else:
        ok = ta >= tb
        out.append(("G3", ta, tb, None, "arm >= BASE (genes)", qb,
                    ("pass" if ok else "fail") if qb else f"below_floor({'pass' if ok else 'fail'})",
                    "no" if ok else "yes", "lower = refute count"))
    # G4 found reference loci, annotated and extra copies separately
    for part in ("annotated", "extra_copy"):
        fa, na_ = _kn(M, arm, f"d.found_{part}")
        fb, nb_ = _kn(M, base, f"d.found_{part}")
        qd = nb_ is not None and nb_ >= FLOORS["d"]
        if fa is None or fb is None:
            note = (M.get(base, {}).get(f"d.found_{part}") or (None, None, None, ""))[3]
            out.append((f"G4.{part}", fa, fb, None, "", None, "not_measured", "", note))
        else:
            ok = 100 * fa >= 99 * fb
            trig = 100 * fa < 98 * fb
            out.append((f"G4.{part}", fa, fb, None, f"arm >= BASE x (1 - {Y})", qd,
                        ("pass" if ok else "fail") if qd else f"below_floor({'pass' if ok else 'fail'})",
                        "yes" if trig else "no", (M[base][f"d.found_{part}"][3])))
    # G5 families
    if species == HUMAN:
        nfam = _val(M, base, "e.compara_families_scored")
        qe = nfam is not None and nfam >= FLOORS["e_human"]
        for meas in ("f", "sens", "prec"):
            va, vb = _val(M, arm, f"e.compara_{meas}"), _val(M, base, f"e.compara_{meas}")
            if va is None or vb is None:
                why = [x[3] for x in (M.get(arm, {}).get(f"e.compara_{meas}"), M.get(base, {}).get(f"e.compara_{meas}"),
                                      M.get(arm, {}).get("e.status"), M.get(arm, {}).get("e.compara_status")) if x]
                out.append((f"G5.compara_{meas}", va, vb, None, "", None if nfam is None else qe, "not_measured", "",
                            why[0] if why else "Compara families / arm families absent"))
                continue
            ok = va >= vb - G5_TOL - 1e-12
            trig = va < vb - G5_REFUTE - 1e-12
            out.append((f"G5.compara_{meas}", va, vb, None, f"arm >= BASE - {G5_TOL}", qe,
                        ("pass" if ok else "fail") if qe else f"below_floor({'pass' if ok else 'fail'})",
                        "yes" if trig else "no", "refute trigger: drop > 0.010"))
    else:
        va, vb = _val(M, arm, "e.liftoff_pair_recall"), _val(M, base, "e.liftoff_pair_recall")
        npair = _kn(M, base, "e.liftoff_pair_recall")[1]
        qe = npair is not None and npair >= FLOORS["e_ape"]
        if va is None or vb is None:
            note = (M.get(base, {}).get("e.liftoff_pair_recall") or M.get(arm, {}).get("e.status")
                    or (None, None, None, "families absent"))[3]
            out.append(("G5.liftoff_pair_recall", va, vb, None, "", None, "not_measured", "", note))
        else:
            ok = va >= vb - G5_TOL - 1e-12
            trig = va < vb - G5_REFUTE - 1e-12
            out.append(("G5.liftoff_pair_recall", va, vb, None, f"arm >= BASE - {G5_TOL}", qe,
                        ("pass" if ok else "fail") if qe else f"below_floor({'pass' if ok else 'fail'})",
                        "yes" if trig else "no", "refute trigger: drop > 0.010"))
    return out


# ================================================================ score
def parse_named(items: list, what: str) -> dict:
    out = collections.OrderedDict()
    for it in items or []:
        k, sep, v = it.partition("=")
        if not sep or not k or not v:
            raise SystemExit(f"--{what} expects NAME=PATH, got {it!r}")
        if k in out:
            raise SystemExit(f"--{what} {k} given twice")
        out[k] = v
    return out


def cmd_score(a) -> int:
    cfg = figlib.load_inputs()
    sub = Substrate(cfg, a.sample, a.contigs, a.drop_contigs)
    if not sub.dev and not a.heldout:
        raise SystemExit(f"{sub.label()} is not a development substrate ({DEV}); held-out substrates are scored only "
                         "after the binary's sha1 is recorded in the prereg (pass --heldout then)")
    arms = parse_named(a.arm, "arm")
    fams = parse_named(a.families, "families")
    if a.base not in arms:
        raise SystemExit(f"--arm {a.base}=<gtf> (the comparator) is required")
    if a.null and a.null not in arms:
        log(f"NULL arm {a.null!r} not given: A1's null part is 'not measured'")
    for k in fams:
        if k not in arms:
            raise SystemExit(f"--families {k}: no such arm")
    gv = gffcompare_version()
    if gv != GFFCOMPARE_VERSION:
        log(f"WARNING gffcompare {gv} != the pre-registered {GFFCOMPARE_VERSION}")
    work = Path(a.work)
    budget = Budget(a.budget_s)
    t0 = time.time()
    ann = load_annotation(sub, work)
    reads = read_pass(sub, ann, work, budget)
    S = Scorer(sub, ann, reads, work)
    subdir = work / "score" / sub.sid / sub.key
    subdir.mkdir(parents=True, exist_ok=True)
    ann_gtf = assembly.annotation_gtf(cfg, sub.sid)
    ref = assembly.restrict_gtf(ann_gtf, subdir / "ref.gtf", sub.cset)
    chains = None
    refd = S.reference_loci(cfg)
    order = [a.base] + [k for k in arms if k != a.base]
    M, GC = {}, {}
    fs_bin = Path(a.fs_bin)
    fs_ok = (fs_bin / "family_score").exists()
    if not fs_ok:
        log(f"WARNING {fs_bin}/family_score absent: (e) not measured")
    fs_sha = sha1_file(fs_bin / "family_score") if fs_ok else "absent"
    scorer_sha = sha1_file(HERE)

    def arm_key(arm):
        fam = families_paths(fams.get(arm), arms[arm])
        return key16(CACHE_VERSION, scorer_sha, fs_sha, sub.label(), fp(arms[arm]),
                     fp(fam["clusters"]) if fam["clusters"] else "-", fp(fam["copies"]) if fam["copies"] else "-",
                     fp(sub.bam), fp(sub.gff), fp(ann_gtf), refd[0], arm == a.null, a.base,
                     fp(arms[a.base]) if arm != a.base else "-")
    base_loci = None
    for arm in order:
        gtf = arms[arm]
        cache = subdir / f"{arm}.metrics.json"
        k = arm_key(arm)
        if cache.exists() and not a.no_cache:
            j = json.loads(cache.read_text())
            if j.get("key") == k:
                M[arm] = collections.OrderedDict((kk, tuple(v)) for kk, v in j["metrics"])
                GC[arm] = {"matched_refs": set(j["matched_refs"])}
                log(f"arm {arm}: cached ({cache})")
                continue
        budget.check(f"arm {arm}")
        log(f"arm {arm}: {gtf}")
        if chains is None:
            chains = {tid: (t["chrom"], t["strand"], t["introns"]) for tid, t in assembly.ref_transcripts(ref).items()
                      if t["introns"]}
        txs, loci = load_arm(gtf, sub.cset)
        if arm == a.base:
            base_loci = loci
        elif base_loci is None:
            base_loci = load_arm(arms[a.base], sub.cset)[1]
        own = S.owners(loci)
        mA, owned = S.metrics_a(loci, own, txs, base_loci if arm != a.base else None)
        m = collections.OrderedDict(mA)
        m.update(S.metrics_b(loci, owned))
        gc = run_gffcompare(arm, gtf, sub, ref, subdir / "gffcompare")
        GC[arm] = gc
        m.update(metrics_c(gc, chains, GC[a.base] if arm != a.base else None, S.multi_gene_transcripts(txs)))
        m.update(S.metrics_d(refd, loci))
        if arm == a.null:
            m["e.status"] = (None, None, None, "NULL: families not run (prereg §2)")
        elif fs_ok:
            m.update(metrics_e(cfg, fs_bin, sub, arm, families_paths(fams.get(arm), gtf), len(loci), subdir / arm))
        else:
            m["e.status"] = (None, None, None, f"not measured: {fs_bin}/family_score absent")
        M[arm] = m
        tmp = cache.with_suffix(".tmp")
        tmp.write_text(json.dumps({"key": k, "metrics": list(m.items()), "matched_refs": sorted(gc["matched_refs"])}))
        tmp.replace(cache)
        del txs, loci, own, owned
    # ---------------------------------------------------------- table
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / f"readthrough_eval.{sub.sid}.{sub.key}"
    fixed = {"a.gene_set": (len(S.genes), None, None, "gene + pseudogene records minus readthrough-described"),
             "a.records_exon_union_not_from_exons": (S.n_span_src, None, None, "records whose exon union is their CDS "
                                                     "or span (Liftoff framework's rule)"),
             "d.mode": (None, None, None, refd[0])}
    rows = []
    base_row = [sub.label(), sub.sid, sub.species, "dev (reported, never in the verdict)" if sub.dev else "held-out"]
    for k, (v, kk, n, note) in fixed.items():
        rows.append(["fixed"] + base_row + ["*", k, v, kk, n, "", "", "", "", "", "", note])
    for arm in order:
        for k, (v, kk, n, note) in M[arm].items():
            rows.append(["metric"] + base_row + [arm, k, v, kk, n, "", "", "", "", "", "", note])
    for arm in order:
        if arm in (a.base, a.null):
            continue
        for c in clauses(M, sub.species, arm, a.base, a.null if a.null in arms else None):
            cl, v, bv, nv, bar, q, res, trig, note = c[:9]
            k, n = (c[9], c[10]) if len(c) > 9 else (None, None)
            if sub.dev:
                note = (note + "; " if note else "") + "development substrate: never in the verdict"
            rows.append(["clause"] + base_row + [arm, cl, v, k, n, bv, nv, bar, q, res, trig, note])
    tmp = Path(str(stem) + ".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write(f"# readthrough_eval score | prereg {PREREG} | {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fo.write(f"# substrate {sub.label()} contigs={','.join(sub.contigs) if len(sub.contigs) <= 30 else len(sub.contigs)}"
                 f" | liftoff framework: {refd[0]} | gffcompare {gv} | family_score {fs_bin}\n")
        fo.write("\t".join(Row.COLS) + "\n")
        for r in rows:
            fo.write("\t".join(fmt(x) for x in r) + "\n")
    tmp.replace(Path(str(stem) + ".tsv"))
    prov = {"prereg": PREREG, "scorer": str(HERE), "scorer_sha1": scorer_sha, "substrate": sub.label(),
            "contigs": sub.contigs, "dev": sub.dev, "bam": fp(sub.bam), "fasta": fp(sub.fasta), "gff": fp(sub.gff),
            "annotation_gtf": fp(ann_gtf), "gffcompare": gv,
            "family_score": fp(fs_bin / "family_score") if (fs_bin / "family_score").exists() else "absent",
            "family_score_sha1": fs_sha,
            "liftoff_framework": refd[0], "arms": {k: fp(v) for k, v in arms.items()},
            "families": {k: families_paths(fams.get(k), v)["prefix"] for k, v in arms.items()},
            "base": a.base, "null": a.null if a.null in arms else None, "wall_s": round(time.time() - t0, 1)}
    Path(str(stem) + ".json").write_text(json.dumps(prov, indent=1))
    log(f"wrote {stem}.tsv ({len(rows)} rows, {time.time() - t0:.0f} s)")
    return 0


# ================================================================ null
def read_junction_table(path: str, default_contig: str | None) -> list:
    """[(contig, start1, end1, strand, rq1)] from the assembler's `<out>.readthrough_junctions.tsv` (prereg Amendment
    1 item 6: contig, donor, acceptor = the 1-based inclusive intron, strand, S, U, V1, rule r|rq1), the design's
    rt_all table (start, end, strand, ..., RQ1) or a NULL list (chrom, intron_start_1b, intron_end_1b, strand)."""
    out = []
    with open(path) as fh:
        hdr = None
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            r = dict(zip(hdr, f))
            c = r.get("chrom") or r.get("contig") or default_contig
            if c is None:
                raise SystemExit(f"{path}: no chrom column; give exactly one contig")
            s = int(r.get("intron_start_1b") or r.get("donor") or r.get("start"))
            e = int(r.get("intron_end_1b") or r.get("acceptor") or r.get("end"))
            q = (str(r.get("RQ1", r.get("rq1", ""))).strip().lower() in ("true", "1", "yes")
                 or str(r.get("rule", "")).strip().lower() == "rq1")
            out.append((c, s, e, r["strand"], q))
    return out


def log2bin(n: int) -> int:
    return int(math.floor(math.log2(n)))


def canonical(fa, contig: str, s: int, e: int, st: str) -> bool:
    d, a = fa.fetch(contig, s - 1, s + 1).upper(), fa.fetch(contig, e - 2, e).upper()
    return (d, a) in CANON[st]


def null_draw(flagged: dict, cand: dict, target: int, rng: random.Random, reads_of: dict | None = None) -> list:
    """Draw candidate junctions matched by log2 bin to the flagged ones (flagged junctions visited in a shuffled
    order, cycling; an exhausted bin falls back to the nearest non-empty bin, the lower on a tie) until the
    alignments carrying them reach `target`: DISTINCT alignments when `reads_of` {junction: read ordinals} is given
    (the last junction may overshoot by its own reads), else the sum of S. flagged / cand: {junction: S}.
    Returns [(junction, S, matched flagged junction)]."""
    if not flagged or target <= 0:
        return []
    covered: set = set()
    order = sorted(flagged)
    rng.shuffle(order)
    pools = collections.defaultdict(list)
    for j in sorted(cand):
        pools[log2bin(cand[j])].append(j)
    out, got, i = [], 0, 0
    while got < target and any(pools.values()):
        j = order[i % len(order)]
        i += 1
        b = log2bin(flagged[j])
        bins = sorted((bb for bb, v in pools.items() if v), key=lambda bb: (abs(bb - b), bb))
        p = pools[bins[0]]
        k = rng.randrange(len(p))
        p[k], p[-1] = p[-1], p[k]
        x = p.pop()
        out.append((x, cand[x], j))
        if reads_of is None:
            got += cand[x]
        else:
            covered.update(reads_of[x])
            got = len(covered)
    return out


def cmd_null(a) -> int:
    import pysam
    cfg = figlib.load_inputs()
    sub = Substrate(cfg, a.sample, a.contigs, a.drop_contigs)
    if not sub.dev and not a.heldout:
        raise SystemExit(f"{sub.label()} is not a development substrate; pass --heldout for a held-out one")
    fl = read_junction_table(a.flags, sub.contigs[0] if len(sub.contigs) == 1 else None)
    rset = collections.defaultdict(set)
    for c, s, e, st, _ in fl:
        if c in sub.cset:
            rset[c].add((s, e, st))
    targets = {}
    if a.target:
        for ln in open(a.target):
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                targets[f[0]] = int(f[1])
    bam, fa = pysam.AlignmentFile(sub.bam), pysam.FastaFile(sub.fasta)
    budget = Budget(a.budget_s)
    rows, summ = [], []
    for c in sub.contigs:
        budget.check(f"null {c}")
        reads_of = collections.defaultdict(lambda: array.array("I"))
        r_reads, ordinal = 0, 0
        for rd in (bam.fetch(c) if c in bam.references else ()):
            if rd.flag & 0x904:
                continue
            x = read_strand_ends(rd)
            if not x:
                continue
            ks = {(s, e, x[0]) for s, e in x[3]}
            for k in ks:
                reads_of[k].append(ordinal)
            ordinal += 1
            if ks & rset[c]:
                r_reads += 1
        S = {k: len(v) for k, v in reads_of.items()}
        flagged = {j: S.get(j, 0) for j in rset[c]}
        missing = sum(1 for v in flagged.values() if v < 2)
        flagged = {j: v for j, v in flagged.items() if v >= 1}
        cand = {j: n for j, n in S.items() if n >= 2 and j not in rset[c] and canonical(fa, c, *j)}
        target = targets.get(c, r_reads)
        draw = null_draw(flagged, cand, target, random.Random(f"{a.seed}:{sub.sid}:{c}"), reads_of)
        achieved = len(set().union(*(reads_of[j] for j, _, _ in draw))) if draw else 0
        for (s, e, st), n, (ms, me, mst) in sorted(draw):
            rows.append((c, s, e, st, n, log2bin(n), f"{ms}-{me}{mst}"))
        summ.append((c, len(rset[c]), missing, r_reads, target, len(draw), sum(n for _, n, _ in draw), achieved,
                     len(cand)))
        log(f"null {c}: R flags {len(rset[c])}, target {target}, drew {len(draw)} junctions, {achieved} primaries")
    with open(a.out, "w") as fo:
        fo.write(f"# NULL junctions (prereg §2) seed {a.seed}, flags {a.flags}\n")
        fo.write("chrom\tintron_start_1b\tintron_end_1b\tstrand\tS\tlog2_bin\tmatched_R_junction\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    with open(a.out + ".summary.tsv", "w") as fo:
        fo.write("contig\tR_flagged\tR_flagged_S_lt2_in_bam\ttarget_primaries\ttarget_used\tnull_junctions\t"
                 "null_sum_S\tnull_distinct_primaries\tcandidates\n")
        for r in summ:
            fo.write("\t".join(map(str, r)) + "\n")
    return 0


# ================================================================ synth (test-only arm)
def synth_gtf(gtf: str, contigs: set, drop: set, out: str) -> dict:
    """Drop transcripts carrying a junction in `drop` {(contig, s1, e1, strand)}; every gene_id that lost a
    transcript is re-split into components of transcripts sharing a junction. Returns counts."""
    lines = collections.defaultdict(list)
    txs, _ = load_arm(gtf, contigs)
    with assembly._open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 1)
            if f[0] not in contigs:
                continue
            t = _attr(line, "transcript_id")
            if t:
                lines[t].append(line)
    dropped = set()
    for t in txs.values():
        if t.contig is None:
            continue
        for a_, b_ in t.introns():
            s1, e1 = a_ + 1, b_ - 1
            if (t.contig, s1, e1, t.strand) in drop or (t.strand == "." and any(
                    (t.contig, s1, e1, x) in drop for x in "+-")):
                dropped.add(t.id)
                break
    by_gene = collections.defaultdict(list)
    for t in sorted(txs.values(), key=lambda t: t.id):
        if t.id not in dropped and t.exons:
            by_gene[t.gene].append(t)
    touched = {txs[t].gene for t in dropped}
    newgene, n_split = {}, 0
    for g, ts in by_gene.items():
        if g not in touched:          # a gene_id that lost nothing keeps its grouping
            for t in ts:
                newgene[t.id] = g
            continue
        par = {t.id: t.id for t in ts}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        jidx = {}
        for t in ts:
            for j in t.introns():
                if j in jidx:
                    ra, rb = find(t.id), find(jidx[j])
                    if ra != rb:
                        par[max(ra, rb)] = min(ra, rb)
                else:
                    jidx[j] = t.id
        comps = collections.OrderedDict()
        for t in ts:
            if t.spliced:
                comps.setdefault(find(t.id), []).append(t)
        for t in ts:
            if t.spliced:
                continue
            home = next((r for r, v in comps.items() if any(L.iv_inter(t.iv, u.iv) > 0 for u in v)), None)
            comps.setdefault(home or t.id, []).append(t)
        for i, (r, v) in enumerate(comps.items()):
            for t in v:
                newgene[t.id] = g if i == 0 else f"{g}.s{i}"
        n_split += max(0, len(comps) - 1)
    with open(out, "w") as fo:
        for t in sorted(newgene, key=lambda x: (txs[x].contig, txs[x].exons[0][0], x)):
            g = txs[t].gene
            for line in lines[t]:
                fo.write(line.replace(f'gene_id "{g}"', f'gene_id "{newgene[t]}"', 1))
    return {"transcripts_in": len(txs), "dropped": len(dropped), "genes_split": n_split}


def cmd_synth(a) -> int:
    cset = {c.strip() for c in a.contigs.split(",") if c.strip()}
    rows = read_junction_table(a.junctions, next(iter(cset)) if len(cset) == 1 else None)
    drop = {(c, s, e, st) for c, s, e, st, q in rows if (q or not a.rq1_only)}
    res = synth_gtf(a.gtf, cset, drop, a.out)
    log(f"synthetic arm {a.out}: {res} (a TEST fixture, never an arm of the prereg)")
    print(json.dumps(res))
    return 0


# ================================================================ verdict (prereg §5, across substrates)
def read_table(path: str) -> list:
    rows = []
    with open(path) as fh:
        hdr = None
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            rows.append(dict(zip(hdr, f)))
    return rows


def verdict(tables: list, npip_cap: str = "unknown") -> dict:
    """{arm: (verdict, reasons, median A1 reduction)} from `score` tables (clause rows of held-out substrates)."""
    by_arm = collections.defaultdict(dict)   # arm -> substrate -> {clause: row}
    species = {}
    for rows in tables:
        for r in rows:
            if r["kind"] != "clause" or r["dev"] != "held-out":
                continue
            by_arm[r["arm"]].setdefault(r["substrate"], {})[r["metric"]] = r
            species[r["substrate"]] = r["species"]
    out = {}
    for arm, subs in by_arm.items():
        vsubs = {s: c for s, c in subs.items() if c.get("A1", {}).get("qualifies") == "yes"}
        human = [s for s in vsubs if species[s] == HUMAN]
        ape = [s for s in vsubs if species[s] != HUMAN]
        reasons = []
        reds = sorted(float(c["A1"]["value"]) for c in vsubs.values() if c["A1"]["value"] not in ("", "NA"))
        h = len(reds) // 2
        med = None if not reds else (reds[h] if len(reds) % 2 else (reds[h - 1] + reds[h]) / 2)
        if not human or not ape:
            out[arm] = ("undecided", ["needs >= 1 human and >= 1 ape verdict substrate"], med)
            continue
        n = len(vsubs)
        a1_fail = sum(1 for c in vsubs.values() if c["A1"]["result"] == "fail")
        low = {g: sum(1 for c in vsubs.values() if c.get(g, {}).get("refute_trigger") == "yes"
                      and c[g]["qualifies"] == "yes") for g in ("G1", "G3")}
        hard = [s for s, c in vsubs.items() if any(r.get("refute_trigger") == "yes" and r["qualifies"] == "yes"
                                                   for k, r in c.items() if k.startswith(("G2", "G4", "G5")))]
        if a1_fail * 2 > n:
            reasons.append(f"A1 fails on {a1_fail}/{n} verdict substrates")
        if hard:
            reasons.append(f"G2/G4 loss > 2Y or G5 drop > 0.010 on {hard}")
        for g, k in low.items():
            if k * 2 > n:
                reasons.append(f"{g} lower on {k}/{n}")
        if reasons:
            out[arm] = ("refute", reasons, med)
            continue
        fails, unmeasured = [], []
        g5_h = g5_a = False
        for s, c in vsubs.items():
            for k, r in c.items():
                if r["result"] == "not_measured":     # "a clause not computed ... caps the verdict" (§8)
                    unmeasured.append(f"{s}:{k}")
                    continue
                if r["qualifies"] != "yes":           # below its power floor: reported, not judged there
                    continue
                if r["result"] != "pass":
                    fails.append(f"{s}:{k}")
                elif k.startswith("G5"):
                    g5_h |= species[s] == HUMAN
                    g5_a |= species[s] != HUMAN
        if fails or unmeasured or not (g5_h and g5_a):
            why = ([f"fails: {fails}"] if fails else []) + ([f"not measured: {unmeasured}"] if unmeasured else [])
            if not (g5_h and g5_a):
                why.append("G5 does not qualify on >= 1 human and >= 1 ape verdict substrate")
            out[arm] = ("keep opt-in", why, med)
        elif npip_cap != "none":
            out[arm] = ("keep opt-in", [f"NPIP cap: {npip_cap}"], med)
        else:
            out[arm] = ("adopt as default (recommendation; the flip is the user's call)", [], med)
    adopt = [k for k, v in out.items() if v[0].startswith("adopt")]
    if len(adopt) > 1:
        best = sorted(adopt, key=lambda k: (-(out[k][2] or 0), 0 if k.upper() == "RQ1" else 1))[0]
        for k in adopt:
            if k != best:
                out[k] = (out[k][0], out[k][1] + [f"not recommended: {best} has the larger median A1 reduction "
                                                   f"(RQ1 on a tie)"], out[k][2])
    return out


def cmd_verdict(a) -> int:
    res = verdict([read_table(p) for p in a.tables], a.npip_cap)
    for arm, (v, why, med) in res.items():
        print(f"{arm}\t{v}\tmedian_A1_reduction={fmt(med)}\t{'; '.join(why)}")
    return 0


# ================================================================ selftest (unit fixtures)
def _write_bam(path: Path, contig: str, length: int, reads: list):
    """reads: [(name, start0, cigar, flag, ts or None)]"""
    import pysam
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": length}]}
    tmp = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as bf:
        for name, s, cig, flag, ts in reads:
            a = pysam.AlignedSegment(bf.header)
            a.query_name, a.reference_id, a.reference_start, a.cigarstring, a.flag = name, 0, s, cig, flag
            qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIS=X])", cig))
            a.query_sequence, a.mapping_quality = "A" * qlen, 60
            if ts:
                a.set_tag("ts", ts, "A")
            bf.write(a)
    pysam.sort("-o", str(path), tmp)
    pysam.index(str(path))
    os.unlink(tmp)


def selftest() -> int:
    import pysam
    ok = 0
    # ---- genome with GT..AG introns
    ctg, n = "c1", 40_000
    seq = list("C" * n)
    for i in range(10700, 10720):          # A-rich: the 20 bases downstream of a + strand 3' end at 10700
        seq[i] = "A"
    def motif(s1, e1, st):
        d, acc = ("GT", "AG") if st == "+" else ("CT", "AC")
        seq[s1 - 1:s1 + 1] = list(d)
        seq[e1 - 2:e1] = list(acc)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # ---- annotation: A (1000-2000, 3000-3500) and B (6000-6200, 7000-7500) on +, non-overlapping spans;
        # C nested antisense; D pseudogene with exons directly under the gene; RT readthrough-described record
        gff = td / "a.gff"
        gff.write_text("\n".join([
            "##gff-version 3",
            f"{ctg}\tR\tgene\t1001\t3500\t.\t+\t.\tID=gene-A;Name=GENEA;gene_biotype=protein_coding",
            f"{ctg}\tR\tmRNA\t1001\t3500\t.\t+\t.\tID=rna-A1;Parent=gene-A",
            f"{ctg}\tR\texon\t1001\t2000\t.\t+\t.\tID=exon-A1-1;Parent=rna-A1",
            f"{ctg}\tR\texon\t3001\t3500\t.\t+\t.\tID=exon-A1-2;Parent=rna-A1",
            f"{ctg}\tR\tgene\t6001\t7500\t.\t+\t.\tID=gene-B;Name=LOC100;gene_biotype=lncRNA",
            f"{ctg}\tR\tlnc_RNA\t6001\t7500\t.\t+\t.\tID=rna-B1;Parent=gene-B",
            f"{ctg}\tR\texon\t6001\t6200\t.\t+\t.\tID=exon-B1-1;Parent=rna-B1",
            f"{ctg}\tR\texon\t7001\t7500\t.\t+\t.\tID=exon-B1-2;Parent=rna-B1",
            f"{ctg}\tR\tgene\t1001\t7500\t.\t+\t.\tID=gene-AB;Name=AB;gene_biotype=protein_coding;"
            f"description=GENEA-LOC100 readthrough",
            f"{ctg}\tR\tmRNA\t1001\t7500\t.\t+\t.\tID=rna-AB;Parent=gene-AB",
            f"{ctg}\tR\texon\t1001\t2000\t.\t+\t.\tID=exon-AB-1;Parent=rna-AB",
            f"{ctg}\tR\texon\t7001\t7500\t.\t+\t.\tID=exon-AB-2;Parent=rna-AB",
            f"{ctg}\tR\tpseudogene\t20001\t21000\t.\t-\t.\tID=gene-D;Name=DP1;gene_biotype=pseudogene",
            f"{ctg}\tR\texon\t20001\t21000\t.\t-\t.\tID=exon-D-1;Parent=gene-D",
            f"{ctg}\tR\tgene\t30001\t30100\t.\t+\t.\tID=gene-E;Name=GENEE;gene_biotype=protein_coding",
            f"{ctg}\tR\tCDS\t30001\t30100\t.\t+\t0\tID=cds-E;Parent=gene-E", ""]))
        ann = parse_gff(str(gff), {ctg})
        by = {r.id: r for r in ann[ctg]}
        assert by["gene-A"].iv == [(1000, 2000), (3000, 3500)] and by["gene-A"].ends == [(1001, 3500)]
        assert by["gene-AB"].rt and not by["gene-A"].rt
        assert by["gene-D"].iv == [(20000, 21000)] and by["gene-D"].ends == [(21000, 20001)]
        assert by["gene-E"].src == "cds" and by["gene-E"].iv == [(30000, 30100)]
        ok += 1
        # ---- genome + BAM: 3 reads of A (ending 3500), 2 read-through reads A->B, 1 secondary, B 2 reads
        motif(2001, 3000, "+")
        motif(2001, 7000, "+")
        motif(6201, 7000, "+")
        fa = td / "g.fa"
        fa.write_text(f">{ctg}\n" + "".join(seq) + "\n")
        pysam.faidx(str(fa))
        reads = [(f"a{i}", 1000, "1000M1000N500M", 0, "+") for i in range(3)]
        reads += [(f"rt{i}", 1000, "1000M5000N500M", 0, "+") for i in range(2)]
        reads += [("sec", 1000, "1000M5000N500M", 256, "+")]
        reads += [(f"b{i}", 6000, "200M800N500M", 0, "+") for i in range(2)]
        reads += [(f"q{i}", 10000, "100M500N100M", 0, "+") for i in range(3)]   # 3' end 10700: internal priming
        reads += [(f"p{i}", 20100, "500M", 16, "-") for i in range(2)]         # unspliced, on D
        bam = td / "r.bam"
        _write_bam(bam, ctg, n, reads)

        class _Sub:   # a minimal substrate (a species with no registry sample: no Liftoff target is looked up)
            sid, species, contigs, cset, dev = "t", "selftest", [ctg], {ctg}, True
            key = ctg
        sub = _Sub()
        sub.bam, sub.fasta, sub.gff = str(bam), str(fa), str(gff)
        rp = read_pass(sub, ann, td / "w", Budget(0))
        sup = {r.id: rp[ctg]["support"][r.idx] for r in ann[ctg]}
        assert sup == {"gene-A": 2, "gene-B": 2, "gene-AB": 2, "gene-D": 2, "gene-E": 0}, sup
        tes = {(t[0], t[1]): t for t in rp[ctg]["tes"]}
        assert tes[("+", 3500)][2] == 3 and tes[("+", 3500)][5] is False     # the 3 A reads
        assert tes[("+", 7500)][2] == 4 and tes[("+", 7500)][5] is False     # 2 RT + 2 B reads
        assert tes[("+", 10700)][2] == 3 and tes[("+", 10700)][5] is True    # internal priming
        tss = {(t[0], t[1]): t[2] for t in rp[ctg]["tss"]}
        assert tss == {("+", 1001): 5, ("+", 6001): 2, ("+", 10001): 3}, tss
        assert rp[ctg]["n_primary"] == 12 and rp[ctg]["n_spliced"] == 10
        ok += 1
        # ---- GTF arms: BASE fuses A and B through the RT transcript; ARM splits them
        def gtf_line(g, t, st, reads_, exons):
            s, e = exons[0][0], exons[-1][1]
            out = [f'{ctg}\tr\ttranscript\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g}"; transcript_id "{t}"; reads "{reads_}";']
            out += [f'{ctg}\tr\texon\t{a}\t{b}\t.\t{st}\t.\tgene_id "{g}"; transcript_id "{t}";' for a, b in exons]
            return out
        base = gtf_line("L1", "T1", "+", 5, [(1001, 2000), (3001, 3500)])
        base += gtf_line("L1", "T2", "+", 2, [(1001, 2000), (7001, 7510)])
        base += gtf_line("L1", "T3", "+", 2, [(6001, 6200), (7001, 7510)])
        base += gtf_line("L2", "T4", "-", 3, [(20001, 21000)])
        arm = gtf_line("L1", "T1", "+", 5, [(1001, 2000), (3001, 3500)])
        arm += gtf_line("L3", "T3", "+", 2, [(6001, 6200), (7001, 7510)])
        arm += gtf_line("L2", "T4", "-", 3, [(20001, 21000)])
        (td / "base.gtf").write_text("\n".join(base) + "\n")
        (td / "arm.gtf").write_text("\n".join(arm) + "\n")
        S = Scorer(sub, ann, rp, td / "w")
        assert {r.id for r in S.universe} == {"gene-A", "gene-B", "gene-D"}
        tx_b, lo_b = load_arm(str(td / "base.gtf"), {ctg})
        tx_a, lo_a = load_arm(str(td / "arm.gtf"), {ctg})
        assert lo_b["L1"].rep.id == "T1" and lo_b["L1"].end1 == 7510
        own_b, own_a = S.owners(lo_b), S.owners(lo_a)
        assert own_b["L1"].id == "gene-A" and own_a["L3"].id == "gene-B" and own_b["L2"].id == "gene-D"
        mb, owned_b = S.metrics_a(lo_b, own_b, tx_b, None)
        ma, owned_a = S.metrics_a(lo_a, own_a, tx_a, lo_b)
        assert mb["a.fused"][0] == 1 and ma["a.fused"][0] == 0, (mb["a.fused"], ma["a.fused"])
        assert mb["a.fused_junction"][0] == 1 and mb["a.rep_fused"][0] == 0 and mb["a.span_cover"][0] == 1
        assert mb["a.fused_any_lncRNA_or_LOC"][0] == 1
        assert mb["a.absorbed_genes"][0] == 1 and ma["a.absorbed_genes"][0] == 0      # B lost to A's locus
        assert ma["a.loci_lost_vs_base"][0] == 0 and ma["a.loci_gained_vs_base"][0] == 0
        assert S.multi_gene_transcripts(tx_b) == 1 and S.multi_gene_transcripts(tx_a) == 0
        ok += 1
        bb, ba = S.metrics_b(lo_b, owned_b), S.metrics_b(lo_a, owned_a)
        # BASE: A owns L1 whose span 3' end is 7510 (not 3500): TES miss; D: - strand, 3' = 20001: hit
        assert bb["b.tes_recovered"][1] == 1 and ba["b.tes_recovered"][1] == 3, (bb["b.tes_recovered"], ba["b.tes_recovered"])
        assert bb["b.tes_recovered_rep"][1] == 2      # A's representative T1 ends at 3500
        assert ba["b.tss_recovered"][1] == 3
        assert bb["b.read_tes_clusters_hit"][1:3] == (1, 2) and ba["b.read_tes_clusters_hit"][1:3] == (2, 2), (bb, ba)
        assert bb["b.read_tss_clusters_hit"][1:3] == (1, 2) and ba["b.read_tss_clusters_hit"][1:3] == (1, 2)
        ok += 1
        refd = S.reference_loci({})
        assert refd[0] == "annotation-only" and {r["id"] for r in refd[1]} == {"gene-A", "gene-B", "gene-AB", "gene-D"}
        db, da = S.metrics_d(refd, lo_b), S.metrics_d(refd, lo_a)
        assert db["d.found_annotated"][1:3] == (4, 4) and da["d.found_annotated"][1:3] == (4, 4)
        assert da["d.reciprocal_one_to_one"][0] == 2 and db["d.reciprocal_one_to_one"][0] == 1, (da, db)
        assert db["d.located"][1:3] == (2, 2)
        assert da["d.found_extra_copy"][0] is None
        ok += 1
        # ---- gffcompare (c) on the fixture
        ref = td / "ref.gtf"
        ref.write_text("\n".join(gtf_line("gA", "rA", "+", 0, [(1001, 2000), (3001, 3500)]) +
                                 gtf_line("gB", "rB", "+", 0, [(6001, 6200), (7001, 7500)])) + "\n")
        gcb = run_gffcompare("BASE", str(td / "base.gtf"), sub, ref, td / "gc")
        gca = run_gffcompare("ARM", str(td / "arm.gtf"), sub, ref, td / "gc")
        assert (gcb["prec_k"], gcb["prec_n"]) == (2, 3) and (gca["prec_k"], gca["prec_n"]) == (2, 2), (gcb, gca)
        chains = {tid: (t["chrom"], t["strand"], t["introns"]) for tid, t in assembly.ref_transcripts(ref).items()}
        mc = metrics_c(gca, chains, gcb, 0)
        assert mc["c.matching_intron_chains"][0] == 2 and mc["c.ref_chains_lost_vs_base"][0] == 0
        ok += 1
        # ---- bar
        M = {"BASE": collections.OrderedDict(), "R": collections.OrderedDict(), "NULL": collections.OrderedDict()}
        def put(arm_, k, v, kk=None, nn=None):
            M[arm_][k] = (v, kk, nn, "")
        for arm_, fz, pk, pn, mic, tes_, fd in (("BASE", 100, 900, 1000, 1000, 800, 1990),
                                                ("R", 88, 901, 1000, 991, 800, 1980), ("NULL", 97, 0, 1, 0, 0, 0)):
            put(arm_, "a.fused", fz)
            put(arm_, "c.chain_precision", pk / pn, pk, pn)
            put(arm_, "c.matching_intron_chains", mic, mic, 2000)
            put(arm_, "b.universe", 5000)
            put(arm_, "b.tes_recovered", tes_ / 5000, tes_, 5000)
            put(arm_, "d.found_annotated", fd / 2000, fd, 2000)
        for arm_, f_ in (("BASE", 0.80), ("R", 0.796)):
            put(arm_, "e.compara_families_scored", 40)
            for meas in ("f", "sens", "prec"):
                put(arm_, f"e.compara_{meas}", f_)
        cl = {c[0]: c for c in clauses(M, HUMAN, "R", "BASE", "NULL")}
        assert cl["A1"][6] == "pass" and abs(cl["A1"][1] - 0.12) < 1e-9 and abs(cl["A1"][3] - 0.03) < 1e-9
        assert cl["G1"][6] == "pass" and cl["G2"][6] == "pass" and cl["G3"][6] == "pass"
        assert cl["G4.annotated"][6] == "pass" and cl["G4.extra_copy"][6] == "not_measured"
        assert cl["G5.compara_f"][6] == "pass"
        M["NULL"]["a.fused"] = (93, None, None, "")
        assert {c[0]: c for c in clauses(M, HUMAN, "R", "BASE", "NULL")}["A1"][6] == "fail"
        M["R"]["c.matching_intron_chains"] = (979, 979, 2000, "")
        cl = {c[0]: c for c in clauses(M, HUMAN, "R", "BASE", None)}
        assert cl["G2"][6] == "fail" and cl["G2"][7] == "yes" and cl["A1"][6] == "not_measured"
        ok += 1
        # verdict aggregation
        def tab(sub_, sp, res):
            return [{"kind": "clause", "dev": "held-out", "arm": "R", "substrate": sub_, "species": sp,
                     "metric": k, "value": v, "qualifies": "yes", "result": r, "refute_trigger": t}
                    for k, v, r, t in res]
        allpass = [("A1", "0.2", "pass", "no"), ("G1", "", "pass", "no"), ("G2", "", "pass", "no"),
                   ("G3", "", "pass", "no"), ("G5.x", "", "pass", "no")]
        v = verdict([tab("h", HUMAN, allpass), tab("g", "gorilla", allpass)], "none")
        assert v["R"][0].startswith("adopt"), v
        assert verdict([tab("h", HUMAN, allpass)], "none")["R"][0] == "undecided"
        assert verdict([tab("h", HUMAN, allpass), tab("g", "gorilla", allpass)], "unknown")["R"][0] == "keep opt-in"
        fail = [("A1", "0.05", "fail", "yes")] + allpass[1:]
        assert verdict([tab("h", HUMAN, fail), tab("g", "gorilla", fail)], "none")["R"][0] == "refute"
        unm = tab("g", "gorilla", allpass) + [{"kind": "clause", "dev": "held-out", "arm": "R", "substrate": "g",
                                               "species": "gorilla", "metric": "G4.extra_copy", "value": "",
                                               "qualifies": "", "result": "not_measured", "refute_trigger": ""}]
        v = verdict([tab("h", HUMAN, allpass), unm], "none")["R"]
        assert v[0] == "keep opt-in" and "G4.extra_copy" in v[1][0], v          # not measured caps
        g1 = [("A1", "0.2", "pass", "no"), ("G1", "", "fail", "yes"), ("G3", "", "pass", "no"), ("G5.x", "", "pass", "no")]
        g3 = [("A1", "0.2", "pass", "no"), ("G1", "", "pass", "no"), ("G3", "", "fail", "yes"), ("G5.x", "", "pass", "no")]
        v = verdict([tab("h", HUMAN, g1), tab("g", "gorilla", g3), tab("h2", HUMAN, allpass)], "none")["R"]
        assert v[0] == "keep opt-in", v      # G1 lower on 1/3 and G3 on 1/3: neither on more than half
        v = verdict([tab("h", HUMAN, g1), tab("g", "gorilla", g1), tab("h2", HUMAN, allpass)], "none")["R"]
        assert v[0] == "refute" and "G1 lower on 2/3" in v[1][0], v
        ok += 1
        # ---- NULL draw: deterministic, bin-matched, reaches the target
        flagged = {("x", 1): 8, ("x", 2): 3}
        cand = {("c", i): (2 + i % 14) for i in range(60)}
        d1 = null_draw(flagged, cand, 11, random.Random(f"{SEED}:t:c1"))
        d2 = null_draw(flagged, cand, 11, random.Random(f"{SEED}:t:c1"))
        assert d1 == d2 and sum(n_ for _, n_, _ in d1) >= 11 and not set(j for j, _, _ in d1) & set(flagged)
        assert all(log2bin(n_) == log2bin(flagged[m]) for _, n_, m in d1)
        assert null_draw({}, cand, 5, random.Random(1)) == []
        (td / "dump.tsv").write_text("contig\tdonor\tacceptor\tstrand\tS\tU\tV1\trule\n"
                                     "c1\t2001\t7000\t+\t2\t40\t4\trq1\nc1\t6201\t7000\t+\t2\t40\t0\tr\n")
        assert read_junction_table(str(td / "dump.tsv"), None) == [("c1", 2001, 7000, "+", True),
                                                                    ("c1", 6201, 7000, "+", False)]
        ro = {j: array.array("I", range(i * 3, i * 3 + cand[j])) for i, j in enumerate(sorted(cand))}
        d3 = null_draw(flagged, cand, 11, random.Random(f"{SEED}:t:c1"), ro)
        assert len(set().union(*(ro[j] for j, _, _ in d3))) >= 11 and d3 == null_draw(
            flagged, cand, 11, random.Random(f"{SEED}:t:c1"), ro)
        ok += 1
        # ---- synth
        drop = {(ctg, 2001, 7000, "+")}
        res = synth_gtf(str(td / "base.gtf"), {ctg}, drop, str(td / "syn.gtf"))
        assert res == {"transcripts_in": 4, "dropped": 1, "genes_split": 1}, res
        _, lo_s = load_arm(str(td / "syn.gtf"), {ctg})
        assert sorted(lo_s) == ["L1", "L1.s1", "L2"], sorted(lo_s)
        ok += 1
    print(f"selftest: {ok} groups passed")
    return 0


# ================================================================ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    def substrate_args(p):
        p.add_argument("--sample", required=True, help="sample id or alias of figures/samples.tsv")
        g = p.add_mutually_exclusive_group()
        g.add_argument("--contigs", help="comma list: the substrate's contigs")
        g.add_argument("--drop-contigs", help="comma list: every annotated contig minus these")
        p.add_argument("--heldout", action="store_true", help="allow a non-development substrate (prereg §8 step 4)")
        p.add_argument("--work", default=str(DEFAULT_WORK))
        p.add_argument("--budget-s", type=float, default=0, help="exit 75 before the next heavy unit after S seconds")
    p = sp.add_parser("score", help="metrics (a)-(e) of every arm + per-clause verdict on one substrate")
    substrate_args(p)
    p.add_argument("--arm", action="append", required=True, help="NAME=GTF (repeat)")
    p.add_argument("--families", action="append", help="NAME=PREFIX (PREFIX.clusters.tsv / PREFIX.copies.tsv)")
    p.add_argument("--base", default="BASE")
    p.add_argument("--null", default="NULL")
    p.add_argument("--out", default=str(DEFAULT_WORK / "tables"))
    p.add_argument("--fs-bin", default=str(DEFAULT_FS_BIN), help="directory holding family_score")
    p.add_argument("--no-cache", action="store_true", help="recompute every arm (per-arm results are cached)")
    p = sp.add_parser("null", help="the NULL arm's junction list")
    substrate_args(p)
    p.add_argument("--flags", required=True, help="arm R's flagged junctions (Rust dump or rt_all table)")
    p.add_argument("--target", help="contig<TAB>alignments: override the per-contig target")
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--out", required=True)
    p = sp.add_parser("synth", help="a synthetic TEST arm (transcripts carrying listed junctions dropped)")
    p.add_argument("--gtf", required=True)
    p.add_argument("--contigs", required=True)
    p.add_argument("--junctions", required=True)
    p.add_argument("--rq1-only", action="store_true", help="drop only rows whose RQ1 column is true")
    p.add_argument("--out", required=True)
    p = sp.add_parser("verdict", help="§5 verdicts from score tables (held-out rows only)")
    p.add_argument("tables", nargs="+")
    p.add_argument("--npip-cap", choices=["none", "capped", "unknown"], default="unknown")
    sp.add_parser("selftest", help="unit fixtures")
    a = ap.parse_args(argv)
    return {"score": cmd_score, "null": cmd_null, "synth": cmd_synth, "verdict": cmd_verdict,
            "selftest": lambda _a: selftest()}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
