# Template-guided copy assembly — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a pre-registered Python stage that assembles a gene's copies with their own reads, guided by an expressed
locus of the same gene (the family relation inside the assembler). Judge it on Liftoff copy recall, on precision (P1, P2),
and on keeping the operating point's lead over StringTie and FLAIR.

**Architecture:** `bench/template_copies/` has six modules:
- `templates.py`: representatives and spliced FASTA;
- `candidates.py`: minimap2 hits, then uncovered, merged candidates;
- `evidence.py`: own primary reads;
- `build.py`: loci from those reads;
- `run_tg.py`: resumable stages `templates` / `project` / `assemble`;
- `score_tg.py`: Liftoff recall, P1, P2, gffcompare chains, and the tool lead through register 1268's scorer.

The pre-registration is committed after the code and before the first real projection. Five runs follow:
- development: human A119b on the shipped base and on `OP_primary`;
- held out: human testis (public, ERR13885926) and gorilla KB3781;
- tool lead: gorilla OR6737 on `OP_primary`.

**Tech Stack:** Python 3 with pysam 0.23.3, minimap2, samtools, gffcompare 0.12.10 and pytest, all installed. It reuses
`bench/guided_pipeline.py` (PAF and projection helpers), `figures/_liftoff.py` (Liftoff loci, read support, coverage,
`located`) and `bench/tool_difference/score_op.py` (register 1268).

**Spec:** `docs/superpowers/specs/2026-10-14-template-copies-design.md` (commit 8d168758, approved by the user 2026-10-14).

## Global Constraints

- **Repo and branch.** Worktree `/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset`, branch `template-copies` (from
  `main` 310dd2b7). Never push; merging is the user's decision.
- **Dependencies.** Python only, with the installed tools. No new dependency, no Rust change, and no change to the
  shipped defaults or to O1's family rule.
- **Tests.** `python3 -m pytest bench/template_copies -q -p no:cacheprovider`. The neighbours must stay green:
  `python3 -m pytest bench/liftoff_pairs bench/soto_subset bench/o1_hier bench/soto -q -p no:cacheprovider`.
- **Real-data jobs.**
  - Projection batches, `assemble`, gffcompare and the tool scorer run under `bash tools/rlock.sh heavy`, one call at a
    time, in the foreground. Each call stays under 10 minutes (Bash timeout 600000).
  - The `templates` stage is a light job (22 s and 0.6 GB on A119b).
  - Never `pkill -f`.
- **Commits.**
  - Path-limited only: `git add <paths>` and `git commit -- <paths>`.
  - Never stage the other session's `bench/hierarchy/*` or `bench/ideal_expression/*`.
  - Every message ends with `Co-Authored-By: <your model name> <noreply@anthropic.com>` and
    `Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z`.
- **Outputs.** `/mnt/linuxdisk/tmp/o1_template_copies/<run>/`, with `<run>` = `A119b_shipped`, `A119b_OP`, `testis`,
  `KB3781` or `OR6737_OP`.
  - `/mnt/linuxdisk/tmp/o1_template_copies/calib/` holds the planning calibration (one 1,000-template batch on A119b:
    89 s, 16.3 GB). Mechanics only; nothing there is scored.
- **Fixed parameters:**
  - identity ≥ 0.80 and aligned ≥ 0.50 (`guided_pipeline.MIN_ID` / `MIN_COV`);
  - coverage ≥ 0.50 (`_liftoff.MATCH_COV`);
  - ≥ 2 own primary reads (`_liftoff.SUPPORT_READS`);
  - `minimap2 -x splice -c -N 50 -p 0.1 -t 4`;
  - batches of 6,000 templates;
  - Liftoff extra copies at identity ≥ 0.95;
  - bootstrap 2,000 resamples, seed 20260925.
- **"0.2 points".** It means −0.2 for gffcompare precisions, which are percentages, and −0.002 for `located` shares,
  which are fractions.
- **Pre-registration.** `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md` is committed before the first real projection
  (Task 5). After that only its `## Results` section grows.
- **Scope.** Human and gorilla are never pooled. The public human testis library is labelled "human testis (public,
  ERR13885926)"; A119b is a testis library too (user, 2026-10-14).

## Rulings on the spec (made while planning; Task 5 copies them into the prereg)

- **R1.** Each emitted locus needs ≥ 2 own reads. The spec's ≥ 2 rule is applied to every component as well as to the
  candidate.
- **R2.** minimap2 also gets two flags:
  - `-c`, which gives the CIGARs; without them there are no projected exons;
  - `-p 0.1`, the guided finder's setting. minimap2's default 0.8 would drop secondary hits below 80% of the best
    score, which are exactly the diverged copies.
- **R3.** `assemble` fails when more than half the candidates lie on contigs absent from the BAM. That guards against
  an index and a BAM that name contigs differently.
- **R4.** Development's operating-point base is `OP_primary`, the base register 1268's B1/B2 judge. `OP_good` is not
  re-run (compute); its pre-design ceiling numbers are reported as context.
- **R5.** The tool lead reads register 1268's scorer's arm table, and `score_tg.lead` evaluates B1 and B2 for the TG
  arm. The scorer prints its own claims only for its four standard arms.
- **R6.** Batches hold 6,000 templates. One 1,000-template batch took 89 s and 16.3 GB on the human index, index load
  included.
- **R7.** P2 compares against base loci whose representative transcript has 2-5 `reads` (from `templates.tsv`). P2
  passes only when at least one locus was added.

## Review Focus

- **Index and BAM contig names differ.** The run must fail loudly, not add zero loci silently. (Task 3 test.)
- **A base locus on a contig missing from the genome FASTA.** The template stage must fail and name the contig. (Task 1
  test.)
- **An interrupted projection batch** (a PAF without its `.done` marker). `assemble` must refuse, and `project` must redo
  only that batch. (Task 3 test.)
- **A base GTF without `reads` attributes**, such as a tool's GTF. Representatives fall back to the longest transcript,
  with no crash. (Task 1 test.)
- **Re-running `assemble`.** Outputs must be byte-identical. (Task 3 test.)

## File structure

| File | Responsibility |
|---|---|
| `bench/template_copies/templates.py` | GTF reading, representatives, gene-level loci, spliced sequences, templates FASTA/TSV |
| `bench/template_copies/candidates.py` | PAF hits passing the thresholds, uncovered hits, merging duplicates |
| `bench/template_copies/evidence.py` | CIGAR blocks; the own primary reads of a candidate |
| `bench/template_copies/build.py` | the copy's loci from its own reads |
| `bench/template_copies/run_tg.py` | the CLI stages `templates`, `project` and `assemble` |
| `bench/template_copies/score_tg.py` | the CLI scores `liftoff`, `chains` and `tools` |
| `bench/template_copies/test_*.py` | one test file per pair of modules, as below |
| `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md` | the pre-registration, then the results |

---

### Task 1: templates and candidates

**Files:**
- Create: `bench/template_copies/templates.py`, `bench/template_copies/candidates.py`
- Test: `bench/template_copies/test_templates_candidates.py`

**Interfaces:**
- Consumes: `bench/guided_pipeline.py`'s `parse_paf`, `tx_exon_blocks`, `MIN_ID` and `MIN_COV`; `figures/_liftoff.py`'s
  `Index`, `iv_inter`, `iv_len` and `MATCH_COV`.
- Produces:
  - `templates.Transcript(tid, gene, contig, strand, exons, reads)`;
  - `read_gtf(path) -> {tid: Transcript}`, `representatives(tx) -> {gene: Transcript}`, `loci(tx) -> [{'contig',
    'locus', 'iv'}]`;
  - `spliced(fasta, t) -> str`, `write_templates(reps, fasta, fa, tsv)`, `read_templates(tsv) -> {tid: Transcript}`,
    `merge(iv)`;
  - `candidates.Hit(template, contig, strand, exons, identity, aligned)`, `hits(paf) -> [Hit]`,
    `uncovered(hits, base_rows) -> [Hit]`, `merge_candidates(hits) -> [(Hit, [Hit])]`.

- [ ] **Step 1: Write the failing tests**

Create `bench/template_copies/test_templates_candidates.py`:

```python
import os
import sys

import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import candidates as C  # noqa: E402
import templates as T  # noqa: E402

GTF = """chr1\tt\ttranscript\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a"; reads "5";
chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a";
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a";
chr1\tt\ttranscript\t101\t450\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b"; reads "5";
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b";
chr1\tt\texon\t301\t450\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b";
chr1\tt\ttranscript\t1001\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "g2.a"; reads "2";
chr1\tt\texon\t1001\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "g2.a";
"""


# ------------------------------------------------------------------------------------------------------- templates

def test_read_gtf_and_representatives(tmp_path):
    p = tmp_path / "base.gtf"
    p.write_text(GTF)
    tx = T.read_gtf(str(p))
    assert tx["g1.a"].exons == [(100, 200), (300, 400)] and tx["g1.a"].reads == 5 and tx["g2.a"].strand == "-"
    reps = T.representatives(tx)
    assert reps["g1"].tid == "g1.b"   # tie on reads (5): the longer exon sum (250 > 200) wins
    assert reps["g2"].tid == "g2.a"
    assert T.loci(tx) == [{"contig": "chr1", "locus": "g1", "iv": [(100, 200), (300, 450)]},
                          {"contig": "chr1", "locus": "g2", "iv": [(1000, 1100)]}]


def test_representative_tie_on_length_takes_smaller_id():
    t1 = T.Transcript("b", "g", "chr1", "+", [(0, 10)], 3)
    t2 = T.Transcript("a", "g", "chr1", "+", [(20, 30)], 3)
    assert T.representatives({"b": t1, "a": t2})["g"].tid == "a"


def test_spliced_and_round_trip(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "A" * 10 + "CCGG" + "T" * 6 + "GGTA" + "A" * 6 + "\n")
    pysam.faidx(str(fa))
    f = pysam.FastaFile(str(fa))
    plus = T.Transcript("p", "g", "chr1", "+", [(10, 14), (20, 24)], 1)
    minus = T.Transcript("m", "g", "chr1", "-", [(10, 14), (20, 24)], 1)
    assert T.spliced(f, plus) == "CCGGGGTA" and T.spliced(f, minus) == "TACCCCGG"
    T.write_templates({"g": plus}, f, str(tmp_path / "t.fa"), str(tmp_path / "t.tsv"))
    assert (tmp_path / "t.fa").read_text() == ">p\nCCGGGGTA\n"
    assert T.read_templates(str(tmp_path / "t.tsv")) == {"p": plus}


# ------------------------------------------------------------------------------------------------------ candidates

def paf_line(q, qlen, qs, qe, strand, ts, te, matches, block, cg):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, "chr1", 100000, ts, te, matches, block, 60, f"cg:Z:{cg}"])) + "\n"


def test_hits_filters_and_projects(tmp_path):
    p = tmp_path / "h.paf"
    p.write_text(paf_line("ok", 200, 0, 200, "+", 1000, 1700, 180, 200, "100M500N100M")
                 + paf_line("low_id", 200, 0, 200, "+", 3000, 3200, 150, 200, "200M")
                 + paf_line("short", 200, 0, 90, "+", 5000, 5090, 90, 90, "90M")
                 + "\t".join(map(str, ["nocg", 200, 0, 200, "+", "chr1", 100000, 7000, 7200, 200, 200, 60])) + "\n")
    hs = C.hits(str(p))
    assert [h.template for h in hs] == ["ok"] and hs[0].exons == [(1000, 1100), (1600, 1700)]
    assert hs[0].identity == 0.9 and hs[0].aligned == 1.0


def hit(t, s, e, ident=0.9, strand="+", contig="chr1"):
    return C.Hit(t, contig, strand, [(s, e)], ident, 1.0)


def test_uncovered_uses_the_half_rule():
    base = [{"contig": "chr1", "locus": "L", "iv": [(0, 600)]}]
    h1, h2 = hit("a", 0, 1000), hit("b", 0, 2000)   # base covers 600/1000 of h1, 600/2000 of h2
    assert C.uncovered([h1, h2], base) == [h2]


def test_merge_candidates():
    a, b = hit("a", 0, 1000, 0.90), hit("b", 400, 1400, 0.95)   # overlap 600 >= 500: one candidate, b kept
    c = hit("c", 900, 1900, 0.99)                               # overlaps b by 500 of 1000: joins (transitive)
    d = hit("d", 0, 1000, 0.99, strand="-")                     # other strand: separate
    e = hit("e", 1800, 3000, 0.80)                               # overlaps c by 100: separate
    groups = C.merge_candidates([a, b, c, d, e])
    kept = [(k.template, sorted(o.template for o in others)) for k, others in groups]
    assert kept == [("c", ["a", "b"]), ("e", []), ("d", [])]


def test_gtf_without_reads_attribute_takes_the_longest(tmp_path):
    p = tmp_path / "st.gtf"
    p.write_text('chr1\tStringTie\ttranscript\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
                 'chr1\tStringTie\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
                 'chr1\tStringTie\texon\t1\t300\t.\t+\t.\tgene_id "G"; transcript_id "G.2";\n')
    reps = T.representatives(T.read_gtf(str(p)))
    assert reps["G"].tid == "G.2" and reps["G"].reads == 0


def test_a_template_on_a_contig_missing_from_the_genome_fails_loudly(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    pysam.faidx(str(fa))
    with pysam.FastaFile(str(fa)) as f, pytest.raises(KeyError, match="chrZ"):
        T.spliced(f, T.Transcript("x", "g", "chrZ", "+", [(0, 10)], 1))
```

- [ ] **Step 2: Run them; they fail**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset && python3 -m pytest bench/template_copies/test_templates_candidates.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'candidates'`.

- [ ] **Step 3: Write the two modules**

Create `bench/template_copies/templates.py`:

```python
#!/usr/bin/env python3
"""Templates of the template-guided copy stage: each assembled locus's representative transcript.

Spec docs/superpowers/specs/2026-10-14-template-copies-design.md, section 1, step 1. A locus is a `gene_id` of the base
GTF. Its representative is the transcript with the most `reads` (the assembler's attribute; 0 when absent), ties broken
by the longer exon sum, then the smaller transcript id. Coordinates are 0-based half-open (GTF lines are 1-based closed).
"""
import collections
import gzip
import re

GENE = re.compile(r'gene_id "([^"]+)"')
TX = re.compile(r'transcript_id "([^"]+)"')
READS = re.compile(r'\breads "([0-9]+)"')
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

Transcript = collections.namedtuple("Transcript", "tid gene contig strand exons reads")


def merge(iv):
    """Sorted union of half-open intervals (overlapping or touching ones join)."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gtf(path):
    """GTF -> {transcript_id: Transcript}. Gene, contig and strand come from a transcript's first line; exons are merged."""
    info, reads, exons = {}, {}, collections.defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "transcript"):
                continue
            t = TX.search(f[8])
            if not t:
                continue
            tid = t.group(1)
            if tid not in info:
                g = GENE.search(f[8])
                info[tid] = (g.group(1) if g else tid, f[0], f[6])
            if f[2] == "exon":
                exons[tid].append((int(f[3]) - 1, int(f[4])))
            else:
                r = READS.search(f[8])
                reads[tid] = int(r.group(1)) if r else 0
    return {tid: Transcript(tid, g, c, s, merge(exons[tid]), reads.get(tid, 0))
            for tid, (g, c, s) in info.items() if exons.get(tid)}


def exon_len(t):
    return sum(e - s for s, e in t.exons)


def representatives(tx):
    """{gene: Transcript}: the most reads, then the longer exon sum, then the smaller transcript id."""
    best = {}
    for t in tx.values():
        b = best.get(t.gene)
        if b is None or (-t.reads, -exon_len(t), t.tid) < (-b.reads, -exon_len(b), b.tid):
            best[t.gene] = t
    return best


def loci(tx):
    """Figure 8 rows of the base loci: {'contig', 'locus', 'iv'}, iv = the union of the gene's transcripts' exons."""
    acc = {}
    for t in tx.values():
        acc.setdefault(t.gene, (t.contig, []))[1].extend(t.exons)
    return [{"contig": c, "locus": g, "iv": merge(ex)} for g, (c, ex) in sorted(acc.items())]


def spliced(fasta, t):
    """A transcript's spliced sequence in transcript orientation (`fasta` = pysam.FastaFile)."""
    seq = "".join(fasta.fetch(t.contig, s, e) for s, e in t.exons).upper()
    return seq.translate(COMP)[::-1] if t.strand == "-" else seq


def iv_str(iv):
    return ",".join(f"{s}-{e}" for s, e in iv)


def write_templates(reps, fasta, fa_path, tsv_path):
    """The templates, sorted by gene: FASTA (`>transcript_id`) and TSV (tid, gene, contig, strand, exons, reads)."""
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("tid\tgene\tcontig\tstrand\texons\treads\n")
        for g in sorted(reps):
            t = reps[g]
            fa.write(f">{t.tid}\n{spliced(fasta, t)}\n")
            tsv.write(f"{t.tid}\t{t.gene}\t{t.contig}\t{t.strand}\t{iv_str(t.exons)}\t{t.reads}\n")


def read_templates(tsv_path):
    """The templates TSV -> {tid: Transcript}."""
    out = {}
    with open(tsv_path) as fh:
        next(fh)
        for line in fh:
            tid, gene, contig, strand, ex, reads = line.rstrip("\n").split("\t")
            out[tid] = Transcript(tid, gene, contig, strand,
                                  [tuple(int(x) for x in b.split("-")) for b in ex.split(",")], int(reads))
    return out
```

Create `bench/template_copies/candidates.py`:

```python
#!/usr/bin/env python3
"""Candidates of the template-guided copy stage: projected template hits at positions no base locus covers.

Spec section 1, steps 2-3. A hit (minimap2 -x splice, PAF with `cg:Z:`) is kept when its identity (matches / block
length) is >= 0.80 and >= 50% of the template is aligned: the guided finder's registered thresholds (`guided_pipeline`
MIN_ID, MIN_COV; figure 8, G1). Its projected exons are its CIGAR's blocks split at `N` (`guided_pipeline.tx_exon_blocks`).
A hit is a candidate when every base locus covers < 50% of its projected exon union (figure 8's rule); the template's own
locus covers its own position, so self-hits drop out. Candidates on one contig and strand whose unions overlap by >= 50% of
the smaller union merge; the highest-identity hit is kept.
"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import guided_pipeline as GP  # noqa: E402
import _liftoff as LF  # noqa: E402

Hit = collections.namedtuple("Hit", "template contig strand exons identity aligned")


def hits(paf_path, min_id=GP.MIN_ID, min_cov=GP.MIN_COV):
    """PAF -> [Hit] passing the identity and aligned-fraction thresholds (hits without a CIGAR are skipped)."""
    out = []
    for h in GP.parse_paf(paf_path):
        if not h["cg"] or not h["bl"] or not h["qlen"]:
            continue
        ident, aligned = h["nm"] / h["bl"], (h["qe"] - h["qs"]) / h["qlen"]
        if ident >= min_id and aligned >= min_cov:
            out.append(Hit(h["q"], h["chrom"], h["strand"], GP.tx_exon_blocks(h), ident, aligned))
    return out


def uncovered(hit_list, base_rows, min_cov=LF.MATCH_COV):
    """The hits whose projected exon union is covered < min_cov by every base locus."""
    idx = LF.Index(base_rows)
    return [h for h in hit_list if idx.best_cov(h.contig, h.exons)[0] < min_cov]


def _rank(h):
    return (-h.identity, -h.aligned, h.template, h.contig, h.exons[0][0])


def merge_candidates(hit_list):
    """Groups of hits on one contig and strand whose exon unions overlap by >= 50% of the smaller union (transitively).
    -> [(kept Hit, [the other Hits])], kept = the highest identity, then aligned fraction, then template id; sorted by
    contig, strand and start."""
    order = sorted(range(len(hit_list)), key=lambda i: (hit_list[i].contig, hit_list[i].strand,
                                                        hit_list[i].exons[0][0], hit_list[i].exons[-1][1]))
    parent = list(range(len(hit_list)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    active = []
    for i in order:
        h = hit_list[i]
        active = [j for j in active if hit_list[j].contig == h.contig and hit_list[j].strand == h.strand
                  and hit_list[j].exons[-1][1] > h.exons[0][0]]
        for j in active:
            o = hit_list[j]
            inter = LF.iv_inter(h.exons, o.exons)
            if inter and inter >= 0.5 * min(LF.iv_len(h.exons), LF.iv_len(o.exons)):
                parent[find(i)] = find(j)
        active.append(i)
    groups = collections.defaultdict(list)
    for i in range(len(hit_list)):
        groups[find(i)].append(hit_list[i])
    out = []
    for g in groups.values():
        g = sorted(g, key=_rank)
        out.append((g[0], g[1:]))
    return sorted(out, key=lambda kg: (kg[0].contig, kg[0].strand, kg[0].exons[0][0], kg[0].template))
```

- [ ] **Step 4: Run the tests; they pass**

Run: `python3 -m pytest bench/template_copies/test_templates_candidates.py -q -p no:cacheprovider`
Expected: `8 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/template_copies/templates.py bench/template_copies/candidates.py bench/template_copies/test_templates_candidates.py
git commit -m "bench: template_copies templates (representatives, spliced FASTA) and candidates (projected hits, uncovered, merged)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/template_copies/templates.py bench/template_copies/candidates.py bench/template_copies/test_templates_candidates.py
```

---

### Task 2: own reads and the copy's loci

**Files:**
- Create: `bench/template_copies/evidence.py`, `bench/template_copies/build.py`
- Test: `bench/template_copies/test_evidence_build.py`

**Interfaces:**
- Consumes: `templates.merge` (Task 1) and `_liftoff.iv_inter`.
- Produces:
  - `evidence.blocks(reference_start, cigartuples) -> (aligned, exon_blocks)`;
  - `evidence.own_reads(bam, contig, exons) -> [(name, exon_blocks)]`;
  - `build.build(reads, proj, min_reads=2) -> [{'transcripts', 'reads', 'exons'}]`.

- [ ] **Step 1: Write the failing tests**

Create `bench/template_copies/test_evidence_build.py`:

```python
import os
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build as B  # noqa: E402
import evidence as E  # noqa: E402


# -------------------------------------------------------------------------------------------------------- evidence

def test_blocks_from_cigar():
    # 5S 10M 2I 5M 3D 4M 100N 6M 2S, from 1000
    al, ex = E.blocks(1000, [(4, 5), (0, 10), (1, 2), (0, 5), (2, 3), (0, 4), (3, 100), (0, 6), (4, 2)])
    assert al == [(1000, 1015), (1018, 1022), (1122, 1128)]
    assert ex == [(1000, 1022), (1122, 1128)]


def write_bam(path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, flag, start, cigar in sorted(reads, key=lambda r: r[2]):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = name, flag, 0, start, 60
            a.cigartuples = cigar
            a.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            out.write(a)
    pysam.index(str(path))


def test_own_reads(tmp_path):
    p = tmp_path / "r.bam"
    write_bam(p, [("prim", 0, 1000, [(0, 100), (3, 500), (0, 100)]),
                  ("second", 256, 1000, [(0, 100)]),
                  ("suppl", 2048, 1000, [(0, 100)]),
                  ("intron_only", 0, 1200, [(0, 100)]),
                  ("prim2", 16, 1550, [(0, 100)])])
    with pysam.AlignmentFile(str(p)) as bam:
        got = E.own_reads(bam, "chr1", [(1000, 1100), (1600, 1700)])
    assert got == [("prim", [(1000, 1100), (1600, 1700)]), ("prim2", [(1550, 1650)])]


# ----------------------------------------------------------------------------------------------------------- build

PROJ = [(100, 200), (300, 400), (500, 600)]


def test_simple_path_is_one_transcript():
    loci = B.build([("r1", [(110, 200), (300, 390)]), ("r2", [(320, 400), (500, 580)])], PROJ)
    assert loci == [{"transcripts": [[(110, 200), (300, 400), (500, 580)]], "reads": 2,
                     "exons": [(110, 200), (300, 400), (500, 580)]}]


def test_branching_keeps_distinct_read_chains():
    reads = [("r1", [(100, 200), (300, 400)]), ("r2", [(100, 200), (500, 600)]), ("r3", [(310, 390)])]
    loci = B.build(reads, PROJ)
    assert loci[0]["transcripts"] == [[(100, 200), (300, 400)], [(100, 200), (500, 600)]] and loci[0]["reads"] == 3


def test_a_dropped_block_breaks_the_link():
    # r1's middle block (700-750) lies outside the projection, so its first and last blocks are not linked.
    reads = [("r1", [(100, 200), (700, 750), (500, 600)]), ("r2", [(120, 200)]), ("r3", [(500, 580)])]
    loci = B.build(reads, PROJ)
    assert [lc["exons"] for lc in loci] == [[(100, 200)], [(500, 600)]] and [lc["reads"] for lc in loci] == [2, 2]


def test_a_component_with_one_read_is_not_emitted():
    loci = B.build([("r1", [(100, 200), (300, 400)]), ("r2", [(500, 560)])], PROJ)
    assert loci == []


def test_single_exon_copy():
    loci = B.build([("r1", [(100, 150)]), ("r2", [(140, 200)])], PROJ)
    assert loci == [{"transcripts": [[(100, 200)]], "reads": 2, "exons": [(100, 200)]}]
```

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/template_copies/test_evidence_build.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'build'`.

- [ ] **Step 3: Write the two modules**

Create `bench/template_copies/evidence.py`:

```python
#!/usr/bin/env python3
"""Own-read evidence of a candidate copy (spec section 1, step 4).

A read is the candidate's own when its PRIMARY alignment (`-F 2308`: not unmapped, secondary or supplementary) has an
aligned base (M / = / X) on the candidate's projected exons: figure 8's read-support rule (C2). Its exon blocks are its
alignment split at `N` (a deletion stays inside an exon), 0-based half-open.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import _liftoff as LF  # noqa: E402
from templates import merge  # noqa: E402

EXCLUDE = 0x4 | 0x100 | 0x800  # 2308
M_OPS, D_OP, N_OP = (0, 7, 8), 2, 3


def blocks(reference_start, cigartuples):
    """-> (aligned blocks: M/=/X, merged; exon blocks: split at N, D inside)."""
    pos = cur = reference_start
    aligned, exons = [], []
    for op, n in cigartuples:
        if op in M_OPS:
            aligned.append((pos, pos + n))
            pos += n
        elif op == D_OP:
            pos += n
        elif op == N_OP:
            if pos > cur:
                exons.append((cur, pos))
            pos += n
            cur = pos
    if pos > cur:
        exons.append((cur, pos))
    return merge(aligned), exons


def own_reads(bam, contig, exons):
    """The candidate's own reads in `bam` (pysam.AlignmentFile) -> [(read name, exon blocks)], sorted by name."""
    out = {}
    for r in bam.fetch(contig, exons[0][0], exons[-1][1]):
        if r.flag & EXCLUDE or r.cigartuples is None:
            continue
        al, ex = blocks(r.reference_start, r.cigartuples)
        if LF.iv_inter(al, exons) > 0:
            out[r.query_name] = ex
    return sorted(out.items())
```

Create `bench/template_copies/build.py`:

```python
#!/usr/bin/env python3
"""The copy's loci from its own reads (spec section 1, step 5).

The projection chooses the region and the reads; the reads give every coordinate. A read's exon blocks that overlap the
projected exons are kept; overlapping kept blocks of different reads merge into one observed exon; two consecutive blocks
of one read (no dropped block between them) link their observed exons. Each connected component is one locus: its
transcript is the path when the component is a simple path (links exactly between neighbouring observed exons), else
each distinct read chain is a transcript and chains contained in another are dropped. A component touched by fewer than
`min_reads` reads is not emitted (the read-support rule applied per locus; plan ruling R1). Nothing is filled from the
template and no boundary tolerance is used.
"""
import bisect
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import _liftoff as LF  # noqa: E402
from templates import merge  # noqa: E402


def _runs(blocks_, proj):
    """A read's kept blocks as runs of consecutive kept blocks: [[block, ...], ...]."""
    runs, prev = [], None
    for i, b in enumerate(blocks_):
        if LF.iv_inter([b], proj) > 0:
            if prev is not None and prev == i - 1:
                runs[-1].append(b)
            else:
                runs.append([b])
            prev = i
    return runs


def _contained(a, b):
    """True when tuple a is a contiguous sub-sequence of tuple b (a != b)."""
    return len(a) < len(b) and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))


def build(reads, proj, min_reads=2):
    """reads = [(name, exon blocks)], proj = projected exons. -> loci, sorted by start:
    [{'transcripts': [[(s, e), ...], ...], 'reads': int, 'exons': [(s, e), ...]}]."""
    runs = [(name, run) for name, bl in reads for run in _runs(bl, proj)]
    obs = merge([b for _, run in runs for b in run])
    starts = [s for s, _ in obs]
    node = lambda b: bisect.bisect_right(starts, b[0]) - 1  # noqa: E731  (the observed exon containing block b)
    parent = list(range(len(obs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edges, chains = set(), []
    for name, run in runs:
        chain = []
        for b in run:
            n = node(b)
            if chain and chain[-1] != n:
                edges.add((chain[-1], n))
                parent[find(n)] = find(chain[-1])
            if not chain or chain[-1] != n:
                chain.append(n)
        chains.append((name, tuple(chain)))
    comp_nodes, comp_reads, comp_chains = (collections.defaultdict(set) for _ in range(3))
    for n in range(len(obs)):
        comp_nodes[find(n)].add(n)
    for name, chain in chains:
        c = find(chain[0])
        comp_reads[c].add(name)
        comp_chains[c].add(chain)
    out = []
    for c, nodes in comp_nodes.items():
        if len(comp_reads[c]) < min_reads:
            continue
        ordered = sorted(nodes)
        cedges = {e for e in edges if find(e[0]) == c}
        if cedges == {(ordered[i], ordered[i + 1]) for i in range(len(ordered) - 1)}:
            txs = [ordered]
        else:
            cs = sorted(comp_chains[c])
            txs = [list(ch) for ch in cs if not any(_contained(ch, o) for o in cs)]
        out.append({"transcripts": [[obs[n] for n in t] for t in txs], "reads": len(comp_reads[c]),
                    "exons": [obs[n] for n in ordered]})
    return sorted(out, key=lambda lc: lc["exons"][0])
```

- [ ] **Step 4: Run the tests; they pass**

Run: `python3 -m pytest bench/template_copies/test_evidence_build.py -q -p no:cacheprovider`
Expected: `7 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/template_copies/evidence.py bench/template_copies/build.py bench/template_copies/test_evidence_build.py
git commit -m "bench: template_copies own primary reads and the copy's loci from them (no template fill, no tolerance)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/template_copies/evidence.py bench/template_copies/build.py bench/template_copies/test_evidence_build.py
```

---

### Task 3: the runner

**Files:**
- Create: `bench/template_copies/run_tg.py`
- Test: `bench/template_copies/test_run_tg.py`

**Interfaces:**
- Consumes: Tasks 1-2.
- Produces: `run_tg.main(argv)` with the stages:
  - `templates`: writes `DIR/templates.fa`, `DIR/templates.tsv` and `DIR/batches/bNNNN.fa`, and prints the counts as
    JSON;
  - `project`: writes `DIR/batches/bNNNN.paf` and `.done`;
  - `assemble`: writes `DIR/tg.gtf`, `DIR/merged.gtf` and `DIR/report.json`, whose keys are `hits`, `candidates`,
    `merged_candidates`, `contig_absent`, `few_reads`, `no_locus`, `with_loci`, `loci` and `transcripts`.
- `TG_<label>_NNNNNN` loci carry the attributes `template`, `template_gene`, `identity`, `reads`, `also_templates` and
  `source "template_copies"`.

- [ ] **Step 1: Write the failing tests**

Create `bench/template_copies/test_run_tg.py`. It runs real minimap2 on a 30 kb toy genome holding gene A and an exact
copy of it:

```python
import json
import os
import random
import sys

import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_tg as R  # noqa: E402

E1, I1, E2 = 300, 500, 300          # exon, intron, exon lengths of gene A
A0, C0 = 1000, 20000                # gene A at 1000; its exact copy at 20000


def toy_genome(path):
    rng = random.Random(3)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    gene = rand(E1) + "GT" + rand(I1 - 4) + "AG" + rand(E2)
    seq = rand(A0) + gene + rand(C0 - A0 - len(gene)) + gene + rand(30000 - C0 - len(gene))
    path.write_text(">chr1\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    pysam.faidx(str(path))
    return seq


def toy_inputs(tmp_path):
    seq = toy_genome(tmp_path / "g.fa")
    a = [(A0, A0 + E1), (A0 + E1 + I1, A0 + E1 + I1 + E2)]
    (tmp_path / "base.gtf").write_text(
        f'chr1\tt\ttranscript\t{a[0][0] + 1}\t{a[1][1]}\t.\t+\t.\tgene_id "A"; transcript_id "A.1"; reads "5";\n'
        + "".join(f'chr1\tt\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "A"; transcript_id "A.1";\n' for s, e in a))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
    reads = [("c1", 0, C0, [(0, E1), (3, I1), (0, E2)]), ("c2", 0, C0 + 50, [(0, E1 - 50), (3, I1), (0, 200)]),
             ("c3", 0, C0 + 100, [(0, E1 - 100), (3, I1), (0, E2)]), ("a1", 0, A0, [(0, E1), (3, I1), (0, E2)])]
    with pysam.AlignmentFile(str(tmp_path / "r.bam"), "wb", header=header) as out:
        for name, flag, start, cig in sorted(reads, key=lambda r: r[2]):
            x = pysam.AlignedSegment(out.header)
            x.query_name, x.flag, x.reference_id, x.reference_start, x.mapping_quality = name, flag, 0, start, 60
            x.cigartuples = cig
            x.query_sequence = "A" * sum(n for op, n in cig if op == 0)
            out.write(x)
    pysam.index(str(tmp_path / "r.bam"))
    return tmp_path


def test_end_to_end_on_a_toy_copy(tmp_path):
    d = toy_inputs(tmp_path)
    out = str(d / "tg")
    assert R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out]) == 0
    assert R.main(["project", "--out", out, "--index", str(d / "g.fa"), "--threads", "1"]) == 0
    assert R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out,
                   "--label", "toy"]) == 0
    rep = json.load(open(os.path.join(out, "report.json")))
    assert rep["merged_candidates"] == 1 and rep["loci"] == 1 and rep["transcripts"] == 1
    exons = [ln.split("\t")[3:5] for ln in open(os.path.join(out, "tg.gtf")) if ln.split("\t")[2] == "exon"]
    assert exons == [[str(C0 + 1), str(C0 + E1)], [str(C0 + E1 + I1 + 1), str(C0 + E1 + I1 + E2)]]
    assert 'template "A.1"' in open(os.path.join(out, "tg.gtf")).read()
    merged = open(os.path.join(out, "merged.gtf")).read()
    assert merged.startswith(open(d / "base.gtf").read()) and "TG_toy_000000" in merged


def test_assemble_refuses_unprojected_batches(tmp_path):
    d = toy_inputs(tmp_path)
    out = str(d / "tg")
    R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out])
    with pytest.raises(SystemExit, match="missing"):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out, "--label", "x"])


def run_until_projected(d):
    out = str(d / "tg")
    R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out])
    R.main(["project", "--out", out, "--index", str(d / "g.fa"), "--threads", "1"])
    return out


def test_contig_names_that_disagree_fail_loudly(tmp_path):
    d = toy_inputs(tmp_path)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrZ", "LN": 1000}]}
    with pysam.AlignmentFile(str(d / "z.bam"), "wb", header=header):
        pass
    pysam.index(str(d / "z.bam"))
    out = run_until_projected(d)
    with pytest.raises(SystemExit, match="name contigs differently"):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "z.bam"), "--out", out, "--label", "t"])


def test_assemble_is_deterministic(tmp_path):
    d = toy_inputs(tmp_path)
    out = run_until_projected(d)
    runs = []
    for _ in range(2):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out, "--label", "t"])
        runs.append((open(os.path.join(out, "tg.gtf")).read(), open(os.path.join(out, "report.json")).read()))
    assert runs[0] == runs[1] and runs[0][0]
```

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/template_copies/test_run_tg.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'run_tg'`.

- [ ] **Step 3: Write the runner**

Create `bench/template_copies/run_tg.py`:

```python
#!/usr/bin/env python3
"""Template-guided copy assembly for one base GTF (spec docs/superpowers/specs/2026-10-14-template-copies-design.md).

    run_tg.py templates --base BASE.gtf --genome FA --out DIR [--batch-size 5000]
    run_tg.py project   --out DIR --index MMI_OR_FA [--batches all|0,3,...] [--threads 4]
    run_tg.py assemble  --base BASE.gtf --bam BAM --out DIR --label NAME

The stages are separate calls, so each stays under the 10-minute limit for one foreground job:
  templates  DIR/templates.fa + .tsv (one representative per locus) and DIR/batches/bNNNN.fa;
  project    minimap2 -x splice -c -N 50 -p 0.1 of each named batch -> DIR/batches/bNNNN.paf + .done
             (a batch with a .done marker is skipped; -c gives the CIGAR the projection needs; -p 0.1 is the guided
             finder's setting, plan ruling R2);
  assemble   every batch must be .done; hits -> candidates -> own reads -> loci; writes DIR/tg.gtf (the added loci),
             DIR/merged.gtf (base + added) and DIR/report.json (counts at every step).
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build as B  # noqa: E402
import candidates as C  # noqa: E402
import evidence as E  # noqa: E402
import templates as T  # noqa: E402

MINIMAP2 = ["minimap2", "-x", "splice", "-c", "-N", "50", "-p", "0.1"]


def stage_templates(base, genome, out, batch_size):
    tx = T.read_gtf(base)
    reps = T.representatives(tx)
    os.makedirs(os.path.join(out, "batches"), exist_ok=True)
    with pysam.FastaFile(genome) as fa:
        T.write_templates(reps, fa, os.path.join(out, "templates.fa"), os.path.join(out, "templates.tsv"))
    n, b, fh = 0, -1, None
    with open(os.path.join(out, "templates.fa")) as src:
        for line in src:
            if line.startswith(">"):
                if n % batch_size == 0:
                    if fh:
                        fh.close()
                    b += 1
                    fh = open(os.path.join(out, "batches", f"b{b:04d}.fa"), "w")
                n += 1
            fh.write(line)
    if fh:
        fh.close()
    rep = dict(transcripts=len(tx), loci=len(reps), templates=n, batches=b + 1, batch_size=batch_size)
    print(json.dumps(rep))
    return rep


def batch_files(out):
    return sorted(glob.glob(os.path.join(out, "batches", "b*.fa")))


def stage_project(out, index, which, threads):
    fas = batch_files(out)
    pick = range(len(fas)) if which == "all" else [int(x) for x in which.split(",")]
    for i in pick:
        fa = fas[i]
        paf, done = fa[:-3] + ".paf", fa[:-3] + ".done"
        if os.path.exists(done):
            continue
        with open(paf + ".tmp", "w") as fh:
            subprocess.run(MINIMAP2 + ["-t", str(threads), index, fa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.replace(paf + ".tmp", paf)
        open(done, "w").close()
        print(f"projected {os.path.basename(fa)}")


def gtf_lines(label, i, kept, others, gene_of, locus):
    gid = f"TG_{label}_{i:06d}"
    attrs = (f'template "{kept.template}"; template_gene "{gene_of.get(kept.template, "")}"; '
             f'identity "{kept.identity:.4f}"; reads "{locus["reads"]}"; also_templates "{len(others)}"; '
             f'source "template_copies";')
    out = []
    for k, ex in enumerate(locus["transcripts"]):
        tid = f"{gid}.{k + 1}"
        head = f"{kept.contig}\ttemplate_copies\t"
        out.append(f'{head}transcript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{kept.strand}\t.\tgene_id "{gid}"; '
                   f'transcript_id "{tid}"; {attrs}\n')
        out += [f'{head}exon\t{s + 1}\t{e}\t.\t{kept.strand}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n'
                for s, e in ex]
    return out


def stage_assemble(base, bam_path, out, label):
    fas = batch_files(out)
    missing = [f for f in fas if not os.path.exists(f[:-3] + ".done")]
    if not fas or missing:
        raise SystemExit(f"assemble needs every batch projected; missing: {[os.path.basename(m) for m in missing]}")
    tx = T.read_gtf(base)
    base_rows = T.loci(tx)
    gene_of = {t.tid: t.gene for t in T.read_templates(os.path.join(out, "templates.tsv")).values()}
    hs = [h for f in fas for h in C.hits(f[:-3] + ".paf")]
    cand = C.uncovered(hs, base_rows)
    groups = C.merge_candidates(cand)
    rep = dict(hits=len(hs), candidates=len(cand), merged_candidates=len(groups), contig_absent=0, few_reads=0,
               no_locus=0, with_loci=0, loci=0, transcripts=0)
    lines = []
    with pysam.AlignmentFile(bam_path) as bam:
        contigs = set(bam.references)
        for kept, others in groups:
            if kept.contig not in contigs:
                rep["contig_absent"] += 1
                continue
            reads = E.own_reads(bam, kept.contig, kept.exons)
            if len(reads) < 2:
                rep["few_reads"] += 1
                continue
            loci = B.build(reads, kept.exons)
            if not loci:
                rep["no_locus"] += 1
                continue
            rep["with_loci"] += 1
            for locus in loci:
                lines += gtf_lines(label, rep["loci"], kept, others, gene_of, locus)
                rep["loci"] += 1
                rep["transcripts"] += len(locus["transcripts"])
    if groups and rep["contig_absent"] > len(groups) / 2:
        raise SystemExit(f"{rep['contig_absent']} of {len(groups)} candidates lie on contigs absent from {bam_path}: "
                         "the index and the BAM name contigs differently")
    with open(os.path.join(out, "tg.gtf"), "w") as fh:
        fh.writelines(lines)
    with open(os.path.join(out, "merged.gtf"), "w") as fh, open(base) as src:
        for line in src:
            fh.write(line)
        fh.writelines(lines)
    with open(os.path.join(out, "report.json"), "w") as fh:
        json.dump(dict(base=base, bam=bam_path, label=label, **rep), fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps(rep))
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("templates")
    p.add_argument("--base", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-size", type=int, default=5000)
    p = sub.add_parser("project")
    p.add_argument("--out", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--batches", default="all")
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("assemble")
    p.add_argument("--base", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--label", required=True)
    a = ap.parse_args(argv)
    if a.stage == "templates":
        stage_templates(a.base, a.genome, a.out, a.batch_size)
    elif a.stage == "project":
        stage_project(a.out, a.index, a.batches, a.threads)
    else:
        stage_assemble(a.base, a.bam, a.out, a.label)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests; they pass**

Run: `python3 -m pytest bench/template_copies/test_run_tg.py -q -p no:cacheprovider`
Expected: `4 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/template_copies/run_tg.py bench/template_copies/test_run_tg.py
git commit -m "bench: template_copies runner (resumable templates / project / assemble; contig-name guard)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/template_copies/run_tg.py bench/template_copies/test_run_tg.py
```

---

### Task 4: the scorer

**Files:**
- Create: `bench/template_copies/score_tg.py`
- Test: `bench/template_copies/test_score_tg.py`

**Interfaces:**
- Consumes: `_liftoff`'s `read_loci`, `read_support`, `denovo_loci`, `Index`, `locus_key`, `MATCH_COV` and
  `SUPPORT_READS`; `templates.read_templates`; `bench/tool_difference/score_op.py` (a subprocess).
- Produces, used by Tasks 6-10:
  - `score_tg.main(argv)` with the commands `liftoff`, `chains` and `tools`, each writing JSON;
  - `liftoff` keys: `recall` {copies, base, merged, reachable, reached}, `p1_located` {base, merged, passes} and `p2`
    {tg_loci, tg_located, tg_ci95, base_low_loci, base_low_located, passes};
  - `chains` keys: {name: {chain_sn, chain_pr}};
  - `tools` keys: `table`, `chain` {arm: [sn, pr]}, `b1` and `b2`.

- [ ] **Step 1: Write the failing tests**

Create `bench/template_copies/test_score_tg.py`:

```python
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_tg as S  # noqa: E402


def lrow(src, cls, s, e, short=False, seqid=1.0):
    return {"source_id": src, "cls": cls, "contig": "chr1", "start0": str(s), "end": str(e), "iv": [(s, e)],
            "short": short, "seqid_f": seqid}


def test_liftoff_scores_by_hand():
    loci = [lrow("S1", "in_place", 0, 1000), lrow("S1", "extra_copy", 5000, 6000), lrow("S1", "extra_copy", 9000, 10000),
            lrow("S2", "in_place", 20000, 21000), lrow("S2", "extra_copy", 30000, 31000, seqid=0.9)]
    sup = {S.LF.locus_key(r): 3 for r in loci}
    base = [{"contig": "chr1", "locus": "b1", "iv": [(0, 1000)]}, {"contig": "chr1", "locus": "b2", "iv": [(40000, 41000)]}]
    tg = [{"contig": "chr1", "locus": "t1", "iv": [(5000, 6000)]}]
    res = S.liftoff_scores(loci, sup, base, base + tg, tg, {"b1": 5, "b2": 9})
    # copies at >= 0.95: S1's two; base covers neither; S1's placement is assembled, so both are reachable; one reached
    assert res["recall"] == {"copies": 2, "base": 0, "merged": 1, "reachable": 2, "reached": 1}
    assert res["p1_located"] == {"base": 0.5, "merged": 2 / 3, "passes": True}
    assert res["p2"]["tg_located"] == 1.0 and res["p2"]["base_low_loci"] == 1 and res["p2"]["base_low_located"] == 1.0
    assert res["p2"]["passes"] is True


TABLE = """[human] genome-wide minus chr20 (sensitivity / precision)
arm                     base        intron         chain    transcript         locus
OP_primary_TG   50.6 /  18.2  66.7 /  56.4  23.9 /  15.5  22.0 /  14.0  36.1 /  23.5
StringTie       59.6 /  14.4  63.5 /  57.8  18.1 /  15.0  16.9 /  12.7  33.7 /  20.1
FLAIR           66.4 /   8.8  67.7 /  30.1  22.4 /   6.2  21.1 /   3.6  31.9 /   6.3
default_primary  48.1 /  17.4  60.8 /  64.4  21.8 /  18.5  20.3 /  16.1  34.5 /  22.8
stratum S1: 926 multi-exon reference transcripts in 304 genes
   OP_primary_TG  matched     0 / 926   precision 0 / 0 = 0.000
"""


def test_parse_score_op():
    arms = S.parse_score_op(TABLE)
    assert arms == {"OP_primary_TG": (23.9, 15.5), "StringTie": (18.1, 15.0), "FLAIR": (22.4, 6.2),
                    "default_primary": (21.8, 18.5)}


def test_parse_stats(tmp_path):
    p = tmp_path / "x.stats"
    p.write_text("#= Summary\n        Intron level:    66.7     |    56.4    |\nIntron chain level:    23.6     |    15.6    |\n")
    assert S.parse_stats(str(p)) == {"intron_sn": 66.7, "intron_pr": 56.4, "intron_chain_sn": 23.6, "intron_chain_pr": 15.6}


def test_interval_is_seeded():
    assert S.interval([1, 0, 1, 1]) == S.interval([1, 0, 1, 1]) and S.interval([]) != S.interval([])  # NaN != NaN


def test_lead():
    arms = S.parse_score_op(TABLE)
    assert S.lead(arms, "OP_primary_TG") == {"b1": True, "b2": True}     # 23.9 > 18.1 and > 22.4; 15.5 > 6.2
    assert S.lead({"OP_primary_TG": (23.9, 15.5)}, "OP_primary_TG") == {"b1": None, "b2": None}
    assert S.lead(dict(arms, FLAIR=(24.0, 6.2)), "OP_primary_TG")["b1"] is False
```

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/template_copies/test_score_tg.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'score_tg'`.

- [ ] **Step 3: Write the scorer**

Create `bench/template_copies/score_tg.py`:

```python
#!/usr/bin/env python3
"""Scores of the template-guided copy stage (spec section 2).

    score_tg.py liftoff --species S --sample ID --base BASE.gtf --merged MERGED.gtf --tg TG.gtf --templates T.tsv --out J
    score_tg.py chains  --ref REF.gtf[.gz] --work DIR --gtf NAME=PATH [--gtf ...] --out J
    score_tg.py tools   --species S --ref REF --exclude CHROM --strata PREFIX --work DIR --arm NAME=GTF ... --tg-arm NAME --out J

liftoff  recall of read-supported Liftoff extra copies (identity >= 0.95), base vs merged, and how many template-reachable
         copies (the copy unassembled, its source placement assembled) were reached; P1's located share, base vs
         merged; P2: the added loci's located share vs base loci whose representative has 2-5 reads, with a 95% interval
         over loci. Located = covered >= 50% by one placed Liftoff locus (figure 8's `located`).
chains   gffcompare 0.12.10 of each GTF against the reference: intron-chain sensitivity / precision (figure 1's rule).
tools    register 1268's scorer (bench/tool_difference/score_op.py) on the arms; B1 (chain sensitivity of the TG arm >
         StringTie's and > FLAIR's) and B2 (its chain precision > FLAIR's).
"""
import argparse
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "figures"))
import _liftoff as LF  # noqa: E402
import templates as T  # noqa: E402

LIFTOFF = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
SC = 0.95
BOOT_N, BOOT_SEED = 2000, 20260925
LEVEL = re.compile(r"^\s*(Base|Exon|Intron|Intron chain|Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)")
ROW = re.compile(r"^(\S+)\s+" + r"\s+".join([r"([\d.]+)\s*/\s*([\d.]+)"] * 5) + r"\s*$")


def placed_index(loci):
    return LF.Index([r for r in loci if r["cls"] in ("in_place", "moved", "extra_copy")])


def located(rows, gidx):
    """Per row (a GTF locus): covered >= 50% by one placed Liftoff locus."""
    return [gidx.best_cov(r["contig"], r["iv"])[0] >= LF.MATCH_COV for r in rows]


def gtf_rows(path, cache):
    return LF.denovo_loci(Path(path), Path(cache))


def interval(hits, n=BOOT_N, seed=BOOT_SEED):
    """95% percentile interval of the mean of 0/1 hits, resampling loci."""
    if not hits:
        return [float("nan"), float("nan")]
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(hits) for _ in hits) / len(hits) for _ in range(n))
    return [means[int(0.025 * n)], means[int(0.975 * n) - 1]]


def liftoff_scores(loci, sup, base_rows, merged_rows, tg_rows, base_reads):
    """loci, sup: figure 8's Liftoff table and the sample's read support; *_rows: GTF loci (`locus` = gene_id);
    base_reads: {gene_id: representative reads}."""
    copies = [r for r in loci if r["cls"] == "extra_copy" and not r["short"] and r["seqid_f"] >= SC - 1e-9
              and sup.get(LF.locus_key(r), 0) >= LF.SUPPORT_READS]
    placement = {r["source_id"]: r for r in loci if r["cls"] in ("in_place", "moved") and not r["short"]}
    bidx, midx = LF.Index(base_rows), LF.Index(merged_rows)
    cov = lambda idx, r: idx.best_cov(r["contig"], r["iv"])[0] >= LF.MATCH_COV  # noqa: E731
    in_base = [cov(bidx, c) for c in copies]
    in_merged = [cov(midx, c) for c in copies]
    reachable = [i for i, c in enumerate(copies) if not in_base[i] and c["source_id"] in placement
                 and cov(bidx, placement[c["source_id"]])]
    gidx = placed_index(loci)
    lb, lm, lt = located(base_rows, gidx), located(merged_rows, gidx), located(tg_rows, gidx)
    low = [h for r, h in zip(base_rows, lb) if 2 <= base_reads.get(r["locus"], 0) <= 5]
    share = lambda xs: sum(xs) / len(xs) if xs else float("nan")  # noqa: E731
    return dict(
        recall=dict(copies=len(copies), base=sum(in_base), merged=sum(in_merged),
                    reachable=len(reachable), reached=sum(in_merged[i] for i in reachable)),
        p1_located=dict(base=share(lb), merged=share(lm), passes=bool(share(lm) >= share(lb) - 0.002)),
        p2=dict(tg_loci=len(lt), tg_located=share(lt), tg_ci95=interval([int(x) for x in lt]),
                base_low_loci=len(low), base_low_located=share(low),
                passes=bool(lt) and share(lt) >= share(low)))


def parse_stats(path):
    out = {}
    for line in open(path):
        m = LEVEL.match(line)
        if m:
            k = m.group(1).lower().replace(" ", "_")
            out[f"{k}_sn"], out[f"{k}_pr"] = float(m.group(2)), float(m.group(3))
    return out


def chains(ref, work, gtfs):
    os.makedirs(work, exist_ok=True)
    out = {}
    for name, path in gtfs:
        pre = os.path.join(work, f"gc_{name}")
        if not os.path.exists(pre + ".stats"):
            subprocess.run(["gffcompare", "-r", ref, "-o", pre, path], check=True, cwd=work,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        st = parse_stats(pre + ".stats")
        out[name] = dict(chain_sn=st["intron_chain_sn"], chain_pr=st["intron_chain_pr"])
    return out


def parse_score_op(text):
    """score_op.py's arm table -> {arm: (chain_sn, chain_pr)}."""
    out = {}
    for line in text.splitlines():
        m = ROW.match(line)
        if m and m.group(1) != "arm":
            v = [float(x) for x in m.groups()[1:]]
            out[m.group(1)] = (v[4], v[5])
    return out


def tools(args):
    cmd = [sys.executable, os.path.join(REPO, "bench", "tool_difference", "score_op.py"), "--species", args.species,
           "--ref", args.ref, "--exclude", args.exclude, "--strata", args.strata, "--work", args.work]
    for a in args.arm:
        cmd += ["--arm", a]
    os.makedirs(args.work, exist_ok=True)
    text = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    arms = parse_score_op(text)
    return dict(table=text, chain={k: list(v) for k, v in arms.items()}, **lead(arms, args.tg_arm))


def lead(arms, tg_arm):
    """B1 / B2 of register 1268 for the TG arm; None when an arm is missing (a cache-filling call with one arm)."""
    if not {tg_arm, "StringTie", "FLAIR"} <= set(arms):
        return dict(b1=None, b2=None)
    tg, st, fl = arms[tg_arm], arms["StringTie"], arms["FLAIR"]
    return dict(b1=bool(tg[0] > st[0] and tg[0] > fl[0]), b2=bool(tg[1] > fl[1]))


def write(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("liftoff")
    for k in ("--species", "--sample", "--base", "--merged", "--tg", "--templates", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--liftoff", default=LIFTOFF)
    p = sub.add_parser("chains")
    p.add_argument("--ref", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--gtf", action="append", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("tools")
    for k in ("--species", "--ref", "--exclude", "--strata", "--work", "--tg-arm", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--arm", action="append", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "liftoff":
        loci = LF.read_loci(Path(a.liftoff) / a.species / "liftoff_loci.tsv")
        sup = LF.read_support(Path(a.liftoff) / "support" / f"{a.sample}.tsv")
        cache = os.path.join(os.path.dirname(os.path.abspath(a.out)), "denovo_cache")
        base_reads = {t.gene: t.reads for t in T.read_templates(a.templates).values()}
        res = liftoff_scores(loci, sup, gtf_rows(a.base, cache), gtf_rows(a.merged, cache), gtf_rows(a.tg, cache),
                             base_reads)
        write(a.out, dict(sample=a.sample, **res))
    elif a.cmd == "chains":
        write(a.out, chains(a.ref, a.work, [g.split("=", 1) for g in a.gtf]))
    else:
        write(a.out, tools(a))
    print(open(a.out).read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run every suite; they pass**

Run: `RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 -m pytest bench/template_copies bench/liftoff_pairs bench/soto_subset bench/o1_hier bench/soto -q -p no:cacheprovider`
Expected: everything passes, with `bench/template_copies` contributing 24 tests. Report the total line.

- [ ] **Step 5: Commit**

```bash
git add bench/template_copies/score_tg.py bench/template_copies/test_score_tg.py
git commit -m "bench: template_copies scorer (Liftoff recall, P1, P2, gffcompare chains, register 1268's tool lead)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/template_copies/score_tg.py bench/template_copies/test_score_tg.py
```

---

### Task 5: templates for every run, and the pre-registration

**Files:**
- Create: `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md`
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_template_copies/<run>/templates.*` and `batches/*.fa`

**Interfaces:**
- Consumes: `run_tg.py templates`.
- Produces: the committed prereg, plus the template and batch counts per run that Tasks 6-9 use.

- [ ] **Step 1: Write the templates of the five runs (light, no projection)**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset
H=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
G=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
O=/mnt/linuxdisk/tmp/o1_template_copies
T="bash tools/rlock.sh light python3 bench/template_copies/run_tg.py templates --batch-size 6000"
RLOCK_TIMEOUT=600 $T --base /mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.gtf --genome $H --out $O/A119b_shipped
RLOCK_TIMEOUT=600 $T --base /mnt/linuxdisk/tmp/asm_gap/op/human.OP_primary.gtf --genome $H --out $O/A119b_OP
RLOCK_TIMEOUT=600 $T --base /mnt/linuxdisk/tmp/rustle_figures/runs/human_testis/human_testis.gtf --genome $H --out $O/testis
RLOCK_TIMEOUT=600 $T --base /mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/gorilla_KB3781.gtf --genome $G --out $O/KB3781
RLOCK_TIMEOUT=600 $T --base /mnt/linuxdisk/tmp/asm_gap/op/gorilla.OP_primary.gtf --genome $G --out $O/OR6737_OP
```

Expected: each prints `{"transcripts": ..., "loci": ..., "templates": ..., "batches": ..., "batch_size": 6000}`. For
A119b_shipped it is 94,685 templates in 16 batches; that count is known from planning. Keep the five JSON lines for
Step 2.

- [ ] **Step 2: Write the prereg**

Create `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md` with this content. Replace `<CODE>` with
`git rev-parse --short HEAD`, which is Task 4's commit, and fill the template table from Step 1's JSON lines.

````markdown
# Template-guided copy assembly (2026-10-14)

**Pre-registered before the first real projection.**
- Spec: `docs/superpowers/specs/2026-10-14-template-copies-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-template-copies.md`.
- Code: `bench/template_copies/` at commit `<CODE>`.
- Every number below comes from that code, run as written in the plan's Tasks 6-9.

## Question

Can copies that have their own reads but no full locus be assembled from those reads, guided by an expressed locus of
the same gene, without lowering precision and without losing the operating point's lead over StringTie and FLAIR? This
is the assembler using the family relation (user, 2026-10-14).

## Known before running (recorded now)

- **Where read-supported Liftoff copy pairs are lost** (identity ≥ 0.95, shipped families):

  | sample | pairs | shared | a copy has no assembled locus | assembled but not in a family | split |
  |---|---|---|---|---|---|
  | human A119b | 514 | 75 | 378 | 45 | 16 |
  | human testis (public, ERR13885926) | 38 | 12 | 24 | 2 | 0 |
  | gorilla KB3781 | 78 | 4 | 71 | 3 | 0 |

- **The operating point helps a little.** Of A119b's 799 read-supported copy loci, the shipped assembler fully covers
  318 and `OP_good` 350 (+43, −11).
- **The ceiling.** 309 extra copies have no full locus under `OP_good`. Their source placement is:
  - in a family: 61;
  - an assembled singleton: 46;
  - not assembled: 202.

  So a template exists for at most 107.
- **Copy shape.** Of A119b's copy loci without a full locus, 360 are multi-exon and 121 single-exon; 217 are partly
  covered and 264 have no locus at all.
- **Calibration (planning).** One 1,000-template batch of A119b's shipped templates was projected on the human index
  (89 s, 16.3 GB; `/mnt/linuxdisk/tmp/o1_template_copies/calib/`) to size the batches. Nothing from it was scored or
  used.
- **Refuted neighbours, not re-proposed:**
  - `--tied-seed` (registers 338, 435, 980);
  - family rescue levers (647, 648);
  - all secondaries seeding loci (1059);
  - SD-region semi-guided nodes (1007/1009);
  - annotated coordinates as de novo vertices (294);
  - the closed loop is expected to be small (1272).
- **Templates per run** (inputs, from `run_tg.py templates`):

  | run | base GTF | templates | batches |
  |---|---|---|---|
  | A119b_shipped | `runs/human_A119b/human_A119b.gtf` | <n> | <b> |
  | A119b_OP | `asm_gap/op/human.OP_primary.gtf` | <n> | <b> |
  | testis | `runs/human_testis/human_testis.gtf` | <n> | <b> |
  | KB3781 | `runs/gorilla_KB3781/gorilla_KB3781.gtf` | <n> | <b> |
  | OR6737_OP | `asm_gap/op/gorilla.OP_primary.gtf` | <n> | <b> |

## Method (fixed)

As the spec's section 1, with the plan's rulings R1-R3:
- identity ≥ 0.80 and ≥ 50% of the template aligned (`minimap2 -x splice -c -N 50 -p 0.1`);
- a candidate is covered < 50% by every base locus, and duplicates overlapping by ≥ 50% of the smaller union merge;
- ≥ 2 own primary reads (`-F 2308`) per candidate and per emitted locus;
- transcripts come only from the own reads' blocks and junctions;
- the run fails when more than half the candidates lie on contigs absent from the BAM.

## Metrics and bars (per run; species never pooled)

- **Recall (reported).** Read-supported Liftoff extra copies (identity ≥ 0.95) covered ≥ 50% by a locus, base vs
  merged; the template-reachable ones and how many were reached.
- **P1 (bar).** Two quantities must not drop: the merged set's gffcompare intron-chain precision against the sample's
  RefSeq must be ≥ the base's − 0.2 (percentage points), and its `located` share ≥ the base's − 0.002.
- **P2 (bar).** The added loci's `located` share must be ≥ that of the base loci whose representative has 2-5 reads,
  with at least one locus added. A 95% interval over the added loci is reported.
- **Tool lead (bar; A119b_OP and OR6737_OP).** Register 1268's scorer, genome-wide minus the development chromosome
  (human chr20, gorilla NC_073244.2). The arm `OP_primary_TG` (merged) must have intron-chain sensitivity > StringTie's
  and > FLAIR's (B1), and intron-chain precision > FLAIR's (B2).

## Readings (fixed now)

| development: A119b_shipped and A119b_OP P1 + P2, tool lead on A119b_OP and OR6737_OP | held out: testis and KB3781 P1 + P2 | reading |
|---|---|---|
| all pass | all pass | Port into the Rust assembler as an opt-in second pass; the recall gain is the reported effect. |
| all pass | any fails | Not adopted; the held-out failures are listed. |
| any fails | — | Not adopted. |

## Rulings (made while planning)

R1-R7 of the plan, copied verbatim.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
````

- [ ] **Step 3: Commit the prereg, before any real projection**

```bash
git add docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md
git commit -m "prereg: template-guided copy assembly (Liftoff recall; P1, P2 precision; tool lead; A119b dev, testis + KB3781 held out)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md
git log --oneline -1 -- docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md
```

---

### Task 6: development run 1, A119b on the shipped base

**Files:**
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_template_copies/A119b_shipped/`

**Interfaces:**
- Consumes: Tasks 3-5.
- Produces: `tg.gtf`, `merged.gtf`, `report.json`, `liftoff.json` and `chains.json`, all read by Task 10.

- [ ] **Step 1: Project every batch, one heavy call each**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset
W=/mnt/linuxdisk/tmp/o1_template_copies/A119b_shipped
X=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy /usr/bin/time -f "%e s %M KB" python3 bench/template_copies/run_tg.py project --out $W --index $X --batches 0 --threads 4
```

- Repeat with `--batches 1`, `2`, … up to the run's last batch, one call each.
- If a call takes under 4 minutes, run two batches per call (`--batches i,j`).
- A call that dies leaves no `.done` for its batch. Rerun the same batch.

- [ ] **Step 2: Assemble**

```bash
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy /usr/bin/time -f "%e s %M KB" python3 bench/template_copies/run_tg.py assemble \
    --base /mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.gtf \
    --bam /mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam --out $W --label A119b_shipped
```

Expected: a JSON line of counts with `contig_absent` near 0. A `SystemExit` about contig names is the R3 guard: stop
and report it.

- [ ] **Step 3: Score**

```bash
B=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.gtf
REF=/mnt/c/Users/jfris/Desktop/benchmark_collapse/ref/A119b.chm13v2.0_RefSeq.gtf.gz
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy python3 bench/template_copies/score_tg.py liftoff --species human \
    --sample human_A119b --base $B --merged $W/merged.gtf --tg $W/tg.gtf --templates $W/templates.tsv --out $W/liftoff.json
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy python3 bench/template_copies/score_tg.py chains --ref $REF --work $W/chains \
    --gtf base=$B --gtf merged=$W/merged.gtf --out $W/chains.json
```

Expected: both JSON files are written. If the `chains` call reaches the 10-minute limit, rerun it: gffcompare's
finished `.stats` files are reused.

---

### Task 7: development run 2, A119b on `OP_primary`, with the tool lead

**Files:**
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_template_copies/A119b_OP/`

**Interfaces:**
- Consumes: Tasks 3-5.
- Produces: as Task 6, plus `tools.json`, all read by Task 10.

- [ ] **Step 1: Project every batch**

As Task 6 Step 1, with `W=/mnt/linuxdisk/tmp/o1_template_copies/A119b_OP` and the same human index.

- [ ] **Step 2: Assemble and score**

```bash
W=/mnt/linuxdisk/tmp/o1_template_copies/A119b_OP
B=/mnt/linuxdisk/tmp/asm_gap/op/human.OP_primary.gtf
BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam
REF=/mnt/c/Users/jfris/Desktop/benchmark_collapse/ref/A119b.chm13v2.0_RefSeq.gtf.gz
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy python3 bench/template_copies/run_tg.py assemble --base $B --bam $BAM --out $W --label A119b_OP
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy python3 bench/template_copies/score_tg.py liftoff --species human \
    --sample human_A119b --base $B --merged $W/merged.gtf --tg $W/tg.gtf --templates $W/templates.tsv --out $W/liftoff.json
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy python3 bench/template_copies/score_tg.py chains --ref $REF --work $W/chains \
    --gtf base=$B --gtf merged=$W/merged.gtf --out $W/chains.json
```

- [ ] **Step 3: The tool lead (register 1268's scorer)**

First fill the gffcompare cache with one arm per call; then run one call with every arm (cached, so fast):

```bash
S=/mnt/c/Users/jfris/Desktop/benchmark_collapse
TOOLS="python3 bench/template_copies/score_tg.py tools --species human --ref $REF --exclude chr20 --strata /mnt/linuxdisk/tmp/asm_gap/op/human.strata --work $W/tools --tg-arm OP_primary_TG"
for arm in OP_primary=$B OP_primary_TG=$W/merged.gtf StringTie=$S/stringtie_A119b/A119b.stringtie.gtf FLAIR=$S/flair_A119b/A119b.flair.isoforms.gtf; do
  RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy $TOOLS --arm $arm --out $W/tools.$(echo $arm | cut -d= -f1).json
done
RLOCK_TIMEOUT=900 bash tools/rlock.sh heavy $TOOLS --arm OP_primary=$B --arm OP_primary_TG=$W/merged.gtf \
    --arm StringTie=$S/stringtie_A119b/A119b.stringtie.gtf --arm FLAIR=$S/flair_A119b/A119b.flair.isoforms.gtf --out $W/tools.json
```

Run the loop's four calls as four separate Bash calls if together they would exceed 10 minutes.

Expected: `tools.json` has `b1` and `b2` as booleans. The `OP_primary` row must reproduce register 1268's human
chain 23.6 / 15.6. If it does not, stop and report it: the scorer or its inputs changed.

---

### Task 8: held-out runs, human testis and gorilla KB3781

**Files:**
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_template_copies/testis/` and `/mnt/linuxdisk/tmp/o1_template_copies/KB3781/`

**Interfaces:**
- Consumes: Tasks 3-5.
- Produces: `report.json`, `liftoff.json` and `chains.json` per run, read by Task 10.

- [ ] **Step 1: testis (human index)**

```bash
W=/mnt/linuxdisk/tmp/o1_template_copies/testis
X=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi
B=/mnt/linuxdisk/tmp/rustle_figures/runs/human_testis/human_testis.gtf
BAM=/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam
REF=/mnt/c/Users/jfris/Desktop/benchmark_collapse/ref/A119b.chm13v2.0_RefSeq.gtf.gz
```

- Run the projection as Task 6 Step 1, one heavy call per batch.
- Then run `assemble` with `--label testis`.
- Then run `score_tg.py liftoff --species human --sample human_testis` and `score_tg.py chains`, exactly as Task 6
  Steps 2-3, with these variables.

- [ ] **Step 2: KB3781 (gorilla index)**

```bash
W=/mnt/linuxdisk/tmp/o1_template_copies/KB3781
X=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
B=/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/gorilla_KB3781.gtf
BAM=/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam
REF=/mnt/c/Users/jfris/Desktop/benchmark_collapse/ref/GGO.GCF_029281585.2_RefSeq.gtf.gz
```

- Run the projection.
- Then run `assemble` with `--label KB3781`.
- Then run `score_tg.py liftoff --species gorilla --sample gorilla_KB3781` and `chains`.
- The gorilla index (13.6 GB) peaks near 17 GB: one heavy call at a time.

---

### Task 9: the tool lead on gorilla OR6737 (`OP_primary`)

**Files:**
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_template_copies/OR6737_OP/`

**Interfaces:**
- Consumes: Tasks 3-5.
- Produces: `report.json`, `liftoff.json`, `chains.json` and `tools.json`, read by Task 10.

- [ ] **Step 1: Project, assemble and score**

```bash
W=/mnt/linuxdisk/tmp/o1_template_copies/OR6737_OP
X=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
B=/mnt/linuxdisk/tmp/asm_gap/op/gorilla.OP_primary.gtf
BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam
REF=/mnt/c/Users/jfris/Desktop/benchmark_collapse/ref/GGO.GCF_029281585.2_RefSeq.gtf.gz
```

- Run the projection.
- Then run `assemble` with `--label OR6737_OP`.
- Then run `score_tg.py liftoff --species gorilla --sample gorilla_OR6737` and `chains`, as Task 6.

- [ ] **Step 2: The tool lead**

Run as Task 7 Step 3, with these values:
- `--species gorilla --exclude NC_073244.2 --strata /mnt/linuxdisk/tmp/asm_gap/op/gorilla.strata`;
- the arms `OP_primary=$B`, `OP_primary_TG=$W/merged.gtf`,
  `StringTie=/mnt/c/Users/jfris/Desktop/benchmark_collapse/stringtie_GGO/GGO.stringtie.gtf` and
  `FLAIR=/mnt/c/Users/jfris/Desktop/benchmark_collapse/flair_GGO/GGO.flair.isoforms.gtf`.

Expected: the `OP_primary` row reproduces register 1268's gorilla chain 29.3 / 30.4. If it does not, stop and report
it.

---

### Task 10: results, readings, register rows, the thesis line

**Files:**
- Modify: `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md` (append under `## Results`)
- Modify: `docs/NEGATIVE_RESULTS_REGISTER.md` (append one row)
- Modify: `docs/THESIS_OBJECTIVES.md` (the O1 row, line 31: one sentence)
- Modify: `docs/superpowers/specs/2026-10-14-template-copies-design.md` (append `## Results`)

**Interfaces:**
- Consumes: the JSON files of Tasks 6-9.
- Produces: the final record.

- [ ] **Step 1: Append the results to the prereg**

Under `## Results`, add `### Results (date)`. Generate every number from the JSON files with a short script, never by
hand. Give one table per metric, with one row per run (A119b_shipped, A119b_OP, testis, KB3781, OR6737_OP):
- `report.json` counts: hits, candidates, merged candidates, contig absent, few reads, no locus, with loci, loci,
  transcripts;
- recall: copies, base, merged, reachable, reached;
- P1: chain precision base / merged, `located` base / merged, PASS or FAIL;
- P2: added loci, their `located` share [95% interval], the base low-support loci and their share, PASS or FAIL;
- the tool lead, A119b_OP and OR6737_OP only: chain sensitivity / precision of `OP_primary`, `OP_primary_TG`,
  StringTie and FLAIR, with B1 and B2.

Format shares to 4 decimals and gffcompare values to 1 decimal, as the scorer prints them.

- [ ] **Step 2: Apply the readings**

Add `### Verdict (date)` with the reading row that applies and one sentence per bar. "Development" means P1 and P2 on
A119b_shipped and A119b_OP plus the tool lead on A119b_OP and OR6737_OP. "Held out" means P1 and P2 on testis and
KB3781.

- [ ] **Step 3: Register row and thesis line**

- Find the last row with `grep -o '^| 1[0-9][0-9][0-9] ' docs/NEGATIVE_RESULTS_REGISTER.md | tail -1` and append row
  N+1 at the end of the table.
- It names the prereg, states the claim ("expressed loci as templates assemble read-supported copies without lowering
  precision or the lead over StringTie/FLAIR"), and gives the verdict with these numbers: recall base → merged per run,
  P2 shares, B1/B2.
- Write `|Δ|`-style pipes inside a cell as `\|`.
- In the O1 row of `docs/THESIS_OBJECTIVES.md` (line 31), append before the cell's closing ` |` one sentence:
  - "Template-guided copy assembly (2026-10-14, docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md, register N+1): <reading>,
    recall <base → merged> in A119b, P2 <share vs base>";
  - with the fail named if it failed.

- [ ] **Step 4: Point the spec at the results; run the suites; commit**

Append to the spec:

```markdown
## Results (date)

<the reading and its one sentence>. Full results: `docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md`; register row N+1.
```

```bash
RLOCK_TIMEOUT=600 bash tools/rlock.sh light python3 -m pytest bench/template_copies bench/liftoff_pairs bench/soto_subset bench/o1_hier bench/soto -q -p no:cacheprovider
git add docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md docs/NEGATIVE_RESULTS_REGISTER.md docs/THESIS_OBJECTIVES.md docs/superpowers/specs/2026-10-14-template-copies-design.md
git commit -m "Template-guided copies: results, verdict, register row, O1 thesis line, spec pointer

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_template_copies_2026-10-14.md docs/NEGATIVE_RESULTS_REGISTER.md docs/THESIS_OBJECTIVES.md docs/superpowers/specs/2026-10-14-template-copies-design.md
```

Expected: every test passes.
