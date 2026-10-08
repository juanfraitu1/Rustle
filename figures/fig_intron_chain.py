"""fig1 — intron-chain sensitivity and precision of Rustle's transcripts vs StringTie, FLAIR and IsoSeq collapse on
the SAME BAMs, and (supplementary) of Rustle's two configurations on every sample of the registry.

Scope: GENOME-WIDE for every sample (assembly.evaluation_contigs: every contig the sample's annotation covers; the
human CHM13 RefSeq annotation leaves out chrM, so chrM is left out for every method alike). The full comparison runs
on every registry sample with all three lab baselines (assembly.benchmark_samples: human_A119b, gorilla_OR6737); the
supplementary table covers every sample (assembly.all_samples), Rustle only, each against its own annotation.

Tables (figures/data/):

  fig1_gffcompare  species, scope, tool, level, sn, pr, matching, n_query, n_ref, n_query_multiexon, n_ref_multiexon,
                   sample
                   gffcompare of each method against the annotation, every method restricted to the same contigs.
                   `species` is the sample key the panels use (gorilla = gorilla_OR6737, human = human_A119b; `sample`
                   is the registry id). sn / pr are FRACTIONS (gffcompare prints percent); `matching` is gffcompare's
                   "Matching intron chains / transcripts / loci" for those three levels and blank for base / exon /
                   intron; n_query / n_ref = query / reference mRNAs, and the *_multiexon columns their multi-exon
                   subsets (the intron-chain denominators: sn = matching / n_ref_multiexon, pr = matching query
                   transcripts / n_query_multiexon).

  fig1_support     species, scope, tool, count_mode, min_reads, n_ref, n_matched, sn, sample
                   intron-chain sensitivity over the DISTINCT multi-exon reference intron chains (unit = contig +
                   intron chain; reference transcripts that differ only in their ends share one chain) supported by
                   >= min_reads primary alignments (samtools -F 2308) whose intron chain EQUALS the reference chain
                   (min_reads 0 = the whole annotation). A chain is matched by a method when gffcompare gives class
                   code '=' to any reference transcript carrying it. count_mode 'reads' counts alignments;
                   'distinct_ends' counts distinct (start, end) alignment spans carrying the chain (the assembler's
                   coordinate-duplicate key, `--keep-coordinate-duplicates` off). The read's strand is not used.

  fig1_paired      species, scope, count_mode, min_reads, n_ref, arm, other, n_both, n_arm_only, n_other_only,
                   n_neither, sn_arm, sn_other, diff, diff_lo, diff_hi, mcnemar_p, sample
                   PAIRED comparison of Rustle (`arm`) with each other method on the SAME reference chains of one
                   stratum (min_reads 1, 2, 5): the 2x2 table of matched-by-both / Rustle-only / other-only / neither,
                   diff = sn_arm - sn_other (a fraction), [diff_lo, diff_hi] = Tango's (1998) asymptotic score 95%
                   interval for a paired difference of proportions, mcnemar_p = exact two-sided McNemar p (binomial
                   test on the discordant chains; written as text so values below 1e-308 keep their exponent).
                   Unadjusted for the four comparisons per stratum; chains of one gene are not independent.

  Mode: every method here is ANNOTATION-FREE (de novo): Rustle assemble (reads + genome), StringTie -L without -G,
  FLAIR collapse without annotation (flair correct skipped), IsoSeq collapse (assembly.MODE_DENOVO_METHODS). The
  figure and every table note say so.

  fig1_guided      (supplementary; ONLY when a sample has an annotation-guided StringTie/FLAIR GTF registered,
                   samples.tsv stringtie_guided_gtf / flair_guided_gtf) the fig1_samples columns plus `mode`
                   (annotation-guided) and `guided_gtf`, for the guided tools only, drawn as fig1g_guided. Never in a
                   panel or table with the annotation-free methods; no Rustle row (Rustle has no annotation-guided
                   transcript assembly; docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md). Until then the
                   figure prints "guided comparison: not available (guided StringTie/FLAIR GTFs not supplied)".

  fig1_samples     (supplementary; every sample, Rustle's two configurations) sample, label, species, tissue, genome,
                   annotation, scope, tool, n_query, n_query_multiexon, n_ref_multiexon, matching, sn, pr, n_ref_ge2,
                   n_matched_ge2, sn_ge2
                   sn / pr / matching: gffcompare intron-chain level against the sample's own annotation (chimpanzee
                   and orangutan: the RefSeq GFF3 converted by gff_to_gtf, the converter of the human and gorilla GTFs);
                   *_ge2: the same sensitivity on the reference chains carried exactly by >= 2 primary alignments
                   (count_mode 'reads', the fig1_support rule). Drawn as figure fig1s_samples.

Heavy steps build() runs (foreground, serial, cached under ${work}; `figs_budget_s` / `fig1_budget_s` bound one call
and `figs_plan=1` lists the units without running them, see assembly.py):
  * assembly.gffcompare() per sample x method (restricted GTFs and gffcompare outputs live in assembly.eval_dir;
    Rustle's genome-wide assemblies come from the run cache and are made there only if missing)
  * one streaming pass over each sample's BAM (pysam, primary alignments of the evaluation contigs) counting the
    reads that carry each reference intron chain, ONE CONTIG AT A TIME (${work}/fig1/<sample>.<scope>.chain_support
    .tsv.parts/<contig>.tsv, each reused while newer than the BAM and the restricted annotation and written from the
    same BAM and annotation), assembled into ${work}/fig1/<sample>.<scope>.chain_support.tsv (the same rows and order
    as one pass over all contigs)
Development only: cfg['fig1_support_contigs_<species>'] (comma list) narrows the exact-chain stratum to those contigs;
the support and paired tables then carry a "provisional" note and the panels say "subset".
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import figlib

try:  # plot() must work without the heavy-side dependencies
    import assembly
    import samples
except Exception:  # noqa: BLE001
    assembly = samples = None

META = {
    "id": "fig1",
    "title": "Intron-chain sensitivity and precision against the annotation: Rustle, StringTie, FLAIR and IsoSeq "
             "collapse on the same alignments, all annotation-free (de novo)",
    "claim": "Annotation-free (de novo) comparison: Rustle assemble (reads + genome), StringTie -L without -G, FLAIR "
             "collapse without annotation (flair correct skipped), IsoSeq collapse; guided comparison: not available "
             "(guided StringTie/FLAIR GTFs not supplied). "
             "Against the RefSeq annotation, with every method's transcripts built from the same alignments and cut "
             "to the same contigs (gorilla OR6737 genome-wide; human A119b chr20-22 in the current tables, "
             "genome-wide after the rebuild), Rustle's intron-chain precision is slightly above StringTie's in "
             "gorilla (35.3% vs 34.6%) and above it in human (17.3% vs 14.8%), and 2.0-7.2 times that of FLAIR and "
             "IsoSeq collapse; its sensitivity is above StringTie's in both species, above FLAIR's in gorilla (27.1% "
             "vs 25.5%) and equal to it in human (23.0% vs 23.1%), and IsoSeq collapse has the highest sensitivity "
             "over the whole annotation (30.8%, 27.7%). On the reference intron chains that >= 2 primary reads carry "
             "exactly (Rustle's minimum read support: independent of every method's output, but matched to Rustle's "
             "design; gorilla 27,447 chains, human 2,677) Rustle's two configurations have the highest sensitivity "
             "point estimates (93.4%, 84.4%; equal in human). Paired on those chains (Tango 95% interval, exact "
             "McNemar test) Rustle's sensitivity is higher than StringTie's, FLAIR's and IsoSeq collapse's in gorilla "
             "and than StringTie's and FLAIR's in human, and not separable from IsoSeq collapse's in human (+1.5 "
             "points, 95% CI -0.4 to +3.5, p = 0.13). On chains carried by >= 1 read IsoSeq collapse has the higher "
             "sensitivity in both species (-10.1 and -13.6 points). Adding the secondary alignments that score within "
             "2% of the read's best alignment score adds 162 / 4 matched chains at a small precision cost (35.6% -> "
             "35.3%, 18.7% -> 17.3%). Supplementary (fig1s_samples, after the rebuild): Rustle's two configurations "
             "on all six samples, each against its own annotation.",
    "tables": ["fig1_gffcompare", "fig1_support", "fig1_paired"],
    # drawn as fig1s_samples / fig1g_guided when present (built by the same `make.py data fig1`; fig1_guided only
    # when an annotation-guided GTF is registered)
    "supplementary_tables": ["fig1_samples", "fig1_guided"],
}

LEVELS = ["base", "exon", "intron", "intron_chain", "transcript", "locus"]
MATCH_KEY = {"intron_chain": "matching_intron_chains", "transcript": "matching_transcripts",
             "locus": "matching_loci"}
SUPPORT_MIN = (0, 1, 2, 5)
PAIRED_MIN = (1, 2, 5)
PAIRED_SHOWN = 2          # the stratum drawn in panels g/h (Rustle's minimum read support)
COUNT_MODES = ("reads", "distinct_ends")
SPECIES = ["gorilla", "human"]
TOOLS = list(figlib.TOOL_ORDER)
BAM_EXCLUDE = 2308  # unmapped | secondary | supplementary  (the `-F 2308` invariant)
Z95 = 1.959963984540054


# ------------------------------------------------------------------------------------------------ helpers (build)
def _natural(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def scope_label(cfg: dict, species: str) -> str:
    return assembly.scope_label(cfg, species)


def support_scope(cfg: dict, species: str):
    """(contigs or None, label, is_subset) of the exact-chain stratum: the evaluation scope, unless the development
    key cfg['fig1_support_contigs_<species>'] narrows it (the tables are then marked provisional)."""
    raw = cfg.get(f"fig1_support_contigs_{species}", "")
    sub = sorted({c.strip() for c in raw.split(",") if c.strip()}, key=_natural)
    if not sub:
        return None, scope_label(cfg, species), False
    ev = assembly.evaluation_contigs(cfg, species)
    if ev is not None:
        sub = [c for c in sub if c in ev]
        if set(sub) == ev:
            return None, scope_label(cfg, species), False
    return set(sub), ",".join(sub) + " only", True


def gffcompare_rows(species: str, scope: str, tool: str, stats: dict) -> list[list]:
    """Tidy rows of one parsed gffcompare .stats (assembly.parse_stats)."""
    rows = []
    for lv in LEVELS:
        if f"{lv}_sn" not in stats:
            continue
        rows.append([species, scope, tool, lv, stats[f"{lv}_sn"] / 100.0, stats[f"{lv}_pr"] / 100.0,
                     stats.get(MATCH_KEY.get(lv, ""), None), stats.get("query_mrnas"), stats.get("ref_mrnas"),
                     stats.get("query_multiexon"), stats.get("ref_multiexon")])
    return rows


def ref_chains(ref_gtf, contigs: set[str] | None = None) -> dict:
    """{(chrom, intron chain): [reference transcript ids]} for multi-exon reference transcripts (0-based half-open
    introns, as assembly.ref_transcripts)."""
    out: dict = {}
    for tid, t in assembly.ref_transcripts(ref_gtf).items():
        if not t["introns"] or (contigs is not None and t["chrom"] not in contigs):
            continue
        out.setdefault((t["chrom"], t["introns"]), []).append(tid)
    return out


def read_intron_chain(cigartuples, ref_start: int):
    """(intron chain, alignment end) of one alignment: every N operation is an intron [start, end) 0-based."""
    pos = ref_start
    introns = []
    for op, n in cigartuples:
        if op == 3:  # N
            introns.append((pos, pos + n))
            pos += n
        elif op in (0, 2, 7, 8):  # M D = X consume the reference
            pos += n
    return tuple(introns), pos


def count_exact_chains(bam_path, wanted: dict[str, set], threads: int = 2) -> dict:
    """Stream the primary alignments of the contigs in `wanted` ({chrom: {chain, ...}}) and return
    {(chrom, chain): (n_alignments, n_distinct_(start,end))} for every wanted chain carried by >= 1 alignment.
    Spans are kept as one int (start << 32 | end) per alignment: a genome-wide pass holds millions of them."""
    import pysam

    n_reads: dict = {}
    spans: dict = {}
    with pysam.AlignmentFile(str(bam_path), "rb", threads=max(1, int(threads))) as bam:
        present = set(bam.references)
        for chrom in sorted(wanted, key=_natural):
            want = wanted[chrom]
            if not want or chrom not in present:
                continue
            for r in bam.fetch(chrom):
                if r.flag & BAM_EXCLUDE:
                    continue
                ct = r.cigartuples
                if not ct:
                    continue
                chain, end = read_intron_chain(ct, r.reference_start)
                if not chain or chain not in want:
                    continue
                key = (chrom, chain)
                n_reads[key] = n_reads.get(key, 0) + 1
                spans.setdefault(key, set()).add((r.reference_start << 32) | end)
    return {k: (n, len(spans[k])) for k, n in n_reads.items()}


def _chain_str(chain) -> str:
    return ",".join(f"{a}-{b}" for a, b in chain)


def _parse_chain(s: str):
    return tuple(tuple(int(x) for x in iv.split("-")) for iv in s.split(",")) if s else ()


def write_chain_support(path: Path, support: dict, *, bam, ref_gtf, contigs: str):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(f"# bam: {bam}\n# ref: {ref_gtf}\n# contigs: {contigs}\n# primary alignments (-F {BAM_EXCLUDE}); "
                 "introns 0-based half-open; only chains of multi-exon reference transcripts are counted\n")
        fh.write("chrom\tintrons\tn_reads\tn_distinct_ends\n")
        for (chrom, chain), (n, d) in sorted(support.items(), key=lambda kv: (_natural(kv[0][0]), kv[0][1])):
            fh.write(f"{chrom}\t{_chain_str(chain)}\t{n}\t{d}\n")
    tmp.replace(path)
    return path


def read_chain_support(path) -> dict:
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom\t"):
                continue
            c, ch, n, d = line.rstrip("\n").split("\t")
            out[(c, _parse_chain(ch))] = (int(n), int(d))
    return out


def chain_support_path(cfg: dict, species: str, label: str) -> Path:
    key = label.replace(" only", "").replace(",", "_")
    if len(key) > 64:
        key = "set_" + __import__("hashlib").md5(key.encode()).hexdigest()[:12]
    return figlib.work_dir(cfg, "fig1") / f"{species}.{key}.chain_support.tsv"  # cache keyed by the contig set


def _part_ok(part: Path, bam, ref_gtf) -> bool:
    """A per-contig part is reused when it is newer than the BAM and the restricted annotation and was counted from
    the same BAM and annotation paths."""
    if not figlib.fresh(part, bam, ref_gtf):
        return False
    head = _cached_header(part)
    return head.get("bam") == str(bam) and head.get("ref") == str(ref_gtf)


def chain_support_todo(cfg: dict, species: str, chains: dict, ref_gtf, label: str) -> list[str]:
    """Contigs whose exact-chain support still has to be counted ([] when the assembled file is current)."""
    bam = samples.get(cfg, species)["bam"]
    path = chain_support_path(cfg, species, label)
    head = _cached_header(path)
    if figlib.fresh(path, bam, ref_gtf) and head.get("bam") == str(bam) and head.get("ref") == str(ref_gtf):
        return []
    parts = Path(str(path) + ".parts")
    return [c for c in sorted({chrom for chrom, _ in chains}, key=_natural)
            if not _part_ok(parts / f"{c}.tsv", bam, ref_gtf)]


def chain_support(cfg: dict, species: str, chains: dict, ref_gtf, label: str, *, force=False, budget=None) -> Path:
    """Cached exact-chain read support of the reference chains (HEAVY: one pass over the evaluation contigs, one
    contig at a time; `budget` (assembly.Budget) is checked before each contig that still has to be counted)."""
    bam = samples.get(cfg, species)["bam"]
    path = chain_support_path(cfg, species, label)
    parts = Path(str(path) + ".parts")
    if force:
        import shutil
        shutil.rmtree(parts, ignore_errors=True)
    elif not chain_support_todo(cfg, species, chains, ref_gtf, label) and path.exists():
        return path
    wanted: dict = {}
    for chrom, chain in chains:
        wanted.setdefault(chrom, set()).add(chain)
    parts.mkdir(parents=True, exist_ok=True)
    for chrom in sorted(wanted, key=_natural):
        part = parts / f"{chrom}.tsv"
        if not force and _part_ok(part, bam, ref_gtf):
            continue
        if budget is not None:
            budget.check(f"exact-chain support of {species} {chrom}")
        support = count_exact_chains(bam, {chrom: wanted[chrom]}, threads=int(cfg.get("threads", "4")))
        write_chain_support(part, support, bam=bam, ref_gtf=ref_gtf, contigs=chrom)
    support: dict = {}
    for chrom in wanted:
        support.update(read_chain_support(parts / f"{chrom}.tsv"))
    write_chain_support(path, support, bam=bam, ref_gtf=ref_gtf, contigs=label.replace(" only", ""))
    return path


def _cached_header(path) -> dict:
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("# "):
                    break
                k, _, v = line[2:].rstrip("\n").partition(": ")
                out[k] = v
    except OSError:
        pass
    return out


def matched_chains(chains: dict, matched_refs: set) -> set:
    """Reference chains (keys of `chains`) with at least one reference transcript matched '=' by the arm."""
    return {key for key, tids in chains.items() if any(t in matched_refs for t in tids)}


def support_rows(species: str, scope: str, tool: str, chains: dict, support: dict, matched: set) -> list[list]:
    """Sensitivity over reference chains with >= k exact-chain reads, for k in SUPPORT_MIN and both count modes."""
    rows = []
    for mi, mode in enumerate(COUNT_MODES):
        for k in SUPPORT_MIN:
            keys = [key for key in chains if support.get(key, (0, 0))[mi] >= k]
            n = len(keys)
            m = sum(1 for key in keys if key in matched)
            rows.append([species, scope, tool, mode, k, n, m, (m / n) if n else None])
    return rows


# ---- paired statistics (pure python: build() and the tests need no scipy)
def mcnemar_log10p(b: int, c: int) -> float:
    """log10 of the exact two-sided McNemar p: 2 * P(X <= min(b, c)), X ~ Binomial(b + c, 1/2), capped at 1."""
    n = b + c
    if n == 0:
        return 0.0
    k = min(b, c)
    lf = math.lgamma(n + 1)
    terms = [lf - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in range(k + 1)]
    mx = max(terms)
    log_tail = mx + math.log(sum(math.exp(t - mx) for t in terms)) - n * math.log(2.0)
    return min(0.0, (log_tail + math.log(2.0)) / math.log(10.0))


def format_p(log10p: float) -> str:
    """A p value as text; below the float range the exponent is kept (e.g. '3.1e-412')."""
    if log10p >= -300:
        p = 10 ** log10p
        return "1" if p >= 0.9995 else f"{p:.3g}"
    e = math.floor(log10p)
    return f"{10 ** (log10p - e):.2f}e{e}"


def parse_p(s: str) -> float:
    """log10 of a p written by format_p."""
    m, _, e = s.partition("e")
    return math.log10(float(m)) + (int(e) if e else 0)


def tango_ci(b: int, c: int, n: int, z: float = Z95):
    """Tango (1998, Stat Med 17:891) asymptotic score interval for a paired difference of proportions
    p1 - p2 = (b - c) / n, with b = pairs positive only under 1, c = only under 2, n = all pairs.
    Score statistic T(d) = (b - c - n d) / sqrt(n (2 q + d (1 - d))), q = the restricted MLE of the 'only 2' cell
    probability; the interval is {d : |T(d)| <= z}. T(0) is McNemar's (b - c) / sqrt(b + c)."""
    if n <= 0:
        return None, None

    def T(d):
        A = 2.0 * n
        B = -b - c + (2.0 * n - b + c) * d
        C = -c * d * (1.0 - d)
        q = (math.sqrt(max(B * B - 4.0 * A * C, 0.0)) - B) / (2.0 * A)
        var = n * (2.0 * q + d * (1.0 - d))
        num = b - c - n * d
        if var <= 1e-300:
            return 0.0 if num == 0 else math.copysign(math.inf, num)
        return num / math.sqrt(var)

    d_hat = (b - c) / n
    eps = 1e-12

    def root(lo, hi, target):  # T is decreasing in d on [-1, 1]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if T(mid) > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    lower = -1.0 if T(-1.0 + eps) <= z else root(-1.0 + eps, d_hat, z)
    upper = 1.0 if T(1.0 - eps) >= -z else root(d_hat, 1.0 - eps, -z)
    return lower, upper


def paired_rows(species: str, scope: str, chains: dict, support: dict, matched: dict) -> list[list]:
    """Rustle vs every other arm on the same chains of each stratum (PAIRED_MIN, both count modes)."""
    rows = []
    if "rustle" not in matched:
        return rows
    A = matched["rustle"]
    for mi, mode in enumerate(COUNT_MODES):
        for k in PAIRED_MIN:
            keys = [key for key in chains if support.get(key, (0, 0))[mi] >= k]
            n = len(keys)
            for other in TOOLS:
                if other == "rustle" or other not in matched:
                    continue
                B = matched[other]
                both = a_only = b_only = 0
                for key in keys:
                    ia, ib = key in A, key in B
                    both += ia and ib
                    a_only += ia and not ib
                    b_only += ib and not ia
                neither = n - both - a_only - b_only
                lo, hi = tango_ci(a_only, b_only, n)
                rows.append([species, scope, mode, k, n, "rustle", other, both, a_only, b_only, neither,
                             (both + a_only) / n if n else None, (both + b_only) / n if n else None,
                             (a_only - b_only) / n if n else None, lo, hi,
                             format_p(mcnemar_log10p(a_only, b_only))])
    return rows


GC_HEADER = ["species", "scope", "tool", "level", "sn", "pr", "matching", "n_query", "n_ref", "n_query_multiexon",
             "n_ref_multiexon", "sample"]
SUPPORT_HEADER = ["species", "scope", "tool", "count_mode", "min_reads", "n_ref", "n_matched", "sn", "sample"]
PAIRED_HEADER = ["species", "scope", "count_mode", "min_reads", "n_ref", "arm", "other", "n_both", "n_arm_only",
                 "n_other_only", "n_neither", "sn_arm", "sn_other", "diff", "diff_lo", "diff_hi", "mcnemar_p", "sample"]
SAMPLES_HEADER = ["sample", "label", "species", "tissue", "genome", "annotation", "scope", "tool", "n_query",
                  "n_query_multiexon", "n_ref_multiexon", "matching", "sn", "pr", "n_ref_ge2", "n_matched_ge2", "sn_ge2"]
RUSTLE_TOOLS = ["rustle", "rustle_primary"]
SAMPLES_MIN_READS = 2   # the supplementary stratum: Rustle's minimum read support (as panels g/h)


def annotation_name(cfg: dict, key: str) -> str:
    row = samples.get(cfg, key)
    if row["annotation_gtf"]:
        return Path(row["annotation_gtf"]).name
    return f"{Path(row['annotation_gff']).name} (GFF3 -> GTF with gff_to_gtf)"


def scope_notes(cfg: dict, keys, tools_of=None) -> list[str]:
    """What 'genome-wide' means for these samples (the contigs left out because nothing annotates them, with the
    transcripts each method has there)."""
    out = []
    for key in keys:
        un = assembly.unannotated_contigs(cfg, key)
        if assembly.is_genome_wide(cfg, key):
            left = ""
            if un:
                tools = (tools_of or {}).get(key, RUSTLE_TOOLS)
                n = {t: sum(assembly.outside_scope(cfg, key, t).values()) for t in tools}
                left = (f"; left out for every method because the annotation has no record there: {', '.join(un)} "
                        f"(transcripts there: " + ", ".join(f"{t} {v:,}" for t, v in n.items()) + ")")
            out.append(f"{key} ({samples.resolve(cfg, key)}): genome-wide = every contig its annotation covers"
                       + (left or " (all contigs)"))
        else:
            out.append(f"provisional: {key} restricted to {assembly.scope_label(cfg, key)} (development key "
                       f"eval_contigs_{samples.resolve(cfg, key)}); the figure build scores every annotated contig")
    return out


class _Runs:
    """gffcompare runs, reference chains and exact-chain support per sample, computed once per build."""

    def __init__(self, cfg, force, budget):
        self.cfg, self.force, self.budget = cfg, force, budget
        self._gc, self._ref, self._sup = {}, {}, {}

    def gc(self, key, tool):
        if (key, tool) not in self._gc:
            if self.force or assembly.gffcompare_state(self.cfg, key, tool):
                self.budget.check(f"gffcompare of {key} {tool}")
            self._gc[(key, tool)] = assembly.gffcompare(self.cfg, key, tool, force=self.force)
        return self._gc[(key, tool)]

    def ref(self, key):
        if key not in self._ref:
            ref_gtf = assembly.ensure_ref(self.cfg, key, force=self.force)
            sub, label, is_subset = support_scope(self.cfg, key)
            self._ref[key] = (ref_gtf, ref_chains(ref_gtf, sub), label, is_subset)
        return self._ref[key]

    def support(self, key):
        if key not in self._sup:
            ref_gtf, chains, label, _ = self.ref(key)
            path = chain_support(self.cfg, key, chains, ref_gtf, label, force=self.force, budget=self.budget)
            self._sup[key] = (path, read_chain_support(path))
        return self._sup[key]


def build(cfg: dict, data_dir: Path, force: bool = False):
    if assembly.plan_only(cfg):
        return plan(cfg)
    budget = assembly.Budget(cfg, "fig1")
    R = _Runs(cfg, force, budget)
    main = assembly.benchmark_samples(cfg)
    everyone = assembly.all_samples(cfg)
    gc_rows, sup_rows, pair_rows, samp_rows, inputs, sub_notes = [], [], [], [], {}, []
    for key in main:   # the full comparison: every sample with the three lab baselines
        sid = samples.resolve(cfg, key)
        scope = scope_label(cfg, key)
        runs = {tool: R.gc(key, tool) for tool in TOOLS}
        for tool in TOOLS:
            gc_rows += [r + [sid] for r in gffcompare_rows(key, scope, tool, assembly.parse_stats(runs[tool]["stats"]))]
            inputs[f"{key} {tool} stats"] = runs[tool]["stats"]
        ref_gtf, chains, label, is_subset = R.ref(key)
        sup_path, support = R.support(key)   # heavy: one BAM pass, one contig at a time
        matched = {tool: matched_chains(chains, assembly.exact_matched_refs(runs[tool]["tmap"])) for tool in TOOLS}
        for tool in TOOLS:
            sup_rows += [r + [sid] for r in support_rows(key, label, tool, chains, support, matched[tool])]
            inputs[f"{key} {tool} tmap"] = runs[tool]["tmap"]
        pair_rows += [r + [sid] for r in paired_rows(key, label, chains, support, matched)]
        if is_subset:
            sub_notes.append(f"provisional: {key} exact-chain stratum restricted to {label} (development subset "
                             f"cfg fig1_support_contigs_{key}); make.py data fig1 without that key counts {scope}")
        row = samples.get(cfg, key)
        inputs[f"{key} bam"] = row["bam"]
        inputs[f"{key} annotation"] = assembly.annotation_gtf(cfg, key)
        inputs[f"{key} restricted annotation"] = ref_gtf
        inputs[f"{key} chain support"] = sup_path
    samp_inputs = {}
    for key in everyone:   # supplementary: Rustle's two configurations on every sample, against its own annotation
        row = samples.get(cfg, key)
        ref_gtf, chains, label, _ = R.ref(key)
        sup_path, support = R.support(key)
        ge2 = [k for k in chains if support.get(k, (0, 0))[0] >= SAMPLES_MIN_READS]
        for tool in RUSTLE_TOOLS:
            run = R.gc(key, tool)
            st = assembly.parse_stats(run["stats"])
            matched = matched_chains(chains, assembly.exact_matched_refs(run["tmap"]))
            m2 = sum(1 for k in ge2 if k in matched)
            samp_rows.append([row["id"], assembly.sample_label(cfg, key), row["species"], row["tissue"], row["genome"],
                              annotation_name(cfg, key), label, tool, st.get("query_mrnas"), st.get("query_multiexon"),
                              st.get("ref_multiexon"), st.get("matching_intron_chains"),
                              st["intron_chain_sn"] / 100.0, st["intron_chain_pr"] / 100.0, len(ge2), m2,
                              (m2 / len(ge2)) if ge2 else None])
            samp_inputs[f"{row['id']} {tool} stats"] = run["stats"]
            samp_inputs[f"{row['id']} {tool} tmap"] = run["tmap"]
        samp_inputs[f"{row['id']} bam"] = row["bam"]
        samp_inputs[f"{row['id']} annotation"] = assembly.annotation_gtf(cfg, key)
        samp_inputs[f"{row['id']} chain support"] = sup_path
    notes = ["mode: " + assembly.MODE_DENOVO_METHODS + "; " + assembly.guided_status(cfg),
             "methods: rustle = pipeline driver default (loci built from primary alignments plus the secondary "
             "alignments whose alignment score is at least 98% of the read's best score anywhere in the genome); "
             "rustle_primary = the same run with --no-seed-secondaries; stringtie/flair/isoseq = the lab's GTFs "
             "(StringTie 3.0.1, FLAIR 3.0.1, IsoSeq collapse 26.2.0): "
             + ", ".join(f"{k} {t} {samples.baseline(cfg, k, t)}" for k in main for t in samples.BASELINES),
             "gffcompare 0.12.10 -r <annotation restricted to the scope> <method restricted to the same contigs>"] \
        + scope_notes(cfg, main, {k: TOOLS for k in main})
    unit = ("unit = distinct multi-exon reference intron chain (contig + chain); support = primary alignments "
            "(-F 2308) whose intron chain equals it; min_reads 0 = whole annotation")
    gen = "figures/fig_intron_chain.py build"
    figlib.write_table("fig1_gffcompare", GC_HEADER, gc_rows, generator=gen, inputs=inputs, notes=notes,
                       data_dir=data_dir)
    figlib.write_table("fig1_support", SUPPORT_HEADER, sup_rows, generator=gen, inputs=inputs,
                       notes=sub_notes + notes + [unit], data_dir=data_dir)
    figlib.write_table("fig1_paired", PAIRED_HEADER, pair_rows, generator=gen, inputs=inputs,
                       notes=sub_notes + notes + [unit, "paired on the same chains: diff = sn_arm - sn_other; "
                                                  "[diff_lo, diff_hi] = Tango asymptotic score 95% interval; "
                                                  "mcnemar_p = exact two-sided McNemar (binomial on the discordant "
                                                  "chains); unadjusted for multiplicity; chains of one gene are not "
                                                  "independent"],
                       data_dir=data_dir)
    figlib.write_table("fig1_samples", SAMPLES_HEADER, samp_rows, generator=gen, inputs=samp_inputs,
                       notes=["supplementary: Rustle's two configurations on every registry sample, each against its "
                              "own annotation; numbers from different samples or species are never pooled",
                              "sn / pr / matching = gffcompare 0.12.10 intron-chain level (fractions; denominators: "
                              "multi-exon reference / query transcripts); *_ge2 = sensitivity on the distinct "
                              f"multi-exon reference chains carried exactly by >= {SAMPLES_MIN_READS} primary alignments "
                              "(-F 2308, count mode 'reads'; the fig1_support rule)",
                              "annotations: human samples CHM13 v2.0 RefSeq; gorilla samples GCF_029281585.2 RefSeq; "
                              "chimpanzee GCF_028858775.2 and orangutan GCF_028885625.2 RefSeq GFF3 converted with "
                              "gff_to_gtf (the converter of the human and gorilla GTFs; the gorilla GFF3 converted "
                              "this way is byte-identical to the gorilla GTF, md5 17052eec)"]
                       + scope_notes(cfg, everyone) + [unit], data_dir=data_dir)
    build_guided(cfg, data_dir, R, unit)


GUIDED_HEADER = ["sample", "label", "species", "tissue", "genome", "annotation", "scope", "mode", "tool", "n_query",
                 "n_query_multiexon", "n_ref_multiexon", "matching", "sn", "pr", "n_ref_ge2", "n_matched_ge2", "sn_ge2",
                 "guided_gtf"]


def build_guided(cfg: dict, data_dir: Path, R, unit: str):
    """fig1_guided: the annotation-guided StringTie/FLAIR runs (only when registered), the fig1_samples measures,
    against each sample's own annotation. Separate table; never mixed with the annotation-free methods."""
    keys = assembly.guided_samples(cfg)
    if not keys:
        print(f"[fig1] {assembly.GUIDED_NA}", file=sys.stderr)
        return
    rows, inputs = [], {}
    for key in keys:
        row = samples.get(cfg, key)
        ref_gtf, chains, label, _ = R.ref(key)
        sup_path, support = R.support(key)
        ge2 = [k for k in chains if support.get(k, (0, 0))[0] >= SAMPLES_MIN_READS]
        for tool in assembly.guided_tools(cfg, key):
            run = R.gc(key, tool)
            st = assembly.parse_stats(run["stats"])
            matched = matched_chains(chains, assembly.exact_matched_refs(run["tmap"]))
            m2 = sum(1 for k in ge2 if k in matched)
            src = samples.baseline(cfg, key, tool)
            rows.append([row["id"], assembly.sample_label(cfg, key), row["species"], row["tissue"], row["genome"],
                         annotation_name(cfg, key), label, assembly.MODE_GUIDED, tool, st.get("query_mrnas"),
                         st.get("query_multiexon"), st.get("ref_multiexon"), st.get("matching_intron_chains"),
                         st["intron_chain_sn"] / 100.0, st["intron_chain_pr"] / 100.0, len(ge2), m2,
                         (m2 / len(ge2)) if ge2 else None, src])
            inputs[f"{row['id']} {tool} gtf"] = src
            inputs[f"{row['id']} {tool} stats"] = run["stats"]
        inputs[f"{row['id']} chain support"] = sup_path
    figlib.write_table("fig1_guided", GUIDED_HEADER, rows, generator="figures/fig_intron_chain.py build",
                       inputs=inputs,
                       notes=["mode: " + assembly.MODE_GUIDED + " (StringTie -G / FLAIR with the annotation, as "
                              "supplied in samples.tsv); " + assembly.GUIDED_CAVEAT,
                              assembly.RUSTLE_NO_GUIDED,
                              "pre-registered: docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md",
                              "sn / pr / matching = gffcompare 0.12.10 intron-chain level against the sample's own "
                              "annotation (the fig1_samples rule); *_ge2 = sensitivity on the reference chains carried "
                              f"exactly by >= {SAMPLES_MIN_READS} primary alignments"]
                       + scope_notes(cfg, keys) + [unit], data_dir=data_dir)


# ---- plan mode (figs_plan=1): the heavy units a build would run, with a cost estimate; runs nothing
GFFCOMPARE_S_PER_MB = 0.1    # human IsoSeq collapse genome-wide: 107 s for its ~1 GB GTF (measured 2026-09-25)
CHAIN_S_PER_M_RECORDS = (10.0, 20.0)   # pysam pass over primary+secondary records; gorilla 10.7 M in ~1.5 min;
#                                        the 96 GB human BAM is I/O-bound (~1.4 KB per record)


def _mb(path) -> float:
    p = Path(path)
    if not p.exists():
        return 0.0
    return p.stat().st_size / 1e6 * (8.0 if str(p).endswith(".gz") else 1.0)


def plan(cfg: dict):
    """Print the heavy units `make.py data fig1` would run now (tab-separated), with estimated seconds."""
    main, everyone = assembly.benchmark_samples(cfg), assembly.all_samples(cfg)
    units, total = [], [0.0, 0.0]
    print(f"# {assembly.guided_status(cfg)}")
    for key in everyone:
        tools = (TOOLS if key in main else RUSTLE_TOOLS) + assembly.guided_tools(cfg, key)
        for tool in tools:
            if tool in assembly.RUSTLE_STAGE:
                state, reason = samples.status(cfg, key, assembly.RUSTLE_STAGE[tool])
                if state not in ("fresh", "adopt"):
                    units.append((key, f"{assembly.RUSTLE_STAGE[tool]} ({state}: {reason})", "make.py runs first", ""))
                    continue
            what = assembly.gffcompare_state(cfg, key, tool)
            if what:
                sec = 10 + GFFCOMPARE_S_PER_MB * _mb(assembly.arm_source(cfg, key, tool, check=False))
                units.append((key, f"{tool}: {what}", f"{sec:.0f}", f"{sec:.0f}"))
                total[0] += sec
                total[1] += sec
        sub, label, _ = support_scope(cfg, key)
        path = chain_support_path(cfg, key, label)
        if not path.exists() or not figlib.fresh(path, samples.get(cfg, key)["bam"]):
            mapped = _idxstats(samples.get(cfg, key)["bam"])
            ev = assembly.evaluation_contigs(cfg, key)
            todo = [c for c in mapped if (ev is None or c in ev) and (sub is None or c in sub)]
            parts = Path(str(path) + ".parts")
            todo = [c for c in todo if not (parts / f"{c}.tsv").exists()]
            n = sum(mapped[c] for c in todo) / 1e6
            lo, hi = (n * r for r in CHAIN_S_PER_M_RECORDS)
            units.append((key, f"exact-chain support: {len(todo)} contig(s), {n:.1f} M records", f"{lo:.0f}",
                          f"{hi:.0f}"))
            total[0] += lo
            total[1] += hi
    print("sample\tunit\test_s_low\test_s_high")
    for u in units:
        print("\t".join(map(str, u)))
    print(f"TOTAL\t{len(units)} unit(s)\t{total[0]:.0f}\t{total[1]:.0f}")


def _idxstats(bam) -> dict:
    import subprocess
    out = subprocess.run(["samtools", "idxstats", str(bam)], capture_output=True, text=True, check=True).stdout
    return {l.split("\t")[0]: int(l.split("\t")[2]) for l in out.splitlines() if l.strip() and l[0] != "*"}


# ------------------------------------------------------------------------------------------------ plotting
def _f(x):
    return float(x) if x not in (None, "") else None


def place_labels(ax, points, texts, *, side="auto", dx_pt=5.0, min_gap_pt=7.5, fontsize=6.3, color=None,
                 leader=True, align=False):
    """Direct labels next to data points (data coords), spread vertically so none overlap (display units).
    side "auto" puts a label right of its point unless it would cross the axes' right edge (then left, if it fits
    there). Labels are ink, never the series colour; a thin leader joins a label that had to move. `align` starts
    every label of a side at the same x (the outermost point), for points that are dodged horizontally."""
    if not points:
        return []
    if side == "auto":  # right of the point unless the label would cross the axes' right edge
        renderer = ax.figure.canvas.get_renderer()
        bbox = ax.get_window_extent(renderer)
        left_edge, right_edge = bbox.x0, bbox.x1
        pt = ax.figure.dpi / 72.0
        groups = {"right": [], "left": []}
        for i, p in enumerate(points):
            probe = ax.text(0, 0, texts[i], fontsize=fontsize)
            width = probe.get_window_extent(renderer).width
            probe.remove()
            xd = ax.transData.transform(p)[0]
            fits_right = xd + dx_pt * pt + width <= right_edge + 2 * pt
            fits_left = xd - dx_pt * pt - width >= left_edge - 2 * pt
            groups["right" if fits_right or not fits_left else "left"].append(i)
        out = []
        for sd, idx in groups.items():
            out += place_labels(ax, [points[i] for i in idx], [texts[i] for i in idx], side=sd, dx_pt=dx_pt,
                                min_gap_pt=min_gap_pt, fontsize=fontsize, color=color, leader=leader, align=align)
        return out
    fig = ax.figure
    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    pt = fig.dpi / 72.0
    disp = [to_disp(p) for p in points]
    order = sorted(range(len(points)), key=lambda i: disp[i][1])
    ys = [disp[i][1] for i in order]
    gap = min_gap_pt * pt
    for _ in range(50):  # relax: push neighbours apart symmetrically
        moved = False
        for j in range(1, len(ys)):
            d = ys[j] - ys[j - 1]
            if d < gap:
                shift = (gap - d) / 2
                ys[j - 1] -= shift
                ys[j] += shift
                moved = True
        if not moved:
            break
    # keep the stack inside the axes' vertical extent (a label above the top would hit the panel title)
    bb = ax.get_window_extent(fig.canvas.get_renderer())
    half = 0.5 * fontsize * pt
    if ys and ys[-1] > bb.y1 - half:
        ys = [y - (ys[-1] - (bb.y1 - half)) for y in ys]
    if ys and ys[0] < bb.y0 + half:
        ys = [y + (bb.y0 + half - ys[0]) for y in ys]
    out = []
    sign = 1 if side == "right" else -1
    x_al = (max if side == "right" else min)(d[0] for d in disp)
    for rank, i in enumerate(order):
        x0, y0 = disp[i]
        tx, ty = to_data(((x_al if align else x0) + sign * dx_pt * pt, ys[rank]))
        t = ax.annotate(texts[i], xy=points[i], xytext=(tx, ty), textcoords="data", fontsize=fontsize,
                        color=color or figlib.INK, ha="left" if side == "right" else "right", va="center",
                        annotation_clip=False,
                        arrowprops=(dict(arrowstyle="-", color=figlib.INK_3, lw=0.4, shrinkA=0, shrinkB=2.5)
                                    if leader and abs(ys[rank] - y0) > 2.0 * pt else None))
        out.append(t)
    return out


def short_label(tool: str) -> str:
    """Direct-label text in dense panels: 'Rustle (prim.)' = Rustle, primaries-only seeding (caption defines it)."""
    return figlib.TOOL_LABEL[tool].replace(" (primaries only)", " (prim.)")


def _placeholder(ax, text):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=6.5, color=figlib.INK_3,
            wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.grid(False)


def _pct(v, decimals=0):
    return f"{100 * v:.{decimals}f}%"


MARKER_PT = 5.0


def _boxes_overlap(a, b, pad=0.0):
    return not (a[2] + pad <= b[0] or b[2] + pad <= a[0] or a[3] + pad <= b[1] or b[3] + pad <= a[1])


def place_scatter_labels(ax, points, texts, *, fontsize=6.0, obstacles=(), other_points=(), marker_pt=MARKER_PT):
    """Direct labels for a scatter: for each point try right, right-below, right-above, left, left-below,
    left-above, below, above (in that order) and keep the first spot that stays inside the axes and overlaps no
    marker, no earlier label and no box in `obstacles` (display coords). A label that left its row gets a leader."""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    pt = fig.dpi / 72.0
    disp = [ax.transData.transform(p) for p in points]
    rad = (marker_pt / 2 + 1.0) * pt
    marks = [(x - rad, y - rad, x + rad, y + rad) for x, y in disp]
    others = [(x - rad, y - rad, x + rad, y + rad) for x, y in (ax.transData.transform(q) for q in other_points)]
    axbb = ax.get_window_extent(r)
    placed = list(obstacles)
    cands = [(6, 0, "left"), (6, -9, "left"), (6, 9, "left"), (-6, 0, "right"), (-6, -9, "right"),
             (-6, 9, "right"), (0, -10, "center"), (0, 10, "center"), (6, -16, "left"), (-6, -16, "right")]
    out = []
    for i, (p, txt) in enumerate(zip(points, texts)):
        probe = ax.text(0, 0, txt, fontsize=fontsize)
        bb = probe.get_window_extent(r)
        probe.remove()
        w, h = bb.width, bb.height
        best = None
        for rank, (dx, dy, ha) in enumerate(cands):
            x0, y0 = disp[i][0] + dx * pt, disp[i][1] + dy * pt
            left = x0 if ha == "left" else x0 - w if ha == "right" else x0 - w / 2
            box = (left, y0 - h / 2, left + w, y0 + h / 2)
            outside = (box[0] < axbb.x0 + pt) + (box[2] > axbb.x1 - pt) + (box[1] < axbb.y0 + pt) + \
                (box[3] > axbb.y1 - pt)
            hits = sum(_boxes_overlap(box, m) for j, m in enumerate(marks) if j != i) + \
                sum(_boxes_overlap(box, m) for m in others) + \
                sum(_boxes_overlap(box, q, pad=1.0 * pt) for q in placed)
            score = (10 * outside + hits, rank)   # first clean spot; otherwise the least-bad one
            if best is None or score < best[0]:
                best = (score, (dx, dy, ha), box)
            if score[0] == 0:
                break
        placed.append(best[2])
        dx, dy, ha = best[1]
        out.append(ax.annotate(txt, xy=p, xytext=(dx, dy), textcoords="offset points", fontsize=fontsize,
                               color=figlib.INK, ha=ha, va="center", annotation_clip=False,
                               arrowprops=(dict(arrowstyle="-", color=figlib.INK_3, lw=0.4, shrinkA=0, shrinkB=3.0)
                                           if dy else None)))
    return out


def _scatter_panel(ax, g: dict, tools: list[str], n_ref_me: int | None):
    """(a/d) intron-chain precision vs sensitivity, one marker per arm. When the two Rustle markers overlap
    (gorilla: 0.2 / 0.3 points apart) the open primaries-only marker is drawn on top and a callout names both
    with their values, so neither is hidden."""
    from matplotlib.ticker import MultipleLocator, PercentFormatter

    xs = {t: _f(g[t]["sn"]) for t in tools}
    ys = {t: _f(g[t]["pr"]) for t in tools}
    for t in tools:  # open (primaries-only) marker above the filled one
        ax.plot([xs[t]], [ys[t]], markersize=MARKER_PT,
                zorder=3 if figlib.TOOL_MARKER_FILLED.get(t, True) else 3.5, **figlib.tool_marker_kwargs(t))
    lo_x, hi_x = min(xs.values()), max(xs.values())
    xpad = max(0.02, (hi_x - lo_x) * 0.25)
    ax.set_xlim(max(0.0, lo_x - xpad * 1.5), hi_x + xpad * 1.5)
    ax.set_ylim(0, max(ys.values()) * 1.18)
    # whole-percent ticks only: a tick at 27.5% printed as "28%" misreads every marker
    ax.xaxis.set_major_locator(MultipleLocator(0.02 if hi_x - lo_x < 0.08 else 0.05))
    ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax.grid(axis="both")
    ax.set_xlabel("Intron-chain sensitivity")
    ax.set_ylabel("Intron-chain precision")
    if n_ref_me:
        ax.set_title(f"n = {n_ref_me:,} multi-exon reference transcripts", fontsize=6.0, color=figlib.INK_2, pad=3)
    pair = [t for t in ("rustle", "rustle_primary") if t in tools]
    overlap = False
    if len(pair) == 2:
        p0, p1 = (ax.transData.transform((xs[t], ys[t])) for t in pair)
        overlap = math.hypot(*(p0 - p1)) / (ax.figure.dpi / 72.0) < 1.1 * MARKER_PT
    obstacles = []
    if overlap:
        obstacles.append(_pair_callout(ax, pair, xs, ys, marker_points=[(xs[t], ys[t]) for t in tools]))
    main = [t for t in tools if not (overlap and t in pair)]
    place_scatter_labels(ax, [(xs[t], ys[t]) for t in main], [short_label(t) for t in main], obstacles=obstacles,
                         other_points=[(xs[t], ys[t]) for t in tools if t not in main])


def _pair_callout(ax, pair, xs, ys, marker_points=()):
    """A framed key for an overlapping pair of markers (glyph + name per arm, then how far apart they are), put
    in the first free spot of a few candidate corners and joined to the pair by a leader. Returns its display bbox
    (an obstacle for the other labels)."""
    from matplotlib.patches import FancyBboxPatch

    fig = ax.figure
    r = fig.canvas.get_renderer()
    pt = fig.dpi / 72.0
    cx, cy = sum(xs[t] for t in pair) / 2, sum(ys[t] for t in pair) / 2
    ax_to = ax.transAxes.inverted().transform
    step = 8.0 * pt / ax.get_window_extent(r).height        # one text line, in axes fraction
    rad = (MARKER_PT / 2 + 1.5) * pt
    marks = [(x - rad, y - rad, x + rad, y + rad) for x, y in (ax.transData.transform(q) for q in marker_points)]
    dsn = abs(xs[pair[0]] - xs[pair[1]]) * 100
    dpr = abs(ys[pair[0]] - ys[pair[1]]) * 100
    spots = ((0.52, 0.62), (0.44, 0.62), (0.52, 0.44), (0.06, 0.62), (0.06, 0.30))
    for i_spot, (x0, y0) in enumerate(spots):
        arts = []
        for k, t in enumerate(["rustle_primary", "rustle"]):
            y = y0 - k * step
            arts += ax.plot([x0], [y], transform=ax.transAxes, markersize=4.2, clip_on=False,
                            **figlib.tool_marker_kwargs(t))
            arts.append(ax.text(x0 + 0.035, y, short_label(t), ha="left", va="center", fontsize=5.8,
                                color=figlib.INK, transform=ax.transAxes))
        arts.append(ax.text(x0 - 0.015, y0 - 2 * step, f"Δ sens {dsn:.1f} / prec {dpr:.1f} pts", ha="left",
                            va="center", fontsize=5.5, color=figlib.INK_3, transform=ax.transAxes))
        boxes = [a.get_window_extent(r) for a in arts]
        bb = (min(b.x0 for b in boxes) - 2.5 * pt, min(b.y0 for b in boxes) - 2.0 * pt,
              max(b.x1 for b in boxes) + 2.5 * pt, max(b.y1 for b in boxes) + 2.0 * pt)
        axbb = ax.get_window_extent(r)
        if bb[2] <= axbb.x1 + 2 * pt and not any(_boxes_overlap(bb, m) for m in marks):
            break
        if i_spot == len(spots) - 1:   # no clean spot: keep the last one rather than an empty frame
            break
        for a in arts:
            a.remove()
    (bx0, by0), (bx1, by1) = ax_to((bb[0], bb[1])), ax_to((bb[2], bb[3]))
    ax.add_patch(FancyBboxPatch((bx0, by0), bx1 - bx0, by1 - by0, boxstyle="round,pad=0,rounding_size=0.02",
                                transform=ax.transAxes, facecolor=figlib.SURFACE, edgecolor=figlib.INK_3, lw=0.4,
                                zorder=1.5, clip_on=False))
    ax.annotate("", xy=(cx, cy), xytext=(0.5 * (bx0 + bx1), by1), xycoords="data", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-", color=figlib.INK_3, lw=0.4, shrinkA=0, shrinkB=4.0))
    return bb


def _bars_panel(ax, g: dict, tools: list[str]):
    """(b/e) matching intron chains per arm; the arm's query transcripts (n) under each name."""
    vals = [int(float(g[t]["matching"])) for t in tools]
    yy = list(range(len(tools)))[::-1]
    for t, y, v in zip(tools, yy, vals):
        ax.barh(y, v, height=0.68, **figlib.tool_bar_kwargs(t))
        ax.text(v, y, f" {v:,}", va="center", ha="left", fontsize=6.0, color=figlib.INK_2)
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{short_label(t)}\nn = {int(float(g[t]['n_query'])):,}" for t in tools], fontsize=6.2)
    ax.set_xlim(0, max(vals) * 1.38)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x")
    ax.xaxis.set_major_formatter(lambda v, _: f"{v / 1000:g}k" if v else "0")
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Matching intron chains")
    ax.set_title("n = transcripts of the method", fontsize=6.0, color=figlib.INK_2, pad=3)


# (c/f) small horizontal dodge per arm, so coinciding points stay visible (the two Rustle arms share a slot: the
# primaries-only arm is a larger open ring AROUND the filled Rustle marker, so both show where they coincide)
SUPPORT_DODGE = {"rustle": -0.15, "rustle_primary": -0.15, "stringtie": -0.05, "flair": 0.05, "isoseq": 0.15}
SUPPORT_MS = {"rustle": 3.4, "rustle_primary": 6.4}   # filled dot inside the open ring; other arms 3.8
SAME_END_PTS = 0.5   # end values of the two Rustle arms closer than this (points) get one label


def _support_panel(ax, s: list[dict], scope: str, species: str):
    """(c/f) sensitivity by exact-chain read support; the >=2 column is Rustle's minimum read support."""
    from matplotlib.ticker import PercentFormatter

    stools = [t for t in TOOLS if any(r["tool"] == t for r in s)]
    if not stools:
        _placeholder(ax, f"{species}: exact-chain read support\nnot computed yet (make.py data fig1)")
        return
    ks = sorted({int(r["min_reads"]) for r in s})
    xpos = {k: i for i, k in enumerate(ks)}
    ends = {}
    for t in stools:
        pts = sorted((int(r["min_reads"]), _f(r["sn"])) for r in s if r["tool"] == t and r["sn"] != "")
        kw = figlib.tool_marker_kwargs(t)
        kw["linestyle"] = "--" if t == "rustle_primary" else "-"
        if t in SUPPORT_MS:
            kw["markeredgewidth"] = 1.0 if t == "rustle_primary" else 0.8
        dx = SUPPORT_DODGE.get(t, 0.0)
        ax.plot([xpos[k] + dx for k, _ in pts], [v for _, v in pts], lw=1.0, markersize=SUPPORT_MS.get(t, 3.8),
                zorder={"rustle_primary": 3.0, "rustle": 3.6}.get(t, 3.3), **kw)
        ends[t] = (xpos[pts[-1][0]] + dx, pts[-1][1])
    labels = [(t, short_label(t)) for t in stools]
    if "rustle" in ends and "rustle_primary" in ends and \
            abs(ends["rustle"][1] - ends["rustle_primary"][1]) * 100 < SAME_END_PTS:
        # the two Rustle arms end at the same value: one label names both (two labels would slide the whole stack
        # one slot off its markers)
        labels = [(t, "Rustle (both)" if t == "rustle" else lab) for t, lab in labels if t != "rustle_primary"]
    ax.set_xticks(range(len(ks)))
    nref = {int(r["min_reads"]): int(r["n_ref"]) for r in s if r["tool"] == stools[0]}
    def _n(v):  # compact reference-chain counts so the stratum labels never collide
        return f"{v / 1000:.1f}k" if v >= 10000 else f"{v:,}"
    ax.set_xticklabels([("all\n" if k == 0 else f"≥{k}\n") + _n(nref[k]) for k in ks], fontsize=5.8)
    if PAIRED_SHOWN in xpos:  # Rustle's minimum read support
        ax.axvspan(xpos[PAIRED_SHOWN] - 0.22, xpos[PAIRED_SHOWN] + 0.22, color="#f1f0eb", zorder=0, lw=0)
        ax.text(xpos[PAIRED_SHOWN], 0.03, "Rustle's\nminimum", ha="center", va="bottom", fontsize=5.5,
                color=figlib.INK_3, linespacing=0.95)
    ax.set_xlim(-0.3, len(ks) - 1 + 1.9)
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    place_labels(ax, [ends[t] for t, _ in labels], [lab for _, lab in labels], side="right", dx_pt=6,
                 fontsize=5.8, min_gap_pt=7.0, align=True)
    sc = next((r["scope"] for r in s), scope)
    ax.set_title("subset: " + sc.replace(",", ", ") if sc != scope else "same chains for every method",
                 fontsize=6.0, color=figlib.INK_3 if sc != scope else figlib.INK_2, pad=3)
    ax.set_xlabel("Primary reads carrying the exact reference chain\n(reference chains in the stratum)", fontsize=6.3)
    ax.set_ylabel("Intron-chain sensitivity")


def _p_text(s: str) -> str:
    lp = parse_p(s)
    if lp >= -0.0002:
        return "1"
    if lp >= -3:
        return f"{10 ** lp:.2g}"
    if lp < -300:
        return "<1e-300"
    e = math.floor(lp)
    return f"{10 ** (lp - e):.1f}e{e}"


def _paired_panel(fig, rect_forest, rect_table, rows: list[dict], species: str, letter: str):
    """(g/h) Rustle minus each other arm on the same >=2-read chains: point = difference in sensitivity,
    bar = Tango 95% interval; the table beside it gives the 2x2 counts and the exact McNemar p."""
    ax = fig.add_axes(rect_forest)
    tab = fig.add_axes(rect_table, sharey=ax)
    others = [t for t in TOOLS if t != "rustle" and any(r["other"] == t for r in rows)]
    if not others:
        _placeholder(ax, f"{species}: paired table\nnot computed yet")
        tab.axis("off")
        return
    r0 = rows[0]
    scope = r0["scope"]
    scope_txt = "genome-wide" if scope == "genome-wide" else scope.replace(",", ", ")
    head_y = rect_forest[1] + rect_forest[3] + 0.018
    fig.text(rect_forest[0] - 0.105, head_y + 0.004, letter, fontsize=9, fontweight="bold", ha="left",
             va="bottom")
    fig.text(rect_forest[0] - 0.080, head_y, f"{figlib.SPECIES_LABEL[species].split(' (')[0]} · {scope_txt}"
             f" · ≥{PAIRED_SHOWN} reads · n = {int(r0['n_ref']):,} chains", fontsize=6.5, ha="left", va="bottom",
             color=figlib.INK_3 if scope.endswith(" only") else figlib.INK)
    n = len(others)
    yy = {t: i for i, t in enumerate(reversed(others))}
    lo_all, hi_all = [0.0], [0.0]
    for t in others:
        r = next(r for r in rows if r["other"] == t)
        d, lo, hi = 100 * _f(r["diff"]), 100 * _f(r["diff_lo"]), 100 * _f(r["diff_hi"])
        lo_all.append(lo)
        hi_all.append(hi)
        ax.plot([lo, hi], [yy[t], yy[t]], color=figlib.INK_2, lw=1.0, solid_capstyle="butt", zorder=2)
        for x in (lo, hi):
            ax.plot([x, x], [yy[t] - 0.14, yy[t] + 0.14], color=figlib.INK_2, lw=0.8, zorder=2)
        ax.plot([d], [yy[t]], markersize=4.2, zorder=3 if figlib.TOOL_MARKER_FILLED.get(t, True) else 3.5,
                **figlib.tool_marker_kwargs(t))
    ax.axvline(0, color=figlib.INK_3, lw=0.6, zorder=1)
    span = max(hi_all) - min(lo_all)
    ax.set_xlim(min(lo_all) - 0.06 * span - 0.8, max(hi_all) + 0.06 * span + 0.8)
    ax.set_ylim(-0.6, n - 0.1 + 0.55)   # the top band holds the table's column headers
    ax.set_yticks([yy[t] for t in others])
    ax.set_yticklabels([f"vs {short_label(t)}" for t in others], fontsize=6.2)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x")
    ax.set_xlabel("Rustle minus method (sensitivity,\npoints; Tango 95% interval)", fontsize=6.3)
    tab.axis("off")
    cols = [("both", 0.20), ("Rustle\nonly", 0.44), ("other\nonly", 0.64), ("McNemar\np", 0.99)]
    for name, x in cols:
        tab.text(x, n - 0.1 + 0.05, name, ha="right", va="bottom", fontsize=5.8, color=figlib.INK_2,
                 transform=tab.get_yaxis_transform(), linespacing=0.95)
    for t in others:
        r = next(r for r in rows if r["other"] == t)
        vals = [f"{int(r['n_both']):,}", f"{int(r['n_arm_only']):,}", f"{int(r['n_other_only']):,}",
                _p_text(r["mcnemar_p"])]
        for (name, x), v in zip(cols, vals):
            tab.text(x, yy[t], v, ha="right", va="center", fontsize=6.0, color=figlib.INK,
                     transform=tab.get_yaxis_transform())


FIG_W = 7.05   # inches; the tight bounding box must stay <= 183 mm
FIG_H = 6.9


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt

    figlib.use_style()
    gc = figlib.read_table("fig1_gffcompare", data_dir)
    try:
        sup = figlib.read_table("fig1_support", data_dir)
    except FileNotFoundError:
        sup = []
    try:
        paired = figlib.read_table("fig1_paired", data_dir)
    except FileNotFoundError:
        paired = []

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.56], hspace=0.86,
                             left=0.085, right=0.915, top=0.915, bottom=0.075)
    letters = iter("abcdef")
    for row, species in enumerate(SPECIES):
        gs = outer[row].subgridspec(1, 3, width_ratios=[1.15, 0.95, 1.1], wspace=0.66)
        g = {r["tool"]: r for r in gc if r["species"] == species and r["level"] == "intron_chain"}
        tools = [t for t in TOOLS if t in g]
        scope = next((r["scope"] for r in gc if r["species"] == species), "")
        scope_txt = "genome-wide" if scope == "genome-wide" else scope.replace(",", ", ")
        y_top = outer[row].get_position(fig).y1
        fig.text(0.01, y_top + 0.045, f"{figlib.SPECIES_LABEL[species]}  ·  {scope_txt}", fontsize=7.5,
                 fontweight="bold", ha="left", va="bottom", color=figlib.INK)
        n_ref_me = next((int(float(g[t]["n_ref_multiexon"])) for t in tools
                         if g[t].get("n_ref_multiexon") not in (None, "")), None)

        ax = fig.add_subplot(gs[0])
        figlib.panel_label(ax, next(letters), x=-0.30, y=1.08)
        if tools:
            _scatter_panel(ax, g, tools, n_ref_me)
        else:
            _placeholder(ax, f"{species}: no gffcompare rows yet\n(make.py data fig1)")

        ax = fig.add_subplot(gs[1])
        figlib.panel_label(ax, next(letters), x=-0.66, y=1.08)
        if tools:
            _bars_panel(ax, g, tools)
        else:
            _placeholder(ax, "no rows yet")

        ax = fig.add_subplot(gs[2])
        figlib.panel_label(ax, next(letters), x=-0.30, y=1.08)
        _support_panel(ax, [r for r in sup if r["species"] == species and r["count_mode"] == "reads"], scope,
                       species)

    # (g/h) paired comparison on the >=2-read stratum, one panel per species (never pooled); explicit geometry
    box = outer[2].get_position(fig)
    fig.text(0.01, box.y1 + 0.058, f"Paired comparison on the ≥{PAIRED_SHOWN}-read stratum (the same reference "
             "chains for every method)", fontsize=7.5, fontweight="bold", ha="left", va="bottom", color=figlib.INK)
    x_left, x_right = 0.125, 0.985
    gap = 0.17   # room for the right panel's row labels, clear of the left panel's table
    half = (x_right - x_left - gap) / 2
    for i, species in enumerate(SPECIES):
        rows = [r for r in paired if r["species"] == species and r["count_mode"] == "reads"
                and int(r["min_reads"]) == PAIRED_SHOWN]
        x0 = x_left + i * (half + gap)
        f_w = half * 0.45
        rect_f = [x0, box.y0, f_w, box.height]
        rect_t = [x0 + f_w + 0.01, box.y0, half - f_w - 0.01, box.height]
        _paired_panel(fig, rect_f, rect_t, rows, species, "gh"[i])

    # the comparison mode, stated on the figure (every panel above compares annotation-free methods only)
    fig.text(0.01, 0.992, figlib.mode_lines("fig1g_guided", (Path(data_dir) / "fig1_guided.tsv").exists()),
             fontsize=6.0, color=figlib.INK_2, ha="left", va="bottom", linespacing=1.3)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "fig1_intron_chain", out_dir)
    plt.close(fig)
    return paths + plot_samples(data_dir, out_dir) + plot_guided(data_dir, out_dir)


GUIDED_TITLES = {"sn": "Intron-chain sensitivity\n(multi-exon reference transcripts)",
                 "pr": "Intron-chain precision\n(multi-exon transcripts of the tool)"}


def plot_guided(data_dir: Path, out_dir: Path) -> list:
    """fig1g_guided: the annotation-guided StringTie/FLAIR runs (fig1_guided), one row per sample; drawn only when
    that table exists (i.e. when the user registered a guided GTF). No annotation-free method and no Rustle row."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import PercentFormatter

    try:
        rows = figlib.read_table("fig1_guided", data_dir)
    except FileNotFoundError:
        return []
    order = list(dict.fromkeys(r["sample"] for r in rows))
    by = {(r["sample"], r["tool"]): r for r in rows}
    tools = [t for t in figlib.GUIDED_TOOL_ORDER if any(r["tool"] == t for r in rows)]
    n = len(order)
    fig_h = 1.35 + 0.34 * n
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.2, fig_h))
    gs = fig.add_gridspec(1, len(SAMPLES_COLS), left=0.25, right=0.975, top=1 - 0.80 / fig_h, bottom=0.62 / fig_h,
                          wspace=0.22)
    yy = {sid: n - 1 - i for i, sid in enumerate(order)}
    for j, (col, title) in enumerate(SAMPLES_COLS):
        ax = fig.add_subplot(gs[j])
        figlib.panel_label(ax, "abc"[j], x=-0.05 if j else -0.95, y=1.0 + 0.30 / (fig_h * 0.7))
        for sid in order:
            for tool in tools:
                r = by.get((sid, tool))
                if not r or r.get(col) in ("", None):
                    continue
                ax.plot([float(r[col])], [yy[sid]], markersize=5.0, clip_on=False, **figlib.tool_marker_kwargs(tool))
                ax.annotate(f"{100 * float(r[col]):.1f}%", (float(r[col]), yy[sid]),
                            xytext=(6, 4 if tool == tools[0] else -6), textcoords="offset points", fontsize=5.6,
                            va="center", ha="left", color=figlib.INK_2)
        ax.set_xlim(0, 1.0)
        ax.set_ylim(-0.6, n - 0.4)
        ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x")
        ax.tick_params(axis="y", length=0)
        ax.set_title(GUIDED_TITLES.get(col, title), fontsize=6.3, pad=4)
        ax.set_yticks([yy[sid] for sid in order])
        if j == 0:
            ax.set_yticklabels([next(r["label"] for r in rows if r["sample"] == sid) for sid in order], fontsize=6.2)
        else:
            ax.tick_params(axis="y", labelleft=False)
    handles = [Line2D([], [], markersize=4.6, **figlib.tool_marker_kwargs(t), label=figlib.TOOL_LABEL[t])
               for t in tools]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 1 - 0.03 / fig_h), ncol=len(tools) or 1,
               fontsize=6.0, frameon=False, handletextpad=0.4, columnspacing=1.6, borderaxespad=0.0)
    fig.text(0.01, 1 - 0.30 / fig_h, "Annotation-guided runs only (never compared with the annotation-free methods "
             "of Fig. 1). The guided tools were given the annotation they are scored against.\nRustle has no "
             "annotation-guided transcript assembly, so it has no row here.", fontsize=5.8, color=figlib.INK_2,
             ha="left", va="top")
    scope = sorted({r["scope"] for r in rows})
    fig.text(0.01, 0.12 / fig_h, "Each sample against its own RefSeq annotation (" + ", ".join(scope) + "). Never "
             "pooled across samples.", fontsize=5.8, color=figlib.INK_2, ha="left", va="bottom")
    figlib.stamp_provisional(fig, ["fig1_guided"], data_dir)
    paths = figlib.save(fig, "fig1g_guided", out_dir)
    plt.close(fig)
    return paths


# ------------------------------------------------------------------------------------------------ supplementary
SAMPLES_COLS = [("sn", "Intron-chain sensitivity\n(every multi-exon reference transcript)"),
                ("pr", "Intron-chain precision\n(every multi-exon Rustle transcript)"),
                ("sn_ge2", "Sensitivity on reference chains\ncarried exactly by ≥ 2 primary reads")]


def plot_samples(data_dir: Path, out_dir: Path) -> list:
    """fig1s_samples: Rustle's two configurations on every sample, each against its own annotation (one row per
    sample, grouped by species; never pooled). Drawn only when figures/data/fig1_samples.tsv exists."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import PercentFormatter

    try:
        rows = figlib.read_table("fig1_samples", data_dir)
    except FileNotFoundError:
        return []
    order = list(dict.fromkeys(r["sample"] for r in rows))
    by = {(r["sample"], r["tool"]): r for r in rows}
    n = len(order)
    fig_h = 1.05 + 0.34 * n
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.2, fig_h))
    gs = fig.add_gridspec(1, len(SAMPLES_COLS), left=0.25, right=0.975, top=1 - 0.62 / fig_h, bottom=0.52 / fig_h,
                          wspace=0.22)
    yy = {sid: n - 1 - i for i, sid in enumerate(order)}
    species_of = {sid: next(r["species"] for r in rows if r["sample"] == sid) for sid in order}
    breaks = [yy[b] + 0.5 for a, b in zip(order, order[1:]) if species_of[a] != species_of[b]]
    for j, (col, title) in enumerate(SAMPLES_COLS):
        ax = fig.add_subplot(gs[j])
        figlib.panel_label(ax, "abc"[j], x=-0.05 if j else -0.95, y=1.0 + 0.30 / (fig_h * 0.7))
        for y in breaks:
            ax.axhline(y, color=figlib.GRID, lw=0.6, zorder=0)
        for sid in order:
            for tool in ("rustle_primary", "rustle"):
                r = by.get((sid, tool))
                if not r or r.get(col) in ("", None):
                    continue
                kw = figlib.tool_marker_kwargs(tool)
                kw["markeredgewidth"] = 1.0 if tool == "rustle_primary" else 0.8
                ax.plot([float(r[col])], [yy[sid]], markersize=6.4 if tool == "rustle_primary" else 3.6,
                        zorder=3 if tool == "rustle" else 2.8, clip_on=False, **kw)
            r = by.get((sid, "rustle"))
            if r and r.get(col) not in ("", None):
                v = float(r[col])
                vs = [float(x[col]) for x in (r, by.get((sid, "rustle_primary"))) if x and x.get(col) not in ("", None)]
                right = max(vs) < 0.8    # print beside the outermost marker, inside the axes
                ax.annotate(f"{100 * v:.1f}%", (max(vs) if right else min(vs), yy[sid]), xytext=(7 if right else -7, 0),
                            textcoords="offset points", fontsize=5.8, va="center", ha="left" if right else "right",
                            color=figlib.INK_2)
        ax.set_xlim(0, 1.0)
        ax.set_ylim(-0.6, n - 0.4)
        ax.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="x", labelsize=6)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x")
        ax.tick_params(axis="y", length=0)
        ax.set_title(title, fontsize=6.3, pad=4)
        if j == 0:
            ax.set_yticks([yy[sid] for sid in order])
            ax.set_yticklabels([f"{by[(sid, 'rustle')]['label']}\n{int(float(by[(sid, 'rustle')]['n_query'])):,} "
                                f"transcripts" if (sid, "rustle") in by else sid for sid in order], fontsize=6.2)
        else:
            ax.set_yticks([yy[sid] for sid in order])
            ax.tick_params(axis="y", labelleft=False)
    handles = [Line2D([], [], markersize=4.2, **figlib.tool_marker_kwargs("rustle"),
                      label="Rustle (default: primary alignments + secondary alignments within 2% of the best score)"),
               Line2D([], [], markersize=5.5,
                      **{**figlib.tool_marker_kwargs("rustle_primary"), "markeredgewidth": 1.0},
                      label="Rustle, primary alignments only")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 1 - 0.03 / fig_h), ncol=2, fontsize=6.0,
               frameon=False, handletextpad=0.4, columnspacing=1.6, borderaxespad=0.0)
    scope = sorted({r["scope"] for r in rows})
    fig.text(0.01, 0.12 / fig_h, "Each sample against its own RefSeq annotation (" + ", ".join(scope) + "); printed: "
             "Rustle (default). Never pooled across samples.", fontsize=5.8, color=figlib.INK_2, ha="left",
             va="bottom")
    figlib.stamp_provisional(fig, META["supplementary_tables"], data_dir)
    paths = figlib.save(fig, "fig1s_samples", out_dir)
    plt.close(fig)
    return paths


# ------------------------------------------------------------------------------------------------ caption numbers
def caption_numbers(data_dir=None) -> str:
    """Every number the caption quotes, read from the tables (`python3 figures/fig_intron_chain.py summary`): refresh
    the caption after `make.py data fig1` without recomputing anything by hand (light)."""
    data_dir = Path(data_dir or figlib.DATA_DIR)
    out = []
    gc = figlib.read_table("fig1_gffcompare", data_dir)
    for n in figlib.table_meta("fig1_gffcompare", data_dir).get("note", []):
        out.append(f"note  {n}")
    for sp in dict.fromkeys(r["species"] for r in gc):
        rows = [r for r in gc if r["species"] == sp and r["level"] == "intron_chain"]
        out.append(f"== {sp} ({rows[0]['scope'] if rows else '?'})")
        for r in rows:
            out.append(f"gffcompare {r['tool']:>15}: sn {100 * float(r['sn']):.1f}%  pr {100 * float(r['pr']):.1f}%  "
                       f"matching {r['matching']}  query {int(float(r['n_query'])):,} (multi-exon "
                       f"{int(float(r['n_query_multiexon'])):,})  ref multi-exon {int(float(r['n_ref_multiexon'])):,}")
    try:
        sup = figlib.read_table("fig1_support", data_dir)
    except FileNotFoundError:
        sup = []
    for r in sup:
        if r["sn"]:
            out.append(f"support {r['species']} {r['count_mode']:>13} >={r['min_reads']} {r['tool']:>15}: "
                       f"{r['n_matched']}/{r['n_ref']} = {100 * float(r['sn']):.1f}%")
    try:
        paired = figlib.read_table("fig1_paired", data_dir)
    except FileNotFoundError:
        paired = []
    for r in paired:
        out.append(f"paired {r['species']} {r['count_mode']:>13} >={r['min_reads']} vs {r['other']:>15}: n {r['n_ref']} "
                   f"both {r['n_both']} rustle-only {r['n_arm_only']} other-only {r['n_other_only']}  diff "
                   f"{100 * float(r['diff']):+.1f} ({100 * float(r['diff_lo']):+.1f} to {100 * float(r['diff_hi']):+.1f})"
                   f"  p {r['mcnemar_p']}")
    try:
        smp = figlib.read_table("fig1_samples", data_dir)
    except FileNotFoundError:
        smp = []
    for r in smp:
        out.append(f"sample {r['sample']:>15} {r['tool']:>15}: sn {100 * float(r['sn']):.1f}%  pr "
                   f"{100 * float(r['pr']):.1f}%  >=2 reads {r['n_matched_ge2']}/{r['n_ref_ge2']}"
                   + (f" = {100 * float(r['sn_ge2']):.1f}%" if r["sn_ge2"] else "")
                   + f"  transcripts {int(float(r['n_query'])):,}  ({r['scope']}; {r['annotation']})")
    try:
        gd = figlib.read_table("fig1_guided", data_dir)
    except FileNotFoundError:
        gd = []
        out.append("guided: not available (guided StringTie/FLAIR GTFs not supplied)")
    for r in gd:
        out.append(f"guided {r['sample']:>15} {r['tool']:>17}: sn {100 * float(r['sn']):.1f}%  pr "
                   f"{100 * float(r['pr']):.1f}%  >=2 reads {r['n_matched_ge2']}/{r['n_ref_ge2']}  (annotation-guided; "
                   "given the annotation it is scored against; no Rustle counterpart)")
    return "\n".join(out)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="fig1: caption numbers from the tables (light)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sm = sub.add_parser("summary", help="print every number the caption quotes, from figures/data")
    sm.add_argument("--data", default=str(figlib.DATA_DIR))
    a = ap.parse_args(argv)
    if a.cmd == "summary":
        print(caption_numbers(Path(a.data)))


if __name__ == "__main__":
    main()
