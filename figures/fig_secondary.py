"""fig_secondary — Figure 3: transcripts assembled in multi-mapping loci (secondary-alignment seeding).

Question (tested, not assumed): Rustle seeds loci with SECONDARY alignments that tie the molecule's genome-wide
best alignment score (AS >= 0.98 x best; the pipeline default since 2026-09-24, `--no-seed-secondaries` turns it
off). StringTie, FLAIR and IsoSeq collapse were run by the lab on the same BAM, secondaries included, each with its
own handling of them (StringTie reads secondary alignments). Does the seeding make Rustle reconstruct annotated
intron chains where reads multi-map, and in particular where no primary alignment reaches the transcript?

Unit: every multi-exon reference transcript on the evaluation contigs (assembly.evaluation_contigs).
Tie fraction of a transcript = n_tied / n_mol, where
  n_mol   = read molecules with a CANDIDATE placement whose ALIGNED BLOCK (a gapless M/=/X run: an intron `N`
            or a deletion `D` is not aligned sequence) overlaps one of the transcript's exons. A candidate
            placement is the molecule's primary record, or a secondary whose AS >= 0.98 x the molecule's
            genome-wide best AS (the pool the assembler seeds from; `is_candidate` mirrors its rule exactly).
            Supplementary (0x800) and unmapped (0x4) records never count; a molecule is counted once per
            transcript however many records it has there.
  n_tied  = those molecules whose genome-wide second-best AS >= 0.98 x best (the molecule has >= 2 candidate
            placements in the genome; from the driver's `as_table` product `<prefix>.molecules.tsv`).
Poor secondaries (minimap2 -N 50 -p 0.1 emits many at < 0.98 x best: 6.16 M of the 10.65 M records streamed on the
gorilla contigs) are candidate placements for no arm, so they do not count as "aligned to" a transcript.
A transcript is matched by an arm when a query transcript of that arm has gffcompare class '=' to it, or to a
reference transcript with the identical intron chain (gffcompare reports one ref_id per query, so duplicated
annotation chains would otherwise be split arbitrarily between the duplicates).

Read-sharing groups: high-tie transcripts (tie fraction > 0.5) that share any tied molecule form one group
(union-find over the whole streamed genome). Near-identical copies (e.g. a tandem array) are all rebuilt from the
SAME reads, so a gain is reported in transcripts AND in groups: the groups are the independent read evidence;
which copy is actually expressed is the copy-assignment question, not this figure's.

Primary support (review B1): `n_mol_primary` = molecules whose PRIMARY record reaches the transcript. A transcript
with n_mol_primary == 0 is reached by no primary alignment, so an assembler of primary alignments has no read
there. Those transcripts are one category ("no primary"), whatever their tie fraction; the tie bins hold only
transcripts with >= 1 primary molecule, where every arm has reads to work with.

Intervals are over CLUSTERS, not transcripts (copies rebuilt from the same reads are not independent): the
cluster of a transcript is its read-sharing group when it has one (tie fraction > 0.5), else its gene (isoforms
share molecules). 95% percentile cluster bootstrap (resampling clusters with replacement, the same resamples for
every arm); where no transcript or every transcript is matched the bootstrap is degenerate, and the interval is
Wilson's with n = the number of clusters.

Tables (figures/data/):
  fig3_ref_tie    one row per EXPRESSED (n_mol >= 1) multi-exon reference transcript: counts, tie fraction, bin,
                  read-sharing group, matched flags per arm
  fig3_bins       species x category (tie bin x primary support) x arm: n transcripts, n clusters, n matched,
                  n clusters with a match, fraction, cluster interval; `plotted` marks panel a's categories (the
                  tie bins with >= 1 primary molecule, and "no primary"); the others are the unstratified bins
                  and the >= 2-primary sensitivity row, kept for the caption
  fig3_gain       species x category: transcripts matched by Rustle but by none of the three baselines (StringTie,
                  FLAIR, IsoSeq collapse) and the reverse (matched by at least one baseline, not by Rustle), and by
                  Rustle but not Rustle primaries-only (and the reverse); in the high-tie categories also the number
                  of independent units behind each count (cluster_of: read-sharing group, else gene)
  fig3_example    one example locus (pick_example: among the transcripts reached by no primary alignment that
                  Rustle rebuilds exactly, the high-tie one with the most molecules): reference / arm exons and read
                  depth split into primary vs candidate-secondary alignments

The per-transcript counting streams the BAM contig by contig (heavy: gorilla genome-wide ~10.7 M records, 3.6 min;
human genome-wide 67.6 M records on the 96 GB BAM, I/O-bound, estimated 17-40 min): it runs as a logged foreground
subprocess (`python3 figures/fig_secondary.py count ...`) cached under work_dir(cfg, fig3)/<sample>/eval_<scope>/ (keyed
by the evaluation scope, like assembly.eval_dir), resumable per contig (`<counts>.parts/`: a part is reused when it is
newer than the BAM, the annotation and the best-AS table and was written by the same COUNT_VERSION); `--budget-s`
stops a call after the first contig that ends past the budget (exit 75; `figs_budget_s` / `fig3_budget_s` pass it).
When the evaluated contigs hold under half of the BAM's records (a development subset), it first collects the names
seen there, so only those rows of the best-AS table are held in memory; genome-wide it loads the whole table (human:
5.8 M reads with more than one alignment).

Scope: every sample with the three lab baselines (assembly.benchmark_samples), GENOME-WIDE: every contig the
sample's annotation covers (human: all but chrM, which the annotation leaves out).

Mode: every method is ANNOTATION-FREE (de novo): Rustle assemble (reads + genome), StringTie -L without -G, FLAIR
collapse without annotation (flair correct skipped), IsoSeq collapse; the figure and the table notes say so.
  fig3_guided_bins  (supplementary; ONLY when a benchmark sample has an annotation-guided StringTie / FLAIR GTF
                  registered, samples.tsv stringtie_guided_gtf / flair_guided_gtf; drawn as fig3g_guided) the fig3_bins
                  columns plus `mode` and `sample`, for the guided tools only, on the same transcripts, categories and
                  clusters (bins_rows_for). Never with the annotation-free methods; no Rustle row (Rustle has no
                  annotation-guided transcript assembly; docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md).
                  Until then the figure prints "guided comparison: not available (guided StringTie/FLAIR GTFs not
                  supplied)".
"""
from __future__ import annotations

import argparse
import bisect
import csv
import heapq
import math
import sys
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import assembly  # noqa: E402
import figlib  # noqa: E402

FIG = "fig3"
TIE_RATIO = 0.98  # the pipeline's RUSTLE_GTF_SECONDARY_AS_RATIO
TOOLS = ["rustle", "rustle_primary", "stringtie", "flair", "isoseq"]
SPECIES = ["gorilla", "human"]

# Tie-fraction bins: (label, lower bound exclusive, upper bound inclusive); the first bin is exactly 0.
# Justification (see captions/fig3.md): the distribution is a spike at 0 (most expressed transcripts have no
# tied molecule at all), a long thin tail of transcripts whose tied molecules are a minority (paralogous
# exons shared by an otherwise unique gene), and a second mode near 1 (every molecule is tied: recent
# duplicates). 0.5 separates "most molecules unique" from "most molecules tied"; > 0.9 isolates the mode at 1.
TIE_BINS = [
    ("0", None, 0.0),
    ("(0, 0.1]", 0.0, 0.1),
    ("(0.1, 0.5]", 0.1, 0.5),
    ("(0.5, 0.9]", 0.5, 0.9),
    ("> 0.9", 0.9, 1.0),
]
HIGH_BINS = ["(0.5, 0.9]", "> 0.9"]  # "multi-mapping" loci: most molecules tied
HIGH_TIE_MIN = 0.5  # lower (exclusive) bound of HIGH_BINS; high-tie transcripts get a read-sharing group

# Primary support of a transcript (n_mol_primary = molecules whose PRIMARY record reaches it). "0" = no primary
# alignment touches the transcript: an assembler of primary alignments has no read there.
PRIMARY_STRATA = {
    ">=1": lambda n: n >= 1,
    "0": lambda n: n == 0,
    ">=2": lambda n: n >= 2,
    "any": lambda n: True,
}
# panel a: the tie bins among transcripts with >= 1 primary molecule, then every transcript with none
PANEL_A = [(label, ">=1") for label, _, _ in TIE_BINS] + [("any", "0")]
# panel b: the high-tie categories of panel a
PANEL_B = [("(0.5, 0.9]", ">=1"), ("> 0.9", ">=1"), ("any", "0")]
BOOT_B = 2000  # cluster-bootstrap resamples

META = {
    "id": FIG,
    "title": "Transcripts assembled in multi-mapping loci",
    "claim": ("Annotation-free (de novo) comparison: Rustle assemble (reads + genome), StringTie -L without -G, FLAIR "
              "collapse without annotation (flair correct skipped), IsoSeq collapse; guided comparison: not available "
              "(guided StringTie/FLAIR GTFs not supplied). "
              "Gorilla OR6737, genome-wide (reference: the RefSeq annotation; a match is an exact intron chain, "
              "gffcompare '='): building loci also from the secondary alignments that score within 2% of the read's "
              "best alignment score lets Rustle rebuild 50 of the 839 multi-exon transcripts that no primary alignment "
              "reaches (23 read-sharing groups or genes; 33 of 794 without NC_073244.2, the contig the seeding default "
              "was decided on); no other method matches any (Rustle with primary alignments only by construction). "
              "Where a primary alignment reaches the transcript the methods are not separated: among transcripts "
              "whose reads are > 90% tied (a second alignment >= 98% of the best score), Rustle reproduces 291 of "
              "1,262 chains and IsoSeq collapse 248, which reaches more read-sharing groups (83 vs 79), with "
              "overlapping 95% intervals (16.8-30.2% vs 14.1-25.7%); against any baseline (StringTie, FLAIR or IsoSeq "
              "collapse) the one-sided matches are 56 Rustle-only vs 62 baseline-only there and 12 vs 57 in (0.5, "
              "0.9]. The 160 high-tie transcripts gained over Rustle with primary alignments only are spread over 66 "
              "read-sharing groups (largest: the CGB-like array LOC129528600-625, 20 gains from 4 reads). Seeding "
              "lowers intron-chain precision (Fig. 1): the 1,074 extra gorilla multi-exon transcripts add 162 "
              "matching chains (15%), the 1,017 extra human ones add 4 (0.4%). Human: chr20-22 in the current "
              "tables (too small to separate the methods); genome-wide after the rebuild."),
    "tables": ["fig3_ref_tie", "fig3_bins", "fig3_gain", "fig3_example"],
    # drawn as fig3g_guided only when an annotation-guided GTF is registered (the same `make.py data fig3`)
    "supplementary_tables": ["fig3_guided_bins"],
}


# ---------------------------------------------------------------- molecule table (driver's as_table product)
def molecule_tied(best: int, second: int) -> bool:
    return best > 0 and second >= 0 and second >= TIE_RATIO * best


def load_molecules(path, keep: set | None = None):
    """{name: best_as} for molecules with >= 2 records, and the set of TIED molecule names (restricted to the
    names in `keep` when given — memory bound for a region subset of a large BAM).

    Single-record molecules need no entry: their one record is the primary and they are never tied."""
    best: dict[str, int] = {}
    tied: set[str] = set()
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 4)
            if int(f[3]) < 2 or (keep is not None and f[0] not in keep):
                continue
            b, s = int(f[1]), int(f[2])
            best[f[0]] = b
            if molecule_tied(b, s):
                tied.add(f[0])
    return best, tied


def molecules_header_bam(path) -> str | None:
    try:
        with open(path) as fh:
            first = fh.readline().rstrip("\n").split("\t")
    except OSError:
        return None
    if not first or first[0] != "#as_table":
        return None
    for kv in first[1:]:
        if kv.startswith("bam="):
            return kv[4:]
    return None


def ensure_molecules(cfg: dict, species: str, *, force=False, budget=None) -> Path:
    """The genome-wide best-AS table the pipeline's `assemble` stage wrote (run cache product `molecules`,
    ${work}/runs/<id>/<id>.molecules.tsv); re-made with `as_table` only if it is absent or from another BAM. It must
    come from the FULL BAM: a region slice's table is not genome-wide."""
    import samples
    path = samples.product(cfg, species, "assemble", "molecules")
    bam_path = samples.get(cfg, species)["bam"]
    bam = str(Path(bam_path).resolve())
    # a deterministic function of the BAM whose header names that BAM: re-made only when absent or from another
    # BAM (`force` does not re-scan it)
    if molecules_header_bam(path) != bam:
        if budget is not None:
            budget.check(f"as_table of {species}")
        figlib.run([str(Path(cfg["bin"]) / "as_table"), "--bam", bam_path, "--out", str(path),
                    "--threads", cfg.get("threads", "4")], log=Path(str(path) + ".as_table.fig3.log"))
    return path


# ---------------------------------------------------------------- reference exon segments
def multi_exon(ref_tx: dict) -> dict:
    return {t: v for t, v in ref_tx.items() if len(v["exons"]) >= 2}


def exon_segments(txs: list[tuple[str, list]]):
    """Elementary segments of one contig: (starts, ends, covering transcript sets); gaps carry an empty set."""
    events: dict[int, list] = {}
    for tid, exons in txs:
        for s, e in exons:
            events.setdefault(s, []).append((1, tid))
            events.setdefault(e, []).append((-1, tid))
    pts = sorted(events)
    active: dict[str, int] = {}
    starts, ends, sets = [], [], []
    interned: dict[frozenset, frozenset] = {}
    for a, b in zip(pts, pts[1:]):
        for sign, tid in events[a]:
            active[tid] = active.get(tid, 0) + sign
            if active[tid] == 0:
                del active[tid]
        fs = frozenset(active)
        starts.append(a)
        ends.append(b)
        sets.append(interned.setdefault(fs, fs))
    return starts, ends, sets


def touched_transcripts(blocks, starts, ends, sets) -> set:
    out: set = set()
    n = len(starts)
    for bs, be in blocks:
        i = bisect.bisect_right(starts, bs) - 1
        if i < 0:
            i = 0
        while i < n and starts[i] < be:
            if ends[i] > bs and sets[i]:
                out |= sets[i]
            i += 1
    return out


# ---------------------------------------------------------------- per-transcript molecule counts (heavy)
PENDING_EXIT = 75    # = assembly.PENDING_EXIT: a budgeted count call stopped with contigs left
COUNT_VERSION = "2"  # bump when the per-contig counting changes: cached parts of another version are recounted
COUNT_HEADER = ["transcript_id", "gene", "chrom", "strand", "start", "end", "n_exons", "n_mol", "n_tied",
                "n_mol_primary", "tie_group", "tie_group_molecules"]


def is_candidate(secondary: bool, as_: int | None, best: int | None) -> bool:
    """The streaming assembler's admission rule (denovo_assemble.rs, `tie_ratio`): a secondary is dropped only
    when it has an AS, the molecule's genome-wide best is known and > 0, and AS < ratio x best."""
    if not secondary:
        return True
    return not (as_ is not None and best is not None and best > 0 and as_ < TIE_RATIO * best)


def _count_contig(bam, chrom, txs, best_as, tied):
    """{tid: (n_mol, n_tied, n_mol_primary)} for one contig; sweep with finalisation so memory stays local."""
    starts, ends, sets = exon_segments(txs)
    tx_end = {tid: exons[-1][1] for tid, exons in txs}
    active: dict[str, list] = {}  # tid -> [candidate names, primary names]
    heap: list = []
    out: dict = {}
    tied_names: dict = {}  # high-tie transcripts only: their tied molecules (for the read-sharing groups)
    st = {"records": 0, "secondary_candidate": 0, "secondary_dropped": 0, "secondary_unknown": 0}

    def finalize(tid):
        cand, prim = active.pop(tid)
        nt = [m for m in cand if m in tied]
        out[tid] = (len(cand), len(nt), len(prim))
        if cand and len(nt) > HIGH_TIE_MIN * len(cand):
            tied_names[tid] = sorted(nt)

    for r in bam.fetch(chrom):
        flag = r.flag
        if flag & 0x804:  # unmapped / supplementary never count
            continue
        st["records"] += 1
        name = r.query_name
        secondary = bool(flag & 0x100)
        if secondary:
            b = best_as.get(name)
            if b is None:
                st["secondary_unknown"] += 1
            a = r.get_tag("AS") if r.has_tag("AS") else None
            if not is_candidate(True, a, b):
                st["secondary_dropped"] += 1
                continue
            st["secondary_candidate"] += 1
        pos = r.reference_start
        while heap and heap[0][0] <= pos:
            finalize(heapq.heappop(heap)[1])
        for tid in touched_transcripts(r.get_blocks(), starts, ends, sets):
            slot = active.get(tid)
            if slot is None:
                slot = active[tid] = [set(), set()]
                heapq.heappush(heap, (tx_end[tid], tid))
            slot[0].add(name)
            if not secondary:
                slot[1].add(name)
    while heap:
        finalize(heapq.heappop(heap)[1])
    return out, tied_names, st


def _part_ok(part: Path, *sources) -> bool:
    """A cached per-contig part is reused when it is newer than every source and has this COUNT_VERSION."""
    if not figlib.fresh(part, *sources):
        return False
    with open(part) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    return first[:2] == ["#version", COUNT_VERSION]


def count_molecules(bam_path: str, ref_gtf: str, molecules_tsv: str, out_tsv: str, contigs: set | None = None,
                    budget_s: float = 0.0):
    """Stream `bam_path` over the contigs holding multi-exon reference transcripts; write COUNT_HEADER rows for
    EVERY multi-exon reference transcript (zeros when no molecule reaches it).

    Resumable: each contig's counts land in `<out>.parts/<contig>.tsv` first, and a re-run skips contigs whose
    part is newer than the BAM, the reference and the molecule table (a killed run loses one contig at most).
    `budget_s` > 0 stops after the first contig that ends past that many seconds (exit PENDING_EXIT; the parts
    written so far are kept, the same command continues)."""
    import time
    import pysam

    t0 = time.time()

    ref = multi_exon(assembly.ref_transcripts(ref_gtf))
    if contigs is not None:
        ref = {t: v for t, v in ref.items() if v["chrom"] in contigs}
    by_chrom: dict[str, list] = {}
    for tid, v in ref.items():
        by_chrom.setdefault(v["chrom"], []).append((tid, v["exons"]))
    parts = Path(str(out_tsv) + ".parts")
    parts.mkdir(parents=True, exist_ok=True)
    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_contigs = set(bam.references)
    todo = [c for c in sorted(by_chrom) if c in bam_contigs and not _part_ok(parts / f"{c}.tsv", bam_path, ref_gtf,
                                                                               molecules_tsv)]
    print(f"[fig3 count] {len(ref)} multi-exon reference transcripts on {len(by_chrom)} contigs; "
          f"{len(todo)} contig(s) to stream", file=sys.stderr)
    if todo:
        # a small subset of the BAM (e.g. a few chromosomes of the human BAM) needs only the names seen there; one
        # cheap pass collects them. Decided on the WHOLE evaluated contig set, not on the contigs left to do, so a
        # genome-wide count resumed over several budgeted calls never adds that extra pass (the lookups, and so the
        # counts, are the same either way: only names seen on the streamed contigs are ever looked up)
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
        keep = None
        evaluated = [c for c in by_chrom if c in bam_contigs]
        if sum(mapped.get(c, 0) for c in evaluated) < 0.5 * max(1, sum(mapped.values())):
            # EVERY name (a primary here can be tied through a secondary on another contig)
            keep = set()
            for c in todo:
                for r in bam.fetch(c):
                    if not r.flag & 0x804:
                        keep.add(r.query_name)
            print(f"[fig3 count] {len(keep)} molecules on the streamed contigs", file=sys.stderr)
        best_as, tied = load_molecules(molecules_tsv, keep)
        print(f"[fig3 count] {len(best_as)} multi-record molecules, {len(tied)} tied", file=sys.stderr)
        for i_chrom, chrom in enumerate(todo):
            if budget_s > 0 and i_chrom > 0 and time.time() - t0 > budget_s:
                print(f"[fig3 count] budget of {budget_s:.0f} s used; {len(todo) - i_chrom} contig(s) left "
                      f"({', '.join(todo[i_chrom:])}); run again to continue", file=sys.stderr)
                sys.exit(PENDING_EXIT)
            out, tnames, st = _count_contig(bam, chrom, by_chrom[chrom], best_as, tied)
            tmp = parts / f"{chrom}.tsv.tmp"
            with open(tmp, "w") as fh:
                fh.write(f"#version\t{COUNT_VERSION}\n")
                fh.write("#stats\t" + "\t".join(f"{k}={v}" for k, v in st.items()) + "\n")
                for tid, c in out.items():
                    fh.write(f"{tid}\t{c[0]}\t{c[1]}\t{c[2]}\n")
                for tid, names in tnames.items():
                    fh.write(f"@tied\t{tid}\t{','.join(names)}\n")
            tmp.replace(parts / f"{chrom}.tsv")
            print(f"[fig3 count] {chrom}: {st}", file=sys.stderr)
    counts: dict[str, tuple] = {}
    totals: dict[str, int] = {}
    tied_of: dict[str, list] = {}
    for chrom in sorted(by_chrom):
        p = parts / f"{chrom}.tsv"
        if not p.exists():
            continue
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if f[0] == "#version":
                continue
            if f[0] == "#stats":
                for kv in f[1:]:
                    k, _, v = kv.partition("=")
                    totals[k] = totals.get(k, 0) + int(v)
                continue
            if f[0] == "@tied":
                tied_of[f[1]] = f[2].split(",") if f[2] else []
                continue
            counts[f[0]] = (int(f[1]), int(f[2]), int(f[3]))
    if totals.get("secondary_unknown"):
        print(f"[fig3 count] WARNING {totals['secondary_unknown']} secondary records have no entry in "
              f"{molecules_tsv} (table from another BAM?); admitted, as the assembler does", file=sys.stderr)
    group = read_sharing_groups(tied_of)
    # distinct tied molecules behind each read-sharing group (a 40-copy group can rest on 4 molecules)
    group_mols: dict[str, set] = {}
    for tid, g in group.items():
        group_mols.setdefault(g, set()).update(tied_of.get(tid, []))
    tmp = Path(str(out_tsv) + ".tmp")
    with open(tmp, "w", newline="") as fh:
        fh.write(f"# fig3 count: bam={bam_path} ref={ref_gtf} molecules={molecules_tsv} "
                 + " ".join(f"{k}={v}" for k, v in totals.items()) + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COUNT_HEADER)
        for tid in sorted(ref, key=lambda t: (ref[t]["chrom"], ref[t]["exons"][0][0], t)):
            v = ref[tid]
            c = counts.get(tid, (0, 0, 0))
            w.writerow([tid, v["gene"], v["chrom"], v["strand"], v["exons"][0][0], v["exons"][-1][1],
                        len(v["exons"]), *c, group.get(tid, ""),
                        len(group_mols[group[tid]]) if tid in group else ""])
    tmp.replace(out_tsv)
    return out_tsv


def read_sharing_groups(tied_of: dict) -> dict:
    """{tid: group} over high-tie transcripts: transcripts that share any tied molecule are one group (union-find;
    the group is named by its smallest transcript id). Copies rebuilt from the SAME reads are not independent
    evidence, so gains are also reported in groups."""
    parent = {t: t for t in tied_of}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    first: dict[str, str] = {}
    for tid in sorted(tied_of):
        for m in tied_of[tid]:
            o = first.setdefault(m, tid)
            if o != tid:
                a, b = find(o), find(tid)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    return {t: find(t) for t in tied_of}


def _counts_current(path) -> bool:
    """The assembled counts file has this module's COUNT_HEADER (else it is re-assembled from the cached parts)."""
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    return line.rstrip("\n").split("\t") == COUNT_HEADER
    except OSError:
        pass
    return False


def read_counts(path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


# ---------------------------------------------------------------- matching, bins, summaries
def chain_key(v: dict):
    return (v["chrom"], v["strand"], v["introns"])


def matched_by_tool(ref_multi: dict, matched_ids: set) -> set:
    """Reference ids matched '=' directly or through an identical reference intron chain."""
    chains = {chain_key(ref_multi[t]) for t in matched_ids if t in ref_multi}
    return {t for t, v in ref_multi.items() if chain_key(v) in chains}


def tie_bin(frac: float) -> str:
    for label, lo, hi in TIE_BINS:
        if lo is None:
            if frac <= hi:
                return label
        elif lo < frac <= hi:
            return label
    raise ValueError(frac)


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


REF_TIE_HEADER = ["species", "transcript_id", "gene", "chrom", "strand", "start", "end", "n_exons", "n_mol",
                  "n_tied", "n_mol_primary", "tie_fraction", "tie_bin", "tie_group", "tie_group_molecules"] + [
                      f"matched_{t}" for t in TOOLS]


def species_rows(species: str, counts_tsv, ref_gtf, matched_ids: dict) -> tuple[list, dict]:
    """fig3_ref_tie rows for one species (EXPRESSED transcripts: n_mol >= 1) + bookkeeping numbers."""
    ref_multi = multi_exon(assembly.ref_transcripts(ref_gtf))
    matched = {tool: matched_by_tool(ref_multi, ids) for tool, ids in matched_ids.items()}
    rows, n_all, n_unexpr = [], 0, 0
    unexpr_matched = {t: 0 for t in TOOLS}
    for c in read_counts(counts_tsv):
        tid = c["transcript_id"]
        if tid not in ref_multi:
            continue
        n_all += 1
        n_mol, n_tied = int(c["n_mol"]), int(c["n_tied"])
        if n_mol == 0:
            # no candidate placement reaches it: outside every bin (an arm can still match it, e.g. from
            # supplementary-only or below-ratio secondary evidence; counted in the notes)
            n_unexpr += 1
            for t in TOOLS:
                unexpr_matched[t] += int(tid in matched[t])
            continue
        frac = n_tied / n_mol
        rows.append([species, tid, c["gene"], c["chrom"], c["strand"], int(c["start"]), int(c["end"]),
                     int(c["n_exons"]), n_mol, n_tied, int(c["n_mol_primary"]), round(frac, 6), tie_bin(frac),
                     c.get("tie_group", ""), c.get("tie_group_molecules", "")]
                    + [int(tid in matched[t]) for t in TOOLS])
    # the streamed-record totals from the count file's header line (records, secondaries kept / dropped)
    with open(counts_tsv) as fh:
        head = fh.readline()
    totals = {k: v for k, _, v in (kv.partition("=") for kv in head.split()[3:]) if k in
              ("records", "secondary_candidate", "secondary_dropped", "secondary_unknown")}
    i_prim, i_frac = REF_TIE_HEADER.index("n_mol_primary"), REF_TIE_HEADER.index("tie_fraction")
    no_prim = [r for r in rows if r[i_prim] == 0]
    info = {"multi_exon_ref_transcripts": n_all, "unexpressed": n_unexpr, "expressed": len(rows),
            **{f"unexpressed_matched_{t}": unexpr_matched[t] for t in TOOLS},
            "no_primary": len(no_prim), "no_primary_tie_below_1": sum(1 for r in no_prim if r[i_frac] < 1),
            **{f"streamed_{k}": v for k, v in totals.items()}}
    return rows, info


def in_category(r: dict, tie: str, primary: str) -> bool:
    """Category membership: `tie` is a TIE_BINS label, "any", or "> 0.5" (every high-tie transcript);
    `primary` a PRIMARY_STRATA key."""
    if tie == "> 0.5":
        ok = float(r["tie_fraction"]) > HIGH_TIE_MIN
    else:
        ok = tie == "any" or r["tie_bin"] == tie
    return ok and PRIMARY_STRATA[primary](int(r["n_mol_primary"]))


def cluster_of(r: dict) -> str:
    """Independent unit of a transcript: its read-sharing group (tie fraction > 0.5), else its gene."""
    if r.get("tie_group"):
        return "group:" + r["tie_group"]
    return "gene:" + r["chrom"] + ":" + (r["gene"] or r["transcript_id"])


def cluster_intervals(n_c, k_by_tool: dict, seed: int, B: int = BOOT_B) -> dict:
    """{tool: (ci_low, ci_high, method)} for the ratio sum(k)/sum(n) over clusters (n_c transcripts, k_c matched).

    95% percentile cluster bootstrap: resample clusters with replacement (the same resamples for every arm, so
    the arms stay paired). Where k = 0 or k = n every resample gives the same ratio; the interval is then
    Wilson's with n = the number of clusters (the independent units)."""
    import numpy as np

    n_c = np.asarray(n_c, dtype=float)
    G = len(n_c)
    out: dict = {}
    if G == 0:
        return {t: (None, None, "") for t in k_by_tool}
    n = n_c.sum()
    boot = {t: [] for t in k_by_tool}
    need = [t for t, k_c in k_by_tool.items() if 0 < sum(k_c) < n]
    if need:
        rng = np.random.default_rng(seed)
        ks = {t: np.asarray(k_by_tool[t], dtype=float) for t in need}
        chunk = max(1, 2_000_000 // G)
        done = 0
        while done < B:
            m = min(chunk, B - done)
            idx = rng.integers(0, G, size=(m, G))
            den = n_c[idx].sum(axis=1)
            for t in need:
                boot[t].append(ks[t][idx].sum(axis=1) / den)
            done += m
    for t, k_c in k_by_tool.items():
        k = sum(k_c)
        if t in need:
            lo, hi = np.percentile(np.concatenate(boot[t]), [2.5, 97.5])
            out[t] = (float(lo), float(hi), "cluster_bootstrap")
        else:
            lo, hi = wilson(0 if k == 0 else G, G)
            out[t] = (lo, hi, "wilson_on_clusters")
    return out


def bins_categories() -> list:
    """fig3_bins categories: every tie bin x (>= 1 primary, no primary, all), the top bin also >= 2 primary,
    and every transcript with no primary molecule. PANEL_A's are plotted."""
    cats = []
    for label, _, _ in TIE_BINS:
        cats += [(label, ">=1"), (label, "0"), (label, "any")]
        if label == TIE_BINS[-1][0]:
            cats.append((label, ">=2"))
    return cats + [("any", "0")]


def gain_categories() -> list:
    """fig3_gain categories: PANEL_A, the unstratified top bin, and every high-tie transcript (panel b's n)."""
    return PANEL_A + [(TIE_BINS[-1][0], "any"), ("> 0.5", "any")]


def _is_high(tie: str, primary: str) -> bool:
    return tie in HIGH_BINS or tie == "> 0.5" or (tie == "any" and primary == "0")


def summarize(ref_rows: list[dict], scopes: dict | None = None) -> tuple[list, list]:
    """(fig3_bins rows, fig3_gain rows) from fig3_ref_tie dict rows; `scopes` = {species: evaluation scope}."""
    bins_rows, gain_rows = [], []
    species_present = [s for s in SPECIES if any(r["species"] == s for r in ref_rows)]
    for sp in species_present:
        scope = (scopes or {}).get(sp, "")
        srows = [r for r in ref_rows if r["species"] == sp]
        bins_rows += bins_rows_for(sp, scope, srows, TOOLS)

        def m(r, t):
            return int(r[f"matched_{t}"]) == 1

        def base(r):  # matched by at least one baseline: StringTie, FLAIR or IsoSeq collapse
            return any(m(r, t) for t in BASELINES)
        for tie, prim in gain_categories():
            sub = [r for r in srows if in_category(r, tie, prim)]
            sets = {
                "rustle_not_baselines": [r for r in sub if m(r, "rustle") and not base(r)],
                "baselines_not_rustle": [r for r in sub if base(r) and not m(r, "rustle")],
                "rustle_and_baseline": [r for r in sub if m(r, "rustle") and base(r)],
                "rustle_not_primary": [r for r in sub if m(r, "rustle") and not m(r, "rustle_primary")],
                "primary_not_rustle": [r for r in sub if m(r, "rustle_primary") and not m(r, "rustle")],
            }
            high = _is_high(tie, prim)

            def groups(rs):  # independent units, as panel a counts them (read-sharing group, else gene)
                return len({cluster_of(r) for r in rs}) if high else None
            gain_rows.append([sp, scope, tie, prim, int((tie, prim) in PANEL_B), len(sub), groups(sub)]
                             + [len(v) for v in sets.values()] + [int(high)] + [groups(v) for v in sets.values()])
    return bins_rows, gain_rows


def bins_rows_for(sp: str, scope: str, srows: list[dict], tools: list[str]) -> list[list]:
    """fig3_bins rows of one species for `tools` (each row dict carries `matched_<tool>`): every category of
    bins_categories(), the cluster units and the paired cluster intervals (one bootstrap per category, seeded by the
    species and category, shared by `tools`). summarize() calls it with the annotation-free TOOLS; the guided path
    with the guided tools only (a separate table)."""
    out = []
    for tie, prim in bins_categories():
        sub = [r for r in srows if in_category(r, tie, prim)]
        groups: dict[str, list] = {}
        for r in sub:
            groups.setdefault(cluster_of(r), []).append(r)
        keys = sorted(groups)
        has_group = {bool(r.get("tie_group")) for r in sub}
        unit = ("read group" if has_group == {True} else "gene" if has_group == {False}
                else "read group or gene" if has_group else "")
        n_c = [len(groups[c]) for c in keys]
        k_by_tool = {t: [sum(int(r[f"matched_{t}"]) for r in groups[c]) for c in keys] for t in tools}
        ci = cluster_intervals(n_c, k_by_tool, zlib.crc32(f"{sp}|{tie}|{prim}".encode()))
        n = len(sub)
        plotted = int((tie, prim) in PANEL_A)
        for tool in tools:
            k = sum(k_by_tool[tool])
            gk = sum(1 for x in k_by_tool[tool] if x)
            lo, hi, method = ci[tool]
            out.append([sp, scope, tie, prim, plotted, tool, n, len(keys), unit, k, gk,
                        (k / n) if n else None, lo, hi, method])
    return out


BINS_HEADER = ["species", "scope", "tie_bin", "primary_support", "plotted", "tool", "n_ref", "n_clusters",
               "cluster_unit", "n_matched", "n_clusters_matched", "fraction", "ci_low", "ci_high", "ci_method"]
BASELINES = ["stringtie", "flair", "isoseq"]  # the trusted baselines; panels c-d set Rustle against their union
GAIN_KEYS = ["rustle_not_baselines", "baselines_not_rustle", "rustle_and_baseline",
             "rustle_not_primary", "primary_not_rustle"]
GAIN_HEADER = (["species", "scope", "tie_bin", "primary_support", "plotted", "n_ref", "n_groups"] + GAIN_KEYS
               + ["high_tie"] + [f"{k}_groups" for k in GAIN_KEYS])


# ---------------------------------------------------------------- example locus
EXAMPLE_HEADER = ["species", "track", "feature_id", "gene", "chrom", "start", "end", "value", "role"]


EXAMPLE_MIN_INTRON = 20  # bp; shorter annotated "introns" are RefSeq indel corrections of the model, not splicing


def pick_example(ref_rows: list[dict], species: str, min_intron: dict | None = None) -> dict | None:
    """The example illustrates the claim: a transcript reached by NO primary alignment that Rustle rebuilds
    exactly (an arm assembling primary alignments has no read there), whose molecules are mostly tied (tie
    fraction > 0.5: the typical case, 794 of gorilla's 839 no-primary transcripts; a low-tie one is reached only
    because the aligner's primary is not the molecule's best placement, a minority mechanism), and whose annotated
    introns are all >= EXAMPLE_MIN_INTRON bp (`min_intron`: {transcript_id: shortest intron}; a 1-2 bp gap in a
    RefSeq model "modified relative to the genomic sequence" is an indel correction, so an intron-chain match there
    shows no splicing; in gorilla no candidate has a shortest intron between 3 and 81 bp, so the cut does not
    choose among real introns). Order: most molecules, then the shortest span (a legible window), then id.
    Fallbacks, in order: drop the intron condition, then the tie condition (same order); the top tie bin, matched by
    Rustle and by none of StringTie, FLAIR, Rustle primaries-only (most molecules, fewest primary, id); then
    (0.5, 0.9]."""
    def m(r, t):
        return r[f"matched_{t}"] == "1"

    def spliced(r):
        return min_intron is None or min_intron.get(r["transcript_id"], 0) >= EXAMPLE_MIN_INTRON
    rows = [r for r in ref_rows if r["species"] == species]
    cand = [r for r in rows if int(r["n_mol_primary"]) == 0 and m(r, "rustle")]
    high = [r for r in cand if float(r["tie_fraction"]) > HIGH_TIE_MIN]
    for pool in ([r for r in high if spliced(r)], high, cand):
        if pool:
            return min(pool, key=lambda r: (-int(r["n_mol"]), int(r["end"]) - int(r["start"]), r["transcript_id"]))
    for label in reversed(HIGH_BINS):
        for strict in (True, False):
            cand = [r for r in rows if r["tie_bin"] == label and m(r, "rustle") and not m(r, "stringtie")
                    and not m(r, "flair") and (not strict or not m(r, "rustle_primary"))]
            if cand:
                return min(cand, key=lambda r: (-int(r["n_mol"]), int(r["n_mol_primary"]), r["transcript_id"]))
    return None


def gtf_transcripts_in(path, chrom: str, lo: int, hi: int) -> dict:
    """{transcript_id: [exons]} of a GTF/GFF (transcript_id attribute or GFF Parent) for the transcripts with an
    exon overlapping [lo, hi); ALL their exons are returned, so a transcript that continues past the window is
    drawn running off its edge rather than looking complete."""
    import re
    tx: dict = {}
    pat = re.compile(r'transcript_id[ =]"?([^";]+)"?')
    with assembly._open(path) as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = pat.search(f[8]) or re.search(r"Parent=([^;]+)", f[8])
            if m:
                tx.setdefault(m.group(1), []).append((int(f[3]) - 1, int(f[4])))
    return {t: sorted(v) for t, v in tx.items() if any(s < hi and e > lo for s, e in v)}


def example_rows(species: str, pick: dict, bam_path: str, molecules_tsv: str, ref_gtf, arm_gtfs: dict,
                 matched_ids: dict, bin_bp: int = 50) -> list:
    """Tracks for the example locus: read depth (primary vs tied-secondary candidate placements) in `bin_bp`
    bins, reference transcripts of the gene, and every arm's transcripts overlapping the gene span."""
    import pysam

    chrom = pick["chrom"]
    ref_all = assembly.ref_transcripts(ref_gtf)
    gene_tx = {t: v for t, v in ref_all.items() if v["gene"] == pick["gene"] and v["chrom"] == chrom}
    glo = min(v["exons"][0][0] for v in gene_tx.values())
    ghi = max(v["exons"][-1][1] for v in gene_tx.values())
    arms = {tool: gtf_transcripts_in(arm_gtfs[tool], chrom, glo, ghi) for tool in TOOLS}
    # the window covers the gene and every arm transcript that matches the example exactly
    spans = [(glo, ghi)] + [(ex[0][0], ex[-1][1]) for tool in TOOLS for q, ex in arms[tool].items()
                            if q in matched_ids.get(tool, set())]
    lo, hi = min(a for a, _ in spans), max(b for _, b in spans)
    pad = max(200, (hi - lo) // 25)
    lo, hi = max(0, lo - pad), hi + pad
    names = set()
    bam = pysam.AlignmentFile(bam_path, "rb")
    for r in bam.fetch(chrom, lo, hi):
        if not r.flag & 0x804:
            names.add(r.query_name)
    best, tied = {}, set()
    with open(molecules_tsv) as fh:
        for line in fh:
            f = line.split("\t", 4)
            if f[0] in names:
                best[f[0]] = int(f[1])
                if molecule_tied(int(f[1]), int(f[2])):
                    tied.add(f[0])
    nb = (hi - lo + bin_bp - 1) // bin_bp
    depth = {"primary_unique": [0] * nb, "primary_tied": [0] * nb, "secondary_tied": [0] * nb,
             "secondary_untied": [0] * nb}
    for r in bam.fetch(chrom, lo, hi):
        flag = r.flag
        if flag & 0x804:
            continue
        name = r.query_name
        if flag & 0x100:
            a = r.get_tag("AS") if r.has_tag("AS") else None
            if not is_candidate(True, a, best.get(name)):
                continue
            track = "secondary_tied" if name in tied else "secondary_untied"
        else:
            track = "primary_tied" if name in tied else "primary_unique"
        # a record counts once per bin: get_blocks() returns one block per =/X CIGAR operation, so a
        # mismatch-dense read has several blocks in one bin
        hit = set()
        for bs, be in r.get_blocks():
            hit.update(range(max(0, (bs - lo) // bin_bp), min(nb, (be - 1 - lo) // bin_bp + 1)))
        for b in hit:
            depth[track][b] += 1
    rows = []
    for track, vals in depth.items():
        for i, v in enumerate(vals):
            rows.append([species, f"depth_{track}", f"bin{i}", pick["gene"], chrom, lo + i * bin_bp,
                         min(hi, lo + (i + 1) * bin_bp), v, "depth"])
    for tid, v in sorted(gene_tx.items(), key=lambda kv: kv[1]["exons"][0][0]):
        role = "example" if tid == pick["transcript_id"] else "ref"
        for s, e in v["exons"]:
            rows.append([species, "ref", tid, v["gene"], chrom, s, e, "", role])
    for tool in TOOLS:
        q = arms[tool]
        for qid, exons in sorted(q.items(), key=lambda kv: (kv[1][0][0], kv[0])):
            role = "match" if qid in matched_ids.get(tool, set()) else "other"
            for s, e in exons:
                rows.append([species, tool, qid, "", chrom, s, e, "", role])
    return rows


def query_ids_matching(tmap_path, ref_ids: set) -> set:
    """Query ids with '=' to one of `ref_ids`."""
    return {r["qry_id"] for r in assembly.read_tmap(tmap_path) if r.get("class_code") == "=" and
            r.get("ref_id") in ref_ids}


# ---------------------------------------------------------------- build
def scope_label(cfg: dict, species: str) -> str:
    """Evaluation scope as printed on the figure: 'genome-wide' (every annotated contig) or the contig list."""
    import re
    if assembly.is_genome_wide(cfg, species):
        return "genome-wide"
    contigs = assembly.evaluation_contigs(cfg, species)
    m = [re.fullmatch(r"(chr)(\d+)", c) for c in contigs]
    if all(m):
        nums = sorted(int(x.group(2)) for x in m)
        if nums == list(range(nums[0], nums[-1] + 1)) and len(nums) > 1:
            return f"chr{nums[0]}-{nums[-1]}"
    return ", ".join(sorted(contigs))


def counts_path(cfg: dict, species: str) -> Path:
    """Per-transcript counts, keyed by the evaluation scope (the same key as the restricted annotation)."""
    wd = figlib.work_dir(cfg, FIG) / species / assembly.eval_dir(cfg, species).name
    wd.mkdir(parents=True, exist_ok=True)
    return wd / "ref_tie_counts.tsv"


COUNT_S_PER_M_RECORDS = (15.0, 35.0)   # gorilla genome-wide 10.7 M records in 3.6 min (2026-09-25); the 96 GB human
#                                        BAM is I/O-bound; plus ~1-2 min per call to load the human best-AS table


def plan(cfg: dict):
    """Print the heavy units `make.py data fig3` would run now, with estimated seconds; run nothing."""
    import samples
    print("sample\tunit\test_s_low\test_s_high")
    tot = [0.0, 0.0]
    print(f"# {assembly.guided_status(cfg)}")
    for key in assembly.benchmark_samples(cfg):
        for tool in TOOLS + assembly.guided_tools(cfg, key):
            if tool in assembly.RUSTLE_STAGE:
                state, reason = samples.status(cfg, key, assembly.RUSTLE_STAGE[tool])
                if state not in ("fresh", "adopt"):
                    print(f"{key}\t{assembly.RUSTLE_STAGE[tool]} ({state}: {reason}): make.py runs first\t\t")
                    continue
            what = assembly.gffcompare_state(cfg, key, tool)
            if what:
                print(f"{key}\t{tool}: {what} (shared with fig. 1)\t10\t120")
        counts = counts_path(cfg, key)
        mol = samples.product(cfg, key, "assemble", "molecules")
        ref = assembly.eval_dir(cfg, key) / "ref.gtf"
        bam = samples.get(cfg, key)["bam"]
        if figlib.fresh(counts, bam, mol, ref) and counts.exists() and _counts_current(counts):
            continue
        import pysam
        with pysam.AlignmentFile(bam, "rb") as fh:
            mapped = {st.contig: st.mapped for st in fh.get_index_statistics()}
        ev = assembly.evaluation_contigs(cfg, key)
        parts = Path(str(counts) + ".parts")
        todo = [c for c in mapped if (ev is None or c in ev) and not
                (ref.exists() and _part_ok(parts / f"{c}.tsv", bam, ref, mol))]
        n = sum(mapped[c] for c in todo) / 1e6
        lo, hi = (n * r for r in COUNT_S_PER_M_RECORDS)
        print(f"{key}\tper-transcript counts: {len(todo)} contig(s), {n:.1f} M records\t{lo:.0f}\t{hi:.0f}")
        tot[0] += lo
        tot[1] += hi
    print(f"TOTAL\t\t{tot[0]:.0f}\t{tot[1]:.0f}")


def build(cfg, data_dir, force):
    """Every sample with the three lab baselines (assembly.benchmark_samples), GENOME-WIDE (every annotated contig).

    Heavy steps, each cached: the methods and gffcompare (assembly.gffcompare; `force` re-restricts and re-runs
    gffcompare but never re-assembles), the genome-wide best-AS table (ensure_molecules), and the per-contig
    counting subprocess (`count`, resumable per contig under work_dir/fig3/<sample>/eval_<scope>/; `figs_budget_s`
    or `fig3_budget_s` bounds one call, exit 75 = run the same command again)."""
    import samples
    if assembly.plan_only(cfg):
        return plan(cfg)
    budget = assembly.Budget(cfg, "fig3")
    all_rows: list = []
    infos: dict = {}
    inputs: dict = {}
    scopes: dict = {}
    example: list = []
    example_note = "no example"
    guided_rows: list = []
    guided_inputs: dict = {}
    print(f"[fig3] {assembly.guided_status(cfg)}", file=sys.stderr)
    skipped = [k for k in assembly.guided_samples(cfg) if k not in assembly.benchmark_samples(cfg)]
    if skipped:
        print(f"[fig3] guided GTFs of {', '.join(skipped)} are not scored here: fig. 3 counts reads per transcript only "
              "on the benchmark samples (figs 1-2 score them)", file=sys.stderr)
    for species in assembly.benchmark_samples(cfg):
        row = samples.get(cfg, species)
        gcs = {}
        for tool in TOOLS:
            if force or assembly.gffcompare_state(cfg, species, tool):
                budget.check(f"gffcompare of {species} {tool}")
            gcs[tool] = assembly.gffcompare(cfg, species, tool, force=force)
        ref_gtf = gcs["rustle"]["ref"]
        mol = ensure_molecules(cfg, species, force=force, budget=budget)
        counts = counts_path(cfg, species)
        wd = counts.parent
        bam = row["bam"]
        if force or not figlib.fresh(counts, bam, mol, ref_gtf) or not _counts_current(counts):
            budget.check(f"per-transcript counts of {species}")
            rem = budget.remaining()
            rc = figlib.run([sys.executable, str(Path(__file__).resolve()), "count", "--bam", bam, "--ref",
                             str(ref_gtf), "--molecules", str(mol), "--out", str(counts)]
                            + (["--budget-s", f"{rem:.0f}"] if rem != math.inf else [])
                            + (["--force"] if force else []), log=wd / "count.log", check=False)
            if rc == PENDING_EXIT:
                raise assembly.Pending(f"fig3: per-transcript counts of {species} are partly done ({wd / 'count.log'})")
            if rc != 0:
                raise RuntimeError(f"fig3 count failed ({rc}) — see {wd / 'count.log'}")
        matched_ids = {tool: assembly.exact_matched_refs(gcs[tool]["tmap"]) for tool in TOOLS}
        rows, info = species_rows(species, counts, ref_gtf, matched_ids)
        gtools = assembly.guided_tools(cfg, species)
        if gtools:   # the separate annotation-guided path: same transcripts, categories and clusters
            for t in gtools:
                if force or assembly.gffcompare_state(cfg, species, t):
                    budget.check(f"gffcompare of {species} {t}")
                gcs[t] = assembly.gffcompare(cfg, species, t, force=force)
            ref_multi = multi_exon(assembly.ref_transcripts(ref_gtf))
            gm = {t: matched_by_tool(ref_multi, assembly.exact_matched_refs(gcs[t]["tmap"])) for t in gtools}
            i_tid = REF_TIE_HEADER.index("transcript_id")
            srows = [{**dict(zip(REF_TIE_HEADER, map(str, r))), **{f"matched_{t}": str(int(r[i_tid] in gm[t]))
                                                                    for t in gtools}} for r in rows]
            guided_rows += [b + [assembly.MODE_GUIDED, row["id"]]
                            for b in bins_rows_for(species, scope_label(cfg, species), srows, gtools)]
            inputs.update({f"{species}_tmap_{t}": gcs[t]["tmap"] for t in gtools})
            guided_inputs.update({f"{species}_{t}_gtf": samples.baseline(cfg, species, t) for t in gtools})
        all_rows += rows
        infos[species] = info
        scopes[species] = scope_label(cfg, species)
        inputs.update({f"{species}_bam": bam, f"{species}_molecules": mol, f"{species}_ref": ref_gtf,
                       f"{species}_counts": counts})
        inputs.update({f"{species}_tmap_{t}": gcs[t]["tmap"] for t in TOOLS})
        if not example and species == "gorilla":
            dict_rows = [dict(zip(REF_TIE_HEADER, map(str, r))) for r in rows]
            shortest = {t: min(b[0] - a[1] for a, b in zip(v["exons"], v["exons"][1:]))
                        for t, v in multi_exon(assembly.ref_transcripts(ref_gtf)).items()}
            pick = pick_example(dict_rows, species, shortest)
            if pick:
                qmatch = {t: query_ids_matching(gcs[t]["tmap"], {pick["transcript_id"]}) for t in TOOLS}
                example = example_rows(species, pick, bam, str(mol), ref_gtf,
                                       {t: gcs[t]["query"] for t in TOOLS}, qmatch)
                example_note = example_description(pick)
    extra = [f"{k}: genome-wide = every contig its annotation covers; not annotated, so left out for every method: "
             f"{', '.join(assembly.unannotated_contigs(cfg, k))}" for k in scopes if assembly.unannotated_contigs(cfg, k)]
    extra.append("mode: " + assembly.MODE_DENOVO_METHODS + "; " + assembly.guided_status(cfg))
    write_all(all_rows, infos, inputs, example, example_note, data_dir, notes_extra=extra, scopes=scopes)
    if guided_rows:
        figlib.write_table("fig3_guided_bins", BINS_HEADER + ["mode", "sample"], guided_rows,
                           generator="figures/fig_secondary.py build", inputs={**inputs, **guided_inputs},
                           notes=["mode: " + assembly.MODE_GUIDED + " (StringTie -G / FLAIR with the annotation, as "
                                  "supplied in samples.tsv); " + assembly.GUIDED_CAVEAT,
                                  assembly.RUSTLE_NO_GUIDED,
                                  "pre-registered: docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md",
                                  "the fig3_bins categories, transcripts, clusters and 95% cluster intervals "
                                  "(bins_rows_for), guided tools only; plotted = the panel-a categories"]
                           + [n for n in extra if not n.startswith("mode:")], data_dir=data_dir)


def example_description(pick: dict) -> str:
    return (f"example = {pick['transcript_id']} ({pick['gene']}), tie fraction {pick['tie_fraction']}, "
            f"n_mol {pick['n_mol']}, n_mol_primary {pick['n_mol_primary']}, read group {pick['tie_group']}, "
            f"rule: pick_example() (no primary, Rustle '=', tie > {HIGH_TIE_MIN}, introns >= {EXAMPLE_MIN_INTRON} bp, "
            "most molecules, shortest span)")


def write_all(all_rows, infos, inputs, example, example_note, data_dir, notes_extra, scopes=None):
    gen = "figures/fig_secondary.py build"
    notes = list(notes_extra) + [
        f"tie ratio {TIE_RATIO}; candidate placement = primary or secondary with AS >= {TIE_RATIO} x genome-wide "
        "best; counted if an aligned block (M/=/X; not N, not D) overlaps an exon; 0x800/0x4 never count",
        "matched = gffcompare '=' to this transcript or to a reference transcript with the identical intron chain",
        "scope: " + ", ".join(f"{sp}={sc}" for sp, sc in (scopes or {}).items()),
    ] + [f"{sp}: " + ", ".join(f"{k}={v}" for k, v in info.items()) for sp, info in infos.items()]
    figlib.write_table("fig3_ref_tie", REF_TIE_HEADER, all_rows, generator=gen, inputs=inputs, notes=notes,
                       data_dir=data_dir)
    dict_rows = [dict(zip(REF_TIE_HEADER, map(str, r))) for r in all_rows]
    bins_rows, gain_rows = summarize(dict_rows, scopes)
    figlib.write_table("fig3_bins", BINS_HEADER, bins_rows, generator=gen, inputs=inputs,
                       notes=notes + [
                           "primary_support: '>=1' / '0' / '>=2' molecules whose PRIMARY record reaches the "
                           "transcript ('any' = unstratified); tie_bin 'any' = every tie bin; plotted = panel a",
                           f"ci = 95% interval over clusters (cluster_unit: read-sharing group if the transcript "
                           f"has one, else gene): percentile cluster bootstrap, {BOOT_B} resamples shared by the "
                           "arms (ci_method cluster_bootstrap); Wilson with n = n_clusters where the bootstrap is "
                           "degenerate (0 or all matched; ci_method wilson_on_clusters)"],
                       data_dir=data_dir)
    figlib.write_table("fig3_gain", GAIN_HEADER, gain_rows, generator=gen, inputs=inputs,
                       notes=notes + ["plotted = panels c-d; tie_bin '> 0.5' = every high-tie transcript (their n); "
                                      "baselines = StringTie, FLAIR and IsoSeq collapse: rustle_not_baselines = matched "
                                      "by Rustle and by none of the three, baselines_not_rustle = matched by at least "
                                      "one of the three and not by Rustle",
                                      "n_groups and *_groups = independent units behind the count (high-tie rows): the "
                                      "read-sharing group, else the gene (cluster_of; the same units as fig3_bins)"],
                       data_dir=data_dir)
    figlib.write_table("fig3_example", EXAMPLE_HEADER, example, generator=gen, inputs=inputs,
                       notes=notes_extra + [example_note], data_dir=data_dir)


# ---------------------------------------------------------------- plot
def _f(x):
    return float(x) if x not in ("", None) else float("nan")


def _contigs_label(scope: str) -> str:
    """'chr20,chr21,chr22' -> 'chr20-22' (fig. 1 writes the contig list)."""
    import re
    parts = [p.strip() for p in scope.split(",") if p.strip()]
    m = [re.fullmatch(r"chr(\d+)", p) for p in parts]
    if len(parts) > 1 and all(m):
        nums = sorted(int(x.group(1)) for x in m)
        if nums == list(range(nums[0], nums[-1] + 1)):
            return f"chr{nums[0]}-{nums[-1]}"
    return scope


def precision_cost(data_dir) -> list[tuple]:
    """[(species, scope, precision rustle, precision rustle_primary, multi-exon queries rustle, multi-exon queries
    primary, matching chains rustle, matching chains primary)] read from fig. 1's table (gffcompare intron-chain
    precision = matching chains / MULTI-EXON query transcripts); [] when that table is absent."""
    try:
        gc = figlib.read_table("fig1_gffcompare", data_dir)
    except FileNotFoundError:
        return []
    out = []
    for sp in SPECIES:
        rr = {r["tool"]: r for r in gc if r["species"] == sp and r["level"] == "intron_chain"}
        if "rustle" in rr and "rustle_primary" in rr:
            a, p = rr["rustle"], rr["rustle_primary"]
            out.append((sp, _contigs_label(a["scope"]), _f(a["pr"]), _f(p["pr"]), int(a["n_query_multiexon"]),
                        int(p["n_query_multiexon"]), int(a["matching"]), int(p["matching"])))
    return out


LETTER_X, TITLE_X = 0.004, 0.028  # panel letters and titles, figure coordinates
XSPAN = 2.75  # panel b: x half-range in units of the largest bar (room for the "N in G groups" labels)


def _fig_dy(fig, inches: float) -> float:
    return inches / fig.get_figheight()


# layout in inches (absolute, so adding the human rows never squeezes the gorilla ones)
A_ABOVE, A_FACETS, A_BELOW = 0.74, 0.85, 0.5      # panel a: titles/headers, facets, ticks + x label
B_ABOVE, B_BARS, B_BELOW = 0.44, 0.98, 0.6        # panels c-d (one per species, stacked): title, bars, x label
NOTE_H = 0.5                                      # precision-cost note under the b stack
C_ABOVE, C_BODY_MIN, C_BELOW = 0.44, 2.0, 0.42    # example: title, depth + tracks, x label
X_LEFT, X_RIGHT = 0.155, 0.985                    # figure fractions of the plotting area
B_RIGHT, C_LEFT = 0.415, 0.595                    # panel b's right edge, the example's left edge


def plot(data_dir, out_dir):
    import matplotlib.pyplot as plt

    bins = figlib.read_table("fig3_bins", data_dir)
    gain = figlib.read_table("fig3_gain", data_dir)
    try:
        ex = figlib.read_table("fig3_example", data_dir)
    except FileNotFoundError:
        ex = []
    species = [s for s in SPECIES if any(r["species"] == s for r in bins)]
    nrow = len(species)
    b_stack = nrow * (B_ABOVE + B_BARS + B_BELOW) + NOTE_H
    bottom_h = max(b_stack, C_ABOVE + C_BODY_MIN + C_BELOW)
    height = 0.05 + nrow * (A_ABOVE + A_FACETS + A_BELOW) + bottom_h + 0.05
    # savefig(bbox_inches="tight") adds 0.1 in of padding on each side: the canvas is narrower by that, so the
    # saved page stays within the 183 mm double column
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.22, height))

    def region(top_in, bottom_in, left=X_LEFT, right=X_RIGHT):
        return fig.add_gridspec(1, 1, left=left, right=right, top=1 - top_in / height,
                                bottom=1 - bottom_in / height)[0]

    # one x range for every species' panel a (so the facets compare across species at a glance)
    vals = [_f(r[k]) for r in bins if r["plotted"] == "1" for k in ("ci_high", "fraction")]
    xmax = min(1.0, math.ceil((max([v for v in vals if v == v] + [0.1]) + 0.12) * 10) / 10)
    letters = iter("abcdefgh")
    y = 0.05
    for sp in species:
        _panel_fraction(fig, region(y + A_ABOVE, y + A_ABOVE + A_FACETS), [r for r in bins if r["species"] == sp],
                        sp, next(letters), xmax)
        y += A_ABOVE + A_FACETS + A_BELOW
    y_bottom = y
    for sp in species:
        ax = fig.add_subplot(region(y + B_ABOVE, y + B_ABOVE + B_BARS, right=B_RIGHT))
        _panel_gain(fig, ax, [r for r in gain if r["species"] == sp], sp, next(letters))
        y += B_ABOVE + B_BARS + B_BELOW
    _precision_note(fig, 1 - (y + 0.02) / height, precision_cost(data_dir))
    if ex:
        _panel_example(fig, region(y_bottom + C_ABOVE, y_bottom + bottom_h - C_BELOW, left=C_LEFT), ex,
                       next(letters))
    # the comparison mode, stated on the figure (every panel above compares annotation-free methods only)
    fig.text(0.01, 1 + 0.06 / height, figlib.mode_lines("fig3g_guided", (Path(data_dir) / "fig3_guided_bins.tsv")
                                                        .exists()),
             fontsize=6.0, color=figlib.INK_2, ha="left", va="bottom", linespacing=1.3)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, f"{FIG}_secondary", out_dir)
    plt.close(fig)
    return paths + plot_guided(data_dir, out_dir)


def plot_guided(data_dir, out_dir) -> list:
    """fig3g_guided: panel a for the annotation-guided StringTie/FLAIR runs only (fig3_guided_bins), one row of facets
    per sample; no annotation-free method and no Rustle row. Drawn only when the table exists."""
    import matplotlib.pyplot as plt

    try:
        rows = figlib.read_table("fig3_guided_bins", data_dir)
    except FileNotFoundError:
        return []
    tools = [t for t in figlib.GUIDED_TOOL_ORDER if any(r["tool"] == t for r in rows)]
    species = list(dict.fromkeys(r["species"] for r in rows))
    height = 0.30 + len(species) * (A_ABOVE + A_FACETS * len(tools) / len(TOOLS) + 0.3 + A_BELOW)
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.22, height))
    vals = [_f(r[k]) for r in rows if r["plotted"] == "1" for k in ("ci_high", "fraction")]
    xmax = min(1.0, math.ceil((max([v for v in vals if v == v] + [0.1]) + 0.12) * 10) / 10)
    y = 0.30
    for sp, letter in zip(species, "abcdefgh"):
        h = A_FACETS * len(tools) / len(TOOLS) + 0.3
        spec = fig.add_gridspec(1, 1, left=X_LEFT, right=X_RIGHT, top=1 - (y + A_ABOVE) / height,
                                bottom=1 - (y + A_ABOVE + h) / height)[0]
        _panel_fraction(fig, spec, [r for r in rows if r["species"] == sp], sp, letter, xmax, tools=tools)
        y += A_ABOVE + h + A_BELOW
    fig.text(0.01, 1 - 0.04 / height, "Annotation-guided runs only (never compared with the annotation-free methods "
             "of Fig. 3). The guided tools were given the annotation they are scored against.\nRustle has no "
             "annotation-guided transcript assembly, so it has no row here.", fontsize=5.8, color=figlib.INK_2,
             ha="left", va="top")
    figlib.stamp_provisional(fig, ["fig3_guided_bins"], data_dir)
    paths = figlib.save(fig, f"{FIG}g_guided", out_dir)
    plt.close(fig)
    return paths


def _category_label(tie: str, prim: str) -> str:
    if prim == "0":
        return "any tie fraction" if tie == "any" else f"tie fraction {tie}"
    return f"tie fraction {tie}"


def _facet_label(tie: str, prim: str) -> str:
    """Facet title (the group header above names the tie fraction)."""
    return "all" if tie == "any" else tie


def _panel_fraction(fig, spec, rows, sp, letter, xmax, tools=None):
    """Dot-and-interval facets, one per panel-a category: the tie bins among transcripts reached by >= 1 primary
    alignment, then the transcripts reached by none. Arms are rows, named once on the left (axis labels). `tools`
    (default: the annotation-free TOOLS) are the rows; the first one carries the category's n."""
    from matplotlib import gridspec
    from matplotlib.lines import Line2D

    TOOLS = tools or globals()["TOOLS"]   # noqa: N806 (the rows of this panel)
    head_tool = TOOLS[0]
    by = {(r["tie_bin"], r["primary_support"], r["tool"]): r for r in rows}
    scope = next((r["scope"] for r in rows if r.get("scope")), "")
    ncat = len(PANEL_A)
    g = gridspec.GridSpecFromSubplotSpec(1, ncat + 1, subplot_spec=spec,
                                         width_ratios=[1] * (ncat - 1) + [0.28, 1.0], wspace=0.13)
    cols = list(range(ncat - 1)) + [ncat]
    ypos = {t: len(TOOLS) - 1 - i for i, t in enumerate(TOOLS)}
    n_expr = sum(int(by[(c, p, head_tool)]["n_ref"]) for c, p in PANEL_A if (c, p, head_tool) in by)
    axes = []
    for j, (tie, prim) in enumerate(PANEL_A):
        ax = fig.add_subplot(g[cols[j]], sharey=axes[0] if axes else None)
        axes.append(ax)
        for i in range(len(TOOLS)):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color="#f4f3ef", zorder=0, linewidth=0)
        head = by.get((tie, prim, head_tool))
        n = int(head["n_ref"]) if head else 0
        high = _is_high(tie, prim)
        for tool in TOOLS:
            r = by.get((tie, prim, tool))
            if not r or not int(r["n_ref"]):
                continue
            y = ypos[tool]
            frac, lo, hi, k = _f(r["fraction"]), _f(r["ci_low"]), _f(r["ci_high"]), int(r["n_matched"])
            # only the primaries-only Rustle arm cannot reach these by construction (StringTie reads secondaries: observed 0)
            structural = prim == "0" and tool == "rustle_primary" and k == 0
            if not structural and lo == lo:
                ax.hlines(y, lo, hi, color=figlib.TOOL_COLOR[tool], linewidth=0.9, zorder=2)
            ax.plot([frac], [y], markersize=3.8, zorder=3, clip_on=False, **figlib.tool_marker_kwargs(tool))
            if high:  # few, clustered transcripts: print the count and the read groups behind it
                txt = f"{k} ({r['n_clusters_matched']})" if k else "0"
                end = hi if (not structural and hi == hi) else frac
                ax.annotate(txt, (end, y), xytext=(4.5, 0), textcoords="offset points", fontsize=5.5,
                            va="center", ha="left", color=figlib.INK_2, zorder=4, annotation_clip=False)
        unit = (head or {}).get("cluster_unit", "")
        n_cl = int(head["n_clusters"]) if head else 0
        unit_txt = {"read group": "groups", "gene": "genes",
                    "read group or gene": "groups or genes"}.get(unit, "units")
        ax.set_title(f"{_facet_label(tie, prim)}\nn = {n:,}\n{n_cl:,} {unit_txt}", fontsize=6, pad=3,
                     color=figlib.INK, linespacing=1.15)
        ax.set_xlim(0, xmax)
        ticks = [t for t in (0, 0.5, 1.0) if t <= xmax + 1e-9]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=6)
        ax.set_ylim(-0.6, len(TOOLS) - 0.4)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        if j == 0:
            ax.set_yticks([ypos[t] for t in TOOLS])
            ax.set_yticklabels([figlib.TOOL_LABEL[t] for t in TOOLS], fontsize=6.5)
        else:
            ax.tick_params(axis="y", labelleft=False)
    # group headers (with rules), panel letter + title, shared x label; placed from the facet positions
    p0, p4, p5 = axes[0].get_position(), axes[ncat - 2].get_position(), axes[ncat - 1].get_position()
    rule_y = p0.y1 + _fig_dy(fig, 0.36)
    for x0, x1 in ((p0.x0, p4.x1), (p5.x0, p5.x1)):
        fig.add_artist(Line2D([x0, x1], [rule_y, rule_y], transform=fig.transFigure, color=figlib.INK_3,
                              linewidth=0.6))
    kw = dict(fontsize=6.5, color=figlib.INK, ha="center", va="bottom")
    fig.text((p0.x0 + p4.x1) / 2, rule_y + _fig_dy(fig, 0.02),
             "Tie fraction: share of reads with a second alignment ≥ 98% of their best score "
             "(transcripts reached by ≥ 1 primary alignment)", **kw)
    fig.text((p5.x0 + p5.x1) / 2, rule_y + _fig_dy(fig, 0.02), "Reached by no primary\nalignment",
             linespacing=1.1, **kw)
    title_y = rule_y + _fig_dy(fig, 0.27)
    fig.text(LETTER_X, title_y, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(TITLE_X, title_y, f"{figlib.SPECIES_LABEL.get(sp, sp)}, {scope}: {n_expr:,} multi-exon RefSeq "
             "transcripts with ≥ 1 read; exact intron-chain match by method", fontsize=7.5, va="bottom", ha="left")
    fig.text((p0.x0 + p5.x1) / 2, p0.y0 - _fig_dy(fig, 0.2),
             "Fraction of the category's transcripts whose intron chain the method reproduces exactly (gffcompare "
             "'=')\nbars: 95% interval, resampling read-sharing groups (groups) or genes; printed: transcripts "
             "matched (groups or genes with a match)",
             fontsize=6.5, ha="center", va="top", color=figlib.INK, linespacing=1.2)


def _panel_gain(fig, ax, rows, sp, letter):
    """Diverging bars per high-tie category: Rustle-only matches (right) vs the other side's only (left)."""
    from matplotlib.transforms import blended_transform_factory

    by = {(r["tie_bin"], r["primary_support"]): r for r in rows}
    plotted = [by[c] for c in PANEL_B if c in by]
    kw_r = figlib.tool_bar_kwargs("rustle")
    kw_p = figlib.tool_bar_kwargs("rustle_primary")
    kw_o = {"color": figlib.INK_3, "edgecolor": figlib.SURFACE, "linewidth": 0.8}
    # the other side: any of the three trusted baselines (their union), then the unseeded Rustle arm
    # (the union is spelled out under the axis: Arial has no set-union glyph)
    comps = [("rustle_not_baselines", "baselines_not_rustle", "Any baseline", kw_o),
             ("rustle_not_primary", "primary_not_rustle", figlib.TOOL_LABEL["rustle_primary"], kw_p)]
    vals = [int(r[k]) for r in plotted for a, b, _, _ in comps for k in (a, b)]
    m = max([1] + vals)
    h = 0.62

    def lab(n, grp, unit):
        if n == 0 or grp in ("", None):
            return f"{n}"
        return f"{n} in {grp} {unit}" + ("" if grp == "1" else "s")

    yt, yl = [], []
    y = 0.0
    for r in reversed(plotted):  # y grows upwards: panel a's order reads top to bottom
        # units: read-sharing groups in the tie bins; read groups or genes where no primary reaches (as in panel a)
        unit = "unit" if r["primary_support"] == "0" else "group"
        for right, left, name, kw in reversed(comps):
            a, b = int(r[right]), int(r[left])
            ax.barh(y, a, height=h, **kw_r)
            ax.barh(y, -b, height=h, **kw)
            ax.text(a + m * 0.06, y, lab(a, r.get(f"{right}_groups"), unit), va="center", ha="left", fontsize=6,
                    color=figlib.INK)
            ax.text(-b - m * 0.06, y, lab(b, r.get(f"{left}_groups"), unit), va="center", ha="right", fontsize=6,
                    color=figlib.INK)
            yt.append(y)
            yl.append(name)
            y += 1
        head = (("reached by no primary alignment" if r["primary_support"] == "0"
                 else f"{_category_label(r['tie_bin'], r['primary_support'])}, ≥ 1 primary alignment")
                + f"  (n = {int(r['n_ref']):,}, {r['n_groups']} "
                + ("groups or genes)" if unit == "unit" else "groups)"))
        # block header across the label column and the bars (x in figure coordinates, aligned with the title)
        ax.text(TITLE_X, y - 0.05, head, ha="left", va="center", fontsize=6, color=figlib.INK, fontweight="bold",
                zorder=5, transform=blended_transform_factory(fig.transFigure, ax.transData), clip_on=False,
                bbox=dict(facecolor=figlib.SURFACE, edgecolor="none", pad=0.6))
        y += 1.0
    ax.axvline(0, color=figlib.INK_2, linewidth=0.6)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=6.5)
    ax.set_ylim(-0.6, y - 0.45)
    ax.set_xlim(-XSPAN * m, XSPAN * m)
    step = next(st for st in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10**4, 10**5) if st >= m / 2.5)
    ticks = [step * i for i in range(-int(m // step), int(m // step) + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(int(round(t)))}" for t in ticks])
    ax.set_xlim(-XSPAN * m, XSPAN * m)
    ax.set_xlabel("Transcripts matched by one side only\n← other side            Rustle →\n"
                  "any baseline = StringTie, FLAIR or IsoSeq collapse", fontsize=6.5, linespacing=1.15)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    allh = by.get(("> 0.5", "any"))
    scope = next((r["scope"] for r in rows if r.get("scope")), "")
    sub = (f"{int(allh['n_ref']):,} transcripts with tie fraction > 0.5 in {allh['n_groups']} read-sharing groups"
           if allh else "")
    pos = ax.get_position()
    ty = pos.y1 + _fig_dy(fig, 0.1)
    fig.text(LETTER_X, ty, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(TITLE_X, ty, f"{figlib.SPECIES_LABEL.get(sp, sp)}, {scope}\n{sub}", fontsize=7, va="bottom",
             ha="left", linespacing=1.2)


def _precision_note(fig, y, cost):
    """The precision cost of seeding, read from fig. 1's table (never hard-coded); `y` = top, figure fraction."""
    if not cost:
        return
    lines = ["Cost of seeding (Fig. 1): gffcompare intron-chain precision, Rustle vs Rustle with",
             "primary alignments only, and what the extra multi-exon transcripts add:"]
    for sp, scope, pr, pp, nq, npq, mc, mpc in cost:
        dq, dm = nq - npq, mc - mpc
        lines.append(f"  {sp.capitalize()}, {scope}: {pr:.1%} vs {pp:.1%}; +{dq:,} transcripts, "
                     f"+{dm:,} matching chains ({dm / dq:.1%})" if dq > 0 else
                     f"  {sp.capitalize()}, {scope}: {pr:.1%} vs {pp:.1%}")
    fig.text(TITLE_X, y, "\n".join(lines), fontsize=6, va="top", ha="left", color=figlib.INK_2, linespacing=1.3)


def _panel_example(fig, spec, ex, letter):
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mticker
    import numpy as np
    from matplotlib import gridspec

    depth = [r for r in ex if r["track"].startswith("depth_")]
    feats = [r for r in ex if not r["track"].startswith("depth_")]
    ref_rows = [r for r in feats if r["track"] == "ref"]
    pick = next((r for r in ref_rows if r["role"] == "example"), ref_rows[0])
    tracks = ["ref"] + TOOLS
    by_track: dict = {t: {} for t in tracks}
    for r in feats:
        by_track[r["track"]].setdefault(r["feature_id"], []).append((int(r["start"]), int(r["end"]), r["role"]))
    lo = min(int(r["start"]) for r in depth)
    hi = max(int(r["end"]) for r in depth)
    g = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=spec, height_ratios=[0.32, 0.32, 1.95], hspace=0.14)
    axp = fig.add_subplot(g[0])
    axs = fig.add_subplot(g[1], sharex=axp, sharey=axp)
    axt = fig.add_subplot(g[2], sharex=axp)
    # read depth in two tracks named on their axes: primary alignments (every arm can use them) and the
    # candidate secondaries within 0.98 x best AS (only the seeded Rustle arm uses them)
    groups = [(axp, ("depth_primary_unique", "depth_primary_tied"), "Primary", figlib.INK_3),
              (axs, ("depth_secondary_tied", "depth_secondary_untied"), "Secondary\n≥ 98% of best",
               figlib.BLUE[650])]
    peak = 1.0
    for ax, keys, name, col in groups:
        acc: dict = {}
        for r in depth:
            if r["track"] in keys:
                acc[int(r["start"])] = acc.get(int(r["start"]), 0.0) + float(r["value"]) 
        xs = sorted(acc)
        v = np.array([acc[x] for x in xs])
        step = (xs[1] - xs[0]) if len(xs) > 1 else 50
        ax.fill_between(np.array(xs) + step / 2, 0, v, step="mid", color=col, linewidth=0)
        peak = max(peak, float(v.max()) if len(v) else 0.0)
        ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=6.5, labelpad=4)
        ax.tick_params(axis="x", labelbottom=False, length=0)
        ax.grid(axis="y", visible=False)
        ax.spines["bottom"].set_color(figlib.GRID)
    axp.set_ylim(0, peak * 1.1)
    axp.set_yticks([0, int(peak)])
    axp.set_yticklabels(["0", f"{int(peak)}"], fontsize=5.5)
    axs.tick_params(axis="y", labelsize=5.5)
    pos = axp.get_position()
    ty = pos.y1 + _fig_dy(fig, 0.1)
    lx = pos.x0 - 1.12 / fig.get_figwidth()  # left of the track labels
    fig.text(lx, ty, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(lx + TITLE_X - LETTER_X, ty,
             f"{pick['species'].capitalize()}, {pick['chrom']}: {pick['gene'].removeprefix('gene-')}, reference "
             f"{pick['feature_id'].removeprefix('rna-')} (star)\nDepth per 50 bp; = its exact intron chain; "
             "faded = another chain", fontsize=7, va="bottom", ha="left", linespacing=1.2)
    ypos = 0.0
    yticks, ylabels = [], []
    band = 0
    for t in tracks:
        items = by_track[t]
        label, color = ("RefSeq", figlib.INK_2) if t == "ref" else (figlib.TOOL_LABEL[t], figlib.TOOL_COLOR[t])
        row_ids = sorted(items, key=lambda k: (0 if any(role in ("match", "example") for _, _, role in items[k])
                                               else 1, min(s for s, _, _ in items[k]), k))
        shown = row_ids[:3]
        top = ypos
        n_rows = max(1, len(shown))
        # one shaded band per track: its transcripts are the rows inside the band, labelled at the first row
        n_more = 0.45 if len(row_ids) > len(shown) else 0.0
        if band % 2 == 0:
            axt.axhspan(top - 0.85 * (n_rows - 1) - 0.45 - n_more, top + 0.45, color="#f4f3ef", zorder=0,
                        linewidth=0)
        band += 1
        yticks.append(top)
        ylabels.append(label)
        if not shown:
            axt.text(lo + (hi - lo) * 0.01, ypos, "no transcript in this window", fontsize=5.5,
                     color=figlib.INK_3, va="center")
        for k in shown:
            exs = sorted((s, e) for s, e, _ in items[k])
            hl = any(role in ("match", "example") for _, _, role in items[k])
            axt.plot([exs[0][0], exs[-1][1]], [ypos, ypos], color=color, linewidth=0.5, zorder=1,
                     alpha=1.0 if hl else 0.5)
            for s, e in exs:
                kw = dict(facecolor=color, alpha=1.0 if hl else 0.45, linewidth=0)
                if t == "rustle_primary":
                    kw.update(hatch=figlib.TOOL_HATCH["rustle_primary"], edgecolor=figlib.SURFACE)
                axt.add_patch(mpatches.Rectangle((s, ypos - 0.3), e - s, 0.6, zorder=2, **kw))
            if hl and t == "ref":  # a drawn star (a path: Arial has no star glyph)
                axt.plot([min(exs[-1][1], hi) + (hi - lo) * 0.016], [ypos], marker="*", markersize=5,
                         color=figlib.INK, markeredgewidth=0, linestyle="none", clip_on=False, zorder=4)
            elif hl:
                axt.text(min(exs[-1][1], hi) + (hi - lo) * 0.006, ypos, "=", fontsize=6,
                         va="center", ha="left", color=figlib.INK)
            ypos -= 0.85
        extra = len(row_ids) - len(shown)
        if extra > 0:  # its own row, so it never sits on a transcript
            axt.text(lo + (hi - lo) * 0.01, ypos + 0.2, f"+{extra} more not drawn", fontsize=5.5,
                     color=figlib.INK_3, va="center", ha="left")
            ypos -= 0.45
        if not shown:
            ypos -= 0.85
        ypos -= 0.25
    axt.set_ylim(ypos + 0.3, 0.55)
    axt.set_xlim(lo, hi)
    axt.set_yticks(yticks)
    axt.set_yticklabels(ylabels, fontsize=6.5)
    axt.tick_params(axis="y", length=0)
    axt.grid(False)
    axt.spines["left"].set_visible(False)
    axt.set_xlabel(f"{pick['chrom']} position (kb)", fontsize=6.5)
    axt.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v / 1000:,.1f}"))
    axt.set_xticks([t for t in mticker.MaxNLocator(5).tick_values(lo, hi) if lo <= t <= hi])


# ---------------------------------------------------------------- caption numbers (light: tables only)
def caption_numbers(data_dir=None) -> str:
    """Every number the caption quotes, read from the tables (`python3 figures/fig_secondary.py summary`), so the
    caption can be refreshed after `make.py data fig3` without recomputing anything by hand."""
    data_dir = data_dir or figlib.DATA_DIR
    bins = figlib.read_table("fig3_bins", data_dir)
    gain = figlib.read_table("fig3_gain", data_dir)
    ref = figlib.read_table("fig3_ref_tie", data_dir)
    out = []
    for note in figlib.table_meta("fig3_ref_tie", data_dir).get("note", []):
        if note.split(":")[0] in SPECIES or note.startswith("scope") or note.startswith("provisional"):
            out.append(f"note  {note}")
    out.append(f"note  {figlib.table_meta('fig3_example', data_dir).get('note', ['?'])[-1]}")
    out += without_contig(ref, "gorilla", "NC_073244.2")
    try:
        gb = figlib.read_table("fig3_guided_bins", data_dir)
    except FileNotFoundError:
        gb = []
        out.append("guided: not available (guided StringTie/FLAIR GTFs not supplied)")
    for r in gb:
        if r["plotted"] == "1":
            out.append(f"guided {r['species']} {r['tool']:>17} tie {r['tie_bin']:>10} primary {r['primary_support']:>3}: "
                       f"{r['n_matched']}/{r['n_ref']} ({r['n_clusters_matched']} units) [{_f(r['ci_low']):.2f},"
                       f"{_f(r['ci_high']):.2f}] (annotation-guided; no Rustle counterpart)")
    for sp in SPECIES:
        b = {(r["tie_bin"], r["primary_support"], r["tool"]): r for r in bins if r["species"] == sp}
        if not b:
            continue
        out.append(f"== {sp} ({next(r['scope'] for r in bins if r['species'] == sp)})")
        for tie, prim in PANEL_A + [("> 0.9", "any"), ("> 0.9", ">=2")]:
            h = b.get((tie, prim, "rustle"))
            if not h or not int(h["n_ref"]):
                continue
            arms = "  ".join(f"{t}={b[(tie, prim, t)]['n_matched']}({b[(tie, prim, t)]['n_clusters_matched']}g)"
                             f"[{_f(b[(tie, prim, t)]['ci_low']):.2f},{_f(b[(tie, prim, t)]['ci_high']):.2f}]"
                             for t in TOOLS)
            out.append(f"bins  tie {tie:>10} primary {prim:>3}: n={h['n_ref']} clusters={h['n_clusters']} "
                       f"({h['cluster_unit']})  {arms}")
        for r in gain:
            if r["species"] != sp:
                continue
            g = {k: (r[k], r.get(f"{k}_groups", "")) for k in GAIN_KEYS}
            out.append(f"gain  tie {r['tie_bin']:>10} primary {r['primary_support']:>3}: n={r['n_ref']} "
                       f"groups={r['n_groups']}  " + "  ".join(f"{k}={v[0]}" + (f"/{v[1]}g" if v[1] else "")
                                                               for k, v in g.items()))
        # read groups behind the seeding gains (Rustle, not Rustle primaries-only; tie > 0.5)
        rows = [r for r in ref if r["species"] == sp and float(r["tie_fraction"]) > HIGH_TIE_MIN]
        gains = [r for r in rows if r["matched_rustle"] == "1" and r["matched_rustle_primary"] == "0"]
        by_g: dict = {}
        for r in gains:
            by_g.setdefault(r["tie_group"], []).append(r)
        for grp, rs in sorted(by_g.items(), key=lambda kv: -len(kv[1])):
            members = [r for r in rows if r["tie_group"] == grp]
            genes = sorted({r["gene"].removeprefix("gene-") for r in members})
            out.append(f"group {grp}: {len(rs)} gains ({sum(1 for r in rs if r['n_mol_primary'] == '0')} with no "
                       f"primary) of {len(members)} high-tie transcripts, {len(genes)} genes {genes[0]}..{genes[-1]}, "
                       f"{members[0]['chrom']}:{min(int(r['start']) for r in members):,}-"
                       f"{max(int(r['end']) for r in members):,}, n_mol {sorted({int(r['n_mol']) for r in members})}, "
                       f"distinct tied molecules in the group {members[0].get('tie_group_molecules', '?')}")
        same = all(b[(tie, prim, "rustle")]["n_matched"] == b[(tie, prim, "rustle_primary")]["n_matched"]
                   for tie, prim in [("0", ">=1")] if (tie, prim, "rustle") in b)
        g0 = next((r for r in gain if r["species"] == sp and r["tie_bin"] == "0" and r["primary_support"] == ">=1"),
                  None)
        if g0:
            out.append(f"control tie 0: rustle {b[('0', '>=1', 'rustle')]['n_matched']} vs rustle_primary "
                       f"{b[('0', '>=1', 'rustle_primary')]['n_matched']} (equal counts {same}; rustle_not_primary "
                       f"{g0['rustle_not_primary']}, primary_not_rustle {g0['primary_not_rustle']})")
    for sp, scope, pr, pp, nq, npq, mc, mpc in precision_cost(data_dir):
        out.append(f"fig1  {sp} {scope}: intron-chain precision rustle {pr} vs rustle_primary {pp} "
                   f"(multi-exon query transcripts {nq} vs {npq}; matching chains {mc} vs {mpc}; marginal "
                   f"+{mc - mpc} of +{nq - npq})")
    return "\n".join(out)


def without_contig(ref: list[dict], species: str, contig: str) -> list[str]:
    """The gorilla evidence without one contig (NC_073244.2: the contig the seeding default was decided on, and the
    one fig. 6d uses): the no-primary category, the seeding gains, and the top tie bin, from fig3_ref_tie."""
    rows = [r for r in ref if r["species"] == species and r["chrom"] != contig]
    if not rows:
        return []

    def m(r, t):
        return r[f"matched_{t}"] == "1"
    nop = [r for r in rows if r["n_mol_primary"] == "0"]
    k_nop = sum(m(r, "rustle") for r in nop)
    others = {t: sum(m(r, t) for r in nop) for t in TOOLS if t != "rustle"}
    high = [r for r in rows if float(r["tie_fraction"]) > HIGH_TIE_MIN]
    gains = [r for r in high if m(r, "rustle") and not m(r, "rustle_primary")]
    lost = [r for r in high if m(r, "rustle_primary") and not m(r, "rustle")]
    top = [r for r in rows if r["tie_bin"] == TIE_BINS[-1][0] and int(r["n_mol_primary"]) >= 1]
    top_k = {t: (sum(m(r, t) for r in top), len({cluster_of(r) for r in top if m(r, t)})) for t in TOOLS}
    return [f"without {contig} ({species}): no primary {k_nop}/{len(nop)} = {k_nop / max(1, len(nop)):.3f} "
            f"(units {len({cluster_of(r) for r in nop if m(r, 'rustle')})}; other arms {others}); seeding gains "
            f"(tie > 0.5) {len(gains)} in {len({cluster_of(r) for r in gains})} read groups on "
            f"{len({r['chrom'] for r in gains})} contigs, {len(lost)} lost; top bin (>= 1 primary) n={len(top)}: "
            + "  ".join(f"{t}={k}({g}g)" for t, (k, g) in top_k.items())]


# ---------------------------------------------------------------- CLI: the heavy counting step
def main(argv=None):
    ap = argparse.ArgumentParser(description="fig3: per-transcript molecule counts (heavy BAM pass); caption numbers")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("count")
    c.add_argument("--bam", required=True)
    c.add_argument("--ref", required=True)
    c.add_argument("--molecules", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--contigs", default="", help="comma list (default: every contig with a reference transcript)")
    c.add_argument("--force", action="store_true", help="discard cached per-contig parts")
    c.add_argument("--budget-s", type=float, default=0.0,
                   help="stop (exit 75) after the first contig that ends past this many seconds; 0 = no limit")
    sm = sub.add_parser("summary", help="print every number the caption quotes, from figures/data (light)")
    sm.add_argument("--data", default=str(figlib.DATA_DIR))
    a = ap.parse_args(argv)
    if a.cmd == "summary":
        print(caption_numbers(Path(a.data)))
        return
    if a.cmd == "count":
        contigs = {x for x in a.contigs.split(",") if x} or None
        if a.force:
            import shutil
            shutil.rmtree(str(a.out) + ".parts", ignore_errors=True)
        count_molecules(a.bam, a.ref, a.molecules, a.out, contigs, budget_s=a.budget_s)


if __name__ == "__main__":
    main()
