"""Figure 5 — copy assignment by copy identity on simulated reads (a-b, one small multiple per sample) and the
hard-locus benchmark on real reads against the lab's other methods (c).

a, b  Per sample, over the simulated reads whose primary alignment has MAPQ 0, per identity band of the read's
      source copy: the fraction assigned and the fraction correct among assigned (95% interval with the SOURCE COPY
      as the unit: reads from one copy are not independent; no interval when fewer than 5 copies are assigned) of
      four readings — the copy-assignment test scored within the read's source family (known only in simulation),
      the default output read as a user would (any assigned result), the union test (`--union-certificate`), and the
      aligner's primary alignment (the pre-registered baseline). Behind the fraction-assigned points, the share of
      the band's MAPQ-0 reads with >= 1 decisive site in their source family (what the test can use). Data:
      figures/_o2.py (same runs as Figure 4). Genome-wide on every sample (docs/PREREG_genome_wide_copy_assignment_
      2026-09-25.md, experiment A), or the development tables (cfg `o2_scope dev`).
c     The hard-locus benchmark: real reads, the lab's StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse on the same
      BAMs, and Rustle's arm = ONE `copy_assign --families --gtf` run (its own transcripts AND the hard set), not the
      assembler of Figures 1-3. For each method, the fraction of the hard molecules whose exact intron chain one of
      its transcripts overlapping a catalog copy carries (`score.py bakeoff-calls` -> `bakeoff-compare`), on all hard
      molecules, the contested stratum and the molecules whose chain is carried by >= 2 reads.
      Genome-wide (experiment B of the pre-registration; human A119b and gorilla OR6737, the samples with lab
      baselines): every multi-copy family of the sample's genome-wide catalog, in shards (figures/_o2.plan_shards);
      one point per family with >= 20 hard molecules plus the pooled fraction. Inset: the chr16 NPIP benchmark
      (docs/archive/2026-09/PREREG_hard_locus_bakeoff_2026-09-09.md), development.
"""
from __future__ import annotations

import collections
import csv
import datetime as dt
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import figlib
import _o2

META = {
    "id": "fig5",
    "title": "Copy assignment by copy identity (simulation) and the hard-locus benchmark (real reads)",
    "claim": ("Development tables. Simulation, human A119b chr16 catalog: scored within the read's source family "
              "(known only in simulation), the copy-assignment test made no wrong call where it assigned: 163 of the "
              "1,263 MAPQ-0 reads, from 21 source copies (95% interval resampling copies 0.85-1; 153 through a "
              "decisive site, 10 with a single candidate copy). It assigns 0.77 and 0.80 of the MAPQ-0 reads from "
              "copies 99-99.5% and 98-99% identical to their most similar copy, and 1 of 858 from copies identical "
              "over the aligned segment. The outputs a user can apply do not reproduce this: the default output is "
              "correct on 60 of 648 assigned (0.09; the aligner's primary 0.50), and the union test assigns none "
              "(1,255 of the 1,263 have an NM-identical twin at another locus). Real reads at the chr16 NPIP family "
              "(a hard set our own gate defines; Rustle's arm = copy_assign --gtf): all hard molecules Rustle 0.858, "
              "IsoSeq collapse 0.732, FLAIR 0.473, StringTie 0.441; contested stratum IsoSeq collapse 0.813 and "
              "Rustle 0.692. The genome-wide tables replace these numbers "
              "(docs/archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md)."),
    "tables": ["fig5_assign_accuracy_bands", "fig5_hard_locus", "fig5_hard_locus_transcripts"],
}
BANDS_TABLE = "fig5_assign_accuracy_bands"
HARD_TABLE = "fig5_hard_locus"
TX_TABLE = "fig5_hard_locus_transcripts"
FAMILY_TABLE = "fig5_hard_locus_families"   # genome-wide only: one row per (sample, family, method)
# the genome-wide family table exists only once experiment B has run; `make.py check` requires it from then on
if (figlib.DATA_DIR / f"{FAMILY_TABLE}.tsv").exists():
    META["tables"].append(FAMILY_TABLE)
# Supplementary Figure 5s (fig5s_margin_rule): Rustle against the alignment-score margin rule of the Eichler lab
# (docs/archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md, Amendment 2, experiment C). Simulation tables come from the
# same runs as a-b; the real-read table from one extra `copy_assign --union-certificate` run per experiment-B shard.
MR_STRATA_TABLE = "fig5s_margin_rule_sim"
MR_ACC_TABLE = "fig5s_margin_rule_accuracy"
MR_REAL_TABLE = "fig5s_margin_rule_real"
for _t in (MR_STRATA_TABLE, MR_ACC_TABLE, MR_REAL_TABLE):
    if (figlib.DATA_DIR / f"{_t}.tsv").exists():
        META["tables"].append(_t)
PROVISIONAL = ("provisional: from the recorded 2026-09-23/24 run (o2sim h16*/g44*, hash-seeded simulation); "
               "make.py data regenerates with the stable seed and current defaults")
PROVISIONAL_HARD = ("provisional: our arm is TWO recorded runs, not the current copy_assign: the transcripts are "
                    "ours.gtf (run `ours`, 2026-09-08 21:46, --gtf, no AS-tied gate) and the hard set is "
                    "ours_final.assignments.tsv (run `ours_final`, 2026-09-09 18:01, AS-tied gate); make.py data "
                    "rebuilds both in one current run")

# ---------------------------------------------------------------- hard-locus arms
# the other methods: (tool, GTF in cfg['hard_locus_dir'], junction fuzz) — the pre-registered setting (fuzz 0 for
# ours/StringTie/FLAIR, 5 for IsoSeq collapse, which runs --max-fuzzy-junction 5). Ours is listed first in every
# bakeoff-compare call because it treats the first arm as "ours".
HARD_TOOLS = [("stringtie", "stringtie_family.gtf", 0), ("flair", "flair_family.gtf", 0),
              ("isoseq", "isoseq_family.gff", 5)]
TOOL_FUZZ = {"rustle": 0, "stringtie": 0, "flair": 0, "isoseq": 5}
RUSTLE_FUZZ = 0
# recorded fallback: two runs (their products are what the pre-registered outcome scored)
RECORDED_GTF, RECORDED_ASSIGN = "ours.gtf", "ours_final.assignments.tsv"
# flags recovered from the recorded runs' params.tsv: `ours` = gtf true; `ours_final` = origin_drop_indels true
# (explicit then, the default since 2026-09-09 18:28). Everything else was, and is, the binary's default.
OUR_FLAGS = ["--gtf", "--origin-drop-indels"]
# the recorded run the rebuilt params.tsv is diffed against: `ours_final_gtf` (2026-09-09 19:12) is the ONE recorded
# process that wrote both a GTF and the gate rows (its assignments are byte-identical to ours_final's, md5 0590d544,
# and its GTF gives byte-identical bakeoff calls to ours.gtf; checked 2026-09-25)
RECORDED_REF_PARAMS = "ours_final_gtf.params.tsv"
# default changes since 2026-09-09 19:12 that params.tsv does not record (git log -- src/bin/copy_assign.rs)
UNRECORDED_DEFAULT_CHANGES = [
    "RUSTLE_JUNCTION_MAJORITY: the --gtf assembly's junction gate tolerates a minority of non-canonical junctions by "
    "default since 2026-09-21 (e44170f2); the recorded GTFs were built with strict canonical junctions",
]
# params.tsv rows that are run OUTCOMES (counts), not settings
COUNT_KEYS = {"origin_rejected", "orphans", "sole_candidates", "placement_assigned", "primary_local_rows",
              "readthrough_explained", "contested_rows", "junction_conflicts", "indel_psv_molecules",
              "indel_psv_columns", "indel_psv_columns_ge10"}
PATH_KEYS = {"families", "copies_fa"}

HARD_STRATA = [("all_molecules", "hard (all gate rows)", "All hard molecules"),
               ("all_molecules", "contested", "Contested"),
               ("chain_ge2_molecules", "hard (all gate rows)", "Chain carried by ≥ 2 reads")]
# panel c's Rustle arm is copy_assign's own transcripts (--families --gtf), not the Figure 1-3 assembler
HARD_LABEL = {"rustle": "Rustle (copy-assign.)", "isoseq": "IsoSeq coll."}
GW_TICK = {"rustle": "Rustle\n(copy-\nassignment)", "isoseq": "IsoSeq\ncollapse"}
# genome-wide experiment B: the samples with the lab's StringTie / FLAIR / IsoSeq collapse (samples.baseline)
HARD_GW_SAMPLES = ["human_A119b", "gorilla_OR6737"]
HARD_GW_MIN_FAMILY = 20      # families with >= this many hard molecules are drawn as points (pre-registered)
SUBSET_PAD = 100_000         # the per-shard subset BAM: copy spans +/- 100 kb (the NPIP hsa16.bam recipe)

METHOD_STYLE = {
    # open markers mean "Rustle (primaries only)" in every figure, so every reading here is a FILLED shape
    "family_certificate": {"color": figlib.TOOL_COLOR["rustle"], "marker": "o"},
    "per_family_table": {"color": figlib.TOOL_COLOR["rustle"], "marker": "P"},
    "union_certificate": {"color": figlib.TOOL_COLOR["rustle"], "marker": "v"},
    "aligner_primary": {"color": figlib.INK_3, "marker": "X"},
}
PSV_BAR = figlib.BLUE[100]   # share of MAPQ-0 reads with >= 1 decisive site in the source family (behind the points)
PSV_LABEL = "Bar: share with ≥ 1 decisive site in the source family*"
BAND_TICK = {"all": "All MAPQ-0 reads", "identical": "100%*", "99.5-100%": "99.5–100%",
             "99-99.5%": "99–99.5%", "98-99%": "98–99%", "<98%": "< 98%"}

# saved with bbox_inches="tight"; the provisional stamp reaches x = 0.995, so the canvas stays under 183 mm
FIG_W = figlib.WIDTH_DOUBLE - 6 * figlib.MM


# ================================================================ data: the NPIP hard-locus benchmark
def catalog_region(copies_tsv) -> str:
    """`chrom:min_start-max_end` over the catalog (one chromosome): the single region the recorded runs swept
    (their logs print chr16:11963320-80438591 for copies16.tsv)."""
    rows = [l.rstrip("\n").split("\t") for l in open(copies_tsv)]
    h = rows[0]
    ci, si, ei = h.index("chrom"), h.index("start"), h.index("end")
    chroms = {r[ci] for r in rows[1:]}
    if len(chroms) != 1:
        raise RuntimeError(f"{copies_tsv}: the hard-locus catalog spans {sorted(chroms)}; expected one chromosome")
    return f"{chroms.pop()}:{min(int(r[si]) for r in rows[1:])}-{max(int(r[ei]) for r in rows[1:])}"


def _read_params(path) -> dict:
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def params_diff(new_path, old_path) -> list[str]:
    """Settings and outcome counts of our rebuilt run that differ from the recorded one (params.tsv rows)."""
    new, old = _read_params(new_path), _read_params(old_path)
    sets, counts = [], []
    for k in list(old) + [k for k in new if k not in old]:
        if k in PATH_KEYS or new.get(k) == old.get(k):
            continue
        (counts if k in COUNT_KEYS else sets).append(f"{k} {old.get(k, 'absent')} -> {new.get(k, 'absent')}")
    return [f"settings that differ from the recorded run ({Path(old_path).name}): " + ("; ".join(sets) or "none"),
            "outcome counts that differ from the recorded run: " + ("; ".join(counts) or "none")]


def ensure_our_arm(cfg: dict, work: Path, *, force=False) -> dict:
    """ONE current `copy_assign` run on the hard-locus BAM that writes both our transcripts (--gtf) and the gate rows
    that define the hard set (HEAVY: a 68-Mb chr16 window, 415,854 records; the 2026-09-25 run took 88 s and
    3.1 GB peak RSS)."""
    hd = Path(cfg["hard_locus_dir"])
    bam, copies, copies_fa = hd / "hsa16.bam", hd / "copies16.tsv", hd / "copies16.fa"
    binary = Path(cfg["bin"]) / "copy_assign"
    d = (Path(work) / "rustle_arm").resolve()
    d.mkdir(parents=True, exist_ok=True)
    out = d / "ours"
    region = catalog_region(copies)
    cmd = [str(binary), "--bam", str(bam), "--fasta", cfg["human_fasta"], "--region", region, "--families",
           str(copies), "--copies-fa", str(copies_fa), *OUR_FLAGS, "--threads", str(cfg.get("threads", "4")),
           "--out", str(out)]
    stamp = d / "ours.cmd"
    gtf, assign, params = Path(f"{out}.gtf"), Path(f"{out}.assignments.tsv"), Path(f"{out}.params.tsv")
    srcs = (bam, copies, copies_fa, binary)
    if force or not (all(figlib.fresh(p, *srcs) for p in (gtf, assign, params)) and _o2._stamp_ok(stamp, " ".join(cmd))):
        # the binary's own defaults only: no RUSTLE_* override from the calling shell leaks into the run
        unset = [x for k in sorted(os.environ) if k.startswith("RUSTLE_") for x in ("-u", k)]
        timer = ["/usr/bin/time", "-f", "TIME elapsed=%e rss_kb=%M"] if os.path.exists("/usr/bin/time") else []
        figlib.run([*timer, "env", *unset, *cmd], log=d / "ours.log", cwd=d)
        stamp.write_text(" ".join(cmd))
    notes = [f"our arm REBUILT: {' '.join(cmd)} (one process writes ours.gtf and ours.assignments.tsv; "
             f"log {d / 'ours.log'})"]
    ref = hd / RECORDED_REF_PARAMS
    if ref.exists():
        notes += params_diff(params, ref)
    notes += [f"default change not recorded in params.tsv: {x}" for x in UNRECORDED_DEFAULT_CHANGES]
    rec_err = hd / "ours_final.err"
    if rec_err.exists():
        m = re.search(r"\]\s+(\S+:\d+-\d+): \d+ mapped reads", rec_err.read_text())
        if m and m.group(1) != region:
            raise RuntimeError(f"swept region {region} differs from the recorded run's {m.group(1)}")
        notes.append(f"swept region {region} = the catalog span, the same region the recorded runs swept")
    when = dt.datetime.fromtimestamp(assign.stat().st_mtime).strftime("%Y-%m-%d")
    return {"gtf": gtf, "assign": assign, "params": params, "notes": notes,
            "label": f"one current copy_assign run ({when}) gives transcripts and hard set"}


def _calls(gtf: Path, bam: Path, copies: Path, tool: str, fuzz: int, work: Path) -> tuple[Path, Path, Path]:
    out = work / f"{tool}_f{fuzz}"
    calls, tx, log = Path(f"{out}.calls.tsv"), Path(f"{out}.tx.tsv"), work / f"{tool}_f{fuzz}.log"
    if not (figlib.fresh(calls, gtf, bam, copies, _o2.SCORE) and figlib.fresh(tx, gtf, bam, copies, _o2.SCORE)
            and figlib.fresh(log, gtf, bam, copies, _o2.SCORE)):
        figlib.run([sys.executable, str(_o2.SCORE), "bakeoff-calls", str(gtf), str(bam), str(copies), "--label", tool,
                    "--fuzz", str(fuzz), "--out", str(out), "--tx-out", str(tx)], log=log)
    return calls, tx, log


def _hard_locus(cfg: dict, work: Path, *, recorded: bool, force=False):
    """Per-tool calls (`score.py bakeoff-calls --tx-out`, light: the 26 copy regions of one contig), then
    `score.py bakeoff-compare --tx-support` on our run's gate rows, without and with --min-mult 2."""
    hd = Path(cfg["hard_locus_dir"])
    bam, copies = hd / "hsa16.bam", hd / "copies16.tsv"
    work = work / ("recorded" if recorded else "current")
    work.mkdir(parents=True, exist_ok=True)
    if recorded:
        arm = {"gtf": hd / RECORDED_GTF, "assign": hd / RECORDED_ASSIGN,
               "label": "recorded runs ours (2026-09-08) + ours_final (2026-09-09)",
               "notes": [PROVISIONAL_HARD,
                         "the one recorded process that wrote both, ours_final_gtf (2026-09-09 19:12), has gate rows "
                         "byte-identical to ours_final (md5 0590d544) and a GTF whose bakeoff calls are "
                         "byte-identical to ours.gtf's (checked 2026-09-25), so the two files are consistent"]}
    else:
        arm = ensure_our_arm(cfg, work, force=force)
    inputs = {"hard_bam": bam, "hard_copies": copies, "hard_gate_rows": arm["assign"], "hard_gtf_rustle": arm["gtf"]}
    if arm.get("params"):
        inputs["hard_params_rustle"] = arm["params"]
    arms = [("rustle", arm["gtf"], RUSTLE_FUZZ)] + [(t, hd / g, f) for t, g, f in HARD_TOOLS]
    calls, txs, logs = [], [], {}
    for tool, gtf, fuzz in arms:
        if tool != "rustle":
            inputs[f"hard_gtf_{tool}"] = gtf
        c, tx, log = _calls(gtf, bam, copies, tool, fuzz, work)
        calls.append(f"{tool}={c}")
        txs += ["--tx-support", f"{tool}={tx}"]
        logs[tool] = log
    rows, tx_rows = [], []
    compare = {}
    for mset, extra in (("all_molecules", []), ("chain_ge2_molecules", ["--min-mult", "2", "--bam", str(bam)])):
        log = work / f"compare_{mset}.txt"
        srcs = (arm["assign"], *[c.split("=", 1)[1] for c in calls], *[t.split("=", 1)[1] for t in txs[1::2]],
                _o2.SCORE)
        if not figlib.fresh(log, *srcs):
            figlib.run([sys.executable, str(_o2.SCORE), "bakeoff-compare", "--assign", str(arm["assign"]), *extra,
                        *txs, *calls], log=log)
        inputs[f"hard_compare_{mset}"] = log
        rows += [["human_A119b", "human", "npip_chr16"] + r[1:] + [arm["label"]] for r in _parse_compare(log, mset)]
        compare[mset] = _parse_tx_support(log)
    for tool, _, fuzz in arms:
        tot, inc, multi = _parse_calls_log(logs[tool])
        a = compare["all_molecules"][tool]
        tx_rows.append(["human_A119b", "human", "npip_chr16", tool, fuzz, tot, inc, multi, a["any"], a["hard"],
                        a["any"] / inc if inc else None, a["hard"] / inc if inc else None,
                        arm["label"] if tool == "rustle" else "recorded (lab)"])
    notes = arm["notes"] + [
        "scope npip_chr16 = the hard-locus benchmark of PREREG_hard_locus_bakeoff_2026-09-09 (development): human "
        "A119b reads (hsa16.bam = samtools view -b -M -L copyregions.bed A119b.t2t.bam, command in its @PG header), "
        "chr16 NPIP family (copies16.tsv, 26 copies)",
        "the hard set is defined by Rustle: the rows of OUR copy_assign run's assignments.tsv (molecules its gate "
        "admits: >= 2 alignments in the swept region with the runner-up alignment score equal to the best, one in a "
        "catalog copy) that also have a primary alignment inside a copy (the scored universe: primaries -F 2308 in "
        "the merged copy intervals, identical for every method); the strata (contested, assigned, ...) are the "
        "copy-assignment test's own verdicts",
        "fraction_carried = derived_one|derived_multi (the method emits a transcript with the molecule's exact intron "
        "chain that overlaps a catalog copy, copy = max raw overlap; for a single-exon molecule the exact chain means "
        "span containment); junction tolerance 0 bp for rustle/stringtie/flair, 5 bp for isoseq (pre-registered)",
        "chain_ge2_molecules = --min-mult 2 --bam: molecules whose exact chain is carried by >= 2 reads; post hoc at "
        "NPIP (PREREG_hard_locus_bakeoff_2026-09-09 section 'Why P5 fails'), pre-registered for the genome-wide run",
        "the Rustle arm is copy_assign --families --gtf (the copy-assignment step's own transcripts), not the "
        "assembler (tools/rustle_pipeline.sh assemble) of Figures 1-3",
        "the 0.846 / 0.574 / 0.502 / 0.755 composite quoted in older notes is bench/isoform_bakeoff.py --summary "
        "(retired, 3,633 molecules) and is NOT what this table reports",
    ]
    tx_notes = [n for n in arm["notes"] if n.startswith("provisional")] + [
        "one row per method on the hard-locus benchmark (same runs as fig5_hard_locus); transcripts_total = "
        "transcripts in the method's GTF (npip: its family GTF; genome: the transcripts overlapping the sample's "
        "copies); transcripts_in_copy = those overlapping a catalog copy (bakeoff-calls); in_copy_spanning_ge2 = "
        "in-copy transcripts overlapping >= 2 copy intervals (conflation)",
        "supported_any / supported_hard = in-copy transcripts carrying the exact intron chain (fuzz as the arm; ends "
        "not compared) of >= 1 molecule of the scored universe / of the hard set (bakeoff-calls --tx-out, "
        "bakeoff-compare --tx-support); an assembler that emits only chains seen in >= N reads scores ~1 on "
        "supported_any by construction",
    ]
    return rows, tx_rows, inputs, notes, tx_notes


def _parse_compare(log: Path, mset: str) -> list[list]:
    rows, tools, on = [], None, False
    for line in open(log):
        if line.startswith("stratum"):
            toks = line.split()
            tools = toks[2:toks.index("(fraction")] if "(fraction" in toks else toks[2:]
            on = True
            continue
        if on:
            if not line.strip():
                break
            name = line[:34].strip()
            rest = line[34:].split()
            n = int(rest[0])
            for t, v in zip(tools, rest[1:]):
                rows.append(["human", mset, name, n, t, float(v) if v != "-" else float("nan")])
    if not rows:
        raise RuntimeError(f"no strata parsed from {log}")
    return rows


_TXS = re.compile(r"^\s+(\S+)\s+in-copy\s+(\d+)\s+any\s+(\d+) \([\d.]+\)\s+hard\s+(\d+) ")


def _parse_tx_support(log: Path) -> dict:
    out, on = {}, False
    for line in open(log):
        if line.startswith("transcript support:"):
            on = True
            continue
        m = _TXS.match(line) if on else None
        if m:
            out[m.group(1)] = {"in_copy": int(m.group(2)), "any": int(m.group(3)), "hard": int(m.group(4))}
    if not out:
        raise RuntimeError(f"no transcript-support section in {log}")
    return out


def _parse_calls_log(log: Path) -> tuple[int, int, int]:
    txt = open(log).read()
    m = re.search(r"transcripts (\d+)\s+in a copy (\d+)", txt)
    k = re.search(r"transcripts spanning >=2 copies (\d+)", txt)
    return int(m.group(1)), int(m.group(2)), int(k.group(1)) if k else 0


# ================================================================ data: the genome-wide hard-locus benchmark
_TID = re.compile(r'transcript_id[ =]"?([^";]*)"?')


def _multi_copy_catalog(cat_tsv, cat_fa, out_prefix: Path) -> tuple[Path, Path]:
    """The catalog restricted to families with >= 2 copies (pre-registration, experiment B)."""
    tsv, fa = Path(f"{out_prefix}.copies.tsv"), Path(f"{out_prefix}.copies.fa")
    if figlib.fresh(tsv, cat_tsv, cat_fa, __file__) and figlib.fresh(fa, cat_tsv, cat_fa, __file__):
        return tsv, fa
    with open(cat_tsv) as fh:
        header = fh.readline()
        rows = [l for l in fh if l.strip()]
    n = collections.Counter(l.split("\t", 1)[0] for l in rows)
    keep = [l for l in rows if n[l.split("\t", 1)[0]] >= 2]
    fams = {l.split("\t", 1)[0] for l in keep}
    recs, order = _o2._fa_records(cat_fa)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    _o2._write_if_changed(tsv, header + "".join(keep))
    _o2._write_if_changed(fa, "".join("".join(recs[k]) for k in order if k[0] in fams))
    return tsv, fa


def _merged_spans(copies_tsv, pad: int, lengths: dict) -> dict:
    """chrom -> merged [lo, hi) of the copy spans +/- pad."""
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(copies_tsv), delimiter="\t"):
        c = r["chrom"]
        by[c].append((max(0, int(r["start"]) - pad), min(lengths.get(c, 1 << 62), int(r["end"]) + pad)))
    out = {}
    for c, v in by.items():
        m = []
        for lo, hi in sorted(v):
            if m and lo <= m[-1][1]:
                m[-1][1] = max(m[-1][1], hi)
            else:
                m.append([lo, hi])
        out[c] = m
    return out


def _bakeoff_copies(shard_tsv: Path, out_tsv: Path) -> dict:
    """The shard's copies with copy_idx renumbered 0..n-1 (bakeoff-calls keys copies by copy_idx, which repeats
    across families); returns {new copy_idx: (family_id, catalog copy_idx)} and writes `out_tsv` (+ .map.tsv)."""
    with open(shard_tsv) as fh:
        header = fh.readline()
        cols = header.rstrip("\n").split("\t")
        ci = cols.index("copy_idx")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    mapping, body = {}, []
    for i, f in enumerate(rows):
        mapping[i] = (f[0], f[ci])
        g = list(f)
        g[ci] = str(i)
        body.append("\t".join(g) + "\n")
    _o2._write_if_changed(out_tsv, header + "".join(body))
    _o2._write_if_changed(Path(f"{out_tsv}.map.tsv"), "copy_idx\tfamily_id\tcatalog_copy_idx\n"
                          + "".join(f"{k}\t{a}\t{b}\n" for k, (a, b) in mapping.items()))
    return mapping


def _split_gtf(gtf, shard_spans: dict, out_name: str, stamp_key: str):
    """Write, into every shard directory, the lines of the transcripts of `gtf` that overlap that shard's copy
    spans (two streaming passes: transcript spans, then lines). bakeoff-calls ignores every transcript that overlaps
    no copy, so the per-shard file changes no call. shard_spans = {shard_dir: {chrom: [[lo, hi], ...]}}."""
    import bisect
    import gzip
    op = (lambda p: gzip.open(p, "rt")) if str(gtf).endswith(".gz") else (lambda p: open(p))
    span = {}
    with op(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "exon":
                continue
            m = _TID.search(f[8])
            if not m:
                continue
            t = m.group(1)
            s, e = int(f[3]) - 1, int(f[4])
            x = span.get(t)
            span[t] = (f[0], s, e) if x is None else (f[0], min(x[1], s), max(x[2], e))
    idx = collections.defaultdict(list)   # chrom -> [(lo, hi, shard_dir)] sorted
    for sd, by in shard_spans.items():
        for c, iv in by.items():
            idx[c] += [(lo, hi, sd) for lo, hi in iv]
    ends = {}
    for c in idx:
        idx[c].sort()
        ends[c] = max(hi - lo for lo, hi, _ in idx[c])
    starts = {c: [x[0] for x in v] for c, v in idx.items()}
    where = collections.defaultdict(set)
    for t, (c, s, e) in span.items():
        v = idx.get(c)
        if not v:
            continue
        j = bisect.bisect_left(starts[c], e)
        k = bisect.bisect_left(starts[c], s - ends[c])
        for lo, hi, sd in v[k:j]:
            if lo < e and hi > s:
                where[t].add(sd)
    outs = {sd: open(Path(sd) / f"{out_name}.tmp", "w") for sd in shard_spans}
    with op(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = _TID.search(line)
            if m and m.group(1) in where:
                for sd in where[m.group(1)]:
                    outs[sd].write(line)
    for sd, fo in outs.items():
        fo.close()
        (Path(sd) / f"{out_name}.tmp").replace(Path(sd) / out_name)
        (Path(sd) / f"{out_name}.key").write_text(stamp_key)


def _tally_shard(sd: Path, tools: list, calls_paths: dict, assign_p: Path, sub_bam: Path, mapping: dict):
    """Per-molecule strata and per-family tallies of one shard, by bakeoff-compare's own definitions: hard = the
    molecules of the scored universe with a row in our assignments (the LAST row per molecule, as bakeoff-compare
    reads it); contested = that row not origin-rejected with >= 2 candidates; carried = derived_one|derived_multi;
    chain_ge2 = the molecule's exact chain (-F 2308, chrom + introns) is carried by >= 2 reads of the subset BAM.
    A molecule's family = the family of its row with primary_local = 1 (its primary overlaps a copy of that family);
    failing that, the family of the last row."""
    import score as sc   # bench/score.py (on sys.path through _o2)
    import lib
    calls = {t: sc.load_calls(p) for t, p in calls_paths.items()}
    allm = set()
    for c in calls.values():
        allm |= set(c)
    last, prim_fam = {}, {}
    for r in csv.DictReader(open(assign_p), delimiter="\t"):
        last[r["read_name"]] = r
        if r.get("primary_local") == "1" and r["read_name"] not in prim_fam:
            prim_fam[r["read_name"]] = r["family_id"]
    chain = {}
    for ln in lib.sam_lines(["-F", "2308", str(sub_bam)]):
        f = ln.split("\t", 6)
        chain.setdefault(f[0], (f[2],) + lib.cigar_introns(int(f[3]) - 1, f[5]))
    mult = collections.Counter(chain.values())
    ge2 = {m for m in allm if mult.get(chain.get(m), 0) >= 2}
    hard = {m for m in allm if m in last}
    contested = {m for m in hard if last[m]["origin_rejected"] == "0" and int(last[m]["n_candidates"]) >= 2}
    carried = lambda t, m: calls[t].get(m, ("derived_none", []))[0] != "derived_none"  # noqa: E731
    strata = {("all_molecules", "hard (all gate rows)"): hard, ("all_molecules", "contested"): contested,
              ("chain_ge2_molecules", "hard (all gate rows)"): hard & ge2,
              ("chain_ge2_molecules", "contested"): contested & ge2}
    pooled = {k: (len(S), {t: sum(carried(t, m) for m in S) for t in tools}) for k, S in strata.items()}
    fam = {m: prim_fam.get(m) or last[m]["family_id"] for m in hard}
    per_fam = collections.defaultdict(lambda: collections.Counter())
    for m in hard:
        f = fam[m]
        per_fam[f]["n_hard"] += 1
        per_fam[f]["n_contested"] += m in contested
        for t in tools:
            per_fam[f][f"{t}_hard"] += carried(t, m)
            per_fam[f][f"{t}_contested"] += carried(t, m) and m in contested
    n_copies = collections.Counter(v[0] for v in mapping.values())
    return pooled, per_fam, n_copies, allm, hard


def _hard_gw_plan(cfg: dict, sid: str):
    """Experiment B's shard plan of one sample (cached; the record pass is bounded per call and raises Pending):
    (work dir, multi-copy catalog tsv / fa, plan, BAM, FASTA, contig lengths). Shared by ensure_hard_gw and
    ensure_margin_real, so both use the same shards."""
    import samples
    row = samples.get(cfg, sid)
    bam, fasta = row["bam"], row["fasta"]
    cat_tsv, cat_fa = _o2.catalog_paths(cfg, sid)
    d = figlib.work_dir(cfg, "fig5") / "hard_gw" / sid
    d.mkdir(parents=True, exist_ok=True)
    mc_tsv, mc_fa = _multi_copy_catalog(cat_tsv, cat_fa, d / "cat.multi")
    header = _o2.bam_contigs(bam)
    lengths = dict(header)
    carrying = {r["chrom"] for r in csv.DictReader(open(mc_tsv), delimiter="\t")}
    contigs = [(c, n) for c, n in header if c in carrying]
    plan = _o2.plan_shards(bam, mc_tsv, mc_fa, contigs, d / "shards", kind="real",
                           budget_s=float(cfg.get("o2_shard_budget_s", _o2.SHARD_BUDGET_S)),
                           deadline=time.time() + float(cfg.get("hard_pass_budget_s", "420")),
                           sample_frac=float(cfg["hard_gw_sample_frac"]) if cfg.get("hard_gw_sample_frac") else None)
    return d, mc_tsv, mc_fa, plan, bam, fasta, lengths


def ensure_hard_gw(cfg: dict, sid: str, budget: _o2.Budget, *, force=False) -> dict:
    """Experiment B on one sample (HEAVY; bounded per call): plan shards of the multi-copy families over the real BAM
    (_o2.plan_shards kind 'real'; the record pass is cached per contig), then per shard: copy_assign --gtf, a subset
    BAM (copy spans +/- 100 kb), the lab GTFs split per shard, bakeoff-calls per method, bakeoff-compare (all and
    --min-mult 2) as a cross-check, and the per-family tallies. Returns pooled / per-family / transcript rows."""
    import samples
    tools_lab = {t: samples.baseline(cfg, sid, t) for t in ("stringtie", "flair", "isoseq")}
    if not all(tools_lab.values()):
        raise SystemExit(f"[fig5] {sid}: no lab baselines in the registry; experiment B needs StringTie, FLAIR and "
                         f"IsoSeq collapse")
    d, mc_tsv, mc_fa, plan, bam, fasta, lengths = _hard_gw_plan(cfg, sid)
    notes = [f"{sid}: {plan['summary']}"]
    binary = Path(cfg["bin"]) / "copy_assign"
    threads = str(cfg.get("threads", "4"))
    shard_spans = {s["dir"]: _merged_spans(s["tsv"], 0, lengths) for s in plan["shards"]}
    # the lab GTFs, split per shard once (one heavy step per method)
    for tool, gtf in tools_lab.items():
        key = f"{gtf}|{os.path.getsize(gtf)}|{int(os.path.getmtime(gtf))}|{plan['summary']}"
        name = f"{tool}.gtf"
        if force or not all(_o2._stamp_ok(Path(sd) / f"{name}.key", key) for sd in shard_spans):
            budget.take(f"{sid}: split the {tool} GTF over {len(shard_spans)} shards")
            _split_gtf(gtf, shard_spans, name, key)
    pooled = collections.defaultdict(lambda: [0, collections.Counter()])
    fam_rows, tx_acc = [], collections.defaultdict(lambda: {"in_copy": set(), "any": set(), "hard": set(),
                                                            "multi": 0})
    tools = ["rustle", "stringtie", "flair", "isoseq"]
    for s in plan["shards"]:
        sd = Path(s["dir"])
        out = sd / "ours"
        cmd = [str(binary), "--bam", bam, "--fasta", fasta, "--regions", s["regions"], "--families", s["tsv"],
               "--copies-fa", s["fa"], *OUR_FLAGS, "--threads", threads, "--out", str(out)]
        _o2._run_copy_assign(cmd, out, [bam, s["tsv"], s["fa"], s["regions"], binary], budget,
                             f"{sid}: copy_assign --gtf shard {s['name']}", force=force)
        sub = sd / "sub.bam"
        bed = sd / "copyregions.bed"
        spans = _merged_spans(s["tsv"], SUBSET_PAD, lengths)
        _o2._write_if_changed(bed, "".join(f"{c}\t{lo}\t{hi}\n" for c in sorted(spans) for lo, hi in spans[c]))
        if force or not (figlib.fresh(sub, bam, bed) and Path(f"{sub}.bai").exists()):
            budget.take(f"{sid}: subset BAM of shard {s['name']}")
            figlib.run(f"samtools view -b -M -L {bed} -o {sub}.tmp {bam} && mv {sub}.tmp {sub} && samtools index {sub}",
                       log=sd / "sub.log")
        copies_bk = sd / "copies.bakeoff.tsv"
        mapping = _bakeoff_copies(Path(s["tsv"]), copies_bk)
        gtfs = {"rustle": Path(f"{out}.gtf"), **{t: sd / f"{t}.gtf" for t in tools_lab}}
        calls, txs, logs = {}, [], {}
        for t in tools:
            c, tx, log = _calls(gtfs[t], sub, copies_bk, t, TOOL_FUZZ[t], sd)
            calls[t], logs[t] = c, log
            txs += ["--tx-support", f"{t}={tx}"]
        p_sh, fam_sh, n_copies, allm, hard = _tally_shard(sd, tools, calls, Path(f"{out}.assignments.tsv"), sub,
                                                          mapping)
        # cross-check against bakeoff-compare on this shard (n exact, fraction to its printed 3 decimals)
        for mset, extra in (("all_molecules", []), ("chain_ge2_molecules", ["--min-mult", "2", "--bam", str(sub)])):
            log = sd / f"compare_{mset}.txt"
            if not figlib.fresh(log, Path(f"{out}.assignments.tsv"), *calls.values(), _o2.SCORE):
                figlib.run([sys.executable, str(_o2.SCORE), "bakeoff-compare", "--assign",
                            str(out) + ".assignments.tsv", *extra, *txs, *[f"{t}={calls[t]}" for t in tools]], log=log)
            for r in _parse_compare(log, mset):
                key = (mset, "contested" if r[2] == "contested" else r[2])
                if key not in p_sh:
                    continue
                n, k = p_sh[key]
                if r[3] != n or (n and abs(k[r[4]] / n - r[5]) > 0.0005 + 1e-9):
                    raise RuntimeError(f"{sd}: per-molecule tally disagrees with bakeoff-compare ({key}, {r[4]}): "
                                       f"{k[r[4]]}/{n} vs {r[5]} of {r[3]}")
        for key, (n, k) in p_sh.items():
            pooled[key][0] += n
            pooled[key][1].update(k)
        for f, c in fam_sh.items():
            for t in tools:
                fam_rows.append([sid, _o2.species_of(cfg, sid), f, n_copies.get(f, 0), c["n_hard"], c["n_contested"],
                                 t, c[f"{t}_hard"], c[f"{t}_contested"]])
        for t in tools:
            tot, inc, multi = _parse_calls_log(logs[t])
            tx_acc[t]["multi"] += multi
            with open(sd / f"{t}_f{TOOL_FUZZ[t]}.tx.tsv") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    tid = (s["name"], r["transcript_id"]) if t == "rustle" else r["transcript_id"]
                    tx_acc[t]["in_copy"].add(tid)
                    ms = {x for x in r["molecules"].split(",") if x}
                    if ms & allm:
                        tx_acc[t]["any"].add(tid)
                    if ms & hard:
                        tx_acc[t]["hard"].add(tid)
    label = f"{len(plan['shards'])} shards of copy_assign --families --gtf (current binary)"
    hard_rows = []
    for (mset, stratum), (n, k) in sorted(pooled.items()):
        for t in tools:
            hard_rows.append([sid, _o2.species_of(cfg, sid), "genome", mset, stratum, n, t,
                              k[t] / n if n else float("nan"), label])
    tx_rows = []
    for t in tools:
        a = tx_acc[t]
        inc = len(a["in_copy"])
        tx_rows.append([sid, _o2.species_of(cfg, sid), "genome", t, TOOL_FUZZ[t], None, inc, a["multi"], len(a["any"]),
                        len(a["hard"]), len(a["any"]) / inc if inc else None, len(a["hard"]) / inc if inc else None,
                        label if t == "rustle" else "recorded (lab)"])
    if plan.get("sampled"):
        notes.append(f"{sid}: SAMPLED read-connected components (fraction {plan['sampled']}, seed 20260925; "
                     f"pre-registered stop rule)")
    return {"hard": hard_rows, "families": fam_rows, "tx": tx_rows, "notes": notes,
            "inputs": {f"{sid}_hard_catalog": mc_tsv, f"{sid}_hard_plan": d / "shards" / "plan.json"}}


# ================================================================ data: Fig. 5s, the alignment-score margin rule
MR_REAL_HEADER = ["sample", "species", "identity_band", "threshold", "stratum", "n_reads", "n_mapq0", "n_margin0"]


def ensure_margin_real(cfg: dict, sid: str, budget: _o2.Budget, *, force=False) -> dict:
    """Experiment C on real reads (Amendment 2, C.3; HEAVY, bounded per call): on experiment B's shards, one
    `copy_assign --union-certificate` run per shard (current defaults, RUSTLE_* unset), then `score.py eichler --real`
    over the reads whose primary alignment overlaps a copy of the shards' families: the margin rule from the sample's
    as_table (genome-wide best / second AS) and the BAM's alignments over the copies, Rustle's answer = the primary's
    copy (MAPQ > 0) or the union test's result (MAPQ 0). Returns the per-read TSV and notes."""
    import samples
    d, mc_tsv, mc_fa, plan, bam, fasta, lengths = _hard_gw_plan(cfg, sid)
    binary = Path(cfg["bin"]) / "copy_assign"
    threads = str(cfg.get("threads", "4"))
    prefixes = []
    for s in plan["shards"]:
        sd = Path(s["dir"])
        out = sd / "u2"
        cmd = [str(binary), "--bam", bam, "--fasta", fasta, "--regions", s["regions"], "--families", s["tsv"],
               "--copies-fa", s["fa"], "--union-certificate", "--threads", threads, "--out", str(out)]
        _o2._run_copy_assign(cmd, out, [bam, s["tsv"], s["fa"], s["regions"], binary], budget,
                             f"{sid}: copy_assign --union-certificate shard {s['name']}", force=force)
        prefixes.append(str(out))
    md = d / "margin_rule"
    md.mkdir(parents=True, exist_ok=True)
    # the catalog actually run: the shards' families (all multi-copy families, or the sampled components)
    cat = md / "catalog.tsv"
    header, body = None, []
    for s in plan["shards"]:
        with open(s["tsv"]) as fh:
            h = fh.readline()
            header = header or h
            body += [l for l in fh if l.strip()]
    _o2._write_if_changed(cat, header + "".join(body))
    molecules = samples.product(cfg, sid, "assemble", "molecules")
    if not molecules.exists():
        raise SystemExit(f"[fig5] {sid}: no as_table ({molecules}); run `make.py runs --sample {sid} --stage assemble`")
    per, log = md / "per_read.tsv", md / "eichler_real.txt"
    srcs = [cat, molecules, bam, _o2.SCORE, *[f"{p}.assignments.tsv" for p in prefixes]]
    if force or not (figlib.fresh(per, *srcs) and figlib.fresh(log, *srcs)):
        budget.take(f"{sid}: margin-rule pass (BAM over the copies + as_table)")
        rc = figlib.run([sys.executable, str(_o2.SCORE), "eichler", "--real", "--bam", bam, "--catalog", str(cat),
                         "--as-table", str(molecules), "--union", *prefixes, "--cache", str(md / "cache"),
                         "--budget-s", str(cfg.get("margin_real_budget_s", "420")), "--per-read", str(per)],
                        log=log, check=False)
        if rc == 75:
            raise _o2.Pending(f"{sid}: the margin-rule pass over the BAM has contigs left ({md / 'cache'})")
        if rc != 0:
            raise RuntimeError(f"score.py eichler --real failed ({rc}); see {log}")
    notes = [f"{sid}: {plan['summary']}; union test = {len(prefixes)} copy_assign --union-certificate runs, one per "
             f"experiment-B shard"]
    if plan.get("sampled"):
        notes.append(f"{sid}: SAMPLED read-connected components (fraction {plan['sampled']}, seed 20260925)")
    return {"per_read": per, "log": log, "notes": notes,
            "inputs": {f"{sid}_margin_rule_catalog": cat, f"{sid}_as_table": molecules}}


def margin_rule_real_rows(per_read, sid: str, species: str) -> list[list]:
    """Strata of the real-read comparison per identity band (of the copy the primary alignment overlaps), T."""
    import score as sc
    c = collections.defaultdict(lambda: [0, 0, 0])
    with open(per_read) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            b = figlib.identity_band(_o2.fnum(r["identity"])) or "unknown"
            for T in _o2.MR_THRESHOLDS:
                s = sc.mr_stratum(r, T, "u")
                for bb in ("all", b):
                    x = c[(bb, T, s)]
                    x[0] += 1
                    x[1] += r["mapq"] == "0"
                    x[2] += r["margin"] == "0"
    out = []
    for b in ["all"] + figlib.IDENTITY_BANDS + ["unknown"]:
        if not any(k[0] == b for k in c):
            continue
        for T in _o2.MR_THRESHOLDS:
            for s in sc.MR_STRATA:
                n, z, m0 = c.get((b, T, s), [0, 0, 0])
                out.append([sid, species, b, T, s, n, z, m0])
    return out


MR_NOTES = [
    "experiment C of docs/archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md (Amendment 2): margin rule MR(T) = assign a "
    "read to its best-scoring alignment iff it has no other alignment or best AS - second AS >= T, over every mapped "
    "non-supplementary alignment of the read in the whole genome (-F 2052; missing AS = 0); T = 10 (headline), 1, 20; "
    "computed by bench/score.py eichler --sim from the same sim.bam as Figs 4-5, not from copy_assign --eichler-margin "
    "(region-local, counts supplementary records as rivals under --no-as-tied-only)",
    "Rustle's answer: primary MAPQ > 0 = the aligner's primary alignment (the copy-assignment step leaves these reads to "
    "the aligner); MAPQ 0 = the copy-assignment result, reading u = the union test (any assigned result of the "
    "--union-certificate run; what a user can apply), reading s = the default run's result for the read's SOURCE family "
    "(score.py reads OWN; known only in simulation); unmapped = unassigned for both",
    "a placement's copy = the catalog copy with the largest raw overlap of the alignment's reference span (ties: first "
    "in catalog order; the Fig. 5 aligner rule); correct = the source copy or a copy at the same locus "
    "(score.same_locus); another catalog copy or outside every catalog copy = wrong (for both rules)",
    "strata: both_same = both assign, same placement (same alignment, or same catalog copy / locus); both_differ = both "
    "assign, different placement; margin_only = the margin rule assigns, Rustle does not; rustle_aligner = the margin "
    "rule discards, Rustle keeps the aligner's primary (MAPQ > 0); rustle_test = the margin rule discards, the "
    "copy-assignment test assigns (MAPQ 0); neither; n_margin0 = reads whose best and second AS are equal",
    "fraction_assigned = assigned / all simulated reads of the row (identity_band 'all' = every band); fraction_correct = "
    "correct / assigned; acc_lo / acc_hi / acc_ci = _o2.copy_level_ci (resampling source copies); rustle_only_* rows: "
    "the reads of that stratum (fraction_assigned = their share of the row's reads); every MAPQ-0 verdict and every "
    "aligner placement was checked read by read against Figs 4-5's own records (score.py reads --per-read, "
    "_o2.per_read)",
]
MR_REAL_NOTES = [
    "real reads (no truth: agreement and coverage, not correctness): the reads whose primary alignment (-F 2308) "
    "overlaps a copy of a family with >= 2 copies in the sample's genome-wide catalog (experiment B's universe); margin "
    "rule from the sample's as_table (genome-wide best / second AS over every non-supplementary alignment) plus the "
    "BAM's alignments over the copies (where the best one lies: a copy, or outside the catalog); Rustle = the "
    "primary's copy (MAPQ > 0) or the union test's result (MAPQ 0; one copy_assign --union-certificate run per "
    "experiment-B shard, current defaults, RUSTLE_* unset); strata as in fig5s_margin_rule_sim (reading u); "
    "identity_band = the band of the copy the primary alignment overlaps",
]


def build_margin_rule_sim(collected: list, data_dir: Path) -> tuple:
    """fig5s_margin_rule_sim + fig5s_margin_rule_accuracy from _o2.collect()'s runs (light: one score.py eichler --sim
    pass over each sample's sim.bam, cached)."""
    srows, arows, notes, inputs = [], [], [], {}
    for runs, logs, reads in collected:
        sid = runs["sample"]
        per, log = _o2.margin_rule_per_read(runs, Path(logs["o2"]).parent)
        rows = _o2.margin_rule_join(per, reads)
        s, a = _o2.margin_rule_tables(rows, sid, runs["species"], runs["catalog_scope"])
        srows += s
        arows += a
        notes += _o2.margin_rule_notes(sid, rows)
        notes += [f"{sid}: {runs['source']}"] + [f"{sid}: {n}" for n in runs.get("notes", [])]
        notes.append(f"{sid}: per-read join checked against _o2.per_read and score.py reads --per-read "
                     f"({len(rows)} simulated reads, {sum(r['mapq'] == '0' for r in rows)} with MAPQ 0)")
        inputs.update({f"{sid}_sim_bam": f"{runs['sim']}.bam", f"{sid}_catalog": runs["catalog"],
                       f"{sid}_assignments_default": f"{runs['o2']}.assignments.tsv",
                       f"{sid}_assignments_union": f"{runs['u2']}.assignments.tsv", f"{sid}_margin_rule": per,
                       f"{sid}_margin_rule_log": log})
    gen = "figures/fig_assign_accuracy.py build_margin_rule_sim()"
    p1 = figlib.write_table(MR_STRATA_TABLE, _o2.MR_STRATA_HEADER, srows, generator=gen, inputs=inputs,
                            notes=MR_NOTES + notes, data_dir=data_dir)
    p2 = figlib.write_table(MR_ACC_TABLE, _o2.MR_ACC_HEADER, arows, generator=gen, inputs=inputs,
                            notes=MR_NOTES + notes, data_dir=data_dir)
    return p1, p2


def margin_rule_cli(argv=None):
    """`python3 figures/fig_assign_accuracy.py margin-rule [--inputs F] [--set KEY=VALUE ...]`: the Fig. 5s simulation
    tables from the runs AS BUILT (runs nothing heavy: no simulation, no copy_assign; one score.py eichler --sim pass
    per sample). `make.py data fig5` builds the same tables together with the rest of Fig. 5."""
    import argparse
    ap = argparse.ArgumentParser(description=margin_rule_cli.__doc__.split("\\n")[0])
    ap.add_argument("--inputs")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--data", default=str(figlib.DATA_DIR), help="where to write the two tables")
    a = ap.parse_args(argv)
    cfg = figlib.load_inputs(a.inputs)
    for kv in a.set:
        k, _, v = kv.partition("=")
        cfg[k] = v
    collected = []
    for sid in _o2.sample_ids(cfg):
        try:
            runs = _o2.runs_as_built(cfg, sid)
        except SystemExit as e:
            print(f"[fig5s] {sid}: skipped ({e})", file=sys.stderr)
            continue
        log_dir = figlib.work_dir(cfg, "o2sim") / (runs["species"] if _o2.o2_scope(cfg) == "dev" else sid) / "score"
        logs = _o2.score_logs(runs, log_dir)
        collected.append((runs, logs, _o2.per_read(runs, logs)))
    for p in build_margin_rule_sim(collected, Path(a.data)):
        print(p)


# ================================================================ data: build
BANDS_NOTES = [
    "reads = simulated reads whose primary alignment has MAPQ 0; identity_band 'all' = every band",
    "family_certificate = the default run's result for the read's SOURCE family (score.py reads OWN; the simulation "
    "picks that result: the test's own accuracy, not available to a user); per_family_table = any assigned result "
    "of the default output (score.py reads ANY; two loci claimed = conflict); union_certificate = any assigned "
    "result of the --union-certificate run (the union test); aligner_primary = the catalog copy the primary "
    "alignment overlaps most (correct = the source copy's locus by score.same_locus, wrong = another catalog copy, "
    "not_scored = no catalog copy)",
    "coverage = fraction assigned = (correct + wrong + conflict) / n_reads; accuracy = fraction correct among "
    "assigned = correct / (correct + wrong + conflict); abstain = a result with no copy assigned; not_scored = no "
    "result (or, for the aligner, a primary outside every catalog copy); verdicts from score.py reads --per-read",
    f"acc_lo / acc_hi: 95% interval with the SOURCE COPY as the unit (reads from one copy are not independent): "
    f"none when fewer than {_o2.MIN_CI_COPIES} copies are assigned (acc_ci too_few_copies); else a bootstrap resampling "
    f"source copies ({_o2.BOOT_B} resamples, seed {_o2.BOOT_SEED}, percentile), or, when every copy has the same ratio "
    f"(e.g. every assigned read correct) and the bootstrap is degenerate, the Wilson interval at the pooled ratio "
    f"with the n_copies_assigned copies as the sample size (acc_ci says which); n_copies = distinct source copies of "
    f"the band's MAPQ-0 reads",
    "n_family_psv (family_certificate rows only) = MAPQ-0 reads whose source family's result has >= 1 decisive site "
    "(a PSV or splice junction covered by the read at which the candidate copies differ; score.py reads --per-read "
    "own_n_decisive >= 1); n_nm_twin (union_certificate rows only) = MAPQ-0 reads with an NM-identical twin "
    "(_o2.twin_state nm_twin: an alignment not overlapping the source copy with the best alignment score and the "
    "source-copy alignment's NM)",
]


def build(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    """Regenerate the tables. Default (make.py data fig5): the simulation and both copy_assign runs of every sample
    (_o2.collect, HEAVY unless cached; bounded per call), the NPIP hard-locus arm rebuilt in one current copy_assign
    run (ensure_our_arm), and — genome scope — experiment B on the samples with lab baselines (ensure_hard_gw,
    HEAVY, bounded per call; cfg `fig5c_samples` narrows it, `none` skips it). Re-run until it completes.

    recorded=True tabulates the recorded runs instead (development / provisional tables): cfg['o2sim_dir'] for
    a/b, and the recorded 2026-09-08/09 products in cfg['hard_locus_dir'] for c."""
    budget = _o2.budget_from(cfg)
    rows, notes, inputs = [], [], {}
    collected = _o2.collect(cfg, force=force, recorded=recorded, budget=budget)
    for runs, logs, reads in collected:
        sid = runs["sample"]
        notes += _o2.run_notes(runs, reads, logs)
        rows += _o2.band_rows(reads, runs["catalog_scope"])
        inputs.update({f"{sid}_sim_bam": f"{runs['sim']}.bam", f"{sid}_catalog": runs["catalog"],
                       f"{sid}_assignments_default": f"{runs['o2']}.assignments.tsv",
                       f"{sid}_assignments_union": f"{runs['u2']}.assignments.tsv",
                       f"{sid}_score_reads_o2": logs["o2"], f"{sid}_score_reads_u2": logs["u2"],
                       f"{sid}_per_read_o2": logs["o2_per_read"], f"{sid}_per_read_u2": logs["u2_per_read"]})
    notes += BANDS_NOTES
    if recorded:
        notes.insert(0, PROVISIONAL)
    hrows, txrows, hinputs, hnotes, txnotes = _hard_locus(cfg, figlib.work_dir(cfg, "fig5") / "hard_locus",
                                                          recorded=recorded, force=force)
    fam_rows = None
    gw = [] if recorded or _o2.o2_scope(cfg) == "dev" else [
        x.strip() for x in cfg.get("fig5c_samples", ",".join(HARD_GW_SAMPLES)).split(",") if x.strip() not in ("", "none")]
    mr_real, mr_notes, mr_inputs = [], [], {}
    if gw:
        fam_rows = []
        for sid in gw:
            try:
                r = ensure_hard_gw(cfg, sid, budget, force=force)
                # Fig. 5s real reads (experiment C): the union test on the same shards (cfg margin_rule_real=0 skips)
                m = (ensure_margin_real(cfg, sid, budget, force=force)
                     if cfg.get("margin_rule_real", "1") != "0" else None)
            except _o2.Pending as e:
                raise SystemExit(f"[fig5] not finished — {e}. Re-run the same command (under flock) to continue.")
            hrows += r["hard"]
            txrows += r["tx"]
            fam_rows += r["families"]
            hnotes += r["notes"]
            hinputs.update(r["inputs"])
            if m:
                mr_real += margin_rule_real_rows(m["per_read"], sid, _o2.species_of(cfg, sid))
                mr_notes += m["notes"]
                mr_inputs.update(m["inputs"], **{f"{sid}_margin_rule_real": m["per_read"],
                                                 f"{sid}_margin_rule_real_log": m["log"]})
        hnotes.append("scope genome = experiment B of PREREG_genome_wide_copy_assignment_2026-09-25: every multi-copy "
                      "family of the sample's genome-wide catalog, copy_assign --families --gtf --origin-drop-indels in "
                      "shards of read-connected families (regions = every contig carrying a copy), each method scored "
                      "on a per-shard subset BAM (copy spans +/- 100 kb) with copy ids renumbered across families; "
                      "fractions pooled over every hard molecule, cross-checked shard by shard against "
                      "bakeoff-compare")
    p1 = figlib.write_table(BANDS_TABLE, _o2.BAND_HEADER, rows, generator="figures/fig_assign_accuracy.py build()",
                            inputs=inputs, notes=notes, data_dir=data_dir)
    p2 = figlib.write_table(HARD_TABLE, ["sample", "species", "scope", "molecule_set", "stratum", "n_molecules",
                                         "tool", "fraction_carried", "rustle_arm"], hrows,
                            generator="figures/fig_assign_accuracy.py build()", inputs=hinputs, notes=hnotes,
                            data_dir=data_dir)
    p3 = figlib.write_table(TX_TABLE, ["sample", "species", "scope", "tool", "fuzz", "transcripts_total",
                                       "transcripts_in_copy", "in_copy_spanning_ge2", "supported_any",
                                       "supported_hard", "frac_supported_any", "frac_supported_hard", "arm"], txrows,
                            generator="figures/fig_assign_accuracy.py build()", inputs=hinputs, notes=txnotes,
                            data_dir=data_dir)
    out = [p1, p2, p3]
    if not recorded:   # Fig. 5s, simulation: the same runs as a-b (the recorded runs are not joined)
        out += list(build_margin_rule_sim(collected, data_dir))
    if mr_real:
        out.append(figlib.write_table(MR_REAL_TABLE, MR_REAL_HEADER, mr_real,
                                      generator="figures/fig_assign_accuracy.py build()", inputs=mr_inputs,
                                      notes=MR_NOTES[:1] + MR_REAL_NOTES + mr_notes, data_dir=data_dir))
    if fam_rows is not None:
        out.append(figlib.write_table(
            FAMILY_TABLE, ["sample", "species", "family_id", "n_copies", "n_hard", "n_contested", "tool",
                           "carried_hard", "carried_contested"], fam_rows,
            generator="figures/fig_assign_accuracy.py build()", inputs=hinputs,
            notes=[f"one row per (sample, family, method); a hard molecule's family = the family of its result with "
                   f"primary_local = 1 (its primary alignment overlaps a copy of that family), else of its last result; "
                   f"the figure draws families with n_hard >= {HARD_GW_MIN_FAMILY}"] + hnotes[-1:],
            data_dir=data_dir))
    return tuple(out)


# ================================================================ plot
def _ypos():
    """Rows of the dot plot, top to bottom: every MAPQ-0 read, then the identity bands (identical first)."""
    groups = ["all"] + list(figlib.IDENTITY_BANDS)
    return groups, {g: (0.0 if g == "all" else 0.4 + i) for i, g in enumerate(groups)}


SHORT = {"family_certificate": "test in source family*", "per_family_table": "default output",
         "union_certificate": "union test", "aligner_primary": "aligner"}


def _band_panel(ax_cov, ax_acc, rows: list[dict], *, left_labels: bool, bottom_labels: bool):
    """Horizontal dot plot: one row per identity band, x = fraction assigned (left axes) or fraction correct among
    assigned (right axes). Both axes share x (0-1) and y (the bands) across every sample's panel."""
    groups, Y = _ypos()
    offs = dict(zip(_o2.METHODS, [-0.3, -0.1, 0.1, 0.3]))
    first = {r["identity_band"]: r for r in rows if r["method"] == _o2.METHODS[0]}
    for g, r in first.items():
        n, k = int(r["n_reads"]), r.get("n_family_psv", "")
        if n and k not in ("", None):
            ax_cov.barh(Y[g], int(k) / n, height=0.86, color=PSV_BAR, linewidth=0, zorder=1)
    none_assigned, no_ci = [], False
    for m in _o2.METHODS:
        st = METHOD_STYLE[m]
        mk = dict(marker=st["marker"], markersize=2.8, markerfacecolor=st["color"], markeredgecolor=st["color"],
                  markeredgewidth=0.5, linestyle="none", color=st["color"], zorder=3)
        tot_asg = 0
        for r in rows:
            if r["method"] != m or int(r["n_reads"]) == 0:
                continue
            y = Y[r["identity_band"]] + offs[m]
            k = int(r["correct"])
            a = k + int(r["wrong"]) + int(r["conflict"])
            if r["identity_band"] == "all":
                tot_asg = a
            ax_cov.plot([_o2.fnum(r["coverage"])], [y], **mk)
            if a > 0:
                lo, hi = _o2.fnum(r["acc_lo"]), _o2.fnum(r["acc_hi"])
                if lo == lo and hi == hi:
                    ax_acc.plot([lo, hi], [y, y], color=st["color"], linewidth=0.7, zorder=2, solid_capstyle="butt")
                else:
                    no_ci = True
                ax_acc.plot([k / a], [y], **mk)
        if tot_asg == 0:
            none_assigned.append(m)
    for ax in (ax_cov, ax_acc):
        ax.set_yticks([Y[g] for g in groups])
        ax.set_ylim(Y[groups[-1]] + 0.5, -0.5)
        ax.set_xlim(-0.04, 1.04)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"] if bottom_labels else [])
        ax.tick_params(axis="both", labelsize=5.4, length=1.8, pad=1.5)
        ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
        ax.grid(axis="y", visible=False)
        ax.axhline(0.5 * (Y["all"] + Y[groups[1]]), color=figlib.GRID, linewidth=0.8, zorder=0)
        for g in groups[1:]:
            if groups.index(g) % 2:
                ax.axhspan(Y[g] - 0.5, Y[g] + 0.5, color="#f7f6f2", zorder=0, linewidth=0)
    ax_cov.set_yticklabels([BAND_TICK[g] for g in groups] if left_labels else [])
    ax_acc.set_yticklabels([])
    # reads (source copies) per band, right of the accuracy axes
    for g in groups if rows else []:
        r = first.get(g)
        n, nc = (int(r["n_reads"]), int(r["n_copies"])) if r else (0, 0)
        ax_acc.text(1.06, Y[g], f"{n:,} ({nc:,})", transform=ax_acc.get_yaxis_transform(), ha="left", va="center",
                    fontsize=4.8, color=figlib.INK_2)
    if bottom_labels:
        ax_cov.set_xlabel("Fraction\nassigned", labelpad=1, fontsize=5.6)
        ax_acc.set_xlabel("Fraction correct\namong assigned", labelpad=1, fontsize=5.6)
    notes = []
    if none_assigned:
        notes.append("0 assigned: " + ", ".join(SHORT[m] for m in none_assigned))
    if no_ci:
        notes.append(f"no interval: < {_o2.MIN_CI_COPIES} copies assigned")
    return notes


def _hbars(fig, ax, items, *, value_fmt, xlabel, header_x, labels=None, fs=5.6):
    """Horizontal bars grouped under headers. items: [(header, [(tool, value, label_suffix)])]. Tool names are the
    y tick labels (direct labels; FLAIR aqua is below 3:1 on white, so colour never carries identity alone).
    Headers start at figure x `header_x` (the left edge of the tick-label column)."""
    import matplotlib.transforms as mtrans

    y, yt, yl, heads = 0.0, [], [], []
    for header, bars in items:
        heads.append((y, header))
        y += 1.0
        for tool, v, suffix in bars:
            if v == v:  # not NaN
                ax.barh(y, v, height=0.72, **figlib.tool_bar_kwargs(tool))
                ax.text(v + 0.02, y, value_fmt(v) + suffix, ha="left", va="center", fontsize=fs, color=figlib.INK)
            yt.append(y)
            yl.append((labels or {}).get(tool, figlib.TOOL_LABEL[tool]))
            y += 1.0
        y += 0.3
    tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
    for yh, header in heads:
        ax.text(header_x, yh + 0.2, header, transform=tr, ha="left", va="center", fontsize=fs + 0.3, color=figlib.INK,
                fontweight="bold", clip_on=False)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=fs, linespacing=0.95)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=5.4)
    ax.set_ylim(y - 0.3, -0.6)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(xlabel, labelpad=1.5, fontsize=5.8)
    ax.spines["left"].set_visible(False)


def _npip_inset(fig, ax, rows: list[dict], header_x: float):
    tools = [t for t in figlib.TOOL_ORDER if t != "rustle_primary" and any(r["tool"] == t for r in rows)]
    items = []
    for mset, stratum, lab in HARD_STRATA:
        sel = {r["tool"]: r for r in rows if r["molecule_set"] == mset and r["stratum"] == stratum}
        n = int(next(iter(sel.values()))["n_molecules"]) if sel else 0
        items.append((f"{lab}, n = {n:,}", [(t, _o2.fnum(sel[t]["fraction_carried"]) if t in sel else float("nan"), "")
                                           for t in tools]))
    _hbars(fig, ax, items, value_fmt=lambda v: f"{v:.2f}", header_x=header_x, labels=HARD_LABEL,
           xlabel="Hard molecules carried", fs=5.2)


def _gw_panel(ax, sid: str, hard: list[dict], fams: list[dict]):
    """Genome-wide experiment B for one sample: per method, one point per family (>= HARD_GW_MIN_FAMILY hard
    molecules; fraction of its hard molecules carried) and the pooled fraction (black bar), all hard molecules."""
    import random
    tools = ["rustle", "stringtie", "flair", "isoseq"]
    pool = {r["tool"]: r for r in hard if r["molecule_set"] == "all_molecules" and r["stratum"] == "hard (all gate rows)"}
    ax.set_xlim(-0.6, len(tools) - 0.4)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels([GW_TICK.get(t, figlib.TOOL_LABEL[t]) for t in tools], fontsize=5.2, linespacing=0.95)
    ax.tick_params(axis="y", labelsize=5.4)
    ax.set_ylabel("Hard molecules carried", fontsize=5.8, labelpad=1)
    if not pool:
        ax.text(0.5, 0.5, "not run yet", transform=ax.transAxes, ha="center", va="center", fontsize=6,
                color=figlib.INK_3, style="italic")
        return 0, 0
    fam = collections.defaultdict(dict)
    for r in fams:
        fam[r["family_id"]][r["tool"]] = r
    keep = [f for f, v in fam.items() if int(next(iter(v.values()))["n_hard"]) >= HARD_GW_MIN_FAMILY]
    rng = random.Random(20260925)
    for i, t in enumerate(tools):
        ys = [int(fam[f][t]["carried_hard"]) / int(fam[f][t]["n_hard"]) for f in keep if t in fam[f]]
        xs = [i + rng.uniform(-0.22, 0.22) for _ in ys]
        ax.scatter(xs, ys, s=3.5, color=figlib.TOOL_COLOR[t], alpha=0.55, linewidths=0, zorder=2)
        if t in pool:
            v = _o2.fnum(pool[t]["fraction_carried"])
            ax.plot([i - 0.32, i + 0.32], [v, v], color=figlib.INK, linewidth=1.2, zorder=3)
            ax.text(i + 0.34, v, f"{v:.2f}", fontsize=5.0, va="center", ha="left", color=figlib.INK)
    n = int(next(iter(pool.values()))["n_molecules"])
    return n, len(keep)


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    bands = figlib.read_table(BANDS_TABLE, data_dir)
    hard = figlib.read_table(HARD_TABLE, data_dir)
    try:
        fams = figlib.read_table(FAMILY_TABLE, data_dir)
    except FileNotFoundError:
        fams = []
    per = collections.defaultdict(list)
    for r in bands:
        per[_o2.table_sample(r)].append(r)
    fig = plt.figure(figsize=(FIG_W, 7.4))
    # ---- a-b: 2 x 3 small multiples (rows = species groups), shared axes
    x_left, cell_w, cell_gap, pair_gap = 0.078, 0.232, 0.085, 0.012
    ax_w = (cell_w - pair_gap) / 2
    y_top, cell_h, row_gap = 0.845, 0.15, 0.125
    letter = iter("ab")
    for i, row in enumerate(_o2.GRID_2x3):
        y0 = y_top - (i + 1) * cell_h - i * row_gap
        fig.text(0.004, y0 + cell_h + 0.036, next(letter), fontsize=9, fontweight="bold", va="bottom", ha="left")
        for j, s in enumerate(row):
            x0 = x_left + j * (cell_w + cell_gap)
            rows = per.get(s, [])
            ax_cov = fig.add_axes([x0, y0, ax_w, cell_h])
            ax_acc = fig.add_axes([x0 + ax_w + pair_gap, y0, ax_w, cell_h])
            scope = _o2.table_scope(rows[0]) if rows else "genome"
            title = _o2.catalog_title(s, scope, short=True).replace("), ", "),\n", 1)
            n_all = next((int(r["n_reads"]) for r in rows if r["identity_band"] == "all"), 0)
            fig.text(x0, y0 + cell_h + 0.008, title + (f"\n{n_all:,} MAPQ-0 reads" if rows else ""), fontsize=5.7,
                     va="bottom", ha="left", linespacing=1.05)
            if rows:
                notes = _band_panel(ax_cov, ax_acc, rows, left_labels=(j == 0), bottom_labels=True)
                if notes:
                    fig.text(x0, y0 - 0.056, "\n".join(notes), fontsize=4.8, color=figlib.INK_2, va="top",
                             ha="left", linespacing=1.05)
                ax_acc.text(1.06, -0.5, "reads\n(copies)", transform=ax_acc.get_yaxis_transform(), ha="left",
                            va="bottom", fontsize=4.8, color=figlib.INK_2)
            else:
                _band_panel(ax_cov, ax_acc, [], left_labels=(j == 0), bottom_labels=True)
                fig.text(x0 + cell_w / 2, y0 + cell_h / 2, "not run yet", ha="center", va="center", fontsize=6,
                         color=figlib.INK_3, style="italic")
    # legend for a-b
    handles = []
    for m in _o2.METHODS:
        st = METHOD_STYLE[m]
        handles.append(mlines.Line2D([], [], marker=st["marker"], markersize=4.0, linestyle="none",
                                     markerfacecolor=st["color"], markeredgecolor=st["color"], markeredgewidth=0.6,
                                     label=_o2.METHOD_LABEL[m]))
    handles.append(mpatches.Patch(color=PSV_BAR, label=PSV_LABEL))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.03, 1.0), ncol=3, fontsize=5.8,
               handletextpad=0.3, columnspacing=1.0, frameon=False, borderaxespad=0.2)
    fig.text(0.03, 0.935, "Rows: identity of the source copy to its most similar directly aligned copy (100%* = over "
             "the aligned segment, ≥ 50% of the shorter copy).\nLines: 95% interval resampling source copies.",
             fontsize=5.4, color=figlib.INK_2, va="top", ha="left", linespacing=1.1)
    # ---- c: genome-wide per-family points (two samples) + NPIP inset
    y_c0, h_c = 0.105, 0.165
    fig.text(0.004, y_c0 + h_c + 0.045, "c", fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(0.03, y_c0 + h_c + 0.045, "Real reads: hard molecules (≥ 2 equal-best alignments, one in a catalog copy) "
             "whose exact intron chain a method's transcript carries", fontsize=6.4, va="bottom", ha="left")
    for j, s in enumerate(HARD_GW_SAMPLES):
        x0 = 0.085 + j * 0.26
        ax = fig.add_axes([x0, y_c0, 0.2, h_c])
        hs = [r for r in hard if _o2.table_sample(r) == s and r.get("scope") == "genome"]
        n, nf = _gw_panel(ax, s, hs, [r for r in fams if r["sample"] == s])
        t = _o2.catalog_title(s, "genome", short=True).replace("), ", "),\n", 1)
        fig.text(x0, y_c0 + h_c + 0.008, t + (f"\n{n:,} hard molecules; {nf:,} families with ≥ {HARD_GW_MIN_FAMILY}"
                                               if n else ""), fontsize=5.7, va="bottom", ha="left", linespacing=1.05)
    npip = [r for r in hard if r.get("scope", "npip_chr16") == "npip_chr16"]
    x_in = 0.735
    ax_in = fig.add_axes([x_in + 0.1, y_c0, 0.12, h_c])
    _npip_inset(fig, ax_in, npip, header_x=x_in)
    fig.text(x_in, y_c0 + h_c + 0.008, "Inset: NPIP family, human chr16\n(26 copies; development)", fontsize=5.7,
             va="bottom", ha="left", linespacing=1.05)
    import textwrap
    foot = ("* Scored within the read's source family, which is known only in simulation. Rustle's arm in c is the "
            "copy-assignment step's own transcripts (copy_assign --gtf), not the assembler of Figs 1–3; the hard set is "
            "defined by that step's gate. Black bars in c: all hard molecules pooled; points: families.")
    fig.text(0.03, 0.004, "\n".join(textwrap.wrap(foot, 150)), fontsize=5.2, color=figlib.INK_2, va="bottom",
             ha="left", linespacing=1.1)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "fig5_assign_accuracy", out_dir)
    plt.close(fig)
    if (Path(data_dir) / f"{MR_STRATA_TABLE}.tsv").exists():
        paths += plot_margin_rule(data_dir, out_dir)
    return paths


# ================================================================ plot: Fig. 5s (supplement), the margin rule
MR_SEG = [
    ("both_same", "Both assign it, to the same place", "#a9a8a2"),
    ("exceptions", "The margin rule assigns it; Rustle places it elsewhere or leaves it unassigned", "#d9a441"),
    ("aligner_ok", "Rustle only: the aligner's placement (MAPQ > 0), correct", figlib.BLUE[550]),
    ("test_ok", "Rustle only: the copy-assignment test (MAPQ 0), correct", figlib.BLUE[250]),
    ("rustle_bad", "Rustle only: wrong", "#b5452f"),
    ("neither", "Neither assigns it", "#e7e6e1"),
]
MR_SEG_REAL = [
    ("both_same", None, "#a9a8a2"),
    ("exceptions", None, "#d9a441"),
    ("aligner", "Rustle only: the aligner's placement (MAPQ > 0)", figlib.BLUE[550]),
    ("test", "Rustle only: the union test (MAPQ 0)", figlib.BLUE[250]),
    ("neither", None, "#e7e6e1"),
]
MR_BAND_TICK = dict(BAND_TICK, all="All reads")
MR_FOOT = ("Margin rule (Eichler lab): a read is assigned to its best-scoring alignment only if no other alignment of the "
           "read, anywhere in the genome, scores within T alignment-score units (a read with no other alignment is "
           "assigned); otherwise it is discarded. Rustle leaves reads with MAPQ > 0 at the aligner's primary alignment "
           "and applies the copy-assignment test to MAPQ-0 reads (union test: one test per read over all its "
           "alignments). * Scored within the read's source family, known only in simulation.")


def _mr_group(rows, **kw) -> dict:
    """{stratum: row} of the strata rows matching every key=value in kw (values compared as strings)."""
    return {r["stratum"]: r for r in rows if all(r[k] == str(v) for k, v in kw.items())}


def _mr_segments(g: dict, real: bool = False) -> dict:
    n = lambda s: int(g[s]["n_reads"]) if s in g else 0  # noqa: E731
    base = {"both_same": n("both_same"), "exceptions": n("both_differ") + n("margin_only"), "neither": n("neither")}
    if real:
        return dict(base, aligner=n("rustle_aligner"), test=n("rustle_test"))
    rc = lambda s: int(g[s]["rustle_correct"]) if s in g else 0  # noqa: E731
    return dict(base, aligner_ok=rc("rustle_aligner"), test_ok=rc("rustle_test"),
                rustle_bad=n("rustle_aligner") - rc("rustle_aligner") + n("rustle_test") - rc("rustle_test"))


def _mr_bar(ax, y, seg: dict, spec, height=0.66):
    total = sum(seg.values())
    left = 0.0
    for key, _, color in spec:
        v = seg.get(key, 0)
        if v <= 0 or total <= 0:
            continue
        ax.barh(y, v / total, left=left, height=height, color=color, edgecolor=figlib.SURFACE, linewidth=0.3)
        left += v / total
    return total


def _mr_ticks(ax, y, acc: list[dict], sample: str, band: str, label: bool, height=0.66):
    """Where the bar's margin-rule part would end at T = 1 and T = 20 (thin ticks; T = 10 is the bar)."""
    for T in (1, 20):
        r = next((x for x in acc if x["sample"] == sample and x["identity_band"] == band and x["method"] == "margin_rule"
                  and x["threshold"] == str(T)), None)
        if r is None:
            continue
        x = _o2.fnum(r["fraction_assigned"])
        ax.plot([x, x], [y - height / 2 - 0.08, y + height / 2 + 0.08], color=figlib.INK, linewidth=0.6, zorder=4)


def _mr_axis(ax, y_max, xlabel=None):
    ax.set_xlim(0, 1)
    ax.set_ylim(y_max + 0.6, -0.6)
    ax.set_yticks([])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.tick_params(axis="x", labelsize=5.2, length=1.8, pad=1.2)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel, labelpad=1.2, fontsize=5.6)


def _mr_panel_a(fig, strata, acc, top, bottom):
    """Simulation, T = 10, union test: one bar per sample over all its simulated reads, and the numbers."""
    import matplotlib.transforms as mtrans
    ys, y, prev = {}, 0.0, None
    for s in _o2.SAMPLE_ORDER:
        sp = _o2.SAMPLE_SPECIES[s]
        if prev is not None and sp != prev:
            y += 0.45
        ys[s] = y
        y += 1.0
        prev = sp
    x0, w = 0.19, 0.26
    ax = fig.add_axes([x0, bottom, w, top - bottom])
    _mr_axis(ax, y - 1.0, "All simulated reads (%)")
    tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
    cols = [("Reads", 0.49), ("Margin rule assigns\n(fraction correct)", 0.572),
            ("Same place\nby Rustle (share)", 0.658), ("Exceptions (correct:\nrule / Rustle)", 0.744),
            ("Rustle only\n(fraction correct)", 0.83), ("Test in source\nfamily*: correct / n", 0.925)]
    for lab, cx in cols:
        ax.text(cx, -0.85, lab, transform=tr, ha="center", va="bottom", fontsize=5.3, linespacing=1.05)
    ax.text(0.008, -0.85, "Sample (copy catalog)", transform=tr, ha="left", va="bottom", fontsize=5.5)
    prev = None
    for s in _o2.SAMPLE_ORDER:
        yy = ys[s]
        sp = _o2.SAMPLE_SPECIES[s]
        if sp != prev:
            ax.text(0.008, yy - 0.6, _o2.SPECIES_TITLE[sp], transform=tr, ha="left", va="center", fontsize=6.0,
                    fontweight="bold")
        prev = sp
        g = _mr_group(strata, sample=s, identity_band="all", threshold=_o2.MR_HEADLINE_T, reading="u")
        scope = next((r["catalog_scope"] for r in strata if r["sample"] == s), "")
        # the accession of a gorilla development contig is named in panel b's title
        cat = "genome-wide catalog" if scope == "genome" else (f"{scope.split(' (')[0]} catalog only" if scope else "")
        ax.text(0.008, yy, f"{_o2.SAMPLE_LABEL[s]}" + (f"\n{cat}" if cat else ""), transform=tr, ha="left",
                va="center", fontsize=5.5, linespacing=1.0)
        if not g:
            ax.text(0.01, yy, "not run yet", ha="left", va="center", fontsize=5.4, color=figlib.INK_3, style="italic")
            for _, cx in cols:
                ax.text(cx, yy, "–", transform=tr, ha="center", va="center", fontsize=5.5, color=figlib.INK_3)
            continue
        n = _mr_bar(ax, yy, _mr_segments(g), MR_SEG)
        _mr_ticks(ax, yy, acc, s, "all", label=(s == _o2.SAMPLE_ORDER[0]))
        mr = next(r for r in acc if r["sample"] == s and r["identity_band"] == "all" and r["method"] == "margin_rule"
                  and r["threshold"] == str(_o2.MR_HEADLINE_T))
        gi = lambda k, f: int(g[k][f]) if k in g else 0  # noqa: E731
        a_mr = int(mr["assigned"])
        same = gi("both_same", "n_reads")
        ex = gi("both_differ", "n_reads") + gi("margin_only", "n_reads")
        ex_mr = gi("both_differ", "margin_rule_correct") + gi("margin_only", "margin_rule_correct")
        ex_ru = gi("both_differ", "rustle_correct")
        only = gi("rustle_aligner", "n_reads") + gi("rustle_test", "n_reads")
        only_ok = gi("rustle_aligner", "rustle_correct") + gi("rustle_test", "rustle_correct")
        gs = _mr_group(strata, sample=s, identity_band="all", threshold=_o2.MR_HEADLINE_T, reading="s")
        t_n = int(gs["rustle_test"]["n_reads"]) if "rustle_test" in gs else 0
        t_ok = int(gs["rustle_test"]["rustle_correct"]) if "rustle_test" in gs else 0
        vals = [f"{n:,}", f"{a_mr:,}\n({_o2.fnum(mr['fraction_correct']):.4f})",
                f"{same:,}\n({same / a_mr:.4f})" if a_mr else "0", f"{ex:,}\n({ex_mr:,} / {ex_ru:,})",
                f"{only:,}\n({only_ok / only:.4f})" if only else "0", f"{t_ok:,} / {t_n:,}"]
        for (_, cx), v in zip(cols, vals):
            ax.text(cx, yy, v, transform=tr, ha="center", va="center", fontsize=5.5, color=figlib.INK_2,
                    linespacing=1.0)
    return ax


def _mr_panel_b(fig, strata, acc, y_top, cell_h, row_gap):
    """Simulation by identity band of the source copy (T = 10, union test), one small multiple per sample."""
    groups = ["all"] + list(figlib.IDENTITY_BANDS)
    x_left, cell_w, gap = 0.1, 0.25, 0.055
    for i, row in enumerate(_o2.GRID_2x3):
        y0 = y_top - (i + 1) * cell_h - i * row_gap
        for j, s in enumerate(row):
            x0 = x_left + j * (cell_w + gap)
            ax = fig.add_axes([x0, y0, cell_w, cell_h])
            _mr_axis(ax, len(groups) - 1, "Simulated reads of the row (%)" if i == 1 else None)
            ax.set_xticklabels(["0", "25", "50", "75", "100"] if i == 1 else [])
            ax.set_yticks(range(len(groups)))
            ax.set_yticklabels([MR_BAND_TICK[g] for g in groups] if j == 0 else [], fontsize=5.2)
            ax.tick_params(axis="y", length=0, pad=1.5)
            scope = next((r["catalog_scope"] for r in strata if r["sample"] == s), "genome")
            fig.text(x0, y0 + cell_h + 0.005, _o2.catalog_title(s, scope, short=True).replace("), ", "),\n", 1),
                     fontsize=5.5, va="bottom", ha="left", linespacing=1.0)
            any_row = False
            for k, b in enumerate(groups):
                g = _mr_group(strata, sample=s, identity_band=b, threshold=_o2.MR_HEADLINE_T, reading="u")
                if not g:
                    continue
                any_row = True
                n = _mr_bar(ax, k, _mr_segments(g), MR_SEG, height=0.7)
                _mr_ticks(ax, k, acc, s, b, label=False, height=0.7)
                ax.text(1.015, k, f"{n:,}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=4.7,
                        color=figlib.INK_2)
            if not any_row:
                ax.text(0.5, 0.5, "not run yet", transform=ax.transAxes, ha="center", va="center", fontsize=5.8,
                        color=figlib.INK_3, style="italic")


def _mr_panel_c(fig, real, top, bottom):
    """Real reads (no truth), human A119b and gorilla OR6737: the strata at T = 1, 10 and 20."""
    import matplotlib.transforms as mtrans
    rows_y, y = [], 0.0
    for s in HARD_GW_SAMPLES:
        for T in _o2.MR_THRESHOLDS:
            rows_y.append((s, T, y))
            y += 1.0
        y += 0.45
    x0, w = 0.19, 0.26
    ax = fig.add_axes([x0, bottom, w, top - bottom])
    _mr_axis(ax, y - 1.45, "Reads whose primary alignment overlaps a catalog copy (%)")
    tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
    cols = [("Reads\ncompared", 0.49), ("Both assign:\nsame place (share)", 0.585),
            ("Margin rule\ndiscards (share)", 0.69), ("Rustle only:\naligner / union test", 0.8),
            ("Margin rule\nonly", 0.905)]
    for lab, cx in cols:
        ax.text(cx, -0.85, lab, transform=tr, ha="center", va="bottom", fontsize=5.5, linespacing=1.05)
    for s, T, yy in rows_y:
        if T == _o2.MR_THRESHOLDS[0]:
            ax.text(0.008, yy - 0.35, f"{_o2.SAMPLE_LABEL[s]}\ngenome-wide catalog", transform=tr, ha="left",
                    va="top", fontsize=5.5, linespacing=1.0)
        ax.text(x0 - 0.006, yy, f"T = {T}", transform=tr, ha="right", va="center", fontsize=5.3)
        g = _mr_group(real, sample=s, identity_band="all", threshold=T)
        if not g:
            ax.text(0.01, yy, "not run yet", ha="left", va="center", fontsize=5.4, color=figlib.INK_3, style="italic")
            continue
        seg = _mr_segments(g, real=True)
        n = _mr_bar(ax, yy, seg, MR_SEG_REAL)
        gi = lambda k: int(g[k]["n_reads"]) if k in g else 0  # noqa: E731
        both = gi("both_same") + gi("both_differ")
        disc = gi("rustle_aligner") + gi("rustle_test") + gi("neither")
        vals = [f"{n:,}", f"{gi('both_same'):,} ({gi('both_same') / both:.4f})" if both else "0",
                f"{disc / n:.3f}" if n else "–", f"{gi('rustle_aligner'):,} / {gi('rustle_test'):,}",
                f"{gi('margin_only'):,}"]
        for (_, cx), v in zip(cols, vals):
            ax.text(cx, yy, v, transform=tr, ha="center", va="center", fontsize=5.4, color=figlib.INK_2)
    return ax


def plot_margin_rule(data_dir: Path, out_dir: Path):
    """Supplementary Figure 5s (fig5s_margin_rule); rendered only once its simulation tables exist."""
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import textwrap
    strata = figlib.read_table(MR_STRATA_TABLE, data_dir)
    acc = figlib.read_table(MR_ACC_TABLE, data_dir)
    try:
        real = figlib.read_table(MR_REAL_TABLE, data_dir)
    except FileNotFoundError:
        real = []
    fig = plt.figure(figsize=(FIG_W, 7.8))
    handles = [mpatches.Patch(color=c, label=l) for _, l, c in MR_SEG]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1.0), ncol=2, fontsize=5.6,
               handlelength=1.0, handleheight=0.8, handletextpad=0.4, columnspacing=1.2, frameon=False,
               title="Each read, under the margin rule at T = 10 and Rustle (union test for MAPQ-0 reads); "
                     "black ticks: where the margin rule's part ends at T = 1 and T = 20",
               title_fontsize=5.8, alignment="left")
    fig.text(0.004, 0.915, "a", fontsize=9, fontweight="bold", va="bottom")
    fig.text(0.03, 0.915, "Simulated reads (truth = source copy), T = 10", fontsize=6.4, va="bottom")
    _mr_panel_a(fig, strata, acc, top=0.87, bottom=0.665)
    fig.text(0.004, 0.61, "b", fontsize=9, fontweight="bold", va="bottom")
    fig.text(0.03, 0.61, "Simulated reads by identity of the source copy to its most similar directly aligned copy "
             "(100%* over the aligned segment), T = 10", fontsize=6.4, va="bottom")
    _mr_panel_b(fig, strata, acc, y_top=0.565, cell_h=0.09, row_gap=0.06)
    fig.text(0.004, 0.262, "c", fontsize=9, fontweight="bold", va="bottom")
    fig.text(0.03, 0.262, "Real reads, no truth: agreement and coverage only", fontsize=6.4, va="bottom")
    real_legend = [mpatches.Patch(color=c, label=l) for k, l, c in MR_SEG_REAL if l]
    fig.legend(handles=real_legend, loc="lower left", bbox_to_anchor=(0.43, 0.257), ncol=2, fontsize=5.5,
               handlelength=1.0, handleheight=0.8, handletextpad=0.4, columnspacing=1.0, frameon=False,
               title=None)
    fig.text(0.43, 0.25, "Real reads are not split by correctness (no truth); grey, orange and light grey as above.",
             fontsize=5.2, color=figlib.INK_2, va="top", ha="left")
    _mr_panel_c(fig, real, top=0.2, bottom=0.085)
    note = [MR_FOOT]
    scopes = {r["catalog_scope"] for r in strata}
    if scopes - {"genome"}:
        note.append("Development tables (one contig's catalog per sample, reads mapped to the whole genome); the "
                    "genome-wide catalogs of all six samples replace them.")
    fig.text(0.008, 0.002, "\n".join(textwrap.wrap(" ".join(note), 185)), fontsize=5.1, color=figlib.INK_2,
             va="bottom", ha="left", linespacing=1.1)
    tables = [MR_STRATA_TABLE, MR_ACC_TABLE] + ([MR_REAL_TABLE] if real else [])
    figlib.stamp_provisional(fig, tables, data_dir)
    paths = figlib.save(fig, "fig5s_margin_rule", out_dir)
    plt.close(fig)
    return paths


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "margin-rule":
        margin_rule_cli(sys.argv[2:])
    else:
        sys.exit("usage: python3 figures/fig_assign_accuracy.py margin-rule [--inputs F] [--set KEY=VALUE ...]")
