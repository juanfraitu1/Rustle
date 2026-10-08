"""_o1 — private helpers for Figure 6 (the default de novo family definition across the paralogue identity spectrum)
and its supplements, and the genome-wide helpers Figure 7 shares (samples, exposure, substrates, annotation caches,
reference families, sharded all-vs-alls).

THE FAMILIES SCORED (user decision 2026-09-25 16:00; docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md, Amendment 1): the
ONE default de novo family definition = the driver's `families` stage: reads -> seeded assembly loci -> one
representative per locus (its "positional exon sum": read-derived exon coordinates, genome bases) -> mcl_families
--from-gtf (exon-sum >= 0.60, MCL 2.8). Its copy table `<id>.fam.copies.tsv` (one copy per member locus; the contract
copy assignment consumes) is what score.py pairs reads. The gw_family_catalog copy catalog is LEGACY and not scored
here. Protein is not part of the default: the translated protein tier (T3) and the protein-homology families appear
only in supplementary figures.

Two scopes. `fig6_scope` in the inputs (default `genome`):

GENOME (every sample, one genome-wide run per object)
    main, panels a-c   human samples (human_A119b, human_testis; Ensembl Compara is human only)
      pairs     Ensembl Compara release 116 human paralogue pairs across chromosomes (`bench/truth.py compara`;
                inputs key `compara_gw`)
      loci      the sample's genome-wide assembly (run-cache stage `assemble`, via assembly.rustle_product: a figure
                build never assembles)
      edges     `bench/score.py spectrum --chrom ALL --skip-t3 --minimap2 tools/mm2_shard.sh` (T1 asm20 / T2 asm20
                -k11 -w5, nucleotide; the per-pair `recovered` flag is T1 or T2)
      families  the default families' copy table (run-cache stage `families`, product `copies`), scored by
                `score.py pairs --chrom ALL --universe <spectrum truth_pairs>`
    supplement fig6s_seeding   every sample, the default families of both seeding configurations (stages `families`
                and `families_primary`; the second is optional)
      Compara   (human samples) Compara pairs at every band, universe = both genes with >= 2 reads whose PRIMARY
                alignment has an aligned block on the gene's annotated exons (`primary_counts_gw`, every gene and
                pseudogene record; `compara_universe`); headline band >= 90%
      Liftoff   (every sample) the Fig. 8 self-lift's (record, extra copy) pairs with both loci read-supported
                (_liftoff.copy_pairs / pair_families; the C2 support table)
      protein-homology families: secondary, only with `fig6s_protein_homology 1` (built by `bench/truth.py
                protein-homology --chrom ALL`; the annotated-mRNA bands of `mrna_pairs_paf`)
    supplement fig6s_protein   the chr16 development spectrum WITH the translated protein tier (T3): a comparator,
                not part of Rustle's rule (never run genome-wide)

DEV (development tables: human A119b chr16; gorilla OR6737 NC_073244.2 for the seeding supplement; `fig6_scope dev`)
    the chr16 restriction of the genome-wide assembly -> score.py spectrum --skip-t3 (main) and with T3 (supplement);
    the driver's families stage on the same GTF (a mcl_families that writes the copy table: `families_bin` in the
    inputs overrides `bin` for this one call); the gorilla contig's protein-homology families for the dev seeding
    supplement (a secondary reference, labelled so).

Heavy units (sharded all-vs-alls, BAM passes) are resumable: a build call stops with exit 75 (assembly.Pending) when
its budget (`fig6_budget_s` / `figs_budget_s`, default 540 s for the genome scope) is used up; the same command
continues. Run-cache products that are missing or stale are never built here: the build names the `make.py runs`
command and leaves the sample out (listed in the table notes; the tables are then marked provisional).
"""
from __future__ import annotations

import ast
import collections
import hashlib
import itertools
import json
import math
import os
import re
import sys
import time
from pathlib import Path

import assembly
import figlib

FIG = "fig6"
HUMAN_CHROM = "chr16"
GORILLA_CONTIG_DEFAULT = "NC_073244.2"   # override with the inputs key `fig6_gorilla_contig` (development contig)
EXPRESSED_MIN_PRIMARY = 2               # a gene is in the sensitivity universe with >= 2 primary reads on its exons
PRIMARY_COUNT_VERSION = "1"             # bump when primary_counts changes: the cached counts are recomputed
COMPARA_BANDS = [">=90", "80-90", "70-80", "60-70", "50-60", "30-50", "<30"]
POOLED = ("60-90", ["80-90", "70-80", "60-70"])   # the claim's pooled range
MRNA_BANDS = [">=90", "80-90", "70-80", "60-70", "<60", "none"]
# score.py spectrum tier unions -> the view names the tables use
SPECTRUM_VIEWS = {"T1": "edge_asm20", "T1+T2": "edge_nt", "T1+T2+T3": "edge_nt_protein"}
# score.py spectrum precision rows -> measure names
SPECTRUM_PRECISION = {"precision T1 (nt)": "edge_asm20", "precision T2 (nt)": "edge_nt_sensitive",
                      "precision T3 (protein)": "edge_protein"}
# the recorded runs (registers 1096-1101); configuration -> recorded clusters in cfg['referee_dir']
RECORDED_GORILLA = {"rustle": "ggo44_GOOD0.98.fam.clusters.tsv", "rustle_primary": "ggo44_P.fam.clusters.tsv"}
RECORDED_CATALOG = "c16_shipped"


def scope(cfg: dict) -> str:
    s = str(cfg.get("fig6_scope", "genome")).strip().lower()
    if s not in ("genome", "dev"):
        raise ValueError(f"fig6_scope {s!r}: expected genome or dev")
    return s


def gorilla_contig(cfg: dict) -> str:
    """The dev-scope gorilla contig: the seeding default's decision contig (development evidence, not held out)."""
    return cfg.get("fig6_gorilla_contig", GORILLA_CONTIG_DEFAULT)


def score_py(cfg: dict) -> str:
    return str(Path(cfg["repo"]) / "bench" / "score.py")


def _bench(cfg: dict):
    """bench/score.py and bench/lib.py as modules (their loaders; imported lazily, plot() never needs them)."""
    b = str(Path(cfg["repo"]) / "bench")
    if b not in sys.path:
        sys.path.insert(0, b)
    import lib  # noqa: E402
    import score  # noqa: E402
    return score, lib


def referee_paths(cfg: dict) -> dict:
    """Dev-scope panel d reference on the gorilla contig (recorded, off-repo products; captions/fig6.md): the
    contig's protein-homology families, gene records, annotated-mRNA PAF, and the RECORDED span-rule gene list of
    registers 1096-1101 (`expressed_recorded`; only checked against primary_counts' span rule)."""
    d = Path(cfg["referee_dir"]) / "ref"
    c = gorilla_contig(cfg)
    return {"referee": d / f"{c}.tsv", "genes": d / f"{c}.genes.gff", "expressed_recorded": d / f"{c}.expressed.tsv",
            "mrna_paf": d / "referee_mrna.paf"}


def _gene_exons(gff: Path, chrom: str) -> tuple[dict, dict]:
    """({gene Name: (start0, end)} as score.load_gene_spans reads it (gene/pseudogene `Name=`, last record wins),
    {gene Name: [exon (start0, end)]}: every exon record whose Parent chain reaches THAT gene record)."""
    name_id: dict = {}
    parent: dict = {}
    ftype: dict = {}
    exons = collections.defaultdict(list)
    spans: dict = {}
    for ln in open(gff):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != chrom:
            continue
        mid = re.search(r"(?:^|;)ID=([^;]+)", f[8])
        mp = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
        if f[2] in ("gene", "pseudogene"):
            mn = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if mn and mid:
                name_id[mn.group(1)] = mid.group(1)
                spans[mn.group(1)] = (int(f[3]) - 1, int(f[4]))
        if mid:
            ftype[mid.group(1)] = f[2]
            if mp:
                parent[mid.group(1)] = mp.group(1).split(",")[0]
        if f[2] == "exon" and mp:
            exons[mp.group(1).split(",")[0]].append((int(f[3]) - 1, int(f[4])))

    def gene_of(x):
        for _ in range(8):
            if x is None or ftype.get(x) in ("gene", "pseudogene"):
                break
            x = parent.get(x)
        return x if ftype.get(x) in ("gene", "pseudogene") else None
    by_gene = collections.defaultdict(list)
    for p, ex in exons.items():
        g = gene_of(p)
        if g:
            by_gene[g].extend(ex)
    return spans, {n: sorted(by_gene.get(i, [])) for n, i in name_id.items()}


def _count_primary(bam, contig: str, s: int, e: int, ivs) -> tuple[int, int]:
    """(reads whose PRIMARY record has an aligned block on one of `ivs`, primary records overlapping [s, e)):
    the rule of primary_counts (pysam get_blocks = M/=/X runs; N and D are not aligned sequence)."""
    n_span, on_exon = 0, set()
    for rec in bam.fetch(contig, s, e):
        if rec.flag & 0x904:
            continue
        n_span += 1
        if ivs and any(bs < xe and be > xs for bs, be in rec.get_blocks() for xs, xe in ivs):
            on_exon.add(rec.query_name)
    return len(on_exon), n_span


def primary_counts(cfg: dict, force: bool = False) -> Path:
    """Dev scope. Per protein-homology gene of the gorilla contig: primary reads on its EXONS and primary records on
    its SPAN.

    `n_primary_exon` = reads whose PRIMARY record (not 0x100/0x800/0x4) has an aligned block (pysam get_blocks:
    M/=/X runs; an intron N or a deletion D is not aligned sequence) overlapping one of the gene's annotated exons
    (every exon record under the gene's `Name=` record in genes.gff) — fig. 3's n_mol_primary rule at gene level.
    `n_primary_span` = primary records overlapping the gene span (`samtools view -c -F 0x904 BAM chrom:span`), the rule
    of the recorded list (spliced-over reads included). No strand condition in either. Light: one region fetch per
    gene (NC_073244.2: 473 k records, ~20 s). Cached under ${work}/fig6/, keyed on the BAM, genes.gff, the reference
    families and PRIMARY_COUNT_VERSION."""
    import pysam

    r = referee_paths(cfg)
    c = gorilla_contig(cfg)
    out = figlib.work_dir(cfg, FIG) / f"{c}.referee_primary_counts.tsv"
    bam_path = cfg["gorilla_bam"]
    if not force and figlib.fresh(out, bam_path, r["genes"], r["referee"]):
        with open(out) as fh:
            if fh.readline().rstrip("\n") == f"#version\t{PRIMARY_COUNT_VERSION}":
                return out
    _, lib = _bench(cfg)
    fam = lib.read_referee(str(r["referee"]))
    spans, gex = _gene_exons(r["genes"], c)
    bam = pysam.AlignmentFile(bam_path, "rb")
    tmp = out.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write(f"#version\t{PRIMARY_COUNT_VERSION}\n")
        fo.write("#rule\tn_primary_exon = primary reads with an aligned block on an exon; n_primary_span = primary "
                 "records overlapping the gene span (-F 0x904)\n")
        fo.write("gene\treferee_family\tstart\tend\tn_exon_records\tn_primary_exon\tn_primary_span\n")
        for g in sorted(fam):
            if g not in spans:
                fo.write(f"{g}\t{fam[g]}\t\t\t0\t0\t0\n")
                continue
            s, e = spans[g]
            ivs = gex.get(g, [])
            n_exon, n_span = _count_primary(bam, c, s, e, ivs)
            fo.write(f"{g}\t{fam[g]}\t{s + 1}\t{e}\t{len(ivs)}\t{n_exon}\t{n_span}\n")
    tmp.replace(out)
    return out


def read_primary_counts(path: Path) -> dict:
    """{gene: {'family', 'exon', 'span'}} from primary_counts' file."""
    out = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("gene\t"):
                continue
            f = ln.rstrip("\n").split("\t")
            out[f[0]] = {"family": f[1], "exon": int(f[5]), "span": int(f[6])}
    return out


def expressed_list(cfg: dict, rule: str, force: bool = False) -> Path:
    """Dev scope: the sensitivity universe's genes under `rule` ('exon' = the figure's rule, 'span' = the recorded
    list's rule): genes with >= EXPRESSED_MIN_PRIMARY primary reads (exon) or records (span). Written in score.py's
    --expressed format (first column; header starting 'Gene')."""
    counts = primary_counts(cfg, force)
    out = counts.with_name(f"{gorilla_contig(cfg)}.expressed.{rule}.tsv")
    if not force and figlib.fresh(out, counts):
        return out
    pc = read_primary_counts(counts)
    with open(out, "w") as fo:
        fo.write(f"Gene Name\tn_primary_{rule}\n")
        for g in sorted(pc):
            if pc[g][rule] >= EXPRESSED_MIN_PRIMARY:
                fo.write(f"{g}\t{pc[g][rule]}\n")
    return out


def chr16_ref_gtf(cfg: dict, force: bool = False) -> Path:
    """The RefSeq annotation restricted to chr16 (gene spans for score.py pairs, symbols for score.py spectrum)."""
    return assembly.restrict_gtf(cfg["human_ref_gtf"], figlib.work_dir(cfg, FIG) / f"{HUMAN_CHROM}.ref.gtf",
                                 {HUMAN_CHROM}, force=force)


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float | None, float | None]:
    """Wilson score 95% interval of k/n."""
    if not n:
        return None, None
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, mid - half), min(1.0, mid + half)


# ---------------------------------------------------------------- dev-scope sources
def recorded_sources(cfg: dict) -> dict:
    """The recorded development runs (registers 1096-1101): the chr16 spectrum WITH T3 (supplement fig6s_protein)
    and the gorilla seeding configurations (supplement fig6s_seeding, development). No recorded run of the default
    families' copy table exists, so a recorded build writes the supplementary tables only."""
    sd, rd = Path(cfg["spectrum_dir"]), Path(cfg["referee_dir"])
    return {
        "source": ("recorded runs: spectrum 2026-09-24 (score.py spectrum's predecessor identity_spectrum.py on the "
                   "2026-09-22 primaries-only assemble-only chr16 GTF), gorilla configurations ggo44_P / "
                   "ggo44_GOOD0.98 2026-09-22/23 (chain.sh region runs on the contig slice)"),
        "spectrum_t3": sd / HUMAN_CHROM,
        "clusters": {arm: rd / f for arm, f in RECORDED_GORILLA.items()},
    }


def families_bin(cfg: dict) -> str:
    """The bin directory of the dev-scope families run: `families_bin` in the inputs (a mcl_families that writes the
    copy table, until `bin` has one), else `bin`."""
    return str(cfg.get("families_bin") or cfg["bin"])


def ensure_sources(cfg: dict, force: bool = False, need_main: bool = True) -> dict:
    """Dev scope: regenerate every prediction with the current binaries (HEAVY; cached under ${work}/fig6/).

    1. the human genome-wide assembly from the run cache (assembly.rustle_product: never assembled here; driver
       default = loci seeded with secondary alignments within 2% of the read's best score) restricted to chr16
       -> `score.py spectrum --skip-t3` (main: nucleotide tiers T1/T2, minimap2 x2) and `score.py spectrum` with the
       translated tier T3 (supplement fig6s_protein; mmseqs)
    2. the DEFAULT de novo families on the same chr16 GTF: `rustle_pipeline.sh families` (mcl_families --from-gtf
       --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units -> chr16_denovo.fam.copies.tsv; ~1-3 min)
    3. gorilla genome-wide assemblies, both seeding configurations (run cache), restricted to the dev contig
       -> `rustle_pipeline.sh families` each (supplement fig6s_seeding, development; ~1 min each)
    `force` recomputes steps 1-3 but never the genome-wide assemblies (`make.py runs` owns those)."""
    w = figlib.work_dir(cfg, FIG)
    threads = cfg.get("threads", "4")
    py = sys.executable or "python3"
    ref16 = chr16_ref_gtf(cfg, force=force)
    out: dict = {"source": "regenerated with the current binaries and pipeline-driver defaults", "clusters": {},
                 "arm_gtf": {}}

    # 1. identity spectrum on the chr16 loci of the default assembly: nucleotide only (main) and with T3 (supplement)
    genome_gtf = assembly.rustle_product(cfg, "human", "rustle")
    gtf = assembly.restrict_gtf(genome_gtf, w / f"{HUMAN_CHROM}.assembled.gtf", {HUMAN_CHROM}, force=force)
    (w / "spectrum").mkdir(exist_ok=True)
    # not keyed on bench/score.py: it is a multi-command file other stages edit (an mtime key would re-run this
    # 2-min step on unrelated edits); the tables fingerprint score.py and lib.py instead
    for tag, extra in (("", ["--mmseqs", cfg.get("mmseqs", "mmseqs")]), ("_nt", ["--skip-t3"])):
        prefix = w / "spectrum" / f"{HUMAN_CHROM}{tag}"
        spec = Path(f"{prefix}.spectrum.tsv")
        if force or not figlib.fresh(spec, gtf, ref16, cfg["compara_chr16"], cfg["human_fasta"]):
            figlib.run([py, score_py(cfg), "spectrum", "--gtf", str(gtf), "--ref", str(ref16),
                        "--fasta", cfg["human_fasta"], "--chrom", HUMAN_CHROM, "--compara", cfg["compara_chr16"],
                        "--out", str(prefix), "--threads", threads, *extra], log=w / f"spectrum{tag}.log")
        out["spectrum" if tag else "spectrum_t3"] = prefix

    # 2. the default de novo families on the same GTF (the driver's families stage writes their copy table)
    if need_main:
        fprefix = w / f"{HUMAN_CHROM}_denovo"
        fgtf = Path(f"{fprefix}.gtf")          # the driver reads PREFIX.gtf: a link to the chr16 restriction
        if not (fgtf.is_symlink() and fgtf.resolve() == gtf.resolve()):
            fgtf.unlink(missing_ok=True)
            fgtf.symlink_to(gtf.resolve())
        copies = Path(f"{fprefix}.fam.copies.tsv")
        fbin = families_bin(cfg)
        if force or not (figlib.fresh(copies, fgtf, Path(fbin) / "mcl_families", cfg["human_fasta"])
                         and _driver_current(cfg, copies)):
            figlib.run(["bash", cfg["driver"], "families", "--bam", cfg["human_chr16_bam"], "--fasta",
                        cfg["human_fasta"], "--out", str(fprefix), "--bin", fbin, "--threads", threads],
                       log=w / f"{HUMAN_CHROM}_denovo.families.driver.log")
            if not copies.exists():
                raise NotBuilt(f"{fbin}/mcl_families does not write the families copy table ({copies}): rebuild it, or "
                               "pass --set families_bin=<a bin dir whose mcl_families --help names <out>.copies.tsv>")
            _driver_stamp(cfg, copies)
        out["families"] = copies
        out["families_bin"] = Path(fbin) / "mcl_families"

    # 3. gorilla: the dev contig (the seeding decision's contig), both seeding configurations (supplement)
    for arm in ("rustle", "rustle_primary"):
        g = assembly.rustle_product(cfg, "gorilla", arm)
        gprefix = w / f"gorilla_{arm}"
        garm = assembly.restrict_gtf(g, Path(f"{gprefix}.gtf"), {gorilla_contig(cfg)}, force=force)
        out["arm_gtf"][arm] = garm
        clusters = Path(f"{gprefix}.fam.clusters.tsv")
        if force or not (figlib.fresh(clusters, garm, Path(cfg["bin"]) / "mcl_families")
                         and _driver_current(cfg, clusters)):
            figlib.run(["bash", cfg["driver"], "families", "--bam", cfg["gorilla_bam"], "--fasta",
                        cfg["gorilla_fasta"], "--out", str(gprefix), "--bin", cfg["bin"], "--threads", threads],
                       log=w / f"gorilla_{arm}.families.driver.log")
            _driver_stamp(cfg, clusters)
        out["clusters"][arm] = clusters
    return out


def _driver_current(cfg: dict, product: Path) -> bool:
    """The driver's CODE (non-comment lines; assembly._code_hash, as for the assemblies) is the one that made
    `product` (its flags, e.g. the families stage's, live there). A product without a stamp is adopted once."""
    stamp = Path(str(product) + ".driver_code")
    code = assembly._code_hash(cfg["driver"])
    if Path(product).exists() and not stamp.exists():
        stamp.write_text(code + "\n")
    return stamp.exists() and stamp.read_text().strip() == code


def _driver_stamp(cfg: dict, product: Path) -> None:
    Path(str(product) + ".driver_code").write_text(assembly._code_hash(cfg["driver"]) + "\n")


# ---------------------------------------------------------------- repo scorers (seconds each)
def score_members(cfg: dict, members: Path, universe: Path, genes_gtf: Path, log: Path, chrom: str = HUMAN_CHROM,
                  compara: str | None = None) -> dict:
    """`score.py pairs` Compara report on a copy table (the default families' `<id>.fam.copies.tsv`, or any file in
    the gw_family_catalog copies contract) over a fixed universe -> parsed dict. chrom 'ALL' = genome-wide."""
    figlib.run([sys.executable or "python3", score_py(cfg), "pairs", "--members", str(members), "--genes",
                str(genes_gtf), "--chrom", chrom, "--truth", f"compara:{compara or cfg['compara_chr16']}",
                "--universe", str(universe)], log=log)
    return parse_pairs_catalog(log)


def score_referee(cfg: dict, clusters: Path, label: str, log: Path, expressed: Path, *, genes=None, truth=None,
                  paf=None, chrom=None) -> dict:
    """`score.py pairs` mRNA-band report -> parsed dict. Defaults: the dev-scope gorilla contig's files."""
    r = referee_paths(cfg)
    figlib.run([sys.executable or "python3", score_py(cfg), "pairs", "--members", str(clusters), "--genes",
                str(genes or r["genes"]), "--chrom", chrom or gorilla_contig(cfg), "--truth",
                f"families:{truth or r['referee']}", "--expressed", str(expressed), "--bands",
                f"paf:{paf or r['mrna_paf']}", "--label", label], log=log)
    return parse_pairs_referee(log)


_CAT_HEAD = re.compile(r"^\[catalog\] .*: (\d+) copies on (\S+), (\d+) families, largest (\d+) genes, (\d+) gene pairs; "
                       r"truth pairs (\d+), with both genes in the catalog (\d+)")
_FRAC = re.compile(r"(\d+)/(\d+)")
_UNI = re.compile(r"recall over the fixed universe \((\d+) truth pairs")


def parse_pairs_catalog(log: Path) -> dict:
    out: dict = {"recall_in_catalog": {}, "recall_universe": {}}
    block = "in"
    for line in open(log):
        s = line.strip()
        m = _CAT_HEAD.match(s)
        if m:
            (out["copies"], out["chrom"], out["families"], out["largest"], out["gene_pairs"], out["truth_pairs"],
             out["truth_in_catalog"]) = (int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4)),
                                         int(m.group(5)), int(m.group(6)), int(m.group(7)))
            continue
        m = _UNI.search(s)
        if m:
            out["universe_pairs"] = int(m.group(1))
            block = "universe"
            continue
        if s.startswith("precision (judgeable"):
            k, n = map(int, _FRAC.search(s).groups())
            out["precision_all"] = (k, n)
        elif s.startswith("precision, multi-exon"):
            k, n = map(int, _FRAC.search(s).groups())
            out["precision_multi_exon"] = (k, n)
        elif s.startswith("recall ") and block == "in":
            band, frac = s.split()[1], s.split()[2]
            out["recall_in_catalog"][band] = tuple(map(int, frac.split("/")))
        elif block == "universe" and s and s.split()[0] in COMPARA_BANDS:
            band, frac = s.split()[0], s.split()[1]
            out["recall_universe"][band] = tuple(map(int, frac.split("/")))
    for key in ("copies", "precision_all", "precision_multi_exon", "universe_pairs"):
        if key not in out:
            raise RuntimeError(f"{log}: could not parse '{key}' from score.py pairs output")
    return out


_REF = re.compile(r"largest\s+(\d+) genes (\{.*?\}) \| prec (\d+)/(\d+)=\S+ \| recall by annotated-mRNA identity: (.*)$")


def parse_pairs_referee(log: Path) -> dict:
    for line in open(log):
        m = _REF.search(line.rstrip("\n"))
        if not m:
            continue
        recall = {}
        for part in m.group(5).split(" · "):
            if not part.strip():
                continue
            band, frac = part.split()
            recall[band] = tuple(map(int, frac.split("/")))
        return {"largest": int(m.group(1)), "largest_composition": ast.literal_eval(m.group(2)),
                "precision": (int(m.group(3)), int(m.group(4))), "recall": recall}
    raise RuntimeError(f"{log}: no mRNA-band line in score.py pairs output")


def parse_spectrum(prefix: Path) -> dict:
    """PREFIX.spectrum.tsv -> {'recall': {band: (n, {tier: fraction})}, 'precision': [(measure, band, n, frac)]}."""
    rec, prec = {}, []
    with open(f"{prefix}.spectrum.tsv") as fh:
        next(fh)
        for line in fh:
            metric, band, n, value = line.rstrip("\n").split("\t")
            if metric == "recall":
                rec[band] = (int(n), ast.literal_eval(value))
            elif metric in SPECTRUM_PRECISION:
                prec.append((SPECTRUM_PRECISION[metric], band, int(n), float(value)))
    return {"recall": rec, "precision": prec}


def widest_tier(spec: dict) -> str:
    """The widest tier union the spectrum run wrote: T1+T2+T3, or T1+T2 under --skip-t3."""
    tiers = set().union(*(set(fr) for _, fr in spec["recall"].values())) if spec["recall"] else set()
    return "T1+T2+T3" if "T1+T2+T3" in tiers else "T1+T2"


def universe_pairs(prefix: Path) -> int:
    with open(f"{prefix}.truth_pairs.tsv") as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith("geneA\t"))


def copy_exon_profile(copies: Path, chrom: str | None) -> tuple[int, int]:
    """(copies on chrom (every contig when None), single-exon copies) of a copy table (families stage or legacy
    catalog; both have `chrom` and `n_exon` columns)."""
    n = single = 0
    with open(copies) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci, ei = header.index("chrom"), header.index("n_exon")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ei and (chrom is None or f[ci] == chrom):
                n += 1
                single += f[ei] == "1"
    return n, single


# ---------------------------------------------------------------- cluster-level attribution
def _pair_groups(groups: dict) -> dict:
    """{frozenset(gene pair): {group ids whose members contain both genes}} (a gene may sit in >1 group)."""
    out = collections.defaultdict(set)
    for c, gs in groups.items():
        for p in itertools.combinations(sorted(gs), 2):
            out[frozenset(p)].add(c)
    return out


def _loo(pairs, p2g, drop) -> int:
    """Pairs still present when group `drop` is removed (a pair survives if another group also holds it)."""
    return sum(1 for k in pairs if p2g[k] - {drop})


def _top(counter: collections.Counter):
    return max(counter, key=lambda c: (counter[c], str(c))) if counter else None


def _components(pairs) -> dict:
    """gene -> component root (min gene name) of the graph whose edges are `pairs` (order-independent)."""
    par: dict = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for k in sorted(tuple(sorted(k)) for k in pairs):
        a, b = find(k[0]), find(k[1])
        if a != b:
            par[max(a, b)] = min(a, b)
    return {g: find(g) for g in par}


def _group_label(genes) -> str:
    """A readable name for a gene group: 'NPIP* (14 genes)' when every gene shares a >= 3-letter name prefix;
    'NPIP (13 of 22 genes NPIP*)' when the most common 4-letter prefix of the named (non-LOC) genes covers at least
    half of them; else '<first gene> + n genes'. The text before ' (' is the short name the plot uses."""
    gs = sorted(g.split(":", 1)[1] if ":" in g else g for g in genes)
    pre = re.match(r"[A-Za-z]*", os.path.commonprefix(gs)).group(0)
    if len(pre) >= 3:
        return f"{pre}* ({len(gs)} genes)"
    named = [g for g in gs if not g.startswith("LOC")]
    if named:
        p4, k = collections.Counter(g[:4] for g in named if len(g) >= 4 and g[:4].isalpha()).most_common(1)[0] \
            if any(len(g) >= 4 and g[:4].isalpha() for g in named) else ("", 0)
        if p4 and 2 * k >= len(named):
            k_all = sum(1 for g in gs if g.startswith(p4))
            return f"{p4} ({k_all} of {len(gs)} genes {p4}*)"
    return f"{gs[0]} + {len(gs) - 1} genes"


def referee_attribution(cfg: dict, clusters: Path, scored: dict, expressed: Path, *, genes=None, truth=None,
                        paf=None, chrom=None) -> dict:
    """Which predicted families and reference families the configuration's recovered / judged pairs come from.

    Reproduces `score.py pairs --bands paf:` with the scorer's own loaders and checks every band total and the
    precision against `scored` (the parsed scorer output); returns per band {k, n, clusters, fams, fams_in_band,
    top, top_pairs, loo_k} and 'precision' {k, n, locus_pairs, gene_pairs, clusters, fams, top, top_pairs, loo_k,
    loo_n}, plus 'cluster_info' {cid: (loci, span label, composition)}. chrom 'ALL' = genome-wide keys
    (CONTIG:NAME); the defaults are the dev-scope gorilla contig's files."""
    score, lib = _bench(cfg)
    r = referee_paths(cfg)
    genes, truth, paf = genes or r["genes"], truth or r["referee"], paf or r["mrna_paf"]
    c = chrom or gorilla_contig(cfg)
    everywhere = c == "ALL"
    if everywhere:
        gene_of = score.span_mapper_all(score.load_gene_spans_all(str(genes), "gff"), lib.genome_gene_key)
        fam = lib.read_families(str(truth))[0]
        expr = score.read_expressed(str(expressed))
    else:
        gene_of = score.span_mapper(score.load_gene_spans(str(genes), c, "gff"))
        fam = lib.read_referee(str(truth))
        expr = set(ln.split("\t")[0] for ln in open(expressed) if not ln.startswith("Gene"))
    ident: dict = {}
    for ln in open(paf):
        f = ln.split("\t")
        if f[0] == f[5]:
            continue
        k = frozenset((f[0], f[5]))
        idn = int(f[9]) / max(1, int(f[10]))
        if k not in ident or idn > ident[k]:
            ident[k] = idn

    def band(k):
        if k not in ident:
            return "none"
        for lo, hi, name in score.BANDS_MRNA:
            if lo <= ident[k] < hi:
                return name
    truth_pairs = {k: band(k) for k in score.family_pairs(fam, expr)}
    cg, _, _ = score.load_members(str(clusters), c, gene_of)
    p2c = _pair_groups(cg)
    out: dict = {}
    for b in MRNA_BANDS:
        tb = [k for k, bb in truth_pairs.items() if bb == b]
        if not tb:
            continue
        hit = [k for k in tb if k in p2c]
        byc = collections.Counter(x for k in hit for x in p2c[k])
        fams = {fam[next(iter(k))] for k in hit}
        top = _top(byc)
        out[b] = {"k": len(hit), "n": len(tb), "clusters": len(byc), "fams": len(fams),
                  "fams_in_band": len({fam[next(iter(k))] for k in tb}), "top": top,
                  "top_pairs": byc[top] if top else 0, "loo_k": _loo(hit, p2c, top) if top else 0}
        if (out[b]["k"], out[b]["n"]) != tuple(scored["recall"].get(b, (None, None))):
            raise RuntimeError(f"{clusters}: band {b} attribution {out[b]['k']}/{out[b]['n']} != scorer "
                               f"{scored['recall'].get(b)}")
    judged = {k for k in p2c if all(g in fam for g in k)}
    tp = {k for k in judged if len({fam[g] for g in k}) == 1}
    byc = collections.Counter(x for k in judged for x in p2c[k])
    top = _top(byc)
    locus_rows = collections.Counter()
    spans: dict = {}
    for ln in open(clusters):
        f = ln.rstrip("\n").split("\t")
        if f[0] == "cluster_id":
            continue
        locus_rows[f[0]] += 1
        s, e = int(f[6]), int(f[7])
        ctg, lo, hi = spans.get(f[0], (f[5], s, e))
        spans[f[0]] = (ctg if ctg == f[5] else "*", min(lo, s), max(hi, e))
    out["precision"] = {"k": len(tp), "n": len(judged), "locus_pairs": sum(v * (v - 1) // 2 for v in locus_rows.values()),
                        "gene_pairs": len(p2c), "clusters": len(byc), "fams": len({fam[next(iter(k))] for k in tp}),
                        "top": top, "top_pairs": byc[top] if top else 0,
                        "loo_k": _loo(tp, p2c, top) if top else 0, "loo_n": _loo(judged, p2c, top) if top else 0}
    if (out["precision"]["k"], out["precision"]["n"]) != tuple(scored["precision"]):
        raise RuntimeError(f"{clusters}: precision attribution {out['precision']['k']}/{out['precision']['n']} != "
                           f"scorer {scored['precision']}")

    def span_label(cid):
        ctg, lo, hi = spans[cid]
        return f"{ctg}:{lo + 1}-{hi}" if ctg != "*" else "several contigs"
    out["cluster_info"] = {cid: (locus_rows[cid], span_label(cid),
                                 dict(collections.Counter(fam.get(g, score.NO_FAMILY) for g in sorted(gs))))
                           for cid, gs in cg.items()}
    out["n_clusters"] = len(locus_rows)
    return out


def gorilla_attribution(cfg: dict, clusters: Path, scored: dict, expressed: Path) -> dict:
    """Dev scope: referee_attribution on the gorilla contig (cluster_info spans as (loci, start0, end, composition)
    for the old table columns)."""
    at = referee_attribution(cfg, clusters, scored, expressed)
    info = {}
    for cid, (loci, lab, comp) in at["cluster_info"].items():
        m = re.match(r".*:(\d+)-(\d+)$", lab)
        info[cid] = (loci, int(m.group(1)) - 1 if m else 0, int(m.group(2)) if m else 0, comp)
    at["cluster_info"] = info
    return at


def clusters_rows(clusters: Path) -> int:
    """Data rows of a clusters.tsv (its header starts with cluster_id)."""
    with open(clusters) as fh:
        return sum(1 for ln in fh if ln.strip() and not ln.startswith("cluster_id"))


def drop_cluster(clusters: Path, drop: str, dst: Path) -> Path:
    """The clusters file without the rows of cluster `drop` (for the leave-largest-family-out rescoring)."""
    with open(clusters) as fi, open(dst, "w") as fo:
        for ln in fi:
            if ln.split("\t", 1)[0] != drop:
                fo.write(ln)
    return dst


def compara_attribution(cfg: dict, members: Path, spectrum_prefix: Path, genes_gtf: Path, spec: dict,
                        scored: dict, *, chrom: str = HUMAN_CHROM, compara: str | None = None, dev_contigs=()) -> dict:
    """For every scored Compara pair: (band, paralogue group, direct alignment, same default family, the families
    holding it, both genes family members, same / cross chromosome, touches a development contig). Checked against
    score.py's band totals. chrom 'ALL' = genome-wide (cross-chromosome pairs kept).

    `members` = the default families' copy table. Paralogue groups = connected components of the scored pairs'
    graph. The direct alignment is the spectrum's per-pair `recovered` flag, which must be the NUCLEOTIDE tiers only
    (T1 or T2: a spectrum run with --skip-t3; checked)."""
    score, lib = _bench(cfg)
    everywhere = chrom == "ALL"
    if everywhere:
        spans_all = score.load_gene_spans_all(str(genes_gtf), "gtf")
        gene_of = score.span_mapper_all(spans_all, lambda c, g: g)
        contigs_of = collections.defaultdict(set)
        for c, gs in spans_all.items():
            for g in gs:
                contigs_of[g].add(c)
    else:
        gene_of = score.span_mapper(score.load_gene_spans(str(genes_gtf), chrom, "gtf"))
        contigs_of = collections.defaultdict(lambda: {chrom})
    if widest_tier(spec) != "T1+T2":
        raise RuntimeError(f"{spectrum_prefix}: the spectrum ran the translated tier; the direct alignment of the main "
                           "figure is nucleotide only (run score.py spectrum with --skip-t3)")
    compara_pairs, judge = lib.load_compara(compara or cfg["compara_chr16"], chrom)

    def band_of(k):
        p = compara_pairs[k][0]
        for lo, hi, n in score.BANDS_COMPARA:
            if lo <= p < hi:
                return n
        return "<30"
    recovered: dict = {}
    for ln in open(f"{spectrum_prefix}.truth_pairs.tsv"):
        f = ln.rstrip("\n").split("\t")
        if f[0] == "geneA" or len(f) < 5:
            continue
        recovered[frozenset((f[0], f[1]))] = f[4] == "1"
    uni = set(recovered) & set(compara_pairs)
    fam_genes, fam_spliced, _ = score.load_members(str(members), chrom, gene_of)
    p2f = _pair_groups(fam_genes)
    mem_genes = set().union(*fam_genes.values()) if fam_genes else set()   # genes with >= 1 family copy
    comp = _components(uni)
    members_of = collections.defaultdict(set)
    for g, root in comp.items():
        members_of[root].add(g)
    dev = set(dev_contigs)
    pairs = []
    for k in uni:
        ga, gb = sorted(k)
        pairs.append({"pair": k, "band": band_of(k), "group": comp[ga], "edge": recovered[k],
                      "family": k in p2f, "families": p2f.get(k, set()), "in_families": k <= mem_genes,
                      "cross": not (contigs_of[ga] & contigs_of[gb]),
                      "dev": bool((contigs_of[ga] | contigs_of[gb]) & dev)})
        if pairs[-1]["family"] and not pairs[-1]["in_families"]:
            raise RuntimeError(f"{members}: family pair {sorted(k)} whose genes are not both family members")
    # checks against the scorers' own totals
    for b in COMPARA_BANDS:
        ps = [p for p in pairs if p["band"] == b]
        if not ps:
            continue
        fam_hits = sum(p["family"] for p in ps)
        if (fam_hits, len(ps)) != tuple(scored["recall_universe"].get(b, (None, None))):
            raise RuntimeError(f"{chrom} band {b}: attribution {fam_hits}/{len(ps)} != score.py pairs "
                               f"{scored['recall_universe'].get(b)}")
        n, fr = spec["recall"][b]
        if round(fr["T1+T2"] * n) != sum(p["edge"] for p in ps) or n != len(ps):
            raise RuntimeError(f"{chrom} band {b}: truth_pairs `recovered` disagrees with spectrum.tsv")
    # within-family judged pairs by family (precision concentration)
    prec = {}
    truth_all = set(compara_pairs)
    for measure, groups in (("families_all_copies", fam_genes), ("families_multi_exon", fam_spliced)):
        p2 = _pair_groups(groups)
        jd = {k for k in p2 if k <= judge}
        tp = {k for k in jd if k in truth_all}
        byf = collections.Counter(x for k in jd for x in p2[k])
        top = _top(byf)
        prec[measure] = {"k": len(tp), "n": len(jd), "families": len(byf), "top": top,
                         "top_pairs": byf[top] if top else 0, "top_genes": len(groups[top]) if top else 0,
                         "loo_k": _loo(tp, p2, top) if top else 0, "loo_n": _loo(jd, p2, top) if top else 0,
                         "top_label": _group_label(groups[top]) if top else ""}
        expect = scored["precision_all" if measure == "families_all_copies" else "precision_multi_exon"]
        if (len(tp), len(jd)) != tuple(expect):
            raise RuntimeError(f"{chrom} {measure}: attribution {len(tp)}/{len(jd)} != score.py pairs {expect}")
    return {"pairs": pairs, "members": dict(members_of), "precision": prec,
            "family_genes": {f: len(g) for f, g in fam_genes.items()}}


def compara_universe_attribution(cfg: dict, members: Path, universe: Path, genes_gtf: Path, compara: str,
                                 scored: dict, bands=(">=90", "80-90")) -> dict:
    """fig6s_seeding, Compara reference: per band, the recovered universe pairs (both genes with copies in one
    family), the family holding the most of them and the count without it; precision's judged pairs likewise.
    Reproduces score.py pairs --chrom ALL --universe with its own loaders and checks every total against `scored`."""
    score, lib = _bench(cfg)
    gene_of = score.span_mapper_all(score.load_gene_spans_all(str(genes_gtf), "gtf"), lambda c, g: g)
    compara_pairs, judge = lib.load_compara(compara, "ALL")
    fam_genes, _, _ = score.load_members(str(members), "ALL", gene_of)
    p2f = _pair_groups(fam_genes)
    uni = set()
    for ln in open(universe):
        f = ln.rstrip("\n").split("\t")
        if f[0] != "geneA" and len(f) >= 2:
            uni.add(frozenset((f[0], f[1])))
    uni &= set(compara_pairs)

    def band_of(k):
        p = compara_pairs[k][0]
        for lo, hi, n in score.BANDS_COMPARA:
            if lo <= p < hi:
                return n
        return "<30"
    out: dict = {}
    for b in bands:
        tb = [k for k in uni if band_of(k) == b]
        hit = [k for k in tb if k in p2f]
        byf = collections.Counter(x for k in hit for x in p2f[k])
        top = _top(byf)
        out[b] = {"k": len(hit), "n": len(tb), "top": top or "", "top_pairs": byf[top] if top else 0,
                  "loo_k": _loo(hit, p2f, top) if top else len(hit)}
        if (len(hit), len(tb)) != tuple(scored["recall_universe"].get(b, (0, 0))):
            raise RuntimeError(f"{members}: band {b} attribution {len(hit)}/{len(tb)} != score.py pairs "
                               f"{scored['recall_universe'].get(b)}")
    jd = {k for k in p2f if k <= judge}
    tp = {k for k in jd if k in compara_pairs}
    byf = collections.Counter(x for k in jd for x in p2f[k])
    top = _top(byf)
    out["precision"] = {"k": len(tp), "n": len(jd), "top": top or "", "top_pairs": byf[top] if top else 0,
                        "loo_k": _loo(tp, p2f, top) if top else len(tp),
                        "loo_n": _loo(jd, p2f, top) if top else len(jd)}
    if (len(tp), len(jd)) != tuple(scored["precision_all"]):
        raise RuntimeError(f"{members}: precision attribution {len(tp)}/{len(jd)} != score.py pairs "
                           f"{scored['precision_all']}")
    return out


# ================================================================ genome-wide (every sample)
# docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md fixes every rule below; the constants are its sections 1-4.
PH_DIR_DEFAULT = "/mnt/linuxdisk/tmp/rustle_figures_dev/truth/gw"      # truth.py protein-homology --chrom ALL PREFIXes
PH_SPECIES_DIR = {"human": "human", "gorilla": "gorilla", "chimpanzee": "chimp", "orangutan": "orangutan"}
COMPARA_GW_DEFAULT = "/mnt/linuxdisk/tmp/rustle_figures_dev/truth/compara/human_e116.tsv"
SPECIES_ORDER = ["human", "gorilla", "chimpanzee", "orangutan"]
GW_BUDGET_S = 540                 # a genome-scope build call stops (exit 75) after about this many seconds
SPECTRUM_SHARD_BP = 3_000_000     # query bases per wrapper shard of the spectrum's T1/T2 (fixed at the first split)
MRNA_SHARD_BP = 5_000_000
MRNA_MM2 = ["-x", "asm20", "-k11", "-w5", "-c", "-X", "-N", "100", "-p", "0.1", "--secondary=yes"]
SPECTRUM_FLAGS = {"T1": "-x asm20", "T2": "-x asm20 -k11 -w5"}    # score.py spectrum's mm2_tier flags
ANNOTATION_VERSION = "1"          # bump when annotation_cache changes
PRIMARY_GW_VERSION = "2"             # 2: every gene / pseudogene record (1: protein-homology genes)

# contig exposure (PREREG section 1); anything else: never used for a family decision
DEVELOPMENT, THRESHOLD, REUSED, SCORED_ONCE, NEVER = ("development", "threshold selection", "reused verdict set",
                                                      "scored once, no decision", "never used for a family decision")
EXPOSURE_CLASSES = [DEVELOPMENT, THRESHOLD, REUSED, SCORED_ONCE, NEVER]
EXPOSURE = {("human", "chr16"): DEVELOPMENT,
            ("human", "chr5"): THRESHOLD, ("human", "chr7"): THRESHOLD, ("human", "chr21"): THRESHOLD,
            ("human", "chr2"): REUSED, ("human", "chr8"): REUSED, ("human", "chr10"): REUSED,
            ("human", "chr6"): SCORED_ONCE,
            ("gorilla", "NC_073244.2"): DEVELOPMENT, ("gorilla", "NC_073234.2"): SCORED_ONCE}
SUBSTRATES = {"S0": "whole genome", "S1": "genome minus development contigs",
              "S2": "genome minus every contig used for a family decision"}


def exposure(species: str, contig: str) -> str:
    return EXPOSURE.get((species, contig), NEVER)


def substrate_drop(species: str, sub: str) -> set:
    """Contigs a substrate leaves out (S0 none; S1 development; S2 every contig used for a family decision)."""
    if sub == "S0":
        return set()
    keep_out = (DEVELOPMENT,) if sub == "S1" else (DEVELOPMENT, THRESHOLD, REUSED)
    return {c for (sp, c), cls in EXPOSURE.items() if sp == species and cls in keep_out}


class NotBuilt(RuntimeError):
    """A run-cache product or reference table this figure needs does not exist yet (the message names the command)."""


def registry(cfg: dict) -> dict:
    import samples
    return samples.registry(cfg)


def samples_by_species(cfg: dict) -> dict:
    """{species: [sample ids in registry order]}, species in SPECIES_ORDER."""
    out: dict = {}
    for sid, row in registry(cfg).items():
        out.setdefault(row["species"], []).append(sid)
    return {sp: out[sp] for sp in SPECIES_ORDER if sp in out} | {sp: v for sp, v in out.items()
                                                                 if sp not in SPECIES_ORDER}


def species_of(cfg: dict, sid: str) -> str:
    return registry(cfg)[sid]["species"]


def species_rep(cfg: dict, species: str) -> str:
    return samples_by_species(cfg)[species][0]


def sample_label(cfg: dict, sid: str) -> str:
    return assembly.sample_label(cfg, sid)


def gw_budget(cfg: dict, fig: str) -> "assembly.Budget":
    """The call budget of a genome-scope build: `<fig>_budget_s` / `figs_budget_s`, else GW_BUDGET_S."""
    b = assembly.Budget(cfg, fig)
    if b.limit <= 0:
        b.limit = float(GW_BUDGET_S)
    return b


def need(budget, est_s: float, what: str):
    """Stop the call (exit 75) unless about `est_s` seconds remain for the next heavy unit."""
    if budget.remaining() < est_s:
        raise assembly.Pending(f"{budget.fig}: {budget.remaining():.0f} s left in this call, next unit ({what}) "
                               f"needs about {est_s:.0f} s")


def stage_product(cfg: dict, sid: str, stage: str, name: str) -> Path:
    """A run-cache product WITHOUT running the stage: fresh -> its path; adoptable -> stamped (no work) and returned;
    otherwise NotBuilt naming the `make.py runs` command (figure builds never run a pipeline stage)."""
    import samples
    if stage not in samples.STAGES:
        raise NotBuilt(f"{sid}: run-cache stage '{stage}' is not defined in figures/samples.py (requested: the "
                       f"driver's families stage on {sid}.primary.gtf, suffix .primary, needs assemble_primary)")
    state, reason = samples.status(cfg, sid, stage)
    if state == "fresh":
        return samples.product(cfg, sid, stage, name)
    if state == "adopt":
        return samples.ensure(cfg, sid, stage)[name]   # writes the stamp only, runs nothing
    raise NotBuilt(f"{sid} {stage} is {state} ({reason}): run `python3 figures/make.py runs --sample {sid} --stage "
                   f"{stage}` (repeat while it exits 75)")


def gw_root(cfg: dict) -> Path:
    """${work}/families_gw: the genome-wide caches figures 6 and 7 share."""
    d = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "families_gw"
    d.mkdir(parents=True, exist_ok=True)
    return d


def species_dir(cfg: dict, species: str) -> Path:
    d = gw_root(cfg) / "species" / species
    d.mkdir(parents=True, exist_ok=True)
    return d


def sample_dir(cfg: dict, sid: str) -> Path:
    d = gw_root(cfg) / "samples" / sid
    d.mkdir(parents=True, exist_ok=True)
    return d


def _key_ok(path: Path, key: str) -> bool:
    try:
        return path.read_text() == key
    except OSError:
        return False


def _md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mm2_shard(cfg: dict) -> str:
    return str(Path(cfg["repo"]) / "tools" / "mm2_shard.sh")


def shard_env(budget, bp: int | None = None, margin_s: float = 45.0) -> dict:
    """Environment of one tools/mm2_shard.sh call: a deadline inside this build call's budget, the shard size."""
    env = {}
    rem = budget.remaining()
    if rem != math.inf:
        env["MM2_SHARD_DEADLINE"] = str(int(time.time() + rem - margin_s))
    if bp:
        env["MM2_SHARD_BP"] = str(int(bp))
    return env


def run_sharded_paf(cfg: dict, out: Path, flags: list, fasta: Path, budget, what: str, log: Path,
                    bp: int | None = None, guided_recipe: bool = False) -> Path:
    """One resumable all-vs-all of `fasta` against itself (tools/mm2_shard.sh; cmp-identical to one minimap2 run).
    Exit 75 of the wrapper -> assembly.Pending (every finished shard is kept); 76 -> a single shard does not fit."""
    if out.exists():
        return out
    need(budget, 120, what)
    threads = str(cfg.get("threads", "4"))
    cmd = ([mm2_shard(cfg), "guided", str(fasta), str(out), threads] if guided_recipe else
           [mm2_shard(cfg), "paf", str(out), *flags, "-t", threads, str(fasta), str(fasta)])
    rc = figlib.run(cmd, log=log, env=shard_env(budget, bp), check=False)
    if rc == 75:
        raise assembly.Pending(f"{what}: all-vs-all shards remain (finished shards kept; see {log})")
    if rc != 0:
        raise RuntimeError(f"{what}: tools/mm2_shard.sh exit {rc} (76 = one shard cannot fit the budget: lower "
                           f"MM2_SHARD_BP or raise the budget); see {log}")
    return out


# ---------------------------------------------------------------- annotation caches (one pass per species)
GENE_TYPES = ("gene", "pseudogene")
# ID / Name / gene_biotype anchored (as _gene_exons); Parent and gene= unanchored, exactly lib.longest_cds's regexes
_ATTR = {"ID": re.compile(r"(?:^|;)ID=([^;]+)"), "Name": re.compile(r"(?:^|;)Name=([^;]+)"),
         "gene_biotype": re.compile(r"gene_biotype=([^;]+)"), "Parent": re.compile(r"Parent=([^;]+)"),
         "gene": re.compile(r"gene=([^;]+)")}


def annotation_gff(cfg: dict, species: str) -> Path:
    p = registry(cfg)[species_rep(cfg, species)]["annotation_gff"]
    if not p or not Path(p).exists():
        raise NotBuilt(f"{species}: no annotation GFF in the sample registry ({p})")
    if str(p).endswith(".gz"):
        raise NotBuilt(f"{species}: annotation GFF {p} is compressed; mcl_families and the scorers read plain text")
    return Path(p)


def annotation_cache(cfg: dict, species: str) -> dict:
    """One pass over the species' RefSeq GFF (contig by contig; the GFF must list each contig in one block):

    genes.tsv       contig, name, type, biotype, start1, end, exons (merged, 0-based half-open 's-e,...'): every
                    gene / pseudogene record; exons = every exon record whose Parent chain reaches the record (the
                    `_gene_exons` rule)
    genes_only.gff  the gene / pseudogene / ncRNA_gene lines verbatim (all family_score and score.py read of a GFF)
    mrna_models.tsv contig, gene (CDS `gene=`), strand, transcript, cds_bases, exons (1-based closed): per gene symbol
                    the transcript with the most CDS bases, first maximum in file order (lib.longest_cds's rule, so
                    the transcript whose protein entered the protein-homology families)
    exons.gtf       (human, whose annotation GTF names the Compara symbols) its exon lines, uncompressed: spectrum
                    --ref, pairs --genes"""
    gff = annotation_gff(cfg, species)
    d = species_dir(cfg, species)
    row = registry(cfg)[species_rep(cfg, species)]
    gtf = row.get("annotation_gtf") if species == "human" else None   # exons.gtf: the Compara panels only (human)
    paths = {"genes": d / "genes.tsv", "genes_gff": d / "genes_only.gff", "mrna_models": d / "mrna_models.tsv"}
    if gtf:
        paths["exons_gtf"] = d / "exons.gtf"
    key = f"v{ANNOTATION_VERSION}\n{figlib.file_fingerprint(gff)}\n{figlib.file_fingerprint(gtf) if gtf else ''}\n"
    kf = d / "annotation.key"
    if _key_ok(kf, key) and all(p.exists() for p in paths.values()):
        return paths
    print(f"[fig6/7] annotation cache {species}: one pass over {gff}", file=sys.stderr)
    tmp = {k: p.with_suffix(p.suffix + ".tmp") for k, p in paths.items()}
    fg, fgff, fm = open(tmp["genes"], "w"), open(tmp["genes_gff"], "w"), open(tmp["mrna_models"], "w")
    fg.write("contig\tname\ttype\tbiotype\tstart1\tend\texons\n")
    fgff.write("##gff-version 3\n")
    fm.write("contig\tgene\tstrand\ttranscript\tcds_bases\texons\n")
    done: set = set()
    st: dict = {}

    def reset():
        st.update(genes=[], ftype={}, parent={}, exons=collections.defaultdict(list), cds=collections.OrderedDict(),
                  tx_gene={}, tx_strand={})

    def flush(contig):
        if contig is None:
            return
        if contig in done:
            raise RuntimeError(f"{gff}: contig {contig} appears in more than one block; sort the GFF by contig")
        done.add(contig)
        ftype, parent, exons = st["ftype"], st["parent"], st["exons"]

        def gene_of(x):
            for _ in range(8):
                if x is None or ftype.get(x) in GENE_TYPES:
                    break
                x = parent.get(x)
            return x if ftype.get(x) in GENE_TYPES else None
        by_gene = collections.defaultdict(list)
        for p, ex in exons.items():
            g = gene_of(p)
            if g:
                by_gene[g].extend(ex)
        for gid, name, typ, bt, s1, e in st["genes"]:
            ex = _merge(by_gene.get(gid, []))
            fg.write(f"{contig}\t{name}\t{typ}\t{bt}\t{s1}\t{e}\t{','.join(f'{a}-{b}' for a, b in ex)}\n")
        best: dict = {}
        for tx, n in st["cds"].items():
            sym = st["tx_gene"].get(tx)
            if sym and (sym not in best or n > best[sym][0]):
                best[sym] = (n, tx)
        for sym, (n, tx) in best.items():
            ex = sorted((a + 1, b) for a, b in exons.get(tx, []))
            if ex:
                fm.write(f"{contig}\t{sym}\t{st['tx_strand'][tx]}\t{tx}\t{n}\t{','.join(f'{a}-{b}' for a, b in ex)}\n")
    reset()
    cur = None
    try:
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[0] != cur:
                    flush(cur)
                    reset()
                    cur = f[0]
                t, a = f[2], f[8]
                mid, mp = _ATTR["ID"].search(a), _ATTR["Parent"].search(a)
                if t in GENE_TYPES or t == "ncRNA_gene":
                    fgff.write(line if line.endswith("\n") else line + "\n")
                if t in GENE_TYPES:
                    mn, mb = _ATTR["Name"].search(a), _ATTR["gene_biotype"].search(a)
                    if mn and mid:
                        st["genes"].append((mid.group(1), mn.group(1), t, mb.group(1) if mb else "", int(f[3]),
                                            int(f[4])))
                if mid:
                    st["ftype"][mid.group(1)] = t
                    if mp:
                        st["parent"][mid.group(1)] = mp.group(1).split(",")[0]
                if t == "exon" and mp:
                    st["exons"][mp.group(1).split(",")[0]].append((int(f[3]) - 1, int(f[4])))
                elif t == "CDS" and mp:
                    tx = mp.group(1)   # longest_cds keys the transcript on the whole Parent value
                    st["cds"][tx] = st["cds"].get(tx, 0) + int(f[4]) - int(f[3]) + 1
                    st["tx_strand"][tx] = f[6]
                    mg = _ATTR["gene"].search(a)
                    if mg:
                        st["tx_gene"][tx] = mg.group(1)
            flush(cur)
    finally:
        fg.close(); fgff.close(); fm.close()
    if gtf:
        with assembly._open(gtf) as fi, open(tmp["exons_gtf"], "w") as fo:
            for line in fi:
                if not line.startswith("#") and line.split("\t", 3)[2:3] == ["exon"]:
                    fo.write(line)
    for k, p in paths.items():
        tmp[k].replace(p)
    kf.write_text(key)
    return paths


def _merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def read_genes(path: Path) -> dict:
    """{(contig, name): (start0, end, [exon (s, e)], type, biotype)} from genes.tsv (last record of a name wins, as
    score.load_gene_spans reads a GFF)."""
    out = {}
    with open(path) as fh:
        next(fh)
        for ln in fh:
            c, n, t, bt, s1, e, ex = ln.rstrip("\n").split("\t")
            ivs = [tuple(map(int, x.split("-"))) for x in ex.split(",") if x]
            out[(c, n)] = (int(s1) - 1, int(e), ivs, t, bt)
    return out


def contig_names(cfg: dict, species: str) -> dict:
    """{contig: chromosome name} from the GFF's region records (`chromosome=`), else the contig itself."""
    d = species_dir(cfg, species)
    p = d / "contig_names.json"
    gff = annotation_gff(cfg, species)
    key = figlib.file_fingerprint(gff)
    if p.exists():
        j = json.loads(p.read_text())
        if j.get("key") == key:
            return j["names"]
    names: dict = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) > 8 and f[2] == "region":
                m = re.search(r"chromosome=([^;\n]+)", f[8])
                if m:
                    names[f[0]] = f"chr{m.group(1)}"
    p.write_text(json.dumps({"key": key, "names": names}))
    return names


# ---------------------------------------------------------------- reference families and pairs
def protein_homology_families(cfg: dict, species: str) -> Path:
    """The species' genome-wide protein-homology families (Gene Name, Family ID, Contig)."""
    prefix = Path(cfg.get("protein_homology_dir", PH_DIR_DEFAULT)) / PH_SPECIES_DIR.get(species, species) / "ph"
    p = Path(f"{prefix}.families.tsv")
    if not p.exists():
        raise NotBuilt(f"{species}: genome-wide protein-homology families absent ({p}): run `flock "
                       f"/mnt/linuxdisk/tmp/rustle_heavy.lock python3 bench/truth.py protein-homology --gff <GFF> "
                       f"--genome <FASTA> --chrom ALL --out {prefix} --budget-s 540` until it exits 0")
    with open(p) as fh:
        head = fh.readline().rstrip("\n").split("\t")
    if head[:3] != ["Gene Name", "Family ID", "Contig"]:
        raise RuntimeError(f"{p}: expected the genome-wide header Gene Name / Family ID / Contig, got {head}")
    return p


def read_ph(path: Path) -> list[tuple[str, str, str]]:
    """[(contig, name, family)] of a genome-wide protein-homology table, file order."""
    out = []
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3:
                out.append((f[2], f[0], f[1]))
    return out


def compara_gw(cfg: dict) -> Path:
    p = Path(cfg.get("compara_gw", COMPARA_GW_DEFAULT))
    if not p.exists():
        raise NotBuilt(f"genome-wide Compara table absent ({p}): run `python3 bench/truth.py compara --out "
                       f"{str(p)[:-4]} --release 116` (network; repeat while it exits 75)")
    return p


def filter_rows(src: Path, dst: Path, contig_col: str | int, drop=(), keep=None) -> Path:
    """Copy a TSV keeping the header and the rows whose contig column is not in `drop` (and in `keep`, if given)."""
    drop, keep = set(drop), (set(keep) if keep is not None else None)
    with open(src) as fi, open(dst, "w") as fo:
        head = fi.readline()
        fo.write(head)
        ci = contig_col if isinstance(contig_col, int) else head.rstrip("\n").split("\t").index(contig_col)
        for ln in fi:
            c = ln.rstrip("\n").split("\t")[ci]
            if c not in drop and (keep is None or c in keep):
                fo.write(ln)
    return dst


# ---------------------------------------------------------------- panel d inputs (per species, per sample)
def mrna_fasta(cfg: dict, species: str) -> Path:
    """The annotated mRNA (mrna_models.tsv) of every protein-homology family member, `>CONTIG:NAME`, sorted."""
    import pysam
    ph = protein_homology_families(cfg, species)
    models = annotation_cache(cfg, species)["mrna_models"]
    fasta = registry(cfg)[species_rep(cfg, species)]["fasta"]
    d = species_dir(cfg, species)
    out, miss = d / "mrna.fa", d / "mrna.missing.tsv"
    key = f"{figlib.file_fingerprint(ph)}\n{figlib.file_fingerprint(models)}\n"
    if _key_ok(d / "mrna.fa.key", key) and out.exists():
        return out
    mod = {}
    with open(models) as fh:
        next(fh)
        for ln in fh:
            c, g, strand, tx, n, ex = ln.rstrip("\n").split("\t")
            mod[(c, g)] = (strand, [tuple(map(int, x.split("-"))) for x in ex.split(",")])
    fa = pysam.FastaFile(fasta)
    _, lib = _bench(cfg)
    genes = sorted({(c, n) for c, n, _ in read_ph(ph)})
    tmp = out.with_suffix(".fa.tmp")
    with open(tmp, "w") as fo, open(miss, "w") as fm:
        fm.write("contig\tname\n")
        for c, n in genes:
            m = mod.get((c, n))
            if not m:
                fm.write(f"{c}\t{n}\n")
                continue
            fo.write(f">{c}:{n}\n{lib.spliced1(fa, c, m[1], m[0])}\n")
    tmp.replace(out)
    (d / "mrna.fa.key").write_text(key)
    return out


def mrna_pairs_paf(cfg: dict, species: str, budget) -> Path:
    """Same-family pairs' best annotated-mRNA alignment: the all-vs-all (MRNA_MM2, sharded) reduced to one 11-column
    record per unordered same-family pair, the record with the highest matches / block length (first maximum).
    score.py's --bands paf: reads columns 1, 6, 10 and 11 and keeps each pair's maximum, and it bands truth pairs
    only, so the reduced file gives it exactly the bands of the full PAF."""
    d = species_dir(cfg, species)
    out = d / "mrna.pairs.paf"
    fa = mrna_fasta(cfg, species)
    ph = protein_homology_families(cfg, species)
    key = f"{_md5(fa)}\n{' '.join(MRNA_MM2)}\n{figlib.file_fingerprint(ph)}\n"
    if _key_ok(d / "mrna.pairs.paf.key", key) and out.exists():
        return out
    full = d / f"mrna.{hashlib.md5(key.encode()).hexdigest()[:10]}.paf"
    run_sharded_paf(cfg, full, MRNA_MM2, fa, budget, f"{species} annotated-mRNA all-vs-all", d / "mrna.mm2.log",
                    bp=MRNA_SHARD_BP)
    fam = {f"{c}:{n}": f for c, n, f in read_ph(ph)}
    best: dict = {}
    with open(full) as fh:
        for ln in fh:
            f = ln.split("\t", 12)
            if f[0] == f[5] or fam.get(f[0]) is None or fam.get(f[0]) != fam.get(f[5]):
                continue
            k = frozenset((f[0], f[5]))
            idn = int(f[9]) / max(1, int(f[10]))
            if k not in best or idn > best[k][0]:
                best[k] = (idn, f[0], f[5], f[9], f[10])
    tmp = out.with_suffix(".paf.tmp")
    with open(tmp, "w") as fo:
        for _, q, t, nm, bl in sorted(best.values(), key=lambda v: (v[1], v[2])):
            fo.write(f"{q}\t0\t0\t0\t+\t{t}\t0\t0\t0\t{nm}\t{bl}\n")
    tmp.replace(out)
    (d / "mrna.pairs.paf.key").write_text(key)
    return out


def primary_counts_gw(cfg: dict, sid: str, budget) -> Path:
    """Per gene / pseudogene record of the sample's species (the annotation cache's genes.tsv), the `primary_counts`
    rule genome-wide, one BAM pass per contig (cached per contig, resumable): n_primary_exon (reads whose primary
    record has an aligned block on the record's annotated exons) and n_primary_span. Header: Gene Name, Contig,
    biotype, n_exon_records, n_primary_exon, n_primary_span (score.py --expressed reads the Contig column as the
    CONTIG:NAME key). Version 2 counts every record (version 1 counted the protein-homology genes only; same rule)."""
    import pysam
    species = species_of(cfg, sid)
    row = registry(cfg)[sid]
    genes_tsv = annotation_cache(cfg, species)["genes"]
    genes = read_genes(genes_tsv)
    d = sample_dir(cfg, sid)
    parts = d / "primary_counts.parts"
    parts.mkdir(exist_ok=True)
    out = d / "primary_counts.tsv"
    key = f"#key\tv{PRIMARY_GW_VERSION} {figlib.file_fingerprint(row['bam'])} {figlib.file_fingerprint(genes_tsv)}\n"
    if out.exists() and open(out).readline() == key:
        return out
    by_contig = collections.defaultdict(list)
    for (c, n) in genes:
        by_contig[c].append(n)
    order = sorted(by_contig, key=assembly._natural)
    bam = pysam.AlignmentFile(row["bam"], "rb")
    mapped = {s_.contig: s_.mapped for s_ in bam.get_index_statistics()}
    rate = 60_000.0   # records per second before the first measured contig (conservative)
    for c in order:
        part = parts / f"{c}.tsv"
        if part.exists() and open(part).readline() == key:
            continue
        need(budget, 20 + mapped.get(c, 0) / rate, f"primary counts {sid} {c} ({mapped.get(c, 0):,} records)")
        t0 = time.time()
        tmp = part.with_suffix(".tsv.tmp")
        with open(tmp, "w") as fo:
            fo.write(key)
            for n in sorted(by_contig[c]):
                s0, e, ivs, _typ, bt = genes[(c, n)]
                n_exon, n_span = _count_primary(bam, c, s0, e, ivs) if c in mapped else (0, 0)
                fo.write(f"{n}\t{c}\t{bt or '-'}\t{len(ivs)}\t{n_exon}\t{n_span}\n")
        tmp.replace(part)
        if mapped.get(c, 0) > 100_000:
            rate = min(rate, mapped[c] / max(1.0, time.time() - t0))
    tmp = out.with_suffix(".tsv.tmp")
    with open(tmp, "w") as fo:
        fo.write(key)
        fo.write("Gene Name\tContig\tbiotype\tn_exon_records\tn_primary_exon\tn_primary_span\n")
        for c in order:
            with open(parts / f"{c}.tsv") as fh:
                next(fh)
                fo.writelines(fh)
    tmp.replace(out)
    return out


def expressed_gw(cfg: dict, sid: str, counts: Path, drop=(), restrict=None) -> tuple[Path, int, int]:
    """The sensitivity universe (>= EXPRESSED_MIN_PRIMARY primary reads on the exons) outside `drop` contigs, among the
    (contig, name) keys in `restrict` when given (e.g. the protein-homology genes): (file in score.py's --expressed
    format with a Contig column, genes in it, genes considered)."""
    tag = ("all" if not drop else "minus_" + "_".join(sorted(drop))) + ("" if restrict is None else ".restricted")
    out = sample_dir(cfg, sid) / f"expressed.{tag}.tsv"
    n_in = n_all = 0
    with open(counts) as fh, open(out, "w") as fo:
        fo.write("Gene Name\tContig\tn_primary_exon\n")
        for ln in fh:
            if ln.startswith("#") or ln.startswith("Gene Name\t"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[1] in drop or (restrict is not None and (f[1], f[0]) not in restrict):
                continue
            n_all += 1
            if int(f[4]) >= EXPRESSED_MIN_PRIMARY:
                n_in += 1
                fo.write(f"{f[0]}\t{f[1]}\t{f[4]}\n")
    return out, n_in, n_all


def compara_universe(cfg: dict, sid: str, counts: Path, drop=()) -> tuple[Path, int, int]:
    """fig6s_seeding, Compara reference: the Compara pairs (every band, cross-chromosome kept; lib.load_compara 'ALL')
    whose two genes both have a record with >= EXPRESSED_MIN_PRIMARY primary reads on its exons, and neither gene a
    record on a `drop` contig. Genes are symbols (Compara's key). Returns (geneA/geneB file for score.py --universe,
    pairs, expressed symbols). Independent of both seeding configurations."""
    _, lib = _bench(cfg)
    compara = compara_gw(cfg)
    expr, dropped = set(), set()
    with open(counts) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("Gene Name\t"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[1] in drop:
                dropped.add(f[0])
            elif int(f[4]) >= EXPRESSED_MIN_PRIMARY:
                expr.add(f[0])
    expr -= dropped
    pairs, _ = lib.load_compara(str(compara), "ALL")
    tag = "all" if not drop else "minus_" + "_".join(sorted(drop))
    out = sample_dir(cfg, sid) / f"compara_universe.{tag}.tsv"
    n = 0
    with open(out, "w") as fo:
        fo.write("geneA\tgeneB\tcompara_pid\n")
        for k, (pid, _sub) in sorted(pairs.items(), key=lambda kv: sorted(kv[0])):
            if k <= expr:
                a, b = sorted(k)
                fo.write(f"{a}\t{b}\t{pid:.1f}\n")
                n += 1
    return out, n, len(expr)


def families_copies(cfg: dict, sid: str, stage: str = "families") -> Path:
    """The default de novo families' copy table of a sample: run-cache stage `families`, product `copies`
    (<id>.fam.copies.tsv); for `families_primary` the sibling <id>.primary.fam.copies.tsv the same driver stage
    writes. NotBuilt when the stage is not fresh or the table is absent (a mcl_families older than the copy table)."""
    if stage == "families":
        p = Path(stage_product(cfg, sid, stage, "copies"))
    else:
        cl = str(stage_product(cfg, sid, stage, "clusters"))
        p = Path(cl[: -len(".clusters.tsv")] + ".copies.tsv")
    if not p.exists():
        raise NotBuilt(f"{sid} {stage}: the families copy table {p} is absent (written by a mcl_families that predates "
                       f"it): rebuild mcl_families, then `python3 figures/make.py runs --sample {sid} --stage {stage} "
                       "--force`")
    return p


# ---------------------------------------------------------------- panels a-c inputs (human samples)
def spectrum_gw(cfg: dict, sid: str, budget) -> Path:
    """`score.py spectrum --chrom ALL --skip-t3` on the sample's genome-wide assembly, its T1/T2 all-vs-alls through
    tools/mm2_shard.sh (resumable: a call that runs out of budget exits 75 and the next call continues)."""
    species = species_of(cfg, sid)
    gtf = assembly.rustle_product(cfg, sid, "rustle")
    ann = annotation_cache(cfg, species)
    if "exons_gtf" not in ann:
        raise NotBuilt(f"{sid}: no annotation GTF in the registry (the spectrum reads gene symbols from it)")
    compara = compara_gw(cfg)
    fasta = registry(cfg)[sid]["fasta"]
    d = sample_dir(cfg, sid) / "spectrum"
    d.mkdir(exist_ok=True)
    prefix = d / "genome"
    threads = str(cfg.get("threads", "4"))
    cmd = [sys.executable or "python3", score_py(cfg), "spectrum", "--gtf", str(gtf), "--ref", str(ann["exons_gtf"]),
           "--fasta", fasta, "--chrom", "ALL", "--compara", str(compara), "--out", str(prefix), "--threads", threads,
           "--minimap2", mm2_shard(cfg), "--skip-t3"]
    key = json.dumps({"cmd": cmd[1:], "inputs": [figlib.file_fingerprint(p) for p in (gtf, ann["exons_gtf"],
                                                                                     compara, fasta)]}) + "\n"
    kf = Path(f"{prefix}.key")
    if _key_ok(kf, key) and Path(f"{prefix}.spectrum.tsv").exists() and Path(f"{prefix}.truth_pairs.tsv").exists():
        return prefix
    need(budget, 240, f"{sid} spectrum")
    log = d / "genome.spectrum.log"
    rc = figlib.run(cmd, log=log, env=shard_env(budget, SPECTRUM_SHARD_BP, margin_s=150), check=False)
    if rc != 0:
        # a wrapper that ran out of budget makes the spectrum fail: tell "shards left" from a real failure
        nodes = Path(f"{prefix}.nodes.fa")
        for i, (tier, flags) in enumerate(SPECTRUM_FLAGS.items()):
            st = figlib.run([mm2_shard(cfg), "status", *flags.split(), "-c", "-X", "--no-long-join", "-N", "50",
                             "-p", "0.1", "--secondary=yes", "-t", threads, str(nodes), str(nodes)],
                            log=d / f"status.{tier}.log", check=False) if nodes.exists() else -1
            if st == 75 or (st == 3 and i > 0):   # shards left, or T2 not started after a complete T1
                raise assembly.Pending(f"{sid} spectrum: the {tier} all-vs-all has shards left (see {log} and "
                                       f"{d / f'status.{tier}.log'})")
            if st != 0:
                break
        raise RuntimeError(f"{sid} spectrum failed (exit {rc}) and no all-vs-all is waiting for shards; see {log}")
    kf.write_text(key)
    return prefix


def chr16_t3_example(cfg: dict) -> Path | None:
    """The dev-scope chr16 spectrum (A119b, WITH the translated tier T3): the comparator of supplement fig6s_protein
    (T3 is not part of Rustle's rule and is never run genome-wide)."""
    p = figlib.work_dir(cfg, FIG) / "spectrum" / HUMAN_CHROM
    return p if Path(f"{p}.spectrum.tsv").exists() else None
