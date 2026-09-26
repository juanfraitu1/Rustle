"""_o1_recovery — private helpers for Figure 7: family recovery by Rustle's two modes (de novo and guided).

A Rustle-internal comparison: two modes of Rustle, not a comparison with other tools.

    de novo   loci assembled from the IsoSeq reads by the pipeline driver (`tools/rustle_pipeline.sh assemble`,
              default seeding), then the driver's families stage (`mcl_families --from-gtf --min-exonic-bp 1
              --min-shared-exon-frac 0.60`)
    guided    loci = the annotated gene and pseudogene bodies: samtools faidx -r REGIONS, minimap2 all-vs-all,
              `mcl_families --paf --gff --min-exonic-bp 1 --min-shared-exon-frac 0.60`

The de novo mode IS Rustle's one default de novo family definition (user decision 2026-09-25): reads -> seeded
assembly loci -> one representative per locus ("positional exon sum") -> families; the families its copy table
feeds to copy assignment. The guided mode is the same family rule on the annotation's gene bodies.

GENOME scope (default, `fig7_scope genome`; docs/PREREG_genome_wide_families_2026-09-25.md, Amendment 1): one
genome-wide run per mode. De novo = the run-cache stage `families` of each sample (never run here: `make.py runs`);
guided = one run per species (it reads no RNA), its all-vs-all with the de novo flags (`-x asm20 -c -X -N 50 -p 0.1
--secondary=yes`) through tools/mm2_shard.sh (`fig7_guided_flags recipe` = the recorded `-x asm20 -c --eqx -P`).
Scored by `family_score --chrom ALL --per-family --pairwise` against EXTERNAL references: Ensembl Compara families
(duplications within primates; human; the headline), Soto et al. 2025 (human; not independent) and the NPIP reference
set (human chr16, an inset); and, for every species, Liftoff copy pairs (the Fig. 8 self-lift; de novo only, see
below). Substrates of _o1.SUBSTRATES, and per contig (a breakdown of the genome-wide run). The protein-homology
families are a SECONDARY reference of the supplementary figure only (`fig7_protein_homology 1`).

DEV scope (`fig7_scope dev`, or `make.py data fig7 --recorded`): the per-chromosome tables of 2026-09-25 (human chr16,
chr2, chr6, chr8, chr10; gorilla NC_073244.2, NC_073234.2): per-contig runs of both modes, Compara families at
Primates restricted to the chromosome (human), Soto 2025, the NPIP reference set; the contig's own protein-homology
families (`bench/truth.py protein-homology --chrom C`) for the supplementary figure. `fig7_dev_cached 1` rescores the
cached per-contig runs of ${work}/fig7/current without re-running them.

References (never pooled across species; how independent each one is of the modes is stated with it):
    compara   Ensembl Compara release 116 human paralogue pairs; a family = a connected component of the pairs whose
              duplication node is at or below Primates, protein-coding genes (a gene-tree clade, no identity threshold)
    soto      Soto et al. 2025 Table S1C, human: a COVER (first family ID kept, as family_score does). NOT independent
              of the family rule: the 0.60 shared-exon threshold was chosen against Soto families (register 903)
    npip_u2   the NPIP reference set (register 990, human chr16 only): Soto NPIP + RefSeq NPIP-named + genes that align
              >= 95% of their own length to a Soto NPIP member; one gene family, an inset
    liftoff   (genome scope, every species) Liftoff's self-lift (record, extra copy) pairs with both loci read-supported
              in the sample, recovered when they lie in one DE NOVO family (copy table, Liftoff's -a); pairwise
              sensitivity only. The guided mode is not scored on it: its loci are the annotation, and an extra copy is
              unannotated by construction
    referee / protein_homology  (supplement only; secondary) longest CDS per gene, pseudogenes and immunoglobulin /
              T-cell-receptor segments excluded, all-vs-all BLASTP, HSP coverage >= 0.30 of the longer protein, MCL
              I = 2.8; fold-level families mostly outside the nucleotide regime (register 1101)

Substrate status of the dev tables (`fig7_human_contigs` / `fig7_gorilla_contigs`: comma-separated `contig[:status]`,
status one of `untouched`, `reused`, `development`; a bare contig takes its default from DEFAULT_STATUS):
    held out, untouched            no rule, threshold or test result was computed on it (human chr6, gorilla NC_073234.2)
    held out, reused verdict set   not used to develop the early family rules, but the held-out verdict set of about 30
                                   guided-mode tests since 2026-09-20 (registers 903-1035; human chr2, chr8, chr10)
    development                    rules were chosen on it (human chr16; gorilla NC_073244.2)

Scoring semantics (family_score, identical for both modes): a locus is labelled with ONE gene (max overlap, first
maximum); predicted families are intersected with the reference genes (genes of reference families with >= 2
members in the substrate), so a predicted member no reference labels is not scored — precision is an upper bound
(register 770/991), the same bound for both modes.
"""
from __future__ import annotations

import collections
import hashlib
import itertools
import json
import re
import sys
from pathlib import Path

import assembly
import figlib
import _o1

FIG = "fig7"
MODES = ["denovo", "guided"]
MODE_LABEL = {"denovo": "de novo", "guided": "guided"}
TRUTH_LABEL = {"compara": "Ensembl Compara families (primates)", "soto": "Soto 2025",
               "npip_u2": "NPIP reference set", "referee": "protein-homology families (secondary)"}
MAIN_DEV_TRUTHS = ("compara", "soto", "npip_u2")   # main figure; `referee` = the supplement's secondary reference

# substrate status (see the module doc); the table's `status` column carries the long label
UNTOUCHED, REUSED, DEVELOPMENT = "held out, untouched", "held out, reused verdict set", "development"
STATUSES = [UNTOUCHED, REUSED, DEVELOPMENT]
STATUS_KEYS = {"untouched": UNTOUCHED, "reused": REUSED, "development": DEVELOPMENT}
DEFAULT_STATUS = {("human", "chr16"): DEVELOPMENT, ("human", "chr2"): REUSED, ("human", "chr6"): UNTOUCHED,
                  ("human", "chr8"): REUSED, ("human", "chr10"): REUSED,
                  ("gorilla", "NC_073244.2"): DEVELOPMENT, ("gorilla", "NC_073234.2"): UNTOUCHED}
HUMAN_DEFAULT = "chr16,chr2,chr6,chr8,chr10"
GORILLA_DEFAULT = "NC_073244.2:development,NC_073234.2:untouched"
_REUSED_NOTE = ("held out, reused verdict set: {n} O1-ledger mention{pl}, but the held-out verdict "
                "set of about 30 pre-registered guided-mode tests since 2026-09-20 (registers 903-1035; "
                "PREREG_heldout_families_2026-09-20 named chr2 with chr6 as its primary held-out pair); the shipped "
                "family rule was kept over every alternative tested there on the GUIDED node set, so the exposure may "
                "favour guided; no de novo arm is recorded on it before this figure")
STATUS_NOTE = {
    ("human", "chr16"): "development: every early O1 decision was scored here (115 ledger mentions)",
    ("human", "chr2"): _REUSED_NOTE.format(n=0, pl="s"),
    ("human", "chr6"): "held out, untouched: 0 O1-ledger mentions (PREREG_heldout_families_2026-09-20 primary held-out "
                       "set with chr2); its 2026-09-20 guided run could not be scored (no Soto family >= 3 members, "
                       "register 903), and no later test used it",
    ("human", "chr8"): _REUSED_NOTE.format(n=0, pl="s"),
    ("human", "chr10"): _REUSED_NOTE.format(n=1, pl=""),
    ("gorilla", "NC_073244.2"): "development here: the pre-registered held-out verdict set of the seeding decision "
                                "(PREREG_locus_read_pool_2026-09-22, registers 1059/1060; 0.98 kept in 1100), whose "
                                "outcome (this contig's protein-homology F) chose the de novo mode's seeding with "
                                "secondary alignments; 29 ledger, 19 register and 50 pre-registration mentions (2026-09-25)",
    ("gorilla", "NC_073234.2"): "held out, untouched: no ledger, register or pre-registration mention before this "
                                "figure (one O2 table row in bench/COPY_ASSIGNMENT_AND_GATE.md); 1,948 gene + "
                                "pseudogene records, the richest unexposed gorilla contig after NC_086017.1",
}
TRUTHS = {"human": ["compara", "soto", "referee"], "gorilla": ["referee"]}
FLAGSHIP = ("NPIP", "TBC1D3", "FAM90A", "AGAP")   # the thesis families + the classic chr10 SD family

# the recorded runs (registers 988-1019, 1100/1101) used for the provisional tables
REC_HUMAN_DENOVO = {"chr16": "/mnt/linuxdisk/tmp/regress/dn16_fam3.clusters.tsv"}
REC_HUMAN_GUIDED = {"chr16": "/mnt/linuxdisk/tmp/regress/chr16_guided.clusters.tsv",
                    "chr2": "/mnt/linuxdisk/tmp/heldout/chr2_fam.clusters.tsv",
                    "chr8": "/mnt/linuxdisk/tmp/heldout/chr8_fam.clusters.tsv",
                    "chr10": "/mnt/linuxdisk/tmp/heldout/chr10_fam.clusters.tsv"}
REC_GORILLA_DENOVO = {"NC_073244.2": "/mnt/linuxdisk/tmp/gw22/sec/ggo44_GOOD0.98.fam.clusters.tsv"}
REC_REFEREE_PREFIX = "/mnt/linuxdisk/tmp/referee/{chrom}_ref"          # cached proteins.faa + blastp.tsv (human)
REC_GORILLA_REFEREE = "/mnt/linuxdisk/tmp/gw22/sec/ref/{contig}.tsv"
REC_GORILLA_GENES = "/mnt/linuxdisk/tmp/gw22/sec/ref/{contig}.genes.gff"
REC_NPIP_U2 = "/mnt/linuxdisk/tmp/union2/U2_truth.tsv"
DEFAULT_HUMAN_GFF = "/mnt/linuxdisk/tmp/regress/chm13.gff"             # uncompressed RefSeq CHM13 v2.0 full GFF
DEFAULT_GUIDED_PAF_DIRS = "/mnt/linuxdisk/tmp/heldout,/mnt/linuxdisk/tmp/regress"


def _status(species: str, contig: str, token: str) -> str:
    token = token.strip()
    default = DEFAULT_STATUS.get((species, contig))
    if token in STATUS_KEYS:
        return STATUS_KEYS[token]
    if token in STATUSES:
        return token
    if token in ("", "held out") and default and (token == "" or default != DEVELOPMENT):
        return default   # bare contig, or the legacy 'held out' of a contig whose default is a held-out kind
    raise ValueError(f"fig7 {species} contig {contig!r}: status {token or '(none)'!r} is not one of "
                     f"{', '.join(STATUS_KEYS)} (write {contig}:untouched, {contig}:reused or {contig}:development)")


def parse_contigs(species: str, spec: str) -> list[tuple[str, str]]:
    """`contig[:status],...` -> [(contig, long status label)]; see the module doc for the syntax."""
    out = []
    for item in spec.split(","):
        c, _, st = item.strip().partition(":")
        if c:
            out.append((c, _status(species, c, st)))
    return out


def gorilla_contigs(cfg: dict) -> list[tuple[str, str]]:
    return parse_contigs("gorilla", cfg.get("fig7_gorilla_contigs") or GORILLA_DEFAULT)


def substrates(cfg: dict) -> list[tuple[str, str, str]]:
    """[(species, contig, status)] in plotting order."""
    hc = parse_contigs("human", cfg.get("fig7_human_contigs") or HUMAN_DEFAULT)
    return [("human", c, s) for c, s in hc] + [("gorilla", c, s) for c, s in gorilla_contigs(cfg)]


def npip_u2_truth(cfg: dict) -> Path | None:
    """The configured NPIP union truth; a configured but absent file is an error (the caption counts its row), and
    `npip_union_truth none` drops it on purpose."""
    v = cfg.get("npip_union_truth", REC_NPIP_U2)
    if not v or v.lower() == "none":
        return None
    if not Path(v).exists():
        raise FileNotFoundError(f"npip_union_truth {v} is absent: rebuild it (register 990, `git show "
                                f"8db314c7^:bench/union_truth_npip.py`) or set `npip_union_truth none`")
    return Path(v)


def bin_path(cfg: dict, name: str) -> str:
    return str(Path(cfg["bin"]) / name)


def bench(cfg: dict) -> Path:
    return Path(cfg["repo"]) / "bench"


# ================================================================ inputs
def gff_slices(src, dsts: dict, *, force=False) -> dict:
    """{contig: dst}: the lines of `src` on each contig (comments dropped), in ONE pass over `src` — per-contig GFFs
    for the scorer and the referee (a 1.7 GB GFF read once, not once per contig)."""
    dsts = {c: Path(d) for c, d in dsts.items()}
    todo = {c: d for c, d in dsts.items() if force or not figlib.fresh(d, src)}
    if todo:
        fhs = {c: open(d.with_suffix(d.suffix + ".tmp"), "w") for c, d in todo.items()}
        try:
            with assembly._open(src) as fi:
                for line in fi:
                    fo = fhs.get(line.split("\t", 1)[0])
                    if fo is not None:
                        fo.write(line)
        finally:
            for fo in fhs.values():
                fo.close()
        for c, d in todo.items():
            d.with_suffix(d.suffix + ".tmp").replace(d)
    return dsts


def gene_regions(gff: Path, contig: str) -> list[str]:
    """`chrom:start-end` of every gene / pseudogene record (1-based, the GFF's own), sorted unique — the guided node
    set (the awk of PREREG_heldout_families §2, sorted as `sort -u` sorts in the C locale)."""
    regs = set()
    for line in open(gff):
        f = line.split("\t", 8)
        if len(f) > 4 and f[0] == contig and f[2] in ("gene", "pseudogene"):
            regs.add(f"{f[0]}:{f[3]}-{f[4]}")
    return sorted(regs, key=lambda s: s.encode())


# ================================================================ the family_score semantics, in Python
def gene_spans(gff, chrom: str) -> list[tuple[int, int, str]]:
    """(start, end, Name) of gene / pseudogene / ncRNA_gene records on `chrom`, as family_score reads them (the file's
    own coordinates), sorted as tuples."""
    out = []
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "pseudogene", "ncRNA_gene"):
            continue
        i = f[8].find("Name=")
        if i < 0:
            continue
        out.append((int(f[3]), int(f[4]), f[8][i + 5:].split(";")[0]))
    out.sort()
    return out


def gene_at(spans, s: int, e: int):
    """The gene with the largest overlap (first maximum, strict >) — family_score's `gene_at`."""
    best = None
    for gs, ge, g in spans:
        if ge < s:
            continue
        if gs > e:
            break
        ov = min(e, ge) - max(s, gs)
        if ov > 0 and (best is None or ov > best[0]):
            best = (ov, g)
    return best[1] if best else None


def truth_table(path) -> list[tuple[str, str]]:
    """(gene, family) in file order, first family per gene, blanks and N/A skipped — family_score's `soto_truth`."""
    out, seen = [], set()
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci_f, ci_g = hdr.index("Family ID"), hdr.index("Gene Name")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            fam = r[ci_f].strip() if ci_f < len(r) else ""
            g = r[ci_g].strip() if ci_g < len(r) else ""
            if fam and fam != "N/A" and g and g not in seen:
                seen.add(g)
                out.append((g, fam))
    return out


def load_clusters(path, chrom: str) -> list[tuple[str, list[tuple[int, int]]]]:
    """cluster_id -> members on `chrom`, first-seen order — family_score's `load_clusters`."""
    order, members = [], {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {h: i for i, h in enumerate(hdr)}
        for need in ("cluster_id", "chrom", "start", "end"):
            if need not in ci:
                raise RuntimeError(f"{path}: no column {need}")
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if len(r) < len(hdr) or r[ci["chrom"]] != chrom:
                continue
            cid = r[ci["cluster_id"]]
            if cid not in members:
                order.append(cid)
                members[cid] = []
            members[cid].append((int(r[ci["start"]]), int(r[ci["end"]])))
    return [(c, members[c]) for c in order]


def family_label(genes) -> str:
    """A readable name for a truth family: the longest common prefix of its named (non-LOC) members, trailing digits
    and dashes stripped; else the first named member."""
    named = sorted(g for g in genes if not g.startswith("LOC"))
    if not named:
        return sorted(genes)[0] if genes else ""
    pre = named[0]
    for g in named[1:]:
        while not g.startswith(pre):
            pre = pre[:-1]
    pre = re.sub(r"[\d\-_]+$", "", pre)
    return pre if len(pre) >= 3 else named[0]


def flagship_of(genes):
    for f in FLAGSHIP:
        if any(g.startswith(f) for g in genes):
            return f
    return ""


def score_arm(clusters, gff_or_spans, truth, chrom: str) -> dict:
    """Pooled bipartite + pairwise + per-family results of one arm (family_score's semantics; see module doc)."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment

    spans = gff_or_spans if isinstance(gff_or_spans, list) else gene_spans(gff_or_spans, chrom)
    fam = truth_table(truth)
    cl = load_clusters(clusters, chrom)
    on_chrom = {g for _, _, g in spans}
    t_order, tset = [], {}
    for g, f in fam:
        if g in on_chrom:
            if f not in tset:
                t_order.append(f)
                tset[f] = set()
            tset[f].add(g)
    t_order = [f for f in t_order if len(tset[f]) >= 2]
    universe = set().union(*(tset[f] for f in t_order)) if t_order else set()
    p_order, pset = [], {}
    for cid, mem in cl:
        gs = {g for g in (gene_at(spans, s, e) for s, e in mem) if g in universe}
        if gs:
            p_order.append(cid)
            pset[cid] = gs
    out = {"truth_families": len(t_order), "truth_genes": sum(len(tset[f]) for f in t_order),
           "clusters_scored": len(p_order), "clusters_total": len(cl), "loci_total": sum(len(m) for _, m in cl)}
    if not t_order or not p_order:   # nothing matched: every truth family is missed (and still listed)
        per = [dict(family=t, n_truth=len(tset[t]), cluster="", n_pred=0, hit=0, sens=0.0, prec=0.0, f=0.0,
                    jaccard=0.0, label=family_label(tset[t]), flagship=flagship_of(tset[t])) for t in t_order]
        tpairs = sum(len(tset[t]) * (len(tset[t]) - 1) // 2 for t in t_order)
        out.update(sens=0.0, prec=0.0, f=0.0, matched=0, pred_members=0, pair_tp=0, truth_pairs=tpairs, pred_pairs=0,
                   pair_sens=0.0, pair_prec=None, exact=0, touched=0, per_family=per)
        return out
    M = np.array([[len(tset[t] & pset[p]) for p in p_order] for t in t_order], dtype=np.int64)
    rows, cols = linear_sum_assignment(-M)
    assign = {int(i): int(j) for i, j in zip(rows, cols)}
    matched = int(sum(M[i, j] for i, j in assign.items()))
    pred_members = int(sum(len(pset[p_order[j]]) for i, j in assign.items() if M[i, j] > 0))
    sens = matched / out["truth_genes"]
    prec = matched / pred_members if pred_members else 0.0
    per = []
    for i, t in enumerate(t_order):
        j = assign.get(i)
        hit = int(M[i, j]) if j is not None else 0
        nt = len(tset[t])
        if hit > 0:
            npred = len(pset[p_order[j]])
            s_, p_ = hit / nt, hit / npred
            per.append(dict(family=t, n_truth=nt, cluster=p_order[j], n_pred=npred, hit=hit, sens=s_, prec=p_,
                            f=2 * s_ * p_ / (s_ + p_), jaccard=hit / (nt + npred - hit),
                            label=family_label(tset[t]), flagship=flagship_of(tset[t])))
        else:
            per.append(dict(family=t, n_truth=nt, cluster="", n_pred=0, hit=0, sens=0.0, prec=0.0, f=0.0,
                            jaccard=0.0, label=family_label(tset[t]), flagship=flagship_of(tset[t])))
    # pairwise: within-family truth pairs vs within-cluster predicted pairs, both over the truth universe
    tpairs = {frozenset(p) for t in t_order for p in itertools.combinations(sorted(tset[t]), 2)}
    ppairs = {frozenset(p) for c in p_order for p in itertools.combinations(sorted(pset[c]), 2)}
    tp = len(tpairs & ppairs)
    out.update(sens=sens, prec=prec, f=2 * sens * prec / (sens + prec) if sens + prec else 0.0, matched=matched,
               pred_members=pred_members, pair_tp=tp, truth_pairs=len(tpairs), pred_pairs=len(ppairs),
               pair_sens=tp / len(tpairs) if tpairs else 0.0, pair_prec=tp / len(ppairs) if ppairs else None,
               exact=sum(1 for r in per if r["f"] == 1.0), touched=sum(1 for r in per if r["hit"] > 0),
               per_family=per)
    return out


_FS = re.compile(r"\|\s*truth (\d+) fams / (\d+) genes \| clusters (\d+) \| sens ([\d.]+) prec ([\d.]+) F ([\d.]+) \| "
                 r"collapsed (\d+) \| no-locus (\d+)")


def family_score(cfg: dict, clusters, gff, truth, chrom: str, label: str, log: Path) -> dict:
    """Run target/release/family_score (the registers' scorer) and parse its one line."""
    figlib.run([bin_path(cfg, "family_score"), "--clusters", str(clusters), "--gff", str(gff), "--soto", str(truth),
                "--chrom", chrom, "--label", label], log=log)
    for line in open(log):
        m = _FS.search(line)
        if m:
            k = ["truth_families", "truth_genes", "clusters_scored", "sens", "prec", "f", "collapsed", "no_locus"]
            v = [int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4)), float(m.group(5)),
                 float(m.group(6)), int(m.group(7)), int(m.group(8))]
            return dict(zip(k, v))
        if "no scoreable truth/prediction overlap" in line:
            return {"empty": True}
    raise RuntimeError(f"{log}: no family_score result line")


def check_against_family_score(mine: dict, fs: dict, what: str):
    """The per-family / pairwise re-derivation must reproduce family_score's pooled numbers exactly."""
    if fs.get("empty"):
        if mine["clusters_scored"]:
            raise RuntimeError(f"{what}: family_score found no overlap, the re-derivation found {mine['clusters_scored']}")
        return
    for k in ("truth_families", "truth_genes", "clusters_scored"):
        if mine[k] != fs[k]:
            raise RuntimeError(f"{what}: {k} {mine[k]} != family_score {fs[k]}")
    for k in ("sens", "prec", "f"):
        if f"{mine[k]:.3f}" != f"{fs[k]:.3f}":
            raise RuntimeError(f"{what}: {k} {mine[k]:.4f} != family_score {fs[k]:.3f}")


# ================================================================ truths
def referee_recorded(cfg: dict, species: str, contig: str, gff: Path, dst: Path) -> Path | None:
    """Provisional truth: the referee families from a RECORDED cache (read-only: human = the cached proteins.faa +
    blastp.tsv of registers 1006/995, re-clustered by truth.protein_referee exactly as score.py referee does with
    reuse_faa; gorilla = the recorded Gene Name / Family ID table)."""
    if species == "gorilla":
        p = Path(REC_GORILLA_REFEREE.format(contig=contig))
        return p if p.exists() else None
    prefix = REC_REFEREE_PREFIX.format(chrom=contig)
    if not (Path(prefix + ".proteins.faa").exists() and Path(prefix + ".blastp.tsv").exists()):
        return None
    if figlib.fresh(dst, gff, prefix + ".blastp.tsv", prefix + ".proteins.faa"):
        return dst
    sys.path.insert(0, str(bench(cfg)))
    import pysam
    import truth as truthlib
    fa = pysam.FastaFile(cfg["human_fasta"])
    fams = truthlib.protein_referee(str(gff), fa, contig, prefix, 2, reuse_faa=True)   # reads only: both cached
    tmp = dst.with_suffix(".tmp")
    with open(tmp, "w") as fh:
        fh.write("Gene Name\tFamily ID\n")
        for fid, mem in fams.items():
            for g in mem:
                fh.write(f"{g}\t{fid}\n")
    tmp.replace(dst)
    return dst


def referee_build(cfg: dict, species: str, contig: str, gff: Path, wdir: Path, force=False) -> Path:
    """`bench/truth.py protein-homology` on the contig's GFF slice (blastp cached at PREFIX.blastp.tsv): the
    contig's own protein-homology families (dev scope; the old per-chromosome path, byte-identical outputs)."""
    prefix = wdir / f"{species}_{contig}_ref"
    out = Path(f"{prefix}.families.tsv")
    if force or not figlib.fresh(out, gff, bench(cfg) / "truth.py", cfg[f"{species}_fasta"]):
        # truth.blastp_all_vs_all caches PREFIX.blastp.tsv on EXISTENCE: a rebuild must not score new proteins
        # against the old all-vs-all
        for stale in [Path(f"{prefix}.blastp.tsv")] + sorted(wdir.glob(f"{prefix.name}_db.*")):
            stale.unlink(missing_ok=True)
        figlib.run([sys.executable or "python3", str(bench(cfg) / "truth.py"), "protein-homology", "--gff", str(gff),
                    "--genome", cfg[f"{species}_fasta"], "--chrom", contig, "--out", str(prefix),
                    "--threads", cfg.get("threads", "4")], log=wdir / f"{species}_{contig}_ref.log")
    return out


# ================================================================ modes (HEAVY; cached under ${work}/fig7/)
def denovo_families(cfg: dict, species: str, contig: str, wdir: Path, force=False) -> Path:
    """Dev scope. De novo: the genome-wide assembly from the run cache (assembly.rustle_product: never assembled
    here), restricted to the contig, through the driver's `families` stage. `force` recomputes the family stage only."""
    genome_gtf = assembly.rustle_product(cfg, species, "rustle")
    prefix = wdir / f"{species}_{contig}.denovo"
    gtf = assembly.restrict_gtf(genome_gtf, Path(f"{prefix}.gtf"), {contig}, force=force)
    clusters = Path(f"{prefix}.fam.clusters.tsv")
    if force or not figlib.fresh(clusters, gtf, Path(cfg["bin"]) / "mcl_families", cfg["driver"],
                                 cfg[f"{species}_fasta"]):
        figlib.run(["bash", cfg["driver"], "families", "--bam", cfg[f"{species}_bam"], "--fasta",
                    cfg[f"{species}_fasta"], "--out", str(prefix), "--bin", cfg["bin"], "--threads",
                    cfg.get("threads", "4")], log=wdir / f"{species}_{contig}.denovo.driver.log")
    return clusters


def _recorded_guided_paf(cfg: dict, contig: str, regions: list[str]):
    """A recorded guided all-vs-all of EXACTLY this node set (same regions file), if one exists: the guided alignment
    depends only on the annotation, the genome and minimap2, so it is reused rather than recomputed (13 min for chr2)."""
    for d in cfg.get("fig7_guided_paf_dirs", DEFAULT_GUIDED_PAF_DIRS).split(","):
        d = Path(d.strip())
        paf, reg = d / f"{contig}.paf", d / f"{contig}.regions"
        if paf.exists() and reg.exists() and [l.strip() for l in open(reg) if l.strip()] == regions:
            return paf
    return None


GUIDED_MM2 = {"recipe": "-x asm20 -c --eqx -P",                                  # the recorded guided recipe
              "denovo": "-x asm20 -c -X -N 50 -p 0.1 --secondary=yes"}               # mcl_families --from-gtf's flags


def guided_families(cfg: dict, species: str, contig: str, gff: Path, full_gff: Path, wdir: Path,
                    force=False, flags: str = "recipe") -> tuple[Path, str]:
    """Guided: annotated gene + pseudogene bodies -> minimap2 all-vs-all -> mcl_families (the documented recipe).
    The node set comes from the contig's GFF slice; mcl_families reads the FULL GFF (exon unions), as the recorded
    guided runs did. flags 'denovo' = the same node set aligned with the de novo mode's minimap2 flags (the bridging
    measurement of the genome-wide pre-registration; prefix `.guided_dnflags`)."""
    prefix = wdir / f"{species}_{contig}.guided" if flags == "recipe" else wdir / f"{species}_{contig}.guided_dnflags"
    regions = gene_regions(gff, contig)
    reg = Path(f"{prefix}.regions")
    if force or not reg.exists() or [l.strip() for l in open(reg) if l.strip()] != regions:
        reg.write_text("\n".join(regions) + "\n")
    paf = _recorded_guided_paf(cfg, contig, regions) if species == "human" and flags == "recipe" else None
    paf_source = f"recorded guided all-vs-all {paf} (identical regions)" if paf else "computed"
    if paf is None:
        paf = Path(f"{prefix}.paf")
        if force or not figlib.fresh(paf, reg, cfg[f"{species}_fasta"]):
            bodies = Path(f"{prefix}.bodies.fa")
            figlib.run(f"samtools faidx {cfg[f'{species}_fasta']} -r {reg} > {bodies}",
                       log=wdir / f"{species}_{contig}.guided.faidx.log")
            figlib.run(f"minimap2 {GUIDED_MM2[flags]} -t {cfg.get('threads', '4')} {bodies} {bodies} > {paf}.tmp "
                       f"&& mv {paf}.tmp {paf}", log=Path(f"{prefix}.mm2.log"))
            bodies.unlink(missing_ok=True)
    clusters = Path(f"{prefix}.clusters.tsv")
    if force or not figlib.fresh(clusters, paf, full_gff, Path(cfg["bin"]) / "mcl_families"):
        figlib.run([bin_path(cfg, "mcl_families"), "--paf", str(paf), "--gff", str(full_gff), "--min-exonic-bp", "1",
                    "--min-shared-exon-frac", "0.60", "--out", str(prefix)],
                   log=wdir / f"{species}_{contig}.guided.mcl.log")
    return clusters, paf_source


# ================================================================ sources
def recorded_sources(cfg: dict, wdir: Path) -> dict:
    """The recorded runs (provisional): human chr16 both modes (2026-09-20: dn16_fam3 = copy_assign --assemble-only on
    chr16, primaries only, and the guided chr16 recipe), human chr2/chr8/chr10 guided (2026-09-20 held-out recipe),
    gorilla NC_073244.2 de novo (2026-09-23 ggo44_GOOD0.98: region run seeded with tied secondaries, the current
    driver default). No recorded de novo run exists on a held-out human chromosome, and no recorded guided run with
    the current rule exists on a gorilla contig: those arms are absent from the provisional tables."""
    gdir = wdir / "gff"
    gdir.mkdir(parents=True, exist_ok=True)
    src = {"source": ("recorded runs: human chr16 de novo dn16_fam3 (2026-09-20: copy_assign --assemble-only "
                      "--assembly-polish full on chr16, loci seeded from PRIMARIES ONLY, the defaults of that date; its "
                      "all-vs-all used the guided recipe's minimap2 flags) + guided chr16_guided "
                      "(2026-09-20), human chr2/chr8/chr10 guided (2026-09-20 held-out recipe), gorilla NC_073244.2 de "
                      "novo ggo44_GOOD0.98 (2026-09-23 chain.sh region run on the contig slice, seeded with tied "
                      "secondaries = today's driver default); human referee families re-clustered read-only from the "
                      "recorded blastp caches (truth.protein_referee, reuse_faa), gorilla referee = recorded table"),
           "arms": {}, "truths": {}, "genes": {}, "extra_inputs": {}}
    hgff = cfg.get("human_ref_gff", DEFAULT_HUMAN_GFF)
    recorded = [(sp, c) for sp, c, _ in substrates(cfg)
                if (REC_HUMAN_DENOVO if sp == "human" else REC_GORILLA_DENOVO).get(c)
                or (sp == "human" and REC_HUMAN_GUIDED.get(c))]
    hslices = gff_slices(hgff, {c: gdir / f"human_{c}.gff" for sp, c in recorded if sp == "human"})
    for species, contig in recorded:
        rec_dn = (REC_HUMAN_DENOVO if species == "human" else REC_GORILLA_DENOVO).get(contig)
        rec_g = REC_HUMAN_GUIDED.get(contig) if species == "human" else None
        genes = hslices[contig] if species == "human" else Path(REC_GORILLA_GENES.format(contig=contig))
        src["genes"][(species, contig)] = genes
        for mode, p in (("denovo", rec_dn), ("guided", rec_g)):
            if p and Path(p).exists():
                src["arms"][(species, contig, mode)] = (Path(p), "recorded")
        ref = referee_recorded(cfg, species, contig, genes, wdir / f"{species}_{contig}.referee.tsv")
        if ref:
            src["truths"][(species, contig, "referee")] = ref
            if species == "human":
                prefix = REC_REFEREE_PREFIX.format(chrom=contig)
                src["extra_inputs"][f"referee_cache_{contig}_faa"] = prefix + ".proteins.faa"
                src["extra_inputs"][f"referee_cache_{contig}_blastp"] = prefix + ".blastp.tsv"
        if species == "human":
            src["truths"][(species, contig, "soto")] = bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
            u2 = npip_u2_truth(cfg) if contig == "chr16" else None
            if u2:
                src["truths"][(species, contig, "npip_u2")] = u2
    return src


def compara_contig_truth(cfg: dict, contig: str, wdir: Path) -> Path | None:
    """Dev scope: the genome-wide Compara families at Primates (compara_families) restricted to one human contig (a
    family keeps its members on the contig; the scorer keeps families with >= 2 of them). None when the Compara table
    or the human annotation cache is absent."""
    try:
        fams, _ = compara_families(cfg)
    except _o1.NotBuilt as e:
        print(f"[fig7] Compara families not available: {e}", file=sys.stderr)
        return None
    dst = wdir / f"human_{contig}_compara.tsv"
    if not figlib.fresh(dst, fams):
        _o1.filter_rows(fams, dst, "Contig", keep={contig})
    return dst


def cached_sources(cfg: dict, wdir: Path) -> dict:
    """Dev scope, `fig7_dev_cached 1`: the per-contig products already in `wdir` (the runs of 2026-09-25), rescored
    against the current references without re-running any mode or the protein-homology BLASTP (their freshness is
    not re-checked; the tables say so). Missing products leave the arm out."""
    gdir = wdir / "gff"
    src = {"source": ("the cached per-contig runs of both modes in ${work}/fig7/current (2026-09-25; not re-run for "
                      "this table: fig7_dev_cached 1), rescored against the current references"),
           "arms": {}, "truths": {}, "genes": {}, "paf_source": {}, "extra_inputs": {}}
    for species, contig, status in substrates(cfg):
        genes = gdir / f"{species}_{contig}.gff"
        if not genes.exists():
            continue
        src["genes"][(species, contig)] = genes
        for mode, p in (("denovo", wdir / f"{species}_{contig}.denovo.fam.clusters.tsv"),
                        ("guided", wdir / f"{species}_{contig}.guided.clusters.tsv")):
            if p.exists():
                src["arms"][(species, contig, mode)] = (p, "cached")
        ref = Path(f"{wdir / f'{species}_{contig}_ref'}.families.tsv")
        if ref.exists():
            src["truths"][(species, contig, "referee")] = ref
        if species == "human":
            cf = compara_contig_truth(cfg, contig, wdir)
            if cf:
                src["truths"][(species, contig, "compara")] = cf
            src["truths"][(species, contig, "soto")] = bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
            u2 = npip_u2_truth(cfg) if contig == "chr16" else None
            if u2:
                src["truths"][(species, contig, "npip_u2")] = u2
        src["paf_source"][(species, contig)] = "cached"
    return src


def ensure_sources(cfg: dict, wdir: Path, force: bool = False) -> dict:
    """Regenerate both modes on every substrate with the current binaries and driver defaults (HEAVY, one step at a
    time; cached). Per substrate: GFF slice; protein-homology families (truth.py protein-homology, blastp); de novo (driver families on
    the restricted genome-wide assembly); guided (bodies -> minimap2 -> mcl_families; a recorded all-vs-all of the
    identical node set is reused)."""
    gdir = wdir / "gff"
    gdir.mkdir(parents=True, exist_ok=True)
    src = {"source": "regenerated with the current binaries and pipeline-driver defaults", "arms": {}, "truths": {},
           "genes": {}, "paf_source": {}, "extra_inputs": {}}
    full = {"human": cfg.get("human_ref_gff", DEFAULT_HUMAN_GFF), "gorilla": cfg["gorilla_ref_gff"]}
    subs = substrates(cfg)
    slices = {sp: gff_slices(full[sp], {c: gdir / f"{sp}_{c}.gff" for s_, c, _ in subs if s_ == sp}, force=force)
              for sp in ("human", "gorilla")}
    for species, contig, status in subs:
        genes = slices[species][contig]
        src["genes"][(species, contig)] = genes
        src["truths"][(species, contig, "referee")] = referee_build(cfg, species, contig, genes, wdir, force=force)
        if species == "human":
            cf = compara_contig_truth(cfg, contig, wdir)
            if cf:
                src["truths"][(species, contig, "compara")] = cf
            src["truths"][(species, contig, "soto")] = bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
            u2 = npip_u2_truth(cfg) if contig == "chr16" else None
            if u2:
                src["truths"][(species, contig, "npip_u2")] = u2
        src["arms"][(species, contig, "denovo")] = (denovo_families(cfg, species, contig, wdir, force=force), "current")
        g, how = guided_families(cfg, species, contig, genes, Path(full[species]), wdir, force=force)
        src["arms"][(species, contig, "guided")] = (g, "current")
        src["paf_source"][(species, contig)] = how
    return src


# ================================================================ genome-wide (every sample)
# docs/PREREG_genome_wide_families_2026-09-25.md sections 2, 3 and 5.
# main references (external): Compara and Soto (human; family_score), Liftoff copy pairs (every species; de novo only,
# scored by fig_family_recovery._liftoff_rows); protein-homology families only in the supplement, on request
GW_TRUTHS = {"human": ["compara", "soto"]}          # every other species: none scored by family_score in the main figure
GW_TRUTH_LABEL = {"protein_homology": "protein-homology families (secondary)",
                  "compara": "Ensembl Compara families (duplications within primates)",
                  "soto": "Soto et al. 2025 (not independent of the exon threshold)",
                  "npip_u2": "NPIP reference set (chr16, development)",
                  "liftoff": "Liftoff copy pairs (record, extra copy; both loci read-supported)"}
GW_TRUTH_SHORT = {"protein_homology": "protein homology", "compara": "Compara, primates", "soto": "Soto 2025",
                  "npip_u2": "NPIP set, chr16", "liftoff": "Liftoff copies"}
# Ensembl species-tree levels of the human lineage, youngest first (hsapiens_paralog_subtype); an unknown value stops
COMPARA_LEVELS = ["Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini",
                  "Primates", "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia", "Amniota",
                  "Tetrapoda", "Sarcopterygii", "Euteleostomi", "Gnathostomata", "Vertebrata", "Chordata", "Bilateria",
                  "Opisthokonta"]
COMPARA_HEADLINE = "Primates"
COMPARA_SWEEP = ["Hominidae", "Catarrhini", "Primates", "Eutheria", "Vertebrata", "Opisthokonta"]
GUIDED_SHARD_BP = 10_000_000       # query bases per wrapper shard of the genome-wide guided all-vs-all
GW_UNIT_VERSION = "1"


def scope(cfg: dict) -> str:
    s = str(cfg.get("fig7_scope", "genome")).strip().lower()
    if s not in ("genome", "dev"):
        raise ValueError(f"fig7_scope {s!r}: expected genome or dev")
    return s


def gw_truths(species: str, protein_homology: bool = False) -> list:
    """family_score references of a species' arms: the main ones, plus the secondary protein-homology families when
    requested (`fig7_protein_homology 1`, supplement only)."""
    return GW_TRUTHS.get(species, []) + (["protein_homology"] if protein_homology else [])


def _write_truth(dst: Path, rows) -> Path:
    """Gene Name / Family ID / Contig (family_score's --chrom ALL contig-keyed truth), written atomically."""
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    with open(tmp, "w") as fo:
        fo.write("Gene Name\tFamily ID\tContig\n")
        for c, g, f in rows:
            fo.write(f"{g}\t{f}\t{c}\n")
    tmp.replace(dst)
    return dst


def _gene_names_by_contig(cfg: dict, species: str) -> dict:
    """{name: [contigs where a gene / pseudogene record carries that Name=]} (genes.tsv, file order)."""
    out = collections.defaultdict(list)
    with open(_o1.annotation_cache(cfg, species)["genes"]) as fh:
        next(fh)
        for ln in fh:
            c, n = ln.split("\t", 2)[:2]
            if c not in out[n]:
                out[n].append(c)
    return out


def soto_gw(cfg: dict) -> Path:
    """Soto 2025 S1C as a contig-keyed truth: each gene's FIRST family (truth_table), on every contig where a RefSeq
    gene record carries that name — exactly what family_score --chrom ALL does with a name-only table, written out so
    that substrates can drop contigs."""
    src = bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
    d = _o1.species_dir(cfg, "human")
    out = d / "soto.families.tsv"
    ann = _o1.annotation_cache(cfg, "human")["genes"]
    if figlib.fresh(out, src, ann):
        return out
    where = _gene_names_by_contig(cfg, "human")
    return _write_truth(out, [(c, g, f) for g, f in truth_table(src) for c in where.get(g, [])])


def npip_gw(cfg: dict) -> Path | None:
    """The NPIP reference set (U2) on human chr16, contig-keyed (an inset: chr16 is the development chromosome)."""
    src = npip_u2_truth(cfg)
    if src is None:
        return None
    out = _o1.species_dir(cfg, "human") / "npip_u2.families.tsv"
    if figlib.fresh(out, src, _o1.annotation_cache(cfg, "human")["genes"]):
        return out
    where = _gene_names_by_contig(cfg, "human")
    return _write_truth(out, [("chr16", g, f) for g, f in truth_table(src) if "chr16" in where.get(g, [])])


def compara_families(cfg: dict, level: str = COMPARA_HEADLINE) -> tuple[Path, dict]:
    """Ensembl Compara families: connected components of the Compara pairs whose duplication node (subtype) is at or
    below `level`, over protein-coding RefSeq CHM13 genes on the chromosome Compara gives (MT -> chrM). Within one
    reconciled gene tree this relation is ultrametric, so the components are the clades below duplication nodes of
    age <= level. Returns (table, stats)."""
    if level not in COMPARA_LEVELS:
        raise ValueError(f"Compara level {level!r} is not one of {COMPARA_LEVELS}")
    compara = _o1.compara_gw(cfg)
    ann = _o1.annotation_cache(cfg, "human")["genes"]
    d = _o1.species_dir(cfg, "human")
    out = d / f"compara.{level.replace(' ', '_')}.families.tsv"
    stats_p = out.with_suffix(".stats.json")
    if figlib.fresh(out, compara, ann) and stats_p.exists():
        return out, json.loads(stats_p.read_text())
    pc = {(c, n) for (c, n), v in _o1.read_genes(ann).items() if v[4] == "protein_coding"}
    lim = COMPARA_LEVELS.index(level)
    par: dict = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    seen_pc, unmatched, rows_used, unknown = set(), set(), 0, collections.Counter()
    with open(compara) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 7 or not f[0] or not f[1] or f[0] == f[1]:
                continue
            if f[4] not in COMPARA_LEVELS:
                unknown[f[4]] += 1
                continue
            a = ("chrM" if f[6] == "MT" else f"chr{f[6]}", f[0])
            b = ("chrM" if f[5] == "MT" else f"chr{f[5]}", f[1])
            for x in (a, b):
                (seen_pc if x in pc else unmatched).add(x)
            if COMPARA_LEVELS.index(f[4]) > lim or a not in pc or b not in pc:
                continue
            rows_used += 1
            ra, rb = find(a), find(b)
            if ra != rb:
                par[max(ra, rb)] = min(ra, rb)
    if unknown:
        raise RuntimeError(f"{compara}: unknown Compara subtype value(s) {dict(unknown)}; extend COMPARA_LEVELS in "
                           "the right species-tree order first")
    comps = collections.defaultdict(list)
    for x in par:
        comps[find(x)].append(x)
    fams = sorted((sorted(v) for v in comps.values() if len(v) >= 2), key=lambda v: v[0])
    _write_truth(out, [(c, g, f"CF{i}") for i, fam in enumerate(fams) for c, g in fam])
    stats = {"level": level, "pair_rows_used": rows_used, "families": len(fams),
             "genes": sum(len(v) for v in fams), "compara_genes_protein_coding_matched": len(seen_pc),
             "compara_genes_unmatched_or_not_protein_coding": len(unmatched)}
    stats_p.write_text(json.dumps(stats))
    return out, stats


def gw_truth_table(cfg: dict, species: str, truth: str) -> Path | None:
    if truth == "protein_homology":
        return _o1.protein_homology_families(cfg, species)
    if truth == "compara":
        return compara_families(cfg)[0]
    if truth == "soto":
        return soto_gw(cfg)
    if truth == "npip_u2":
        return npip_gw(cfg)
    raise ValueError(truth)


def guided_gw(cfg: dict, species: str, budget) -> Path:
    """Guided mode genome-wide, one run per species: every gene / pseudogene record (genes_only.gff), bodies by
    samtools faidx -r, the resumable all-vs-all (tools/mm2_shard.sh; flags `fig7_guided_flags`, default the de novo
    mode's), mcl_families --paf --gff <full GFF> --min-exonic-bp 1 --min-shared-exon-frac 0.60."""
    flags = str(cfg.get("fig7_guided_flags", "denovo")).strip()
    if flags not in GUIDED_MM2:
        raise ValueError(f"fig7_guided_flags {flags!r}: expected one of {list(GUIDED_MM2)}")
    gff = _o1.annotation_gff(cfg, species)
    ann = _o1.annotation_cache(cfg, species)
    fasta = _o1.registry(cfg)[_o1.species_rep(cfg, species)]["fasta"]
    d = _o1.species_dir(cfg, species) / "guided"
    d.mkdir(exist_ok=True)
    regs = set()
    with open(ann["genes_gff"]) as fh:
        for line in fh:
            f = line.split("\t", 5)
            if len(f) > 4 and f[2] in ("gene", "pseudogene"):
                regs.add(f"{f[0]}:{f[3]}-{f[4]}")
    regions = sorted(regs, key=lambda x: x.encode())
    reg = d / "regions.txt"
    text = "\n".join(regions) + "\n"
    if not reg.exists() or reg.read_text() != text:
        reg.write_text(text)
    prefix = d / f"guided.{flags}"
    clusters = Path(f"{prefix}.clusters.tsv")
    mcl = Path(cfg["bin"]) / "mcl_families"
    key = json.dumps({"regions": hashlib.md5(text.encode()).hexdigest(), "gff": figlib.file_fingerprint(gff),
                      "fasta": figlib.file_fingerprint(fasta), "flags": GUIDED_MM2[flags],
                      "mcl_families": figlib.file_fingerprint(mcl)}) + "\n"
    kf = Path(f"{prefix}.key")
    if _o1._key_ok(kf, key) and clusters.exists():
        return clusters
    paf = Path(f"{prefix}.paf")
    bodies = d / "bodies.fa"
    bkey = d / "bodies.fa.key"
    btext = f"{hashlib.md5(text.encode()).hexdigest()} {figlib.file_fingerprint(fasta)}\n"
    if not paf.exists() and not (bodies.exists() and _o1._key_ok(bkey, btext)):
        _o1.need(budget, 240, f"{species} guided bodies (samtools faidx, {len(regions):,} regions)")
        figlib.run(f"samtools faidx {fasta} -r {reg} > {bodies}.tmp && mv {bodies}.tmp {bodies}",
                   log=d / "faidx.log")
        bkey.write_text(btext)
    _o1.run_sharded_paf(cfg, paf, GUIDED_MM2[flags].split(), bodies, budget, f"{species} guided all-vs-all",
                        d / f"guided.{flags}.mm2.log", bp=GUIDED_SHARD_BP, guided_recipe=flags == "recipe")
    _o1.need(budget, 300, f"{species} guided mcl_families")
    figlib.run([str(mcl), "--paf", str(paf), "--gff", str(gff), "--min-exonic-bp", "1", "--min-shared-exon-frac",
                "0.60", "--out", str(prefix)], log=d / f"guided.{flags}.mcl.log")
    kf.write_text(key)
    bodies.unlink(missing_ok=True)
    return clusters


_FS_PAIR = re.compile(r"\| pairwise \| truth pairs (\d+) \| predicted pairs (\d+) \| tp (\d+)")


def fs_score(cfg: dict, clusters: Path, genes_gff: Path, truth: Path, label: str, wdir: Path, keep=None,
             drop=()) -> dict:
    """`family_score --chrom ALL --per-family --pairwise` on the clusters and the contig-keyed truth, both restricted
    to the kept contigs (every contig but `drop`, or only `keep`). Pooled numbers are recomputed exactly from the
    per-family rows and checked against the printed 3-decimal line."""
    wdir.mkdir(parents=True, exist_ok=True)
    tag = hashlib.md5(json.dumps([str(clusters), str(truth), sorted(keep) if keep is not None else None,
                                  sorted(drop)]).encode()).hexdigest()[:12]
    c_in, t_in = Path(clusters), Path(truth)
    if keep is not None or drop:
        c_in = _o1.filter_rows(clusters, wdir / f"{tag}.clusters.tsv", "chrom", drop, keep)
        t_in = _o1.filter_rows(truth, wdir / f"{tag}.truth.tsv", "Contig", drop, keep)
    pf, log = wdir / f"{tag}.per_family.tsv", wdir / f"{tag}.family_score.txt"
    figlib.run([bin_path(cfg, "family_score"), "--clusters", str(c_in), "--gff", str(genes_gff), "--soto", str(t_in),
                "--chrom", "ALL", "--label", label, "--per-family", str(pf), "--pairwise"], log=log)
    pooled = pair = None
    for line in open(log):
        m = _FS.search(line)
        if m:
            pooled = m
        m = _FS_PAIR.search(line)
        if m:
            pair = m
    per = []
    if pf.exists():
        with open(pf) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for ln in fh:
                r = dict(zip(hdr, ln.rstrip("\n").split("\t")))
                per.append(r)
    out = {"truth_families": len(per), "truth_genes": sum(int(r["n_truth"]) for r in per),
           "per_family": per, "log": str(log)}
    if pooled is None:
        out.update(empty=True, sens=0.0, prec=None, f=0.0, clusters_scored=0, collapsed=None, no_locus=None,
                   matched=0, pred_members=0, pair_truth=None, pair_pred=None, pair_tp=None, pair_sens=None,
                   pair_prec=None, exact=0, partial=0, missed=len(per))
        return out
    matched = sum(int(r["hit"]) for r in per)
    pred = sum(int(r["n_pred"]) for r in per if int(r["hit"]) > 0)
    tg = out["truth_genes"]
    sens = matched / tg if tg else 0.0
    prec = matched / pred if pred else 0.0
    f = 2 * sens * prec / (sens + prec) if sens + prec else 0.0
    for got, i in ((sens, 4), (prec, 5), (f, 6)):
        if abs(got - float(pooled.group(i))) > 0.0006:
            raise RuntimeError(f"{log}: recomputed {got:.4f} != family_score {pooled.group(i)}")
    if int(pooled.group(1)) != len(per) or int(pooled.group(2)) != tg:
        raise RuntimeError(f"{log}: per-family rows ({len(per)} families, {tg} genes) != the pooled line")
    fv = [float(r["f"]) for r in per]
    out.update(empty=False, sens=sens, prec=prec, f=f, clusters_scored=int(pooled.group(3)),
               collapsed=int(pooled.group(7)), no_locus=int(pooled.group(8)), matched=matched, pred_members=pred,
               exact=sum(1 for x in fv if x >= 1.0 - 1e-9), partial=sum(1 for x in fv if 0 < x < 1.0 - 1e-9),
               missed=sum(1 for x in fv if x <= 0))
    if pair:
        pt, pp, tp = (int(pair.group(i)) for i in (1, 2, 3))
        out.update(pair_truth=pt, pair_pred=pp, pair_tp=tp, pair_sens=tp / pt if pt else None,
                   pair_prec=tp / pp if pp else None)
    else:
        out.update(pair_truth=None, pair_pred=None, pair_tp=None, pair_sens=None, pair_prec=None)
    return out


def clusters_profile(clusters: Path) -> tuple[int, int]:
    """(predicted families, loci in them) of a clusters.tsv."""
    ids, n = set(), 0
    with open(clusters) as fh:
        next(fh)
        for ln in fh:
            ids.add(ln.split("\t", 1)[0])
            n += 1
    return len(ids), n


def annotation_contigs(cfg: dict, species: str) -> list:
    """Contigs with >= 1 gene record, natural order (the per-contig breakdown)."""
    cs = set()
    with open(_o1.annotation_cache(cfg, species)["genes"]) as fh:
        next(fh)
        for ln in fh:
            cs.add(ln.split("\t", 1)[0])
    return sorted(cs, key=assembly._natural)
