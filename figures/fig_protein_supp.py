"""Supplementary figure S-P — the manual extra-sensitive step (protein attachment; `tools/protein_attach.py`, driver
stage `families-protein`, never in `all`): what it adds to the default de novo families, scored against EXTERNAL
references (docs/PREREG_protein_attach_2026-09-25.md, section 3).

Protein is not part of the default family definition. The default de novo families (driver `families`: seeded
assembly loci -> one representative per locus -> mcl_families --from-gtf, exon-sum >= 0.60, MCL 2.8) are read and
never modified; the step can only add a locus the RNA rule left out to ONE existing family, and never merges two.

References (section 3.2): Ensembl Compara release 116 paralogues at any duplication age (primary, human) and within
primates (secondary, human); Liftoff copies (the Fig. 8 self-lift; certifies, never falsifies; every species); the
annotation's protein-homology families only as a secondary reference labelled circular (a protein step agrees with a
protein-built truth partly by construction; registers 1029, 1031).

Substrates: the development tables use the de novo families of the Fig. 7 runs (the genome-wide assembly of each
species restricted to one chromosome; `${work}/fig7/current/<species>_<contig>.denovo.fam.*`, cfg
`figS_protein_families_dir`): human chr16 (development), human chr6 (held out, never used), gorilla chr20
(NC_073244.2; the contig the seeding rule was chosen on) and chr10 (NC_073234.2; held out, never used). The step is
light there and runs from here (outputs in cfg `figS_protein_work`, default `${work}/figS_protein`). Genome-wide rows
appear per sample once the driver's `families-protein` stage has written `${work}/runs/<id>/<id>.fam_protein.params.tsv`
(this module never runs a genome-wide stage).

    python3 figures/fig_protein_supp.py data [--inputs F] [--force]   build the tables (dev runs are light)
    python3 figures/fig_protein_supp.py plot                          render figures/out/figS_protein.*
    python3 figures/fig_protein_supp.py summary                       the numbers the caption and PREREG quote
"""
from __future__ import annotations

import collections
import csv
import math
import os
import sys
from pathlib import Path

import figlib

META = {
    "id": "figS_protein",
    "title": "Supplementary: the manual extra-sensitive step (protein attachment) scored against external references",
    "claim": ("Protein is not in the default family definition. The optional step attaches a locus the RNA rule "
              "left out to one existing family only when its protein is at least as close to the family as the "
              "family's own loosest member (and >= 60% identical); scored against Ensembl Compara and Liftoff, "
              "never against protein-built families except as a labelled secondary reference."),
    "tables": ["figS_protein_substrates", "figS_protein_refs"],
}
T_SUB, T_REF = "figS_protein_substrates", "figS_protein_refs"
DEV_SUBSTRATES = [("human", "chr16", "development"), ("human", "chr6", "held out, never used"),
                  ("gorilla", "NC_073244.2", "development (seeding rule)"),
                  ("gorilla", "NC_073234.2", "held out, never used")]
CONTIG_LABEL = {"chr16": "chr16", "chr6": "chr6", "NC_073244.2": "chr20 (NC_073244.2)",
                "NC_073234.2": "chr10 (NC_073234.2)"}
REFS = ["compara_any", "compara_primates", "liftoff", "protein_homology"]
REF_LABEL = {"compara_any": "Ensembl Compara, paralogues of any age",
             "compara_primates": "Ensembl Compara, duplications within primates",
             "liftoff": "Liftoff copies (self-lift, ≥ 95% identity)",
             "protein_homology": "protein-homology families (secondary; circular for a protein step)"}
REF_ROLE = {"compara_any": "primary", "compara_primates": "secondary", "liftoff": "certifies only",
            "protein_homology": "secondary, circular"}
VERDICT_ORDER = ["attached", "ambiguous", "below_family_identity", "below_scope_floor", "below_coverage",
                 "family_uncalibrated"]
VERDICT_LABEL = {"attached": "attached", "ambiguous": "passes for ≥ 2 families (merge hint, not attached)",
                 "below_family_identity": "less identical than the family's loosest member",
                 "below_scope_floor": "below 60% protein identity",
                 "below_coverage": "hit covers < 30% of the longer protein",
                 "family_uncalibrated": "family's members share no protein homology"}
VERDICT_COLOR = {"attached": figlib.BLUE[650], "ambiguous": "#8a6fb0", "below_family_identity": figlib.BLUE[300],
                 "below_scope_floor": figlib.BLUE[100], "below_coverage": "#bdbcb6", "family_uncalibrated": "#dddcd6"}
FIG_W = figlib.WIDTH_DOUBLE - 4 * figlib.MM
SUB_HEAD = ["species", "substrate", "scope", "exposure", "loci", "families", "member_loci", "families_calibrated",
            "loci_protein_ok", "loci_te_majority", "loci_orf_50_99", "candidates", "attached", "attached_families",
            "attached_nt_corroborated", "null_attached", "ambiguous", "below_family_identity", "below_scope_floor",
            "below_coverage", "family_uncalibrated", "no_hit", "overlaps_member", "proposed_merges"]
REF_HEAD = ["species", "substrate", "scope", "exposure", "reference", "role", "status", "members_judgeable",
            "members_true", "baseline_precision", "baseline_lo", "baseline_hi", "att_judgeable", "att_true",
            "att_false", "att_same_gene", "att_unjudgeable", "att_precision", "att_lo", "att_hi", "missing",
            "missing_protein", "missing_base_hit", "missing_attached", "sensitivity", "merges", "merges_judgeable",
            "merges_true"]
PROVISIONAL = ("provisional: development scope (the Fig. 7 de novo families of four chromosomes); the genome-wide rows "
               "of every sample appear once the driver's families-protein stage has run")


def _tools():
    p = str(figlib.REPO / "tools")
    if p not in sys.path:
        sys.path.insert(0, p)
    import protein_attach
    return protein_attach


# ================================================================ statistics
def wilson(k: int, n: int, z: float = 1.959964) -> tuple:
    if not n:
        return (None, None, None)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, c - h), min(1.0, c + h))


# ================================================================ substrates and running the step (dev scope)
def work(cfg: dict) -> Path:
    w = cfg.get("figS_protein_work")
    d = Path(w) if w else figlib.work_dir(cfg, "figS_protein")
    d.mkdir(parents=True, exist_ok=True)
    return d


def substrates(cfg: dict) -> list[dict]:
    """Every substrate with its families prefix, its protein-step prefix, FASTA and exposure."""
    fam_dir = Path(cfg.get("figS_protein_families_dir") or
                   Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "fig7" / "current")
    out = []
    for sp, contig, exposure in DEV_SUBSTRATES:
        fam = fam_dir / f"{sp}_{contig}.denovo.fam"
        out.append({"species": sp, "substrate": contig, "scope": "dev", "exposure": exposure, "fam": str(fam),
                    "prot": str(work(cfg) / f"{sp}_{contig}"), "fasta": cfg.get(f"{sp}_fasta"), "contigs": {contig}})
    try:
        import samples
        reg = samples.registry(cfg)
    except Exception as e:  # noqa: BLE001  (no registry on this machine: dev rows only)
        print(f"[figS_protein] no sample registry ({e}); development rows only", file=sys.stderr)
        reg = {}
    for sid, row in reg.items():
        pre = Path(cfg.get("work", "/mnt/linuxdisk/tmp/rustle_figures")) / "runs" / sid / sid
        if Path(f"{pre}.fam_protein.params.tsv").exists():
            exp = "held out (closed-loop verdict sample)" if sid in ("human_testis", "gorilla_KB3781", "chimp_PTR",
                                                                   "orangutan_PPY") else "development sample"
            out.append({"species": row["species"], "substrate": sid, "scope": "genome", "exposure": exp,
                        "fam": f"{pre}.fam", "prot": f"{pre}.fam_protein", "fasta": row["fasta"], "contigs": None})
    return out


def run_dev(s: dict, force: bool = False):
    """The step on one development substrate (light: seconds to a minute; BLASTP of a few thousand proteins)."""
    tool = figlib.REPO / "tools" / "protein_attach.py"
    params = Path(s["prot"] + ".params.tsv")
    srcs = [s["fam"] + ".clusters.tsv", s["fam"] + ".loci.gff3", s["fam"] + ".loci.tsv", tool]
    if not Path(s["fam"] + ".clusters.tsv").exists():
        raise FileNotFoundError(f"{s['fam']}.clusters.tsv absent: run `make.py data fig7` (development scope) first")
    if force or not figlib.fresh(params, *srcs):
        figlib.run([sys.executable or "python3", str(tool), "--fam", s["fam"], "--fasta", s["fasta"],
                    "--out", s["prot"], "--threads", "4"], log=Path(s["prot"] + ".log"))


# ================================================================ gene labels and references
def read_genes(path: Path, contigs) -> dict:
    """{contig: [(start0, end, name, biotype, [exon (s, e) 0-based])]} from the annotation cache genes.tsv."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        next(fh)
        for ln in fh:
            c, n, t, bt, s1, e, ex = ln.rstrip("\n").split("\t")
            if contigs is not None and c not in contigs:
                continue
            ivs = [tuple(map(int, x.split("-"))) for x in ex.split(",") if x]
            if ivs:
                out[c].append((int(s1) - 1, int(e), n, bt, ivs))
    return out


class GeneLabeler:
    """A locus takes the gene whose merged exons cover the most of its representative's exonic bases, if >= 0.50
    (PREREG 3.1; register 637's large-gene attractor is why the cut is on the LOCUS's bases)."""
    B = 1_000_000

    def __init__(self, genes: dict):
        self.bins = collections.defaultdict(list)
        for c, gl in genes.items():
            for g in gl:
                for k in range(g[0] // self.B, g[1] // self.B + 1):
                    self.bins[(c, k)].append(g)

    def label(self, locus: dict):
        ex = [(a - 1, b) for a, b in locus["exons"]]
        tot = sum(b - a for a, b in ex)
        if not tot:
            return None
        best, bcov = None, 0.0
        seen = set()
        for k in range(ex[0][0] // self.B, ex[-1][1] // self.B + 1):
            for g in self.bins.get((locus["contig"], k), ()):
                if id(g) in seen or g[1] <= ex[0][0] or ex[-1][1] <= g[0]:
                    continue
                seen.add(id(g))
                ov = sum(max(0, min(b, y) - max(a, x)) for a, b in ex for x, y in g[4])
                if ov / tot > bcov:
                    best, bcov = g, ov / tot
        return (locus["contig"], best[2]) if best is not None and bcov >= 0.5 else None


def compara_relation(cfg: dict, level: str, wdir: Path) -> tuple[set, dict]:
    """(universe, {gene: group}) of Compara at `level` (components of the pairs at or below it), over protein-coding
    CHM13 genes (the construction of _o1_recovery.compara_families; universe = protein-coding genes in the table)."""
    import _o1
    import _o1_recovery as rec
    compara = _o1.compara_gw(cfg)
    ann = _o1.annotation_cache(cfg, "human")["genes"]
    pc = {(c, n) for (c, n), v in _o1.read_genes(ann).items() if v[4] == "protein_coding"}
    lim = rec.COMPARA_LEVELS.index(level)
    par: dict = {}

    def find(x):
        par.setdefault(x, x)
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    uni = set()
    with open(compara) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 7 or not f[0] or not f[1] or f[0] == f[1] or f[4] not in rec.COMPARA_LEVELS:
                continue
            a = ("chrM" if f[6] == "MT" else f"chr{f[6]}", f[0])
            b = ("chrM" if f[5] == "MT" else f"chr{f[5]}", f[1])
            for x in (a, b):
                if x in pc:
                    uni.add(x)
                    find(x)
            if rec.COMPARA_LEVELS.index(f[4]) <= lim and a in pc and b in pc:
                ra, rb = find(a), find(b)
                if ra != rb:
                    par[max(ra, rb)] = min(ra, rb)
    return uni, {x: find(x) for x in uni}


def protein_homology_relation(fam_file: Path, contig: str | None) -> tuple[set, dict]:
    """(universe, {gene: family}) of a protein-homology families table (per-contig two-column or genome-wide
    three-column); universe = the genes of its .proteins.faa when present (a gene in no family is a singleton)."""
    grp, uni = {}, set()
    with open(fam_file) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            key = (f[2], f[0]) if len(head) >= 3 else (contig, f[0])
            grp[key] = f[1]
            uni.add(key)
    faa = Path(str(fam_file).replace(".families.tsv", ".proteins.faa"))
    if faa.exists() and contig is not None:
        for ln in open(faa):
            if ln.startswith(">"):
                uni.add((contig, ln[1:].strip()))
    return uni, {g: grp.get(g, ("single",) + g) for g in uni}


def liftoff_labels(cfg: dict, species: str, loci: list[dict]):
    """{locus idx: {(source_id, placement key)}} for loci >= 0.50 covered by a Liftoff placement, or None when the
    species' self-lift is not merged yet."""
    import _liftoff
    try:
        path = _liftoff.loci_path(cfg, species)
    except Exception as e:  # noqa: BLE001
        print(f"[figS_protein] Liftoff {species}: {e}", file=sys.stderr)
        return None
    if path is None:
        return None
    rows = _liftoff.read_loci(path)
    bins = collections.defaultdict(list)
    for r in rows:
        s0, e0 = int(r["start0"]), int(r["end"])
        for k in range(s0 // 1_000_000, e0 // 1_000_000 + 1):
            bins[(r["contig"], k)].append(r)
    out = {}
    for d in loci:
        ex = [(a - 1, b) for a, b in d["exons"]]
        tot = sum(b - a for a, b in ex)
        labs = set()
        for k in range(ex[0][0] // 1_000_000, ex[-1][1] // 1_000_000 + 1):
            for r in bins.get((d["contig"], k), ()):
                ov = sum(max(0, min(b, y) - max(a, x)) for a, b in ex for x, y in r["iv"])
                if tot and ov / tot >= 0.5:
                    labs.add((r["source_id"], (r["contig"], r["start0"], r["end"])))
        if labs:
            out[d["idx"]] = labs
    return out


# ================================================================ scoring one substrate against one reference
def score(labels: dict, fam_of: dict, members_by_fam: dict, overlapping: set, status: dict, cand: dict,
          hits: dict, merges: list) -> dict:
    """labels: {locus idx: set of (group, unit)}; related = same group, different unit; same gene = shared unit."""
    def related(x, y):
        lx, ly = labels.get(x), labels.get(y)
        return bool(lx and ly and any(g1 == g2 and u1 != u2 for g1, u1 in lx for g2, u2 in ly))

    def units(x):
        return {u for _, u in labels.get(x, ())}
    out = collections.Counter()
    # baseline: the RNA families' own members (PREREG 3.4)
    for F, mem in members_by_fam.items():
        for x in mem:
            if x not in labels:
                continue
            others = [y for y in mem if y != x and y in labels and not (units(x) & units(y))]
            if others:
                out["members_judgeable"] += 1
                out["members_true"] += any(related(x, y) for y in others)
    held_units = set()
    for x in fam_of:
        held_units |= units(x)
    group_fams = collections.defaultdict(set)          # group -> {(family, unit)}
    for x, F in fam_of.items():
        for g, u in labels.get(x, ()):
            group_fams[g].add((F, u))
    # attachments (PREREG 3.3)
    for u, (v, F) in cand.items():
        if v != "attached":
            continue
        mem = members_by_fam[F]
        if any(units(u) & units(m) for m in mem):
            out["att_same_gene"] += 1
        elif u in labels and any(m in labels for m in mem):
            out["att_judgeable"] += 1
            if any(related(u, m) for m in mem):
                out["att_true"] += 1
            else:
                out["att_false"] += 1
        else:
            out["att_unjudgeable"] += 1
    # missing members (PREREG 3.5)
    for u in status:
        if u in fam_of or u in overlapping or u not in labels or (units(u) & held_units):
            continue
        rf = {F for g, un in labels[u] for F, uu in group_fams.get(g, ()) if uu != un}
        if not rf:
            continue
        out["missing"] += 1
        out["missing_protein"] += status[u] == "ok"
        out["missing_base_hit"] += any(hits.get((u, F)) for F in rf)
        v = cand.get(u)
        out["missing_attached"] += bool(v and v[0] == "attached" and v[1] in rf)
    for A, B in merges:
        ma = [x for x in members_by_fam.get(A, ()) if x in labels]
        mb = [x for x in members_by_fam.get(B, ()) if x in labels]
        out["merges"] += 1
        if ma and mb:
            out["merges_judgeable"] += 1
            out["merges_true"] += any(related(x, y) for x in ma for y in mb)
    return out


def load_step(s: dict) -> dict:
    """The protein step's outputs of one substrate, keyed by locus idx."""
    pa = _tools()
    loci = pa.read_loci(s["fam"] + ".loci.gff3")
    by_id = {d["id"]: d["idx"] for d in loci}
    fam_of_span, fold, _ = pa.read_families(s["fam"] + ".clusters.tsv", s["fam"] + ".loci.tsv")
    fam_of = {}
    for d in loci:
        k = (d["contig"], d["start"], d["end"])
        k = k if k in fam_of_span else fold.get(k)
        if k in fam_of_span:
            fam_of[d["idx"]] = fam_of_span[k]
    members_by_fam = collections.defaultdict(list)
    for i, F in fam_of.items():
        members_by_fam[F].append(i)
    mindex = pa.ExonIndex([loci[i] for i in fam_of])
    overlapping = {d["idx"] for d in loci if d["idx"] not in fam_of and mindex.hits(d)}
    status = {int(r["idx"]): r["status"] for r in csv.DictReader(open(s["prot"] + ".orfs.tsv"), delimiter="\t")}
    cand = {by_id[r["locus"]]: (r["verdict"], r["best_family"])
            for r in csv.DictReader(open(s["prot"] + ".candidates.tsv"), delimiter="\t")}
    hits = {}
    for r in csv.DictReader(open(s["prot"] + ".candidate_hits.tsv"), delimiter="\t"):
        hits[(by_id[r["locus"]], r["family"])] = r["base_identity"] != "NA"
    merges = [(r["family_a"], r["family_b"]) for r in csv.DictReader(open(s["prot"] + ".merges.tsv"), delimiter="\t")
              if r["via"] == "member_pair"]
    params = dict(ln.rstrip("\n").split("\t", 1) for ln in open(s["prot"] + ".params.tsv"))
    return {"loci": loci, "fam_of": fam_of, "members_by_fam": members_by_fam, "overlapping": overlapping,
            "status": status, "cand": cand, "hits": hits, "merges": merges, "params": params}


def gene_labels(cfg: dict, species: str, loci: list[dict], contigs) -> dict:
    import _o1
    genes = read_genes(_o1.annotation_cache(cfg, species)["genes"], contigs)
    lab = GeneLabeler(genes)
    return {d["idx"]: g for d in loci if (g := lab.label(d)) is not None}


def ref_labels(cfg: dict, s: dict, ref: str, glab: dict, loci: list[dict], cache: dict):
    """{locus idx: {(group, unit)}} for one reference, or (None, reason) when it is not available here."""
    if ref.startswith("compara"):
        if s["species"] != "human":
            return None, "Compara export in hand is human only"
        level = "Opisthokonta" if ref == "compara_any" else "Primates"
        if level not in cache:
            try:
                cache[level] = compara_relation(cfg, level, work(cfg))
            except Exception as e:  # noqa: BLE001
                return None, f"not available ({e})"
        uni, grp = cache[level]
        return {i: {(grp[g], g)} for i, g in glab.items() if g in uni}, "ok"
    if ref == "liftoff":
        lab = liftoff_labels(cfg, s["species"], loci)
        if lab is None:
            return None, "not available (Liftoff self-lift not finished for this species)"
        return lab, "ok"
    if ref == "protein_homology":
        if s["scope"] == "dev":
            fam_dir = Path(s["fam"]).parent
            p = fam_dir / f"{s['species']}_{s['substrate']}_ref.families.tsv"
            contig = s["substrate"]
        else:
            import _o1
            try:
                p, contig = _o1.protein_homology_families(cfg, s["species"]), None
            except Exception as e:  # noqa: BLE001
                return None, f"not available ({str(e).split(':')[0]})"
        if not p.exists():
            return None, "not available (no protein-homology families for this substrate)"
        uni, grp = protein_homology_relation(p, contig)
        # gene labels are (contig, RefSeq name); protein-homology ids are the gene symbols
        return {i: {(grp[g], g)} for i, g in glab.items() if g in uni}, "ok"
    raise ValueError(ref)


# ================================================================ build
def build(cfg: dict, data_dir: Path, force: bool = False):
    sub_rows, ref_rows, inputs, cache = [], [], {}, {}
    subs = substrates(cfg)
    for s in subs:
        if s["scope"] == "dev":
            if not s["fasta"]:
                print(f"[figS_protein] {s['species']}: no FASTA key; skipped", file=sys.stderr)
                continue
            run_dev(s, force)
        st = load_step(s)
        p = st["params"]
        inputs[f"{s['species']}_{s['substrate']}"] = s["prot"] + ".params.tsv"
        sub_rows.append([s["species"], s["substrate"], s["scope"], s["exposure"], p["loci"], p["families"],
                         p["member_loci"], p["families_calibrated"], p["loci_protein_ok"], p["loci_te_majority"],
                         p["loci_orf_50_99"], p["candidates"], p["verdict_attached"], p["attached_families"],
                         p["attached_nt_corroborated"], p["null_attached"], p["verdict_ambiguous"],
                         p["verdict_below_family_identity"], p["verdict_below_scope_floor"],
                         p["verdict_below_coverage"], p["verdict_family_uncalibrated"], p["verdict_no_hit"],
                         p["verdict_overlaps_member"], p["proposed_merges_member_pairs"]])
        glab = gene_labels(cfg, s["species"], st["loci"], s["contigs"])
        for ref in REFS:
            lab, why = ref_labels(cfg, s, ref, glab, st["loci"], cache)
            base = [s["species"], s["substrate"], s["scope"], s["exposure"], ref, REF_ROLE[ref]]
            if lab is None:
                ref_rows.append(base + [why] + [None] * (len(REF_HEAD) - 7))
                continue
            o = score(lab, st["fam_of"], st["members_by_fam"], st["overlapping"], st["status"], st["cand"],
                      st["hits"], st["merges"])
            bp, blo, bhi = wilson(o["members_true"], o["members_judgeable"])
            ap, alo, ahi = wilson(o["att_true"], o["att_judgeable"])
            ref_rows.append(base + ["ok", o["members_judgeable"], o["members_true"], bp, blo, bhi,
                                    o["att_judgeable"], o["att_true"], o["att_false"], o["att_same_gene"],
                                    o["att_unjudgeable"], ap, alo, ahi, o["missing"], o["missing_protein"],
                                    o["missing_base_hit"], o["missing_attached"],
                                    (o["missing_attached"] / o["missing"]) if o["missing"] else None,
                                    o["merges"], o["merges_judgeable"], o["merges_true"]])
    notes = [PROVISIONAL] if all(r[2] == "dev" for r in sub_rows) else []
    gen = "figures/fig_protein_supp.py (tools/protein_attach.py; docs/PREREG_protein_attach_2026-09-25.md)"
    figlib.write_table(T_SUB, SUB_HEAD, sub_rows, generator=gen, inputs=inputs, notes=notes, data_dir=data_dir)
    figlib.write_table(T_REF, REF_HEAD, ref_rows, generator=gen, inputs=inputs, notes=notes, data_dir=data_dir)


# ================================================================ plot
def _num(x):
    return None if x in (None, "", "NA") else float(x)


def _row_label(r: dict) -> str:
    sp = "Human" if r["species"] == "human" else r["species"].capitalize()
    sub = CONTIG_LABEL.get(r["substrate"], r["substrate"])
    return f"{sp} {sub}\n{r['exposure']}"


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import textwrap

    subs = figlib.read_table(T_SUB, data_dir)
    refs = figlib.read_table(T_REF, data_dir)
    fig = plt.figure(figsize=(FIG_W, 5.6))
    n = len(subs)
    ys = list(range(n))

    def frame(ax, xlabel, title, letter, x_letter):
        ax.set_xlabel(xlabel, labelpad=1.5, fontsize=6.0)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
        ax.set_title(title, fontsize=6.4, loc="left")
        figlib.panel_label(ax, letter, x=x_letter)

    # a: what the step adds, per substrate
    ax = fig.add_axes([0.15, 0.665, 0.25, 0.28])
    top = max([int(r["attached"]) for r in subs] + [int(r["null_attached"]) for r in subs] + [1])
    for y, r in zip(ys, subs):
        att, cand, null, nf = int(r["attached"]), int(r["candidates"]), int(r["null_attached"]), int(r["attached_families"])
        ax.barh(y, att, height=0.55, color=figlib.BLUE[650])
        ax.plot([null], [y], marker="o", markerfacecolor=figlib.SURFACE, markeredgecolor=figlib.INK, markersize=3.5,
                linestyle="none", clip_on=False)
        ax.text(top * 1.12, y, f"{att} of {cand:,} candidates attached\n({nf} famil{'y' if nf == 1 else 'ies'}); "
                f"null {null}; {r['proposed_merges']} merge proposals", va="center", fontsize=5.2, color=figlib.INK_2,
                linespacing=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels([_row_label(r) for r in subs], fontsize=5.6, linespacing=1.0)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0, top / 0.22)
    ax.set_xticks(range(0, top + 1, max(1, top // 3)))
    frame(ax, "Loci added to an existing family", "What the step adds (circle: shuffled null)", "a", -0.62)

    # b: fate of the candidates that hit a member's protein
    ax = fig.add_axes([0.15, 0.285, 0.25, 0.28])
    for y, r in zip(ys, subs):
        tot = sum(int(r[k]) for k in VERDICT_ORDER)
        left = 0.0
        for k in VERDICT_ORDER:
            v = int(r[k])
            if tot and v:
                ax.barh(y, v / tot, left=left, height=0.55, color=VERDICT_COLOR[k], edgecolor=figlib.SURFACE,
                        linewidth=0.5)
                if v / tot > 0.08:
                    ax.text(left + v / tot / 2, y, str(v), ha="center", va="center", fontsize=5.2,
                            color="white" if k in ("attached", "ambiguous") else figlib.INK)
                left += v / tot
        ax.text(1.02, y, f"n = {tot}", va="center", fontsize=5.2, color=figlib.INK_2)
    ax.set_yticks(ys)
    ax.set_yticklabels([_row_label(r) for r in subs], fontsize=5.6, linespacing=1.0)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "50", "100"])
    frame(ax, "Candidates with a protein hit to a family member (%)", "Why candidates were not attached", "b", -0.62)
    fig.legend(handles=[mpatches.Patch(color=VERDICT_COLOR[k], label=VERDICT_LABEL[k]) for k in VERDICT_ORDER],
               loc="upper left", bbox_to_anchor=(0.005, 0.205), ncol=2, fontsize=5.2, handlelength=1.0,
               handleheight=0.7, columnspacing=1.0, frameon=False)

    # c: precision against external references, the RNA members (baseline) vs the attached loci
    ax = fig.add_axes([0.67, 0.665, 0.25, 0.28])
    rows = [r for r in refs if r["status"] == "ok" and r["reference"] != "liftoff"]
    short = {"compara_any": "Compara, any age", "compara_primates": "Compara, primates",
             "protein_homology": "protein homology*"}
    for y, r in enumerate(rows):
        c = figlib.INK_3 if r["reference"] == "protein_homology" else figlib.BLUE[650]
        bp, blo, bhi = _num(r["baseline_precision"]), _num(r["baseline_lo"]), _num(r["baseline_hi"])
        if bp is not None:
            ax.plot([blo, bhi], [y - 0.14, y - 0.14], color=c, linewidth=0.8)
            ax.plot([bp], [y - 0.14], marker="o", markerfacecolor=figlib.SURFACE, markeredgecolor=c, markersize=3.3,
                    linestyle="none", clip_on=False)
        ax.text(1.03, y - 0.14, f"{r['members_true']}/{r['members_judgeable']}", va="center", fontsize=5.0,
                color=figlib.INK_2)
        ap, alo, ahi = _num(r["att_precision"]), _num(r["att_lo"]), _num(r["att_hi"])
        if ap is not None:
            ax.plot([alo, ahi], [y + 0.14, y + 0.14], color=c, linewidth=0.8)
            ax.plot([ap], [y + 0.14], marker="o", color=c, markersize=3.3, linestyle="none", clip_on=False)
            ax.text(1.03, y + 0.14, f"{r['att_true']}/{r['att_judgeable']}", va="center", fontsize=5.0,
                    color=figlib.INK_2)
        else:
            ax.text(1.03, y + 0.14, "none", va="center", fontsize=5.0, color=figlib.INK_3, style="italic")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{_row_label(r).splitlines()[0]}: {short[r['reference']]}" for r in rows], fontsize=5.3)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    frame(ax, "Precision over judgeable loci (95% Wilson interval)",
          "Open: the families' own members; filled: attached loci", "c", -0.98)

    # d: missing members (human, Compara)
    ax = fig.add_axes([0.67, 0.285, 0.25, 0.28])
    mall = [r for r in refs if r["status"] == "ok" and r["reference"].startswith("compara")]
    mrows = [r for r in mall if int(r["missing"] or 0) > 0]
    zero = [f"{_row_label(r).splitlines()[0]}, {short[r['reference']]}" for r in mall if int(r["missing"] or 0) == 0]
    keys = [("missing", "missing members"), ("missing_protein", "with a protein (>= 100 aa)"),
            ("missing_base_hit", "with a hit to the related family"), ("missing_attached", "attached to it")]
    shades = [figlib.BLUE[250], figlib.BLUE[350], figlib.BLUE[500], figlib.BLUE[700]]
    h = 0.8 / len(keys)
    for y, r in enumerate(mrows):
        for j, (k, _) in enumerate(keys):
            v = int(r[k] or 0)
            yy = y - 0.4 + h * (j + 0.5)
            ax.barh(yy, v, height=h * 0.9, color=shades[j])
            ax.text(v, yy, f" {v}", va="center", fontsize=4.8, color=figlib.INK_2)
    ax.set_yticks(range(len(mrows)))
    ax.set_yticklabels([f"{_row_label(r).splitlines()[0]}: {short[r['reference']]}" for r in mrows], fontsize=5.3)
    ax.set_ylim(len(mrows) - 0.5, -0.5)
    frame(ax, "Loci", "Missing members (human; Ensembl Compara)", "d", -0.98)

    ax.legend(handles=[mpatches.Patch(color=c, label=l) for c, (_, l) in zip(shades, keys)], fontsize=5.0,
              loc="upper left", bbox_to_anchor=(-0.02, -0.16), ncol=2, frameon=False, handlelength=0.9,
              columnspacing=1.0)

    by_status = collections.defaultdict(list)
    for r in refs:
        if r["status"] != "ok" and REF_LABEL[r["reference"]] not in by_status[r["status"]]:
            by_status[r["status"]].append(REF_LABEL[r["reference"]])
    miss = [f"Not shown, {' and '.join(v)}: {k.replace('not available (', '').rstrip(')')}" for k, v in
            sorted(by_status.items())]
    note = ("Protein is not part of the default family definition: this optional, manually run step adds a locus the "
            "RNA rule left out to one existing family and never merges two. Missing member: an unattached locus "
            "whose gene no family holds and which Compara relates to a member's gene. * Protein-homology families "
            "are built from protein homology, so agreement with them is partly by construction (secondary only). "
            + (f"No missing member at the primates level ({'; '.join(zero)}). " if zero else "")
            + " ".join(f"{m}." for m in miss))
    fig.text(0.005, 0.125, "\n".join(textwrap.wrap(note, 200)), fontsize=5.0, color=figlib.INK_2, va="top",
             linespacing=1.1)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "figS_protein", out_dir)
    plt.close(fig)
    return paths


# ================================================================ summary / CLI
def summary(data_dir: Path = figlib.DATA_DIR):
    for r in figlib.read_table(T_SUB, data_dir):
        print("\t".join(f"{k}={r[k]}" for k in SUB_HEAD))
    for r in figlib.read_table(T_REF, data_dir):
        print("\t".join(f"{k}={r[k]}" for k in REF_HEAD if r.get(k) not in (None, "")))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["data", "plot", "summary"])
    ap.add_argument("--inputs")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--data", default=str(figlib.DATA_DIR))
    ap.add_argument("--out", default=str(figlib.OUT_DIR))
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    a = ap.parse_args()
    if a.cmd == "data":
        cfg = figlib.load_inputs(a.inputs)
        for kv in a.set:
            k, _, v = kv.partition("=")
            cfg[k] = v
        build(cfg, Path(a.data), a.force)
    elif a.cmd == "plot":
        figlib.use_style()
        plot(Path(a.data), Path(a.out))
    else:
        summary(Path(a.data))
