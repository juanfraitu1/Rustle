"""The gene A: from an annotation + genome (any species), synthetic, or a saved template file; and the decoys.

Annotation formats: GFF3 (RefSeq / CAT / Liftoff: gene-like records with `Name=`, transcript-like records with
`Parent=` to the gene, exons with `Parent=` to the transcript — or straight to the gene for pseudogenes) and GTF
(`gene_id` / `gene_name` / `transcript_id`), plain or `.gz`. The canonical transcript of a gene = most exons, tie:
longest span. The template is the genomic span from its first to its last exon in TRANSCRIPT orientation.
"""
import gzip
import json
import os
import re

from .model import Exon, GeneModel, rc, synthetic_model

GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene"}
TX_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "lncRNA", "pseudogenic_transcript", "primary_transcript", "tRNA",
            "rRNA", "snRNA", "snoRNA", "miRNA", "misc_RNA", "V_gene_segment", "C_gene_segment", "J_gene_segment",
            "D_gene_segment", "unconfirmed_transcript"}


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _attr(col9, key):
    m = re.search(r'(?:^|;)\s*' + re.escape(key) + r'=([^;]+)', col9)
    return m.group(1) if m else None


def _gtf_attr(col9, key):
    m = re.search(key + r' "([^"]+)"', col9)
    return m.group(1) if m else None


def load_genes(annotation, chrom=None, gene=None):
    """{gene_name: {"chrom", "strand", "transcripts": {tx_id: [(start1, end1), ...]}}} for one chromosome (or, with
    `gene`, for that gene wherever it is; a name scan finds its chromosome first so a 100 MB GFF is parsed once)."""
    is_gtf = ".gtf" in os.path.basename(str(annotation)).lower()
    if gene is not None and chrom is None:
        chrom = _find_gene_chrom(annotation, gene, is_gtf)
        if chrom is None:
            raise ValueError(f"gene {gene!r} not found in {annotation}")
    genes, id2name, tx2gene, tx_exons, tx_strand = {}, {}, {}, {}, {}
    with _open(annotation) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or (chrom is not None and f[0] != chrom):
                continue
            typ = f[2]
            if is_gtf:
                if typ != "exon":
                    continue
                g = _gtf_attr(f[8], "gene_name") or _gtf_attr(f[8], "gene_id")
                t = _gtf_attr(f[8], "transcript_id")
                if not g or not t:
                    continue
                genes.setdefault(g, {"chrom": f[0], "strand": f[6], "transcripts": {}})
                genes[g]["transcripts"].setdefault(t, []).append((int(f[3]), int(f[4])))
                continue
            if typ in GENE_TYPES:
                gid, name = _attr(f[8], "ID"), _attr(f[8], "Name") or _attr(f[8], "gene") or _attr(f[8], "ID")
                if gid and name:
                    id2name[gid] = name
                    genes.setdefault(name, {"chrom": f[0], "strand": f[6], "transcripts": {}})
            elif typ in TX_TYPES:
                tid, par = _attr(f[8], "ID"), _attr(f[8], "Parent")
                if tid and par:
                    par = par.split(",")[0]
                    if par in id2name:
                        tx2gene[tid] = id2name[par]; tx_strand[tid] = f[6]
                    elif par in tx2gene:                      # nested transcript records (rare)
                        tx2gene[tid] = tx2gene[par]; tx_strand[tid] = f[6]
            elif typ == "exon":
                par = _attr(f[8], "Parent")
                if not par:
                    continue
                par = par.split(",")[0]
                if par in tx2gene:
                    tx_exons.setdefault(par, []).append((int(f[3]), int(f[4])))
                elif par in id2name:                          # exon parented straight to the gene (pseudogenes)
                    tx2gene.setdefault(par, id2name[par]); tx_strand.setdefault(par, f[6])
                    tx_exons.setdefault(par, []).append((int(f[3]), int(f[4])))
    for t, ex in tx_exons.items():
        g = tx2gene[t]
        genes[g]["transcripts"][t] = sorted(set(ex))
        genes[g]["strand"] = tx_strand.get(t, genes[g]["strand"])
    for g in list(genes):
        if not genes[g]["transcripts"]:
            del genes[g]
        else:
            for t in genes[g]["transcripts"]:
                genes[g]["transcripts"][t].sort()
    if gene is not None:
        if gene not in genes:
            raise ValueError(f"gene {gene!r} has no exons on {chrom} in {annotation}")
        return {gene: genes[gene]}
    return genes


def _find_gene_chrom(annotation, gene, is_gtf):
    pat = (f'gene_name "{gene}"', f'gene_id "{gene}"') if is_gtf else (f"Name={gene};", f"Name={gene}\n", f"gene={gene};")
    with _open(annotation) as fh:
        for line in fh:
            if any(p in line for p in pat):
                f = line.split("\t", 3)
                if len(f) > 2 and (is_gtf or f[2] in GENE_TYPES or f[2] == "exon"):
                    return f[0]
    return None


def canonical_transcript(gene_rec):
    """(tx_id, exons) with the most exons, tie: the longest span."""
    best = None
    for t, ex in gene_rec["transcripts"].items():
        key = (len(ex), ex[-1][1] - ex[0][0], t)
        if best is None or key > best[0]:
            best = (key, t, ex)
    return best[1], best[2]


def model_from_gene(fa, name, rec, tx=None, intron_cap=None):
    """GeneModel of gene `name` (pysam FastaFile `fa`), in transcript orientation, with provenance."""
    tid, exons = (tx, rec["transcripts"][tx]) if tx else canonical_transcript(rec)
    chrom, strand = rec["chrom"], rec["strand"]
    gs, ge = exons[0][0], exons[-1][1]
    seq = fa.fetch(chrom, gs - 1, ge).upper()
    if strand == "-":
        seq = rc(seq)
        local = sorted((ge - e, ge - s + 1) for s, e in exons)
    else:
        local = [(s - gs, e - gs + 1) for s, e in exons]
    # overlapping/abutting annotated exons would break the model: merge them
    merged = []
    for s, e in local:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    m = GeneModel(seq, [Exon(s, e, f"e{i + 1}") for i, (s, e) in enumerate(merged)],
                  provenance={"source": "annotation", "gene": name, "transcript": tid, "chrom": chrom, "start": gs, "end": ge,
                              "strand": strand, "n_exons_annotated": len(exons), "motifs": None})
    m.provenance["motifs"] = ["".join(x) for x in m.junction_motifs()]
    if intron_cap:
        from .ops import op_intron_resize
        import random
        rng = random.Random(0)
        for k, (s, e) in enumerate(m.introns()):
            if e - s > intron_cap:
                op_intron_resize(m, {"intron": k + 1, "length": int(intron_cap)}, rng)
        m.provenance["intron_cap"] = int(intron_cap)
    return m


def pick_genes(genes, rng, n, min_exons=3, max_span=60000, exclude=(), require_canonical=True, fa=None, max_exons=None):
    """`n` gene names drawn at random among those with >= min_exons exons on their canonical transcript, span <=
    max_span, not in `exclude`, and (with a FASTA) every intron GT…AG."""
    cands = []
    for g, rec in genes.items():
        if g in exclude or g.startswith("LOC") and rec["strand"] not in "+-":
            continue
        t, ex = canonical_transcript(rec)
        if len(ex) < min_exons or (max_exons and len(ex) > max_exons) or ex[-1][1] - ex[0][0] > max_span:
            continue
        cands.append(g)
    cands.sort()
    out = []
    order = list(cands)
    rng.shuffle(order)
    for g in order:
        if len(out) >= n:
            break
        if require_canonical and fa is not None:
            m = model_from_gene(fa, g, genes[g])
            if not m.canonical() or "N" in m.seq:
                continue
        out.append(g)
    if len(out) < n:
        raise ValueError(f"only {len(out)} eligible genes (need {n}); relax min_exons/max_span or pick another chromosome")
    return out


def resolve_template(spec, rng):
    """GeneModel from a template spec: source annotation | synthetic | file."""
    src = spec.get("source", "annotation")
    if src == "file":
        return GeneModel.load(spec["path"])
    if src == "synthetic":
        m = synthetic_model(spec.get("exons", [200, 150, 300, 120, 500]), spec.get("introns", [1500, 2500, 900, 1800]), rng,
                            spec.get("gc", 0.41))
        m.provenance["motifs"] = ["".join(x) for x in m.junction_motifs()]
        return m
    if src != "annotation":
        raise ValueError(f"unknown template source {src!r}")
    import pysam
    fa = pysam.FastaFile(spec["genome"])
    ann = spec["annotation"]
    if spec.get("gene"):
        genes = load_genes(ann, spec.get("chrom"), spec["gene"])
        name = spec["gene"]
    elif spec.get("transcript"):
        genes = load_genes(ann, spec.get("chrom"))
        hits = [g for g, r in genes.items() if spec["transcript"] in r["transcripts"]]
        if not hits:
            raise ValueError(f"transcript {spec['transcript']} not found")
        name = hits[0]
    else:
        if not spec.get("chrom"):
            raise ValueError("a random template needs 'chrom'")
        genes = load_genes(ann, spec["chrom"])
        name = pick_genes(genes, rng, 1, spec.get("min_exons", 4), spec.get("max_span", 60000), fa=fa,
                          max_exons=spec.get("max_exons"))[0]
    m = model_from_gene(fa, name, genes[name], spec.get("transcript"), spec.get("intron_cap"))
    m.provenance["genome"] = spec["genome"]; m.provenance["annotation"] = ann
    return m


def resolve_decoys(spec, template, rng):
    """[GeneModel] of `n` unrelated genes (same annotation/genome as the template unless overridden; `chrom` defaults to
    the template's), or the models of a `file` list."""
    if not spec or not spec.get("n", 0) and not spec.get("path"):
        return []
    if spec.get("source") == "file" or spec.get("path"):
        with open(spec["path"]) as fh:
            return [GeneModel.from_json(d) for d in json.load(fh)]
    prov = template.provenance
    ann = spec.get("annotation") or prov.get("annotation")
    genome = spec.get("genome") or prov.get("genome")
    chrom = spec.get("chrom") or prov.get("chrom")
    if not (ann and genome and chrom):
        raise ValueError("decoys need annotation, genome and chrom (a synthetic template has none: give them explicitly)")
    import pysam
    fa = pysam.FastaFile(genome)
    genes = load_genes(ann, chrom)
    names = pick_genes(genes, rng, int(spec["n"]), spec.get("min_exons", 3), spec.get("max_span", 30000),
                       exclude={prov.get("gene")}, fa=fa, max_exons=spec.get("max_exons"))
    out = []
    for g in names:
        m = model_from_gene(fa, g, genes[g], intron_cap=spec.get("intron_cap"))
        m.provenance["genome"] = genome; m.provenance["annotation"] = ann
        out.append(m)
    return out
