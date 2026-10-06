"""The artificial genome: background contigs + planted copies -> FASTA, truth GTF/GFF3, copies table, manifest.

Planting REPLACES background bases (contig lengths stay what the background gives), so every coordinate is known
before the FASTA is written. A copy on strand `-` is planted as the reverse complement of its model; its exons map to
genomic coordinates accordingly. Copies with `in_reference: false` go on their own contig(s) that `genome.ref.fa`
omits (the reads come from them; the aligner never sees them; `genome.truth.fa` has everything).

Products (DIR/): genome.truth.fa, genome.ref.fa (+ .fai), truth.gtf, truth.gff3, truth.families.tsv, copies.tsv,
copies.fa, manifest.json, template.json, decoys.json.
"""
import json
import os
import random
import subprocess

from .model import GeneModel, rc, random_seq, seeded
from .ops import apply_ops
from .template import resolve_decoys, resolve_template

DEFAULT_SPACING = 20000
DEFAULT_START = 20000


class Planted:
    """One planted sequence: a copy of A, or a decoy."""

    def __init__(self, cid, kind, model, contig, pos, strand, in_reference=True, expression=None, isoforms=None, family=None,
                 explicit_contig=False):
        self.id, self.kind, self.model, self.contig, self.pos, self.strand = cid, kind, model, contig, int(pos), strand
        self.in_reference, self.expression, self.isoforms = bool(in_reference), expression, isoforms or []
        self.family = family or ("A" if kind == "copy" else cid)
        self.explicit_contig = bool(explicit_contig)

    @property
    def length(self):
        return len(self.model.seq)

    @property
    def end(self):
        return self.pos + self.length

    def planted_seq(self):
        return rc(self.model.seq) if self.strand == "-" else self.model.seq

    def genomic(self, s, e):
        """model-local [s,e) -> contig 0-based half-open."""
        if self.strand == "-":
            return self.pos + self.length - e, self.pos + self.length - s
        return self.pos + s, self.pos + e

    def exon_intervals(self, skip=()):
        """[(gstart, gend, label)] of an isoform's exons, in genomic order."""
        iv = [(*self.genomic(e.start, e.end), e.label) for e in self.model.rna_exons(skip)]
        return sorted(iv)

    def dna_exon_intervals(self):
        return sorted((*self.genomic(e.start, e.end), e.label, e.in_rna, e.inverted) for e in self.model.exons)

    def junctions(self, skip=()):
        """set of genomic 0-based half-open introns of an isoform."""
        return {self.genomic(s, e) for s, e in self.model.chain_junctions(skip)}

    def isoform_list(self):
        """[(name, skip tuple, weight)]: iso0 = the full chain, then the spec's."""
        out = [("iso0", (), 1.0)]
        for k, iso in enumerate(self.isoforms):
            out.append((f"iso{k + 1}", tuple(iso.get("skip", [])), float(iso.get("weight", 0.5))))
        return out


def build(spec, out_dir, log=print):
    """Resolve the spec, build everything, write the products. Returns (planted list, manifest dict)."""
    os.makedirs(out_dir, exist_ok=True)
    seed = int(spec.get("seed", 1))
    rng_t = seeded(seed, "template")
    template = resolve_template(spec["template"], rng_t)
    template.save(os.path.join(out_dir, "template.json"))
    log(f"template: {template.summary()}  [{template.provenance.get('gene') or template.provenance.get('source')}]")
    decoys = resolve_decoys(spec.get("decoys"), template, seeded(seed, "decoys"))
    with open(os.path.join(out_dir, "decoys.json"), "w") as fh:
        json.dump([d.to_json() for d in decoys], fh)

    # ---- copies
    planted, models = [], {}
    for k, c in enumerate(spec["copies"]):
        cid = c.get("id", f"A{k + 1}" if k else "A")
        m = template.copy()
        m.ops = []
        strand = c.get("strand", "+")
        ops = list(c.get("ops", []))
        if any(o.get("op") == "invert" and o.get("whole") for o in ops):
            strand = "-" if strand == "+" else "+"
        apply_ops(m, ops, seeded(seed, "copy", cid), donors=models)
        models[cid] = m
        planted.append(Planted(cid, "copy", m, c.get("contig", "sim"), c.get("pos", -1), strand, c.get("in_reference", True),
                               c.get("expression"), c.get("isoforms"), c.get("family"), "contig" in c))
    for k, d in enumerate(decoys):
        dspec = (spec.get("decoys") or {})
        planted.append(Planted(f"D{k + 1}", "decoy", d, dspec.get("contig", "sim"), -1, "+" if k % 2 == 0 else "-", True,
                               dspec.get("expression"), None))

    # ---- layout
    layout = spec.get("layout", {})
    spacing, start = int(layout.get("spacing", DEFAULT_SPACING)), int(layout.get("start", DEFAULT_START))
    last_end = {}
    for p in planted:
        if not p.in_reference and not p.explicit_contig:
            p.contig = f"absent_{p.id}"          # an absent copy gets its own contig unless the spec placed it
        if p.pos < 0:
            p.pos = last_end.get(p.contig, start - spacing) + spacing
        last_end[p.contig] = max(last_end.get(p.contig, 0), p.end)
    by_contig = {}
    for p in planted:
        by_contig.setdefault(p.contig, []).append(p)
    for ctg, ps in by_contig.items():
        ps.sort(key=lambda p: p.pos)
        for a, b in zip(ps, ps[1:]):
            if b.pos < a.end:
                raise ValueError(f"{a.id} [{a.pos},{a.end}) overlaps {b.id} [{b.pos},{b.end}) on {ctg}")
        absent = {p.in_reference for p in ps}
        if len(absent) > 1:
            raise ValueError(f"contig {ctg} mixes reference-present and reference-absent copies; give the absent ones their own contig")

    # ---- background per contig
    bg = spec.get("background", {"source": "random", "length": 300000})
    contigs = {}
    realised_bg = {}
    for i, (ctg, ps) in enumerate(sorted(by_contig.items())):
        need = max(p.end for p in ps) + spacing
        seq, desc = _background(bg, i, need, seeded(seed, "background", ctg), planted)
        contigs[ctg] = seq
        realised_bg[ctg] = desc
    for ctg, ps in by_contig.items():
        seq = contigs[ctg]
        for p in ps:
            s = p.planted_seq()
            seq = seq[:p.pos] + s + seq[p.end:]
        contigs[ctg] = seq

    # ---- products
    _write_fasta(os.path.join(out_dir, "genome.truth.fa"), contigs)
    ref = {c: s for c, s in contigs.items() if all(p.in_reference for p in by_contig[c])}
    _write_fasta(os.path.join(out_dir, "genome.ref.fa"), ref)
    for fn in ("genome.truth.fa", "genome.ref.fa"):
        subprocess.run(["samtools", "faidx", os.path.join(out_dir, fn)], check=True)
    write_truth(planted, out_dir)
    manifest = {
        "spec": spec, "seed": seed, "template": template.provenance, "template_motifs": template.provenance.get("motifs"),
        "background": realised_bg, "contigs": {c: len(s) for c, s in contigs.items()},
        "copies": [{"id": p.id, "kind": p.kind, "family": p.family, "contig": p.contig, "pos": p.pos, "end": p.end, "strand": p.strand,
                    "in_reference": p.in_reference, "expression": p.expression, "isoforms": p.isoforms, "n_exons_dna": len(p.model.exons),
                    "n_exons_rna": len(p.model.rna_exons()), "ops": p.model.ops,
                    "exons": [{"label": e.label, "start": e.start, "end": e.end, "in_rna": e.in_rna, "inverted": e.inverted}
                              for e in p.model.exons],
                    "provenance": p.model.provenance} for p in planted],
    }
    with open(os.path.join(out_dir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    log(f"genome: {', '.join(f'{c} {len(s):,} bp' for c, s in contigs.items())}; reference omits "
        f"{[c for c in contigs if c not in ref] or 'nothing'}")
    for p in planted:
        log(f"  {p.id:5s} {p.kind:5s} {p.contig}:{p.pos + 1}-{p.end} {p.strand} exons {len(p.model.exons)} "
            f"(rna {len(p.model.rna_exons())}) ops {[o['op'] for o in p.model.ops]}")
    return planted, manifest


def _background(bg, idx, need, rng, planted):
    src = bg.get("source", "random")
    if src == "random":
        L = max(int(bg.get("length", 300000)), need)
        return random_seq(L, rng, bg.get("gc", 0.41)), {"source": "random", "length": L}
    if src != "fasta":
        raise ValueError(f"unknown background source {src!r}")
    import pysam
    fa = pysam.FastaFile(bg["path"])
    chrom, rng_ = bg["region"].split(":")
    s, e = (int(x.replace(",", "")) for x in rng_.split("-"))
    L0 = e - s + 1
    # contig idx takes the idx-th consecutive slice of the same length after the region (idx 0 = the region itself)
    s = s + idx * L0
    e = max(s + L0 - 1, s + need - 1)
    if e > fa.get_reference_length(chrom):
        raise ValueError(f"background {chrom}:{s}-{e} runs past the end of {chrom}")
    # the real template / decoy loci must not sit inside the slice (they would be unplanned extra copies)
    for p in planted:
        pv = p.model.provenance
        if pv.get("chrom") == chrom and pv.get("start") and not (pv["end"] < s or pv["start"] > e):
            raise ValueError(f"background {chrom}:{s}-{e} contains the real locus of {p.id} ({pv['chrom']}:{pv['start']}-{pv['end']}); "
                             f"pick another region")
    seq = fa.fetch(chrom, s - 1, e).upper()
    return seq, {"source": "fasta", "path": bg["path"], "region": f"{chrom}:{s}-{e}"}


def _write_fasta(path, contigs):
    with open(path, "w") as fh:
        for c, s in contigs.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def write_truth(planted, out_dir):
    """truth.gtf (one transcript per isoform), truth.gff3 (gene + mRNA + exon, `Name=` on genes and `gene=` on exons:
    the attributes mcl_families --gff reads), truth.families.tsv (family_score's `Gene Name / Family ID / Contig`),
    copies.tsv, copies.fa (each planted sequence in transcript orientation)."""
    gtf = open(os.path.join(out_dir, "truth.gtf"), "w")
    gff = open(os.path.join(out_dir, "truth.gff3"), "w")
    fam = open(os.path.join(out_dir, "truth.families.tsv"), "w")
    cop = open(os.path.join(out_dir, "copies.tsv"), "w")
    cfa = open(os.path.join(out_dir, "copies.fa"), "w")
    gff.write("##gff-version 3\n")
    fam.write("Gene Name\tFamily ID\tContig\n")
    cop.write("copy_id\tkind\tfamily\tcontig\tstart\tend\tstrand\tin_reference\texpression\tn_exons_dna\tn_exons_rna\tlength\tops\n")
    for p in planted:
        gs, ge = p.pos + 1, p.end
        gff.write(f"{p.contig}\tfamsim\tgene\t{gs}\t{ge}\t.\t{p.strand}\t.\tID=gene-{p.id};Name={p.id};gene={p.id};"
                  f"gene_biotype={'protein_coding' if p.kind == 'copy' else 'decoy'};family={p.family}\n")
        for name, skip, w in p.isoform_list():
            iv = p.exon_intervals(skip)
            if not iv:
                continue
            tid = f"{p.id}.{name}"
            gtf.write(f"{p.contig}\tfamsim\ttranscript\t{iv[0][0] + 1}\t{iv[-1][1]}\t.\t{p.strand}\t.\t"
                      f'gene_id "{p.id}"; transcript_id "{tid}"; family "{p.family}"; kind "{p.kind}";\n')
            gff.write(f"{p.contig}\tfamsim\tmRNA\t{iv[0][0] + 1}\t{iv[-1][1]}\t.\t{p.strand}\t.\tID=rna-{tid};Parent=gene-{p.id};Name={tid};gene={p.id}\n")
            order = iv if p.strand == "+" else iv[::-1]
            for k, (s, e, label) in enumerate(order):
                gtf.write(f"{p.contig}\tfamsim\texon\t{s + 1}\t{e}\t.\t{p.strand}\t.\t"
                          f'gene_id "{p.id}"; transcript_id "{tid}"; exon_number "{k + 1}"; exon_label "{label}";\n')
                gff.write(f"{p.contig}\tfamsim\texon\t{s + 1}\t{e}\t.\t{p.strand}\t.\tID=exon-{tid}-{k + 1};Parent=rna-{tid};gene={p.id};exon_label={label}\n")
        fam.write(f"{p.id}\t{p.family}\t{p.contig}\n")
        cop.write(f"{p.id}\t{p.kind}\t{p.family}\t{p.contig}\t{gs}\t{ge}\t{p.strand}\t{int(p.in_reference)}\t"
                  f"{'' if p.expression is None else p.expression}\t{len(p.model.exons)}\t{len(p.model.rna_exons())}\t{p.length}\t"
                  f"{json.dumps([o['op'] for o in p.model.ops])}\n")
        cfa.write(f">{p.id} {p.kind} {p.contig}:{gs}-{ge} {p.strand} family={p.family}\n{p.model.seq}\n")
    for fh in (gtf, gff, fam, cop, cfa):
        fh.close()


def load_planted(out_dir):
    """Rebuild the Planted list from DIR/manifest.json + copies' models (template.json / decoys.json are not enough:
    the models carry the ops, so they are stored in the manifest's copies as exons + the sequence in copies.fa)."""
    with open(os.path.join(out_dir, "manifest.json")) as fh:
        man = json.load(fh)
    seqs = {}
    name = None
    for line in open(os.path.join(out_dir, "copies.fa")):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    out = []
    for c in man["copies"]:
        m = GeneModel(seqs[c["id"]], [(e["start"], e["end"], e["label"], e["in_rna"], e["inverted"]) for e in c["exons"]],
                      c["ops"], c.get("provenance"))
        out.append(Planted(c["id"], c["kind"], m, c["contig"], c["pos"], c["strand"], c["in_reference"], c["expression"],
                           c.get("isoforms"), c.get("family")))
    return out, man
