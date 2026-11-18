"""PREREG r1797 / r1798 (docs/PREREG_famwise_sim_2026-11-18.md): the read model of a real library, measured on genes
outside the families.

  read_profile.py BAM GENES.tsv GFF3.gz GENOME.fa FAM.json OUT.json [--cap 50] [--seed 2026] [--chroms chr1,...,chr22]

BAM        coordinate-sorted, indexed, aligned with --eqx.
GENES.tsv  gene_id gene_name gene_biotype source chrom start0 end strand n_transcripts exon_blocks (with a header).
GFF3.gz    bgzip + tabix: gene / transcript (ID, Parent = gene id) / exon (Parent = transcript id).
FAM.json   {"sets": {NAME: [{"family": .., "copies": [{"gene_id": ..}]}]}}: every gene in it is a member and left out.
Genes: protein_coding, one transcript with >= 2 exons, on --chroms, not a member. Their reads: primary, MAPQ >= 60, on the
gene's strand, within 1000 bases of its span. Per gene a seeded sample of at most --cap reads; a substitution is the
individual's allele, not an error, when >= 3 reads and >= 20% of the coverage carry it among ALL the gene's reads.
OUT.json   classes (bounds of the transcript length classes), end_pairs (per class: the fractions of the transcript missing
           at 5' and 3', one pair per read), rates (per read: substitution, insertion and deletion events per aligned base),
           ins_len / del_len / ins_hp_frac / del_hp_frac and kept_intron_per_intron as bench/o3_graph/realworld.py reads
           them; pooled sub_error / ins_event / del_event, read_length, full5_frac and class_reads for the record."""
import argparse, bisect, collections, csv, importlib.util, json, os, random
import pysam

CLASSES = [1000, 2000, 3000, 5000]
MIN_MAPQ, FLANK, FULL = 60, 1000, 50
_spec = importlib.util.spec_from_file_location("o3_graph_calibrate", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "o3_graph", "calibrate.py"))
C = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(C)

def klass(tx_len):
    """number of bounds <= the spliced transcript length: 0 (< 1 kb) .. 4 (>= 5 kb)"""
    return bisect.bisect_right(CLASSES, tx_len)

def load_members(fam_json):
    with open(fam_json) as fh:
        return {c["gene_id"] for fams in json.load(fh)["sets"].values() for f in fams for c in f["copies"]}

def transcripts(gff, g):
    """exons (0-based half-open, sorted) of every transcript of the gene"""
    tids, ex = set(), collections.defaultdict(list)
    for row in gff.fetch(g["chrom"], int(g["start0"]), int(g["end"])):
        f = row.split("\t")
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "transcript" and at.get("Parent") == g["gene_id"]:
            tids.add(at["ID"])
        elif f[2] == "exon":
            ex[at.get("Parent")].append((int(f[3]) - 1, int(f[4])))
    return [sorted(ex[t]) for t in sorted(tids) if ex.get(t)]

def gene_set(genes_tsv, gff, members, chroms):
    """-> [(gene row, exons)]: protein-coding, one transcript of >= 2 exons, on `chroms`, not a member"""
    out = []
    with open(genes_tsv) as fh:
        for g in csv.DictReader(fh, delimiter="\t"):
            if g["gene_biotype"] != "protein_coding" or g["n_transcripts"] != "1" or g["chrom"] not in chroms or g["gene_id"] in members:
                continue
            tx = transcripts(gff, g)
            if len(tx) == 1 and len(tx[0]) >= 2:
                out.append((g, tx[0]))
    return out

def gene_reads(bam, fa, g, exons, cap, seed):
    """the gene's sampled reads [{name, ch, qlen, ev, sub, ends, tx_len, kept, introns}]; sub = substitutions that are not
    alleles, judged on the pileup of every qualifying read (a sample of 50 cannot tell an allele from an error)"""
    ch, gs, ge, cache = g["chrom"], int(g["start0"]), int(g["end"]), {}

    def ref_base(p):
        if p not in cache:
            cache[p] = fa.fetch(ch, max(0, p), p + 1).upper() if p >= 0 else "N"
        return cache[p]
    reads = []
    for r in bam.fetch(ch, gs, ge):
        if (r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < MIN_MAPQ or r.query_sequence is None
                or ("-" if r.is_reverse else "+") != g["strand"] or r.reference_start < gs - FLANK or r.reference_end > ge + FLANK):
            continue
        reads.append(dict(name=r.query_name, ch=ch, qlen=r.query_length,
                          ev=C.cigar_events(r.cigartuples, r.reference_start, r.query_sequence, ref_base)))
    counts = collections.Counter((ch, p, b) for rd in reads for p, b in rd["ev"]["mismatch"])
    pos, cov = sorted({p for (_, p, _), n in counts.items() if n >= C.MIN_VAR_READS}), collections.Counter()
    for rd in reads:
        for bs, be in rd["ev"]["blocks"]:
            for p in pos[bisect.bisect_left(pos, bs):bisect.bisect_left(pos, be)]:
                cov[(ch, p)] += 1
    kind = C.classify_mismatches(counts, cov)
    out, tx_len = [], sum(e - s for s, e in exons)
    for rd in reads if len(reads) <= cap else random.Random(f"{seed}:{g['gene_id']}").sample(reads, cap):
        blocks = rd["ev"]["blocks"]
        ends = C.transcript_ends(blocks, exons, g["strand"])
        if ends is not None:
            out.append(dict(rd, sub=sum(kind[(ch, p, b)] == "error" for p, b in rd["ev"]["mismatch"]), ends=ends, tx_len=tx_len,
                            kept=C.retained_introns(blocks, exons), introns=len(exons) - 1))
    return out

def profile(rows, source, genes):
    es = C.error_stats(rows)              # its sub_error would judge alleles on the sample alone: replaced below
    pairs = {str(k): [] for k in range(len(CLASSES) + 1)}
    for rd in rows:
        pairs[str(klass(rd["tx_len"]))].append([round(x / rd["tx_len"], 4) for x in rd["ends"]])
    return dict(source=source, genes=genes, reads=len(rows), classes=CLASSES, end_pairs=pairs,
                rates=[[round(n / max(1, rd["ev"]["aligned"]), 6) for n in (rd["sub"], len(rd["ev"]["ins"]), len(rd["ev"]["del"]))] for rd in rows],
                **{k: es[k] for k in ("ins_len", "del_len", "ins_hp_frac", "del_hp_frac", "ins_event", "del_event")},
                kept_intron_per_intron=sum(rd["kept"] for rd in rows) / max(1, sum(rd["introns"] for rd in rows)),
                sub_error=sum(rd["sub"] for rd in rows) / max(1, es["aligned_bases"]), read_length=C._quantiles([rd["qlen"] for rd in rows]),
                full5_frac=sum(rd["ends"][0] <= FULL for rd in rows) / max(1, len(rows)), class_reads={k: len(v) for k, v in pairs.items()})

def summary(rows, prof):
    lines = [f"genes {prof['genes']}  reads {prof['reads']}  events per aligned base: sub {prof['sub_error']:.6f}  ins {prof['ins_event']:.6f}"
             f"  del {prof['del_event']:.6f}  kept introns per intron {prof['kept_intron_per_intron']:.5f}"]
    names = ["< 1 kb", "1-2 kb", "2-3 kb", "3-5 kb", ">= 5 kb"]
    for k, name in enumerate(names):
        v = [rd for rd in rows if klass(rd["tx_len"]) == k]
        both = f"{sum(max(rd['ends']) <= FULL for rd in v) / len(v):.3f}" if v else "-"
        lines.append(f"class {k} ({name}): reads {len(v)}  both ends within {FULL} bases {both}")
    tot = sorted(sum(r) for r in prof["rates"])
    if tot:
        p50, p90, p99 = (tot[min(len(tot) - 1, int(q * len(tot)))] for q in (0.5, 0.9, 0.99))
        lines.append(f"events per aligned base, per read: p50 {p50:.6f}  p90 {p90:.6f}  p99 {p99:.6f}  reads without an event {tot.count(0) / len(tot):.3f}")
    return "\n".join(lines)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for x in ("bam", "genes", "gff", "genome", "fam", "out"):
        ap.add_argument(x)
    ap.add_argument("--cap", type=int, default=50); ap.add_argument("--seed", default="2026")
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 23)))
    a = ap.parse_args()
    rows, genes = [], 0
    with pysam.TabixFile(a.gff) as gff, pysam.AlignmentFile(a.bam) as bam, pysam.FastaFile(a.genome) as fa:
        for g, exons in gene_set(a.genes, gff, load_members(a.fam), a.chroms.split(",")):
            got = gene_reads(bam, fa, g, exons, a.cap, a.seed)
            rows += got
            genes += bool(got)
    prof = profile(rows, a.bam, genes)
    with open(a.out, "w") as fh:
        json.dump(prof, fh)
    print(summary(rows, prof))

if __name__ == "__main__":
    main()
