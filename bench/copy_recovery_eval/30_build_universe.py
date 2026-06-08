#!/usr/bin/env python3
"""U3 CLI: build the paralog universe from the RefSeq annotation + genome.
Steps: (1) parse annotation -> gene_tx + per-transcript exons; (2) write a
representative-sequence FASTA per gene (longest transcript) via gffread; (3)
self-align with minimap2; (4) keep pairs >= identity/coverage; (5) build_universe.
Outputs universe.tsv and exons.tsv (transcript_id,chrom,exon_starts,exon_ends).

Real-GFF inspection findings (GGO_tx.gff, GFF3 format):
  - Transcript features: mRNA / lnc_RNA / transcript (+ snRNA/snoRNA/tRNA etc.)
  - Transcript ID: transcript_id=<accession>  (no rna- prefix; this is canonical)
  - Gene ID:       gene=<gene_name>            (also in Parent=gene-<name>)
  - Exon parent:   transcript_id=<accession>  (same clean accession on every exon line)
  ID=rna-<accession> also present on transcript lines but has unwanted rna- prefix.
  The parser uses transcript_id= exclusively for both transcript and exon features."""
import argparse, csv, subprocess, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import lib_eval


# Feature types treated as transcripts in the NCBI GFF3 annotation.
_TRANSCRIPT_FEATS = {
    "mRNA", "lnc_RNA", "transcript", "snRNA", "snoRNA", "tRNA",
    "miRNA", "rRNA", "ncRNA", "primary_transcript",
}


def parse_annotation(gff):
    """Return (gene_tx: {gene:[tx]}, exons: {tx:(chrom,[(s,e)])}, tx_gene:{tx:gene}).

    Adapted for NCBI GFF3 (GGO_tx.gff):
      - transcript_id= is present on BOTH transcript-level and exon-level features.
      - gene= holds the gene symbol / locus tag on transcript-level features.
      - We fall back to ID= (stripping a leading 'rna-' prefix) and then to
        Parent= (stripping a leading 'gene-' prefix) if the preferred keys are absent.
    """
    gene_tx, exons, tx_gene = {}, {}, {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            chrom, feat = c[0], c[2]
            try:
                start, end = int(c[3]) - 1, int(c[4])
            except ValueError:
                continue
            attrs = c[8]

            # Parse key=value; (GFF3 style).  The GFF may also contain GTF-style
            # fields inside some attributes — but GGO_tx.gff is pure GFF3.
            kv = {}
            for field in attrs.split(";"):
                field = field.strip()
                if not field:
                    continue
                if "=" in field:
                    k, v = field.split("=", 1)
                    kv[k.strip()] = v.strip()

            if feat in _TRANSCRIPT_FEATS:
                # transcript_id= is the canonical clean accession (no rna- prefix).
                tx = kv.get("transcript_id")
                if not tx:
                    raw_id = kv.get("ID", "")
                    tx = raw_id[4:] if raw_id.startswith("rna-") else raw_id or None

                # gene= holds the gene name on transcript-level features.
                gene = kv.get("gene")
                if not gene:
                    raw_par = kv.get("Parent", "")
                    gene = raw_par[5:] if raw_par.startswith("gene-") else raw_par or None

                if tx and gene:
                    gene_tx.setdefault(gene, []).append(tx)
                    tx_gene[tx] = gene
                    exons.setdefault(tx, (chrom, []))

            elif feat == "exon":
                # transcript_id= on exon lines is the same clean accession.
                # Fallback: Parent=rna-<accession> → strip the rna- prefix.
                # Do NOT fall back to bare gene-level Parents (pseudo=true exons
                # have Parent=gene-XXX with no transcript_id= — skip them).
                tx = kv.get("transcript_id")
                if not tx:
                    raw_par = kv.get("Parent", "")
                    tx = raw_par[4:] if raw_par.startswith("rna-") else None

                if tx:
                    ch, lst = exons.setdefault(tx, (chrom, []))
                    lst.append((start, end))

    return gene_tx, exons, tx_gene


def write_exons_tsv(exons, path):
    with open(path, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["transcript_id", "chrom", "exon_starts", "exon_ends"])
        for tx in sorted(exons):
            chrom, lst = exons[tx]
            lst = sorted(lst)
            w.writerow([tx, chrom,
                        ",".join(str(s) for s, _ in lst),
                        ",".join(str(e) for _, e in lst)])


def gene_paralog_pairs(gene_fasta, min_identity, min_cov_frac, minimap2):
    """All-vs-all minimap2 of per-gene representative seqs; return passing gene pairs."""
    paf = subprocess.run(
        [minimap2, "-x", "asm20", "-X", "--no-long-join", "-c", gene_fasta, gene_fasta],
        capture_output=True, text=True, check=True).stdout
    pairs = set()
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, qlen, tname = f[0], int(f[1]), f[5]
        nmatch, alen = int(f[9]), int(f[10])
        if q == tname:
            continue
        ident = nmatch / alen if alen else 0.0
        cov = alen / qlen if qlen else 0.0
        if ident >= min_identity and cov >= min_cov_frac:
            pairs.add(tuple(sorted((q, tname))))
    return sorted(pairs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--min-identity", type=float, default=0.90)
    ap.add_argument("--min-cov-frac", type=float, default=0.50)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--gffread", default="gffread")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sys.stderr.write("[U3] Parsing annotation...\n")
    gene_tx, exons, tx_gene = parse_annotation(args.ref_gff)
    sys.stderr.write(f"[U3] {len(tx_gene)} transcripts in {len(gene_tx)} genes\n")

    write_exons_tsv(exons, os.path.join(args.outdir, "exons.tsv"))
    sys.stderr.write(f"[U3] exons.tsv written ({len(exons)} transcripts)\n")

    # Pick the longest representative transcript per gene (by total exon length).
    tx_len = {tx: sum(e - s for s, e in lst) for tx, (_, lst) in exons.items()}
    rep_tx = {}
    for gene, txs in gene_tx.items():
        best = max((t for t in txs if t in tx_len), key=lambda t: tx_len[t], default=None)
        if best:
            rep_tx[best] = gene

    sys.stderr.write(f"[U3] Extracting {len(rep_tx)} representative sequences via gffread...\n")
    tx_fa = os.path.join(args.outdir, "tx.fa")
    subprocess.run([args.gffread, "-w", tx_fa, "-g", args.genome, args.ref_gff], check=True)

    gene_fa = os.path.join(args.outdir, "gene_rep.fa")
    kept = 0
    with open(tx_fa) as f, open(gene_fa, "w") as o:
        cur_keep = False
        for line in f:
            if line.startswith(">"):
                # gffread emits the GFF3 ID= value as the FASTA header, which for
                # NCBI GFF3 is "rna-<accession>" (with the rna- prefix).
                # Our rep_tx dict is keyed by the clean accession (transcript_id=).
                # Accept both forms: strip a leading "rna-" if present.
                raw = line[1:].split()[0]
                name = raw[4:] if raw.startswith("rna-") else raw
                cur_keep = name in rep_tx
                if cur_keep:
                    o.write(f">{rep_tx[name]}\n")
                    kept += 1
            elif cur_keep:
                o.write(line)
    sys.stderr.write(f"[U3] gene_rep.fa written ({kept} sequences)\n")

    sys.stderr.write("[U3] Running minimap2 all-vs-all self-alignment...\n")
    pairs = gene_paralog_pairs(gene_fa, args.min_identity, args.min_cov_frac, args.minimap2)
    sys.stderr.write(f"[U3] {len(pairs)} paralog gene pairs passing thresholds\n")

    uni = lib_eval.build_universe(gene_tx, pairs)
    with open(os.path.join(args.outdir, "universe.tsv"), "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["transcript_id", "gene_id", "family_id", "n_family_copies"])
        for r in uni:
            w.writerow([r["transcript_id"], r["gene_id"], r["family_id"], r["n_family_copies"]])
    n_fam = len({r["family_id"] for r in uni})
    sys.stderr.write(f"[U3] {len(uni)} transcripts in {n_fam} paralog families\n")
    sys.stderr.write(f"[U3] universe.tsv written to {args.outdir}\n")


if __name__ == "__main__":
    main()
