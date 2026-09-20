#!/usr/bin/env bash
# QuicK-mer2 paralog-specific copy number (parCN) for the gorilla T2T genome — the GOLD-STANDARD copy-vs-allele
# resolver (Soto-2025 / Eichler lab). Settles the O4 candidates' COPY-vs-ALLELE that the read-free structural
# test (bench/copy_vs_allele_structural.py) left AMBIGUOUS (the MHC-class: region duplicated AND het-range).
#
# STATUS (verified 2026-06-30 on the dev box): QuicK-mer2 installs + builds + runs (20 Mb slice, 74 s). But the
# full pipeline needs a CLUSTER, for two independent reasons:
#   (1) `search` (genome reference) needs ~40 GB RAM (-s 1G used 10.7 GB on a 20 Mb slice; the 3.5 Gb genome
#        needs -s ~4G). The dev box has 19 GB. ONE-TIME, then the .qm is reusable forever.
#   (2) `count`/`est` need gorilla DNA WGS reads (we have only RNA IsoSeq on disk). The mGorGor1/"Kamilah" HiFi
#        DNA is public (T2T-Primates project, Makova/Eichler 2024) — set $DNA below to that run.
#
# Run on a >=48 GB node:  bash quickmer2_parcn.sh
set -euo pipefail
QM=${QM:-/home/juanfra/QuicK-mer2/quicKmer2}          # built binary
GENOME=${GENOME:-/home/juanfra/winloci_scratch/GGO.fasta}   # the T2T reference (GCF_029281585.2)
DNA=${DNA:?set DNA to the gorilla WGS reads — a BAM/CRAM/fastq of Kamilah/mGorGor1 HiFi-DNA (SRA, T2T-Primates)}
OUT=${OUT:-/home/juanfra/winloci_scratch/qm2}
THREADS=${THREADS:-16}
SAM=/home/juanfra/miniforge3/bin/samtools
mkdir -p "$OUT"; cd "$OUT"

# ---- 1. SEARCH: build the single-copy k-mer reference (read-free, ~40 GB RAM, ONE-TIME, reusable) ----
# -e 2 masks k-mers within edit distance 2 (the standard, most specific); -s 4G sizes the hash for 3.5 Gb;
# -w 1000 = 1000-kmer windows. Produces GGO.fasta.qm (index), .bed (windows), .qgc (GC).
if [ ! -s "$GENOME.qm" ]; then
  echo "[1/4] search (genome reference; ~40 GB RAM, hours) ..."
  "$QM" search -k 30 -e 2 -t "$THREADS" -s 4G "$GENOME"
fi

# ---- 2. COUNT: k-mer depth from the DNA reads (needs DNA WGS) ----
# Stream reads as fasta into quicKmer2 count; this is the read-depth over the single-copy k-mers.
echo "[2/4] count (DNA k-mer depth) ..."
case "$DNA" in
  *.bam|*.cram) "$SAM" fasta -@ "$THREADS" "$DNA" | "$QM" count "$GENOME" /dev/stdin sample.qm ;;
  *)            "$QM" count "$GENOME" "$DNA" sample.qm ;;
esac

# ---- 3. EST: GC-normalize depth into COPY NUMBER per window ----
echo "[3/4] est (GC normalization -> copy number) ..."
"$QM" est "$GENOME" sample.qm sample.CN.bed   # col4 = estimated diploid-scaled copy number per window

# ---- 4. CALL copy-vs-allele for the O4 candidates from per-window CN ----
# CN ~2 over a candidate locus => single-copy gene => the divergent variant is an ALLELE/het.
# CN >=3 (and the candidate's divergent k-mers carry the extra depth) => a real COPY (paralog).
echo "[4/4] intersect candidate loci with parCN ..."
/home/juanfra/miniforge3/bin/python - "$OUT/sample.CN.bed" <<'PY'
import sys, csv, json, os
cnbed = sys.argv[1]
# load per-window CN
wins = []
for ln in open(cnbed):
    f = ln.rstrip().split("\t")
    try: wins.append((f[0], int(f[1]), int(f[2]), float(f[3])))
    except (ValueError, IndexError): continue
def cn_over(chrom, s, e):
    vals = [w[3] for w in wins if w[0] == chrom and not (w[2] <= s or w[1] >= e)]
    return sum(vals)/len(vals) if vals else None
struct = "bench/copy_vs_allele_structural.tsv"
rows = list(csv.DictReader(open(struct), delimiter="\t")) if os.path.exists(struct) else []
print("cid\tchrom\tstart\tend\tstructural_call\tparCN\tparCN_call")
copy=allele=0
for r in rows:
    cn = cn_over(r["chrom"], int(r["start"]), int(r["end"]))
    call = "n/a" if cn is None else ("COPY" if cn >= 2.5 else "ALLELE/HET" if cn < 2.5 else "?")
    copy += call=="COPY"; allele += call=="ALLELE/HET"
    print(f"{r['cid']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{r['call']}\t{cn if cn is not None else ''}\t{call}")
print(f"# parCN resolves: {copy} COPY (CN>=2.5), {allele} ALLELE/HET (CN~2) — settles the structural AMBIGUOUS set", file=sys.stderr)
PY
echo "DONE. sample.CN.bed = per-window parCN; the table above = per-candidate copy-vs-allele from DNA depth."
