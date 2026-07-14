# Illustrative example — the gorilla GSTM multi-copy family (real data)

A worked, corroboratable example on real gorilla Iso-Seq: the GST-Mu cluster on NC_073224.2. Figure:
`bench/slides/gstm_real_copies.png` (regenerate with `bench/make_gstm_copies_fig.py`).

## Reproduce

```bash
# in /home/juanfra/winloci_scratch (GGO_mm.bam + GGO.fasta)
copy_assign --bam GGO_mm.bam --fasta GGO.fasta \
  --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 \
  --dump-psv --gtf --out gstm_vg
python bench/make_gstm_copies_fig.py    # reads gstm_vg.psv_{reads,copies,cols}.tsv
```

## What it corroborates (for a skeptic)

1. **The 3 called copies ARE annotated paralogs.** copy 0/1/2 land exactly on the annotated GST-Mu genes
   **GSTM3 / GSTM5 / GSTM1** (positions from the copy tids vs `annot_intervals.tsv`). Anyone can check these
   are real, textbook recent gene duplications.
2. **They are genuinely similar** — grouped into ONE family only because their exon-sum identity clears the
   homology gate (≥ 80%). The clean pair GSTM5↔GSTM1 measures **83%**; GSTM3 is the divergent member (agrees
   with the others at only ~15% of PSV sites — consistent with GSTM3 being the outlier Mu paralog).
3. **The reads sort into 3 clean copy-blocks by their PSV alleles** (the SDA / Vollger read×PSV matrix, on
   real data): 411 PSV columns distinguish the copies; the read block below each copy-consensus matches it.
4. **The assignment is verifiable.** Where minimap2 maps a read *uniquely* (MAPQ > 0), our PSV-based copy call
   agrees with it **1341/1341 = 100%** — so on the MAPQ-0 reads where the aligner gives up, the method is
   doing what the aligner would if it could. 2654/2673 reads assigned; **16 certified TIED** (abstain, never
   1/k). Every copy carries hundreds of **private SUNs**, so a single read over a SUN deterministically pins
   its copy.

## The one-line takeaway

The variation graph is not a schematic: each PSV column is a bubble, each copy is a path, and the 2673 real
reads visibly thread onto the 3 copy-paths — assigned when a distinguishing bubble is significant, tied when
none is.
