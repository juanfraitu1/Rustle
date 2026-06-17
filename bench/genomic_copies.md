# Genomic self-alignment — completing the copy roster (pseudogene + unannotated copies)

The DNA/protein tier (mmseqs2 on translated CDS) catches ancient *coding* families but is blind to two
copy classes: **pseudogene copies** (no CDS → no protein) and **unannotated copies** (no annotation at
all). This pass finds them by genomic self-alignment: map every gene-rep transcript back to the genome
(`minimap2 -cx splice -N 20 -p 0.4 --secondary=yes`) and classify each strong homology hit (≥50% of
the transcript aligned at ≥80% identity) by the annotated feature it lands on.

## Result
PAF hit classes: self=23,123 · annotated_gene(paralog)=18,495 · **pseudogene_copy=2,825** · **unannotated=1,949**.

Merged to distinct loci: **1,978 NEW copy loci** the annotation + protein-clustering missed:
- **pseudogene copies: 1,313** (annotated as pseudogenes — the protein-clustering blind spot)
- **unannotated copies: 665** (no annotation at all — the annotation blind spot)

## RNA overlay — what transcribes (the functional readout)
- new copies are **mostly silent**: pseudogene copies 18% transcribed (239/1,313), unannotated 16%
  (104/665) — vs **72%** for coding-family copies. Exactly the expected biology: duplicated/pseudogenized
  copies are largely dark, with a transcribed minority.
- the transcribed minority is the interesting part — e.g. an **SDHA retrocopy** (335 reads, 94.8% id),
  and several high-identity LOC copies expressed at 180–390 reads. These are real transcribed copies
  invisible to both the annotation and the protein tier.

## The complete picture (three tiers + RNA overlay)
| copy class | source | copies | transcribed |
|---|---|---|---|
| annotated coding-family | mmseqs2 protein clusters | 14,545 | 72% |
| pseudogene copy | genomic self-align → annotated pseudogene | 1,313 | 18% |
| unannotated copy | genomic self-align → outside all annotation | 665 | 16% |

**DNA enumerates every copy** (annotated coding + pseudogene + unannotated); **RNA shows which are
actually transcribed**, copy-by-copy. The silent vs transcribed split *is* the genome-vs-transcriptome
result, and copy-assignment (the identifiability theorem) then says which of the transcribed copies are
individually resolvable.

## Honest scope
- Copies found by mapping ANNOTATED transcripts → genome, so a fully novel family with NO annotated
  member anywhere would still be missed (would need de-novo transcript assembly as the query).
- ORF-intactness of the new copies (true pseudogene vs intact-but-unannotated gene) is the natural next
  annotation — read the ORF off the transcribed copies' reads.
- Identity/coverage thresholds (≥0.8 id, ≥0.5 qcov) set the recent/divergence horizon for "a copy."

## Reproduce
- annotated intervals: extract `gene`/`pseudogene` features (chrom,start,end,biotype,gene) from
  `/tmp/gw/ref_*.gff3` → `annot_intervals.tsv`
- `MINIFORGE python bench/build_genomic_copies.py` (minimap2 self-align + classify + RNA overlay)
- `python3 bench/genomic_copies_fig.py`
