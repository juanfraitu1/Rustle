# From reads to results — the four objectives (real-data examples)

Two figures for the "how does it work / does it work" conversation:

- `bench/slides/objectives_flow.png` (`make_objectives_flow.py`) — the **flow**: reads → variation graph
  → result, one row per objective, real numbers.
- `bench/slides/gstm_real_copies.png` (`make_gstm_copies_fig.py`) — the **raw evidence** behind O2: the
  2673 GSTM reads sorting into 3 copy-blocks by their PSV alleles (see `GSTM_REAL_COPIES.md`).

## What each row shows (and how to reproduce)

All from `GGO_mm.bam` + `GGO.fasta` in `/home/juanfra/winloci_scratch`.

**O1 — define the family + count copies (GSTM).** Reads at 3 mutually-homologous loci → one family
(γ-quasi-clique of the transcribed-homology graph) → copy number **χ_H = 3**, and the 3 copies are the
annotated paralogs GSTM3 / GSTM5 / GSTM1.
```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --region NC_073224.2:129160000-129230000 \
            --homology-primary --min-copies 2 --dump-psv --out gstm_vg     # families.tsv → n_copies=3
```

**O2 — assign each read to a copy (GSTM).** A read threads its alleles through the PSV bubbles; if a
distinguishing bubble is significant it is **assigned** (example reads span 414 / 347 / 328 decisive PSVs,
p ≈ 0), and a read spanning **0** distinguishing bubbles is certified **TIED** (abstain, never 1/k).
2654 assigned · 16 tied · **unique-mapper agreement 1341/1341 = 100%** (same run as O1; `assignments.tsv`).

**O3 — allele-specific junctions (ABCC4).** On single molecules, the **T allele never uses** the junction
at donor 109,947,543 (0/173), while the **C allele uses it 38%** (45/117): ΔPSI = 0.38, q = 3×10⁻¹⁸ —
allele-linked splicing with no phasing.
```
asj --bam GGO_mm.bam --region NC_073238.2:109943482-109971515 --out abcc4   # asj.tsv; verified in asj_calls_verified.tsv
```

**O4 — a copy not in the reference (cross-species check, human MAGEA).** When the reads resolve more copies
than the reference annotates (χ_H or depth > reference), the surplus is a reference-absent / collapsed copy.
Concrete corroboration: on human testis Iso-Seq the same method resolves **11 MAGEA copies** against the
2-copy ape baseline, matching the T2T annotation (`HUMAN_CROSSSPECIES.md`). Detect-and-flag; separating a
real extra *copy* from a heterozygous *allele* needs DNA / parental copy number.

## The one-line story

The pipeline is one thing seen four ways: reads become a **variation graph** (spine + PSV bubbles); O1 reads
it as *which loci are one family and how many copies*, O2 as *which copy each read is*, O3 links *allele to
junction* on the same molecule, and O4 flags *a copy the reference doesn't show*.
