# Milestone — Reference-absent gene-family copies from gorilla IsoSeq

**One line:** expressed, endogenous, divergent gene-family *second haplotypes* that are **absent from the
T2T gorilla reference** at their locus, found at the RNA level — strongest (protein-confirmed) in the MHC
and supported genome-wide in the recombination (PRDM9) and zinc-finger families. The RNA-level analog of
Eichler's novel TBC1D3 groups and Soto's paralogs missing from GRCh38. **Caveat up front:** from RNA
alone "reference-absent *copy*" vs "hyperdivergent *allele*" is not separable — these are strong
*candidates* pending DNA copy-number (parCN); see Honest caveats.

## The question (advisor's standing interest)

"Find copies of a gene that are *not in the genome*." A copy expressed in the individual but absent
from the reference assembly can never align correctly, so an alignment-keyed method is blind to it. The
question: can we recover such copies from gorilla IsoSeq against the (near-complete) T2T reference?

## Approach — two populations, and where they actually live

A copy absent from the reference surfaces in one of two read piles, and we tested both:

1. **Unmapped reads (divergent-absent).** Tested directly (`unmapped_rescue_poc.md`): the unmapped pile
   is tiny (0.13 %), 79 % already-present (99.7 % identity), residual = expressed transposable elements,
   leaving **one** genuine hit (`cl0`). On a *complete* T2T reference this route is **dry**. Its own
   conclusion told us where to look: a divergent/CNV copy with residual homology gets **force-mapped
   onto its paralog**, not flagged unmapped — so the gene-family copies live in the **mapped** pile.

2. **Mapped reads (collapsed / CNV-absent) — the yield.** A reference copy that is really ≥2 collapsed
   copies (or carries an extra CNV paralog the reference lacks) shows, among its **primary** alignments
   (the paralog-bleed firewall), a coherent **second haplotype**: ≥12 balanced-frequency (0.20–0.60)
   alt columns co-segregating on ≥5 reads (`hidden_copy.rs::detect_hidden_copy`). A diploid het is 1–2
   columns; sequencing error never reaches the balanced band.

## The pipeline (detect → validate → promote → lock)

Run directly on the copy loci (19 s for 357 family loci; the de-novo + POA path was 100× slower and was
abandoned). Each flag is then put through a funnel that each removes a distinct false-positive class:

- **RNA-editing filter** — A→I (read as A>G) makes co-segregating columns; a genuine copy has a
  *diverse* substitution spectrum, editing is A>G/T>C-dominated. (Load-bearing: removed 5/30 family
  flags, incl. a 67 %-A>G GSTM signal.)
- **Assembly + genome-divergence** — POA-assemble the alt-haplotype reads into the copy's consensus,
  map to the whole T2T genome. The decisive signal is **divergence in a band**: the consensus maps back
  to its own locus at **3–20 % divergence** — too divergent to be a het allele (<0.5 %) or the copy
  itself, too conserved to be a chimera (>20 %). (The >20 % band is consensus chimera/repeat — the
  DSFAM213 failure mode.)
- **Protein BLAST** — blastx vs the 22,614-protein gorilla proteome: names the paralog and **rules out
  contamination** (endogenous gorilla protein hit, e ≤ 1e-57 — the check `cl0` failed and these pass).
- **Multi-mapping** (genome-wide) — a real copy of a *multi-locus* family maps to its own divergent
  locus **plus** the reference paralogs; an in-place second haplotype maps to one locus only.

## Result

| tier | n | what | confidence |
|---|---|---|---|
| **MHC copy candidates (protein-confirmed)** | **4** | class I (Gogo-B\*0201), DRB1, DQ-α, DQ-β — 3.9–20.4 % divergent, ORFs 97–362 aa, endogenous (86–98 % gorilla-protein), clustered in the MHC (NC_073229.2 49–50 Mb) | high protein homology; copy-vs-allele needs DNA |
| **Dispersed-paralog copies** | **15** | genome-wide; map to own divergent locus + paralog loci — **PRDM9** (13.5 %, distinct member of the recombination-hotspot family), **ZNF208**, ZDHHC13, large multi-locus families | medium |
| single-locus / family-CNV | 66 | het / tandem-CNV / absent-retrocopy / depth-doubled — RNA cannot separate | needs DNA parCN |

**Funnel:** 12,793 expressed loci → 1,015 flags (8 %) → 73 paralog candidates → 15 dispersed-paralog
copies; plus the 4 family-scoped MHC. The copies land **exactly where biology predicts** — immune
(MHC), recombination (PRDM9), and zinc-finger / segmental-duplication families — the most polymorphic,
copy-number-variable, reference-divergent regions of the genome.

## Why it matters (thesis through-line)

This is the **identifiability thesis turned outward**: the reference cannot *model* these copies (one
gene where the reads imply more), and we detect that from the reads while **abstaining from placing the
copy** (DAZ3/Canzar discipline — detect + flag the evidence, never fabricate a sequence). It complements
the within-reference copy-assignment work: there we assign reads to the copies the reference has; here we
flag copies the reference lacks.

## Measured specificity — false-positive BOUND (`bench/o4_fp_bound.py`)

Quantified (no BAM re-scan; intersect the genome-wide flags with single-copy genes). On **17,334
single-copy genes** (unique-name protein_coding, non-LOC — definitionally cannot host a reference-absent
*copy*), the raw hidden-copy flag still fires on **828/11,206 = 7.39%** of overlapping loci — essentially
the genome-wide background rate (7.93%). **The raw flag is a non-specific SCREEN, not a specific copy
detector**: ~7.4% is an upper bound on its false-positive rate, dominated by **heterozygosity** (a
diploid gene with enough balanced SNP columns trips the ≥12-column threshold). This is the het-vs-copy
wall, measured. *Implication:* the 1,015 genome-wide flags are a candidate screen, not a copy catalog;
specificity comes from the DOWNSTREAM filters (divergence >3% ≫ het ~0.5%; RNA-editing filter; protein
BLAST), and only the **4 protein-BLAST-confirmed MHC candidates** survive an external check. Do not
assert the flagged pool as a copy set.

## Honest caveats (the loose ends)

- **Het vs copy (the central one) — now quantified at ~7.4% raw-flag FP (above).** For the
  hyperpolymorphic MHC and the 66 single-locus candidates,
  "reference-absent paralog *copy*" vs "hyperdivergent *allele*" is genuinely blurry from RNA alone. The
  resolver is **DNA copy-number (parCN)** — the concrete experimental ask. The finding holds regardless:
  expressed sequence 4–20 % divergent from the reference's gene, endogenous and coding.
- **No independent validation yet.** No DNA, no second individual, no assembled full copy locus. The 4
  MHC are protein-confirmed *endogenous* but not *cross-validated* as copies.
- **Genome-wide is a candidate-generator, not a clean catalog.** Precision drops off the
  family/segdup-scoped regions; the 15 dispersed-paralog are multi-mapping-supported candidates, the 66
  single-locus need DNA.
- **Intronless/retrocopy detection is not possible here** — an absent retrocopy's reads map to the
  parent and splice out its introns, so they look like a normal spliced transcript; the intronless
  signature only appears where the retrocopy's own locus is in the reference (absent by definition).

## Reproduce

`hidden_copy_scan.py` · `hidden_copy_validate.py` · `promote_hidden_copies.py` ·
`promote_evidence_plus.py` · `hidden_copy_scan_genomewide.py` · `promote_genomewide.py`. Consensus
sequences + evidence: `/home/juanfra/winloci_scratch/refabsent/`. Full technical writeup:
`bench/reference_absent_catalog.md`.
