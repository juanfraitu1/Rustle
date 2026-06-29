# Can we find UNANNOTATED copies (not even a LOC)? — yes, at three levels of "unknown"

The advisor's question: not reference-absent copies (those mapped to *annotated* MHC/etc. loci), but
expressed multi-copy-family members at loci with **no annotation at all** — not even an uncharacterised
`LOC######`. Checkable by intersecting the de-novo expressed loci with the *entire* gorilla annotation
(gene/lnc_RNA/pseudogene/mRNA/tRNA/…, 30,191 spans, 1,510 Mb) and asking which unannotated ones are
paralog copies. **Answer: yes.**

## The funnel

- **12,793 de-novo expressed loci → 291 unannotated** (2.3%; overlap *zero* annotation).
- Of the 291: **42 are paralog copies** (transcript maps to ≥2 loci at ≥80% identity). Copy-count
  2–8 loci (31 are 2–4, 11 are 5–8) — **gene-family-like, none >8 = no dispersed-TE contamination**.
- blastx vs the gorilla proteome: **203/291 unannotated transcripts hit a protein** (most "unannotated
  loci" are coding paralogs the annotation simply missed at that spot); 88 have no protein hit.

## Three levels of "unknown"

1. **Unannotated copies of KNOWN families (the robust win).** ~35–36 paralog loci whose paralog *and*
   protein are annotated — a copy of a known gene sitting at a locus with **no annotation, not even a
   LOC**. Concrete: **ROBO2** — a 94%-identity copy family at **5 unannotated loci** (likely
   retrocopies / segmental duplications of the axon-guidance gene); ZDHHC13 paralog (53%); several
   `LOC#####` copies (75–81%). Plus **6 read-conflict-VALIDATED unannotated copies inside the de-tie
   families** — including **DSFAM817**, the difficult-case flagship win: *we assigned reads to a copy
   that isn't even annotated.*
2. **Loci whose paralogs are also unannotated, but the protein family is known.** e.g. **ERVFRD-1**
   (endogenous-retrovirus env / syncytin-2 family, 31% — a retroelement family, not a classic gene),
   and copies of annotated `LOC` genes at unannotated loci.
3. **LITERALLY unknown (the rarest).** **NC_073240.2:23.57 Mb** — a **2-copy** family with **no
   annotation, no protein homology, no LOC** — genuinely uncharacterised. The broader pool of **88
   unannotated loci with no protein hit** are the candidates for this class (novel genes / lncRNA /
   to-be-vetted).

## Honest caveats

- "Unannotated copy of a known family" is the solid, common result (40+ loci, gene-family copy-counts,
  protein-confirmed, 6 read-conflict-validated). This robustly answers "can we find unannotated copies."
- "Literally nothing known" (no annotation *and* no protein *and* no LOC) is **rare** — 1 clear
  multi-copy case here, plus the 88-no-protein-hit pool that needs ORF/coherence vetting (novel gene vs
  lncRNA vs assembly/mapping artifact). Some (ROBO2 ×5 @94%) are likely processed pseudogenes/segdups;
  ERVFRD-1 is a retroelement — both are "copies" but not novel *gene* families.
- The annotation is gorilla RefSeq (fairly complete: 34k genes incl. 7k pseudogenes, 10k lncRNAs), so
  "unannotated" here is a genuine gap, not a sparse-annotation artifact.

## Bottom line for the advisor

**Yes — we find expressed paralog copies at loci the gorilla annotation does not cover at all (no LOC):
42 genome-wide (ROBO2 ×5, ZDHHC13, LOC-copies…) + 6 read-conflict-validated (incl. the DSFAM817 win
copy).** Most are unannotated copies of families whose protein is known elsewhere; a small residue
(1 clear multi-copy + an 88-locus no-protein pool) is the "literally unknown" frontier worth vetting.

## Vetting the 88 no-protein-hit pool (`vet_unannotated_novel.py`)

Separating novel gene families from lncRNA/noise (ORF, exon count, paralog count, 4-mer complexity):
- **4 NOVEL multi-copy gene-FAMILY candidates** (multi-copy + ORF ≥100 aa + spliced + complex, no
  annotation, no gorilla-protein hit): **NC_073236.2:139.89 Mb — a 5-copy family, 183-aa ORF, spliced,
  complexity 0.96** (the flagship); NC_086018.1:20.93 Mb (2 cp, 157 aa); NC_086018.1:30.31 Mb (2 cp,
  142 aa, 3 ex); NC_073240.2:23.57 Mb (2 cp, 117 aa).
- 50 novel **single-copy** gene candidates (ORF + spliced, 1 locus — novel genes, not families).
- 15 **lncRNA** candidates (spliced, no/low ORF). 19 ambiguous. 0 low-complexity/noise (complexity
  filter clean).

**Caveat — "truly novel" needs a cross-species check.** "No protein hit" is vs the *gorilla* proteome
only; these 4 could still be homologs of a protein in another species (= gorilla-unannotated, not
universally novel). The airtight confirmation is blastp/diamond vs SwissProt/nr — the next step. Reads
are modest (3–11/locus); the 183-aa ORF is well above the ~21-aa random expectation, so coding potential
is real; complexity 0.96–1.0 + spliced + multi-copy argue against artifact.

Artifacts: `unannotated_loci.json` · `unannotated_copies.json` · `unannotated_novel_vetted.json` ·
`unannot_blastx.tsv` · `unannot_tx.{fa,paf}` (in `/home/juanfra/winloci_scratch/refabsent/`).
