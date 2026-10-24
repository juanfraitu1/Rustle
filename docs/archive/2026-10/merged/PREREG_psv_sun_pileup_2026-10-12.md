# PREREG — SUNs vs our PSVs (M1) and a per-copy pileup check of PSVs (M2), gorilla KB3781 (2026-10-12)

**Why.** The advisor (a) calls PSVs "SUNs" and (b) argues pileup-style variant evidence would tell real variants from
errors better than our alignment-derived PSVs. M1 measures how our O2 PSVs relate to SUNs; M2 measures whether a
per-copy pileup AFTER assignment finds places where this animal does not carry the reference PSV allele, with
KB3781's own other haplotype as truth. Descriptive study; bars below say what each outcome means.

## Substrate and run (fixed before any output is read)

- Reads: KB3781 fibroblast IsoSeq `fibroblasts/GCA_029281585.2_flnc_mm.bam`; genome `winloci_scratch/GGO.fasta`
  (mGorGor1 v2 primary, RefSeq names); catalog = the 2026-10-10 de novo run
  `/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/` (`fam.copies.tsv`, `fam.copies.fa`, 326 families).
- **Subset S (56 families):** rows of that run's `assign.families.tsv` with `psv_cols >= 1`, `n_reads >= 20`,
  `n_copies <= 20` (chosen for resources — the full assign is 9.4 min / 5.3 GB; no outcome of this study is known).
- **Binary:** snapshot of `rustle_target_m2/release/copy_assign` (md5 `8b0026b9b82de0878743814be1b01297`, the binary
  of the 10-10 run), same arguments as that run's driver call plus `--dump-psv --threads 2`, families filtered to S,
  regions = the 10-10 `regions.txt` lines overlapping an S copy. Run under `tools/rlock.sh light` with a 4-GB
  virtual-memory cap (no interference with other sessions: no cargo build, own output dir `/mnt/linuxdisk/tmp/psv_sun/`).
- **Sanity S1:** the subset run's (read, family) rows agree with the 10-10 run on status and copy for >= 95% of shared
  rows (differences allowed only for reads shared with families outside S). Fail -> stop and diagnose.
- **Gate V (coordinates):** `psv_col_pos` gives one copy's genome position per column; the other copies' positions
  are derived by aligning their `fam.copies.fa` spliced sequences to that copy's (`minimap2 -cx asm20 --cs`) and
  mapping offsets through each copy's exons and strand. A derived position is used only if the genome base there
  (complemented on `-`) equals the dumped allele. Gate: >= 95% of (copy, column) pairs with a non-gap allele verify;
  else stop and fix before reading M1/M2.

## M1 — how many of our decisive PSVs are SUNs

- **Private allele:** at a decisive column, an allele carried by exactly one copy of the family (gaps excluded).
- **SUN** (our operational form of Sudmant 2010): a private allele with >= 1 overlapping genomic 30-mer (taken from
  that copy's genome position, canonical) that occurs exactly once in GGO.fasta. Counted with a purpose-built k-mer
  scan (candidate k-mers only; no genome-wide k-mer database — disk).
- **Reported (per family and totals, with n):** decisive columns; columns with >= 1 private allele; SUN alleles /
  private alleles; **G = private alleles with NO genome-unique 30-mer / private alleles** (the allele+context exists
  elsewhere in the genome = an outside-the-catalog sequence our catalog-relative PSV definition cannot see);
  **R = O2-assigned reads that carry >= 1 SUN allele of their assigned copy / O2-assigned reads** (what a SUN-only
  method could assign; n ~ 529 rows on 10-10, so R is indicative only).
- **Meaning:** G >= 0.05 -> recommend annotating O2's PSVs with genome-uniqueness (a catalog-private allele is often
  not genome-private); G < 0.05 -> the catalog-relative definition loses little. R is descriptive (how much the
  non-private PSVs add over SUN-only evidence).

## M2 — per-copy pileup check of PSVs

- For copy X and decisive column j (verified position g): the per-copy pileup = (a) reads PLACED on X by the aligner —
  primary alignments (`-F 2308`) with MAPQ >= 20 covering g (`samtools mpileup -q 20 -Q 0 --ff 2308 -B`, per copy
  window, bases in genome orientation; ref-skips and deletions not counted) — plus (b) reads O2 ASSIGNED to X in the
  dump (only 529 such rows in S on 10-10, hence (a); the default AS-tied gate never certifies uncontested reads).
  Pools are reported separately; flags use their sum.
- **Flagged:** >= 10 covering reads, a single non-X allele in >= 3 reads and >= 0.20 of them (thresholds of the
  RNA allele caller, register row 1194). **Control:** >= 10 covering reads, X's allele in >= 0.95.
- **Truth:** X's locus (genomic span +-2 kb) aligned to the OTHER haplotype's chromosome (`rna_allele/chrmap.tsv`
  `B_name`; the primary chromosome is identical to the same-name haplotype) with `minimap2 -cx asm5 --cs`; counterpart
  = one alignment covering >= 95% of the query. Classes: **B_ALT** (other haplotype's base = the reads' alt allele),
  **B_OTHER** (differs from X's reference allele, not the alt), **B_SAME** (= X's reference allele), **NO_COUNTERPART**.
  One minimap2 call per chromosome, each as a light job.
- **Bar (pileup check carries real signal):** among flagged pairs with a counterpart, B_ALT >= 0.50 AND the flagged
  B_ALT+B_OTHER rate >= 5x the control rate. Otherwise flags are mostly not the animal's own sequence (mis-assignment,
  editing, error).
- **Negative control:** copies on chrX (male, no second X: `B_name` empty) — reported flag rate per covered pair vs
  autosomes; het-driven flags cannot exist there.
- Also reported: transitions share among flags (A<->G editing; O2 already filters within-copy A<->G columns), and the
  share of O2-assigned reads crossing a flagged column.

## Not claimed

One animal, one tissue, 56 families; SUN = our operationalization (30-mer, primary assembly only); M2 truth tests only
the other haplotype (the primary IS one haplotype). Code `bench/psv_sun/`; outputs `/mnt/linuxdisk/tmp/psv_sun/`.

## Amendment 1 (2026-10-12, after gate V failed; before any M1/M2 output existed)

- **S1 PASS** (26,721 / 26,821 = 99.6% of the 10-10 rows identical). **Gate V FAIL**: 118,169 / 216,593 (54.6%)
  (copy, column) pairs verify (`gate_v.json`), after the mapper was fixed to use every alignment piece on either strand.
- **Diagnosis:** (i) 10 of 56 families contain copies whose spliced transcripts are REVERSE-COMPLEMENT homologs of the
  anchor's (inconsistent strand calls on 1-2-exon de novo copies). Their dumped alleles match the homologous base 26%
  and its complement 29% (chance level): O2's column frame pairs NON-HOMOLOGOUS bases there. Median decisive-column
  density 723 / kb vs 9.8 / kb in the 46 strand-consistent families. (ii) In strand-consistent families 37% of pairs
  are unmapped: no local alignment covers the anchor's position (non-homologous segments forced into the global frame,
  or divergence beyond `asm20`). Reported as a separate finding; not fixed here (O2 code is another session's).
- **Amended population:** a decisive column is **homologous** iff every copy carrying a non-gap allele has a VERIFIED
  position. M1 and M2 are computed over homologous columns only, and the share of O2's decisive columns that are
  homologous is reported (overall and by family strand-consistency). Gate V is replaced by this restriction; no other
  definition or bar changes.
