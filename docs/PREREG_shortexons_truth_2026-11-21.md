# Row 1878: does the short-exon repair pay on complete truth and outside primates? (pre-registered 2026-11-21)

Lead from row 1876 (`docs/PREREG_whats_left_2026-11-21.md`, post hoc look `r1876_whats_left/probe.py`): among the wrong chains of
today's `final.gtf` on the simulated sets, 19 of 187 (fly), 121 of 787 (mouse) and 77 of 1,170 (zebrafish) have a new intron that
skips exactly one true exon; 100% / 92% / 83% of those exons are 15 bp or shorter. That is the aligner dropping a short exon, the
case the opt-in stage `shortexons` repairs (rows 1746, 1753, 1754, 1757; `docs/SHORT_EXONS_2026-11-16.md`). The stage was judged
on eight primate samples against RefSeq, under the defaults of 2026-11-16, and never on a complete truth. Its open problem then:
the new chains were 3-28% RefSeq-correct and RefSeq cannot say which are real. This row runs the stage as it is, at today's
defaults. Nothing is built or tuned. Workspace `r1878_shortexons_truth`. Registered before any run.

## Runs
Per dataset one job, driver = `git show HEAD:tools/rustle_pipeline.sh` at this commit (snapshot), binaries `~/rustle_target/release`,
`RUSTLE_PTC=0 RUSTLE_ONEREAD_TIER=0`:
- `base`: `assemble --bam B --fasta G --out base/x`
- `se`: `shortexons --bam B --fasta G --out se/x` (-> `se/x.shortexons.gtf`, `se/x.shortexons.final.gtf`)

| group | datasets (truth) |
|---|---|
| complete truth | fly, mouse, zebrafish simulated (rows 1779 / 1800 truths) |
| non-primate real | mouse ES (mouse contigs), heart, adrenal, cortex (GENCODE vM27) |
| primate real | OR6737 gorilla testis, KB3781 gorilla fibroblast (`apetools/GGO.refseq.gtf`), chimp testis (`cellline_asm/PTR.refseq.gtf`), WTC11 r2 human iPS (`ptr_ppa_rules/gffcmp/hsa.ref.gtf`) |
A119b is left out (the stage's work directory peaks near 180 GB there); HG002 likewise. OR6737 and KB3781 were development
samples of the stage, chimp testis was a check sample of row 1757; the seven non-primate datasets and WTC11 never saw it.

## Measures (multi-exon intron chains, `final.gtf` of `se` against `final.gtf` of `base`)
- chains, correct, precision of both; NEW = chains only in `se`, LOST = chains only in `base`, each with its correct count.
- On the complete truths: wrong chains whose new intron skips exactly one true exon of <= 15 bp, in `base` and in `se`.
- Stage output: exons fired, reads rewritten, wall time.

## Reading, fixed now
- The stage PAYS on a dataset when `se` has at least as many correct chains as `base` and at least `base`'s precision, and is
  strictly better on one of the two.
- **CANDIDATE** (for a default, the user's call) = pays on all three complete-truth sets, on at least 3 of the 4 non-primate real
  libraries and on at least 3 of the 4 primate ones, and no dataset loses more than 0.1% of its correct chains or 0.001 of precision.
- **TRUTH-ONLY** = pays on the three complete-truth sets and misses the real-library counts: then the misses on real libraries are
  reported with their new / lost tables, and the question "are the unannotated new chains real" stays open; no proposal.
- Anything else = no candidate; the stage stays opt-in as it is.

## Limits, stated in advance
- Simulated reads come from annotated transcripts only: a simulated set cannot contain an unannotated skip form, so the stage's
  known cost (an annotated skip form losing its reads to the exon form, row 1757) can show only where both forms are in the truth.
- The stage re-runs the whole assembly on a rewritten BAM; chains can change anywhere reads were rewritten.
- Real-library truths are incomplete; the stage's exons are unannotated by construction when they are new.
