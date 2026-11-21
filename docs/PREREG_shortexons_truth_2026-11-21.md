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

## Result (2026-11-21): TRUTH-ONLY as registered (primates 2 of 4); pays on 9 of 11, the two gorillas within 4 chains

Jobs 4494570-4494580, every arm exit 0; `r1878_shortexons_truth/score.txt` (`score.py`), driver `f0b9fbd3`, binaries of
2026-11-20. `base` equals today's default numbers of row 1849 on every dataset.

`final.gtf`, chains / correct / precision; NEW and LOST with their correct counts:
| dataset | base | with the repair | new | lost | pays |
|---|---|---|---|---|---|
| fly simulated | 6,940 / 6,753 / 0.9731 | 6,938 / 6,762 / 0.9746 | 11 (10) | 13 (1) | yes |
| mouse simulated | 25,864 / 25,077 / 0.9696 | 25,826 / 25,141 / 0.9735 | 76 (72) | 114 (8) | yes |
| zebrafish simulated | 13,300 / 12,129 / 0.9120 | 13,271 / 12,147 / 0.9153 | 30 (20) | 59 (2) | yes |
| mouse ES | 41,679 / 17,074 / 0.4097 | 41,669 / 17,093 / 0.4102 | 122 (24) | 132 (5) | yes |
| mouse heart | 21,176 / 12,923 / 0.6103 | 21,183 / 12,950 / 0.6113 | 87 (30) | 80 (3) | yes |
| mouse adrenal | 23,042 / 14,163 / 0.6147 | 23,046 / 14,184 / 0.6155 | 64 (26) | 60 (5) | yes |
| mouse cortex | 29,787 / 15,491 / 0.5201 | 29,803 / 15,526 / 0.5210 | 277 (61) | 261 (26) | yes |
| OR6737 gorilla testis | 57,830 / 24,315 / 0.4205 | 57,859 / 24,311 / 0.4202 | 163 (13) | 134 (17) | no |
| KB3781 gorilla fibroblast | 53,153 / 23,407 / 0.4404 | 53,175 / 23,406 / 0.4402 | 128 (12) | 106 (13) | no |
| chimpanzee testis | 46,969 / 20,309 / 0.4324 | 46,976 / 20,317 / 0.4325 | 111 (18) | 104 (10) | yes |
| WTC11 r2 human iPS | 31,336 / 20,016 / 0.6388 | 31,345 / 20,054 / 0.6398 | 82 (46) | 73 (8) | yes |

- **Verdict as registered: TRUTH-ONLY.** Complete truth 3 of 3, non-primate real 4 of 4, primate 2 of 4 (bar 3 of 4). No
  dataset breaks the harm bounds (largest loss: OR6737, 4 correct chains = 0.016%, precision -0.0003).
- **The question row 1757 left open is answered on complete truth:** the new chains are true. Mouse simulated 72 of 76, fly 10
  of 11, zebrafish 20 of 30; the lost chains are wrong ones (8 of 114, 1 of 13, 2 of 59 correct). Wrong chains that skip one
  true exon of <= 15 bp: fly 19 -> 11, mouse 108 -> 31, zebrafish 62 -> 17. `x.gtf` gives the same table.
- On the real libraries outside gorilla the new chains are 3-9 times as often annotated as the lost ones (WTC11 46 of 82
  against 8 of 73).
- **The two gorillas (post hoc, not registered).** Against the second gorilla annotation (CAT, `r1735_ggo_cat/cat.gff3`
  converted with gffread) both pay: OR6737 new 38 correct, lost 13 (+25 correct chains, +0.0003 precision); KB3781 new 35, lost 7
  (+28, +0.0004). The miss is specific to the gorilla RefSeq models. That is one annotation against another, not a verdict.
- Stage output: 32-567 exons fired per dataset, 1,504-76,728 reads rewritten. Wall time on 48 threads against the plain assembly
  on 24: 135-769 s against 71-281 s on nine datasets; fly 1,351 s (1,870 contigs, one temporary BAM each) and KB3781 1,232 s.
  The stage assembles a second time; folded into `assemble` the cost would be the rewrite alone.

**Reading (mine).** This is a repair, not a filter: a wrong chain is exchanged for its true one, so both axes rise, and the
effect is as large outside primates as inside (it follows the aligner, not the species). It is small: +0.04% to +0.26% correct
chains, +0.0001 to +0.004 precision. By the registered reading it is not a default candidate; whether it becomes one against
that reading, and whether it is folded into `assemble` (one assembly on the rewritten reads) or stays a second pass, is the
user's call. Not judged here: A119b and HG002 (left out for the work directory size), the family stages on rewritten reads.
