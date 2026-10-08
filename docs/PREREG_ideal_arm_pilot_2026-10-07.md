# Pre-registration: the ideal arm of the expression ladder, pilot on the 2,334 Soto genes (2026-10-07)

**Status: written before any simulated read exists and before any number of this study. Descriptive pilot, no verdict class.**
Authorisation: the user's message of 2026-10-07 21:49 PDT, "can we try the pilot of the ideal arm?", after the question of 21:44 ("to prove that the problem is expression, can we show Soto's way with our RNA mode but using a dataset that has all their proposed genes as expressed with simulated reads ...").
It covers the ideal arm only (simulate, map, assemble, families, score) as bounded foreground heavy calls. It does not cover the multi-tissue and testis arms, which need their own pre-registration after the survey of 2026-10-07 (wf_ee69e34a-6d5).

## 1. Question

If every Soto gene is expressed, how many of them does our RNA mode (de novo assembly and the families stage, current default) recover as one-to-one nodes, by biotype, and which losses remain that expression cannot explain? The real-library baseline is the genome-wide A119b run of section 6w4/6w5 (`bench/SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md`): 765 of 1,078 Soto genes with a family id matched (71.0%): protein_coding 79.3%, transcribed_unprocessed_pseudogene 85.9%, unprocessed_pseudogene 46.6%, processed_pseudogene 42.1%, lncRNA 85.4%; 33.7% of the matched genes sit on a locus that covers two or more annotated genes. That baseline is from an older pipeline version and a different library, so a difference to it is **not** a like-for-like contrast; the like-for-like contrast (same harness, realistic profiles) is the later arms.

## 2. Design (fixed now)

| item | value |
|---|---|
| Genome and annotation | T2T-CHM13 v2.0 (`winloci_data/chm13v2.0.fa`) and the CAT/Liftoff v2.0 slim GFF (`winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz`); all 2,334 Soto gene ids (CHM13_G.../LOFF_G...) are gene ids in it (checked) |
| Genes | the unique `Gene ID` values of `bench/soto/soto_famCN_S1C.tsv` (2,334) |
| Windows | gene span +- 25 kb, merged; transcripts = every transcript of every gene (any biotype) overlapping a window, so expressed neighbours are present as in the earlier ideal runs |
| Reads | 10 per simulated transcript (ideal: equal depth for every transcript), the model of `bench/sim.py` (error .001, indel error/3, no truncation, 0 to 30 bp end jitter; ends must vary), canonical-intron recipe of `sim_windows.py`, seed 20261007, one replicate |
| Mapping | minimap2 with the command, `-N 50` and the whole-genome splice index of `bench/ideal_expression/run.sh` (read-disjoint parts, one part per call) |
| Pipeline | the driver's `assemble` and `families` stages with the HEAD release binaries (`/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release`), current defaults, nothing tuned, sha1 of every binary recorded in the run log |
| Scorer | `bench/ideal_expression/score_soto_nodes.py`, which wraps the matching rules of `bench/soto_family_locus_fidelity.py` (one-to-one exonic matching, fused loci reported separately) |
| Scored population | the Soto genes that were simulated; the family-id population of the baseline (1,078) is reported separately |

## 3. Outputs

Per biotype and overall: the one-to-one matched share, the fused-locus share (reported separately), median exonic Jaccard, |d5|, |d3| and width ratio of the matched genes, and a mutually exclusive ordered loss class per gene: matched; matched_but_fused; no_locus; wrong_extent (lowest 5% of the matched genes' Jaccard, the cut written to the output from the data). A table against the baseline above with the paired caveat of section 1. The estimand for each biotype is the **algorithm-limited share** = 1 - (matched share of the ideal arm), and, against the old baseline, the **gain** = ideal - baseline; no threshold turns either into a class.

## 4. Predictions (written before any outcome; directional, from the earlier ideal runs and the A119b natural experiment)

- P1. The ideal arm matches a larger share than the baseline in every biotype.
- P2. The gain is largest for unprocessed and processed pseudogenes and smallest for protein-coding genes; the gap between transcribed and untranscribed unprocessed pseudogenes (85.9% against 46.6% in the baseline) shrinks in the ideal arm, where both are expressed.
- P3. The ideal arm does not reach 100% in any biotype; the largest remaining class is matched_but_fused (the earlier NPIP ideal run lost 14 of 46 member nodes to a swallowed neighbour), then wrong_extent.
- P4. The fused-locus share of matched genes is at least that of the baseline, because expressed neighbours are simulated for every gene.
A failed prediction is reported as failed.

## 5. Gates (each prints PASS or FAIL and no outcome statistic; a table is INVALID unless all pass)

G1 manifest: reads equal the sum of per-transcript counts, every read name parses to a listed transcript, every Soto gene is simulated or listed with a reason. G2 read ends vary (the 5' and 3' offsets are not constant). G3 determinism: the simulation is byte-identical under `PYTHONHASHSEED` 0 and 1; the scorer's tables are identical under two hash seeds. G4 mapping sanity from the BAM alone, before the pipeline: the share of reads with a primary alignment, and the share of primary alignments that land on their source gene (a measure of paralog ambiguity that cannot be removed by any pipeline). G5 binaries: the stamps match HEAD.

## 6. Limits stated now

The reads are simulated from the annotation, so the truth for node extent is the annotation itself (circular for extent, independent for presence and fusion). Near-identical paralogs produce MAPQ-0 reads in the simulation as in real data, and the truth label stays with the source copy. Retrocopies, annotation errors and truly unexpressed copies are not represented. The comparison with A119b crosses pipeline versions and libraries. One replicate, one depth (10 reads per transcript); the node boundary is the k-th most extreme read end with k fixed at 2, so extent depends on depth, and a depth ladder belongs to the full study. No copy-number step and no family-level scoring in the pilot. Both halves of the Soto data are already spent for copy-number rules; this pilot tunes nothing.

## 7. Amendments

*(none yet; code hashes and the run calendar are appended here, below this line, before the first heavy call.)*
