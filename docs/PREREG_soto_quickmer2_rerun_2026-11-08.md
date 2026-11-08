# PREREG row 1588: Soto 2025's QuicK-mer2 parCN pipeline, run end to end as released, on four 1KGP genomes

Written 2026-11-08, before any QuicK-mer2 output exists. Advisor's point: show the published method can be followed from
the authors' instructions, start to finish. Register 1167 (2026-09-29) stopped here on memory (index needs ~52 GB; the
laptop had 18 GB). Fritz nodes have 256 GB.

## What is run (as released, `github.com/mydennislab/HSD_brain_evolution`, `section_I&II/QM2_pipeline/`)

- Their `quickmer2.smk` with their config: `distance=1 window=500`, reference **CHM13 draft v1.0**
  (`chm13.draft_v1.0.fasta`, human-pangenomics S3), their control BED `control/chm13.draft_v1.0.control.bed`
  (6,222 autosomal regions). QuicK-mer2 = Kidd Lab master `a81d4e8` (the repository their readme links; their lab's
  installed build is not public).
- One edit, recorded: the CRAM decoding reference is hard-coded to their cluster
  (`/share/dennislab/.../GRCh38_full_analysis_set_plus_decoy_hla.fa`); it points to the same file downloaded from the
  1KGP FTP.
- Second edit, recorded 2026-11-08 after the first `search` attempt and before any output: the released `search` line
  puts `-e 1 -w 500` after the FASTA (`search -t N -c control.bed ref.fasta -e 1 -w 500`). Public QuicK-mer2 reads its
  file arguments from the end of argv (`QuicKmer.c`: `argv[argc-1]`), so it took `500` as the genome and segfaulted (job
  4452893; reproduced on chr21). The same options are moved before the FASTA; chr21 then indexes in 26 s. Their lab build
  may parse differently. This is a finding about following their instructions, not a change of method.
- Third edit, recorded 2026-11-08 after the first `count` attempt (job 4452959) and before any CN output: the CRAM
  decode option `--input-fmt-option required_fields=0x202` is written as `514`, the same value in decimal. htslib 1.25
  parses this option as a decimal integer (`[E::get_int_kmg] Unrecognised suffix 'x202'`), so samtools printed its
  usage and `count` received no reads. Nothing else in the Snakefile changes.
- Environment fix, recorded 2026-11-08 after the first complete run (job 4452964) and before scoring any CN value
  that depended on it: QuicK-mer2's GC smoother `smooth_GC_mrsfast.py` uses `numpy.float` (removed in numpy 1.24) and
  crashed. `est` does not check that pipe and wrote CN of about 0 everywhere (chr1 mean 1e-5). The first scoring of
  that broken output (`qm2_compare.py`: Bar A 0.000, Spearman -0.07) was a broken run, not a result. Fix: the helper runs under
  python 3.10 / numpy 1.23.5 / pandas 1.5 / matplotlib (`qm2/envpy`, first on PATH); no QuicK-mer2 code changed.
  `count` output is kept; only `est` and bigBed rerun.
- Samples: one 1KGP high-coverage genome per population, the first listed in
  `1000G_2504_high_coverage.sequence.index`: NA18488 (YRI), NA18530 (CHB), NA12718 (CEU), HG02285 (PEL).
- Per-gene parCN: Soto's `genotype_cn_parallel.py` is not public (register 1167). Its public sibling is an unweighted
  mean of window CNs; we take the mean CN of the windows overlapping the gene's `SD98_v1.0` coordinates (S1E column
  12), per sample. This is the one step we must reconstruct.

## Comparison with Soto's Table S1E (2,504 genomes; we have 4)

- **Bar A (method reproduces, fixed genes):** among S1E `parCN class` = Fixed genes with a window, >= 90% of
  gene x sample values lie within 0.5 of S1E `Median parCN`.
- **Bar B (rank agreement):** Spearman between our 4-sample median and S1E `Median parCN`, all genes with a value,
  >= 0.80.
- Reported, no bar: Nearly-Fixed and Polymorphic agreement, per-sample results, genes with no window, wall time and
  peak memory of `search` and `count`.

PASS = both bars. A FAIL is read as: the released pipeline does not reproduce S1E from the authors' instructions
(their per-gene aggregator, their QuicK-mer2 build or the sample set being the named candidates).
Work dir: `/anvme/workspace/v138da13-rustle/qm2/`.

## Outcome (2026-11-08): FAIL as registered (Bar A PASS, Bar B FAIL)

Scored with `bench/soto/qm2_compare.py` on the corrected `est` output (job 4453019; chr1 mean CN 2.02). All 1,833 S1E
genes have windows.
- Bar A: Fixed genes, 0.982 of 1,360 gene x sample values within 0.5 of S1E Median parCN (bar 0.90): PASS.
- Bar B: Spearman of our 4-sample median vs S1E Median parCN 0.760 (bar 0.80): FAIL.
- Reported: Nearly-Fixed 0.711 (3,912 values), Polymorphic 0.566 (2,060).
- Post hoc, not part of the verdict: Pearson 0.786. Median absolute difference 0.046 (Fixed), 0.136 (Nearly-Fixed),
  0.234 (Polymorphic). 66% of S1E medians lie in 1.5-2.5, so with 4 genomes against 2,504, ranks inside that band are
  mostly noise.

Reading: the released pipeline does run end to end from the authors' instructions, after four mechanical fixes (listed
above, none a change of method). It reproduces S1E values on Fixed genes. The registered rank-agreement bar fails.
Outputs: `/anvme/workspace/v138da13-rustle/qm2/run/output/{CN,bigBed}/`, per gene `qm2/r1588_per_gene.tsv`.
