# PREREG — does O2 recover the right copy? Truth from the read's own unique anchor (2026-10-13)

**Why.** The user asked for examples that prove O2 works. Existing evidence:
- O2's calibration (it never assigns undecidable reads; YAG Arm A, chr16 simulation 157/157);
- haplotype phasing on the assembly's own animal (ASJ, T1 discordance .001 held-out).

Real-read ACCURACY has never been measured with a truth independent of O2's own evidence. Simulated reads come from
the reference copies, and real reads have had no truth. The user chose three truth sources, run in this order, each
pre-registered before it runs:
1. anchor trim (this document);
2. copy-specific junctions;
3. split-half on the reference animal itself.

Arms 2 and 3 are added as amendments before they run.

## Arm 1 — anchor trim (gorilla = dev; human = held-out)

**Truth construction** (`bench/o2_anchor/anchor_trim.py`, 7 tests; run before this prereg, counts only):
- **Read selection.** Per catalog family with >= 2 copies (gorilla_KB3781 default catalog, 326 families / 918 copies),
  take primary MAPQ-60 reads of the fibroblast BAM (`-F 2308 -q 60`). The reads come from the assembly's own animal.
  Each read's reference span must touch exactly one family; truth copy = the copy it overlaps most. At most 20 reads per
  copy, the lowest md5 of the name.
- **Trimming.** Each read, in its original orientation, is aligned to the genomic spans of its family's copies
  (`minimap2 -c -x splice:hq -uf`, secondaries kept). It is trimmed to the interval covered by both the truth copy's best
  hit and the best sibling's best hit. The trimmed read is kept iff that interval is >= 300 bp and >= 100 bp were
  removed. The truth is therefore the removed, copy-unique part of the read. O2 sees only the part shared with a sibling.
- **Feasibility (gorilla).** 10,395 reads picked, 4,048 trimmed, from 240 families.
- **Re-alignment.** The trimmed reads are re-aligned genome-wide with the original BAM's settings (`-ax splice:hq -uf
  --eqx -Y -N 50 -p 0.1 --secondary=yes`, GGO.splice.mmi).
- **Primary MAPQ after trimming.** 3,152 stay MAPQ 60 (the shared part still places them), 744 are MAPQ 1-59 and 152
  are MAPQ 0.

**O2 run.** `copy_assign` exactly as the pipeline's assign stage:
- inputs: `--regions`, `--families`, `--copies-fa`, plus `--skip-poa-diagnostic --region-threads 4 --tied-prefilter`
  and `--dump-star`;
- BAM: the original records in every catalog copy's ±50 kb window (what `copy_assign` reads) merged with all records
  of the trimmed reads (names `<read>|trim`);
- families: `--only-families` = every family a trimmed read's record touches, run in batches if needed (rows are per
  read x family, so batching does not change them).

The original full reads stay in the BAM; each trimmed read is an extra molecule.

**Scoring** (`bench/o2_anchor/score.py`; locus-level: a copy is right iff its catalog span overlaps the truth copy's
span):
- **T0 truth gate.** Among trimmed reads still MAPQ 60, the aligner's primary is on the truth copy for >= 0.99. If not,
  the truth construction is suspect: stop, nothing below is read.
- **A1 precision.** Over the contested trimmed reads (`contested` = 1 in the truth family's row), count every `assigned`
  row, in the truth family and in any other family (a foreign assigned row is wrong).
  - PASS iff wrong / assigned <= 0.01 with >= 30 assigned. UNDERPOWERED if < 30 assigned.
  - The 95% upper bound is reported.
- **A2 calibration.** Some contested trimmed reads are NOT distinguishable: the truth copy is not a candidate, or some
  rival differs from it at < 2 read-star columns the read covers (YAG Arm A's definition, O2's m >= 2 step). Of
  those, <= 0.05 assigned.
- **Reported, no bar:**
  - recall over distinguishable contested reads;
  - the aligner primary's accuracy on the same contested reads, and on those O2 assigns (O2 against a coin toss);
  - per family;
  - an examples file: every wrong call, and up to 20 correct calls with their columns.

**Held-out: human A119b** (testis, catalog human_A119b, CHM13; the donor is not the reference man). Same code and bars
(H0-H2), run only after the gorilla results are written below. A donor mismatch can lower recall there. Its effect on
precision is the question.

Outputs `/mnt/linuxdisk/tmp/o2_anchor/`; the merged BAM is on C: (`/mnt/c/Users/jfris/Desktop/data/o2_anchor/`).

## Results

(pending)
