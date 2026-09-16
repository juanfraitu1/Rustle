# Chr20 assembler comparison: ours vs StringTie vs FLAIR (gffcompare + SQANTI3)

Date: 2026-09-15. Purpose: a standard, defensible comparison of our multi-copy-family tool's de novo
transcript-assembly output against StringTie and FLAIR, on an **ordinary** human chromosome that was not
chosen or tuned for multi-copy-family detection — the advisor is skeptical that this tool is a "decent
assembler" in general, so this benchmark asks the question directly, at defaults, with no thumb on the
scale for our tool.

## Substrate

- **BAM**: `/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam` (T2T-CHM13-aligned
  human testis IsoSeq, genome-wide), restricted to **chr20** (66,210,255 bp).
  `samtools idxstats` on chr20: 188,864 mapped records; `samtools flagstat` on the chr20-only BAM: 25,341
  primary / 163,244 secondary / 279 supplementary (100% mapped, unpaired long reads). chr20 was chosen
  because it is moderate-sized and, per this project's own working ledger, only lightly touched elsewhere
  — genuinely not cherry-picked for multi-copy-family work.
- **Genome FASTA**: `/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0.fa`, chr20 extracted
  (66,210,255 bp, matches the BAM's `@SQ` line exactly) and re-indexed.
- **Reference annotation**: `/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz`,
  chr20 subset (99,049 GFF3 rows), converted to GTF (4,574 transcripts / 43,173 exons / 32,511 CDS rows).
  gffcompare's own reference parse reports 4,563 reference mRNAs in 1,140 loci (4,286 multi-exon) — the
  11-transcript gap from gffread's 4,574 is gffcompare's own standard collapsing of redundant/duplicate
  reference transcripts (identical intron chain + contained boundaries) during its internal annotation
  build, not a data-prep error.
- **Tools**: `copy_assign` (this repo, release build, commit-current), StringTie
  (`tools/stringtie/stringtie`, long-read mode `-L`), FLAIR 3.0.0 (conda env `flair`). SQANTI3 (conda env
  `sqanti3`, checkout at `/mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3/`).

All working data lives under `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/` (not in git — large
intermediate/output files). The scripts below are committed to the repo under `bench/`.

## A genuine reference-conversion bug found and fixed (worth flagging)

The RefSeq GFF3 lists a feature's leftmost **exon** before its own **gene/transcript** record whenever
they share the same start coordinate (confirmed present in the un-filtered genome-wide
`chm13v2.0_RefSeq_full.gff.gz`, not introduced by our chr20 filtering). `gffread`'s single-pass GFF3
parser needs to see the parent before the child to propagate the true gene ID into GTF's `gene_id`; left
in original order, **every** converted transcript's `gene_id` silently collapsed to its own
`transcript_id` (verified: the same two lines convert correctly in a 2-line isolated extract, but wrong
inside the full 99,049-row file). Fixed with a stable resort (chrom, start, feature-rank
[gene/pseudogene=0, mRNA/transcript/other=1, exon/CDS/UTR=2], original order) before conversion — see
`bench/prep_chr20_ref.sh` step 4. This doesn't change gffcompare's own stats (which cluster reference
transcripts by genomic overlap, not by `gene_id`), but it would have broken any per-gene reporting.

## A genuine FLAIR 3.0.0 packaging bug found and worked around

1. **`flair correct` requires an annotation.** The installed FLAIR 3.0.0 (`flair_brookslab-3.0.0`) hard-
   requires `-f/--gtf`, `--junction_tab`, or `--junction_bed` for its `correct` subcommand
   (`FlairInputDataError: No junctions from GTF or junctionsBed to correct with. Exiting...`); fully
   unguided splice correction, which the project's existing `bakeoff_flair.sh` (gorilla substrate) relies
   on, no longer exists in this version. Supplying our chr20 reference GTF here would make FLAIR's arm
   *guided* while ours/StringTie stay unguided, so instead `bench/bakeoff_chr20_flair.sh` skips `correct`
   and feeds the raw `flair align` bed straight into `flair collapse` — unguided at the cost of skipping
   the splice-site-correction step. Documented as a deviation, not hidden.
2. **`filter_transcriptome_align.py` fails under the wrong Python.** `flair collapse` internally calls
   `flair.set_unix_path()`, which prepends the flair *package directory* itself
   (`site-packages/flair/`, where `filter_transcriptome_align.py` physically lives) to `PATH`. That
   script's shebang is `#!/usr/bin/env python3`; since the package dir has no `python3` binary of its own,
   `env` searches the rest of `PATH` — and on this machine `/home/linuxbrew/.linuxbrew/bin` (no `flair`
   package) sits ahead of the flair conda env's own `bin/` in the inherited shell `PATH`, so the script
   executed under the wrong Python and died with `ModuleNotFoundError: No module named 'flair'` even
   though flair is fully installed. Fixed non-invasively in `bench/bakeoff_chr20_flair.sh` by exporting
   `/home/juanfra/miniforge3/envs/flair/bin` to the front of `PATH` before invoking any flair subcommand
   (still ahead of linuxbrew after flair's own self-prepend). No conda environment files were modified.

## Commands run

```bash
# 1. substrate prep (chr20 BAM/FASTA/reference GTF)
bash bench/prep_chr20_ref.sh

# 2. build check
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign

# 3. the three assemblers, one at a time, foreground
bash bench/bakeoff_chr20_ours.sh
bash bench/bakeoff_chr20_stringtie.sh
bash bench/bakeoff_chr20_flair.sh

# 4. score: gffcompare (x3 + multicopy=false subset) and SQANTI3 QC (x3)
bash bench/chr20_score.sh
```

Exact per-tool invocations (as run by the scripts above):

- **Ours**: `copy_assign --gtf --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out ours`
  (no `--families` — pure de novo detection, no catalog).
- **StringTie**: `stringtie -L -p 4 -o st.gtf chr20.bam`.
- **FLAIR**: `samtools fastq -F 2308 chr20.bam > reads.fq` → `flair align -g chr20.fa -r reads.fq -o flair
  --threads 4` → (correct skipped, see above) → `flair collapse -g chr20.fa -q flair.bed -r reads.fq -o
  flair --threads 4 --generate_map`.
- **gffcompare**: `gffcompare -r chr20_ref.gtf -o <label> <tool.gtf>`, once per tool plus once for the
  `multicopy=false`-only subset of `ours.gtf`.
- **SQANTI3**: `sqanti3_qc.py --isoforms <tool.gtf> --refGTF chr20_ref.gtf --refFasta chr20.fa -o <label>
  -d <dir> --report skip -t 4` (`--report skip`: only the `_classification.txt` structural-category table
  was needed, so the R/HTML report rendering path was skipped; nothing else was skipped — ORF prediction
  ran at defaults).

## gffcompare results (whole chr20)

| Tool | Query mRNAs (loci) | Base Sn/Pr | Exon Sn/Pr | Intron Sn/Pr | Intron-chain Sn/Pr | Transcript Sn/Pr | Locus Sn/Pr |
|---|---|---|---|---|---|---|---|
| **Ours** | 976 (456 loci) | 11.6 / 69.3 | 18.0 / 71.0 | 18.6 / 86.1 | 8.0 / 44.6 | 7.6 / 35.6 | 18.2 / 46.1 |
| **StringTie** | 712 (359 loci) | 11.9 / 81.4 | 19.6 / 79.6 | 20.8 / 86.9 | 7.7 / 47.4 | 7.3 / 47.1 | 18.9 / 60.7 |
| **FLAIR** | 820 (335 loci) | 10.4 / 79.8 | 18.3 / 68.1 | 19.4 / 74.2 | 6.2 / 35.1 | 5.8 / 32.3 | 14.3 / 49.0 |

Reference: 4,563 mRNAs in 1,140 loci (4,286 multi-exon), same for all three comparisons.

Full stats blocks (verbatim gffcompare output):

### Ours

```
#     Query mRNAs :     976 in     456 loci  (773 multi-exon transcripts)
#            (170 multi-transcript loci, ~2.1 transcripts per locus)
# Reference mRNAs :    4563 in    1140 loci  (4286 multi-exon)
# Super-loci w/ reference transcripts:      314
#-----------------| Sensitivity | Precision  |
        Base level:    11.6     |    69.3    |
        Exon level:    18.0     |    71.0    |
      Intron level:    18.6     |    86.1    |
Intron chain level:     8.0     |    44.6    |
  Transcript level:     7.6     |    35.6    |
       Locus level:    18.2     |    46.1    |

     Matching intron chains:     345
       Matching transcripts:     347
              Matching loci:     208

          Missed exons:    7353/9619	( 76.4%)
           Novel exons:     256/2806	(  9.1%)
        Missed introns:    6089/8563	( 71.1%)
         Novel introns:      36/1849	(  1.9%)
           Missed loci:     821/1140	( 72.0%)
            Novel loci:     123/456	( 27.0%)

 Total union super-loci across all input datasets: 437
976 out of 976 consensus transcripts written (0 discarded as redundant)
```

### StringTie

```
#     Query mRNAs :     712 in     359 loci  (698 multi-exon transcripts)
#            (151 multi-transcript loci, ~2.0 transcripts per locus)
# Reference mRNAs :    4563 in    1140 loci  (4286 multi-exon)
# Super-loci w/ reference transcripts:      325
#-----------------| Sensitivity | Precision  |
        Base level:    11.9     |    81.4    |
        Exon level:    19.6     |    79.6    |
      Intron level:    20.8     |    86.9    |
Intron chain level:     7.7     |    47.4    |
  Transcript level:     7.3     |    47.1    |
       Locus level:    18.9     |    60.7    |

     Matching intron chains:     331
       Matching transcripts:     335
              Matching loci:     215

          Missed exons:    7086/9619	( 73.7%)
           Novel exons:     171/2320	(  7.4%)
        Missed introns:    5802/8563	( 67.8%)
         Novel introns:      51/2050	(  2.5%)
           Missed loci:     812/1140	( 71.2%)
            Novel loci:      19/359	(  5.3%)

 Total union super-loci across all input datasets: 344
712 out of 712 consensus transcripts written (0 discarded as redundant)
```

### FLAIR

```
#     Query mRNAs :     820 in     335 loci  (753 multi-exon transcripts)
#            (147 multi-transcript loci, ~2.4 transcripts per locus)
# Reference mRNAs :    4563 in    1140 loci  (4286 multi-exon)
# Super-loci w/ reference transcripts:      261
#-----------------| Sensitivity | Precision  |
        Base level:    10.4     |    79.8    |
        Exon level:    18.3     |    68.1    |
      Intron level:    19.4     |    74.2    |
Intron chain level:     6.2     |    35.1    |
  Transcript level:     5.8     |    32.3    |
       Locus level:    14.3     |    49.0    |

     Matching intron chains:     264
       Matching transcripts:     265
              Matching loci:     163

          Missed exons:    7335/9619	( 76.3%)
           Novel exons:     454/2907	( 15.6%)
        Missed introns:    5995/8563	( 70.0%)
         Novel introns:     292/2241	( 13.0%)
           Missed loci:     864/1140	( 75.8%)
            Novel loci:      61/335	( 18.2%)

 Total union super-loci across all input datasets: 322
820 out of 820 consensus transcripts written (0 discarded as redundant)
```

## `multicopy` stratification (ours only)

`copy_assign --gtf` tags every emitted transcript `multicopy "true"`/`"false"` (`src/bin/copy_assign.rs:
193-198, 3575-3588`) based on **its own de novo family/copy detector**, run with no external catalog. On
chr20, that detector found **0 co-located families** (`ours/ours.stderr.log`: `"chr20:1-66210255: 91
mapped reads -> 0 families"`, `"refine: 0 co-located families -> 0 homology-gated"`). Consequently **all
976/976 emitted transcripts are `multicopy "false"`** — chr20 has essentially no multi-copy families in
the detector's de novo sense, which is exactly the "ordinary chromosome" property this benchmark wanted.
The `multicopy="false"`-only subset is therefore **byte-identical** to the whole-chromosome `ours` set;
we materialized it explicitly (`bench/chr20_score.sh` step 2) and reran gffcompare on it rather than just
asserting the identity — the resulting `.stats` block is character-for-character the same as the whole-set
table above. This is an honest degenerate result, not a shortcut: on this substrate, there is no
"hard-locus" subset to separate out, so the whole-chromosome numbers already ARE the "ordinary assembler"
numbers for our tool.

## SQANTI3 structural category classification

SQANTI3 ran successfully end-to-end for all three tools (`--isoforms <gtf> --refGTF chr20_ref.gtf
--refFasta chr20.fa --report skip -t 4`), no dependency issues encountered (cDNA_Cupcake is vendored in
the SQANTI3 checkout under `src/utilities/cupcake`, not an external pip dependency, and R/Rscript were
present and not needed once `--report skip` bypassed the HTML/PDF rendering step). StringTie's arm logged
13 transcripts discarded for unknown strand (`.` strand — StringTie sometimes emits this for single-exon
transcripts with no splice evidence to call strand from); 699/712 were classified.

| Category | Ours (n=976) | StringTie (n=699) | FLAIR (n=820) |
|---|---|---|---|
| full-splice_match (FSM) | 352 | 332 | 271 |
| incomplete-splice_match (ISM) | 217 | 79 | 58 |
| novel_in_catalog (NIC) | 77 | 84 | 62 |
| novel_not_in_catalog (NNC) | 166 | 177 | 335 |
| genic | 31 | 1 | 11 |
| genic_intron | 16 | 1 | 8 |
| antisense | 80 | 12 | 47 |
| intergenic | 31 | 9 | 16 |
| fusion | 6 | 4 | 12 |

(`ours` and `flair` classify 100% of their emitted transcripts; StringTie classifies 699/712, 13 dropped
for unknown strand — see above. The `multicopy=false` subset for `ours` is not shown as a separate SQANTI3
row: it is the identical 976-transcript set, per the stratification section above.)

## Honest summary

On chr20 — an ordinary chromosome where our tool's own de novo family detector finds nothing to grab onto
(0 multi-copy families, so every transcript is emitted through the same plain FLAIR-style intron-chain-
collapse-plus-gate assembly path FLAIR/StringTie also use) — all three tools land in a similar, modest
regime versus the RefSeq annotation, and StringTie is the strongest of the three on most measures.
StringTie has higher precision than ours at every single level gffcompare reports (35.6–86.1% for ours vs
47.1–86.9% for StringTie, level by level) and higher sensitivity at base/exon/intron/locus level too (e.g.
exon 18.0% vs 19.6%, intron 18.6% vs 20.8%), plus by far the fewest novel loci (5.3%, vs 27.0% for ours
and 18.2% for FLAIR) — consistent with its long-standing engineering for clean, conservative long-read
transcript models. Ours has a narrow sensitivity edge over StringTie specifically at intron-chain and transcript
level (8.0%/7.6% vs StringTie's 7.7%/7.3%; FLAIR trails both at 6.2%/5.8%), the only levels where ours
leads any tool, and it does so while giving up precision at those same two levels (44.6%/35.6% vs
StringTie's 47.4%/47.1%) and emitting the most novel loci of the three (27.0%) — a more permissive,
lower-precision profile than StringTie, closer to FLAIR's permissiveness (FLAIR has the most NNC
transcripts and novel introns of the three, and the lowest intron-chain/transcript-level precision).
None of the three is a strong instrument on this ordinary chromosome in absolute terms — all sit under
21% sensitivity at every level and under 87% precision at any level — which itself is informative: this
is what an "ordinary" long-read assembly problem looks like at IsoSeq depth on one chromosome without
short-read or annotation guidance, and our tool's numbers here are in the same ballpark as two
purpose-built, widely used long-read assemblers, not degraded relative to them, though StringTie is
plainly the better plain assembler of the two "general" alternatives on this substrate. The main honest
takeaway for the advisor: ours trades StringTie's precision (and most of its sensitivity) for a thin edge
in intron-chain/transcript recall and a larger novel-locus count — the expected shape for a tool whose
primary design target is copy-level recall inside multi-copy families, tested here with that machinery
idle.

## Follow-up: does TSS/TES boundary snapping help on THIS substrate? (2026-09-15, negative)

`RUSTLE_TSS_SNAP` (start/end quantile → sharp-peak snapping) and `RUSTLE_TES_EXTEND` (3'-peak sequence
extension) are real, opt-in, off-by-default mechanisms in `denovo_assemble.rs`/`denovo_pipeline.rs`. Both
were previously validated only on the Soto multi-copy-family benchmark (chr1/7/15/16), where they were
found to move almost nothing (2/43 copy boundaries, paired p=0.69) and kept opt-in "for absence of benefit
rather than demonstrated harm." That verdict had never been tested against a general/ordinary-locus,
gffcompare-style target — chr20 is exactly that test.

**Command**: `RUSTLE_TSS_SNAP=1 RUSTLE_TES_EXTEND=1` + the same `bakeoff_chr20_ours.sh` invocation, output to
`ours_tessnap/`, scored against the same `chr20_ref.gtf`.

**Real effect on the GTF**: 48 of 976 transcripts got a different start/end coordinate (confirmed by diffing
transcript lines) — the flags are genuinely reachable and active on this substrate, not a no-op.

**Real effect on gffcompare**: none, at the level gffcompare measures.

| | Matching intron chains | Matching transcripts | Matching loci | Base Sn/Pr | Transcript Sn/Pr |
|---|---|---|---|---|---|
| Baseline (flags off) | 345 | 347 | 208 | 11.6 / 69.3 | 7.6 / 35.6 |
| TSS_SNAP + TES_EXTEND | 345 | 347 | 208 | 11.3 / 69.0 | 7.6 / 35.6 |

Every transcript/intron-chain/locus match count is IDENTICAL; base/exon-level sensitivity and precision move
by ≤0.3 percentage points, in the negative direction. **Conclusion: the prior "no benefit" verdict
generalizes from the multi-copy Soto substrate to this ordinary chromosome — freshly re-derived here, not
assumed.** Not worth enabling for general assembly either. Left off by default; no code change.

## Files

- `bench/prep_chr20_ref.sh` — chr20 BAM/FASTA/reference-GTF extraction (incl. the GFF3 resort fix).
- `bench/bakeoff_chr20_ours.sh` — our tool, pure de novo (`--gtf`, no `--families`).
- `bench/bakeoff_chr20_stringtie.sh` — StringTie, long-read mode, no annotation.
- `bench/bakeoff_chr20_flair.sh` — FLAIR, unguided (incl. the FLAIR 3.0.0 workarounds above).
- `bench/chr20_score.sh` — gffcompare (x3 + multicopy=false subset) + SQANTI3 QC (x3) driver.
- Large outputs (BAMs, GTFs, SQANTI3 intermediates) live under
  `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/` and are NOT in git.
