# Pre-registration: family layer (FW), singleton confirmation across copies (FAMILY_LAYER_SPEC Part B)

**Written 2026-10-26, before any FW run or FW number.** Branch `assembler-stringtie-gaps`. A default flip is the user's call.

## Seen before this file
Part A stratified bakeoff (human A119b only): multi-copy recall / precision 20.7 / 12.6 (shipped pooled run), single-copy
22.1 / 19.8. The ledger (`docs/LEDGER_GENERAL.md` "How StringTie chooses single-read chains") measured genome-blind
single-read admission (`--min-chain-reads 1`, `--read-isoform-k 2`): more chains, lower precision (register 848 / 855).
No family-restricted number has been computed.

## Rule (annotation-free; `bench/stgaps/fw_layer.py`)
Inputs: BASE = the shipped `assemble` GTF; MIN1 = the same command with `RUSTLE_MIN_CHAIN_READS=1`; the sample's O1
families (`<sample>.fam.copies.tsv`, copy loci `locus_start..locus_end`); the BAM; the genome.
1. **Copy chains.** Per copy locus, primary spliced reads (`-F 0x904`) overlapping it, keyed by exact intron chain;
   molecules de-duplicated on (chain, start, end).
2. **Copy map.** Per family with >= 2 copies, the copy with most reads (`n_reads`) is the representative; every other
   copy locus is aligned to it (`minimap2 -x asm20 -c --cs`, best hit to its own representative). A junction lifts when
   its last exon base and next exon base both fall on aligned (match or mismatch) columns; its family junction is the
   sorted pair of lifted positions. The representative's junctions lift to themselves.
3. **Candidates.** MIN1 multi-exon transcripts whose intron chain is absent from BASE (same contig) and that overlap a
   family copy locus. Own support = molecules at that copy with exactly the candidate's chain.
4. **Admission (must-observe).** A candidate with own support == 1 is **confirmed** when every junction lifts and the
   family chain's molecules summed over all copies of the family are >= 2 (the pass-1 floor). Confirmed candidates are
   appended to BASE unchanged (MIN1 coordinates). Nothing is emitted where no read carries the whole chain; no structure
   is transferred between copies.
FW = BASE + confirmed. BASE lines are never altered or removed.

## Substrates
Dev (reported, not judged): human A119b genome-wide (`A119b.t2t.bam`, CHM13 v2.0, `fam.copies.tsv` of
`runs/human_A119b`). **Held-out: gorilla KB3781 fibroblast genome-wide** (`fibroblasts/GCA_029281585.2_flnc_mm.bam`,
`runs/gorilla_KB3781` families), strata from `fw/ggo_reftx.tsv` / `ggo_reftx.paf` / `ggo_ref.gtf` (GGO RefSeq
GCF_029281585.2, Amendment-1 definition, built before this file, never scored).

## Clauses (each substrate; `stratify.py` strata; BASE and MIN1 run by the same binary)
- **F1 gain:** multi-copy recall FW − BASE >= 0.5 points.
- **F2 precision:** multi-copy precision FW >= BASE − 0.5 points.
- **F3 the confirmation carries information:** among candidates with own support == 1, the gffcompare `=` share of
  confirmed >= 1.5 × the `=` share of unconfirmed.
- **F4 single-copy untouched:** every BASE line is in FW, and single-copy recall and precision |FW − BASE| <= 0.1 point.
**ADOPT as opt-in** = F1-F4 on held-out gorilla KB3781. Otherwise a register row.
Descriptive: confirmed candidates with no reference gene (possible unannotated copies, an O3 caution, not gain);
family size of confirmed candidates.
