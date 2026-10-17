# Real fusion genes as their own category (Yoo 2025), 2026-10-17

**User (2026-10-17):** "one of our papers warn about fusion loci that are real, lets ensure to keep those in a special
category to not confuse them with spurious readthroughs or discard them if they are still not annotated"; then "I think it
was Yoo, 2025, lets implement that".

**The paper.** Yoo et al. 2025 (complete ape genomes, Nature; vault `Papers/Yoo_2025_Complete_Ape_Genomes.md`) count fusion
genes of segmental-duplication expansions as gene copies:
- gorilla MAPKBP1 + JMJD7-PLA2G4B + SPTBN5, 8 units on chr1;
- chimpanzee/bonobo PLA2G4B-JMJD7, 6/3 copies;
- orangutan HERC2-GOLGA, 5/11 copies.

A copy is real when Iso-Seq supports it with more than 3 reads and it carries a valid ORF (TransDecoder). Human RefSeq
annotates JMJD7-PLA2G4B as a readthrough, and gorilla RefSeq annotates no fusion at all, only the two genes (LOC129529234
JMJD7 and LOC129529226 PLA2G4B at chr1 15.67 Mb). An annotation therefore cannot decide.

## What was built

Branch `fusion-category` (Rust) and `o3-fast` (O3 label). Defaults are unchanged except where noted.

1. **`<out>.fusions.tsv`** (with `--bridge-regroup f1|f1v2`): every F1v2 bridge transcript gets a class from read and sequence
   evidence alone.
   - **`fusion_gene`:** at least 4 primary reads on every bridge junction, each junction canonical (GT-AG, GC-AG or AT-AC),
     and the longest complete ORF (at least 100 aa, TransDecoder's default minimum) is a fusion protein across every bridge
     junction.
     - "A fusion protein" means each partner gives at least one whole codon and is read in its own dominant open frame.
     - Before the junction, the ORF's open run is a longest stop-free run of that part; after it, its run to its stop is a
       longest stop-free run of the rest.
     - There is no length threshold. An ORF that only runs on a few junk codons past the junction, out of the partner's frame
       (review I1), is a readthrough.
   - **`readthrough`:** supported, but no fusion protein (no fusion ORF, or one that leaves a partner's frame).
   - **`unsupported`:** fewer reads, or a non-canonical junction.
   - **`contested`** (a column, not a class): some bridge junction is carried by more secondary than primary alignments, so
     which copy the reads came from is O2's call.
   - The GTF is unchanged. Bridges stay their own `<g>.fus<k>` loci tagged `fusion_of`, outside the families input, as F1v2
     already did.
2. **`RUSTLE_READTHROUGH_PROTECT`** (on by default whenever `RUSTLE_READTHROUGH_JUNCTIONS` runs a statistic arm, `r`, `rq1`,
   `r2` or `r3`).
   - A junction the rule flags keeps its reads when it passes the fusion_gene test: at least 4 primary reads, canonical, and
     its modal molecule codes a fusion protein across it.
     - The modal molecule is the intron chain most of its reads share, as that chain's longest read. One aberrant read no
       longer decides (review I5).
     - The test is the same ORF test, run on the genome along that read's exon blocks.
   - Kept junctions are listed in `<out>.readthrough_protected.tsv`, with the tier columns of r2/r3.
   - The ALL dump gains a last `protected` column. The tier counts and `params.tsv` count only removed junctions (review I3).
   - `off` reproduces the registered arms byte for byte. The variable is read only when the filter is on.
   - `readthrough_eval.py null --exclude <readthrough_protected.tsv>` keeps the kept junctions out of the NULL draw.
   - A `list:` arm removes what it lists. A read that also carries another, unprotected flagged junction still leaves.
3. **O3 `detect.py --fusions <out>.fusions.tsv`:** `calls.tsv` gains a `fusion` column (fusion_gene or readthrough) when at
   least half of a call's members' primary alignments splice that class's junction.
   - This is a label, not a filter. A missing copy of a fusion gene, such as a ninth JMJD7-PLA2G4B unit in another gorilla,
     has reads that splice its sister copies' junctions. Removing them would blind O3 to exactly Yoo's case.

## Validation

**Yoo's own example, gorilla KB3781 chr1 15-30 Mb (fibroblast FLNC), the positive control:**
- **All 8 JMJD7-PLA2G4B copies are assembled as whole loci, not bridges.**
  - Each locus holds JMJD7-only, PLA2G4B-only and spanning transcripts of 24-26 exons.
  - F1v2 finds one bridge in the whole region: LYST -> LOC129529029, classed `readthrough` (77-aa ORF, not across).
- **The ORF test sees the fusion protein, in frame on both sides.**
  - `fusion_gene`: copies 1 and 2 code 1,012 / 965 / 913-aa ORFs across the JMJD7 -> PLA2G4B junction. 1,012 aa is the
    human JMJD7-PLA2G4B fusion protein's length. One transcript of the 19.02 Mb copy codes 706 aa.
  - `readthrough`: the other 19.0-19.3 Mb transcripts end their ORF in JMJD7 (309-319 aa), and at 29.08 Mb the ORF is not
    across in frame (543-639 aa). This matches Yoo's "2 with full-length valid start/stop".
  - Every copy is `contested`.
- **⚠ The first rule failed this control.**
  - That rule required primary >= secondary alignments (the §6n9 test).
  - Every copy's junction is carried by more secondary than primary alignments: 14 / 28 at copy 1, 24 / 429 at 19.02 Mb,
    45 / 203 at 29.08 Mb. Each read also aligns to the 7 sister copies.
  - The old rule would have called all of Yoo's real copies unsupported. The secondary majority is now the `contested` column
    (the copy of origin is O2's question, the user's 10-11 rule), never the class.
- **r3 flags none of the 8 JMJD7 -> PLA2G4B junctions.**
  - Its 8 flags are a copy 2 -> copy 3 tandem readthrough (110 kb), junctions inside RAB4A, B3GALNT2 and NID1, and the
    LYST bridge.
  - None is a fusion gene, so the protection keeps nothing here.

**Human A119b chr16 (testis):**
- **38 bridges:** 3 fusion_gene, 31 readthrough, 4 unsupported; 3 contested (all of them unsupported).
- **The bridges that pass or fail the fusion test:**

| bridge | junction | primary reads | ORF across | annotation |
|---|---|---|---|---|
| CREBBP internal | 3,798,252-3,801,016 (-) | 199 | 779 aa | one gene: F1v2 split CREBBP and the bridge is its own coding mRNA |
| MSTRG.2130 -> DCTN5 | 23,842,869-23,918,762 (+) | 12 | 237 aa | no fusion annotated |
| PALB2 -> NDUFAB1 | 23,863,133-23,879,395 (-) | 8 | 852 aa, not in frame: **readthrough** | no fusion annotated (`fusion_gene` before the in-frame test) |
| FBXO31 -> C16orf95 | 93,386,708-93,404,772 (-) | 10 | 573 aa | CAT AC010531.1, an annotated readthrough protein |

- **Under r3 without protection,** the reads of MSTRG.2130-DCTN5 and FBXO31-C16orf95 are removed.
- **With protection, 36 of r3's 1,024 flagged junctions keep their reads,** both of those included.
  - The first version, with the longest read and no frame test, kept 123. Most of them were junctions inside one coding gene,
    e.g. RBL2.
  - r3 now removes 3,394 alignments instead of 3,751 (-9.5%). Tier A 783 + B 205 = 988 removed.
- **`RUSTLE_READTHROUGH_PROTECT=off` is byte-identical to the registered r3 arm** on chr16. Its readthrough_junctions TSV, ALL
  dump, GTF, families GTF and `params.tsv` are all identical.
- **The families input (review I4):** 9 of the 36 kept junctions sit in families-input loci that r3 without protection
  splits into 2-10 pieces. Without the filter, which is the default, they are one locus anyway. Whether any of them joins
  two families is unmeasured.

**The known real fused transcripts PKD1P6-NPIPP1 and PDXDC2P-NPIPB14P** are CAT lncRNAs and not bridges. F1v2 keeps them as
loci, so the class never meets them. Under the ORF test they would be `readthrough`: real molecules with no fusion protein.

## What is not claimed

- **The class only judges F1v2 bridges (review I6).** A fusion gene that is its own locus, like all 8 gorilla JMJD7-PLA2G4B
  copies, is kept but not labelled, so the O3 label cannot fire for it either.
  - No read signal marks it as two genes. V1 (reads starting at PLA2G4B's own promoter inside the intron) is 0 at every
    copy.
  - A census of junctions joining two independently transcribed units (U >= 1 and V1 >= 1) would miss all 8 copies. It
    would also flag 2,548 junctions on chr16.
  - Labelling whole loci needs the annotation's gene boundaries. `bench/fusion/yoo_jmjd7_control.py` does that for JMJD7;
    a general guided labeller is the open step.
- **The ORF test is annotation-free and lenient** (TransDecoder's 100 aa; Yoo also required a blastp hit covering >= 75% of a
  known protein). The r3 protection therefore also keeps coding junctions inside one gene.
- **Its effect on r3's pre-registered metrics** (r1117-r1120) is unmeasured; `RUSTLE_READTHROUGH_PROTECT=off` is the
  registered arm.
- **The O3 label is exercised by unit tests only.** The simulated worlds carry no fusions.

## Files

- Rust:
  - `src/family.rs`, module `bridge_regroup`: `fusion_class`, `fusion_contested`, `fusion_frame`, `orf_across`,
    `classify_bridge`, `classify_bridges` (one contig's sequence at a time, review I2).
  - `src/family.rs`, module `denovo_assemble`: `ReadthroughStats::with_protect`, `rt_modal_molecule`, `rt_fusion_orf`,
    `ReadthroughSwitch::from_values`, `ReadthroughFlags::write_protected_tsv`.
  - `src/genome.rs`: `contig_seq`.
  - `src/copy_assign_bin.rs`: outputs and params rows.
- Python:
  - `bench/o3_graph/detect.py`: `read_fusions` / `fusion_carriers` / `fusion_label` (`o3-fast`, 4b504017, ab0b5f18).
  - `bench/mechanism/readthrough_eval.py`: `null_candidates` and `null --exclude`.
- Runs:
  - `/mnt/linuxdisk/tmp/fusion_ggo/` (gorilla; positive control `bench/fusion/yoo_jmjd7_control.py GTF GFF BAM FASTA`).
  - `/mnt/linuxdisk/tmp/fusion_chr16/` (human).

## Review (fresh reviewer, whole branch, 2026-10-17)

There were no critical findings, 6 Important and 8 Minor.
- **Fixed with tests:**
  - I1: the in-frame test.
  - I2: one contig at a time.
  - I3: tier counts only for removed junctions, the `protected` column, and `null --exclude`.
  - I5: the modal molecule.
  - M1: contested per junction.
  - M2: a fusion_gene junction is never demoted.
  - M3: the variable is read only with the filter on, and `on` is accepted.
  - M5: the chains are kept only when the protection can act.
  - M6: `--fusions` is read at startup, with region fetch.
  - M7: help text, wording and tier columns.
  - M8: a minus-strand ORF test.
- **I4:** measured, see Validation.
- **I6:** not fixed, see What is not claimed.
- **M4:** documented.
