# Family-wise assembly as a layer on the single-copy assembler — spec (2026-10-26)

User (2026-10-26): *"a layer on top of the stringtie and flair like assembly so we can tell that the assembler works for
single copy mode and then extrapolate to multi copy."*

## Part A — single-copy mode, stratified bakeoff (descriptive; definitions frozen here before any stratified number)

**Copy status of a reference gene is annotation-only (no reads, no assembler):** every RefSeq transcript's spliced
sequence is mapped to its own genome with `minimap2 -x splice -N 50 --secondary=yes` (human: the prebuilt
`chm13v2.0.splice.mmi`). A hit counts when it covers >= 80% of the transcript at >= 95% identity (matches / aligned
columns) and does NOT overlap the transcript's own locus. A gene is
- **multi-copy** if any of its transcripts has >= 1 counted hit;
- **single-copy** otherwise.
Arms: ours (shipped `--assemble-only` polish, the 09-22 `benchmark_collapse/assembler_*` GTFs), StringTie, FLAIR,
IsoSeq (the lab's runs). Metric per stratum: reference multi-exon transcripts recovered (gffcompare `=`, ref side) and
query precision (`=` / multi-exon queries assigned to that stratum's genes or to no gene).
Reading: the single-copy stratum shows the assembler works as a StringTie/FLAIR-class tool; the multi-copy gap (ours
single-copy recall minus ours multi-copy recall, and vs the tools) is what the layer targets.

## Part B — the family layer (FW), on top of the unchanged assembler

Input: the base assembler's per-copy chains (exact intron chains of primary reads, the pass-1 groups), O1 families
(`<sample>.fam.copies.tsv`), the genome. Steps:
1. **Copy map.** Each family's copy loci are aligned to the family's best-supported copy (minimap2, genomic); a junction
   is lifted when both its donor and acceptor lift exactly. Junctions that lift to the same pair form one family junction.
2. **Family chain support.** A chain at copy c, rewritten in family junctions, is supported by the reads carrying that
   family chain at ANY copy.
3. **Admission (must-observe).** A chain is admitted at copy c when it has >= 1 own read at c (whole chain observed at c)
   AND family support >= the pass-1 floor (2). Nothing is ever emitted at a copy without its own read for the whole chain:
   no structure is transferred between copies (r647, r1285).
4. Admitted chains enter the unchanged polish; output tags `family_support`, `family_copies` (copies where the chain is
   read-observed).
**Single-copy invariance (by construction, checked):** a one-copy family has family support == own support, so FW's
output equals the base assembler's byte for byte on every single-copy locus.

Tests (pre-registered separately before any FW run): single-copy byte-identity; multi-copy stratum recall gain vs
precision; NPIP 26-copy complete chains / junctions (current 20/26, 197/209; ceiling 22/26, 205/209); held-out
gorilla fibroblast families. Constraints: O1-only (no quantification), no facility location / bipartite matching, no
human+gorilla pooling.

**Amendment 1 (2026-10-26, after one table was seen).** "Identity = matches / aligned columns" used PAF col10/col11, but a
spliced PAF's col11 counts intron bases, so only intronless hits (retrocopies) passed: 107 multi-copy genes of 37,330,
which is a broken stratum, not a finding. Identity is now `1 - dv` (minimap2's divergence estimate). The broken table
(ours / ours_cur / StringTie recall single 22.1 / 22.0 / 18.1, multi 33.3 / 33.3 / 26.5 on 411 transcripts) is void.
