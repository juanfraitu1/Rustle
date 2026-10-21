# Pre-registration: are the multi-copy gains real copies, and are those loci real duplications? (2026-10-21)

**User (2026-10-21), after register 1366:** "lests check the multi-copy gains, also can we prove they are real multicopies? My
advisor is only familiar with NPIP and other mainstays and if he sees that O2 doesnt apply there he might think all other
'families' are artifacts".

Written before the instrument runs on these genes and before any SD overlap is computed.

## Population (fixed from existing tables)

- **Source:** `figures/data/fig3_ref_tie.tsv` (genome-wide; gorilla OR6737 on GGO, human A119b on CHM13 scored against RefSeq,
  as in register 1366).
- **High-tie strata:** expressed multi-exon reference transcripts with tie fraction > 0.5 at 0.98 × best AS, or with no primary
  molecule.
- **Gains (G):** high-tie transcripts that Rustle with secondary seeding matches and Rustle primaries-only does not. Gorilla has
  162 transcripts in 126 genes; human has 284 in 203 genes.
- **High-tie loci (H):** every gene of a high-tie transcript. This is the wider population for Q2.

## Q1: are the gains built from that copy's own reads? (O2's certificate, register 1269's instrument unchanged)

- **Instrument:** `bench/tool_difference/s1_o2.py` (sha1 6b6831e3...), with the gained genes as its gene list.
  - Every good secondary record (AS >= 0.98 × the molecule's best) overlapping a gained gene's exons (copy B) is paired with the
    molecule's primary record (copy A).
  - It certifies on the PSV columns the read spans (Binomial, p < 1e-3, strict maximum, abstain otherwise).
- **Gene classes:**
  - **supported:** >= 2 records certified to B;
  - **contradicted:** 0 to B and >= 2 to A;
  - **undecidable:** otherwise.
- **Reported:** per species, gained transcripts by their gene's class, and records by verdict.
- **Expected, from register 1269:** mostly contradicted or undecidable.
- **Reading, fixed now:**
  - the gains are "real copies" only if at least half of the gained transcripts lie in supported genes;
  - "shadows" if most lie in contradicted genes;
  - "unresolvable" if most lie in undecidable genes, meaning identical copies where the read cannot say which copy produced it.

## Q2: are these loci real duplications in the DNA? (independent of O2, O1 and the RNA reads)

A gene counts as **in a segmental duplication (SD)** when its span overlaps an SD interval of >= 90% identity:
- **human:** the published CHM13 SD track (`CHM13_SDs.SD90.bed`, SDs >= 90% identity; we did not produce it);
- **gorilla:** our genome-wide DNA SD scan of the GGO T2T (`sd_cache/.../genome.sd.pairs.tsv`, both sides of every pair with
  identity >= 0.90), which is sequence-only and blind to the annotation and the reads.

**Also reported:**
- each gene's location class:
  - pseudoautosomal (X/Y PAR: two placements by construction in a male reference);
  - chrY ampliconic;
  - acrocentric short arm;
  - other.
- the named mainstays among the genes: NPIP, TBC1D3, GOLGA6/8, SMN, LRPAP1, and the YAGs (TSPY, RBMY, DAZ, HSFY, BPY2, CDY).

**Bar, fixed now:** "the high-tie loci are real duplications" when >= 90% of H's genes lie in an SD of >= 90% identity, in each
species. The same fraction is reported for G.

## What this can and cannot show

- **Q2** shows whether the places where reads multi-map are duplicated sequence in a T2T genome. That is real, independent of
  our family definition.
- **Q2 cannot show** whether every copy is expressed or present in the sample's own genome: OR6737 and A119b are not the reference
  individuals.
- **Q1** shows whether the extra transcripts are expressed copies or paralog shadows.

## Rules

- Descriptive with fixed readings; nothing adopted or tuned.
- Heavy steps in the foreground, one at a time, each under 10 minutes. The human BAM is 96 GB; if its primary-record scan cannot
  finish in time, an amendment says how it is split.
- The metric-auditor checks the result before it is reported. Never pool human and gorilla. Push only when asked.

## Amendment 1 (2026-10-21, before the human run)

- **Instrument version:** the instrument in the tree is register 1270's version (3d4553af, sha1 a79218c0). Next to register
  1269's (da774afe, sha1 6b6831e3) it only adds the Eichler AS-margin report and the `--stratum` label; the certificate is
  byte-identical.
- **Human primaries scan, split:** the 96 GB human BAM cannot be scanned for primaries within one 10-minute call. The scan is
  split by chromosome group into 5 foreground calls, using two new options on `s1_o2.py` (the default path is unchanged):
  - `--names-only` writes the molecule names and stops;
  - `--primaries FILE` reads the primary records from a file instead of scanning.
- **Check:** on gorilla, the split path reproduces the default run's `records.tsv`, `genes.tsv` and summary byte for byte.
