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

## Results (run once; `/mnt/linuxdisk/tmp/multicopy_gains/`)

### Q1: the gains are mostly NOT that copy's own expression

| | gorilla | human |
|---|---|---|
| gained transcripts | 162 | 284 |
| in **supported** genes (>= 2 reads certified to the copy) | 20 (12%) | 44 (15%) |
| in **contradicted** genes (reads belong to the other copy: shadows) | 73 (45%) | **200 (70%)** |
| in **undecidable** genes (identical over the reads) | 69 (43%) | 40 (14%) |
| genes: supported / contradicted / undecidable | 7 / 63 / 63 | 34 / 134 / 37 |
| good secondary records, tied: to this copy / other / abstain | 4 / 544 / 1,076 | 45 / 2,073 / 13,707 |
| good secondary records, strictly worse: to this copy / other / abstain | 46 / 2,975 / 3,107 | 400 / 33,871 / 8,964 |

**Reading (fixed rules):**
- Neither species reaches "real copies" (half of the gains in supported genes).
- **Human: shadows**, with 70% of the gains in contradicted genes.
- **Gorilla:** no class holds a majority (45% contradicted, 43% undecidable). The gains there are shadows or unresolvable, not real
  copies.

**Mainstays among the gained genes:**
- **Human, supported:** TBC1D3D, TBC1D3H, TBC1D3I, TBC1D3K, NPIPB5, GOLGA8F, SMN2, RBMY1A1, RBMY1D.
- **Human, contradicted:** NPIPB15, GOLGA6A, TSPY2 / 3 / 8 / 9 / 10, CDY1B.
- **Human, undecidable:** TBC1D3B, GOLGA6B.
- **Gorilla:** three TBC1D3G-like genes, all contradicted.

### Q2: the loci where reads multi-map are real duplications in the DNA

| genes in a segmental duplication (>= 90% identity) | gorilla (our GGO SD scan, stage-2 identity) | human (published CHM13 SD90) |
|---|---|---|
| **high-tie genes, span overlap (registered bar >= 90%)** | **1,286 / 1,423 (90.4%): holds** | **1,919 / 1,936 (99.1%): holds** |
| gained genes, span overlap | 104 / 133 (78.2%) | 202 / 205 (98.5%) |
| high-tie genes, >= 50% of the span inside an SD (descriptive) | 81.5% | 97.4% |
| high-tie genes with the scan's coarse identity estimate >= 90 (descriptive) | 82.9% | – |

- **Location of the high-tie genes:**
  - human: autosomes 1,444, chrY 337, chrX 118, PAR 37;
  - gorilla: autosomes 1,218, chrX 122, chrY 71, PAR 12.
- **Gorilla gained genes outside any SD** include a run of ~15 adjacent LOC129528xxx genes (undecidable, probably a tandem array
  the scan's 1-kb segment floor misses) and two PAR genes.

**Reading.** In human the loci where long reads multi-map are, almost without exception, segmental duplications in an independent
published annotation. In gorilla, by our DNA-only scan, 90% are (weaker on stricter cuts). The families are real duplicated DNA,
not artefacts of our family definition. What is NOT supported is that the extra transcripts secondaries build there are expressed
copies. Most are rebuilt from another copy's reads (human) or cannot be told apart (gorilla); 12-15% are supported, among them
several of the advisor's mainstays in human (TBC1D3 × 4, NPIPB5, SMN2, GOLGA8F, RBMY1 × 2).

Audit: below (Q1 holds with caveats; Q2 human holds with caveats; Q2 gorilla does NOT hold as stated; the reading above is withdrawn).

## Audit (metric-auditor, 2026-10-21): Q1 holds with caveats; Q2 human holds with caveats; Q2 gorilla does not hold as stated

Every number reproduces to the digit. The populations rebuild identically. The certificate is unchanged since register 1269.
Corrections:

- **The gorilla SD bar rested on a filter that removes nothing.**
  - Every row of our scan has stage-2 identity >= 90 by construction (blocks are kept at >= 90), so "90.4%" means "overlaps any
    scanned pair".
  - With the scan's own identity definition (est_id) it is 82.9%: **FAIL**.
  - Register 1278 kept SEDEF as gorilla's SD source. On SEDEF SD90, as a **post-hoc** check, 1,397 / 1,423 = **98.2%**.
- **Q2 is close to true by construction.** A read ties at >= 0.98 × best only where two loci share near-identical sequence, and
  Q2 had no null.
  - Base rate, genes with no tied reads: 15.3% in an SD (human), 23.6% (gorilla, SEDEF).
  - With the satellite-filtered CHM13 track (`SD90.nosat`), human falls to **77.3%**: chrY 18.4%, acrocentric chromosomes 68.8%,
    other autosomes 97.6%.
  - Of the 1,936 human high-tie genes, 69 lie in the X/Y PARs (allelic, not paralogous) and 371 are RefSeq pseudogenes.
  - "Not artefacts of our family definition" overreaches: no O1 family enters Q2. Its population is read ties.
- **"Supported" mainstays are certified 97-100% to the other copy.** Reads certified to this copy vs the other:
  - TBC1D3D 5 / 931, TBC1D3H 2 / 743, TBC1D3I 11 / 869, TBC1D3K 4 / 674;
  - NPIPB5 11 / 354, GOLGA8F 2 / 498, SMN2 5 / 593;
  - RBMY1A1 4 / 138, RBMY1D 3 / 211.
  
  Only 1 of the 34 supported human genes has as many reads certified to its own copy as to the other. The mainstay lists also
  matched on gene symbol only: counted by RefSeq description, human mainstay-family gains are 11 supported / 43 contradicted / 3
  undecidable.
- **"Contradicted = shadow" mostly restates minimap2's ranking.** Of the calls to the other copy, 94% (human) / 85% (gorilla) come
  from secondaries scoring strictly below the read's best, which minimap2 already placed elsewhere (outside O2's tied-only scope).
  - On TIED records alone, the gained transcripts split supported / contradicted / undecidable as human 9 / 59 / 216 (76%
    undecidable) and gorilla 2 / 25 / 135 (83% undecidable).
  - The instrument never tests reads whose primary lies at the gene itself.
- **Minor:**
  - The location split was computed by gene name, not coordinates (by coordinates: human PAR 69 / chrX 102 / chrY 320).
  - The acrocentric class was not computed.
  - Human is scored on RefSeq, not CAT.

**What stands:**
- **Q2:** the places where reads multi-map are strongly enriched for segmental duplications, 4-6× the base rate. Human: 97.6% on
  ordinary autosomes in the published track, and less for chrY / acrocentric satellite-adjacent regions. Gorilla: 98.2% in SEDEF
  (post hoc). That is expected when reads tie only where sequence is duplicated. It does not validate any O1 family.
- **Q1:**
  - The extra transcripts secondaries build in these loci are not shown to be the copies' own expression.
  - Where a read is truly tied, its copy is mostly undecidable (76-83% of gains).
  - Where it is not tied, minimap2 already put it elsewhere.
  - No mainstay gene is shown to be an expressed copy by this test.

**Correction (2026-10-21, from the audit of `PREREG_o1_families_sd_pairs`):** `CHM13_SDs.SD90.bed` is not a published track. It is a
SEDEF run on CHM13 v2.0 supplied by the user (fracMatch >= 0.90, >= 1 kb). Every "published" above should read "external SEDEF run".
SEDEF is independent of our code but not independent evidence (registers 650, 996).
