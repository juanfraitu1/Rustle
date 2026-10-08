# Pre-registration: loci in the Liftoff framework (Figure 8)

**Written 2026-09-25, before any comparison between a Rustle output and a Liftoff output exists.** User decision of
2026-09-25 12:15 (recorded in the figure-phase notes): every locus comparison is made in the Liftoff framework.

1. For each genome (human T2T-CHM13 v2.0, gorilla mGorGor1, chimpanzee mPanTro3, orangutan mPonPyg2), `liftoff
   -copies` lifts the genome's own RefSeq annotation onto the same genome. This is the **annotation-guided locus
   baseline**: every annotated record placed again, plus the extra copies Liftoff finds.
2. Every locus comparison uses Liftoff's own matching criteria, coverage and sequence identity per locus, and reports
   them as Liftoff reports its `coverage` and `sequence_ID` attributes.
3. Rustle's guided-mode loci are compared with Liftoff's (both start from the annotation and the genome; like for
   like). Rustle's de novo loci and copy catalog are scored against Liftoff's loci as a **reference, not a
   competitor**. Rustle's missing-copy flags are compared with Liftoff's extra copies.

Not chosen by the user: lifting the human annotation onto the apes.

## 0. What existed and what was seen before this file

- **One pilot, run by the main session (2026-09-25 12:18), cost only.** CHM13 chr21 records of type `gene`
  (Liftoff's default; no pseudogenes) lifted from chr21 onto chr21 with `-copies -sc 0.95 -p 4`: 36 s wall, 0.49 GB
  peak; the log printed 776 genes lifted, 257 extra copies, 0 unmapped
  (`/mnt/linuxdisk/tmp/rustle_figures_dev/liftoff_pilot/`). Its target was one chromosome, so it is not a unit of
  the design below and none of its numbers is reported.
- **Earlier Liftoff use** (read, not re-run): `docs/o1_ledger.md` §6iv (Liftoff v1.6.3 lifting a 2,334-gene CHM13
  v1.0 annotation onto v2.0; `-copies` found 68 extra loci; known paralogues are lifted individually, with no family
  notion) and §6iw (a Liftoff-style own-sequence realignment used as a family-edge signal; not adopted).
- No Rustle genome-wide catalog, de novo family table or missing-copy flag table exists yet for any of the six
  samples (`${work}/runs/<sample>/` holds assemblies only; checked with `ls`). No read-support count on any Liftoff
  locus exists.
- Annotation metadata read to write this file (feature-type counts of the four GFFs): gene + pseudogene records are
  human 41,561 + 17,002; gorilla 34,114 + 7,079; chimpanzee 34,730 + 7,085; orangutan 30,310 + 7,207.

Anything measured after this file is written goes into the amendments at the end, with the date.

## 1. The Liftoff baseline (one per species; independent of the RNA sample)

**Tool.** Liftoff v1.6.3 (conda env `liftoff`), with that environment's minimap2 2.24 passed explicitly (`-m`).
Rustle's own steps use minimap2 2.31; the two are never mixed inside one step.

**Input.** Target = reference = the species' genome FASTA (the registry's `fasta`, `figures/samples.tsv`);
annotation = the species' RefSeq GFF3 (the registry's `annotation_gff`: human `chm13v2.0_RefSeq_full` uncompressed,
`GGO_genomic.gff`, `PTR_genomic.gff`, `PPY_genomic.gff`).

**Records lifted.** Parent types `gene` and `pseudogene` (Liftoff's default is `gene` alone; `-f` adds
`pseudogene`), so the set equals the guided node set of Figure 7 (every gene and pseudogene record). Only these
records and their descendants are passed to Liftoff; other GFF features (cDNA_match, match, enhancer,
biological_region, region, ...) have no gene or pseudogene ancestor and cannot be lifted as genes, so dropping them
does not change what is lifted.

**Parameters (pre-registered; all Liftoff defaults except `-copies`, `-sc`, `-f`, `-p`).**

| option | value | meaning |
|---|---|---|
| `-a` | 0.5 (default) | a record is mapped only if ≥ 50% of its exon/CDS bases align (`coverage`) |
| `-s` | 0.5 (default) | ... and its exon/CDS columns are ≥ 50% identical (`sequence_ID`) |
| `-copies` | on | search for extra copies |
| `-sc` | **0.95** | an extra copy needs exon/CDS `sequence_ID` ≥ 0.95 to its source record |
| `-overlap` | 0.1 (default) | two placed records may overlap by ≤ 10% of the shorter one; an extra copy may overlap none on the same strand |
| `-d`, `-flank` | 2.0, 0 (defaults) | |
| `-mm2_options` | default: `-a --end-bonus 5 --eqx -N 50 -p 0.5` | |
| `-exclude_partial` | off (default) | mappings below `-a`/`-s` are written with `partial_mapping=True` / `low_identity=True` |
| `-p` | 4 | minimap2 threads |

**What `-copies` searches (from the v1.6.3 source, `run_liftoff.map_extra_copies`, `align_features`).** After the
annotated records are placed, every record's whole gene body (introns included; flank 0) is aligned again against the
**whole target genome**, whatever `-chroms` says, with minimap2 keeping the primary alignment and up to 50 secondary
alignments that score ≥ 0.5 of the primary (`-N 50 -p 0.5`). Only **end-to-end** alignments (the whole body aligned,
no clipping) become candidate copies. A candidate is lifted onto the record's exons/CDS and accepted when its exon/CDS
`sequence_ID` ≥ `-sc`; its `coverage` is not thresholded. A candidate that overlaps, on the same contig and strand, any
already placed record (its own source included) by ≥ 1 bp is re-placed or dropped. Consequences: copies on other
chromosomes are found; a record with more than 51 end-to-end copies is truncated at 51; an extra copy is by definition
**not annotated** (it overlaps no annotated record on its strand).

**Why `-sc 0.95`.** The default 1.0 admits only exact copies. The copies this thesis is about are the ones whose reads
are ambiguous or nearly so: the seeding rule admits alignments within 2% of the best score, the copy-assignment
figures report identity bands from 100% down to below 98%, and Soto et al. 2025 define recent segmental duplications
at ≥ 98% identity. 0.95 contains that whole band with margin. **Sensitivity rows:** every number that uses extra
copies is also reported with the extra copies restricted to `sequence_ID` ≥ 0.98, ≥ 0.99 and = 1.00 (the default).
These rows filter the `-sc 0.95` run; a separate run at a higher `-sc` could differ only through Liftoff's overlap
resolution between competing copies (said in the caption).

**Execution (the 10-minute rule).** One Liftoff call on the whole genome does not fit the machine rule (≤ 10 min, ≤ 20
GB per call). The run is split by the **reference annotation**: one call per contig, lifting that contig's records onto
the **whole** genome with `-copies`, so cross-contig copies are found exactly as in one run. A contig whose call would
exceed the budget is split into blocks of records at positions no record spans. The minimap2 index of the genome is
built once per species with Liftoff's own command (`minimap2 -d genome.mmi -a --end-bonus 5 --eqx -N 50 -p 0.5`) and
reused by every call (human reuses the index built for §6iv, header checked: k15 w10, the same options).

**Merge (emulating one run's overlap rule across calls).** Inside a call, Liftoff resolves overlaps itself. Across
calls, the merge applies Liftoff's rule for what one run would have seen:
- M1: an extra copy is dropped when it overlaps, on the same contig and strand, by ≥ 1 bp, the placed locus of any
  annotated record from another call (in one run that record would have been placed first, and the copy has a single
  alignment, so it could not be re-placed);
- M2: extra copies from different calls that overlap on the same strand: the one with the higher `sequence_ID` is
  kept, then the longer one, then the one whose source record comes first in GFF order;
- M3: annotated records from different calls placed onto overlapping positions on the same strand (a record lifted
  onto another record's copy) are counted and listed, and kept.

**Validation V-L1 (before the genome-wide runs; pre-registered bar).** Mini-genome = CHM13 chr20 + chr21 + chr22
(one FASTA; its own records): one Liftoff run versus three per-contig calls plus the merge, same parameters. Compared
as sets of (source record, contig, start, end, strand): annotated placements must be identical, and extra copies
must agree with Jaccard ≥ 0.98. If V-L1 fails, the per-contig design is reported as an approximation on the figure,
and a single run per species is requested on the cluster.

## 2. The Liftoff locus table

Per species, one row per placed parent feature:
- **annotated locus** (copy tag `_0`): class **in place** when it lies on its own record's contig with reciprocal
  span overlap ≥ 0.5; **moved** otherwise (placed onto another copy); **partial** when Liftoff marks it
  `partial_mapping` or `low_identity`; **unmapped** records come from the unmapped list.
- **extra copy** (tag `_k`, k ≥ 1).
- Columns: species, source record ID and Name, record type, `gene_biotype`, source contig/start/end/strand, placed
  contig/start/end/strand, class, `coverage`, `sequence_ID`, `extra_copy_number`, exon union (merged exon children;
  CDS if the record has no exon; else the span), exonic length.
- **Strata** (for display and every per-stratum number): protein-coding; lncRNA; pseudogene (record type
  `pseudogene`, or a biotype containing "pseudogene"); other; and, crossing these, **short** = exon union < 200 bp
  (Iso-Seq reads and the catalog's 200-bp representative floor cannot see these; they are reported, never scored for
  Rustle).
- **Reference set for scoring Rustle** = in place + moved + extra copies, exon union ≥ 200 bp. Partial and unmapped
  records are excluded (counted).

## 3. The matching rule (every comparison)

Liftoff's `coverage` of a record is the fraction of the record's exon/CDS bases that are aligned at the placed locus;
its `sequence_ID` is the fraction of identical columns over those bases, with unaligned and inserted bases counting
against it. Two loci compared **in the same genome** share their bases, so a coordinate match has no mismatch and no
insertion: Liftoff's `sequence_ID` of a locus G placed on locus R's exons equals its `coverage`. The two criteria
therefore reduce to one number:

    cov(G | R) = |X(G) ∩ X(R)| / |X(G)|,   X = exon union, same contig, strand ignored

- **G is found by R** iff cov(G | R) ≥ 0.5 (`-a`) for some R. **R lies at a Liftoff locus** iff cov(R | G) ≥ 0.5 for
  some G. Strand agreement is reported, not required (Rustle's single-exon loci may carry no strand).
- Exon unions: de novo locus = union of the exons of every assembled transcript of the `gene_id`; catalog copy = its
  `exons` column; guided candidate = the exon blocks of its transcript hit, else its span; missing-copy flag = the span
  of the flagged annotated record, and the span of the consensus's best other hit (`other_locus`).
- Where a Rustle object carries its own alignment identity (guided candidates, flag consensus hits), that identity is
  reported beside Liftoff's `sequence_ID`, never used as a gate.

## 4. Comparisons

Samples: the six of `figures/samples.tsv`; each is scored against its species' Liftoff table. Species and samples are
never pooled.

**C1. Guided mode versus Liftoff (per species; like for like: annotation + genome, no reads).**
- *Annotated loci.* Rustle's guided node set is every gene and pseudogene record, taken as given
  (`_o1_recovery.gene_regions`; `mcl_families` folds bodies whose exon unions overlap). Liftoff re-places them.
  Reported: in place / moved / partial / unmapped per species. No Rustle claim: guided mode does not move records.
- *New loci.* Rustle's guided mode, as defined in `docs/seeded_family_definition.md` §0★★ clause 1, adds candidate
  loci found by homology search from the annotated records (chain-first construction, §6js; implemented in
  `bench/guided_pipeline.py`: `transcript_hits`, `gene_body_chains`, `chain_first`). Figure 7's guided recipe does not
  run this search. Here it is run genome-wide with **every** gene and pseudogene record as a seed, unchanged
  otherwise: each record's unit (its longest NM_/NR_ transcript, else any transcript, else its exons, else its span)
  aligned with `-x splice`, and its CDS envelope (else exon envelope, else span) with `-x asm20`, both
  `-c -N 100 -p 0.1` against the species' splice index; transcript hits need identity ≥ 0.80 and ≥ 50% of the unit;
  gene-body chains need identity ≥ 0.80 and ≥ 50% of min(query, extrapolated span); hits overlapping any seed's span
  are blocked; loci are built chain-first. The candidates are then compared with Liftoff's extra copies by §3, both
  ways, at the finder's own identity (≥ 0.80) and restricted to identity ≥ 0.95 (the `-sc` row). Liftoff extra copies
  that lie inside an annotated record's span on the other strand are reported separately (the finder blocks every
  hit inside a seed's span; Liftoff blocks only same-strand overlaps).

**C2. De novo loci versus Liftoff (per sample; reference, not competitor).** Rustle's de novo loci = `gene_id` groups
of the sample's genome-wide assembly (`make.py runs` stage `assemble`, product `gtf`). Reference loci are restricted
to those with read support in the sample: **≥ 2 reads whose primary alignment (`-F 2308`) has an aligned block on the
locus's exon union** (the rule of Fig. 6d). Reported per sample: the fraction of read-supported reference loci found
by a de novo locus, for annotated loci (in place, moved) and extra copies separately, by stratum; and the fraction of
de novo loci that lie at a Liftoff locus (the rest are outside every Liftoff locus, which is not an error).

**C3. Copy catalog versus Liftoff (per sample).** The same as C2 with the sample's genome-wide copy catalog
(`catalog` stage, `copies`) as R, plus a pair test: for every (source record, extra copy) pair of Liftoff where both
loci are found by catalog copies, are the two in one catalog family (any matched copy of one and any matched copy of
the other share a `family_id`)?

**C4. Missing-copy flags versus Liftoff's extra copies (per sample).** From the `flag` stage table
(`missing_copy.tsv`; one row per annotated locus with ≥ 10 reads): (i) the fraction of fired loci, and of
reference-absent candidates, whose annotated record has ≥ 1 Liftoff extra copy, beside the same fraction among all
scanned loci (enrichment); (ii) among those, the fraction whose consensus best other hit (`other_locus`) covers
≥ 50% of the span of one of that record's extra copies (the flag's own home search found the copy Liftoff found);
(iii) the reference-absent candidates whose record has an extra copy that `other_locus` does not cover, listed for
review. Liftoff's extra copies are in the reference; the flag looks for copies the reference lacks. So (ii)-(iii)
measure how often an unannotated copy present in the reference could explain a flag.

## 5. Claims, bars and what changes them

| id | claim | bar | if it fails |
|---|---|---|---|
| L1 | Liftoff's self-lift is a usable baseline | per species, ≥ 99.0% of gene + pseudogene records placed in place | the species' figure says so; records not in place are excluded from C1-C4 and listed |
| G1 | Rustle's guided candidate search finds what `-copies` finds | per species, ≥ 0.90 of Liftoff's extra copies (exon union ≥ 200 bp, not inside another record's span) found by a candidate at identity ≥ 0.80 | Liftoff's end-to-end body search reaches copies the transcript/CDS search misses; the missed copies are listed by stratum |
| D1 | de novo finds unannotated copies as well as annotated ones | per sample, sensitivity on read-supported extra copies ≥ sensitivity on read-supported in-place annotated loci − 0.10 | read placement absorbs unannotated copies into annotated ones; report by `sequence_ID` band |
| F1 | the catalog's families contain Liftoff's copy relation | per sample, ≥ 0.90 of (source, extra copy) pairs with both loci in the catalog share a family | the catalog splits near-identical copies; list the pairs with their catalog copies |

Descriptive only (no bar, no base rate known before looking): C1 Rustle-only candidates; C2 fraction of de novo loci
at a Liftoff locus; C3 catalog sensitivities; all of C4.

**Intervals.** Every per-sample proportion gets a 95% interval from 2,000 bootstrap resamples of **source records**
(an annotated locus and its extra copies are one resampling unit; random seed 20260925).

**Exposure.** Nothing here was scored before. The family and catalog rules were developed on human chr16 and gorilla
chr20 (NC_073244.2); the guided finder's thresholds (0.80 / 0.50) on NPIP, TBC1D3 and AMY (human). D1 and F1 are
also reported for the genome minus those contigs. Human A119b and human testis share the human Liftoff table; the
two gorilla samples share the gorilla table.

## 6. Stop rules

- A Liftoff call, a guided-search mapping shard or a read-support contig pass that cannot fit 10 min / 20 GB is split
  (blocks of records, query shards, contig parts). If a single unit cannot fit, it stops with a message for the
  cluster; nothing is dropped silently.
- If V-L1 fails its bar, §1 applies (approximation stated on the figure; single run requested).
- Numbers are new register rows; the §6iv/§6iw rows are not overwritten.

## Amendments

**Amendment 1 (2026-09-25, 13:45, before any Liftoff output of this design or any flag table existed).** Reading
`missing_copy::verdict` showed that the flag already has a verdict for the case Liftoff's extra copies describe:
`unannotated_paralogue` = the consensus's best hit elsewhere in the reference (≥ 80% of the consensus aligned) is
more identical than its host. C4 gains a fourth descriptive item: (iv) of the fired loci with verdict
`unannotated_paralogue`, the fraction whose `other_locus` covers ≥ 50% of the span of (a) an extra copy of the same
record, (b) any Liftoff extra copy, (c) the placed locus of another annotated record. No bar. The flag table's loci
are the annotation's gene/pseudogene records (`missing_copy::load_loci` on the GFF: `locus` = the record ID), so a flag
row maps to its Liftoff source record by ID.

Also fixed here, from the code: the cost pilot below is the first shard of the human run itself (records of chr1,
block 1 of 3: 1,829 records, lifted onto the whole genome), not a separate experiment.

**Amendment 2 (2026-09-25, 14:20, after the cost pilot; before V-L1 and before any Rustle comparison).**
- *Cost pilot (human shard chr1.b1, 1,829 records, 51.8 Mb of gene bodies, whole-genome target, under the heavy
  lock):* 538 s wall, 9.35 GB peak. Each of Liftoff's two minimap2 passes took 261 s at 1.2 CPUs of 4 (one or a few
  query sequences dominate); the Python lifting took about 16 s. A probe of the longest human record alone (RBFOX1,
  2.48 Mb) took 99 s including a 30 s index load, one alignment, so body length alone does not set the cost; the
  records with many long secondary alignments (segmental duplications) do. Counts seen in this shard (Liftoff only,
  no comparison): 1,828 records placed in place, 1 unmapped, 205 extra copies before the cross-call merge.
- *Exact halving of the cost.* Liftoff's `-copies` step re-extracts the same record sequences into the same file and
  re-runs the same minimap2 command against the same index. Liftoff's `-m` is given a wrapper
  (`_liftoff.MM2_CACHE_WRAPPER`) that returns the first pass's SAM when the query's md5, every argument and the
  index's size and mtime are identical. Checked before use: the chr21 pilot re-run through the wrapper gives a
  lifted GFF3 identical to the pilot's (0 differing lines outside the command-line header) and the same SAM (`cmp`).
- *Shards that do not finish.* A Liftoff call killed at the hard limit (590 s) is split in two at the record-free cut
  nearest its middle and re-run (labels `.s1`/`.s2`); a single record or a shard with no record-free cut stops for
  the cluster (§6). V-L1's arm C (blocks) covers this kind of split.
- *Cost projection:* about 5 min per shard with the wrapper; human 36 shards (about 3 h of calls), each ape 25-30
  shards (about 2-2.5 h each), plus one index build per ape (about 5 min, about 12 GB).

**Amendment 3 (2026-09-25, 14:40): V-L1 result. The bar is missed narrowly; the pre-registered consequence applies.**
Mini-genome CHM13 chr20 + chr21 + chr22, 3,759 records, same parameters, run through the wrapper of amendment 2
(`${work}/liftoff/vl1/vl1_report.json`):

| arm | annotated placements | extra copies | bar |
|---|---|---|---|
| A: one run | 3,750 placed | 442 | reference |
| B: one call per contig + merge | 3,751 placed; 11 records placed at a different position than in A | 445; 439 shared with A (Jaccard 0.980) | FAIL (placements not identical; 0.9799 < 0.98) |
| C: blocks of <= 600 records + merge | 3,751; 21 records placed differently | 445; 438 shared (Jaccard 0.976) | FAIL |

Every difference lies in chr21:3.1-5.1 Mb or chr22:4.8-5.6 Mb, the acrocentric short arms, where records sit in arrays of
identical copies (rDNA units and LOC124908xxx copies): which of several identical positions a record or a copy
receives depends on which other records compete for them in the same call. Elsewhere the three arms agree exactly.
As pre-registered: the genome-wide baseline is run per contig (blocks of records) and **reported on the figure and in
the table notes as an approximation of one Liftoff run**, with these numbers; a single run per species is requested
for the cluster (it does not fit the 10-minute rule here: the human run is about 3 h of mapping). No rule is changed
after this result. Consequence for the comparisons: extra copies inside identical-copy arrays can move between
equivalent positions; the scoring counts a locus as found by coverage, so a moved copy can change a found/not-found
call only where the array copies are not all covered alike.

**Amendment 4 (2026-09-25, 15:10): cost of the guided candidate search (C1), before any candidate was compared.**
Human query preparation: 58,563 units (146 Mb) and 58,563 envelopes (1.57 Gb), 43 s, 2.4 GB. First units shard
(20,000 units, 50 Mb, `-x splice -c -N 100 -p 0.1 -t 4` against `target.splice.mmi`): 381 s, **20.5 GB peak**, above
the 20 GB rule. Shards are therefore halved (10,000 units; 25 Mb of envelopes per call), which changes no alignment
(each query is mapped independently; only the batch held in memory shrinks). Projection: about 72 calls per
species, 2-5 min each, about 3-5 h per species; the envelope (`-x asm20`) shards are unmeasured. The search stays
optional (`fig8_guided_finder=1`); until it runs, panel b and claim G1 say "not run". The pilot's hits were deleted
unread.

**Amendment 5 (2026-09-25 17:10, before any comparison between a Rustle family and a Liftoff locus existed; user
decision of 16:00).** Seen when written: no Liftoff self-lift is merged for any species (human 1 of 36 shards,
gorilla, chimpanzee and orangutan 0), no `fig8_*` table exists, no families-stage output exists for any sample, and
the two genome-wide legacy catalogs (chimp_PTR, human_testis) were never compared with anything.
- The user made the driver's `families` stage the ONE default de novo family definition (reads → seeded assembly
  loci → one representative per locus → `mcl_families --from-gtf`, exon-sum ≥ 0.60, MCL 2.8), whose copy table
  (`<id>.fam.copies.tsv`, one copy per member locus with its representative's exons;
  `docs/PREREG_families_copy_table_2026-09-25.md`) copy assignment consumes; `gw_family_catalog` (`catalog`) is
  legacy.
- **C3 and F1 therefore apply to the default families' copy table** in place of the catalog: C3's sensitivities use
  the copies' `exons` column (as for the catalog), and F1's pair test asks whether the covering copies of the source
  locus and of the extra copy share a `family_id` of the default families. Bar unchanged (≥ 0.90 of the pairs whose
  two loci are both covered). The copy table holds only loci that are members of a family, so C3's sensitivity is
  "read-supported Liftoff loci covered by a family member", not "by any Rustle locus" (C2 keeps that).
- The legacy catalog is scored by the same code when its product exists, as labelled secondary rows (`catalog`);
  it carries no claim.
- The same (record, extra copy) pairs, restricted to pairs whose two loci are both read-supported in the sample (the
  C2 rule), are used as a family reference by Figures 6s and 7 (defined in
  `docs/PREREG_genome_wide_families_2026-09-25.md`, Amendment 1, items 3–4): there the denominator is every
  read-supported pair, not only the covered ones.
