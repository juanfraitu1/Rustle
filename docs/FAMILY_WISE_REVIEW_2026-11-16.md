# Family-wise assembly review, and the `names` / `eqclasses` / `psi` stages (2026-11-16)

Descriptive session, no registered bars and no register rows. Tables, scripts and the running summary are in the
workspace: `/anvme/workspace/v138da13-rustle/famwise_look_2026-11-16/` (`SUMMARY.md`, `copies.tsv`, `kb/`, `eq/`).
Assemblies used are the stored default runs, made before transfrag extension became the default.

## Current state
- Family-wise assembly has nothing copy-specific to borrow in the families looked at: either the primary alignment is
  already on the right copy, or the copies cannot be told apart from the reads. The honest product at tied copies is
  one structure set, the set of compatible copies and a family-level count.
- That product now exists as the opt-in stage `eqclasses` (tie classes, per-locus count ranges, EM for quantification
  only). Its tie bands are read on minimap2's `ms` score since 2026-11-16 (user; section 4a), which removes the
  gene / processed-pseudogene ties of AS in the exact band. The bands are still each wrong somewhere; neither is validated.
- `--names ANNOT` / stage `names` copies annotation names onto the assembled GTF as attributes, independent of `--guide`.
- Opt-in stage `psi` (section 8): PSI of each intron within its bunch, with the range tied reads allow; at multi-copy
  family loci 0.224 of the bunches of gorilla OR6737 are not measurable (range > 0.1) against 0.0025 elsewhere.
- Not worth a feature (counted, below): short introns, non-canonical junctions, junction repair.

## 1. Where primary-only tools lose copies (A119b on CHM13)
`copy_table.py` (workspace): per annotated copy, read records on its exons and strand classed by the molecule's
genome-wide best / second AS, and models in ours (r1666 `g_off` default arm), StringTie 3.0.1 `-L` and FLAIR collapse.
Both tools ran on the full BAM; FLAIR's own filter keeps one line per primary alignment.

| family | read support at the copy | copies | ours | StringTie | FLAIR |
|---|---|---|---|---|---|
| TSPY array (chrY) | >= 1 read best only here | 8 | 8 | 8 | 6 |
| TSPY array | ties only | 24 | 24 | 10 | 7 |
| TSPY array | secondaries within 0.98 only | 2 | 2 | 0 | 0 |
| GAGE (chrX) | >= 1 read best only here | 10 | 10 | 9 | 9 |
| GAGE | ties only | 12 | 12 | 7 | 6 |
| LILR (control, no ties) | all reads unique | 10 | 5 | 4 | 5 |

- Ten copies have no primary record (TSPY 8, GAGE 2): StringTie 0, FLAIR 0, ours 10.
- TSPY: 274 molecules with a primary in the array over 34 copies; ours emits 237 transcripts = 11 distinct structures,
  six of them on 25-44 copies each. One structure set repeated on every tied copy, not 34 copies resolved.
- GAGE rows are not independent (GAGE12F spans GAGE12E; two 117 bp GAGE12B fragments). The `clean` column fails for
  ours there on a longer terminal exon. TSPY2 (15 unique primaries, no model in ours) is unexplained.
- CHM13 with HG002's Y is not this individual: "best only here" is not proven origin.

## 2. The same with per-read origin (KB3781 on its own genome)
Four tied pairs chosen from the score table (MAGED4, NSMCE1-like, SPIN2, FAM156). StringTie and FLAIR on a sub-BAM of the
copies +- 100 kb. Origin = copies whose loci on the animal's haplotypes (mat, pat) hold the read's best hits, NM slack 0.

- 520 reads have a single origin copy by this truth; the primary is on it for all 520. All three tools model every copy.
- MAGED4: all 1,089 judged reads match both copies equally; primaries split 655 / 687 by chance.
- FAM156: 175-195 secondaries per copy sit within 0.98 of the best score (inside the default seeding pool) although
  each read's origin is unique.
- Corrections after reading the haplotype hits:
  - The 110 NSMCE1-like reads with origin "elsewhere" come from the maternal second copy, which the slack-0 copy loci
    left unclaimed (both reference copies' best maternal hit is the first copy). By position they are placed correctly.
  - The 141 reads "equal on both copies" come from the maternal first copy, which is 3 mismatches from each reference
    copy (the reference holds the paternal alleles on that chromosome); 116 of their primaries sit on the second copy.
    So the primary is not always right even on the animal's own genome.
  - The X-linked pairs have usable hits on the maternal index only.
  - The truth compares mismatch counts only: SPIN2 has 130 primaries that prefer one copy by score yet count as equal.

## 3. Stage `names` (`tools/gtf_names.py`, `tools/test_gtf_names.py`)
`--names ANNOT.gff3|gtf[.gz]` on `assemble` / `all`, or `rustle_pipeline.sh names --names ANNOT --out PREFIX` alone
(no BAM or FASTA). Writes `PREFIX.named.gtf` and `PREFIX.final.named.gtf`: the same lines with `ref_gene_name`,
`ref_gene_id`, `ref_transcript_id` (exact intron chain only), `ref_match` (`exact` | `overlap`) and `ref_genes` (when a
model touches several genes) appended. Ids, coordinates and every other product are unchanged; nothing reads the named
files. Without an exact chain the name is the same-strand gene sharing most exon bases, so fragments and readthrough
models take the name of the gene they overlap most; ends are not compared.
KB3781 vs RefSeq: 72,460 transcripts, 25,844 exact, 44,466 by overlap, 2,150 outside (34 s).

## 4. Stage `eqclasses` (`tools/eq_classes.py`, `tools/test_eq_classes.py`)
`rustle_pipeline.sh eqclasses --bam B --fasta G --out PREFIX [--names ANNOT]`, opt-in, never part of `all`. One BAM pass
over the loci `families` reads, using the score table.
- A read's placements = loci of its records whose AS ties its best; the primary record is always a placement; a
  placement on no locus is `outside`. Two bands side by side (user, 2026-11-16): exact ties, and AS >= 0.98 x best
  (`RUSTLE_EQ_RATIO`; columns `_w`).
- `PREFIX.eq.classes.tsv`: band, reads, the set of loci. `PREFIX.eq.loci.tsv` per locus and band: `placed` (reads
  placed only there: "aligns best only here", not origin), `max` (placed + every tied read), `em` (counts-only EM;
  quantification only, no read is assigned), tie group with its size, contigs and placed reads. Loci mostly inside the
  exclude BED are `repeat` and join no group. Where `group_placed` is small the EM is the starting split: read the range.

| sample | band | reads in classes | classes | tie groups | largest group |
|---|---|---|---|---|---|
| KB3781 | exact | 159,218 | 724 | 244 | 34 loci |
| KB3781 | 0.98 | 225,502 | 1,212 | 344 | 34 loci |
| A119b | exact | 185,228 | 27,613 | 778 | 83 loci |
| A119b | 0.98 | 467,055 | 73,495 | 971 | 144 loci |

KB3781 5 min, 0.3 GB; A119b 1 h 53 min, 1.3 GB (64 processes). A119b: 39 of 70 TSPY-named loci in one exact group
(46 wide); 20 of 22 GAGE loci in one exact group (21 wide).

Check on the four KB3781 pairs (seen data = development only; origin by position, unjudged reads not counted):

| copy | origin | exact [placed, max], EM | 0.98 [placed, max], EM |
|---|---|---|---|
| NSMCE1-like LOC129527614 | 253 | [136, 215], 158 | [66, 496], 286 |
| NSMCE1-like LOC101123970 | 221 | [344, 423], 401 | [63, 493], 273 |
| FAM156A | 137 | [175, 189], 182 | [0, 384], 192 |
| FAM156B (LOC101125357) | 160 | [195, 209], 202 | [0, 384], 192 |
| MAGED4 (each copy) | undecidable | [9, 1332] / [10, 1333], 636 / 706 | [0, 1342], 671 |

- Exact ties put the NSMCE1 origin outside the range of both copies (the maternal-allele reads above); the 0.98 band
  contains it but turns FAM156, which exact ties get right, into [0, 384].
- AS ties join spliced genes with their processed pseudogenes, because AS charges for introns: the largest KB3781 class
  is 80,611 reads between FTH1 and LOC115931247 (GAPDH: 11,225). EM resolves them from the placed reads (FTH1 93,543 of
  93,715) but the pseudogene's range is [28, 80,639]. A band on `ms` does not tie them: section 4a.
- The 14,747-read class between a chrX and a chrY locus is the pseudo-autosomal region (the reference carries a Y).
- No registered accuracy test exists. One would need tie groups other than these four pairs and a stated scoring of
  unjudged and equal reads.

## 4a. Bands on `ms` (default since 2026-11-16, user decision)
`eq_classes.py --score ms|AS`, default `ms` (minimap2's max-scoring segment, the score it ranks records on); the driver
stage passes nothing, so it runs on ms. Sections 4's tables are the AS run (`--score AS` reproduces them, except that a
read whose primary is not its best record is now always tie-capable). The score table must have its ms columns (fields
7-9); a table without them is refused. Outputs: workspace `eq/kb.ms.ms.*`, `eq/a1.ms.*`.

| sample | band | reads in classes | classes | tie groups | largest group |
|---|---|---|---|---|---|
| KB3781 | ms exact | 21,453 | 407 | 116 | 34 loci |
| KB3781 | ms 0.98 | 142,997 | 884 | 240 | 34 loci |
| A119b | ms exact | 117,172 | 19,869 | 678 | 133 loci |
| A119b | ms 0.98 | 423,701 | 77,093 | 867 | 150 loci (22 contigs) |

A119b on the ratio band: 2 h 41 min, 1.4 GB.

KB3781, same four pairs and the pseudogene classes (AS values in section 4):
- FTH1 exact [93,711, 93,717], its pseudogene LOC115931247 [26, 30]; GAPDH [14,398, 14,398]. The AS classes of 80,611
  and 11,225 reads are gone from the exact band.
- NSMCE1-like, FAM156 exact ranges are identical to AS ([136, 215] / [344, 423]; [175, 189] / [195, 209]): these copies
  differ by mismatches, not splice shape, so the NSMCE1 origin (253 / 221) is still outside the exact range. Not fixed.
- MAGED4 [12, 1,328] / [14, 1,330]: undecidable on either score. The pseudo-autosomal class persists (14,750 reads).
- The 0.98 band on ms brings gene / pseudogene ties back where the pseudogene is near-identical: its largest class is
  54,506 reads between EEF1A1 and LOC101142679 (EEF1A1 [25,723, 82,717], EM 82,717), then SEPTIN7 9,313. NSMCE1-like
  becomes [6, 553] each. A ratio of 0.98 is a window of ~40 on a 2,000-score read; the project's standing ms near-tie
  definition is an absolute window of 10 (O2, row 1508). Which wide band to keep is open; no third band was built.
- Primary below its best record: 1,421 molecules on ms vs 415,293 on AS (KB3781); 26,745 vs 230,174 (A119b).

## 4b. Wide band = within 10 of the best `ms` (default since 2026-11-16, user decision)
The 0.98 ratio is replaced as the default wide band by an absolute window: `eq_classes.py --window W` (default 10),
driver `RUSTLE_EQ_WINDOW=W`; a ratio band is still available with `--ratio R` / `RUSTLE_EQ_RATIO=R` (overrides the
window). Outputs: workspace `eq/kb.w10.*`, `eq/a1.w10.*`; code snapshot `eq/tools_w10/`.

| sample | wide band | reads in classes | classes | tie groups | largest group |
|---|---|---|---|---|---|
| KB3781 | ms within 10 | 51,477 | 746 | 209 | 34 loci |
| KB3781 | ms 0.98 (4a) | 142,997 | 884 | 240 | 34 loci |
| A119b | ms 0.98 (4a) | 423,701 | 77,093 | 867 | 150 loci |
| A119b | ms within 10 | 269,608 | 51,665 | 898 | 145 loci |

A119b on the window: 49 min, 0.8 GB (64 processes).

KB3781 (seen data, development only):
- Gene / pseudogene ties of the 0.98 band are gone: EEF1A1 [82,497, 82,716] (was [25,723, 82,717]), SEPTIN7
  [13,666, 13,818], FTH1 [92,900, 93,745], GAPDH [14,397, 14,398].
- NSMCE1-like [132, 420] / [139, 427], EM 272 / 287: the origin (253 / 221) is inside both ranges, which the exact band
  misses, and the range is narrower than 0.98 gave ([6, 553]).
- FAM156 is still [0, 384] each (exact: [175, 189] / [195, 209], right): its copies sit within 10 ms of each other, so
  no window that catches the NSMCE1 allele reads (3 mismatches) keeps FAM156 apart. Read both bands.
- MAGED4 [0, 1,341] / [1, 1,342]: undecidable, as before. Largest wide classes: pseudo-autosomal 14,750, then a
  same-chromosome pair of 7,075 reads (NC_073239.2 35.40 / 35.52 Mb).
 that closed three ideas (A119b default assembly vs RefSeq on CHM13)
- Short introns: 623 of 349,884 distinct RefSeq introns are < 70 bp. 1-3 bp (90) are frameshift-correction gaps that
  reads run through; 51-69 bp (284) are called as introns by minimap2; 21-50 bp (230) appear as a deletion in 94 reads
  genome-wide. Same conclusion as register row 1085. No feature.
- Non-canonical junctions: the assembly emits none (all GT-AG, GC-AG or AT-AC); RefSeq has 1,871 others (0.5%,
  including the 1-3 bp gaps). No feature.
- Junction repair: of 131,211 wrong chains with an unannotated intron, 4,423 have every such intron within 10 bp of an
  annotated one; moving them gives an annotated chain for 2,571, already emitted in all but 203. About 200 chains gained.
  The rest: 126,788 wrong chains with a canonical unannotated intron far from any annotated one, 49,591 with annotated
  introns only (combination or fragment).

## 6. Judgement recorded, not measured
- Short-read RNA would not help copy assignment (shorter reads tie more) and adds little to assembly beyond what the
  chain score does. Short-read DNA of the same individual could separate alleles from copy differences for samples
  aligned to another individual's genome; test on HG002 (public Illumina, haplotype truth) before acquiring ape data.
- Read-level EM stays closed (row 1648); EM is used here for counts only.

## 7. Environment
Compute nodes have no `/usr/bin/time`. GNU time 1.10 is installed with `pixi global install time`
(`~/.pixi/bin/time`); the driver's `all` falls back to it for the per-stage peak memory.

## 8. Stage `psi` (`tools/psi_bunches.py`, `tools/psi_family.py`, tests beside them)
`rustle_pipeline.sh psi --bam B --fasta G --out PREFIX`, opt-in, never part of `all`, nothing reads its tables. Needs
`PREFIX.gtf` and the score table. Measurement of one sample: no test between samples (no replicates; MntJULiP or
LeafCutter would be the tools), intron retention is not seen.

- Unit: a bunch = the introns sharing a donor, or an acceptor (MntJULiP). Counts come from the CIGAR N operations of the
  primary alignment and of the secondaries tied in the band. Per band (exact, and `_w` = wide, the bands and switches of
  section 4): `placed`, `max`, `psi_placed`, `psi_lo`, `psi_hi`. A PSI is quoted with its range; a wide range means not
  measurable, and `placed` is not origin.
- OR6737, whole sample (58 s, 0.5 GB, 64 threads; workspace `psi_2026-11-16/or.psi.tsv`): 62,370 bunches; range wider than
  0.1 in 409 (exact band) / 782 (wide band). Of the wide-band 782, 633 are among the 2,821 bunches at multi-copy family
  loci (0.224) and 149 among the 59,549 others (0.0025).
- Between copies: after stage `exonmap`, `psi` also writes `PREFIX.fam.psi.tsv`, the psi rows keyed by exon columns
  (family, anchor, column of the anchor, column of the other end), so one junction is read across the copies of a family.
  `differs` = the copies' wide ranges are disjoint with >= 20 certain reads each; a count, not a test. A copy's range is
  that of its deepest bunch under the key.
- OR6737: 3,652 of 148,190 psi rows get a key, in 307 families (only the representative's exons are in the columns);
  2,416 keys, 505 seen in >= 2 copies, 212 of them differ (55 families). By the gap between the ranges: > 0.1 75 keys,
  > 0.2 44, > 0.5 16 (side script, not in the tool). Example MCL123: one junction at 0.01 in the copy at 117 Mb and 0.48
  in the copy at 125 Mb of NC_073224.2.
- Not done: any check of the differing junctions against the annotation or in a second sample.
