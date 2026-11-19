# PREREG r1820: do undecidable reads depict the copies they cannot be placed on? (2026-11-19)

Written on PSU (Roar) while Fritz is in maintenance; workspace `~/ws/r1820_undecidable_copies` (claimed on PSU after the
register tail on origin/main was 1819; Fritz sessions told). Readings C1, C1b, C2, C3 and the class thresholds accepted
by the user on 2026-11-19 ("yes, looks good"), before any read of either sample was looked at for this row. KB3781
development may resize the `shadow` threshold as written below; any other change after this point is an amendment,
dated, in a section of its own.

## Why
Advisor, as relayed on 2026-11-18 / 19: secondary alignments of reads that cannot be placed depict pseudogenes and other
copies, and should be used when forming transcripts. A gene's copies, silent or not, could then be found from the
multi-mapping reads of the expressed copy. This row tests that claim on its own terms: as a statement about the
STRUCTURE of a copy, not about its expression.

What bears on it (do not re-run):
- Placement other than minimap2's primary is closed (rows 1642-1651, 1696, 1777-1785); near ties are 5.7% of family
  reads on HG002 / CHM13 (1790), of which 44% in the 14 tie-heavy families fit several HG002 copies equally (1791).
- Even exact provenance changes RefSeq-exact chains in tied regions by x1.0-1.25 (1794, 1798, 1809).
- Loci held by secondaries alone are mostly shadows when read as expression: KB3781 strict truth 0 real / 155 shadow
  (1671); the shadow flag failed held out (1672); 1673 = the real ones are rDNA scaffold copies. Wider seeding bands add
  chains mostly at shadow loci (1767). Unspliced secondaries lift single-exon pseudogene copies over the floor (1518,
  1520). Read-free projection of family models onto members: 1519 (failed, bench prototype). This row differs from
  1519's `tied` tier (a family-internal secondary carrying the exact chain): 1519 judged projected models as
  expression (SQANTI against CAT, artifact share) on an A119b panel; here the reads' records feed the assembler, the
  unit is a copy reached ONLY by undecidable reads, the judge is structure, and projection is the control (PROJ, C1b).
- Annotated structure is not origin (1767), and the RefSeq-exact measure credits false copies (1796 / 1798): CAT and
  Liftoff chains at paralogs and unprocessed pseudogenes are largely the parent's structure lifted onto the copy, so a
  parent read's secondary lands that same chain almost by construction. C1 alone cannot say the reads add anything;
  C1b compares against a read-free projection.
- The `famwise` layer already reports a family transcript as a set of copies with a status; on held-out HG002 its
  `tied_only` claims are 0.84 shared or false as expression claims (B2 of 1796).

So: if undecidable records are counted as expression they are known to be wrong. The open question is whether they put
the right STRUCTURE at the right copy, so that a separate output class ("structure at a copy, no expression claim") is
worth having.

## Questions
- **Q1 (where).** Where do the non-primary records of undecidable reads land: on which biotypes, in which families?
- **Q2 (what kind of copy).** For each landing copy, do the reads that CAN be placed say it is expressed (own evidence),
  silent although distinguishable (shadow), or indistinguishable (cannot be judged)?
- **Q3 (the claim).** Do undecidable records, added to the primary alignments, make the assembler emit the annotated
  intron chain of copies that the primaries alone do not reach; how much of that do today's defaults already emit; and
  how many chains do they add that match nothing annotated at the copy?

## Inputs
- **HG002 Kinnex on CHM13** (`~/reads/HSA/HG002_kinnex.t2t.bam`, on PSU). Truth: `bench/o2_truth/neutral_truth.py` with
  `NT_REPAIR=1` against HG002 T2T v1.1 mat / pat (`~/ws/r1776_hg002_room/hg002v1.1.{mat,pat}.fasta.gz`; indices rebuilt on PSU with
  the same minimap2 2.30-r1287 and command as on Fritz: `minimap2 -x splice:hq -d`), then `bench/o2_truth/tie_split.py` (row 1791).
- **KB3781 on its own reference** (`~/reads/GGO/GGO_KB3781_fib_mm.bam`, `~/genomes/GGO/GGO.fa`). Truth: KB3781 mat / pat
  (`~/ws/o2contested/haps/KB3781.{mat,pat}.mmi`, copied from Fritz, built with
  minimap2 2.30-r1287).
- **Biotypes.** CHM13: `$DENOVO_CAT_TSV` (gene level, on PSU). Gorilla: `$DENOVO_GGO_GFF` (NCBI) and `~/ws/r1735_ggo_cat/cat.gff3`.
- **Annotated intron chains.** CHM13: `~/ws/r1797_famsim_mimic/annot/chm13v2.0_CAT_Liftoff.slim.gff3.gz` (234,903
  transcripts, 16,000 pseudogene transcripts by `transcript_biotype`). Gorilla: `$DENOVO_GGO_GFF` and
  `~/ws/r1735_ggo_cat/cat.gff3`.
- **Family membership.** The truth families of `bench/denovo_eval/famtruth.py`; per-read family as in rows 1790 / 1791.

## Samples and the judging set
- KB3781 develops: every script is written and debugged there, then frozen (hash, commit) before HG002 is opened.
- HG002 judges, on a quarter nobody has looked at: `NT_SELECT_PCT=25 NT_SELECT_OFFSET=50` (hash residues 50-74;
  0-24 are seen, 25-49 were spent by row 1796). Checked 2026-11-19 in the repository and the workspace on PSU: no run has
  selected residues >= 50. The parked cross-copy test (rows 1802-1804, not run) has a ruling that its C may take HG002
  reads of every residue at its own call loci; if it runs after this row, those loci are no longer unseen there, and
  that is said in its record.
- Rule for both: one alignment, no re-alignment; scores re-read from the BAM; `ms` first, AS beside.

## Definitions, fixed now
- **Near-tied read:** its two best places are within 10 `ms` (records on one contig and strand with overlapping spans are
  one place; the truth's merge).
- **Undecidable read (U):** near-tied, truth status `tie`, tie split class SEVERAL COPIES. ONE PLACE ties, `outside`
  and `no near-perfect hit` are counted and reported apart, never pooled with U.
- **Decidable read (Dd):** truth status `decidable` (one origin copy), near-tied or not.
- **U records**, two versions:
  - oracle `U*`: a U read's records on the loci of its truth origin set, other than its primary's place;
  - practical `Ums`: the non-primary records of every read whose best two places are tied EXACTLY on `ms` (no truth used).
    Its overlap with `U*` is printed.
- **Landing copy:** a CHM13 (or GGO) place holding at least one `U*` record, merged as the truth merges loci; annotated
  by the gene it overlaps on the same strand (largest overlap), biotype grouped as protein_coding / unprocessed
  pseudogene (unprocessed, transcribed_unprocessed, unitary, transcribed_unitary, polymorphic) / processed pseudogene
  (processed, transcribed_processed, translated_processed) / lncRNA / other / none.
- **Class of a landing copy, by decidable reads** (Q2):
  - `own`: at least 2 Dd reads have their origin at the copy;
  - `shadow`: 0 Dd reads have their origin there, AND at least 5 Dd reads with origin elsewhere in the family have a
    record there at ANY score (`-N 50 -p 0.1` keeps them; the copy is distinguishable from where those reads came, and
    none came from it). A record-within-10-`ms` version is printed beside (it is near-tied by definition and expected
    nearly empty: about 600 decidable near ties per HG002 quarter, row 1785). The threshold 5 is checked on KB3781 class
    counts before freezing; if `shadow` holds fewer than 30 copies there, the threshold is resized on KB3781 only and
    the new value written here before HG002 is opened;
  - `unjudgeable`: neither (no read could separate this copy). Copies with exactly 1 Dd origin read are `own_1`, printed
    apart.
- **Depicted:** a landing copy carries an emitted multi-exon transcript whose intron chain equals an annotated chain of
  the gene there (exact), or equals a contiguous run of >= 2 of its introns (`depicted_part`, printed beside). For
  processed pseudogenes (no introns): a single-exon model covering >= 80% of the annotated gene's exonic span, printed,
  no bar.
- **Unmatched chain:** an emitted multi-exon chain on a landing copy whose introns match no annotated chain at the copy
  in either sense above.

## Arms (one selected BAM per sample, built by a new `bench/famwise_eval/undecidable_bam.py` on the pattern of
`origin_bam.py`)
- **P**: primary alignments only (`assemble --seed-pool primary`).
- **P+U\***: P plus the `U*` records, kept as secondaries (`--seed-pool all`).
- **P+Ums**: P plus the `Ums` records (`--seed-pool all`).
- **D**: the selected BAM at the defaults of `all` and every stage but `assemble` alone: `--seed-pool good` (seeding band
  AS >= 0.98 x best). `assemble` alone seeds with primaries since row 1812, which would make D equal P.
- **PROJ** (read-free control, no assembly): P's multi-exon chains on the primary's place of each U read projected onto
  that read's landing copies (minimap2 `splice:hq`, cov >= 0.9, id >= 0.85, as row 1519); a projected chain is
  "emitted" at the copy for the depicted test.
- Rustle read on `PREFIX.gtf` (bars) and `PREFIX.final.gtf` (beside). StringTie `-L` on P and P+U\* as the trusted
  baseline, beside, no bar. Regions = merged spans of every landing copy and of its U reads' primaries (no genome-wide
  assembly).
- **Panel first** (working rules): the arms run first on 10 families of KB3781 in hash order, base rerun beside; arm wall
  <= 2 x P's wall, whole row wall limit 1 h 30 on one `basic` node; hitting the limit = FAIL.

## Readings (accepted 2026-11-19)
Read on held-out HG002, landing copies of classes `shadow` + `unjudgeable` (the copies the claim is about); KB3781 beside,
no bar. Intervals Wilson 95%.
- **C1 (Q3, does it depict).** Added depicted copies = depicted under P+U\* and not under P.
  The advisor's claim HOLDS when added depicted copies >= 20 AND among the chains P+U\* adds on these copies the depicted
  share (exact) has its lower bound >= 0.5. It FAILS when added depicted copies < 20, or the depicted share's upper bound
  < 0.5. Otherwise: between, reported.
- **C1b (do the reads add anything over projection).** On the same copies: added depicted copies under P+U\* that PROJ
  does not depict, and P+U\*'s depicted share beside PROJ's. "The reads add structure" is claimed only when C1 holds
  AND the copies depicted by P+U\* but not by PROJ number >= 10. C1 holding with C1b not met reads: the copies are
  depicted, but no better than by projecting the parent.
- **C2 (is it new against the defaults).** Of C1's added depicted copies, the share D does not already depict. Printed
  with its interval; "new over defaults" is claimed only when that share's lower bound >= 0.3.
- **C3 (is it reachable without truth).** P+Ums reaches >= 0.5 of C1's added depicted copies, and its depicted share on
  the same copies is within 0.1 of P+U\*'s. Only read when C1 holds.
- **Cost, printed, no bar:** unmatched chains added by P+U\* and by P+Ums on `shadow` copies and on `own` copies (a
  chain added at an `own` copy changes an expressed gene's output).
- **Descriptive, no bar:** Q1 table (U reads, U\* records and landing copies by biotype and family, `ms` and AS
  versions); Q2 table (classes x biotypes); the same for ONE PLACE ties; processed-pseudogene single-exon coverage;
  landing copies with biotype `none` (unannotated copies; candidate O3 input, not scored).

If C1 holds, what it makes possible (not decided here, the user's call): an opt-in `famwise` output class, e.g.
`fw_status "depicted_copy"`, for a chain at a copy reached only by undecidable reads, never counted as expression.

## Prediction (mine, before any count)
- Q1: most U\* records land on protein_coding paralogs of the 14 tie-heavy families; pseudogene biotypes 10-30% of
  landing copies.
- Q2: `unjudgeable` is the largest class (SEVERAL COPIES reads come from copies that are near identical), `shadow`
  second.
- C1: holds on count (>= 20 added depicted copies) with depicted share 0.5-0.7.
- C2: fails to show "new": D already depicts >= 70% of them (its seeding band admits the same tied records; row 1767).

## Time to answer (estimate)
Truth + tie split per sample about 1 h (row 1644's recipe: 56 min for HG002); BAM selection and the four arms on the
landing regions: minutes each (KB3781 whole-sample `all` 477-482 s); PROJ and scorer under 30 min. KB3781 development
about half a day of work; the HG002 judging run about 2 h wall from submission.

## Open before registering
1. Fetch from Fritz: HG002 T2T v1.1 mat / pat (FASTA or `.mmi`, with the minimap2 version and command that built
   them; public HPRC FASTAs are the fallback, rebuilt with the same command), KB3781 mat / pat, the CHM13
   transcript-level annotation behind "RefSeq-exact" in rows 1794-1809 (the file `sim_chains.py` / `gtf_origins.py` were
   given), `r1797_famsim_mimic/annot/`.
2. Verify the residue range 50-74 is unseen.
3. The user accepts or changes C1, C1b, C2, C3 and the class thresholds (2 Dd origin reads for `own`, 5 for `shadow`).
4. Write `undecidable_bam.py` and the scorer with tests on KB3781; freeze; then claim the row and commit this file as
   the pre-registration.

## Amendment 1 (2026-11-19, before any count of this row)
KB3781 holds too few undecidable reads to size anything: the selection of rows 1775 / 1796 has 186 molecules in exact
`ms` tie classes of 100,810 (row 1796's `eq.log`; seen before this row, not a count of it). Development therefore uses
two sets: KB3781 (mechanics: every script runs end to end there) and the SEEN HG002 quarter (hash residues 0-24; selected
BAM `r1796_famwise_origin/hg_seen/sel.bam`, repaired truth `r1785_truth_repair/hg/reads.tsv`, tie split of row 1791).
The `shadow` threshold is checked and, if needed, resized on the seen quarter, by the rule of the definitions (>= 30
`shadow` copies), and frozen before residues 50-74 are selected. Bars are still read on residues 50-74 only.
Extra inputs to fetch from Fritz: `r1790_family_ties/` (`family_ties.py`), `r1791_tie_split/` (`tie_split.tsv`).

## Amendment 2 (2026-11-19, on the seen HG002 quarter, before residues 50-74 are selected)
The rule of the definitions applied: at 5, `shadow` held 28 copies (< 30), so the threshold was resized on the seen
quarter: `shadow` = 0 Dd origin reads and >= 1 Dd read with origin elsewhere holding a record at the copy (any score).
Seen-quarter counts at thresholds 1 / 2 / 3-5: 30 / 29 / 28; the copies split into 18 reached by no Dd read and 28
reached by >= 5, so the threshold moves 2 copies. Chosen by the user ("yes, lets lower it"). Frozen at 1.

## Side table: tandem arrays on A119b (added 2026-11-19 on the user's word, descriptive, no truth, no bar)
Why: the user asked about TSPY. HG002 is a cell line with 14 primary reads on the whole TSPY array (A119b testis: 786),
and the truth keeps at most 12 loci per read, so arrays cannot be judged by origin here.
- Arrays: CAT / Liftoff genes (`chm13v2.0_CAT_Liftoff.slim.gff3.gz`) whose names share a stem (name minus a trailing
  number with an optional letter and `P`: TSPY10, TSPY6P -> TSPY), groups protein_coding / pseudogenes / lncRNA, on one
  contig with neighbours <= 200 kb apart, >= 5 units (`undecidable.py arrays`; 227 arrays, 2,264 units on CHM13).
- Input: `~/reads/HSA/A119b.t2t.bam` sliced at the arrays (+-5 kb). P = primaries; PUms = P + the slice's non-primary
  records at the best ms of reads whose two best places in the slice tie exactly (a read's records elsewhere are not
  in the slice); D = the slice at `--seed-pool good` (its best-AS table is the slice's); PROJ = P's multi-exon chains
  projected onto the other units of the same array. No U* arm (no truth).
- Printed per unit and per family x biotype group: units with annotated chains, units depicted (exact) per arm,
  unmatched chains per arm. The pseudogene units of TSPY (TSPY6P, 7P, 9P...) are the cases the user asked about.

## Freeze (2026-11-19, before residues 50-74 are selected)
Snapshot `~/ws/r1820_undecidable_copies/snap/` (tools, bench, release binaries copied); hashes in `snap/FROZEN.sha256`.

## Result (2026-11-19, job 57494600 on PSU, 37 min; workspace `judge/` and `arrays/`)
Integrity: the text above this section, minus its final blank line, hashes to the frozen value in `snap/FROZEN.sha256`
(068bd12a..., written 22:17 EST, before job 57494600 started). Correction to that text: "Fritz sessions told" was not
true when written; the Fritz session is to be told by the user before row 1820 is pushed.

**C1 FAILS: adding the undecidable reads' records gives the annotated chain to almost no copy that the primaries do not
already reach. C1b, C2 not met; C3 not read.** The share bar alone decides it (upper bound 0.30 < 0.5); the count bar
asked for 20 of 35 claim copies.
- Of the 35 claim copies, 30 carry an annotated multi-exon chain; the exact chain is emitted at 5 by P, 6 by P+U*, 7 by
  P+Ums, 8 by D, 6 by PROJ (contained part: 6 / 9 / 10 / 12 / 9). StringTie `-L` (beside, no bar): 6 by P, 6 by P+U*;
  it adds no chain from the secondaries.
- Judging quarter (HG002 residues 50-74): 136,262 reads selected; truth kinds one 18,321, several (U) 2,985, one place
  1,227, tie unsplit 589, outside 1,459. U* records 4,021 (Ums 3,904, both 2,589).
- Q1: 116 landing copies: protein_coding 82, lncRNA 13, unprocessed pseudogene 10, processed pseudogene 1, none 8,
  other 2. U reads land on protein-coding paralogs (65 of the 116 copies have own evidence from decidable reads), rarely
  on pseudogenes (my prediction 10-30% pseudogenes: 9.5%).
- Q2: own 65, own_1 16, shadow 23 (6 under the within-10-ms rule), unjudgeable 12 -> 35 claim copies (prediction:
  unjudgeable largest; wrong, shadow is).
- C1 (`PREFIX.gtf`): added depicted copies 1 (bar >= 20); depicted share of the chains P+U* adds 2/20 = 0.10
  [0.03, 0.30]. C1b: 0 of them missed by projection. C2: D already depicts it (0/1). C3: P+Ums reaches it (1/1).
- Cost: P+U* adds 82 unmatched chains on own copies and 6 on shadow copies; P+Ums 19 and 5.
- Beside, `PREFIX.final.gtf`: same verdicts (added depicted 1; share 2/13; unmatched on own copies 16).
- Development sets (no bar) agreed: seen quarter 2 added depicted copies, share 0.13; KB3781 0.

Side table, A119b tandem arrays (descriptive; `arrays/arrays_summary.tsv`, `arrays_depicted.tsv`): 227 arrays, 1.16 M
reads in the slices. Units with an annotated chain, depicted by P / P+Ums / D / PROJ:
- all protein_coding 1,166 units: 570 / 593 / 630 / 185; all pseudogene units 364: 27 / 32 / 47 / 15; lncRNA 155:
  36 / 36 / 37 / 0. Unmatched chains P / P+Ums / D: protein_coding 4,528 / 4,666 / 6,184; pseudogenes 460 / 492 / 949.
- TSPY coding units 38: 11 / 25 / 38 / 38 (unmatched 4 / 58 / 138 / 165); TSPY pseudogene units 14: 3 / 4 / 6 / 6
  (unmatched 6 / 16 / 48 / 36). FAM90A pseudogene units 35: 0 / 4 / 16 / 0 (unmatched 0 / 8 / 59). USP17L 14 coding and
  REXO1L 52 processed-pseudogene units: none depicted by any arm.
- Reading (mine): inside near-identical arrays, tied records do spread the expressed unit's structure onto other units
  (TSPY coding 11 -> 25 with exact ties, 38 at the default band), at a cost of 4 -> 58 / 138 unmatched chains; pseudogene
  units gain little (TSPY 3 -> 4 / 6). No truth here: a unit "depicted" in an array of identical units says the
  structure fits it, not that it is expressed.

Registered items not delivered, or changed:
- Q1 per family: not done (the family table of row 1790 was not fetched; `classes.tsv` carries each copy's CAT gene name
  instead). Q1 in `ms` and AS versions: not split (the U set is the truth's, no score enters it).
- The ONE PLACE table: not made (1,227 such reads counted in `inputs.tsv`, not classed).
- Processed-pseudogene single-exon line: computed in `depicted.tsv` (`*_single`), 1 such copy, depicted by no arm.
- Panel first on 10 KB3781 families: replaced by the whole KB3781 selection (100,810 molecules; arms 30-60 s each).
- StringTie on P and P+U*: run after the job, on the login node, same BAMs (`judge/stringtie/`).
