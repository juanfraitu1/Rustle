# Multi-copy gene families: what the reads can decide (session of 2026-11-18 / 19, on PSU Roar)

Work done on PSU while Fritz was in maintenance. Goal: `GOAL_multicopy_families_2026-11-19.md`. Registered rows: 1820, 1821,
1822 (`NEGATIVE_RESULTS_REGISTER.md`). Development record with every number: `MULTICOPY_DEV_LEDGER_2026-11-19.md` (a copy of the PSU workspace
`goal_multicopy_dev/LEDGER.md`, steps 1-11).

## The question
The advisor's position: the provenance of multi-mapping reads is the main problem, and an assembler that borrows information
across the copies of a family should gain several-fold. The user's view after this session: no several-fold gain is in the
reads (even every secondary alignment adds few transcripts); the better test is whether the members of known families are
complete, and why they are not.

## What was established

### 1. Reads that cannot be placed do not reveal hidden copies (row 1820)
Reads that fit several of HG002's own copies equally (SEVERAL COPIES ties) land on 116 copies on CHM13: 82 protein-coding
paralogs, 11 pseudogenes. Adding their other alignments to the assembly gives the annotated structure to 1 copy that
primaries do not reach (bar 20), 2 of 20 added chains match the annotation, and 82 unmatched chains appear at copies with
their own reads. **C1 FAILS.** Side table (A119b testis arrays, no truth): exact ties spread structure inside near-identical
arrays (TSPY coding units 11 -> 25 of 38) at a cost in unmatched chains; pseudogene units gain little.

### 2. A measure for copy claims
`bench/famwise_eval/copyset_score.py`: every transcript is a claim about which copies its reads come from. A claim is right
(EQUAL) when its copies are the copies its supporting reads fit equally on the individual's own haplotypes; NARROW = it
claims one copy where the reads fit several; WIDE = it claims copies the reads do not fit. `--match contained` lets a read
fragment support the transcript that holds it. Precision, a molecule-weighted share, resolved units (true (copy set, chain)
units held by one claim with exactly their copies), and new / lost chains against a base.
Finding: the per-copy output's main error is NARROW (a transcript at one copy whose reads fit several).

### 3. A per-family resolution call (row 1821, judged on the last unseen HG002 quarter)
`resolution_call.py`: keep a family transcript of the `famwise` layer as one copy set only when >= half of its own supporting
reads tie across its copies (within 10 ms); otherwise report the copies one by one.
| output | precision | resolved units |
| primaries, per copy | 0.888 | 1,019 |
| primaries -> famwise | 0.909 | 1,034 |
| primaries -> famwise -> call | 0.909 | 1,044 |
| StringTie | 0.861 | 919 |
R1 (call over the uncalled layer, +1%) **FAILS narrowly**; R2 (over the per-copy output) PASSES (+2.5%, precision up);
R3 (copy sets hold 0.94 of their reads) PASSES. Grouping copies is what helps; the call adds little over the layer.

### 4. Secondary alignments as support only, never as seeds (built, judged in row 1821)
`RUSTLE_SEED_SPLICED=support` (opt-in; `src/family.rs`; refused off the streaming `--assemble-only` path; default output
byte-identical: SIRV md5 and a whole HG002 quarter): a spliced chain seeds only when a primary record carries it; secondaries
still add read support. On the judging quarter, against primaries alone:
| | new right | new NARROW | new wrong copy | lost (judged) |
| default good pool | 51 | 298 | 60 | 30 |
| support only | 23 | 88 | **15** | 26 |
About half the right new chains are kept and wrong-copy chains fall 4x.

### 5. Tandem arrays (simulation, `famsim`: 10 TSPY3 units)
- Identical units: minimap2 scatters reads evenly (10% on their own unit) whatever the true expression; every assembly puts
  the transcript on every unit, silent ones included; the spread over units cannot be recovered. The call makes one copy set
  over all units: the least wrong output (against the literal source unit nothing is right).
- From 0.2% exonic divergence primaries place every read; the default pool still puts chains on all silent units;
  support-only removes them (0.1%: 7 -> 3; 0.2%: 7 -> 0).
- Real testis arrays have no origin truth here (A119b has no assembled haplotypes). On A119b the call groups MAGEA, AMY,
  TTTY, RBMY, GAGE, CT45A / CT47A into copy sets of up to 8 units (descriptive).

### 6. Borrowing across paralogs
- As structure (`borrow.py`, labelled `borrowed_structure`): paralog transcripts projected onto a copy where read evidence
  matches within 5 bp. Development: 5 more copy-set units with structure; read as expression only 0.17 right, so the class
  must stay labelled. Judged on KB3781 (row 1822, B1): **FAILS**, 2 -> 2 of the 2 copy-set units KB3781 has.
- As shared isoforms (the user's reading: isoforms several family members produce; `family_isoform.py`): 8% of isoforms
  match an isoform assembled at another place; one-read chains promoted because a paralog has the same isoform are 17 of 19
  right on the copy measure (KB3781, descriptive). Promising, too few to judge.
- Ceiling on development: of 102 units no arm recovers, 4 have their junctions emitted at a paralog.

### 7. Are the members of known families complete? (the user's test; StringTie as the reference tool)
`member_completeness.py`: members = 1,243 genes of 365 truth families (Compara recent paralogues); a member is found when a
transcript maps to it; each miss is classed by the reads.
| sample | members with own reads | StringTie | Rustle primaries | Rustle defaults | missed by all: no reads / own reads / tied only / truly theirs but placed elsewhere |
| HG002 quarter (truth) | 338 | 290 | 268 | 272 | 890 / 42 / 4 / **0** |
| HG002 all reads | 409 | 348 | 312 | 332 | 626 / 44 / <= 98 / - |
| A119b testis | 835 | 721 | 646 | 702 | 289 / 70 / <= 66 / - |
- Missing members are mostly silent: no read reaches them.
- No member is missing because its reads were placed on another copy (0 on the truth quarter).
- Seeding tied secondaries at the members nobody finds restores 3 (HG002) and 7 (A119b); the extra members the default pool
  finds through secondaries are mostly shadows (11 of 15 on truth have no read from them).
- **The actionable gap is the assembler's floor**: StringTie finds 24-75 more expressed members than Rustle. Of the 24 on the
  truth quarter, 19 get no Rustle transcript (2-36 own reads, low-support chains: row 1536's class), 5 are fused into a
  neighbour. `RUSTLE_READ_ISOFORM_K=2` closed 9 of 22 on development but found no extra member on KB3781 (row 1822 K1 FAILS;
  KB3781 leaves little room: primaries already find 66 of 69).

## In one paragraph
The reads do not hide a several-fold gain. Reads that fit several copies are best reported as a copy set, which raises
precision and resolution a few percent over per-copy output; secondaries should support chains that primaries carry, not
seed new ones (4x fewer wrong-copy chains); identical tandem units can only be reported as one set; borrowing structure
across paralogs adds little and must stay labelled. Family members are missing because they are silent or because the
assembler's read floor drops low-support chains that StringTie keeps, not because their reads went to another copy.

## Code (all pushed; opt-in, defaults unchanged)
- `src/family.rs`, `src/copy_assign_bin.rs`: `RUSTLE_SEED_SPLICED=support` (+ unit test; driver comment).
- `bench/famwise_eval/`: `undecidable.py` (row 1820), `copyset_score.py`, `resolution_call.py`, `support_filter.py`,
  `array_score.py`, `sim_truth.py`, `borrow.py`, `family_isoform.py`, `member_completeness.py`, `member_bed.py` (+ tests).

## Data spent
HG002 on CHM13, hash quarters: 0-24 seen (rows 1776-1791), 25-49 row 1796, 50-74 row 1820, 75-99 row 1821. No unseen HG002
quarter remains. KB3781: development in rows 1796 / 1820, judged in row 1822.

## Open (see PARKED.md item 3)
1. The member gap: why Rustle drops low-support members that StringTie keeps (row 1536 class), judged on a sample with room.
2. Shared isoforms as a promotion rule for one-read chains: judge on a sample with many one-read chains in families (needs
   new bars; HG002 only by reuse; KB3781 with the Fritz gene map, `DENOVO_GGO_ENS_MAP`, gives a larger family truth).
3. Tandem arrays on real data: needs a testis Iso-Seq sample from an individual with assembled haplotypes.
4. Whether `RUSTLE_SEED_SPLICED=support` becomes a default: the user's call.

## Where the files are
PSU workspace (`/storage/group/kdm16/default/jxi21/fritz/ws`): `r1820_undecidable_copies`, `r1821_resolution_call`,
`r1822_borrow_isoformk`, `goal_multicopy_dev` (LEDGER.md, arrays_sim, members_full). To Fritz:
`rsync -a psu:makova_group/fritz/ws/{r1820_undecidable_copies,r1821_resolution_call,r1822_borrow_isoformk,goal_multicopy_dev} ~/ws/`
(the snapshots and `goal_multicopy_dev/target` hold binaries; exclude them if space is short).

## Added later the same day
- Correlation clustering (Strainy port) against the HG002 truth: `CORRELATION_CLUSTERING_2026-11-19.md` (register 1823).
- The oracle ceiling of correcting wrongly chosen primaries: +0.2% resolved units, no family member
  (`PRIMARY_CORRECTION_CEILING_2026-11-19.md`, register 1824).
- Secondary-only places, real or shadow: a few real (mostly pseudogenes), each dominated by another copy's reads
  (`SHADOW_LOCI_2026-11-19.md`, register 1825).
- Tandem arrays: how StringTie and Rustle assemble them, read artifacts, member size by folding reads onto the period
  (`TANDEM_ARRAYS_2026-11-19.md`, register 1826).
