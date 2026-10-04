# O3 candidate copies in the pipeline (`o3_candidates`) — design

Date 2026-10-02. Branch `machine2/soto-evidence` (== `main`). Status: implemented 2026-10-02; the pre-registered acceptance
(Amendment 12) FAILED, so the driver's `candidates` stage is OPT-IN (ruling R14, §9b; `docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`).

## 1. Goal

Make the reference-absent-copy chain validated in `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md` (Amendments 7-11) a stage
of `tools/rustle_pipeline.sh` between O1 (`families`) and O2 (`assign`), **without IsoCon**: per family, reads -> read clusters -> one
consensus per cluster -> flag / link / merge -> candidate copies, each represented by the exon union of its component -> O2 run over the
reference plus the candidates. Only minimap2 stays external, as today.

## 2. What the measurements fix (not open for redesign)

| decision | measured basis |
|---|---|
| Candidates come from the family's **read net**, not from piles at loci | Amendment 11: 31/53 deleted copies have no read whose primary lands on a surviving copy (5,257 reads unmapped); `missing_copy_flag` reached 14.7% of the deleted copies' reads vs IsoCon's 74% |
| Consensus = the **reads' own** consensus, not a reference patched at sites | Amendment 11: 6/14 fired loci gave a host-backbone hybrid that labels as the survivor |
| Link rule: output within **delta = 0.00958** of a reference locus = that locus' allele; delta = p99 allelic divergence at single-copy genes | Amendment 7 (WORKS, HELP); on the Y, delta_Y from single-copy Y genes (YAG prereg) |
| Merge rule: components of pairwise `de` <= delta over >= 50% of the shorter | Amendment 8 (M1 WORKS = truth ceiling; identical at delta/2 and 2 x delta) |
| Flag floor: a candidate needs **>= 2 clusters** | Amendment 9 post hoc (LR 8.2 vs 2.75), pre-registered and confirmed in Amendment 10 (D1-D3) and the YAG runs |
| Representative = **exon union** of the component | `bench/rna_allele/rep_choice.py`: longest member keeps 73% of the copy's reads, most-supported 71%, union 100% |
| The per-locus screen stays | `missing_copy_flag` keeps verdicts / editing / Ig / contamination screens and the DNA-depth expectation; it is not the consensus source |

## 3. Scope

In: the `o3_candidates` binary; the `candidates` driver stage with augmentation and patch realignment; `assign` reading the default O1
output; the chromosome-aware overlap fix in `copy_assign`; `flag` corroboration column. Out: POA-graph consensus (v2, §10), the
hierarchy work, YAG-specific delta (the stage takes `--delta`; the Y value is the caller's business).

## 4. Pipeline placement and data flow

```
assemble -> families (O1: P.fam.copies.{tsv,fa,regions}) -> candidates (new; opt-in, R14) -> assign (O2, augmented) -> flag (O3 screen)
```

```
P.fam.copies.tsv + P.fam.copies.fa + BAM + FASTA + splice index
   | o3_candidates
   v
P.cand.candidates.tsv   one row per candidate copy (family, id, n_clusters, n_reads, flag, union length, nearest locus, d)
P.cand.contigs.fa       one union sequence per candidate: >cand_<family>_<k>
P.cand.clusters.tsv     one row per cluster (candidate, cluster id, n_reads, consensus length, linked_to, d)
P.cand.reads.tsv        read -> cluster (diagnostics; not read by O2)
P.cand.nets.fa          the net reads of families WITH a flagged candidate (input of the patch realignment)
   | driver: augmentation + patch realignment
   v
P.aug.fa = FASTA + P.cand.contigs.fa          P.aug.copies.{tsv,fa} = P.fam.copies.* + one row per flagged candidate
P.aug.regions.txt = copy hulls +- 5 kb, and cand_*:0-len
P.aug.bam = P.cand.nets.fa realigned to P.aug.fa (minimap2 splice:hq, the pipeline's own flags)
   | copy_assign --families P.aug.copies.tsv --copies-fa P.aug.copies.fa --regions P.aug.regions.txt --bam P.aug.bam   (candidate families)
   | copy_assign --families P.fam.copies.tsv ... --bam BAM                                                               (other families)
   v
P.assign.* (tables concatenated; families are disjoint between the two runs)
```

## 5. `o3_candidates` (new binary, `src/bin/o3_candidates.rs` + `src/rustle/vg_family/o3_candidates.rs`)

CLI: `o3_candidates --bam B --fasta G --copies P.fam.copies.tsv --copies-fa P.fam.copies.fa --index G.splice.mmi --out P.cand
[--delta 0.00958] [--max-reads 1000] [--min-cluster 3] [--min-clusters 2] [--threads 4] [--families F1,F2]`.
Environment: `RUSTLE_MINIMAP2` (binary path, as `mcl_families`), `RUSTLE_CACHE_DIR` (as the other binaries). Deterministic: seeded
sampling (seed 1), stable sort orders, no threads in anything that affects output (threads only inside minimap2).

### 5.1 Net (per family)

1. Copies = rows of `--copies` with the family id; intervals = `locus_start..locus_end` when present else `start..end`.
2. BAM pass A (indexed, per interval): every record overlapping a copy, primary or secondary (`aligned_read_from_record`,
   `denovo_assemble.rs:1327`), supplementary excluded. Primary records give the read's sequence (orientation = as sequenced: reverse
   complement when the record is reverse). Secondary records only name the read.
3. BAM pass B (one sequential pass over the whole BAM, once for all families): (a) the sequence of every read named only by a secondary
   record; (b) the ATTRIBUTION SET of §5.2 (prereg Amendments 13b / 13c): every unmapped record, and every read in no net of this run
   whose PRIMARY record is poorly placed (`de > 0.02` or MAPQ 0; ruling R18: "no record on a family copy" means in no net of THIS run,
   pass A's scope under `--families`), each >= 300 bp (ruling R19: the floor holds for both classes), streamed as sequenced to a FASTA
   (decoded once, never held in memory). (2026-10-02: unmapped records only, attributed by a k-mer index; retired by Amendment 13.)
4. Cap: if a family's net exceeds `--max-reads`, sample that many (seed 1, after sorting names).
5. Output: `P.cand.nets.fa` holds the WHOLE net (before the cap of 4.) of each family that ends with a flagged candidate (§5.6), each read
   once (ruling R9: the first family in `--copies` order keeps a read two nets share) — the input of the patch realignment (§7). The
   `cand` cache entry stores the products as written: no other family's net is kept.

### 5.2 Attribution of unmapped and poorly placed reads (prereg Amendment 13b, ruling R16)

The attribution set of §5.1.3b is aligned once, `MM2_ATTRIB = -x map-ont -c -N 5 -p 0.5`, against the run's net reads (every read pass A
put in a net, as `>{family}|{read}`) and every family's `--copies-fa` records. A read joins the family of its best hit (most matches;
the first on a tie) iff that hit covers >= 50% of the READ (`(qe - qs) / qlen`) and its `de` <= 0.20: the family definition's own edge
rule (identity >= 0.80 over >= 50%) in read space. A best hit that fails either joins the read to nothing (a lesser hit never stands
in); a read given to a family outside this run (`--families`) joins no net. Joiners enter their net before the cap (§5.1.4). The
chain's genome check (§5.6.2, identity x coverage 0.999) remains the guard against reads of foreign genes pulled in this way. Under
`--families` the attribution set and the targets are the run's (R18), so batched runs may attribute one read in two batches.
(2026-10-02's rule — canonical 31-mers of `--copies-fa`, >= 30% of a read's k-mers and >= 2 x the runner-up — attributed 0 of 5,312
unmapped reads on the held-out, Amendment 12; Amendment 13 retired it. The A13 acceptance: 557 unmapped reads joined a net, 555 of them
the right family by label; 1,171 poorly placed reads, 927 the right family; 5 batches.)

### 5.3 Read clustering (greedy template clustering, per family)

Reads sorted by length descending. For each read: candidate clusters = those whose template shares >= 50% of the read's 31-mer sketch
(minimizers, w = 5); the read is aligned to each candidate's template with the poasta global aligner used by `family_graph.rs:23-62`
(`PoastaAligner`, `AlignmentType::Global`, `GapAffine` costs as there); it joins the first cluster (in creation order) whose alignment
covers >= 50% of the shorter sequence with gap-compressed divergence `de` <= delta; otherwise it founds a new cluster with itself as
template. `de` = (mismatches + gap openings) / aligned columns, computed from the alignment's operations — the same quantity minimap2's
`de` tag reports, so the thresholds carry over from Amendment 8. The read's alignment to its template is kept for §5.4.
As implemented (the minimap2 engine of §9b): the net's reads are aligned all-vs-all (`MM2_AVA = -x asm20 -c --cs --dual=no -N 100
-p 0.1 --secondary=yes`) and clustered by union-find over the pairs whose best hit (most matches) covers >= 50% of the shorter read
with `de` <= delta (`cluster_reads`: Amendment 8's merge rule applied to reads); the same all-vs-all is what §5.4's template rule reads.

### 5.4 Consensus per cluster (template and vote)

Clusters with < `--min-cluster` reads are dropped (IsoCon's floor of 3).

**Template** (prereg Amendment 13 as corrected by 13d / 13e, rulings R20 / R21; it replaces "the longest read", which let an
intron-retaining read drag the union, Amendment 12): the medoid of the cluster under a structural distance read from the net's
all-vs-all (§5.3): d(m, p) = the bases of indels >= 20 bp (insertions, deletions and `~` introns alike) in m's best alignment to p PLUS
p's terminal bases (>= 20 bp at either end) that the alignment leaves uncovered. A member is eligible when it is aligned to >= min(0.5 x
(n - 1), 50) of the n - 1 other members (50 = half of the all-vs-all's `-N 100`); the template is the eligible member with the lowest
mean d over its aligned partners, ties -> longest -> smallest name; with no eligible member, the longest member that has an aligned
partner, else the longest member (a member with no aligned partner is never chosen while another has one). A cluster of one is its own
template.

**Vote** (rulings R2 / R5, Amendment 13): the members are aligned to the template with the splice preset (`MM2_MEMBERS = -x splice:hq
-uf -c --cs -N 5 -p 0.5`; asm20 cut alignments at exon skips, as R6 found for the union) and the template is polished by column majority
over them: a column with >= 3 covering members takes the majority base; before a column at most one insertion is made, the insertions
>= 20 bp first (structure: the most frequent one with >= 3 carriers, whatever its share), and only without one the < 20 bp majority
insertion (>= 50% of >= 3 covering members); a deletion < 20 bp carried by >= 50% of >= 3 covering members is applied, a deletion >= 20
bp never (the consensus is the exon union of the cluster's isoforms); columns with < 3 covering members keep the template; template
ends covered by fewer than 2 members are trimmed.

**Refinement** (R5): one pass — the members are re-aligned to the consensus and those that do not fit (`de` > delta or < 50% of the
shorter covered) are split off, as one new cluster polished on its own structural template when they are >= `--min-cluster`; a kept set
that still holds its template is re-polished on it, and a kept set whose template was split off is re-templated by the same structural
rule over its own pairs and re-polished (Amendment 13).

### 5.5 Cluster merge by the significance test

Pairs of clusters whose consensus sketches share >= 50% are aligned (poasta, as §5.3); k = distinguishing columns (mismatches; indels
in homopolymer runs excluded, as `read_conflict.rs` does). The smaller cluster B (n_B reads) is a real variant of A (n_A) when the
number of B's reads carrying B's base at all k columns is unlikely under error: p = P(X >= n_B), X ~ Binomial(n_A + n_B, eps^k) with
eps = 0.001 (the per-column error proxy of `read_conflict.rs:77`) and alpha = the assignment gate's alpha (the same constant
`read_conflict.rs` uses). If p >= alpha the clusters merge (union of reads, re-polished on A's template); if k = 0 they merge. Iterate
until no pair merges. This is IsoCon's statistical test in our own code, and the same test O2 uses to de-tie. (Amendment 13: a merged
cluster is re-polished on the structural template of all its members (§5.4), not on A's; a merge whose re-polished consensus is empty is
undone — the clusters stay separate, nothing is dropped — and that pair is not tested again.)

### 5.6 Flag, link, merge, floor (the chain, Amendments 7-9)

1. minimap2 (`-c -x splice:hq -uf -N 20 --eqx`) of every consensus against `--index` -> PAF (cached, kind `cand`).
2. Flag: best hit identity x coverage < 0.999 (identity = matches / block length; coverage = aligned query span / length).
3. Link: d = 1 - matches / length of the best hit; d <= delta -> the cluster is an allele of that locus (`linked_to` = locus, not a
   candidate). Else new copy.
4. Merge: among a family's new-copy consensus sequences, components of pairs whose poasta alignment covers >= 50% of the shorter with
   `de` <= delta (Amendment 8's rule; in-process instead of minimap2 asm20).
5. Floor: a component is a **flagged candidate** iff it holds >= `--min-clusters` clusters. Unflagged components are written with
   `flag = 0` and do not enter the augmentation.

### 5.7 Exon-union representative

For a component: backbone = the longest consensus; the others in decreasing length are each aligned to the CURRENT union (poasta
global); every insertion >= 20 bp relative to the union is spliced in at its aligned position, unaligned prefixes/suffixes >= 20 bp are
appended at the respective end; insertions < 20 bp are ignored (they are errors or microindels, not exons). The union is the
representative written to `P.cand.contigs.fa` as `>cand_<family>_<k>` (k = component order within the family), with its members in
`P.cand.clusters.tsv`.

### 5.8 Caching

`run_cache` kind `cand` (new arm in the `match kind` at `run_cache.rs:320`, required `candidates.tsv`, `contigs.fa`). Key = header
`rustle o3 candidates v1`, `cmd`, `exe_fingerprint()`, `file_fingerprint(bam, bai, fasta, fai, copies.tsv, copies.fa, index)`,
`env_fingerprint`. The minimap2 PAF of §5.6 is a second entry (kind `paf`, replayed as the other binaries do).

## 6. `copy_assign`: chromosome-aware overlap

`AlignedRead` (`copy_split.rs:177`) gains `chrom: String`, filled by `aligned_read_from_record`. `best_overlap_copy`
(`copy_assign_pipeline.rs:1588`) and the mosaic attribution in `assign_family_detailed_once` (`:2223`) require `read.chrom ==
copy.chrom` before comparing positions. Expected effect on existing behaviour: none where a family is on one chromosome (every fixture);
cross-chromosome families stop attributing a read to a copy on another chromosome by numeric coincidence. The 843-test suite must pass
byte-identical except where a test is added for the cross-chromosome case.

## 7. Driver (`tools/rustle_pipeline.sh`)

- `assign` reads `P.fam.copies.tsv` / `.fa` (the default O1 output) and derives `P.regions.txt` from `P.fam.copies.regions` (second
  column; the first is the family id — `parse_region` takes the first token). The legacy `P.cat.*` path stays behind `--legacy-catalog`.
- New stage `candidates`, OPT-IN since ruling R14 (§9b): naming the stage runs it; `all` runs it between `families` and `assign` only
  with `--candidates`, and `assign` and `flag` use its products only with `--candidates` (`--no-candidates`, the default, says so
  explicitly; `--candidates` with `--legacy-catalog` exits 2); `--delta`, `--cand-max-reads` pass through. As first written here it
  was a default stage of `all`. It runs `o3_candidates`, then the augmentation (§4): `P.aug.fa` (+ `samtools faidx`), `P.aug.copies.tsv` rows
  `family_id copy_idx=<next> tid=cand_<f>_<k> chrom=cand_<f>_<k> start=0 end=<len> n_exon=1 strand=+ n_reads=0 exons=0-<len> ...
  source=o3_candidate`, `P.aug.copies.fa` entries `>{fid}|{idx}|cand_<f>_<k>:0-<len>|+|nexon=1`, `P.aug.regions.txt` with
  `cand_<f>_<k>:0-<len>`, and the patch realignment of `P.cand.nets.fa` to `P.aug.fa` with the same minimap2 flags the pipeline's BAM
  was made with (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`), sorted and indexed as `P.aug.bam`.
- `assign --candidates` with flagged candidates: two `copy_assign` runs (candidate families on `P.aug.*`, the rest on the original
  inputs; `--families` restricted by `--only-families` / `--skip-families`), outputs concatenated with one header. Ruling R13
  (§9b): O2 assigns AS-tied molecules only, so a read the realignment places uniquely on a candidate has no row; its placement is
  `P.aug.bam`'s.
- `flag --candidates` passes `--candidates P.cand.candidates.tsv`: `missing_copy_flag` adds a column `o3_candidate` (candidate id or
  `-`) to `P.flag.missing_copy.tsv` when a flagged candidate's nearest locus is the row's locus — the two O3 sources corroborate each
  other.

## 8. Error handling

- Missing or empty `--copies`: exit 2 with the message naming the file. A family with no reads: written with 0 clusters, no error.
- minimap2 absent: exit 2 naming `RUSTLE_MINIMAP2`. A failed minimap2 run: the stage fails; nothing is committed to the cache.
- A net above `--max-reads` is sampled, never truncated silently: `P.cand.candidates.tsv` carries `n_net` and `n_used`.
- The augmentation refuses to write `P.aug.copies.tsv` if a candidate name collides with a FASTA sequence name.

## 9. Testing

Unit (in the module): clustering separates two synthetic copies at 2% divergence with 0.2% random read error and merges them at 0.3%
(delta 0.00958); the consensus recovers the true sequence from 10 reads with errors (incl. a homopolymer indel); the significance test
keeps a 5-read variant at k = 3 and merges a 2-read variant at k = 1; the union of two isoforms sharing exons contains each exon once and
both exclusive exons; flag / link / merge on a hand-written PAF reproduce Amendment 7's `contigs.tsv` logic (same d, same linking).
Fixture: a small BAM/FASTA (`tests/fixtures/o3_candidates/`: a 2-copy family with the second copy deleted from the FASTA) where the
binary flags one candidate and `copy_assign` on the augmented inputs assigns the deleted copy's reads to it. `copy_assign` suite:
byte-identical outputs on every existing fixture after §6, plus one cross-chromosome fixture.
**Acceptance (prereg Amendment 12, written before the run):** the 53-family held-out of Amendment 7 with `o3_candidates` in IsoCon's
place: D right >= 80% of 12,787 and false moves <= 5% -> adopt; the union representatives keep >= 95% of the components' reads (the
`rep_choice.py` measure); wall time <= 2 x IsoCon's (~20 min for 53 families at the 1,000-read cap).
**Re-run acceptance (prereg Amendment 13 + 13b-13e, written before the A13 run):** A13-1 = D right >= 0.80 x C, C = IsoCon's right D
reads over the truth-free attainable D reads (ruling R17), and false moves <= 5%; A13-2 = A12-2; A13-3 = A12-3. All three PASSED on
2026-10-03 (`docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md`).

## 9b. Plan rulings (2026-10-02, recorded here so the spec and the plan agree)

- Alignment engine for read-vs-template, consensus votes, cluster/consensus all-vs-all and the union: **batched minimap2** calls through
  `run_cache` (as `mcl_families` does), not poasta. Poasta's exact affine search runs ~100 ms per 3 kb pair, which puts the 53-family
  acceptance at hours; minimap2's `de`, coverage and `cs` are the quantities Amendments 7-8 validated. §5.3, §5.5 and §5.7 read with
  this substitution; the rules are unchanged.
- §6: the chromosome reaches `best_overlap_copy` as a parallel `read_chroms: Option<&[String]>` slice on `assign_family_detailed_once`
  instead of a field on `AlignedRead` (41 literal constructors); same observable behaviour, byte-identical when `None`.

- Pre-flight rulings R1/R2 (2026-10-02, before Task 1): **R1** the flag floor is `--min-support 6` reads over a component's clusters,
  not `>= 2 clusters` — the stage's clusters merge a copy's isoforms at delta (Amendment 8's rule applied to reads), so the cluster count is
  not IsoCon's transcript count; 6 = 2 x IsoCon's 3-read transcript minimum. The >= 2-cluster count is reported beside in the acceptance.
  **R2** §5.4 consensus: indels >= 20 bp are structure — an insertion >= 20 bp carried by >= 3 members is inserted whatever its share, a
  deletion >= 20 bp is never applied — so a cluster's consensus is the exon union of its reads' isoforms (the representative decision
  carried down one level); indels < 20 bp follow the 50% majority.

- Rulings made during the implementation (2026-10-02):
  - **R3** (§6): `detect_and_assign` hands a supplied CROSS-chromosome family, besides the reads on its first chromosome within the
    family span, the reads on each of its other chromosomes that overlap that chromosome's copy hull, each read matched only to the
    copies of its own chromosome. Byte-identical for single-chromosome families; for existing `~xchrom~` families a behaviour change
    (their reads on other chromosomes were silently never assigned), and the precondition for O2 over candidate contigs (§4, §7).
  - **R6** (§5.7): the member-vs-union alignment is minimap2 `-x splice:hq -uf -c --cs -N 5 -p 0.5` (`MM2_UNION`), not `asm20`:
    measured on Amendment 8's 50 multi-member components (540 real IsoCon contigs), `asm20` cut the alignment at an exon skip and
    inserted the far side again, duplicating 11.7% of the unions' bases; `splice:hq` duplicated 0% and contained all 540 members.
  - **R13** (§7): O2's scope is AS-tied molecules (assign-or-abstain, user 2026-09-09). A deleted copy's reads that realign uniquely
    to its candidate (the fixture: 60/60, MAPQ 60) are placed by the aligner and never enter the certificate; the candidate family is
    assigned as a family that includes its candidate copies, under the unchanged AS-tied gate. Amendment 12 scores placement (arm M);
    the `--no-as-tied-only` route is not taken.
  - **R14** (after Amendment 12): A12-1 and A12-2 FAILED (D right 5,995 vs the bar 10,230; union representatives 90.9% vs 95%;
    `docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`), so the stage does not replace IsoCon yet: the driver's `candidates` stage is
    OPT-IN (`--candidates`; default off; `all` skips it) until a new prereg (A13) passes. The code ships, inert by default.
  - **R15** (final review): in `copy_assign`, a sweep bound to no family skips the §6gz tie-outside registration only when its region
    holds a read window of a cross-chromosome (`~xchrom~`) family (every candidate family is one); a catalog without cross-chromosome
    families registers exactly as before 2026-10-02 (A/B against a b9c412f9 build on the human chr16 O2 simulation, the O2 and
    `--union-certificate` commands of `figures/_o2.py`: every table byte-identical).

- Amendment 13 and its rulings (2026-10-03, each written into the prereg before the A13 run):
  - **Amendment 13** (`e3e9d4bf`): the net attribution by alignment (§5.1-§5.2) and the structurally central template (§5.4), with the
    consensus details the reviews named (splice preset for the votes, the insertion vote by size class, the refinement re-template, the
    empty-merge fallback, §5.5); the chain's rules (delta, the merge rule, `--min-support 6`, 0.98, the 1,000-read cap, R13) unchanged.
  - **R16** (Amendment 13b, `d69f0e02`): the attribution rule is the family's own edge rule in read space (the best hit covers >= 50% of
    the READ, `de` <= 0.20), applied to the unmapped AND the poorly placed un-netted reads, against the families' net reads plus the
    copies, with map-ont (§5.2). Measured at the attribution step only: no truth-free rule reaches the ~5,200 reads IsoCon received by
    label in Amendment 8.
  - **R17** (Amendment 13b): A13-1's comparator is re-registered to C = IsoCon's right D reads over the truth-free attainable D reads (a
    record on a surviving copy of their family in `R.bam`, or attributed by the rule); A12-1's bar (10,230) is reported beside, not
    decided on.
  - **R18** (Amendment 13c, `065b2b46`): under `--families`, "no record on any family copy" (the poorly placed selection) is judged
    against THIS run's families (pass A's scope), and the attribution targets are this run's nets plus every copy; batched runs may
    attribute a read in more than one batch.
  - **R19** (Amendment 13c): the 300-bp floor applies to the poorly placed reads as well.
  - **R20** (Amendment 13d, `f31f663e`): the template is the medoid under the structural distance (big indels + the partner's uncovered
    terminal bases), mean over aligned partners, eligibility by aligned fraction, ties longest then name (§5.4); Amendment 13's "lowest
    total indel bases" picked fragments.
  - **R21** (Amendment 13e, `ff869c40`): eligibility = aligned to >= min(0.5 x (n - 1), 50) other members (the all-vs-all's `-N 100` made
    50% unattainable in 400+-read clusters); a member with no aligned partner is never chosen while another has a mean.

## 10. Open items (deferred, named)

- POA-graph consensus (`build_poa_graph`, `family_graph.rs:18`) instead of template-and-vote, if the acceptance run shows consensus
  errors (visible as clusters that fail to link at delta although their reads are the survivor's).
- The k-mer attribution thresholds (30%, 2 x) are retired (Amendment 13). The alignment rule (R16) reaches 1,482 of the held-out's
  17,286 deleted-copy reads by attribution; poorly placed reads joined a net of the wrong family in 244 of 1,171 joins (batched, R18).
  A library with long unmapped reads of another kind would revisit the read-coverage floor.
- The two-run O2 split is v1; a single run over a merged BAM is the cleaner end state once BAM patching is worth its cost.
