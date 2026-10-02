# O3 candidate copies in the pipeline (`o3_candidates`) — design

Date 2026-10-02. Branch `machine2/soto-evidence` (== `main`). Status: spec for review; no code yet.

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
assemble -> families (O1: P.fam.copies.{tsv,fa,regions}) -> candidates (new) -> assign (O2, augmented) -> flag (O3 screen)
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
   record; (b) every unmapped record with length >= 300 bp, attributed by the k-mer index of §5.2.
4. Cap: if a family's net exceeds `--max-reads`, sample that many (seed 1, after sorting names).
5. Output: `P.cand.nets.fa` holds the nets of families that end with a flagged candidate (§5.6); all nets are kept in the cache entry.

### 5.2 Unmapped-read attribution

Index: canonical 31-mers of every sequence in `--copies-fa`, each k-mer -> set of family ids (k-mers shared by > 8 families dropped as
repeats). A read is attributed to family F when >= 30% of its 31-mers hit F and F's hits are >= 2 x the runner-up's; otherwise it is
not in any net. Reads attributed to a family join its net before the cap.

### 5.3 Read clustering (greedy template clustering, per family)

Reads sorted by length descending. For each read: candidate clusters = those whose template shares >= 50% of the read's 31-mer sketch
(minimizers, w = 5); the read is aligned to each candidate's template with the poasta global aligner used by `family_graph.rs:23-62`
(`PoastaAligner`, `AlignmentType::Global`, `GapAffine` costs as there); it joins the first cluster (in creation order) whose alignment
covers >= 50% of the shorter sequence with gap-compressed divergence `de` <= delta; otherwise it founds a new cluster with itself as
template. `de` = (mismatches + gap openings) / aligned columns, computed from the alignment's operations — the same quantity minimap2's
`de` tag reports, so the thresholds carry over from Amendment 8. The read's alignment to its template is kept for §5.4.

### 5.4 Consensus per cluster (template and vote)

Clusters with < `--min-cluster` reads are dropped (IsoCon's floor of 3). Consensus = the template polished by column majority over the
cluster's alignments: a column with >= 3 covering reads takes the majority base; an insertion relative to the template present in >= 50%
of the reads covering that position is inserted (the majority insertion sequence); a deletion in >= 50% is applied; columns with < 3
reads keep the template. Template ends covered by fewer than 2 reads are trimmed. One polishing pass (HiFi).

### 5.5 Cluster merge by the significance test

Pairs of clusters whose consensus sketches share >= 50% are aligned (poasta, as §5.3); k = distinguishing columns (mismatches; indels
in homopolymer runs excluded, as `read_conflict.rs` does). The smaller cluster B (n_B reads) is a real variant of A (n_A) when the
number of B's reads carrying B's base at all k columns is unlikely under error: p = P(X >= n_B), X ~ Binomial(n_A + n_B, eps^k) with
eps = 0.001 (the per-column error proxy of `read_conflict.rs:77`) and alpha = the assignment gate's alpha (the same constant
`read_conflict.rs` uses). If p >= alpha the clusters merge (union of reads, re-polished on A's template); if k = 0 they merge. Iterate
until no pair merges. This is IsoCon's statistical test in our own code, and the same test O2 uses to de-tie.

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
- New stage `candidates` (in `all` between `families` and `assign`; `--no-candidates` skips it; `--delta`, `--cand-max-reads` pass
  through). It runs `o3_candidates`, then the augmentation (§4): `P.aug.fa` (+ `samtools faidx`), `P.aug.copies.tsv` rows
  `family_id copy_idx=<next> tid=cand_<f>_<k> chrom=cand_<f>_<k> start=0 end=<len> n_exon=1 strand=+ n_reads=0 exons=0-<len> ...
  source=o3_candidate`, `P.aug.copies.fa` entries `>{fid}|{idx}|cand_<f>_<k>:0-<len>|+|nexon=1`, `P.aug.regions.txt` with
  `cand_<f>_<k>:0-<len>`, and the patch realignment of `P.cand.nets.fa` to `P.aug.fa` with the same minimap2 flags the pipeline's BAM
  was made with (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`), sorted and indexed as `P.aug.bam`.
- `assign` with candidates: two `copy_assign` runs (candidate families on `P.aug.*`, the rest on the original inputs; `--families`
  restricted by a family list file each binary already accepts or gains), outputs concatenated with one header.
- `flag` gains `--candidates P.cand.candidates.tsv`: `missing_copy_flag` adds a column `o3_candidate` (candidate id or `-`) to
  `P.flag.missing_copy.tsv` when a flagged candidate's nearest locus is the row's locus — the two O3 sources corroborate each other.

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

## 9b. Plan rulings (2026-10-02, recorded here so the spec and the plan agree)

- Alignment engine for read-vs-template, consensus votes, cluster/consensus all-vs-all and the union: **batched minimap2** calls through
  `run_cache` (as `mcl_families` does), not poasta. Poasta's exact affine search runs ~100 ms per 3 kb pair, which puts the 53-family
  acceptance at hours; minimap2's `de`, coverage and `cs` are the quantities Amendments 7-8 validated. §5.3, §5.5 and §5.7 read with
  this substitution; the rules are unchanged.
- §6: the chromosome reaches `best_overlap_copy` as a parallel `read_chroms: Option<&[String]>` slice on `assign_family_detailed_once`
  instead of a field on `AlignedRead` (41 literal constructors); same observable behaviour, byte-identical when `None`.

## 10. Open items (deferred, named)

- POA-graph consensus (`build_poa_graph`, `family_graph.rs:18`) instead of template-and-vote, if the acceptance run shows consensus
  errors (visible as clusters that fail to link at delta although their reads are the survivor's).
- Unmapped-read attribution thresholds (30%, 2 x) are set from the deletion tests' unmapped reads (median 69 bp, no evidence); a
  prereg on a library with long unmapped reads would revisit them.
- The two-run O2 split is v1; a single run over a merged BAM is the cleaner end state once BAM patching is worth its cost.
