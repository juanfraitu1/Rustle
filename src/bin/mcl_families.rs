//! `mcl_families` — define multi-copy gene families by MCL clustering of the ANNOTATION, by sequence,
//! then corroborate them with RNA.
//!
//! ⭐ **DNA PROPOSES, RNA DISPOSES.** Sequence clustering generates candidate families; read support
//! decides which are real. Clustering alone is OrthoFinder; the corroborated fraction is the
//! contribution. Measured genome-wide on gorilla (§6de): a repeat clique is a PERFECT clique
//! (median density **1.000**) and a real family is not (**0.700**) — at identical median size 4.0.
//!
//! ⚠ **NO GENE SYMBOL ENTERS THE DEFINITION.** Genome-wide 95.2% of members across 1,021
//! product-defined families are `LOC*`-named; a `gene=NPIP*` grep recovers 1 of 44 NPIP members.
//!
//! INPUT is the all-vs-all PAF of annotated gene sequences, e.g.
//! `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -K 100M genes.fa genes.fa`
//! whose FASTA headers are `CONTIG:START-END` in **GFF 1-based coordinates, verbatim**.

use anyhow::{Context, Result};
use clap::Parser;
use rustle::vg_family::annotation_families::{
    build_clusters, fold_parts_into_loci, graph_from_paf_loci, loci_from_exon_blocks, mcl, sd_blocks, Cluster, CoreStatus,
    GeneKey, GraphParams, SdPairs,
};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

/// Which transcript of a `--from-gtf` locus is its representative (`--representative`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Representative {
    /// the transcript with the most `reads`; ties to the longer span, then the last `transcript_id` (the shipped rule)
    MostReads,
    /// the transcript with the most junctions; ties to the most `reads`, then the longer span, then the last `transcript_id`
    MostJunctions,
}

#[derive(Parser, Debug)]
#[command(about = "Multi-copy gene families by MCL over annotation sequence, corroborated by RNA")]
struct Args {
    /// All-vs-all PAF of the locus sequences (`minimap2 -x asm20 -c -X`). Not needed with `--from-gtf`.
    #[arg(long, default_value = "")]
    paf: String,

    /// One-command de novo family stage: derive the loci from an assembled GTF (`copy_assign --assemble-only`
    /// output; locus = `gene_id` group, span = min/max over its transcripts, representative = the transcript with
    /// most `reads` (tie: longer span; see `--representative`), exons = the representative's), write `<out>.loci.gff3` and
    /// `<out>.loci.fa` (genomic spans, `--fasta` required), run the all-vs-all (`minimap2 -x asm20 -c -X -N 50
    /// -p 0.1 --secondary=yes`) to `<out>.loci.paf`, then proceed as with `--paf <out>.loci.paf --gff
    /// <out>.loci.gff3`. Replaces the former scratch step `loci_from_gtf.py` + a hand-run minimap2.
    #[arg(long)]
    from_gtf: Option<String>,

    /// With `--from-gtf` only: which transcript of a locus is its representative — the locus's exons in
    /// `<out>.loci.gff3` and the copy in `<out>.copies.tsv`. `most-reads` (default; every product before
    /// 2026-10-04): the transcript with the most `reads`, ties to the longer span, then the last `transcript_id`.
    /// `most-junctions` (opt-in; pre-registered test, `docs/PREREG_locus_representative_rule_2026-10-04.md`): the
    /// transcript with the most junctions — gaps of >= 50 bp between consecutive exons, the junction floor of the
    /// strict "found" rule (`bench/copy_support.py`) — ties to the most `reads`, then the longer span, then the
    /// last `transcript_id`. In a 5'-truncated library the most-read transcript can be a 3' fragment (NPIPA9: 1
    /// junction of the locus's 38). Locus spans, so `<out>.loci.paf`, do not depend on the rule. Driver:
    /// `RUSTLE_REPRESENTATIVE`.
    #[arg(long, value_enum, default_value_t = Representative::MostReads)]
    representative: Representative,

    /// minimap2 threads for `--from-gtf`.
    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// GFF supplying each gene's EXON-UNION length. ⚠ Without it the coverage denominator falls back to
    /// the genomic span, which removed 62.8% of eligible genes in the pilot (§6dc) — the run warns loudly.
    #[arg(long)]
    gff: Option<String>,
    /// PREREG (identity-weighted density, 09-10): also write `<out>.pairs.tsv` (`cluster_id a b identity`,
    /// one row per WITHIN-cluster homology edge — `identity_gap.py`'s own input format), reusing
    /// `g.idents`/`g.edges` directly so a downstream metric is scored on the SAME aggregated identity
    /// (Σnmatch/Σblocklen per pair, the deferred-exonic path) that decided admission, not a re-derived one.
    /// Default off; output-only, existing files unchanged.
    #[arg(long, default_value_t = false)]
    dump_pairs: bool,

    /// MCL inflation. ⭐ With a size-safe prune (§6ec) the anchored families are STABLE from I=2.0 to 4.0
    /// (NPIP 44/44, both tandem halves intact); §6dd's "cliff at 3.6" was the old prune emptying NPIP's
    /// columns. 2.8 is kept as the historical default, not a tuned value.
    #[arg(long, default_value_t = 2.8)]
    inflation: f64,

    #[arg(long, default_value_t = 0.70)]
    min_identity: f64,

    /// Coverage floor on the LONGER sequence. ⚠ Not the shorter: a ~300 bp Alu covers most of a fragment
    /// and almost none of a real gene, and every adjudicated NPIP false merge was that shape (§6cr).
    #[arg(long, default_value_t = 0.30)]
    min_cov_longer: f64,

    /// ⭐ §6x4 CONTAINMENT ESCAPE — **default `0.70` since 2026-09-29** (the user's decision; shipped opt-in in
    /// f2144faf, register 1006/1014: held-out F up on 2 of 3 chromosomes, precision up-or-equal everywhere,
    /// insensitive to C over 0.40-0.80). `--min-cov-shorter 0` = OFF, byte-identical to every catalog built
    /// before the flip. When `> 0.0`, a pair ALSO passes the coverage gate if the alignment covers at least
    /// this fraction of the SHORTER gene's exonic length, even when `cov_longer` fails, and its edge weight is
    /// that coverage. ⚠ Known regressions (register 1007/1009): NPIP in GUIDED mode, Soto F .833 -> .800
    /// (sensitivity .750 -> .700); semi-guided SD-region nodes, precision .973 -> .833 — an SD region has no
    /// gene boundary, so never use it with `--from-genome-sd` nodes.
    ///
    /// Motivation (§6x3/r1002): of the 1,026 chr16 de novo loci that fail admission, 679 are evicted by a
    /// LONGER partner while aligning along a median 1.00 of their own length at passing identity and
    /// `alen` — nothing is wrong with them except the denominator. Splitting the giant cannot fix it
    /// (r1001: median 4.83x shrink required, a binary split gives 2x) and boundary pull-in is negative
    /// (§6w6).
    ///
    /// ⚠⚠ **Never enable without the exon conjunct** (`--min-exonic-bp 1 --min-shared-exon-frac`).
    /// §6x4 measured the difference on chr16: guarded at 0.90 admits 216 of the 679 at a largest-component
    /// cost of 2.1x baseline; UNGUARDED the same norm runs 5.7-16.3x — which is register 913's refuted
    /// `min(la,lb)` hub failure. The guard is the result, not the norm.
    #[arg(long, default_value_t = 0.70)]
    min_cov_shorter: f64,

    #[arg(long, default_value_t = 300)]
    min_bp: u64,

    /// Smallest cluster reported. **2** (user decision 2026-09-05, §6ex): two loci sharing a duplicated core
    /// is the minimum object the definition names. Pre-registered on Soto's slice, 3 → 2 recovered 18 of the
    /// 362 members (19 of the 43 annotated misses sat in size-2 clusters) with the [0.90,1) band precision
    /// unchanged (0.949 → 0.949) and no new pair in the Soto-silent 0.80–0.90 band. Density carries no
    /// signal on a pair; for pairs the SD-core certificate carries the whole weight. `--min-size 3`
    /// reproduces every catalog built before this date.
    #[arg(long, default_value_t = 2)]
    min_size: usize,
    /// MCL prune threshold: after inflation, matrix entries below this are dropped. ⚠ An ABSOLUTE
    /// threshold interacts with cluster SIZE: in a near-uniform clique of n nodes every entry is ~1/n
    /// and inflation maps it to (1/n)^I, so the whole column empties once n > prune^(-1/I) — ≈61
    /// nodes at I=2.8 with the old 1e-5 (§6ec: the anchored 84+22-copy tandem array dissolved genome-wide).
    /// ⭐ Default 1e-9 (§6ed, user decision 2026-09-04): safe to n ≈ 1,635 at I=2.8. Every catalog before
    /// §6ed was built at 1e-5; pass `--prune 1e-5` to reproduce them byte-for-byte.
    #[arg(long, default_value_t = 1e-9)]
    prune: f64,

    /// ⭐ Charge `cov_longer`'s NUMERATOR in exonic bases too. Without it the numerator is aligned
    /// GENOMIC span and the denominator is exon-union length — different units, so intronic repeat
    /// homology satisfies an exonic floor. Default OFF ⟹ byte-identical to every prior run.
    #[arg(long, default_value_t = false)]
    exonic_overlap: bool,

    /// ⭐ Reject pairs whose annotation intervals overlap on the same contig (the `q == t` guard only
    /// catches identical headers, so a nested gene aligns to its host as a "paralog"). Default OFF.
    #[arg(long, default_value_t = false)]
    reject_overlapping: bool,

    /// ⭐ Minimum ABSOLUTE exonic bases an edge must rest on (additive; `cov_longer` is unchanged).
    /// 0 = off. ⚠ Prefer this over --exonic-overlap: replacing the coverage measure shatters real
    /// families, because a segmental duplication copies introns too.
    #[arg(long, default_value_t = 0)]
    min_exonic_bp: u64,
    /// ⭐ A LOCUS is the node (§6ee): annotation records whose EXON-UNIONS overlap on a contig (a lncRNA
    /// model over a gene's exons, two models of one transcription unit, an antisense model over the same
    /// exons) are ONE node — a gene inside another's INTRON stays a separate locus — represented by the
    /// record with the greatest exon-union length, and an edge admitted for ANY record of the locus is the
    /// locus's edge; records between two annotations of one locus are skipped. Represented
    /// by the record with the greatest exon-union length; PAF records of the folded-away annotations are
    /// skipped. Writes `<out>.loci.tsv` (annotation -> representative). Measured on NPIP: 13 overlapping
    /// copy pairs, 607/1,221 O2 ties were the same locus twice. Default OFF ⟹ byte-identical
    #[arg(long, default_value_t = true)]
    merge_overlapping_loci: bool,
    /// Escape hatch: annotation RECORDS as nodes (every catalog before §6er was built this way).
    #[arg(long, default_value_t = false)]
    no_merge_overlapping_loci: bool,
    /// Locus evidence = ATTRIBUTION (every model's admitted edge is the locus's edge). Default OFF =
    /// representative-only. ⚠ Attribution reconstructs the duplication BLOCK (LCR16a+LCR16u, §6eg).
    #[arg(long, default_value_t = false)]
    locus_attribute_edges: bool,
    /// SEDEF pairs (BED: chr1 s1 e1 chr2 s2 e2 …, 0-based half-open) for `--core-refine`.
    #[arg(long)]
    sedef: Option<String>,
    /// Measurement scaffolding (§6iz proposal #4, `docs/o1_ledger.md`): dump the PRE-MCL homology graph
    /// (`node<TAB>weight` self-row per node, then `u<TAB>v<TAB>weight` per edge, nodes as
    /// `chrom:start-end`) to this path before partitioning, so a controlled ablation can feed the SAME
    /// graph into a different partitioner (e.g. `gamma_refine`) instead of re-deriving it. Not consumed by
    /// any shipped pipeline; default off, no effect on any existing output.
    #[arg(long)]
    dump_graph: Option<String>,
    /// ⭐ §6fo: derive the SD-like pairs for the core refinement (and `blocks.tsv`) from the input `--paf`
    /// itself — every alignment between two annotated loci is one pair of genomic intervals — instead of a
    /// SEDEF bed. No external SD caller; the intergenic SD context is not available this way. Ignored when
    /// `--sedef` is given.
    #[arg(long, default_value_t = false)]
    core_from_paf: bool,
    /// §6ft polish 2: the core majority counts the locus itself (shared with ≥ half of the family, depth + 1 ≥ n/2)
    /// instead of ≥ half of the other members. Default OFF (pre-registered test pending).
    #[arg(long, default_value_t = false)]
    core_majority_inclusive: bool,
    /// §6ft polish 1: a member whose reads leave no block in its chain keeps its annotated model as its unit
    /// (`source gff_fallback`, `n_reads 0`) — membership does not depend on expression. Default ON.
    #[arg(long, default_value_t = true)]
    units_keep_unexpressed: bool,
    /// Escape hatch: skip members with no read in their chain (the row set before §6ft).
    #[arg(long, default_value_t = false)]
    no_units_keep_unexpressed: bool,
    /// ⭐ DUPLICON-FIRST refinement (§6eh; pre-registered adj/core/PREREG.md): within each cluster, a
    /// member's CORE is the part of it linked by SEDEF pairs to ≥ half the other members. Clusters whose
    /// members lack SD depth (median max-depth < half) are UNTOUCHED (old ZNF/OR families have no SEDEF
    /// pairs). In SD-evidenced clusters: core ≥ span/2 ⟹ kept (full/partial copy); else core ≥ half the
    /// cluster's median core ⟹ kept and TRIMMED to the core hull (a chimeric model: ABCC1+NPIP, SORL1+NPIP);
    /// else DROPPED (EIF3C carries 808 bp of NPIP's 23 kb LCR16a core). Writes `<out>.cores.tsv` and
    /// `<out>.refined.clusters.tsv`; `<out>.clusters.tsv` is byte-identical. Default OFF
    #[arg(long, default_value_t = false)]
    core_refine: bool,
    /// Escape hatch: do NOT run the core refinement even though `--sedef` is given (§6er: it runs whenever
    /// a SEDEF bed is supplied).
    #[arg(long, default_value_t = false)]
    no_core_refine: bool,
    /// ⭐ O1-10b (§6el): emit per-locus UNITS in the `copy_assign --families/--copies-fa` contract:
    /// `<out>.units.tsv` / `<out>.units.fa` / `<out>.units.regions`. A unit = the member's locus (its core hull
    /// under `--core-refine`, else its span) with the READ-SUPPORTED exon chain: primaries with an aligned block
    /// inside the locus, cut at introns >50 kb whose junction has <3 reads (the shipped mis-chain rule), exon
    /// chain = bases covered by >= `--min-reads` reads, strand = majority transcript strand (`ts` tag, else
    /// flag) — a base is exonic iff covered by >= `--min-reads` aligned blocks AND by more blocks than reads
    /// that splice over it (§6en: pre-mRNA reads must not glue exons across an intron the majority splices);
    /// a locus with < `--min-reads` reads keeps its GFF exons (reported). Needs `--bam` and `--fasta`.
    /// Measured (§6el): with these units O2's junction-anchored agreement went 7/11 -> 5/5 and the control
    /// 52/52 -> 57/57 — the GFF model was the cause of O2's confident wrong calls. Default OFF
    ///
    /// ⭐ WITH `--from-gtf` (the de novo families stage; user decision 2026-09-25: copy assignment consumes the
    /// SAME families) the unit is the LOCUS REPRESENTATIVE and no `--bam` is needed: `<out>.copies.tsv` /
    /// `<out>.copies.fa` / `<out>.copies.regions` / `<out>.copies.merged.tsv`, one copy per member of every
    /// cluster of `clusters.tsv` (`family_id` = `cluster_id`), in the `gw_family_catalog` copies contract
    /// (columns 1-11 identical: `copy_assign --families/--copies-fa` and `bench/sim.py copies` read it as they read
    /// the legacy catalog). Copy = the representative transcript (`tid` = its `transcript_id`, exons = its exons,
    /// sequence = its spliced exon sum, `n_reads` = its `reads`), source `locus_rep`; see `write_locus_rep_copies`.
    /// The read-chain units above are the annotation mode's and are not written. Without `--emit-units` a
    /// `--from-gtf` run is byte-identical to before.
    #[arg(long, default_value_t = false)]
    emit_units: bool,
    /// Escape hatch: do NOT emit units even though `--bam` and `--fasta` are given (§6er: units are emitted
    /// whenever both are supplied).
    #[arg(long, default_value_t = false)]
    no_emit_units: bool,
    /// With `--from-gtf` only; default OFF (every product byte-identical without it). The CONTAINER of each family
    /// member's extra pieces (`rustle::vg_family::family_container`, the port of the frozen `bench/family_container.py`,
    /// `docs/PREREG_fusion_container_sim_2026-09-28.md` §1 + Amendment 1): after the families are written, every
    /// clustered locus (with the records `loci.tsv` folds into it) gets its exon blocks = the union of the exons of
    /// ALL transcripts of its gene_ids in the `--from-gtf` GTF; a block is `core` iff one aligned CIGAR column of
    /// `<out>.loci.paf` joins it to an exon base of another member of the same family, else `accessory`, and an
    /// accessory block records every OTHER family it aligns to the same way. Writes `<out>.container.tsv` (one row per
    /// locus x block), `<out>.container_relations.tsv` (one row per directed family relation, with `reciprocal`) and
    /// `<out>.container_summary.tsv` (counts), byte for byte the script's `blocks` / `relations` / `summary` tables.
    /// It never changes a family. Driver: `RUSTLE_FAMILY_CONTAINER=1`.
    #[arg(long, default_value_t = false)]
    emit_container: bool,
    /// Optional RepeatMasker `.out` (curated library) — adds `rep_frac` (interspersed-repeat fraction of the
    /// unit's exon chain) to the unit table.
    #[arg(long)]
    rmsk: Option<String>,
    /// Genome FASTA (for `--emit-units` sequences).
    #[arg(long)]
    fasta: Option<String>,

    /// Optional RNA BAM. Without it `corroborated` is reported as `NA` — ⚠ which is NOT 0.000, the
    /// repeat-clique signature. The two must never be conflated.
    #[arg(long)]
    bam: Option<String>,

    /// A member is corroborated when it carries at least this many primary reads with an ALIGNED BLOCK
    /// inside it (`-F 2308`; a read spliced OVER a locus is not evidence for it).
    #[arg(long, default_value_t = 3)]
    min_reads: usize,

    /// ⭐ Fold overlapping annotation records into loci AFTER clustering, inside each cluster (§6ey), instead
    /// of before (the default since §6ef). Fold-first with representative-only evidence loses every record that
    /// overlaps a DIFFERENT family's record on exon bases (Soto: 19 of 43 annotated misses — ANAPC1P1 folded
    /// under CD8B's locus, PMS2P4 under SPDYE21, FAM72A under SRGAP2, LRRC37A under ARL17B); attribution rebuilds
    /// the duplication block. With this flag records are the graph's nodes, MCL runs on them, and two records
    /// become one locus only if they overlap on exon bases AND share a cluster. **Default ON (user decision 2026-09-05, §6ey)**;
    /// `--no-fold-within-clusters` restores fold-first. Implies `--no-merge-overlapping-loci` for the graph;
    /// `--locus-attribute-edges` is ignored.
    #[arg(long, default_value_t = true)]
    fold_within_clusters: bool,
    /// Escape hatch: the behaviour before 2026-09-05 (`--fold-within-clusters` off).
    #[arg(long, default_value_t = false)]
    no_fold_within_clusters: bool,

    /// ⭐ A gene/pseudogene record with NO exon children counts as one exon spanning the record (§6ey). RefSeq
    /// (CHM13) leaves 160 of 747 records in Soto's neighbourhoods without exon features (NF1P, CNTNAP3P, PMS2P
    /// pseudogenes); without a block they have no exonic denominator and no exonic overlap, so no edge can
    /// reach them (20 Soto members exist only as such records). The gorilla annotation has none. **Default ON
    /// (user decision 2026-09-05, §6ey)**; `--no-exonless-span` restores the old behaviour.
    #[arg(long, default_value_t = true)]
    exonless_span: bool,
    /// Escape hatch: the behaviour before 2026-09-05 (`--exonless-span` off).
    #[arg(long, default_value_t = false)]
    no_exonless_span: bool,

    /// ⭐ The pair's alignment must cover ≥ `--min-exonic-bp` exon bases on BOTH records (§6ey). Records are
    /// aligned as genomic spans, so a pseudogene inside another family's gene aligns to that family's paralogs
    /// on the host's exons alone (Soto: PMS2P7 inside SPDYE8 → 26 false PMS2P×SPDYE pairs under
    /// `--fold-within-clusters`). Homologous copies share exon bases on both sides. **Default ON (user decision
    /// 2026-09-05, §6ey)**; `--no-exonic-both-sides` restores the one-sided rule.
    #[arg(long, default_value_t = true)]
    exonic_both_sides: bool,
    /// Escape hatch: the behaviour before 2026-09-05 (`--exonic-both-sides` off).
    #[arg(long, default_value_t = false)]
    no_exonic_both_sides: bool,

    /// ⭐⭐ §6ks: the pair's best record must cover this FRACTION of the smaller gene's exonic length with
    /// shared exon-to-exon evidence (not merely >= 1 bp, `--exonic-both-sides`'s structural floor). Measured
    /// against Soto et al. 2025's family calls: pairs both definitions agree on share a median 52-89% on their
    /// best record; pairs this project alone joined share a median 11-25%, over 1,000 of them <5% — a
    /// co-duplicated neighbour riding one shared base of flanking sequence, not the two genes' own homology.
    /// Implies `--exonic-both-sides` (default on, above). **Default ON at 0.30 (user decision 2026-09-14,
    /// §6ks)**: T=0.30 was fixed from development (chr1/chr15/17) and held out on FRESH Soto families never
    /// used to pick it (chr5/chr7/chr21) — bipartite F (universe) 0.831 -> 0.881, pairwise precision (universe)
    /// 0.815 -> 1.000, zero Soto-verified true pairs lost on any of the 11 held-out families; TBC1D3's 9-copy
    /// family stays whole (reported, not a discovery test there). `--min-shared-exon-frac 0.0` restores the
    /// behaviour before 2026-09-14 byte-for-byte — every catalog in this repo built before that date used 0.0.
    #[arg(long, default_value_t = 0.30)]
    min_shared_exon_frac: f64,

    /// ⭐ Units of one family that share EXON bases are one locus (§6fb): the longest exon union represents them,
    /// the others go to `<out>.units.merged.tsv`. A base cannot belong to two copies (MCL108: a 1.16-Mb read-followed
    /// unit with two units nested in its exons, 13,000 reads counted three times, all K = 0 ties). Default ON.
    #[arg(long, default_value_t = true)]
    merge_overlapping_units: bool,
    /// Escape hatch: keep every unit (the catalogs before 2026-09-05 §6fb, byte-identical).
    #[arg(long, default_value_t = false)]
    no_merge_overlapping_units: bool,

    /// ⭐ Units FOLLOW THE READS (§6eu): extend a unit's read-supported exon chain beyond its core hull to every
    /// block supported by the reads that overlap the hull, WITHIN the member's annotated locus span — the same
    /// `>= --min-reads` rule and the same giant-intron cut, only the hull is no longer a clip. OFF, the
    /// chain is clipped to the window: that emitted the ZNF569-like unit as an 809-bp 5' fragment without the
    /// gene's annotated 3.3-kb 3' exon, and O2 then assigned 190 MAPQ-60 reads of that locus to ZNF875 at
    /// p <= 1e-133 (adj/worst2). **Default ON (user decision 2026-09-05)**; `--no-units-follow-reads` clips to the hull.
    #[arg(long, default_value_t = true)]
    units_follow_reads: bool,
    /// Escape hatch: the behaviour before 2026-09-05 (`--units-follow-reads` off).
    #[arg(long, default_value_t = false)]
    no_units_follow_reads: bool,
    /// ⭐ L1 (`docs/O1_O2_LOOSE_ENDS.md`): emit cluster members DROPPED by the core rule (core = 0) as units too,
    /// with `member_status = dropped`. Family membership is the flag; O2's candidate set is every locus of the
    /// cluster — MCL clustered the locus by homology, so it competes for the family's reads whatever its core
    /// status (NPIP's ABCC1-region records: 332 reads with no candidate, surfacing as a false "missing copy"
    /// signal, §6fg). Default ON (2026-09-05).
    #[arg(long, default_value_t = true)]
    units_include_dropped: bool,
    /// Escape hatch: the row set before L1 (dropped members are not units). Columns `member_status`,
    /// `locus_start`, `locus_end` are still appended; the previous columns are byte-identical.
    #[arg(long, default_value_t = false)]
    no_units_include_dropped: bool,
    /// ⭐ Emit the CONJOINED READ-THROUGH as its own unit (`PREREG_readthrough_object_2026-09-07`, md5 1a51fa3b).
    /// Of the 71 NPIP reads rejected for running past their locus, **48 are real read-throughs** — canonical
    /// splice sites, the intron shared by ≥ 3 molecules, the far end inside another catalog unit (40 of them one
    /// event: NPIP unit 2 → MCL27:0 through a single 15-kb intron). The catalog has no object for them, so they
    /// arrive as rejections. A read-through unit is emitted for an ordered pair of emitted units (A, B) on one
    /// contig when ≥ `--min-reads` distinct primaries have a block in A's chain and share ONE intron whose far
    /// end lands in B's chain, and that intron is canonical (`GT..AG`, `GC..AG`, `AT..AC` on the unit's strand).
    /// `member_status = readthrough`, `source = readthrough`, chain = A's ∪ B's. It is a UNIT, not a member:
    /// `bench/o1_eval.py` counts it with the candidates. Default OFF ⟹ byte-identical.
    /// ⛔ The rejected alternative was `minimap2 --junc-bed` from the annotation: it puts the annotation inside
    /// the read layer, tilts near-ties toward better-annotated copies, and suppresses what O3 looks for.
    #[arg(long, default_value_t = false)]
    emit_readthrough_units: bool,
    /// Escape hatch: emit read-through units WITHOUT the §6fw guard (the unguarded set: 46 units on the three
    /// gorilla contigs instead of 40). The guard is the user's objection turned into a rule — a read-through is
    /// admitted only when the two units share a strand and their donor/acceptor flanks are NOT linked by a
    /// duplication pair, because a cross-copy mis-chain needs the two flanks to be copies of one another.
    #[arg(long, default_value_t = false)]
    no_readthrough_guard: bool,
    /// ⭐ CODING CORE (§6ga addendum): a member must additionally preserve the reading frame the family
    /// preserves — its longest ORF must reach **half the family's best**. Annotation-free and threshold-free
    /// in absolute terms: the comparison is to the family, which is the same majority the core rule takes,
    /// one level up. Failing members are emitted with `member_status = noncoding` (a candidate, not a
    /// member), so they stay visible and are never conflated with the core rule's `dropped`.
    /// ⛔ The rejected alternative was reading the annotation's CDS records: a RefSeq pseudogene has none by
    /// construction, so that test re-derives the annotation's own class label. Default OFF ⟹ byte-identical.
    #[arg(long, default_value_t = false)]
    coding_core: bool,
    /// ⭐ §6ge: forbid two units of DIFFERENT families from claiming the same exon bases. The §6fb merge folds
    /// overlapping units only WITHIN a family, so a cross-family overlap survives — and it manufactures
    /// "read-throughs" out of ordinary introns: `MCL1:3`'s chain overran 2 kb into `MCL27:0`, whose FIRST EXON
    /// is that overlap, so every read splicing across `LOC129527585`'s own first intron looked like a molecule
    /// joining two loci. **32 of 42 guarded read-throughs were annotated introns of one gene.**
    /// Contested bases go to the unit whose own annotated member span contains them; if that does not decide,
    /// to the unit with more reads. Default OFF ⟹ byte-identical.
    /// ⛔ Default stays OFF. The 2026-09-08 flip (PREREG md5 92c2c1e0) was **reverted**: on human chr16+18 —
    /// the one substrate the rule was not designed against — it LOSES `NPIPA1` and `NPIPA6`, taking NPIP
    /// sensitivity 26/26 → 24/26 (P4 failed). The rule is right about the gorilla read-throughs and wrong as
    /// a universal default; turning it on is a per-run decision until the loss is understood.
    #[arg(long, default_value_t = false)]
    no_cross_family_exon_overlap: bool,
    /// Escape hatch: permit two units of different families to claim the same exon bases, reproducing every
    /// catalog built before 2026-09-08. ⚠ With overlap allowed, an ordinary intron of one gene is reported as
    /// a read-through out of its neighbour — 32 of 42 guarded read-throughs were exactly that (register 754).
    #[arg(long, default_value_t = false)]
    allow_cross_family_exon_overlap: bool,

    #[arg(long)]
    out: String,
}

/// Exon-union length per gene, keyed by the gene's own GFF 1-based `(contig, start, end)`.
///
/// ⚠ Reads `exon` features and joins them to their gene by the `gene=` attribute (a Name), which is why
/// the Name -> ID map is built first. Reports the join rate: **a no-op result is the signature of a failed
/// join** (§6dd), so a silent fallback to the span must never be possible.
fn exonic_blocks(gff: &str, exonless_span: bool) -> Result<BTreeMap<GeneKey, Vec<(u64, u64)>>> {
    let text = std::fs::read_to_string(gff).with_context(|| format!("reading {gff}"))?;
    let attr = |a: &str, k: &str| -> Option<String> {
        a.split(';').find_map(|f| f.strip_prefix(k).map(|v| v.to_string()))
    };
    let mut span_of_name: BTreeMap<String, GeneKey> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || !matches!(f[2], "gene" | "pseudogene") {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        if let Some(n) = attr(f[8], "Name=") {
            span_of_name.insert(n, (f[0].to_string(), s, e)); // ⚠ GFF 1-based, verbatim
        }
    }
    let mut blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let Some(name) = attr(f[8], "gene=") else { continue };
        let Some(g) = span_of_name.get(&name) else { continue };
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        if f[0] == g.0 {
            blocks.entry(g.clone()).or_default().push((s - 1, e));
        }
    }
    if exonless_span {
        // a record with no exon children is one exon: its own span (0-based half-open)
        for g in span_of_name.values() {
            blocks.entry(g.clone()).or_insert_with(|| vec![(g.1 - 1, g.2)]);
        }
    }
    let mut out: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
    for (g, mut v) in blocks {
        v.sort_unstable();
        let (mut merged, mut cur) = (Vec::new(), v[0]);
        for &(s, e) in &v[1..] {
            if s <= cur.1 {
                cur.1 = cur.1.max(e);
            } else {
                merged.push(cur);
                cur = (s, e);
            }
        }
        merged.push(cur);
        out.insert(g, merged);
    }
    Ok(out)
}

/// Exon-union LENGTH per gene, derived from the same merge as [`exonic_blocks`] so the numerator and the
/// denominator can never disagree about what an exon is.
/// GFF gene strand per gene key (for `--emit-units` fallbacks and a tie-break when reads carry no strand).
fn gene_strands(gff: &str) -> Result<BTreeMap<GeneKey, char>> {
    let text = std::fs::read_to_string(gff).with_context(|| format!("reading {gff}"))?;
    let mut m = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || !matches!(f[2], "gene" | "pseudogene") {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
        m.insert((f[0].to_string(), s, e), f[6].chars().next().unwrap_or('+'));
    }
    Ok(m)
}

/// Aligned blocks and introns of a read, 0-based half-open on the reference (`D` extends a block; `N` closes it).
fn blocks_and_introns(br: &rustle::vg_family::denovo_assemble::BamRead) -> (Vec<(u64, u64)>, Vec<(u64, u64)>) {
    (br.read.exon_blocks(), rustle::vg_family::copy_split::intron_chain_of(&br.read))
}

/// The mis-chain rule shared by the chain and the extent (§6el; one shipped constant pair: 50 kb / `min_reads`):
/// a record is split at every intron longer than `GIANT_BP` whose junction fewer than `min_reads` records
/// support; the segment returned is the one with an aligned base in `[lo, hi)`, `None` if no segment has one.
const GIANT_BP: u64 = 50_000;
fn kept_segment(
    blocks: &[(u64, u64)],
    introns: &[(u64, u64)],
    lo: u64,
    hi: u64,
    support: &BTreeMap<(u64, u64), usize>,
    min_reads: usize,
) -> Option<Vec<(u64, u64)>> {
    let mut segs: Vec<Vec<(u64, u64)>> = vec![vec![*blocks.first()?]];
    for (k, intr) in introns.iter().enumerate() {
        if intr.1 - intr.0 > GIANT_BP && support.get(intr).copied().unwrap_or(0) < min_reads {
            segs.push(Vec::new());
        }
        if k + 1 < blocks.len() {
            segs.last_mut().unwrap().push(blocks[k + 1]);
        }
    }
    segs.into_iter().find(|sg| sg.iter().any(|&(s, e)| e > lo && s < hi))
}

/// ⭐ L2: the read-supported EXTENT of a unit — the union, over every PRIMARY record (`-F 2308`, the project
/// invariant) with an aligned block inside the unit's exon chain, of the record's kept segment under the
/// mis-chain rule above. Two definitions were refuted on the way (3-contig, unit median 21 kb): every record's
/// reference span gave a median extent of 548 kb (records spliced over the locus through giant introns), and
/// every record's kept segment 142 kb (MAPQ-0 secondaries with 400–800 blocks hopping through a ZNF array).
/// Primaries with a block in the chain are the molecules O2 aligns for this copy. This is O2's alignment
/// target (`copies.tsv` `locus_start`/`locus_end`), replacing the padding rule O2 invented in §6fd. `None`
/// when no primary record has a block in the chain.
fn read_extent(reads: &[rustle::vg_family::denovo_assemble::BamRead], chain: &[(u64, u64)], min_reads: usize) -> Option<(u64, u64)> {
    let (lo, hi) = (chain.first()?.0, chain.last()?.1);
    let parsed: Vec<(Vec<(u64, u64)>, Vec<(u64, u64)>)> =
        reads.iter().filter(|br| !br.is_supplementary && !br.is_secondary).map(blocks_and_introns).collect();
    let mut support: BTreeMap<(u64, u64), usize> = BTreeMap::new();
    for (_, introns) in &parsed {
        for &i in introns {
            *support.entry(i).or_insert(0) += 1;
        }
    }
    parsed
        .iter()
        .filter_map(|(b, i)| kept_segment(b, i, lo, hi, &support, min_reads))
        .filter(|sg| sg.iter().any(|&(bs, be)| chain.iter().any(|&(s, e)| be > s && bs < e)))
        .map(|sg| (sg[0].0, sg.last().unwrap().1))
        .reduce(|(a, b), (c, d)| (a.min(c), b.max(d)))
}

/// ⭐ L2: a locus never contains another catalog unit — of ANY family (§6fm). `units[i] = (chain start, chain end,
/// extent)` on ONE contig; the extent is clipped to the nearest chain ends of the other units that lie entirely outside its own
/// span (units overlapping the span — nested, interleaved — do not clip). Without this, in a tandem array a
/// single read-through molecule extends copy 1's locus over copy 2's chain, a copy-2 read then aligns identically
/// inside both targets and ties with zero decisive columns (MCL4: 226 assigned → tied, arm A0 of PREREG L1/L2).
fn clip_extents_to_neighbours(units: &[(u64, u64, (u64, u64))]) -> Vec<(u64, u64)> {
    units
        .iter()
        .enumerate()
        .map(|(i, &(s, e, (a, b)))| {
            let (mut lo, mut hi) = (a.min(s), b.max(e));
            for (j, &(s2, e2, _)) in units.iter().enumerate() {
                if j == i {
                    continue;
                }
                if e2 <= s {
                    lo = lo.max(e2);
                } else if s2 >= e {
                    hi = hi.min(s2);
                }
            }
            (lo, hi)
        })
        .collect()
}

/// ⭐ The read-supported exon chain of one locus `[lo, hi)` (§6el rule; one shipped constant pair: 50 kb / 3).
/// Returns `(blocks, strand, n_reads)`; `blocks` empty when fewer than `min_reads` reads or no base reaches it.
fn read_chain(
    reads: &[rustle::vg_family::denovo_assemble::BamRead],
    lo: u64,
    hi: u64,
    min_reads: usize,
    follow: Option<(u64, u64)>,
) -> (Vec<(u64, u64)>, Option<char>, usize) {
    let parsed: Vec<_> = reads
        .iter()
        .filter(|br| !br.is_supplementary && !br.is_secondary)
        .map(|br| {
            let (b, i) = blocks_and_introns(br);
            let strand = match br.ts {
                Some(t) => {
                    if br.reverse {
                        if t == '+' { '-' } else { '+' }
                    } else {
                        t
                    }
                }
                None => {
                    if br.reverse { '-' } else { '+' }
                }
            };
            (b, i, strand)
        })
        .filter(|(b, _, _)| b.iter().any(|&(s, e)| e > lo && s < hi))
        .collect();
    let n_reads = parsed.len();
    if n_reads < min_reads || hi <= lo {
        return (Vec::new(), None, n_reads);
    }
    let mut support: BTreeMap<(u64, u64), usize> = BTreeMap::new();
    for (_, introns, _) in &parsed {
        for &i in introns {
            *support.entry(i).or_insert(0) += 1;
        }
    }
    // the segment of each read that overlaps the window, after the mis-chain cut (split at giant,
    // unsupported introns); `None` = the read's kept segment does not touch the window.
    let mut kept: Vec<Option<Vec<(u64, u64)>>> = Vec::with_capacity(parsed.len());
    let mut strands: BTreeMap<char, usize> = BTreeMap::new();
    for (blocks, introns, strand) in &parsed {
        kept.push(kept_segment(blocks, introns, lo, hi, &support, min_reads));
        *strands.entry(*strand).or_insert(0) += 1;
    }
    // the coverage window: the locus window itself, or (units follow the reads, `follow = Some(bound)`) the
    // extent of every kept segment CLAMPED to `bound` = the member's annotated locus span — a block outside
    // the window but inside the annotation is then judged by the same support rule instead of being clipped.
    // ⚠ Unbounded following ENGULFS neighbours: read-through molecules (≥3) chained MCL7's 32-kb kept-full
    // unit into a 133-kb unit over 8 genes (CDR2 among them) and lifted the family's "unit reads" 659 → 9,098
    // (adj/worst2, rna_units_v3_unbounded). The annotation bounds the locus; the reads shape it inside.
    let (wlo, whi) = match follow {
        Some((blo, bhi)) => {
            let (a, b) = kept.iter().flatten().flatten().fold((lo, hi), |(a, b), &(s, e)| (a.min(s), b.max(e)));
            (a.max(blo.min(lo)), b.min(bhi.max(hi)))
        }
        None => (lo, hi),
    };
    let mut cov = vec![0u32; (whi - wlo) as usize];
    let mut spliced = vec![0u32; (whi - wlo) as usize]; // reads whose intron spans the base (they vote AGAINST exon)
    for sg in kept.iter().flatten() {
        for &(s, e) in sg {
            for x in s.max(wlo)..e.min(whi) {
                cov[(x - wlo) as usize] += 1;
            }
        }
        // introns INSIDE the kept segment splice over their bases
        for w in sg.windows(2) {
            for x in w[0].1.max(wlo)..w[1].0.min(whi) {
                spliced[(x - wlo) as usize] += 1;
            }
        }
    }
    let mut blocks: Vec<(u64, u64)> = Vec::new();
    for (k, &c) in cov.iter().enumerate() {
        // exonic iff covered by >= min_reads blocks AND by more blocks than reads splicing over it
        if c as usize >= min_reads && c > spliced[k] {
            let x = wlo + k as u64;
            match blocks.last_mut() {
                Some(b) if b.1 == x => b.1 = x + 1,
                _ => blocks.push((x, x + 1)),
            }
        }
    }
    // majority strand; ties -> '+' (python `Counter.most_common` order-dependence removed on purpose)
    let strand = strands.iter().max_by_key(|(c, n)| (**n, if **c == '+' { 1 } else { 0 })).map(|(c, _)| *c);
    (blocks, strand, n_reads)
}

/// Introns of a read that LEAVE the emitted chain: the read must have an aligned block inside `chain`
/// (`-F 2308` is the caller's job) and the intron's far end must lie beyond the chain's last exon.
/// These are the candidate read-through junctions; the far end is resolved against the catalog later.
fn leaving_introns(blocks: &[(u64, u64)], introns: &[(u64, u64)], chain: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if chain.is_empty() {
        return Vec::new();
    }
    let (clo, chi) = (chain[0].0, chain.last().unwrap().1);
    let _ = clo;
    if !blocks.iter().any(|&(bs, be)| chain.iter().any(|&(s, e)| be > s && bs < e)) {
        return Vec::new();
    }
    // only the DOWNSTREAM direction: the mirrored case is found when the other unit is the source.
    // ⭐ §6ge: the DONOR must be this chain's own splice donor — the intron must start at the end of one of
    // its exons. Without that, a read merely ANCHORED in the chain and splicing anywhere downstream counted,
    // so an ordinary intron of the NEXT gene was reported as a read-through out of this one (32 of 42 guarded
    // junctions were annotated introns). The acceptor is checked separately by `readthrough_target`.
    introns
        .iter()
        .copied()
        .filter(|&(_, e)| e > chi)
        .filter(|&(s, _)| chain.iter().any(|&(_, ce)| ce == s))
        .collect()
}

/// The emitted unit whose chain CONTAINS the first base after the intron (`intron.1`), i.e. where the
/// read-through lands. `units` is `(contig, exon chain)` per candidate unit; the source unit never matches.
fn readthrough_target(intron: (u64, u64), contig: &str, src: usize, units: &[(String, Vec<(u64, u64)>)]) -> Option<usize> {
    units.iter().enumerate().position(|(k, (c, ch))| {
        k != src && c == contig && ch.iter().any(|&(s, e)| intron.1 >= s && intron.1 < e)
    })
}

/// A canonical intron on `strand`: `GT..AG`, `GC..AG` or `AT..AC` read 5'→3' on the transcript.
/// On `-` the genomic dinucleotides are the reverse complements (`CT..AC`, `CT..GC`, `GT..AT`).
fn canonical_intron(genome: &rustle::genome::GenomeIndex, contig: &str, s: u64, e: u64, strand: char) -> bool {
    if e < s + 4 {
        return false;
    }
    let (Some(d), Some(a)) = (genome.fetch_sequence(contig, s, s + 2), genome.fetch_sequence(contig, e - 2, e)) else {
        return false;
    };
    let up = |v: &[u8]| -> String { String::from_utf8_lossy(v).to_uppercase() };
    let (d, a) = (up(&d), up(&a));
    if strand == '-' {
        matches!((d.as_str(), a.as_str()), ("CT", "AC") | ("CT", "GC") | ("GT", "AT"))
    } else {
        matches!((d.as_str(), a.as_str()), ("GT", "AG") | ("GC", "AG") | ("AT", "AC"))
    }
}

use rustle::vg_family::denovo_assemble::longest_orf;

fn lengths_from_blocks(b: &BTreeMap<GeneKey, Vec<(u64, u64)>>) -> BTreeMap<GeneKey, u64> {
    b.iter()
        .map(|(g, v)| (g.clone(), v.iter().map(|&(s, e)| e - s).sum::<u64>().max(1)))
        .collect()
}

/// Does this read have at least one ALIGNED BLOCK inside `[start-1, end)`?
///
/// ⚠ **NOT span overlap.** `N` in an RNA CIGAR is an intron, spliced OUT: a read that splices straight
/// OVER a locus contributes no aligned base and is no evidence for it. Counting span overlap instead
/// inflated a headline 3.4x and produced a retracted mechanism (ledger §6cm) — 71.6% of the reads
/// "supporting" one locus were merely passing through it. Supplementary records are excluded by the
/// caller so one molecule is never two witnesses.
fn has_block_in(br: &rustle::vg_family::denovo_assemble::BamRead, m: &GeneKey) -> bool {
    let (lo, hi) = (m.1.saturating_sub(1), m.2);
    let mut p = br.read.ref_start;
    let mut cur: Option<(u64, u64)> = None;
    for &(op, n) in &br.read.cigar {
        match op {
            'M' | '=' | 'X' | 'D' => {
                cur = Some((cur.map_or(p, |c| c.0), p + n));
                p += n;
            }
            'N' => {
                if let Some((s, e)) = cur.take() {
                    if s < hi && lo < e {
                        return true;
                    }
                }
                p += n;
            }
            _ => {}
        }
    }
    matches!(cur, Some((s, e)) if s < hi && lo < e)
}

/// A RepeatMasker `.out` (curated library) as sorted interspersed-repeat intervals per contig (0-based half-open).
fn load_rmsk(path: &str) -> Result<BTreeMap<String, Vec<(u64, u64)>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut m: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 11 || f[0].parse::<u64>().is_err() {
            continue;
        }
        let class = f[10].split('/').next().unwrap_or("");
        if !matches!(class, "LINE" | "SINE" | "LTR" | "Retroposon" | "DNA" | "RC" | "Unknown") {
            continue;
        }
        if let (Ok(a), Ok(b)) = (f[5].parse::<u64>(), f[6].parse::<u64>()) {
            m.entry(f[4].to_string()).or_default().push((a - 1, b));
        }
    }
    for v in m.values_mut() {
        v.sort_unstable();
    }
    Ok(m)
}

/// Interspersed-repeat fraction of an exon chain (`None` when the contig has no RepeatMasker interval at all).
fn rep_frac_in(rmsk: &BTreeMap<String, Vec<(u64, u64)>>, chrom: &str, exons: &[(u64, u64)]) -> Option<f64> {
    let v = rmsk.get(chrom)?;
    let (mut tot, mut inter) = (0u64, 0u64);
    for &(s, e) in exons {
        tot += e - s;
        let i = v.partition_point(|x| x.1 <= s);
        for &(a, b) in &v[i..] {
            if a >= e {
                break;
            }
            inter += b.min(e).saturating_sub(a.max(s));
        }
    }
    Some(inter as f64 / tot.max(1) as f64)
}

/// Header of `<out>.copies.tsv` (`--from-gtf --emit-units`). Columns 1-11 are `gw_family_catalog`'s `copies.tsv`
/// header, in its order, so `copy_assign --families` (parsed by name), `bench/sim.py copies` (positional 1-9) and
/// every reader of the legacy catalog read it unchanged; the rest are appended.
const COPIES_HEADER: &str = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\
     \tsource\tgene_id\tcore_hull\tsd_depth\tcore_bp\trep_frac\tmember_status\tlocus_start\tlocus_end";

/// Counts of one `write_locus_rep_copies` call (the params certificate and the log).
#[derive(Debug, Default, Clone, PartialEq)]
struct RepCopyStats {
    /// rows written to `copies.tsv`
    copies: usize,
    /// families with >= 1 copy / with >= 2 copies
    families: usize,
    multi_copy_families: usize,
    /// members folded into an exon-overlapping copy of the same family (`copies.merged.tsv`)
    merged: usize,
    /// members skipped because the core rule dropped them and `--no-units-include-dropped` was given
    skipped_dropped: usize,
    dropped_emitted: usize,
    noncoding: usize,
    /// copies whose representative has `reads 0` (or no `reads` attribute)
    unexpressed: usize,
    /// representatives with strand `.`, written as `+` (the copies contract has no unstranded copy)
    unstranded: usize,
    /// loci sharing one `(chrom, start, end)` with another `gene_id` (one graph node; the representative with more
    /// reads is the copy)
    key_collisions: usize,
    /// representatives whose exons overlapped or abutted and were coalesced into one block
    coalesced: usize,
}

/// Genome bases at `exons` (0-based half-open, ascending), concatenated, reverse-complemented on `-`: a copy's
/// spliced sequence in transcription orientation, as `gw_family_catalog` writes `copies.fa`.
fn spliced_exon_sum(genome: &rustle::genome::GenomeIndex, chrom: &str, exons: &[(u64, u64)], strand: char) -> Result<Vec<u8>> {
    let mut seq: Vec<u8> = Vec::new();
    for &(s, e) in exons {
        let part = genome
            .fetch_sequence(chrom, s, e)
            .with_context(|| format!("copies: {chrom}:{s}-{e} is not in --fasta"))?;
        seq.extend_from_slice(&part);
    }
    if strand == '-' {
        seq = rustle::vg_family::seq_utils::revcomp_keep_case(&seq);
    }
    Ok(seq)
}

/// ⭐ `--from-gtf --emit-units` (user decision 2026-09-25 16:00: ONE default de novo family definition, and copy
/// assignment consumes the SAME families): write the families as a COPY TABLE in the `copy_assign --families /
/// --copies-fa` contract. One copy per member locus of every reported cluster (the rows of `clusters.tsv`,
/// `family_id` = its `cluster_id`); the copy IS the locus representative — its exons (the positional exon sum),
/// its spliced sequence, its `transcript_id` (`tid`, which joins back to the assembled GTF), its `reads`
/// (`n_reads`) and its strand. No BAM is read: the de novo loci already are read-derived.
///
/// `max_family_identity` = the identity of the best families-stage edge (genomic-span `-x asm20` alignment, the
/// one the admission rule scored) from this locus to another member, `NA` when it has no direct edge. ⚠ The legacy
/// catalog's column of that name is an exon-sum alignment identity: same role, different alignment.
/// `locus_start`/`locus_end` = the de novo locus span (all its transcripts), clipped at the exon-chain ends of the
/// neighbouring copies on the contig (the L2 rule, `clip_extents_to_neighbours`).
/// Also applied as for read-chain units: `--no-units-include-dropped`, `--coding-core`, and the §6fb merge of
/// copies of one family that share exon bases (a no-op after `--fold-within-clusters`, kept as the guarantee that
/// no base belongs to two copies of one family). Writes `<out>.copies.tsv/.fa/.regions/.merged.tsv`.
#[allow(clippy::too_many_arguments)]
fn write_locus_rep_copies(
    out: &str,
    fasta: &str,
    clusters: &[Cluster],
    g: &rustle::vg_family::annotation_families::HomologyGraph,
    core_records: &[Vec<rustle::vg_family::annotation_families::CoreRecord>],
    loci: &[GtfLocus],
    rmsk: Option<&BTreeMap<String, Vec<(u64, u64)>>>,
    include_dropped: bool,
    merge_overlapping: bool,
    coding_core: bool,
) -> Result<RepCopyStats> {
    struct RepCopy {
        member: GeneKey,
        gene_id: String,
        tid: String,
        reads: u64,
        strand: char,
        exons: Vec<(u64, u64)>,
        seq: Vec<u8>,
        ident: Option<f64>,
        hull_col: String,
        sd_depth: String,
        core_bp: String,
        status: &'static str,
        orf: usize,
        locus: (u64, u64),
    }
    let mut st = RepCopyStats::default();
    // representative per graph node (`loci.gff3` gene line = the node key, GFF 1-based)
    let mut rep_of: BTreeMap<GeneKey, &GtfLocus> = BTreeMap::new();
    for l in loci {
        let k: GeneKey = (l.chrom.clone(), l.start, l.end);
        match rep_of.get(&k) {
            Some(prev) => {
                st.key_collisions += 1;
                if l.rep_reads > prev.rep_reads {
                    rep_of.insert(k, l);
                }
            }
            None => {
                rep_of.insert(k, l);
            }
        }
    }
    let contigs: std::collections::HashSet<String> =
        clusters.iter().flat_map(|c| c.members.iter().map(|m| m.0.clone())).collect();
    // ⚠ `from_fasta_contigs` with an EMPTY set loads the whole genome: no family, no genome
    let genome = if contigs.is_empty() {
        rustle::genome::GenomeIndex::empty()
    } else {
        rustle::genome::GenomeIndex::from_fasta_contigs(fasta, &contigs)?
    };
    let node_idx: BTreeMap<&GeneKey, usize> = g.genes.iter().enumerate().map(|(k, gk)| (gk, k)).collect();
    let mut staged: Vec<(String, Vec<RepCopy>, Vec<Option<usize>>)> = Vec::new();
    for (i, c) in clusters.iter().enumerate() {
        let fid = format!("MCL{i}");
        let mut pending: Vec<RepCopy> = Vec::new();
        for (mi, m) in c.members.iter().enumerate() {
            let rec = core_records.get(i).and_then(|v| v.get(mi));
            let status: &'static str = match rec.map(|r| r.status) {
                Some(CoreStatus::Dropped) => "dropped",
                Some(CoreStatus::KeptTrimmed) => "kept_trimmed",
                Some(_) => "kept_full",
                None => "ungated",
            };
            if status == "dropped" && !include_dropped {
                st.skipped_dropped += 1;
                continue;
            }
            let l = rep_of.get(m).with_context(|| {
                format!("copies: family member {}:{}-{} is not a locus of the --from-gtf GTF (a failed join)", m.0, m.1, m.2)
            })?;
            // the representative's exons, 0-based half-open; overlapping or abutting blocks coalesced
            let mut exons: Vec<(u64, u64)> = Vec::with_capacity(l.rep_exons.len());
            let mut coalesced = false;
            for &(_, a, b) in &l.rep_exons {
                let (s, e) = (a.saturating_sub(1), b);
                match exons.last_mut() {
                    Some(p) if s <= p.1 => {
                        p.1 = p.1.max(e);
                        coalesced = true;
                    }
                    _ => exons.push((s, e)),
                }
            }
            anyhow::ensure!(!exons.is_empty(), "copies: locus {} has a representative without exons", l.gene_id);
            if coalesced {
                st.coalesced += 1;
            }
            let strand = match l.strand.as_str() {
                "+" => '+',
                "-" => '-',
                _ => {
                    st.unstranded += 1;
                    '+'
                }
            };
            let seq = spliced_exon_sum(&genome, &m.0, &exons, strand)?;
            let ident = node_idx.get(m).and_then(|&a| {
                c.members
                    .iter()
                    .filter_map(|o| node_idx.get(o).copied())
                    .filter(|&b| b != a)
                    .filter_map(|b| g.idents.get(&(a.min(b), a.max(b))).copied())
                    .reduce(f64::max)
            });
            let hull_col = match rec.and_then(|r| r.hull) {
                Some((a, b)) => format!("{}-{}", a.saturating_sub(1), b),
                None => "NA".to_string(),
            };
            let (sd_depth, core_bp) =
                rec.map(|r| (r.max_depth.to_string(), r.core_bp.to_string())).unwrap_or_else(|| ("NA".into(), "NA".into()));
            let (us, ue) = (exons[0].0, exons.last().unwrap().1);
            if status == "dropped" {
                st.dropped_emitted += 1;
            }
            pending.push(RepCopy {
                member: m.clone(),
                gene_id: l.gene_id.clone(),
                tid: l.rep.clone(),
                reads: l.rep_reads,
                strand,
                orf: if coding_core { longest_orf(&seq) } else { 0 },
                exons,
                seq,
                ident,
                hull_col,
                sd_depth,
                core_bp,
                status,
                locus: (m.1.saturating_sub(1).min(us), m.2.max(ue)),
            });
        }
        // copies in genomic order within the family (the catalog's order)
        pending.sort_by(|a, b| {
            (a.member.0.as_str(), a.exons[0].0, a.exons.last().unwrap().1, a.tid.as_str())
                .cmp(&(b.member.0.as_str(), b.exons[0].0, b.exons.last().unwrap().1, b.tid.as_str()))
        });
        // CODING CORE: the family's best ORF is the reference (see `--coding-core`)
        if coding_core {
            let best_orf = pending.iter().filter(|u| u.status != "dropped").map(|u| u.orf).max().unwrap_or(0);
            if best_orf > 0 {
                for u in pending.iter_mut() {
                    if u.status != "dropped" && u.orf * 2 < best_orf {
                        u.status = "noncoding";
                        st.noncoding += 1;
                    }
                }
            }
        }
        // §6fb: copies of one family that share exon bases are one copy (kept before dropped, then the longest
        // exon union represents them)
        let n = pending.len();
        let mut merged_into: Vec<Option<usize>> = vec![None; n];
        if merge_overlapping {
            let mut parent: Vec<usize> = (0..n).collect();
            fn find(p: &mut [usize], mut x: usize) -> usize {
                while p[x] != x {
                    p[x] = p[p[x]];
                    x = p[x];
                }
                x
            }
            for a in 0..n {
                for b in (a + 1)..n {
                    if pending[a].member.0 != pending[b].member.0 {
                        continue;
                    }
                    let share = pending[a].exons.iter().any(|&(s1, e1)| pending[b].exons.iter().any(|&(s2, e2)| s1 < e2 && s2 < e1));
                    if share {
                        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                        if ra != rb {
                            parent[ra.max(rb)] = ra.min(rb);
                        }
                    }
                }
            }
            let rank = |u: &RepCopy| (u.status != "dropped", u.exons.iter().map(|(s, e)| e - s).sum::<u64>());
            let mut rep_of_root: BTreeMap<usize, usize> = BTreeMap::new();
            for k in 0..n {
                let r = find(&mut parent, k);
                let e = rep_of_root.entry(r).or_insert(k);
                if rank(&pending[k]) > rank(&pending[*e]) {
                    *e = k;
                }
            }
            for k in 0..n {
                let rep = rep_of_root[&find(&mut parent, k)];
                if rep != k {
                    merged_into[k] = Some(rep);
                }
            }
        }
        staged.push((fid, pending, merged_into));
    }
    // L2: every copy's locus extent stops at the exon-chain ends of the copies (of any family) around it
    let mut clipped: Vec<Vec<(u64, u64)>> = staged.iter().map(|(_, p, _)| p.iter().map(|u| u.locus).collect()).collect();
    {
        let mut by_ctg: BTreeMap<&str, Vec<(usize, usize)>> = BTreeMap::new();
        for (fi, (_, pending, merged_into)) in staged.iter().enumerate() {
            for (k, u) in pending.iter().enumerate() {
                if merged_into[k].is_none() {
                    by_ctg.entry(u.member.0.as_str()).or_default().push((fi, k));
                }
            }
        }
        for ks in by_ctg.values() {
            let spans: Vec<(u64, u64, (u64, u64))> = ks
                .iter()
                .map(|&(fi, k)| {
                    let u = &staged[fi].1[k];
                    (u.exons[0].0, u.exons.last().unwrap().1, u.locus)
                })
                .collect();
            for (&(fi, k), c) in ks.iter().zip(clip_extents_to_neighbours(&spans)) {
                clipped[fi][k] = c;
            }
        }
    }
    let mut ct = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.tsv"))?);
    let mut cf = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.fa"))?);
    let mut cr = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.regions"))?);
    let mut cm = std::io::BufWriter::new(std::fs::File::create(format!("{out}.copies.merged.tsv"))?);
    writeln!(ct, "{COPIES_HEADER}")?;
    writeln!(cm, "family_id\tmerged_member\tmerged_tid\tinto_member\tinto_tid")?;
    for (fi, (fid, pending, merged_into)) in staged.iter().enumerate() {
        let mut idx = 0usize;
        let mut hulls: BTreeMap<&str, (u64, u64)> = BTreeMap::new();
        for (k, u) in pending.iter().enumerate() {
            if let Some(rep) = merged_into[k] {
                let r = &pending[rep];
                writeln!(
                    cm,
                    "{fid}\t{}:{}-{}\t{}\t{}:{}-{}\t{}",
                    u.member.0, u.member.1, u.member.2, u.tid, r.member.0, r.member.1, r.member.2, r.tid
                )?;
                st.merged += 1;
                continue;
            }
            let chrom = u.member.0.as_str();
            let (us, ue) = (u.exons[0].0, u.exons.last().unwrap().1);
            let rep_col = rmsk
                .and_then(|r| rep_frac_in(r, chrom, &u.exons))
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "NA".into());
            writeln!(
                ct,
                "{fid}\t{idx}\t{}\t{chrom}\t{us}\t{ue}\t{}\t{}\t{}\t{}\t{}\tlocus_rep\t{}\t{}\t{}\t{}\t{rep_col}\t{}\t{}\t{}",
                u.tid,
                u.exons.len(),
                u.strand,
                u.reads,
                u.exons.iter().map(|(s, e)| format!("{s}-{e}")).collect::<Vec<_>>().join(","),
                u.ident.map(|v| format!("{v:.6}")).unwrap_or_else(|| "NA".into()),
                u.gene_id,
                u.hull_col,
                u.sd_depth,
                u.core_bp,
                u.status,
                clipped[fi][k].0,
                clipped[fi][k].1
            )?;
            writeln!(cf, ">{fid}|{idx}|{chrom}:{us}-{ue}|{}|nexon={}", u.strand, u.exons.len())?;
            cf.write_all(&u.seq)?;
            writeln!(cf)?;
            let h = hulls.entry(chrom).or_insert((us, ue));
            h.0 = h.0.min(us);
            h.1 = h.1.max(ue);
            if u.reads == 0 {
                st.unexpressed += 1;
            }
            st.copies += 1;
            idx += 1;
        }
        for (ctg, (a, b)) in hulls {
            writeln!(cr, "{fid}\t{ctg}:{}-{}", a.saturating_sub(5_000).max(1), b + 5_000)?;
        }
        if idx >= 1 {
            st.families += 1;
        }
        if idx >= 2 {
            st.multi_copy_families += 1;
        }
    }
    for w in [&mut ct, &mut cf, &mut cr, &mut cm] {
        w.flush()?;
    }
    Ok(st)
}

fn main() -> Result<()> {
    let mut args = Args::parse();
    anyhow::ensure!(
        !args.emit_container || args.from_gtf.is_some(),
        "--emit-container needs --from-gtf (the container's exon blocks are the assembled GTF's transcripts)"
    );
    anyhow::ensure!(
        args.representative == Representative::MostReads || args.from_gtf.is_some(),
        "--representative most-junctions needs --from-gtf (it picks the representative transcript of each de novo locus)"
    );
    // `--from-gtf`: the de novo loci (their representatives become the copy table under `--emit-units`)
    let mut gtf_loci_list: Option<Vec<GtfLocus>> = None;
    if let Some(gtf) = args.from_gtf.clone() {
        let fasta = args.fasta.clone().context("--from-gtf needs --fasta (the genome the GTF was assembled on)")?;
        let (gff3, fa, paf, loci) = loci_from_gtf(&gtf, &fasta, &args.out, args.threads, args.representative)?;
        args.gff = Some(gff3);
        args.paf = paf;
        gtf_loci_list = Some(loci);
        eprintln!("[mcl_families] --from-gtf: loci in {fa}");
    }
    anyhow::ensure!(!args.paf.is_empty(), "--paf is required unless --from-gtf is given");
    // §6er (S2): the unit is the catalog row. Stages engage on their inputs; escape hatches reproduce the
    // record-level catalogs (`--no-merge-overlapping-loci`, `--no-core-refine`, `--no-emit-units`).
    if args.no_merge_overlapping_loci {
        args.merge_overlapping_loci = false;
    }
    if (args.sedef.is_some() || args.core_from_paf) && !args.no_core_refine {
        args.core_refine = true;
    }
    if args.no_fold_within_clusters {
        args.fold_within_clusters = false;
    }
    if args.no_exonless_span {
        args.exonless_span = false;
    }
    if args.no_exonic_both_sides {
        args.exonic_both_sides = false;
    }
    if args.no_units_follow_reads {
        args.units_follow_reads = false;
    }
    if args.no_units_include_dropped {
        args.units_include_dropped = false;
    }
    if args.bam.is_some() && args.fasta.is_some() && !args.no_emit_units {
        args.emit_units = true;
    }
    let p = GraphParams {
        min_identity: args.min_identity,
        min_cov_longer: args.min_cov_longer,
        min_cov_shorter: args.min_cov_shorter,
        min_bp: args.min_bp,
        exonic_overlap: args.exonic_overlap,
        reject_overlapping: args.reject_overlapping,
        min_exonic_bp: args.min_exonic_bp,
        exonic_both_sides: args.exonic_both_sides,
        min_shared_exon_frac: args.min_shared_exon_frac,
    };

    let blocks = match &args.gff {
        Some(g) => exonic_blocks(g, args.exonless_span)?,
        None => {
            eprintln!(
                "[mcl_families] WARNING: no --gff, so the coverage denominator is the GENOMIC SPAN. The \
                 median gene is ~23% exon, and on the pilot this removed 62.8% of eligible genes (§6dc)."
            );
            BTreeMap::new()
        }
    };
    let exonic = lengths_from_blocks(&blocks);
    eprintln!("[mcl_families] exon-union lengths for {} genes", exonic.len());

    let paf = std::fs::read_to_string(&args.paf).with_context(|| format!("reading {}", args.paf))?;
    let loci = if args.merge_overlapping_loci && !args.fold_within_clusters {
        let mut m = loci_from_exon_blocks(&blocks);
        m.attribute_edges = args.locus_attribute_edges;
        eprintln!(
            "[mcl_families] merge-overlapping-loci (EXON-union overlap): {} annotation record(s) folded \
             into {} multi-record loci over {} GFF genes",
            m.n_merged(),
            m.n_multi,
            blocks.len()
        );
        let mut lf = std::fs::File::create(format!("{}.loci.tsv", args.out))?;
        writeln!(lf, "annotation\trepresentative")?;
        for (k, r) in &m.rep_of {
            if k != r {
                writeln!(lf, "{}:{}-{}\t{}:{}-{}", k.0, k.1, k.2, r.0, r.1, r.2)?;
            }
        }
        Some(m)
    } else {
        None
    };
    let g = graph_from_paf_loci(&paf, &exonic, &blocks, &p, loci.as_ref());
    if let Some(path) = &args.dump_graph {
        let mut f = std::fs::File::create(path)?;
        for (&(i, j), &w) in &g.edges {
            let (ci, si, ei) = &g.genes[i];
            let (cj, sj, ej) = &g.genes[j];
            writeln!(f, "{ci}:{si}-{ei}\t{cj}:{sj}-{ej}\t{w:.6}")?;
        }
        eprintln!(
            "[mcl_families] --dump-graph: {} node(s), {} edge(s) -> {path}",
            g.n_nodes(),
            g.n_edges()
        );
    }
    if let Some(m) = &loci {
        eprintln!(
            "[mcl_families] merge-overlapping-loci ({}): {} PAF record(s) skipped ({} loci with >1 record)",
            if m.attribute_edges { "attribution" } else { "representative-only" },
            g.same_locus_records,
            m.n_multi
        );
    }
    // ⭐ The join rate, always. A failed coordinate join returns a byte-identical graph and reads as
    // "the fix is inert" — it must be visible, not inferred.
    eprintln!(
        "[mcl_families] graph: {} nodes, {} edges (identity>={}, cov_longer>={}, >={}bp); \
         {} node(s) had NO exon-union length and fell back to the span",
        g.n_nodes(),
        g.n_edges(),
        p.min_identity,
        p.min_cov_longer,
        p.min_bp,
        g.missing_exonic
    );
    if p.exonic_overlap {
        eprintln!(
            "[mcl_families] exonic-overlap numerator: {} edge(s) joined exon blocks, {} fell back to \
             the span numerator",
            g.exonic_overlap_joined, g.exonic_overlap_missing
        );
        if g.exonic_overlap_joined == 0 {
            eprintln!(
                "[mcl_families] WARNING: 0 edges joined exon blocks — --exonic-overlap did NOTHING. \
                 That is the §6dd coordinate bug, not an inert clause."
            );
        }
    }
    if p.min_exonic_bp > 0 {
        eprintln!(
            "[mcl_families] min-exonic-bp={}: {} pair(s) dropped for resting on no exonic evidence",
            p.min_exonic_bp, g.rejected_no_exonic
        );
    }
    if p.min_shared_exon_frac > 0.0 {
        eprintln!(
            "[mcl_families] min-shared-exon-frac={}: {} pair(s) dropped for sharing too small a fraction of exons",
            p.min_shared_exon_frac, g.rejected_low_shared_exon
        );
    }
    if p.reject_overlapping {
        eprintln!(
            "[mcl_families] reject-overlapping: {} pair(s) dropped as overlapping annotation intervals",
            g.rejected_overlapping
        );
    }
    if !exonic.is_empty() && g.missing_exonic * 2 > g.n_nodes() {
        eprintln!(
            "[mcl_families] WARNING: {}/{} nodes fell back to the span — this is what a FAILED \
             COORDINATE JOIN looks like. The headers must be GFF 1-based verbatim.",
            g.missing_exonic,
            g.n_nodes()
        );
    }

    let parts = mcl(&g, args.inflation, 100, args.prune);
    // ⭐ §6ey: fold overlapping records into loci INSIDE each cluster (sequence decided the cluster, coordinates
    // decide the locus). `parts` becomes representative-only; the fold map is written like the fold-first one.
    let (parts, loci) = if args.fold_within_clusters {
        let (folded, m) = fold_parts_into_loci(&g, &parts, &blocks);
        eprintln!(
            "[mcl_families] fold-within-clusters: {} annotation record(s) folded into {} multi-record loci inside their clusters",
            m.n_merged(),
            m.n_multi
        );
        let mut lf = std::fs::File::create(format!("{}.loci.tsv", args.out))?;
        writeln!(lf, "annotation\trepresentative")?;
        for (k, r) in &m.rep_of {
            if k != r {
                writeln!(lf, "{}:{}-{}\t{}:{}-{}", k.0, k.1, k.2, r.0, r.1, r.2)?;
            }
        }
        (folded, Some(m))
    } else {
        (parts, loci)
    };

    // RNA corroboration, when a BAM is supplied.
    let corr_set: Option<BTreeSet<GeneKey>> = match &args.bam {
        None => None,
        Some(bam) => {
            let members: BTreeSet<GeneKey> = parts
                .iter()
                .filter(|p| p.len() >= args.min_size)
                .flat_map(|p| p.iter().map(|&i| g.genes[i].clone()))
                .collect();
            eprintln!("[mcl_families] corroborating {} member(s) against {bam}", members.len());
            let mut ok = BTreeSet::new();
            for m in &members {
                // `start` is GFF 1-based; the BAM query wants 0-based half-open.
                let (_, reads) = rustle::vg_family::denovo_assemble::reads_in_region(
                    bam,
                    &m.0,
                    m.1.saturating_sub(1),
                    m.2,
                    1,
                )
                .with_context(|| format!("reading {}:{}-{}", m.0, m.1, m.2))?;
                let n = reads.iter().filter(|br| !br.is_supplementary && has_block_in(br, m)).count();
                if n >= args.min_reads {
                    ok.insert(m.clone());
                }
            }
            Some(ok)
        }
    };
    let pred = corr_set.as_ref().map(|s| {
        let f = move |k: &GeneKey| s.contains(k);
        Box::new(f) as Box<dyn Fn(&GeneKey) -> bool>
    });
    let clusters: Vec<Cluster> =
        build_clusters(&g, &parts, args.min_size, pred.as_ref().map(|b| b.as_ref()));

    let mut ch = std::fs::File::create(format!("{}.clusters.tsv", args.out))?;
    writeln!(ch, "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend")?;
    for (i, c) in clusters.iter().enumerate() {
        let corr = c.corroborated.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".into());
        for m in &c.members {
            writeln!(
                ch,
                "MCL{i}\t{}\t{:.4}\t{:.4}\t{corr}\t{}\t{}\t{}",
                c.members.len(),
                c.density,
                c.frac_in,
                m.0,
                m.1,
                m.2
            )?;
        }
    }

    // ⭐ PREREG (identity-weighted density, 09-10): `--dump-pairs` — every WITHIN-cluster edge's identity,
    // straight from `g.idents` (the same aggregated Σnmatch/Σblocklen value the admission rule computed),
    // so a downstream density-weighting or identity_gap.py run scores exactly the shipped graph.
    if args.dump_pairs {
        let node_idx: std::collections::BTreeMap<&GeneKey, usize> =
            g.genes.iter().enumerate().map(|(k, gk)| (gk, k)).collect();
        let mut pf = std::fs::File::create(format!("{}.pairs.tsv", args.out))?;
        writeln!(pf, "cluster_id\ta\tb\tidentity")?;
        let mut n_pairs = 0usize;
        for (i, c) in clusters.iter().enumerate() {
            for j in 0..c.members.len() {
                for k in (j + 1)..c.members.len() {
                    let (Some(&a), Some(&b)) = (node_idx.get(&c.members[j]), node_idx.get(&c.members[k])) else {
                        continue;
                    };
                    let key = if a < b { (a, b) } else { (b, a) };
                    if let Some(&identity) = g.idents.get(&key) {
                        writeln!(
                            pf,
                            "MCL{i}\t{}:{}-{}\t{}:{}-{}\t{identity:.4}",
                            c.members[j].0, c.members[j].1, c.members[j].2,
                            c.members[k].0, c.members[k].1, c.members[k].2
                        )?;
                        n_pairs += 1;
                    }
                }
            }
        }
        eprintln!("[mcl_families] --dump-pairs: {n_pairs} within-cluster edges -> {}.pairs.tsv", args.out);
    }

    // ⭐ Duplicon-first core refinement (§6eh). Post-MCL, per cluster; clusters.tsv above is untouched.
    let mut core_stats = (0usize, 0usize, 0usize, 0usize, 0usize); // gated clusters, kept-full, trimmed, dropped, untouched clusters
    let mut core_records: Vec<Vec<rustle::vg_family::annotation_families::CoreRecord>> = Vec::new();
    let mut sd_pairs: Option<SdPairs> = None; // kept for the §6fw read-through guard
    if args.core_refine {
        let sd = match args.sedef.as_ref() {
            Some(bed) => {
                let text = std::fs::read_to_string(bed).with_context(|| format!("reading {bed}"))?;
                let sd = SdPairs::from_bed_str(&text);
                eprintln!("[mcl_families] core-refine: {} SEDEF pair(s) loaded from {bed}", sd.n_pairs());
                sd
            }
            None => {
                let text = std::fs::read_to_string(&args.paf).with_context(|| format!("reading {}", args.paf))?;
                let sd = SdPairs::from_paf_str(&text);
                eprintln!("[mcl_families] core-refine: {} pair(s) derived from the input PAF (--core-from-paf)", sd.n_pairs());
                sd
            }
        };
        sd_pairs = Some(sd.clone());
        let mut cf = std::fs::File::create(format!("{}.cores.tsv", args.out))?;
        writeln!(cf, "cluster_id\tmember\tgate\tmax_depth\tcore_bp\tspan\tmedian_core\tstatus\tcore_hull")?;
        let mut rf = std::fs::File::create(format!("{}.refined.clusters.tsv", args.out))?;
        writeln!(rf, "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\tstatus")?;
        for (i, c) in clusters.iter().enumerate() {
            let recs = rustle::vg_family::annotation_families::refine_cluster_cores_with(&c.members, &sd, args.core_majority_inclusive);
            core_records.push(recs.clone());
            let gate = recs.first().map_or(false, |r| r.gate_passed);
            if gate {
                core_stats.0 += 1;
            } else {
                core_stats.4 += 1;
            }
            let corr = c.corroborated.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".into());
            let kept: Vec<&rustle::vg_family::annotation_families::CoreRecord> =
                recs.iter().filter(|r| r.status != CoreStatus::Dropped).collect();
            for r in &recs {
                let st = match r.status {
                    CoreStatus::Untouched => "untouched",
                    CoreStatus::KeptFull => "kept_full",
                    CoreStatus::KeptTrimmed => "kept_trimmed",
                    CoreStatus::Dropped => "dropped",
                };
                match r.status {
                    CoreStatus::KeptFull => core_stats.1 += 1,
                    CoreStatus::KeptTrimmed => core_stats.2 += 1,
                    CoreStatus::Dropped => core_stats.3 += 1,
                    CoreStatus::Untouched => {}
                }
                let hull = r.hull.map(|(a, b)| format!("{a}-{b}")).unwrap_or_else(|| "NA".into());
                writeln!(
                    cf,
                    "MCL{i}\t{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{st}\t{hull}",
                    r.member.0, r.member.1, r.member.2, gate, r.max_depth, r.core_bp, r.span, r.median_core
                )?;
                if r.status != CoreStatus::Dropped {
                    let (s0, e0) = match (r.status, r.hull) {
                        (CoreStatus::KeptTrimmed, Some(h)) => h,
                        _ => (r.member.1, r.member.2),
                    };
                    writeln!(
                        rf,
                        "MCL{i}\t{}\t{:.4}\t{:.4}\t{corr}\t{}\t{s0}\t{e0}\t{st}",
                        kept.len(),
                        c.density,
                        c.frac_in,
                        r.member.0
                    )?;
                }
            }
        }
        eprintln!(
            "[mcl_families] core-refine: {} cluster(s) SD-evidenced, {} untouched; members kept-full {}, \
             trimmed {}, dropped {}",
            core_stats.0, core_stats.4, core_stats.1, core_stats.2, core_stats.3
        );
        // ⭐ Duplication blocks (user request 2026-09-05): union-find over every member's core hull, two hulls
        // linked when ONE SEDEF pair overlaps both. Written to `<out>.blocks.tsv`; no existing table changes.
        // A block shared by several clusters is the SEDEF object (LCR16a + LCR16u); the cluster is the family.
        let mut owners: Vec<(usize, GeneKey, (u64, u64))> = Vec::new();
        for (ci, recs) in core_records.iter().enumerate() {
            for r in recs {
                if let Some(h) = r.hull {
                    owners.push((ci, r.member.clone(), h));
                }
            }
        }
        let hulls: Vec<(String, u64, u64)> = owners.iter().map(|(_, m, h)| (m.0.clone(), h.0, h.1)).collect();
        let (blocks, links) = sd_blocks(&hulls, &sd);
        // direct SD partners per cluster: clusters holding a hull joined to one of this cluster's hulls by ONE pair
        let mut direct: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for &(i, j) in &links {
            let (ci, cj) = (owners[i].0, owners[j].0);
            if ci != cj {
                direct.entry(ci).or_default().insert(cj);
                direct.entry(cj).or_default().insert(ci);
            }
        }
        let mut clusters_of_block: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for ((ci, _, _), &b) in owners.iter().zip(blocks.iter()) {
            clusters_of_block.entry(b).or_default().insert(*ci);
        }
        let mut bf = std::fs::File::create(format!("{}.blocks.tsv", args.out))?;
        writeln!(bf, "cluster_id\tmember\tcore_hull\tsd_block\tblock_n_hulls\tblock_clusters\tdirect_sd_partner_clusters")?;
        let n_of_block: BTreeMap<usize, usize> = blocks.iter().fold(BTreeMap::new(), |mut m, &b| {
            *m.entry(b).or_insert(0) += 1;
            m
        });
        for ((ci, m, h), &b) in owners.iter().zip(blocks.iter()) {
            let cl: Vec<String> = clusters_of_block[&b].iter().map(|c| format!("MCL{c}")).collect();
            let dp: Vec<String> = direct.get(ci).map(|s| s.iter().map(|c| format!("MCL{c}")).collect()).unwrap_or_default();
            writeln!(
                bf,
                "MCL{ci}\t{}:{}-{}\t{}-{}\tSDB{b}\t{}\t{}\t{}",
                m.0, m.1, m.2, h.0, h.1, n_of_block[&b], cl.join(","),
                if dp.is_empty() { "-".to_string() } else { dp.join(",") }
            )?;
        }
        let shared: Vec<(usize, &BTreeSet<usize>)> =
            clusters_of_block.iter().filter(|(_, cs)| cs.len() > 1).map(|(b, cs)| (*b, cs)).collect();
        eprintln!(
            "[mcl_families] duplication blocks: {} hull(s) in {} block(s); {} block(s) shared by >1 cluster (e.g. {})",
            owners.len(),
            clusters_of_block.len(),
            shared.len(),
            shared
                .iter()
                .take(3)
                .map(|(b, cs)| format!("SDB{b}={{{}}}", cs.iter().map(|c| format!("MCL{c}")).collect::<Vec<_>>().join(",")))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }

    // ⭐ O1-10b: per-locus units = locus (core hull if refined) + read-supported exon chain (§6el).
    let mut unit_stats = (0usize, 0usize, 0usize, 0usize); // read-chain, gff-fallback, skipped (dropped / no exons), merged into another unit
    struct PendingUnit {
        member: GeneKey,
        exons: Vec<(u64, u64)>,
        strand: char,
        source: String,
        n_reads: usize,
        seq: Vec<u8>,
        hull_col: String,
        sd_depth: String,
        core_bp: String,
        nearest_col: String,
        rep_col: String,
        /// L1: `kept_full` / `kept_trimmed` / `dropped` / `ungated` (the core rule's verdict on the member).
        status: &'static str,
        /// L2: the read-supported locus extent (0-based half-open) — the chain's extent unioned with the
        /// reference span of every BAM record overlapping the locus region. O2's alignment target.
        locus: (u64, u64),
        /// `--coding-core`: longest ORF of `seq`, in bases. 0 when the flag is off.
        orf: usize,
    }
    let mut units_dropped_emitted = 0usize;
    let mut units_unexpressed = 0usize;
    let mut readthrough_units = 0usize;
    let mut noncoding_units = 0usize; // --coding-core: members demoted for not preserving the family's frame
    let mut rt_rejected = (0usize, 0usize, 0usize); // (opposite strand, duplicate flanks, donor not ours)
    // ⭐ `--from-gtf --emit-units`: the copy table of the de novo families (`write_locus_rep_copies`); the read-chain
    // units below are the ANNOTATION mode's (`--paf --gff --bam`), unchanged.
    let mut rep_copy_stats: Option<RepCopyStats> = None;
    if args.emit_units && gtf_loci_list.is_some() {
        anyhow::ensure!(
            !args.emit_readthrough_units,
            "--emit-readthrough-units is not available with --from-gtf: a read-through between two de novo loci is \
             an assembled transcript of its own, not a unit to add"
        );
        anyhow::ensure!(
            !(args.no_cross_family_exon_overlap && !args.allow_cross_family_exon_overlap),
            "--no-cross-family-exon-overlap is not available with --from-gtf (a copy is the locus representative \
             as assembled; it is never trimmed)"
        );
        let fasta = args.fasta.as_ref().ok_or_else(|| anyhow::anyhow!("--from-gtf --emit-units needs --fasta"))?;
        let rmsk = match &args.rmsk {
            Some(path) => Some(load_rmsk(path)?),
            None => None,
        };
        let s = write_locus_rep_copies(
            &args.out,
            fasta,
            &clusters,
            &g,
            &core_records,
            gtf_loci_list.as_deref().unwrap_or(&[]),
            rmsk.as_ref(),
            args.units_include_dropped,
            args.merge_overlapping_units && !args.no_merge_overlapping_units,
            args.coding_core,
        )?;
        unit_stats.2 = s.skipped_dropped;
        unit_stats.3 = s.merged;
        units_dropped_emitted = s.dropped_emitted;
        units_unexpressed = s.unexpressed;
        noncoding_units = s.noncoding;
        eprintln!(
            "[mcl_families] copies (--from-gtf --emit-units): {} cop(ies) = locus representatives in {} famil(ies) \
             ({} with >= 2 copies) -> {}.copies.tsv/.fa/.regions; {} merged into an exon-overlapping copy of the same \
             family, {} dropped member(s) skipped, {} unstranded representative(s) written as +, {} (chrom,start,end) \
             collision(s), {} representative(s) with coalesced exons",
            s.copies, s.families, s.multi_copy_families, args.out, s.merged, s.skipped_dropped, s.unstranded,
            s.key_collisions, s.coalesced
        );
        rep_copy_stats = Some(s);
    } else if args.emit_units {
        let bam = args.bam.as_ref().ok_or_else(|| anyhow::anyhow!("--emit-units needs --bam"))?;
        let fasta = args.fasta.as_ref().ok_or_else(|| anyhow::anyhow!("--emit-units needs --fasta"))?;
        let genome = rustle::genome::GenomeIndex::from_fasta(fasta)?;
        let strands = match &args.gff {
            Some(g) => gene_strands(g)?,
            None => BTreeMap::new(),
        };
        let mut ut = std::fs::File::create(format!("{}.units.tsv", args.out))?;
        let mut um = std::fs::File::create(format!("{}.units.merged.tsv", args.out))?;
        writeln!(um, "cluster_id\tmerged_unit_member\tinto_member")?;
        let mut uf = std::fs::File::create(format!("{}.units.fa", args.out))?;
        let mut ur = std::fs::File::create(format!("{}.units.regions", args.out))?;
        writeln!(ut, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tsource\tcore_hull\tsd_depth\tcore_bp\tnearest_ident\trep_frac\tmember_status\tlocus_start\tlocus_end")?;
        // curated repeats (optional --rmsk): per contig, sorted interspersed intervals
        let rmsk: BTreeMap<String, Vec<(u64, u64)>> = match &args.rmsk {
            Some(path) => load_rmsk(path)?,
            None => BTreeMap::new(),
        };
        let rep_frac = |chrom: &str, exons: &[(u64, u64)]| -> Option<f64> { rep_frac_in(&rmsk, chrom, exons) };
        let node_idx: BTreeMap<&GeneKey, usize> = g.genes.iter().enumerate().map(|(k, gk)| (gk, k)).collect();
        // every family's units are staged first: the L2 clipping (below) needs EVERY unit on a contig, whatever
        // its family — a locus never contains another catalog unit
        let mut staged: Vec<(String, Vec<PendingUnit>, Vec<Option<usize>>)> = Vec::new();
        // (family index, pending index, intron, molecule count) — read-through candidates, resolved after staging
        let mut rt_evidence: Vec<(usize, usize, (u64, u64), usize)> = Vec::new();
        for (i, c) in clusters.iter().enumerate() {
            let fid = format!("MCL{i}");
            let mut pending: Vec<PendingUnit> = Vec::new();
            for (mi, m) in c.members.iter().enumerate() {
                // locus = core hull (trimmed) / member span; dropped members are not units
                // the member's SEDEF core hull (0-based half-open), for `copy_assign --psv-genomic` (§6ep)
                let hull_col = match core_records.get(i).and_then(|v| v.get(mi)).and_then(|r| r.hull) {
                    Some((a, b)) => format!("{}-{}", a.saturating_sub(1), b),
                    None => "NA".to_string(),
                };
                let status: &'static str = match core_records.get(i).and_then(|v| v.get(mi)).map(|r| r.status) {
                    Some(CoreStatus::Dropped) => "dropped",
                    Some(CoreStatus::KeptTrimmed) => "kept_trimmed",
                    Some(_) => "kept_full",
                    None => "ungated",
                };
                let (lo, hi) = match core_records.get(i).and_then(|v| v.get(mi)) {
                    Some(r) if r.status == CoreStatus::Dropped => {
                        if !args.units_include_dropped {
                            unit_stats.2 += 1;
                            continue;
                        }
                        (m.1.saturating_sub(1), m.2) // L1: a dropped member's locus is its annotated span
                    }
                    Some(r) if r.status == CoreStatus::KeptTrimmed => match r.hull {
                        Some((a, b)) => (a.saturating_sub(1), b),
                        None => (m.1.saturating_sub(1), m.2),
                    },
                    _ => (m.1.saturating_sub(1), m.2),
                };
                let (_, reads) = rustle::vg_family::denovo_assemble::reads_in_region(bam, &m.0, lo, hi, 1)
                    .with_context(|| format!("reading {}:{}-{}", m.0, lo, hi))?;
                let (chain, rstrand, n_reads) = read_chain(&reads, lo, hi, args.min_reads, args.units_follow_reads.then(|| (m.1.saturating_sub(1), m.2)));
                let gff_strand = strands.get(m).copied().unwrap_or('+');
                let (exons, strand, source): (Vec<(u64, u64)>, char, &str) = if !chain.is_empty() {
                    (chain, rstrand.unwrap_or(gff_strand), "read_chain")
                } else {
                    let ex: Vec<(u64, u64)> = blocks
                        .get(m)
                        .map(|v| v.iter().filter(|&&(s, e)| e > lo && s < hi).map(|&(s, e)| (s.max(lo), e.min(hi))).collect())
                        .unwrap_or_default();
                    if ex.is_empty() {
                        unit_stats.2 += 1;
                        continue;
                    }
                    (ex, gff_strand, "gff_fallback")
                };
                // reads with an aligned block inside the EMITTED chain (the copies.tsv contract: every copy has
                // ≥1 read inside it; counting reads in the hull aborted 4 sweep families on chain-less units, §6es)
                let n_in_chain = reads
                    .iter()
                    .filter(|br| !br.is_supplementary && !br.is_secondary)
                    .filter(|br| {
                        let (bl, _) = blocks_and_introns(br);
                        bl.iter().any(|&(bs, be)| exons.iter().any(|&(s, e)| be > s && bs < e))
                    })
                    .count();
                let keep_unexpressed = args.units_keep_unexpressed && !args.no_units_keep_unexpressed;
                if n_in_chain == 0 && !(keep_unexpressed && source == "gff_fallback") {
                    unit_stats.2 += 1;
                    continue;
                }
                if n_in_chain == 0 {
                    units_unexpressed += 1;
                }
                // read-through evidence: introns of primaries with a block in this chain that leave it
                // downstream. Molecules are counted by name so a split record cannot vote twice.
                if args.emit_readthrough_units {
                    let mut seen: BTreeMap<(u64, u64), std::collections::BTreeSet<String>> = BTreeMap::new();
                    for br in reads.iter().filter(|br| !br.is_supplementary && !br.is_secondary) {
                        let (bl, itr) = blocks_and_introns(br);
                        for j in leaving_introns(&bl, &itr, &exons) {
                            seen.entry(j).or_default().insert(br.name.clone());
                        }
                    }
                    for (j, names) in seen {
                        rt_evidence.push((staged.len(), pending.len(), j, names.len()));
                    }
                }
                let n_reads = n_in_chain;
                let mut seq: Vec<u8> = Vec::new();
                for &(s, e) in &exons {
                    let Some(part) = genome.fetch_sequence(&m.0, s, e) else {
                        anyhow::bail!("--emit-units: {}:{s}-{e} is not in --fasta", m.0)
                    };
                    seq.extend_from_slice(&part);
                }
                if strand == '-' {
                    seq = rustle::vg_family::seq_utils::revcomp_keep_case(&seq);
                }
                let (sd_depth, core_bp) = core_records
                    .get(i)
                    .and_then(|v| v.get(mi))
                    .map(|r| (r.max_depth.to_string(), r.core_bp.to_string()))
                    .unwrap_or_else(|| ("NA".into(), "NA".into()));
                let nearest = node_idx.get(m).map(|&a| {
                    c.members
                        .iter()
                        .filter_map(|o| node_idx.get(o).copied())
                        .filter(|&b| b != a)
                        .filter_map(|b| g.idents.get(&(a.min(b), a.max(b))).copied())
                        .fold(0.0f64, f64::max)
                });
                let nearest_col = nearest.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".into());
                let rep_col = rep_frac(&m.0, &exons).map(|v| format!("{v:.3}")).unwrap_or_else(|| "NA".into());
                let locus = {
                    let (us, ue) = (exons[0].0, exons.last().unwrap().1);
                    match read_extent(&reads, &exons, args.min_reads) {
                        Some((a, b)) => (a.min(us), b.max(ue)),
                        None => (us, ue),
                    }
                };
                if status == "dropped" {
                    units_dropped_emitted += 1;
                }
                pending.push(PendingUnit {
                    orf: if args.coding_core { longest_orf(&seq) } else { 0 },
                    member: m.clone(),
                    exons,
                    strand,
                    source: source.to_string(),
                    n_reads,
                    seq,
                    hull_col,
                    sd_depth,
                    core_bp,
                    nearest_col,
                    rep_col,
                    status,
                    locus,
                });
            }
            // ⭐ CODING CORE: a member keeps its status only if its longest ORF reaches half the family's
            // best — the family, not an absolute aa cutoff, is the reference. Failing members become
            // `noncoding` candidates. Dropped members are left alone: the core rule already spoke.
            if args.coding_core {
                let best_orf = pending.iter().filter(|u| u.status != "dropped").map(|u| u.orf).max().unwrap_or(0);
                if best_orf > 0 {
                    for u in pending.iter_mut() {
                        if u.status != "dropped" && u.orf * 2 < best_orf {
                            u.status = "noncoding";
                            noncoding_units += 1;
                        }
                    }
                }
            }
            // ⭐ Units of one family that share EXON bases are one locus (§6fb): the read-followed chain of a large
            // record can cover records nested in it (MCL108: a 1.16-Mb unit with a 143-bp and a 2.7-kb unit inside
            // its exons, 13,000 reads counted three times, every one of them a K = 0 tie under read-star). A base
            // cannot belong to two copies. Representative = the longest exon union; the others are recorded in
            // `<out>.units.merged.tsv`. `--no-merge-overlapping-units` keeps every unit (byte-identical to before).
            let merged_into: Vec<Option<usize>> = if args.merge_overlapping_units && !args.no_merge_overlapping_units {
                let n = pending.len();
                let mut parent: Vec<usize> = (0..n).collect();
                fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
                    while p[x] != x {
                        p[x] = p[p[x]];
                        x = p[x];
                    }
                    x
                }
                for a in 0..n {
                    for b in (a + 1)..n {
                        if pending[a].member.0 != pending[b].member.0 {
                            continue;
                        }
                        let share = pending[a].exons.iter().any(|&(s1, e1)| pending[b].exons.iter().any(|&(s2, e2)| s1 < e2 && s2 < e1));
                        if share {
                            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                            if ra != rb {
                                parent[ra.max(rb)] = ra.min(rb);
                            }
                        }
                    }
                }
                // representative = kept before dropped (L1), then the longest exon union
                let rank = |u: &PendingUnit| (u.status != "dropped", u.exons.iter().map(|(s, e)| e - s).sum::<u64>());
                let mut rep_of_root: BTreeMap<usize, usize> = BTreeMap::new();
                for k in 0..n {
                    let r = find(&mut parent, k);
                    let e = rep_of_root.entry(r).or_insert(k);
                    if rank(&pending[k]) > rank(&pending[*e]) {
                        *e = k;
                    }
                }
                (0..n)
                    .map(|k| {
                        let r = find(&mut parent, k);
                        let rep = rep_of_root[&r];
                        if rep == k { None } else { Some(rep) }
                    })
                    .collect()
            } else {
                vec![None; pending.len()]
            };
            staged.push((fid, pending, merged_into));
        }
        // ⭐ §6ge: cross-family exon overlap. Two units of different families claiming the same bases turn an
        // ordinary intron of one into a "read-through" of the other. Contested bases are given to the unit
        // whose own annotated member span contains them, else to the unit with more reads; the loser's exons
        // are trimmed. Within one family the §6fb merge already handles overlap, so same-family pairs are skipped.
        let mut cross_trimmed = 0usize;
        if args.no_cross_family_exon_overlap && !args.allow_cross_family_exon_overlap {
            let mut by_ctg: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
            for (fi, (_, pending, _)) in staged.iter().enumerate() {
                for (k, u) in pending.iter().enumerate() {
                    by_ctg.entry(u.member.0.clone()).or_default().push((fi, k));
                }
            }
            for ks in by_ctg.values() {
                for a in 0..ks.len() {
                    for b in (a + 1)..ks.len() {
                        let ((fa, ka), (fb, kb)) = (ks[a], ks[b]);
                        if fa == fb {
                            continue; // §6fb owns same-family overlap
                        }
                        let (ea, eb) = (staged[fa].1[ka].exons.clone(), staged[fb].1[kb].exons.clone());
                        let shared: Vec<(u64, u64)> = ea
                            .iter()
                            .flat_map(|&(s1, e1)| eb.iter().filter_map(move |&(s2, e2)| {
                                let (lo, hi) = (s1.max(s2), e1.min(e2));
                                (hi > lo).then_some((lo, hi))
                            }))
                            .collect();
                        if shared.is_empty() {
                            continue;
                        }
                        // owner: whose ANNOTATED member span holds more of the contested bases; ties -> reads
                        let inside = |u: &PendingUnit| -> u64 {
                            let (ms, me) = (u.member.1.saturating_sub(1), u.member.2);
                            shared.iter().map(|&(s, e)| e.min(me).saturating_sub(s.max(ms))).sum()
                        };
                        let (ia, ib) = (inside(&staged[fa].1[ka]), inside(&staged[fb].1[kb]));
                        let loser = if ia != ib {
                            if ia < ib { (fa, ka) } else { (fb, kb) }
                        } else if staged[fa].1[ka].n_reads <= staged[fb].1[kb].n_reads {
                            (fa, ka)
                        } else {
                            (fb, kb)
                        };
                        let u = &mut staged[loser.0].1[loser.1];
                        let mut out: Vec<(u64, u64)> = Vec::new();
                        for &(s, e) in &u.exons {
                            let mut segs = vec![(s, e)];
                            for &(cs, ce) in &shared {
                                segs = segs
                                    .into_iter()
                                    .flat_map(|(a0, b0)| {
                                        let mut v = Vec::new();
                                        if a0 < cs.min(b0) { v.push((a0, cs.min(b0))); }
                                        if ce.max(a0) < b0 { v.push((ce.max(a0), b0)); }
                                        if ce <= a0 || cs >= b0 { v.clear(); v.push((a0, b0)); }
                                        v
                                    })
                                    .filter(|&(a0, b0)| b0 > a0)
                                    .collect();
                            }
                            out.extend(segs);
                        }
                        out.sort_unstable();
                        if out != u.exons {
                            cross_trimmed += 1;
                        }
                        u.exons = out;
                    }
                }
            }
            // a unit trimmed to nothing is no longer a unit
            for (_, pending, _) in staged.iter_mut() {
                pending.retain(|u| !u.exons.is_empty());
            }
            eprintln!("[mcl_families] no-cross-family-exon-overlap: {cross_trimmed} unit(s) trimmed");
        }
        // ⭐ L2: clip every emitted unit's extent at the chain ends of its neighbours on the contig — units of
        // EVERY family (§6fm: MCL1971's 951-kb extent contained MCL42's unit and other genes; reads of those
        // loci aligned perfectly inside the target and were accepted as MCL1971's sole candidates at 3 %
        // divergence from its own unit). Nested units (overlapping spans) do not clip.
        let clipped: Vec<Vec<(u64, u64)>> = {
            let mut by_ctg: BTreeMap<&str, Vec<(usize, usize)>> = BTreeMap::new();
            for (fi, (_, pending, merged_into)) in staged.iter().enumerate() {
                for (k, u) in pending.iter().enumerate() {
                    if merged_into[k].is_none() {
                        by_ctg.entry(u.member.0.as_str()).or_default().push((fi, k));
                    }
                }
            }
            let mut out: Vec<Vec<(u64, u64)>> = staged.iter().map(|(_, p, _)| p.iter().map(|u| u.locus).collect()).collect();
            for ks in by_ctg.values() {
                let spans: Vec<(u64, u64, (u64, u64))> = ks
                    .iter()
                    .map(|&(fi, k)| {
                        let u = &staged[fi].1[k];
                        (u.exons[0].0, u.exons.last().unwrap().1, u.locus)
                    })
                    .collect();
                for (&(fi, k), c) in ks.iter().zip(clip_extents_to_neighbours(&spans)) {
                    out[fi][k] = c;
                }
            }
            out
        };
        // ⭐ The conjoined read-through as its own unit (PREREG_readthrough_object, md5 1a51fa3b). Resolved
        // only now: the far end of a junction must land in an EMITTED chain, which is known after staging.
        // Rows are (family index, source pending index, target family, target pending, intron, molecules).
        let readthrough_rows: Vec<(usize, usize, usize, usize, (u64, u64), usize)> = if args.emit_readthrough_units {
            let mut flat: Vec<(String, Vec<(u64, u64)>)> = Vec::new();
            let mut where_of: Vec<(usize, usize)> = Vec::new();
            for (fi, (_, pending, merged_into)) in staged.iter().enumerate() {
                for (k, u) in pending.iter().enumerate() {
                    if merged_into[k].is_none() {
                        flat.push((u.member.0.clone(), u.exons.clone()));
                        where_of.push((fi, k));
                    }
                }
            }
            let pos_of: BTreeMap<(usize, usize), usize> = where_of.iter().enumerate().map(|(a, &b)| (b, a)).collect();
            let mut best: BTreeMap<(usize, usize, usize, usize), ((u64, u64), usize)> = BTreeMap::new();
            for &(fi, k, j, n) in &rt_evidence {
                if n < args.min_reads {
                    continue;
                }
                let Some(&src) = pos_of.get(&(fi, k)) else { continue };
                let u = &staged[fi].1[k];
                // ⭐ §6ge: re-check the DONOR against the FINAL chain. The evidence above was gathered during
                // staging, before the cross-family trim, so a junction whose donor was trimmed away belongs to
                // the neighbouring gene and is its ordinary intron, not a read-through out of this unit.
                if !u.exons.iter().any(|&(_, ce)| ce == j.0) {
                    rt_rejected.2 += 1;
                    continue;
                }
                let Some(tgt) = readthrough_target(j, &u.member.0, src, &flat) else { continue };
                if !canonical_intron(&genome, &u.member.0, j.0, j.1, u.strand) {
                    continue;
                }
                let (tfi, tk) = where_of[tgt];
                // ⭐ §6fw guard (default ON): the two units must share a strand — a transcript cannot join
                // opposite strands — and their donor/acceptor flanks must NOT be duplicates of one another,
                // since that is exactly what a cross-copy mis-chain needs in order to jump. PAD is the 2 kb
                // window the adjudication used. Without SD pairs only the strand half of the guard can run.
                if !args.no_readthrough_guard {
                    const PAD: u64 = 2_000;
                    if staged[tfi].1[tk].strand != u.strand {
                        rt_rejected.0 += 1;
                        continue;
                    }
                    if let Some(sd) = sd_pairs.as_ref() {
                        let donor = (j.0.saturating_sub(PAD), j.0 + PAD);
                        let acceptor = (j.1.saturating_sub(PAD), j.1 + PAD);
                        if sd.links(&u.member.0, donor, acceptor) {
                            rt_rejected.1 += 1;
                            continue;
                        }
                    }
                }
                let e = best.entry((fi, k, tfi, tk)).or_insert((j, 0));
                if n > e.1 {
                    *e = (j, n);
                }
            }
            best.into_iter().map(|((a, b, c, d), (j, n))| (a, b, c, d, j, n)).collect()
        } else {
            Vec::new()
        };
        readthrough_units = readthrough_rows.len();
        // side file: which two units each read-through joins, and through which intron (the 19-column
        // units.tsv contract is untouched, so nothing downstream has to change to read this)
        if args.emit_readthrough_units {
            let mut rf = std::fs::File::create(format!("{}.readthrough.tsv", args.out))?;
            writeln!(rf, "from_family\tchrom\tfrom_unit\tto_family\tto_unit\tintron_start\tintron_end\tintron_bp\tn_molecules")?;
            for &(sfi, sk, tfi, tk, j, n) in &readthrough_rows {
                let (a, b) = (&staged[sfi].1[sk], &staged[tfi].1[tk]);
                writeln!(
                    rf, "{}\t{}\t{}:{}-{}\t{}\t{}:{}-{}\t{}\t{}\t{}\t{n}",
                    staged[sfi].0, a.member.0, a.member.0, a.exons[0].0, a.exons.last().unwrap().1,
                    staged[tfi].0, b.member.0, b.exons[0].0, b.exons.last().unwrap().1, j.0, j.1, j.1 - j.0
                )?;
            }
        }
        for (fi, (fid, pending, merged_into)) in staged.iter().enumerate() {
            let mut idx = 0usize;
            let mut hulls: BTreeMap<String, (u64, u64)> = BTreeMap::new();
            for (k, u) in pending.iter().enumerate() {
                if let Some(rep) = merged_into[k] {
                    let r = &pending[rep];
                    writeln!(um, "{fid}\t{}:{}-{}\t{}:{}-{}", u.member.0, u.member.1, u.member.2, r.member.0, r.member.1, r.member.2)?;
                    unit_stats.3 += 1;
                    continue;
                }
                if u.source == "read_chain" {
                    unit_stats.0 += 1;
                } else {
                    unit_stats.1 += 1;
                }
                let (us, ue) = (u.exons[0].0, u.exons.last().unwrap().1);
                let (m, exons, strand, source) = (&u.member, &u.exons, u.strand, u.source.as_str());
                writeln!(
                    ut,
                    "{fid}\t{idx}\tMCL_{}_{us}\t{}\t{us}\t{ue}\t{}\t{strand}\t{}\t{}\t{source}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    m.0,
                    m.0,
                    exons.len(),
                    u.n_reads,
                    exons.iter().map(|(s, e)| format!("{s}-{e}")).collect::<Vec<_>>().join(","),
                    u.hull_col, u.sd_depth, u.core_bp, u.nearest_col, u.rep_col, u.status, clipped[fi][k].0, clipped[fi][k].1
                )?;
                writeln!(uf, ">{fid}|{idx}|{}:{us}-{ue}|{strand}|nexon={}", m.0, exons.len())?;
                uf.write_all(&u.seq)?;
                writeln!(uf)?;
                let h = hulls.entry(m.0.clone()).or_insert((us, ue));
                h.0 = h.0.min(us);
                h.1 = h.1.max(ue);
                idx += 1;
            }
            for &(sfi, sk, tfi, tk, _j, n) in readthrough_rows.iter().filter(|r| r.0 == fi) {
                let (a, b) = (&staged[sfi].1[sk], &staged[tfi].1[tk]);
                // the two chains may share bases (a unit's read-followed chain can reach into its neighbour),
                // and the copies.tsv contract requires ascending disjoint blocks: coalesce after sorting
                let mut exons: Vec<(u64, u64)> = a.exons.iter().chain(b.exons.iter()).copied().collect();
                exons.sort_unstable();
                exons = exons.into_iter().fold(Vec::new(), |mut acc: Vec<(u64, u64)>, (s, e)| {
                    match acc.last_mut() {
                        Some(l) if s <= l.1 => l.1 = l.1.max(e),
                        _ => acc.push((s, e)),
                    }
                    acc
                });
                let (us, ue) = (exons[0].0, exons.last().unwrap().1);
                let ctg = &a.member.0;
                let mut seq: Vec<u8> = Vec::new();
                for &(s, e) in &exons {
                    let Some(part) = genome.fetch_sequence(ctg, s, e) else {
                        anyhow::bail!("--emit-readthrough-units: {ctg}:{s}-{e} is not in --fasta")
                    };
                    seq.extend_from_slice(&part);
                }
                if a.strand == '-' {
                    seq = rustle::vg_family::seq_utils::revcomp_keep_case(&seq);
                }
                writeln!(
                    ut,
                    "{fid}\t{idx}\tMCL_{ctg}_{us}\t{ctg}\t{us}\t{ue}\t{}\t{}\t{n}\t{}\treadthrough\tNA\tNA\tNA\tNA\tNA\treadthrough\t{us}\t{ue}",
                    exons.len(),
                    a.strand,
                    exons.iter().map(|(s, e)| format!("{s}-{e}")).collect::<Vec<_>>().join(","),
                )?;
                writeln!(uf, ">{fid}|{idx}|{ctg}:{us}-{ue}|{}|nexon={}|readthrough", a.strand, exons.len())?;
                uf.write_all(&seq)?;
                writeln!(uf)?;
                let h = hulls.entry(ctg.clone()).or_insert((us, ue));
                h.0 = h.0.min(us);
                h.1 = h.1.max(ue);
                idx += 1;
            }
            for (ctg, (a, b)) in hulls {
                writeln!(ur, "{fid}\t{ctg}:{}-{}", a.saturating_sub(5_000).max(1), b + 5_000)?;
            }
        }
        eprintln!(
            "[mcl_families] emit-units: {} read-chain unit(s), {} GFF-fallback unit(s), {} member(s) without a unit \
             (dropped, no exon inside the locus, or no read inside the chain), {} unit(s) merged into an overlapping unit of the same family",
            unit_stats.0, unit_stats.1, unit_stats.2, unit_stats.3
        );
        if args.coding_core {
            eprintln!("[mcl_families] coding-core: {noncoding_units} member(s) demoted to `noncoding` (longest ORF below half the family's best)");
        }
        if args.emit_readthrough_units {
            eprintln!(
                "[mcl_families] emit-readthrough-units: {readthrough_units} conjoined read-through unit(s); \
                 guard rejected {} on opposite strands, {} whose flanks are linked by a duplication pair, \
                 {} whose donor is not this unit's own exon end",
                rt_rejected.0, rt_rejected.1, rt_rejected.2
            );
        }
    }

    // ⭐ `--emit-container` (with `--from-gtf`): the container of every member's extra pieces, read from the products
    // written above exactly as the frozen post-processor reads them (clusters.tsv, loci.gff3, the fold table this run
    // wrote, the PAF, the GTF); it changes nothing already written.
    let mut container_counts: Option<[(&str, i64); 3]> = None;
    if args.emit_container {
        use rustle::vg_family::family_container as fc;
        let gtf = args.from_gtf.as_deref().expect("--emit-container is checked to come with --from-gtf");
        let gff3 = args.gff.as_deref().expect("--from-gtf sets --gff to <out>.loci.gff3");
        let clusters_path = format!("{}.clusters.tsv", args.out);
        // the fold table exists iff this run wrote it (never a stale one from an earlier run)
        let loci_tsv = loci.as_ref().map(|_| format!("{}.loci.tsv", args.out));
        let mut lt = match &loci_tsv {
            Some(p) => Some(fc::open_text(p)?),
            None => None,
        };
        let c = fc::run(
            &mut *fc::open_text(gtf)?,
            &mut *fc::open_text(&clusters_path)?,
            &mut *fc::open_text(gff3)?,
            lt.as_deref_mut().map(|r| r as &mut dyn std::io::BufRead),
            &mut paf.as_bytes(),
            &args.paf,
        )
        .context("--emit-container")?;
        c.write(&args.out)?;
        eprintln!("{}", c.summary_line());
        eprintln!(
            "[mcl_families] container: {} block(s), {} accessory, {} directed family relation(s) -> {}.container.tsv, \
             {}.container_relations.tsv, {}.container_summary.tsv",
            c.count("blocks"),
            c.count("accessory_blocks"),
            c.count("family_relations_directed"),
            args.out,
            args.out,
            args.out
        );
        container_counts = Some([
            ("container_blocks", c.count("blocks")),
            ("container_accessory_blocks", c.count("accessory_blocks")),
            ("container_family_relations", c.count("family_relations_directed")),
        ]);
    }

    // Params certificate: a flag with no certificate row makes two arms indistinguishable.
    let mut ph = std::fs::File::create(format!("{}.params.tsv", args.out))?;
    for (k, v) in [
        ("paf".to_string(), args.paf.clone()),
        ("gff".to_string(), args.gff.clone().unwrap_or_else(|| "<unset>".into())),
        ("exonic_denominator".to_string(), (!exonic.is_empty()).to_string()),
        ("exonic_lengths_loaded".to_string(), exonic.len().to_string()),
        ("nodes_fell_back_to_span".to_string(), g.missing_exonic.to_string()),
        ("exonic_overlap".to_string(), args.exonic_overlap.to_string()),
        ("exonic_overlap_joined".to_string(), g.exonic_overlap_joined.to_string()),
        ("exonic_overlap_missing".to_string(), g.exonic_overlap_missing.to_string()),
        ("reject_overlapping".to_string(), args.reject_overlapping.to_string()),
        ("rejected_overlapping".to_string(), g.rejected_overlapping.to_string()),
        ("min_exonic_bp".to_string(), args.min_exonic_bp.to_string()),
        ("rejected_no_exonic".to_string(), g.rejected_no_exonic.to_string()),
        ("min_shared_exon_frac".to_string(), args.min_shared_exon_frac.to_string()),
        ("rejected_low_shared_exon".to_string(), g.rejected_low_shared_exon.to_string()),
        ("merge_overlapping_loci".to_string(), args.merge_overlapping_loci.to_string()),
        ("locus_attribute_edges".to_string(), args.locus_attribute_edges.to_string()),
        ("annotations_folded_into_loci".to_string(), loci.as_ref().map_or(0, |m| m.n_merged()).to_string()),
        ("paf_records_same_locus_skipped".to_string(), g.same_locus_records.to_string()),
        ("core_refine".to_string(), args.core_refine.to_string()),
        ("sedef".to_string(), args.sedef.clone().unwrap_or_else(|| "<unset>".into())),
        ("core_from_paf".to_string(), (args.core_from_paf && args.sedef.is_none()).to_string()),
        ("core_clusters_gated".to_string(), core_stats.0.to_string()),
        ("core_members_kept_full".to_string(), core_stats.1.to_string()),
        ("core_members_trimmed".to_string(), core_stats.2.to_string()),
        ("core_members_dropped".to_string(), core_stats.3.to_string()),
        ("emit_units".to_string(), args.emit_units.to_string()),
        ("units_follow_reads".to_string(), args.units_follow_reads.to_string()),
        ("fold_within_clusters".to_string(), args.fold_within_clusters.to_string()),
        ("exonless_span".to_string(), args.exonless_span.to_string()),
        ("exonic_both_sides".to_string(), args.exonic_both_sides.to_string()),
        ("merge_overlapping_units".to_string(), (args.merge_overlapping_units && !args.no_merge_overlapping_units).to_string()),
        ("units_merged".to_string(), unit_stats.3.to_string()),
        ("units_include_dropped".to_string(), args.units_include_dropped.to_string()),
        ("emit_readthrough_units".to_string(), args.emit_readthrough_units.to_string()),
        ("readthrough_units".to_string(), readthrough_units.to_string()),
        ("readthrough_guard".to_string(), (!args.no_readthrough_guard).to_string()),
        ("coding_core".to_string(), args.coding_core.to_string()),
        ("no_cross_family_exon_overlap".to_string(), (args.no_cross_family_exon_overlap && !args.allow_cross_family_exon_overlap).to_string()),
        ("noncoding_units".to_string(), noncoding_units.to_string()),
        ("readthrough_rejected_strand".to_string(), rt_rejected.0.to_string()),
        ("readthrough_rejected_duplicate_flanks".to_string(), rt_rejected.1.to_string()),
        ("readthrough_rejected_foreign_donor".to_string(), rt_rejected.2.to_string()),
        ("units_dropped_emitted".to_string(), units_dropped_emitted.to_string()),
        ("units_keep_unexpressed".to_string(), (args.units_keep_unexpressed && !args.no_units_keep_unexpressed).to_string()),
        ("units_unexpressed".to_string(), units_unexpressed.to_string()),
        ("core_majority_inclusive".to_string(), args.core_majority_inclusive.to_string()),
        ("rmsk".to_string(), args.rmsk.clone().unwrap_or_else(|| "NA".into())),
        ("units_read_chain".to_string(), unit_stats.0.to_string()),
        ("units_gff_fallback".to_string(), unit_stats.1.to_string()),
        ("inflation".to_string(), args.inflation.to_string()),
        ("prune".to_string(), format!("{:e}", args.prune)),
        ("min_identity".to_string(), p.min_identity.to_string()),
        ("min_cov_longer".to_string(), p.min_cov_longer.to_string()),
        ("min_bp".to_string(), p.min_bp.to_string()),
        ("min_size".to_string(), args.min_size.to_string()),
        ("bam".to_string(), args.bam.clone().unwrap_or_else(|| "<unset>".into())),
        ("corroboration_min_reads".to_string(), args.min_reads.to_string()),
        ("dump_pairs".to_string(), args.dump_pairs.to_string()),
        ("n_nodes".to_string(), g.n_nodes().to_string()),
        ("n_edges".to_string(), g.n_edges().to_string()),
        ("n_clusters".to_string(), clusters.len().to_string()),
        // §6x4, appended LAST so every existing positional reader of params.tsv is unaffected (r936).
        ("min_cov_shorter".to_string(), p.min_cov_shorter.to_string()),
        ("admitted_by_containment".to_string(), g.admitted_by_containment.to_string()),
    ] {
        writeln!(ph, "{k}\t{v}")?;
    }
    // `--from-gtf --emit-units` only (appended LAST, and only then, so every other run's params.tsv is unchanged)
    if let Some(s) = &rep_copy_stats {
        for (k, v) in [
            ("copies_from_gtf", "true".to_string()),
            ("copies_written", s.copies.to_string()),
            ("copies_families", s.families.to_string()),
            ("copies_multi_copy_families", s.multi_copy_families.to_string()),
            ("copies_unstranded_as_plus", s.unstranded.to_string()),
            ("copies_gene_key_collisions", s.key_collisions.to_string()),
            ("copies_exons_coalesced", s.coalesced.to_string()),
        ] {
            writeln!(ph, "{k}\t{v}")?;
        }
    }
    // `--emit-container` only (appended LAST, and only then, so every other run's params.tsv is unchanged)
    if let Some(cc) = container_counts {
        writeln!(ph, "emit_container\ttrue")?;
        for (k, v) in cc {
            writeln!(ph, "{k}\t{v}")?;
        }
    }
    // `--representative most-junctions` only (appended LAST, and only then: a default run's params.tsv is unchanged)
    if args.representative == Representative::MostJunctions {
        writeln!(ph, "representative\tmost-junctions")?;
    }

    let members: usize = clusters.iter().map(|c| c.members.len()).sum();
    let largest = clusters.iter().map(|c| c.members.len()).max().unwrap_or(0);
    let zero = clusters.iter().filter(|c| c.corroborated == Some(0.0)).count();
    eprintln!(
        "[mcl_families] {} cluster(s) >= {} members, {members} members, largest {largest}",
        clusters.len(),
        args.min_size
    );
    if corr_set.is_some() {
        eprintln!(
            "[mcl_families] ZERO-corroboration clusters (the repeat-clique signature): {zero}/{}",
            clusters.len()
        );
    }
    eprintln!("[mcl_families] wrote {}.clusters.tsv + {}.params.tsv", args.out, args.out);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustle::vg_family::copy_split::AlignedRead;
    use rustle::vg_family::denovo_assemble::BamRead;

    /// `--min-cov-shorter` defaults to 0.70 (2026-09-29); `0` is the explicit OFF that `graph_from_paf` treats as
    /// "no escape" (`p.min_cov_shorter > 0.0`), the edge weights of every catalog built before the flip.
    #[test]
    fn min_cov_shorter_defaults_to_0_70_and_zero_turns_it_off() {
        let a = Args::try_parse_from(["mcl_families", "--out", "o"]).expect("parse");
        assert_eq!(a.min_cov_shorter, 0.70);
        let a = Args::try_parse_from(["mcl_families", "--out", "o", "--min-cov-shorter", "0"]).expect("parse");
        assert_eq!(a.min_cov_shorter, 0.0);
        assert_eq!(GraphParams::default().min_cov_shorter, 0.0, "the library default stays off: only the CLI flips");
    }

    fn br(start: u64, cigar: &[(char, u64)]) -> BamRead {
        let n: u64 = cigar.iter().filter(|(o, _)| matches!(o, 'M' | '=' | 'X' | 'I' | 'S')).map(|(_, l)| l).sum();
        BamRead {
            chrom: "c".into(),
            read: AlignedRead { ref_start: start, cigar: cigar.to_vec(), seq: vec![b'A'; n as usize], qual: vec![] },
            mapq: 60,
            name: String::new(),
            as_score: 0,
            de: 0.0,
            is_supplementary: false,
            is_secondary: false,
            reverse: false,
            ts: Some('+'),
        }
    }

    /// The window clips the chain unless the units follow the reads: three reads share an exon at 100-200
    /// that lies OUTSIDE the locus window 900-1200 and splice into 1000-1100 inside it.
    #[test]
    fn read_chain_follows_reads_beyond_the_window_only_when_asked() {
        let reads: Vec<BamRead> = (0..3).map(|_| br(100, &[('M', 100), ('N', 800), ('M', 100)])).collect();
        let (clipped, _, n) = read_chain(&reads, 900, 1200, 3, None);
        assert_eq!((clipped, n), (vec![(1000, 1100)], 3));
        let (followed, strand, _) = read_chain(&reads, 900, 1200, 3, Some((0, 2000)));
        assert_eq!(followed, vec![(100, 200), (1000, 1100)]);
        assert_eq!(strand, Some('+'));
        // the annotated locus span bounds the following: a block outside it stays clipped (no engulfment)
        let (bounded, _, _) = read_chain(&reads, 900, 1200, 3, Some((150, 2000)));
        assert_eq!(bounded, vec![(150, 200), (1000, 1100)]);
        // a block supported by fewer than min_reads reads is not followed
        let mut reads2 = reads;
        reads2.push(br(5000, &[('M', 50), ('N', 200), ('M', 50)]));
        reads2[3].read.ref_start = 1150; // 1150-1200 (in window) N 1400-1450 (outside, 1 read)
        let (followed2, _, _) = read_chain(&reads2, 900, 1200, 3, Some((0, 2000)));
        assert_eq!(followed2, vec![(100, 200), (1000, 1100)]);
    }

    /// L2: a locus is clipped at the chain ends of its family neighbours, never inside its own span, and a
    /// nested or interleaved unit does not clip.
    #[test]
    fn clip_extents_stop_at_the_neighbouring_units_of_the_family() {
        let units = [(1000, 2000, (100, 5000)), (3000, 4000, (1500, 9000)), (1200, 1300, (900, 1400)), (8000, 8500, (7000, 9500))];
        // unit 0 stops at unit 1's chain (3000); unit 1 starts at unit 0's chain end (2000) and stops at unit 3's
        // start (8000); the nested unit 2 is not clipped by its host (they overlap) and stops at 3000 → 1400 is
        // its own bound; unit 3's own extent (7000) is inside every bound
        assert_eq!(clip_extents_to_neighbours(&units), vec![(100, 3000), (2000, 8000), (900, 1400), (7000, 9500)]);
        assert_eq!(clip_extents_to_neighbours(&[(10, 20, (5, 50))]), vec![(5, 50)]);
    }

    /// The ORF scan is frame-aware, takes the LONGEST ATG..stop across the three forward frames, is
    /// case-insensitive, and reports 0 when no complete ORF exists. A frameshift shortens it, which is the
    /// whole point of `--coding-core` (the rule compares this length to the family's best, never to a
    /// fixed number of amino acids).
    #[test]
    fn longest_orf_takes_the_best_complete_reading_frame() {
        assert_eq!(longest_orf(b""), 0);
        assert_eq!(longest_orf(b"ATGAAACCC"), 0, "no stop codon: not a complete ORF");
        assert_eq!(longest_orf(b"ATGAAATAA"), 9);
        assert_eq!(longest_orf(b"atgaaataa"), 9, "case-insensitive");
        // frame 1: a leading base pushes the same ORF into another frame
        assert_eq!(longest_orf(b"CATGAAATAA"), 9);
        // a single-base insertion after the start breaks the frame, so the ORF collapses
        assert_eq!(longest_orf(b"ATGAAACCCTAA"), 12);
        assert_eq!(longest_orf(b"ATGAAAGCCCTAA"), 0, "frameshift removes the in-frame stop");
        // the longer of two ORFs wins
        assert_eq!(longest_orf(b"ATGTAAGGGATGAAACCCTAA"), 12);
    }

    /// A read-through junction is an intron of a read that HAS a block in the chain and whose far end
    /// leaves the chain downstream. A read spliced over the chain with no block in it contributes nothing,
    /// and an ordinary intron inside the chain is not a candidate.
    #[test]
    fn leaving_introns_are_downstream_junctions_of_reads_anchored_in_the_chain() {
        let chain = [(1000u64, 1100u64), (4000, 4120)];
        // no block in the chain: nothing, however far the intron reaches
        let r = br(100, &[('M', 50), ('N', 90_000), ('M', 50)]);
        let (b, i) = blocks_and_introns(&r);
        assert!(leaving_introns(&b, &i, &chain).is_empty());
        // anchored, intron entirely inside the chain span: not a read-through
        let r = br(1000, &[('M', 100), ('N', 2900), ('M', 120)]);
        let (b, i) = blocks_and_introns(&r);
        assert!(leaving_introns(&b, &i, &chain).is_empty());
        // anchored, intron leaving downstream FROM THIS CHAIN'S OWN DONOR (4120 is an exon end): a candidate
        let r = br(1000, &[('M', 100), ('N', 2900), ('M', 120), ('N', 15_000), ('M', 200)]);
        let (b, i) = blocks_and_introns(&r);
        assert_eq!(leaving_introns(&b, &i, &chain), vec![(4120, 19_120)]);
        // ⭐ §6ge: anchored, but the intron starts BEYOND the chain — the donor belongs to the next gene, so
        // this is that gene's ordinary intron and not a read-through out of this unit
        let r = br(1000, &[('M', 100), ('N', 2900), ('M', 5000), ('N', 15_000), ('M', 200)]);
        let (b, i) = blocks_and_introns(&r);
        assert!(leaving_introns(&b, &i, &chain).is_empty(), "donor is not this chain's exon end");
        // an empty chain is never a source
        assert!(leaving_introns(&b, &i, &[]).is_empty());
    }

    /// `--from-gtf` loci: a locus is a `gene_id` group spanning all its transcripts; its representative is the
    /// transcript with the most reads, ties to the longer span, then to the LAST transcript_id; a gene without
    /// exons is not a locus; the representative's exons come back sorted.
    #[test]
    fn gtf_loci_takes_the_most_read_transcript_ties_to_the_longer_span_then_the_last_id() {
        let gtf = "\
c1\tr\ttranscript\t101\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; reads \"5\";
c1\tr\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
c1\tr\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
c1\tr\ttranscript\t101\t450\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; reads \"5\";
c1\tr\texon\t301\t450\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
c1\tr\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
c1\tr\ttranscript\t51\t90\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T0\"; reads \"1\";
c1\tr\texon\t51\t90\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T0\";
c1\tr\ttranscript\t1001\t1100\t.\t-\t.\tgene_id \"G2\"; transcript_id \"Ta\"; reads \"3\";
c1\tr\texon\t1001\t1100\t.\t-\t.\tgene_id \"G2\"; transcript_id \"Ta\";
c1\tr\ttranscript\t1001\t1100\t.\t-\t.\tgene_id \"G2\"; transcript_id \"Tb\"; reads \"3\";
c1\tr\texon\t1001\t1100\t.\t-\t.\tgene_id \"G2\"; transcript_id \"Tb\";
c1\tr\ttranscript\t5001\t5100\t.\t+\t.\tgene_id \"G3\"; transcript_id \"Tx\"; reads \"9\";
";
        let loci = gtf_loci(std::io::Cursor::new(gtf), Representative::MostReads).unwrap();
        assert_eq!(loci.len(), 2, "G3 has no exon: not a locus");
        assert_eq!(
            loci[0],
            GtfLocus {
                gene_id: "G1".into(),
                chrom: "c1".into(),
                start: 51,
                end: 450,
                rep: "T2".into(),
                rep_reads: 5,
                strand: "+".into(),
                rep_exons: vec![("c1".into(), 101, 200), ("c1".into(), 301, 450)],
            }
        );
        assert_eq!((loci[1].rep.as_str(), loci[1].strand.as_str(), loci[1].rep_reads), ("Tb", "-", 3));
    }

    /// GTF lines of one assembled transcript on `c1`, `+`: the `transcript` line with its `reads`, then one `exon` line per
    /// `(start, end)` (GFF 1-based closed) in the order given.
    fn gtf_tx(gene: &str, tid: &str, reads: u64, exons: &[(u64, u64)]) -> String {
        let (s, e) = (exons.iter().map(|x| x.0).min().unwrap(), exons.iter().map(|x| x.1).max().unwrap());
        let mut out = format!("c1\tr\ttranscript\t{s}\t{e}\t.\t+\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; reads \"{reads}\";\n");
        for (a, b) in exons {
            out += &format!("c1\tr\texon\t{a}\t{b}\t.\t+\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\";\n");
        }
        out
    }

    /// `(gene_id, representative transcript_id)` of every locus of `gtf` under `rule`, in locus order.
    fn reps_of(gtf: &str, rule: Representative) -> Vec<(String, String)> {
        gtf_loci(std::io::Cursor::new(gtf.to_string()), rule).unwrap().into_iter().map(|l| (l.gene_id, l.rep)).collect()
    }

    /// `--representative`: `most-reads` is the default (every product before 2026-10-04), `most-junctions` the opt-in arm,
    /// anything else is refused at parse time.
    #[test]
    fn representative_defaults_to_most_reads_and_accepts_exactly_the_two_rules() {
        let parse = |v: &[&str]| Args::try_parse_from(["mcl_families", "--out", "o"].into_iter().chain(v.iter().copied()));
        assert_eq!(parse(&[]).expect("parse").representative, Representative::MostReads);
        assert_eq!(parse(&["--representative", "most-reads"]).expect("parse").representative, Representative::MostReads);
        assert_eq!(
            parse(&["--representative", "most-junctions"]).expect("parse").representative,
            Representative::MostJunctions
        );
        assert!(parse(&["--representative", "longest-span"]).is_err());
    }

    /// A junction is a gap of >= 50 bp between consecutive exons in coordinate order. GFF 1-based closed: after an exon
    /// ending at 200, an exon starting at 251 leaves bases 201-250, a 50-bp gap (the boundary); 250 leaves 49.
    #[test]
    fn junction_count_is_the_gaps_of_at_least_50_bp_between_consecutive_exons() {
        let ex = |v: &[(u64, u64)]| -> Vec<(String, u64, u64)> { v.iter().map(|&(a, b)| ("c1".to_string(), a, b)).collect() };
        assert_eq!(junction_count(&ex(&[])), 0);
        assert_eq!(junction_count(&ex(&[(101, 200)])), 0, "one exon: no junction");
        assert_eq!(junction_count(&ex(&[(101, 200), (251, 300)])), 1, "a 50-bp gap is a junction");
        assert_eq!(junction_count(&ex(&[(101, 200), (250, 300)])), 0, "a 49-bp gap is not");
        assert_eq!(junction_count(&ex(&[(101, 200), (231, 300)])), 0, "a 30-bp gap is not");
        assert_eq!(junction_count(&ex(&[(101, 200), (201, 300)])), 0, "abutting exons leave no gap");
        assert_eq!(junction_count(&ex(&[(101, 200), (150, 300)])), 0, "overlapping exons leave no gap");
        // coordinate order, not listing order
        assert_eq!(junction_count(&ex(&[(501, 600), (101, 200), (301, 400)])), 2);
        // a short gap is skipped, the two long ones are counted
        assert_eq!(junction_count(&ex(&[(101, 200), (231, 300), (501, 600), (1001, 1100)])), 2);
    }

    /// ⭐ `most-junctions` (5'-truncated libraries, where the most-read transcript is a 3' fragment): a 1-junction transcript
    /// with 10 reads represents its locus under `most-reads` and loses to the 5-junction transcript with 2 reads under
    /// `most-junctions`. Only the representative moves: the locus, its span and every locus with one transcript, or with
    /// no spliced transcript, are the same under both rules.
    #[test]
    fn most_junctions_beats_reads_where_the_most_read_transcript_is_a_fragment() {
        let full = [(101, 200), (301, 400), (501, 600), (701, 800), (901, 1000), (1101, 1200)];
        let gtf = gtf_tx("G1", "Tfrag", 10, &full[4..])
            + &gtf_tx("G1", "Tfull", 2, &full)
            + &gtf_tx("G2", "Tshort", 5, &[(5001, 5100)])
            + &gtf_tx("G2", "Tlong", 2, &[(5001, 5300)]);
        let by_reads = gtf_loci(std::io::Cursor::new(&gtf), Representative::MostReads).unwrap();
        let by_junctions = gtf_loci(std::io::Cursor::new(&gtf), Representative::MostJunctions).unwrap();
        let shape = |l: &GtfLocus| (l.rep.clone(), l.rep_reads, l.rep_exons.len());
        assert_eq!(shape(&by_reads[0]), ("Tfrag".to_string(), 10, 2));
        assert_eq!(shape(&by_junctions[0]), ("Tfull".to_string(), 2, 6));
        assert_eq!((by_junctions[0].start, by_junctions[0].end), (101, 1200));
        assert_eq!(by_reads[0].end, 1200, "the span is the locus's, not the representative's");
        // no spliced transcript in G2: 0 junctions each, so `most-junctions` falls back to reads (then span, then id)
        assert_eq!(by_reads[1], by_junctions[1]);
        assert_eq!(by_junctions[1].rep, "Tshort");
        assert_eq!(by_reads.len(), by_junctions.len());
    }

    /// `most-junctions` ties go to the most `reads`, then the longer span, then the LAST `transcript_id` in sorted order
    /// (never file order). Every locus holds a decoy `Tz_<gene>` that wins on reads, span and id but has 1 junction against
    /// the candidates' 2, so each winner has to beat it on junctions first and the lower keys only among the tied.
    #[test]
    fn most_junctions_ties_go_to_reads_then_span_then_the_last_id() {
        let decoy = |g: &str| gtf_tx(g, &format!("Tz_{g}"), 100, &[(1, 100), (901, 1000)]);
        let ex3 = |a: u64, z: u64| vec![(a, 200), (301, 400), (501, z)]; // 2 junctions, span z - a
        let gtf = [
            // G1: Ta has more reads but a shorter span and the earlier id than Tb: reads decide
            gtf_tx("G1", "Ta", 9, &ex3(101, 600)),
            gtf_tx("G1", "Tb", 4, &ex3(51, 650)),
            decoy("G1"),
            // G2: equal reads, so the longer span decides (Tc, the earlier id)
            gtf_tx("G2", "Tc", 5, &ex3(101, 700)),
            gtf_tx("G2", "Td", 5, &ex3(101, 600)),
            decoy("G2"),
            // G3: everything ties; the sorted-last id wins, whatever order the file lists them in
            gtf_tx("G3", "Tf", 5, &ex3(101, 600)),
            gtf_tx("G3", "Th", 5, &ex3(101, 600)),
            gtf_tx("G3", "Te", 5, &ex3(101, 600)),
            decoy("G3"),
        ]
        .concat();
        let want: Vec<(String, String)> = [("G1", "Ta"), ("G2", "Tc"), ("G3", "Th")].iter().map(|&(g, t)| (g.into(), t.into())).collect();
        assert_eq!(reps_of(&gtf, Representative::MostJunctions), want);
        // `most-reads` takes the decoy everywhere: the junction count is the one thing it does not look at
        assert!(reps_of(&gtf, Representative::MostReads).iter().all(|(_, t)| t.starts_with("Tz_")));
    }

    /// A gap under 50 bp is not a junction: a transcript whose exons are split by 30-bp gaps has none, so a one-junction
    /// transcript with fewer reads represents the locus under `most-junctions`.
    #[test]
    fn most_junctions_does_not_count_a_30_bp_gap() {
        let gtf = gtf_tx("G1", "Tgap", 9, &[(101, 200), (231, 330), (361, 460)]) // 30-bp gaps: 0 junctions
            + &gtf_tx("G1", "Tone", 1, &[(101, 200), (401, 500)]); // one 200-bp intron: 1 junction
        assert_eq!(reps_of(&gtf, Representative::MostJunctions), vec![("G1".to_string(), "Tone".to_string())]);
        assert_eq!(reps_of(&gtf, Representative::MostReads), vec![("G1".to_string(), "Tgap".to_string())]);
    }

    /// ⭐ The `--from-gtf --emit-units` copy table is read by `copy_assign --families/--copies-fa` exactly as the
    /// legacy catalog is: header columns 1-11 are the catalog's; every row passes `parse_copies_tsv`; the FASTA
    /// passes `parse_copies_fa` and `to_colocated`; the sequence is the representative's spliced exon sum
    /// (reverse-complemented on `-`); copies are in genomic order; `max_family_identity` is the best direct
    /// edge; an exon-overlapping member of the same family is merged; `.` strand is written as `+`; the locus
    /// extent is the de novo locus span clipped at the neighbouring copy's chain.
    #[test]
    fn locus_rep_copies_are_in_the_catalog_contract() {
        use rustle::vg_family::annotation_families::HomologyGraph;
        use rustle::vg_family::catalog_input::{group_families, parse_copies_fa, parse_copies_tsv, to_colocated};
        let dir = tempfile::tempdir().unwrap();
        let mut x: u64 = 12345;
        let mut rnd = |n: usize| -> String {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[(x >> 62) as usize] as char
                })
                .collect()
        };
        let (c1, c2) = (rnd(4000), rnd(2000));
        let fasta = dir.path().join("g.fa");
        std::fs::write(&fasta, format!(">c1\n{c1}\n>c2 description\n{c2}\n")).unwrap();
        let locus = |g: &str, chrom: &str, s: u64, e: u64, rep: &str, reads: u64, st: &str, ex: &[(u64, u64)]| GtfLocus {
            gene_id: g.into(),
            chrom: chrom.into(),
            start: s,
            end: e,
            rep: rep.into(),
            rep_reads: reads,
            strand: st.into(),
            rep_exons: ex.iter().map(|&(a, b)| (chrom.to_string(), a, b)).collect(),
        };
        let loci = vec![
            locus("G1", "c1", 101, 900, "T1", 7, "+", &[(101, 200), (301, 400)]),
            locus("G2", "c1", 2001, 2900, "T2", 4, "-", &[(2001, 2100), (2301, 2500)]),
            locus("G3", "c2", 101, 600, "T3", 0, ".", &[(101, 600)]),
            locus("G5", "c1", 151, 260, "T5", 9, "+", &[(151, 260)]), // shares exon bases with G1
            locus("G4", "c1", 3001, 3500, "T4", 2, "+", &[(3001, 3500)]), // no family
        ];
        let key = |l: &GtfLocus| -> GeneKey { (l.chrom.clone(), l.start, l.end) };
        let mut g = HomologyGraph::default();
        g.genes = vec![key(&loci[0]), key(&loci[1]), key(&loci[2]), key(&loci[3])];
        g.idents.insert((0, 1), 0.95);
        g.idents.insert((0, 2), 0.90);
        g.idents.insert((0, 3), 0.99);
        let clusters = vec![Cluster {
            members: vec![key(&loci[1]), key(&loci[0]), key(&loci[2]), key(&loci[3])],
            density: 1.0,
            frac_in: 1.0,
            corroborated: None,
        }];
        let out = dir.path().join("fam").to_string_lossy().to_string();
        let st = write_locus_rep_copies(&out, &fasta.to_string_lossy(), &clusters, &g, &[], &loci, None, true, true, false)
            .unwrap();
        assert_eq!((st.copies, st.families, st.multi_copy_families, st.merged, st.unstranded), (3, 1, 1, 1, 1));
        let tsv = std::fs::read_to_string(format!("{out}.copies.tsv")).unwrap();
        let header: Vec<&str> = tsv.lines().next().unwrap().split('\t').collect();
        assert_eq!(
            header[..11],
            ["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "exons", "max_family_identity"]
        );
        let copies = parse_copies_tsv(&tsv).unwrap();
        let got: Vec<(&str, usize, &str, u64, u64, char, u32)> =
            copies.iter().map(|c| (c.tid.as_str(), c.copy_idx, c.chrom.as_str(), c.start, c.end, c.strand, c.n_reads)).collect();
        assert_eq!(
            got,
            vec![("T1", 0, "c1", 100, 400, '+', 7), ("T2", 1, "c1", 2000, 2500, '-', 4), ("T3", 2, "c2", 100, 600, '+', 0)]
        );
        assert_eq!(copies[0].exons, vec![(100, 200), (300, 400)]);
        // the L2 extent: the locus span, clipped at the neighbouring copy's chain on the contig
        assert_eq!(copies.iter().map(|c| c.locus.unwrap()).collect::<Vec<_>>(), vec![(100, 900), (2000, 2900), (100, 600)]);
        let ident: Vec<&str> = tsv.lines().skip(1).map(|l| l.split('\t').nth(10).unwrap()).collect();
        assert_eq!(ident, vec!["0.990000", "0.950000", "0.900000"], "best direct edge, merged member's included");
        let merged = std::fs::read_to_string(format!("{out}.copies.merged.tsv")).unwrap();
        assert_eq!(merged.lines().nth(1).unwrap(), "MCL0\tc1:151-260\tT5\tc1:101-900\tT1");
        let fa = std::fs::read_to_string(format!("{out}.copies.fa")).unwrap();
        let seqs = parse_copies_fa(&fa).unwrap();
        let rc = |s: &str| -> String {
            s.bytes().rev().map(|b| match b { b'A' => 'T', b'C' => 'G', b'G' => 'C', _ => 'A' }).collect()
        };
        assert_eq!(String::from_utf8(seqs[&("MCL0".to_string(), 0)].seq.clone()).unwrap(), format!("{}{}", &c1[100..200], &c1[300..400]));
        assert_eq!(String::from_utf8(seqs[&("MCL0".to_string(), 1)].seq.clone()).unwrap(), rc(&format!("{}{}", &c1[2000..2100], &c1[2300..2500])));
        let genome = rustle::genome::GenomeIndex::from_fasta(&fasta.to_string_lossy()).unwrap();
        let fams = group_families(copies).unwrap();
        assert_eq!(fams.len(), 1);
        let (cf, _) = to_colocated(&fams[0], Some(&seqs), &genome).unwrap();
        assert_eq!(cf.copies.len(), 3);
        let regions = std::fs::read_to_string(format!("{out}.copies.regions")).unwrap();
        assert_eq!(regions, "MCL0\tc1:1-7500\nMCL0\tc2:1-5600\n");
    }

    /// The target is the unit whose chain contains the first base AFTER the intron; the source never
    /// matches itself, and a junction landing between units resolves to nothing.
    #[test]
    fn readthrough_target_is_the_unit_holding_the_first_base_after_the_intron() {
        let units = vec![
            ("c1".to_string(), vec![(1000u64, 1100u64), (4000, 4120)]),
            ("c1".to_string(), vec![(19_100, 19_400)]),
            ("c2".to_string(), vec![(19_100, 19_400)]),
        ];
        assert_eq!(readthrough_target((4120, 19_120), "c1", 0, &units), Some(1));
        // the same coordinates on another contig do not match
        assert_eq!(readthrough_target((4120, 19_120), "c3", 0, &units), None);
        // landing in a gap
        assert_eq!(readthrough_target((4120, 18_000), "c1", 0, &units), None);
        // the source unit is never its own target
        assert_eq!(readthrough_target((1000, 4010), "c1", 0, &units), None);
    }

    /// L2: the extent spans the kept segment of every PRIMARY record with a block in the chain, through
    /// ordinary introns; secondaries and supplementaries never count; a record spliced over the chain with no
    /// block in it contributes nothing; a giant unsupported intron is cut (mis-chain rule).
    #[test]
    fn read_extent_is_the_union_of_kept_segments_of_primaries_with_a_block_in_the_chain() {
        let chain = [(1000u64, 1100u64)];
        assert_eq!(read_extent(&[], &chain, 3), None);
        assert_eq!(read_extent(&[br(1000, &[('M', 100)])], &[], 3), None);
        let mut reads = vec![br(1000, &[('M', 100)]), br(500, &[('M', 50), ('N', 3000), ('M', 20)])];
        // intron 550-3550 spans the chain with no block in it: the first read alone counts
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 1100)));
        // a primary WITH a block in the chain extends it through its 3-kb intron
        reads[1].read.ref_start = 1050; // 1050-1100 N 4100-4120
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 4120)));
        // ...but not when it is a secondary or a supplementary (-F 2308)
        reads[1].is_secondary = true;
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 1100)));
        reads[1].is_secondary = false;
        let mut supp = br(900, &[('S', 40), ('M', 150)]);
        supp.is_supplementary = true;
        reads.push(supp);
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 4120)));
        // a giant (> 50 kb) intron supported by one record is cut; the same intron in three records is kept
        reads.push(br(1090, &[('M', 20), ('N', 80_000), ('M', 30)]));
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 4120)));
        reads.push(br(1090, &[('M', 20), ('N', 80_000), ('M', 30)]));
        reads.push(br(1090, &[('M', 20), ('N', 80_000), ('M', 30)]));
        assert_eq!(read_extent(&reads, &chain, 3), Some((1000, 81_140)));
    }
}

/// One de novo locus of an assembled GTF (`--from-gtf`): a `gene_id` group, its span (GFF 1-based, min/max over the
/// exons of all its transcripts) and its REPRESENTATIVE, picked by [`Representative`] (`--representative`). `MostReads`,
/// the default: the transcript with the most `reads`, ties to the longer span, then to the lexicographically last
/// `transcript_id` (the order `max_by_key` has always resolved ties in). `MostJunctions`: the transcript with the most
/// [`junction_count`] junctions, then the same key (reads, longer span, last `transcript_id`). The span does not depend on
/// the rule. The representative's exons are the locus's exons in `loci.gff3` and, with `--emit-units`, the copy's exons in
/// `copies.tsv`: the "positional exon sum" (read-derived coordinates, genome bases).
#[derive(Clone, Debug, PartialEq)]
struct GtfLocus {
    gene_id: String,
    chrom: String,
    start: u64,
    end: u64,
    rep: String,
    rep_reads: u64,
    /// The representative's strand column, verbatim (`.` when the transcript line had none).
    strand: String,
    /// The representative's exons, GFF 1-based closed, sorted by start (stable).
    rep_exons: Vec<(String, u64, u64)>,
}

/// The junction floor of `--representative most-junctions`, in bp: the gap between two consecutive exons that makes a
/// junction. The floor of the strict "found" rule (`bench/copy_support.py` `MIN_INTRON`), not a new constant.
const MIN_JUNCTION_GAP: u64 = 50;

/// Junctions of one transcript, for `--representative most-junctions`: the gaps of at least [`MIN_JUNCTION_GAP`] bp between
/// consecutive exons in coordinate order. Exons are GFF 1-based closed, so the gap after `(_, e1)` and before `(s2, _)` is
/// the `s2 - e1 - 1` bases between them; abutting or overlapping exons leave none.
fn junction_count(exons: &[(String, u64, u64)]) -> usize {
    let mut iv: Vec<(u64, u64)> = exons.iter().map(|x| (x.1, x.2)).collect();
    iv.sort_unstable();
    iv.windows(2).filter(|w| w[1].0.saturating_sub(w[0].1 + 1) >= MIN_JUNCTION_GAP).count()
}

/// Parse an assembled GTF into its loci, representatives picked by `rule`, in first-appearance order of `gene_id` (genes
/// without exons are left out).
/// ⚠ `loci.gff3`, and through it the whole family graph, is written from exactly this list: any change here must
/// be cmp-checked on the families stage products.
fn gtf_loci<R: std::io::BufRead>(reader: R, rule: Representative) -> Result<Vec<GtfLocus>> {
    use std::collections::{BTreeMap, HashMap, HashSet};
    fn attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
        let pat = format!("{key} \"");
        let i = s.find(&pat)? + pat.len();
        let j = s[i..].find('"')? + i;
        Some(&s[i..j])
    }
    let mut exons: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
    let mut gene_of: HashMap<String, String> = HashMap::new();
    let mut strand: HashMap<String, String> = HashMap::new();
    let mut reads: HashMap<String, u64> = HashMap::new();
    let mut gene_order: Vec<String> = Vec::new();
    // genes already in `gene_order` (was a scan of every `gene_of` value per transcript line: O(T^2) on a
    // whole-genome GTF); same first-appearance order for any GTF whose transcript_ids do not switch gene
    let mut seen_genes: HashSet<String> = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 {
            continue;
        }
        let Some(t) = attr(r[8], "transcript_id") else { continue };
        if r[2] == "transcript" {
            let g = attr(r[8], "gene_id").unwrap_or(t).to_string();
            if seen_genes.insert(g.clone()) {
                gene_order.push(g.clone());
            }
            gene_of.insert(t.to_string(), g);
            strand.insert(t.to_string(), r[6].to_string());
            reads.insert(t.to_string(), attr(r[8], "reads").and_then(|v| v.parse().ok()).unwrap_or(0));
        } else if r[2] == "exon" {
            exons.entry(t.to_string()).or_default().push((r[0].to_string(), r[3].parse()?, r[4].parse()?));
        }
    }
    let mut txs_of: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (t, g) in &gene_of {
        txs_of.entry(g.clone()).or_default().push(t.clone());
    }
    let mut out = Vec::new();
    for g in &gene_order {
        let Some(ts) = txs_of.get(g) else { continue };
        let all: Vec<&(String, u64, u64)> = ts.iter().flat_map(|t| exons.get(t).into_iter().flatten()).collect();
        if all.is_empty() {
            continue;
        }
        let chrom = all[0].0.clone();
        let (s, e) = (all.iter().map(|x| x.1).min().unwrap(), all.iter().map(|x| x.2).max().unwrap());
        let span_of = |t: &String| exons.get(t).map(|v| v.iter().map(|x| x.2).max().unwrap() - v.iter().map(|x| x.1).min().unwrap()).unwrap_or(0);
        let mut sorted_ts = ts.clone();
        sorted_ts.sort();
        // `max_by_key` keeps the LAST maximum: over the sorted ids, a full tie goes to the last `transcript_id` under both
        // rules (`most-junctions` only puts the junction count in front of the `most-reads` key)
        let rep = match rule {
            Representative::MostReads => sorted_ts.iter().max_by_key(|t| (reads.get(*t).copied().unwrap_or(0), span_of(t))),
            Representative::MostJunctions => sorted_ts.iter().max_by_key(|t| {
                (exons.get(*t).map_or(0, |v| junction_count(v)), reads.get(*t).copied().unwrap_or(0), span_of(t))
            }),
        }
        .unwrap()
        .clone();
        let st = strand.get(&rep).cloned().unwrap_or_else(|| ".".into());
        let mut ex = exons.get(&rep).cloned().unwrap_or_default();
        ex.sort_by_key(|x| x.1);
        out.push(GtfLocus {
            gene_id: g.clone(),
            chrom,
            start: s,
            end: e,
            rep_reads: reads.get(&rep).copied().unwrap_or(0),
            rep,
            strand: st,
            rep_exons: ex,
        });
    }
    Ok(out)
}

/// Write `<out>.loci.fa`: one record `>CONTIG:START-END` + the genome's forward strand over that span, per locus span.
/// With `hash`, also the [`ContentHash`](rustle::vg_family::run_cache::ContentHash) of every byte written (the PAF
/// cache key), taken as the bytes go out so the file is never read back.
fn write_loci_fa(
    genome: &rustle::genome::GenomeIndex,
    spans: &[(String, u64, u64)],
    fa_path: &str,
    fasta: &str,
    hash: bool,
) -> Result<Option<rustle::vg_family::run_cache::ContentHash>> {
    use rustle::vg_family::run_cache::HashingWriter;
    let mut fa = HashingWriter::new(std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(fa_path)?), hash);
    for (c, s, e) in spans {
        let seq = genome.fetch_sequence(c, s - 1, *e).with_context(|| format!("{c}:{s}-{e} not in {fasta}"))?;
        writeln!(fa, ">{c}:{s}-{e}")?;
        fa.write_all(&seq)?;
        writeln!(fa)?;
    }
    fa.flush()?;
    Ok(fa.hash)
}

/// The families PAF cache key (`paf/<fnv(key)>/key.tsv`): the minimap2 command line and build, and the content
/// hash + byte length of the loci FASTA it aligns. v2 (2026-09-28): the hash is the word-wise 128-bit
/// `ContentHash` taken while the FASTA is written; v1 was a byte-wise FNV-1a 64 of the FASTA read back from disk.
fn families_paf_key(cmd: &str, minimap2_version: &str, loci_fa: &rustle::vg_family::run_cache::ContentHash) -> String {
    format!(
        "rustle families paf v2\ncmd\t{cmd}\nminimap2\t{minimap2_version}\nquery_hash\tcontent128:{}\nquery_bytes\t{}\n",
        loci_fa.hex(),
        loci_fa.len()
    )
}

/// `--from-gtf`: the de novo locus set of an assembled GTF, as the family stage consumes it (see the flag doc).
/// Returns `(loci.gff3, loci.fa, loci.paf)` paths and the loci themselves (for `--emit-units`' copy table).
fn loci_from_gtf(
    gtf: &str,
    fasta: &str,
    out: &str,
    threads: usize,
    rule: Representative,
) -> Result<(String, String, String, Vec<GtfLocus>)> {
    use std::collections::HashSet;
    use std::io::Write;
    let f = std::fs::File::open(gtf).with_context(|| format!("opening {gtf}"))?;
    let loci = gtf_loci(std::io::BufReader::new(f), rule)?;
    let gff3 = format!("{out}.loci.gff3");
    let fa_path = format!("{out}.loci.fa");
    let paf = format!("{out}.loci.paf");
    let mut g3 = std::fs::File::create(&gff3)?;
    writeln!(g3, "##gff-version 3")?;
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    for l in &loci {
        let (chrom, s, e, st, g) = (&l.chrom, l.start, l.end, &l.strand, &l.gene_id);
        writeln!(g3, "{chrom}\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{g};Name={g}")?;
        for (_, a, b) in &l.rep_exons {
            writeln!(g3, "{chrom}\t.\texon\t{a}\t{b}\t.\t{st}\t.\tParent=gene-{g};gene={g}")?;
        }
        spans.push((chrom.clone(), s, e));
    }
    use rustle::vg_family::run_cache as rc;
    let root = rc::cache_root();
    let contigs: HashSet<String> = spans.iter().map(|x| x.0.clone()).collect();
    let genome = rustle::genome::GenomeIndex::from_fasta_contigs(fasta, &contigs)?;
    let fa_hash = write_loci_fa(&genome, &spans, &fa_path, fasta, root.is_some())?;
    eprintln!("[mcl_families] --from-gtf: {} loci -> all-vs-all", spans.len());
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let mm_args: Vec<String> = ["-x", "asm20", "-c", "-X", "-N", "50", "-p", "0.1", "--secondary=yes", "-t"]
        .iter()
        .map(|s| s.to_string())
        .chain(std::iter::once(threads.to_string()))
        .collect();
    // PAF cache (`RUSTLE_CACHE_DIR`, see `rustle::vg_family::run_cache`): keyed by EVERY byte of the loci FASTA
    // (hashed while it was written above, never re-read), the command line and the minimap2 build; a hit hard-links
    // the cached PAF to `<out>.loci.paf` instead of re-aligning (and instead of copying it). The entry is pinned, so a
    // write through the link is a miss on the next run, never a stale replay.
    let paf_entry = root.zip(fa_hash).map(|(root, h)| {
        let key = families_paf_key(&format!("{mm2} {}", mm_args.join(" ")), &rc::minimap2_version(&mm2), &h);
        rc::Entry::new(&root, "paf", key).pinned()
    });
    let paf_path = std::path::Path::new(&paf);
    if let Some(e) = paf_entry.as_ref().filter(|e| e.is_hit()) {
        // a replay that fails (another run replacing the entry) falls through to running minimap2
        if let Ok(linked) = e.replay("out.paf", paf_path) {
            eprintln!(
                "[cache] all-vs-all PAF replayed from {} ({}; minimap2 skipped)",
                e.dir.display(),
                if linked { "hard link" } else { "copy" }
            );
            return Ok((gff3, fa_path, paf, loci));
        }
    }
    // never truncate in place: `<out>.loci.paf` may be a hard link to a cache entry from an earlier run
    match std::fs::remove_file(paf_path) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
        Err(err) => return Err(err).with_context(|| format!("removing {paf}")),
    }
    let out_paf = std::fs::File::create(&paf)?;
    let status = std::process::Command::new(&mm2)
        .args(&mm_args)
        .arg(&fa_path)
        .arg(&fa_path)
        .stdout(out_paf)
        .stderr(std::process::Stdio::null())
        .status()
        .with_context(|| format!("running {mm2}"))?;
    anyhow::ensure!(status.success(), "minimap2 all-vs-all failed");
    if let Some(e) = paf_entry.as_ref() {
        let stored = e.staging().and_then(|st| {
            e.stage_link(&st, "out.paf", paf_path)?;
            e.commit(&st)
        });
        if let Err(err) = stored {
            eprintln!("[cache] could not store the PAF ({err:#}); continuing");
        }
    }
    Ok((gff3, fa_path, paf, loci))
}

#[cfg(test)]
mod paf_cache_key_tests {
    use super::*;
    use rustle::vg_family::run_cache as rc;

    /// The families PAF key covers every byte of the loci FASTA (hashed as it is written, equal to a re-read of the
    /// file): the same loci FASTA written again (a new mtime, a later run) keeps the key, one changed base inside a
    /// locus span changes it at the same file size, and an entry committed under one key is a miss for the other.
    #[test]
    fn families_paf_key_changes_with_any_loci_fasta_byte_and_only_then() {
        let dir = std::env::temp_dir().join(format!("rustle_fam_paf_key_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let g = dir.join("g.fa");
        let spans = vec![("c1".to_string(), 3u64, 12u64), ("c2".to_string(), 1, 8)];
        let key_of = |genome_text: &str, name: &str| -> (String, Vec<u8>) {
            std::fs::write(&g, genome_text).unwrap();
            let genome = rustle::genome::GenomeIndex::from_fasta(g.to_str().unwrap()).unwrap();
            let fa = dir.join(name);
            let h = write_loci_fa(&genome, &spans, fa.to_str().unwrap(), "g.fa", true).unwrap().unwrap();
            let bytes = std::fs::read(&fa).unwrap();
            assert_eq!(rc::ContentHash::of_file(&fa).unwrap().hex(), h.hex(), "in-stream hash == hash of the file");
            assert_eq!(h.len(), bytes.len() as u64);
            // without a cache nothing is hashed and the file is the same
            assert!(write_loci_fa(&genome, &spans, fa.to_str().unwrap(), "g.fa", false).unwrap().is_none());
            assert_eq!(std::fs::read(&fa).unwrap(), bytes);
            (families_paf_key("minimap2 -x asm20 -t 4", "2.28-r1209", &h), bytes)
        };
        let (k1, b1) = key_of(">c1\nACGTACGTACGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(b1, b">c1:3-12\nGTACGTACGT\n>c2:1-8\nGGGGCCCC\n");
        std::thread::sleep(std::time::Duration::from_millis(20));
        let (k1b, b1b) = key_of(">c1\nACGTACGTACGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!((&k1, &b1), (&k1b, &b1b), "rewritten, same bytes: same key");
        let (k2, b2) = key_of(">c1\nACGTACGTTCGTAC\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(b1.len(), b2.len());
        assert_ne!(k1, k2, "one base inside a span, same size: another key");
        let (k3, _) = key_of(">c1\nACGTACGTACGTAT\n>c2\nGGGGCCCCAA\n", "a.loci.fa");
        assert_eq!(k1, k3, "a base outside every span leaves the loci FASTA, and so the key, unchanged");
        let root = dir.join("cache");
        let e1 = rc::Entry::new(&root, "paf", k1.clone()).pinned();
        let st = e1.staging().unwrap();
        std::fs::write(dir.join("p.paf"), b"c1:3-12\t10\n").unwrap();
        e1.stage_link(&st, "out.paf", &dir.join("p.paf")).unwrap();
        e1.commit(&st).unwrap();
        assert!(rc::Entry::new(&root, "paf", k1).pinned().is_hit_verify(true));
        assert!(!rc::Entry::new(&root, "paf", k2).pinned().is_hit_verify(false), "a changed loci FASTA is a miss");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
