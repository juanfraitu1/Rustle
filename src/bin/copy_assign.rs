//! De-novo multi-copy family DETECTION + per-read copy ASSIGNMENT CLI.
//!
//! Annotation-free, read-coherence pipeline: scan a region of a coordinate-sorted BAM, assemble de-novo
//! transcripts, detect co-located paralog families, and assign each read — including the hard multimappers
//! minimap2 leaves at MAPQ 0. A family IS one VARIATION GRAPH (copies = PATHS, PSV columns = BUBBLES); each
//! read is THREADED through it and scored to its maximum-likelihood copy-path — by its PSV bases + copy-
//! specific junctions — then significance-gated (assign-or-abstain, never 1/k). `--phase` emits the
//! materialized GFA with the reads threaded through it as walks (the Canzar shared-evidence flip, visualized).
//!
//! Writes `<out>.families.tsv` (per-family roster + two-pass + unique-mapper agreement stats) and
//! `<out>.assignments.tsv` (per-read copy assignment). A `.bai` next to the BAM makes the region read fast.

use anyhow::{Context, Result};
use clap::Parser;
use lru::LruCache;
use std::collections::HashSet;
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use rustle::genome::GenomeIndex;
use rustle::vg_family::denovo_assemble::longest_orf;
use rustle::vg_family::absent_copy::DnaNeedsRecord;
use rustle::vg_family::linearize::LinearizeCertificate;
use rustle::vg_family::copy_assign::{AssignParams, AssignStatus};
use rustle::vg_family::em_copy_assign::em_assign_family;
use rustle::vg_family::denovo_assemble::{
    assemble_gate, assemble_gate_census, merge_fuzzy_skeletons, pass1_skeletons_widened, reads_in_region,
    tied_secondary_reads_in_region,
    BamIndexCache, BamRead, GATE_MIN_READS,
};
use rustle::vg_family::catalog_input::{
    group_families, parse_copies_fa, parse_copies_tsv, to_colocated, CatalogFamily, SeqIndex,
};
use rustle::vg_family::denovo_pipeline::{
    catalog_overlaps, detect_and_assign, ColocatedFamily, DenovoConfig, FallbackEdge, FamilyAssignment,
    OverlapKind,
};
use rustle::vg_family::family_detect::collapse_loci_groups;
use rustle::vg_family::read_conflict::{as_evidence, AsEvidence};
use rustle::vg_family::readonly_copy_number::{chi_h_with_junctions, depth_cn};

/// Read one GTF attribute out of an attribute string (`key "value";`). Used by `--productivity` to recover
/// the family and copy it already wrote, rather than threading them separately.
fn re_attr(attrs: &str, key: &str) -> Option<String> {
    let pat = format!("{key} \"");
    let i = attrs.find(&pat)? + pat.len();
    let j = attrs[i..].find('"')? + i;
    Some(attrs[i..j].to_string())
}

/// One assembled isoform (FLAIR-style intron-chain collapse), kept for the optional `--gtf` emit. `gene_tid`
/// is the locus this isoform collapses into (shared-junction gene); a family copy is its own gene, so a
/// family-copy gene_tid matches a `copy_tid` in the assignment output and is tagged `multicopy` in the GTF.
struct TranscriptRec {
    tid: String,
    gene_tid: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
    introns: Vec<(u64, u64)>,
    n_reads: u32,
}

/// The expensive, independent per-region work (BAM read + `detect_and_assign`, which holds the dominant
/// poasta alignment). Computed possibly in parallel across regions (any contig), then drained SERIALLY in the
/// original region order so `CAFAM` ids and every output row stay byte-identical to the serial sweep.
///
/// The heavy read SEQUENCES (`BamRead`) are dropped inside the worker — the output stage needs only the read
/// NAMES (everything else lives in `fams`) — so collecting all regions' results out-of-order is lightweight.
struct RegionWork {
    contig: String,
    lo: u64,
    hi: u64,
    read_names: Vec<String>,
    /// Primary-alignment MAPQ per read, parallel to `read_names`. Kept (the read sequences are dropped) so the
    /// output stage can compute the tie-break invariance certificate: `mapq > 0` = a unique mapper whose copy
    /// support survives any primary/secondary relabeling. See `anchored_support`.
    read_mapqs: Vec<u8>,
    /// Per-RECORD placement `(ref_start, ref_end, flags)`, parallel to `read_names`, with `flags` bit0 =
    /// secondary and bit1 = supplementary. `read_names` alone names the MOLECULE; this names the
    /// alignment RECORD a family's claim rests on, which is what distinguishes a molecule that SPANS two
    /// loci (one record, claimed by both families) from one placed independently at two disjoint loci
    /// (two records) — the whole of the cross-family reconciliation rule. Kept for the same reason
    /// `read_mapqs` is: the heavy `BamRead` sequences are dropped in the worker.
    read_spans: Vec<(u64, u64, u8)>,
    /// Aligned blocks per record (0-based half-open; `N` closes a block, `D` extends it), parallel to
    /// `read_names`. Kept so the output stage can say whether a read has an aligned BASE inside a copy —
    /// a read spliced OVER a copy is no evidence for it (§6es hygiene; ledger §6cm).
    read_blocks: Vec<Vec<(u64, u64)>>,
    /// Each record's OWN chromosome, parallel to `read_names`. For every region this binary ever swept
    /// before cross-chromosome family support, every record here shares `contig` — so this used to be
    /// redundant with the struct's own `contig` field and nothing read it. A cross-chromosome family's
    /// `RegionWork` pools records from several chromosomes into one work unit, and `contig` (kept for
    /// logging/back-compat) is then just one of them — `xfam_pass1`'s cross-family reconciliation needs
    /// each record's real chromosome, not the work unit's label, so it reads this instead.
    read_chrom: Vec<String>,
    /// Genomic strand per read, parallel to `read_names` — `ts:A` (transcript strand relative to the READ)
    /// flipped by alignment orientation (`BamRead::ts`'s own doc), falling back to the read's own FLAG 0x10
    /// when minimap2 emitted no `ts` (an unspliced read has no junction motif to read it from). Kept for
    /// `--rescue-singletons` (§6hw), whose emitted transcripts otherwise have no strand source of their own
    /// once `BamRead` is dropped.
    read_strand: Vec<char>,
    /// Alignment-score evidence per read, parallel to `read_names`. Reported, never decisive — see
    /// `read_conflict::AsEvidence` for why raw AS is length-confounded and `de` makes the call.
    as_ev: Vec<AsEvidence>,
    n_mapped: usize,
    fams: Vec<FamilyAssignment>,
    fallback: Vec<FallbackEdge>,
    dna_needs: Vec<DnaNeedsRecord>,
    /// Augment-and-linearize certificates (Task 4), one per Stage-2-admitted reference-absent copy:
    /// `(family_id, certificate, (chrom, start, end))`. Empty unless the opt-in is on (`--linearize` or
    /// `--linearize-gate`); written to `<out>.linearize.tsv` (Task 5) when set. Under `--linearize-gate`,
    /// a non-LINEARIZES verdict also demotes that candidate out of `admitted` in `detect_and_assign` (it
    /// appears here as a certificate row but not as an admitted copy).
    linearize_certs: Vec<(String, LinearizeCertificate, (String, u64, u64))>,
    transcripts: Vec<TranscriptRec>, // FLAIR-style isoforms for the --gtf emit (empty unless --gtf)
    /// `--gtf-copy-set`: the primaries the AS-tied gate dropped (unique / clear-best mappers) as
    /// (chrom, start, end, intron chain) — the EVIDENCE that places an isoform at a copy (§6hn)
    uniq_reads: Vec<(String, u64, u64, Vec<(u64, u64)>)>,
    /// O3: this family's raw (uncorrected) missing-copy pair statistics. Empty unless `--flag-missing-copies`.
    o3_raw_pairs: Vec<rustle::vg_family::missing_copy_flag_pass::RawPair>,
    /// O3: candidate orphan-read loci outside every unit of this family. Empty unless `--flag-missing-copies`.
    o3_orphan_loci: Vec<rustle::vg_family::missing_copy_flag_pass::OrphanLocus>,
    /// Read-seeded copy discovery: candidate new copies clustered from AS-tied reads' out-of-catalog
    /// placements. Empty unless `--discover-copies`. Report only (Task 4 drains this to
    /// `<out>.discovered_copies.tsv`) -- never feeds back into this run's own catalog or assignments.
    discovered: Vec<rustle::vg_family::copy_discovery::DiscoveredCopy>,
    /// `--union-certificate`: what the union pass did in this region (its side-file rows + counts). Default
    /// (empty) unless the flag is on -- the pass runs INSIDE the worker because it needs the read sequences.
    union: rustle::vg_family::denovo_pipeline::UnionSummary,
}

#[derive(Parser, Debug)]
#[command(
    name = "copy_assign",
    about = "De-novo multi-copy family detection + per-read copy assignment (PSV + copy-specific junctions)"
)]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast indexed region query).
    #[arg(long)]
    bam: String,
    /// Genome FASTA (with a `.fai` so only the needed contig is loaded).
    #[arg(long)]
    fasta: String,
    /// A single region to scan as `chrom:start-end`.
    #[arg(long)]
    region: Option<String>,
    /// A regions FILE (one `chrom:start-end` per line; extra whitespace-separated columns ignored) for a
    /// genome-wide sweep. Regions are grouped by contig so each contig is loaded ONCE.
    #[arg(long)]
    regions: Option<String>,
    /// §6zb: sweep EVERY contig with mapped reads (from the BAM index), each as one region, in one process —
    /// meant for `--assemble-only`, whose streaming pass-1 keeps memory at O(distinct chains) per contig, so
    /// no batching script is needed. Mutually exclusive with `--region`/`--regions`.
    #[arg(long, default_value_t = false)]
    genome_wide: bool,
    /// §6zb: force the historical read-materialising path under `--assemble-only` (every record kept as
    /// `PrimaryRead` + `BamRead`, the O2 AS-tied gate run). Same transcripts; the only output difference is
    /// the O2-only `matched_reads` attribute (tied reads per transcript), which the streaming path reports as
    /// 0 because it never evaluates ties. Costs ~2.4 GB per million records and ~4× the wall time.
    #[arg(long, default_value_t = false)]
    materialize_reads: bool,
    /// Output prefix; writes `<out>.families.tsv` and `<out>.assignments.tsv`.
    #[arg(long)]
    out: String,
    /// ⭐ §6m5 / PREREG `docs/PREREG_assembler_widening_2026-09-18.md`: READ-ISOFORM WIDENING for the
    /// `--gtf` assembly. The assembler groups reads by EXACT intron chain and filters on that chain's own
    /// read count, so a junction carried by many reads spread over many chains produces no surviving group
    /// and the junction disappears (measured on NPIP: dropped junctions have a median largest-chain of 2
    /// reads vs 32 for kept ones; five dropped junctions carry 255-287 reads over 117-134 chains). With
    /// `k > 0` a chain is ALSO admitted when every one of its junctions has >= k reads supporting it,
    /// counted per junction over all spliced reads in the region. Chains are never concatenated. Port of
    /// `shared_definition::widen_with_read_isoforms`, which ships at k=5. **Default 0 = OFF, byte-identical.**
    #[arg(long, default_value_t = 0)]
    read_isoform_k: u32,

    /// ⭐ ASSEMBLER-ONLY MODE (§6p6). Skip family detection, homology refinement and copy assignment
    /// entirely, and run ONLY the assembly path: reads -> intron-chain skeletons -> `assemble_gate` ->
    /// `collapse_loci_groups` -> GTF. This is the "define loci and cluster isoforms into transcripts"
    /// product, with none of the all-vs-all comparison work.
    ///
    /// It implies `--gtf` (the assembly IS the output) and leaves `<out>.families.tsv` /
    /// `<out>.assignments.tsv` empty by construction — there is no assignment in this mode, and an empty
    /// file is the honest record of that rather than a missing one.
    ///
    /// Composes with the assembly knobs: `--read-isoform-k`, `RUSTLE_JUNCTION_MAJORITY` (default ON since
    /// 2026-09-21; `=0` restores the old strict-canonicity behaviour), `RUSTLE_GTF_SECONDARY`,
    /// `RUSTLE_GATE_CENSUS`.
    #[arg(long, default_value_t = false)]
    assemble_only: bool,

    /// ⭐ §6za JUNCTION MODE for the `--assemble-only` transcript product
    /// (`docs/PREREG_assembly_precision_levers_2026-09-23.md`): `strict` = every junction of a transcript
    /// must be canonical (GT-AG / GC-AG / AT-AC) on one strand; `majority` = the §6m8 family-recovery rule
    /// that tolerates a minority of short non-canonical junctions when a canonical majority fixes the
    /// strand. The majority rule was adopted for FAMILY recovery (NPIPB12); measured on the transcript
    /// product it admits chains that are annotation-exact at 1.2% (human chr20-22, 9.7% of the output)
    /// and 0.0% (gorilla NC_073244.2), against 16-37% for the rest (register 1069). Default `strict`
    /// under `--assemble-only`; ignored otherwise (the family paths keep `RUSTLE_JUNCTION_MAJORITY`'s
    /// default). An explicit `RUSTLE_JUNCTION_MAJORITY` in the environment always wins.
    #[arg(long, default_value = "strict")]
    assembly_junctions: String,

    /// ⭐ ASSEMBLY POLISH (§6p8). Post-assembly precision filters over the emitted GTF, using only the
    /// `reads "N"` attribute — no reference and no annotation, so this is legal in de novo mode.
    ///
    /// * `none` (default) — byte-identical to the unfiltered emit.
    /// * `mono` — the MONO-EXONIC SUPPORT FLOOR. A single-exon transcript carries no junction evidence at
    ///   all, so it must reach the upper quartile of the multi-exon read support in the same run
    ///   (`--polish-mono-quantile`). Self-tuning: the threshold is read off this run's own distribution.
    ///   Measured cost on chr20 and held-out chr11: **zero matching intron chains, zero sensitivity**;
    ///   gain: transcript precision 35.6 → 43.6 (chr20) and 34.8 → 41.5 (chr11).
    /// * `full` — `mono` plus the SUPPORT-AWARE ISM COLLAPSE: drop a transcript whose intron chain is a
    ///   contiguous sub-chain of another's on the same contig/strand, unless it carries at least as much
    ///   read support as its container. Adds ~6 more precision points but costs chains (345 → 337 on
    ///   chr20, 715 → 690 on chr11), so it is a deliberate recall/precision trade, not a free win.
    #[arg(long, default_value = "none", value_parser = ["none", "mono", "full"])]
    assembly_polish: String,

    /// Quantile of the multi-exon `reads` distribution used as the mono-exonic support floor under
    /// `--assembly-polish`. 0.75 reproduces the validated setting.
    #[arg(long, default_value_t = 0.75)]
    polish_mono_quantile: f64,

    /// §6q2 3'-ANCHORED ISM for `--assembly-polish full`. Collapse a sub-chain fragment ONLY when it
    /// shares the container's 3' end — i.e. its chain is a suffix of the container's on `+` and a prefix
    /// on `-`. §6p4 localised the artifact as 5' truncation specifically (TBC1D3's 5'UTR 75.3% covered
    /// against 100% of its 3'UTR), so a fragment that is 3'-anchored is the truncation signature, while a
    /// fragment sitting in the middle of a chain, or sharing only the 5' end, is a distinct short isoform
    /// the collapse should not absorb.
    ///
    /// ⛔**REFUTED as a default (register row 861):** 19-21 of 30 cells against 28/30 with the
    /// unrestricted collapse. Mid-chain and 5'-anchored fragments are junk at a similar rate to
    /// 3'-anchored ones, so restricting the collapse to the truncation signature keeps the rest.
    #[arg(long, default_value_t = false)]
    polish_ism_3p: bool,

    /// ⭐ §6q6 FUZZY JUNCTION TOLERANCE for `--assembly-polish` — the mechanism `isoseq collapse` calls
    /// `--max-fuzzy-junction` (its default is 5). Two intron chains count as the SAME chain when they have
    /// the same number of junctions and every corresponding donor/acceptor is within N bp.
    ///
    /// ⛔**REFUTED as a default at isoseq's own setting (register row 866).** 0 and 2 bp are 28/30 cells
    /// against StringTie; **5 bp collapses to 15/30, costing 81 matching intron chains** over six
    /// chromosomes for +0.24 points of transcript precision. Our pass-1 already enforces canonical GT-AG
    /// motifs, so a sub-5 bp junction difference that survives into our GTF is a real tandem splice site
    /// rather than alignment noise — on chr20 the modal offset among the pairs 5 bp merges is **3 bp, the
    /// NAGNAG signature** — and the merge deletes whichever variant has less read support, which is not
    /// always the annotated one. `isoseq collapse` needs the tolerance because it collapses raw
    /// alignments without a motif constraint; we do not.
    ///
    /// Duplicates are merged into the best-supported member, whose `reads` attribute is left untouched
    /// (the polish never rewrites support, only removes rows). By default the tolerance applies ONLY to
    /// that near-duplicate merge; `--polish-fuzzy-ism` extends it to the ISM sub-chain test as well.
    #[arg(long, default_value_t = 0)]
    polish_fuzzy_junction: i64,

    /// Extend `--polish-fuzzy-junction` to the ISM sub-chain containment test, which is what `isoseq
    /// collapse` does. Measured as a near no-op on top of the merge — the two halves cost the same
    /// (register row 866), so the damage is the merge itself, not the containment test.
    #[arg(long, default_value_t = false)]
    polish_fuzzy_ism: bool,

    /// ⭐ §6r5 TPM IN THE GTF. Add `cov` and `TPM` to every transcript line, from the `reads` support the
    /// assembler already records. **Count-based, NOT length-normalised**: a long read is one molecule, so
    /// `TPM_i = reads_i / sum(reads) * 1e6`. Measured on chr20 against StringTie's own TPM over the 507
    /// intron chains both tools call: count-based **rho 0.879**, length-normalised only **0.714** — so the
    /// short-read convention is the wrong one here and is not offered. Against FLAIR's isoform counts,
    /// rho 0.755 (388 shared chains). `cov` is reads per kb of spliced length, for orientation only.
    ///
    /// Default off: with it unset the GTF is byte-identical to a run without this flag.
    #[arg(long, default_value_t = false)]
    gtf_tpm: bool,

    /// §6r2 ABSOLUTE FLOOR EXEMPTION for `--polish-isoform-fraction`: never drop a transcript carrying at
    /// least this many reads for being a minor fraction of its locus. Distinct from
    /// `--polish-fraction-exempt`, whose bar is the run's `--polish-mono-quantile` of multi-exon support
    /// (≈12 reads on a deep library) and which therefore recovers nothing here. 0 = off.
    #[arg(long, default_value_t = 0)]
    polish_fraction_min_reads: u64,

    /// ⭐ §6za RETAINED-INTRON FILTER for `--assembly-polish` (`docs/PREREG_assembly_precision_levers_2026-09-23.md`).
    /// Drop a transcript when a junction of ANOTHER transcript at its locus (same `gene_id`, same strand)
    /// lies strictly inside one of its exons and that junction's support — the reads of every surviving
    /// transcript at the locus that carries it — is at least this many times the transcript's own reads.
    /// The simulation of r1070 showed the aligner reading through a short exon consistently across the
    /// reads of one molecule, which yields exactly such a chain at ≥ 2 reads; genuine retained-intron
    /// isoforms carry a substantial share of their locus and survive at a high ratio. Reads-only, so legal
    /// de novo. Held-out gorilla at 10: intron-chain precision 33.1 → 34.3, matching chains −0.42%, the
    /// dropped set 3.4% annotation-exact against 34.2% for the kept set. **Default 10 (2026-09-23,
    /// user's call)**; 0 = off (the pre-2026-09-23 output).
    #[arg(long, default_value_t = 10.0)]
    polish_retained_ratio: f64,

    /// ⭐ §6q4 FRACTION EXEMPTION for `--assembly-polish`. Exempt a transcript from
    /// `--polish-isoform-fraction` when its OWN read support reaches the run's `--polish-mono-quantile`
    /// level of multi-exon support. The fraction test is purely relative, so at a deep locus it discards
    /// isoforms that carry plenty of absolute evidence merely because the dominant isoform is deeper
    /// still; the exemption makes the filter say "a minor flow AND thinly supported" instead of "a minor
    /// flow". It reuses the existing quantile, so it adds no new constant.
    #[arg(long, default_value_t = false)]
    polish_fraction_exempt: bool,

    /// ⭐ §6q3 ISM SUPPORT RATIO for `--assembly-polish full`. Keep a sub-chain fragment when its read
    /// support reaches this fraction of its container's. 1.0 (the default) demands parity; a genuine
    /// shorter isoform typically carries a substantial share of its locus, while a 5'-truncation artifact
    /// carries a small one, so the ratio separates them **independently of library depth** — unlike
    /// `--polish-ism-escape`, whose absolute bar rises with coverage.
    #[arg(long, default_value_t = 1.0)]
    polish_ism_ratio: f64,

    /// §6q1 ISM ABSOLUTE ESCAPE for `--assembly-polish full`. Keep a sub-chain fragment when its own read
    /// support reaches the run's `--polish-mono-quantile` level of multi-exon support, even if its
    /// container is deeper still. The ISM collapse gets harsher as coverage grows, and the deepest
    /// chromosome measured (chr5, 519,887 records) is the only one where it pushes matching intron chains
    /// below StringTie's (473 vs 476).
    ///
    /// ⚠**Not a good default, and measured as such (register row 858):** turning it on recovers chr5's
    /// chains (473 → 476, a tie) but costs chr11 its precision lead (50.3 → 49.9 intron-chain, 50.1 → 49.7
    /// transcript), for a net 27/30 against 28/30 with it off. It is a deliberate recall/precision dial
    /// with both endpoints measured, not an improvement. Reuses the existing quantile, so it adds no new
    /// constant.
    #[arg(long, default_value_t = false)]
    polish_ism_escape: bool,

    /// ⭐ §6q0 SHADOW RULE for `--assembly-polish`: a single-exon transcript that lies in a spliced gene's
    /// shadow is that gene's unspliced, intronic or antisense signal, not an independent single-exon gene.
    /// Unlike the ISM host rule it needs neither containment nor a support test. Drop a single-exon
    /// transcript when it
    ///   * overlaps any EXON of a multi-exon transcript, on EITHER strand (a mono read pile has no splice
    ///     motif, so its strand label carries no evidence and antisense overlap is not informative of
    ///     independence), or
    ///   * overlaps the SPAN of a multi-exon transcript on the SAME strand (its own gene's introns).
    ///
    /// Measured on chr20 (development): of 175 single-exon predictions that match no reference transcript
    /// it removes a large share, and of the single-exon predictions that survive the read floor AND match
    /// a reference transcript it removes **none** — the one chr20 survivor has only an anti-strand SPAN
    /// overlap, which is deliberately not a criterion.
    #[arg(long, default_value_t = false)]
    polish_mono_shadow: bool,

    /// ⭐ §6p9 LOCUS ISOFORM FRACTION for `--assembly-polish`. Drop a transcript whose read support is
    /// below this fraction of the BEST-supported transcript at the same `gene_id`; the dominant isoform of
    /// a locus is never dropped. This is StringTie's `-f` criterion (fraction of the locus maximum), and it
    /// is the lever that closes the class-`j` gap — on held-out chr11 we emitted 589 "novel junction
    /// combination" transcripts against StringTie's 481, which is the whole of our precision deficit there.
    ///
    /// Distinct from `--min-isoform-fraction`, which is a fraction of the locus TOTAL and only tags
    /// `low_confidence`; this one removes rows. Default 0.0 = off.
    #[arg(long, default_value_t = 0.0)]
    polish_isoform_fraction: f64,

    /// Minimum copies for a co-located family. Two-copy homologous families are the majority and were
    /// invisible to assignment at the old default of 3; lowering it to 2 changes default family detection
    /// on its own, independently of `--homology-primary`.
    #[arg(long, default_value_t = 2)]
    min_copies: usize,
    /// poasta→minimap2 fallback length cap (bp) for the INTRON-RETENTION PSV discovery
    /// (`discover_intron_psvs`, opt-in via `RUSTLE_INTRON_PSV=1`): above this a copy pair's genomic span
    /// uses minimap2, not poasta (which is exact but O(n^2)). Threaded to the use site via `RUSTLE_POA_CAP`
    /// (an env var, not a struct field — the const sits many call frames below `main`; see the module doc at
    /// its use site). Default 20000 matches the prior hard-coded constant exactly, so leaving this flag
    /// unset is byte-identical to before it existed.
    #[arg(long, default_value_t = 20_000)]
    poa_cap: usize,
    /// Per-family multimapper read-pool cap (`o2_materialize::READ_CAP`, `MaterializeConfig::read_cap`).
    /// EXPOSED HERE FOR AUDITABILITY ONLY: `o2_materialize` is a Rust byte-parity port of the Python
    /// genome-wide-catalog materializer (`bench/o2_vg_visualization.py::materialize_family`) that no
    /// `src/bin/*.rs` binary — including this one — imports, so this flag is currently a NO-OP in
    /// `copy_assign` (parses so `RUSTLE_READ_CAP`/CLI usage never hard-errors; a non-default value warns at
    /// startup rather than silently doing nothing). Default 6000 matches the constant.
    #[arg(long, default_value_t = 6_000)]
    read_cap: usize,
    /// Co-located window (bp): copies must cluster within this span.
    #[arg(long, default_value_t = 5_000_000)]
    win: u64,
    /// BAM-reading threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// REGION-PARALLEL sweep: process this many regions (independent families) concurrently, ACROSS contigs.
    /// Each family is an independent unit (`detect_and_assign` is pure; `BamIndexCache`/genome are read-only),
    /// so the heavy per-family poasta alignments of N families — even on different chromosomes — overlap on N
    /// cores; a bounded LRU cache of loaded contig genomes (capacity ≈ N) avoids reloads and caps resident
    /// chromosomes. Output is collected and `CAFAM` ids assigned in the SAME (serial) order afterward, so the
    /// result is BYTE-IDENTICAL to the serial run. Peak memory ≈ N regions' reads + ≈N contig genomes (raise N
    /// for speed on a many-core box, lower it if memory-bound; the documented genome-wide OOM is why this is
    /// opt-in). `1` (default) = the exact serial path. The pool is sized to N and shared with the inner
    /// per-copy alignment parallelism. Speedup ceiling = the single heaviest family (already internally
    /// parallel); on a full genome-wide sweep this approaches ~N×.
    #[arg(long, default_value_t = 1)]
    region_threads: usize,
    /// FLAIR-LIKE ASSEMBLY emit: also write `<out>.gtf` — every de-novo isoform in the swept regions as a
    /// transcript+exon GTF (IGV-loadable), grouped into genes by shared junctions. Multi-copy family copies
    /// appear as separate genes at their own loci, tagged `family_id`/`copy_index`/`multicopy "true"`; all
    /// other loci are ordinary isoforms (`multicopy "false"`). Annotation-free (intron-chain collapse + the
    /// canonical-junction gate; no short-read junction correction). Independent of the assignment outputs;
    /// default off. Pair with `bench/igv_tracks.py` for the copy-coloured reads.
    #[arg(long, default_value_t = false)]
    gtf: bool,
    /// ⭐ §6hn/§6ho (PREREG 95409846): with `--gtf`, emit the GTF O2 BELIEVES. Family isoforms are grouped
    /// across copies by lifting their intron chains through copy-to-copy alignments (minimap2 asm20 on the
    /// copy spans); a group's addresses are the copies with EVIDENCE (a unique mapper whose primary lies
    /// there, or a certificate-assigned read); transcripts at copies without evidence are dropped (phantoms),
    /// a certified read whose copy has no transcript gets a LIFTED placement, and an isoform with no evidence
    /// anywhere is emitted ONCE with `copies "A,B[,outside]"` = the AS-tied placement set of its reads.
    /// Non-family and single-intron transcripts pass through. ⭐ DEFAULT ON (user, 2026-09-09 §6hp);
    /// `--no-gtf-copy-set` restores the aligner-placed GTF byte-for-byte.
    #[arg(long, default_value_t = true)]
    gtf_copy_set: bool,
    /// Escape for the 2026-09-09 default: the aligner-placed GTF with the old copy attributes (pre-§6hp).
    #[arg(long, default_value_t = false)]
    no_gtf_copy_set: bool,
    /// ⭐ (`docs/OPEN_ITEMS_2026-09-09.md`, 09-10, user request): with `--gtf`, write
    /// `<out>.read_provenance.tsv` — one row per AS-tied alignment record, naming exactly which transcript
    /// it contributed to (gate-passed or rescued) or exactly why it did not (excluded_ambiguous,
    /// excluded_tied, excluded_chain_at_unassigned_locus, excluded_supplementary, ...). "Complete, not
    /// necessarily good": every record this region's O2 saw gets exactly one row. Default off.
    #[arg(long, default_value_t = false)]
    read_provenance: bool,
    /// ⭐ B2 (`docs/OPEN_ITEMS_2026-09-09.md`, 09-10): with `--gtf`, emit a low-support (`support "N"`,
    /// N < 3) transcript for any certificate-ASSIGNED read whose own exon chain does not otherwise reach
    /// `assemble_gate`'s min-reads floor — a real, resolvable read the GTF currently drops on the floor
    /// entirely (measured: 220/290 such contested-but-uncarried molecules are singletons, §6hh row 798).
    /// Default off, byte-identical when unset.
    #[arg(long, default_value_t = false)]
    rescue_singletons: bool,
    /// O3 (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`): port of
    /// `bench/missing_copy_flag_pass.py`'s missing-copy detector. Requires `--families`. Adds `missing_copy_flag`/`missing_copy_class`/
    /// `missing_copy_rate_per_kb`/`missing_copy_p`/`missing_copy_n_rejected` columns to `<out>.family_join.tsv` and writes
    /// `<out>.missing_copy_loci.tsv`. Default off, byte-identical when unset.
    #[arg(long, default_value_t = false)]
    flag_missing_copies: bool,
    /// O3: Bonferroni alpha for the genome-wide missing-copy threshold (`alpha / n_pairs`).
    #[arg(long, default_value_t = 0.001)]
    missing_copy_alpha: f64,
    /// O3: cap on reads realigned per side (test/control) per candidate copy, for wall-clock control.
    #[arg(long, default_value_t = 500)]
    missing_copy_max_reads: usize,
    /// ⭐ A3 (`docs/OPEN_ITEMS_2026-09-09.md`, 09-10): append `sibling_identity` and `n_cols_vs_sibling` to
    /// `.assignments.tsv` — the whole-family PSV identity between an assigned copy and the competitor
    /// governing its p_value, and how many distinguishing positions this read spans against it. Reporting
    /// only; changes no decision. Default off, byte-identical schema when unset.
    #[arg(long, default_value_t = false)]
    sibling_report: bool,
    /// ⭐ §6u6: append the EICHLER-STYLE AS-margin call beside ours, for the head-to-head the advisor
    /// asks for. The rule he cites: assign a multi-mapper to its best alignment iff no rival scores
    /// within T units of it, else discard as ambiguous. Two columns are appended:
    /// `eichler_call` (`assign`/`discard`) and `eichler_same_copy` (1/0/NA — whether OUR assigned copy
    /// is the read's best-AS placement, decidable only where both rules assign).
    ///
    /// ⚠ The two rules barely share a subject and the counts are NOT comparable without saying so:
    /// under the default AS-tied gate every surviving read has margin 0, so Eichler discards 100% of
    /// them by construction. Run with `--no-as-tied-only` for the honest comparison — measured on the
    /// Y ampliconic genes, 69.8% of reads have margin 0 and his rule discards 83.7% of the
    /// multi-mapping population, which is precisely O2's subject (`docs/EICHLER_COMPARISON_2026-09-21.md`).
    ///
    /// ⚠ T is a convention, not a constant: 4,687 assignments at T=1 vs 1,588 at T=20 on that substrate.
    /// Reporting only; changes no decision. Default off, byte-identical schema when unset.
    #[arg(long)]
    eichler_margin: Option<i32>,
    /// ⭐ §6hq/PREREG fd894558: StringTie-style per-locus relative-depth demotion. Per `gene_tid`
    /// (the locus `collapse_loci_groups` already assigns), `isoform_fraction = n_reads / max(n_reads in the
    /// locus)`; below this floor a transcript is tagged `low_confidence "true"` and, under `--gtf-copy-set`,
    /// excluded from evidence/grouping (it never seeds or extends a copy's span) — StringTie's flow
    /// decomposition reports the identical readthrough at copy 4 as a minor isoform (coverage 6.1 vs the
    /// dominant 50.5); we collapsed by exact chain with one flat `min_reads` floor, so a 2-read and a
    /// 50-read chain at the same locus were both just "a transcript." Default `0.0` = off, byte-identical.
    #[arg(long, default_value_t = 0.0)]
    min_isoform_fraction: f64,
    /// ⭐ §6hr/PREREG 35e290a9: the locus-boundary outlier test. Per `gene_tid` group, bucket members by
    /// their outer boundary (50 bp); a transcript in the FARTHEST bucket is `boundary_low_confidence` iff
    /// (a) that bucket sits > `--min-boundary-gap` bp past the next-farthest bucket, AND (b) the farthest
    /// bucket's read total is < this fraction of the group's total reads. Narrower than
    /// `--min-isoform-fraction` (which over-flagged ordinary heterogeneity, row 808): BOTH depth AND a clear
    /// positional gap are required, matching the copy-4 readthrough (18.4 kb gap, 2/41 = 4.9 % reads) without
    /// catching routine alternative-TSS/TES variants that cluster closely. ⭐ DEFAULT `0.10` (user,
    /// 2026-09-10 §6hs): measured clean on both species (human 1.03 %, gorilla 0.5 % flagged, 0 false
    /// positives on hand inspection at the paired `--min-boundary-gap` default). `--min-boundary-fraction 0`
    /// is the escape (byte-identical to pre-§6hs).
    #[arg(long, default_value_t = 0.10)]
    min_boundary_fraction: f64,
    /// Minimum gap (bp) from a group's farthest boundary bucket to its next-farthest, before
    /// `--min-boundary-fraction` even considers flagging it. ⭐ 5000, not the originally pre-registered 1000
    /// (PREREG 35e290a9 outcome): at 1000, 2 of 9 human MCL0 flags were ordinary smooth-tail heterogeneity
    /// (a continuum of alternative termini whose last two points happened to sit > 1000 bp apart by chance,
    /// not an isolated outlier) — copies 9 and 25, gaps 1400/1150 bp. At 5000 both drop out and the
    /// remaining 7 are all clean, isolated single-or-few-transcript outliers 10.5–40 kb from a well-
    /// supported majority cluster.
    #[arg(long, default_value_t = 5000)]
    min_boundary_gap: u64,
    /// Lift tolerance (bp) per intron boundary when matching an isoform across copies (`--gtf-copy-set`).
    #[arg(long, default_value_t = 5)]
    gtf_lift_tol: u64,
    /// ⭐ PRODUCTIVITY (§6gp), with `--gtf`: every transcript gains `orf_aa` — the longest ORF over its
    /// strand-oriented exon-sum — and `productive`, true when that ORF reaches **half the best ORF among the
    /// transcripts assigned to the same family**. The bar is relative on purpose: an absolute amino-acid cut
    /// discarded 26 % of protein-coding units when it was tried on the core rule (§6gb, register 746).
    /// Also writes `<out>.productivity.tsv`, one row per copy. flair predicts productivity per ISOFORM; doing
    /// it per COPY is the part flair cannot reach, because it cannot attribute an isoform to a copy.
    /// Default off.
    #[arg(long, default_value_t = false)]
    productivity: bool,
    /// poasta memory threshold (bp) for POA homology confirmation. A candidate family pair whose larger
    /// transcript exceeds this is confirmed via the linear-memory longest-common-substring FALLBACK instead of
    /// poasta (which OOMs on long sequences); those edges are recorded in `<out>.fallback.tsv`. Lower it (e.g.
    /// 8000) on dense/large-gene regions to keep poasta off the big operands. Default 20000 matches the python.
    #[arg(long, default_value_t = 20_000)]
    max_poa_len: usize,
    /// Recover COLLAPSED copies: feed AS-tied SECONDARY reads (a copy whose reads minimap2 flagged secondary
    /// because it picked a sibling as primary) into the rescue, so the starved copy can clear the support gate.
    /// Additive to the rescue only; default OFF (primary-only, byte-identical).
    #[arg(long, default_value_t = false)]
    recover_copies: bool,
    /// AS-tie ratio for --recover-copies: a secondary counts only if its AS >= ratio * the read's best AS
    /// (1.0 = exact tie; 0.98 admits a 2% margin). Guards against homology-shadow spillover.
    #[arg(long, default_value_t = 0.98)]
    as_ratio: f64,

    /// ⭐ O2 SCOPE (user, 2026-09-09; DEFAULT ON by user decision the same day): copy assignment runs on
    /// **AS-TIED MULTIMAPPERS ONLY** — the molecules the aligner could not place (≥2 placements in the region
    /// with runner-up AS ≥ `--as-tie-ratio` × best). Everything else is dropped **before the certificate**:
    /// no read-star, no assignment, no row. The advisor's definition: "O2 takes all the ambiguous reads —
    /// same AS, a coin toss which is best — and infers their copy; uniquely mapped transcripts are
    /// irrelevant to O2." Measured before the gate existed, 4,969 of gorilla MCL1's 6,142 "assigned"
    /// molecules were single-placement uncontested reads (§6gv, PREREG_as_tied_only_2026-09-09.md).
    ///
    /// The `--families` contract check and the `--gtf` assembly still see every read: a GTF is judged on the
    /// hard transcripts but is not forbidden the easy ones (ADVISOR_QUESTIONS Part 0h). ⚠ Region-local:
    /// placements on other contigs are invisible, so a molecule tied only across contigs is (conservatively)
    /// dropped. This flag is the ESCAPE: it restores the pre-2026-09-09 behaviour byte-for-byte.
    #[arg(long)]
    no_as_tied_only: bool,

    /// Tie width for the AS-tied gate. 1.0 = exact tie (the runner-up scores exactly the best, i.e. the
    /// aligner's primary pick was a coin toss); 0.98 admits a 2 % margin. Reported at both widths before
    /// any default is proposed.
    #[arg(long, default_value_t = 1.0)]
    as_tie_ratio: f64,

    /// ⭐ §6hd (user 2026-09-09): widen the gate to ALIGNER SELF-DISAGREEMENT — admit a molecule whose PRIMARY
    /// record sits in one supplied-family unit while its best-AS record sits in a DIFFERENT one, even with no
    /// AS tie. minimap2's chaining stage (which picks the primary) and its base-level scoring (AS) disagree
    /// about the copy: by construction the aligner cannot reliably tell, the "coin toss" of O2's definition
    /// in a stronger form. Found while auditing the origin-drop-indels assignments: a real 58 bp insertion
    /// costs ~50–70 AS in gap penalties, so the TRUE copy (the primary) scores below a wrong copy; when the
    /// wrong copy's AS is unique there is no tie and the gate skipped the read as a "clear best". Human MCL0:
    /// 2,278 such molecules, 64 % carrying ≥ 50 bp of insertion in the primary placement, and `--as-tie-ratio`
    /// does not reach them (0.98 admits 42 %). Default OFF; a gate-only column `aligner_disagreement` marks them.
    ///
    /// ⚠⚠ A7 (`docs/OPEN_ITEMS_2026-09-09.md`, register row 791): REFUTED for the reason it was built. PREREG
    /// 58a5496b predicted this recovers the §6hc blind spot as real assignments; measured 33% assigned
    /// (<40% predicted), median margin 13.8 (<20 predicted) — most of the widening lands on PSV dead heats
    /// (margin 0, 830 of the admitted set) or a copy worse than the primary by AS, not the SV-shape mechanism
    /// it was named for. Kept as a flag (not deleted) only so the PREREG's own negative result stays
    /// reproducible on request; do not enable it expecting the row-306/791 mechanism to fire.
    #[arg(long)]
    admit_aligner_disagreement: bool,
    /// Also dump the per-read PSV GENOTYPE MATRIX — `<out>.psv_reads.tsv` (each read's base at every PSV column
    /// + its assignment), `<out>.psv_copies.tsv` (each copy's PSV alleles), `<out>.psv_cols.tsv` (column →
    /// genome position). The raw per-molecule evidence behind each assignment, for the proof visualization.
    #[arg(long, default_value_t = false)]
    dump_psv: bool,
    /// Cluster AS-tied reads' out-of-catalog placements into candidate new copies, written to
    /// `<out>.discovered_copies.tsv`. Report only -- never mutates the input catalog or this run's own
    /// assignments (two-pass: inspect the report, append accepted rows to the catalog by hand, re-run).
    /// Default off; unset, output is byte-identical to a run without this flag.
    #[arg(long, default_value_t = false)]
    discover_copies: bool,
    /// ⭐ UNION CERTIFICATE (register rows 1092/1093, 2026-09-24). With `--families`, every AS-tied molecule
    /// whose tied placements touch copies of >= 2 supplied families and/or >= 1 locus outside every supplied
    /// unit gets ONE certificate over the UNION of those candidates — all copies of every family touched,
    /// plus one pseudo-copy per outside locus built from `--fasta` over the placement's aligned blocks — and
    /// the verdict is applied to every family's row for it: the winning family's row `assigned` (with the
    /// union's evidence), every other family's row `tied`; an outside winner or a union abstention ties /
    /// abstains every row. Cures the per-family table's cross-family double claims (904 foreign `assigned`
    /// rows on the chr16 truth simulation) and the §6gz demotion of correct votes whose outside tie partner
    /// was merely never scored. A molecule with a tied placement in another region (its primary record not
    /// loaded in the worker) is counted and left as today. Requires `--families` and the AS-tied gate.
    /// Writes `<out>.union_certificate.tsv` (+ a `union_certificate` row in `params.tsv`); default off =
    /// every existing output byte-identical. What the pass changes is `fa.assignments` (every status emit
    /// site) and the touched families' `assigned_*`/`resolvable_*` counters; what it does NOT recompute:
    /// `quant.tsv` `abundance`/`ci` and the `--prior abundance` weighting (the per-family EM, pre-union),
    /// `n_hard` (counts `best_copy` regardless of status, so a union-tied row still counts at its per-family
    /// best copy), `uniq`/`uniq_agree`, `--dump-star` proofs and `readthrough_into` (registered by the
    /// per-family pass). Read `n_soft` and the statuses under this flag, not those columns.
    #[arg(long, default_value_t = false)]
    union_certificate: bool,
    /// One-flag IGV bundle: implies `--dump-psv` (the PSV genotype matrix), so a subsequent
    /// `bench/igv_tracks.py --assignments <out>.assignments.tsv --bam <bam> --regions <regions> --out <out>`
    /// emits `<out>.tagged.bam` (reads coloured by assigned copy), `<out>.copies.bed`, and `<out>.psv.vcf`
    /// (PSVs as an IGV variant track, copies as samples) — everything IGV needs to SEE each assignment.
    #[arg(long, default_value_t = false)]
    igv: bool,
    /// Decisive-margin τ: the minimum log-likelihood-ratio over the runner-up copy to call a read ASSIGNED
    /// (else AMBIGUOUS); the single calibrated knob that replaces the vote-count (min_psv, margin) integers.
    /// τ = ln((1−p)/p) for a target per-read misassignment p (τ=6.9 default ≈ p 1e-3, the Eichler AS≥10
    /// analog; τ=2.0 ≈ p 0.12). The identifiability gate (n_decisive≥1) is independent of τ and always applied.
    #[arg(long, default_value_t = 6.9)]
    margin: f64,
    /// Per-base PSV error rate used in the likelihood when a read carries no per-base quality (HiFi ~0.003).
    #[arg(long, default_value_t = 0.003)]
    error_rate: f64,
    /// IsoCon significance level α for the DEFAULT gate: assign iff the per-read certificate p < α/(n−1)
    /// (Bonferroni over the n−1 competitors) AND the read is the strict MLE. α is the FAMILY-WIDE
    /// misassignment rate over assigned reads (1e-3 ≈ the τ=6.9 / Eichler AS≥10 precision point). Ignored
    /// when `--margin-gate` is set.
    #[arg(long, default_value_t = 1e-3)]
    alpha: f64,
    /// Use the LEGACY τ-margin gate (the `--margin` knob) instead of the IsoCon significance gate. For
    /// reproducing pre-significance-gate numbers and the gate A/B comparison.
    #[arg(long, default_value_t = false)]
    margin_gate: bool,
    /// Disable the RNA-editing filter (Clair3-RNA): by default, A↔G PSV columns showing within-copy
    /// heterogeneity are flagged as A-to-I editing sites and downweighted in the certificate so an edited
    /// base cannot fake copy-support. This reverts to trusting every PSV column uniformly.
    #[arg(long, default_value_t = false)]
    no_editing_filter: bool,
    /// ⭐ O2-8a (§6ej): abstain on junction/PSV CONFLICT — a read whose own splice junctions fit another copy
    /// strictly better than its PSV-best copy is `ambiguous`, never `assigned` (on NPIP's LCR16a cores 4 of 11
    /// junction-anchored assignments contradicted the junction at min_p 1e-16..1e-35). Writes
    /// `<out>.conflicts.tsv`. Likelihoods untouched; default OFF ⟹ byte-identical
    #[arg(long, default_value_t = false)]
    junction_conflict_abstain: bool,

    /// ⭐ §6eu / register 689: the read-support PSV filter (`read_supported_columns`: keep a copy-vs-copy column
    /// only if two alleles each reach two reads). It cannot tell a mis-assembled base from an UNEXPRESSED
    /// paralogue — both are monomorphic pileups — and on ZSCAN5 it kept 4 of 216 columns and produced 11
    /// confident wrong calls. **Default OFF (user decision 2026-09-05)**; on the 35-family sweep OFF raised
    /// assignments 18 % and MAPQ-60 agreement 95.3 → 96.0 %. `--psv-read-filter` turns it back on;
    /// `RUSTLE_PSV_READFILTER=0|1`, when set, overrides the flag (the pre-2026-09-05 escape).
    #[arg(long, default_value_t = false)]
    psv_read_filter: bool,

    /// ⭐ O2-9 / D3 (PREREG adj/d3, register 691): pool every BAM record of a molecule into ONE observation
    /// vector — "the read is the star" (§6fa): the molecule's sequence aligned to every copy's unit, its own
    /// columns, an origin certificate — instead of scoring each BAM record and abstaining on contradiction
    /// (which masked wrong column positions, rows 700–702). **Default ON (user decision 2026-09-05)**;
    /// `--no-molecule-observations` restores the record-level path byte-for-byte.
    #[arg(long, default_value_t = true)]
    molecule_observations: bool,
    /// Escape hatch: the record-level observation path (every catalog before 2026-09-05 §6fa).
    #[arg(long, default_value_t = false)]
    no_molecule_observations: bool,
    /// ⭐ §6fc: the origin certificate counts SUBSTITUTIONS only (indels = isoform structure, not origin). Higher
    /// yield (NPIP 175 → 442 assigned) at 1 wrong of 20 audited anchors and 98.6 % MAPQ-60 agreement; the
    /// default counts every edit and made no wrong anchor call. Opt-in.
    #[arg(long, default_value_t = false)]
    origin_substitutions_only: bool,
    /// ⭐ §6hc (traced from the `best_by_psv` fix, §6ha): in the GENOMIC origin certificate, drop indels
    /// from the edit count — test substitutions + unaligned bases only, keeping the "must explain the
    /// whole read" requirement `--origin-substitutions-only` drops. 20 of 40 human MCL0 molecules newly
    /// origin-rejected by the `best_by_psv` fix share one best copy with X=2–5 substitutions (inside error
    /// tolerance) but I=57–66 — a recurring ~60 bp indel, not evidence of the wrong copy; measured 32/40
    /// flip REJECT→pass, the other 8 (dominated by 14–327 unaligned bases) correctly stay rejected. Default
    /// off (byte-identical); takes a back seat to `--origin-substitutions-only` if both are set.
    /// ⭐ DEFAULT ON (user, 2026-09-09 step 2, §6hf): sequence-confirmed on both mechanisms, excision 97.4 % /
    /// 100 %, indel PSV columns measured separately and kept off. `--no-origin-drop-indels` is the escape.
    #[arg(long, default_value_t = true)]
    origin_drop_indels: bool,
    /// Escape for the 2026-09-09 default: the origin certificate counts I+D bases again (pre-§6hf, byte-identical).
    #[arg(long, default_value_t = false)]
    no_origin_drop_indels: bool,
    /// ⭐ A6 (`docs/OPEN_ITEMS_2026-09-09.md`, register row 790, 09-10): with `--gtf-copy-set`, name an
    /// undecided isoform's outside tie partner by its own locus (`outside:chrom:start-end`) instead of the
    /// bare `outside` bareword — e.g. the EIF3C/EIF3CL class, 8,944 human molecules, previously unnamed and
    /// unactionable. Default off, byte-identical `copies_undecided` schema when unset.
    #[arg(long, default_value_t = false)]
    name_outside_tie: bool,
    /// ⭐ A2 (`docs/OPEN_ITEMS_2026-09-09.md`, register row 799): before calling a read `origin_rejected`
    /// (foreign to the family), also test every other candidate, not only the PSV-duel winner. A read that
    /// fails the certificate at its best copy but is origin-consistent with another candidate is reported
    /// `Ambiguous` instead of `origin_rejected` — both abstain, only the stronger label changes. Default
    /// off (byte-identical); this WIDENS which reads get the milder label, moving the `origin_rejected`
    /// count and the contested denominator that's built from it (measured: 27 human reads).
    #[arg(long, default_value_t = false)]
    origin_consistency_check: bool,
    /// ⭐ PREREG 021446fb: indel PSV columns in read-star — I/D events ≥ `--indel-psv-min-len` bp, clustered
    /// within 20 bp along the read, one column per cluster, the same per-column error as substitution
    /// columns. Default off (byte-identical).
    ///
    /// ⚠⚠ A7 (`docs/OPEN_ITEMS_2026-09-09.md`, register rows 793-794): REFUTED for the reason it was built,
    /// AND carries a known artifact class it does not exclude. Measured: 0/397 convertible human `tied`
    /// molecules convert, 0/23 gorilla MCL1, 0/5 MCL7 — reference copies this similar do not differ by an
    /// indel inside one read's span, so indels strengthen a margin but never break a genuine tie. Separately,
    /// columns within 20bp of a read END contradict the substitution-best copy 34% of the time (an aligner
    /// end-gap placement artifact, not real copy evidence) — this flag does not guard against that class.
    /// Kept as a flag only so the PREREG's negative result stays reproducible; do not enable it expecting it
    /// to resolve ties.
    #[arg(long, default_value_t = false)]
    indel_psv: bool,
    #[arg(long, default_value_t = 3)]
    indel_psv_min_len: u32,
    /// ⭐ §6fc: splice junctions as pairwise evidence in read-star (opt-in: +4 points of assignment at −0.3 points
    /// of MAPQ-60 agreement on the paired 35 families).
    #[arg(long, default_value_t = false)]
    read_star_junctions: bool,
    /// ⭐ §6fd: read-star against each candidate's GENOMIC locus (unit extent padded by the family's longest
    /// read), splice-aware; the origin certificate counts X + I + D + unaligned read bases over the read length.
    /// **Default ON (user decision 2026-09-05; PREREG adj/gstar held: 0 wrong anchors, 100 % placement
    /// agreement on the paired 35)**. `--read-star-unit` restores the spliced-unit form of §6fc (10× faster).
    #[arg(long, default_value_t = true)]
    read_star_genomic: bool,
    /// Escape hatch: read-star against the spliced UNITS (§6fc), the default before 2026-09-05 §6fd.
    #[arg(long, default_value_t = false)]
    read_star_unit: bool,
    /// Escape hatch (L2): ignore the catalog's `locus_start`/`locus_end` columns and pad each unit's extent by
    /// the family's longest molecule (the §6fd rule). Without the columns both arms are byte-identical.
    #[arg(long, default_value_t = false)]
    read_star_pad_locus: bool,
    /// Escape hatch (L3, §6fi): report every single-candidate molecule as `tied` (the §6fa rule) instead of
    /// `assigned` with `sole_candidate = 1` when its only candidate passes the origin certificate.
    #[arg(long, default_value_t = false)]
    no_sole_candidate: bool,
    /// ⭐ §6ha (user 2026-09-09; PREREG_best_by_psv_score_2026-09-09.md): escape hatch — pick read-star's BEST
    /// candidate (`bk`) by alignment (aligned bases − substitutions over the whole read), the pre-2026-09-09
    /// choice, instead of by PSV score (Σ ±1 over the columns the pairwise certificate itself uses). The
    /// alignment choice can prefer a copy that aligns more of the read at the cost of more substitutions over
    /// one that matches the read exactly wherever both are covered, ties it against the true best, and reports
    /// `Tied` — human MCL0: 12 of 343 in-family `tied` molecules had a unique PSV-perfect candidate whose LLR
    /// margin toward it was −48…−62. Required for `--no-as-tied-only` to reproduce the pre-2026-09-09 output.
    #[arg(long, default_value_t = false)]
    best_by_alignment: bool,
    /// ⭐ PREREG 819c1615: the read-star best = the candidate whose WORST pairwise duel (the certificate's own
    /// LLR over the columns both carry) is best — a candidate that beats every rival wins; twins fall through
    /// to the column-count score. Fixes the non-pairwise `psv_score` (§6he: 65 of 1,143 human contested rows
    /// had a best that loses a duel). Governs even with `--best-by-alignment`. ⭐ DEFAULT ON (user, 2026-09-09
    /// §6hk: +32 human assignments, 0 regressions, excision 21/23); `--no-best-by-duel` is the escape.
    #[arg(long, default_value_t = true)]
    best_by_duel: bool,
    /// Escape for the 2026-09-09 default: the read-star best by the column-count score again (pre-§6hk, byte-identical).
    #[arg(long, default_value_t = false)]
    no_best_by_duel: bool,
    /// L6: write `<out>.star_reads.tsv` — per molecule its read-star proof: the read positions of its columns,
    /// its base there and every candidate's base (the assignment-proof figure's source). Default OFF.
    #[arg(long, default_value_t = false)]
    dump_star: bool,
    /// Escape hatch (§6fm): count every read-star hit for a candidate even when it does not overlap the
    /// candidate's unit inside the target (the §6fd behaviour).
    #[arg(long, default_value_t = false)]
    no_read_star_hit_in_unit: bool,
    /// Escape hatch (§6fp): certify against the genomic locus alone (the §6fd form) instead of the better of
    /// locus and expressed chain per candidate.
    #[arg(long, default_value_t = false)]
    read_star_genomic_only: bool,
    /// Escape hatch (§6fq): do NOT assign uncontested molecules (primary MAPQ ≥ 60) to their placement; run the
    /// certificate machinery on every molecule as before.
    #[arg(long, default_value_t = false)]
    no_placement_assign: bool,
    /// §6fq variant: the placement OVERRIDES the machinery's certified call on uncontested molecules (the pure
    /// "any assembler would" rule). Default is certified-first with the placement as the fallback.
    #[arg(long, default_value_t = false)]
    placement_first: bool,
    /// Escape hatch (§6ft): no read-through certificate — every unaligned base counts against the candidate.
    #[arg(long, default_value_t = false)]
    no_readthrough_certificate: bool,
    /// ⭐ O2-8c (§6eo): discover PSV columns on the GENOMIC alignment of the copies' spans (exons + introns,
    /// reverse-complement retry for inverted duplications) instead of their spliced sequences. Read-chain units
    /// of unequal exon composition sent the spliced star projection to min_p 3e-270 on a wrong call (register
    /// 683); SEDEF core hulls are collinear genomic segments. Default OFF ⟹ byte-identical
    #[arg(long, default_value_t = false)]
    psv_genomic: bool,
    /// εⱼ used for an editing-flagged PSV column in the certificate (the rate a base shows the other allele
    /// by editing rather than sequencing error). Default 0.2.
    #[arg(long, default_value_t = 0.2)]
    edit_rate: f64,
    /// IsoCon-style iterative copy pruning: after assignment, repeatedly merge copies that have no read with
    /// significant evidence distinguishing them from their nearest neighbor, reassigning reads until all
    /// surviving copies are defensible. Default OFF (byte-identical baseline).
    #[arg(long, default_value_t = false)]
    iterative_prune: bool,
    /// RE-EXPRESS THE COPY ASSIGNMENT IN PHASING VOCABULARY (PS/HP tags + a GFA). This runs NO phasing
    /// algorithm of its own: it is a RELABELING of the assignment this binary already computed —
    /// phase set = family, haplotype = the assigned copy, unphased (-1) = abstained (ambiguous/tied).
    /// Measured: over 6 historical runs / 119,524 read-rows the multiset of
    /// `(read_name, family, haplotype)` equals `(read_name, family, assigned_copy-if-assigned-else -1)`
    /// EXACTLY, symmetric difference 0. (Compare per read NAME and you get a spurious mismatch — read
    /// names repeat within a family; the identity is only visible as a multiset.)
    /// Writes `<out>.phase_blocks.tsv` (one PHASE SET per family), `<out>.phased_haplotypes.tsv` (each
    /// copy's `pos:allele` PSV string), `<out>.phased_reads.tsv` (read → copy, with `n_psv_spanned` =
    /// `n_decisive` and `margin` = the assignment's log-LR margin), and `<out>.phase.gfa`.
    /// Block ≙ PS tag, haplotype ≙ HP tag, so downstream phasing tooling can consume it.
    /// ⚠ HISTORICAL NOTE: this help previously claimed a "min-path-cover over the PSV graph". No such
    /// computation exists in the O2 path — see `docs/copy_assignment_definition.md` §9.4.
    #[arg(long, default_value_t = false)]
    phase: bool,

    /// Optional gene annotation (GFF3/GTF/BED) to tag in-genome copies annotated vs unannotated in the
    /// --phase copy graph. Without it, in-genome copies are tagged `annotation-unknown`.
    #[arg(long)]
    gff: Option<String>,

    /// Skip the POA homology DIAGNOSTIC pass. It is the dominant per-region cost (the poasta all-pairs
    /// alignment over candidate rep pairs — ~85% of wall-clock on dense families) but is purely diagnostic:
    /// families come from the de-tie conflict graph, so the emitted families/assignments/abundance are
    /// BYTE-IDENTICAL with or without it. Only the `.fallback.tsv` report and the POA edge counts in the log
    /// are lost. STRONGLY recommended for genome-wide sweeps (measured ~6.8× faster on the heaviest family).
    /// Equivalent to setting `RUSTLE_SKIP_POA_DIAGNOSTIC=1`.
    #[arg(long, default_value_t = false)]
    skip_poa_diagnostic: bool,

    /// Discover reference-ABSENT (collapsed) copies from reads and re-thread the abstain pool against them
    /// (two-stage freeze; default OFF = byte-identical). Candidates failing the admission gate are written
    /// to `<out>.dna_needs.tsv`.
    #[arg(long, default_value_t = false)]
    absent_copies: bool,

    /// EXPERIMENT-ONLY override of the reference-absent admission gate's cluster floor
    /// (`absent_copy.rs` gate 1, `n_clusters >= min_clusters`). Unset = the shipping default 3, and
    /// output is byte-identical to a build without this flag.
    ///
    /// Gate 1 is an identifiability claim, NOT a tuning knob: below three co-varying clusters a
    /// second COPY is indistinguishable from a heterozygous ALLELE without DNA copy-number data.
    /// Lowering it in production manufactures false positives. The one legitimate use is a
    /// removal-recovery ablation (V4b), where a known copy is DELETED from the assembly, so copy
    /// status holds BY CONSTRUCTION and the gate is re-asking a question the design already
    /// answered — and where it is unreachable anyway (deleting one copy of a 3-copy family leaves
    /// at most 2 clusters). Such a run is interpretable ONLY beside the identical INTACT-assembly
    /// control at the same value: if the control also admits a copy, the recovery is an artefact.
    ///
    /// A value that is not a positive integer falls back to 3 (a typo must not disable the gate).
    #[arg(long, value_name = "N")]
    absent_min_clusters: Option<usize>,

    /// Opt-in augment-and-linearize REPORT (requires --absent-copies; no-op otherwise): compute a
    /// `LinearizeCertificate` for every Stage-2-admitted reference-absent candidate and write it to
    /// `<out>.linearize.tsv`. This costs one minimap2 realign-pool subprocess per admitted candidate, so it
    /// is OFF by default — plain `--absent-copies` keeps its prior admission cost and emits no certificate
    /// file (byte-identical to the pre-feature `--absent-copies` output). `--linearize-gate` implies this.
    #[arg(long, default_value_t = false)]
    linearize: bool,

    /// Opt-in augment-and-linearize GATE (implies --linearize; requires --absent-copies; no-op otherwise): a
    /// Stage-2 candidate whose `LinearizeCertificate` verdict is NOT `Linearizes` (its MAPQ-0 read pool does
    /// not land on it distinguishably more often than on a dinucleotide-shuffled decoy) is DEMOTED — written
    /// to `<out>.dna_needs.tsv` instead of admitted as a copy. Because it turns the certificate into an
    /// admission decision it also enables the certificate computation + `<out>.linearize.tsv` report (as if
    /// `--linearize` were set). Default OFF: admission is unchanged (byte-identical to plain --absent-copies).
    #[arg(long, default_value_t = false)]
    linearize_gate: bool,

    /// Emit `<out>.posterior.tsv`: per read, the soft per-copy POSTERIOR and the consistent ZONE (the genomic
    /// region of the copies it is compatible with) — the Bayesian complement to the hard assign/abstain, so an
    /// unassignable (Tied) read is localized to a zone with a distribution instead of a bare flag. The prior is
    /// uniform by default; set `RUSTLE_POSTERIOR_PRIOR=abundance` to weight by the EM copy abundance. Default off.
    #[arg(long, default_value_t = false)]
    posterior: bool,

    /// Run the EM soft-relaxation (Vollger 2019 PSV correlation-clustering, maximum-likelihood soft
    /// version) over each family's PSV evidence and emit `<out>.em.tsv` (per-read soft posterior +
    /// K-frontier label) and `<out>.em_abundance.tsv` (per-copy recovered abundance). This is a
    /// PSV-only reduction of the gate likelihood (editing-filtered; junctions and per-base quality
    /// are NOT used), so per-read labels may differ from the hard `.assignments.tsv` gate on
    /// junction/quality-resolvable reads; the abundance estimate is robust to these. `.em_abundance.tsv`
    /// is the convergent, gate-likelihood EM (the estimator the consistency theorem describes); it
    /// uses `error_rate` and a convergence gate, so its `pi_hat` can differ from the legacy
    /// `.quant.tsv` (`soft_quantify_em`, fixed error 0.01 / 100 iters) -- prefer `.em_abundance.tsv`
    /// for the theorem's estimator. Additive: leaves `.assignments.tsv`/`.families.tsv`/`.quant.tsv`
    /// byte-identical. Default off.
    #[arg(long, default_value_t = false)]
    em: bool,
    /// Max E/M sweeps for `--em`.
    #[arg(long, default_value_t = 500)]
    em_max_iter: usize,
    /// EM convergence tolerance (absolute+relative on the observed-data log-likelihood) for `--em`.
    #[arg(long, default_value_t = 1e-6)]
    em_eps: f64,

    /// Reference-free RNA single-copy expression floor (lambda_global): the genome-wide median
    /// n_reads over single-copy transcripts, precomputed by
    /// `bench/rna_copy_number_depth.py::global_single_copy_anchor` -- an RNA-only quantity, NOT
    /// genomic. When given, enables the `depth_cn` column of `<out>.famcn_readonly.tsv` (else
    /// `NA`); this file is ALWAYS written (additive; independent of every other output).
    #[arg(long)]
    lambda_global: Option<f64>,

    /// Read lambda_global from a `gw_family_catalog --single-copy-baseline` `<prefix>.lambda_global.tsv`
    /// instead of passing the scalar by hand. `--lambda-global <f64>` (if given) takes precedence. This makes
    /// the copy-number normalizer an in-binary quantity rather than an external-script number.
    #[arg(long)]
    lambda_file: Option<String>,

    /// VG re-align supplement (opt-in): for every co-located family, re-align each poor-fit/candidate
    /// read (low MAPQ, heavy clipping, or high divergence — `vg_realign::is_candidate`) to the family's
    /// copy-paths and record the decision (`reassigned` / `rejected` / `novel-candidate`) to
    /// `<out>.vg_realign.tsv`. FEEDS BACK into the assignment: BOTH the correction leg (re-thread + reassign
    /// among existing copies) AND the admission leg (admit novel-read pools as NEW copies — genome-touching,
    /// widens the roster). Default off; when off, every output is byte-identical. Use `--vg-realign-correct`
    /// for the correction leg WITHOUT the FP-risk admission.
    #[arg(long, default_value_t = false)]
    vg_realign: bool,

    /// VG re-align CORRECTION leg ONLY (opt-in): re-thread hard reads through the family copy-paths and
    /// correct their assignments among the EXISTING copies, WITHOUT admitting novel copies (no genome touch,
    /// no roster widening). The safe VG-native assignment leg. Off by default; byte-identical when off.
    #[arg(long, default_value_t = false)]
    vg_realign_correct: bool,

    /// Define family MEMBERSHIP by E_r transcript homology instead of the E_c read-conflict graph. The
    /// conflict graph links two copies only when reads map ambiguously between them, so a copy whose reads
    /// all map uniquely is dropped from its family and its reads come back `tied` — not because they are
    /// unassignable, but because their true copy was never admitted. Conflict, PSVs, and chi(H) remain
    /// within-family. Admitting a dropped copy enlarges the copy set, so the Bonferroni certificate
    /// alpha/(K-1) tightens and existing assignments shift; this is why the mode is opt-in. Requires
    /// minimap2 (honors RUSTLE_MINIMAP2); aborts rather than falling back to the conflict graph.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,

    /// DISABLE the mutual-homology + distinct-locus family gate. Refinement is ON BY DEFAULT: each co-located
    /// family must have its copies MUTUALLY HOMOLOGOUS (the shared E_r primary tier -- `-k 11 -w 5` at
    /// `sensitive_identity` 0.60 by default since X.4, `-x asm20` at 0.80 under `RUSTLE_ER_SENSITIVE_ONLY=0`
    /// -- cov-of-shorter>=0.50)
    /// across >= 2 distinct loci.
    /// ⚠ This is NOT "the same criterion `gw_family_catalog` refines by" and the two paths do NOT
    /// automatically agree — that claim was false from D1 (2026-08-09) onward and is corrected here
    /// (O-4). `gw_family_catalog`'s DEFAULT homology catalog does not call refine at all
    /// (`refine_enabled` = `refine_flag || !o1_homology`), so it applies no such gate; and where refine
    /// DOES run it additionally UNIONS a genomic-span tier in (see `additive_genomic_tier`), which the
    /// catalog's own E_r site never does. Both facts are certified per run: `additive_genomic_tier` in
    /// `<prefix>.rule.tsv` and `n_edges_genomic_tier_added` in `<prefix>.refine.params.tsv`.
    /// Without this gate the conflict oracle admits large-gene mis-chains (PBX1) and
    /// repeat-bridges as families (`bench/GW_CATALOG_FP_AUDIT.md`). `--no-refine` assigns the raw families and
    /// needs no minimap2.
    #[arg(long, default_value_t = false)]
    no_refine: bool,

    /// Keep unspliced readthrough transcripts as candidate copies. A single-exon de-novo transcript that
    /// engulfs >= 5 distinct splice junctions (each with >= 2 reads) is intronic pileup / unspliced pre-mRNA,
    /// not an mRNA, and is dropped by default. Validated on 15 such transcripts (minimum 14 engulfed
    /// junctions) against 260 expressed intronless genes (maximum 4), including the EEF1A1 retrocopy whose
    /// spliced parent cross-maps onto it. Pass this to disable the filter and reproduce the old behaviour.
    #[arg(long, default_value_t = false)]
    keep_readthrough: bool,

    /// EXPERIMENTAL, OFF BY DEFAULT. Admit a single-rep locus whose reads are ambiguously placed (MAPQ 0) at a
    /// rate incompatible with a unique locus as a multi-copy family with `n_copies = chi(H)`, reads certified
    /// tied, no per-copy sequence materialised.
    ///
    /// ⚠ The instrument detects unresolvable PARALOGY, not collapse. It fires on EEF1A1, whose MAPQ-0 reads
    /// align to processed pseudogenes on other chromosomes, and reports chi(H) = 7 for a one-copy locus. A copy
    /// genuinely absent from the reference would pile reads on at HIGH mapq, giving depth excess and no
    /// ambiguity -- which is why SDA detects collapses by depth. Do not use for copy number.
    #[arg(long, default_value_t = false)]
    collapse_gate: bool,

    /// NO-OP in copy_assign: this flag parses here but has no consumer in this binary. The
    /// collapse-enumeration gate (re-admitting near-identical <2-RNA-loci families as K0_COLLAPSED
    /// copy NUMBER, `<out>.collapsed.tsv`) runs only in `gw_family_catalog` -- use `--collapse-enumerate`
    /// there instead. Kept here only so `RUSTLE_COLLAPSE_ENUMERATE`/CLI parsing doesn't hard-error;
    /// NOT wired into copy_assign's per-read assignment path (that would violate the
    /// COPY-NUMBER-ONLY contract of this feature).
    #[arg(long, default_value_t = false)]
    collapse_enumerate: bool,

    /// Background per-read ambiguity rate for the collapse test (fraction of PRIMARY reads at MAPQ 0
    /// genome-wide). Must be a genome-wide quantity: a region-local estimate is degenerate, since in a
    /// collapsed window the reads outside the assembled rep are precisely the ambiguous ones. Default is the
    /// value measured on GGO_mm.bam (5785 / 4404440 = 0.001313). Recompute per sample with:
    ///   `echo $(( $(samtools view -c -F 2308 b.bam) - $(samtools view -c -F 2308 -q 1 b.bam) ))`
    #[arg(long)]
    eps_amb: Option<f64>,

    /// Apply the assemble gate's `min_reads` per ISOFORM (the pre-fix behaviour) instead of per LOCUS.
    /// Diagnostic: use to isolate the effect of junction-incidence pooling.
    #[arg(long, default_value_t = false)]
    no_pool_locus_support: bool,

    /// §6z6 (`docs/PREREG_primary_dedupe_2026-09-22.md`): keep PRIMARY reads that share
    /// (chrom, start, end, intron chain) with another read of the region instead of collapsing them to one.
    /// The historical key (default) was written to drop a record fetched twice from two overlapping copy
    /// windows, but it also collapses DISTINCT molecules with identical coordinates — 25-32% of primary
    /// records on A119b, ~23% of its 2-read chains pushed below pass-1's floor. With this flag a record is
    /// dropped only if it overlaps a window already fetched for the region (a true double-fetch).
    /// ⛔ Not the default: the shipped polish was fitted on de-duplicated counts and loses more matching
    /// chains than the raw assembly gains (human dev −147, gorilla held-out −1,096 at +3.2 pts precision).
    /// Byte-identical when unset.
    #[arg(long, default_value_t = false)]
    keep_coordinate_duplicates: bool,

    /// Seed candidate loci from AS-tied secondary reads that share an intron chain, even with no primary
    /// (recovers covered-but-tied K=0 copies as detected-but-unassignable). Implies fetching tied secondaries.
    #[arg(long, default_value_t = false)]
    tied_seed: bool,

    /// CONSUME the O1 catalog: a `gw_family_catalog` `<out>.copies.tsv`. Those rows ARE the copy set for the
    /// swept regions — this binary then detects NOTHING. Reads still come from `--bam`; only the copy set is
    /// supplied.
    ///
    /// # Why
    /// O1 and O2 share one node type, one edge engine and one admission primitive BY FUNCTION CALL and
    /// nothing BY FILE, so each binary re-derived its own families and the two tables had no join key
    /// (`GWFAM{i}` vs `CAFAM{i}`, assigned independently). Measured at defaults, they built different
    /// objects: the GSTM catalog has 4 copies where `copy_assign` found 0 families on 6031 reads. With
    /// `--families` the two copy sets agree BY CONSTRUCTION, and every emitted row carries the catalog's own
    /// `family_id` (so `family_id` here IS `GWFAM{i}`, not `CAFAM{i}`) plus a `<out>.family_join.tsv`
    /// mapping each assigned copy back to its `(family_id, copy_idx, tid)` catalog row.
    ///
    /// # What is switched OFF (documented, not silent)
    /// Family CONSTRUCTION, in full: pass-1 skeletons, the assemble gate, locus collapse, the membership
    /// oracle (E_c conflict graph / E_r homology), the POA diagnostic, co-location, the REFINE gate
    /// (`--no-refine` becomes moot), and the thin-locus RESCUE leg (it ADDS copies below the assembly read
    /// floor). Anything else that would change the roster is REFUSED rather than silently applied:
    /// `--absent-copies`, `--vg-realign` (admission leg), `--iterative-prune`, `--collapse-gate`,
    /// `--tied-seed`, `--recover-copies`. `--vg-realign-correct` is allowed (it only re-threads reads among
    /// the given copies). `--min-copies`/`--win` are not applied: the catalog already decided membership.
    ///
    /// # Contract (all loud, none silent)
    /// Every supplied copy must (a) be named by the `copies.tsv` header columns, (b) fall inside some
    /// `--region`/`--regions` entry on its OWN chromosome (c) have a sequence (see `--copies-fa`), and
    /// (d) have at least one overlapping read in the BAM. A violation aborts the run.
    ///
    /// A CROSS-CHROMOSOME family (RABL2's 5 contigs) is not truncated to whichever copies happen to fall
    /// in one region (2026-09-15): its reads are gathered directly from every one of its copies' own
    /// (chromosome, span) windows and pooled before assignment, via a synthetic `~xchrom~<family_id>`
    /// sweep key (see `catalog_input::group_families`/`load_supplied_families`). ⚠ KNOWN LIMITATION: the
    /// deeper PSV/mosaic certificate (`assign_family_detailed_once`, `best_overlap_copy`) still compares
    /// bare numeric positions with no chromosome field at all (`AlignedRead` carries none) — for a
    /// cross-chromosome family whose copies happen to sit at OVERLAPPING numeric coordinates on different
    /// chromosomes, that layer can attribute a read to the wrong copy. Safe whenever a family's per-
    /// chromosome coordinate ranges do not numerically coincide; not a general guarantee.
    #[arg(long)]
    families: Option<String>,

    /// The `gw_family_catalog` `<out>.copies.fa` beside `--families`. When given, each copy's spliced
    /// sequence is the catalog's OWN emitted bytes (checked against its `copies.tsv` row: chrom, span,
    /// strand, exon count) — so O2 assigns against exactly the sequence O1 defined the family with, with no
    /// reconstruction step that could differ. Without it the sequence is rebuilt from `--fasta` at the
    /// catalog's exon coordinates through the SAME `build_spliced_seq` the catalog used, and the strand it
    /// derives from the junction motifs must match the strand the catalog recorded (a mismatch means the
    /// FASTA is not the assembly the catalog was built against, and aborts). Prefer `--copies-fa`.
    #[arg(long)]
    copies_fa: Option<String>,
}

fn status_str(s: AssignStatus) -> &'static str {
    match s {
        AssignStatus::Assigned => "assigned",
        AssignStatus::Ambiguous => "ambiguous",
        AssignStatus::Tied => "tied",
    }
}

// ---- CROSS-FAMILY RECONCILIATION (RUSTLE_XFAM_RECONCILE) -------------------------------------------
//
// THE DEFECT. Every family is genotyped INDEPENDENTLY, and a molecule whose records fall in two
// families' read pools is put through both significance gates with no communication between them. On
// `mec` (12 supplied families, 2 regions, one NPIP-region BAM) that produces 79,175 assignment rows over
// 77,372 distinct molecules; 1,587 molecules are genotyped in >=2 families and 517 come back `assigned`
// in >=2 families at once, for 519 assigned-copy pairs. 210 of those pairs name copy intervals that are
// DISJOINT (median separation 16,341 bp, min 4,797 bp) — and a single molecule cannot originate from two
// disjoint loci, so at least one member of each such pair is wrong. Nothing downstream reconciled them.
// MEASURED split of the 519 (report arm, in-binary and authoritative): 309 `shared_locus`,
// 99 `readthrough_span`, 111 `cross_family_contradiction`; abstaining demotes 221/53,715 = 0.0041 of
// assigned rows, 110 molecules, and 221/1,035 = 0.2135 of CONTESTED rows — the third rate being the one
// that matters, since a rule that only fires on multi-placement molecules is selection-biased by
// construction. ⚠ ONE BAM, 2 regions, 6 family pairs, 170/210 disjoint pairs involving GWFAM111 alone
// and all 309 overlaps from 2 copy pairs: these COUNTS are not claimed to transport, only the mechanism.
//
// THE RULE, in three strata, and only the third one abstains:
//
//   `shared_locus`             the two claimed copy intervals OVERLAP. Two nominally different families
//                              claim one locus. That is an O1 PARTITION artifact (the same thing
//                              `catalog_overlaps` already warns about, below), not an O2 assignment
//                              error: the molecule really does come from that one interval, and both
//                              families name it. Reported, never demoted — demoting it would charge an
//                              O1 defect to O2's abstention rate and gut real copies of their support
//                              (MEASURED on `mec`, double-claimed/assigned: GWFAM113 copy0 141/141,
//                              GWFAM111 copy1 168/179, GWFAM112 copy12 141/153, GWFAM96 copy11
//                              168/234 — the same two flagged copy pairs `catalog_overlaps` reports).
//   `readthrough_span`         the two copies are disjoint but ONE alignment record is claimed by both —
//                              a molecule whose N-gap spans both loci. "A molecule cannot come from two
//                              disjoint loci" does not apply to a molecule that demonstrably spans both,
//                              so this stratum is reported, never demoted.
//   `cross_family_contradiction`  the two copies are disjoint AND the two claims rest on DIFFERENT
//                              alignment records. Two independent placements naming disjoint loci: the
//                              molecule has one origin, so at least one claim is false and there is no
//                              admissible way to tell which. This is exactly O2's assign-or-abstain
//                              contract at cross-family scope, so under `abstain` the molecule is
//                              demoted to `Ambiguous` in EVERY family it is assigned in.
//
// WHY NO ARBITRATION (a "keep the better one" arm was never built, deliberately). At cross-family scope
// none of the certificate fields is comparable (all MEASURED on `mec`, 519 contested pairs over 12
// families): the median assigned margin per family spans 9.9 to 10,745.2 (1,085x), `p_value` is 0.0 on
// BOTH sides of 434/519 = 0.8362 of contested pairs and is gated against a family-size-dependent
// `alpha/(n-1)`, and the margin winner is simply the larger-`n_decisive` side in 385/519 = 0.7418 (with
// 0 n_decisive ties). AS is byte-identical on both rows in 519/519, so an AS tie-break IS minimap2's
// primary flag in disguise — the defect that retired `uniq_agree` — and the primary flag itself agrees
// with the margin winner in only 61/109 = 0.5596 of the pairs where exactly one side is primary
// (chance, and underpowered). Abstention is
// therefore the whole of the intervention: it demotes ALL of a contradicting molecule's claims, because
// choosing which side to keep IS the arbitration being refused.
//
// WHAT THIS DOES NOT FIX (stated here and repeated above the `.quant.tsv` write). A demotion changes
// STATUS only. `n_reads_hard` counts `fa.assignments` by argmax `best_copy` with NO status filter, and
// `abundance`/`ci95` come from `soft_quantify_em` inside the per-family pipeline whose `obs_for_em` is
// populated regardless of status — so `.quant.tsv` is byte-identical between `report` and `abstain` and
// the row inflation there is untouched. Removing a molecule from a family's EM is a two-pass
// architecture change and a separate decision. Anyone claiming this fixes the abundance double-count is
// wrong.
//
// DEFAULT OFF. Under `off` (or unset) the detector does not run, nothing new is written beyond the
// unconditional `<out>.params.tsv`, `demote` is empty and every status emit site returns `status_str`
// verbatim — the OFF arm is byte-identical BY CONSTRUCTION, not by inspection.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum XfamMode {
    /// unset / `off`: the detector does not run and no side file is written (byte-identical arm).
    Off,
    /// `report`: write `<out>.xfam_conflicts.tsv` + a stderr summary; change no existing byte.
    Report,
    /// `abstain`: additionally demote every claim of a `cross_family_contradiction` molecule to `Ambiguous`.
    Abstain,
}

impl XfamMode {
    /// The EFFECTIVE value recorded in `<out>.params.tsv`. An unrecognized value is an error, not a
    /// silent `off`: a typo'd flag that quietly disables the pass is how an ON arm gets reported as OFF.
    fn from_env() -> Result<Self> {
        match std::env::var("RUSTLE_XFAM_RECONCILE").ok().as_deref() {
            None | Some("") | Some("off") => Ok(XfamMode::Off),
            Some("report") => Ok(XfamMode::Report),
            Some("abstain") => Ok(XfamMode::Abstain),
            Some(other) => anyhow::bail!(
                "RUSTLE_XFAM_RECONCILE={other:?} is not one of off|report|abstain"
            ),
        }
    }
    fn as_str(self) -> &'static str {
        match self {
            XfamMode::Off => "off",
            XfamMode::Report => "report",
            XfamMode::Abstain => "abstain",
        }
    }
}

/// A REGION-INDEPENDENT alignment-record key: `(contig, ref_start, ref_end, flags)` with
/// `flags` bit0 = secondary, bit1 = supplementary.
///
/// `fa.assignments` already carries the region-global `bam_reads` index (`idx_map[r.read_index]` in
/// `denovo_pipeline`), so the record a family's claim rests on is known exactly — but two REGIONS have
/// independent index spaces, so the index alone cannot be compared across them. This key can. Two records
/// at identical contig/start/end/flags are the same placement; collapsing them is correct and
/// conservative (it can only move a pair OUT of `cross_family_contradiction`).
type RecKey = (String, u64, u64, u8);

/// One family's `assigned` claim on one molecule: which (region, family) made it, which copy interval it
/// names, and which alignment record it rests on.
struct XfamClaim<'a> {
    g: usize,
    f: usize,
    fid: &'a str,
    copy: usize,
    span: (String, u64, u64),
    rec: RecKey,
}

/// One contested pair, as written to `<out>.xfam_conflicts.tsv`.
struct XfamConflict {
    read_name: String,
    stratum: &'static str,
    fid_a: String,
    copy_a: usize,
    span_a: (String, u64, u64),
    fid_b: String,
    copy_b: usize,
    span_b: (String, u64, u64),
    same_record: bool,
    sep_bp: u64,
    demoted: bool,
}

/// PASS 1: read-only sweep over every region's assignments, classifying every molecule claimed
/// `assigned` by two or more families.
///
/// Deterministic by construction: `BTreeMap`/`BTreeSet` keyed on `(read_name, region ordinal, family
/// ordinal)` with insertion-ordered inner vectors, and no `HashMap` iteration anywhere — so the output
/// does not depend on region-thread scheduling.
fn xfam_pass1(
    works: &[RegionWork],
    named_families: bool,
) -> (Vec<XfamConflict>, std::collections::BTreeSet<(String, usize, usize)>) {
    use std::collections::{BTreeMap, BTreeSet};
    // The family ids the DRAIN will mint, recomputed here with the same counter so the side file names
    // families exactly as `.assignments.tsv` does. Declared before `by_read` so it outlives the borrow.
    let mut gfam = 0usize;
    let mut fids: Vec<Vec<String>> = Vec::with_capacity(works.len());
    for work in works.iter() {
        let mut here = Vec::with_capacity(work.fams.len());
        for fa in &work.fams {
            here.push(if named_families { fa.family_id.clone() } else { format!("CAFAM{gfam}") });
            gfam += 1;
        }
        fids.push(here);
    }
    let mut by_read: BTreeMap<&str, Vec<XfamClaim>> = BTreeMap::new();
    for (g, work) in works.iter().enumerate() {
        for (f, fa) in work.fams.iter().enumerate() {
            for (ri, a) in &fa.assignments {
                if !matches!(a.status, AssignStatus::Assigned) {
                    continue;
                }
                let Some(span) = fa.copy_spans.get(a.best_copy).cloned() else { continue };
                let (rs, re, fl) = work.read_spans.get(*ri).copied().unwrap_or((0, 0, 0));
                by_read.entry(work.read_names[*ri].as_str()).or_default().push(XfamClaim {
                    g,
                    f,
                    fid: fids[g][f].as_str(),
                    copy: a.best_copy,
                    span,
                    rec: (work.read_chrom.get(*ri).cloned().unwrap_or_else(|| work.contig.clone()), rs, re, fl),
                });
            }
        }
    }
    let mut conflicts: Vec<XfamConflict> = Vec::new();
    let mut demote: BTreeSet<(String, usize, usize)> = BTreeSet::new();
    for (name, claims) in by_read.iter() {
        // ">= 2 DISTINCT (region, family)" is the contest condition. There is at most one row per
        // (read, family) today, but the guard is stated on the key, not on the count.
        let distinct: BTreeSet<(usize, usize)> = claims.iter().map(|c| (c.g, c.f)).collect();
        if distinct.len() < 2 {
            continue;
        }
        let mut contradicts = false;
        let mut here: Vec<XfamConflict> = Vec::new();
        for i in 0..claims.len() {
            for j in (i + 1)..claims.len() {
                let (a, b) = (&claims[i], &claims[j]);
                if (a.g, a.f) == (b.g, b.f) {
                    continue; // same family: the intra-family reduction already owns this case
                }
                let overlaps =
                    a.span.0 == b.span.0 && a.span.1 < b.span.2 && b.span.1 < a.span.2;
                let same_record = a.rec == b.rec;
                let stratum = if overlaps {
                    "shared_locus"
                } else if same_record {
                    "readthrough_span"
                } else {
                    contradicts = true;
                    "cross_family_contradiction"
                };
                let sep_bp = if overlaps || a.span.0 != b.span.0 {
                    0
                } else {
                    a.span.1.max(b.span.1) - a.span.2.min(b.span.2)
                };
                here.push(XfamConflict {
                    read_name: (*name).to_string(),
                    stratum,
                    fid_a: a.fid.to_string(),
                    copy_a: a.copy,
                    span_a: a.span.clone(),
                    fid_b: b.fid.to_string(),
                    copy_b: b.copy,
                    span_b: b.span.clone(),
                    same_record,
                    sep_bp,
                    demoted: false,
                });
            }
        }
        if contradicts {
            // CONTRADICTION IS A PROPERTY OF THE MOLECULE, so every claim it makes abstains — not only
            // the two in the contradicting pair. Abstaining one side would be arbitration, which is
            // exactly what is refused here.
            for c in claims {
                demote.insert(((*name).to_string(), c.g, c.f));
            }
            for c in here.iter_mut() {
                c.demoted = true;
            }
        }
        conflicts.extend(here);
    }
    (conflicts, demote)
}

/// GFA W-line SampleId sanitizer (a walk id must be a whitespace-free GFA token).
fn sanitize_gfa_id(s: &str) -> String {
    s.chars().map(|c| if c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-' | ':') { c } else { '_' }).collect()
}

/// Per-copy status shared by `build_copy_graph` (v1, PSV-bubble graph) and `build_exon_graph` (v2,
/// exon presence/absence graph) — the single source of truth for the in-genome/absent axis so the two
/// graph builders can never disagree about a copy's status.
///
/// A copy is ABSENT iff at least one read assigned to it is `discovery_coupled` — i.e. the copy exists only
/// because absent-copy discovery (`--absent-copies`) admitted it. This is the EXACT v1 rule, unchanged:
/// `copy_map_identity` feeds the `MI:f:` tag ONLY and must NEVER drive the absent/ST status (a copy admitted
/// via `--absent-copies` carries a remap identity but is not absent unless a `discovery_coupled` read pins
/// it). (The `collapsed_copies`/`rescued_copies` fields are diagnostic COUNTS, NOT per-`copy_tids` absence
/// markers — `collapsed_copies` can exceed `n_copies` — so they must NOT drive absence either.) An absent
/// copy is `AbsentCollapsed` when its genomic span OVERLAPS a non-absent (in-genome) copy of the same
/// family — a hidden CO-LOCATED haplotype — else `AbsentDivergent` (dispersed, or spans unavailable).
/// Non-absent in-genome copies are tagged by `annotation_status`: `InGenomeAnnotated`/`InGenomeUnannotated`
/// when `--gff` was given, else `AnnotationUnknown` (we never claim "unannotated" unchecked).
fn copy_status(
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ci: usize,
    ann: Option<&[(String, u64, u64)]>,
) -> rustle::vg_family::copy_graph::CopyStatus {
    use rustle::vg_family::copy_graph::CopyStatus;
    let n = fa.copy_tids.len();
    let is_coupled =
        |k: usize| -> bool { fa.assignments.iter().any(|(_, a)| a.discovery_coupled && a.best_copy == k) };
    if !is_coupled(ci) {
        return annotation_status(fa, ci, ann);
    }
    // Half-open genomic-span overlap of two copies (same chrom): distinguishes an absent copy that is a
    // hidden CO-LOCATED haplotype (AbsentCollapsed) from a dispersed one (AbsentDivergent).
    let span_overlap = |i: usize, j: usize| match (fa.copy_spans.get(i), fa.copy_spans.get(j)) {
        (Some((ci, si, ei)), Some((cj, sj, ej))) => ci == cj && si < ej && sj < ei,
        _ => false,
    };
    // AbsentCollapsed iff this absent copy's span overlaps some IN-GENOME (non-coupled) copy — a hidden
    // co-located haplotype; else it is dispersed / unlocalized => AbsentDivergent.
    let collapsed = (0..n).any(|k| k != ci && !is_coupled(k) && span_overlap(ci, k));
    if collapsed { CopyStatus::AbsentCollapsed } else { CopyStatus::AbsentDivergent }
}

/// Build one family's copy-graph (`--phase`): the REFERENCE walk + every copy as a tagged, corroborable
/// PATH + every read as a threaded WALK, over the family's usable PSV columns (those with both a genome
/// position and a fetchable reference base). Pure w.r.t. the caller except for `ref_base` (injected so this
/// is unit-testable without genome I/O — see `build_copy_graph_maps_family_to_graph`). Per-copy status is
/// `copy_status` (shared with `build_exon_graph` — see its doc for the absence rule).
/// `eff` is the EFFECTIVE per-row status, parallel to `fa.assignments` — `a.status` verbatim unless the
/// cross-family reconciliation demoted that row (see `XfamMode`). Passed in rather than read off
/// `a.status` so the graph's `Assigned` filters cannot disagree with `.assignments.tsv`.
fn build_copy_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ref_base: impl Fn(&str, u64) -> Option<u8>,
    bam_reads: &[String],
    ann: Option<&[(String, u64, u64)]>,
    eff: &[AssignStatus],
) -> rustle::vg_family::copy_graph::CopyGraph {
    use rustle::vg_family::copy_graph::*;
    // usable columns (both a genome position and a reference base), remembering the original index.
    let mut cols: Vec<PsvColumn> = Vec::new();
    let mut keep: Vec<usize> = Vec::new();
    for (j, p) in fa.psv_col_pos.iter().enumerate() {
        if let Some(pos) = p {
            if let Some(rb) = ref_base(&fa.chrom, *pos) {
                cols.push(PsvColumn { col: j, genome_pos: Some(*pos), ref_allele: Some(rb) });
                keep.push(j);
            }
        }
    }
    let sel = |row: &Vec<Option<u8>>| keep.iter().map(|&j| row.get(j).copied().flatten()).collect::<Vec<_>>();
    let backbone = vec![b"NNNNNNNNNN".to_vec(); cols.len() + 1];

    let n = fa.copy_tids.len();
    let copies = (0..n)
        .map(|ci| {
            let status = copy_status(fa, ci, ann);
            let reads = fa
                .assignments
                .iter()
                .zip(eff.iter())
                .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                .count() as u32;
            CopyPath {
                id: format!("{}_copy{}", fid, ci),
                alleles: fa.copy_psv_alleles.get(ci).map(|r| sel(r)).unwrap_or_default(),
                status,
                corrob: Corrob { reads: Some(reads), suns: None, map_identity: fa.copy_map_identity.get(ci).copied().flatten() },
            }
        })
        .collect();

    let reads = fa
        .assignments
        .iter()
        .zip(fa.read_psv_obs.iter())
        .enumerate()
        .map(|(k, ((ri, a), obs))| {
            let st = eff.get(k).copied().unwrap_or(a.status);
            ReadWalk {
                name: sanitize_gfa_id(&format!("{}_{}", fid, bam_reads[*ri])),
                obs: sel(obs),
                assigned_copy: if matches!(st, AssignStatus::Assigned) { Some(a.best_copy) } else { None },
                cert: Some(ReadCert { p_value: a.p_value, min_p_value: a.min_p_value, status: st }),
            }
        })
        .collect();

    CopyGraph { family: fid.to_string(), columns: cols, backbone, copies, reads }
}

/// Build one family's exon presence/absence graph (`--phase` v2): reconstructs each copy's genomic exon
/// chain from `copy_introns` + `copy_spans` (donor/acceptor intron chain -> exon intervals) and folds it
/// through `ExonGraph::from_copies` — the exon-level sibling of `build_copy_graph`'s PSV-bubble graph.
/// Per-copy status reuses `copy_status` (DRY with `build_copy_graph`), threading the same `--gff`
/// annotation overlay (`ann`) so in-genome copies come back `InGenomeAnnotated`/`InGenomeUnannotated`
/// under `--gff`, exactly as v1's `.phase.gfa`. The per-exon reference sequence is fetched later, by the
/// caller, via `ExonGraph::to_gfa`'s own `exon_seq` closure at write time (this builder lays out intervals
/// only, never sequence — no genome I/O here).
/// `eff` — see `build_copy_graph`.
fn build_exon_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ann: Option<&[(String, u64, u64)]>,
    eff: &[AssignStatus],
) -> rustle::vg_family::copy_graph::ExonGraph {
    use rustle::vg_family::copy_graph::*;
    let n = fa.copy_tids.len();
    let copies: Vec<(String, CopyStatus, Corrob, String, Vec<(u64, u64)>)> = (0..n)
        .map(|ci| {
            // Guard: `.get(ci)` never panics on a length-0/short fixture — missing span/introns => an
            // empty-exon copy (still walks zero classes; from_copies skips zero-length intervals anyway).
            let (chrom, start, end) =
                fa.copy_spans.get(ci).cloned().unwrap_or_else(|| (fa.chrom.clone(), 0, 0));
            let introns = fa.copy_introns.get(ci).cloned().unwrap_or_default();
            // genomic exons from the intron chain + outer span bounds:
            let mut exons = Vec::with_capacity(introns.len() + 1);
            let mut prev = start;
            for (d, a) in &introns {
                exons.push((prev, *d));
                prev = *a;
            }
            exons.push((prev, end));
            let reads = fa
                .assignments
                .iter()
                .zip(eff.iter())
                .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                .count() as u32;
            let status = copy_status(fa, ci, ann);
            let corrob = Corrob {
                reads: Some(reads),
                suns: None,
                map_identity: fa.copy_map_identity.get(ci).copied().flatten(),
            };
            (format!("{}_copy{}", fid, ci), status, corrob, chrom, exons)
        })
        .collect();
    ExonGraph::from_copies(fid, &copies)
}

/// In-genome annotation axis: overlap of copy `ci`'s span with any annotated interval.
/// `None` intervals => AnnotationUnknown (we never claim "unannotated" unchecked).
fn annotation_status(
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ci: usize,
    ann: Option<&[(String, u64, u64)]>,
) -> rustle::vg_family::copy_graph::CopyStatus {
    use rustle::vg_family::copy_graph::CopyStatus;
    let Some(ann) = ann else { return CopyStatus::AnnotationUnknown };
    let Some((c, s, e)) = fa.copy_spans.get(ci) else { return CopyStatus::AnnotationUnknown };
    let hit = ann.iter().any(|(ac, as_, ae)| ac == c && *as_ < *e && *s < *ae);
    if hit { CopyStatus::InGenomeAnnotated } else { CopyStatus::InGenomeUnannotated }
}

/// Parse a gene annotation file (`--gff`) into `(chrom, start0, end)` intervals: accepts BED (0-based, cols
/// 0/1/2) and GFF3/GTF (1-based, cols 0/3/4 -> start converted to 0-based). Lines starting with `#` and blank
/// lines are skipped; malformed lines (wrong column count or non-numeric coords) are skipped, never panic.
///
/// Format is decided from the FILE EXTENSION first — a BED6/BED12 line whose NAME column (col 3) is numeric
/// (`chr1\t1000\t2000\t5\t0\t+`) would otherwise be silently misread as GFF, yielding the wrong interval.
/// `.bed` => BED; `.gff`/`.gff2`/`.gff3`/`.gtf` => GFF/GTF. An unknown extension falls back to the per-line
/// numeric heuristic (GFF cols 3/4 preferred, else BED cols 1/2) for best-effort on mislabeled files.
fn parse_annotation(path: &str) -> anyhow::Result<Vec<(String, u64, u64)>> {
    enum Fmt { Bed, Gff, Auto }
    let lower = path.to_ascii_lowercase();
    let fmt = if lower.ends_with(".bed") {
        Fmt::Bed
    } else if lower.ends_with(".gff") || lower.ends_with(".gff2") || lower.ends_with(".gff3") || lower.ends_with(".gtf") {
        Fmt::Gff
    } else {
        Fmt::Auto
    };
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        let bed = |f: &[&str], out: &mut Vec<(String, u64, u64)>| {
            if f.len() < 3 { return; }
            let (Ok(s), Ok(e)) = (f[1].parse::<u64>(), f[2].parse::<u64>()) else { return }; // BED 0-based
            out.push((f[0].to_string(), s, e));
        };
        let gff = |f: &[&str], out: &mut Vec<(String, u64, u64)>| {
            if f.len() < 5 { return; }
            let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { return }; // GFF/GTF 1-based
            out.push((f[0].to_string(), s.saturating_sub(1), e));
        };
        match fmt {
            Fmt::Bed => bed(&f, &mut out),
            Fmt::Gff => gff(&f, &mut out),
            // unknown extension: prefer the GFF shape (5+ cols, numeric 3/4), else BED (3+ cols, numeric 1/2).
            Fmt::Auto => {
                let before = out.len();
                gff(&f, &mut out);
                if out.len() == before { bed(&f, &mut out); }
            }
        }
    }
    Ok(out)
}

fn parse_region(s: &str) -> Result<(String, u64, u64)> {
    let tok = s.split_whitespace().next().context("empty region")?;
    let (chrom, range) = tok.split_once(':').context("region must be chrom:start-end")?;
    let (lo_s, hi_s) = range.split_once('-').context("region must be chrom:start-end")?;
    Ok((chrom.to_string(), lo_s.parse().context("bad region start")?, hi_s.parse().context("bad region end")?))
}

/// Every `(contig, lo, hi)` from `--region`/`--regions` becomes its own independently-swept `RegionWork`
/// (own BAM query, own certificate). Two overlapping windows on the same contig would each independently
/// see and report the physical alignment records in the overlap — silently duplicating rows in
/// `--read-provenance` (and, unaudited, possibly other per-record outputs). Half-open `[lo, hi)`: touching
/// (one ends exactly where the next starts) is NOT an overlap.
fn validate_no_overlapping_regions(by_contig: &std::collections::BTreeMap<String, Vec<(u64, u64)>>) -> Result<()> {
    for (contig, windows) in by_contig {
        let mut sorted = windows.clone();
        sorted.sort_unstable();
        for w in sorted.windows(2) {
            let (a_lo, a_hi) = w[0];
            let (b_lo, b_hi) = w[1];
            if b_lo < a_hi {
                anyhow::bail!(
                    "overlapping --regions on {contig}: {contig}:{a_lo}-{a_hi} and {contig}:{b_lo}-{b_hi} \
                     would each independently see records in the overlap and duplicate them in \
                     --read-provenance (and possibly other per-record outputs); merge them into one region"
                );
            }
        }
    }
    Ok(())
}

/// A swept region, keyed exactly as the sweep iterates it.
type RegionKey = (String, u64, u64);
/// `--families`: the supplied catalog families BOUND to the region that will assign them.
type RegionFamilies = std::collections::BTreeMap<RegionKey, Vec<CatalogFamily>>;
/// Flank loaded around each supplied copy when gathering a dispersed family's reads (§6dh). Comfortably
/// above any long read, so a read reaching into a copy from outside its span is still collected.
const COPY_READ_PAD: u64 = 50_000;
/// Padded read-fetch windows for a swept region: `(chrom, lo, hi)` rather than bare `(lo, hi)` because a
/// CROSS-CHROMOSOME family's windows are not all on the region key's own contig — see `load_supplied_families`.
/// For every other region this is a redundant per-window copy of the key's own contig (harmless: `compute`
/// fetches each window from its own tagged chrom either way, and for those windows that is always `contig`).
type RegionWindows = std::collections::BTreeMap<RegionKey, Vec<(String, u64, u64)>>;
/// `--families`: catalog `tid` -> `(catalog family_id, catalog copy_idx)`. The JOIN KEY. Built from the
/// supplied table (never from the assignment output), so `<out>.family_join.tsv` reports the catalog's own
/// identity for a copy rather than an index this binary re-derived.
type CatalogIndex = std::collections::HashMap<String, (String, usize)>;

fn build_catalog_index(rf: &RegionFamilies) -> CatalogIndex {
    let mut ix = CatalogIndex::new();
    for fams in rf.values() {
        for f in fams {
            for c in &f.copies {
                ix.insert(c.tid.clone(), (c.family_id.clone(), c.copy_idx));
            }
        }
    }
    ix
}

/// Reference end (0-based, exclusive) of an aligned read. Local mirror of
/// `copy_assign_pipeline::read_ref_end` (which is `pub(crate)`), used only for the "every supplied copy has
/// reads" contract check below.
/// Aligned blocks of a read (0-based half-open): `M`/`=`/`X`/`D` extend, `N` closes.
fn aligned_blocks_local(read: &rustle::vg_family::copy_split::AlignedRead) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut p = read.ref_start;
    let mut cur: Option<(u64, u64)> = None;
    for &(op, n) in &read.cigar {
        match op {
            'M' | '=' | 'X' | 'D' => {
                cur = Some((cur.map_or(p, |c| c.0), p + n));
                p += n;
            }
            'N' => {
                if let Some(c) = cur.take() {
                    out.push(c);
                }
                p += n;
            }
            _ => {}
        }
    }
    if let Some(c) = cur {
        out.push(c);
    }
    out
}

fn read_ref_end_local(read: &rustle::vg_family::copy_split::AlignedRead) -> u64 {
    read.ref_start
        + read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X' | 'D' | 'N')).map(|(_, n)| n).sum::<u64>()
}

/// Reference-block overlap between an alignment and `[s, e)`, mirroring pysam's `get_blocks()` (used by
/// `bench/missing_copy_flag_pass.py`'s `truth` dict, lines 88-90: `sum(min(b1,e)-max(b0,s) for b0,b1 in
/// a.get_blocks() if b1>s and b0<e)`): only `M`/`=`/`X` runs count as aligned reference bases; `D`/`N`
/// advance the reference position without contributing overlap. A large intron (`N`) that merely SPANS
/// a target window contributes zero here, unlike `read_ref_end_local`'s span (`ref_start..ref_end`),
/// which would wrongly count the whole intron as covering it -- Task 7's reproduction-gate fix.
fn block_overlap(read: &rustle::vg_family::copy_split::AlignedRead, s: u64, e: u64) -> u64 {
    let mut pos = read.ref_start;
    let mut total = 0u64;
    for &(op, n) in &read.cigar {
        match op {
            'M' | '=' | 'X' => {
                let blk_end = pos + n;
                let lo = pos.max(s);
                let hi = blk_end.min(e);
                if hi > lo {
                    total += hi - lo;
                }
                pos = blk_end;
            }
            'D' | 'N' => pos += n,
            _ => {}
        }
    }
    total
}

/// The O3 detector's "truth" table: for every primary-aligned read, the family's own candidate copy
/// whose span the read's PRIMARY alignment overlaps MOST (by [`block_overlap`]), independent of O2's own
/// (possibly MAPQ-0-ambiguous) `catalog_copy_idx` call -- mirrors `bench/missing_copy_flag_pass.py`'s `truth` dict
/// (lines 84-90: `for a in BAM.fetch(...): ... if o > 0 and (name not in truth or o > truth[name][1]):
/// truth[name] = (i, o)`). Ties keep the FIRST-seen candidate (Python's strict `>` compare over
/// `cp.items()`'s insertion order, i.e. `copy_spans`' iteration order here). A candidate copy with no
/// entry in `catalog_index` (or no `catalog_index` at all) is silently skipped for that read.
///
/// Returns `(cf, cidx)` -- the TRUE catalog family id AND copy idx, not bare `cidx` alone (final
/// whole-branch-review fix round, follow-up): `copy_spans`/`copy_tids` here are `fa`'s ENTIRE co-located
/// copy set, which can span more than one true catalog family (Fix 1's own finding). Bare `cidx` is not
/// unique across that whole set -- two copies from DIFFERENT catalog families sharing a bare index would
/// make the caller's `*truth_cidx == cidx` comparison succeed for the wrong reason whenever a read's true
/// best-overlap copy is one of that colliding pair and O2's own `assignment.best_copy` names the other.
/// Confirmed live: Task 8's real data has 3 local families/arm with exactly this bare-cidx collision.
///
/// Extracted out of `main()`'s per-family loop so the tie-break is directly unit-testable.
fn best_overlap_truth_copy<'a>(
    bam_reads: &'a [BamRead],
    copy_spans: &[(String, u64, u64)],
    copy_tids: &[String],
    catalog_index: Option<&CatalogIndex>,
) -> std::collections::HashMap<&'a str, ((String, String), u64)> {
    let mut truth_copy: std::collections::HashMap<&str, ((String, String), u64)> = std::collections::HashMap::new();
    for br in bam_reads.iter().filter(|br| !br.is_secondary && !br.is_supplementary) {
        let end = read_ref_end_local(&br.read); // cheap span bound, to skip non-overlapping copies fast
        for (ci, (c, s, e)) in copy_spans.iter().enumerate() {
            if br.chrom != *c || end <= *s || br.read.ref_start >= *e {
                continue;
            }
            let ov = block_overlap(&br.read, *s, *e);
            if ov == 0 {
                continue;
            }
            let Some(tid) = copy_tids.get(ci) else { continue };
            let Some((cf, cidx)) = catalog_index.and_then(|ix| ix.get(tid)) else { continue };
            match truth_copy.get(br.name.as_str()) {
                Some((_, best_ov)) if *best_ov >= ov => {}
                _ => {
                    truth_copy.insert(br.name.as_str(), ((cf.clone(), cidx.to_string()), ov));
                }
            }
        }
    }
    truth_copy
}

/// `--discover-copies`: the candidate new copies ONE family's own AS-tied reads point at.
///
/// ⚠ Fix (final whole-branch review, Critical): the pre-fix code inlined this in the `fams.iter()`
/// closure and passed the WHOLE region's tied-read list to EVERY family, so one identical out-of-catalog
/// site -- identical read-name list and all -- was reported under 3-8 different `family_id` values in real
/// output. The design spec's own step 3 says clustering happens "across all reads IN A FAMILY". The read
/// set a family actually considered is exactly `fa.assignments`, whose `usize` is the region-global
/// `bam_reads` index (see the `RecKey` doc comment above), so this restricts `tied` to those names first.
///
/// ⚠ Second half of the same fix: `existing_copies` used to come from a `colocated.iter().find(|cf|
/// cf.family_id == fa.family_id)` join, which silently yielded an EMPTY exclusion list (`unwrap_or_default`)
/// whenever no `ColocatedFamily` carried that id -- and `fa` can bundle more than one true catalog family
/// (see the O3 block's own Fix 1 comment in `main`), so the join can and does fail on real input, making
/// every catalog copy invisible to `inside_any_copy`. `fa.copy_spans` / `fa.copy_tids` are parallel arrays
/// already assembled on `FamilyAssignment` -- the copies THIS family was assigned against, no join needed.
///
/// Residual, deliberately unchanged: the reported `family_id` is `fa.family_id`, the LOCAL co-located
/// group's id. When one `fa` bundles several true catalog families that label is one group, not one catalog
/// family -- but then `existing_copies` spans all of them too, so the exclusion stays conservative (a
/// candidate is only reported when it is outside EVERY copy the group holds).
fn discover_copies_for_family(
    fa: &FamilyAssignment,
    bam_reads: &[BamRead],
    tied: &[(String, Vec<rustle::vg_family::copy_discovery::TiePlacement>)],
) -> Vec<rustle::vg_family::copy_discovery::DiscoveredCopy> {
    let considered: std::collections::HashSet<&str> = fa
        .assignments
        .iter()
        .filter_map(|&(ri, _)| bam_reads.get(ri).map(|br| br.name.as_str()))
        .collect();
    let mine: Vec<(String, Vec<rustle::vg_family::copy_discovery::TiePlacement>)> =
        tied.iter().filter(|(name, _)| considered.contains(name.as_str())).cloned().collect();
    let existing_copies: Vec<(String, u64, u64, String)> = fa
        .copy_spans
        .iter()
        .zip(fa.copy_tids.iter())
        .map(|((chrom, start, end), tid)| (chrom.clone(), *start, *end, tid.clone()))
        .collect();
    rustle::vg_family::copy_discovery::cluster_tie_partners(
        &mine,
        &fa.family_id,
        &existing_copies,
        rustle::vg_family::copy_discovery::TIE_PARTNER_MERGE_DISTANCE_BP,
        rustle::vg_family::copy_discovery::TIE_PARTNER_MIN_SUPPORT,
    )
}

/// Load, VALIDATE and region-bind the `--families` catalog (see the flag's help for the contract).
///
/// Returns `(None, None, None)` when `--families` was not given — the historical path, untouched.
///
/// Everything here is a hard error. The one thing this function must never do is drop a supplied copy:
/// a copy silently missing from O2's roster is indistinguishable, downstream, from a copy O2 legitimately
/// found no evidence for, and that ambiguity is the exact defect class this contract exists to remove.
fn load_supplied_families(
    args: &Args,
    by_contig: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> Result<(Option<RegionFamilies>, Option<SeqIndex>, Option<RegionWindows>)> {
    let Some(path) = args.families.as_deref() else {
        if args.copies_fa.is_some() {
            anyhow::bail!("--copies-fa is only meaningful with --families (it supplies the copies' sequences)");
        }
        if args.flag_missing_copies {
            anyhow::bail!("--flag-missing-copies requires --families (it tests catalog copies for a missing sibling)");
        }
        if args.union_certificate {
            anyhow::bail!("--union-certificate requires --families (the union is over the supplied families' copies)");
        }
        return Ok((None, None, None));
    };
    // Roster-CHANGING legs are refused, not silently applied: with --families the copy set must be exactly
    // the supplied one, and each of these adds or removes copies. --vg-realign-correct is deliberately NOT
    // here: it re-threads reads among the GIVEN copies and never touches the roster.
    for (on, flag, why) in [
        (args.absent_copies, "--absent-copies", "admits reference-absent copies"),
        (args.vg_realign, "--vg-realign", "admits novel read pools as new copies"),
        (args.iterative_prune, "--iterative-prune", "merges/drops copies"),
        (args.collapse_gate, "--collapse-gate", "admits collapsed loci as extra copies"),
        (args.tied_seed, "--tied-seed", "seeds additional loci as copies"),
        (args.recover_copies, "--recover-copies", "feeds tied secondaries into copy rescue"),
    ] {
        if on {
            anyhow::bail!(
                "--families is incompatible with {flag}: it {why}, so the assigned copy set would no longer \
                 be the supplied catalog. Drop {flag}, or run without --families."
            );
        }
    }
    if args.no_refine {
        eprintln!(
            "[copy_assign] NOTE: --no-refine is moot under --families (no family construction runs at all)."
        );
    }
    if args.homology_primary {
        eprintln!(
            "[copy_assign] NOTE: --homology-primary is moot under --families (the membership oracle is the \
             supplied catalog, not E_c or E_r)."
        );
    }
    let text = std::fs::read_to_string(path).with_context(|| format!("reading --families {path}"))?;
    let fams = group_families(parse_copies_tsv(&text).with_context(|| format!("parsing --families {path}"))?)?;
    let seqs = match args.copies_fa.as_deref() {
        Some(p) => {
            let t = std::fs::read_to_string(p).with_context(|| format!("reading --copies-fa {p}"))?;
            Some(parse_copies_fa(&t).with_context(|| format!("parsing --copies-fa {p}"))?)
        }
        None => None,
    };
    // Cross-chromosome families (2026-09-15) never bind to a single swept region — see `cross_chrom` below
    // and its own containment check, per chromosome. Only single-chromosome families go through the
    // region-binding contract that follows, unchanged from every catalog built before this date.
    let (same_chrom, cross_chrom): (Vec<CatalogFamily>, Vec<CatalogFamily>) =
        fams.into_iter().partition(|f| !f.is_cross_chrom());
    // Which single real region (if exactly one) contains `[start, end)` on `chrom`. Shared by the
    // same-chromosome path (checked once for the whole family span) and the cross-chromosome path
    // (checked once per chromosome the family touches) — same containment contract either way.
    let contained_in = |chrom: &str, start: u64, end: u64| -> Vec<(u64, u64)> {
        by_contig
            .get(chrom)
            .map(|rs| rs.iter().copied().filter(|&(lo, hi)| start >= lo && end <= hi).collect())
            .unwrap_or_default()
    };
    // Bind each family to the ONE swept region that contains it. Containment (not overlap) is required:
    // a family straddling a region boundary would be assigned against the reads of only part of its own
    // span, which is the truncation this mode exists to prevent.
    let mut bound: RegionFamilies = RegionFamilies::new();
    // Per-family, per-chromosome clip bounds: the (lo, hi) of the one real region each chromosome's span
    // falls inside, keyed by the region key the family is BOUND to (same-chrom: its one real key;
    // cross-chrom: its one synthetic key) so the window builder below can look them up uniformly.
    let mut clip_bounds: std::collections::BTreeMap<RegionKey, std::collections::BTreeMap<String, (u64, u64)>> =
        std::collections::BTreeMap::new();
    for f in same_chrom {
        let hits = contained_in(&f.chrom, f.start, f.end);
        match hits.len() {
            0 => anyhow::bail!(
                "--families: {} ({}:{}-{}) lies outside every --region/--regions entry, so its reads would \
                 never be read. Add a region containing it, or remove it from the catalog.",
                f.family_id, f.chrom, f.start, f.end
            ),
            1 => {
                let key: RegionKey = (f.chrom.clone(), hits[0].0, hits[0].1);
                clip_bounds.entry(key.clone()).or_default().insert(f.chrom.clone(), hits[0]);
                bound.entry(key).or_default().push(f);
            }
            n => anyhow::bail!(
                "--families: {} ({}:{}-{}) is contained in {n} different swept regions, so which reads it \
                 would be assigned against is ambiguous. De-duplicate the region list.",
                f.family_id, f.chrom, f.start, f.end
            ),
        }
    }
    // CROSS-CHROMOSOME FAMILIES (2026-09-15): each of this family's per-chromosome spans (`chrom_spans`)
    // must individually sit inside some supplied region on ITS OWN chromosome — the same containment
    // contract as above, just checked once per chromosome instead of once for the whole (meaningless,
    // multi-chromosome) family span. There is no single real region to bind such a family to, so it is
    // bound instead to a SYNTHETIC region key (`~xchrom~<family_id>`, guaranteed not to collide with a
    // real contig name — no FASTA/BAM contig starts with `~`) that `main` adds to the swept list
    // alongside the real ones. `compute` fetches each of that key's windows from the window's OWN tagged
    // chromosome (see `RegionWindows`'s doc) rather than the key's, so the same assignment logic that
    // runs per real region pools this family's reads across every chromosome it touches and compares them
    // against its FULL copy set — never a truncated one.
    let mut n_cross_chrom = 0usize;
    for f in cross_chrom {
        let key: RegionKey = (format!("~xchrom~{}", f.family_id), 0, 0);
        let mut per_chrom_bounds: std::collections::BTreeMap<String, (u64, u64)> = std::collections::BTreeMap::new();
        for (chrom, (start, end)) in f.chrom_spans() {
            let hits = contained_in(&chrom, start, end);
            match hits.len() {
                0 => anyhow::bail!(
                    "--families: {} ({chrom}:{start}-{end}, one of its {} chromosomes) lies outside every \
                     --region/--regions entry, so its reads on {chrom} would never be read. Add a region \
                     containing it, or remove it from the catalog.",
                    f.family_id,
                    f.chrom_spans().len()
                ),
                1 => {
                    per_chrom_bounds.insert(chrom, hits[0]);
                }
                n => anyhow::bail!(
                    "--families: {} ({chrom}:{start}-{end}) is contained in {n} different swept regions on \
                     {chrom}, so which reads it would be assigned against is ambiguous. De-duplicate the \
                     region list.",
                    f.family_id
                ),
            }
        }
        eprintln!(
            "[copy_assign] --families: {} spans {} chromosomes ({}) — {} copies will be assigned together \
             via the cross-chromosome pass (key {:?}), not truncated to one region.",
            f.family_id,
            f.chrom_spans().len(),
            f.chrom_spans().keys().cloned().collect::<Vec<_>>().join(","),
            f.copies.len(),
            key.0
        );
        n_cross_chrom += 1;
        clip_bounds.insert(key.clone(), per_chrom_bounds);
        bound.entry(key).or_default().push(f);
    }
    // ⭐ DISPERSED-FAMILY READ WINDOWS (§6dh). A family binds to the ONE region containing its whole
    // span, but a genuinely DISPERSED family (NPIP: 38 copies over 89.5 Mb) makes that region enormous
    // and loading it whole costs 254,726 primaries to assign copies occupying a few hundred kb — it OOMs.
    // A read that overlaps NO copy can never be assigned to one, so the region's reads are gathered from
    // the union of the copies' own neighbourhoods instead. The anti-truncation guarantee is preserved:
    // every copy's reads are still loaded in full. Each window carries its OWN chrom (`c.chrom`, not the
    // region key's) so a cross-chrom family's windows on different chromosomes are never merged together,
    // and is clipped to the bounds of the one real region THAT chromosome's span was found inside (never
    // the key's own bounds, which for a cross-chrom family's synthetic key are meaningless placeholders).
    let mut windows: RegionWindows = std::collections::BTreeMap::new();
    for (k, fs) in &bound {
        let bounds_for = clip_bounds.get(k);
        let mut w: Vec<(String, u64, u64)> = fs
            .iter()
            .flat_map(|f| f.copies.iter())
            .filter_map(|c| {
                let (rlo, rhi) = *bounds_for?.get(&c.chrom)?;
                let (lo, hi) = (c.start.saturating_sub(COPY_READ_PAD).max(rlo), (c.end + COPY_READ_PAD).min(rhi));
                (lo < hi).then_some((c.chrom.clone(), lo, hi))
            })
            .collect();
        w.sort_unstable();
        let mut merged: Vec<(String, u64, u64)> = Vec::with_capacity(w.len());
        for (chrom, lo, hi) in w {
            match merged.last_mut() {
                Some(last) if last.0 == chrom && lo <= last.2 => last.2 = last.2.max(hi),
                _ => merged.push((chrom, lo, hi)),
            }
        }
        windows.insert(k.clone(), merged);
    }
    let n_fam: usize = bound.values().map(|v| v.len()).sum();
    let n_copy: usize = bound.values().flatten().map(|f| f.copies.len()).sum();
    eprintln!(
        "[copy_assign] --families {path}: {n_fam} famil{} / {n_copy} copies bound to {} region(s) ({} \
         cross-chromosome); sequences from {}",
        if n_fam == 1 { "y" } else { "ies" },
        bound.len(),
        n_cross_chrom,
        if seqs.is_some() { "--copies-fa (the catalog's own bytes)" } else { "--fasta (rebuilt at the catalog's exon coordinates)" },
    );
    Ok((Some(bound), seqs, Some(windows)))
}

/// Alignment-score evidence for every read in the region, indexed like the region's `bam_reads`.
///
/// Each `BamRead` is ONE alignment record, so a multimapper contributes several records under one name.
/// We group by name, so `best`/`second` are that READ's top two placements anywhere in the region — the
/// familiar "AS of the best hit vs the next best hit". Purely reported; `de` decides.
/// `exclude_supplementary` is true under the AS-tied gate (a supplementary is another segment of the read,
/// not an alternative placement) and false on the `--no-as-tied-only` escape, which must reproduce the
/// pre-2026-09-09 columns byte-for-byte.
fn as_evidence_per_read(bam_reads: &[BamRead], exclude_supplementary: bool) -> Vec<AsEvidence> {
    let aligned_len = |br: &BamRead| -> u32 {
        br.read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X')).map(|(_, n)| *n).sum::<u64>() as u32
    };
    // ⚠ A SUPPLEMENTARY record is another SEGMENT of the same read (a split/chimeric alignment), not an
    // alternative placement of it, so it can never be a tie partner: a primary + supplementary with equal
    // AS is one MAPQ-60 read in two pieces, and counting it as a tie let one such molecule through the
    // AS-tied gate and into `placement_assign` (§6gz addendum). Only primary + secondary records vote.
    let mut by_name: std::collections::HashMap<&str, Vec<(i32, u32)>> = std::collections::HashMap::new();
    for br in bam_reads.iter().filter(|br| !(exclude_supplementary && br.is_supplementary)) {
        by_name.entry(br.name.as_str()).or_default().push((br.as_score, aligned_len(br)));
    }
    bam_reads
        .iter()
        .map(|br| {
            // a molecule seen ONLY through supplementary records has no placement evidence of its own
            let fallback = vec![(br.as_score, aligned_len(br))];
            let placements = by_name.get(br.name.as_str()).unwrap_or(&fallback);
            as_evidence(placements).expect("read is its own placement, so the slice is non-empty")
        })
        .collect()
}

/// ⭐ O2 SCOPE (user, 2026-09-09). An **AS-TIED MULTIMAPPER**: ≥2 placements in the region whose runner-up
/// alignment score reaches `ratio` × the best. This is the population copy assignment exists for — the reads
/// where the aligner's primary/secondary label is a coin toss. A single-placement read, or one with a clear
/// best, had nothing for O2 to decide. ⚠ Region-local: placements on other contigs are not counted, so this
/// means "tied among the placements O2 was handed", not genome-wide multi-mapping.
static GATE_MOL_ALL: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static GATE_MOL_TIED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static GATE_REC_ALL: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static GATE_REC_TIED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static GATE_MOL_OUTSIDE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
/// `--gtf-copy-set`: per AS-tied molecule, the CATALOG copy indices at its tied placements and whether a tied
/// placement lies outside every unit (the copy SET an undecided isoform is emitted with).
static TIE_SET: std::sync::OnceLock<Mutex<std::collections::HashMap<String, (std::collections::BTreeSet<String>, bool)>>> = std::sync::OnceLock::new();
fn tie_set_of(name: &str) -> Option<(std::collections::BTreeSet<String>, bool)> {
    TIE_SET.get_or_init(Default::default).lock().unwrap().get(name).cloned()
}

/// One minimap2 fragment between two copy spans as blocks `(q_lo, q_hi, t_lo, sign)`: `q ∈ [q_lo, q_hi)` maps
/// to `t_lo + (q − q_lo)` (sign +) or `t_lo + (q_hi − 1 − q)` (sign −). The inverse of a block is
/// `(t_lo, t_lo + n, q_lo, sign)` under the same formula. Span-relative, 0-based (§6hn).
struct LiftBlocks {
    blocks: Vec<(u64, u64, u64, i8)>,
}
impl LiftBlocks {
    fn from_cigar(qs: u64, qe: u64, ts: u64, minus: bool, cg: &str) -> LiftBlocks {
        let mut blocks = Vec::new();
        let (mut q, mut t) = (if minus { qe } else { qs }, ts);
        let mut num: u64 = 0;
        for ch in cg.bytes() {
            if ch.is_ascii_digit() {
                num = num * 10 + (ch - b'0') as u64;
                continue;
            }
            match ch {
                b'=' | b'X' | b'M' => {
                    if minus {
                        blocks.push((q - num, q, t, -1i8));
                        q -= num;
                    } else {
                        blocks.push((q, q + num, t, 1i8));
                        q += num;
                    }
                    t += num;
                }
                b'I' => { if minus { q -= num } else { q += num } }
                b'D' | b'N' => t += num,
                _ => {}
            }
            num = 0;
        }
        blocks.sort_unstable();
        LiftBlocks { blocks }
    }
    fn inverse(&self) -> LiftBlocks {
        let mut blocks: Vec<(u64, u64, u64, i8)> = self.blocks.iter().map(|&(ql, qh, tl, s)| (tl, tl + (qh - ql), ql, s)).collect();
        blocks.sort_unstable();
        LiftBlocks { blocks }
    }
    /// `(mapped position, distance to the nearest aligned base)`; `None` when the fragment has no block.
    fn map(&self, q: u64) -> Option<(u64, u64)> {
        let at = |&(ql, qh, tl, s): &(u64, u64, u64, i8), x: u64| -> u64 { if s > 0 { tl + (x - ql) } else { tl + (qh - 1 - x) } };
        let i = self.blocks.partition_point(|b| b.0 <= q);
        if i > 0 {
            let b = self.blocks[i - 1];
            if q < b.1 {
                return Some((at(&b, q), 0));
            }
        }
        let mut best: Option<(u64, u64)> = None;
        for j in [i.wrapping_sub(1), i] {
            if let Some(b) = self.blocks.get(j) {
                for edge in [b.0, b.1 - 1] {
                    let d = edge.abs_diff(q);
                    if best.map_or(true, |(_, bd)| d < bd) {
                        let base = at(b, edge) as i64 + (q as i64 - edge as i64) * (b.3 as i64);
                        best = Some((base.max(0) as u64, d));
                    }
                }
            }
        }
        best
    }
}
/// All-vs-all copy-span lifts: `(a, b)` → the fragments mapping span `a` (relative) into span `b` (relative).
fn copy_span_lifts(spans: &[(String, u64, u64)], gi: &GenomeIndex, tag: &str) -> std::collections::HashMap<(usize, usize), Vec<LiftBlocks>> {
    let mut out: std::collections::HashMap<(usize, usize), Vec<LiftBlocks>> = std::collections::HashMap::new();
    let fa = format!("{tag}.copyset_spans.fa");
    {
        let Ok(mut fh) = std::fs::File::create(&fa) else { return out };
        for (i, (c, s0, e0)) in spans.iter().enumerate() {
            let seq = gi.fetch_sequence(c, *s0, *e0).unwrap_or_default();
            let _ = writeln!(fh, ">{i}");
            let _ = fh.write_all(&seq);
            let _ = fh.write_all(b"\n");
        }
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let res = std::process::Command::new(&mm2)
        .args(["-c", "-x", "asm20", "--eqx", "-X", "-t", "4"])
        .arg(&fa)
        .arg(&fa)
        .stderr(std::process::Stdio::null())
        .output();
    let _ = std::fs::remove_file(&fa);
    let Ok(res) = res else { return out };
    for line in String::from_utf8_lossy(&res.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (Ok(qi), Ok(qs), Ok(qe), Ok(ti), Ok(ts)) = (f[0].parse::<usize>(), f[2].parse::<u64>(), f[3].parse::<u64>(), f[5].parse::<usize>(), f[7].parse::<u64>()) else { continue };
        if qi == ti {
            continue;
        }
        let Some(cg) = f[12..].iter().copied().find_map(|t| t.strip_prefix("cg:Z:")) else { continue };
        let fwd = LiftBlocks::from_cigar(qs, qe, ts, f[4] == "-", cg);
        out.entry((ti, qi)).or_default().push(fwd.inverse());
        out.entry((qi, ti)).or_default().push(fwd);
    }
    out
}
static GATE_MOL_DISAGREE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
/// §6hd: molecules admitted to the gate by aligner self-disagreement rather than an AS tie.
static DISAGREE: std::sync::OnceLock<std::sync::Mutex<std::collections::HashSet<String>>> = std::sync::OnceLock::new();
fn register_disagreement(n: &str) { DISAGREE.get_or_init(Default::default).lock().unwrap().insert(n.to_string()); }
fn is_disagreement(n: &str) -> bool { DISAGREE.get().map_or(false, |m| m.lock().unwrap().contains(n)) }

fn as_tied(ev: &AsEvidence, ratio: f64) -> bool {
    match ev.second {
        Some(second) => (second as f64) >= ratio * (ev.best as f64),
        None => false,
    }
}

/// `NA` for an absent runner-up (single-placement read); otherwise the formatted value.
fn opt_i32(v: Option<i32>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| x.to_string())
}
fn opt_f32(v: Option<f32>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| format!("{x:.3}"))
}
/// `NA` for an absent distance (`--discover-copies`: a candidate on a chromosome this family has no copy
/// on), otherwise the value. Same convention as [`opt_i32`].
fn opt_u64(v: Option<u64>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| x.to_string())
}
/// `--read-provenance`: a record's own intron chain as `d1-a1,d2-a2,...`, or `none` for an unspliced record.
fn fmt_chain(chain: &[(u64, u64)]) -> String {
    if chain.is_empty() {
        "none".to_string()
    } else {
        chain.iter().map(|(d, a)| format!("{d}-{a}")).collect::<Vec<_>>().join(",")
    }
}

/// One assignment-table row (resolved while the region's reads are in scope).
struct AssignRow {
    read_name: String,
    family_id: String,
    assigned_copy: usize,
    status: &'static str,
    n_decisive: usize,
    margin: f64,
    p_value: f64,
    min_p_value: f64,
    /// Reported alignment-score evidence (see `as_evidence_per_read`). Never feeds the decision.
    as_ev: AsEvidence,
    junction_conflict: bool,
    origin_rejected: bool,
    n_candidates: usize,
    /// The read has an aligned BASE inside a copy of its family (§6es hygiene): rows with `false` are reads
    /// gathered from the copies' neighbourhoods that overlap no copy; report O2 on `in_copy == true`.
    in_copy: bool,
    /// ⭐ register 734: the molecule has a PRIMARY alignment overlapping a copy of this family. `in_copy`
    /// fires on any aligned block, so a genome-wide multimapper that only visits as a SECONDARY satisfies it
    /// and inflates every per-family rate — on DAZ, 16,257 of 18,192 rows were such visitors and the assigned
    /// fraction read 3.7 % instead of 34.7 %. **Every per-family read rate must be taken over this column.**
    primary_local: bool,
    /// §6fq: the molecule's primary MAPQ < 60 — the aligner could not place it; the certificate machinery
    /// decided this row. `false` = assigned to its placement (the certificate only reported).
    contested: bool,
    /// §6ft: the partner copy (catalog idx) whose alignment explained this molecule's unaligned tail, `cut` when
    /// only the mis-chain cut explained it, `-` otherwise.
    readthrough_into: String,
    /// The catalog `copy_idx` of `assigned_copy` under `--families` (copy_assign SORTS copies and reports its
    /// own index; `family_join.tsv` carries the same map). `NA` without a catalog.
    catalog_copy_idx: String,
    /// A3 (`docs/OPEN_ITEMS_2026-09-09.md`): whole-family PSV identity between `assigned_copy` and its
    /// nearest sibling (the competitor governing `p_value`). Emitted only with `--sibling-report`.
    sibling_identity: f64,
    /// A3: how many distinguishing PSV/junction positions this read spans against that nearest sibling.
    n_cols_vs_nearest_sibling: usize,
}
/// One family-table row.
struct FamilyRow {
    family_id: String,
    chrom: String,
    n_copies: usize,
    n_reads: usize,
    psv_cols: usize,
    resolvable_psv: usize,
    resolvable_j: usize,
    junction_only: usize,
    assigned_j: usize,
    uniq_agree: usize,
    uniq: usize,
    collapsed_copies: usize,
    rescued_copies: usize,
}
/// One per-copy soft-quantification row.
struct QuantRow {
    family_id: String,
    copy_index: usize,
    copy_tid: String,
    /// Genomic span of the copy. Emitted so a catalog can be audited for the same-locus artifact:
    /// two copies of one family whose spans overlap are one locus admitted twice, not two copies.
    copy_chrom: String,
    copy_start: u64,
    copy_end: u64,
    abundance: f64,
    ci: f64,
    n_hard: usize,
    /// L4: Σ posterior over molecules the certificate did not reject (a sole candidate counts 1, a K = 0 tie
    /// 1/k, an uncertified pair its softmax). Under the read-star this is the abundance's numerator.
    n_soft: f64,
    /// Tie-break invariance certificate: reads assigned to this copy that map UNIQUELY (`mapq > 0`), so their
    /// support survives any primary/secondary relabeling of tied reads. See `anchored_support`.
    anchored: usize,
    /// `anchored >= GATE_MIN_READS`: the copy exists under every tie-break (adversarially invariant) via unique
    /// mappers. FALSE = not guaranteed by unique mappers alone (may still be junction-defined, e.g. DAZ2).
    ///
    /// ⚠ A5 (`docs/OPEN_ITEMS_2026-09-09.md`, register row 786): NEAR-VACUOUS under the default AS-tied gate.
    /// `anchored` counts exactly the unique (MAPQ>0) mappers the gate removes BEFORE the certificate ever
    /// runs — measured on gorilla MCL1, 52/80 copies "true" by this column collapse to 2/80 once the gate is
    /// on, because the population this column certifies over has mostly left the certificate's business
    /// entirely. Kept (not removed) because `junction_invariant` — the OTHER half of the reported OR,
    /// `n_inv` below — is NOT vacuous the same way (copy-specific splice structure survives the gate); do
    /// not read a high `tie_invariant` count under `--no-as-tied-only` as evidence of anything under the
    /// default gate.
    tie_invariant: bool,
    /// `copy_junction_support >= GATE_MIN_READS`: the copy is pinned by >= 3 reads carrying a copy-specific
    /// JUNCTION (identifies it by splice structure regardless of the primary label) — the DAZ2-rescue mechanism.
    /// A copy is invariant overall if `tie_invariant || junction_invariant`.
    junction_invariant: bool,
}
/// One family-confirmed gene-conversion row.
struct MosaicRow {
    family_id: String,
    copy_a: usize,
    copy_b: usize,
    bp_lo: u64,
    bp_hi: u64,
    n_reads: usize,
    dispersion: u64,
    confirmed: bool,
}
/// One reference-free per-family copy-number row (Task R1: `chi_h` + `depth_cn`, no genome/assembly).
struct FamCnRow {
    family_id: String,
    chrom: String,
    n_copies: usize,
    n_reads: usize,
    chi_h: usize,
    depth_cn: f64, // NaN when --lambda-global was not given
    regime: &'static str,
}
/// One COPY-level historical gene-conversion row (a copy that is a mosaic of two others).
struct CopyConvRow {
    family_id: String,
    copy_c: String,
    copy_a: usize,
    copy_b: usize,
    bp_lo: u64,
    bp_hi: u64,
    n_decisive: usize,
}

/// Read the `lambda_global` scalar from a `lambda_global.tsv` (header + one data row, first column). `None` if
/// missing, unreadable, or the value is `NA`.
fn read_lambda_file(path: &str) -> Option<f64> {
    let text = std::fs::read_to_string(path).ok()?;
    let data = text.lines().nth(1)?; // skip header
    data.split('\t').next()?.trim().parse::<f64>().ok()
}

/// Resolve lambda: explicit scalar > file > none.
/// Anchored (tie-break-invariant) support for a copy: assigned reads that map UNIQUELY (`mapq > 0`) and whose
/// `best_copy == ci`. minimap2 sets `mapq = 0` exactly when the primary is tied/arbitrary, so a `mapq > 0` read
/// is never a candidate for primary/secondary relabeling — its support for `ci` is fixed under every tie-break.
/// A copy with `anchored >= GATE_MIN_READS` therefore exists under EVERY relabeling (adversarially invariant).
/// `best_copies` and `mapqs` are parallel over the family's assignments. Reproduces, in-binary, the experiment's
/// `samtools view -c -F 2308 -q 1 <copy_span>` bound.
fn anchored_support(best_copies: &[usize], mapqs: &[u8], ci: usize) -> usize {
    best_copies.iter().zip(mapqs.iter()).filter(|(&bc, &mq)| bc == ci && mq > 0).count()
}

fn resolve_lambda(explicit: Option<f64>, from_file: Option<f64>) -> Option<f64> {
    explicit.or(from_file)
}

fn verdict_str(v: rustle::vg_family::linearize::Verdict) -> &'static str {
    use rustle::vg_family::linearize::Verdict::*;
    match v {
        Linearizes => "LINEARIZES",
        Not => "NOT",
        Undetermined => "UNDETERMINED",
    }
}

/// One `<out>.linearize.tsv` row for a single Stage-2-admitted reference-absent candidate's
/// augment-and-linearize certificate (Task 4's `linearize_certs`). `NA` for NaN fracs/perm_p
/// (the `n_pool < min_pool` short-circuit in `linearize_certificate`).
fn linearize_tsv_row(fam: &str, loc: (&str, u64, u64), c: &LinearizeCertificate) -> String {
    let f = |x: f64| if x.is_nan() { "NA".to_string() } else { format!("{:.3}", x) };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        fam,
        loc.0,
        loc.1,
        loc.2,
        c.n_pool,
        f(c.linearized_frac_real),
        f(c.mean_frac_decoy),
        f(c.delta),
        if c.perm_p.is_nan() { "NA".to_string() } else { format!("{:.4}", c.perm_p) },
        verdict_str(c.verdict)
    )
}


/// §6p8 assembly polish: drop low-evidence transcripts from an emitted GTF using only the `reads "N"`
/// attribute. `mode` is "none" (no-op), "mono" (mono-exonic support floor) or "full" (floor + the
/// support-aware ISM collapse). Returns (ism_dropped, mono_dropped, floor) for the log line.
///
/// Validated in `bench/ASSEMBLY_POLISH.md` against `docs/PREREG_assembly_polish_2026-09-19.md`: the mono
/// floor costs zero matching intron chains on both chr20 and the held-out chr11; the ISM collapse trades
/// chains for precision.
fn polish_gtf_lines(
    lines: &mut Vec<String>,
    mode: &str,
    mono_quantile: f64,
    isoform_fraction: f64,
    mono_shadow: bool,
    ism_absolute_escape: bool,
    ism_3p_anchored: bool,
    ism_ratio: f64,
    fraction_exempt: bool,
    fuzzy: i64,
    fuzzy_ism: bool,
    fraction_min_reads: u64,
    retained_ratio: f64,
) -> (usize, usize, usize, u64, usize) {
    use std::collections::{HashMap, HashSet};
    if mode == "none" {
        return (0, 0, 0, 0, 0);
    }
    // exons per transcript, in genomic order, plus the transcript's read support
    let mut exons: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    let mut key: HashMap<String, (String, String)> = HashMap::new(); // tid -> (contig, strand)
    let mut reads: HashMap<String, u64> = HashMap::new();
    let mut gene: HashMap<String, String> = HashMap::new();
    for line in lines.iter() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let Some(tid) = re_attr(f[8], "transcript_id") else { continue };
        if let Some(r) = re_attr(f[8], "reads").and_then(|v| v.parse::<u64>().ok()) {
            let e = reads.entry(tid.clone()).or_insert(0);
            *e = (*e).max(r);
        }
        if let Some(g) = re_attr(f[8], "gene_id") {
            gene.entry(tid.clone()).or_insert(g);
        }
        if f[2] != "exon" {
            continue;
        }
        let (Ok(a), Ok(b)) = (f[3].parse::<i64>(), f[4].parse::<i64>()) else { continue };
        exons.entry(tid.clone()).or_default().push((a - 1, b));
        key.entry(tid).or_insert_with(|| (f[0].to_string(), f[6].to_string()));
    }
    let mut chain: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    let mut span: HashMap<String, (i64, i64)> = HashMap::new();
    for (tid, ex) in exons.iter_mut() {
        ex.sort_unstable();
        chain.insert(tid.clone(), (0..ex.len().saturating_sub(1)).map(|i| (ex[i].1, ex[i + 1].0)).collect());
        span.insert(tid.clone(), (ex[0].0, ex[ex.len() - 1].1));
    }

    // §6q1 the run's own "well-supported transcript" level: the mono_quantile-th percentile of multi-exon
    // read support. It is the mono-exonic floor AND the ISM pass's absolute escape, so a run with deep
    // coverage does not lose well-supported short isoforms merely because their containers are deeper.
    let mut multi_all: Vec<u64> = chain
        .iter()
        .filter(|(_, c)| !c.is_empty())
        .map(|(t, _)| reads.get(t).copied().unwrap_or(0))
        .collect();
    multi_all.sort_unstable();
    let support_level = if multi_all.is_empty() || mono_quantile <= 0.0 {
        0
    } else {
        multi_all[((mono_quantile * multi_all.len() as f64) as usize).min(multi_all.len() - 1)]
    };

    // §6q6 fuzzy junction comparison (isoseq's `--max-fuzzy-junction`): equality and sub-chain
    // containment up to a per-junction tolerance.
    let near = |a: (i64, i64), b: (i64, i64)| (a.0 - b.0).abs() <= fuzzy && (a.1 - b.1).abs() <= fuzzy;
    let chain_eq = |x: &[(i64, i64)], y: &[(i64, i64)]| -> bool {
        x.len() == y.len() && x.iter().zip(y.iter()).all(|(&a, &b)| near(a, b))
    };

    let mut drop: HashSet<String> = HashSet::new();

    // §6q6 pass 0: merge near-duplicate chains into their best-supported member. Runs before the ISM
    // collapse so a wobbled duplicate cannot act as a container, and before the mono floor so it cannot
    // shift the support quantile.
    if fuzzy > 0 {
        let mut buckets: HashMap<(&str, &str, usize), Vec<&String>> = HashMap::new();
        for (t, c) in chain.iter() {
            if c.is_empty() {
                continue;
            }
            if let Some(k) = key.get(t) {
                buckets.entry((k.0.as_str(), k.1.as_str(), c.len())).or_default().push(t);
            }
        }
        let mut keys: Vec<_> = buckets.keys().copied().collect();
        keys.sort();
        for bk in keys {
            let mut ts = buckets.remove(&bk).unwrap();
            // best-supported first, ties by id, so the kept representative is deterministic
            ts.sort_by(|a, b| {
                reads.get(*b).copied().unwrap_or(0).cmp(&reads.get(*a).copied().unwrap_or(0)).then_with(|| a.cmp(b))
            });
            for i in 0..ts.len() {
                if drop.contains(ts[i]) {
                    continue;
                }
                for j in (i + 1)..ts.len() {
                    if !drop.contains(ts[j]) && chain_eq(&chain[ts[i]], &chain[ts[j]]) {
                        drop.insert(ts[j].clone());
                    }
                }
            }
        }
    }
    let n_fuzzy = drop.len();

    if mode == "full" {
        // a fragment survives if it carries at least as much support as the chain that contains it, OR if
        // it is itself a well-supported transcript by this run's own standard
        let supported = |frag: &str, cont: &str| -> bool {
            let rf = reads.get(frag).copied().unwrap_or(0);
            if ism_absolute_escape && support_level > 0 && rf >= support_level {
                return true;
            }
            let rc = reads.get(cont).copied().unwrap_or(0);
            rc > 0 && (rf as f64) >= ism_ratio * rc as f64
        };
        let mut groups: HashMap<(String, String), Vec<String>> = HashMap::new();
        for (tid, k) in key.iter() {
            groups.entry(k.clone()).or_default().push(tid.clone());
        }
        for (_, tids) in groups.iter() {
            // deterministic: longest chain first, ties broken by transcript id. Both the container
            // scan and the mono-exonic host search depend on this order, so it must not come from a
            // HashMap's iteration order.
            let mut multi: Vec<&String> = tids.iter().filter(|t| !chain[*t].is_empty()).collect();
            multi.sort_by(|a, b| chain[*b].len().cmp(&chain[*a].len()).then_with(|| a.cmp(b)));
            // §6zb: the container scan is O(m²) over every multi-exon transcript of the contig+strand (chr1:
            // ~30k ⇒ ~10⁹ pair tests, 240 s of a 306 s run). A contiguous sub-chain's FIRST junction is one of
            // the container's junctions, so indexing candidates by first junction visits exactly the pairs
            // the full scan would drop, in the same `multi` order — byte-identical at tolerance 0. A fuzzy ISM
            // (tolerance > 0) cannot use exact-junction buckets and keeps the full scan.
            let exact_index = !(fuzzy_ism && fuzzy > 0);
            let mut by_first: HashMap<(i64, i64), Vec<&String>> = HashMap::new();
            if exact_index {
                for y in multi.iter() {
                    by_first.entry(chain[*y][0]).or_default().push(y);
                }
            }
            for x in multi.iter() {
                if drop.contains(*x) {
                    continue;
                }
                let cx = &chain[*x];
                let candidates: Vec<&String> = if exact_index {
                    cx.iter().flat_map(|j| by_first.get(j).into_iter().flatten().copied()).collect()
                } else {
                    multi.clone()
                };
                for y in candidates.iter() {
                    if x == y || drop.contains(*y) {
                        continue;
                    }
                    let cy = &chain[*y];
                    if cy.len() >= cx.len() {
                        continue;
                    }
                    // a contiguous sub-chain; under `ism_3p_anchored` only the one flush with the
                    // container's 3' end (suffix on '+', prefix on '-') counts as a truncation
                    let offsets: Vec<usize> = if ism_3p_anchored {
                        match key.get(*y).map(|k| k.1.as_str()) {
                            Some("-") => vec![0],
                            _ => vec![cx.len() - cy.len()],
                        }
                    } else {
                        (0..=cx.len() - cy.len()).collect()
                    };
                    let tol = if fuzzy_ism { fuzzy } else { 0 };
                    let sub = offsets.iter().any(|&k| {
                        cx[k..k + cy.len()]
                            .iter()
                            .zip(cy.iter())
                            .all(|(&a, &b)| (a.0 - b.0).abs() <= tol && (a.1 - b.1).abs() <= tol)
                    });
                    if sub && !supported(y, x) {
                        drop.insert((*y).clone());
                    }
                }
            }
            // §6zb: the mono-exonic host search was a linear scan of `multi` per single-exon transcript (O(mono ×
            // multi), the other half of chr1's 240 s). Host = the FIRST surviving multi in `multi` order whose span
            // contains the mono's span. Offline: monos by end descending, multis inserted by end descending into
            // a prefix-min Fenwick over their (sorted) starts holding their `multi` rank; the query is the min
            // rank among starts ≤ mono.start. Same host as the scan, so byte-identical.
            let alive: Vec<(usize, (i64, i64))> =
                multi.iter().enumerate().filter(|(_, m)| !drop.contains(**m)).map(|(r, m)| (r, span[*m])).collect();
            let mut starts: Vec<i64> = alive.iter().map(|(_, s)| s.0).collect();
            starts.sort_unstable();
            starts.dedup();
            let nfen = starts.len();
            let mut fen: Vec<usize> = vec![usize::MAX; nfen + 1];
            let mut by_end: Vec<(i64, usize, i64)> = alive.iter().map(|(r, s)| (s.1, *r, s.0)).collect();
            by_end.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
            let mut tids: Vec<&String> = tids.iter().collect();
            tids.sort();
            let mut monos: Vec<(&String, (i64, i64))> =
                tids.iter().filter(|t| chain[**t].is_empty() && !drop.contains(**t)).map(|t| (*t, span[*t])).collect();
            monos.sort_by(|a, b| b.1 .1.cmp(&a.1 .1).then_with(|| a.0.cmp(b.0)));
            let mut p = 0usize;
            for (t, s) in monos {
                while p < by_end.len() && by_end[p].0 >= s.1 {
                    let mut i = starts.partition_point(|&x| x < by_end[p].2) + 1;
                    while i <= nfen {
                        fen[i] = fen[i].min(by_end[p].1);
                        i += i & i.wrapping_neg();
                    }
                    p += 1;
                }
                let mut i = starts.partition_point(|&x| x <= s.0);
                let mut best = usize::MAX;
                while i > 0 {
                    best = best.min(fen[i]);
                    i -= i & i.wrapping_neg();
                }
                if best != usize::MAX {
                    let h = multi[best];
                    if !supported(t, h) {
                        drop.insert(t.clone());
                    }
                }
            }
        }
    }
    let n_ism = drop.len() - n_fuzzy;

    // mono-exonic support floor, read off this run's own multi-exon distribution
    let mut multi_reads: Vec<u64> = chain
        .iter()
        .filter(|(t, c)| !c.is_empty() && !drop.contains(*t))
        .map(|(t, _)| reads.get(t).copied().unwrap_or(0))
        .collect();
    multi_reads.sort_unstable();
    // computed over the SURVIVORS of the ISM pass, which is why it can differ from `support_level`
    let floor = if multi_reads.is_empty() || mono_quantile <= 0.0 {
        0
    } else {
        multi_reads[((mono_quantile * multi_reads.len() as f64) as usize).min(multi_reads.len() - 1)]
    };
    if floor > 0 {
        for (t, c) in chain.iter() {
            if c.is_empty() && !drop.contains(t) && reads.get(t).copied().unwrap_or(0) < floor {
                drop.insert(t.clone());
            }
        }
    }
    // §6q0 same-strand shadow: a single-exon transcript on a spliced gene's own footprint is that gene's
    // unspliced signal, not a gene. Anti-strand overlap is deliberately NOT a criterion.
    if mono_shadow {
        // multi-exon EXONS keyed by contig alone (either strand), and multi-exon SPANS keyed by
        // contig+strand (same strand only)
        let mut exons_any: HashMap<&str, Vec<(i64, i64)>> = HashMap::new();
        let mut spans_same: HashMap<(&str, &str), Vec<(i64, i64)>> = HashMap::new();
        for (t, c) in chain.iter() {
            if c.is_empty() || drop.contains(t) {
                continue;
            }
            let Some(k) = key.get(t) else { continue };
            spans_same.entry((k.0.as_str(), k.1.as_str())).or_default().push(span[t]);
            if let Some(ex) = exons.get(t) {
                exons_any.entry(k.0.as_str()).or_default().extend(ex.iter().copied());
            }
        }
        // merge each interval list so the overlap probe is a single sorted scan
        let merge = |v: &mut Vec<(i64, i64)>| {
            v.sort_unstable();
            let mut out: Vec<(i64, i64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match out.last_mut() {
                    Some(last) if a <= last.1 => last.1 = last.1.max(b),
                    _ => out.push((a, b)),
                }
            }
            *v = out;
        };
        for v in exons_any.values_mut() {
            merge(v);
        }
        for v in spans_same.values_mut() {
            merge(v);
        }
        // merged and disjoint: the only candidate is the last interval starting at or before `e`
        let hits = |v: &Vec<(i64, i64)>, s: i64, e: i64| -> bool {
            let hi = v.partition_point(|&(a, _)| a <= e);
            hi > 0 && v[hi - 1].1 >= s
        };
        let mut mono: Vec<&String> =
            chain.iter().filter(|(t, c)| c.is_empty() && !drop.contains(*t)).map(|(t, _)| t).collect();
        mono.sort();
        for t in mono {
            let Some(k) = key.get(t) else { continue };
            let (s, e) = span[t];
            let exon_hit = exons_any.get(k.0.as_str()).is_some_and(|v| hits(v, s, e));
            let span_hit = spans_same.get(&(k.0.as_str(), k.1.as_str())).is_some_and(|v| hits(v, s, e));
            if exon_hit || span_hit {
                drop.insert(t.clone());
            }
        }
    }
    let n_mono = drop.len() - n_ism - n_fuzzy;

    // §6p9 locus isoform fraction: a transcript far below the best-supported isoform of its own locus is
    // a minor-flow artifact. The locus dominant is never dropped, so no locus is ever emptied.
    if isoform_fraction > 0.0 {
        let mut best: HashMap<&str, u64> = HashMap::new();
        for (t, g) in gene.iter() {
            if drop.contains(t) {
                continue;
            }
            let r = reads.get(t).copied().unwrap_or(0);
            let e = best.entry(g.as_str()).or_insert(0);
            *e = (*e).max(r);
        }
        let mut candidates: Vec<&String> = gene.keys().filter(|t| !drop.contains(*t)).collect();
        candidates.sort();
        for t in candidates {
            let Some(g) = gene.get(t) else { continue };
            let b = best.get(g.as_str()).copied().unwrap_or(0);
            if b == 0 {
                continue;
            }
            let r = reads.get(t).copied().unwrap_or(0);
            if fraction_exempt && support_level > 0 && r >= support_level {
                continue;
            }
            if fraction_min_reads > 0 && r >= fraction_min_reads {
                continue;
            }
            if r < b && (r as f64) < isoform_fraction * b as f64 {
                drop.insert(t.clone());
            }
        }
    }
    let n_frac = drop.len() - n_ism - n_mono - n_fuzzy;

    // §6za retained-intron filter, over the survivors of every step above, in one pass (support and the
    // candidate junction sets are fixed before any drop, so the result does not depend on order).
    let mut n_ret = 0usize;
    if retained_ratio > 0.0 {
        let mut support: HashMap<(&str, &str, (i64, i64)), u64> = HashMap::new();
        let mut gj: HashMap<(&str, &str), Vec<(i64, i64)>> = HashMap::new();
        for (t, c) in chain.iter() {
            if drop.contains(t) {
                continue;
            }
            let Some(k) = key.get(t) else { continue };
            let r = reads.get(t).copied().unwrap_or(0);
            for &j in c.iter() {
                *support.entry((k.0.as_str(), k.1.as_str(), j)).or_insert(0) += r;
            }
            if let Some(g) = gene.get(t) {
                gj.entry((g.as_str(), k.1.as_str())).or_default().extend(c.iter().copied());
            }
        }
        for v in gj.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        let mut cands: Vec<&String> = chain.keys().filter(|t| !drop.contains(*t)).collect();
        cands.sort();
        let mut newly: HashSet<String> = HashSet::new();
        for t in cands {
            let (Some(g), Some(k), Some(ex)) = (gene.get(t), key.get(t), exons.get(t)) else { continue };
            let Some(js) = gj.get(&(g.as_str(), k.1.as_str())) else { continue };
            let own: HashSet<(i64, i64)> = chain[t].iter().copied().collect();
            let rt = reads.get(t).copied().unwrap_or(0);
            for &j in js.iter() {
                if own.contains(&j) {
                    continue;
                }
                let inside = ex.iter().any(|&(a, b)| a < j.0 && j.1 < b);
                if inside && (support.get(&(k.0.as_str(), k.1.as_str(), j)).copied().unwrap_or(0) as f64) >= retained_ratio * rt as f64 {
                    newly.insert(t.clone());
                    break;
                }
            }
        }
        n_ret = newly.len();
        drop.extend(newly);
    }

    if !drop.is_empty() {
        lines.retain(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 9 {
                return true;
            }
            match re_attr(f[8], "transcript_id") {
                Some(t) => !drop.contains(&t),
                None => true,
            }
        });
    }
    (n_ism + n_fuzzy, n_mono, n_frac, floor, n_ret)
}


/// §6r5: add `cov` and `TPM` to every transcript line of an emitted GTF, from its `reads` attribute.
///
/// Count-based by design — `TPM_i = reads_i / sum_j(reads_j) * 1e6` — because a long read is one molecule,
/// so dividing by transcript length would down-weight long transcripts that were sequenced end to end.
/// Measured against StringTie's TPM on chr20 (507 shared intron chains): count-based rho 0.879,
/// length-normalised 0.714. Returns the number of transcript lines annotated.
fn annotate_tpm(lines: &mut [String]) -> usize {
    let reads_of = |line: &str| -> Option<f64> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "transcript" {
            return None;
        }
        re_attr(f[8], "reads").and_then(|v| v.parse::<f64>().ok())
    };
    let total: f64 = lines.iter().filter_map(|l| reads_of(l)).sum();
    if total <= 0.0 {
        return 0;
    }
    // spliced length per transcript id, from its exon lines
    let mut len: std::collections::HashMap<String, i64> = std::collections::HashMap::new();
    for line in lines.iter() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let (Some(t), Ok(a), Ok(b)) = (re_attr(f[8], "transcript_id"), f[3].parse::<i64>(), f[4].parse::<i64>())
        else {
            continue;
        };
        *len.entry(t).or_insert(0) += b - a + 1;
    }
    let mut n = 0;
    for line in lines.iter_mut() {
        let Some(r) = reads_of(line) else { continue };
        let tid = {
            let f: Vec<&str> = line.split('\t').collect();
            re_attr(f[8], "transcript_id")
        };
        let l = tid.and_then(|t| len.get(&t).copied()).unwrap_or(0).max(1) as f64;
        line.push_str(&format!(" cov \"{:.6}\"; TPM \"{:.6}\";", r / (l / 1000.0), r / total * 1e6));
        n += 1;
    }
    n
}

fn main() -> Result<()> {
    let mut args = Args::parse();
    // §6p6: --assemble-only IS the assembly product, so it implies --gtf. Setting it here means every
    // existing `if args.gtf` gate fires unchanged rather than each one needing a second condition.
    if args.assemble_only {
        args.gtf = true;
        // §6r8: assembly never reads a read's bases or qualities (only its CIGAR/introns), and O2 — the
        // only consumer — is skipped in this mode. Dropping them at parse time is where the memory is:
        // peak RSS was 8.4-11.8 GB for one chromosome, which capped concurrency at 2 and OOM-killed 4.
        rustle::vg_family::denovo_assemble::SKIP_READ_SEQUENCE
            .store(true, std::sync::atomic::Ordering::Relaxed);
        eprintln!(
            "[copy_assign] ASSEMBLE-ONLY: family detection, homology refinement and copy assignment are \
             SKIPPED; running loci + isoform assembly only. `<out>.families.tsv`/`.assignments.tsv` will be \
             empty by construction."
        );
        // §6za: the transcript product uses strict canonical junctions by default; the family paths are
        // untouched because this runs only under --assemble-only, and an explicit env value wins.
        match args.assembly_junctions.as_str() {
            "strict" => {
                if std::env::var_os("RUSTLE_JUNCTION_MAJORITY").is_none() {
                    std::env::set_var("RUSTLE_JUNCTION_MAJORITY", "0");
                }
            }
            "majority" => {}
            other => anyhow::bail!("--assembly-junctions must be `strict` or `majority`, got `{other}`"),
        }
        eprintln!(
            "[copy_assign] ASSEMBLE-ONLY junctions: {} (RUSTLE_JUNCTION_MAJORITY={})",
            args.assembly_junctions,
            std::env::var("RUSTLE_JUNCTION_MAJORITY").unwrap_or_else(|_| "unset (majority)".into())
        );
    }
    // §6eu: the pipeline reads RUSTLE_PSV_READFILTER; an explicit env value wins, else the flag decides.
    if std::env::var_os("RUSTLE_PSV_READFILTER").is_none() {
        std::env::set_var("RUSTLE_PSV_READFILTER", if args.psv_read_filter { "1" } else { "0" });
    }
    if args.igv {
        args.dump_psv = true; // --igv is a bundle: the PSV matrix feeds bench/igv_tracks.py -> tagged BAM + PSV VCF
    }
    // --gff: parsed ONCE before the sweep into the annotation axis intervals (None => every in-genome copy in
    // the --phase copy graph stays AnnotationUnknown, byte-identical to the no-flag path).
    let annotation: Option<Vec<(String, u64, u64)>> =
        args.gff.as_deref().map(parse_annotation).transpose().context("parsing --gff")?;
    // Minor (final whole-branch review): `bench/missing_copy_flag_pass.py`'s reference requires `--gff` to reach its
    // `annotated_no_unit` orphan-locus class at all (no gene intervals -> every locus is either
    // `other_family` or `unannotated`) -- without it, `--flag-missing-copies` silently never produces an
    // `annotated_no_unit` row rather than erroring, which is easy to miss on a real run.
    if args.flag_missing_copies && args.gff.is_none() {
        eprintln!("[copy_assign] WARNING: --flag-missing-copies without --gff can never classify an orphan locus as annotated_no_unit (no gene intervals to test against)");
    }

    // Augment-and-linearize is opt-in: the report (`--linearize`) or the gate (`--linearize-gate`, which
    // implies the report) turns it on. When false, the certificate is skipped entirely inside
    // `detect_and_assign` (no minimap2 realign-pool per candidate) and no `<out>.linearize.tsv` is written,
    // so plain `--absent-copies` is byte-identical to the pre-feature `--absent-copies` path.
    let do_linearize = args.linearize || args.linearize_gate;

    // collect the regions, then group by contig so each contig loads ONCE (memory-bounded sweep).
    let regions: Vec<(String, u64, u64)> = match (&args.region, &args.regions, args.genome_wide) {
        (Some(r), None, false) => vec![parse_region(r)?],
        (None, Some(f), false) => std::fs::read_to_string(f)
            .with_context(|| format!("reading {f}"))?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(parse_region)
            .collect::<Result<_>>()?,
        (None, None, true) => {
            // §6zb: every reference sequence with >= 1 mapped record in the .bai, in header order
            use noodles_csi::binning_index::ReferenceSequence as _;
            let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(&args.bam)?;
            let header = reader.read_header()?;
            let index = noodles_bam::bai::read(format!("{}.bai", args.bam)).context("--genome-wide needs a .bai")?;
            let mut v = Vec::new();
            for (i, (name, rs)) in header.reference_sequences().iter().enumerate() {
                let mapped = index.reference_sequences().get(i).and_then(|r| r.metadata()).map(|m| m.mapped_record_count()).unwrap_or(0);
                if mapped > 0 {
                    v.push((String::from_utf8_lossy(name).to_string(), 0, usize::from(rs.length()) as u64));
                }
            }
            eprintln!("[copy_assign] --genome-wide: {} contig(s) with mapped reads", v.len());
            v
        }
        _ => anyhow::bail!("provide exactly one of --region, --regions or --genome-wide"),
    };
    let mut by_contig: std::collections::BTreeMap<String, Vec<(u64, u64)>> = std::collections::BTreeMap::new();
    for (c, lo, hi) in regions {
        by_contig.entry(c).or_default().push((lo, hi));
    }
    // Each top-level (contig, lo, hi) becomes its own independent `RegionWork` — its own BAM query, its own
    // certificate, its own `read_provenance` rows. Two overlapping top-level windows would each independently
    // fetch and report the same physical alignment records in the overlap zone (duplicate provenance rows;
    // `--read-provenance`'s doc comment promises exactly one row per record). Validated at the boundary
    // (`--region`/`--regions` is user input), before any read is touched, rather than silently double-counted.
    validate_no_overlapping_regions(&by_contig)?;
    eprintln!(
        "[copy_assign] sweeping {} region(s) over {} contig(s)",
        by_contig.values().map(|v| v.len()).sum::<usize>(),
        by_contig.len()
    );

    // ---- O1 -> O2 FILE CONTRACT (`--families`) -------------------------------------------------------
    // The catalog is loaded, validated and BOUND TO REGIONS here, before a single read is touched, so a
    // malformed or unassignable roster fails in the first second rather than after an hour of alignment.
    // `region_families` is keyed by the exact `(contig, lo, hi)` triple the sweep iterates, so each region's
    // worker gets its own families with no re-derivation and no ambiguity about which region owns a family.
    let (region_families, catalog_seqs, region_windows) = load_supplied_families(&args, &by_contig)?;
    let catalog_index: Option<CatalogIndex> = region_families.as_ref().map(build_catalog_index);

    // O3: built ONCE, read-only across every parallel region worker (same pattern as genome_cache/
    // bam_cache below) -- only when the flag is set, so the unset path pays nothing.
    let o3_all_units_by_chrom: std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>> =
        if args.flag_missing_copies {
            let mut m: std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>> = std::collections::BTreeMap::new();
            if let Some(rf) = &region_families {
                for fams in rf.values() {
                    for f in fams {
                        for c in &f.copies {
                            m.entry(c.chrom.clone()).or_default().push((c.start, c.end, c.family_id.clone(), c.copy_idx.to_string()));
                        }
                    }
                }
            }
            m
        } else {
            std::collections::BTreeMap::new()
        };
    // Reuses `annotation` (parsed once above, unconditionally, for the AnnotationUnknown axis) instead of
    // re-reading/re-parsing --gff a second time -- same (chrom, start, end) rows the O3 orphan-locus
    // classifier needs, just bucketed by chrom. `annotation` is borrowed, not consumed, here.
    let o3_genes_by_chrom: std::collections::BTreeMap<String, Vec<(u64, u64)>> = if args.flag_missing_copies {
        let mut m: std::collections::BTreeMap<String, Vec<(u64, u64)>> = std::collections::BTreeMap::new();
        if let Some(ann) = &annotation {
            for (chrom, s, e) in ann {
                m.entry(chrom.clone()).or_default().push((*s, *e));
            }
        }
        m
    } else {
        std::collections::BTreeMap::new()
    };
    // Fix 2 (Task 6, carried forward from Task 5's review): gated behind the same flag as its sibling
    // indices `o3_all_units_by_chrom`/`o3_genes_by_chrom` above, for literal consistency with the stated
    // design constraint ("no allocation happens on the unset path") -- a 3-field POD with no side effect,
    // so this changes nothing observable, just tidiness.
    let o3_params = if args.flag_missing_copies {
        rustle::vg_family::missing_copy_flag_pass::O3Params {
            alpha: args.missing_copy_alpha,
            max_reads: args.missing_copy_max_reads,
            min_reads: 3,
        }
    } else {
        rustle::vg_family::missing_copy_flag_pass::O3Params::default()
    };

    let lambda = resolve_lambda(args.lambda_global, args.lambda_file.as_deref().and_then(read_lambda_file));
    let mut cfg = DenovoConfig::from_env();
    cfg.detect.len_cap = args.max_poa_len; // poasta memory threshold: above it, the bounded LCS fallback
    cfg.vg_realign = args.vg_realign || args.vg_realign_correct; // VG correction leg (re-thread hard reads)
    cfg.vg_realign_admit = args.vg_realign; // roster-widening admission: only via the full --vg-realign
    cfg.homology_primary = args.homology_primary; // E_r membership; off => the E_c path is untouched
    cfg.refine = !args.no_refine; // mutual-homology family gate (matches gw_family_catalog); on by default
    cfg.filter_readthrough = !args.keep_readthrough; // unspliced pre-mRNA spans are not copies
    cfg.gate.pool_locus_support = !args.no_pool_locus_support;
    cfg.tied_seed = args.tied_seed; // seed from AS-tied secondaries
    cfg.collapse_gate = args.collapse_gate; // experimental; detects paralogy, not collapse (see module header)
    cfg.collapse_enumerate = args.collapse_enumerate || cfg.collapse_enumerate; // CLI OR env (RUSTLE_COLLAPSE_ENUMERATE)
    if cfg.collapse_enumerate {
        eprintln!(
            "[copy_assign] WARNING: --collapse-enumerate / RUSTLE_COLLAPSE_ENUMERATE is a no-op in copy_assign \
             (no consumer in this binary); run gw_family_catalog --collapse-enumerate instead."
        );
    }
    if let Some(e) = args.eps_amb {
        cfg.eps_amb = Some(e);
    }
    let params = AssignParams {
        margin: args.margin,
        error_rate: args.error_rate,
        alpha: args.alpha,
        use_margin_gate: args.margin_gate,
        rna_editing_filter: !args.no_editing_filter,
        edit_rate: args.edit_rate,
        iterative_prune: args.iterative_prune,
        junction_conflict_abstain: args.junction_conflict_abstain,
        molecule_pool: args.molecule_observations && !args.no_molecule_observations,
        origin_subst_only: args.origin_substitutions_only,
        origin_drop_indels: args.origin_drop_indels && !args.no_origin_drop_indels,
        indel_psv: args.indel_psv,
        indel_psv_min_len: args.indel_psv_min_len,
        read_star_junctions: args.read_star_junctions,
        read_star_genomic: args.read_star_genomic && !args.read_star_unit,
        read_star_catalog_locus: !args.read_star_pad_locus,
        sole_candidate: !args.no_sole_candidate,
        best_by_psv: !args.best_by_alignment,
        best_by_duel: args.best_by_duel && !args.no_best_by_duel,
        dump_star: args.dump_star,
        read_star_hit_in_unit: !args.no_read_star_hit_in_unit,
        read_star_two_form: !args.read_star_genomic_only,
        read_star_readthrough: !args.no_readthrough_certificate,
        origin_consistency_check: args.origin_consistency_check,
        ..AssignParams::default()
    };
    eprintln!("[copy_assign] decisive-margin tau={} error_rate={}", args.margin, args.error_rate);
    let mut family_rows: Vec<FamilyRow> = Vec::new();
    let mut placement_assigned_total = 0usize; // §6fq: uncontested molecules assigned to their placement
    let mut primary_local_rows = 0usize; // register 734: the enforced per-family denominator
    let mut assign_rows: Vec<AssignRow> = Vec::new();
    let mut posterior_lines: Vec<String> = Vec::new();
    // EM-abundance prior for the posterior (else uniform).
    let prior_abundance = std::env::var("RUSTLE_POSTERIOR_PRIOR").ok().as_deref() == Some("abundance");
    // Cross-family reconciliation mode (see `XfamMode`). Parsed HERE, before any read is touched, so an
    // unrecognized value fails in the first second rather than silently running as `off`.
    let xfam_mode = XfamMode::from_env()?;
    if args.union_certificate && args.no_as_tied_only {
        anyhow::bail!(
            "--union-certificate is defined on AS-tied molecules (their tied placements are the candidates) and \
             needs the AS-tied gate; drop --no-as-tied-only."
        );
    }
    // locus from a de-novo tid `DN_<chrom>_<start>_<n>` (chrom may contain `_`, so split from the right).
    fn parse_locus(tid: &str) -> Option<(String, u64)> {
        let rest = tid.strip_prefix("DN_")?;
        let parts: Vec<&str> = rest.rsplitn(3, '_').collect(); // [n, start, chrom]
        Some((parts.get(2)?.to_string(), parts.get(1)?.parse().ok()?))
    }
    let mut quant_rows: Vec<QuantRow> = Vec::new();
    // O3 Phase 2 (Task 6): accumulated across the WHOLE serial drain (every region, every family) -- the
    // genome-wide Bonferroni flag threshold in `finalize_flags` can only be computed once every region has
    // drained, so nothing downstream of `compute()` can act on these until the loop below finishes.
    let mut o3_all_raw_pairs: Vec<rustle::vg_family::missing_copy_flag_pass::RawPair> = Vec::new();
    let mut o3_all_orphan_loci: Vec<rustle::vg_family::missing_copy_flag_pass::OrphanLocus> = Vec::new();
    // `--discover-copies`: read-seeded candidate copies found while scanning each region, accumulated the
    // same way as the O3 vectors above -- `RegionWork.discovered` is already gated on `args.discover_copies`
    // at the `compute()` call site, so this just drains whatever each region produced.
    let mut all_discovered: Vec<rustle::vg_family::copy_discovery::DiscoveredCopy> = Vec::new();
    // `--union-certificate`: every region's union rows + counts, drained in region order (side file + summary).
    let mut union_all = rustle::vg_family::denovo_pipeline::UnionSummary::default();
    // `--families`: one row per ASSIGNED copy, naming the catalog row it came from. The explicit join
    // between `<out>.quant.tsv` and the O1 `copies.tsv`, and the place a copy that failed to survive
    // assignment would be visible as a missing row.
    //
    // O3 (Task 6): carries its own join key (`family_id`/`copy_idx`) alongside the pre-formatted `line`,
    // so the 5 `o3_*` columns can be appended at write time -- the genome-wide flag threshold is only known
    // AFTER every region has drained, long after each row's `line` string was built.
    struct JoinRow {
        line: String,
        family_id: String,
        copy_idx: String,
    }
    let mut join_rows: Vec<JoinRow> = Vec::new();
    let mut mosaic_rows: Vec<MosaicRow> = Vec::new();
    let mut famcn_rows: Vec<FamCnRow> = Vec::new(); // reference-free chi_H + depth_cn (always emitted)
    let mut copyconv_rows: Vec<CopyConvRow> = Vec::new();
    let mut psv_read_lines: Vec<String> = Vec::new(); // --dump-psv: per-read genotype (alleles at every PSV col)
    let mut psv_copy_lines: Vec<String> = Vec::new(); // --dump-psv: per-copy PSV alleles
    let mut psv_col_lines: Vec<String> = Vec::new(); // --dump-psv: PSV column -> genome position
    let mut em_lines: Vec<String> = Vec::new();           // --em: per-read soft posterior + K-frontier label
    let mut em_abundance_lines: Vec<String> = Vec::new();  // --em: per-copy recovered abundance
    let mut phase_block_lines: Vec<String> = Vec::new();  // --phase: one phase set (PS) per family
    let mut phased_hap_lines: Vec<String> = Vec::new();   // --phase: each haplotype's PSV variant string
    let mut phased_read_lines: Vec<String> = Vec::new();  // --phase: read -> haplotype (HP) haplotag
    // --phase: a self-contained variation graph (GFA) of the phasing — PSV columns = BUBBLES
    // (one segment per allele), copies = PATHS through the bubbles. Loadable in Bandage/vg.
    let mut gfa_segs: HashSet<String> = HashSet::new();        // dedup'd S-lines (shared allele = shared node = bubble anchor)
    let mut gfa_links: HashSet<String> = HashSet::new();        // dedup'd full "L\t..." strings (copy_graph emits complete lines)
    let mut gfa_paths: Vec<String> = Vec::new();
    // VG read-threading (the Canzar flip, materialized): each read WALKS the PSV-bubble nodes for the alleles
    // it observes, REUSING a copy's node wherever their alleles agree — so multimapping reads become shared
    // threaded evidence through the one family graph. W-lines (GFA 1.1) + a Bandage node/path colour CSV.
    let mut gfa_walks: Vec<String> = Vec::new();
    let mut gfa_colors: Vec<String> = Vec::new(); // "name,colour" for Bandage (copies distinct, reads by assigned copy)
    let mut legend_rows: Vec<String> = Vec::new(); // "status\tcolour" (de-duplicated at write time)
    // --phase v2: one exon presence/absence graph per family (built during the drain, where `fa` is in
    // scope; sequence-free — `to_gfa` fetches reference bases lazily at write time via `genome_for`).
    let mut exon_graphs: Vec<rustle::vg_family::copy_graph::ExonGraph> = Vec::new();
    let mut fallback_all: Vec<FallbackEdge> = Vec::new(); // family edges confirmed via the LCS fallback
    let mut dna_needs_rows: Vec<DnaNeedsRecord> = Vec::new(); // --absent-copies: candidates needing DNA validation
    let mut prov_rows: Vec<String> = Vec::new(); // --read-provenance: one row per AS-tied alignment record
    // --absent-copies + opt-in --linearize/--linearize-gate: linearize certificates, one per Stage-2-admitted
    // candidate (Task 4), written to `<out>.linearize.tsv` below (Task 5) when `do_linearize`. Empty otherwise
    // (the cert is skipped in `detect_and_assign`). `--linearize-gate` also uses the verdict to gate admission
    // itself, so a demoted candidate shows up here but not in `fams`.
    let mut linearize_certs_all: Vec<(String, LinearizeCertificate, (String, u64, u64))> = Vec::new();
    let mut vg_realign_lines: Vec<String> = Vec::new(); // --vg-realign: per-read re-align decisions (report-only)
    let mut gfam = 0usize; // global family counter (unique ids across regions)
    let mut gtf_lines: Vec<String> = Vec::new(); // --gtf: FLAIR-style isoform GTF (transcript + exon rows)
    // --productivity: (attribute string, transcript id, ORF in aa) — the `productive` call needs the family's
    // best ORF, which is only known after every region is drained, so it is a second pass over the GTF below
    let mut prod_rows: Vec<(String, String, String, usize)> = Vec::new(); // family, copy, transcript, ORF aa

    // `--skip-poa-diagnostic` is read by `detect_and_assign` via this env var (it is purely diagnostic and
    // does not change the emitted families/assignments — see the flag's help).
    if args.skip_poa_diagnostic {
        std::env::set_var("RUSTLE_SKIP_POA_DIAGNOSTIC", "1");
    }
    // `--poa-cap` is read by `discover_intron_psvs` (copy_assign_pipeline.rs) via this env var, the same
    // "flag -> env var -> deep read" idiom as `--skip-poa-diagnostic` above: the const it replaces lives many
    // call frames below `main` (through `assign_family`/`assign_family_detailed`/`find_weak_copies`), so a
    // signature thread-through would touch dozens of call sites for an opt-in (`RUSTLE_INTRON_PSV=1`) code
    // path. Always set (not gated on non-default) so the resolved value is unambiguous; the default 20000
    // reproduces the prior hard-coded constant exactly.
    std::env::set_var("RUSTLE_POA_CAP", args.poa_cap.to_string());
    if args.psv_genomic {
        std::env::set_var("RUSTLE_PSV_GENOMIC", "1");
    }
    // `--absent-min-clusters` reaches `absent_copy::AbsentCopyParams::from_env` through the same
    // "flag -> env var -> deep read" idiom (the params are built inside `detect_and_assign`'s
    // candidate loop, many frames below `main`). Set ONLY when the flag is given, so (a) an
    // unset run is byte-identical and (b) a caller who exported RUSTLE_ABSENT_MIN_CLUSTERS
    // directly is not silently clobbered by the flag's default.
    if let Some(n) = args.absent_min_clusters {
        std::env::set_var(rustle::vg_family::absent_copy::MIN_CLUSTERS_ENV, n.to_string());
        if !args.absent_copies {
            eprintln!(
                "[copy_assign] WARNING: --absent-min-clusters={n} has no effect without --absent-copies"
            );
        }
    }
    // `--read-cap` is a NO-OP in copy_assign (see the flag's help): `o2_materialize::READ_CAP` has no
    // consumer in any `src/bin/*.rs` binary. Warn rather than silently ignore a non-default value.
    if args.read_cap != 6_000 {
        eprintln!(
            "[copy_assign] WARNING: --read-cap={} has no consumer in this binary (o2_materialize's READ_CAP \
             is not on copy_assign's execution path — it backs a Rust byte-parity port of the Python \
             genome-wide-catalog materializer that no shipped binary imports); the value is ignored.",
            args.read_cap
        );
    }
    let timing = std::env::var_os("RUSTLE_TIMING").is_some();
    // Parse the BAM index + header ONCE and reuse across every region (the per-region path re-parses the
    // multi-MB `.bai` otherwise). None => no usable index; fall back to the per-region open (which scans).
    let bam_cache = BamIndexCache::open(&args.bam).ok();
    // Region-parallel pool (opt-in via --region-threads > 1). Sized to region_threads; the inner per-copy
    // poasta parallelism (discover_psvs) composes on the SAME pool, so total concurrency is bounded to N.
    let region_pool = if args.region_threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.region_threads)
                .build()
                .context("building region thread pool")?,
        )
    } else {
        None
    };
    // FLAT region list across ALL contigs, in the deterministic by_contig order (sorted contig, file-order
    // ranges). Out-of-order parallel processing over this flat list lets the globally-heaviest families —
    // which live on DIFFERENT contigs — overlap, while the serial drain below (in this same order) keeps
    // CAFAM ids + every row byte-identical to the serial sweep.
    // Cross-chromosome families (2026-09-15) are bound to a SYNTHETIC key (`~xchrom~<family_id>`, never a
    // real contig from `by_contig`/`--regions`) that exists only in `region_families`'s keys — append them
    // here so the sweep actually visits them. `compute` recognizes such a key purely by its `region_windows`
    // entry tagging every window with the REAL chromosome to fetch from (never the key's own placeholder).
    let flat: Vec<(String, u64, u64)> = by_contig
        .iter()
        .flat_map(|(c, ranges)| ranges.iter().map(move |&(lo, hi)| (c.clone(), lo, hi)))
        .chain(
            region_families
                .iter()
                .flat_map(|rf| rf.keys())
                .filter(|k| k.0.starts_with("~xchrom~"))
                .cloned(),
        )
        .collect();
    // Bounded LRU cache of loaded contig genomes — so a worker on any contig reuses an already-loaded genome
    // instead of reloading, and at most ~capacity contig sequences are resident (the memory bound; Arc keeps
    // a genome alive while an evicting worker still uses it). Capacity tracks the concurrency.
    let genome_cap = NonZeroUsize::new((args.region_threads + 1).max(2)).unwrap();
    let genome_cache: Arc<Mutex<LruCache<String, Arc<GenomeIndex>>>> =
        Arc::new(Mutex::new(LruCache::new(genome_cap)));
    let genome_for = |contig: &str| -> Result<Arc<GenomeIndex>> {
        if let Some(g) = genome_cache.lock().unwrap().get(contig).cloned() {
            return Ok(g);
        }
        // load OUTSIDE the lock (a chromosome load is seconds; never block other workers on it). A rare
        // double-load on a concurrent miss is harmless — the second insert just wins.
        let contigs: HashSet<String> = std::iter::once(contig.to_string()).collect();
        let g = Arc::new(
            GenomeIndex::from_fasta_contigs(&args.fasta, &contigs)
                .with_context(|| format!("loading {} for {contig}", args.fasta))?,
        );
        genome_cache.lock().unwrap().put(contig.to_string(), g.clone());
        Ok(g)
    };
    // A region's genome, generalized to however many chromosomes its windows actually touch. The
    // single-element case (every region before cross-chromosome families existed) delegates to the
    // cached `genome_for` above and is therefore byte-for-byte the same load as before; a cross-chromosome
    // family's multi-element set builds one combined, uncached `GenomeIndex` instead (cross-chromosome
    // families are rare enough that a dedicated per-contig cache slot for them is not worth the complexity).
    let genome_for_multi = |contigs: &std::collections::BTreeSet<String>| -> Result<Arc<GenomeIndex>> {
        match contigs.len() {
            1 => genome_for(contigs.iter().next().expect("len == 1")),
            _ => {
                let wanted: HashSet<String> = contigs.iter().cloned().collect();
                Ok(Arc::new(
                    GenomeIndex::from_fasta_contigs(&args.fasta, &wanted).with_context(|| {
                        format!("loading {} for cross-chromosome contigs {:?}", args.fasta, contigs)
                    })?,
                ))
            }
        }
    };
    // The expensive, INDEPENDENT per-region work: BAM read + detect_and_assign (the dominant poasta alignment
    // lives here). Pure w.r.t. the read-only genome/bam_cache. The heavy read SEQUENCES are dropped here —
    // only the read NAMES + computed `fams` are returned — so collecting every region's result is lightweight.
    let compute = |contig: &String, lo: u64, hi: u64| -> Result<RegionWork> {
        // §6dh: on the --families path, gather from the supplied copies' own neighbourhoods rather than
        // the whole bound region — a dispersed family's hull can be tens of Mb while its copies occupy a
        // few hundred kb, and a read overlapping no copy can never be assigned to one. Cross-chromosome
        // families (2026-09-15) push windows tagged with a chromosome OTHER than `contig` (their synthetic
        // key's own "chromosome" is a placeholder, not a real one) — already clipped to the one real region
        // each window's own chromosome was found inside (`load_supplied_families`), so no clipping happens
        // here any more.
        let wins: Vec<(String, u64, u64)> = region_windows
            .as_ref()
            .and_then(|w| w.get(&(contig.clone(), lo, hi)))
            .cloned()
            .unwrap_or_else(|| vec![(contig.clone(), lo, hi)]);
        // The genome this region's assignment needs: every chromosome any window actually reads from —
        // for every region before cross-chromosome families existed this is the single-element set
        // `{contig}`, so `genome_for_multi` delegates to the ORIGINAL cached single-contig `genome_for`
        // and behaviour is unchanged; a cross-chromosome family's windows pull in its other chromosomes.
        let win_contigs: std::collections::BTreeSet<String> = wins.iter().map(|(c, _, _)| c.clone()).collect();
        let genome = genome_for_multi(&win_contigs)?;
        let t_read = std::time::Instant::now();
        // §6zb streaming pass-1 (`docs/PREREG_streaming_assembly_2026-09-23.md`): under --assemble-only with
        // no read-level extras, reduce every record on arrival and never materialise the reads. Every other
        // configuration takes the historical path below unchanged.
        let streaming = args.assemble_only
            && !args.materialize_reads
            && args.read_isoform_k == 0
            && !rustle::vg_family::denovo_assemble::footprint_nodes_enabled()
            && !(args.recover_copies || args.tied_seed)
            // GOOD seeding (r1060/r1100) streams too once a genome-wide best-AS table is loaded: the
            // streaming reader applies `AS >= ratio x table best` itself (`stream_pass1_region`); without a
            // table the ratio needs the region's buffered records to know a local best, as before.
            && (rustle::vg_family::denovo_assemble::gtf_secondary_as_ratio() <= 0.0
                || rustle::vg_family::denovo_assemble::global_best_as().is_some());
        let mut streamed: Option<Vec<rustle::vg_family::denovo_assemble::Skeleton>> = None;
        let mut n_mapped_streamed = 0usize;
        if streaming {
            let mut acc = rustle::vg_family::denovo_assemble::Pass1Acc::new(1, None);
            let mut fetched: Vec<(String, u64, u64)> = Vec::new();
            for (wchrom, wlo, whi) in &wins {
                n_mapped_streamed += rustle::vg_family::denovo_assemble::stream_pass1_region(
                    &args.bam, wchrom, *wlo, *whi,
                    rustle::vg_family::denovo_assemble::gtf_secondary_enabled(),
                    !args.keep_coordinate_duplicates,
                    if args.keep_coordinate_duplicates { &fetched } else { &[] },
                    &mut acc,
                )
                .with_context(|| format!("streaming {wchrom}:{wlo}-{whi}"))?;
                fetched.push((wchrom.clone(), *wlo, *whi));
            }
            streamed = Some(acc.finish(cfg.pass1_min_reads, 0, None));
        }
        let (primary, mut bam_reads) = if streaming { (Vec::new(), Vec::new()) } else {
            let mut pr: Vec<_> = Vec::new();
            let mut br: Vec<_> = Vec::new();
            // Windows already fetched for this region. A read spanning a window boundary is returned by
            // both queries (every reader yields each record overlapping `[lo, hi)`), so a record from
            // window i that overlaps an EARLIER window j was necessarily fetched by j and is skipped here.
            //
            // ⚠ §6z6 (`docs/PREREG_primary_dedupe_2026-09-22.md`): this used to key on
            // `(chrom, ref_start, ref_end, intron_chain)` — PrimaryRead has no name — which also collapsed
            // DISTINCT molecules with identical coordinates. On A119b that key dropped 25-32% of primary
            // records per chromosome and pushed ~23% of the 2-read chains below pass-1's floor (chr20:
            // 3,732 of 16,166; traced at ZBTB21, two MAPQ-60 reads with the exact RefSeq chain). With a
            // single window — the whole `--assemble-only` path — nothing is dropped now.
            //
            // ⛔ MEASURED (same prereg, OUTCOME): the raw assembly does gain chains (+81 on human
            // chr20/21/22), but the SHIPPED POLISH was fitted on the de-duplicated counts and, fed the
            // true counts, removes more matching chains than the fix adds (polished −147 human dev,
            // −1,096 gorilla held-out, with gorilla precision +3.2 pts). So the historical key stays the
            // DEFAULT (byte-identical) and the window rule is opt-in via `--keep-coordinate-duplicates`.
            let mut seen = std::collections::HashSet::new();
            let mut fetched: Vec<(String, u64, u64)> = Vec::new();
            for (wchrom, wlo, whi) in &wins {
                let (wlo, whi) = (*wlo, *whi);
                let (p, b) = match &bam_cache {
                    Some(c) => c.reads_in_region(&args.bam, wchrom, wlo, whi),
                    None => reads_in_region(&args.bam, wchrom, wlo, whi, args.threads),
                }
                .with_context(|| format!("reading {wchrom}:{wlo}-{whi}"))?;
                for x in p {
                    let keep = if args.keep_coordinate_duplicates {
                        !fetched.iter().any(|(c, l, h)| c == &x.chrom && x.ref_start < *h && x.ref_end > *l)
                    } else {
                        // historical key: PrimaryRead has no name; (chrom, span, intron chain) stands in
                        seen.insert((x.chrom.clone(), x.ref_start, x.ref_end, x.introns.clone()))
                    };
                    if keep {
                        pr.push(x);
                    }
                }
                for x in b {
                    br.push(x);
                }
                fetched.push((wchrom.clone(), wlo, whi));
            }
            // ⚠ Must include `chrom`, not just `(name, ref_start)`: before cross-chromosome families
            // (2026-09-15) every record `compute` ever saw shared one contig, so `ref_start` alone was
            // already a sufficient tiebreaker. A cross-chromosome family's windows span several real
            // contigs in one call, and a read placed at the SAME offset on two of them (a real, distinct
            // alignment record each) would otherwise collide onto one key and the second record would be
            // silently dropped as a "duplicate" — exactly the kind of silent truncation this whole feature
            // exists to avoid.
            let mut bseen = std::collections::HashSet::new();
            br.retain(|x: &rustle::vg_family::denovo_assemble::BamRead| {
                bseen.insert((x.name.clone(), x.chrom.clone(), x.read.ref_start))
            });
            (pr, br)
        };
        if timing && wins.len() > 1 {
            eprintln!(
                "[timing] {contig}:{lo}-{hi} gathered from {} copy window(s) ({:.1} Mb across {} \
                 chromosome(s), of {:.1} Mb hull)",
                wins.len(),
                wins.iter().map(|(_, a, b)| (b - a) as f64).sum::<f64>() / 1e6,
                win_contigs.len(),
                (hi - lo) as f64 / 1e6
            );
        }
        if timing {
            eprintln!(
                "[timing] reads_in_region {contig}:{lo}-{hi} ({} reads): {:.1}s",
                bam_reads.len(),
                t_read.elapsed().as_secs_f64()
            );
        }
        let extra = if args.recover_copies || args.tied_seed {
            tied_secondary_reads_in_region(&args.bam, contig, lo, hi, args.as_ratio).unwrap_or_default()
        } else {
            Vec::new()
        };
        // `--families`: materialize THIS region's supplied catalog families into the copy set, and enforce
        // the last two contract clauses (a sequence exists for every copy; every copy has reads here).
        // Done inside the worker because it needs this region's genome and this region's BAM slice.
        let supplied: Option<Vec<ColocatedFamily>> = match &region_families {
            None => None,
            Some(rf) => {
                let mine = rf.get(&(contig.clone(), lo, hi)).map(|v| v.as_slice()).unwrap_or(&[]);
                let mut v: Vec<ColocatedFamily> = Vec::with_capacity(mine.len());
                for f in mine {
                    let (cf, _src) = to_colocated(f, catalog_seqs.as_ref(), &genome)?;
                    for c in &cf.copies {
                        // A supplied copy with no read here cannot be assigned anything, and dropping it
                        // silently would understate K and loosen the Bonferroni certificate for the rest.
                        let n = bam_reads
                            .iter()
                            .filter(|br| {
                                br.chrom == c.chrom
                                    && br.read.ref_start < c.end
                                    && read_ref_end_local(&br.read) > c.start
                            })
                            .count();
                        // §6ft: a catalog copy the catalog itself marks unexpressed (`n_reads 0`, an annotated model
                        // kept as the unit) or a partner may legitimately have no read here — it stays a target
                        let catalog_zero = f.copies.iter().any(|cc| cc.tid == c.tid && (cc.n_reads == 0 || cc.partner));
                        if n == 0 && !catalog_zero {
                            anyhow::bail!(
                                "--families: {} copy {} ({}:{}-{}) has NO reads in {contig}:{lo}-{hi} of \
                                 {}. It cannot be assigned, and silently dropping it would understate the \
                                 family's copy count. Check the BAM is the one the catalog was built from \
                                 (subset BAMs are the recurring trap here).",
                                cf.family_id, c.tid, c.chrom, c.start, c.end, args.bam
                            );
                        }
                    }
                    v.push(cf);
                }
                Some(v)
            }
        };
        // ⭐ THE AS-TIED GATE (default on). AS evidence needs every placement of a molecule, so it is
        // computed on the full record set; then non-tied molecules leave BEFORE the certificate. A unique
        // mapper is not O2's business and is never assigned. `--no-as-tied-only` skips this block.
        let mut uniq_reads: Vec<(String, u64, u64, Vec<(u64, u64)>)> = Vec::new();
        if !args.no_as_tied_only && !streaming {
            let ev = as_evidence_per_read(&bam_reads, !args.no_as_tied_only);
            let mut tied: std::collections::HashSet<&str> = std::collections::HashSet::new();
            let mut all: std::collections::HashSet<&str> = std::collections::HashSet::new();
            for (br, e) in bam_reads.iter().zip(ev.iter()) {
                all.insert(br.name.as_str());
                if as_tied(e, args.as_tie_ratio) {
                    tied.insert(br.name.as_str());
                }
            }
            let (n_all, n_tied, n_rec) = (all.len(), tied.len(), bam_reads.len());
            let mut tied_owned: std::collections::HashSet<String> = tied.into_iter().map(|s| s.to_string()).collect();
            // ⭐ §6gz: which tied molecules have a tied placement OUTSIDE every supplied family UNIT?
            // ⚠ The test is against the UNIT SPAN (`start`/`end`, the read-supported exon chain), NOT the
            // padded read-star locus: the locus is exactly what swallowed EIF3C into NPIP copy 16 (its locus
            // 28,982,252–29,053,456 contains EIF3C; its unit starts at 29,016,080), and a first form of this
            // detector that used the locus flagged 541 of the 4,706 leaks — it defined "inside" by the
            // swallowing target. A tied record overlapping no unit is a competitor O2 will never score ⟹ the
            // molecule is registered and can never be `Assigned`. Only meaningful with --families.
            let mut n_outside = 0usize;
            if let Some(sup) = supplied.as_deref() {
                let targets: Vec<(String, u64, u64)> = sup
                    .iter()
                    .flat_map(|f| f.copies.iter())
                    .map(|c| (c.chrom.clone(), c.start, c.end))
                    .collect();
                // parallel to `targets`: the catalog copy index of each unit (`--gtf-copy-set`)
                let target_idx: Vec<String> = sup
                    .iter()
                    .flat_map(|f| f.copies.iter())
                    .map(|c| catalog_index.as_ref().and_then(|ix| ix.get(&c.tid)).map(|(_, i)| i.to_string()).unwrap_or_default())
                    .collect();
                let best_as: std::collections::HashMap<&str, i32> = bam_reads
                    .iter()
                    .filter(|br| tied_owned.contains(&br.name) && !br.is_supplementary)
                    .fold(std::collections::HashMap::new(), |mut m, br| {
                        let e = m.entry(br.name.as_str()).or_insert(br.as_score);
                        *e = (*e).max(br.as_score);
                        m
                    });
                let mut flagged: std::collections::HashSet<&str> = std::collections::HashSet::new();
                for br in bam_reads.iter().filter(|br| tied_owned.contains(&br.name) && !br.is_supplementary) {
                    if br.as_score < best_as[br.name.as_str()] {
                        continue; // not one of the tied placements
                    }
                    let (s0, e0) = (br.read.ref_start, read_ref_end_local(&br.read));
                    let inside = targets.iter().any(|(c, a, b)| *c == br.chrom && s0 < *b && e0 > *a);
                    if !inside {
                        flagged.insert(br.name.as_str());
                        // A6: name the outside placement's own locus, not just the fact it exists.
                        // `s0` is 0-based (`AlignedRead::ref_start`); every other coordinate this binary
                        // emits into a GTF attribute is 1-based (the exon/transcript rows below, `+ 1`), so
                        // the registered start needs the same `+ 1` or the printed `outside:chrom:start-end`
                        // token is off by one relative to the file it sits in.
                        rustle::vg_family::copy_assign_pipeline::register_tie_outside_locus(&br.name, &br.chrom, s0 + 1, e0);
                    }
                    if (args.gtf_copy_set && !args.no_gtf_copy_set) {
                        // the copy SET of an undecided isoform (§6hn): catalog indices at the tied placements
                        let hit = targets.iter().position(|(c, a, b)| *c == br.chrom && s0 < *b && e0 > *a);
                        let mut reg = TIE_SET.get_or_init(Default::default).lock().unwrap();
                        let e = reg.entry(br.name.clone()).or_insert_with(|| (std::collections::BTreeSet::new(), false));
                        match hit.and_then(|i| target_idx.get(i)) {
                            Some(idx) => { e.0.insert(idx.clone()); }
                            None => e.1 = true,
                        }
                    }
                }
                n_outside = flagged.len();
                for n in flagged {
                    rustle::vg_family::copy_assign_pipeline::register_tie_outside(n);
                }
                // ⭐ §6hd: aligner self-disagreement. Per molecule: the unit index of its PRIMARY record and of
                // its best-AS record(s), both by unit-span overlap (same `targets` as above). Disagree ⟹ admit.
                if args.admit_aligner_disagreement {
                    let unit_of = |br: &rustle::vg_family::denovo_assemble::BamRead| -> Option<usize> {
                        let (s0, e0) = (br.read.ref_start, read_ref_end_local(&br.read));
                        targets.iter().position(|(c, a, b)| *c == br.chrom && s0 < *b && e0 > *a)
                    };
                    let mut prim: std::collections::HashMap<&str, Option<usize>> = std::collections::HashMap::new();
                    let mut best: std::collections::HashMap<&str, (i32, Vec<Option<usize>>)> = std::collections::HashMap::new();
                    for br in bam_reads.iter().filter(|br| !br.is_supplementary) {
                        let u = unit_of(br);
                        if !br.is_secondary {
                            prim.insert(br.name.as_str(), u);
                        }
                        let e = best.entry(br.name.as_str()).or_insert((br.as_score, Vec::new()));
                        if br.as_score > e.0 {
                            *e = (br.as_score, vec![u]);
                        } else if br.as_score == e.0 {
                            e.1.push(u);
                        }
                    }
                    let mut n_dis = 0usize;
                    for (name, pu) in &prim {
                        if tied_owned.contains(*name) {
                            continue; // already admitted by the AS tie
                        }
                        let Some(pu) = pu else { continue }; // primary outside every unit: not this family's
                        let Some((_, bus)) = best.get(name) else { continue };
                        // disagreement = the primary's unit is NOT among the best-AS units, and some best-AS
                        // record IS inside a unit (a best-AS placement outside every unit is the §6gz case)
                        if !bus.contains(&Some(*pu)) && bus.iter().any(|u| u.is_some()) {
                            tied_owned.insert(name.to_string());
                            register_disagreement(name);
                            n_dis += 1;
                        }
                    }
                    GATE_MOL_DISAGREE.fetch_add(n_dis, std::sync::atomic::Ordering::Relaxed);
                }
            }
            GATE_MOL_OUTSIDE.fetch_add(n_outside, std::sync::atomic::Ordering::Relaxed);
            if (args.gtf_copy_set && !args.no_gtf_copy_set) {
                for br in bam_reads.iter().filter(|br| !tied_owned.contains(&br.name) && !br.is_secondary && !br.is_supplementary) {
                    let bl = aligned_blocks_local(&br.read);
                    let chain: Vec<(u64, u64)> = bl.windows(2).map(|w| (w[0].1, w[1].0)).filter(|&(a, b)| b > a).collect();
                    uniq_reads.push((br.chrom.clone(), br.read.ref_start, read_ref_end_local(&br.read), chain));
                }
            }
            bam_reads.retain(|br| tied_owned.contains(&br.name));
            GATE_MOL_ALL.fetch_add(n_all, std::sync::atomic::Ordering::Relaxed);
            GATE_MOL_TIED.fetch_add(n_tied, std::sync::atomic::Ordering::Relaxed);
            GATE_REC_ALL.fetch_add(n_rec, std::sync::atomic::Ordering::Relaxed);
            GATE_REC_TIED.fetch_add(bam_reads.len(), std::sync::atomic::Ordering::Relaxed);
        }
        let t_da = std::time::Instant::now();
        // §6p6 --assemble-only: the assembly path below needs `primary` and nothing detect_and_assign
        // produces, so skip it outright. Empty results keep every downstream writer on its normal path.
        let (mut fams, fallback, dna_needs, linearize_certs) = if args.assemble_only {
            (Vec::new(), Vec::new(), Vec::new(), Vec::new())
        } else {
            detect_and_assign(
                &primary, &bam_reads, &genome, &cfg, args.win, args.min_copies, &params, &extra,
                args.absent_copies, do_linearize, args.linearize_gate, &args.fasta,
                supplied.as_deref(),
            )
        };
        if timing && !args.assemble_only {
            eprintln!("[timing] detect_and_assign {contig}:{lo}-{hi}: {:.1}s", t_da.elapsed().as_secs_f64());
        }
        // ⭐ --union-certificate: one certificate over the union of each cross-family / outside-tied
        // molecule's candidates, applied to every family's row IN PLACE (`fa.assignments`), so the four
        // status emit sites in the drain agree by construction. It must run HERE: the read sequences the
        // certificate aligns exist only inside this worker (`RegionWork` drops them), which is why
        // `xfam_pass1` -- which runs later, without them -- can only report or demote, never re-score.
        let union = match (args.union_certificate && !args.assemble_only, supplied.as_deref()) {
            (true, Some(sup)) => {
                let t_u = std::time::Instant::now();
                let label = |fid: &str, tid: &str, ci: usize| -> String {
                    match catalog_index.as_ref().and_then(|ix| ix.get(tid)) {
                        Some((cf, idx)) => format!("{cf}:{idx}"),
                        None => format!("{fid}:#{ci}"),
                    }
                };
                let s = rustle::vg_family::denovo_pipeline::union_certificate_pass(
                    &mut fams, sup, &bam_reads, &genome, &params, &label,
                );
                eprintln!(
                    "[union] {contig}:{lo}-{hi}: {} molecule(s) in scope, {} scored in {} group(s): assigned to a \
                     family {}, to an outside locus {}, tied {}, ambiguous {}, no result {}; {} left as today (a \
                     tie partner in another region: no primary record loaded here); {} outside pseudo-cop(y/ies) \
                     unbuildable; {} row(s) added ({:.1}s)",
                    s.n_in_scope, s.n_scored(), s.n_groups, s.n_assigned_family, s.n_assigned_outside, s.n_tied,
                    s.n_ambiguous, s.n_no_result, s.n_other_region, s.n_pseudo_unbuildable, s.n_rows_added,
                    t_u.elapsed().as_secs_f64()
                );
                s
            }
            _ => rustle::vg_family::denovo_pipeline::UnionSummary::default(),
        };
        // FLAIR-style isoform assembly for the optional GTF (intron-chain collapse -> gate -> gene grouping).
        // Recomputed here only under --gtf (cheap: pass1/gate are ~0s); independent of the assignment.
        let transcripts: Vec<TranscriptRec> = if args.gtf {
            // §6m5 / PREREG 6d586b2d: read-isoform widening. `--read-isoform-k 0` (the default) is the
            // explicit no-op, so this line is byte-identical to the previous `pass1_skeletons` call.
            let skeletons = match streamed.take() {
                Some(s) => s,
                None => pass1_skeletons_widened(&primary, cfg.pass1_min_reads, 1, None, args.read_isoform_k),
            };
            // Opt-in (RUSTLE_JUNCTION_FUZZ_BP, default off): merge skeletons whose intron chains match in
            // count and differ only by a pre-registered per-junction tolerance --
            // docs/PREREG_junction_fuzz_2026-09-15.md, docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md.
            // This IS the "--gtf pure de novo path" call site the design targets (traced 2026-09-15: the
            // `detect_and_assign`/`pass1_skeletons_robust` skeletons feed only the multi-copy family/O1
            // oracle below and are never read by this block). Zero effect when unset (tolerance 0 is
            // `merge_fuzzy_skeletons`'s own explicit no-op), so every existing catalog stays byte-identical.
            // ⚠ MEASURED NET-NEGATIVE at the pre-registered 672bp tolerance on real chr20 data: matching
            // intron chains 345 -> 284, matching transcripts 347 -> 286 vs baseline
            // (`bench/CHR20_ASSEMBLER_COMPARISON.md`'s "fuzzy junction merging" follow-up). Left off by
            // default for this reason, not merely because it is untested -- do not enable it without
            // re-registering a new tolerance via a fresh, separate experiment.
            let fuzz_bp: u64 = std::env::var("RUSTLE_JUNCTION_FUZZ_BP")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0);
            let skeletons = if fuzz_bp > 0 { merge_fuzzy_skeletons(skeletons, fuzz_bp) } else { skeletons };
            // §6m6 follow-up: localise where pass-1 skeletons die before the GTF. `RUSTLE_GATE_CENSUS=1`
            // only PRINTS — the transcripts are the same objects either way.
            let iso = if matches!(std::env::var("RUSTLE_GATE_CENSUS"), Ok(v) if v != "0" && !v.is_empty()) {
                let use_rs = matches!(std::env::var("RUSTLE_READ_STRAND"), Ok(v) if v != "0" && !v.is_empty());
                let margin: f64 = std::env::var("RUSTLE_READ_STRAND_MARGIN")
                    .ok().and_then(|v| v.parse().ok()).unwrap_or(0.90);
                let (iso, c) = assemble_gate_census(&skeletons, &genome, &cfg.gate, use_rs, margin);
                eprintln!(
                    "[gate-census] {} skeletons -> kept {} | rejected: reads {} span {} seq(motif/coords) {} len {}",
                    c.total(), c.kept, c.rej_reads, c.rej_span, c.rej_seq, c.rej_len
                );
                iso
            } else {
                assemble_gate(&skeletons, &genome, &cfg.gate)
            };
            let groups = collapse_loci_groups(&iso);
            iso.iter()
                .enumerate()
                .map(|(i, t)| TranscriptRec {
                    tid: t.tid.clone(),
                    gene_tid: iso[groups[i]].tid.clone(),
                    chrom: t.chrom.clone(),
                    start: t.start,
                    end: t.end,
                    strand: t.strand,
                    introns: t.introns.clone(),
                    n_reads: t.n_reads,
                })
                .collect()
        } else {
            Vec::new()
        };
        let read_names: Vec<String> = bam_reads.iter().map(|r| r.name.clone()).collect();
        let read_chrom: Vec<String> = bam_reads.iter().map(|r| r.chrom.clone()).collect();
        let read_mapqs: Vec<u8> = bam_reads.iter().map(|r| r.mapq).collect();
        let read_spans: Vec<(u64, u64, u8)> = bam_reads
            .iter()
            .map(|r| {
                (
                    r.read.ref_start,
                    read_ref_end_local(&r.read),
                    (r.is_secondary as u8) | ((r.is_supplementary as u8) << 1),
                )
            })
            .collect();
        let read_blocks: Vec<Vec<(u64, u64)>> = bam_reads.iter().map(|r| aligned_blocks_local(&r.read)).collect();
        let read_strand: Vec<char> = bam_reads
            .iter()
            .map(|r| match (r.ts, r.reverse) {
                (Some('+'), rev) => if rev { '-' } else { '+' },
                (Some('-'), rev) => if rev { '+' } else { '-' },
                (_, rev) => if rev { '-' } else { '+' },
            })
            .collect();
        let as_ev = as_evidence_per_read(&bam_reads, !args.no_as_tied_only);
        let n_mapped = if streaming { n_mapped_streamed } else { bam_reads.len() };
        // O3 (Task 5, Phase 1 of the flag-pass wiring): per-family raw pair statistics + orphan-locus scan.
        // Only runs when the flag is set -- the `if` guards every allocation and every minimap2 realign
        // call, so the unset path is untouched (byte-identical `RegionWork` in every other field).
        let (o3_raw_pairs, o3_orphan_loci): (Vec<_>, Vec<_>) = if args.flag_missing_copies {
            let mut pairs = Vec::new();
            let mut loci = Vec::new();
            // Fix 1 (Task 6, revised per review): `bench/missing_copy_flag_pass.py`'s outer loop (lines 76-82)
            // rebuilds its `A`/`targets` lookup FRESH, per family, from that ONE family's own
            // `A.assignments.tsv` -- each `fam_*` sweep is a separate, single-family `copy_assign`
            // invocation. There is no cross-family lookback anywhere in the Python: a read that actually
            // belongs to a neighboring family is let into the locus cluster and reclassified afterward by
            // the genome-wide `other_family`/`annotated_no_unit`/`unannotated` classifier
            // (`bench/missing_copy_flag_pass.py:148-149`, ported as `classify_orphan_locus` below), not excluded
            // upstream. Pooling every family in `fams` into one map (the previous version of this fix) is
            // therefore not a faithful port, and has its own bug on top: when two co-located families both
            // claim the same read (measured as real on this codebase's own data -- see
            // `copy_assign.rs:771-831`'s `xfam_pass1` docs), a plain `HashMap::insert` lets whichever
            // family is iterated last silently overwrite the other's verdict, with no reconciliation (that
            // reconciliation is `xfam_pass1`'s job, and it runs later, in the serial drain, after this
            // `compute()` closure has already finished). So this is scoped to `fa.assignments` ONLY --
            // built fresh inside the `for fa in &fams` loop, per family, exactly like the Python.
            for fa in &fams {
                let fa_verdict: std::collections::HashMap<&str, (bool, usize)> = fa
                    .assignments
                    .iter()
                    .filter_map(|&(ri, ref a)| {
                        bam_reads.get(ri).map(|br| (br.name.as_str(), (a.origin_rejected, a.n_candidates)))
                    })
                    .collect();
                // Resolve catalog_copy_idx -> (chrom, start, end, locus_extent): fa.copy_spans is indexed by
                // the SWEEP's own internal position (`ci`); catalog_copy_idx is the CATALOG's own separate
                // namespace (see `cat_idx_of`, further down in this same file) -- rebuild the mapping via
                // fa.copy_tids + catalog_index, the same pattern the assignment-row loop already uses.
                // The locus extent (§L2, `register_locus_extent`/`locus_extent_of`, populated from the
                // catalog's `locus_start`/`locus_end` columns by `catalog_input.rs`) is threaded through
                // here so `detect_missing_copy_pairs` can realign against the SAME target window
                // `bench/missing_copy_flag_pass.py`'s `detector()` does -- Task 7's reproduction-gate fix: the bare
                // copy span alone under-/over-sizes the window whenever it differs from the L2 extent.
                //
                // Fix 1 (final whole-branch review, confirmed live on Task 8's real output -- 8 rows/arm):
                // `fa` is a LOCALLY co-located physical family, which can group copies from MORE THAN ONE
                // true catalog family (`cf`). Keying `copy_span_by_catalog_idx`/`rejected_by_idx`/
                // `accepted_by_idx` by BARE `cidx` across the whole `fa` (the pre-fix design) had two bugs:
                // (a) two catalog copies from DIFFERENT `cf`s that happen to share a bare index silently
                // overwrote each other's window/read-pool entry (pool-mixing corruption), and (b) the pair
                // detector was called with `&fa.family_id` (the LOCAL id) as `RawPair.family_id`, a
                // different namespace than `JoinRow.family_id` (always the TRUE `cf`, read at that row's
                // own construction site further down) -- so the later `o3_flags.get(&(r.family_id,
                // r.copy_idx))` lookup silently missed for every such row, indistinguishable from "never
                // tested". Fixed by partitioning `fa`'s own copies by their TRUE `cf` BEFORE building any
                // lookup map: within one `cf`'s own partition a bare `cidx` is safe (catalog copy indices
                // are unique within a single catalog family), and `detect_missing_copy_pairs` is called
                // once per `cf` with that `cf` itself as `family_id`, exactly matching `JoinRow.family_id`.
                let mut copy_span_by_cf: std::collections::HashMap<
                    String,
                    std::collections::HashMap<String, (String, u64, u64, Option<(u64, u64)>)>,
                > = std::collections::HashMap::new();
                for (ci, tid) in fa.copy_tids.iter().enumerate() {
                    if let Some((cf, cidx)) = catalog_index.as_ref().and_then(|ix| ix.get(tid)) {
                        if let Some((chrom, s, e)) = fa.copy_spans.get(ci) {
                            let locus = rustle::vg_family::copy_assign_pipeline::locus_extent_of(tid);
                            copy_span_by_cf
                                .entry(cf.clone())
                                .or_default()
                                .insert(cidx.to_string(), (chrom.clone(), *s, *e, locus));
                        }
                    }
                }
                // Task 7 fix (reproduction-gate finding): `bench/missing_copy_flag_pass.py`'s `truth` dict gates
                // EVERY read entering the detector on its PRIMARY alignment record physically overlapping
                // at least one of this family's own copy spans (lines 83-96: `for i, r in cp.items(): for a
                // in BAM.fetch(c, s, e): ... if a.flag & 2308: continue ...`). `catalog_copy_idx`
                // (`assignment.best_copy` here) is O2's INFERRED origin call, which for a MAPQ-0
                // multimapper can point to a copy the read's primary never physically touched -- Python
                // only trusts a read as "this copy's own abandoned reads" evidence when the primary
                // actually landed there. Measured without this gate: MCL106 copy 0 tested 25 origin-rejected
                // reads against Python's 9 (all 16 extra had their primary elsewhere), and genome-wide on
                // the 76-family sweep_v13 substrate the omission inflated 173 Python-tested pairs to 272
                // and 50 Python flags to 144, including `structural`-class flags Python's own detector
                // never produces on this substrate.
                //
                // `truth_copy[name] = (cidx, overlap)`: the candidate copy whose span the read's PRIMARY
                // overlaps MOST (ties keep the first-seen candidate, matching Python's strict `>` compare
                // over `cp.items()`'s insertion order). Built once per family from `bam_reads` (already
                // resident for the region) instead of a second BAM fetch. Extracted as `best_overlap_truth_copy`
                // (below) so the tie-break is directly unit-testable.
                let truth_copy = best_overlap_truth_copy(&bam_reads, &fa.copy_spans, &fa.copy_tids, catalog_index.as_ref());
                // Group this family's bam_reads by best-candidate catalog_copy_idx, split into rejected
                // (origin_rejected==true) and accepted (this family's own certificate-passed reads at that
                // copy), NOW bucketed by (cf, cidx) rather than bare cidx (Fix 1 above -- `cf` here is the
                // TRUE catalog family of `assignment.best_copy`'s own tid, captured instead of discarded).
                // `Assignment` (src/rustle/vg_family/copy_assign.rs:101-146) carries `best_copy: usize`
                // (an index into `copy_tids`/`copy_spans`, the same namespace `ci` uses above -- `.get()`,
                // not direct indexing, since no invariant here guarantees every family's assignments stay
                // in range), `status: AssignStatus` (Assigned/Ambiguous/Tied) and `origin_rejected: bool`.
                let mut rejected_by_cf: std::collections::HashMap<String, std::collections::HashMap<String, Vec<(String, Vec<u8>)>>> =
                    std::collections::HashMap::new();
                let mut accepted_by_cf: std::collections::HashMap<String, std::collections::HashMap<String, Vec<(String, Vec<u8>)>>> =
                    std::collections::HashMap::new();
                for &(read_i, ref assignment) in &fa.assignments {
                    let Some(br) = bam_reads.get(read_i) else { continue };
                    let Some(tid) = fa.copy_tids.get(assignment.best_copy) else { continue };
                    let Some((cf, cidx)) = catalog_index.as_ref().and_then(|ix| ix.get(tid)) else { continue };
                    let cidx = cidx.to_string();
                    // Python's `n in truth` gate: rejected reads need a primary overlapping ANY of this
                    // family's copies (grouped by O2's own best-copy call, which may differ from the
                    // truth-copy); control/accepted reads need the primary's OWN best-overlap copy to BE
                    // the candidate under test (`truth[n][0] == y` in the Python) AND not be rejected --
                    // that is `ctl_names`'s full condition (`bench/missing_copy_flag_pass.py:109`): `truth[n][0] ==
                    // y and A[n]['origin_rejected'] != '1' and A[n]['catalog_copy_idx'] == y`. Note there is
                    // NO status/verdict condition -- Python's own inline comment says so explicitly: "any
                    // MAPQ, any status: NPIP copies have few MAPQ-60 reads". A prior version of this code
                    // additionally required `assignment.status == AssignStatus::Assigned`, which is NOT in
                    // Python and is stricter: on `sweep_v13` it collapsed MCL1_073242 copy 15's true
                    // 31-read control pool to 1 read (forcing `Untestable` instead of Python's real
                    // verdict) and affected 8 (family, copy) groups genome-wide (docs/o1_ledger.md §6ib).
                    // Follow-up fix (final whole-branch review, round 2): compare BOTH `truth_cf` and
                    // `truth_cidx` against this assignment's own `(cf, cidx)` -- bare `cidx` equality alone
                    // (the pre-fix comparison) can succeed for the WRONG reason when two catalog copies from
                    // DIFFERENT families share a bare index inside the same `fa` (see `best_overlap_truth_copy`'s
                    // own doc comment).
                    let Some(((truth_cf, truth_cidx), _)) = truth_copy.get(br.name.as_str()) else { continue };
                    let entry = (br.name.clone(), br.read.seq.clone());
                    if assignment.origin_rejected {
                        rejected_by_cf.entry(cf.clone()).or_default().entry(cidx).or_default().push(entry);
                    } else if truth_cf == cf && *truth_cidx == cidx {
                        accepted_by_cf.entry(cf.clone()).or_default().entry(cidx).or_default().push(entry);
                    }
                }
                // Fix 1: one `detect_missing_copy_pairs` call PER real catalog family (`cf`) present among
                // this `fa`'s own copies, passing that TRUE `cf` as `family_id` (not `fa.family_id`, the
                // local co-located id) -- makes `RawPair.family_id` always equal `JoinRow.family_id`, closing
                // the namespace mismatch above. `cfs` sorted for run-to-run determinism (this diff's
                // BTreeMap-everywhere discipline elsewhere).
                let mut cfs: Vec<&String> = copy_span_by_cf.keys().collect();
                cfs.sort();
                let empty_rej: std::collections::HashMap<String, Vec<(String, Vec<u8>)>> = std::collections::HashMap::new();
                let empty_acc: std::collections::HashMap<String, Vec<(String, Vec<u8>)>> = std::collections::HashMap::new();
                for cf in cfs {
                    let spans = copy_span_by_cf.get(cf).unwrap();
                    let rej = rejected_by_cf.get(cf).unwrap_or(&empty_rej);
                    let acc = accepted_by_cf.get(cf).unwrap_or(&empty_acc);
                    // Fix 3 (Task 6, carried forward from Task 5's review): iterating a HashMap's `.keys()`
                    // is nondeterministic order -- sort by `copy_idx` so `o3_raw_pairs` (and therefore its
                    // `family_join.tsv`/`missing_copy_loci.tsv` row order) is stable run-to-run.
                    let mut inputs: Vec<rustle::vg_family::missing_copy_flag_pass::PairInput> = spans
                        .keys()
                        .map(|cidx| rustle::vg_family::missing_copy_flag_pass::PairInput {
                            copy_idx: cidx.clone(),
                            // CatalogCopy::partner is not threaded through FamilyAssignment yet -- default
                            // false never OVER-claims a partner exclusion (see the design doc's is_partner note).
                            is_partner: false,
                            rejected: rej.get(cidx).cloned().unwrap_or_default(),
                            accepted: acc.get(cidx).cloned().unwrap_or_default(),
                        })
                        .collect();
                    inputs.sort_by(|a, b| a.copy_idx.cmp(&b.copy_idx));
                    pairs.extend(rustle::vg_family::missing_copy_flag_pass::detect_missing_copy_pairs(
                        cf, spans, &genome, &inputs, &o3_params,
                    ));
                }
                // Orphan-locus scan: bam_reads whose primary lands outside every unit of this family,
                // clustered by proximity (<=5kb gap, matching bench/missing_copy_flag_pass.py), classified via the
                // shared genome-wide indices built once above.
                let fam_units: Vec<(String, u64, u64)> = fa.copy_spans.clone();
                let mut outside: Vec<&BamRead> = bam_reads
                    .iter()
                    .filter(|br| !br.is_secondary && !br.is_supplementary)
                    .filter(|br| {
                        // Fix 3 (final whole-branch review): this used to be a genomic-SPAN overlap check
                        // (`ref_start < e && s < ref_end`), the exact bug class Task 7 already fixed at the
                        // pair-detector's truth-overlap gate via `block_overlap()` -- a read whose intron (an
                        // `N` CIGAR op) merely SPANS a unit without any aligned block actually landing inside
                        // it was wrongly counted "inside", undercounting orphans. Now uses the same
                        // aligned-block overlap `block_overlap()` (M/=/X runs only) the truth gate uses.
                        !fam_units.iter().any(|(c, s, e)| br.chrom == *c && block_overlap(&br.read, *s, *e) > 0)
                    })
                    // Fix 1 (Task 6, revised): restrict to reads that are demonstrably ORPHANED for `fa`
                    // specifically -- present in `fa`'s own assignments (i.e. `fa` actually considered this
                    // read as a candidate) with that assignment rejected or with no candidate copy at all.
                    // A read absent from `fa.assignments` is excluded either way, matching the Python: if
                    // `fa` never considered the read, that says nothing about whether it is an orphan FOR
                    // `fa`'s purposes. See the `fa_verdict` comment above for why this must not be pooled
                    // across `fams`.
                    .filter(|br| {
                        fa_verdict
                            .get(br.name.as_str())
                            .map_or(false, |&(rejected, n_cand)| rejected || n_cand == 0)
                    })
                    .collect();
                outside.sort_by(|a, b| (a.chrom.as_str(), a.read.ref_start).cmp(&(b.chrom.as_str(), b.read.ref_start)));
                // Follow-up fix (final whole-branch review, round 2): `classify_orphan_locus`'s exclusion
                // test also used to take `&fa.family_id` -- the LOCAL co-located group's own (possibly
                // arbitrary) label, not necessarily any real catalog family id. Since `fa` can bundle more
                // than one true catalog family (Fix 1's own finding), the correct "self" set for this
                // exclusion is EVERY true catalog family id present among `fa`'s own copies -- exactly
                // `copy_span_by_cf`'s key set, already built above.
                let own_family_ids: std::collections::HashSet<String> = copy_span_by_cf.keys().cloned().collect();
                // Minor (final whole-branch review, "n_orphans hardcoded 0" parked in Task 6): a TRUE orphan
                // (`n_candidates==0`, nowhere to go at all) is a strict subset of this cluster's members
                // (which also include merely-rejected-but-had-a-candidate reads, per the `fa_verdict` filter
                // above) -- `fa_verdict` already carries `n_candidates` per read name, matching
                // `bench/missing_copy_flag_pass.py`'s own `n_orph = sum(1 for a in reads if
                // A[a.query_name].get('n_candidates','1')=='0')`.
                let mut clusters: Vec<(String, u64, u64, usize, usize)> = Vec::new();
                for br in &outside {
                    let end = read_ref_end_local(&br.read);
                    let is_true_orphan = fa_verdict.get(br.name.as_str()).map_or(false, |&(_, n_cand)| n_cand == 0);
                    if let Some(last) = clusters.last_mut() {
                        if last.0 == br.chrom && br.read.ref_start.saturating_sub(last.2) <= 5000 {
                            last.2 = last.2.max(end);
                            last.3 += 1;
                            if is_true_orphan {
                                last.4 += 1;
                            }
                            continue;
                        }
                    }
                    clusters.push((br.chrom.clone(), br.read.ref_start, end, 1, if is_true_orphan { 1 } else { 0 }));
                }
                for (chrom, start, end, n_reads, n_orphans) in clusters {
                    if n_reads < 3 {
                        continue;
                    }
                    let (class, n_genes, other_units) = rustle::vg_family::missing_copy_flag_pass::classify_orphan_locus(
                        &chrom, start, end, &own_family_ids, &o3_all_units_by_chrom, &o3_genes_by_chrom,
                    );
                    loci.push(rustle::vg_family::missing_copy_flag_pass::OrphanLocus {
                        chrom, start, end, n_reads, n_orphans, class,
                        n_genes_overlapping: n_genes, other_family_units: other_units,
                    });
                }
            }
            (pairs, loci)
        } else {
            (Vec::new(), Vec::new())
        };
        // Read-seeded copy discovery (opt-in, --discover-copies): cluster AS-tied reads' out-of-catalog
        // placements into candidate new copies. Gated the same way as the O3 block above -- empty Vec, no
        // allocation, when the flag is unset.
        let discovered: Vec<rustle::vg_family::copy_discovery::DiscoveredCopy> = if args.discover_copies {
            // The region's AS-tied reads are extracted ONCE; `discover_copies_for_family` then restricts
            // them, per family, to the reads that family actually considered (`fa.assignments`) before
            // clustering. Pooling them across families is the cross-family attribution bug the final
            // whole-branch review caught -- see that function's own doc comment.
            let tied = rustle::vg_family::copy_discovery::tie_partner_placements(&bam_reads);
            fams.iter().flat_map(|fa| discover_copies_for_family(fa, &bam_reads, &tied)).collect()
        } else {
            Vec::new()
        };
        Ok(RegionWork { contig: contig.clone(), lo, hi, read_names, read_chrom, read_mapqs, read_spans, read_blocks, read_strand, as_ev, n_mapped, fams, fallback, dna_needs, linearize_certs, transcripts, uniq_reads, o3_raw_pairs, o3_orphan_loci, discovered, union })
    };
    // Compute all regions (out-of-order across contigs when region_threads > 1), collected in the flat order.
    let works: Vec<RegionWork> = match &region_pool {
        Some(pool) => pool.install(|| {
            flat.par_iter().map(|(c, lo, hi)| compute(c, *lo, *hi)).collect::<Result<Vec<_>>>()
        })?,
        None => flat.iter().map(|(c, lo, hi)| compute(c, *lo, *hi)).collect::<Result<Vec<_>>>()?,
    };
    // PASS 1 (read-only): cross-family reconciliation. It must run BEFORE the drain, not as a post-pass
    // over `assign_rows`, because a molecule's status is emitted from `fa.assignments` at FOUR sites
    // inside the drain (`.assignments.tsv`, `.posterior.tsv`, `.psv_reads.tsv`, `.phased_reads.tsv`) plus
    // the two `--phase` graph builders' `Assigned` filters. A post-pass would fix only the first and
    // leave the rest disagreeing with it — `status_consistency_across_outputs` is the regression that
    // catches exactly that. Under `off` this does not run at all.
    let (xfam_conflicts, xfam_demote) = if xfam_mode == XfamMode::Off {
        (Vec::new(), std::collections::BTreeSet::new())
    } else {
        xfam_pass1(&works, region_families.is_some())
    };
    // The status a row is EMITTED with. Under `off`/`report` (and for every non-demoted key) this is
    // `a.status` verbatim, so those arms are byte-identical BY CONSTRUCTION. The `Abstain` test comes
    // FIRST so the other two modes never even build the lookup key.
    let eff_astatus = |read_name: &str, g: usize, f: usize, a: &rustle::vg_family::copy_assign::Assignment| -> AssignStatus {
        if xfam_mode == XfamMode::Abstain
            && matches!(a.status, AssignStatus::Assigned)
            && xfam_demote.contains(&(read_name.to_string(), g, f))
        {
            // Ambiguous, NEVER Tied: `Tied` is reserved for the K=0 identifiability wall
            // (`min_p >= alpha/(n-1)`), and the warning below exists precisely so catalog artifacts do
            // not masquerade as that wall. The intra-family record contradiction already demotes to
            // Ambiguous; the cross-family scope fix inherits the same bucket.
            AssignStatus::Ambiguous
        } else {
            a.status
        }
    };
    let eff_status = |read_name: &str, g: usize, f: usize, a: &rustle::vg_family::copy_assign::Assignment| -> &'static str {
        status_str(eff_astatus(read_name, g, f, a))
    };
    // SERIAL drain (PASS 2) in the original region order — every row push + the `gfam` id counter is
    // exactly the serial path, so the output is byte-identical.
    {
        for (gwork, work) in works.into_iter().enumerate() {
            let RegionWork { contig, lo, hi, read_names, read_chrom: _, read_mapqs, read_spans, read_blocks, read_strand, as_ev, n_mapped, fams, fallback, dna_needs, linearize_certs, transcripts, uniq_reads, o3_raw_pairs, o3_orphan_loci, discovered, union } = work;
            // O3 Phase 2 (Task 6): fold this region's raw pair stats + orphan loci into the genome-wide
            // vectors. Nothing is written here -- the Bonferroni threshold in `finalize_flags` needs every
            // region's pairs first, so `family_join.tsv`/`missing_copy_loci.tsv` are written once, after
            // this whole drain loop finishes.
            o3_all_raw_pairs.extend(o3_raw_pairs);
            o3_all_orphan_loci.extend(o3_orphan_loci);
            all_discovered.extend(discovered);
            union_all.absorb(union);
            let contig = &contig;
            let bam_reads = &read_names; // output stage indexes read NAMES (sequences were dropped)
            fallback_all.extend(fallback);
            dna_needs_rows.extend(dna_needs);
            linearize_certs_all.extend(linearize_certs);
            // Best mapq over each MOLECULE's records in this region (only its primary record can be >0).
            // See the `mqs` construction below: the tie-break invariance certificate is a molecule property.
            let mol_mapq: std::collections::HashMap<&str, u8> = {
                let mut m: std::collections::HashMap<&str, u8> = std::collections::HashMap::new();
                for (n, &q) in bam_reads.iter().zip(read_mapqs.iter()) {
                    let e = m.entry(n.as_str()).or_insert(0);
                    *e = (*e).max(q);
                }
                m
            };
            // ⭐ §6fq (user, 2026-09-06): the certificate machinery is for the CONTESTED molecules — the ones the
            // aligner could not place (primary MAPQ < 60). An uncontested molecule is assigned to its placement
            // (the copy its primary's blocks overlap most), as any assembler would use it; the certificate is
            // still computed for it and reported (`origin_rejected`), never applied. One sensitivity over every
            // read; abstention only among the contested. `--no-placement-assign` = the machinery on every read.
            let readthroughs = rustle::vg_family::copy_assign_pipeline::take_readthroughs();
            // ⭐ §6gz: under the AS-tied gate every molecule that reaches this point is tied by ALIGNMENT SCORE,
            // and a MAPQ of 60 is the aligner's chaining-stage opinion, not a guarantee — one human read
            // carried a MAPQ-60 primary at AS 1323 with three secondaries at AS 1384, and placement put it at
            // the primary's copy after the certificate had rejected every candidate. No tied molecule is ever
            // placed by its primary; placement exists only on the escape path.
            let placement_assign = args.molecule_observations && !args.no_molecule_observations && !args.no_placement_assign && args.no_as_tied_only;
            let mut placement_assigned = 0usize;
            let mut fams = fams;
            // the molecule's PRIMARY record (its highest-MAPQ record): the row's `ri` is the read-star
            // representative, which can be a secondary record at another copy
            let mol_primary: std::collections::HashMap<&str, usize> = {
                let mut m: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
                for (i, n) in bam_reads.iter().enumerate() {
                    let e = m.entry(n.as_str()).or_insert(i);
                    if read_mapqs[i] > read_mapqs[*e] {
                        *e = i;
                    }
                }
                m
            };
            if placement_assign {
                for fa in fams.iter_mut() {
                    for (ri, a) in fa.assignments.iter_mut() {
                        let mq = mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]);
                        if mq < 60 || rustle::vg_family::copy_assign_pipeline::is_tie_outside(&bam_reads[*ri]) {
                            continue; // §6gz: a competitor O2 never scored forbids placement too
                        }
                        // certified first: a molecule the machinery already assigned keeps that call (it can
                        // correct a placement: 4 % of MAPQ-60 simulated reads sit at the wrong copy, §6fq);
                        // the placement is the fallback when the machinery abstains or ties
                        if !args.placement_first && a.status == rustle::vg_family::copy_assign::AssignStatus::Assigned {
                            continue;
                        }
                        let pri = mol_primary.get(bam_reads[*ri].as_str()).copied().unwrap_or(*ri);
                        let Some(bl) = read_blocks.get(pri) else { continue };
                        let mut best: Option<(usize, u64)> = None;
                        for (ci, (c, s0, e0)) in fa.copy_spans.iter().enumerate() {
                            if c != contig || fa.copy_tids.get(ci).map_or(false, |t| rustle::vg_family::copy_assign_pipeline::is_partner(t)) {
                                continue; // §6ft: never place a molecule at a partner
                            }
                            let o: u64 = bl.iter().map(|&(bs, be)| be.min(*e0).saturating_sub(bs.max(*s0))).sum();
                            if o > 0 && best.map_or(true, |(_, bo)| o > bo) {
                                best = Some((ci, o));
                            }
                        }
                        let Some((pc, _)) = best else { continue };
                        if std::env::var_os("RUSTLE_STAR_DEBUG").is_some() && placement_assigned < 5 {
                            eprintln!("[placement] read {} mapq {mq} blocks {:?} -> copy {pc} span {:?} (was best_copy {} status {:?}); spans {:?}", bam_reads[*ri], &bl[..bl.len().min(3)], fa.copy_spans.get(pc), a.best_copy, a.status, &fa.copy_spans[..fa.copy_spans.len().min(3)]);
                        }
                        a.status = rustle::vg_family::copy_assign::AssignStatus::Assigned;
                        a.best_copy = pc;
                        a.resolvable = true;
                        let mut one = vec![0.0f64; fa.copy_spans.len()];
                        one[pc] = 1.0;
                        a.posterior = one;
                        placement_assigned += 1;
                    }
                }
            }
            placement_assigned_total += placement_assigned;
            let fams = fams;
            // --gtf: gene_tid (a copy's own locus) -> (family id, copy index), filled as fids are assigned below.
            let mut copy_gene: std::collections::HashMap<String, (String, usize)> = std::collections::HashMap::new();
            for (fwork, fa) in fams.iter().enumerate() {
                // JOIN KEY. Without `--families` this binary invents its own id (`CAFAM{i}`), which is
                // precisely why the O1 and O2 tables could not be joined. With `--families` the family
                // KEEPS the catalog's own `GWFAM{i}` — no id is minted, so `family_id` means the same
                // thing in both tables. `gfam` still advances so the two modes cannot alias.
                let fid = if region_families.is_some() { fa.family_id.clone() } else { format!("CAFAM{gfam}") };
                gfam += 1;
                if args.gtf {
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        copy_gene.insert(tid.clone(), (fid.clone(), ci));
                    }
                }
                // ⭐ register 734: the molecules that BELONG here — those with a PRIMARY (non-secondary,
                // non-supplementary) alignment overlapping a copy of THIS family. `in_copy` is not enough:
                // it fires on any aligned block, so a genome-wide multimapper visiting as a secondary counts
                // as a family read and inflates every per-family rate. On DAZ that inflated the denominator
                // 9.4x (18,192 rows, 1,935 of them local) and turned 34.7% assigned into a reported 3.7%.
                let primary_local: std::collections::HashSet<&str> = bam_reads
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| read_spans.get(*i).map_or(false, |sp| sp.2 == 0))
                    .filter(|(i, _)| {
                        read_blocks.get(*i).map_or(false, |bl| {
                            bl.iter().any(|&(bs, be)| fa.copy_spans.iter().any(|(c, s0, e0)| c == contig && be > *s0 && bs < *e0))
                        })
                    })
                    .map(|(_, n)| n.as_str())
                    .collect();
                let cat_idx_of = |ci: usize| -> String {
                    match (&catalog_index, fa.copy_tids.get(ci)) {
                        (Some(ix), Some(tid)) => ix.get(tid).map(|(_, i)| i.to_string()).unwrap_or_else(|| "NA".into()),
                        _ => "NA".into(),
                    }
                };
                for (ri, a) in &fa.assignments {
                    let in_copy = read_blocks.get(*ri).map_or(false, |bl| {
                        bl.iter().any(|&(bs, be)| fa.copy_spans.iter().any(|(c, s0, e0)| c == contig && be > *s0 && bs < *e0))
                    });
                    assign_rows.push(AssignRow {
                        read_name: bam_reads[*ri].clone(),
                        family_id: fid.clone(),
                        assigned_copy: a.best_copy,
                        status: eff_status(&bam_reads[*ri], gwork, fwork, a),
                        n_decisive: a.n_decisive,
                        margin: a.log_lr_margin,
                        p_value: a.p_value,
                        min_p_value: a.min_p_value,
                        as_ev: as_ev[*ri],
                        junction_conflict: a.junction_conflict,
                        origin_rejected: a.origin_rejected,
                        n_candidates: a.n_candidates,
                        in_copy,
                        primary_local: primary_local.contains(bam_reads[*ri].as_str()),
                        contested: mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]) < 60,
                        readthrough_into: match readthroughs.get(bam_reads[*ri].as_str()) {
                            Some(&(pc, _)) if pc == usize::MAX => "cut".to_string(),
                            Some(&(pc, _)) => cat_idx_of(pc),
                            None => "-".to_string(),
                        },
                        catalog_copy_idx: cat_idx_of(a.best_copy),
                        sibling_identity: a.sibling_identity,
                        n_cols_vs_nearest_sibling: a.n_cols_vs_nearest_sibling,
                    });
                }
                // --vg-realign (report-only): the re-align supplement's per-read decisions for this family.
                // Empty unless --vg-realign was passed (cfg.vg_realign gates run_family_realign itself).
                for r in &fa.realign_records {
                    vg_realign_lines.push(format!(
                        "{}\t{}\t{}\t{}\t{:.6}\t{}",
                        r.read_name, fid, r.action, r.target_copy, r.id_best, r.linear_copy
                    ));
                }
                // soft per-copy POSTERIOR + consistent ZONE (opt-in): localize even the unassignable reads.
                if args.posterior {
                    const FLOOR: f64 = 0.01; // a copy is in the consistent ZONE if its posterior exceeds this
                    let loci: Vec<Option<(String, u64)>> =
                        fa.copy_tids.iter().map(|t| parse_locus(t)).collect();
                    for (ri, a) in &fa.assignments {
                        if a.posterior.len() != fa.copy_tids.len() {
                            continue; // posterior frame must line up with the copy roster (e.g. post-freeze)
                        }
                        // apply the prior (uniform = posterior as-is; else weight by EM abundance), renormalize.
                        let mut post: Vec<f64> = a.posterior.clone();
                        if prior_abundance {
                            for (c, x) in post.iter_mut().enumerate() {
                                *x *= fa.copy_abundance.get(c).copied().unwrap_or(0.0).max(1e-9);
                            }
                            let z: f64 = post.iter().sum();
                            if z > 0.0 {
                                for x in &mut post {
                                    *x /= z;
                                }
                            }
                        }
                        // consistent zone = copies above the floor; its genomic extent + the posterior string.
                        let mut idx: Vec<usize> = (0..post.len()).filter(|&c| post[c] > FLOOR).collect();
                        idx.sort_by(|&a2, &b2| post[b2].partial_cmp(&post[a2]).unwrap());
                        let zone: Vec<u64> = idx.iter().filter_map(|&c| loci[c].as_ref().map(|l| l.1)).collect();
                        let chrom = idx
                            .iter()
                            .find_map(|&c| loci[c].as_ref().map(|l| l.0.clone()))
                            .unwrap_or_default();
                        let (zs, ze) = (zone.iter().min().copied().unwrap_or(0), zone.iter().max().copied().unwrap_or(0));
                        let pstr = idx
                            .iter()
                            .map(|&c| format!("{}:{:.3}", c, post[c]))
                            .collect::<Vec<_>>()
                            .join(",");
                        posterior_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            bam_reads[*ri], fid, eff_status(&bam_reads[*ri], gwork, fwork, a), idx.len(), chrom, zs, ze, pstr
                        ));
                    }
                }
                // soft per-copy abundance (+ the hard read count for comparison) + the tie-break invariance
                // certificate (anchored = unique-mapper support; see anchored_support).
                let bcs: Vec<usize> = fa.assignments.iter().map(|(_, a)| a.best_copy).collect();
                // MOLECULE-level mapq: only a molecule's PRIMARY record carries mapq>0, and after the
                // record->molecule reduction the surviving record need not be the primary. Max over the
                // molecule's records keeps `anchored`/`tie_invariant` counting the same molecules as before.
                let mqs: Vec<u8> = fa
                    .assignments
                    .iter()
                    .map(|(ri, _)| mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]))
                    .collect();
                for (ci, tid) in fa.copy_tids.iter().enumerate() {
                    let anchored = anchored_support(&bcs, &mqs, ci);
                    quant_rows.push(QuantRow {
                        family_id: fid.clone(),
                        copy_index: ci,
                        copy_tid: tid.clone(),
                        copy_chrom: fa.copy_spans.get(ci).map(|s| s.0.clone()).unwrap_or_default(),
                        copy_start: fa.copy_spans.get(ci).map_or(0, |s| s.1),
                        copy_end: fa.copy_spans.get(ci).map_or(0, |s| s.2),
                        abundance: fa.copy_abundance.get(ci).copied().unwrap_or(0.0),
                        ci: fa.copy_abundance_ci.get(ci).copied().unwrap_or(0.0),
                        // an ORPHAN (no candidate, §6fg) is nobody's hard read (it carried copy 0's index by default)
                        n_hard: fa.assignments.iter().filter(|(_, a)| a.best_copy == ci && !(a.origin_rejected && a.n_candidates == 0)).count(),
                        n_soft: fa.assignments.iter().filter(|(_, a)| !a.origin_rejected).map(|(_, a)| a.posterior.get(ci).copied().unwrap_or(0.0)).sum(),
                        anchored,
                        tie_invariant: anchored as u32 >= GATE_MIN_READS,
                        junction_invariant: fa.copy_junction_support.get(ci).copied().unwrap_or(0) as u32
                            >= GATE_MIN_READS,
                    });
                    // `--families`: name the catalog row this copy IS. Looked up by the catalog `tid`
                    // (carried through `DenovoTranscript::tid` untouched), never by position, so the join
                    // does not depend on the copy ordering surviving the assignment stage.
                    if let Some(ix) = &catalog_index {
                        let (cf, cidx) = match ix.get(tid) {
                            Some((f, i)) => (f.clone(), i.to_string()),
                            // Unreachable while the roster-changing legs are refused; reported rather than
                            // hidden, because a copy O2 invented is exactly what --families must never do.
                            None => ("NOT_IN_CATALOG".to_string(), "NA".to_string()),
                        };
                        join_rows.push(JoinRow {
                            line: format!(
                                "{fid}\t{ci}\t{tid}\t{cf}\t{cidx}\t{}\t{}\t{}\t{}",
                                fa.copy_spans.get(ci).map(|s| s.0.clone()).unwrap_or_default(),
                                fa.copy_spans.get(ci).map_or(0, |s| s.1),
                                fa.copy_spans.get(ci).map_or(0, |s| s.2),
                                fa.assignments.iter().filter(|(_, a)| a.best_copy == ci).count(),
                            ),
                            family_id: cf.clone(),
                            copy_idx: cidx.clone(),
                        });
                    }
                }
                // gene-conversion: report per-read candidate switches (RT-switch-like) vs recurrence-confirmed
                // events, so the discriminator's anti-artifact gate is visible on real data.
                eprintln!(
                    "[mosaic] {fid}: {} reads showed a candidate copy-switch -> {} recurrence-confirmed conversion event(s)",
                    fa.mosaic_reads,
                    fa.conversions.iter().filter(|e| e.confirmed).count()
                );
                // gene-conversion events (per-molecule PSV-path switches that recur across reads)
                for ev in &fa.conversions {
                    mosaic_rows.push(MosaicRow {
                        family_id: fid.clone(),
                        copy_a: ev.copy_a,
                        copy_b: ev.copy_b,
                        bp_lo: ev.breakpoint_ref.0,
                        bp_hi: ev.breakpoint_ref.1,
                        n_reads: ev.n_supporting_reads,
                        dispersion: ev.breakpoint_dispersion,
                        confirmed: ev.confirmed,
                    });
                }
                // copy-level historical conversions (a copy whose PSV vector is a mosaic of two others)
                for cv in &fa.copy_conversions {
                    copyconv_rows.push(CopyConvRow {
                        family_id: fid.clone(),
                        copy_c: fa.copy_tids.get(cv.copy_c).cloned().unwrap_or_else(|| cv.copy_c.to_string()),
                        copy_a: cv.copy_a,
                        copy_b: cv.copy_b,
                        bp_lo: cv.breakpoint.0,
                        bp_hi: cv.breakpoint.1,
                        n_decisive: cv.n_decisive,
                    });
                }
                // raw per-molecule PSV genotype evidence (the assignment-proof matrix)
                if args.dump_psv {
                    let allele_str = |v: &Vec<Option<u8>>| -> String {
                        v.iter().map(|o| o.map(|b| b as char).unwrap_or('.')).collect()
                    };
                    for ((ri, a), obs) in fa.assignments.iter().zip(fa.read_psv_obs.iter()) {
                        psv_read_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
                            bam_reads[*ri], fid, a.best_copy, eff_status(&bam_reads[*ri], gwork, fwork, a),
                            a.log_lr_margin, a.n_decisive, allele_str(obs)
                        ));
                    }
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let alleles = fa.copy_psv_alleles.get(ci).map(allele_str).unwrap_or_default();
                        let strand = fa.copy_strand.get(ci).copied().unwrap_or('+');
                        psv_copy_lines.push(format!("{}\t{}\t{}\t{}\t{}", fid, ci, tid, alleles, strand));
                    }
                    for (col, pos) in fa.psv_col_pos.iter().enumerate() {
                        psv_col_lines.push(format!("{}\t{}\t{}", fid, col, pos.map(|x| x as i64).unwrap_or(-1)));
                    }
                }
                // EM soft-relaxation (opt-in): re-runs the family's PSV evidence through the maximum-
                // likelihood EM engine (Task 1's exact gate likelihood, Task 4's binary wiring) for a soft
                // posterior + recovered abundance, alongside (not instead of) the hard PSV+junction
                // assignment above. Fully gated: with `--em` absent this block never runs, so the hard
                // outputs are untouched.
                if args.em {
                    let em_result = em_assign_family(
                        &fa.read_psv_obs,
                        &fa.copy_psv_alleles,
                        &fa.read_junctions,
                        &fa.copy_junctions,
                        &params,
                        args.em_eps,
                        args.em_max_iter,
                    );
                    for (row_idx, (ri, _)) in fa.assignments.iter().enumerate() {
                        if row_idx >= em_result.posteriors.len() {
                            continue; // posterior frame must line up with the read roster (e.g. post-freeze)
                        }
                        let post = &em_result.posteriors[row_idx];
                        let argmax = post
                            .iter()
                            .enumerate()
                            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                            .map(|(k, _)| k)
                            .unwrap_or(0);
                        let label = match em_result.labels[row_idx] {
                            rustle::vg_family::em_copy_assign::EmLabel::Certified => "Certified",
                            rustle::vg_family::em_copy_assign::EmLabel::SoftZone => "SoftZone",
                        };
                        let post_str = post
                            .iter()
                            .enumerate()
                            .map(|(k, p)| format!("{}:{:.4}", k, p))
                            .collect::<Vec<_>>()
                            .join(";");
                        em_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            bam_reads[*ri], fid, argmax, label, post_str, em_result.n_iter
                        ));
                    }
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let n_reads_soft: f64 = em_result.posteriors.iter().map(|p| p.get(ci).copied().unwrap_or(0.0)).sum();
                        em_abundance_lines.push(format!(
                            "{}\t{}\t{:.4}\t{:.2}",
                            fid,
                            tid,
                            em_result.abundances.get(ci).copied().unwrap_or(0.0),
                            n_reads_soft
                        ));
                    }
                }
                // FACULTATIVE phasing: phase set = family; haplotypes = its copies; read->haplotype =
                // the PSV assignment. A read is phased iff it clears the decisive-margin gate (Assigned);
                // Ambiguous/Tied reads are emitted with haplotype = -1 (unphaseable = K-frontier).
                if args.phase {
                    // Effective per-row status, built ONCE so the phase-block counters, the per-copy
                    // support counts, the read haplotypes and the two graph builders' `Assigned` filters
                    // all agree row-for-row with `.assignments.tsv`. Element-wise equal to `a.status`
                    // under `off`/`report`.
                    let eff: Vec<AssignStatus> = fa
                        .assignments
                        .iter()
                        .map(|(ri, a)| eff_astatus(&bam_reads[*ri], gwork, fwork, a))
                        .collect();
                    let n_phased = eff.iter().filter(|s| matches!(s, AssignStatus::Assigned)).count();
                    phase_block_lines.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        fid, fa.chrom, fa.n_copies, fa.psv_cols, n_phased,
                        fa.assignments.len() - n_phased
                    ));
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let mut vs = String::new();
                        if let Some(alleles) = fa.copy_psv_alleles.get(ci) {
                            for (col, a) in alleles.iter().enumerate() {
                                if let (Some(b), Some(Some(pos))) = (a, fa.psv_col_pos.get(col)) {
                                    if !vs.is_empty() {
                                        vs.push(';');
                                    }
                                    vs.push_str(&format!("{}:{}", pos, *b as char));
                                }
                            }
                        }
                        let n_sup = fa
                            .assignments
                            .iter()
                            .zip(eff.iter())
                            .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                            .count();
                        phased_hap_lines.push(format!("{}\t{}\t{}\t{}\t{}", fid, ci, tid, n_sup, vs));
                    }
                    for ((ri, a), es) in fa.assignments.iter().zip(eff.iter()) {
                        let hap: i64 = if matches!(es, AssignStatus::Assigned) {
                            a.best_copy as i64
                        } else {
                            -1
                        };
                        phased_read_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{:.3}\t{}",
                            bam_reads[*ri], fid, hap, a.n_decisive, a.log_lr_margin, status_str(*es)
                        ));
                    }

                    // Materialize this family's copy-graph (REFERENCE walk + tagged copy paths + read
                    // walks over the shared PSV-bubble nodes) and fold its GFA lines / Bandage colours /
                    // status legend into the region's accumulators. Replaces the inline emitter above.
                    let ref_base = |chrom: &str, pos: u64| {
                        genome_for(chrom).ok()
                            .and_then(|g| g.fetch_sequence(chrom, pos, pos + 1))
                            .and_then(|v| v.first().copied())
                    };
                    let cg = build_copy_graph(&fid, fa, ref_base, bam_reads, annotation.as_deref(), &eff);
                    let gl = cg.gfa_lines();
                    for s in gl.segs { gfa_segs.insert(s); }
                    for l in gl.links { gfa_links.insert(l); }
                    gfa_paths.extend(gl.paths);
                    gfa_walks.extend(gl.walks);
                    for row in cg.colours_csv().lines() { gfa_colors.push(row.to_string()); }
                    legend_rows.extend(cg.legend_tsv().lines().map(|s| s.to_string()));

                    // v2: this family's exon presence/absence graph (copies = walks over shared exon
                    // classes; a copy-specific exon reads as a visible arm). Only built when the family
                    // carries an intron chain (no `copy_introns` => no exon structure to graph).
                    // Sequence-free at this point; folded into `<out>.exon.gfa` at write time below, where
                    // `genome_for` fetches each exon's bases. Same `--gff` annotation overlay as v1.
                    if !fa.copy_introns.is_empty() {
                        exon_graphs.push(build_exon_graph(&fid, fa, annotation.as_deref(), &eff));
                    }
                }
                // reference-free per-family copy number (Task R1): chi_H (PSV conflict-structure
                // lower bound) always computed; depth_cn only when --lambda-global was given.
                famcn_rows.push(FamCnRow {
                    family_id: fid.clone(),
                    chrom: fa.chrom.clone(),
                    n_copies: fa.n_copies,
                    n_reads: fa.n_reads,
                    chi_h: chi_h_with_junctions(&fa.copy_psv_alleles, &fa.copy_junctions),
                    depth_cn: lambda.map(|lam| depth_cn(fa.n_reads, lam)).unwrap_or(f64::NAN),
                    regime: if fa.collapsed_copies > 0 { "reference_collapsed" } else { "reference_resolved" },
                });
                family_rows.push(FamilyRow {
                    family_id: fid,
                    chrom: fa.chrom.clone(),
                    n_copies: fa.n_copies,
                    n_reads: fa.n_reads,
                    psv_cols: fa.psv_cols,
                    resolvable_psv: fa.resolvable_psv,
                    resolvable_j: fa.resolvable_j,
                    junction_only: fa.junction_only,
                    assigned_j: fa.assigned_j,
                    uniq_agree: fa.uniq_agree,
                    uniq: fa.uniq,
                    collapsed_copies: fa.collapsed_copies,
                    rescued_copies: fa.rescued_copies,
                });
            }
            eprintln!("[copy_assign]   {contig}:{lo}-{hi}: {} mapped reads -> {} families", n_mapped, fams.len());
            // --gtf: emit every isoform of this region (transcript + exon rows), tagging family-copy genes.
            //
            // ⭐ §6gl. Three fixes over the first form, all here:
            //  (1) the copy tag matched the isoform's `gene_tid` against a map of catalog copy TIDs. Those are
            //      different namespaces — isoforms are `DN_<contig>_<pos>_<n>`, catalog copies
            //      `MCL_<contig>_<pos>` — so the lookup NEVER fired with `--families` and every transcript was
            //      emitted `multicopy "false"` (register 756). It is a positional match now.
            //  (2) the isoform now carries the COPY ITS OWN READS were assigned to, from the certificate, not
            //      merely the copy it sits inside.
            //  (3) `copy_status` distinguishes three cases a blank used to conflate: `assigned` (≥1 read
            //      certificate-assigned), `undecidable` (reads were adjudicated and none carried a
            //      certificate — the identifiability wall) and `unadjudicated` (no matched read reached the
            //      assignment at all, i.e. the isoform lies outside the swept copies). On NPIP those are
            //      108 / 328 / 376 of 817 — reporting them as one blank hid the difference between "we cannot
            //      tell" and "we did not look".
            let read_chain: Vec<Vec<(u64, u64)>> = read_blocks
                .iter()
                .map(|bl| bl.windows(2).map(|w| (w[0].1, w[1].0)).filter(|&(a, b)| b > a).collect())
                .collect();
            let verdict: std::collections::HashMap<&str, &AssignRow> =
                assign_rows.iter().map(|r| (r.read_name.as_str(), r)).collect();
            // ⚠ §6gl: `DenovoTranscript::tid` is `DN_<contig>_<start>_<n_exon>`, which COLLIDES — two distinct
            // isoforms sharing a start and an exon count get the same id. On NPIP 53 ids covered 122 of the
            // 886 transcript rows, and any consumer keyed on `transcript_id` (IGV, gffcompare, bedtools, our
            // own join) silently merges them into one impossible model. Disambiguated HERE, in the GTF only,
            // so no catalog's tids move; the assembler's own id scheme is left for a separate change.
            let mut tid_seen: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
            // B2/read-provenance (below) needs to look up "which gate-passed transcript did this AS-tied
            // RECORD contribute to", keyed by the record's OWN chain — filled alongside `uniq_tid` so the
            // answer is the disambiguated id, not the raw (colliding) `t.tid`. A spliced chain is matched by
            // exact intron equality (unambiguous); an EMPTY chain is not (every unspliced record on the
            // contig has one) so unspliced transcripts are looked up by span containment instead — the
            // fix for register 757, which the isoform-vote code below (`ch.is_empty() && ...`) already
            // applies but this record-classification path had not, and it had recurred there.
            let mut chain_uniq_tid: std::collections::HashMap<Vec<(u64, u64)>, String> = std::collections::HashMap::new();
            let mut unspliced_gate_passed: Vec<(u64, u64, String)> = Vec::new();
            let mut prod_genome: Option<std::sync::Arc<GenomeIndex>> = None;
            // --gtf-copy-set: held-back family transcripts, placed by evidence after the loop
            struct PendingTx<'a> { fw: usize, cidx: String, tline: String, elines: Vec<String>, t: &'a TranscriptRec, uniq: std::collections::BTreeMap<String, usize>, asg: std::collections::BTreeMap<String, usize>, abst: Vec<String>, uniq_tid: String }
            let mut pending: Vec<PendingTx> = Vec::new();
            // unique mappers by (chrom, intron chain) -> their primary spans (evidence for a copy, §6hn)
            let mut uniq_by_chain: std::collections::HashMap<(&str, &[(u64, u64)]), Vec<(u64, u64)>> = std::collections::HashMap::new();
            for (c, s0, e0, ch) in &uniq_reads {
                uniq_by_chain.entry((c.as_str(), ch.as_slice())).or_default().push((*s0, *e0));
            }
            let sweep_to_catalog = |fw: usize, ci: usize| -> String {
                match (&catalog_index, fams[fw].copy_tids.get(ci)) {
                    (Some(ix), Some(tid)) => ix.get(tid).map(|(_, i)| i.to_string()).unwrap_or_else(|| ci.to_string()),
                    _ => ci.to_string(),
                }
            };
            // ⭐ PREREG fd894558: per-locus read depth, keyed on the SAME `gene_tid` `collapse_loci_groups`
            // already assigned (threaded into every `TranscriptRec` unconditionally, not only under
            // --gtf-copy-set). SUM, not max: exact-chain collapse FRAGMENTS one locus's reads across many
            // near-identical chains (measured: one locus's 41 reads split 9 ways, largest chain only 12), so
            // no single chain ever reaches the locus's true depth the way a splice-graph "gene" would — the
            // max badly undercounts it (row 807). The sum is the closer proxy for StringTie's per-locus
            // coverage; it over-counts only where two OVERLAPPING junction-groups both touch this locus,
            // which the two-form/one-region sweep does not produce in practice (measured, not assumed).
            let mut group_total_reads: std::collections::HashMap<&str, u64> = std::collections::HashMap::new();
            for t in &transcripts {
                *group_total_reads.entry(t.gene_tid.as_str()).or_insert(0) += t.n_reads as u64;
            }
            // ⭐ PREREG 35e290a9: the boundary-outlier test, precomputed once per region by transcript INDEX
            // (a `DenovoTranscript::tid` can collide before the `uniq_tid` disambiguation below, so the key
            // must be the position in `transcripts`, not the tid string).
            let mut boundary_far: std::collections::HashMap<usize, (bool, bool, u64, u64)> = std::collections::HashMap::new();
            if args.min_boundary_fraction > 0.0 {
                let mut by_locus: std::collections::HashMap<&str, Vec<usize>> = std::collections::HashMap::new();
                for (i, t) in transcripts.iter().enumerate() {
                    by_locus.entry(t.gene_tid.as_str()).or_default().push(i);
                }
                let bucket = |x: u64| x / 50;
                for idxs in by_locus.values() {
                    if idxs.len() < 2 {
                        continue;
                    }
                    let total: u64 = idxs.iter().map(|&i| transcripts[i].n_reads as u64).sum();
                    // RIGHT (largest `end`): the farthest bucket vs. the next-farthest.
                    let mut by_end: std::collections::BTreeMap<u64, (u64, Vec<usize>)> = std::collections::BTreeMap::new();
                    for &i in idxs {
                        let e = by_end.entry(bucket(transcripts[i].end)).or_insert((0, Vec::new()));
                        e.0 += transcripts[i].n_reads as u64;
                        e.1.push(i);
                    }
                    if by_end.len() >= 2 {
                        let mut v: Vec<(u64, u64, Vec<usize>)> = by_end.into_iter().map(|(b, (r, idx))| (b * 50, r, idx)).collect();
                        v.sort_unstable_by_key(|&(pos, _, _)| pos);
                        let (far_pos, far_reads, far_idx) = v.pop().unwrap();
                        let second_pos = v.last().unwrap().0;
                        let gap = far_pos.saturating_sub(second_pos);
                        if gap > args.min_boundary_gap && (far_reads as f64 / total.max(1) as f64) < args.min_boundary_fraction {
                            for i in far_idx {
                                let e = boundary_far.entry(i).or_insert((false, false, 0, 0));
                                e.1 = true;
                                e.3 = gap;
                            }
                        }
                    }
                    // LEFT (smallest `start`): the farthest bucket vs. the next-farthest.
                    let mut by_start: std::collections::BTreeMap<u64, (u64, Vec<usize>)> = std::collections::BTreeMap::new();
                    for &i in idxs {
                        let e = by_start.entry(bucket(transcripts[i].start)).or_insert((0, Vec::new()));
                        e.0 += transcripts[i].n_reads as u64;
                        e.1.push(i);
                    }
                    if by_start.len() >= 2 {
                        let v: Vec<(u64, u64, Vec<usize>)> = by_start.into_iter().map(|(b, (r, idx))| (b * 50, r, idx)).collect();
                        let (far_pos, far_reads, far_idx) = v[0].clone();
                        let second_pos = v[1].0;
                        let gap = second_pos.saturating_sub(far_pos);
                        if gap > args.min_boundary_gap && (far_reads as f64 / total.max(1) as f64) < args.min_boundary_fraction {
                            for i in far_idx {
                                let e = boundary_far.entry(i).or_insert((false, false, 0, 0));
                                e.0 = true;
                                e.2 = gap;
                            }
                        }
                    }
                }
            }
            for (ti, t) in transcripts.iter().enumerate() {
                let n = tid_seen.entry(t.tid.as_str()).or_insert(0);
                *n += 1;
                let uniq_tid = if *n == 1 { t.tid.clone() } else { format!("{}.{}", t.tid, *n) };
                if t.introns.is_empty() {
                    unspliced_gate_passed.push((t.start, t.end, uniq_tid.clone()));
                } else {
                    chain_uniq_tid.entry(t.introns.clone()).or_insert_with(|| uniq_tid.clone());
                }
                let isoform_fraction = t.n_reads as f64 / (*group_total_reads.get(t.gene_tid.as_str()).unwrap_or(&1)).max(1) as f64;
                let depth_low = args.min_isoform_fraction > 0.0 && isoform_fraction < args.min_isoform_fraction;
                let (b_left, b_right, gap_left, gap_right) = boundary_far.get(&ti).copied().unwrap_or((false, false, 0, 0));
                let boundary_low = b_left || b_right;
                let low_confidence = depth_low || boundary_low;
                // (1) positional: the catalog copy this isoform overlaps most
                let best = fams.iter().enumerate().flat_map(|(fw, fa)| {
                    fa.copy_spans.iter().enumerate().map(move |(ci, (c, s0, e0))| (fw, ci, c, *s0, *e0))
                }).filter(|(_, _, c, s0, e0)| *c == &t.chrom && t.end > *s0 && t.start < *e0)
                  .max_by_key(|(_, _, _, s0, e0)| t.end.min(*e0).saturating_sub(t.start.max(*s0)));
                let (fam_attr, multicopy) = match best {
                    Some((fw, ci, _, _, _)) => {
                        let fid = if region_families.is_some() { fams[fw].family_id.clone() } else { String::new() };
                        // report the CATALOG copy index, the same namespace `assigned_copy` uses below — the
                        // sweep's own index differs (sweep 4 == catalog 24 on NPIP), and printing the two
                        // schemes on one line reads as a disagreement when they in fact agree
                        let idx = match (&catalog_index, fams[fw].copy_tids.get(ci)) {
                            (Some(ix), Some(tid)) => ix.get(tid).map(|(_, i)| i.to_string()).unwrap_or_else(|| ci.to_string()),
                            _ => ci.to_string(),
                        };
                        (format!(" family_id \"{fid}\"; copy_index \"{idx}\";"), "true")
                    }
                    None => (String::new(), "false"),
                };
                // (2)+(3) the copy its own reads were assigned to, and how sure that is
                let mut votes: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
                let (mut seen, mut matched) = (0usize, 0usize);
                let mut matched_ri: Vec<usize> = Vec::new();
                for (ri, ch) in read_chain.iter().enumerate() {
                    if read_spans.get(ri).map_or(true, |sp| sp.2 != 0) {
                        continue; // primaries only
                    }
                    let same = if t.introns.is_empty() {
                        // an unspliced isoform's chain is EMPTY and would otherwise match every unspliced read
                        // in the region (register 757, which recurred): require containment instead
                        ch.is_empty() && read_spans[ri].0 < t.end && t.start < read_spans[ri].1
                    } else {
                        ch.as_slice() == t.introns.as_slice()
                    };
                    if !same {
                        continue;
                    }
                    matched += 1;
                    if let Some(r) = verdict.get(bam_reads[ri].as_str()) {
                        seen += 1;
                        if r.status == "assigned" && !r.origin_rejected {
                            *votes.entry(r.catalog_copy_idx.as_str()).or_insert(0) += 1;
                        }
                    }
                    if (args.gtf_copy_set && !args.no_gtf_copy_set) {
                        matched_ri.push(ri);
                    }
                }
                let total: usize = votes.values().sum();
                let copy_attr = match votes.iter().max_by_key(|(_, n)| **n) {
                    Some((c, n)) => format!(
                        " assigned_copy \"{c}\"; copy_votes \"{n}/{total}\"; copy_purity \"{:.3}\"; copy_status \"assigned\";",
                        *n as f64 / total.max(1) as f64
                    ),
                    None if seen > 0 => format!(" copy_status \"undecidable\"; adjudicated_reads \"{seen}\";"),
                    None => format!(" copy_status \"unadjudicated\"; matched_reads \"{matched}\";"),
                };
                let fam_attr = format!("{fam_attr}{copy_attr}");
                let fam_attr = if args.min_isoform_fraction > 0.0 || args.min_boundary_fraction > 0.0 {
                    let reason = match (depth_low, boundary_low) {
                        (true, true) => "both",
                        (true, false) => "depth",
                        (false, true) => "boundary",
                        (false, false) => "none",
                    };
                    let mut extra = String::new();
                    if args.min_isoform_fraction > 0.0 {
                        extra.push_str(&format!(" isoform_fraction \"{isoform_fraction:.3}\";"));
                    }
                    if args.min_boundary_fraction > 0.0 {
                        extra.push_str(&format!(" boundary_gap_left \"{gap_left}\"; boundary_gap_right \"{gap_right}\";"));
                    }
                    format!("{fam_attr}{extra} low_confidence \"{low_confidence}\"; low_confidence_reason \"{reason}\";")
                } else {
                    fam_attr
                };
                // §6gp: longest ORF over the strand-oriented exon-sum. Sequence comes from the genome at the
                // transcript's own exons, so this needs no annotation and no CDS.
                let orf_aa = if args.productivity {
                    let gi = match prod_genome.as_ref() {
                        Some(g) => g.clone(),
                        None => {
                            let g = genome_for(&contig)?;
                            prod_genome = Some(g.clone());
                            g
                        }
                    };
                    let mut seq: Vec<u8> = Vec::new();
                    let mut prev = t.start;
                    let mut ex: Vec<(u64, u64)> = Vec::new();
                    for &(d, a) in &t.introns {
                        ex.push((prev, d));
                        prev = a;
                    }
                    ex.push((prev, t.end));
                    for &(es, ee) in &ex {
                        if let Some(part) = gi.fetch_sequence(&t.chrom, es, ee) {
                            seq.extend_from_slice(&part);
                        }
                    }
                    if t.strand == '-' {
                        seq.reverse();
                        for b in seq.iter_mut() {
                            *b = match *b {
                                b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
                                b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c', x => x,
                            };
                        }
                    }
                    Some(longest_orf(&seq) / 3)
                } else {
                    None
                };
                if let Some(aa) = orf_aa {
                    let fid = re_attr(&fam_attr, "family_id").unwrap_or_default();
                    let cp = re_attr(&fam_attr, "assigned_copy")
                        .or_else(|| re_attr(&fam_attr, "copy_index"))
                        .unwrap_or_default();
                    prod_rows.push((fid, cp, uniq_tid.clone(), aa));
                }
                let fam_attr = match orf_aa {
                    Some(aa) => format!("{fam_attr} orf_aa \"{aa}\";"),
                    None => fam_attr,
                };
                let gs = t.start + 1; // GTF is 1-based, end-inclusive (our coords are 0-based half-open)
                let tline = format!(
                    "{}\trustle\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; reads \"{}\"; multicopy \"{}\";{}",
                    t.chrom, gs, t.end, t.strand, t.gene_tid, uniq_tid, t.n_reads, multicopy, fam_attr
                );
                // exons = the gene span minus the introns (the read's spliced structure)
                let mut prev = t.start;
                let mut exons: Vec<(u64, u64)> = Vec::new();
                for &(d, a) in &t.introns {
                    exons.push((prev, d));
                    prev = a;
                }
                exons.push((prev, t.end));
                let elines: Vec<String> = exons.iter().enumerate().map(|(k, (es, ee))| format!(
                    "{}\trustle\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";",
                    t.chrom, es + 1, ee, t.strand, t.gene_tid, uniq_tid, k + 1
                )).collect();
                // --gtf-copy-set: multi-intron family transcripts are held back and placed by evidence below
                if (args.gtf_copy_set && !args.no_gtf_copy_set) && t.introns.len() >= 2 && best.is_some() && !low_confidence {
                    let (fw, ci, _, _, _) = best.unwrap();
                    let (mut uniq, mut asg, mut abst) = (std::collections::BTreeMap::new(), std::collections::BTreeMap::new(), Vec::new());
                    // unique mappers (gate-dropped primaries) with this chain: evidence at the copy their primary lies in
                    for &(s0, e0) in uniq_by_chain.get(&(t.chrom.as_str(), t.introns.as_slice())).map(|v| v.as_slice()).unwrap_or(&[]) {
                        if let Some(pc) = fams[fw].copy_spans.iter().position(|(c, a, b)| *c == t.chrom && s0 < *b && e0 > *a) {
                            *uniq.entry(sweep_to_catalog(fw, pc)).or_insert(0usize) += 1;
                        }
                    }
                    for &ri in &matched_ri {
                        match verdict.get(bam_reads[ri].as_str()) {
                            Some(r) if r.status == "assigned" && !r.origin_rejected => { *asg.entry(r.catalog_copy_idx.clone()).or_insert(0usize) += 1; }
                            Some(_) => abst.push(bam_reads[ri].clone()),
                            None => {} // an AS-tied read without a row (no catalog placement): no evidence
                        }
                    }
                    pending.push(PendingTx { fw, cidx: sweep_to_catalog(fw, ci), tline, elines, t, uniq, asg, abst, uniq_tid: uniq_tid.clone() });
                } else {
                    gtf_lines.push(tline);
                    gtf_lines.extend(elines);
                }
            }
            if (args.gtf_copy_set && !args.no_gtf_copy_set) && !pending.is_empty() {
                // ⭐ §6hn: group the held transcripts across copies by LIFT, then place each group by evidence.
                let tol = args.gtf_lift_tol;
                let gi = match prod_genome.as_ref() { Some(g) => g.clone(), None => { let g = genome_for(&contig)?; prod_genome = Some(g.clone()); g } };
                let mut lifts_by_fam: std::collections::HashMap<usize, std::collections::HashMap<(usize, usize), Vec<LiftBlocks>>> = std::collections::HashMap::new();
                let fam_ids: std::collections::BTreeSet<usize> = pending.iter().map(|p| p.fw).collect();
                for &fw in &fam_ids {
                    lifts_by_fam.insert(fw, copy_span_lifts(&fams[fw].copy_spans, &gi, &format!("{}.{}", args.out, fw)));
                }
                let strand_of: std::collections::HashMap<String, char> = region_families.as_ref().map(|rf| rf.values().flatten().flat_map(|f| f.copies.iter()).map(|c| (c.copy_idx.to_string(), c.strand)).collect()).unwrap_or_default();
                let sweep_ci = |fw: usize, cidx: &str| -> Option<usize> { (0..fams[fw].copy_spans.len()).find(|&ci| sweep_to_catalog(fw, ci) == cidx) };
                let lift_pos = |fw: usize, g: u64, a: usize, b: usize| -> Option<(u64, u64)> {
                    let (_, sa, _) = &fams[fw].copy_spans[a];
                    let (_, sb, _) = &fams[fw].copy_spans[b];
                    let rel = g.checked_sub(*sa)?;
                    lifts_by_fam.get(&fw)?.get(&(a, b))?.iter().filter_map(|l| l.map(rel)).min_by_key(|&(_, d)| d).map(|(t, d)| (t + sb, d))
                };
                fn find(p: &mut Vec<usize>, mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
                let n = pending.len();
                let mut parent: Vec<usize> = (0..n).collect();
                for i in 0..n {
                    for j in (i + 1)..n {
                        let (p, q) = (&pending[i], &pending[j]);
                        if p.fw != q.fw || p.cidx == q.cidx || p.t.strand != q.t.strand || p.t.introns.len() != q.t.introns.len() { continue; }
                        let (Some(a), Some(b)) = (sweep_ci(p.fw, &p.cidx), sweep_ci(q.fw, &q.cidx)) else { continue };
                        // same isoform: every boundary lifts onto the other's within `tol`. The lift may extrapolate
                        // past an aligned block's end (a last exon beyond the alignment): agreement of the
                        // extrapolated coordinate is what counts here (the measured rule, `bench/gtf_copy_set.py`);
                        // the distance-to-block guard applies to PLACEMENT lifts below, not to grouping.
                        let ok = p.t.introns.iter().zip(q.t.introns.iter()).all(|(&(x0, x1), &(y0, y1))| {
                            matches!((lift_pos(p.fw, x0, a, b), lift_pos(p.fw, x1, a, b)), (Some((l0, _)), Some((l1, _))) if l0.abs_diff(y0) <= tol && l1.abs_diff(y1) <= tol)
                        });
                        if let Some(dbg) = std::env::var_os("RUSTLE_COPYSET_DEBUG_TID") {
                            if p.uniq_tid == dbg.to_string_lossy() || q.uniq_tid == dbg.to_string_lossy() {
                                let detail: Vec<String> = p.t.introns.iter().zip(q.t.introns.iter()).map(|(&(x0, x1), &(y0, y1))| format!("{}-{}=>{:?}/{:?} vs {}-{}", x0, x1, lift_pos(p.fw, x0, a, b), lift_pos(p.fw, x1, a, b), y0, y1)).collect();
                                eprintln!("[copyset] pair {}@{} vs {}@{}: ok={ok} {}", p.uniq_tid, p.cidx, q.uniq_tid, q.cidx, detail.join(" | "));
                            }
                        }
                        if ok { let (ri, rj) = (find(&mut parent, i), find(&mut parent, j)); parent[ri] = rj; }
                    }
                }
                let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
                for i in 0..n { let r = find(&mut parent, i); groups.entry(r).or_default().push(i); }
                let (mut n_kept, mut n_drop, mut n_lift, mut n_lift_fail, mut n_und) = (0usize, 0usize, 0usize, 0usize, 0usize);
                // ⭐ B1 (`docs/OPEN_ITEMS_2026-09-09.md`): a lift failure used to be an anonymous count — a
                // certificate-assigned read's evidence at a copy silently had no transcript anywhere, which
                // contradicts "the GTF O2 believes". Named here (source transcript, source copy, target
                // copy) and printed explicitly below, so it is auditable instead of a bare "N lifts failed".
                let mut lift_fail_detail: Vec<(String, String, String)> = Vec::new();
                let fmt_map = |m: &std::collections::BTreeMap<String, usize>| -> String { let mut v: Vec<(&String, &usize)> = m.iter().collect(); v.sort_by_key(|(k, _)| k.parse::<i64>().unwrap_or(i64::MAX)); v.iter().map(|(k, c)| format!("{k}:{c}")).collect::<Vec<_>>().join(",") };
                // `starts_with("outside")` (not `==`) so this sorts correctly whether A6's
                // `--name-outside-tie` is on (`outside:chrom:start-end`) or off (bare `outside`).
                let fmt_set = |s: &std::collections::BTreeSet<String>| -> String { let mut v: Vec<&String> = s.iter().collect(); v.sort_by_key(|k| (k.starts_with("outside"), k.parse::<i64>().unwrap_or(i64::MAX))); v.iter().map(|k| k.as_str()).collect::<Vec<_>>().join(",") };
                let copyset_debug = std::env::var_os("RUSTLE_COPYSET_DEBUG").is_some();
                if copyset_debug {
                    for &fw in &fam_ids {
                        let l = &lifts_by_fam[&fw];
                        eprintln!("[copyset] family {fw}: {} copy spans, {} lift pairs ({} fragments); spans: {}", fams[fw].copy_spans.len(), l.len(), l.values().map(|v| v.len()).sum::<usize>(),
                            fams[fw].copy_spans.iter().enumerate().map(|(ci, (c, s0, e0))| format!("{}={}:{}-{}", sweep_to_catalog(fw, ci), c, s0, e0)).collect::<Vec<_>>().join(" "));
                    }
                }
                for (_, members) in groups {
                    let fw = pending[members[0]].fw;
                    if copyset_debug {
                        eprintln!("[copyset] group: {}", members.iter().map(|&i| format!("{}@{}", pending[i].uniq_tid, pending[i].cidx)).collect::<Vec<_>>().join(" "));
                    }
                    let (mut uniq, mut asg) = (std::collections::BTreeMap::<String, usize>::new(), std::collections::BTreeMap::<String, usize>::new());
                    let mut und: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
                    let mut n_abst = 0usize;
                    let mut outside_bare = false; // ran out of locus info for >=1 outside tie
                    let mut outside_raw: Vec<(String, u64, u64)> = Vec::new(); // pre-merge, A6
                    for &i in &members {
                        for (k, c) in &pending[i].uniq { *uniq.entry(k.clone()).or_insert(0) += c; }
                        for (k, c) in &pending[i].asg { *asg.entry(k.clone()).or_insert(0) += c; }
                        for name in &pending[i].abst {
                            n_abst += 1;
                            if let Some((set, outside)) = tie_set_of(name) {
                                und.extend(set);
                                if outside {
                                    // A6 (`docs/OPEN_ITEMS_2026-09-09.md`, register row 790): name the outside
                                    // placement's own locus instead of a bare "outside" — default off so the
                                    // existing `copies_undecided` schema stays byte-identical.
                                    let loci = if args.name_outside_tie {
                                        rustle::vg_family::copy_assign_pipeline::tie_outside_loci(name)
                                    } else {
                                        std::collections::BTreeSet::new()
                                    };
                                    if loci.is_empty() {
                                        outside_bare = true;
                                    } else {
                                        for l in loci {
                                            if let Some((c, rest)) = l.split_once(':') {
                                                if let Some((s, e)) = rest.split_once('-') {
                                                    if let (Ok(s), Ok(e)) = (s.parse::<u64>(), e.parse::<u64>()) {
                                                        outside_raw.push((c.to_string(), s, e));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if outside_bare {
                        und.insert("outside".to_string());
                    }
                    if !outside_raw.is_empty() {
                        // Merge overlapping/near-adjacent (<=1kb gap) outside loci per chromosome — many
                        // near-identical secondary-alignment coordinates for the SAME outside locus otherwise
                        // produce hundreds of one-bp-apart entries, which defeats the point of naming it.
                        outside_raw.sort_unstable();
                        let mut merged: Vec<(String, u64, u64)> = Vec::new();
                        for (c, s, e) in outside_raw {
                            if let Some(last) = merged.last_mut() {
                                if last.0 == c && s <= last.2 + 1000 {
                                    last.2 = last.2.max(e);
                                    continue;
                                }
                            }
                            merged.push((c, s, e));
                        }
                        for (c, s, e) in merged {
                            und.insert(format!("outside:{c}:{s}-{e}"));
                        }
                    }
                    let evidence: std::collections::BTreeSet<String> = uniq.keys().chain(asg.keys()).cloned().collect();
                    let base_attr = format!(" evidence_unique \"{}\"; evidence_assigned \"{}\"; reads_undecided \"{}\";", fmt_map(&uniq), fmt_map(&asg), n_abst);
                    let rep = *members.iter().max_by_key(|&&i| pending[i].uniq.values().sum::<usize>() + pending[i].asg.values().sum::<usize>() + pending[i].abst.len()).unwrap();
                    if evidence.is_empty() {
                        let p = &pending[rep];
                        gtf_lines.push(format!("{}{} copies \"{}\"; placed_by \"aligner_primary\"; copy_status_final \"undecidable\";", p.tline, base_attr, fmt_set(&und)));
                        gtf_lines.extend(p.elines.iter().cloned());
                        n_und += 1;
                        n_drop += members.len() - 1;
                        continue;
                    }
                    let und_only: std::collections::BTreeSet<String> = und.difference(&evidence).cloned().collect();
                    let mut have: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
                    for &i in &members {
                        let p = &pending[i];
                        if evidence.contains(&p.cidx) {
                            have.insert(p.cidx.clone());
                            let by = if uniq.contains_key(&p.cidx) { "unique_mapper" } else { "assigned_read" };
                            gtf_lines.push(format!("{}{} copies \"{}\"; copies_undecided \"{}\"; placed_by \"{by}\";", p.tline, base_attr, fmt_set(&evidence), fmt_set(&und_only)));
                            gtf_lines.extend(p.elines.iter().cloned());
                            n_kept += 1;
                        } else {
                            n_drop += 1; // a phantom: nothing at this copy tells it from the evidence copies
                        }
                    }
                    let p = &pending[rep];
                    for c in evidence.difference(&have) {
                        let (Some(a), Some(b)) = (sweep_ci(fw, &p.cidx), sweep_ci(fw, c)) else {
                            n_lift_fail += 1;
                            lift_fail_detail.push((p.uniq_tid.clone(), p.cidx.clone(), c.clone()));
                            continue;
                        };
                        let mut prev = p.t.start;
                        let mut exons: Vec<(u64, u64)> = Vec::new();
                        for &(d, aa) in &p.t.introns { exons.push((prev, d)); prev = aa; }
                        exons.push((prev, p.t.end));
                        let mut lifted: Vec<(u64, u64)> = Vec::new();
                        let mut ok = true;
                        for &(es, ee) in &exons {
                            match (lift_pos(fw, es, a, b), lift_pos(fw, ee, a, b)) {
                                (Some((l0, d0)), Some((l1, d1))) if d0 <= tol && d1 <= tol => lifted.push((l0.min(l1), l0.max(l1))),
                                _ => { ok = false; break; }
                            }
                        }
                        lifted.sort_unstable();
                        if !ok || lifted.windows(2).any(|w| w[0].1 > w[1].0) {
                            n_lift_fail += 1;
                            lift_fail_detail.push((p.uniq_tid.clone(), p.cidx.clone(), c.clone()));
                            continue;
                        }
                        let strand = strand_of.get(c).copied().unwrap_or(p.t.strand);
                        let (chrom, _, _) = &fams[fw].copy_spans[b];
                        let tid = format!("{}_lift{}", p.uniq_tid, c);
                        let gene = format!("{}_copy{}", p.t.gene_tid, c);
                        let attrs = p.tline.splitn(9, '\t').nth(8).unwrap_or("").to_string();
                        let attrs = attrs.replace(&format!("transcript_id \"{}\"", p.uniq_tid), &format!("transcript_id \"{tid}\""))
                            .replace(&format!("gene_id \"{}\"", p.t.gene_tid), &format!("gene_id \"{gene}\""));
                        let attrs = match attrs.find("copy_index \"") {
                            Some(i) => match attrs[i + 12..].find('"') { Some(j) => format!("{}copy_index \"{c}{}", &attrs[..i], &attrs[i + 12 + j..]), None => attrs },
                            None => attrs,
                        };
                        gtf_lines.push(format!("{}\trustle\ttranscript\t{}\t{}\t.\t{}\t.\t{}{} copies \"{}\"; copies_undecided \"{}\"; placed_by \"assigned_read\"; lifted_from \"{}\";",
                            chrom, lifted[0].0 + 1, lifted[lifted.len() - 1].1, strand, attrs, base_attr, fmt_set(&evidence), fmt_set(&und_only), p.uniq_tid));
                        let order: Vec<(u64, u64)> = if strand == '+' { lifted.clone() } else { lifted.iter().rev().cloned().collect() };
                        for (k, (es, ee)) in order.iter().enumerate() {
                            gtf_lines.push(format!("{}\trustle\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; exon_number \"{}\";", chrom, es + 1, ee, strand, k + 1));
                        }
                        n_lift += 1;
                    }
                }
                eprintln!("[copy_assign]   ⭐ --gtf-copy-set {contig}:{lo}-{hi}: {} family isoforms placed by evidence, {} phantoms/duplicates dropped, {} lifted placements ({} lifts failed), {} undecided isoforms emitted once with a copy set",
                    n_kept, n_drop, n_lift, n_lift_fail, n_und);
                // B1: name every lift failure — which transcript, which source copy, which target copy has
                // evidence but no transcript anywhere. Never quote "N lifts failed" without this list beside it.
                for (tid, from, to) in &lift_fail_detail {
                    eprintln!("[copy_assign]     lift_failed: transcript {tid} (source copy {from}) has evidence at copy {to} but could not be placed there");
                }
            }
            // ⭐ B2 (`docs/OPEN_ITEMS_2026-09-09.md`, 09-10) + read-provenance (09-10, same day, user request):
            // EVIDENCE-BACKED SINGLETON rescue, now paired with a COMPLETE per-record accounting.
            // `assemble_gate` drops any exon chain with < GATE_MIN_READS (3) reads, so a certificate-
            // assigned read whose own chain never reaches that floor has NO transcript in the emitted GTF
            // at all — measured: 220 of 290 contested molecules O2 assigned but the GTF does not carry are
            // singletons (row 798, §6hh), and the flagship case is copy 22's SV-carrying isoform. The
            // certificate is a STRONGER claim than "3 reads agree" (a per-read significance test against
            // every other candidate copy), so a read that clears it deserves a transcript even at support
            // 1 — tagged, not silently equal to a normal gate-passed model.
            //
            // ⚠ CORRECTNESS FIX over the first form of this rescue: that version grouped by a record's own
            // chain and labelled it with the CERTIFICATE'S assigned copy without checking the record's own
            // genomic position actually overlaps that copy's span. In a secondary-alignment-heavy region a
            // molecule's records sit at DIFFERENT loci (primary at A, secondary at B); the fix requires the
            // rescued record to be the one AT the assigned copy's own span (`fams[..].copy_spans`), else it
            // is excluded with an explicit `chain_at_unassigned_locus` reason rather than silently mislabelled.
            //
            // `--read-provenance`: emit ONE row per AS-tied alignment record covering EVERY exclusion, not
            // only the rescued ones — "complete, not necessarily good": every record this region's O2 saw
            // gets exactly one row saying which transcript it became, or exactly why it did not.
            if args.gtf && (args.rescue_singletons || args.read_provenance) {
                // Gate-passed lookup, position-aware for the empty (unspliced) chain — see the comment at
                // `chain_uniq_tid`'s declaration above. Never match an unspliced record by chain alone.
                let gate_passed_tid_for = |chain: &[(u64, u64)], s0: u64, e0: u64| -> Option<&str> {
                    if chain.is_empty() {
                        unspliced_gate_passed.iter().find(|(ts, te, _)| s0 < *te && *ts < e0).map(|(_, _, tid)| tid.as_str())
                    } else {
                        chain_uniq_tid.get(chain).map(|s| s.as_str())
                    }
                };
                let copy_span_of = |fid: &str, ci: usize| -> Option<(String, u64, u64)> {
                    fams.iter().find(|f| f.family_id == fid).and_then(|f| f.copy_spans.get(ci)).cloned()
                };
                enum Rec<'a> {
                    ExcludedNoSpan,
                    ExcludedSupplementary,
                    ContributesGatePassed(&'a str),
                    ExcludedNoCertificate,
                    ExcludedNotAssigned(&'a str),
                    ExcludedOffLocus(String),
                    Eligible { chain: Vec<(u64, u64)>, s0: u64, e0: u64, strand: char, catalog_idx: &'a str },
                }
                struct Resc { starts: Vec<u64>, ends: Vec<u64>, introns: Vec<(u64, u64)>, catalog_idx: String, n: usize, fwd: u32, rev: u32 }
                let mut groups: std::collections::BTreeMap<Vec<(u64, u64)>, Resc> = std::collections::BTreeMap::new();
                let mut classified: Vec<(usize, Rec)> = Vec::with_capacity(bam_reads.len());
                for (ri, name) in bam_reads.iter().enumerate() {
                    let is_supplementary = read_spans.get(ri).map_or(false, |&(_, _, f)| f & 2 != 0);
                    let chain = read_chain.get(ri).cloned().unwrap_or_default();
                    let rec = if read_blocks.get(ri).map_or(true, |b| b.is_empty()) {
                        Rec::ExcludedNoSpan
                    } else if is_supplementary {
                        Rec::ExcludedSupplementary
                    } else {
                        let blocks = &read_blocks[ri];
                        let (s0, e0) = (blocks.first().unwrap().0, blocks.last().unwrap().1);
                        if let Some(tid) = gate_passed_tid_for(&chain, s0, e0) {
                            Rec::ContributesGatePassed(tid)
                        } else if let Some(row) = verdict.get(name.as_str()) {
                            if row.status != "assigned" {
                                Rec::ExcludedNotAssigned(row.status)
                            } else {
                                let span = copy_span_of(&row.family_id, row.assigned_copy);
                                let overlaps = span.as_ref().is_some_and(|(c, s, e)| c == contig && s0 < *e && e0 > *s);
                                if !overlaps {
                                    Rec::ExcludedOffLocus(row.catalog_copy_idx.clone())
                                } else {
                                    let strand = read_strand.get(ri).copied().unwrap_or('+');
                                    Rec::Eligible { chain, s0, e0, strand, catalog_idx: row.catalog_copy_idx.as_str() }
                                }
                            }
                        } else {
                            Rec::ExcludedNoCertificate
                        }
                    };
                    classified.push((ri, rec));
                }
                if args.rescue_singletons {
                    for (_, rec) in &classified {
                        if let Rec::Eligible { chain, s0, e0, strand, catalog_idx } = rec {
                            // `chain` is guaranteed to have found no gate-passed match here (that is exactly
                            // the condition classification checked before ever producing `Rec::Eligible`),
                            // via `gate_passed_tid_for`'s position-aware lookup — not raw chain membership.
                            let e = groups.entry(chain.clone()).or_insert_with(|| Resc {
                                starts: Vec::new(), ends: Vec::new(), introns: chain.clone(), catalog_idx: catalog_idx.to_string(), n: 0, fwd: 0, rev: 0,
                            });
                            e.starts.push(*s0);
                            e.ends.push(*e0);
                            e.n += 1;
                            if *strand == '-' { e.rev += 1 } else { e.fwd += 1 }
                        }
                    }
                }
                let mut n_rescued = 0usize;
                let mut rescue_tid_of: std::collections::HashMap<Vec<(u64, u64)>, String> = std::collections::HashMap::new();
                if args.rescue_singletons {
                    for (i, (chain, g)) in groups.iter().enumerate() {
                        let s0 = *g.starts.iter().min().unwrap();
                        let e0 = *g.ends.iter().max().unwrap();
                        let tid = format!("RESCUE_{contig}_{s0}_{i}");
                        let gene = format!("RESCUE_{contig}_{s0}");
                        // Strand: majority vote of the group's own reads (`read_strand`, `ts` flipped by
                        // alignment orientation, falling back to the read's own FLAG 0x10 when minimap2
                        // emitted no `ts`) — a read-specific call, not a flat placeholder. A tie resolves to
                        // `'+'`, matching `majority_read_strand`'s own convention elsewhere in this codebase.
                        let strand = if g.rev > g.fwd { '-' } else { '+' };
                        gtf_lines.push(format!(
                            "{contig}\trustle\ttranscript\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; copies \"{}\"; placed_by \"assigned_read_singleton\"; support \"{}\"; low_confidence \"true\"; low_confidence_reason \"singleton_rescued\";",
                            s0 + 1, e0, g.catalog_idx, g.n
                        ));
                        let mut exons = Vec::new();
                        let mut prev = s0;
                        for &(d, a) in &g.introns {
                            exons.push((prev, d));
                            prev = a;
                        }
                        exons.push((prev, e0));
                        let order: Vec<(u64, u64)> = if strand == '+' { exons.clone() } else { exons.iter().rev().cloned().collect() };
                        for (k, &(es, ee)) in order.iter().enumerate() {
                            gtf_lines.push(format!("{contig}\trustle\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{tid}\"; exon_number \"{}\";", es + 1, ee, k + 1));
                        }
                        rescue_tid_of.insert(chain.clone(), tid);
                        n_rescued += 1;
                    }
                    if n_rescued > 0 {
                        eprintln!("[copy_assign]   ⭐ --rescue-singletons {contig}:{lo}-{hi}: {n_rescued} certificate-assigned read(s) emitted as low-support transcripts the min-reads gate dropped");
                    }
                }
                if args.read_provenance {
                    for (ri, rec) in &classified {
                        let name = &bam_reads[*ri];
                        let (chain_str, tid, reason) = match rec {
                            Rec::ExcludedNoSpan => ("NA".to_string(), "NA".to_string(), "excluded_no_aligned_span"),
                            Rec::ExcludedSupplementary => ("NA".to_string(), "NA".to_string(), "excluded_supplementary"),
                            Rec::ContributesGatePassed(tid) => (fmt_chain(&read_chain[*ri]), tid.to_string(), "contributed_gate_passed"),
                            Rec::ExcludedNoCertificate => ("NA".to_string(), "NA".to_string(), "excluded_no_certificate_row"),
                            Rec::ExcludedNotAssigned(st) => (fmt_chain(&read_chain[*ri]), "NA".to_string(),
                                match *st { "ambiguous" => "excluded_ambiguous", "tied" => "excluded_tied", _ => "excluded_unassigned" }),
                            Rec::ExcludedOffLocus(cidx) => (fmt_chain(&read_chain[*ri]), format!("NA(assigned_copy={cidx})"), "excluded_chain_at_unassigned_locus"),
                            Rec::Eligible { chain, .. } => {
                                if !args.rescue_singletons {
                                    (fmt_chain(chain), "NA".to_string(), "eligible_for_rescue_flag_off")
                                } else if let Some(t) = rescue_tid_of.get(chain) {
                                    (fmt_chain(chain), t.clone(), "contributed_rescued_singleton")
                                } else {
                                    // defensive only: every `Eligible` chain is grouped and emitted above,
                                    // so `rescue_tid_of` always has an entry here in practice.
                                    (fmt_chain(chain), "NA".to_string(), "excluded_unexpected_no_rescue_tid")
                                }
                            }
                        };
                        prov_rows.push(format!("{name}\t{contig}\t{chain_str}\t{tid}\t{reason}"));
                    }
                }
            }
        } // for work in works (serial drain, region order)
    } // serial-drain block

    if args.gtf {
        // §6gp: the `productive` call is RELATIVE to the family's best ORF, which is only known once every
        // region has been drained — so it is stamped here, in a second pass over the finished GTF lines.
        // An absolute amino-acid cut was tried on the core rule and discarded 26 % of protein-coding units
        // (§6gb, register 746), which is why the bar is the family's own best rather than a constant.
        if args.productivity {
            // ⚠ The bar is half the family's MEDIAN ORF, not half its best. "Half the best" was tried first
            // and failed the same way it failed for the core rule (§6gb, register 746): one 2,578-aa outlier
            // put the bar at 1,289 aa and called 0 of copy 27's 39 isoforms productive, when every one of
            // them carries ~700 aa. A single long transcript must not define the family's standard.
            let mut all: std::collections::HashMap<&str, Vec<usize>> = std::collections::HashMap::new();
            for (fid, _, _, aa) in &prod_rows {
                all.entry(fid.as_str()).or_default().push(*aa);
            }
            let best: std::collections::HashMap<&str, usize> = all
                .into_iter()
                .map(|(f, mut v)| {
                    v.sort_unstable();
                    (f, v[v.len() / 2])
                })
                .collect();
            let by_tid: std::collections::HashMap<&str, (&str, &str, usize)> = prod_rows
                .iter()
                .map(|(f, c, t, aa)| (t.as_str(), (f.as_str(), c.as_str(), *aa)))
                .collect();
            for line in gtf_lines.iter_mut() {
                if !line.contains("\ttranscript\t") {
                    continue;
                }
                let Some(tid) = re_attr(line, "transcript_id") else { continue };
                let Some(&(fid, _, aa)) = by_tid.get(tid.as_str()) else { continue };
                let bar = best.get(fid).copied().unwrap_or(0);
                let prod = bar > 0 && aa * 2 >= bar;
                line.push_str(&format!(" productive \"{prod}\"; family_median_orf_aa \"{bar}\";"));
            }
            let mut ph = std::fs::File::create(format!("{}.productivity.tsv", args.out))?;
            writeln!(ph, "family_id\tcopy\tisoforms\tproductive\tmedian_orf_aa\tmax_orf_aa\tfamily_median_orf_aa")?;
            let mut per: std::collections::BTreeMap<(&str, &str), Vec<usize>> = std::collections::BTreeMap::new();
            for (f, c, _, aa) in &prod_rows {
                per.entry((f.as_str(), c.as_str())).or_default().push(*aa);
            }
            for ((f, c), mut v) in per {
                v.sort_unstable();
                let bar = best.get(f).copied().unwrap_or(0);
                let n_prod = v.iter().filter(|&&aa| bar > 0 && aa * 2 >= bar).count();
                let f = if f.is_empty() { "NA" } else { f };
                let c = if c.is_empty() { "NA" } else { c };
                writeln!(ph, "{f}\t{c}\t{}\t{n_prod}\t{}\t{}\t{bar}", v.len(), v[v.len() / 2], v[v.len() - 1])?;
            }
            eprintln!("[copy_assign] wrote {}.productivity.tsv ({} isoform(s) with an ORF, bar = half the family median ORF)",
                args.out, prod_rows.len());
        }
        if args.assembly_polish != "none" {
            // §6zb: the polish runs PER CONTIG. Its mono-exonic floor is a quantile of the run's own
            // multi-exon support, validated per chromosome (§6p8-§6q4) and reproduced by the per-contig
            // sweep; a `--genome-wide` run in one process must not pool that quantile across contigs. For a
            // single-region run this is exactly the former single call (one contig), byte-identical.
            let before = gtf_lines.iter().filter(|l| l.contains("\ttranscript\t")).count();
            let mut contigs: Vec<String> = Vec::new();
            for l in gtf_lines.iter() {
                if let Some(c) = l.split('\t').next() {
                    if !l.starts_with('#') && contigs.last().map_or(true, |p| p != c) && !contigs.iter().any(|p| p == c) {
                        contigs.push(c.to_string());
                    }
                }
            }
            let (mut n_ism, mut n_mono, mut n_frac, mut n_ret) = (0usize, 0usize, 0usize, 0usize);
            let mut floors: Vec<u64> = Vec::new();
            let mut polished: Vec<String> = Vec::with_capacity(gtf_lines.len());
            for c in &contigs {
                let mut part: Vec<String> = gtf_lines.iter().filter(|l| !l.starts_with('#') && l.split('\t').next() == Some(c.as_str())).cloned().collect();
                let (a, b, d, floor, e) = polish_gtf_lines(
                    &mut part,
                    &args.assembly_polish,
                    args.polish_mono_quantile,
                    args.polish_isoform_fraction,
                    args.polish_mono_shadow,
                    args.polish_ism_escape,
                    args.polish_ism_3p,
                    args.polish_ism_ratio,
                    args.polish_fraction_exempt,
                    args.polish_fuzzy_junction,
                    args.polish_fuzzy_ism,
                    args.polish_fraction_min_reads,
                    args.polish_retained_ratio,
                );
                n_ism += a; n_mono += b; n_frac += d; n_ret += e; floors.push(floor);
                polished.extend(part);
            }
            let comments: Vec<String> = gtf_lines.iter().filter(|l| l.starts_with('#')).cloned().collect();
            gtf_lines = comments.into_iter().chain(polished).collect();
            let after = gtf_lines.iter().filter(|l| l.contains("\ttranscript\t")).count();
            let floor_s = if floors.len() == 1 { floors[0].to_string() } else { format!("{:?} (per contig)", floors) };
            eprintln!(
                "[copy_assign] ⭐ ASSEMBLY POLISH ({}): {before} transcripts -> ISM dropped {n_ism} -> \
                 mono floor {floor_s} reads dropped {n_mono} -> isoform fraction {} dropped {n_frac} -> \
                 retained-intron ratio {} dropped {n_ret} -> {after} kept",
                args.assembly_polish, args.polish_isoform_fraction, args.polish_retained_ratio
            );
        }
        if args.gtf_tpm {
            let n = annotate_tpm(&mut gtf_lines);
            eprintln!("[copy_assign] ⭐ TPM: annotated {n} transcript lines with count-based `TPM` and `cov` (§6r5)");
        }
        let mut gh = std::fs::File::create(format!("{}.gtf", args.out))?;
        for line in &gtf_lines {
            writeln!(gh, "{line}")?;
        }
        eprintln!("[copy_assign] wrote {}.gtf ({} GTF rows = FLAIR-style isoforms; family copies tagged multicopy)",
            args.out, gtf_lines.len());
    }
    let mut fh = std::fs::File::create(format!("{}.families.tsv", args.out))?;
    writeln!(
        fh,
        "family_id\tchrom\tn_copies\trescued_copies\tcollapsed_copies\tn_reads\tpsv_cols\tresolvable_psv\tresolvable_j\tjunction_only\tassigned_j\tuniq_agree\tuniq"
    )?;
    for r in &family_rows {
        writeln!(
            fh,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family_id, r.chrom, r.n_copies, r.rescued_copies, r.collapsed_copies, r.n_reads, r.psv_cols,
            r.resolvable_psv, r.resolvable_j, r.junction_only, r.assigned_j, r.uniq_agree, r.uniq
        )?;
    }
    if args.read_provenance {
        let mut pvh = std::fs::File::create(format!("{}.read_provenance.tsv", args.out))?;
        // Scope note (read before quoting): covers every alignment RECORD of every AS-TIED molecule this
        // region's O2 saw (`--no-as-tied-only` disables the AS-tied gate entirely, and this file with it —
        // an uncontested unique mapper is not O2's business and gets its transcript the ordinary way). This
        // is deliberately the secondary-alignment-heavy population, not literally every FLNC read in the BAM.
        writeln!(pvh, "read_name\tchrom\tchain\ttranscript_id\treason")?;
        for r in &prov_rows {
            writeln!(pvh, "{r}")?;
        }
        eprintln!("[copy_assign] wrote {}.read_provenance.tsv ({} record(s))", args.out, prov_rows.len());
    }
    let mut ah = std::fs::File::create(format!("{}.assignments.tsv", args.out))?;
    // `tie_outside_catalog` (§6gz) exists only under the gate, so `--no-as-tied-only` stays byte-identical
    // to the pre-2026-09-09 schema.
    let hdr = "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin\tp_value\tmin_p_value\tas_best\tas_second\tas_margin\tas_per_base_best\tas_per_base_2nd\tin_copy\tcatalog_copy_idx\torigin_rejected\tn_candidates\tsole_candidate\tcontested\treadthrough_into\tprimary_local";
    let sibling_hdr = if args.sibling_report { "\tsibling_identity\tn_cols_vs_sibling" } else { "" };
    // §6u6: appended LAST so every existing column keeps its position.
    let eichler_hdr = if args.eichler_margin.is_some() { "\teichler_call\teichler_same_copy" } else { "" };
    if args.no_as_tied_only { writeln!(ah, "{hdr}{sibling_hdr}{eichler_hdr}")?; } else { writeln!(ah, "{hdr}\ttie_outside_catalog\taligner_disagreement{sibling_hdr}{eichler_hdr}")?; }
    for r in &assign_rows {
        // L3: a CONTESTED molecule assigned with exactly one candidate is a sole candidate (§6fi); an uncontested
        // one is assigned to its placement (§6fq) whatever its candidate count
        let sole = (r.status == "assigned" && r.n_candidates == 1 && r.contested) as u8;
        let status = r.status;
        let outside = if args.no_as_tied_only {
            String::new()
        } else {
            format!("\t{}\t{}", rustle::vg_family::copy_assign_pipeline::is_tie_outside(&r.read_name) as u8, is_disagreement(&r.read_name) as u8)
        };
        let sibling = if args.sibling_report {
            format!("\t{:.4}\t{}", r.sibling_identity, r.n_cols_vs_nearest_sibling)
        } else {
            String::new()
        };
        // §6u6: Eichler's rule. A read with NO rival placement (margin None) has nothing within T, so
        // he assigns it; otherwise he needs a margin of at least T. `eichler_same_copy` is decidable
        // only where BOTH rules assign — his pick is the read's best-AS placement, which is exactly
        // what `primary_local` marks for our assigned copy.
        let eichler = match args.eichler_margin {
            None => String::new(),
            Some(t) => {
                let assigns = r.as_ev.margin().is_none_or(|m| m >= t);
                let same = if assigns && r.status == "assigned" {
                    (r.primary_local as u8).to_string()
                } else {
                    "NA".to_string()
                };
                format!("\t{}\t{}", if assigns { "assign" } else { "discard" }, same)
            }
        };
        writeln!(
            ah,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3e}\t{:.3e}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}{}{}",
            r.read_name, r.family_id, r.assigned_copy, status, r.n_decisive, r.margin, r.p_value, r.min_p_value,
            r.as_ev.best, opt_i32(r.as_ev.second), opt_i32(r.as_ev.margin()),
            r.as_ev.best_per_base, opt_f32(r.as_ev.second_per_base), r.in_copy, r.catalog_copy_idx, r.origin_rejected as u8, r.n_candidates, sole, r.contested as u8, r.readthrough_into, r.primary_local as u8, outside, sibling, eichler
        )?;
    }
    let indel_stats = rustle::vg_family::copy_assign_pipeline::take_indel_stats();
    {
        // §6es hygiene: reads with an aligned base inside a copy. ⚠ Kept for continuity only — it counts
        // secondary-only visitors, so it is NOT the denominator to quote (register 734).
        let inc: Vec<&AssignRow> = assign_rows.iter().filter(|r| r.in_copy).collect();
        let cnt = |st: &str| inc.iter().filter(|r| r.status == st).count();
        eprintln!(
            "[copy_assign] reads with an aligned base inside a copy: {} of {} rows — assigned {} / tied {} / ambiguous {}",
            inc.len(),
            assign_rows.len(),
            cnt("assigned"),
            cnt("tied"),
            cnt("ambiguous")
        );
        if !args.no_as_tied_only {
            let (ma, mt, ra, rt) = (
                GATE_MOL_ALL.load(std::sync::atomic::Ordering::Relaxed),
                GATE_MOL_TIED.load(std::sync::atomic::Ordering::Relaxed),
                GATE_REC_ALL.load(std::sync::atomic::Ordering::Relaxed),
                GATE_REC_TIED.load(std::sync::atomic::Ordering::Relaxed),
            );
            eprintln!(
                "[copy_assign] ⭐ AS-TIED GATE (ratio {:.2}): {} of {} molecules in the swept regions are \
                 AS-tied and entered the certificate ({} of {} records); {} unique/clear-best molecules were \
                 skipped before it and are NOT assigned (`--no-as-tied-only` restores them)",
                args.as_tie_ratio, mt, ma, rt, ra, ma.saturating_sub(mt)
            );
            if args.admit_aligner_disagreement {
                eprintln!(
                    "[copy_assign]   ⭐ §6hd ALIGNER-DISAGREEMENT: {} additional molecules admitted whose PRIMARY unit ≠ best-AS unit (no AS tie; `aligner_disagreement` column)",
                    GATE_MOL_DISAGREE.load(std::sync::atomic::Ordering::Relaxed)
                );
            }
            if args.indel_psv {
                eprintln!(
                    "[copy_assign]   ⭐ INDEL-PSV (PREREG 021446fb, floor {} bp): {} molecules gained {} indel columns ({} with an event ≥ 10 bp)",
                    args.indel_psv_min_len, indel_stats.0, indel_stats.1, indel_stats.2
                );
            }
            let mo = GATE_MOL_OUTSIDE.load(std::sync::atomic::Ordering::Relaxed);
            eprintln!(
                "[copy_assign]   ⚠ {} of those tied molecules have a tied placement OUTSIDE every supplied family \
                 target (§6gz) — their competitor was never scored, so they can be `tied` but never `assigned` \
                 (`tie_outside_catalog` column)",
                mo
            );
        }
        // ⭐⭐ O2 SCOPE (user, 2026-09-09): the population copy assignment EXISTS for — AS-tied multimappers,
        // where the aligner's primary/secondary pick was a coin toss. Reported at BOTH tie widths, always,
        // with or without `--as-tied-only`, because the rates above are otherwise read as if every molecule
        // posed a question. One row per MOLECULE (a molecule's records share their AS evidence).
        {
            let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
            let mols: Vec<&AssignRow> =
                assign_rows.iter().filter(|r| seen.insert(r.read_name.as_str())).collect();
            // Under the gate only the gated width is meaningful: molecules outside it never reached the
            // certificate, so a wider decomposition would re-count the same rows. Both widths print only
            // with `--no-as-tied-only`.
            let widths: Vec<f64> = if args.no_as_tied_only { vec![1.0, 0.98] } else { vec![args.as_tie_ratio] };
            for ratio in widths {
                let el: Vec<&&AssignRow> = mols.iter().filter(|r| as_tied(&r.as_ev, ratio) || is_disagreement(&r.read_name)).collect();
                // ⚠⚠ The AS-tied set is NOT yet O2's subject: it still holds molecules the catalog cannot
                // explain (origin-rejected — O3's material) and molecules with a single candidate locus
                // (nothing to choose). The two arms fail this in OPPOSITE ways — gorilla MCL1 is 95.6 %
                // origin-rejected with 0 single-candidate, human MCL0 is 62.0 % single-candidate (§6gv) —
                // so the decomposition is printed, never a single pooled rate.
                let rej = el.iter().filter(|r| r.origin_rejected).count();
                let one = el.iter().filter(|r| !r.origin_rejected && r.n_candidates < 2).count();
                let con: Vec<&&&AssignRow> =
                    el.iter().filter(|r| !r.origin_rejected && r.n_candidates >= 2).collect();
                let c = |st: &str| con.iter().filter(|r| r.status == st).count();
                let pc = |n: usize| if con.is_empty() { 0.0 } else { 100.0 * n as f64 / con.len() as f64 };
                eprintln!(
                    "[copy_assign] AS-TIED @ratio {:.2}: {} of {} molecules — origin-rejected {} (O3's) / \
                     single-candidate {} (nothing to choose) / CONTESTED {}",
                    ratio, el.len(), mols.len(), rej, one, con.len()
                );
                eprintln!(
                    "[copy_assign]   ⭐ over the CONTESTED set (O2's actual subject): \
                     assigned {} ({:.1}%) / tied {} ({:.1}%) / ambiguous {} ({:.1}%)",
                    c("assigned"), pc(c("assigned")), c("tied"), pc(c("tied")),
                    c("ambiguous"), pc(c("ambiguous"))
                );
            }
        }
        // ⭐ register 734: THE denominator. Molecules with a PRIMARY alignment inside a copy of their family.
        let loc: Vec<&AssignRow> = assign_rows.iter().filter(|r| r.primary_local).collect();
        let lcnt = |st: &str| loc.iter().filter(|r| r.status == st).count();
        let pct = |n: usize| if loc.is_empty() { 0.0 } else { 100.0 * n as f64 / loc.len() as f64 };
        eprintln!(
            "[copy_assign] ⭐ RATES ARE OVER THIS SET — molecules with a PRIMARY alignment in a copy: {} of {} rows \
             ({} secondary-only visitors excluded) — assigned {} ({:.1}%) / tied {} ({:.1}%) / ambiguous {} ({:.1}%)",
            loc.len(),
            assign_rows.len(),
            assign_rows.len() - loc.len(),
            lcnt("assigned"), pct(lcnt("assigned")),
            lcnt("tied"), pct(lcnt("tied")),
            lcnt("ambiguous"), pct(lcnt("ambiguous"))
        );
        primary_local_rows = loc.len();
    }

    if args.junction_conflict_abstain {
        let mut cf = std::fs::File::create(format!("{}.conflicts.tsv", args.out))?;
        writeln!(cf, "read_name\tfamily_id\tpsv_best_copy\tstatus\tn_decisive\tmin_p_value")?;
        let mut n = 0usize;
        for r in assign_rows.iter().filter(|r| r.junction_conflict) {
            writeln!(cf, "{}\t{}\t{}\t{}\t{}\t{:.3e}", r.read_name, r.family_id, r.assigned_copy, r.status, r.n_decisive, r.min_p_value)?;
            n += 1;
        }
        eprintln!("[copy_assign] junction-conflict-abstain: {n} read(s) whose splice junctions contradict their PSV-best copy -> ambiguous ({}.conflicts.tsv)", args.out);
    }

    // Reference-free per-family copy number (Task R1, additive; needs no flag): chi_H (PSV
    // conflict-structure lower bound, always computed) + depth_cn (read-depth leg, only when
    // --lambda-global was supplied -- else "NA"). famcn_readonly = the max of the two lower
    // bounds, so it recovers Tier-3 collapsed copies chi_H alone misses.
    let mut cnh = std::fs::File::create(format!("{}.famcn_readonly.tsv", args.out))?;
    writeln!(cnh, "family_id\tchrom\tn_copies\tn_reads\tchi_H\tdepth_cn\tregime\tfamcn_readonly")?;
    for r in &famcn_rows {
        if r.depth_cn.is_finite() {
            let famcn = (r.chi_h as f64).max(r.depth_cn);
            writeln!(
                cnh,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{:.3}",
                r.family_id, r.chrom, r.n_copies, r.n_reads, r.chi_h, r.depth_cn, r.regime, famcn
            )?;
        } else {
            writeln!(
                cnh,
                "{}\t{}\t{}\t{}\t{}\tNA\t{}\t{}",
                r.family_id, r.chrom, r.n_copies, r.n_reads, r.chi_h, r.regime, r.chi_h
            )?;
        }
    }
    eprintln!("[copy_assign] wrote {}.famcn_readonly.tsv ({} families; depth_cn={})",
        args.out, famcn_rows.len(), if lambda.is_some() { "on" } else { "NA (pass --lambda-global or --lambda-file)" });

    // per-read posterior + consistent zone (opt-in via --posterior).
    if args.posterior {
        let mut ph = std::fs::File::create(format!("{}.posterior.tsv", args.out))?;
        writeln!(ph, "read_name\tfamily_id\tstatus\tn_consistent\tzone_chrom\tzone_start\tzone_end\tposterior")?;
        for line in &posterior_lines {
            writeln!(ph, "{line}")?;
        }
        eprintln!("[copy_assign] wrote {}.posterior.tsv ({} reads, prior={})",
            args.out, posterior_lines.len(), if prior_abundance { "abundance" } else { "uniform" });
    }

    // EM soft-relaxation outputs (opt-in via --em): per-read soft posterior + K-frontier label, and the
    // recovered per-copy abundance. Only written under --em; the hard outputs above are unaffected either way.
    if args.em {
        let mut eh = std::fs::File::create(format!("{}.em.tsv", args.out))?;
        writeln!(eh, "read_name\tfamily_id\targmax_copy\tlabel\tposterior\tn_iter")?;
        for l in &em_lines {
            writeln!(eh, "{l}")?;
        }
        let mut eah = std::fs::File::create(format!("{}.em_abundance.tsv", args.out))?;
        writeln!(eah, "family_id\tcopy_id\tpi_hat\tn_reads_soft")?;
        for l in &em_abundance_lines {
            writeln!(eah, "{l}")?;
        }
        eprintln!(
            "[copy_assign] wrote {}.em.tsv ({} reads) + {}.em_abundance.tsv",
            args.out, em_lines.len(), args.out
        );
    }

    // soft per-copy quantification: family/copy, EM abundance ± 95% CI half-width, + the hard read count for
    // comparison. The EM uses partial PSV evidence (the benchmark: beats hard at sparse PSVs; uniform at K=0).
    //
    // SCOPE NOTE (`RUSTLE_XFAM_RECONCILE`): this file is DELIBERATELY unmoved by a cross-family abstention.
    // `n_reads_hard` below counts `fa.assignments` by argmax `best_copy` with NO status filter, and
    // `abundance`/`ci95` come from `soft_quantify_em` inside the per-family pipeline, whose `obs_for_em` is
    // populated regardless of status — so a demoted molecule still counts here, in both arms, and
    // `.quant.tsv` is byte-identical between `report` and `abstain` (`quant_is_unmoved_by_demotion` pins
    // it). The row inflation a double-assigned molecule causes here is therefore NOT fixed by that flag;
    // removing a molecule from a family's EM is a two-pass architecture change and a separate decision.
    let mut qh = std::fs::File::create(format!("{}.quant.tsv", args.out))?;
    writeln!(qh, "family_id\tcopy_index\tcopy_tid\tcopy_chrom\tcopy_start\tcopy_end\tabundance\tci95_halfwidth\tn_reads_hard\tanchored_reads\ttie_invariant\tjunction_invariant\tn_reads_soft")?;
    for r in &quant_rows {
        writeln!(qh, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.2}", r.family_id, r.copy_index, r.copy_tid,
            r.copy_chrom, r.copy_start, r.copy_end, r.abundance, r.ci, r.n_hard, r.anchored, r.tie_invariant,
            r.junction_invariant, r.n_soft)?;
    }
    // A copy is invariant to the arbitrary primary/secondary label if it is pinned by unique mappers OR by a
    // copy-specific junction (splice structure identifies it regardless of the label). Report the OR bottom line.
    let n_inv = quant_rows.iter().filter(|r| r.tie_invariant || r.junction_invariant).count();
    let n_inv_junction_only = quant_rows.iter().filter(|r| !r.tie_invariant && r.junction_invariant).count();
    eprintln!(
        "[copy_assign] tie-break invariance: {}/{} copies invariant (>= {} unique-mapper OR copy-specific-junction reads; FALSE = existence leans on the arbitrary primary label). \
         ⚠ the unique-mapper half is near-vacuous under the default AS-tied gate (A5, register 786) — {n_inv_junction_only} of those {n_inv} are invariant ONLY via copy-specific junctions",
        n_inv, quant_rows.len(), GATE_MIN_READS
    );

    // `--families`: the explicit O1<->O2 JOIN, plus the closing half of the no-silent-drop contract. The
    // input side was validated before any read was touched; this checks the OUTPUT side — that every
    // supplied catalog copy actually came back out as an assigned copy. It is the only place a copy lost
    // inside the assignment stage could be seen, so it is an ERROR, not a log line.
    if let Some(ix) = &catalog_index {
        // O3 Phase 2 (Task 6): the genome-wide Bonferroni flag pass runs exactly ONCE, here, only after
        // every region has drained into `o3_all_raw_pairs` -- `finalize_flags`'s threshold is
        // `alpha / n_pairs_with_a_p_value` over the WHOLE run, so it cannot be computed per-region or
        // per-family. Keyed by `(family_id, copy_idx)`, the same join key `JoinRow` now carries.
        let o3_flags: std::collections::HashMap<(String, String), rustle::vg_family::missing_copy_flag_pass::FlaggedPair> =
            if args.flag_missing_copies {
                rustle::vg_family::missing_copy_flag_pass::finalize_flags(&o3_all_raw_pairs, args.missing_copy_alpha)
                    .into_iter()
                    .map(|fp| ((fp.pair.family_id.clone(), fp.pair.copy_idx.clone()), fp))
                    .collect()
            } else {
                std::collections::HashMap::new()
            };
        let mut jh = std::fs::File::create(format!("{}.family_join.tsv", args.out))?;
        let header = "family_id\tcopy_index\tcopy_tid\tcatalog_family_id\tcatalog_copy_idx\tchrom\tstart\tend\tn_reads_hard";
        if args.flag_missing_copies {
            writeln!(jh, "{header}\to3_flag\to3_class\to3_rate_per_kb\to3_p\to3_n_rejected")?;
        } else {
            writeln!(jh, "{header}")?;
        }
        for r in &join_rows {
            if args.flag_missing_copies {
                match o3_flags.get(&(r.family_id.clone(), r.copy_idx.clone())) {
                    Some(fp) => {
                        let flag_str = match fp.flag {
                            rustle::vg_family::missing_copy_flag_pass::Flag::MissingCopy => "missing_copy",
                            rustle::vg_family::missing_copy_flag_pass::Flag::Untestable => "untestable",
                            rustle::vg_family::missing_copy_flag_pass::Flag::NoFlag => "none",
                        };
                        let class_str = match fp.pair.class {
                            rustle::vg_family::missing_copy_flag_pass::Class::Divergent => "divergent",
                            rustle::vg_family::missing_copy_flag_pass::Class::Structural => "structural",
                        };
                        let rate = if fp.pair.covered_kb > 0.0 { fp.pair.n_sites as f64 / fp.pair.covered_kb } else { 0.0 };
                        let p_str = fp.pair.p_uncorrected.map_or("NA".to_string(), |p| format!("{p:.3e}"));
                        writeln!(jh, "{}\t{flag_str}\t{class_str}\t{rate:.2}\t{p_str}\t{}", r.line, fp.pair.n_rejected)?;
                    }
                    // Fix 2 (final whole-branch review): a genuine `o3_flags` lookup miss (never reached
                    // `detect_missing_copy_pairs` at all -- skipped for <3 rejected reads, or, before Fix 1,
                    // lost to the namespace mismatch) is NOT the same thing as "tested and found clean"
                    // (`none`). Writing the same `none` literal for both conflated them; `not_tested` is a
                    // distinct token so a reader (and `bench/o3_cross_individual_diff.py`) can tell "we have
                    // no information" from "we looked and it was negative".
                    None => writeln!(jh, "{}\tnot_tested\tNA\t0.00\tNA\t0", r.line)?,
                }
            } else {
                writeln!(jh, "{}", r.line)?;
            }
        }
        let emitted: HashSet<&str> = join_rows
            .iter()
            .filter_map(|r| r.line.split('\t').nth(2))
            .collect();
        let missing: Vec<&String> = ix.keys().filter(|t| !emitted.contains(t.as_str())).collect();
        if !missing.is_empty() {
            anyhow::bail!(
                "--families: {} of {} supplied copies did not come back as assigned copies ({}). A supplied \
                 copy must never be dropped; this is a bug in the assignment path, not a filter.",
                missing.len(),
                ix.len(),
                missing.iter().take(10).map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
            );
        }
        eprintln!(
            "[copy_assign] wrote {}.family_join.tsv ({} copies, all {} supplied catalog copies present)",
            args.out,
            join_rows.len(),
            ix.len()
        );
    }
    // O3 Phase 2 (Task 6): candidate unannotated/reference-absent loci found while scanning for orphan
    // read clusters outside every family's own units. Independent of `--families`/`catalog_index` -- an
    // orphan locus is about reads with nowhere to go, not about the catalog join -- so this sits OUTSIDE
    // the `if let Some(ix) = &catalog_index` block above, gated only on the flag itself.
    if args.flag_missing_copies {
        let mut lh = std::fs::File::create(format!("{}.missing_copy_loci.tsv", args.out))?;
        writeln!(lh, "chrom\tstart\tend\tn_reads\tn_orphans\tclass\tn_genes_overlapping\tother_family_units")?;
        for l in &o3_all_orphan_loci {
            let class_str = match l.class {
                rustle::vg_family::missing_copy_flag_pass::LocusClass::OtherFamily => "other_family",
                rustle::vg_family::missing_copy_flag_pass::LocusClass::AnnotatedNoUnit => "annotated_no_unit",
                rustle::vg_family::missing_copy_flag_pass::LocusClass::Unannotated => "unannotated",
            };
            writeln!(
                lh, "{}\t{}\t{}\t{}\t{}\t{class_str}\t{}\t{}",
                l.chrom, l.start, l.end, l.n_reads, l.n_orphans, l.n_genes_overlapping,
                if l.other_family_units.is_empty() { "-".to_string() } else { l.other_family_units.join(";") },
            )?;
        }
        eprintln!("[copy_assign] wrote {}.missing_copy_loci.tsv ({} loci)", args.out, o3_all_orphan_loci.len());
    }

    // `--discover-copies`: read-seeded candidate copies accumulated across every region above. Report only
    // -- never mutates the input catalog or this run's own assignment output. Independent of `--families`/
    // `catalog_index`, same as the O3 orphan-loci block above.
    if args.discover_copies {
        let mut dh = std::fs::File::create(format!("{}.discovered_copies.tsv", args.out))?;
        writeln!(dh, "family_id\tchrom\tstart\tend\tstrand\tn_supporting_reads\tread_names\tnearest_copy_tid\tnearest_copy_distance")?;
        for d in &all_discovered {
            writeln!(
                dh, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                d.family_id, d.chrom, d.start, d.end, d.strand, d.n_supporting_reads,
                d.read_names.join(","), d.nearest_copy_tid,
                // `NA` when the family has no copy on this chromosome at all (`nearest_copy_tid == "NA"`),
                // the same convention `opt_i32`/`opt_f32` use above -- never a `u64::MAX` sentinel printed
                // verbatim as 18446744073709551615.
                opt_u64(d.nearest_copy_distance)
            )?;
        }
        eprintln!(
            "[copy_assign] --discover-copies: {} candidate cop{} -> {}.discovered_copies.tsv",
            all_discovered.len(), if all_discovered.len() == 1 { "y" } else { "ies" }, args.out
        );
    }

    // FACULTATIVE long-read phasing output (dependency-free): phase set (PS) per family, each haplotype's
    // PSV variant string, and read -> haplotype (HP) haplotag. Only written under --phase.
    if args.phase {
        let mut pb = std::fs::File::create(format!("{}.phase_blocks.tsv", args.out))?;
        writeln!(pb, "block_id\tchrom\tn_haplotypes\tn_psv_sites\tn_reads_phased\tn_unphased")?;
        for l in &phase_block_lines {
            writeln!(pb, "{}", l)?;
        }
        let mut ph = std::fs::File::create(format!("{}.phased_haplotypes.tsv", args.out))?;
        writeln!(ph, "block_id\thaplotype\tcopy_tid\tn_support_reads\tvariants")?;
        for l in &phased_hap_lines {
            writeln!(ph, "{}", l)?;
        }
        let mut pr = std::fs::File::create(format!("{}.phased_reads.tsv", args.out))?;
        writeln!(pr, "read_name\tblock_id\thaplotype\tn_psv_spanned\tmargin\tstatus")?;
        for l in &phased_read_lines {
            writeln!(pr, "{}", l)?;
        }
        // self-contained variation graph of the phasing (copies = paths, PSVs = bubbles)
        let mut gf = std::fs::File::create(format!("{}.phase.gfa", args.out))?;
        writeln!(gf, "H\tVN:Z:1.1")?; // GFA 1.1 (W-lines carry the read walks)
        let mut segs: Vec<&String> = gfa_segs.iter().collect();
        segs.sort();
        for s in segs {
            writeln!(gf, "{}", s)?;
        }
        let mut links: Vec<&String> = gfa_links.iter().collect();
        links.sort();
        for l in links {
            writeln!(gf, "{}", l)?;
        }
        for p in &gfa_paths {
            writeln!(gf, "{}", p)?;
        }
        // reads threaded through the family graph (W-lines) — the shared-evidence flip made visible.
        for w in &gfa_walks {
            writeln!(gf, "{}", w)?;
        }
        // Bandage colour CSV: copies distinct, reads coloured by assigned copy (grey = tied/K=0, unresolvable).
        {
            let mut cf = std::fs::File::create(format!("{}.phase.gfa.colours.csv", args.out))?;
            writeln!(cf, "Name,Colour")?;
            for c in &gfa_colors {
                writeln!(cf, "{}", c)?;
            }
        }
        // Legend: status -> colour, de-duplicated across every family's copy-graph (first-seen order, which
        // is deterministic given the serial region/family drain order above).
        {
            let mut lf = std::fs::File::create(format!("{}.phase.gfa.legend.tsv", args.out))?;
            writeln!(lf, "status\tcolour")?;
            let mut seen: HashSet<&String> = HashSet::new();
            for r in &legend_rows {
                if seen.insert(r) {
                    writeln!(lf, "{}", r)?;
                }
            }
        }
        let n_phased = phased_read_lines.iter().filter(|l| !l.contains("\t-1\t")).count();
        eprintln!(
            "[copy_assign] phasing: {} blocks, {} haplotypes, {}/{} reads phased -> {}.phased_*.tsv + {}.phase.gfa ({} bubble-nodes, {} copy-paths, {} read-walks; Bandage colours -> {}.phase.gfa.colours.csv, legend -> {}.phase.gfa.legend.tsv)",
            phase_block_lines.len(), phased_hap_lines.len(), n_phased, phased_read_lines.len(),
            args.out, args.out, gfa_segs.len(), gfa_paths.len(), gfa_walks.len(), args.out, args.out
        );

        // v2: fold every family's exon presence/absence graph (built during the drain above) into
        // <out>.exon.gfa. Each exon's reference bases are fetched HERE, lazily, via genome_for (the
        // builder itself only lays out intervals) — a missing/uncovered stretch falls back to an N-run,
        // counted rather than silently faked (never claim sequence we didn't fetch).
        let n_seq_fallback = std::cell::Cell::new(0usize);
        let exon_seq = |ec: &rustle::vg_family::copy_graph::ExonClass| -> Vec<u8> {
            genome_for(&ec.chrom)
                .ok()
                .and_then(|g| g.fetch_sequence(&ec.chrom, ec.start, ec.end))
                .unwrap_or_else(|| {
                    n_seq_fallback.set(n_seq_fallback.get() + 1);
                    vec![b'N'; (ec.end - ec.start) as usize]
                })
        };
        let mut eg_file = std::fs::File::create(format!("{}.exon.gfa", args.out))?;
        writeln!(eg_file, "H\tVN:Z:1.1")?; // one header; each family's S/L/P folded below (fids are unique, no dedup needed)
        for fam_eg in &exon_graphs {
            for line in fam_eg.to_gfa(&exon_seq).lines().skip(1) {
                // skip the per-family embedded "H\tVN:Z:1.1" line — the file header above already covers it
                writeln!(eg_file, "{}", line)?;
            }
        }
        {
            let mut cf = std::fs::File::create(format!("{}.exon.gfa.colours.csv", args.out))?;
            writeln!(cf, "Name,Colour")?;
            for fam_eg in &exon_graphs {
                for row in fam_eg.colours_csv().lines() {
                    writeln!(cf, "{}", row)?;
                }
            }
        }
        {
            let mut lf = std::fs::File::create(format!("{}.exon.gfa.legend.tsv", args.out))?;
            writeln!(lf, "status\tcolour")?;
            let mut seen: HashSet<String> = HashSet::new();
            for fam_eg in &exon_graphs {
                for row in fam_eg.legend_tsv().lines() {
                    if seen.insert(row.to_string()) {
                        writeln!(lf, "{}", row)?;
                    }
                }
            }
        }
        eprintln!(
            "[copy_assign] exon graph: {} families -> {}.exon.gfa ({} exon-node(s) fell back to an N-run, no fetchable reference sequence); Bandage colours -> {}.exon.gfa.colours.csv, legend -> {}.exon.gfa.legend.tsv",
            exon_graphs.len(), args.out, n_seq_fallback.get(), args.out, args.out
        );
    }

    // gene-conversion events: per-molecule PSV-path switches confirmed by RECURRENCE across reads (vs one-off
    // chimeras). Only written when something was found. The enriched per-molecule multimapper signal.
    if !mosaic_rows.is_empty() {
        let mut mh = std::fs::File::create(format!("{}.mosaic.tsv", args.out))?;
        writeln!(mh, "family_id\tcopy_a\tcopy_b\tbreakpoint_lo\tbreakpoint_hi\tn_reads\tdispersion\tconfirmed")?;
        for r in &mosaic_rows {
            writeln!(
                mh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.family_id, r.copy_a, r.copy_b, r.bp_lo, r.bp_hi, r.n_reads, r.dispersion, r.confirmed as u8
            )?;
        }
        let conf = mosaic_rows.iter().filter(|r| r.confirmed).count();
        eprintln!(
            "[copy_assign] {} read-level gene-conversion event(s) ({} confirmed by recurrence) -> {}.mosaic.tsv",
            mosaic_rows.len(), conf, args.out
        );
    }

    // copy-level historical gene conversions (a de-novo copy whose PSV-allele vector is a mosaic of two
    // others) -- the APOBEC3/RFPL signal, baked into the copy sequence. Written only when found.
    if !copyconv_rows.is_empty() {
        let mut ch = std::fs::File::create(format!("{}.copy_conversions.tsv", args.out))?;
        writeln!(ch, "family_id\tconverted_copy\tdonor_a\tdonor_b\tbreakpoint_lo\tbreakpoint_hi\tn_decisive")?;
        for r in &copyconv_rows {
            writeln!(
                ch, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.family_id, r.copy_c, r.copy_a, r.copy_b, r.bp_lo, r.bp_hi, r.n_decisive
            )?;
        }
        eprintln!(
            "[copy_assign] {} COPY-level historical gene conversion(s) -> {}.copy_conversions.tsv",
            copyconv_rows.len(), args.out
        );
    }

    // ⭐ L6 --dump-star: each molecule's read-star proof (its own columns, its bases, every candidate's bases).
    if args.dump_star {
        let proofs = rustle::vg_family::copy_assign_pipeline::take_star_proofs();
        let mut sh = std::fs::File::create(format!("{}.star_reads.tsv", args.out))?;
        writeln!(sh, "read_name\tfamily_id\tstatus\tassigned_copy\tcatalog_copy_idx\tn_candidates\tcandidates\tn_cols\tcolumns")?;
        let mut n = 0usize;
        for r in &assign_rows {
            let Some(pf) = proofs.get(&r.read_name) else { continue };
            let ch = |o: Option<u8>| o.map(|b| b as char).unwrap_or('.');
            let cols: Vec<String> = pf
                .cols
                .iter()
                .enumerate()
                .map(|(j, &pos)| {
                    let cands: String = pf.alleles.iter().map(|al| ch(al.get(j).copied().flatten())).collect();
                    format!("{pos}:{}:{cands}", ch(pf.obs.get(j).copied().flatten()))
                })
                .collect();
            writeln!(
                sh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.read_name, r.family_id, r.status, r.assigned_copy, r.catalog_copy_idx, r.n_candidates,
                pf.cand.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","), pf.cols.len(), cols.join(",")
            )?;
            n += 1;
        }
        eprintln!("[copy_assign] --dump-star: {n} molecule proofs -> {}.star_reads.tsv (columns = read_pos:read_base:candidate_bases in `candidates` order, `.` = uncovered)", args.out);
    }
    // --dump-psv: the raw per-molecule PSV genotype matrix (the assignment-proof evidence). reads × PSV columns
    // (each read's base + its assignment), the per-copy alleles, and the column→genome map — for the figure.
    if args.dump_psv {
        let mut rh = std::fs::File::create(format!("{}.psv_reads.tsv", args.out))?;
        writeln!(rh, "read_name\tfamily_id\tassigned_copy\tstatus\tmargin\tn_decisive\talleles")?;
        for l in &psv_read_lines {
            writeln!(rh, "{l}")?;
        }
        let mut ch = std::fs::File::create(format!("{}.psv_copies.tsv", args.out))?;
        writeln!(ch, "family_id\tcopy_index\tcopy_tid\talleles\tstrand")?;
        for l in &psv_copy_lines {
            writeln!(ch, "{l}")?;
        }
        let mut lh = std::fs::File::create(format!("{}.psv_cols.tsv", args.out))?;
        writeln!(lh, "family_id\tcol_index\tgenome_pos")?;
        for l in &psv_col_lines {
            writeln!(lh, "{l}")?;
        }
        eprintln!(
            "[copy_assign] dumped PSV genotype matrix: {} read rows, {} copy rows -> {}.psv_reads.tsv/.psv_copies.tsv/.psv_cols.tsv",
            psv_read_lines.len(), psv_copy_lines.len(), args.out
        );
        if args.igv {
            eprintln!(
                "[copy_assign] --igv: now run  python bench/igv_tracks.py --assignments {0}.assignments.tsv \
--bam <bam> --regions <regions> --out {0}  -> {0}.tagged.bam + {0}.copies.bed + {0}.psv.vcf (load in IGV)",
                args.out
            );
        }
    }

    // document every family edge confirmed via the large-sequence LCS fallback (poasta memory threshold
    // exceeded), so edges resting on the approximate metric are auditable. Only written when the fallback ran.
    if !fallback_all.is_empty() {
        let mut sh = std::fs::File::create(format!("{}.fallback.tsv", args.out))?;
        writeln!(sh, "chrom\ttid_a\tstart_a\tend_a\tlen_a\ttid_b\tstart_b\tend_b\tlen_b")?;
        for s in &fallback_all {
            writeln!(
                sh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.chrom, s.tid_a, s.start_a, s.end_a, s.len_a, s.tid_b, s.start_b, s.end_b, s.len_b
            )?;
        }
        eprintln!(
            "[copy_assign] {} family edge(s) confirmed via the large-seq LCS fallback (transcript > --max-poa-len={}); wrote {}.fallback.tsv",
            fallback_all.len(),
            args.max_poa_len,
            args.out
        );
    }

    // --absent-copies: surface candidates that failed the admission gate and need DNA-level validation.
    // Only written when --absent-copies is set so an OFF run produces exactly the same output files.
    if args.absent_copies {
        let mut dh = std::fs::File::create(format!("{}.dna_needs.tsv", args.out))?;
        writeln!(dh, "chrom\tstart\tend\tn_clusters\tread_count\treason")?;
        for r in &dna_needs_rows {
            writeln!(dh, "{}\t{}\t{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.n_clusters, r.read_count, r.reason)?;
        }
        eprintln!(
            "[copy_assign] {} DNA-needs candidate(s) -> {}.dna_needs.tsv",
            dna_needs_rows.len(),
            args.out
        );

        // Augment-and-linearize certificates (Task 4/5): one row per Stage-2-admitted candidate, report-only
        // (`--linearize`) or gating admission too (`--linearize-gate`, which implies the report). ONLY written
        // when the opt-in is set (`do_linearize`) -- plain `--absent-copies` skips the certificate entirely
        // (no minimap2 per candidate) and emits no `.linearize.tsv`, byte-identical to the pre-feature path.
        if do_linearize {
            let mut lh = std::fs::File::create(format!("{}.linearize.tsv", args.out))?;
            writeln!(
                lh,
                "family_id\tchrom\tstart\tend\tn_pool\tlinearized_frac_real\tmean_frac_decoy\tdelta\tperm_p\tverdict"
            )?;
            for (fam, cert, (chrom, start, end)) in &linearize_certs_all {
                writeln!(lh, "{}", linearize_tsv_row(fam, (chrom, *start, *end), cert))?;
            }
            eprintln!(
                "[copy_assign] {} linearize certificate(s) -> {}.linearize.tsv{}",
                linearize_certs_all.len(),
                args.out,
                if args.linearize_gate { " (--linearize-gate: non-LINEARIZES candidates demoted to .dna_needs.tsv)" } else { "" }
            );
        }
    }

    // --vg-realign: the re-align supplement's per-family/per-read decisions (report-only — not fed back
    // into the assignment). Only written when --vg-realign is set so an OFF run produces exactly the same
    // output files (cfg.vg_realign OFF also means fa.realign_records is always empty, so this is belt-
    // and-suspenders with the flag check).
    if args.vg_realign || args.vg_realign_correct {
        let mut vh = std::fs::File::create(format!("{}.vg_realign.tsv", args.out))?;
        writeln!(vh, "read_name\tfamily_id\taction\ttarget_copy\tid_best\tlinear_copy")?;
        for l in &vg_realign_lines {
            writeln!(vh, "{l}")?;
        }
        eprintln!(
            "[copy_assign] {} vg-realign decision(s) -> {}.vg_realign.tsv",
            vg_realign_lines.len(),
            args.out
        );
    }

    let (uniq, agree): (usize, usize) = family_rows.iter().fold((0, 0), |(u, g), f| (u + f.uniq, g + f.uniq_agree));
    eprintln!(
        "[copy_assign] {} families, {} read assignments",
        family_rows.len(),
        assign_rows.len()
    );

    // ---- <out>.union_certificate.tsv (--union-certificate only) ------------------------------------
    // The union verdict per touched molecule. A NEW file, like `xfam_conflicts.tsv`: the OFF arm's outputs
    // stay byte-identical, and the ON arm announces itself by this file plus its `params.tsv` row.
    if args.union_certificate {
        let mut uh = std::fs::File::create(format!("{}.union_certificate.tsv", args.out))?;
        writeln!(uh, "read_name\tn_candidates\tcandidates\twinner\tn_decisive\tmargin\tp_value\tverdict")?;
        for r in &union_all.rows {
            writeln!(
                uh,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3e}\t{}",
                r.read_name, r.n_candidates, r.candidates, r.winner, r.n_decisive, r.margin, r.p_value, r.verdict
            )?;
        }
        eprintln!(
            "[union] --union-certificate: {} molecule(s) in scope over every region, {} scored in {} group(s) \
             (assigned to a family {}, to an outside locus {}, tied {}, ambiguous {}, no result {}), {} left as \
             today (tie partner in another region), {} pseudo-copies unbuildable, {} row(s) added -> \
             {}.union_certificate.tsv",
            union_all.n_in_scope, union_all.n_scored(), union_all.n_groups, union_all.n_assigned_family,
            union_all.n_assigned_outside, union_all.n_tied, union_all.n_ambiguous, union_all.n_no_result,
            union_all.n_other_region, union_all.n_pseudo_unbuildable, union_all.n_rows_added, args.out
        );
    }

    // ---- <out>.xfam_conflicts.tsv (report/abstain only) -------------------------------------------
    // The NEW information goes to a NEW file: adding a column to any existing output would break the OFF
    // arm's byte-identity, which is the gate this change is judged on.
    if xfam_mode != XfamMode::Off {
        let mut xh = std::fs::File::create(format!("{}.xfam_conflicts.tsv", args.out))?;
        writeln!(
            xh,
            "read_name\tstratum\tfamily_a\tcopy_a\tchrom_a\tstart_a\tend_a\tfamily_b\tcopy_b\tchrom_b\tstart_b\tend_b\tsame_record\tsep_bp\tdemoted"
        )?;
        for c in &xfam_conflicts {
            writeln!(
                xh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.read_name, c.stratum, c.fid_a, c.copy_a, c.span_a.0, c.span_a.1, c.span_a.2,
                c.fid_b, c.copy_b, c.span_b.0, c.span_b.1, c.span_b.2, c.same_record, c.sep_bp, c.demoted
            )?;
        }
        // Every rate with its denominator (a rate without one is how four headlines died here).
        let n_pairs = xfam_conflicts.len();
        let n_shared = xfam_conflicts.iter().filter(|c| c.stratum == "shared_locus").count();
        let n_rt = xfam_conflicts.iter().filter(|c| c.stratum == "readthrough_span").count();
        let n_contra = xfam_conflicts.iter().filter(|c| c.stratum == "cross_family_contradiction").count();
        let contested_reads: std::collections::BTreeSet<&str> =
            xfam_conflicts.iter().map(|c| c.read_name.as_str()).collect();
        let contested_rows: std::collections::BTreeSet<(&str, &str)> = xfam_conflicts
            .iter()
            .flat_map(|c| [(c.read_name.as_str(), c.fid_a.as_str()), (c.read_name.as_str(), c.fid_b.as_str())])
            .collect();
        // Denominator = the assigned rows the OFF arm would have emitted. A demotion only ever turns
        // "assigned" into "ambiguous", so that is the post-demotion count plus the demotions — quoting the
        // POST-demotion count would shrink the denominator by exactly the numerator.
        let n_assigned = assign_rows.iter().filter(|r| r.status == "assigned").count()
            + if xfam_mode == XfamMode::Abstain { xfam_demote.len() } else { 0 };
        let distinct_reads: std::collections::BTreeSet<&str> =
            assign_rows.iter().map(|r| r.read_name.as_str()).collect();
        let demoted_mols: std::collections::BTreeSet<&str> =
            xfam_demote.iter().map(|(n, _, _)| n.as_str()).collect();
        // Distinct molecules with at least one assigned row IN THE OFF ARM — the honest denominator for
        // "how many molecules did this demote": the post-demotion set is smaller by the numerator.
        let distinct_assigned: std::collections::BTreeSet<&str> = assign_rows
            .iter()
            .filter(|r| r.status == "assigned")
            .map(|r| r.read_name.as_str())
            .chain(demoted_mols.iter().copied())
            .collect();
        eprintln!(
            "[xfam] RUSTLE_XFAM_RECONCILE={}: {n_pairs} contested assigned-copy pair(s) over \
             {}/{} distinct molecules — {n_shared} shared_locus, {n_rt} readthrough_span, \
             {n_contra} cross_family_contradiction (the only demoting stratum). \
             Contested (read,family) rows: {}. -> {}.xfam_conflicts.tsv",
            xfam_mode.as_str(), contested_reads.len(), distinct_reads.len(), contested_rows.len(), args.out
        );
        if xfam_mode == XfamMode::Abstain {
            let n_dem = xfam_demote.len();
            eprintln!(
                "[xfam] abstain: demoted {n_dem}/{n_assigned} assigned rows ({:.4}), \
                 {}/{} distinct assigned molecules ({:.4}), {n_dem}/{} CONTESTED rows ({:.4}). \
                 The third rate is the one that matters: a rule that only ever fires on \
                 multi-placement molecules is selection-biased by construction, so the run-wide rate \
                 alone hides where all of it lands. Demotion is a STATUS change only: no row is added \
                 or deleted, and .quant.tsv (n_reads_hard, abundance, ci95) is unmoved by design.",
                if n_assigned > 0 { n_dem as f64 / n_assigned as f64 } else { 0.0 },
                demoted_mols.len(), distinct_assigned.len(),
                if !distinct_assigned.is_empty() { demoted_mols.len() as f64 / distinct_assigned.len() as f64 } else { 0.0 },
                contested_rows.len(),
                if !contested_rows.is_empty() { n_dem as f64 / contested_rows.len() as f64 } else { 0.0 }
            );
        }
    }

    // ---- <out>.params.tsv (ALWAYS) ----------------------------------------------------------------
    // The M2 landing spot this binary lacked: it writes 24 output files and no params certificate, so an
    // ON and an OFF arm of ANY env-driven knob were previously indistinguishable from their outputs.
    // A NEW file changes no existing byte, so writing it unconditionally is compatible with the OFF gate.
    //
    // ⚠ DELIBERATE OMISSION, do not "fix": there is NO `discover_copies` row here. This table is written
    // unconditionally, so adding a row for `--discover-copies` would make `params.tsv` differ between the
    // flag-off and flag-on arms -- breaking that flag's own byte-identity contract (and the regression that
    // pins it, `discover_copies_off_by_default_is_byte_identical` in tests/copy_assign_families.rs, which
    // compares `params.tsv` explicitly). The flag announces itself by the PRESENCE of its own additive
    // output file, `<out>.discovered_copies.tsv`, instead.
    {
        let mut ph = std::fs::File::create(format!("{}.params.tsv", args.out))?;
        writeln!(ph, "key\tvalue")?;
        let mut row = |k: &str, v: String| -> Result<()> {
            writeln!(ph, "{k}\t{v}")?;
            Ok(())
        };
        row("xfam_reconcile", xfam_mode.as_str().to_string())?;
        if args.union_certificate {
            row("union_certificate", "true".to_string())?;
        }
        row("posterior_prior", if prior_abundance { "abundance".into() } else { "uniform".to_string() })?;
        row("margin", format!("{}", args.margin))?;
        row("error_rate", format!("{}", args.error_rate))?;
        row("alpha", format!("{}", args.alpha))?;
        row("margin_gate", format!("{}", args.margin_gate))?;
        row("rna_editing_filter", format!("{}", !args.no_editing_filter))?;
        row("junction_conflict_abstain", format!("{}", args.junction_conflict_abstain))?;
        row("psv_genomic", format!("{}", args.psv_genomic))?;
        row("psv_read_filter", std::env::var("RUSTLE_PSV_READFILTER").unwrap_or_else(|_| "unset".into()))?;
        row("molecule_observations", format!("{}", args.molecule_observations && !args.no_molecule_observations))?;
        row("origin_rejected", format!("{}", assign_rows.iter().filter(|r| r.origin_rejected).count()))?;
        row("orphans", format!("{}", assign_rows.iter().filter(|r| r.origin_rejected && r.n_candidates == 0).count()))?;
        row("origin_substitutions_only", format!("{}", args.origin_substitutions_only))?;
        row("origin_drop_indels", format!("{}", args.origin_drop_indels && !args.no_origin_drop_indels))?;
        row("best_by_duel", format!("{}", args.best_by_duel && !args.no_best_by_duel))?;
        row("gtf_copy_set", format!("{}", (args.gtf_copy_set && !args.no_gtf_copy_set)))?;
        row("min_isoform_fraction", format!("{}", args.min_isoform_fraction))?;
        row("min_boundary_fraction", format!("{}", args.min_boundary_fraction))?;
        row("min_boundary_gap", format!("{}", args.min_boundary_gap))?;
        row("indel_psv", format!("{}", args.indel_psv))?;
        row("indel_psv_min_len", format!("{}", args.indel_psv_min_len))?;
        row("indel_psv_molecules", format!("{}", indel_stats.0))?;
        row("indel_psv_columns", format!("{}", indel_stats.1))?;
        row("indel_psv_columns_ge10", format!("{}", indel_stats.2))?;
        row("admit_aligner_disagreement", format!("{}", args.admit_aligner_disagreement))?;
        row("read_star_junctions", format!("{}", args.read_star_junctions))?;
        row("read_star_genomic", format!("{}", args.read_star_genomic && !args.read_star_unit))?;
        row("read_star_catalog_locus", format!("{}", !args.read_star_pad_locus))?;
        row("sole_candidate", format!("{}", !args.no_sole_candidate))?;
        row("sole_candidates", format!("{}", assign_rows.iter().filter(|r| r.status == "assigned" && r.n_candidates == 1 && r.contested).count()))?;
        row("placement_assign", format!("{}", args.molecule_observations && !args.no_molecule_observations && !args.no_placement_assign))?;
        row("placement_assigned", format!("{}", placement_assigned_total))?;
        row("primary_local_rows", format!("{primary_local_rows}"))?;
        row("placement_first", format!("{}", args.placement_first))?;
        row("readthrough_certificate", format!("{}", !args.no_readthrough_certificate))?;
        row("readthrough_explained", format!("{}", assign_rows.iter().filter(|r| r.readthrough_into != "-").count()))?;
        row("contested_rows", format!("{}", assign_rows.iter().filter(|r| r.contested).count()))?;
        row("dump_star", format!("{}", args.dump_star))?;
        row("read_star_hit_in_unit", format!("{}", !args.no_read_star_hit_in_unit))?;
        row("read_star_two_form", format!("{}", !args.read_star_genomic_only))?;
        row("junction_conflicts", format!("{}", assign_rows.iter().filter(|r| r.junction_conflict).count()))?;
        row("edit_rate", format!("{}", args.edit_rate))?;
        row("iterative_prune", format!("{}", args.iterative_prune))?;
        row("families", args.families.clone().unwrap_or_else(|| "NONE".to_string()))?;
        row("copies_fa", args.copies_fa.clone().unwrap_or_else(|| "NONE".to_string()))?;
        row("dump_psv", format!("{}", args.dump_psv))?;
        row("phase", format!("{}", args.phase))?;
        row("posterior", format!("{}", args.posterior))?;
        row("em", format!("{}", args.em))?;
        row("gtf", format!("{}", args.gtf))?;
        eprintln!("[copy_assign] wrote {}.params.tsv (run certificate)", args.out);
    }

    // Same-locus artifact: two copies of ONE family whose genomic spans OVERLAP are one locus admitted
    // twice, not two copies. Such a family reports min_p == 1 for every read, so it abstains wholesale and
    // its reads masquerade as the K=0 identifiability wall. Warn loudly rather than fail — the catalog is
    // still emitted, but its abstention must not be read as biology. `bench/artifact_audit.py` audits this.
    {
        let catalog: Vec<(String, String, u64, u64)> = quant_rows
            .iter()
            .map(|r| (r.family_id.clone(), r.copy_chrom.clone(), r.copy_start, r.copy_end))
            .collect();
        let flagged = catalog_overlaps(&catalog);
        if !flagged.is_empty() {
            // Two distinct meanings, and they should not be described the same way (verified on RFPL/r4,
            // bench/CONTAINMENT_COVERAGE_FLOOR.md):
            //  - DuplicateLocus (recip ~ 1): one locus admitted twice. Every read scores min_p == 1, so the
            //    family abstains wholesale and its reads masquerade as the K=0 wall.
            //  - Containment (recip << 1): a shorter transcript nested/staggered inside a longer one. On
            //    low-coverage regions this is usually a fragment or a chimeric readthrough, so the COPY COUNT
            //    is inflated — NOT the min_p == 1 masquerade. It cannot be pruned without also deleting
            //    genuine overlapping tandem paralogs (they occupy the same feature cell), so it is reported,
            //    not removed.
            let n_dup = flagged.iter().filter(|f| f.3 == OverlapKind::DuplicateLocus).count();
            let n_contain = flagged.iter().filter(|f| f.3 == OverlapKind::Containment).count();
            let n_shared = flagged.iter().filter(|f| f.3 == OverlapKind::SharedAcrossFamilies).count();
            eprintln!(
                "[copy_assign] WARNING: {} copy pair(s) share genomic sequence \
                 ({n_dup} DuplicateLocus, {n_contain} Containment, {n_shared} SharedAcrossFamilies). \
                 DuplicateLocus = one locus twice, its reads abstain at min_p == 1 (not the K=0 wall). \
                 Containment = a fragment/readthrough nested in a real copy, inflating the copy count on \
                 low-coverage regions — reported, not pruned (it shares its feature cell with real \
                 overlapping paralogs).",
                flagged.len()
            );
            // QUANTIFY the structural warning. A flagged pair is a STRUCTURE ("these two copies share
            // sequence"); what makes it consequential is how many MOLECULES it actually double-claims —
            // the same molecules the cross-family reconciliation's `shared_locus` stratum sees, and the
            // reason that stratum is reported rather than demoted (demoting it would strip these copies
            // of most of their hard support and charge an O1 partition defect to O2's abstention rate).
            // `catalog` is built from `quant_rows` in order, so index i IS quant row i.
            let mut qidx: std::collections::BTreeMap<(&str, usize), usize> = std::collections::BTreeMap::new();
            for (i, r) in quant_rows.iter().enumerate() {
                qidx.insert((r.family_id.as_str(), r.copy_index), i);
            }
            let mut claims: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
            for r in assign_rows.iter().filter(|r| r.status == "assigned") {
                if let Some(&i) = qidx.get(&(r.family_id.as_str(), r.assigned_copy)) {
                    claims.entry(r.read_name.as_str()).or_default().push(i);
                }
            }
            let mut double: std::collections::BTreeMap<(usize, usize), usize> = std::collections::BTreeMap::new();
            for v in claims.values().filter(|v| v.len() >= 2) {
                for a in 0..v.len() {
                    for b in (a + 1)..v.len() {
                        let (lo2, hi2) = (v[a].min(v[b]), v[a].max(v[b]));
                        *double.entry((lo2, hi2)).or_insert(0) += 1;
                    }
                }
            }
            let n_double: usize = flagged
                .iter()
                .map(|&(i, j, _, _)| double.get(&(i.min(j), i.max(j))).copied().unwrap_or(0))
                .sum();
            eprintln!(
                "[copy_assign]   ... those pairs double-claim {n_double} assigned molecule(s) in total \
                 (one molecule counted once per flagged pair it is assigned to both sides of)."
            );
            for &(i, j, recip, kind) in flagged.iter().take(10) {
                let n = double.get(&(i.min(j), i.max(j))).copied().unwrap_or(0);
                eprintln!(
                    "[copy_assign]   {kind:?} recip={recip:.2}  {}/{}:{}-{}  vs  {}/{}-{}  double_claimed_molecules={n}",
                    catalog[i].0, catalog[i].1, catalog[i].2, catalog[i].3, catalog[j].0, catalog[j].2, catalog[j].3
                );
            }
        }
    }

    if uniq > 0 {
        eprintln!(
            "[copy_assign] genome-wide unique-mapper agreement: {agree}/{uniq} ({:.1}%)",
            100.0 * agree as f64 / uniq as f64
        );
    }
    eprintln!("[copy_assign] wrote {0}.families.tsv + {0}.assignments.tsv + {0}.quant.tsv", args.out);
    Ok(())
}

#[cfg(test)]
mod tests {
    /// §6p8: build a tiny GTF and check both polish passes. Layout on chr1/+:
    ///  - LONG   2 introns, 9 reads   (the container)
    ///  - SHORT  1 intron  = LONG's first intron, 2 reads   -> ISM-dropped (2 < 9)
    ///  - STRONG 1 intron  = LONG's first intron, 9 reads   -> kept (support ties the container)
    ///  - MONOHI single exon inside LONG's span, 9 reads    -> kept
    ///  - MONOLO single exon inside LONG's span, 1 read     -> dropped by the ISM pass (host = LONG)
    ///  - FREE   single exon at its own locus, 1 read       -> dropped only by the mono floor
    #[test]
    fn polish_drops_unsupported_fragments_and_bare_mono_loci() {
        fn gtf(tid: &str, reads: u64, exons: &[(i64, i64)]) -> Vec<String> {
            let at = format!("gene_id \"g_{tid}\"; transcript_id \"{tid}\";");
            let mut v = vec![format!(
                "chr1\trustle\ttranscript\t{}\t{}\t.\t+\t.\t{at} reads \"{reads}\";",
                exons[0].0, exons[exons.len() - 1].1
            )];
            for (k, (s, e)) in exons.iter().enumerate() {
                v.push(format!("chr1\trustle\texon\t{s}\t{e}\t.\t+\t.\t{at} exon_number \"{}\";", k + 1));
            }
            v
        }
        let build = || {
            let mut l = Vec::new();
            l.extend(gtf("LONG", 9, &[(100, 200), (300, 400), (500, 600)]));
            l.extend(gtf("SHORT", 2, &[(100, 200), (300, 400)]));
            l.extend(gtf("STRONG", 9, &[(100, 200), (300, 400)]));
            l.extend(gtf("MONOHI", 9, &[(310, 390)]));
            l.extend(gtf("MONOLO", 1, &[(320, 380)]));
            l.extend(gtf("FREE", 1, &[(9000, 9500)]));
            l
        };
        let tids = |l: &[String]| -> Vec<String> {
            l.iter()
                .filter(|x| x.contains("\ttranscript\t"))
                .filter_map(|x| re_attr(x.split('\t').nth(8).unwrap(), "transcript_id"))
                .collect()
        };

        // none is a no-op
        let mut l = build();
        assert_eq!(polish_gtf_lines(&mut l, "none", 0.75, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0), (0, 0, 0, 0, 0));
        assert_eq!(l, build());

        // mono: floor = p75 of {9, 2, 9} = 9, so both 1-read mono transcripts go, MONOHI stays
        let mut l = build();
        let (ism, mono, _, floor, _) = polish_gtf_lines(&mut l, "mono", 0.75, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!((ism, floor), (0, 9));
        assert_eq!(mono, 2);
        assert_eq!(tids(&l), vec!["LONG", "SHORT", "STRONG", "MONOHI"]);

        // full: SHORT is an unsupported sub-chain, MONOLO an unsupported mono inside LONG
        let mut l = build();
        let (ism, mono, _, _, _) = polish_gtf_lines(&mut l, "full", 0.75, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(ism, 2);
        assert_eq!(mono, 1); // FREE has no host, so only the floor removes it
        assert_eq!(tids(&l), vec!["LONG", "STRONG", "MONOHI"]);

        // quantile 0 disables the floor entirely
        let mut l = build();
        let (_, mono, _, floor, _) = polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!((mono, floor), (0, 0));
        assert!(tids(&l).contains(&"FREE".to_string()));
    }

    /// §6p9: the locus isoform fraction removes minor flows and never empties a locus. Two transcripts
    /// share `gene_id "g_LOCUS"`: BIG with 100 reads and TINY with 1 (1% of the locus best).
    #[test]
    fn polish_isoform_fraction_drops_minor_flows_but_keeps_the_dominant() {
        fn gtf(tid: &str, gene: &str, reads: u64, exons: &[(i64, i64)]) -> Vec<String> {
            let at = format!("gene_id \"{gene}\"; transcript_id \"{tid}\";");
            let mut v = vec![format!(
                "chr1\trustle\ttranscript\t{}\t{}\t.\t+\t.\t{at} reads \"{reads}\";",
                exons[0].0, exons[exons.len() - 1].1
            )];
            for (k, (s, e)) in exons.iter().enumerate() {
                v.push(format!("chr1\trustle\texon\t{s}\t{e}\t.\t+\t.\t{at} exon_number \"{}\";", k + 1));
            }
            v
        }
        let build = || {
            let mut l = Vec::new();
            l.extend(gtf("BIG", "g_LOCUS", 100, &[(100, 200), (300, 400), (500, 600)]));
            l.extend(gtf("TINY", "g_LOCUS", 1, &[(100, 200), (350, 400), (500, 600)]));
            l.extend(gtf("SOLO", "g_OTHER", 1, &[(9000, 9100), (9300, 9400)]));
            l
        };
        let tids = |l: &[String]| -> Vec<String> {
            l.iter()
                .filter(|x| x.contains("\ttranscript\t"))
                .filter_map(|x| re_attr(x.split('\t').nth(8).unwrap(), "transcript_id"))
                .collect()
        };
        // TINY is 1% of BIG, so F = 0.02 removes it; SOLO is its own locus's dominant and survives
        let mut l = build();
        let (_, _, frac, _, _) = polish_gtf_lines(&mut l, "full", 0.0, 0.02, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(frac, 1);
        assert_eq!(tids(&l), vec!["BIG", "SOLO"]);
        // F below TINY's share keeps everything
        let mut l = build();
        let (_, _, frac, _, _) = polish_gtf_lines(&mut l, "full", 0.0, 0.005, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(frac, 0);
        assert_eq!(tids(&l), vec!["BIG", "TINY", "SOLO"]);
        // even a huge F never empties a locus: the dominant of each gene_id survives
        let mut l = build();
        polish_gtf_lines(&mut l, "full", 0.0, 0.99, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["BIG", "SOLO"]);
    }

    /// §6q6: fuzzy junction tolerance merges near-duplicate chains into the best-supported member, and
    /// leaves chains that differ by more than the tolerance alone.
    #[test]
    fn polish_fuzzy_junction_merges_near_duplicates() {
        fn gtf(tid: &str, reads: u64, exons: &[(i64, i64)]) -> Vec<String> {
            let at = format!("gene_id \"g_{tid}\"; transcript_id \"{tid}\";");
            let mut v = vec![format!(
                "chr1\trustle\ttranscript\t{}\t{}\t.\t+\t.\t{at} reads \"{reads}\";",
                exons[0].0, exons[exons.len() - 1].1
            )];
            for (k, (s, e)) in exons.iter().enumerate() {
                v.push(format!("chr1\trustle\texon\t{s}\t{e}\t.\t+\t.\t{at} exon_number \"{}\";", k + 1));
            }
            v
        }
        let tids = |l: &[String]| -> Vec<String> {
            l.iter()
                .filter(|x| x.contains("\ttranscript\t"))
                .filter_map(|x| re_attr(x.split('\t').nth(8).unwrap(), "transcript_id"))
                .collect()
        };
        // WOBBLE's acceptor sits 3 bp from BEST's (the NAGNAG distance); FAR's is 9 bp away
        let build = || {
            let mut l = Vec::new();
            l.extend(gtf("BEST", 30, &[(100, 200), (300, 400)]));
            l.extend(gtf("WOBBLE", 5, &[(100, 200), (303, 400)]));
            l.extend(gtf("FAR", 5, &[(100, 200), (309, 400)]));
            l
        };
        // tolerance 0 (the default): nothing merges
        let mut l = build();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l).len(), 3);
        // tolerance 5: WOBBLE folds into BEST (the better-supported member survives), FAR does not
        let mut l = build();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 5, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["BEST", "FAR"]);
        // tolerance 10: FAR folds in too
        let mut l = build();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 10, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["BEST"]);
    }

    /// §6q0/§6q1: the shadow rule drops single-exon transcripts in a spliced gene's shadow (exon overlap
    /// on EITHER strand, or same-strand span overlap) but not one whose only overlap is an anti-strand
    /// SPAN; and the ISM escape keeps a well-supported fragment its container would otherwise absorb.
    #[test]
    #[test]
    /// §6za: a chain whose exon contains another transcript's junction carrying >= ratio x its reads is
    /// dropped; the same chain survives when the spanning junction is not dominant enough, and at ratio 0.
    fn polish_retained_intron_drops_the_read_through_and_keeps_the_minor_isoform() {
        let mk = |tid: &str, reads: u32, exons: &[(u32, u32)]| -> Vec<String> {
            let mut v = vec![format!(
                "c\tr\ttranscript\t{}\t{}\t.\t+\t.\tgene_id \"G\"; transcript_id \"{tid}\"; reads \"{reads}\";",
                exons[0].0, exons[exons.len() - 1].1
            )];
            for (a, b) in exons {
                v.push(format!("c\tr\texon\t{a}\t{b}\t.\t+\t.\tgene_id \"G\"; transcript_id \"{tid}\";"));
            }
            v
        };
        // SPLICED: 100-200, 300-400, 500-600, 700-800 with 50 reads; RETAIN: 100-200, 300-600, 700-800 with 4
        // reads — its middle exon contains SPLICED's 401-499 intron, and its chain is NOT a contiguous
        // sub-chain of SPLICED's (it skips one junction), so the ISM pass leaves it to this rule.
        let base = || {
            let mut l = mk("SPLICED", 50, &[(100, 200), (300, 400), (500, 600), (700, 800)]);
            l.extend(mk("RETAIN", 4, &[(100, 200), (300, 600), (700, 800)]));
            l
        };
        let mut l = base();
        let (_, _, _, _, n_ret) = polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 10.0);
        assert_eq!(n_ret, 1, "50 >= 10 x 4: the read-through chain is dropped");
        assert!(!l.iter().any(|x| x.contains("\"RETAIN\"")) && l.iter().any(|x| x.contains("\"SPLICED\"")));
        let mut l = base();
        let (_, _, _, _, n_ret) = polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 20.0);
        assert_eq!(n_ret, 0, "50 < 20 x 4: a minor isoform with real share survives");
        let mut l = base();
        let (_, _, _, _, n_ret) = polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(n_ret, 0, "ratio 0 is off");
    }

    #[test]
    fn polish_shadow_and_ism_escape() {
        fn gtf(tid: &str, strand: &str, reads: u64, exons: &[(i64, i64)]) -> Vec<String> {
            let at = format!("gene_id \"g_{tid}\"; transcript_id \"{tid}\";");
            let mut v = vec![format!(
                "chr1\trustle\ttranscript\t{}\t{}\t.\t{strand}\t.\t{at} reads \"{reads}\";",
                exons[0].0, exons[exons.len() - 1].1
            )];
            for (k, (s, e)) in exons.iter().enumerate() {
                v.push(format!("chr1\trustle\texon\t{s}\t{e}\t.\t{strand}\t.\t{at} exon_number \"{}\";", k + 1));
            }
            v
        }
        let tids = |l: &[String]| -> Vec<String> {
            l.iter()
                .filter(|x| x.contains("\ttranscript\t"))
                .filter_map(|x| re_attr(x.split('\t').nth(8).unwrap(), "transcript_id"))
                .collect()
        };
        // PLUS spans 100..600 on '+' with exons 100-200/300-400/500-600 (introns 201-299, 401-499)
        let shadow = || {
            let mut l = Vec::new();
            l.extend(gtf("PLUS", "+", 20, &[(100, 200), (300, 400), (500, 600)]));
            l.extend(gtf("SAMEEX", "+", 20, &[(150, 190)]));  // same-strand exon  -> drop
            l.extend(gtf("ANTIEX", "-", 20, &[(150, 190)]));  // anti-strand exon  -> drop
            l.extend(gtf("SAMEIN", "+", 20, &[(220, 280)]));  // same-strand intron (span) -> drop
            l.extend(gtf("ANTIIN", "-", 20, &[(220, 280)]));  // anti-strand SPAN only -> KEEP
            l.extend(gtf("FAR", "+", 20, &[(9000, 9500)]));   // no overlap -> KEEP
            l
        };
        let mut l = shadow();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, true, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["PLUS", "ANTIIN", "FAR"]);
        // shadow off leaves them all
        let mut l = shadow();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l).len(), 6);

        // ISM escape: FRAG's chain is a sub-chain of DEEP's; 10 reads is below DEEP's 100 but reaches the
        // support level of the multi-exon set {100, 10, 10}, so at quantile 0 it escapes
        let ism = || {
            let mut l = Vec::new();
            l.extend(gtf("DEEP", "+", 100, &[(100, 200), (300, 400), (500, 600)]));
            l.extend(gtf("FRAG", "+", 10, &[(100, 200), (300, 400)]));
            l.extend(gtf("OTHER", "-", 10, &[(9000, 9100), (9300, 9400)]));
            l
        };
        let mut l = ism();
        polish_gtf_lines(&mut l, "full", 0.0, 0.0, false, false, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["DEEP", "OTHER"], "without the escape the fragment is absorbed");
        let mut l = ism();
        polish_gtf_lines(&mut l, "full", 0.10, 0.0, false, true, false, 1.0, false, 0, false, 0, 0.0);
        assert_eq!(tids(&l), vec!["DEEP", "FRAG", "OTHER"], "with the escape a well-supported fragment survives");
    }

    #[test]
    fn lift_blocks_map_both_strands_and_inverse() {
        // '+': query 0..10 aligned to target 100..110 with a 2-bp query insertion after 4 and a 3-bp deletion after 7
        let f = LiftBlocks::from_cigar(0, 12, 100, false, "4=2I3=3D3=");
        assert_eq!(f.map(2), Some((102, 0)));
        assert_eq!(f.map(6), Some((104, 0)), "after the insertion the query is 2 ahead");
        assert_eq!(f.map(9), Some((110, 0)), "after the deletion the target is 3 ahead");
        assert_eq!(f.map(4), Some((104, 1)), "inside the insertion: nearest edge (q=3 -> t=103) plus the offset, distance 1");
        let inv = f.inverse();
        assert_eq!(inv.map(102), Some((2, 0)));
        assert_eq!(inv.map(110), Some((9, 0)));
        // '-': query span [0,10) reverse-complemented onto target 200..210: query 9 <-> target 200, query 0 <-> target 209
        let r = LiftBlocks::from_cigar(0, 10, 200, true, "10=");
        assert_eq!(r.map(9), Some((200, 0)));
        assert_eq!(r.map(0), Some((209, 0)));
        let rinv = r.inverse();
        assert_eq!(rinv.map(200), Some((9, 0)));
        assert_eq!(rinv.map(209), Some((0, 0)));
        assert_eq!(rinv.map(205), Some((4, 0)));
        // round trip on the '+' case
        for q in [0u64, 3, 7, 9, 11] {
            if let Some((t, 0)) = f.map(q) { assert_eq!(inv.map(t), Some((q, 0)), "round trip at {q}"); }
        }
    }

    use super::*;

    #[test]
    fn block_overlap_ignores_an_intron_that_merely_spans_the_window() {
        // ref_start=0, "10M5000N10M": aligned blocks are [0,10) and [5010,5020); the intron covers
        // [10,5010) with no M/=/X inside it. A window fully inside the intron (e.g. [100,200)) must
        // score 0 overlap, even though it lies strictly between the read's ref_start and ref_end --
        // this is the exact bug `read_ref_end_local`-based span checks were vulnerable to.
        let read = rustle::vg_family::copy_split::AlignedRead {
            ref_start: 0,
            cigar: vec![('M', 10), ('N', 5000), ('M', 10)],
            seq: vec![],
            qual: vec![],
        };
        assert_eq!(block_overlap(&read, 100, 200), 0, "window sits entirely inside the spliced-out intron");
        // sanity: a window over the trailing M block still scores correctly.
        assert_eq!(block_overlap(&read, 5010, 5020), 10, "window exactly covers the second aligned block");
    }

    #[test]
    fn block_overlap_counts_the_real_match_run_inside_the_window() {
        // ref_start=100, "50M": a simple hand-computable case -- window [110,130) is fully inside the
        // aligned block [100,150), so overlap is the window's own width, 20.
        let read = rustle::vg_family::copy_split::AlignedRead {
            ref_start: 100,
            cigar: vec![('M', 50)],
            seq: vec![],
            qual: vec![],
        };
        assert_eq!(block_overlap(&read, 110, 130), 20);
        // partial overlap at the trailing edge: window [140,160) vs aligned block ending at 150 -> 10.
        assert_eq!(block_overlap(&read, 140, 160), 10);
        // window entirely outside the block -> 0.
        assert_eq!(block_overlap(&read, 200, 210), 0);
    }

    #[test]
    fn best_overlap_truth_copy_picks_the_copy_with_more_overlap() {
        // A read whose primary alignment ("M",100 from ref_start 0) overlaps copy A's span [0,80) by 80bp
        // and copy B's span [60,100) by only 40bp -- the truth copy must be A ("best overlap wins").
        let read = rustle::vg_family::copy_split::AlignedRead { ref_start: 0, cigar: vec![('M', 100)], seq: vec![], qual: vec![] };
        let br = BamRead {
            chrom: "chr1".to_string(), read, mapq: 0, name: "r1".to_string(), as_score: 0, de: 0.0,
            is_supplementary: false, is_secondary: false, reverse: false, ts: None,
        };
        let copy_spans = vec![("chr1".to_string(), 0u64, 80u64), ("chr1".to_string(), 60u64, 100u64)];
        let copy_tids = vec!["tidA".to_string(), "tidB".to_string()];
        let mut catalog_index: CatalogIndex = std::collections::HashMap::new();
        catalog_index.insert("tidA".to_string(), ("famA".to_string(), 0usize));
        catalog_index.insert("tidB".to_string(), ("famB".to_string(), 1usize));
        let truth = best_overlap_truth_copy(std::slice::from_ref(&br), &copy_spans, &copy_tids, Some(&catalog_index));
        assert_eq!(
            truth.get("r1"),
            Some(&(("famA".to_string(), "0".to_string()), 80)),
            "copy A (80bp overlap) beats copy B (40bp); return value now carries (cf, cidx), not bare cidx"
        );
    }

    #[test]
    fn best_overlap_truth_copy_ties_keep_the_first_seen_candidate() {
        // Two candidate copies with EQUAL overlap (50bp each): the first one in `copy_spans`' iteration
        // order wins, matching Python's strict `>` compare over `cp.items()`'s insertion order -- a later
        // equal-overlap candidate never displaces it.
        let read = rustle::vg_family::copy_split::AlignedRead { ref_start: 0, cigar: vec![('M', 100)], seq: vec![], qual: vec![] };
        let br = BamRead {
            chrom: "chr1".to_string(), read, mapq: 0, name: "r1".to_string(), as_score: 0, de: 0.0,
            is_supplementary: false, is_secondary: false, reverse: false, ts: None,
        };
        let copy_spans = vec![("chr1".to_string(), 0u64, 50u64), ("chr1".to_string(), 50u64, 100u64)];
        let copy_tids = vec!["tidA".to_string(), "tidB".to_string()];
        let mut catalog_index: CatalogIndex = std::collections::HashMap::new();
        catalog_index.insert("tidA".to_string(), ("famA".to_string(), 0usize));
        catalog_index.insert("tidB".to_string(), ("famB".to_string(), 1usize));
        let truth = best_overlap_truth_copy(std::slice::from_ref(&br), &copy_spans, &copy_tids, Some(&catalog_index));
        assert_eq!(
            truth.get("r1"),
            Some(&(("famA".to_string(), "0".to_string()), 50)),
            "copy A (seen first) wins the tie over copy B; return value now carries (cf, cidx), not bare cidx"
        );
    }

    #[test]
    fn discover_copies_is_attributed_only_to_the_family_that_considered_the_read() {
        // ⚠ THE CROSS-FAMILY POOLING BUG (final whole-branch review, Critical). Two families, A and B,
        // in the SAME region. One AS-tied read (`tied`) is in family A's `assignments` and NOT in B's.
        // Its two max-AS placements are 1000-1100 (inside A's own catalog copy) and 5000-5100 (out of
        // catalog) -- so A should report the 5000-5100 site and B should report NOTHING AT ALL, because
        // that read was never B's to reason about. Pre-fix, the whole region's tied list was handed to
        // every family, so the identical site with the identical read list came out under both ids.
        use rustle::vg_family::copy_assign::Assignment;
        use rustle::vg_family::copy_discovery::tie_partner_placements;
        let mk = |name: &str, start: u64, as_score: i32| BamRead {
            chrom: "chr1".to_string(),
            read: rustle::vg_family::copy_split::AlignedRead {
                ref_start: start, cigar: vec![('M', 100)], seq: vec![], qual: vec![],
            },
            mapq: 0, name: name.to_string(), as_score, de: 0.0,
            is_supplementary: false, is_secondary: start != 1000, reverse: false, ts: None,
        };
        // Two reads, each AS-tied across two placements; only `tied_a` belongs to family A.
        let bam_reads = vec![
            mk("tied_a", 1000, 200), mk("tied_a", 5000, 200),
            mk("tied_b", 2000, 300), mk("tied_b", 7000, 300),
        ];
        let assign = || Assignment {
            best_copy: 0, log_lr_margin: 0.0, n_decisive: 0, resolvable: false,
            status: AssignStatus::Tied, p_value: 1.0, min_p_value: 1.0, discovery_coupled: false,
            junction_conflict: false, origin_rejected: false, n_candidates: 0,
            posterior: vec![1.0], sibling_identity: 1.0, n_cols_vs_nearest_sibling: 0,
        };
        let mut fam_a = FamilyAssignment::empty();
        fam_a.family_id = "FAM_A".to_string();
        fam_a.copy_tids = vec!["tidA".to_string()];
        fam_a.copy_spans = vec![("chr1".to_string(), 900, 1200)]; // contains the 1000-1100 placement
        fam_a.assignments = vec![(0, assign()), (1, assign())]; // indices of tied_a's two records
        let mut fam_b = FamilyAssignment::empty();
        fam_b.family_id = "FAM_B".to_string();
        fam_b.copy_tids = vec!["tidB".to_string()];
        fam_b.copy_spans = vec![("chr1".to_string(), 1900, 2200)]; // contains the 2000-2100 placement
        fam_b.assignments = vec![(2, assign()), (3, assign())]; // indices of tied_b's two records

        let tied = tie_partner_placements(&bam_reads);
        assert_eq!(tied.len(), 2, "both reads are AS-tied across two placements each");

        // min_support is 2, so a single read cannot clear it -- give each family's own read a second,
        // co-located supporter that the OTHER family still never considered.
        let bam_reads = {
            let mut v = bam_reads;
            v.push(mk("tied_a2", 1000, 200));
            v.push(mk("tied_a2", 5050, 200));
            v.push(mk("tied_b2", 2000, 300));
            v.push(mk("tied_b2", 7050, 300));
            v
        };
        fam_a.assignments.extend([(4, assign()), (5, assign())]);
        fam_b.assignments.extend([(6, assign()), (7, assign())]);
        let tied = tie_partner_placements(&bam_reads);

        let a = discover_copies_for_family(&fam_a, &bam_reads, &tied);
        let b = discover_copies_for_family(&fam_b, &bam_reads, &tied);

        assert_eq!(a.len(), 1, "family A must report exactly its own out-of-catalog site: {a:#?}");
        assert_eq!(a[0].family_id, "FAM_A");
        assert_eq!((a[0].chrom.as_str(), a[0].start, a[0].end), ("chr1", 5000, 5150));
        assert_eq!(a[0].read_names, vec!["tied_a".to_string(), "tied_a2".to_string()]);
        assert_eq!(a[0].nearest_copy_tid, "tidA");
        // and the 7000 site, which only family B's reads support, must NOT appear under A:
        assert!(!a.iter().any(|d| d.start >= 7000), "family A must not inherit family B's reads: {a:#?}");

        assert_eq!(b.len(), 1, "family B likewise reports only its own site: {b:#?}");
        assert_eq!(b[0].family_id, "FAM_B");
        assert_eq!((b[0].chrom.as_str(), b[0].start, b[0].end), ("chr1", 7000, 7150));
        assert!(!b.iter().any(|d| d.start == 5000), "family B must not inherit family A's reads: {b:#?}");

        // The second half of the same fix: `existing_copies` comes from `fa.copy_spans`/`copy_tids`, so a
        // family with NO catalog copy at the tie's own position still excludes its own copies -- here,
        // stripping A's copy set makes the 1000-1100 placement surface as a candidate too (proof the
        // exclusion list is really being read from the family, not silently empty).
        let mut fam_a_no_copies = fam_a.clone();
        fam_a_no_copies.copy_spans.clear();
        fam_a_no_copies.copy_tids.clear();
        let a2 = discover_copies_for_family(&fam_a_no_copies, &bam_reads, &tied);
        assert_eq!(a2.len(), 2, "with no catalog copies both of A's tied sites are out-of-catalog: {a2:#?}");
        assert_eq!(a2[0].nearest_copy_tid, "NA");
        assert_eq!(opt_u64(a2[0].nearest_copy_distance), "NA", "no copy at all -> NA distance, not a sentinel");
    }

    #[test]
    fn overlapping_regions_on_the_same_contig_are_rejected() {
        let mut by_contig = std::collections::BTreeMap::new();
        by_contig.insert("chr1".to_string(), vec![(100, 200), (150, 250)]);
        let err = validate_no_overlapping_regions(&by_contig).unwrap_err();
        assert!(err.to_string().contains("chr1:100-200"), "{err}");
        assert!(err.to_string().contains("chr1:150-250"), "{err}");
    }

    #[test]
    fn touching_regions_are_not_an_overlap() {
        let mut by_contig = std::collections::BTreeMap::new();
        by_contig.insert("chr1".to_string(), vec![(100, 200), (200, 300)]);
        assert!(validate_no_overlapping_regions(&by_contig).is_ok());
    }

    #[test]
    fn disjoint_regions_on_different_contigs_are_fine_even_if_the_coordinates_overlap() {
        let mut by_contig = std::collections::BTreeMap::new();
        by_contig.insert("chr1".to_string(), vec![(100, 200)]);
        by_contig.insert("chr2".to_string(), vec![(100, 200)]);
        assert!(validate_no_overlapping_regions(&by_contig).is_ok());
    }

    #[test]
    fn linearize_tsv_row_formats() {
        use rustle::vg_family::linearize::{LinearizeCertificate, Verdict};
        let c = LinearizeCertificate { n_pool: 40, linearized_frac_real: 0.82, mean_frac_decoy: 0.01,
            delta: 0.81, perm_p: 0.05, verdict: Verdict::Linearizes };
        let row = linearize_tsv_row("GWFAM1", ("chr9", 100, 200), &c);
        assert_eq!(row, "GWFAM1\tchr9\t100\t200\t40\t0.820\t0.010\t0.810\t0.0500\tLINEARIZES");
    }

    #[test]
    fn read_lambda_file_parses_the_scalar() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_test_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\n25.5\t1234\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), Some(25.5));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn read_lambda_file_none_on_na_or_missing() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_na_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\nNA\t0\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), None);
        assert_eq!(read_lambda_file("/nonexistent/path.tsv"), None);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn resolve_lambda_precedence_explicit_over_file() {
        assert_eq!(resolve_lambda(Some(30.0), Some(25.0)), Some(30.0));
        assert_eq!(resolve_lambda(None, Some(25.0)), Some(25.0));
        assert_eq!(resolve_lambda(None, None), None);
    }

    #[test]
    fn anchored_support_counts_only_unique_mappers_of_the_copy() {
        // reads: (best_copy, mapq) = (0,60) (0,0) (1,60) (0,60). Copy 0's unique-mapper support is the two
        // mapq>0 reads (the mapq=0 read is a tie-break-arbitrary primary and does not count).
        let bcs = [0usize, 0, 1, 0];
        let mqs = [60u8, 0, 60, 60];
        assert_eq!(anchored_support(&bcs, &mqs, 0), 2, "copy 0: two mapq>0 reads, the mapq=0 excluded");
        assert_eq!(anchored_support(&bcs, &mqs, 1), 1, "copy 1: one mapq>0 read");
        assert_eq!(anchored_support(&bcs, &mqs, 2), 0, "copy 2: no reads");
    }

    #[test]
    fn anchored_support_is_zero_when_every_primary_is_tied() {
        // The TSPY case: every assigned read is mapq=0 -> no copy has any anchored support -> none invariant.
        let bcs = [0usize, 1, 2, 0, 1];
        let mqs = [0u8, 0, 0, 0, 0];
        for ci in 0..3 {
            assert_eq!(anchored_support(&bcs, &mqs, ci), 0);
        }
    }

    #[test]
    fn tie_invariant_threshold_is_the_locus_gate() {
        use rustle::vg_family::denovo_assemble::GATE_MIN_READS;
        // The certificate boolean is anchored >= GATE_MIN_READS (=3): 2 -> false, 3 -> true.
        assert!(!(2u32 >= GATE_MIN_READS));
        assert!(3u32 >= GATE_MIN_READS);
    }

    #[test]
    fn build_copy_graph_maps_family_to_graph() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100), Some(200)];
        fa.copy_psv_alleles = vec![vec![Some(b'A'), Some(b'A')], vec![Some(b'A'), Some(b'G')]];
        fa.read_psv_obs = vec![];
        fa.assignments = vec![];
        // injected reference: A at both positions
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &[], None, &eff);
        assert_eq!(g.columns.len(), 2);
        assert_eq!(g.copies.len(), 2);
        assert_eq!(g.columns[0].ref_allele, Some(b'A'));
        // graph renders the reference walk + a divergent copy
        let gfa = g.to_gfa();
        assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
        assert!(gfa.contains("CAFAM0_c1_G+"));
    }

    // A discovery_coupled Assignment for read `ri` pinned to copy `best_copy` (status Assigned) — the ONLY
    // signal that makes a copy absent. Mirrors the default Assignment (copy_assign.rs) with the flag flipped.
    fn coupled_assignment(ri: usize, best_copy: usize) -> (usize, rustle::vg_family::copy_assign::Assignment) {
        (ri, rustle::vg_family::copy_assign::Assignment {
            best_copy,
            log_lr_margin: 10.0,
            n_decisive: 1,
            resolvable: true,
            status: AssignStatus::Assigned,
            p_value: 0.0,
            min_p_value: 0.0,
            discovery_coupled: true,
            junction_conflict: false,
            origin_rejected: false,
            n_candidates: 0,
            posterior: vec![],
            sibling_identity: 1.0,
            n_cols_vs_nearest_sibling: 0,
        })
    }

    #[test]
    fn build_copy_graph_coupled_copy_colocated_is_absent_collapsed() {
        // Absence is driven by a discovery_coupled read, NOT the collapsed/rescued counts. copy1 has a
        // coupled read AND its span overlaps copy0's span => AbsentCollapsed (hidden co-located haplotype).
        // copy0 (no coupled read) stays AnnotationUnknown (no _ABSENT).
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        // copy0 = reference allele A; copy1 = divergent G (the absent, co-located copy).
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // overlapping spans on the same chrom => co-located.
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr1".into(), 1500, 2500)];
        fa.read_psv_obs = vec![vec![Some(b'G')]];
        fa.assignments = vec![coupled_assignment(0, 1)]; // one read, coupled to copy1
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &["read0".to_string()], None, &eff);
        assert_eq!(g.copies.len(), 2);
        assert_eq!(g.copies[1].status, CopyStatus::AbsentCollapsed);
        assert!(g.copies[1].status.is_absent());
        let gfa = g.to_gfa();
        assert!(gfa.contains("_copy1_ABSENT"), "coupled co-located copy must render _ABSENT:\n{}", gfa);
        assert!(gfa.contains("ST:Z:absent-collapsed"), "expected ST:Z:absent-collapsed in:\n{}", gfa);
        // copy 0 is NOT absent (no coupled read) => AnnotationUnknown, no _ABSENT.
        assert_eq!(g.copies[0].status, CopyStatus::AnnotationUnknown);
        assert!(!g.copies[0].status.is_absent());
        assert!(!gfa.contains("_copy0_ABSENT"), "non-coupled copy must NOT be _ABSENT:\n{}", gfa);
    }

    #[test]
    fn build_copy_graph_coupled_copy_dispersed_is_absent_divergent() {
        // copy1 has a coupled read but its span is DISJOINT from copy0's (different chrom) => AbsentDivergent
        // (dispersed, no overlapping in-genome copy).
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // disjoint spans (different chrom) => dispersed.
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr9".into(), 1000, 2000)];
        fa.read_psv_obs = vec![vec![Some(b'G')]];
        fa.assignments = vec![coupled_assignment(0, 1)];
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &["read0".to_string()], None, &eff);
        assert_eq!(g.copies[1].status, CopyStatus::AbsentDivergent);
        let gfa = g.to_gfa();
        assert!(gfa.contains("_copy1_ABSENT"), "dispersed coupled copy must render _ABSENT:\n{}", gfa);
        assert!(gfa.contains("ST:Z:absent-divergent"), "expected ST:Z:absent-divergent in:\n{}", gfa);
        assert_eq!(g.copies[0].status, CopyStatus::AnnotationUnknown);
    }

    #[test]
    fn build_copy_graph_gstm_regression_no_coupled_reads_no_absent() {
        // REGRESSION for the GSTM bug: collapsed_copies (9) >> n_copies (3) with NO discovery_coupled reads.
        // The old code let absent_tail_start underflow to 0 and mislabeled ALL 3 in-genome copies _ABSENT.
        // Correct behavior: no coupled read => NO copy is absent; all are AnnotationUnknown, no _ABSENT.
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 3;
        fa.copy_tids = vec!["c0".into(), "c1".into(), "c2".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![Some(b'G')]];
        fa.read_psv_obs = vec![];
        fa.assignments = vec![]; // NO discovery_coupled reads (no --absent-copies)
        fa.collapsed_copies = 9;  // diagnostic count, far exceeds n_copies — must NOT force absence.
        fa.rescued_copies = 0;
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &[], None, &eff);
        assert_eq!(g.copies.len(), 3);
        for ci in 0..3 {
            assert_eq!(g.copies[ci].status, CopyStatus::AnnotationUnknown, "copy {} must be in-genome", ci);
            assert!(!g.copies[ci].status.is_absent(), "copy {} must NOT be absent", ci);
        }
        let gfa = g.to_gfa();
        assert!(!gfa.contains("_ABSENT"), "no copy may render _ABSENT when no read is discovery_coupled:\n{}", gfa);
        assert!(!gfa.contains("absent-collapsed"), "no absent-collapsed tag expected:\n{}", gfa);
    }

    #[test]
    fn build_copy_graph_fills_mi_from_copy_map_identity() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // copy1 has a remap identity but NO discovery_coupled read: copy_map_identity feeds the MI tag
        // ONLY and must NOT flip absent status — so copy1 stays IN-GENOME (AnnotationUnknown, no --gff).
        fa.copy_map_identity = vec![None, Some(0.952)];
        fa.assignments = vec![]; // no discovery_coupled reads => no absence
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, |_c, _p| Some(b'A'), &[], None, &eff);
        // copy1 is in-genome, NOT absent (the LOCK on the copy_map_identity-drives-absence bug).
        assert_eq!(g.copies[1].status, CopyStatus::AnnotationUnknown, "copy_map_identity alone must NOT make a copy absent");
        assert!(!g.copies[1].status.is_absent());
        let gfa = g.to_gfa();
        // P-line name is EXACTLY CAFAM0_copy1 (the trailing tab excludes CAFAM0_copy1_ABSENT).
        let c1 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy1\t")).expect("copy1 must be P\\tCAFAM0_copy1 (NOT _ABSENT)");
        assert!(!gfa.contains("CAFAM0_copy1_ABSENT"), "copy_map_identity must not render _ABSENT:\n{}", gfa);
        assert!(c1.contains("MI:f:0.952"), "copy1 MI missing: {}", c1);
        assert!(c1.contains("ST:Z:annotation-unknown"), "copy1 must carry a non-absent ST:Z: status: {}", c1);
        let c0 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy0\t")).unwrap();
        assert!(!c0.contains("MI:f:"), "copy0 must omit MI: {}", c0);
    }

    #[test]
    fn build_exon_graph_makes_copy_specific_arm() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.copy_spans = vec![("chr1".into(), 0, 400), ("chr1".into(), 0, 400)];
        // copy0 introns skip 100-300 (exons 0-100, 300-400); copy1 has an extra exon (exons 0-100,150-250,300-400)
        fa.copy_introns = vec![ vec![(100,300)], vec![(100,150),(250,300)] ];
        fa.copy_map_identity = vec![None, Some(0.95)];
        // copy1 is reference-ABSENT via a discovery_coupled read (the v1 absence mechanism — NOT
        // copy_map_identity, which feeds only the MI tag); its span overlaps copy0 => AbsentCollapsed.
        fa.assignments = vec![coupled_assignment(0, 1)];
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_exon_graph("CAFAM0", &fa, None, &eff);
        // copy1 walks one more class than copy0
        assert!(g.copies[1].exon_nodes.len() > g.copies[0].exon_nodes.len());
        let gfa = g.to_gfa(|ec| vec![b'A'; (ec.end-ec.start) as usize]);
        assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
        assert!(gfa.lines().any(|l| l.starts_with("P\tCAFAM0_copy1_ABSENT")));
    }

    #[test]
    fn annotation_axis_from_intervals() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr1".into(), 5000, 6000)];
        fa.psv_col_pos = vec![Some(1500)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'A')]];
        // annotation covers only copy0's span
        let ann = vec![("chr1".to_string(), 900u64, 2100u64)];
        let s0 = annotation_status(&fa, 0, Some(&ann));
        let s1 = annotation_status(&fa, 1, Some(&ann));
        assert_eq!(s0, CopyStatus::InGenomeAnnotated);
        assert_eq!(s1, CopyStatus::InGenomeUnannotated);
        assert_eq!(annotation_status(&fa, 1, None), CopyStatus::AnnotationUnknown);
    }

    #[test]
    fn parse_annotation_bed_with_numeric_name() {
        // Regression: a BED6 line whose NAME column (col 3) is numeric must NOT be misread as GFF (which would
        // yield (4,5) from cols 3/4). Extension-first dispatch (`.bed` => BED cols 0/1/2) keeps it correct.
        let dir = std::env::temp_dir();
        let bed = dir.join(format!("rustle_task8_ann_{}.bed", std::process::id()));
        std::fs::write(&bed, "chr1\t1000\t2000\t5\t0\t+\n").unwrap();
        let got = parse_annotation(bed.to_str().unwrap()).unwrap();
        std::fs::remove_file(&bed).ok();
        assert_eq!(got, vec![("chr1".to_string(), 1000u64, 2000u64)], "BED numeric-name line must NOT parse as (4,5)");

        // A GFF/GTF file still parses cols 3/4 with the 1-based -> 0-based start conversion.
        let gff = dir.join(format!("rustle_task8_ann_{}.gff3", std::process::id()));
        std::fs::write(&gff, "##gff-version 3\nchr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=g1\n").unwrap();
        let got = parse_annotation(gff.to_str().unwrap()).unwrap();
        std::fs::remove_file(&gff).ok();
        assert_eq!(got, vec![("chr1".to_string(), 999u64, 2000u64)], "GFF cols 3/4, 1-based start -> 0-based");
    }
}
