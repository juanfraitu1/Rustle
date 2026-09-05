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
    build_clusters, fold_parts_into_loci, graph_from_paf_loci, loci_from_exon_blocks, mcl, sd_blocks, refine_cluster_cores, Cluster, CoreStatus,
    GeneKey, GraphParams, SdPairs,
};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

#[derive(Parser, Debug)]
#[command(about = "Multi-copy gene families by MCL over annotation sequence, corroborated by RNA")]
struct Args {
    /// All-vs-all PAF of the annotated gene sequences (headers `CONTIG:START-END`, GFF 1-based).
    #[arg(long)]
    paf: String,

    /// GFF supplying each gene's EXON-UNION length. ⚠ Without it the coverage denominator falls back to
    /// the genomic span, which removed 62.8% of eligible genes in the pilot (§6dc) — the run warns loudly.
    #[arg(long)]
    gff: Option<String>,

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
    #[arg(long, default_value_t = false)]
    emit_units: bool,
    /// Escape hatch: do NOT emit units even though `--bam` and `--fasta` are given (§6er: units are emitted
    /// whenever both are supplied).
    #[arg(long, default_value_t = false)]
    no_emit_units: bool,
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
    let (mut blocks, mut introns): (Vec<(u64, u64)>, Vec<(u64, u64)>) = (Vec::new(), Vec::new());
    let mut p = br.read.ref_start;
    let mut cur: Option<(u64, u64)> = None;
    for &(op, n) in &br.read.cigar {
        match op {
            'M' | '=' | 'X' | 'D' => {
                cur = Some((cur.map_or(p, |c| c.0), p + n));
                p += n;
            }
            'N' => {
                if let Some(c) = cur.take() {
                    blocks.push(c);
                }
                introns.push((p, p + n));
                p += n;
            }
            _ => {}
        }
    }
    if let Some(c) = cur {
        blocks.push(c);
    }
    (blocks, introns)
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
    const GIANT_BP: u64 = 50_000;
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
        let mut segs: Vec<Vec<(u64, u64)>> = vec![vec![blocks[0]]];
        for (k, intr) in introns.iter().enumerate() {
            if intr.1 - intr.0 > GIANT_BP && support.get(intr).copied().unwrap_or(0) < min_reads {
                segs.push(Vec::new());
            }
            if k + 1 < blocks.len() {
                segs.last_mut().unwrap().push(blocks[k + 1]);
            }
        }
        kept.push(segs.into_iter().find(|sg| sg.iter().any(|&(s, e)| e > lo && s < hi)));
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

fn main() -> Result<()> {
    let mut args = Args::parse();
    // §6er (S2): the unit is the catalog row. Stages engage on their inputs; escape hatches reproduce the
    // record-level catalogs (`--no-merge-overlapping-loci`, `--no-core-refine`, `--no-emit-units`).
    if args.no_merge_overlapping_loci {
        args.merge_overlapping_loci = false;
    }
    if args.sedef.is_some() && !args.no_core_refine {
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
    if args.bam.is_some() && args.fasta.is_some() && !args.no_emit_units {
        args.emit_units = true;
    }
    let p = GraphParams {
        min_identity: args.min_identity,
        min_cov_longer: args.min_cov_longer,
        min_bp: args.min_bp,
        exonic_overlap: args.exonic_overlap,
        reject_overlapping: args.reject_overlapping,
        min_exonic_bp: args.min_exonic_bp,
        exonic_both_sides: args.exonic_both_sides,
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

    // ⭐ Duplicon-first core refinement (§6eh). Post-MCL, per cluster; clusters.tsv above is untouched.
    let mut core_stats = (0usize, 0usize, 0usize, 0usize, 0usize); // gated clusters, kept-full, trimmed, dropped, untouched clusters
    let mut core_records: Vec<Vec<rustle::vg_family::annotation_families::CoreRecord>> = Vec::new();
    if args.core_refine {
        let bed = args.sedef.as_ref().ok_or_else(|| anyhow::anyhow!("--core-refine needs --sedef <bed>"))?;
        let text = std::fs::read_to_string(bed).with_context(|| format!("reading {bed}"))?;
        let sd = SdPairs::from_bed_str(&text);
        eprintln!("[mcl_families] core-refine: {} SEDEF pair(s) loaded from {bed}", sd.n_pairs());
        let mut cf = std::fs::File::create(format!("{}.cores.tsv", args.out))?;
        writeln!(cf, "cluster_id\tmember\tgate\tmax_depth\tcore_bp\tspan\tmedian_core\tstatus\tcore_hull")?;
        let mut rf = std::fs::File::create(format!("{}.refined.clusters.tsv", args.out))?;
        writeln!(rf, "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\tstatus")?;
        for (i, c) in clusters.iter().enumerate() {
            let recs = refine_cluster_cores(&c.members, &sd);
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
    }
    if args.emit_units {
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
        writeln!(ut, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tsource\tcore_hull\tsd_depth\tcore_bp\tnearest_ident\trep_frac")?;
        // curated repeats (optional --rmsk): per contig, sorted interspersed intervals
        let rmsk: BTreeMap<String, Vec<(u64, u64)>> = match &args.rmsk {
            Some(path) => {
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
                m
            }
            None => BTreeMap::new(),
        };
        let rep_frac = |chrom: &str, exons: &[(u64, u64)]| -> Option<f64> {
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
        };
        let node_idx: BTreeMap<&GeneKey, usize> = g.genes.iter().enumerate().map(|(k, gk)| (gk, k)).collect();
        for (i, c) in clusters.iter().enumerate() {
            let fid = format!("MCL{i}");
            let mut idx = 0usize;
            let mut hulls: BTreeMap<String, (u64, u64)> = BTreeMap::new();
            let mut pending: Vec<PendingUnit> = Vec::new();
            for (mi, m) in c.members.iter().enumerate() {
                // locus = core hull (trimmed) / member span; dropped members are not units
                // the member's SEDEF core hull (0-based half-open), for `copy_assign --psv-genomic` (§6ep)
                let hull_col = match core_records.get(i).and_then(|v| v.get(mi)).and_then(|r| r.hull) {
                    Some((a, b)) => format!("{}-{}", a.saturating_sub(1), b),
                    None => "NA".to_string(),
                };
                let (lo, hi) = match core_records.get(i).and_then(|v| v.get(mi)) {
                    Some(r) if r.status == CoreStatus::Dropped => {
                        unit_stats.2 += 1;
                        continue;
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
                if n_in_chain == 0 {
                    unit_stats.2 += 1;
                    continue;
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
                    seq.reverse();
                    for b in seq.iter_mut() {
                        *b = match *b {
                            b'A' => b'T',
                            b'T' => b'A',
                            b'C' => b'G',
                            b'G' => b'C',
                            b'a' => b't',
                            b't' => b'a',
                            b'c' => b'g',
                            b'g' => b'c',
                            x => x,
                        };
                    }
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
                pending.push(PendingUnit {
                    member: m.clone(),
                    exons,
                    strand,
                    source: source.to_string(),
                    n_reads: n_reads.max(1),
                    seq,
                    hull_col,
                    sd_depth,
                    core_bp,
                    nearest_col,
                    rep_col,
                });
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
                let exlen = |u: &PendingUnit| u.exons.iter().map(|(s, e)| e - s).sum::<u64>();
                let mut rep_of_root: BTreeMap<usize, usize> = BTreeMap::new();
                for k in 0..n {
                    let r = find(&mut parent, k);
                    let e = rep_of_root.entry(r).or_insert(k);
                    if exlen(&pending[k]) > exlen(&pending[*e]) {
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
                    "{fid}\t{idx}\tMCL_{}_{us}\t{}\t{us}\t{ue}\t{}\t{strand}\t{}\t{}\t{source}\t{}\t{}\t{}\t{}\t{}",
                    m.0,
                    m.0,
                    exons.len(),
                    u.n_reads,
                    exons.iter().map(|(s, e)| format!("{s}-{e}")).collect::<Vec<_>>().join(","),
                    u.hull_col, u.sd_depth, u.core_bp, u.nearest_col, u.rep_col
                )?;
                writeln!(uf, ">{fid}|{idx}|{}:{us}-{ue}|{strand}|nexon={}", m.0, exons.len())?;
                uf.write_all(&u.seq)?;
                writeln!(uf)?;
                let h = hulls.entry(m.0.clone()).or_insert((us, ue));
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
        ("merge_overlapping_loci".to_string(), args.merge_overlapping_loci.to_string()),
        ("locus_attribute_edges".to_string(), args.locus_attribute_edges.to_string()),
        ("annotations_folded_into_loci".to_string(), loci.as_ref().map_or(0, |m| m.n_merged()).to_string()),
        ("paf_records_same_locus_skipped".to_string(), g.same_locus_records.to_string()),
        ("core_refine".to_string(), args.core_refine.to_string()),
        ("sedef".to_string(), args.sedef.clone().unwrap_or_else(|| "<unset>".into())),
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
        ("n_nodes".to_string(), g.n_nodes().to_string()),
        ("n_edges".to_string(), g.n_edges().to_string()),
        ("n_clusters".to_string(), clusters.len().to_string()),
    ] {
        writeln!(ph, "{k}\t{v}")?;
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
}
