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
    build_clusters, graph_from_paf, mcl, Cluster, GeneKey, GraphParams,
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

    /// Smallest cluster reported. 3 is the floor at which density carries signal.
    #[arg(long, default_value_t = 3)]
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

    /// Optional RNA BAM. Without it `corroborated` is reported as `NA` — ⚠ which is NOT 0.000, the
    /// repeat-clique signature. The two must never be conflated.
    #[arg(long)]
    bam: Option<String>,

    /// A member is corroborated when it carries at least this many primary reads with an ALIGNED BLOCK
    /// inside it (`-F 2308`; a read spliced OVER a locus is not evidence for it).
    #[arg(long, default_value_t = 3)]
    min_reads: usize,

    #[arg(long)]
    out: String,
}

/// Exon-union length per gene, keyed by the gene's own GFF 1-based `(contig, start, end)`.
///
/// ⚠ Reads `exon` features and joins them to their gene by the `gene=` attribute (a Name), which is why
/// the Name -> ID map is built first. Reports the join rate: **a no-op result is the signature of a failed
/// join** (§6dd), so a silent fallback to the span must never be possible.
fn exonic_blocks(gff: &str) -> Result<BTreeMap<GeneKey, Vec<(u64, u64)>>> {
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
    let args = Args::parse();
    let p = GraphParams {
        min_identity: args.min_identity,
        min_cov_longer: args.min_cov_longer,
        min_bp: args.min_bp,
        exonic_overlap: args.exonic_overlap,
        reject_overlapping: args.reject_overlapping,
        min_exonic_bp: args.min_exonic_bp,
    };

    let blocks = match &args.gff {
        Some(g) => exonic_blocks(g)?,
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
    let g = graph_from_paf(&paf, &exonic, &blocks, &p);
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
