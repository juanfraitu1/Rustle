//! ASJ verification — deterministic confound control for the `asj` calls (mode-aware).
//!
//! Adds per call: `canonical` (the junction's splice dinucleotides — GT-AG/CT-AC means a REAL splice site,
//! not a mapping-induced gap), `frac_mq0` (MAPQ-0 fraction at the anchor/donor), `dist` (anchor→junction bp),
//! `transversion` (genetic vs RNA-edit), and `high_confidence`. The frac_mq0 reading INVERTS by mode:
//!   - het ASJ (single/multisnp): high frac_mq0 = collapsed-paralog masquerade -> low is good. HC = transversion
//!     AND uniquely-mapped AND canonical.
//!   - copy-specific junctions (copy): multimapping is EXPECTED; the artifact-ruler is a canonical junction
//!     (a real splice site, per-molecule). HC = canonical.
//! Reads `<prefix>.asj.tsv`, writes `<prefix>.asj_verified.tsv`.

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::{HashMap, HashSet};
use std::io::Write;

use rustle::genome::GenomeIndex;
use rustle::vg_family::allele_specific_junctions::{anchor_junction_dist, is_transversion};
use rustle::vg_family::denovo_assemble::frac_mq0_at;

#[derive(Parser, Debug)]
#[command(name = "asj_verify", about = "Confound control for ASJ / copy-specific-junction calls")]
struct Args {
    /// Coordinate-sorted BAM with a `.bai` (for frac_mq0).
    #[arg(long)]
    bam: String,
    /// Genome FASTA with a `.fai` (for canonical splice dinucleotides).
    #[arg(long)]
    fasta: String,
    /// The `<prefix>.asj.tsv` produced by `asj`.
    #[arg(long)]
    calls: String,
    /// frac_mq0 threshold below which a het anchor counts as uniquely-mapped.
    #[arg(long, default_value_t = 0.10)]
    mq0: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let text = std::fs::read_to_string(&args.calls).with_context(|| format!("reading {}", args.calls))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().context("empty calls file")?.split('\t').collect();
    let idx = |name: &str| header.iter().position(|h| *h == name);
    let single = header.contains(&"anchor_pos"); // single-anchor mode has the anchor columns
    let rows: Vec<Vec<String>> = lines.map(|l| l.split('\t').map(String::from).collect()).collect();

    let chroms: HashSet<String> = rows.iter().filter_map(|r| r.first().cloned()).collect();
    let genome = GenomeIndex::from_fasta_contigs(&args.fasta, &chroms).context("loading FASTA contigs")?;
    let get = |r: &[String], name: &str| idx(name).and_then(|i| r.get(i)).cloned().unwrap_or_default();

    let mut mq0cache: HashMap<(String, u64), f64> = HashMap::new();
    let out_path = args.calls.replace(".asj.tsv", ".asj_verified.tsv");
    let mut fh = std::fs::File::create(&out_path)?;
    writeln!(fh, "{}\tcanonical\tfrac_mq0\tdist\ttransversion\thigh_confidence", header.join("\t"))?;

    let mut n_hc = 0usize;
    for r in &rows {
        let chrom = r.first().cloned().unwrap_or_default();
        let donor: u64 = get(r, "donor").parse().unwrap_or(0);
        let acceptor: u64 = get(r, "acceptor").parse().unwrap_or(0);
        let canonical = if genome.is_canonical_junction(&chrom, donor, acceptor, '+') {
            "+"
        } else if genome.is_canonical_junction(&chrom, donor, acceptor, '-') {
            "-"
        } else {
            "no"
        };
        // frac_mq0 at the anchor (single) or the junction donor (multi/copy)
        let pos = if single { get(r, "anchor_pos").parse().unwrap_or(donor) } else { donor };
        let fm = *mq0cache
            .entry((chrom.clone(), pos))
            .or_insert_with(|| frac_mq0_at(&args.bam, &chrom, pos).unwrap_or(0.0));
        let (tv, dist) = if single {
            let ax = get(r, "ax").bytes().next().unwrap_or(b'N');
            let ay = get(r, "ay").bytes().next().unwrap_or(b'N');
            let anc: u64 = get(r, "anchor_pos").parse().unwrap_or(0);
            (is_transversion(ax, ay), anchor_junction_dist(anc, donor, acceptor))
        } else {
            (false, 0)
        };
        let hc = if single {
            tv && fm < args.mq0 && canonical != "no"
        } else {
            canonical != "no" // copy/multi: a real (canonical) splice site is the artifact-ruler
        };
        if hc {
            n_hc += 1;
        }
        writeln!(fh, "{}\t{}\t{:.3}\t{}\t{}\t{}", r.join("\t"), canonical, fm, dist, tv as u8, hc as u8)?;
    }
    eprintln!(
        "[asj_verify] mode={} | {} calls -> {} high-confidence ({}) -> {}",
        if single { "het" } else { "copy/multi" },
        rows.len(),
        n_hc,
        if single { "transversion + uniquely-mapped + canonical" } else { "canonical splice site" },
        out_path
    );
    Ok(())
}
