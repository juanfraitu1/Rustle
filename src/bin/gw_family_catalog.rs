//! gw_family_catalog — the GENOME-WIDE de-tie read-conflict multi-copy-family catalog (interest I / O1).
//!
//! Replaces the shipped `denovo_families.tsv` (built by the arbitrary `core_recip≥0.13` similarity
//! threshold over coordinate windows, which over-merges: DNFAM0 = 728 members chr1→chrY) with the
//! PRINCIPLED, threshold-free family definition run AT SCALE: a family is a connected component of loci
//! among which reads are genuinely confused (the de-tie read-conflict graph), with the same-strand +
//! disjoint-loci fixes applied. Writes `<out>.families.tsv` + `<out>.copies.tsv`.

use anyhow::Result;
use clap::Parser;
use std::io::Write;

use rustle::vg_family::denovo_pipeline::{
    detect_conflict_catalog_genome_wide, detect_conflict_catalog_genome_wide_xchrom,
    detect_homology_catalog_genome_wide, refine_families_exon_sum, DenovoConfig, RefineParams,
};
use rustle::vg_family::family_detect::DenovoTranscript;

#[derive(Parser, Debug)]
#[command(about = "Genome-wide de-tie read-conflict multi-copy-family catalog (no similarity threshold).")]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast per-contig region query).
    #[arg(long)]
    bam: String,
    /// Genome FASTA (with a `.fai`).
    #[arg(long)]
    fasta: String,
    /// Output prefix; writes `<out>.families.tsv` + `<out>.copies.tsv`.
    #[arg(long)]
    out: String,
    /// Minimum copies for a family (real paralog PAIRS are common → default 2).
    #[arg(long, default_value_t = 2)]
    min_copies: usize,
    /// Co-located window (bp): a family's copies must cluster within this span.
    #[arg(long, default_value_t = 5_000_000)]
    win: u64,
    /// BAM-reading threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Capture CROSS-CHROMOSOME paralog families (RABL2A/B, DAZ-class) — a read confused between copies on
    /// different chromosomes forms a conflict edge; families are NOT restricted to one chromosome. Without
    /// this flag, only co-located same-chrom tandem arrays are emitted (the `--win`-bounded default).
    #[arg(long, default_value_t = false)]
    cross_chrom: bool,
    /// REFINE the raw conflict-graph catalog by exon-sum (FLNC) homology (requires minimap2 on PATH). A real
    /// multi-copy family must (i) have its copies MUTUALLY HOMOLOGOUS full-length (their spliced exon-sum
    /// sequences all-vs-all align, asm20 identity>=0.80 cov-of-shorter>=0.50) and (ii) span >=2 spatially-
    /// DISTINCT loci (disjoint genomic spans). This removes the Alu-SINE-bridged cross-chrom over-merges and
    /// the same-locus isoform/fragment contamination that read-conflict alone admits. Emits one row per
    /// distinct locus. Strongly recommended with `--cross-chrom`.
    #[arg(long, default_value_t = false)]
    refine: bool,
    /// With `--refine`, align the GENOMIC span (introns INCLUDED) instead of the exon-sum. Captures intron
    /// divergence (separates copies that are exon-identical but intron-divergent), at the cost of being
    /// stricter on older paralogs whose introns have diverged. Default off (exon-sum).
    #[arg(long, default_value_t = false)]
    refine_introns: bool,
    /// With `--refine`, also add the PROTEIN divergent tier (longest-ORF -> mmseqs, fident>=0.50): recovers
    /// synonymous-divergent CODING paralogs (RABL2B-class) that nucleotide seeds cannot reach. Batched into one
    /// mmseqs run. Additive only. Needs `mmseqs`. Default off (the sensitive nucleotide tier is already on).
    #[arg(long, default_value_t = false)]
    protein_tail: bool,
    /// With `--refine`, DISABLE the default sensitive nucleotide tier (`-k11 -w5`), reverting to the pure
    /// asm20 (id>=0.80) baseline. For reproducing the exact asm20-only catalog.
    #[arg(long, default_value_t = false)]
    no_sensitive: bool,
    /// POA-CORE COMPLETION (cross-chrom only): after the read-conflict families are built, attach
    /// loosely-related paralogs at NEW loci whose contiguous POA core to a family member clears `t_core` —
    /// reaching divergent copies the conflict graph (confusable-only, ~>=87% identity) misses. Bounded +
    /// seeded by the conflict families. Default off.
    #[arg(long, default_value_t = false)]
    complete_core: bool,
    /// HOMOLOGY-PRIMARY (E_r) family definition: build families from exon-sum nucleotide homology
    /// (gamma-quasi-clique), not the read-conflict graph. Captures ancient paralogs + K=0 collapses the
    /// conflict path misses. Ships alongside --cross-chrom.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,
    /// E_r nucleotide identity floor (sensitive tier). Default from RefineParams (~0.60). `0.98` = Soto SD98 mode.
    #[arg(long)]
    min_identity: Option<f64>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut cfg = DenovoConfig::default();
    cfg.complete_poa_core = args.complete_core;
    // unify to `Vec<Vec<DenovoTranscript>>` (each = a family's copies) for a single emit path. The
    // cross-chrom path emits same-chromosome families (incl. inverted duplications and distant paralogs)
    // and cross-chromosome families together; the default path is the tight same-strand tandem-array view.
    let mut refine_params = RefineParams { threads: args.threads, ..Default::default() };
    if let Some(mi) = args.min_identity {
        refine_params.min_identity = mi;
        refine_params.sensitive_identity = mi.min(refine_params.sensitive_identity);
    }
    let raw: Vec<Vec<DenovoTranscript>> = if args.homology_primary {
        detect_homology_catalog_genome_wide(
            &args.bam, &args.fasta, args.threads, args.min_copies, &cfg, &refine_params, 0.20,
        )?
    } else if args.cross_chrom {
        detect_conflict_catalog_genome_wide_xchrom(
            &args.bam, &args.fasta, args.threads, args.min_copies, &cfg,
        )?
    } else {
        let catalog = detect_conflict_catalog_genome_wide(
            &args.bam, &args.fasta, args.threads, args.win, args.min_copies, &cfg,
        )?;
        catalog.into_iter().map(|c| c.copies).collect()
    };
    // Optional exon-sum (FLNC) homology + distinct-locus refinement (the principled membership criterion).
    let fams: Vec<Vec<DenovoTranscript>> = if args.refine {
        let n_raw = raw.len();
        let params = RefineParams {
            threads: args.threads,
            include_introns: args.refine_introns,
            intron_fasta: if args.refine_introns { Some(args.fasta.clone()) } else { None },
            nucleotide_sensitive: !args.no_sensitive,
            protein_tail: args.protein_tail,
            ..Default::default()
        };
        let refined = refine_families_exon_sum(raw, &params)?;
        eprintln!(
            "[gw-catalog] {}{}{} refine: {} raw families -> {} refined (homology component AND >= 2 distinct loci)",
            if args.refine_introns { "genomic(intron-inclusive)" } else { "exon-sum" },
            if args.no_sensitive { "" } else { " +sensitive" },
            if args.protein_tail { " +protein" } else { "" },
            n_raw,
            refined.len()
        );
        refined
    } else {
        raw
    };

    let mut fh = std::fs::File::create(format!("{}.families.tsv", args.out))?;
    let mut ch = std::fs::File::create(format!("{}.copies.tsv", args.out))?;
    // The exon-sum (spliced, FLNC-derived) sequence of every copy, in transcription orientation. This is
    // the substrate for the ANNOTATION-FREE family validation: all-vs-all align a family's copies and a
    // copy is confirmed iff its spliced sequence aligns full-length to a sibling (independent of both the
    // conflict-graph family definition AND of RefSeq gene annotation → de-circularising).
    let mut sh = std::fs::File::create(format!("{}.copies.fa", args.out))?;
    writeln!(fh, "family_id\tn_copies\tn_chroms\tchroms\tcross_chrom\tavg_reads")?;
    writeln!(ch, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads")?;

    use std::collections::BTreeSet;
    // deterministic order: by the first copy's (chrom, start).
    let mut fams = fams;
    fams.sort_by(|a, b| {
        let ka = a.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        let kb = b.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        ka.cmp(&kb)
    });
    let mut n_xchrom = 0usize;
    for (fi, copies) in fams.iter().enumerate() {
        let fid = format!("GWFAM{fi}");
        let chroms: BTreeSet<&str> = copies.iter().map(|c| c.chrom.as_str()).collect();
        let cross = chroms.len() > 1;
        if cross {
            n_xchrom += 1;
        }
        let avg_reads = if copies.is_empty() {
            0.0
        } else {
            copies.iter().map(|c| c.n_reads as f64).sum::<f64>() / copies.len() as f64
        };
        writeln!(
            fh,
            "{fid}\t{}\t{}\t{}\t{}\t{:.1}",
            copies.len(),
            chroms.len(),
            chroms.iter().cloned().collect::<Vec<_>>().join(","),
            cross,
            avg_reads
        )?;
        // copies sorted by (chrom, start) for readability.
        let mut sorted = copies.clone();
        sorted.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
        for (ci, c) in sorted.iter().enumerate() {
            writeln!(
                ch,
                "{fid}\t{ci}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.tid,
                c.chrom,
                c.start,
                c.end,
                c.introns.len() + 1,
                c.strand,
                c.n_reads
            )?;
            // exon-sum FASTA record: header carries family/copy/locus; sequence is the spliced consensus.
            if !c.seq.is_empty() {
                writeln!(
                    sh,
                    ">{fid}|{ci}|{}:{}-{}|{}|nexon={}",
                    c.chrom,
                    c.start,
                    c.end,
                    c.strand,
                    c.introns.len() + 1
                )?;
                sh.write_all(&c.seq)?;
                writeln!(sh)?;
            }
        }
    }
    eprintln!(
        "[gw-catalog] wrote {} families ({} cross-chromosome) -> {}.families.tsv + {}.copies.tsv",
        fams.len(),
        n_xchrom,
        args.out,
        args.out
    );
    Ok(())
}
