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
    detect_single_copy_baseline_genome_wide,
    detect_conflict_catalog_genome_wide, detect_conflict_catalog_genome_wide_xchrom,
    detect_homology_catalog_genome_wide, family_protein_coheres, homology_refine_params,
    refine_families_exon_sum, DenovoConfig, RefineParams,
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
    /// DISABLE homology refinement and emit the RAW read-conflict catalog. Refinement is ON BY DEFAULT: the raw
    /// conflict graph admits repeat-bridge (unrelated genes sharing an intronic repeat) and large-gene
    /// mis-chain false positives (~9% of families genome-wide, `bench/GW_CATALOG_FP_AUDIT.md`); the homology +
    /// distinct-locus gate removes them. Use `--no-refine` for the raw catalog (e.g. reproducing the pre-fix
    /// numbers, or when minimap2 is unavailable). The `--refine` flag is now a no-op kept for back-compat.
    #[arg(long, default_value_t = false)]
    no_refine: bool,
    /// With `--refine`, align the GENOMIC span (introns INCLUDED) instead of the exon-sum. Captures intron
    /// divergence (separates copies that are exon-identical but intron-divergent), at the cost of being
    /// stricter on older paralogs whose introns have diverged. Default off (exon-sum).
    #[arg(long, default_value_t = false)]
    refine_introns: bool,
    /// Add the PROTEIN divergent tier (longest-ORF -> mmseqs, fident>=0.50): recovers synonymous-divergent
    /// CODING paralogs (RABL2B-class) that nucleotide seeds cannot reach. Batched into one mmseqs run.
    /// Additive only. Needs `mmseqs`. Default off (the sensitive nucleotide tier is already on). Applies in
    /// TWO places: with `--refine`, as an extra within-family membership edge; with `--homology-primary`,
    /// as a THIRD genome-wide E_r DEFINITION-edge tier (alongside asm20/sensitive) that groups coding
    /// paralogs past the ~0.65 nucleotide-identity floor where reads map all-primary.
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
    /// Enumerate genomic copy number (famCN) per family via Liftoff-style genome projection (spec §7):
    /// project each family's best-supported consensus onto the genome (minimap2) to recover K=0 collapses
    /// the RNA read-conflict/homology path merges into one locus. Writes `<out>.famcn.tsv`. Defaults ON
    /// when `--min-identity 0.98` (Soto SD98 mode), since that mode targets exactly this near-identical
    /// segdup regime.
    #[arg(long, default_value_t = false)]
    enumerate_copies: bool,

    /// Emit the single-copy baseline instead of the family catalog: `<out>.single_copy.tsv` (one row per
    /// single-copy chi(H)=1 locus) + `<out>.lambda_global.tsv` (the genome-wide median n_reads = lambda_global,
    /// the copy-number normalizer). A TABLE, not a GTF -- this is copy-number calibration, not an isoform catalog.
    #[arg(long, default_value_t = false)]
    single_copy_baseline: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut cfg = DenovoConfig::default();
    cfg.complete_poa_core = args.complete_core;

    if args.single_copy_baseline {
        use rustle::vg_family::single_copy::lambda_global;
        let loci = detect_single_copy_baseline_genome_wide(
            &args.bam, &args.fasta, args.threads, args.win, args.min_copies, !args.no_refine, &cfg,
        )?;
        let mut sc = std::fs::File::create(format!("{}.single_copy.tsv", args.out))?;
        writeln!(sc, "chrom\tstart\tend\tstrand\tn_reads\tn_exons\tchi_h\tn_psv")?;
        for l in &loci {
            writeln!(sc, "{}\t{}\t{}\t{}\t{}\t{}\t1\t0", l.chrom, l.start, l.end, l.strand, l.n_reads, l.n_exons)?;
        }
        let lam = lambda_global(&loci);
        let lam_str = lam.map(|x| format!("{x}")).unwrap_or_else(|| "NA".into());
        let mut lf = std::fs::File::create(format!("{}.lambda_global.tsv", args.out))?;
        writeln!(lf, "lambda_global\tn_single_copy_loci")?;
        writeln!(lf, "{}\t{}", lam_str, loci.len())?;
        eprintln!(
            "[gw-catalog] single-copy baseline: {} loci -> lambda_global={} ({}.single_copy.tsv, {}.lambda_global.tsv)",
            loci.len(), lam_str, args.out, args.out
        );
        return Ok(());
    }
    // unify to `Vec<Vec<DenovoTranscript>>` (each = a family's copies) for a single emit path. The
    // cross-chrom path emits same-chromosome families (incl. inverted duplications and distant paralogs)
    // and cross-chromosome families together; the default path is the tight same-strand tandem-array view.
    let mut refine_params = homology_refine_params(args.min_identity, args.threads);
    // `--protein-tail` also promotes protein homology to a genome-wide E_r DEFINITION edge in the
    // homology-primary catalog (in addition to feeding the `--refine` block below): it recovers coding
    // paralogs that have diverged past the nucleotide seeds' ~0.65 identity floor.
    refine_params.protein_tail = args.protein_tail;
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
    // Exon-sum (FLNC) homology + distinct-locus refinement (the principled membership criterion). ON BY DEFAULT
    // (removes repeat-bridge + large-gene mis-chain FPs); `--no-refine` opts out to the raw conflict catalog.
    let refine = !args.no_refine;
    let fams: Vec<Vec<DenovoTranscript>> = if refine {
        let n_raw = raw.len();
        let params = RefineParams {
            threads: args.threads,
            include_introns: args.refine_introns,
            // Always supply the genome: the additive genomic-span tier (recovers near-identical segdups whose
            // partial transcript models fail the exon-sum coverage floor) needs it even without --refine-introns.
            intron_fasta: Some(args.fasta.clone()),
            nucleotide_sensitive: !args.no_sensitive,
            protein_tail: args.protein_tail,
            ..Default::default()
        };
        let refined = refine_families_exon_sum(raw, &params, None)?;
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
    writeln!(fh, "family_id\tn_copies\tn_chroms\tchroms\tcross_chrom\tavg_reads\tprotein_coheres")?;
    writeln!(ch, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads")?;

    use std::collections::BTreeSet;
    // deterministic order: by the first copy's (chrom, start).
    let mut fams = fams;
    fams.sort_by(|a, b| {
        let ka = a.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        let kb = b.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        ka.cmp(&kb)
    });
    // spec §6: orthogonal protein-coherence QC flag — does each family's ORFs form a connected
    // protein-homology graph? `None` (emitted as "NA") when `--protein-tail` is off (no mmseqs call).
    let qc_params = &refine_params;
    let coheres = family_protein_coheres(&fams, qc_params);
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
        let coheres_str = match coheres[fi] {
            Some(true) => "true",
            Some(false) => "false",
            None => "NA",
        };
        writeln!(
            fh,
            "{fid}\t{}\t{}\t{}\t{}\t{:.1}\t{coheres_str}",
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

    // famCN / totalCN via genome projection (spec §7): a family's RNA-observed copies are a LOWER bound
    // on its true genomic copy number when copies collapse onto one locus (K=0). Project each family's
    // best-supported consensus (most-read-supported copy's exon-sum) back onto the genome in ONE batched
    // minimap2 pass (one genome-index load for ALL families, not one per family -- avoids re-indexing the
    // whole genome hundreds/thousands of times) and count additional disjoint loci beyond the already-
    // known copy loci, bucketed by identity: totalCN at the >=0.80 asm20 floor, famCN at the >=0.98 Soto
    // SD98 near-identical floor.
    // totalCN cov>=0.80 (batch floor): full-length divergent copies without partial-domain fragment
    // inflation; famCN keeps cov>=0.90 (near-identical, Soto SD98). A coverage sweep on real families
    // showed cov>=0.50 inflates totalCN with partial/domain-fragment hits (GWFAM18: 16 vs 12=truth at
    // cov>=0.70), while cov>=0.90 is too strict and drops divergent full-length copies (GSTM 4 vs 19
    // truth); cov>=0.80 is the sweet spot (verified: GSTM 20 vs 19, RABL2 6 vs 5, GWFAM18 12 vs 11,
    // GWFAM21 22 vs 22).
    let enumerate = (args.enumerate_copies || args.min_identity == Some(0.98)) && args.homology_primary;
    if enumerate {
        use std::collections::HashMap;
        let consensuses: Vec<(String, Vec<u8>)> = fams
            .iter()
            .enumerate()
            .map(|(fi, copies)| {
                let cons = copies.iter().max_by_key(|c| c.n_reads).map(|c| c.seq.clone()).unwrap_or_default();
                (format!("GWFAM{fi}"), cons)
            })
            .collect();
        let known: HashMap<String, Vec<(String, u64, u64)>> = fams
            .iter()
            .enumerate()
            .map(|(fi, copies)| {
                let loci: Vec<(String, u64, u64)> = copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect();
                (format!("GWFAM{fi}"), loci)
            })
            .collect();
        let proj_by_fam = match rustle::vg_family::genome_projection::project_families_batch(
            &consensuses, &args.fasta, &known, 0.80, 0.80, &refine_params.minimap2, args.threads,
        ) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("[gw-catalog] batch famCN projection failed ({e}); famCN/totalCN fall back to n_rna_copies");
                HashMap::new()
            }
        };
        let mut ff = std::fs::File::create(format!("{}.famcn.tsv", args.out))?;
        writeln!(ff, "family_id\tn_rna_copies\tfamCN\ttotalCN\tprojection_loci")?;
        for (fi, copies) in fams.iter().enumerate() {
            let fid = format!("GWFAM{fi}");
            let n_rna = copies.len();
            let proj = proj_by_fam.get(&fid).cloned().unwrap_or_default();
            // totalCN: divergent copies admitted (id>=0.80; cov>=0.80 already enforced by the batch call
            // above -- full-length divergent copies without partial/domain-fragment inflation).
            let n_total_loci = proj.iter().filter(|p| p.identity >= 0.80).count();
            // famCN: Soto SD98 near-identical FULL-LENGTH copies only (id>=0.98 AND cov>=0.90) -- excludes
            // fragment/partial hits that would otherwise inflate the Soto-comparable metric. cov>=0.90 is
            // a subset of the returned cov>=0.80 loci.
            let n_fam_loci = proj.iter().filter(|p| p.identity >= 0.98 && p.cov >= 0.90).count();
            let total_cn = n_rna + n_total_loci;
            let fam_cn = n_rna + n_fam_loci;
            // projection_loci lists all totalCN-contributing loci (id>=0.80, cov>=0.80).
            let loci = proj
                .iter()
                .filter(|p| p.identity >= 0.80)
                .map(|p| format!("{}:{}-{}@{:.4}", p.chrom, p.start, p.end, p.identity))
                .collect::<Vec<_>>()
                .join(";");
            writeln!(ff, "{fid}\t{n_rna}\t{fam_cn}\t{total_cn}\t{loci}")?;
        }
        eprintln!("[gw-catalog] famCN/totalCN batch projection -> {}.famcn.tsv", args.out);
    }
    Ok(())
}
