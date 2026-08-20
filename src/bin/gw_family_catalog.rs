//! gw_family_catalog — the GENOME-WIDE multi-copy-gene-family catalog (interest I / O1).
//!
//! DEFINITION: a family is a HOMOLOGY component — a γ-quasi-clique of the transcribed-homology graph E_r
//! spanning ≥2 physically-distinct loci — enforced by `refine_families_exon_sum` (the shared E_r primary
//! tier, `-k 11 -w 5` @ 0.60 by default since X.4,
//! cov-of-shorter≥0.50 + sensitive tier + ≥2-distinct-loci gate; default ON). The read-conflict graph is
//! NOT the family definition: it SEEDS the catalog (the loci among which reads are genuinely confused = where
//! copies must be co-resolved) and is the within-family copy-number oracle (χ_H). Homology decides
//! membership; conflict counts copies. (`--homology-primary` makes E_r the membership proposer directly,
//! in-engine via minimap2 — the cleanest form; the default seeds with the read-conflict graph, then the
//! homology gate enforces the definition.) This replaces the shipped `denovo_families.tsv` (built by an
//! arbitrary `core_recip≥0.13` similarity threshold that over-merges: DNFAM0 = 728 members chr1→chrY).
//! Writes `<out>.families.tsv` + `<out>.copies.tsv`.

use anyhow::Result;
use clap::Parser;
use std::io::Write;

use rustle::vg_family::denovo_pipeline::{
    detect_single_copy_baseline_genome_wide,
    detect_conflict_catalog_genome_wide, detect_conflict_catalog_genome_wide_xchrom,
    detect_homology_catalog_genome_wide, families_from_reps_certified, family_protein_coheres,
    homology_refine_params, refine_families_exon_sum, write_joint_rna_dna_certificate, DenovoConfig,
    certificates_for_families, FamilyCertificate, RefineParams, Substrate};
use rustle::vg_family::family_detect::DenovoTranscript;

#[derive(Parser, Debug)]
#[command(about = "Genome-wide de-tie read-conflict multi-copy-family catalog (no similarity threshold).")]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast per-contig region query). Required
    /// unless `--from-genome` is given.
    #[arg(long)]
    bam: Option<String>,
    /// GENOME-ONLY mode: discover SD families from the genome alone (no reads, no annotation). Argument is
    /// a BED of search windows (e.g. `bench/soto/80_fams.chr.bed`). Reps are duplicated genomic loci found by
    /// self-alignment; grouping uses the SAME `homology_blocks` core as `--homology-primary`. Mutually
    /// exclusive with `--bam`.
    #[arg(long)]
    from_genome: Option<String>,
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
    /// DEPRECATED — build the O1 catalog from the READ-CONFLICT graph (E_c) instead of sequence homology.
    ///
    /// The name is misleading and is kept only to reproduce historical catalogs. This flag never was a
    /// "chromosome setting": the DEFAULT homology mode has no chromosome restriction either, and reports
    /// `cross_chrom` per family as an OBSERVATION. What this flag actually selects is the conflict-graph
    /// catalog, where a read confused between two copies forms the family edge.
    ///
    /// Two reasons not to use it. It loses on every measurement: 65.5% against homology's 81.8% overall,
    /// and 14% against 57% on DISPERSED paralogs — the very case the name promises — where it fails
    /// silently. And it makes E_c participate in O1, so "a family is E_r, and O1 is independent of O2" is
    /// false of this mode alone. E_c remains correct and necessary for O2 (copy assignment); that path is
    /// untouched by this flag.
    #[arg(long, default_value_t = false)]
    cross_chrom: bool,
    /// REFINE the raw conflict-graph catalog by exon-sum (FLNC) homology (requires minimap2 on PATH). A real
    /// multi-copy family must (i) have its copies MUTUALLY HOMOLOGOUS full-length (their spliced exon-sum
    /// sequences all-vs-all align, asm20 identity>=0.80 cov-of-shorter>=0.50) and (ii) span >=2 spatially-
    /// DISTINCT loci (disjoint genomic spans). This removes the Alu-SINE-bridged cross-chrom over-merges and
    /// the same-locus isoform/fragment contamination that read-conflict alone admits. Emits one row per
    /// distinct locus. ON BY DEFAULT for the legacy conflict catalogs (`--window-catalog`/`--cross-chrom`).
    ///
    /// On the DEFAULT (homology, E_r) catalog this flag is OPT-IN and changes the emitted OBJECT: the
    /// definition is γ-quasi-clique(E_r) (`docs/seeded_family_definition.md` §1), and refine appends a
    /// SECOND clustering stage — its own primary + genomic-span edge tiers (⚠ the tier itself is now the
    /// SAME object as E_r's since X.4; what still differs is the clustering), no core-coverage
    /// denominator, no stub guard, clustered by CONNECTED COMPONENTS. Use it only to reproduce the
    /// pre-2026-08-09 default catalog. At shipped defaults it is a measured no-op (0/50 files change on
    /// the 25-region panel); under `--homology-genomic-span` it is not (7f/28c on vs 6f/29c off).
    ///
    /// ⚠ **REJECTED on the O1 homology catalog since 2026-08-20** — that path has ONE object and refine
    /// is not it. Meaningful only on `--window-catalog` / `--cross-chrom`, where it is already the
    /// default (so this flag is a no-op there too, and `--no-refine` is the lever).
    #[arg(long, default_value_t = false)]
    refine: bool,
    /// DISABLE homology refinement and emit the RAW read-conflict catalog. Refinement is ON BY DEFAULT: the raw
    /// conflict graph admits repeat-bridge (unrelated genes sharing an intronic repeat) and large-gene
    /// mis-chain false positives (~9% of families genome-wide, `bench/GW_CATALOG_FP_AUDIT.md`); the homology +
    /// distinct-locus gate removes them. Use `--no-refine` for the raw catalog (e.g. reproducing the pre-fix
    /// numbers, or when minimap2 is unavailable). It wins over `--refine` and applies to every catalog.
    ///
    /// ⚠ On the DEFAULT homology catalog refine is already OFF since 2026-08-09, so `--no-refine` is a
    /// no-op there; it is the conflict catalogs (`--window-catalog`/`--cross-chrom`) that it still turns off.
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
    /// [DEFAULT, accepted for compatibility] HOMOLOGY-PRIMARY (E_r) family definition: build families from
    /// nucleotide homology between locus representatives (gamma-quasi-clique), not the read-conflict graph.
    /// Captures ancient paralogs + the K=0 collapses the conflict path misses.
    ///
    /// ALREADY CROSS-CHROMOSOME. Chromosome never enters the grouping — membership is decided by sequence
    /// alone — so dispersed families (MAGEA, TSPY, RFPL, GSTM) are found without any extra flag. The
    /// `cross_chrom` column in families.tsv reports whether a family turned out to span chromosomes; it is
    /// an output, not an input. This is now the default, so passing this flag is a no-op.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,
    /// Build the O1 catalog from the read-conflict graph restricted to a `--win`-bounded window on ONE
    /// chromosome (the historical pre-2026-07 default). Superseded by the homology mode; kept addressable
    /// so the old behaviour can still be reproduced.
    #[arg(long, default_value_t = false)]
    window_catalog: bool,
    /// E_r nucleotide identity floor (sensitive tier). Default from RefineParams (~0.60). `0.98` = Soto SD98 mode.
    #[arg(long)]
    min_identity: Option<f64>,
    /// With `--homology-primary`: compute the E_r family edge on each locus's GENOMIC SPAN instead of its
    /// assembled exon-sum. Fixes family FRAGMENTATION — copies of one family whose transcripts assemble
    /// DISJOINT exon subsets share almost no assembled sequence, so no edge forms and the family splits
    /// (measured: 75% of should-link pairs had no exon-sum alignment at all). RNA still decides WHICH loci
    /// are expressed; the reference supplies complete sequence for grouping. Needs no read depth/copy number.
    /// Default off (exon-sum) — existing catalogs are byte-identical unless this is set.
    #[arg(long, default_value_t = false)]
    homology_genomic_span: bool,
    /// RNA-only precision guard: accept a nucleotide E_r edge only when at least one passing alignment is
    /// forward (`+`) between transcript-oriented representatives. Reverse-only matches are consistent with
    /// inverted repeats rather than homologous transcripts. Incompatible with `--from-genome`, whose DNA
    /// reps are deliberately stored in reference orientation and may be true inverted duplications.
    ///
    /// DEFAULT ON since 2026-08-19 — this flag is retained only so an explicit request is still accepted
    /// (and so it can still conflict-check against `--from-genome`). Use `--no-rna-forward-only` to opt out.
    #[arg(long, default_value_t = false)]
    rna_forward_only: bool,
    /// Opt OUT of the forward-only E_r guard and accept reverse-only nucleotide matches, restoring the
    /// pre-2026-08-19 behaviour. Reverse-only edges are enriched ~46× for overlapping ANTISENSE pairs
    /// (a gene and its own antisense partner), so expect lower precision.
    #[arg(long, default_value_t = false)]
    no_rna_forward_only: bool,
    /// Compare the RNA exon-sum and DNA genomic-span homology graphs on the SAME emitted RNA loci. Writes
    /// `<out>.joint_edges.tsv`, `.joint_families.tsv`, and `.joint_rule.tsv` with RNA+DNA corroborated,
    /// RNA-only, and DNA-only edges. Reporting-only: it never silently unions/intersects family membership.
    #[arg(long, default_value_t = false)]
    joint_dna_rna: bool,
    /// Enumerate genomic copy number (famCN) per family via Liftoff-style genome projection (spec §7):
    /// project each family's best-supported consensus onto the genome (minimap2) to recover K=0 collapses
    /// the RNA read-conflict/homology path merges into one locus. Writes `<out>.famcn.tsv`. Defaults ON
    /// when `--min-identity 0.98` (Soto SD98 mode), since that mode targets exactly this near-identical
    /// segdup regime.
    #[arg(long, default_value_t = false)]
    enumerate_copies: bool,
    /// Re-admit near-identical families that collapse to <2 RNA loci as K0_COLLAPSED, reporting
    /// genome-projected copy NUMBER (needs DNA parCN for per-read resolution). Writes <out>.collapsed.tsv.
    /// Default off; when off, all existing output is byte-identical.
    #[arg(long, default_value_t = false)]
    collapse_enumerate: bool,
    /// Re-admit EXON-IDENTICAL (0-PSV) but heavily-EXPRESSED families that collapse to <2 RNA loci as a
    /// K0_COLLAPSED_EXPRESSED copy-number class (projection >=2 loci, each read-supported; needs DNA parCN for
    /// per-read resolution). Writes <out>.expressed_collapsed.tsv. Requires --homology-primary (the drop
    /// point lives on that path; a silent no-op without it). Default off; byte-identical when off.
    #[arg(long, default_value_t = false)]
    collapse_expressed: bool,

    /// Emit the single-copy baseline instead of the family catalog: `<out>.single_copy.tsv` (one row per
    /// single-copy chi(H)=1 locus) + `<out>.lambda_global.tsv` (the genome-wide median n_reads = lambda_global,
    /// the copy-number normalizer). A TABLE, not a GTF -- this is copy-number calibration, not an isoform catalog.
    #[arg(long, default_value_t = false)]
    single_copy_baseline: bool,

    /// Generalized projection: project EVERY resolved copy consensus (not one per family) at id>=0.98 +
    /// >=3-read support to localize members a single-consensus projection misses. Writes <out>.allproj.tsv
    /// (DNA-localized parCN, NOT added to copies.tsv). Requires --homology-primary. Default off (byte-identical).
    #[arg(long, default_value_t = false)]
    project_all_families: bool,

    /// DNA-family fallback (DNA edge oracle): for RNA-orphan expressed loci (0 homology edges -> no RNA
    /// family), project the transcript onto the genome at a LOOSER identity floor (--dna-family-min-identity,
    /// default 0.90 vs collapse-expressed's 0.98) to recover DIVERGENT genomic paralogs whose expressed
    /// sequence is non-homologous (the DNA-family != RNA-family case). Writes <out>.dna_family.tsv (copy
    /// NUMBER only; per-read resolution needs DNA parCN). Requires --homology-primary. Default off; byte-identical when off.
    #[arg(long, default_value_t = false)]
    dna_family_fallback: bool,
    /// Projection identity floor for --dna-family-fallback (default 0.90).
    #[arg(long)]
    dna_family_min_identity: Option<f64>,

    /// QUERY the emitted catalog: print the family (the connected component / γ-quasi-clique block)
    /// containing `chrom:start-end`. Repeatable. Writes `<out>.seed.tsv` in addition to stdout.
    ///
    /// ⚠ **This is a projection, not a second pipeline.** The seed is NOT an input to the definition:
    /// it does not filter the BAM, it does not restrict the window set, and it does not enter node
    /// construction or `E_r`. The catalog is built exactly as it would be without `--seed`, and the
    /// flag then looks up which emitted block the seed lands in (max overlapping bp). Every output
    /// file other than `<out>.seed.tsv` is byte-identical with and without the flag — that identity
    /// is what makes the answer a property of the partition rather than of the query.
    ///
    /// Coordinates are 1-based inclusive (GFF/samtools convention). Seeds that overlap no emitted
    /// copy ABSTAIN; they are not snapped to a nearest locus. See
    /// `src/rustle/vg_family/seed_projection.rs` for why a seed-DERIVED node set (the
    /// `bench/crossspecies` probe) is not seed-invariant and is not what this queries.
    #[arg(long)]
    seed: Vec<String>,
}

/// `--seed`: project each seed onto the EMITTED catalog and report the block containing it.
///
/// The view built here must match `emit_catalog`'s rows exactly — families in emitted order (so the
/// index is `GWFAM{i}`) and copies sorted by `(chrom, start)` (so the index is `copy_idx`) —
/// otherwise the ids printed would not name rows the user can look up in `copies.tsv`.
fn project_seeds(
    specs: &[String],
    fams: &[Vec<DenovoTranscript>],
    out: &str,
    min_copies: usize,
) -> Result<()> {
    use rustle::vg_family::seed_projection::{
        format_seed_rows, parse_seed, project_seed, SeedHit, SEED_TSV_HEADER,
    };
    if specs.is_empty() {
        return Ok(());
    }
    let view: Vec<Vec<(String, u64, u64)>> = fams
        .iter()
        .map(|copies| {
            let mut v: Vec<(String, u64, u64)> =
                copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect();
            v.sort_by(|a, b| (a.0.as_str(), a.1).cmp(&(b.0.as_str(), b.1)));
            v
        })
        .collect();
    // Parse EVERY seed before writing anything, so a typo in the third seed does not leave a
    // half-written .seed.tsv beside a good catalog.
    let seeds = specs.iter().map(|s| parse_seed(s)).collect::<Result<Vec<_>>>()?;
    let mut fh = std::fs::File::create(format!("{out}.seed.tsv"))?;
    writeln!(fh, "{SEED_TSV_HEADER}")?;
    println!("{SEED_TSV_HEADER}");
    let (mut n_hit, mut n_abstain) = (0usize, 0usize);
    for seed in &seeds {
        let hit = project_seed(seed, &view);
        match hit {
            SeedHit::Hit { .. } => n_hit += 1,
            SeedHit::Abstain => n_abstain += 1,
        }
        for row in format_seed_rows(seed, &view, &hit) {
            writeln!(fh, "{row}")?;
            println!("{row}");
        }
    }
    eprintln!(
        "[gw-catalog] --seed: {n_hit} projected onto emitted families, {n_abstain} abstained -> {out}.seed.tsv"
    );
    if n_abstain > 0 {
        eprintln!(
            "[gw-catalog]   ABSTAIN_NO_OVERLAP = no EMITTED copy overlaps the seed. At --min-copies \
             {min_copies} a block below that size is not written, so this does NOT distinguish \"the \
             seed's locus is a singleton\" from \"there is no transcribed locus there\"; re-run with \
             --min-copies 1 to separate the two."
        );
    }
    Ok(())
}

/// Parse a BED of search windows (`chrom\tstart\tend`, extra columns ignored) for `--from-genome`.
fn read_windows_bed(path: &str) -> Result<Vec<(String, u64, u64)>> {
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 3 {
            out.push((f[0].to_string(), f[1].parse()?, f[2].parse()?));
        }
    }
    Ok(out)
}

/// Write the family catalog (`<out>.families.tsv` + `.copies.tsv` + `.copies.fa`). Shared by the RNA
/// catalog path and the DNA `--from-genome` path: both hand it the same `Vec<Vec<DenovoTranscript>>`
/// (each inner vec = a family's copies). Sorts families by their first copy's (chrom, start) and returns
/// the sorted vec so the RNA path's downstream sections (collapse/famcn/project-all) see the same order.
/// Genomic exon blocks of a copy as `start-end,start-end,...` (half-open, genomic order, ascending).
///
/// `copies.tsv` carried only `n_exon`, so every downstream consumer that needed the actual exons had to
/// re-derive them from an annotation — which defeats the point, since the loci are annotation-free. It also
/// silently forced read-depth work onto the transcript SPAN: measuring famCN over `start..end` averages in
/// every intron, and on a 208 kb NBPF model that measures intronic repeat content rather than the gene
/// (measured: within-family famCN spread 60 over spans vs 4.4 over exons).
///
/// An unspliced copy (no introns) yields the single block `start-end`.
fn exon_blocks(c: &DenovoTranscript) -> String {
    let mut out = String::new();
    let mut prev = c.start;
    for &(d, a) in &c.introns {
        // Defensive: a malformed chain (donor before the running cursor, or acceptor <= donor) would
        // otherwise emit a reversed block that a consumer would read as a valid interval.
        if d > prev && a >= d {
            if !out.is_empty() {
                out.push(',');
            }
            out.push_str(&format!("{prev}-{d}"));
            prev = a;
        }
    }
    if c.end > prev {
        if !out.is_empty() {
            out.push(',');
        }
        out.push_str(&format!("{prev}-{}", c.end));
    }
    out
}

fn emit_catalog(
    out: &str,
    fams: Vec<Vec<DenovoTranscript>>,
    certs: Option<Vec<FamilyCertificate>>,
    refine_params: &RefineParams,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    use std::collections::BTreeSet;
    let mut fh = std::fs::File::create(format!("{out}.families.tsv"))?;
    let mut ch = std::fs::File::create(format!("{out}.copies.tsv"))?;
    // The exon-sum (spliced, FLNC-derived) sequence of every copy, in transcription orientation (for the
    // RNA path); for the DNA path it is the genomic locus sequence. Annotation-free family-validation substrate.
    let mut sh = std::fs::File::create(format!("{out}.copies.fa"))?;
    // `n_edges`/`density`/`lambda`/`cut_certified` are APPENDED last, so existing header-keyed readers
    // keep working unchanged. They are a REPORT on the emitted family, never an input to it:
    //   lambda        = minimum number of E_r edges whose removal splits this family (0 = n<2)
    //   cut_certified = lambda >= 2, i.e. "no single alignment record's loss can split this family"
    // A 2-copy family necessarily has lambda = 1, so `cut_certified=false` is NOT a defect flag — see
    // docs/seeded_family_definition.md §1★.5.
    writeln!(
        fh,
        "family_id\tn_copies\tn_chroms\tchroms\tcross_chrom\tavg_reads\tprotein_coheres\
         \tn_edges\tdensity\tlambda\tcut_certified"
    )?;
    // `exons` is APPENDED last so existing header-keyed readers keep working unchanged.
    writeln!(ch, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons")?;

    // The certificate is positional, so it MUST travel with its family through the sort below — sorting
    // the two lists independently would silently attach each λ to the wrong family.
    if let Some(c) = &certs {
        anyhow::ensure!(
            c.len() == fams.len(),
            "certificate/family count mismatch: {} certificates for {} families",
            c.len(),
            fams.len()
        );
    }
    let mut paired: Vec<(Vec<DenovoTranscript>, Option<FamilyCertificate>)> = match certs {
        Some(c) => fams.into_iter().zip(c.into_iter().map(Some)).collect(),
        None => fams.into_iter().map(|f| (f, None)).collect(),
    };
    // deterministic order: by the first copy's (chrom, start).
    paired.sort_by(|a, b| {
        let ka = a.0.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        let kb = b.0.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        ka.cmp(&kb)
    });
    let cert_of: Vec<Option<FamilyCertificate>> = paired.iter().map(|p| p.1).collect();
    let fams: Vec<Vec<DenovoTranscript>> = paired.into_iter().map(|p| p.0).collect();
    // orthogonal protein-coherence QC flag; "NA" when --protein-tail is off (no mmseqs call).
    let coheres = family_protein_coheres(&fams, refine_params);
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
        // Structural certificate; "NA" on paths that do not carry one (never a silent 0 — a missing
        // certificate and a genuinely disconnected family must not print the same thing).
        let (n_edges_s, density_s, lambda_s, certified_s) = match cert_of[fi] {
            Some(c) => (
                c.n_edges.to_string(),
                format!("{:.4}", c.density),
                c.lambda.to_string(),
                (c.lambda >= 2).to_string(),
            ),
            None => ("NA".into(), "NA".into(), "NA".into(), "NA".into()),
        };
        writeln!(
            fh,
            "{fid}\t{}\t{}\t{}\t{}\t{:.1}\t{coheres_str}\t{n_edges_s}\t{density_s}\t{lambda_s}\t{certified_s}",
            copies.len(),
            chroms.len(),
            chroms.iter().cloned().collect::<Vec<_>>().join(","),
            cross,
            avg_reads
        )?;
        let mut sorted = copies.clone();
        sorted.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
        for (ci, c) in sorted.iter().enumerate() {
            writeln!(
                ch,
                "{fid}\t{ci}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.tid,
                c.chrom,
                c.start,
                c.end,
                c.introns.len() + 1,
                c.strand,
                c.n_reads,
                exon_blocks(c)
            )?;
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
        "[gw-catalog] wrote {} families ({} cross-chromosome) -> {out}.families.tsv + {out}.copies.tsv",
        fams.len(),
        n_xchrom,
    );
    Ok(fams)
}

/// Should the second, exon-sum `refine_families_exon_sum` stage run?
///
/// ⚠ **This is the D1 fix: the shipped catalog used to be TWO definitions.** `refine` ran
/// unconditionally on the homology (E_r) path while the DNA `--from-genome` path returned before it, so
/// RNA emitted `refine(γ-QC(E_r))` and DNA emitted `γ-QC(E_r)`. Refine is not a pass-through — it
/// recomputes edges on its own core and genomic-span SUBSTRATES (no core-coverage denominator,
/// no stub guard) and clusters by CONNECTED COMPONENTS
/// (`denovo_pipeline::homology_components`) instead of the γ-quasi-clique rule the definition names.
/// None of that appears in `docs/seeded_family_definition.md` §1.
///
/// Resolution: refine is **opt-in on the homology (default) catalog** and stays **on by default for the
/// legacy conflict catalogs** (`--window-catalog` / `--cross-chrom`), which is what it was written for
/// and where it was measured to bite (APOBEC3 1f/2c → 0f/0c, SHARP 2f/4c → 1f/2c).
///
/// Cost at shipped defaults: **zero** — 0 of 50 output files changed over the 25-region gorilla panel
/// (strict single-copy controls 0/7, family panel 6 families / 24 copies, SHARP 1f/2c), and one whole
/// per-family minimap2 stage disappears. It stops being free the moment E_r's edge substrate stops
/// matching refine's: under `--homology-genomic-span` the same panel reads 7f/28c with refine vs 6f/29c
/// without, refine having split two recovered MAGEA copies into a spurious 2-copy family and dropped
/// GSTM4 from the GSTM family.
///
/// ⭐ **2026-08-20: refine can no longer be opted INTO on the O1 path — there is ONE default path.**
/// It was opt-in from 2026-08-09 so the pre-fix catalog stayed reproducible. That option is what let the
/// shipped 494-family catalog be built with `refine(γ-QC(E_r))` while the default moved to `γ-QC(E_r)`,
/// and for the discrepancy to go unnoticed for six weeks (`docs/o1_catalog_provenance.md`). A flag that
/// changes what the catalog *is* is not a convenience, it is a provenance hazard: the emitted object
/// must not depend on an invocation detail. Passing `--refine` on the O1 path is now an ERROR that says
/// so, rather than silently producing a non-definitional object.
///
/// Refine keeps its real home: the legacy conflict catalogs, where it was written to run and where it
/// was measured to help.
///
/// `--no-refine` still wins over everything (the documented escape hatch when minimap2 is unavailable).
fn refine_enabled(o1_homology: bool, no_refine_flag: bool) -> bool {
    if no_refine_flag {
        return false;
    }
    // The O1 homology catalog emits γ-quasi-clique(E_r) and nothing else. Refine is for the legacy
    // conflict catalogs only — there is no argument that turns it on here.
    !o1_homology
}

fn main() -> Result<()> {
    let args = Args::parse();
    // O1 mode. Homology (E_r) is THE mode and the default; the two conflict-graph catalogs are legacy and
    // must now be asked for by name. Derived once here rather than read off `args.homology_primary`, which
    // is a no-op compatibility flag and would be true even when a legacy catalog was requested.
    let o1_homology = !args.cross_chrom && !args.window_catalog;
    if args.cross_chrom && args.window_catalog {
        anyhow::bail!("--cross-chrom and --window-catalog select different legacy catalogs; pass at most one");
    }
    // ONE DEFAULT PATH (2026-08-20). Refine re-clusters by CONNECTED COMPONENTS over its own substrates,
    // which is not the γ-quasi-clique(E_r) object the definition names, so it cannot be requested on the
    // O1 catalog at all. Erroring beats ignoring: a silently-ignored flag is how the shipped catalog's
    // provenance was lost in the first place.
    if o1_homology && (args.refine || args.refine_introns) {
        anyhow::bail!(
            "--refine/--refine-introns cannot be used with the O1 homology catalog: refine clusters by \
             connected components over its own substrates, which is NOT the γ-quasi-clique(E_r) object \
             docs/seeded_family_definition.md §1 defines. The O1 catalog has ONE path and emits ONE \
             object. Refine remains available on the legacy conflict catalogs (--window-catalog / \
             --cross-chrom), where it is on by default. See docs/o1_catalog_provenance.md for why the \
             opt-in was removed."
        );
    }
    if args.joint_dna_rna && !o1_homology {
        anyhow::bail!(
            "--joint-dna-rna is a certificate for the default RNA homology catalog; it is not defined for legacy read-conflict catalogs"
        );
    }
    if args.joint_dna_rna && args.single_copy_baseline {
        anyhow::bail!("--joint-dna-rna requires a family catalog and is incompatible with --single-copy-baseline");
    }
    if args.joint_dna_rna
        && matches!(std::env::var("RUSTLE_SHARED_EXON"), Ok(v) if v != "0" && !v.is_empty())
    {
        anyhow::bail!(
            "--joint-dna-rna currently requires the standard nucleotide E_r tiers; RUSTLE_SHARED_EXON is a different edge definition"
        );
    }
    if args.cross_chrom {
        eprintln!(
            "[gw-catalog] WARNING: --cross-chrom builds O1 from the READ-CONFLICT graph (E_c), not sequence\n\
             [gw-catalog]          homology. It is retained only to reproduce historical catalogs.\n\
             [gw-catalog]          * it is NOT what makes families cross-chromosome — the default homology\n\
             [gw-catalog]            mode has no chromosome restriction either, and finds dispersed families\n\
             [gw-catalog]            BETTER (57% vs 14% on dispersed paralogs; 81.8% vs 65.5% overall);\n\
             [gw-catalog]          * it makes E_c participate in O1, so O1 is no longer independent of O2.\n\
             [gw-catalog]          Drop the flag to use the homology (E_r) catalog."
        );
    }
    if args.window_catalog {
        eprintln!(
            "[gw-catalog] WARNING: --window-catalog is the pre-2026-07 conflict catalog, bounded to a \
             --win window on one chromosome. Superseded by the default homology mode."
        );
    }
    let project_all = (args.project_all_families
        || std::env::var("RUSTLE_PROJECT_ALL_FAMILIES").ok().as_deref() == Some("1"))
        && o1_homology;
    let mut cfg = DenovoConfig::from_env();
    cfg.complete_poa_core = args.complete_core;
    cfg.collapse_enumerate = args.collapse_enumerate || cfg.collapse_enumerate;
    cfg.collapse_expressed = args.collapse_expressed || cfg.collapse_expressed;
    cfg.dna_family_fallback = args.dna_family_fallback || cfg.dna_family_fallback;
    if let Some(x) = args.dna_family_min_identity { cfg.dna_family_min_identity = x; }

    // --from-genome: read-free/annotation-free DNA family catalog. Discovers duplicated genomic loci by
    // self-alignment, then groups them with the SAME homology_blocks core the RNA --homology-primary path
    // uses (via families_from_reps). Returns before any BAM-consuming path.
    if let Some(win_bed) = args.from_genome.as_deref() {
        use rustle::vg_family::from_genome::{genome_reps, GenomeRepParams};
        if args.bam.is_some() {
            anyhow::bail!("--from-genome and --bam are mutually exclusive (genome-only mode takes no reads)");
        }
        // The guard is default-ON, but it is meaningless on reference-oriented DNA reps, so
        // --from-genome silently runs without it. Only an EXPLICIT request is an error.
        if args.rna_forward_only {
            anyhow::bail!(
                "--rna-forward-only is defined only for transcript-oriented RNA representatives; \
                 it cannot be applied to reference-oriented --from-genome DNA loci"
            );
        }
        if args.joint_dna_rna {
            anyhow::bail!(
                "--joint-dna-rna must be run on the RNA catalog with --bam; --from-genome has no RNA arm"
            );
        }
        let mut refine_params = homology_refine_params(args.min_identity, args.threads);
        // ⚠ These reps are stored in REFERENCE orientation, where a '-' record can be a REAL inverted
        // segmental duplication. DECLARING the substrate makes every orientation-sensitive option inert
        // here by construction — see `RefineParams::forward_only_active`. Flipping a default can no
        // longer reach this path (that exact mistake was caught by
        // `genome_mode_grouping_keeps_an_inverted_duplication` on 2026-08-19).
        refine_params.substrate = Substrate::ReferenceOriented;
        refine_params.require_forward_alignment = false;
        refine_params.protein_tail = args.protein_tail;
        // DNA-mode GROUPING knobs (anti over/under-merge levers). Env-only and applied ONLY on this branch,
        // so the RNA path stays byte-identical.
        //   RUSTLE_GENOME_MIN_COVERAGE — homology-edge coverage floor (fraction of the SHORTER locus the
        //     alignment must cover). Default 0.50. This is the documented repeat/duplicon-bridge defense:
        //     raising it rejects "two families share one duplicated block" edges that fuse distinct families.
        //   RUSTLE_GENOME_GAMMA — γ-quasi-clique density for family blocks. Default 0.20. Raising it demands
        //     denser within-family connectivity, splitting blocks held together by a few bridge edges.
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_COVERAGE") {
            if let Ok(x) = v.parse::<f64>() { refine_params.min_coverage = x; }
        }
        let gamma: f64 = std::env::var("RUSTLE_GENOME_GAMMA")
            .ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.20);
        eprintln!(
            "[gw-catalog-genome] grouping: min_identity={:.2} min_coverage={:.2} gamma={:.2}",
            refine_params.min_identity, refine_params.min_coverage, gamma
        );
        let windows = read_windows_bed(win_bed)?;
        let mut gp = GenomeRepParams::from_env();
        gp.threads = args.threads;
        gp.minimap2 = refine_params.minimap2.clone();
        // `--min-identity` sets BOTH floors: the SD locus-discovery floor here AND the family-grouping floor
        // (via refine_params above). Env RUSTLE_GENOME_MIN_IDENTITY (already applied by from_env) is overridden
        // only when the flag is given, so `--from-genome --min-identity 0.98` = SD98 discovery + grouping.
        if let Some(mi) = args.min_identity { gp.min_identity = mi; }
        let reps = genome_reps(&args.fasta, &windows, &gp)?;
        eprintln!("[gw-catalog-genome] {} windows -> {} duplicated-locus reps", windows.len(), reps.len());
        let (dna_fams, dna_certs) = families_from_reps_certified(
            reps, &refine_params, gamma, args.min_copies, 0,
        )?;
        let dna_fams = emit_catalog(&args.out, dna_fams, Some(dna_certs), &refine_params)?;
        project_seeds(&args.seed, &dna_fams, &args.out, args.min_copies)?;
        return Ok(());
    }

    // BAM is required for every path below (single-copy baseline, conflict/homology catalogs, project-all).
    let bam = args
        .bam
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("--bam is required unless --from-genome is given"))?;

    if args.single_copy_baseline {
        use rustle::vg_family::single_copy::lambda_global;
        let loci = detect_single_copy_baseline_genome_wide(
            bam, &args.fasta, args.threads, args.win, args.min_copies, !args.no_refine, &cfg,
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
    // Same knob as the --from-genome path, honoured here too: the E_r coverage floor (fraction of the
    // SHORTER exon-sum the alignment must cover, default 0.50) is the duplicon-bridge defence, and it is the
    // criterion that decides borderline members. Exposed so a member's exclusion can be attributed to
    // coverage rather than inferred -- e.g. NPIPB12, whose 13 homology edges all reach >= 0.89 identity but
    // top out at 0.25 coverage (bench/soto/merge_quality_analysis.md §24).
    if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_COVERAGE") {
        if let Ok(x) = v.parse::<f64>() { refine_params.min_coverage = x; }
    }
    // `--protein-tail` also promotes protein homology to a genome-wide E_r DEFINITION edge in the
    // homology-primary catalog (in addition to feeding the `--refine` block below): it recovers coding
    // paralogs that have diverged past the nucleotide seeds' ~0.65 identity floor.
    refine_params.protein_tail = args.protein_tail;
    // DEFAULT ON (2026-08-19): forward-only unless explicitly opted out. `--rna-forward-only`
    // remains accepted so an explicit request still conflict-checks against `--from-genome`.
    // The substrate declaration is what ACTIVATES it: on a reference-oriented substrate the flag is
    // inert by construction (`RefineParams::forward_only_active`).
    refine_params.substrate = Substrate::TranscriptOriented;
    refine_params.require_forward_alignment = !args.no_rna_forward_only;
    // `--homology-genomic-span`: E_r edges on the genomic span of each RNA-detected locus (anti-fragmentation).
    // Needs the genome path; harmless when the flag is off (field stays false, edges stay exon-sum).
    refine_params.homology_genomic_span = args.homology_genomic_span;
    if args.homology_genomic_span {
        refine_params.intron_fasta = Some(args.fasta.clone());
        if !o1_homology {
            eprintln!("[gw-catalog] note: --homology-genomic-span affects the E_r edge, which drives family membership only in the (default) homology catalog — it is inert for the legacy conflict catalogs");
        }
    }
    let (raw, raw_certs, collapsed, expressed, dna_families): (
        Vec<Vec<DenovoTranscript>>,
        Vec<FamilyCertificate>,
        Vec<rustle::vg_family::collapse_enumerate::CollapsedFamily>,
        Vec<rustle::vg_family::collapse_enumerate::ExpressedCollapsedFamily>,
        Vec<rustle::vg_family::collapse_enumerate::ExpressedCollapsedFamily>,
    ) = if o1_homology {
        detect_homology_catalog_genome_wide(
            bam,
            &args.fasta,
            args.threads,
            args.min_copies,
            &cfg,
            &refine_params,
            // γ-quasi-clique density floor. Same knob the --from-genome path reads, so the two E_r
            // callers cannot drift apart. Default 0.20; 0.40 is the high-precision setting, which splits
            // low-density blocks instead of admitting them as one family.
            std::env::var("RUSTLE_GENOME_GAMMA").ok().and_then(|v| v.parse().ok()).unwrap_or(0.20),
        )?
    } else if args.cross_chrom {
        // CONFLICT catalogs build no E_r graph, so there is no λ to certify: an EMPTY certificate vector
        // makes the column print "NA" rather than a fabricated 0.
        (
            detect_conflict_catalog_genome_wide_xchrom(
                bam, &args.fasta, args.threads, args.min_copies, &cfg,
            )?,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    } else {
        let catalog = detect_conflict_catalog_genome_wide(
            bam, &args.fasta, args.threads, args.win, args.min_copies, &cfg,
        )?;
        (
            catalog.into_iter().map(|c| c.copies).collect(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    };
    // Exon-sum (FLNC) homology + distinct-locus refinement. ON BY DEFAULT for the legacy CONFLICT catalogs
    // (--window-catalog / --cross-chrom), where the raw read-conflict graph admits repeat-bridge and
    // large-gene mis-chain FPs that this stage removes. OPT-IN (`--refine`) on the default HOMOLOGY catalog:
    // there the definition is already γ-quasi-clique(E_r), and refine would append a second, undocumented
    // clustering stage over a different edge set. See `refine_enabled` for the measurements.
    let refine = refine_enabled(o1_homology, args.no_refine);
    if o1_homology && refine {
        eprintln!(
            "[gw-catalog] NOTE: --refine on the homology catalog appends a SECOND clustering stage \
             (core + genomic-span substrates, connected components) on top of γ-quasi-clique(E_r). \
             The emitted object is then no longer the definition in docs/seeded_family_definition.md §1."
        );
    }
    // The λ certificate describes the γ-quasi-clique(E_r) families. `--refine` re-clusters them over a
    // DIFFERENT edge set, so the pre-refine certificates would no longer describe the emitted rows and
    // are never carried across — a wrong λ is worse than an absent one. Since 2026-08-19 they are
    // RECOMPUTED over the refined rows instead of dropped, so the column is a number rather than "NA"
    // and it describes what is actually written.
    let mut certs: Option<Vec<FamilyCertificate>> = if raw_certs.is_empty() { None } else { Some(raw_certs) };
    let fams: Vec<Vec<DenovoTranscript>> = if refine {
        certs = None;
        let n_raw = raw.len();
        let params = RefineParams {
            threads: args.threads,
            include_introns: args.refine_introns,
            // Always supply the genome: the additive genomic-span tier (recovers near-identical segdups whose
            // partial transcript models fail the exon-sum coverage floor) needs it even without --refine-introns.
            intron_fasta: Some(args.fasta.clone()),
            nucleotide_sensitive: !args.no_sensitive,
            protein_tail: args.protein_tail,
            // Inherit ALL THREE E_r thresholds from `refine_params` rather than taking struct defaults --
            // otherwise `--min-identity` and RUSTLE_GENOME_MIN_COVERAGE are silently ignored on this path,
            // which is the code path `--refine` actually runs.
            //
            // ⚠ `sensitive_identity` was missing here. `homology_refine_params` deliberately sets BOTH
            // floors from `--min-identity` (see `homology_refine_params_min_identity_sets_both_tiers`), but
            // this struct then rebuilt the sensitive floor from `Default`, so `--min-identity 0.93` gave
            // run 1 a floor of 0.93 while run 3 silently kept 0.60. The user asked for one floor and got
            // two, and the looser one is the effective floor of the union.
            min_identity: refine_params.min_identity,
            min_coverage: refine_params.min_coverage,
            sensitive_identity: refine_params.sensitive_identity,
            require_forward_alignment: refine_params.require_forward_alignment,
            substrate: refine_params.substrate,
            ..Default::default()
        };
        let refined = refine_families_exon_sum(raw, &params, None, cfg.conflict.min_reads)?;
        // Rebuild the certificate over exactly the rows that will be emitted. On failure the column
        // falls back to "NA" — an absent certificate is honest, a stale one is not.
        certs = match certificates_for_families(&refined, &params) {
            Ok(c) => Some(c),
            Err(e) => {
                eprintln!("[gw-catalog] refine: λ certificate not recomputed ({e}); column stays NA");
                None
            }
        };
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

    // Sort + write the family catalog (families.tsv/copies.tsv/copies.fa). Shared by the RNA path
    // and the DNA --from-genome path; returns the sorted families for the downstream sections.
    let fams = emit_catalog(&args.out, fams, certs, &refine_params)?;
    if args.joint_dna_rna {
        write_joint_rna_dna_certificate(&args.out, &fams, &args.fasta, &refine_params)?;
    }
    // `--seed` is a QUERY over what was just emitted: it runs AFTER the catalog exists and cannot
    // influence it. No-op (and no `.seed.tsv`) when the flag is absent.
    project_seeds(&args.seed, &fams, &args.out, args.min_copies)?;

    // K=0-collapsed families re-admitted by Task 3 (`--collapse-enumerate`): copy-NUMBER only, kept fully
    // separate from families.tsv/copies.tsv so the OFF path (flag off or nothing collapsed) writes no file
    // and byte-identical-OFF holds.
    if cfg.collapse_enumerate && !collapsed.is_empty() {
        let mut cf = std::fs::File::create(format!("{}.collapsed.tsv", args.out))?;
        writeln!(cf, "family_id\tchrom\tstart\tend\tfamCN\tn_alt_reads\talt_frac\tstatus\tprojection_loci")?;
        for (i, fam) in collapsed.iter().enumerate() {
            writeln!(cf, "{}", rustle::vg_family::collapse_enumerate::format_collapsed_row(&format!("GWFAMc{i}"), fam))?;
        }
        eprintln!("[gw-catalog] wrote {} K=0-collapsed families -> {}.collapsed.tsv", collapsed.len(), args.out);
    }

    // K0_COLLAPSED_EXPRESSED families re-admitted by Task 3 (`--collapse-expressed`): same isolation as
    // `--collapse-enumerate` above -- OFF path (flag off or nothing expressed-collapsed) writes no file.
    if cfg.collapse_expressed && !expressed.is_empty() {
        let mut ef = std::fs::File::create(format!("{}.expressed_collapsed.tsv", args.out))?;
        writeln!(ef, "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci")?;
        for (i, fam) in expressed.iter().enumerate() {
            writeln!(ef, "{}", rustle::vg_family::collapse_enumerate::format_expressed_collapsed_row(&format!("GWFAMe{i}"), fam))?;
        }
        eprintln!("[gw-catalog] collapse-expressed: {} K0_COLLAPSED_EXPRESSED families -> {}.expressed_collapsed.tsv", expressed.len(), args.out);
    }

    // DNA-family fallback (`--dna-family-fallback`): RNA-orphan loci recovered by the DNA edge oracle. Same
    // isolation contract -- OFF path (flag off or nothing recovered) writes no file (byte-identical).
    if cfg.dna_family_fallback && !dna_families.is_empty() {
        let mut df = std::fs::File::create(format!("{}.dna_family.tsv", args.out))?;
        writeln!(df, "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci")?;
        for (i, fam) in dna_families.iter().enumerate() {
            writeln!(df, "{}", rustle::vg_family::collapse_enumerate::format_dna_family_row(&format!("GWFAMdna{i}"), fam))?;
        }
        eprintln!("[gw-catalog] dna-family-fallback: {} DNA_FAMILY_RNA_NONHOMOLOGOUS loci -> {}.dna_family.tsv", dna_families.len(), args.out);
    }

    // famCN / totalCN via genome projection (spec §7) -- the VG copy-number leg: LAND the family variation
    // graph's consensus PATH on the genome and count near-identical landing sites. A family's RNA-observed
    // copies are a LOWER bound on famCN when copies collapse onto one locus (K=0). Project each family's
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
    let enumerate = (args.enumerate_copies || args.min_identity == Some(0.98)) && o1_homology;
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

    // Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not just
    // the best-supported one per family) to localize members a single-consensus projection misses. OPTIONAL,
    // additive, DNA-localized parCN leg -- dedup'd (reciprocal-overlap >=0.50) + read-support gated (>=3
    // primary reads over the locus), written to its own file so the RNA-split catalog (families.tsv/
    // copies.tsv) and the famCN/totalCN batch projection (famcn.tsv) are untouched.
    if project_all {
        use rustle::vg_family::project_all::{CopyIn, all_copy_consensuses, known_from_fams, dedup_overlapping, overlaps_any, format_allproj_row};
        // Build (fid, copies) with the SAME fid the catalog uses.
        let fam_copies: Vec<(String, Vec<CopyIn>)> = fams.iter().enumerate().map(|(fi, copies)| {
            (format!("GWFAM{fi}"), copies.iter().map(|c| CopyIn { seq: c.seq.clone(), chrom: c.chrom.clone(), start: c.start, end: c.end }).collect())
        }).collect();
        let consensuses = all_copy_consensuses(&fam_copies);
        let known = known_from_fams(&fam_copies);
        let proj = rustle::vg_family::genome_projection::project_families_batch(
            &consensuses, &args.fasta, &known, 0.98, 0.90, &refine_params.minimap2, args.threads,
        ).unwrap_or_default();
        let spans_by_fam: std::collections::HashMap<&String, &Vec<(String,u64,u64)>> = known.iter().collect();
        let mut rows: Vec<String> = Vec::new();
        for (fid, locs) in &proj {
            for l in dedup_overlapping(locs.clone()) {
                // read-support gate: >=3 primary reads over the locus
                let n_support = rustle::vg_family::denovo_assemble::reads_in_region(bam, &l.chrom, l.start, l.end, args.threads)
                    .map(|(p, _)| p.len()).unwrap_or(0);
                if n_support < 3 { continue; }
                let overlaps = spans_by_fam.get(fid).map(|s| overlaps_any(&l.chrom, l.start, l.end, s)).unwrap_or(false);
                rows.push(format_allproj_row(fid, &l, n_support, overlaps));
            }
        }
        rows.sort(); // `proj` is a HashMap (random iteration order); sort so allproj.tsv is reproducible run-to-run
        if !rows.is_empty() {
            let mut af = std::fs::File::create(format!("{}.allproj.tsv", args.out))?;
            writeln!(af, "family_id\tchrom\tstart\tend\tidentity\tn_support_reads\toverlaps_existing_copy")?;
            for r in &rows { writeln!(af, "{r}")?; }
            eprintln!("[gw-catalog] project-all-families: {} loci (id>=0.98, >=3 reads) -> {}.allproj.tsv", rows.len(), args.out);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(start: u64, end: u64, introns: Vec<(u64, u64)>) -> DenovoTranscript {
        DenovoTranscript {
            tid: "t".into(),
            chrom: "chr1".into(),
            start,
            end,
            n_reads: 5,
            strand: '+',
            introns,
            seq: Vec::new(),
            distinguishing_uniq: 0,
            core_bp: 0,
            stub: false,
            tes: None,
        }
    }

    #[test]
    fn exon_blocks_splits_on_the_intron_chain() {
        // two introns -> three exons; blocks are half-open and in ascending genomic order
        let t = tx(1000, 2000, vec![(1100, 1300), (1500, 1800)]);
        assert_eq!(exon_blocks(&t), "1000-1100,1300-1500,1800-2000");
        // the block count must agree with the n_exon column emitted beside it
        assert_eq!(exon_blocks(&t).split(',').count(), t.introns.len() + 1);
    }

    #[test]
    fn exon_blocks_unspliced_copy_is_one_block() {
        assert_eq!(exon_blocks(&tx(500, 900, vec![])), "500-900");
    }

    #[test]
    fn exon_blocks_rejects_a_malformed_chain_rather_than_emitting_a_reversed_block() {
        // acceptor before donor, and a donor behind the running cursor: both must be skipped, never
        // written as "1500-1200" which a consumer would read as a valid interval.
        let t = tx(1000, 2000, vec![(1500, 1200), (900, 950)]);
        let blocks = exon_blocks(&t);
        assert_eq!(blocks, "1000-2000", "malformed introns must not produce reversed blocks");
        for b in blocks.split(',') {
            let (s, e) = b.split_once('-').unwrap();
            assert!(s.parse::<u64>().unwrap() < e.parse::<u64>().unwrap(), "reversed block {b}");
        }
    }

    // ---- D1: the shipped catalog must be ONE definition, not two -------------------------------
    //
    // `refine_families_exon_sum` used to run unconditionally on the homology (E_r) path, which meant
    // the default catalog emitted `refine(gamma-QC(E_r))` while the DNA `--from-genome` path emitted
    // `gamma-QC(E_r)`. Refine is not a pass-through: it recomputes edges on its own core and
    // genomic-span@0.80 tiers, WITHOUT the core-coverage denominator or the stub guard, and clusters by
    // CONNECTED COMPONENTS (denovo_pipeline::homology_components) rather than by the gamma-quasi-clique
    // rule the definition names (docs/seeded_family_definition.md §1). None of those tiers appear in the
    // definition. So on the homology path refine is now OPT-IN; on the legacy conflict catalogs
    // (--window-catalog / --cross-chrom), where it was written and where it was MEASURED to bite
    // (APOBEC3 1f/2c -> 0f/0c, SHARP 2f/4c -> 1f/2c), it stays on by default.
    //
    // Measured at shipped defaults when this flipped: 0 of 50 output files changed across the 25-region
    // gorilla panel (strict single-copy controls 0/7, family panel 6 families / 24 copies, SHARP 1f/2c).
    // The flip is NOT free once E_r's substrate moves: under --homology-genomic-span the same panel goes
    // 7f/28c (refine on) -> 6f/29c (refine off), refine having split two recovered MAGEA copies into a
    // spurious family and dropped GSTM4.

    #[test]
    fn the_o1_catalog_has_exactly_one_path() {
        // ⭐ 2026-08-20. The O1 homology catalog emits γ-QC(E_r) and NOTHING else. There is no argument
        // that turns refine on here — the opt-in was removed because it is what let the shipped
        // 494-family catalog be built with a different object than the default produced, unnoticed for
        // six weeks (docs/o1_catalog_provenance.md). `--refine` on this path is rejected in main().
        assert!(!refine_enabled(true, false), "default O1 run must be γ-QC(E_r) alone");
        assert!(!refine_enabled(true, true), "--no-refine changes nothing: refine was never on here");
    }

    #[test]
    fn refine_stays_on_by_default_on_the_conflict_catalogs() {
        // --window-catalog / --cross-chrom: refine is the documented FP filter for the raw conflict
        // graph (repeat-bridge + large-gene mis-chain), and is what it was written for. Removing the
        // O1 opt-in must not disturb its real home.
        assert!(refine_enabled(false, false));
    }

    #[test]
    fn no_refine_still_wins_on_the_conflict_catalogs() {
        // --no-refine remains the documented escape hatch (e.g. minimap2 absent).
        assert!(!refine_enabled(false, true));
    }

    #[test]
    fn exon_blocks_sum_of_lengths_is_the_spliced_length() {
        let t = tx(0, 1000, vec![(200, 400), (600, 700)]);
        let total: u64 = exon_blocks(&t)
            .split(',')
            .map(|b| {
                let (s, e) = b.split_once('-').unwrap();
                e.parse::<u64>().unwrap() - s.parse::<u64>().unwrap()
            })
            .sum();
        assert_eq!(total, 1000 - (400 - 200) - (700 - 600));
    }
}
